import tempfile
import time
import json
from abc import ABC, abstractmethod
from datetime import datetime
import os
import shutil
import subprocess
import gzip
import csv
import glob
import math
import re
from typing import Optional, Dict
import concurrent.futures
import threading
import hashlib
import uuid

from ..task import Task
from ...accession_utils import canonicalize_accession, canonicalize_accessions
from ...database import DBManager
from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE, RAW_PROFILE, resolve_profile_selector
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
)

class OrthoFinderTask(Task):
    '''A class that handles the OrthoFinder task'''
    _PROTEOME_SUFFIXES = (".faa", ".fasta", ".faa.gz", ".fasta.gz")

    @classmethod
    def default_thread_count(cls, registry_required_threads: int, daemon_max_threads: int) -> int:
        return min(max(int(daemon_max_threads or 1), 1), 24)

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=64):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.input_dir = self.data.get("input_dir")
        self.out_dir = self.data.get("out_dir", )
        self.force = self.data.get("force", False)
        # Resolve accessions later once the DB connection is available.
        self.accessions = normalize_accessions(self.data.get("accessions") or [])
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        self.library_id = self.data.get("library_id", None)
        self.library_name = self.data.get("library_name", None)
        self.check_for_previous_run_folders = self.data.get("check_for_previous_run_folders", True)
        self.mcl_inflation = self.data.get("mcl_inflation", None)

    def _clean_up(self, paths):
        fail = False
        for path in paths:
            try:
                if os.path.exists(path):
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
            except OSError as e:
                self.log(f"Error cleaning up {path}: {e}", "ERROR")
                fail = True
        return not fail

    def _clean_up_and_error(self, orthofinder_id, paths, error_message, error_context=None):
        if not self.db_manager.orthofinder.delete_results(orthofinder_id):
            self.error("Error deleting OrthoFinder results after failure.")
        if not self._clean_up(paths):
            self.error("Error cleaning up files after failure.")
        return self.handle_exception(error_message, error_context)

    def _iter_proteome_filenames(self, path):
        for filename in sorted(os.listdir(path)):
            full_path = os.path.join(path, filename)
            if filename.startswith(".") or not os.path.isfile(full_path):
                continue
            if filename.lower().endswith(self._PROTEOME_SUFFIXES):
                yield filename

    def get_accessions_from_proteomes(self, path):
        # From a file path attempt to get the accession from the file name - for creating the run hash
        accessions = []
        for filename in self._iter_proteome_filenames(path):
            stem = filename[:-3] if filename.endswith(".gz") else filename
            stem = os.path.splitext(stem)[0]
            match = re.search(r"(GC[AF]_?\d+\.\d+)", stem, re.IGNORECASE)
            if match:
                accession = canonicalize_accession(match.group(1))
            else:
                accession = canonicalize_accession(stem.split("_", 1)[0])
            if accession:
                accessions.append(accession)
        return accessions

    def _resolve_orthofinder_executable(self, raw_path):
        candidate = str(raw_path or "").strip()
        if not candidate:
            return None
        expanded = os.path.expanduser(candidate)
        if os.path.isdir(expanded):
            binary = os.path.join(expanded, "orthofinder")
            if os.path.isfile(binary) and os.access(binary, os.X_OK):
                return binary
        if os.path.isfile(expanded) and os.access(expanded, os.X_OK):
            return expanded
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
        if os.path.dirname(expanded):
            binary = os.path.join(expanded, "orthofinder")
            if os.path.isfile(binary) and os.access(binary, os.X_OK):
                return binary
        return None

    def _select_proteome_filename(self, location, accession):
        canonical = canonicalize_accession(accession) or str(accession)
        compact = canonical.replace("_", "").lower()
        candidates = []
        for filename in self._iter_proteome_filenames(location):
            lower = filename.lower()
            if "busco_input" in lower:
                continue
            score = 0
            normalized_name = lower.replace("_", "")
            if compact and compact in normalized_name:
                score += 40
            if "protein" in lower:
                score += 20
            if lower.endswith(".faa"):
                score += 8
            elif lower.endswith(".faa.gz"):
                score += 6
            elif lower.endswith(".fasta"):
                score += 4
            else:
                score += 2
            candidates.append((score, filename))
        if not candidates:
            return None
        candidates.sort(key=lambda item: (-item[0], item[1]))
        return candidates[0][1]

    def _build_staged_proteome_filename(self, accession, source_filename, selected_profile=None):
        accession_token = canonicalize_accession(accession) or str(accession)
        original = os.path.basename(str(source_filename))
        stem = original[:-3] if original.endswith(".gz") else original
        suffix = ".faa"
        lowered = stem.lower()
        if lowered.endswith(".fasta"):
            suffix = ".fasta"
            stem = stem[:-6]
        elif lowered.endswith(".faa"):
            suffix = ".faa"
            stem = stem[:-4]
        else:
            stem = os.path.splitext(stem)[0]
        label = str(selected_profile or stem).strip() or stem or "proteome"
        safe_label = re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("._-") or "proteome"
        return f"{accession_token}_{safe_label}{suffix}"

    def _read_logged_species_used_filenames(self, log_path):
        if not log_path or not os.path.exists(log_path):
            return []
        filenames = []
        in_species_block = False
        try:
            with open(log_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not in_species_block:
                        if line.startswith("Species used:"):
                            in_species_block = True
                        continue
                    if not line:
                        if filenames:
                            break
                        continue
                    match = re.match(r"^\d+:\s+(.+)$", line)
                    if match:
                        filenames.append(os.path.basename(match.group(1).strip()))
                        continue
                    if filenames:
                        break
            return sorted(filenames)
        except OSError as exc:
            self.log(f"Error reading OrthoFinder log file {log_path}: {exc}", "WARNING")
            return []

    def _read_logged_command_line(self, log_path):
        if not log_path or not os.path.exists(log_path):
            return None
        try:
            with open(log_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if line.startswith("Command Line:"):
                        return line.split("Command Line:", 1)[1].strip()
        except OSError as exc:
            self.log(f"Error reading OrthoFinder command line from {log_path}: {exc}", "WARNING")
        return None

    @staticmethod
    def _parse_logged_mcl_inflation(command_line):
        text = str(command_line or "").strip()
        if not text:
            return None
        match = re.search(r"(?:^|\s)-I\s+([0-9]*\.?[0-9]+)(?:\s|$)", text)
        if not match:
            return None
        try:
            return float(match.group(1))
        except ValueError:
            return None

    def _has_reusable_results(self, results_folder):
        if not results_folder or not os.path.isdir(results_folder):
            return False
        return os.path.isdir(os.path.join(results_folder, "Orthogroup_Sequences"))

    def _expected_analysis_input_filenames(self, genome_info, location_idx):
        expected = []
        if self.accessions:
            for accession in self.accessions:
                if accession not in genome_info:
                    raise FileNotFoundError(f"Accession {accession} not found in genome information.")
                location = genome_info[accession][location_idx]
                selected_profile = None
                if self.proteome_profile or self.prefer_proteome_profile:
                    selected_profile, _, proteome_path = self._resolve_proteome_input(accession, location)
                    source_filename = os.path.basename(proteome_path)
                else:
                    source_filename = self._select_proteome_filename(location, accession)
                    if not source_filename:
                        raise FileNotFoundError(f"No valid proteome files found for accession {accession}.")
                expected.append(
                    self._build_staged_proteome_filename(
                        accession,
                        source_filename,
                        selected_profile=selected_profile,
                    )
                )
        elif self.input_dir:
            expected.extend(list(self._iter_proteome_filenames(self.input_dir)))
        return sorted(expected)

    def _lookup_proteome_profile_row(self, accession, profile_name):
        requested = self.db_manager.proteomes.resolve_selector_profile_name(str(accession), profile_name) or profile_name
        if requested == RAW_PROFILE:
            raw_id = self.db_manager.proteomes.ensure_raw_profile(str(accession), is_default=False)
            if raw_id is None:
                return requested, None
            return requested, self.db_manager.proteomes.get(int(raw_id))
        return requested, self.db_manager.proteomes.get_profile(str(accession), str(requested))

    def _resolve_requested_proteome_profile(self, accession):
        if self.proteome_profile:
            requested = self.db_manager.proteomes.resolve_selector_profile_name(str(accession), self.proteome_profile)
            if requested:
                return requested
            if self.proteome_profile == DEFAULT_CLEAN_PROFILE:
                default_cleaned = self.db_manager.proteomes.get_default_cleaned_profile_name(str(accession))
                if default_cleaned:
                    return default_cleaned
            return self.proteome_profile
        if self.prefer_proteome_profile:
            preferred, row = self._lookup_proteome_profile_row(str(accession), self.prefer_proteome_profile)
            if row is not None:
                return preferred
        default_profile = self.db_manager.proteomes.get_default_profile_name(str(accession))
        if default_profile:
            return default_profile
        default_cleaned = self.db_manager.proteomes.get_default_cleaned_profile_name(str(accession))
        if default_cleaned:
            return default_cleaned
        return RAW_PROFILE

    def _resolve_proteome_input(self, accession, genome_path):
        requested_profile = self._resolve_requested_proteome_profile(str(accession))
        resolved_profile, row = self._lookup_proteome_profile_row(str(accession), requested_profile)
        if row is None:
            raise FileNotFoundError(
                f"Proteome profile '{requested_profile}' does not exist for accession '{accession}'."
            )
        path = self.db_manager.proteomes.resolve_path(row)
        if not path or not os.path.exists(path):
            raise FileNotFoundError(
                f"Proteome profile '{resolved_profile}' has no readable artifact for accession '{accession}'."
            )
        return resolved_profile, int(row[0]), path

    def get_folder_creation_date(self, folder_path):
        # Get the creation date of the folder
        return datetime.fromtimestamp(os.path.getctime(folder_path)).strftime("%Y-%m-%d %H:%M:%S")

    def generate_orthofinder_id_hash_code(self, accessions):
        return hash(tuple(accessions))
        # # Create a number based on a list of accessions that is reproducable in other instances of the program
        # val = 0
        # for accession in accessions:
        #     for char in accession:
        #         val += ord(char)
        #     # Use a prime number to reduce collisions
        #     val *= 7
        # print(val)
        # return f"orthofinder_{val}"

    def _create_accessions_list_file(self, accessions, filepath):
        try:
            with open(filepath, "w") as f:
                for accession in accessions:
                    f.write(f"{accession}\n")
            return True
        except OSError as e:
            self.log(f"Error creating accessions list file: {e}", "ERROR")
            return False

    def _create_consolidated_tree_file(self, tree_dir, output_name):
        if not os.path.isdir(tree_dir):
            return False
        output_path = os.path.join(tree_dir, output_name)
        if os.path.exists(output_path):
            return True

        tree_files = sorted(
            path for path in glob.glob(os.path.join(tree_dir, "*.txt"))
            if os.path.basename(path) != output_name
        )
        if not tree_files:
            return False

        try:
            with open(output_path, "w") as out_handle:
                written = 0
                for tree_path in tree_files:
                    base = os.path.basename(tree_path)
                    og_name = base[:-9] if base.endswith("_tree.txt") else os.path.splitext(base)[0]
                    with open(tree_path, "r") as tree_handle:
                        tree_newick = tree_handle.read().strip()
                    if not tree_newick:
                        continue
                    out_handle.write(f"{og_name}:{tree_newick}\n")
                    written += 1
            if written == 0:
                os.remove(output_path)
                return False
            self.log(
                f"Built fallback {output_name} from {written} per-orthogroup tree files in {tree_dir}.",
                "DEBUG",
            )
            return True
        except OSError as exc:
            self.log(f"Failed to create consolidated tree file {output_path}: {exc}", "WARNING")
            if os.path.exists(output_path):
                try:
                    os.remove(output_path)
                except OSError as cleanup_exc:
                    self.log(f"Failed to remove incomplete consolidated tree file {output_path}: {cleanup_exc}", "WARNING")
            return False

    def _ensure_core_set_result_layout(self, results_folder):
        resolved_dir = os.path.join(results_folder, "Resolved_Gene_Trees")
        if not self._create_consolidated_tree_file(resolved_dir, "Resolved_Gene_Trees.txt"):
            self.log(
                f"Resolved gene tree consolidation skipped for {results_folder}; "
                "core-set analysis will rely on any per-orthogroup resolved trees already present.",
                "DEBUG",
            )
        return True

    def run(self):
        '''
        ORTHOFINDER Emms D.M. & Kelly S. (2019)

        There are two modes this task can be run in:
        a. An existing input_dir is provided with no accessions - program assumes proteomes are already present.
        b. Accessions are provided - proteomes will be staged into a clean per-run input directory.
        1. Create output location
        2. Copy proteomes into output location
        3. Run OrthoFinder orthofinder -f "$INPUT_DIR" -o "$OUTPUT_DIR"
        4. Log into database with location and library reference

        '''
        # Expand any accession variables now that the DB connection is available.
        if self.accessions:
            try:
                self.accessions = self.selector_accessions()
            except Exception as exc:  # boundary: selector resolution failure becomes this task error.
                return self.handle_exception(f"Failed to resolve accessions: {exc}", {"accessions": self.accessions})

        # Check provided with everything
        if not self.input_dir and not self.accessions:
            return self.handle_exception("Input directory is not specified or does not exist.", {"input_dir": self.input_dir})
        if not self.out_dir:
            return self.handle_exception("Output directory is not specified or does not exist.", {"out_dir": self.out_dir})

        if not self.library_id:
            if not self.library_name:
                return self.handle_exception("Either library_id or library_name must be provided.", {"library_id": self.library_id, "library_name": self.library_name})
            # Get library_id from library_name
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
            if not self.library_id:
                return self.handle_exception(f"Library name {self.library_name} not found in database.", {"library_name": self.library_name})

        if self.library_id and not self.library_name:
            self.library_name = self.db_manager.libraries.get_name(self.library_id)
            if not self.library_name:
                return self.handle_exception(f"Library ID {self.library_id} not found in database.", {"library_id": self.library_id})

        # Check the outpath exists, if not create it
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        if self.input_dir and not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)

        orthofinder_binaries_path = self.db_manager.env.get('ORTHOFINDER_BINARIES_PATH')
        orthofinder_executable = self._resolve_orthofinder_executable(orthofinder_binaries_path)
        if not orthofinder_executable:
            return self.handle_exception("OrthoFinder binaries path is not set.", {"orthofinder_binaries_path": orthofinder_binaries_path})

        input_file_count = sum(1 for _ in self._iter_proteome_filenames(self.input_dir)) if self.input_dir else 0
        usable_input_count = input_file_count if not self.accessions else 0

        # We need the input dir to contain at least two proteome files OR supply at least two proteome accessions
        if usable_input_count + len(self.accessions) < 2:
            return self.handle_exception("At least two proteome files are required either by passing an accession in the database or by providing them in the input directory.", {"input_dir": self.input_dir, "accessions": self.accessions})
        
        # Get accessions from proteomes in the input_directory
        input_dir_accessions = self.get_accessions_from_proteomes(self.input_dir) if self.input_dir and not self.accessions else []
        if self.accessions and input_file_count:
            self.log(
                "Input directory already contains proteome files, but selector accessions were supplied; "
                "PhyloODB will stage a clean per-run input directory and ignore the pre-existing shared input contents.",
                "WARNING",
            )

        # If there are accessions provided check that they  are in the database and have protein information available
        genome_info = self.db_manager.genomes.get_many(self.accessions)

        # Convert genome_info to a dict where the first item in each sublist is the key
        if isinstance(genome_info, list):
            genome_info = {row[0]: row for row in genome_info}

        genome_col_map = self.db_manager.genomes.get_column_map()
        protein_idx = genome_col_map.get("protein")
        location_idx = genome_col_map.get("location")
        if protein_idx is None or location_idx is None:
            return self.handle_exception(
                "Genome table is missing required columns.",
                {"columns": list(genome_col_map.keys())},
            )

        for accession in self.accessions:
            if accession not in genome_info:
                return self.handle_exception(f"Accession {accession} not found in genome information.", {"accession": accession})
            if not bool(genome_info[accession][protein_idx]):
                return self.handle_exception(f"Accession {accession} has no protein information available.", {"accession": accession})
            if not bool(genome_info[accession][location_idx]):
                return self.handle_exception(f"Accession {accession} has no location information available.", {"accession": accession})

        # First, detect if an equivalent run already exists for this accession set
        clean_up_paths = []
        rows = self.db_manager.orthofinder.assert_results_exist(
            self.accessions,
            mcl_inflation=self.mcl_inflation,
        )
        previous_run_id = rows[0] if rows else None
        previous_run_location = rows[1] if rows else None
        if previous_run_id and not self._has_reusable_results(previous_run_location):
            self.log(
                f"Ignoring incomplete cached OrthoFinder run id={previous_run_id}: {previous_run_location}",
                "WARNING",
            )
            previous_run_id = None
            previous_run_location = None
        try:
            expected_input_filenames = self._expected_analysis_input_filenames(genome_info, location_idx)
        except Exception as exc:  # boundary: required input/profile resolution failure becomes this task error.
            return self.handle_exception(
                "Failed to resolve the expected staged proteome filenames for this OrthoFinder run.",
                {"error": str(exc), "accessions": self.accessions},
            )

        # If we have not found the previous run in the db we can scan the folders for a matching run
        # Look for a folder within the out_dir with a name *_{library_id} and check the accessions list file
        if not previous_run_id and self.check_for_previous_run_folders:
            for folder in os.listdir(self.out_dir):
                if folder.endswith(f"_{self.library_name}"):
                    self.log(f"Scanning folder for existing OrthoFinder run: {folder}", "DEBUG")
                    folder_path = os.path.join(self.out_dir, folder)
                    prior_results_folder = glob.glob(os.path.join(folder_path, "Results*"))
                    accession_list_file = os.path.join(prior_results_folder[0], "accession_list.txt") if prior_results_folder else None
                    if accession_list_file and os.path.exists(accession_list_file):
                        try:
                            with open(accession_list_file, "r") as f:
                                file_accessions = [line.strip() for line in f.readlines()]
                            expected_accessions = list(dict.fromkeys(canonicalize_accessions(self.accessions + input_dir_accessions)))
                            if sorted(canonicalize_accessions(file_accessions)) == sorted(expected_accessions):
                                results_folder = prior_results_folder[0]
                                log_file = os.path.join(results_folder, "Log.txt")
                                logged_input_filenames = self._read_logged_species_used_filenames(log_file)
                                if not logged_input_filenames:
                                    self.log(
                                        f"Skipping folder {folder_path}: no parsable 'Species used' block found in {log_file}.",
                                        "DEBUG",
                                    )
                                    continue
                                if logged_input_filenames != expected_input_filenames:
                                    self.log(
                                        f"Skipping folder {folder_path}: logged staged input filenames do not match the current request.",
                                        "DEBUG",
                                    )
                                    continue
                                command_line = self._read_logged_command_line(log_file)
                                logged_mcl_inflation = self._parse_logged_mcl_inflation(command_line)
                                requested_inflation = None if self.mcl_inflation in (None, "") else float(self.mcl_inflation)
                                if logged_mcl_inflation != requested_inflation:
                                    self.log(
                                        f"Skipping folder {folder_path}: logged MCL inflation ({logged_mcl_inflation}) does not match the current request ({requested_inflation}).",
                                        "DEBUG",
                                    )
                                    continue
                                if not self._has_reusable_results(results_folder):
                                    self.log(
                                        f"Skipping folder {folder_path}: missing Orthogroup_Sequences in {results_folder}.",
                                        "DEBUG",
                                    )
                                    continue
                                self.log(f"Matching accession list and staged input filenames found in folder: {folder_path}", "DEBUG")
                                # Add to the database if not already present
                                old_orthofinder_id = self.db_manager.orthofinder.add_results(
                                    library_id=self.library_id,
                                    datetime="2025-09-19 14:41:04",
                                    location=results_folder,
                                    mcl_inflation=logged_mcl_inflation,
                                    command_line=command_line,
                                )
                                if not old_orthofinder_id:
                                    return self.handle_exception("Error adding existing OrthoFinder run found by scanning folders to database.", {"folder_path": folder_path})
                                if not self.db_manager.orthofinder.add_accessions(old_orthofinder_id, file_accessions):
                                    return self.handle_exception("Error adding accessions for existing OrthoFinder run found by scanning folders to database.", {"folder_path": folder_path, "accessions": file_accessions})
                                previous_run_location = results_folder
                                self.log(
                                    f"Found existing OrthoFinder run by scanning folders: id={old_orthofinder_id}, location={previous_run_location}",
                                    "DEBUG",
                                )
                                previous_run_id = old_orthofinder_id
                                break
                            else:
                                self.log(f"Accession list in folder {folder_path} does not match current accessions.", "DEBUG")
                        except (OSError, UnicodeError, ValueError) as e:
                            self.log(f"Error reading accession list file {accession_list_file}: {e}", "ERROR")
                            continue

        if previous_run_id and not self.force:
            # Reuse existing run; update timestamp implicitly when used, location should already be set
            self.log(f"Existing OrthoFinder run found (id={previous_run_id}); skipping a new run.", "INFO")
            # Nothing to add to DB; ensure inputs are cleaned if any were staged (none at this point)
            return True


        # Create a new OrthoFinder_Results row to obtain a unique run id
        requested_inflation = None if self.mcl_inflation in (None, "") else float(self.mcl_inflation)
        orthofinder_id = self.db_manager.orthofinder.add_results(
            library_id=self.library_id,
            datetime=None,
            mcl_inflation=requested_inflation,
            command_line=None,
        )
        self.log(f"Generated OrthoFinder ID: {orthofinder_id}", "DEBUG")
        if not orthofinder_id:
            return self.handle_exception("Error creating OrthoFinder results entry in database.", {"library_id": self.library_id})

        run_folder = f"{self.out_dir}/{orthofinder_id}_{self.library_name}"
        analysis_input_dir = self.input_dir
        if self.accessions:
            analysis_input_dir = os.path.join(self.out_dir, f".orthofinder_input_{orthofinder_id}")
            try:
                os.makedirs(analysis_input_dir, exist_ok=True)
            except OSError as exc:
                return self._clean_up_and_error(
                    orthofinder_id,
                    clean_up_paths,
                    "Failed to create the OrthoFinder staging directory.",
                    {"analysis_input_dir": analysis_input_dir, "error": str(exc)},
                )
            clean_up_paths.append(analysis_input_dir)

        # If a prior run exists but we are forcing, clear previous artifacts
        if previous_run_id and self.force:
            try:
                if previous_run_location and os.path.exists(previous_run_location):
                    shutil.rmtree(previous_run_location)
                if not self.db_manager.orthofinder.delete_results(previous_run_id):
                    return self._clean_up_and_error(orthofinder_id, clean_up_paths, "Error deleting previous OrthoFinder results from database.", {"previous_run_id": previous_run_id})
            except OSError as e:
                return self._clean_up_and_error(orthofinder_id, clean_up_paths, f"Error deleting previous output directory: {e}", {"run_folder": previous_run_location or run_folder})
        
        # For proteomes not already in the input directory we need to copy them there 
        for accession in self.accessions:
            location = genome_info[accession][location_idx]
            selected_profile = None
            if self.proteome_profile or self.prefer_proteome_profile:
                try:
                    selected_profile, _, proteome_path = self._resolve_proteome_input(accession, location)
                except Exception as exc:  # boundary: required proteome profile resolution failure becomes this task error.
                    return self._clean_up_and_error(
                        orthofinder_id,
                        clean_up_paths,
                        f"Failed to resolve proteome input for accession {accession}.",
                        {"accession": accession, "error": str(exc)},
                    )
                file = os.path.basename(proteome_path)
                self.log(
                    f"Using proteome profile '{selected_profile}' for accession {accession}: {proteome_path}",
                    "DEBUG",
                )
            else:
                file = self._select_proteome_filename(location, accession)
                if not file:
                    return self._clean_up_and_error(orthofinder_id, clean_up_paths, f"No valid proteome files found for accession {accession}.", {"accession": accession})
                proteome_path = os.path.join(location, file)
            staged_name = self._build_staged_proteome_filename(accession, file, selected_profile=selected_profile)
            staged_path = os.path.join(analysis_input_dir, staged_name)
            if file.endswith(".faa") or file.endswith(".fasta"):
                if not os.path.exists(proteome_path):
                    return self._clean_up_and_error(orthofinder_id, clean_up_paths, f"Proteome file for accession {accession} no longer exists.", {"accession": accession, "path": proteome_path})
                try:
                    shutil.copyfile(proteome_path, staged_path)
                except OSError as e:
                    return self._clean_up_and_error(orthofinder_id, clean_up_paths, f"Error copying proteome file for accession {accession}.", {"error": str(e)})
                clean_up_paths.append(staged_path)
            elif file.endswith(".faa.gz") or file.endswith(".fasta.gz"):
                # Decompress in its location, then copy decompressed version to input_dir
                decompressed_path = os.path.join(location, file[:-3])  # Remove .gz
                try:
                    with gzip.open(proteome_path, 'rb') as f_in, open(decompressed_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                except (OSError, EOFError, gzip.BadGzipFile) as e:
                    return self._clean_up_and_error(orthofinder_id, clean_up_paths, f"Error decompressing proteome file for accession {accession}.", {"error": str(e)})
                try:
                    shutil.copyfile(decompressed_path, staged_path)
                except OSError as e:
                    return self._clean_up_and_error(orthofinder_id, clean_up_paths, f"Error copying decompressed proteome file for accession {accession}.", {"error": str(e)})
                clean_up_paths.append(staged_path)
                # Optionally clean up decompressed file after run
                clean_up_paths.append(decompressed_path)
            
        # Now rename the accessions and sequences for analysis
        # SKIPPED but may come back to.

        # We are now ready to run orthofinder
        command = [
            orthofinder_executable,
            "-f", analysis_input_dir,
            "-t", str(self.REQUIRED_THREADS),
            "-o", run_folder,
        ]
        if requested_inflation is not None:
            command.extend(["-I", str(requested_inflation)])
        if not self.db_manager.orthofinder.update_run_metadata(
            orthofinder_id,
            mcl_inflation=requested_inflation,
            command_line=" ".join(command),
        ):
            return self._clean_up_and_error(
                orthofinder_id,
                clean_up_paths,
                "Failed to persist OrthoFinder run metadata before execution.",
            )

        self.log(f"Running OrthoFinder with command: {' '.join(command)}", "DEBUG")

        self.log(
            f"Starting OrthoFinder run for {len(all_accessions) if 'all_accessions' in locals() else len(self.accessions) + len(input_dir_accessions)} "
            f"proteomes in library '{self.library_name}'.",
            "INFO",
        )

        try:
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError as exc:
            return self._clean_up_and_error(
                orthofinder_id,
                clean_up_paths,
                f"Failed to launch OrthoFinder: {exc}",
                {"command": command, "error": str(exc)},
            )

        if result.returncode != 0:
            self.error(f"OrthoFinder failed: {result.stderr}")
            if not self._clean_up(clean_up_paths):
                    self.error("Error cleaning up files after OrthoFinder run.")
            return self.handle_exception("OrthoFinder command failed.", {"returncode": result.returncode, "stderr": result.stderr})
        self.log("OrthoFinder completed successfully.", "INFO")
        # time_stamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        # Now need to enter results into database
        # Persist location for this specific run folder, and then record accessions used
        all_accessions = list(dict.fromkeys(canonicalize_accessions(list(self.accessions) + list(input_dir_accessions))))
        # Find the OrthoFinder results folder (should be run_folder/Results*/ using glob)
        results_folders = glob.glob(os.path.join(run_folder, "Results*"))
        if not results_folders:
            return self._clean_up_and_error(orthofinder_id, clean_up_paths, "OrthoFinder results folder not found.", {"run_folder": run_folder})
        results_folder = results_folders[0]

        orthogroup_sequences_dir = os.path.join(results_folder, "Orthogroup_Sequences")
        if not os.path.isdir(orthogroup_sequences_dir):
            return self._clean_up_and_error(
                orthofinder_id,
                clean_up_paths,
                "OrthoFinder completed but did not produce Orthogroup_Sequences.",
                {"results_folder": results_folder},
            )

        if not self._ensure_core_set_result_layout(results_folder):
            return self._clean_up_and_error(orthofinder_id, clean_up_paths, "Failed to normalize OrthoFinder result layout.", {"results_folder": results_folder})

        if not self.db_manager.orthofinder.update_location(orthofinder_id, results_folder) or not self.db_manager.orthofinder.add_accessions(orthofinder_id, all_accessions):
            return self._clean_up_and_error(orthofinder_id, clean_up_paths, "Error adding OrthoFinder results to database")
        try:
            self.db_manager.artifacts.register(
                owner_type="orthofinder_run",
                owner_id=orthofinder_id,
                artifact_type="orthofinder_results_dir",
                path=results_folder,
                is_dir=True,
                format="directory",
                metadata={
                    "library_id": self.library_id,
                    "library_name": self.library_name,
                    "mcl_inflation": requested_inflation,
                    "command_line": " ".join(command),
                },
            )
        except Exception as exc:  # boundary: optional artifact catalog metadata; OrthoFinder results already persisted.
            self.log(f"Failed to register OrthoFinder artifact: {exc}", "WARNING")
        
        if not self._create_accessions_list_file(all_accessions, os.path.join(results_folder, "accession_list.txt")):
            return self._clean_up_and_error(orthofinder_id, clean_up_paths, "Error creating accession list file.", {"run_folder": run_folder})

        self.log(f"Recorded OrthoFinder results as run ID {orthofinder_id}.", "INFO")

        if not self._clean_up(clean_up_paths):
            self.error("Error cleaning up files after OrthoFinder run.")

        return True
