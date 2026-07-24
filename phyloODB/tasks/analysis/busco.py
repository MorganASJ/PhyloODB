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
from ...database import DBManager
from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE, RAW_PROFILE, is_staged_busco_input_path, resolve_profile_selector
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
)

class DownloadBuscoLibraryTask(Task):
    '''A Task that downloads the BUSCO library'''
    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.lineage = self.data.get("lineage", None)
        self.libraries_dir = self.data.get("libraries_dir", None)
        self.busco_path = self.data.get("busco_path", None)
        self.parent_library_name = self.data.get("parent_library_name", None)
        self.coverage = self.data.get("coverage", None)
        self.size = self.data.get("size", None)
        self.debug_skip_dl = self.data.get("debug_skip_dl", False)

        
    def run(self):
        # Check we have everything we need
        if not self.busco_path:
            self.busco_path = self.db_manager.env.get("BUSCO_BINARIES_PATH")
            if not self.busco_path:
                return self.handle_exception("BUSCO binaries path is not specified.", {"busco_path": self.busco_path})
        if not self.lineage:
            return self.handle_exception("Lineage is not specified.", {"lineage": self.lineage})
        if not self.libraries_dir:
            self.libraries_dir = self.db_manager.storage.get_root_base("libraries")
            if not self.libraries_dir:
                return self.handle_exception("Lineage directory is not specified.", {"libraries_dir": self.libraries_dir})

        # If the folder for this lineage already exists delete it
        lineage_dir = os.path.join(self.libraries_dir, f"lineages/{self.lineage}")
        if os.path.exists(lineage_dir) and not self.debug_skip_dl:
            self.log(f"BUSCO library {self.lineage} already exists. Deleting existing directory to redownload.", "WARNING")
            try:
                shutil.rmtree(lineage_dir)
            except OSError as e:
                # self.error(f"Failed to delete existing BUSCO library directory {lineage_dir}: {e}")
                return self.handle_exception(f"Failed to delete existing BUSCO library directory {lineage_dir}.", e)

        command = [
            "busco",
            "--download", self.lineage,
            "--download_path", self.libraries_dir,
        ]

        self.log(f"Downloading BUSCO library {self.lineage}", "INFO")

        if not self.debug_skip_dl:

            try:
                result = subprocess.run(
                    command,
                    check=True,
                    capture_output=True,
                    text=True
                )
                self.log(f"Downloaded BUSCO library {self.lineage} successfully.")

            except FileNotFoundError as e:
                # self.error(f"BUSCO not found. Ensure BUSCO is installed at the set location: {self.busco_path}")
                self.handle_exception("BUSCO not found at file path set in environment variables.", e)
                return "ERROR"
            except PermissionError as e:
                # self.error(f"Permission error: {e}")
                return self.handle_exception("Permission error when trying to run BUSCO.", e)
            except subprocess.CalledProcessError as e:
                err = f"BUSCO failed with exit code {e.returncode}"
                # self.error(f"stderr: {e.stderr}")
                # self.handle_exception("BUSCO command failed.", e)
                # Optional crude hints
                if "download" in e.stderr.lower():
                    err += ("\nLikely a network or server issue.")
                elif "permission" in e.stderr.lower():
                    err += ("\nLikely a write permission issue.")
                # self.error(err)
                return self.handle_exception(err, e)
        
        # Check lineage_dir now exists
        if not os.path.exists(lineage_dir):
            # self.error(f"BUSCO library directory not found after download: {lineage_dir}")
            return self.handle_exception(f"BUSCO library directory not found after download", {"lineage_dir": lineage_dir})
        
        # TEMPORARY CODE FOR LIBRARY UNTIL TAXONOMY MODULE IS INSTALLED:
        # if "metazoa" in self.lineage:
        #     self.coverage = 33208

        # Only get missing values from config file
        config_file = os.path.join(lineage_dir, "dataset.cfg")
        # if (self.coverage is None or self.size is None):
        if not os.path.exists(config_file):
            return self.handle_exception(f"Config file not found: {config_file}")
        odb_version = None
        with open(config_file) as f:
            for line in f:
                if self.size is None and line.startswith("number_of_BUSCOs="):
                    self.size = int(line.split("=")[1])
                if self.coverage is None and line.startswith("ncbi_taxid="):
                    self.coverage = int(line.split("=")[1])
                if line.startswith("OrthoDB_version="):
                    odb_version = line.split("=")[1].strip()
                if self.size is not None and self.coverage is not None and odb_version is not None:
                    break

        if not self.size or not self.coverage or not odb_version:
            return self.handle_exception("Missing coverage or size information. Could not retrieve from the config file.", {"coverage": self.coverage, "size": self.size, "odb_version": odb_version})

        self.log(f"BUSCO library {self.lineage} info: coverage={self.coverage}, size={self.size}, odb_version={odb_version}", "DEBUG")

        library_id = self.db_manager.libraries.add(
            library_name=self.lineage,
            taxid=self.coverage,
            size=self.size,
            location=lineage_dir,
            parent_id=None,
            ref_accessions=None,
            odb_version=odb_version,
        )

        if not library_id:
            return self.handle_exception("Failed to create or update BUSCO library record.", {"lineage": self.lineage})

        self.log(f"Registered BUSCO library {self.lineage} (ID {library_id}) in database.", "INFO")
        try:
            self.db_manager.artifacts.register(
                owner_type="library",
                owner_id=library_id,
                artifact_type="library_root",
                path=lineage_dir,
                is_dir=True,
                format="directory",
                metadata={"library_id": library_id, "library_name": self.lineage},
            )
            self.db_manager.artifacts.register(
                owner_type="library",
                owner_id=library_id,
                artifact_type="library_dataset_cfg",
                path=config_file,
                format="cfg",
                metadata={"library_id": library_id, "library_name": self.lineage},
            )
        except Exception as exc:  # boundary: artifact persistence is required for registered BUSCO libraries
            return self.handle_exception(
                f"Failed to register lineage library artifacts for {self.lineage}: {exc}",
                {"lineage": self.lineage, "library_id": library_id, "lineage_dir": lineage_dir},
            )

        # Get BUSCO family information and submit to database (BUSCO_Description table)
        links_to_ODB_files = glob.glob(os.path.join(lineage_dir, "links_*.txt"))
        if links_to_ODB_files:
            odb_links = links_to_ODB_files[0]
            try:
                with open(odb_links, newline="") as f:
                    reader = csv.reader(f, delimiter="\t")
                    descriptions = []
                    for row in reader:
                        if len(row) == 3:
                           descriptions.append((row[0], library_id, row[1], row[2]))  # (busco_id, description, link)
                           self.log(f"Found BUSCO description: {row[0]} - {row[1]}", "DEBUG")
                        else:
                            return self.handle_exception("Invalid format in BUSCO links to ODB file.", {"file": links_to_ODB_files[0]})
            except (OSError, UnicodeError, csv.Error) as e:
                return self.handle_exception("Failed to read BUSCO links to ODB file.", e)
            if descriptions and not self.db_manager.libraries.add_busco_descriptions(descriptions):
                return self.handle_exception("Failed to add BUSCO descriptions to database.", {"library_id": library_id, "descriptions": descriptions})
        else:
            self.handle_exception("BUSCO links to ODB file not found.", {"lineage": self.lineage})

        return True

class BatchBuscoTask(Task):
    '''A Task that runs BUSCO on multiple assemblies over subtasks'''
    @classmethod
    def default_thread_count(cls, registry_required_threads: int, daemon_max_threads: int) -> int:
        return 1

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.lineage = self.data.get("lineage")
        self.library = self.data.get("library", self.lineage)
        fmt = str(self.data.get("format") or "auto").strip().lower()
        if fmt in ("genome", "nucleotide", "geno"):
            self.format = "genome"
        elif fmt in ("protein", "proteome", "prot"):
            self.format = "protein"
        elif fmt in ("auto",):
            self.format = "auto"
        else:
            self.format = fmt
        self.accessions = self.data.get("accessions")
        self.selector_requested_accessions = self.data.get("selector_requested_accessions") or []
        self.selector_skipped_accessions = self.data.get("selector_skipped_accessions") or []
        self.output_path = self.data.get("output_path", None)
        self.force = self.data.get("force", False)
        self.pipeline = str(self.data.get("pipeline") or "auto").strip().lower()
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        self.augustus_evalue = self.data.get("augustus_evalue")
        self.augustus_limit = self.data.get("augustus_limit")
        self.augustus_long = self.data.get("augustus_long")
        self.augustus_species = self.data.get("augustus_species")
        self.augustus_parameters = self.data.get("augustus_parameters")
        self.metaeuk_parameters = self.data.get("metaeuk_parameters")
        self.metaeuk_rerun_parameters = self.data.get("metaeuk_rerun_parameters")
        self.miniprot_parameters = self.data.get("miniprot_parameters")
        self.keep_miniprot_ref_file_raw = self.data.get("keep_miniprot_ref_file", None)
        self.keep_miniprot_ref_file = bool(self.keep_miniprot_ref_file_raw) if self.keep_miniprot_ref_file_raw is not None else False
        self.max_concurrent = max(1, int(self.data.get("max_concurrent", 1) or 1))
        self.busco_lib_wait_seconds = self.data.get("busco_lib_wait_seconds", 0)
        self.busco_lib_retries = self.data.get("busco_lib_retries", 0)
        self.stage = checkpoint if checkpoint is not None else 0

    def _busco_missing(self, *, respect_force: bool = True):
        """Get accessions that are missing BUSCO results.

        When ``respect_force`` is true, ``force`` makes every accession eligible
        for queueing so BUSCO reruns are scheduled even if prior results exist.
        Completion checks must ignore that override and inspect the real BUSCO
        state instead.
        """
        if respect_force and self.force:
            # Force rerun: treat all accessions as missing so they re-queue.
            return list(self.accessions or [])
        def _has_expected_run_folder(result_dir: Optional[str], pipeline_name: str) -> bool:
            if not result_dir or not os.path.isdir(result_dir):
                return False
            run_dirs = glob.glob(os.path.join(result_dir, "run_*"))
            if not run_dirs:
                return False
            expected = f"run_{pipeline_name}_{self.lineage}"
            legacy = f"run_{self.lineage}"
            for run_dir in run_dirs:
                base = os.path.basename(run_dir)
                if base.startswith(expected):
                    return True
                # Legacy pre-pipeline naming is treated as miniprot.
                if pipeline_name == "miniprot" and base.startswith(legacy):
                    return True
            return False
        # Use the library_id to scope BUSCO results correctly
        library_id = self.db_manager.libraries.get_id(self.library) if self.library else None
        requested_pipeline = str(self.pipeline or "auto").strip().lower()
        requested_format = str(self.format or "auto").strip().lower()
        explicit_profile = str(self.proteome_profile or "").strip() or None
        effective_requested_format = requested_format
        if explicit_profile and effective_requested_format == "auto":
            # Proteome-profile selectors only apply to protein-mode BUSCO runs.
            effective_requested_format = "protein"
        if library_id is not None and (
            requested_pipeline != "auto"
            or requested_format != "auto"
            or explicit_profile is not None
        ):
            missing = []
            for acc in self.accessions or []:
                clauses = ["r.accession = ?", "r.library_id = ?"]
                params = [str(acc), int(library_id)]
                if requested_pipeline != "auto":
                    clauses.append("LOWER(r.pipeline) = ?")
                    params.append(requested_pipeline)
                if effective_requested_format in {"protein", "genome"}:
                    clauses.append("LOWER(r.input_mode) = ?")
                    params.append(effective_requested_format)
                sql = (
                    "SELECT r.result_dir, pp.profile_name "
                    "FROM BUSCO_Runs r "
                    "LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id "
                    f"WHERE {' AND '.join(clauses)} "
                    "ORDER BY r.run_id DESC"
                )
                self.db_manager.cursor.execute(sql, tuple(params))
                rows = self.db_manager.cursor.fetchall() or []
                exists = False
                for row in rows:
                    result_dir = row[0] if row else None
                    row_profile = row[1] if row and len(row) > 1 else None
                    if effective_requested_format == "protein" and explicit_profile:
                        if not self.db_manager.proteomes.profile_matches_selector(
                            str(acc),
                            str(row_profile) if row_profile is not None else None,
                            explicit_profile,
                        ):
                            continue
                    if requested_pipeline == "auto":
                        if result_dir and os.path.isdir(str(result_dir)):
                            exists = True
                            break
                    elif _has_expected_run_folder(str(result_dir) if result_dir is not None else None, requested_pipeline):
                        exists = True
                        break
                if not exists:
                    missing.append(acc)
            return missing
        present = self.db_manager.busco.get_processed_accessions(library_id)
        missing = [acc for acc in self.accessions if acc not in present]
        return missing
    
    def queue_busco(self):
        if self.selector_skipped_accessions:
            self.log(
                "Skipped accessions before BUSCO queueing because BatchBuscoTask only targets "
                f"assemblies with local genome status >= 1: {self.selector_skipped_accessions}",
                "WARNING",
            )
        missing = self._busco_missing()
        self.log(f"BUSCO missing: {missing}", "DEBUG")
        queued = False
        child_threads = max(1, (int(self.REQUIRED_THREADS or 1) + self.max_concurrent - 1) // self.max_concurrent)
        for acc in missing:
            self.queue_subtask(
                job_type=4,
                status="P",
                priority=1,
                data={
                    "lineage": self.lineage,
                    "library": self.library,
                    "format": self.format,
                    "accession": acc,
                    "output_path": self.output_path,
                    "force": self.force,
                    "pipeline": self.pipeline,
                    "proteome_profile": self.proteome_profile,
                    "prefer_proteome_profile": self.prefer_proteome_profile,
                    "augustus_evalue": self.augustus_evalue,
                    "augustus_limit": self.augustus_limit,
                    "augustus_long": self.augustus_long,
                    "augustus_species": self.augustus_species,
                    "augustus_parameters": self.augustus_parameters,
                    "metaeuk_parameters": self.metaeuk_parameters,
                    "metaeuk_rerun_parameters": self.metaeuk_rerun_parameters,
                    "miniprot_parameters": self.miniprot_parameters,
                    "keep_miniprot_ref_file": self.keep_miniprot_ref_file,
                    "required_threads": child_threads,
                },
            )
            queued = True
        return queued

    def busco_done(self):
        # Check that all accessions have BUSCO results
        missing = self._busco_missing(respect_force=False)
        return (len(missing) == 0)

    def busco_incomplete_message(self):
        missing = self._busco_missing(respect_force=False)
        summary = f"BUSCO phase incomplete. Missing BUSCO results for {len(missing)} accession(s): {', '.join(missing)}"
        stack = ""
        return (summary, stack)

    def run(self):
        # Queue BUSCO tasks for each accession
        if not self.accessions or not isinstance(self.accessions, list):
            return self.handle_exception("Accessions list is not specified or invalid.", {"accessions": self.accessions}) 
        
        # Check lineage is valid
        if not self.lineage:
            return self.handle_exception("Lineage is not specified.", {"lineage": self.lineage})
        
        if not self.db_manager.busco.is_lineage_downloaded(self.lineage):
            return self.handle_exception(f"BUSCO lineage '{self.lineage}' is not downloaded.", {"lineage": self.lineage})
    

        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=self.queue_busco,
            done_fn=self.busco_done,
            wait_seconds=self.busco_lib_wait_seconds,
            retry_key=None,
            max_retries=self.busco_lib_retries,
            incomplete_message_fn=self.busco_incomplete_message,
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            # Suspended; daemon will resume when subtasks complete
            return False
        # continue

        return True

class BuscoTask(Task):
    '''A Task that runs BUSCO on an assembly'''
    @classmethod
    def default_thread_count(cls, registry_required_threads: int, daemon_max_threads: int) -> int:
        daemon_max_threads = max(int(daemon_max_threads or 1), 1)
        if daemon_max_threads >= 16:
            return 8
        if daemon_max_threads >= 8:
            return 4
        if daemon_max_threads >= 4:
            return 2
        return 1

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.lineage = self.data.get("lineage")
        self.library = self.data.get("library", self.lineage)
        fmt = (str(self.data.get("format") or "auto").strip().lower())
        if fmt in ("genome", "nucleotide", "geno"):
            self.format = "genome"
        elif fmt in ("protein", "proteome", "prot"):
            self.format = "protein"
        elif fmt in ("auto",):
            self.format = "auto"
        else:
            self.format = fmt
        self.accession = self.data.get("accession")
        self.output_path = self.data.get("output_path")
        self.force = self.data.get("force", False)
        self.pipeline = str(self.data.get("pipeline") or "auto").strip().lower()
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        self.expected_proteome_profile_id = self.data.get("expected_proteome_profile_id")
        self.expected_proteome_checksum = self.data.get("expected_proteome_checksum")
        self.augustus_evalue = self.data.get("augustus_evalue")
        self.augustus_limit = self.data.get("augustus_limit")
        self.augustus_long = self.data.get("augustus_long")
        self.augustus_species = self.data.get("augustus_species")
        self.augustus_parameters = self.data.get("augustus_parameters")
        self.metaeuk_parameters = self.data.get("metaeuk_parameters")
        self.metaeuk_rerun_parameters = self.data.get("metaeuk_rerun_parameters")
        self.miniprot_parameters = self.data.get("miniprot_parameters")
        self.busco_lib_wait_seconds = self.data.get("busco_lib_wait_seconds", 0)
        self.busco_lib_retries = self.data.get("busco_lib_retries", 0)
        self.keep_miniprot_ref_file_raw = self.data.get("keep_miniprot_ref_file", None)
        self.keep_miniprot_ref_file = bool(self.keep_miniprot_ref_file_raw) if self.keep_miniprot_ref_file_raw is not None else False
        self.stage = checkpoint if checkpoint is not None else 0

    def _resolve_requested_proteome_profile(self, accession: str) -> str:
        if self.proteome_profile:
            if self.proteome_profile == DEFAULT_CLEAN_PROFILE:
                default_cleaned = self.db_manager.proteomes.get_default_cleaned_profile_name(accession)
                if default_cleaned:
                    return default_cleaned
            return self.proteome_profile
        default_profile = self.db_manager.proteomes.get_default_profile_name(accession)
        if default_profile:
            return default_profile
        default_cleaned = self.db_manager.proteomes.get_default_cleaned_profile_name(accession)
        if default_cleaned:
            return default_cleaned
        return RAW_PROFILE

    def _resolve_proteome_profile_row(self, accession: str, profile_name: str):
        if profile_name == RAW_PROFILE:
            raw_id = self.db_manager.proteomes.ensure_raw_profile(accession, is_default=False)
            if raw_id is None:
                return None
            return self.db_manager.proteomes.get(int(raw_id))
        return self.db_manager.proteomes.get_profile(accession, profile_name)

    def _iter_proteome_profile_candidates(self, accession: str):
        seen: set[str] = set()

        def add(profile_name: Optional[str], *, required: bool = False):
            token = str(profile_name or "").strip()
            if not token or token in seen:
                return
            seen.add(token)
            yield token, required

        if self.proteome_profile:
            yield from add(self._resolve_requested_proteome_profile(accession), required=True)
            return

        if self.prefer_proteome_profile:
            preferred = self.db_manager.proteomes.resolve_selector_profile_name(accession, self.prefer_proteome_profile)
            yield from add(preferred or self.prefer_proteome_profile)

        yield from add(self.db_manager.proteomes.get_default_profile_name(accession))
        yield from add(self.db_manager.proteomes.get_default_cleaned_profile_name(accession))
        yield from add(RAW_PROFILE)

        for row in self.db_manager.proteomes.list_profiles(accessions=[accession]):
            profile_name = str(row[2] or "").strip() if row and len(row) > 2 else ""
            yield from add(profile_name)

    def _resolve_proteome_input(self, accession: str, genome_path: str) -> tuple[str, Optional[int], str]:
        attempted: list[str] = []
        for profile_name, required in self._iter_proteome_profile_candidates(accession):
            attempted.append(profile_name)
            row = self._resolve_proteome_profile_row(accession, profile_name)
            if row is None:
                if required:
                    raise FileNotFoundError(
                        f"Proteome profile '{profile_name}' does not exist for accession '{accession}'."
                    )
                continue
            path = self.db_manager.proteomes.resolve_path(row)
            if path and os.path.exists(path):
                if attempted[:-1]:
                    self.log(
                        f"Proteome profile fallback for {accession}: using '{profile_name}' after skipping unreadable profiles {attempted[:-1]}.",
                        "WARNING",
                    )
                return profile_name, int(row[0]), path
            if required:
                raise FileNotFoundError(
                    f"Proteome profile '{profile_name}' has no readable artifact for accession '{accession}'."
                )

        attempted_display = ", ".join(f"'{name}'" for name in attempted) or "<none>"
        raise FileNotFoundError(
            f"No readable proteome profile artifact found for accession '{accession}' (attempted: {attempted_display})."
        )

    def _parse_parameters_csv(self, value: Optional[str]) -> list[str]:
        if value is None:
            return []
        raw = str(value).strip()
        if not raw:
            return []
        return [token.strip() for token in raw.split(",") if token.strip()]

    def _resolve_pipeline_name(self, effective_format: str) -> str:
        requested = str(self.pipeline or "auto").strip().lower()
        if requested in {"miniprot", "metaeuk", "augustus"}:
            return requested
        env_default = str(self.db_manager.env.get("DEFAULT_BUSCO_PIPELINE") or "miniprot").strip().lower()
        if env_default in {"miniprot", "metaeuk", "augustus"}:
            if effective_format == "protein" and env_default in {"miniprot", "metaeuk", "augustus"}:
                # In protein mode BUSCO may ignore genome-pipeline selectors; still retain chosen label for traceability.
                return env_default
            return env_default
        return "miniprot"

    def _resolve_pipeline_args(self, pipeline_name: str) -> tuple[list[str], dict, dict]:
        args: list[str] = []
        effective: dict = {"pipeline": pipeline_name}
        source: dict = {"pipeline": "task" if str(self.pipeline or "").strip().lower() not in {"", "auto"} else "env"}

        if pipeline_name == "augustus":
            args.append("--augustus")
            evalue = self.augustus_evalue if self.augustus_evalue is not None else self.db_manager.env.get("BUSCO_AUGUSTUS_EVALUE")
            limit = self.augustus_limit if self.augustus_limit is not None else self.db_manager.env.get("BUSCO_AUGUSTUS_LIMIT")
            long_mode = self.augustus_long if self.augustus_long is not None else self.db_manager.env.get("BUSCO_AUGUSTUS_LONG")
            species = self.augustus_species if self.augustus_species is not None else self.db_manager.env.get("BUSCO_AUGUSTUS_SPECIES")
            parameters = self.augustus_parameters if self.augustus_parameters is not None else self.db_manager.env.get("BUSCO_AUGUSTUS_PARAMETERS")
            if evalue is not None:
                args.extend(["-e", str(evalue)])
            if limit is not None:
                args.extend(["--limit", str(limit)])
            if bool(long_mode):
                args.append("--long")
            if species:
                args.extend(["--augustus_species", str(species)])
            extra = self._parse_parameters_csv(parameters)
            if extra:
                args.extend(["--augustus_parameters", ",".join(extra)])
            effective.update(
                {
                    "evalue": evalue,
                    "limit": limit,
                    "long": bool(long_mode),
                    "species": species,
                    "parameters": extra,
                }
            )
            source.update(
                {
                    "evalue": "task" if self.augustus_evalue is not None else "env",
                    "limit": "task" if self.augustus_limit is not None else "env",
                    "long": "task" if self.augustus_long is not None else "env",
                    "species": "task" if self.augustus_species is not None else "env",
                    "parameters": "task" if self.augustus_parameters is not None else "env",
                }
            )
        elif pipeline_name == "metaeuk":
            args.append("--metaeuk")
            params = self.metaeuk_parameters if self.metaeuk_parameters is not None else self.db_manager.env.get("BUSCO_METAEUK_PARAMETERS")
            rerun = self.metaeuk_rerun_parameters if self.metaeuk_rerun_parameters is not None else self.db_manager.env.get("BUSCO_METAEUK_RERUN_PARAMETERS")
            pvals = self._parse_parameters_csv(params)
            rvals = self._parse_parameters_csv(rerun)
            if pvals:
                args.extend(["--metaeuk_parameters", ",".join(pvals)])
            if rvals:
                args.extend(["--metaeuk_rerun_parameters", ",".join(rvals)])
            effective.update({"parameters": pvals, "rerun_parameters": rvals})
            source.update(
                {
                    "parameters": "task" if self.metaeuk_parameters is not None else "env",
                    "rerun_parameters": "task" if self.metaeuk_rerun_parameters is not None else "env",
                }
            )
        else:
            args.append("--miniprot")
            params = self.miniprot_parameters if self.miniprot_parameters is not None else self.db_manager.env.get("BUSCO_MINIPROT_PARAMETERS")
            pvals = self._parse_parameters_csv(params)
            if pvals:
                # BUSCO accepts miniprot parameters through generic extra args in recent versions.
                args.extend(pvals)
            effective.update({"parameters": pvals})
            source.update({"parameters": "task" if self.miniprot_parameters is not None else "env"})
        return args, effective, source

    def _resolve_effective_format(self, genome_path: str) -> str:
        """Resolve BUSCO mode for this accession, including 'auto' mode."""
        requested = (self.format or "auto").strip().lower()
        if requested in {"protein", "genome"}:
            return requested
        if requested not in {"auto", "nucleotide"}:
            raise ValueError("Format must be one of 'auto', 'protein', 'genome' (or nucleotide aliases).")

        has_fna = False
        for fname in os.listdir(genome_path):
            low = fname.lower()
            if low.endswith(".fna") or low.endswith(".fna.gz"):
                has_fna = True
        has_faa = (
            self.db_manager.proteomes.get_default_profile(self.accession) is not None
            or self.db_manager.proteomes.ensure_raw_profile(self.accession, is_default=False) is not None
        )

        if has_fna and not has_faa:
            self.log(f"Auto BUSCO format for {self.accession}: genome (only nucleotide file present).", "INFO")
            return "genome"
        if has_faa and not has_fna:
            self.log(f"Auto BUSCO format for {self.accession}: protein (only proteome file present).", "INFO")
            return "protein"

        default_fmt_raw = str(self.db_manager.env.get("DEFAULT_BUSCO_FORMAT") or "protein").strip().lower()
        if default_fmt_raw in ("nucleotide", "geno"):
            default_fmt = "genome"
        elif default_fmt_raw in ("protein", "proteome", "prot"):
            default_fmt = "protein"
        elif default_fmt_raw in ("genome",):
            default_fmt = "genome"
        else:
            default_fmt = "protein"
            self.log(
                f"Invalid DEFAULT_BUSCO_FORMAT='{default_fmt_raw}'. Falling back to 'protein'.",
                "WARNING",
            )

        if has_fna and has_faa:
            self.log(
                f"Auto BUSCO format for {self.accession}: {default_fmt} (both nucleotide and proteome present; using DEFAULT_BUSCO_FORMAT).",
                "INFO",
            )
        else:
            self.log(
                f"Auto BUSCO format for {self.accession}: {default_fmt} (no BUSCO input file detected yet; using DEFAULT_BUSCO_FORMAT).",
                "WARNING",
            )
        return default_fmt

    def get_json_file(self, results_dir):
        # BUSCO outputs a file named short_summary_<run_name>.json in the results directory
        if not os.path.exists(results_dir):
            return None

        json_files = []
        for fname in os.listdir(results_dir):
            if fname.startswith("short_summary") and fname.endswith(".json"):
                json_files.append(os.path.join(results_dir, fname))

        if not json_files:
            self.log(f"BUSCO results JSON file does not exist: {results_dir}/short_summary*.json", "ERROR")
            return None
        if len(json_files) > 1:
            self.log(f"Multiple BUSCO results JSON files found: {json_files}", "ERROR")
            return False

        return json_files[0]

    def _pipeline_run_dir_exists(self, results_dir: str, pipeline_name: str, lineage: str) -> bool:
        if not results_dir or not os.path.isdir(results_dir):
            return False
        prefix = f"run_{pipeline_name}_{lineage}"
        legacy_prefix = f"run_{lineage}"
        for run_dir in glob.glob(os.path.join(results_dir, "run_*")):
            base = os.path.basename(run_dir)
            if base.startswith(prefix):
                return True
            # Legacy naming (pre pipeline-aware folders) is assumed miniprot.
            if pipeline_name == "miniprot" and base.startswith(legacy_prefix):
                return True
        return False

    def _find_existing_pipeline_results_dir(
        self,
        output_path: str,
        lineage: str,
        pipeline_name: str,
        *,
        accession: Optional[str] = None,
        library_id: Optional[int] = None,
        input_mode: Optional[str] = None,
        proteome_profile_id: Optional[int] = None,
    ) -> Optional[str]:
        """Find an existing results root for the requested run identity.

        We first consult BUSCO_Runs metadata (accession/library/pipeline/mode/profile)
        so profile-specific protein runs do not collide. A filesystem scan is used as
        a legacy fallback when strict metadata matching is unavailable.
        """
        mode = str(input_mode or "").strip().lower() or None

        # Strict match against completed BUSCO run metadata.
        if accession and library_id is not None:
            clauses = [
                "accession = ?",
                "library_id = ?",
                "LOWER(COALESCE(status, 'completed')) = 'completed'",
                "LOWER(pipeline) = ?",
            ]
            params = [str(accession), int(library_id), str(pipeline_name).strip().lower()]
            if mode in {"protein", "genome"}:
                clauses.append("LOWER(input_mode) = ?")
                params.append(mode)
            if mode == "protein":
                if proteome_profile_id is not None:
                    clauses.append("proteome_profile_id = ?")
                    params.append(int(proteome_profile_id))
                else:
                    clauses.append("proteome_profile_id IS NULL")

            sql = f"SELECT result_dir FROM BUSCO_Runs WHERE {' AND '.join(clauses)} ORDER BY completed_at DESC, run_id DESC"
            self.db_manager.cursor.execute(sql, tuple(params))
            for row in self.db_manager.cursor.fetchall() or []:
                candidate = str(row[0]) if row and row[0] is not None else ""
                if not candidate or not os.path.isdir(candidate):
                    continue
                if not self.get_json_file(candidate):
                    continue
                if self._pipeline_run_dir_exists(candidate, pipeline_name, lineage):
                    return candidate

        # Legacy fallback: scan lineage folders on disk.
        # Do not use this fallback for explicit protein-profile runs; directory names
        # are not profile-scoped and can cause false skips across profiles.
        if mode == "protein" and proteome_profile_id is not None:
            return None
        if not output_path or not os.path.isdir(output_path):
            return None
        candidates = sorted(glob.glob(os.path.join(output_path, f"{lineage}_results*")))
        for candidate in candidates:
            if not os.path.isdir(candidate):
                continue
            if not self.get_json_file(candidate):
                continue
            if self._pipeline_run_dir_exists(candidate, pipeline_name, lineage):
                return candidate
        return None

    def _find_existing_pipeline_run(
        self,
        output_path: str,
        lineage: str,
        pipeline_name: str,
        *,
        accession: Optional[str] = None,
        library_id: Optional[int] = None,
        input_mode: Optional[str] = None,
        proteome_profile_id: Optional[int] = None,
    ) -> tuple[Optional[str], Optional[int]]:
        mode = str(input_mode or "").strip().lower() or None

        if accession and library_id is not None:
            clauses = [
                "accession = ?",
                "library_id = ?",
                "LOWER(COALESCE(status, 'completed')) = 'completed'",
                "LOWER(pipeline) = ?",
            ]
            params = [str(accession), int(library_id), str(pipeline_name).strip().lower()]
            if mode in {"protein", "genome"}:
                clauses.append("LOWER(input_mode) = ?")
                params.append(mode)
            if mode == "protein":
                if proteome_profile_id is not None:
                    clauses.append("proteome_profile_id = ?")
                    params.append(int(proteome_profile_id))
                else:
                    clauses.append("proteome_profile_id IS NULL")

            sql = (
                f"SELECT run_id, result_dir FROM BUSCO_Runs WHERE {' AND '.join(clauses)} "
                "ORDER BY completed_at DESC, run_id DESC"
            )
            self.db_manager.cursor.execute(sql, tuple(params))
            for row in self.db_manager.cursor.fetchall() or []:
                run_id = int(row[0]) if row and row[0] is not None else None
                candidate = str(row[1]) if row and len(row) > 1 and row[1] is not None else ""
                if not candidate or not os.path.isdir(candidate):
                    continue
                if not self.get_json_file(candidate):
                    continue
                if self._pipeline_run_dir_exists(candidate, pipeline_name, lineage):
                    return candidate, run_id

        result_dir = self._find_existing_pipeline_results_dir(
            output_path,
            lineage,
            pipeline_name,
            accession=accession,
            library_id=library_id,
            input_mode=input_mode,
            proteome_profile_id=proteome_profile_id,
        )
        return result_dir, None

    def _cleanup_staged_busco_input(self, staged_input_path: Optional[str], *, preserve: bool = False) -> None:
        if preserve or not staged_input_path:
            return
        try:
            if os.path.exists(staged_input_path):
                os.remove(staged_input_path)
        except OSError as exc:
            self.log(f"Failed to remove staged BUSCO input {staged_input_path}: {exc}", "WARNING")

    def _cleanup_stale_staged_busco_inputs(self, genome_path: str) -> None:
        if not genome_path or not os.path.isdir(genome_path):
            return
        removed = 0
        for fname in sorted(os.listdir(genome_path)):
            if not is_staged_busco_input_path(fname):
                continue
            path = os.path.join(genome_path, fname)
            if not os.path.isfile(path):
                continue
            try:
                os.remove(path)
                removed += 1
            except FileNotFoundError:
                continue
            except OSError as exc:
                self.log(f"Failed to remove stale staged BUSCO input {path}: {exc}", "WARNING")
        if removed:
            self.log(f"Removed {removed} stale staged BUSCO input file(s) from {genome_path}.", "DEBUG")
    
    def get_full_table_tsv(self, results_dir):
        # BUSCO outputs a file named full_table_<run_name>.tsv in the results directory
        if not os.path.exists(results_dir):
            return None
        
        #Open the run_[lineage] folder (run_*)
        lineage_dir = glob.glob(os.path.join(results_dir, "run_*"))
        if not lineage_dir:
            return None
        results_dir = lineage_dir[0]

        tsv_files = []
        for fname in os.listdir(results_dir):
            if fname.startswith("full_table") and fname.endswith(".tsv"):
                tsv_files.append(os.path.join(results_dir, fname))

        if not tsv_files:
            self.log(f"BUSCO full table TSV file does not exist: {results_dir}/full_table*.tsv", "ERROR")
            return None
        if len(tsv_files) > 1:
            self.log(f"Multiple BUSCO full table TSV files found: {tsv_files}", "ERROR")
            return False

        return tsv_files[0]

    def _get_run_dir(self, results_dir: str) -> Optional[str]:
        if not results_dir or not os.path.exists(results_dir):
            return None
        run_dirs = sorted(glob.glob(os.path.join(results_dir, "run_*")))
        if not run_dirs:
            return None
        return run_dirs[0]

    def _status_sequence_dirs(self, status: int) -> list[str]:
        if status == 1:
            return ["single_copy_busco_sequences"]
        if status == 2:
            return ["multi_copy_busco_sequences"]
        if status == 3:
            return ["fragmented_busco_sequences"]
        return []

    def _sequence_kind_for_path(self, path: str) -> Optional[str]:
        low = str(path).lower()
        if low.endswith(".faa") or low.endswith(".faa.gz"):
            return "prot"
        if low.endswith(".fna") or low.endswith(".fna.gz"):
            return "nucl"
        return None

    def _find_busco_sequence_files(self, results_dir, status, family_id):
        run_dir = self._get_run_dir(results_dir)
        if not run_dir:
            return []
        found = []
        for subdir in self._status_sequence_dirs(int(status)):
            base = os.path.join(run_dir, "busco_sequences", subdir, f"{family_id}")
            for suffix in (".faa", ".faa.gz", ".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz"):
                candidate = f"{base}{suffix}"
                if os.path.exists(candidate):
                    found.append(candidate)
        return found

    def _get_busco_location(self, results_dir, status, family_id, preferred_sequence_kind: str = "prot"):
        files = self._find_busco_sequence_files(results_dir, status, family_id)
        if not files:
            return None
        preferred = next((path for path in files if self._sequence_kind_for_path(path) == preferred_sequence_kind), None)
        return preferred or files[0]

    def _register_busco_run_artifacts(self, run_id: Optional[int], results_dir: str, json_file: Optional[str], family_table: Optional[str]) -> None:
        if run_id is None:
            return
        try:
            self.db_manager.busco.register_run_artifact(run_id, "busco_result_root", results_dir, is_dir=True, format="directory")
            run_dir = self._get_run_dir(results_dir)
            if run_dir:
                self.db_manager.busco.register_run_artifact(run_id, "busco_run_dir", run_dir, is_dir=True, format="directory")
                seq_dir = os.path.join(run_dir, "busco_sequences")
                if os.path.isdir(seq_dir):
                    self.db_manager.busco.register_run_artifact(run_id, "busco_sequences_dir", seq_dir, is_dir=True, format="directory")
            if json_file and os.path.isfile(json_file):
                self.db_manager.busco.register_run_artifact(run_id, "busco_summary_json", json_file, format="json")
            if family_table and os.path.isfile(family_table):
                self.db_manager.busco.register_run_artifact(run_id, "busco_full_table_tsv", family_table, format="tsv")
        except Exception as exc:  # boundary: BUSCO run artifacts are helpful metadata, not required results
            self.log(f"Failed to register BUSCO run artifacts for run_id={run_id}: {exc}", "WARNING")

    def _cleanup_miniprot_ref_file(self, results_dir):
        if self.keep_miniprot_ref_file:
            return
        run_dirs = glob.glob(os.path.join(results_dir, "run_*"))
        if not run_dirs:
            return
        removed = 0
        for run_dir in run_dirs:
            ref_path = os.path.join(run_dir, "miniprot_output", "ref.mpi")
            if not os.path.isfile(ref_path):
                continue
            try:
                os.remove(ref_path)
                removed += 1
            except OSError as exc:
                self.log(f"Failed to remove miniprot ref file {ref_path}: {exc}", "WARNING")
        if removed:
            self.log(f"Removed {removed} miniprot ref file(s) from {results_dir}.", "DEBUG")

    def _rename_run_dir_with_pipeline(self, results_dir: str, pipeline_name: str, lineage: str, run_id: Optional[int]) -> None:
        run_dirs = sorted(glob.glob(os.path.join(results_dir, "run_*")))
        if not run_dirs:
            return
        src = run_dirs[0]
        target_base = f"run_{pipeline_name}_{lineage}"
        target = os.path.join(results_dir, target_base)
        if os.path.abspath(src) == os.path.abspath(target):
            return
        if os.path.exists(target):
            suffix = str(run_id) if run_id is not None else uuid.uuid4().hex[:8]
            target = os.path.join(results_dir, f"{target_base}__{suffix}")
        try:
            os.rename(src, target)
        except OSError as exc:
            self.log(f"Failed to rename BUSCO run folder {src} -> {target}: {exc}", "WARNING")


    def run(self):
        # Phase 10: Ensure BUSCO lineage library is available (download if needed)
        if not self.library:
            return self.handle_exception("Library is not specified.", {"library": self.library})
        if not self.lineage:
            return self.handle_exception("Lineage is not specified.", {"lineage": self.lineage})
        if not self.accession:
            return self.handle_exception("Accession is not specified.", {"accession": self.accession})
        self.accession = self.resolve_assembly_accession(self.accession)
        if self.format not in ("auto", "protein", "genome"):
            return self.handle_exception(ValueError("Format must be one of 'auto', 'protein', or 'genome' (or nucleotide aliases)."))

        # Resolve lineage directory root
        env = self.db_manager.env.get_many(["BUSCO_LINEAGE_DIR", "BUSCO_BINARIES_PATH"]) or {}
        lineages_root = env.get("BUSCO_LINEAGE_DIR")
        if not lineages_root:
            libs_dir = self.db_manager.storage.get_root_base("libraries")
            if not libs_dir:
                return self.handle_exception("Missing libraries storage root or BUSCO_LINEAGE_DIR configuration.", {})
            lineages_root = os.path.join(libs_dir, "lineages")
        lineage_dir = os.path.join(lineages_root, self.lineage)

        # Subtask to download -------------------------------

        def queue_busco_library_download():
            libraries_dir = self.db_manager.storage.get_root_base("libraries", fallback=os.path.dirname(lineages_root))
            busco_path = env.get("BUSCO_BINARIES_PATH")
            self.queue_subtask(
                job_type=6,
                status="P",
                priority=1,
                data={
                    "lineage": self.lineage,
                    "libraries_dir": libraries_dir,
                    "busco_path": busco_path,
                },
            )
            return True

        def busco_library_ready():
            return os.path.exists(lineage_dir)

        # If the lineage folder doesn't exist, queue a download subtask and suspend until ready
        if not busco_library_ready():
            outcome = self.manage_subtasks(
                stage=1,
                queue_fn=queue_busco_library_download,
                done_fn=busco_library_ready,
                wait_seconds=self.busco_lib_wait_seconds,
                retry_key=None,
                max_retries=self.busco_lib_retries,
                incomplete_message_fn=lambda: (f"BUSCO lineage '{self.lineage}' not yet available", ""),
                retry_incomplete=False,
            )
            if outcome == "ERROR":
                return "ERROR"
            if outcome is False:
                # Suspended; daemon will resume when subtask completes
                return False
            # continue
            
        # -------------------------------

        # Checkpoint 1: Get location of genome
        genome_path = self.db_manager.genomes.resolve_path(self.accession)
        if not self.output_path:
            self.output_path = genome_path
        run_suffix = uuid.uuid4().hex[:10]
        results_folder = f"{self.lineage}_results__{run_suffix}"
        results_dir = os.path.join(self.output_path, results_folder)
        legacy_results_dir = os.path.join(self.output_path, f"{self.lineage}_results")
        result = None
        if not genome_path:
            return self.handle_exception("No path to genome. Does the genome exist and is it downloaded?", {"accession": self.accession})

        try:
            effective_format = self._resolve_effective_format(genome_path)
        except (OSError, ValueError) as exc:
            return self.handle_exception(exc, {"format": self.format, "accession": self.accession})

        library_id = self.db_manager.libraries.get_id(self.library)
        if library_id is None:
            return self.handle_exception("BUSCO library is not registered.", {"library": self.library})
        if self.keep_miniprot_ref_file_raw is None:
            self.keep_miniprot_ref_file = bool(self.db_manager.env.get("BUSCO_MINIPROT_KEEP_REF_FILE"))
        pipeline_name = self._resolve_pipeline_name(effective_format)
        pipeline_args, pipeline_effective, pipeline_source = self._resolve_pipeline_args(pipeline_name)
        if effective_format != "genome":
            pipeline_args = []
        run_id = None
        created_new_run = False

        def fail_run() -> None:
            if run_id is None or not created_new_run:
                return
            try:
                self.db_manager.busco.update_run(run_id, status="failed", completed=True)
            except Exception as exc:  # boundary: preserve primary BUSCO failure while logging failed status update
                self.log(f"Failed to mark BUSCO run {run_id} as failed: {exc}", "WARNING")

        def create_failed_preflight_run(*, proteome_profile_id: Optional[int] = None) -> None:
            nonlocal run_id, created_new_run
            if run_id is not None:
                return
            run_id = self.db_manager.busco.create_run(
                accession=self.accession,
                library_id=library_id,
                lineage_name=self.lineage,
                input_mode=effective_format,
                pipeline=pipeline_name,
                pipeline_params_effective=pipeline_effective,
                pipeline_params_source=pipeline_source,
                busco_cli_args=[],
                result_dir=results_dir,
                proteome_profile_id=proteome_profile_id,
                status="running",
            )
            created_new_run = run_id is not None
            fail_run()

        selected_profile_name = None
        selected_profile_id = None
        proteome_input_path = None
        if effective_format == "protein":
            try:
                selected_profile_name, selected_profile_id, proteome_input_path = self._resolve_proteome_input(
                    self.accession,
                    genome_path,
                )
                selected_row = self.db_manager.proteomes.get(int(selected_profile_id))
                selected_checksum = selected_row[8] if selected_row and len(selected_row) > 8 else None
                if (
                    self.expected_proteome_profile_id is not None
                    and int(selected_profile_id) != int(self.expected_proteome_profile_id)
                ):
                    raise ValueError(
                        f"Pinned proteome profile changed for accession '{self.accession}': "
                        f"expected id {self.expected_proteome_profile_id}, observed {selected_profile_id}."
                    )
                if (
                    self.expected_proteome_checksum is not None
                    and str(selected_checksum) != str(self.expected_proteome_checksum)
                ):
                    raise ValueError(
                        f"Pinned proteome checksum changed for accession '{self.accession}': "
                        f"expected {self.expected_proteome_checksum}, observed {selected_checksum}."
                    )
            except Exception as exc:  # boundary: preflight profile resolution failure becomes this task's failure
                create_failed_preflight_run()
                return self.handle_exception(exc, {"accession": self.accession, "proteome_profile": self.proteome_profile})

        if self.force:
            if library_id is not None:
                # Purge stale BUSCO + paralog filtering rows before re-running.
                self.db_manager.busco.delete_records(self.accession, library_id)
                self.db_manager.filtering.delete_paralog_records(
                    self.accession,
                    busco_library_id=library_id,
                )
            else:
                self.log(
                    f"Force requested for {self.accession}, but library_id not found for '{self.library}'. Skipping DB purge.",
                    "WARNING",
                )
        
        # Decide whether we can skip execution for this specific lineage/pipeline.
        existing_pipeline_results_dir, existing_run_id = self._find_existing_pipeline_run(
            self.output_path,
            self.lineage,
            pipeline_name,
            accession=self.accession,
            library_id=library_id,
            input_mode=effective_format,
            proteome_profile_id=selected_profile_id,
        )
        should_run = bool(self.force or not existing_pipeline_results_dir)
        staged_busco_input_path = None

        if should_run:
            run_id = self.db_manager.busco.create_run(
                accession=self.accession,
                library_id=library_id,
                lineage_name=self.lineage,
                input_mode=effective_format,
                pipeline=pipeline_name,
                pipeline_params_effective=pipeline_effective,
                pipeline_params_source=pipeline_source,
                busco_cli_args=[],
                result_dir=results_dir,
                proteome_profile_id=selected_profile_id,
                status="running",
            )
            created_new_run = run_id is not None
        elif existing_run_id is not None:
            run_id = int(existing_run_id)
        else:
            run_id = self.db_manager.busco.create_run(
                accession=self.accession,
                library_id=library_id,
                lineage_name=self.lineage,
                input_mode=effective_format,
                pipeline=pipeline_name,
                pipeline_params_effective=pipeline_effective,
                pipeline_params_source=pipeline_source,
                busco_cli_args=[],
                result_dir=existing_pipeline_results_dir,
                proteome_profile_id=selected_profile_id,
                status="completed",
            )
            created_new_run = run_id is not None

        if should_run:
            self._cleanup_stale_staged_busco_inputs(genome_path)
            try:
                # Find the input sequence file based on requested format
                seq_file = None
                if effective_format == "protein":
                    if proteome_input_path and os.path.exists(proteome_input_path):
                        if str(proteome_input_path).endswith(".gz"):
                            extracted_path = os.path.join(genome_path, f".busco_input_{selected_profile_name or 'protein'}.faa")
                            self._cleanup_staged_busco_input(extracted_path)
                            self.log(f"Extracting {proteome_input_path} to {extracted_path}", "DEBUG")
                            try:
                                with gzip.open(proteome_input_path, 'rb') as f_in, open(extracted_path, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                    seq_file = extracted_path
                                    staged_busco_input_path = extracted_path
                            except (OSError, EOFError, gzip.BadGzipFile) as e:
                                fail_run()
                                return self.handle_exception(
                                    "Failed to extract compressed proteome profile.",
                                    {"file": proteome_input_path, "error": str(e)},
                                )
                        else:
                            seq_file = proteome_input_path
                else:
                    target_ext = ".fna"
                    for fname in os.listdir(genome_path):
                        if fname.endswith(target_ext):
                            seq_file = os.path.join(genome_path, fname)
                            break
                        gz_ext = f"{target_ext}.gz"
                        if fname.endswith(gz_ext):
                            gz_path = os.path.join(genome_path, fname)
                            extracted_path = os.path.join(genome_path, fname[:-3])  # strip .gz
                            self.log(f"Extracting {gz_path} to {extracted_path}", "DEBUG")
                            try:
                                with gzip.open(gz_path, 'rb') as f_in, open(extracted_path, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                    seq_file = extracted_path
                                    break
                            except (OSError, EOFError, gzip.BadGzipFile) as e:
                                fail_run()
                                return self.handle_exception("Failed to extract compressed genome file.", {"file": gz_path, "error": str(e)})

                if not seq_file or not os.path.exists(seq_file):
                    missing = ".fna/.fna.gz" if effective_format == "genome" else ".faa/.faa.gz"
                    fail_run()
                    return self.handle_exception(f"Genome {missing} file does not exist.", {"genome_path": genome_path})

                vars = self.db_manager.env.get_many(['BUSCO_BINARIES_PATH'])
                busco_binaries_path = vars.get('BUSCO_BINARIES_PATH')

                command = [
                    f"{busco_binaries_path}",
                    "-i", seq_file,
                    "-o", results_folder,
                    "--out_path", self.output_path,
                    "-l", f"{lineage_dir}",
                    "-m", effective_format,
                    "--offline",
                    "--force",
                    "-c", str(self.REQUIRED_THREADS)
                ]
                if pipeline_args:
                    command.extend(pipeline_args)
                if created_new_run and run_id is not None:
                    try:
                        with self.db_manager.transaction(operation=f"record BUSCO command for run {run_id}"):
                            self.db_manager.busco.update_run(run_id, result_dir=results_dir, status="running", completed=False)
                            self.db_manager.cursor.execute(
                                "UPDATE BUSCO_Runs SET busco_cli_args_json = ? WHERE run_id = ?",
                                (json.dumps(command), int(run_id)),
                            )
                    except Exception as exc:  # boundary: command metadata is required for a newly created BUSCO run
                        fail_run()
                        return self.handle_exception(
                            "Failed to record BUSCO command metadata.",
                            {"run_id": run_id, "accession": self.accession, "error": str(exc)},
                        )

                self.log(
                    f"Using command: {' '.join(command)}",
                    "DEBUG",
                )
                self.log(
                    f"Running BUSCO for {self.accession} "
                    f"(format={effective_format}, pipeline={pipeline_name}, lineage={self.lineage}).",
                    "INFO",
                )
                # Execute the command
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode != 0:
                    fail_run()
                    self.error(f"BUSCO failed: {result.stderr}")
                    return self.handle_exception("BUSCO command failed.", {"returncode": result.returncode, "stderr": result.stderr})
                self._rename_run_dir_with_pipeline(results_dir, pipeline_name, self.lineage, run_id)
            finally:
                self._cleanup_staged_busco_input(staged_busco_input_path)
        else:
            # We have skipped rerunning BUSCO
            results_dir = existing_pipeline_results_dir or legacy_results_dir
            self.log(f"Reusing existing BUSCO results for {self.accession} (pipeline={pipeline_name}) at {results_dir}.", "INFO")
            # self.checkpoint(2, {"done": True, "out_path": f"{self.output_path}/{self.lineage}_results", "stdout": result.stdout})

        # Add information about BUSCO run to the database

        # Find the BUSCO results JSON file
        json_file = self.get_json_file(results_dir)

        if json_file is None:
            fail_run()
            return self.handle_exception("BUSCO results JSON file not found.", {"results_dir": results_dir})
        if json_file is False:
            fail_run()
            return self.handle_exception("Multiple BUSCO results JSON files found.", {"results_dir": results_dir})

        try:
            with open(json_file, "r") as f:
                busco_results = json.load(f).get("results", {})
        except (OSError, UnicodeError, json.JSONDecodeError) as e:
            fail_run()
            return self.handle_exception("Failed to read BUSCO results JSON file.", {"error": str(e)})

        # If we did actually run BUSCO for the first time we will use the current date.
        if result:
            modified_time = None
        else:
            # Otherwise we will retrieve from the JSON file modified stat
            modified_time = os.path.getmtime(json_file)
            modified_time = datetime.fromtimestamp(modified_time).strftime("%Y-%m-%d %H:%M:%S")

        if not self.db_manager.busco.add_results(self.accession, library_id, busco_results, modified_time):
            fail_run()
            return self.handle_exception("Failed to add BUSCO results to database.", {"accession": self.accession, "library_id": library_id, "results": busco_results, "modified_time": modified_time})
        if run_id is not None:
            self.db_manager.busco.update_run(
                run_id,
                status="completed",
                result_dir=results_dir,
                counts=busco_results,
                proteome_profile_id=selected_profile_id,
                completed=True,
            )
            self.db_manager.busco.refresh_auto_primary_runs_for_accession(
                self.accession,
                library_id,
                updated_by="busco-task",
                policy="auto_best",
            )
            # Keep legacy folder path available for tasks that still resolve <lineage>_results.
            if os.path.abspath(results_dir) != os.path.abspath(legacy_results_dir):
                try:
                    if os.path.lexists(legacy_results_dir):
                        if os.path.islink(legacy_results_dir) or os.path.isfile(legacy_results_dir):
                            os.remove(legacy_results_dir)
                        else:
                            # Preserve pre-existing real results folder to avoid deleting prior pipeline outputs.
                            backup = f"{legacy_results_dir}__legacy_{datetime.now().strftime('%Y%m%d%H%M%S')}"
                            shutil.move(legacy_results_dir, backup)
                            self.log(
                                f"Moved existing legacy BUSCO results folder to {backup} before updating symlink.",
                                "WARNING",
                            )
                    if os.path.isdir(results_dir):
                        os.symlink(results_dir, legacy_results_dir)
                except OSError as exc:
                    self.log(f"Failed to update legacy BUSCO results symlink {legacy_results_dir}: {exc}", "WARNING")
        
        # Update genome status to indicate BUSCO complete (2)
        if not self.db_manager.genomes.update_status(self.accession, 2):
            fail_run()
            return self.handle_exception(f"Failed to update genome status for {self.accession}")
        
        # Now need to extract the per BUSCO family results and add to the BUSCO_Family_Data table
        family_table = self.get_full_table_tsv(results_dir)
        if family_table is None:
            fail_run()
            return self.handle_exception("BUSCO full table TSV file not found.", {"results_dir": results_dir})
        if family_table is False:
            fail_run()
            return self.handle_exception("Multiple BUSCO full table TSV files found.", {"results_dir": results_dir})
        self._register_busco_run_artifacts(run_id, results_dir, json_file, family_table)

        try:
            with open(family_table, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                family_data = []
                family_locations = []
                family_artifacts = []
                for row in reader:
                    if row[0].startswith("#"):
                        continue  # Skip header or comment lines
                    if len(row) < 2:
                        self.log(f"Skipping invalid BUSCO family data row: {row}", "WARNING")
                        continue
                    family_id = row[0]
                    status_str = row[1]
                    status_map = {
                        "Complete": 1,
                        "Duplicated": 2,
                        "Fragmented": 3,
                        "Missing": 4
                    }
                    status = status_map.get(status_str, 0)
                    if len(row) < 5:
                        # self.log(f"Skipping BUSCO family data row with insufficient columns: {row}", "WARNING")
                        continue
                    sequence = row[2] if row[2] != "" else None
                    # BUSCO full_table columns (v5/v6): Busco id, Status, Sequence, Gene Start, Gene End, Strand, Score, Length, OrthoDB url, Description
                    score = None
                    length = None
                    if len(row) >= 7 and row[6] != "":
                        try:
                            score = float(row[6])
                        except (TypeError, ValueError, OverflowError):
                            score = None
                    if len(row) >= 8 and row[7] != "":
                        try:
                            length = int(row[7])
                        except (TypeError, ValueError):
                            length = None
                    # Fallback for older/short rows that may have been parsed differently
                    if score is None and len(row) >= 4 and row[3] != "":
                        try:
                            score = float(row[3])
                        except (TypeError, ValueError, OverflowError):
                            score = None
                    if length is None and len(row) >= 5 and row[4] != "":
                        try:
                            length = int(row[4])
                        except (TypeError, ValueError):
                            length = None
                    # Detect probable mis-parse (coordinate mistaken as score) and drop it
                    if score is not None and score > 1e5:
                        self.log(
                            f"Discarding anomalous BUSCO score {score} for {family_id} ({self.accession}); likely a coordinate column.",
                            "WARNING",
                        )
                        score = None
                    location = self._get_busco_location(results_dir, status, family_id)
                    family_data.append((family_id, library_id, self.accession, status, sequence, score, length))
                    family_locations.append((family_id, library_id, self.accession, location))
                    if run_id is not None:
                        for seq_path in self._find_busco_sequence_files(results_dir, status, family_id):
                            family_artifacts.append(
                                (
                                    family_id,
                                    seq_path,
                                    self._sequence_kind_for_path(seq_path),
                                    {
                                        "status": status,
                                        "accession": self.accession,
                                        "library_id": library_id,
                                    },
                                )
                            )
                if run_id is not None:
                    if not self.db_manager.busco.add_run_family_data(run_id, family_data):
                        fail_run()
                        return self.handle_exception("Failed to add BUSCO run family data.", {"accession": self.accession, "library_id": library_id, "run_id": run_id})
                    if not self.db_manager.busco.add_run_family_locations(run_id, list(set(family_locations))):
                        fail_run()
                        return self.handle_exception("Failed to add BUSCO run family locations.", {"accession": self.accession, "library_id": library_id, "run_id": run_id})
                    for family_id, seq_path, sequence_kind, metadata in family_artifacts:
                        self.db_manager.busco.register_family_artifact(
                            run_id=run_id,
                            family_id=family_id,
                            library_id=library_id,
                            accession=self.accession,
                            path=seq_path,
                            sequence_kind=sequence_kind,
                            format="fasta",
                            metadata=metadata,
                        )
        except Exception as e:  # boundary: BUSCO family ingestion failure is converted to this task error
            fail_run()
            return self.handle_exception("Failed to read BUSCO family data TSV file.", {"error": str(e), "file": family_table})
        
        summary_fields = {
            "Complete": busco_results.get("Complete percentage"),
            "Single": busco_results.get("Single copy percentage"),
            "Duplicated": busco_results.get("Multi copy percentage"),
            "Fragmented": busco_results.get("Fragmented percentage"),
            "Missing": busco_results.get("Missing percentage"),
        }
        summary_text = ", ".join(
            f"{label}={value}"
            for label, value in summary_fields.items()
            if value is not None
        )
        if summary_text:
            self.log(f"BUSCO completed for {self.accession}: {summary_text}.", "INFO")
        else:
            self.log(f"BUSCO completed for {self.accession}.", "INFO")
        self.log(f"BUSCO results: {busco_results}", "DEBUG")

        self._cleanup_miniprot_ref_file(results_dir)

        return True
