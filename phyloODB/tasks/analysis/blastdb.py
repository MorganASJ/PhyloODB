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
from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE, RAW_PROFILE, resolve_profile_selector
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
)


def _cache_blastdb_prefix(task: Task, *parts: str) -> str:
    cache_root = task.db_manager.storage.ensure_cache_root()
    path = os.path.join(cache_root, *[str(part) for part in parts])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def _sanitize_profile_label(label: Optional[str]) -> str:
    token = str(label or "").strip() or RAW_PROFILE
    return re.sub(r"[^A-Za-z0-9._-]+", "_", token).strip("._-") or RAW_PROFILE


def resolve_proteome_profile_input(
    manager: DBManager,
    accession: str,
    *,
    proteome_profile: Optional[str] = None,
    prefer_proteome_profile: Optional[str] = None,
) -> tuple[str, object, str]:
    accession_token = str(accession)
    requested_profile = str(proteome_profile or "").strip() or None
    preferred_profile = str(prefer_proteome_profile or "").strip() or None

    if requested_profile:
        resolved = manager.proteomes.resolve_selector_profile_name(accession_token, requested_profile)
        if resolved:
            requested_profile = resolved
        elif requested_profile == DEFAULT_CLEAN_PROFILE:
            default_cleaned = manager.proteomes.get_default_cleaned_profile_name(accession_token)
            if default_cleaned:
                requested_profile = default_cleaned
        profile_name = requested_profile
    else:
        profile_name = None
        if preferred_profile:
            preferred_resolved = manager.proteomes.resolve_selector_profile_name(accession_token, preferred_profile) or preferred_profile
            preferred_row = manager.proteomes.get_profile(accession_token, preferred_resolved)
            if preferred_row is not None:
                profile_name = preferred_resolved
        if profile_name is None:
            default_profile = manager.proteomes.get_default_profile_name(accession_token)
            if default_profile:
                profile_name = default_profile
            else:
                default_cleaned = manager.proteomes.get_default_cleaned_profile_name(accession_token)
                if default_cleaned:
                    profile_name = default_cleaned
        if profile_name is None:
            profile_name = RAW_PROFILE

    if profile_name == RAW_PROFILE:
        raw_id = manager.proteomes.ensure_raw_profile(accession_token, is_default=False)
        if raw_id is None:
            raise FileNotFoundError(f"Proteome profile '{profile_name}' does not exist for accession '{accession_token}'.")
        row = manager.proteomes.get(int(raw_id))
    else:
        row = manager.proteomes.get_profile(accession_token, profile_name)
    if row is None:
        raise FileNotFoundError(f"Proteome profile '{profile_name}' does not exist for accession '{accession_token}'.")
    path = manager.proteomes.resolve_path(row)
    if not path or not os.path.exists(path):
        raise FileNotFoundError(
            f"Proteome profile '{profile_name}' has no readable artifact for accession '{accession_token}'."
        )
    return str(profile_name), row, str(path)


def build_proteome_blastdb_prefix(task: Task, accession: str, library_id: Optional[int], selected_profile: Optional[str]) -> str:
    library_token = f"library_{library_id}" if library_id is not None else "default"
    acc_token = str(accession or "proteome").strip() or "proteome"
    profile_token = _sanitize_profile_label(selected_profile)
    return _cache_blastdb_prefix(task, "proteome-blastdb", library_token, f"{acc_token}_{profile_token}_blastdb")

class CreateProteomeBlastDB(Task):
    '''A class that handles the creation of a proteome BLAST database'''
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=4):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.accession = self.data.get("accession", None)
        self.faa = self.data.get("faa", None) # The faa file that will become the db
        self.out_path = self.data.get("out_path", None) # the path that the db will be created at
        self.force = self.data.get("force", False)
        self.library_name = self.data.get("library_name", None)
        self.library_id = self.data.get("library_id", None)
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        self.selected_profile: Optional[str] = None

    def _materialize_faa_path(self, path: str) -> str:
        if not str(path).lower().endswith(".gz"):
            return str(path)
        decompressed_path = str(path)[:-3]
        if os.path.exists(decompressed_path):
            return decompressed_path
        try:
            with gzip.open(path, "rb") as f_in, open(decompressed_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        except (OSError, EOFError, gzip.BadGzipFile) as exc:
            raise RuntimeError(str(exc)) from exc
        return decompressed_path
        
    def run(self):
        self.accession = self.resolve_assembly_accession(self.accession)
        
        # if a direct path is not given to a proteome check the accession has one and get it
        if not self.faa:
            genome_info = self.db_manager.genomes.get(self.accession)
            if not genome_info:
                return self.handle_exception(f"Accession {self.accession} not found in database.", {"accession": self.accession})
            genome_col_map = self.db_manager.genomes.get_column_map()
            location_idx = genome_col_map.get("location")
            if location_idx is None:
                return self.handle_exception(
                    "Genome table is missing location column.",
                    {"columns": list(genome_col_map.keys())},
                )
            genome_dir = genome_info[location_idx]
            if not genome_dir:
                return self.handle_exception(
                    f"Accession {self.accession} has no location information available.",
                    {"accession": self.accession},
                )
            try:
                selected_profile, _profile_row, profile_path = resolve_proteome_profile_input(
                    self.db_manager,
                    str(self.accession),
                    proteome_profile=self.proteome_profile,
                    prefer_proteome_profile=self.prefer_proteome_profile,
                )
                self.selected_profile = str(selected_profile)
                self.faa = self._materialize_faa_path(str(profile_path))
            except Exception as exc:  # boundary: convert missing/invalid required proteome input into task error state.
                return self.handle_exception(
                    f"Accession {self.accession} has no protein information available for the requested proteome profile.",
                    {
                        "accession": self.accession,
                        "proteome_profile": self.proteome_profile,
                        "prefer_proteome_profile": self.prefer_proteome_profile,
                        "error": str(exc),
                    },
                )
        else:
            self.faa = self._materialize_faa_path(str(self.faa))

        if not os.path.exists(self.faa):
            return self.handle_exception(f"FASTA file for accession {self.accession} does not exist.", {"faa": self.faa})

        if self.library_name and not self.library_id:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
            if not self.library_id:
                return self.handle_exception(f"Library '{self.library_name}' not found in database.", {"library_name": self.library_name})

        # If the output path is not given (expected) then create one based on the faa file location
        if not self.out_path:
            if not self.selected_profile and self.accession:
                try:
                    self.selected_profile, _profile_row, _profile_path = resolve_proteome_profile_input(
                        self.db_manager,
                        str(self.accession),
                        proteome_profile=self.proteome_profile,
                        prefer_proteome_profile=self.prefer_proteome_profile,
                    )
                except (FileNotFoundError, OSError, ValueError):
                    self.selected_profile = os.path.splitext(os.path.basename(self.faa))[0]
            self.out_path = build_proteome_blastdb_prefix(self, str(self.accession), self.library_id, self.selected_profile)
        else:
            out_dir = os.path.dirname(self.out_path)
            if not os.path.exists(out_dir):
                return self.handle_exception("Output directory for BLAST database does not exist.", {"out_dir": out_dir})

        if not self.force:
            # Check if BLAST DB already exists for this accession
            rows = self.db_manager.filtering.get_blast_dbs(accession=self.accession, library_id=self.library_id)
            expected_files = [f"{self.out_path}.{ext}" for ext in ["phr", "pin", "psq"]]
            for row in rows or []:
                row_location = str(row[3]) if row and len(row) > 3 and row[3] is not None else None
                if row_location == self.out_path and all(os.path.exists(path) for path in expected_files):
                    self.log(
                        f"BLAST database already exists for accession {self.accession} using proteome profile '{self.selected_profile or RAW_PROFILE}', skipping creation."
                    )
                    return True

            # If out_path is given, check if the BLAST DB files already exist there
            if self.out_path:
                db_files = [f"{self.out_path}.{ext}" for ext in ["phr", "pin", "psq"]]
                if all(os.path.exists(f) for f in db_files):
                    self.log(f"BLAST database files already exist at {self.out_path}, skipping creation.")
                    return True

        # Create BLAST database using makeblastdb
        makeblastdb_path = self.db_manager.env.get('MAKEBLASTDB_PATH')
        if not makeblastdb_path:
            return self.handle_exception("makeblastdb path is not set in environment variables.", {"MAKEBLASTDB_PATH": makeblastdb_path})

        command = [
            makeblastdb_path,
            "-in", self.faa,
            "-dbtype", "prot",
            "-out", self.out_path
        ]

        self.log(f"Creating BLAST database for accession {self.accession} with command: {' '.join(command)}", "DEBUG")

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode != 0:
            self.error(f"makeblastdb failed: {result.stderr}")
            return self.handle_exception("makeblastdb command failed.", {"returncode": result.returncode, "stderr": result.stderr})

        # Record the BLAST database in the database
        blastdb_location = self.out_path
        blastdb_id = self.db_manager.filtering.add_proteome_blastdb(
            accession=self.accession,
            location=blastdb_location,
            library_id=self.library_id
        )
        if not blastdb_id:
            return self.handle_exception("Failed to record BLAST database in database.", {"accession": self.accession, "location": blastdb_location})

        self.db_manager.artifacts.register(
            owner_type="blast_db",
            owner_id=int(blastdb_id),
            artifact_type="blast_db_prefix",
            path=self.out_path,
            format="prot",
            metadata={"source": "proteome", "accession": self.accession, "library_id": self.library_id},
        )

        self.log(f"BLAST database created and recorded for accession {self.accession} at {blastdb_location}", "INFO")
        return True

class ConstructBuscoBlastDB(Task):
    """Build a BLAST DB from BUSCO sequences for a set of accessions."""
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=4):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.accessions = list(dict.fromkeys(normalize_accessions(self.data.get("accessions", []))))
        self.busco_library_id = self.data.get("busco_library_id")
        self.target_library_id = self.data.get("target_library_id")
        self.family_ids = list(dict.fromkeys(self.data.get("family_ids", [])))
        self.output_path = self.data.get("output_path")
        self.use_paralog_filtered = bool(self.data.get("use_paralog_filtered", False))
        self.force = bool(self.data.get("force", False))
        self.id_mode = (self.data.get("id_mode") or "legacy").lower()
        self.id_map_path = self.data.get("id_map_path")
        self.db_type = (self.data.get("db_type") or "").lower() or None
        self.deduplicate_family_per_accession = bool(self.data.get("deduplicate_family_per_accession", True))

    def run(self):
        if not self.accessions:
            return self.handle_exception("No accessions provided to build BUSCO BLAST DB.", {"accessions": self.accessions})
        if not self.busco_library_id:
            return self.handle_exception("busco_library_id is required for BUSCO sequence lookup.", {"busco_library_id": self.busco_library_id})

        if not self.output_path:
            target_token = f"target_{self.target_library_id}" if self.target_library_id is not None else "unscoped"
            busco_token = f"busco_{self.busco_library_id}"
            digest = hashlib.sha1(
                "\t".join(
                    [
                        target_token,
                        busco_token,
                        ",".join(sorted(self.accessions)),
                        ",".join(sorted(str(f) for f in self.family_ids)),
                        str(self.use_paralog_filtered),
                        str(self.id_mode),
                        str(self.db_type or ""),
                    ]
                ).encode("utf-8")
            ).hexdigest()[:16]
            self.output_path = _cache_blastdb_prefix(self, "busco-blastdb", target_token, busco_token, digest, "ref_buscos")
        out_dir = os.path.dirname(self.output_path)
        try:
            os.makedirs(out_dir, exist_ok=True)
        except OSError as exc:
            return self.handle_exception("Failed to create output directory.", {"out_dir": out_dir, "error": str(exc)})

        fasta_path = f"{self.output_path}.faa"
        makeblastdb_path = self.db_manager.env.get('MAKEBLASTDB_PATH')
        if not makeblastdb_path:
            return self.handle_exception("MAKEBLASTDB_PATH not configured.", {})

        clean_lookup = {}
        if self.use_paralog_filtered:
            try:
                rows = self.db_manager.filtering.get_paralog_results(target_library_id=self.target_library_id)
                for fam, lib_id, targ_lib, acc, clean in rows:
                    if lib_id == self.busco_library_id and clean:
                        clean_lookup.setdefault(acc, set()).add(fam)
            except Exception as exc:  # boundary: optional paralog filter can be omitted without invalidating BLAST DB construction.
                self.log("Paralog filter lookup failed; proceeding without filter.", "WARNING")
                clean_lookup = {}

        family_filter = set(self.family_ids) if self.family_ids else None
        id_map_path = self.id_map_path
        if self.id_mode != "legacy" and not id_map_path:
            id_map_path = f"{self.output_path}.id_map.tsv"
        if self.deduplicate_family_per_accession:
            self.log("BUSCO BLAST DB: deduplicating to one sequence per accession/family by highest BUSCO bitscore.", "INFO")
        acc_taxids = {}
        for acc in self.accessions:
            try:
                g = self.db_manager.genomes.get(acc)
                acc_taxids[acc] = g[1] if g and len(g) > 1 else None
            except Exception as exc:  # boundary: optional taxonomy labels may be absent without invalidating BLAST DB construction.
                self.log(f"Taxid lookup failed for {acc}; using empty taxid in ID map: {exc}", "WARNING")
                acc_taxids[acc] = None
        map_rows = []
        detected_db_type = None
        written = 0

        def _stable_internal_seq_id(accession, family_id, taxid, orig_header, sequence):
            taxid_label = taxid if taxid is not None else "NA"
            # Stable per biological sequence and BUSCO context, independent of write order.
            digest_source = f"{accession}\t{family_id}\t{orig_header or ''}\t{sequence or ''}"
            digest = hashlib.sha1(digest_source.encode("utf-8")).hexdigest()[:16]
            acc_label = re.sub(r"[^A-Za-z0-9_.-]", "_", str(accession))
            return f"t{taxid_label}|f{family_id}|a{acc_label}|h{digest}"

        def _iter_fasta_records(path):
            header = None
            seq_lines = []
            with open(path) as fin:
                for raw in fin:
                    line = raw.rstrip("\n")
                    if line.startswith(">"):
                        if header is not None:
                            yield header, "".join(seq_lines), seq_lines
                        header = line[1:].strip()
                        seq_lines = []
                    else:
                        seq_lines.append(line.strip())
                if header is not None:
                    yield header, "".join(seq_lines), seq_lines

        with open(fasta_path, "w") as fout:
            for acc in sorted(self.accessions):
                rows = self.db_manager.busco.get_family_results_for_library(
                    library_id=self.busco_library_id,
                    accessions=[acc],
                    status=[1, 2],
                )
                if not rows:
                    self.log(f"No BUSCO family rows for {acc} (busco_library_id={self.busco_library_id})", "WARNING")
                    continue
                allowed_fams = clean_lookup.get(acc) if clean_lookup else None
                selected_rows = []
                if self.deduplicate_family_per_accession:
                    # Keep one BUSCO sequence per (accession, family), selecting the row
                    # with the highest BUSCO bitscore (length as deterministic tie-breaker).
                    best_by_family = {}
                    for row in rows:
                        family_id, lib_id, accession, status, seq_id, bitscore, length = row
                        if family_filter and family_id not in family_filter:
                            continue
                        if allowed_fams is not None and family_id not in allowed_fams:
                            continue
                        try:
                            score = float(bitscore) if bitscore is not None else float("-inf")
                        except (TypeError, ValueError):
                            score = float("-inf")
                        try:
                            seq_len = int(length) if length is not None else -1
                        except (TypeError, ValueError):
                            seq_len = -1
                        existing = best_by_family.get(family_id)
                        if existing is None:
                            best_by_family[family_id] = (row, score, seq_len)
                            continue
                        _prev_row, prev_score, prev_len = existing
                        if score > prev_score or (score == prev_score and seq_len > prev_len):
                            best_by_family[family_id] = (row, score, seq_len)
                    selected_rows = [entry[0] for entry in best_by_family.values()]
                else:
                    selected_rows = rows
                selected_rows = sorted(selected_rows, key=lambda r: (str(r[0]), str(r[1]), str(r[2]), str(r[4])))

                for family_id, lib_id, accession, status, seq_id, bitscore, length in selected_rows:
                    if family_filter and family_id not in family_filter:
                        continue
                    if allowed_fams is not None and family_id not in allowed_fams:
                        continue
                    requested_kind = self.db_type if self.db_type in ("prot", "nucl") else None
                    loc = self.db_manager.busco.get_family_location(
                        family_id,
                        lib_id,
                        accession,
                        sequence_kind=requested_kind,
                    )
                    if not loc or not os.path.exists(loc):
                        continue
                    try:
                        for orig_header, seq_str, seq_lines in _iter_fasta_records(loc):
                            if self.id_mode == "legacy":
                                seq_id = f"{accession}|{family_id}"
                            else:
                                taxid = acc_taxids.get(accession)
                                seq_id = _stable_internal_seq_id(accession, family_id, taxid, orig_header, seq_str)
                            fout.write(f">{seq_id}\n")
                            for seq_line in seq_lines:
                                fout.write(f"{seq_line}\n")
                            if id_map_path:
                                map_rows.append(
                                    (seq_id, accession, family_id, acc_taxids.get(accession), orig_header)
                                )
                        written += 1
                        if detected_db_type is None:
                            if loc.endswith(".fna") or loc.endswith(".fna.gz"):
                                detected_db_type = "nucl"
                            else:
                                detected_db_type = "prot"
                    except OSError as exc:
                        self.log(f"Failed to read BUSCO sequence {loc}: {exc}", "WARNING")

        if written == 0:
            return self.handle_exception("No BUSCO sequences written for BLAST DB.", {"accessions": self.accessions})

        db_type = self.db_type or detected_db_type or "prot"
        if db_type not in ("prot", "nucl"):
            db_type = "prot"
        cmd = [makeblastdb_path, "-in", fasta_path, "-dbtype", db_type, "-out", self.output_path]
        self.log(f"Building BUSCO BLAST DB at {self.output_path} from {written} sequences.", "INFO")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                return self.handle_exception("makeblastdb failed.", {"cmd": " ".join(cmd), "stderr": result.stderr})
        except (OSError, subprocess.SubprocessError) as exc:
            return self.handle_exception("Failed to run makeblastdb.", {"error": str(exc)})

        if id_map_path and map_rows:
            try:
                with open(id_map_path, "w") as fh:
                    fh.write("\t".join(["sseqid", "accession", "family_id", "taxid", "orig_header"]) + "\n")
                    for row in map_rows:
                        fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")
            except (OSError, UnicodeError) as exc:
                return self.handle_exception("Failed to write BUSCO BLAST ID map.", {"path": id_map_path, "error": str(exc)})

        self.data["blastdb_path"] = self.output_path
        if id_map_path:
            self.data["id_map_path"] = id_map_path
        self.data["db_type"] = db_type
        try:
            self.db_manager.tasks.update_data(self.task_id, data=self.data)
        except Exception as exc:  # boundary: persist enrichment only; produced BLAST DB remains usable.
            self.log(f"Failed to persist BUSCO BLAST DB task metadata: {exc}", "WARNING")
        self.db_manager.artifacts.register(
            owner_type="blast_db",
            owner_id=self.task_id,
            artifact_type="blast_db_prefix",
            path=self.output_path,
            format=db_type,
            metadata={"id_map_path": id_map_path, "source": "busco"},
        )
        if id_map_path:
            self.db_manager.artifacts.register(
                owner_type="blast_db",
                owner_id=self.task_id,
                artifact_type="blast_db_id_map",
                path=id_map_path,
                format="tsv",
                metadata={"db_type": db_type},
            )
        self.log(f"BUSCO BLAST DB created at {self.output_path}", "INFO")
        return True
