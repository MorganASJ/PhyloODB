import glob
import os
import sys
import time
import random
import threading
import shutil
import gzip
import json
import hashlib
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from ..task import Task
from ..utilities import NCBIHelper, FnaDownloadError, FaaDownloadError, GffDownloadError


PHYLOODB_MD5_MANIFEST = "phyloodb_md5checksums.txt"
NCBI_MD5_MANIFEST = "md5checksums.txt"


def _parse_ncbi_md5_file(path: str) -> Dict[str, str]:
    checksums: Dict[str, str] = {}
    if not path or not os.path.isfile(path):
        return checksums
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            digest = parts[0].strip().lower()
            filename = os.path.basename(parts[-1].lstrip("./"))
            if len(digest) == 32 and all(ch in "0123456789abcdef" for ch in digest) and filename:
                checksums[filename] = digest
    return checksums


def _md5_file(path: str) -> Optional[str]:
    if not path or not os.path.isfile(path):
        return None
    hasher = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _validate_gzip_file(path: str) -> None:
    with gzip.open(path, "rb") as handle:
        for _chunk in iter(lambda: handle.read(1024 * 1024), b""):
            pass


def _write_phyloodb_md5_manifest(folder: str, paths: Iterable[Optional[str]]) -> Optional[str]:
    """Write a PhyloODB checksum manifest in the same line format as NCBI."""

    if not folder:
        return None
    os.makedirs(folder, exist_ok=True)
    entries = []
    seen = set()
    for raw_path in paths:
        if not raw_path:
            continue
        path = os.path.abspath(str(raw_path))
        if path in seen or not os.path.isfile(path):
            continue
        seen.add(path)
        if path.lower().endswith(".gz"):
            _validate_gzip_file(path)
        entries.append((os.path.basename(path), _md5_file(path)))
    if not entries:
        return None
    manifest_path = os.path.join(folder, PHYLOODB_MD5_MANIFEST)
    with open(manifest_path, "w", encoding="utf-8") as handle:
        for filename, digest in sorted(entries):
            if digest:
                handle.write(f"{digest}  ./{filename}\n")
    return manifest_path


def _register_genome_file_artifacts(task: Task, accession: str, folder: str) -> None:
    if not folder or not os.path.isdir(folder):
        return
    existing_rows = task.db_manager.artifacts.find(owner_type="genome", owner_id=accession)
    existing_paths = set()
    for row in existing_rows:
        if not row or len(row) < 4:
            continue
        resolved = task.db_manager.artifacts.resolve_path(row)
        if resolved:
            existing_paths.add((str(row[3]), os.path.abspath(str(resolved))))
    preferred = {
        "genome_fna": None,
        "genome_faa": None,
        "genome_gff": None,
        "genome_faa_archive": None,
        "genome_ncbi_md5checksums": None,
        "genome_phyloodb_md5checksums": None,
    }
    for fname in sorted(os.listdir(folder)):
        path = os.path.join(folder, fname)
        low = fname.lower()
        if preferred["genome_ncbi_md5checksums"] is None and low == NCBI_MD5_MANIFEST:
            preferred["genome_ncbi_md5checksums"] = path
        elif preferred["genome_phyloodb_md5checksums"] is None and low == PHYLOODB_MD5_MANIFEST:
            preferred["genome_phyloodb_md5checksums"] = path
        elif preferred["genome_fna"] is None and low.endswith((".fna", ".fna.gz")):
            preferred["genome_fna"] = path
        elif preferred["genome_faa"] is None and low.endswith((".faa", ".faa.gz")) and ".archive" not in low:
            preferred["genome_faa"] = path
        elif preferred["genome_gff"] is None and low.endswith((".gff", ".gff.gz", ".gff3", ".gff3.gz")):
            preferred["genome_gff"] = path
        elif preferred["genome_faa_archive"] is None and low.endswith((".faa.archive", ".faa.archive.gz")):
            preferred["genome_faa_archive"] = path
    md5_manifest = preferred.get("genome_ncbi_md5checksums") or preferred.get("genome_phyloodb_md5checksums")
    expected_md5 = _parse_ncbi_md5_file(md5_manifest) if md5_manifest else {}
    checksum_source = "NCBI md5checksums.txt" if preferred.get("genome_ncbi_md5checksums") else "PhyloODB local import checksum manifest"
    for artifact_type, path in preferred.items():
        if path and os.path.exists(path):
            abs_path = os.path.abspath(path)
            if (artifact_type, abs_path) in existing_paths:
                continue
            basename = os.path.basename(path)
            checksum = (
                _md5_file(path)
                if artifact_type in {"genome_ncbi_md5checksums", "genome_phyloodb_md5checksums"}
                else expected_md5.get(basename)
            )
            metadata = {"accession": accession}
            if artifact_type == "genome_ncbi_md5checksums":
                metadata["source"] = "NCBI assembly FTP md5checksums.txt"
            elif artifact_type == "genome_phyloodb_md5checksums":
                metadata["source"] = "PhyloODB local import checksum manifest"
            elif checksum:
                metadata["checksum_source"] = checksum_source
            task.db_manager.artifacts.register(
                owner_type="genome",
                owner_id=accession,
                artifact_type=artifact_type,
                path=path,
                format="directory" if os.path.isdir(path) else os.path.splitext(path)[1].lstrip(".") or None,
                checksum=checksum,
                size_bytes=os.path.getsize(path) if os.path.isfile(path) else None,
                metadata=metadata,
            )


class UpdateAssemblyInformation(Task):
    '''A class that adds the metadata for a set of assemblies to the assembly table in the database.'''
    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads)
        from ...selector_utils import normalize_accessions

        self.taxid = self.data.get("taxid")
        self.released_after = self.data.get("after")
        self.released_before = self.data.get("before")
        self.level_filter = self.data.get("level")
        self.debug_path = self.data.get("debug_path")
        if self.taxid is not None:
            try:
                self.taxid = int(self.taxid)
            except (TypeError, ValueError) as exc:
                raise ValueError("'taxid' must be an integer if provided.") from exc
        raw_accessions = normalize_accessions(self.data.get("accessions") or [])
        if not self.taxid and not raw_accessions:
            raise ValueError("Data must contain either a 'taxid' or 'accessions' key.")
        if self.taxid and raw_accessions:
            raise ValueError("Data cannot contain both 'taxid' and 'accessions' keys. Please choose one.")
        self.accessions: List[str] = raw_accessions
        self.force_update = bool(self.data.get("force_update", False))

    # @staticmethod
    # def get_interface_prompts():
    
    # @staticmethod
    # def validate_data(data):
    #     if not isinstance(data, dict):
    #         raise ValueError("Data must be a dictionary.")
    #     if "taxid" not in data and "accessions" not in data:
    #         raise ValueError("Data must contain either a 'taxid' or 'accessions' key.")
    #     if "taxid" in data and "accessions" in data:Eoraptor676$
    
    #         raise ValueError("Data cannot contain both 'taxid' and 'accessions' keys. Please choose one.")
    #     if "taxid" in data and not isinstance(data["taxid"], int):
    #         raise ValueError("'taxid' must be an integer.")
    #     if "accessions" in data:
    #         if not isinstance(data["accessions"], list) or not all(isinstance(acc, str) for acc in data["accessions"]):
    #             raise ValueError("'accessions' must be a list of strings.")
    #     return True

    def run(self):
        '''Runs the task'''
        email = self.db_manager.env.get('EMAIL')
        if not email:
            return self.handle_exception("Email environment variable not set in the database.", {})
        api_key = self.db_manager.env.get('NCBI_API_KEY')
        # Initialize the NCBI helper
        ncbi_helper = NCBIHelper(email, db_manager=self.db_manager, api_key=api_key)
        # Expand any accession variables now that the DB connection is available.
        if self.accessions:
            try:
                self.accessions = self.selector_accessions()
            except ValueError as exc:
                return self.handle_exception(f"Failed to resolve accessions: {exc}", {"accessions": self.accessions})

        progress_state = {"last_step": -1}

        def _log_progress(done, total):
            if not total:
                return
            if total < 20:
                if done >= total and progress_state["last_step"] != 100:
                    self.log(f"Fetched metadata: 100% ({done}/{total} assemblies)", "INFO")
                    progress_state["last_step"] = 100
                return
            pct = int((done / total) * 100)
            step = (pct // 10) * 10
            if step >= 10 and step != progress_state["last_step"]:
                self.log(f"Fetched metadata: {step}% ({done}/{total} assemblies)", "INFO")
                progress_state["last_step"] = step

        debug_timestamp = None
        raw_debug_path = None
        full_debug_path = None
        if self.debug_path:
            try:
                os.makedirs(self.debug_path, exist_ok=True)
                debug_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                raw_debug_path = os.path.join(
                    self.debug_path,
                    f"assembly_esummary_raw_{self.task_id}_{debug_timestamp}.jsonl",
                )
                full_debug_path = os.path.join(
                    self.debug_path,
                    f"assembly_esummary_full_{self.task_id}_{debug_timestamp}.jsonl",
                )
            except OSError as exc:
                self.log(f"Failed to prepare debug path {self.debug_path}: {exc}", "WARNING")
                raw_debug_path = None
                full_debug_path = None

        hidden_rows = []
        if self.accessions:
            placeholders = ",".join("?" for _ in self.accessions)
            self.db_manager.cursor.execute(
                f"SELECT accession, status, reason FROM Hidden_Genomes WHERE accession IN ({placeholders})",
                tuple(self.accessions),
            )
            hidden_rows = self.db_manager.cursor.fetchall() or []

        if self.taxid:
            self.log(f"Fetching assembly metadata for taxid {self.taxid}.", "INFO")
            self.log(f"TaxID provided: {self.taxid}", "DEBUG")
            # Get the assemblies
            assembly_dataset, tax_info_dataset, taxonomy_update_dict, genome_dataset = ncbi_helper.fetch_assemblies_v2(
                self.taxid,
                progress_cb=_log_progress,
                debug_raw_path=raw_debug_path,
                debug_full_path=full_debug_path,
            )
        else:
            self.log(f"Fetching assembly metadata for {len(self.accessions)} requested accessions.", "INFO")
            # downloaded_only_flag = bool(self.data.get("downloaded_only", False))
            # released_after = self.data.get("after")
            # released_before = self.data.get("before")
            # level_filter = self.data.get("level")
            # protein_only_flag = bool(self.data.get("protein_only", False))
            # try:
            #     self.accessions = self.prepare_selectors(
            #         taxid=None,
            #         use_rule_selection=False,
            #         require_candidates=True,
            #         downloaded_only=downloaded_only_flag,
            #         released_after=released_after,
            #         released_before=released_before,
            #         level=level_filter,
            #         protein_only=protein_only_flag,
            #     )
            # except ValueError as exc:
            #     return self.handle_exception(
            #         str(exc),
            #         {
            #             "selectors": self.selector_accessions(),
            #             "after": released_after,
            #             "before": released_before,
            #             "level": level_filter,
            #             "protein_only": protein_only_flag,
            #         },
            #     )
            # self.log(f"Accessions provided: {len(self.accessions)}")
            # Get the assemblies
            assembly_dataset, tax_info_dataset, taxonomy_update_dict, genome_dataset = ncbi_helper.fetch_assemblies_v2(
                self.accessions,
                accessions_only=True,
                progress_cb=_log_progress,
                debug_raw_path=raw_debug_path,
                debug_full_path=full_debug_path,
            )

        def _passes_filters(asm, genome):
            """Filter by release date and assembly level if provided."""
            rel = asm.get("release_date")
            if self.released_after:
                if not rel or str(rel) < str(self.released_after):
                    return False
            if self.released_before:
                if not rel or str(rel) > str(self.released_before):
                    return False
            if self.level_filter:
                lvl = (genome.get("assembly_level") or "").strip().lower()
                if lvl != str(self.level_filter).strip().lower():
                    return False
            return True

        if assembly_dataset:
            filtered = [
                (asm, tax, genome)
                for asm, tax, genome in zip(assembly_dataset, tax_info_dataset, genome_dataset)
                if _passes_filters(asm, genome)
            ]
            if filtered:
                assembly_dataset, tax_info_dataset, genome_dataset = map(list, zip(*filtered))
            else:
                assembly_dataset, tax_info_dataset, genome_dataset = [], [], []

        if hidden_rows:
            suppressed_msg = ", ".join(
                [
                    f"{acc} ({reason or status or 'suppressed'})"
                    for acc, status, reason in hidden_rows
                ]
            )
            self.log(f"Suppressed accessions skipped: {suppressed_msg}", "WARNING")

        if not assembly_dataset:
            if self.accessions:
                hidden_set = {row[0] for row in hidden_rows}
                missing = [acc for acc in self.accessions if acc not in hidden_set]
                if hidden_rows and not missing:
                    self.log(
                        f"All requested accessions are suppressed/hidden; skipping metadata insert: {', '.join(sorted(hidden_set))}",
                        "WARNING",
                    )
                    return True
            return self.handle_exception(
                "No assemblies matched the provided selectors/filters.",
                {
                    "taxid": self.taxid,
                    "accessions": self.accessions,
                    "after": self.released_after,
                    "before": self.released_before,
                    "level": self.level_filter,
                },
            )

        if self.debug_path:
            try:
                os.makedirs(self.debug_path, exist_ok=True)
                timestamp = debug_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
                meta_path = os.path.join(self.debug_path, f"assembly_metadata_{self.task_id}_{timestamp}.jsonl")
                with open(meta_path, "w") as handle:
                    for asm, tax, genome in zip(assembly_dataset, tax_info_dataset, genome_dataset):
                        record = {"assembly": asm, "taxonomy": tax, "genome": genome}
                        handle.write(json.dumps(record, ensure_ascii=True) + "\n")
                tax_update_path = os.path.join(self.debug_path, f"taxonomy_updates_{self.task_id}_{timestamp}.json")
                with open(tax_update_path, "w") as handle:
                    json.dump(taxonomy_update_dict, handle, ensure_ascii=True, indent=2)
                if self.accessions:
                    placeholders = ",".join("?" for _ in self.accessions)
                    self.db_manager.cursor.execute(
                        f"SELECT accession, status, reason FROM Hidden_Genomes WHERE accession IN ({placeholders})",
                        tuple(self.accessions),
                    )
                    hidden_rows = self.db_manager.cursor.fetchall() or []
                    if hidden_rows:
                        hidden_path = os.path.join(self.debug_path, f"hidden_genomes_{self.task_id}_{timestamp}.jsonl")
                        with open(hidden_path, "w") as handle:
                            for acc, status, reason in hidden_rows:
                                handle.write(
                                    json.dumps(
                                        {"accession": acc, "status": status, "reason": reason},
                                        ensure_ascii=True,
                                    )
                                    + "\n"
                                )
                self.log(f"Wrote metadata debug files to {self.debug_path}", "DEBUG")
            except (OSError, TypeError, ValueError) as exc:
                self.log(f"Failed to write metadata debug files: {exc}", "WARNING")

        
    
        # # For testing
        # time.sleep(2)
        # # Random 1/4 chance of failing here for testing
        # test = random.randint(1, 4)
        # if test == 1:
        #     raise Exception("An error occurred while fetching assemblies.")
        # if test == 2:
        #     raise Exception("A different error occurred while fetching assemblies.")
        # else:
        #     pass

        # The code below would run in normal operation (remove the raise above to enable)
        # Check all same number of records
        if not (len(assembly_dataset) == len(tax_info_dataset) == len(genome_dataset)):
            return self.handle_exception("Mismatch in number of records fetched.", {
                "assembly_count": len(assembly_dataset),
                "tax_info_count": len(tax_info_dataset),
                "genome_count": len(genome_dataset)
            })
            

        # Update the taxonomy table so genome foreign keys are satisfied
        if taxonomy_update_dict:
            if not self.db_manager.genomes.insert_taxonomy_information(taxonomy_update_dict):
                return self.handle_exception("Failed to insert taxonomy information.", taxonomy_update_dict)
            

        # Update the assembly and genome tables
        for i in range(len(assembly_dataset)):
            # Update the database with the assembly data
            if self.force_update:
                if not self.db_manager.genomes.upsert_assembly(assembly_dataset[i]):
                    return self.handle_exception("Failed to upsert assembly information.", assembly_dataset[i])
            else:
                if not self.db_manager.genomes.insert_assembly(assembly_dataset[i]):
                    return self.handle_exception("Failed to insert assembly information.", assembly_dataset[i])
            # Update the database with the genome data
            if self.force_update:
                if not self.db_manager.genomes.upsert(genome_dataset[i]):
                    return self.handle_exception("Failed to upsert genome information.", genome_dataset[i])
            else:
                if not self.db_manager.genomes.insert(genome_dataset[i]):
                    return self.handle_exception("Failed to insert genome information.", genome_dataset[i])
            self.log(f"Added information for taxon: {tax_info_dataset[i]['tax_name']}", "DEBUG")

        self.log(f"Registered {len(assembly_dataset)} assemblies in the database.", "INFO")
        return True

class BatchImportLocalAssemblyTask(Task):
    '''A task that takes a path to a folder and imports all locally stored FNA or FAA files into the database with metadata if provided.'''
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.assembly_dir = self.data.get("assembly_dir")
        self.accessions_for_import = list(self.data.get("accessions", []) or [])
        self.accession_to_payload_dict = dict(self.data.get("accession_to_payload_dict", {}) or {})
        self.clean_isoforms = self.payload_bool("clean_isoforms", self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.skip_clean_isoforms = self.payload_bool("skip_clean_isoforms", not self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.clean_skip_gff = self.payload_bool("clean_skip_gff", not self.env_bool("DEFAULT_PROTEOME_USE_GFF", True))
        self.clean_skip_cdhit = self.payload_bool("clean_skip_cdhit", not self.env_bool("DEFAULT_PROTEOME_USE_CDHIT", False))
        self.clean_gff_priority = self.payload_bool("clean_gff_priority", self.env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False))
        self.clean_cdhit_identity = self.data.get("clean_cdhit_identity", self.env_float("DEFAULT_PROTEOME_CDHIT_IDENTITY", None))
        self.clean_max_concurrent = int(self.data.get("clean_max_concurrent", self.env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1)) or 1)
        self.clean_threads_per_job = int(self.data.get("clean_threads_per_job", self.env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1)) or 1)
        self.stage = checkpoint if checkpoint is not None else 0
        if not self.assembly_dir:
            raise ValueError("'assembly_dir' is required.")

    def _infer_file_kind(self, file_name: str):
        """Return 'faa' or 'fna' based on filename, or None if unknown."""
        name = file_name.lower()
        if name.endswith(".faa") or name.endswith(".faa.gz"):
            return "faa"
        if (
            name.endswith(".fna")
            or name.endswith(".fna.gz")
            or name.endswith(".fasta")
            or name.endswith(".fasta.gz")
            or name.endswith(".fa")
            or name.endswith(".fa.gz")
        ):
            return "fna"
        return None

    def _get_missing_accessions(self):
        if not self.accessions_for_import:
            return []

        rows = self.db_manager.genomes.get_many(self.accessions_for_import)
        if not rows:
            present = []
        else:
            present = [r[0] for r in rows]
        missing = [acc for acc in self.accessions_for_import if acc not in present]
        return missing
        
    def _queue_import_tasks(self):
        # Check which accessions are already present in the database.
        # Queue a subtask for missing accession.
        missing = self._get_missing_accessions()
        if not missing:
            return False

        queued = False
        for acc in missing:
            payload = self.accession_to_payload_dict.get(acc)
            if not payload:
                self.log(f"Skipping accession {acc}; no payload data recorded.", "WARNING")
                continue
            self.log(f"Queueing import subtask for accession {acc}", "DEBUG")
            self.queue_subtask(job_type=7, status="P", priority=1, data=payload)
            queued = True
        return queued

    def _import_done(self):
        missing = self._get_missing_accessions()
        for acc in missing:
            self.log(f"Missing accession: {acc}", "DEBUG")
        return len(missing) == 0
    
    def _import_incomplete_message(self):
        missing = self._get_missing_accessions()
        if not missing:
            return None
        summary = f"Unable to import {len(missing)} assemblies. Missing files: {', '.join(missing)}"
        stack = ""
        return (summary, stack)

    def run(self):

        if self.stage == 0:
            # Get the .txt in the assembly_dir
            keys = glob.glob(os.path.join(self.assembly_dir, "*.txt"))
            if not keys or len(keys) > 1:
                return self.handle_exception("There should be exactly one .txt file in the assembly directory.")
            key_file = keys[0]

            # txt is tab delimited with file name and the accession
            with open(key_file, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if not line.strip():
                    continue
                file_name, accession, taxid = line.strip().split()
                if self.accessions_for_import and accession not in self.accessions_for_import:
                    continue
                file_kind = self._infer_file_kind(file_name)
                if not file_kind:
                    return self.handle_exception(
                        f"Unrecognized file extension for {file_name}; expected .faa/.fna/.fasta/.fa",
                        {"line": line.strip()},
                    )
                data = {
                    "taxid": int(taxid),
                    "accession": accession,
                    "copy_to_genome_dir": True,
                }
                payload = self.accession_to_payload_dict.get(accession)
                if not payload:
                    payload = data
                else:
                    if int(taxid) != int(payload.get("taxid")):
                        return self.handle_exception(
                            f"Conflicting taxid for accession {accession}: {payload.get('taxid')} vs {taxid}",
                            {"line": line.strip()},
                        )
                file_path = os.path.join(self.assembly_dir, file_name)
                if payload.get(file_kind):
                    self.log(
                        f"Duplicate {file_kind.upper()} for accession {accession}: {payload.get(file_kind)} and {file_path}",
                        "WARNING",
                    )
                    return self.handle_exception(
                        f"Duplicate {file_kind} entry for accession {accession}",
                        {"line": line.strip()},
                    )
                payload[file_kind] = file_path
                payload["clean_isoforms"] = self.clean_isoforms
                payload["skip_clean_isoforms"] = self.skip_clean_isoforms
                payload["clean_skip_gff"] = self.clean_skip_gff
                payload["clean_skip_cdhit"] = self.clean_skip_cdhit
                payload["clean_gff_priority"] = self.clean_gff_priority
                payload["clean_max_concurrent"] = self.clean_max_concurrent
                payload["clean_threads_per_job"] = self.clean_threads_per_job
                self.accession_to_payload_dict[accession] = payload

            self.accessions_for_import = list(self.accession_to_payload_dict.keys())

            self.checkpoint(0, {
                "accessions": self.accessions_for_import,
                "accession_to_payload_dict": self.accession_to_payload_dict
            })

        # queue sub-tasks for import
        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=self._queue_import_tasks,
            done_fn=self._import_done,
            wait_seconds=int(self.data.get("import_wait_seconds", 5)),
            # Use internal retry counter; treat provided value as the max
            retry_key=None,
            max_retries=int(self.data.get("import_retries", 1)),
            incomplete_message_fn=self._import_incomplete_message,
            retry_incomplete=True,  # do not retry incomplete outcomes; immediately record message
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        if self.clean_isoforms and not self.skip_clean_isoforms and self.accessions_for_import:
            def queue_clean_subtask():
                self.queue_subtask(
                    job_type=31,
                    status="P",
                    priority=1,
                    data={
                        "accessions": self.accessions_for_import,
                        "downloaded_only": True,
                        "skip_gff": self.clean_skip_gff,
                        "skip_cdhit": self.clean_skip_cdhit,
                        "gff_priority": self.clean_gff_priority,
                        "cdhit_identity": self.clean_cdhit_identity,
                        "profile_name": "clean_default",
                        "input_profile": "raw",
                        "set_default": True,
                        "max_concurrent": self.clean_max_concurrent,
                        "threads_per_job": self.clean_threads_per_job,
                        "parent_threads_budget": self.REQUIRED_THREADS,
                    },
                )
                self.log(f"Queued prepare-proteome subtask for {len(self.accessions_for_import)} imported assemblies.", "INFO")
                return True

            outcome = self.manage_subtasks(
                stage=2,
                queue_fn=queue_clean_subtask,
                done_fn=None,
                wait_seconds=0,
                retry_key=None,
                max_retries=0,
                incomplete_message_fn=lambda: ("Import prepare-proteome subtask did not complete.", ""),
                retry_incomplete=False,
            )
            if outcome == "ERROR":
                return "ERROR"
            if outcome is False:
                return False

        self.log(f"Successfully imported {len(self.accessions_for_import)} assemblies.", "INFO")
        return True

class ImportLocalAssemblyTask(Task):
    '''A task that takes a path to a locally stored FNA or FAA and imports it into the database with metadata if provided.'''
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.stage = checkpoint if checkpoint is not None else 0
        self.fna = self.data.get("fna", None)
        self.faa = self.data.get("faa", None)
        self.gff = self.data.get("gff", None)
        self.others = self.data.get("path_to_others", None)
        self.accession = self.data.get("accession", None)
        self.metadata = self.data.get("metadata", {})
        self.taxid = self.data.get("taxid", None)
        self.taxon_name = self.data.get("taxon_name", None)
        self.genus = self.data.get("genus", None)
        self.species = self.data.get("species", None)
        self.location = self.data.get("location", None)
        self.copy_to_genome_dir = self.data.get("copy_to_genome_dir", True)
        self.clean_isoforms = self.payload_bool("clean_isoforms", self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.skip_clean_isoforms = self.payload_bool("skip_clean_isoforms", not self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.clean_skip_gff = self.payload_bool("clean_skip_gff", not self.env_bool("DEFAULT_PROTEOME_USE_GFF", True))
        self.clean_skip_cdhit = self.payload_bool("clean_skip_cdhit", not self.env_bool("DEFAULT_PROTEOME_USE_CDHIT", False))
        self.clean_gff_priority = self.payload_bool("clean_gff_priority", self.env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False))
        self.clean_cdhit_identity = self.data.get("clean_cdhit_identity", self.env_float("DEFAULT_PROTEOME_CDHIT_IDENTITY", None))
        self.clean_max_concurrent = int(self.data.get("clean_max_concurrent", self.env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1)) or 1)
        self.clean_threads_per_job = int(self.data.get("clean_threads_per_job", self.env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1)) or 1)
        if not self.fna and not self.faa:
            raise ValueError("Data must contain either a 'fna' or 'faa' file path.")
        has_taxon_lookup = bool(self.taxon_name or (self.genus and self.species))
        if not self.taxid and not has_taxon_lookup:
            raise ValueError("Data must contain either a 'taxid' or taxon lookup fields.")
        if self.taxid and has_taxon_lookup:
            raise ValueError("Data cannot contain both 'taxid' and taxon lookup fields.")
        if self.location and self.copy_to_genome_dir:
            raise ValueError("Data cannot contain both 'location' and 'copy_to_genome_dir' set to True.")
        
    def import_data(self):
        '''Imports the data into the database'''

        assembly_data = {
            "accession": self.accession,
            "origin": "local",
            "uid": self.metadata.get("uid", None),
            "assembly_method": self.metadata.get("assembly_method", None),
            "assembly_type": self.metadata.get("assembly_type", None),
            "assembly_status": self.metadata.get("assembly_status", None),
            "release_date": self.metadata.get("release_date", None),
            "warnings": self.metadata.get("warnings", None),
            "bioproject_accession": self.metadata.get("bioproject_accession", None),
            "biosample_accession": self.metadata.get("biosample_accession", None),
            "comments": self.metadata.get("comments", None),
            "diploid_role": self.metadata.get("diploid_role", None),
            "refseq_category": self.metadata.get("refseq_category", None),
            "sequencing_tech": self.metadata.get("sequencing_tech", None),
            "submitter": self.metadata.get("submitter", None),
            "contig_l50": self.metadata.get("contig_l50", None),
            "contig_n50": self.metadata.get("contig_n50", None),
            "gc_count": self.metadata.get("gc_count", None),
            "gc_percent": self.metadata.get("gc_percent", None),
            "genome_coverage": self.metadata.get("genome_coverage", None),
            "number_of_component_sequences": self.metadata.get("number_of_component_sequences", None),
            "number_of_contigs": self.metadata.get("number_of_contigs", None),
            "number_of_organelles": self.metadata.get("number_of_organelles", None),
            "number_of_scaffolds": self.metadata.get("number_of_scaffolds", None),
            "scaffold_l50": self.metadata.get("scaffold_l50", None),
            "scaffold_n50": self.metadata.get("scaffold_n50", None),
            "total_number_of_chromosomes": self.metadata.get("total_number_of_chromosomes", None),
            "total_sequence_length": self.metadata.get("total_sequence_length", None),
            "total_ungapped_length": self.metadata.get("total_ungapped_length", None),
        }
        assembly_payload = {"accession": self.accession, "origin": "local"}
        assembly_payload.update({key: value for key, value in assembly_data.items() if key != "accession" and value is not None})

        if not self.db_manager.genomes.upsert_assembly(assembly_payload):
            return self.handle_exception("Failed to insert assembly information.", {assembly_data})

        genome_data = {
            "accession": self.accession,
            "taxid": self.taxid,
            "assembly_level": self.metadata.get("assembly_level", None),
            "assembly_name": self.metadata.get("assembly_name", None),
            "comments": self.metadata.get("comments", None),
        }
        genome_payload = {"accession": self.accession, "taxid": self.taxid}
        genome_payload.update({key: value for key, value in genome_data.items() if key not in {"accession", "taxid"} and value is not None})

        if not self.db_manager.genomes.upsert(genome_payload):
            return self.handle_exception("Failed to insert genome information.", {genome_data})
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if not self.db_manager.genomes.update_status(self.accession, 1, (current_time, self.location, 1 if self.faa else 0)):
            return self.handle_exception(
                "Failed to update genome status for imported assembly.",
                {"accession": self.accession, "location": self.location},
            )
        
        self.log(f"Imported local assembly {self.accession} (taxid {self.taxid}).", "INFO")
        return True

    def _write_taxid_file(self) -> bool:
        if not self.location or not self.taxid:
            return True
        try:
            os.makedirs(self.location, exist_ok=True)
            Path(os.path.join(self.location, "taxid")).write_text(f"{int(self.taxid)}\n", encoding="utf-8")
            return True
        except (OSError, TypeError, ValueError) as exc:
            return self.handle_exception("Failed to write taxid file for discovered custom assembly.", {"location": self.location, "taxid": self.taxid, "error": str(exc)})

    def run(self):
        '''Runs the task'''
        clean_targets: List[str] = []

        if self.stage < 1:
            if self.fna:
                if not os.path.isfile(self.fna):
                    return self.handle_exception(f"FNA file not found: {self.fna}", {"fna": self.fna})
            if self.faa:
                if not os.path.isfile(self.faa):
                    return self.handle_exception(f"FAA file not found: {self.faa}", {"faa": self.faa})
            if not self.accession:
                # Generate a unique accession
                self.accession = f"LOCAL_{int(time.time())}_{random.randint(1000,9999)}"
                self.data["accession"] = self.accession
                self.log(f"No accession provided; generated local accession {self.accession}.", "INFO")

            if not self.taxid:
                # We need to look up the taxid from the taxon name
                # First try the local db
                if self.taxon_name:
                    tax_info = self.db_manager.genomes.get_taxid_for_genus_species(self.taxon_name)
                elif self.genus and self.species:
                    tax_info = self.db_manager.genomes.get_taxid_for_genus_species(f"{self.genus} {self.species}")
                if tax_info:
                    self.taxid = tax_info.get("taxid")
                    self.log(f"Found taxid {self.taxid} for taxon name {self.taxon_name} in local database.", "DEBUG")
                else:
                    #TODO: Implement a remote lookup if not found locally
                    return self.handle_exception("Taxid not found for provided taxon name.", {"taxon_name": self.taxon_name, "genus": self.genus, "species": self.species})

            if self.copy_to_genome_dir:
                # Copy files from their current location to genome_dir/accession_folder/
                genomes_path = self.db_manager.storage.get_root_base("genomes")
                if not genomes_path:
                    return self.handle_exception("Genomes storage root is not configured.", {})
                target_dir = os.path.join(genomes_path, self.accession)
                self.location = target_dir
                os.makedirs(target_dir, exist_ok=True)
                # Copy FNA file
                if self.fna and os.path.isfile(self.fna):
                    target_fna = os.path.join(target_dir, os.path.basename(self.fna))
                    shutil.copy2(self.fna, target_fna)
                    self.fna = target_fna
                    self.log(f"Copied FNA to {target_fna}", "DEBUG")
                # Copy FAA file
                if self.faa and os.path.isfile(self.faa):
                    target_faa = os.path.join(target_dir, os.path.basename(self.faa))
                    shutil.copy2(self.faa, target_faa)
                    self.faa = target_faa
                    self.log(f"Copied FAA to {target_faa}", "DEBUG")
                # Copy GFF file
                if self.gff and os.path.isfile(self.gff):
                    target_gff = os.path.join(target_dir, os.path.basename(self.gff))
                    shutil.copy2(self.gff, target_gff)
                    self.gff = target_gff
                    self.log(f"Copied GFF to {target_gff}", "DEBUG")
                # Copy other files
                if self.others:
                    new_others = []
                    if isinstance(self.others, list):
                        for other in self.others:
                            if os.path.isfile(other):
                                target_other = os.path.join(target_dir, os.path.basename(other))
                                shutil.copy2(other, target_other)
                                new_others.append(target_other)
                                self.log(f"Copied auxiliary file to {target_other}", "DEBUG")
                            else:
                                self.log(f"Other file not found: {other}", "WARNING")
                            self.others = new_others
                    elif os.path.isfile(self.others):
                        target_other = os.path.join(target_dir, os.path.basename(self.others))
                        shutil.copy2(self.others, target_other)
                        self.others = [target_other]
                        self.log(f"Copied auxiliary file to {target_other}", "DEBUG")
                    else:
                        self.log(f"Other files not found or invalid: {self.others}", "WARNING")
            elif not self.location:
                anchor = self.fna or self.faa or self.gff
                if anchor:
                    self.location = os.path.dirname(os.path.abspath(anchor))

            try:
                manifest = _write_phyloodb_md5_manifest(
                    self.location,
                    [self.fna, self.faa, self.gff],
                )
                if manifest:
                    self.log(f"Wrote local checksum manifest to {manifest}", "DEBUG")
            except (OSError, EOFError, gzip.BadGzipFile, ValueError) as exc:
                return self.handle_exception(
                    f"Failed to validate/write local checksum manifest for {self.accession}: {exc}",
                    {"accession": self.accession, "location": self.location},
                )

            result = self.import_data()
            if result is not True:
                return result
            wrote_taxid = self._write_taxid_file()
            if wrote_taxid is not True:
                return wrote_taxid
            try:
                _register_genome_file_artifacts(self, self.accession, self.location)
            except Exception as exc:  # boundary: convert artifact registration failure into this task error
                return self.handle_exception(
                    f"Failed to register imported genome artifacts for {self.accession}: {exc}",
                    {"accession": self.accession, "location": self.location},
                )
            if self.clean_isoforms and not self.skip_clean_isoforms and self.faa:
                clean_targets = [self.accession]
            self.data["_clean_targets"] = clean_targets
        else:
            clean_targets = list(self.data.get("_clean_targets", []) or [])

        if not (self.clean_isoforms and not self.skip_clean_isoforms and clean_targets):
            return True

        def queue_clean_subtask():
            self.queue_subtask(
                job_type=31,
                status="P",
                priority=1,
                data={
                    "accessions": clean_targets,
                    "downloaded_only": True,
                    "skip_gff": self.clean_skip_gff,
                    "skip_cdhit": self.clean_skip_cdhit,
                    "gff_priority": self.clean_gff_priority,
                    "cdhit_identity": self.clean_cdhit_identity,
                    "profile_name": "clean_default",
                    "input_profile": "raw",
                    "set_default": True,
                    "max_concurrent": 1,
                    "threads_per_job": self.clean_threads_per_job,
                    "parent_threads_budget": self.REQUIRED_THREADS,
                },
            )
            self.log(f"Queued prepare-proteome subtask for imported accession {self.accession}.", "INFO")
            return True

        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=queue_clean_subtask,
            done_fn=None,
            wait_seconds=0,
            retry_key=None,
            max_retries=0,
            incomplete_message_fn=lambda: ("Isoform cleaning did not complete.", ""),
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False
        return True

class ExampleTask(Task):
    """An example task that does nothing"""
    def __init__(self, db_path, task_id, required_threads=1):
        super().__init__(db_path, task_id, data=None, required_threads=required_threads)

    def run(self):
        self.db_manager.connect()
        self.update_on_start()
        time.sleep(random.uniform(1, 4))
        self.update_on_complete()
        self.db_manager.close()

class DownloadAssembliesTask(Task):
    """A task that downloads assemblies from NCBI"""
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.accessions: List[str] = []
        self.stage = checkpoint if checkpoint is not None else 0
        self.taxid = self.data.get("taxid", None)
        if self.taxid is not None:
            try:
                self.taxid = int(self.taxid)
            except (TypeError, ValueError) as exc:
                raise ValueError("'taxid' must be an integer if provided.") from exc
        self.protein = self.data.get("protein", False)
        self.max_concurrent = self.data.get("max_concurrent", 1)
        self.force_redownload = self.data.get("force_redownload", False)
        self.download_retries = int(self.data.get("download_retries", 0) or 0)
        self.clean_isoforms = self.payload_bool("clean_isoforms", self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.skip_clean_isoforms = self.payload_bool("skip_clean_isoforms", not self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.clean_skip_gff = self.payload_bool("clean_skip_gff", not self.env_bool("DEFAULT_PROTEOME_USE_GFF", True))
        self.clean_skip_cdhit = self.payload_bool("clean_skip_cdhit", not self.env_bool("DEFAULT_PROTEOME_USE_CDHIT", False))
        self.clean_gff_priority = self.payload_bool("clean_gff_priority", self.env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False))
        self.clean_cdhit_identity = self.data.get("clean_cdhit_identity", self.env_float("DEFAULT_PROTEOME_CDHIT_IDENTITY", None))
        self.clean_max_concurrent = int(self.data.get("clean_max_concurrent", self.env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1)) or 1)
        self.clean_threads_per_job = int(self.data.get("clean_threads_per_job", self.env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1)) or 1)
        self.rule_quantity = self.data.get("quantity")
        provided_rank = self.data.get("rank")
        if self.rule_quantity is not None and provided_rank is None:
            self.rule_rank = "species"
        else:
            self.rule_rank = provided_rank

    def run(self):
        '''Runs the task'''
        clean_targets: List[str] = []

        if self.stage < 1:
            # Initialize the NCBI helper
            vars = self.db_manager.env.get_many(['EMAIL'])
            email = vars.get('EMAIL')
            genomes_path = self.db_manager.storage.get_root_base("genomes")
            if not email:
                return self.handle_exception("Email environment variable not set in the database.", vars)
            if not genomes_path:
                return self.handle_exception("Genomes storage root is not configured.", vars)
            downloaded_only_flag = self.data.get("downloaded_only")
            released_after = self.data.get("after")
            released_before = self.data.get("before")
            level_filter = self.data.get("level")
            primary_only = self.data.get("primary_only")
            use_busco = self.data.get("use_busco")
            min_complete = self.data.get("min_completeness")
            min_sc = self.data.get("min_single_copy_complete")

            try:
                self.accessions = self.prepare_selectors(
                    taxid=self.taxid,
                    rule_quantity=self.rule_quantity,
                    rule_rank=self.rule_rank,
                    downloaded_only=downloaded_only_flag,
                    released_after=released_after,
                    released_before=released_before,
                    level=level_filter,
                    primary_only=primary_only,
                    use_busco=use_busco,
                    min_completeness=min_complete,
                    min_single_copy_complete=min_sc,
                )
            except ValueError as exc:
                return self.handle_exception(
                    str(exc),
                    {
                        "taxid": self.taxid,
                        "quantity": self.rule_quantity,
                        "rank": self.rule_rank,
                        "after": released_after,
                        "before": released_before,
                        "level": level_filter,
                        "primary_only": primary_only,
                        "use_busco": use_busco,
                        "min_completeness": min_complete,
                        "min_single_copy_complete": min_sc,
                    },
                )

            self.log(
                f"Downloading {len(self.accessions)} assemblies to {genomes_path} "
                f"(protein={'yes' if self.protein else 'no'}, force_redownload={'yes' if self.force_redownload else 'no'}).",
                "INFO",
            )
            ncbi_helper = NCBIHelper(email, db_manager=self.db_manager)
            status_by_accession: Dict[str, int] = {}
            protein_by_accession: Dict[str, bool] = {}

            import concurrent.futures

            existing_paths = {
                accession: self.db_manager.genomes.resolve_path(accession)
                for accession in self.accessions
            }

            def download_accession(accession, target_path, ncbi_helper, protein_required):
                def _validate_gzip(path: str) -> bool:
                    if not path or not os.path.isfile(path):
                        return False
                    try:
                        with gzip.open(path, "rb") as handle:
                            for _chunk in iter(lambda: handle.read(1024 * 1024), b""):
                                pass
                        return True
                    except (OSError, EOFError, gzip.BadGzipFile) as exc:
                        self.log(f"Gzip validation failed for {path}: {exc}", "ERROR")
                        return False

                # Rename the current thread
                threading.current_thread().name = f"DL_{accession}"
                path = os.path.abspath(str(target_path))
                # Check for existing files in target directory
                self.log(f"Initiating download for {accession}", "DEBUG")
                dir_exists = os.path.isdir(path)
                fna_exists = any(fname.endswith(".fna.gz") for fname in os.listdir(path)) if dir_exists else False
                faa_exists = any(fname.endswith(".faa.gz") for fname in os.listdir(path)) if dir_exists else False

                if dir_exists and (fna_exists or faa_exists) and not self.force_redownload:
                    invalid_existing = []
                    for fname in os.listdir(path):
                        if fname.endswith(".gz"):
                            fpath = os.path.join(path, fname)
                            if not _validate_gzip(fpath):
                                invalid_existing.append(fpath)
                    if invalid_existing:
                        self.log(
                            f"Existing download for {accession} has invalid gzip files; redownloading: "
                            f"{', '.join(invalid_existing)}",
                            "WARNING",
                        )
                        for bad in invalid_existing:
                            try:
                                os.remove(bad)
                            except OSError as exc:
                                self.error(f"Failed to remove invalid gzip {bad}: {exc}")
                                return (accession, 4, False, False, path)
                        fna_exists = any(fname.endswith(".fna.gz") for fname in os.listdir(path)) if os.path.isdir(path) else False
                        faa_exists = any(fname.endswith(".faa.gz") for fname in os.listdir(path)) if os.path.isdir(path) else False
                        if fna_exists or faa_exists:
                            # Some files were still valid, but redownload the accession as a coherent set.
                            for fname in os.listdir(path):
                                try:
                                    file_path = os.path.join(path, fname)
                                    if os.path.isfile(file_path) or os.path.islink(file_path):
                                        os.remove(file_path)
                                    elif os.path.isdir(file_path):
                                        shutil.rmtree(file_path)
                                except OSError as exc:
                                    self.error(f"Failed to clean existing download {file_path}: {exc}")
                                    return (accession, 4, False, False, path)
                    else:
                        # Do not re-download; validate success requirements
                        if protein_required and not faa_exists:
                            self.log(f"Protein required but FAA missing for {accession}; not downloading due to existing folder.", "ERROR")
                            return (accession, 2, False, True, path)
                        # Consider success based on what is present
                        self.log(f"Download already exists for {accession} at {path}; skipping re-download.", "DEBUG")
                        return (accession, 0, faa_exists, True, path)

                if dir_exists and self.force_redownload:
                    self.log(f"Force redownload enabled; deleting existing files for {accession}.", "INFO")
                    # Delete everything in the genome folder
                    for fname in os.listdir(path):
                        try:
                            file_path = os.path.join(path, fname)
                            if os.path.isfile(file_path) or os.path.islink(file_path):
                                os.remove(file_path)
                            elif os.path.isdir(file_path):
                                shutil.rmtree(file_path)
                        except OSError as e:
                            self.error(f"Failed to delete {file_path}: {e}")
                            return (accession, 4, False, False, path)  # General error, no protein

                attempts = self.download_retries + 1
                for attempt in range(1, attempts + 1):
                    try:
                        ncbi_helper.download_assembly(accession=accession, location=path, uid=None, protein=True)
                        self.log(f"Successfully downloaded assembly for {accession} to {path}.", "DEBUG")
                        # Validate gzip files
                        invalid = []
                        for fname in os.listdir(path):
                            if fname.endswith(".gz"):
                                fpath = os.path.join(path, fname)
                                if not _validate_gzip(fpath):
                                    invalid.append(fpath)
                        if invalid:
                            for bad in invalid:
                                try:
                                    os.remove(bad)
                                except OSError:
                                    pass
                            raise FnaDownloadError(f"Validation failed for files: {', '.join(invalid)}")

                        protein_present = any(fname.endswith(".faa.gz") for fname in os.listdir(path)) if os.path.isdir(path) else False
                        return (accession, 0, protein_present, False, path)  # Success
                    except FnaDownloadError as e:
                        self.error(f"Failed to download FNA for {accession} (attempt {attempt}/{attempts}): {e}")
                        if attempt == attempts:
                            return (accession, 1, False, False, path)
                    except FaaDownloadError as e:
                        if protein_required:
                            self.error(
                                f"Failed to download FAA for {accession} (attempt {attempt}/{attempts}). "
                                f"Protein output is required, so this accession will fail. {e}"
                            )
                        else:
                            self.log(
                                f"No FAA file found for {accession}; continuing because protein output was not requested.",
                                "WARNING",
                            )
                        if attempt == attempts:
                            if protein_required:
                                return (accession, 2, False, False, path)
                            os.makedirs(path, exist_ok=True)
                            return (accession, 0, False, False, path)
                    except GffDownloadError as e:
                        self.error(f"Failed to download GFF for {accession} (attempt {attempt}/{attempts}): {e}")
                        if attempt == attempts:
                            protein_present = any(fname.endswith(".faa.gz") for fname in os.listdir(path)) if os.path.isdir(path) else False
                            return (accession, 0, protein_present, False, path)
                    except Exception as e:  # boundary: isolate one download attempt in a per-accession batch
                        self.error(f"Failed to download assembly for {accession} (attempt {attempt}/{attempts}): {e}")
                        if attempt == attempts:
                            return (accession, 4, False, False, path)
                    # If here, retry loop continues

            with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_concurrent, thread_name_prefix="DL_") as executor:
                futures = {
                    executor.submit(
                        download_accession,
                        accession,
                        existing_paths.get(accession) or os.path.join(genomes_path, accession),
                        ncbi_helper,
                        self.protein,
                    ): accession
                    for accession in self.accessions
                }
                for future in concurrent.futures.as_completed(futures):
                    accession, result, protein_present, reused_existing, target_path = future.result()
                    status_by_accession[accession] = result
                    protein_by_accession[accession] = bool(protein_present)
                    if result == 0:
                        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        if reused_existing:
                            self.log(f"Using existing download for {accession} at {current_time}", "DEBUG")
                        else:
                            self.log(f"Download completed for {accession} at {current_time}", "DEBUG")
                        if not self.db_manager.genomes.update_status(accession, 1, (current_time, target_path, protein_present)):
                            self.error(f"Failed to update genome status for {accession}")
                            self.log(f"Failed to update genome status for {accession}", "ERROR")
                            status_by_accession[accession] = 5
                        else:
                            try:
                                _register_genome_file_artifacts(self, accession, target_path)
                            except Exception as exc:  # boundary: collect one failed artifact registration without hiding batch status
                                self.collect_batch_failure(accession, "register downloaded genome artifacts", exc)
                                status_by_accession[accession] = 5

            if not (status_by_accession and all(s == 0 for s in status_by_accession.values())):
                error_types = {1: "FNA", 2: "FAA", 3: "GFF", 4: "General", 5: "DB update"}
                failures = [
                    (acc, error_types.get(status_by_accession.get(acc), "Unknown"))
                    for acc in self.accessions
                    if status_by_accession.get(acc) not in (None, 0)
                ]
                recorded = {failure.item for failure in self._batch_failures}
                for accession, error_type in failures:
                    if accession not in recorded:
                        self.collect_batch_failure(
                            accession,
                            "download assembly",
                            RuntimeError(error_type),
                        )
                return self.fail_if_batch_failures("Assembly download batch failed")

            self.log(f"Downloaded {len(self.accessions)}/{len(self.accessions)} assemblies successfully.", "INFO")

            if self.clean_isoforms and not self.skip_clean_isoforms:
                clean_targets = [acc for acc in self.accessions if protein_by_accession.get(acc)]
            self.data["_clean_targets"] = clean_targets
        else:
            clean_targets = list(self.data.get("_clean_targets", []) or [])

        if not (self.clean_isoforms and not self.skip_clean_isoforms and clean_targets):
            return True

        def queue_clean_subtask():
            self.queue_subtask(
                job_type=31,
                status="P",
                priority=1,
                data={
                    "accessions": clean_targets,
                    "downloaded_only": True,
                    "skip_gff": self.clean_skip_gff,
                    "skip_cdhit": self.clean_skip_cdhit,
                    "gff_priority": self.clean_gff_priority,
                    "cdhit_identity": self.clean_cdhit_identity,
                    "profile_name": "clean_default",
                    "input_profile": "raw",
                    "set_default": True,
                    "max_concurrent": self.clean_max_concurrent,
                    "threads_per_job": self.clean_threads_per_job,
                    "parent_threads_budget": self.REQUIRED_THREADS,
                },
            )
            self.log(f"Queued prepare-proteome subtask for {len(clean_targets)} downloaded assemblies.", "INFO")
            return True

        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=queue_clean_subtask,
            done_fn=None,
            wait_seconds=0,
            retry_key=None,
            max_retries=0,
            incomplete_message_fn=lambda: ("Isoform cleaning did not complete.", ""),
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False
        return True
