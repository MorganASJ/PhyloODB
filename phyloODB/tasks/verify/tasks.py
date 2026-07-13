import os
import glob
import gzip
import shutil
import csv
import json
import re
import hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE, is_staged_busco_input_path, staged_busco_input_profile_name
from ..task import Task
from ..reporting import resolve_report_run_dir
from ...selector_utils import normalize_accessions


def _md5_file(path: str) -> Optional[str]:
    if not path or not os.path.isfile(path):
        return None
    hasher = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _parse_md5_manifest(path: str) -> Dict[str, str]:
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


def _default_verify_report_dir(task: Task) -> Path:
    cached = getattr(task, "_verify_report_dir", None)
    if cached is not None:
        return cached
    out = resolve_report_run_dir(
        task,
        namespace="verify-reports",
        explicit_root=task.data.get("report_root"),
        run_label=task.data.get("library_name") or task.data.get("library_id") or task.data.get("taxid"),
        cache_attr="_verify_report_dir",
    )
    setattr(task, "_verify_report_dir", out)
    return out


def _resolve_verify_report_path(task: Task, default_name: str) -> Path:
    report = task.data.get("report")
    if report:
        return Path(report)
    return _default_verify_report_dir(task) / default_name


def _write_tsv_report(task: Task, default_name: str, headers: List[str], rows: List[Dict[str, Any]]) -> None:
    out_path = _resolve_verify_report_path(task, default_name)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(headers)]
    for row in rows:
        lines.append("\t".join(str(row.get(col, "")) for col in headers))
    out_path.write_text("\n".join(lines), encoding="utf-8")
    task.log(f"Wrote verify report to {out_path}", "INFO")


def _action_tokens(actions: List[str] | str | None) -> List[str]:
    if not actions:
        return []
    if isinstance(actions, str):
        return [token for token in actions.split(";") if token]
    return [str(token) for token in actions if token]


def _run_verify_workers(task: Task, items: List[Any], worker_fn, item_label: str) -> Tuple[List[Any], List[Tuple[Any, Exception]]]:
    if not items:
        return [], []
    worker_count = min(max(1, int(task.REQUIRED_THREADS or 1)), len(items))
    results: List[Any] = [None] * len(items)
    errors: List[Tuple[Any, Exception]] = []
    if worker_count <= 1:
        for idx, item in enumerate(items):
            try:
                results[idx] = worker_fn(item)
            except Exception as exc:  # boundary: one sequential verify worker failure is aggregated.
                errors.append((item, exc))
        return results, errors
    task.log(f"Running {len(items)} {item_label} verification item(s) with {worker_count} worker(s).", "INFO")
    with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix=f"verify-{task.task_id}") as executor:
        future_map = {
            executor.submit(worker_fn, item): (idx, item)
            for idx, item in enumerate(items)
        }
        for future in as_completed(future_map):
            idx, item = future_map[future]
            try:
                results[idx] = future.result()
            except Exception as exc:  # boundary: one parallel verify worker failure is aggregated.
                errors.append((item, exc))
    return results, errors


def _parallel_error_message(prefix: str, errors: List[Tuple[Any, Exception]]) -> str:
    preview = "; ".join(f"{item}: {exc}" for item, exc in errors[:5])
    if len(errors) > 5:
        preview += f"; ... and {len(errors) - 5} more"
    return f"{prefix} for {len(errors)} item(s): {preview}"


def _sync_profile_state(task: Task, accession: str, folder: str) -> dict[str, Any]:
    if not folder or not os.path.isdir(folder):
        return {
            "profiles": {},
            "ready_profiles": [],
            "stale_profiles": [],
            "default_profile": None,
            "has_protein": False,
            "removed_staged_inputs": [],
        }
    sync = task.db_manager.proteomes.sync_profiles_from_filesystem(str(accession), str(folder), set_default=True)
    removed_staged_inputs = _cleanup_staged_busco_inputs(task, accession, folder)
    ready_profiles = set(sync.get("ready_profiles") or [])
    task.db_manager.genomes.set_protein(str(accession), bool(sync.get("has_protein")))
    task.db_manager.genomes.set_isoforms_cleaned(
        str(accession),
        bool(task.db_manager.proteomes.get_default_cleaned_profile_name(str(accession)) or any(name != "raw" for name in ready_profiles)),
    )
    sync["removed_staged_inputs"] = removed_staged_inputs
    return sync


def _cleanup_staged_busco_inputs(task: Task, accession: str, folder: str) -> List[str]:
    if not folder or not os.path.isdir(folder):
        return []
    running = task.db_manager.cursor.execute(
        """
        SELECT COUNT(*)
        FROM BUSCO_Runs
        WHERE accession = ?
          AND COALESCE(status, '') = 'running'
        """,
        (str(accession),),
    ).fetchone()
    if running and int(running[0] or 0) > 0:
        return []
    removed: List[str] = []
    for fname in sorted(os.listdir(folder)):
        path = os.path.join(folder, fname)
        if not os.path.isfile(path) or not is_staged_busco_input_path(fname):
            continue
        try:
            os.remove(path)
            removed.append(path)
        except FileNotFoundError:
            continue
        except OSError as exc:
            task.log(f"{accession}: failed to remove staged BUSCO input {path}: {exc}", "WARNING")
    if removed:
        task.log(f"{accession}: removed {len(removed)} stale staged BUSCO input file(s).", "INFO")
    return removed


def _inspect_profile_state(task: Task, accession: str, folder: str) -> dict[str, Any]:
    ready_names: set[str] = set()
    stale_names: set[str] = set()
    default_profile: Optional[str] = None
    has_protein = False

    if folder and os.path.isdir(folder):
        for fname in sorted(os.listdir(folder)):
            low = fname.lower()
            if (
                low.endswith((".faa", ".faa.gz"))
                and ".archive" not in low
                and not is_staged_busco_input_path(fname)
            ):
                ready_names.add("raw")
                has_protein = True
                break

        profiles_dir = os.path.join(folder, "proteome_profiles")
        if os.path.isdir(profiles_dir):
            for fname in sorted(os.listdir(profiles_dir)):
                path = os.path.join(profiles_dir, fname)
                if not os.path.isfile(path):
                    continue
                low = fname.lower()
                for suffix in (".faa.gz", ".faa"):
                    if low.endswith(suffix):
                        profile_name = fname[: -len(suffix)].strip()
                        if profile_name:
                            ready_names.add(profile_name)
                        break

    for row in task.db_manager.proteomes.list_profiles(accessions=[str(accession)]):
        if not row or row[2] is None:
            continue
        profile_name = str(row[2])
        status = str(row[6] or "ready").strip().lower()
        resolved = task.db_manager.proteomes.resolve_path(row)
        exists = bool(resolved and os.path.exists(resolved))
        if status == "ready" and exists:
            ready_names.add(profile_name)
            if row[9]:
                default_profile = profile_name
        else:
            stale_names.add(profile_name)

    stale_names.difference_update(ready_names)

    if default_profile is None:
        default_from_db = task.db_manager.proteomes.get_default_profile_name(str(accession))
        if default_from_db and default_from_db in ready_names:
            default_profile = str(default_from_db)
        elif DEFAULT_CLEAN_PROFILE in ready_names:
            default_profile = DEFAULT_CLEAN_PROFILE
        else:
            nonraw_ready = sorted(name for name in ready_names if name != "raw")
            if nonraw_ready:
                default_profile = nonraw_ready[0]
            elif "raw" in ready_names:
                default_profile = "raw"

    clean_profile_ready = any(name != "raw" for name in ready_names)
    if clean_profile_ready:
        profile_status = "ready"
    elif "raw" in ready_names:
        profile_status = "raw_only"
    elif stale_names:
        profile_status = "stale"
    else:
        profile_status = "missing"

    return {
        "profiles": {},
        "ready_profiles": sorted(ready_names),
        "stale_profiles": sorted(stale_names),
        "default_profile": default_profile,
        "has_protein": has_protein,
        "clean_profile_ready": clean_profile_ready,
        "status": profile_status,
    }


class VerifyAssemblyTask(Task):
    """Verify assembly folders: existence, gzip integrity, optional tidy/discover/reacquire."""

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads=required_threads)
        self.stage = checkpoint if checkpoint is not None else 0
        self.repair = bool(self.data.get("repair", False))
        self.reacquire = bool(self.data.get("reacquire", False))
        self.discover = bool(self.data.get("discover", False))
        self.discover_protein = bool(self.data.get("discover_protein", False))
        self.tidy = bool(self.data.get("tidy", False))
        self.organise = bool(self.data.get("organise", False))
        self.organise_check_only = bool(self.data.get("organise_check_only", False))
        self.report = self.data.get("report")
        self.report_root = self.data.get("report_root")
        self.all = bool(self.data.get("all", False))
        self.downloaded_only = bool(self.data.get("downloaded_only", not bool(self.discover)))
        self.clean_isoforms = self.payload_bool("clean_isoforms", False)
        self.skip_clean_isoforms = self.payload_bool("skip_clean_isoforms", not self.env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True))
        self.clean_skip_gff = self.payload_bool("clean_skip_gff", not self.env_bool("DEFAULT_PROTEOME_USE_GFF", True))
        self.clean_skip_cdhit = self.payload_bool("clean_skip_cdhit", not self.env_bool("DEFAULT_PROTEOME_USE_CDHIT", False))
        self.clean_gff_priority = self.payload_bool("clean_gff_priority", self.env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False))
        self.clean_revert_from_archive = bool(self.data.get("clean_revert_from_archive", False))
        self.clean_cdhit_identity = self.data.get("clean_cdhit_identity", self.env_float("DEFAULT_PROTEOME_CDHIT_IDENTITY", None))
        self.clean_max_concurrent = int(self.data.get("clean_max_concurrent", self.env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1)) or 1)
        self.clean_threads_per_job = int(self.data.get("clean_threads_per_job", self.env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1)) or 1)

    def _assembly_join(self, items: List[str]) -> str:
        return ";".join(sorted(items)) if items else ""

    def _mark_assembly_checks_skipped(self, row: Dict[str, str]) -> None:
        for key in (
            "nuc_gz_valid",
            "nuc_gz_invalid",
            "nuc_ok",
            "nuc_missing_or_invalid",
            "prot_gz_valid",
            "prot_gz_invalid",
            "prot_ok",
            "prot_missing_or_invalid",
        ):
            row[key] = "skipped"

    def _finalize_assembly_row(
        self,
        row: Dict[str, str],
        actions: List[str],
        reacquire_reasons: List[str],
        *,
        status_before: Optional[int],
        protein_before: Optional[int],
        isoforms_cleaned_before: bool,
        artifact_statuses: Optional[Dict[str, str]] = None,
        profile_state: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        row["actions"] = ";".join(actions) if actions else ""
        row["reacquire_reason"] = self._assembly_join(reacquire_reasons)
        accession = row["accession"]
        status_after = self.db_manager.genomes.get_status(accession)
        genome_after = self.db_manager.genomes.get(accession)
        col_map = self.db_manager.genomes.get_column_map()
        protein_after = None
        isoforms_cleaned_after = False
        if genome_after:
            if "protein" in col_map:
                try:
                    protein_after = genome_after[col_map["protein"]]
                except (TypeError, IndexError):
                    protein_after = None
            if "isoforms_cleaned" in col_map:
                try:
                    isoforms_cleaned_after = bool(genome_after[col_map["isoforms_cleaned"]])
                except (TypeError, IndexError, ValueError):
                    isoforms_cleaned_after = False
        assembly_usable = row.get("nuc_ok") == "yes"
        proteome_usable = row.get("prot_ok") == "yes"
        artifact_statuses = artifact_statuses or {}
        profile_state = profile_state or {}
        ready_profiles = [str(name) for name in (profile_state.get("ready_profiles") or [])]
        stale_profiles = [str(name) for name in (profile_state.get("stale_profiles") or [])]
        default_profile = str(profile_state.get("default_profile") or "")
        profile_status = str(profile_state.get("status") or "missing")
        clean_profile_ready = bool(profile_state.get("clean_profile_ready"))
        critical_reasons = []
        if row.get("folder_exists") != "yes":
            critical_reasons.append("folder_missing")
        if not assembly_usable:
            critical_reasons.append("genome_fna_missing_or_invalid")
        if not proteome_usable:
            critical_reasons.append("proteome_missing_or_invalid")
        if isoforms_cleaned_after and not clean_profile_ready:
            critical_reasons.append("clean_proteome_profile_missing_or_invalid")
        if reacquire_reasons:
            for reason in reacquire_reasons:
                if reason not in critical_reasons:
                    critical_reasons.append(reason)
        artifact_warnings = [
            f"{artifact_type}:{artifact_status}"
            for artifact_type, artifact_status in artifact_statuses.items()
            if artifact_status != "ready"
        ]
        if profile_status == "stale":
            artifact_warnings.append("proteome_profiles:stale")
        summary_row = {
            "accession": accession,
            "folder": row.get("folder", ""),
            "folder_exists": row.get("folder_exists", ""),
            "status_before": "" if status_before is None else str(status_before),
            "status_after": "" if status_after is None else str(status_after),
            "status_changed": "yes" if status_before != status_after else "no",
            "assembly_usable": "yes" if assembly_usable else "no",
            "proteome_usable": "yes" if proteome_usable else "no",
            "protein_before": "" if protein_before is None else str(protein_before),
            "protein_after": "" if protein_after is None else str(protein_after),
            "isoforms_cleaned_before": "yes" if isoforms_cleaned_before else "no",
            "isoforms_cleaned_after": "yes" if isoforms_cleaned_after else "no",
            "assembly_stale": "yes" if not assembly_usable else "no",
            "proteome_stale": "yes" if not proteome_usable else "no",
            "fna_artifact_status": artifact_statuses.get("genome_fna", "n/a"),
            "faa_artifact_status": artifact_statuses.get("genome_faa", "n/a"),
            "gff_artifact_status": artifact_statuses.get("genome_gff", "n/a"),
            "default_proteome_profile": default_profile,
            "ready_proteome_profiles": self._assembly_join(ready_profiles),
            "stale_proteome_profiles": self._assembly_join(stale_profiles),
            "proteome_profile_status": profile_status,
            "deactivated": "yes" if status_after is not None and int(status_after) <= 0 else "no",
            "critical_reasons": self._assembly_join(critical_reasons),
            "artifact_warnings": self._assembly_join(artifact_warnings),
            "reacquire_reason": row.get("reacquire_reason", ""),
            "actions": row.get("actions", ""),
        }
        action_tokens = _action_tokens(actions)
        changed = any(
            (
                summary_row["status_changed"] == "yes",
                summary_row["protein_before"] != summary_row["protein_after"],
                summary_row["isoforms_cleaned_before"] != summary_row["isoforms_cleaned_after"],
                row.get("discover_marked") == "yes",
                row.get("reacquire_queued") == "yes",
                row.get("organise_action") == "queued",
                any(
                    token.startswith(
                        (
                            "removed_",
                            "compressed_",
                            "discover_marked",
                            "reacquire_queued",
                            "organise_queued",
                        )
                    )
                    for token in action_tokens
                ),
            )
        )
        return {"report_row": row, "summary_row": summary_row if changed else None}

    def _run_assembly_worker(self, accession: str, genome_dir: str, split_isolated: bool) -> Dict[str, Any]:
        worker = VerifyAssemblyTask(
            self.db_manager.get_path(),
            self.task_id,
            0,
            data=json.dumps(self.data),
            required_threads=1,
        )
        try:
            return worker._process_assembly_accession(accession, genome_dir, split_isolated)
        finally:
            worker.db_manager.close()

    def _find_gff_files(self, folder: str) -> List[str]:
        matches = []
        for fname in sorted(os.listdir(folder)):
            low = fname.lower()
            if low.endswith((".gff", ".gff.gz", ".gff3", ".gff3.gz")):
                matches.append(os.path.join(folder, fname))
        return matches

    def _strip_gff_ext(self, path: str) -> str:
        base = os.path.basename(path)
        lower = base.lower()
        for suffix in (".gff.gz", ".gff3.gz", ".gff", ".gff3"):
            if lower.endswith(suffix):
                return base[: -len(suffix)]
        return os.path.splitext(base)[0]

    def _reconcile_owner_artifacts(self, owner_type: str, owner_id: str, expected: Dict[str, Optional[str]], *, metadata: Optional[dict] = None, repair: bool = False) -> Dict[str, str]:
        current = self.db_manager.artifacts.find(owner_type=owner_type, owner_id=owner_id)
        by_type = {str(row[3]): row for row in current}
        statuses: Dict[str, str] = {}
        ncbi_manifest = expected.get("genome_ncbi_md5checksums")
        phyloodb_manifest = expected.get("genome_phyloodb_md5checksums")
        manifest_path = ncbi_manifest or phyloodb_manifest
        manifest_checksums = _parse_md5_manifest(manifest_path) if manifest_path else {}
        checksum_source = "NCBI md5checksums.txt" if ncbi_manifest else "PhyloODB local import checksum manifest"
        for artifact_type, path in expected.items():
            row = by_type.get(artifact_type)
            if path and os.path.exists(path):
                artifact_metadata = dict(metadata or {})
                if artifact_type == "genome_ncbi_md5checksums":
                    artifact_metadata["source"] = "NCBI assembly FTP md5checksums.txt"
                elif artifact_type == "genome_phyloodb_md5checksums":
                    artifact_metadata["source"] = "PhyloODB local import checksum manifest"
                elif manifest_checksums.get(os.path.basename(path)):
                    artifact_metadata["checksum_source"] = checksum_source
                expected_checksum = (
                    _md5_file(path)
                    if artifact_type in {"genome_ncbi_md5checksums", "genome_phyloodb_md5checksums"}
                    else manifest_checksums.get(os.path.basename(path))
                )
                if repair:
                    self.db_manager.artifacts.register(
                        owner_type=owner_type,
                        owner_id=owner_id,
                        artifact_type=artifact_type,
                        path=path,
                        is_dir=os.path.isdir(path),
                        format="directory" if os.path.isdir(path) else os.path.splitext(path)[1].lstrip(".") or None,
                        checksum=expected_checksum,
                        size_bytes=os.path.getsize(path) if os.path.isfile(path) else None,
                        metadata=artifact_metadata,
                    )
                    current = self.db_manager.artifacts.find(owner_type=owner_type, owner_id=owner_id)
                    by_type = {str(updated[3]): updated for updated in current}
                    row = by_type.get(artifact_type)
                checksum_status = "ready"
                if row and len(row) > 12 and row[12] and os.path.isfile(path):
                    try:
                        actual = _md5_file(path)
                    except OSError as exc:
                        self.log(f"{owner_id}: failed to checksum {artifact_type} at {path}: {exc}", "WARNING")
                        actual = None
                    if actual and str(actual).lower() != str(row[12]).lower():
                        checksum_status = "checksum_mismatch"
                        if repair:
                            self.db_manager.artifacts.set_status(int(row[0]), "stale", checksum=str(row[12]))
                statuses[artifact_type] = checksum_status
            elif row:
                if repair:
                    self.db_manager.artifacts.set_status(int(row[0]), "stale")
                statuses[artifact_type] = "stale"
            else:
                statuses[artifact_type] = "missing"
        return statuses

    def _append_comment(self, accession: str, note: str) -> None:
        """Append a timestamped note to Genome.comments."""
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        try:
            self.db_manager.cursor.execute("SELECT comments FROM Genome WHERE accession = ?", (accession,))
            row = self.db_manager.cursor.fetchone()
            existing = row[0] if row else ""
            new_comment = f"[{ts}] {note}"
            combined = f"{existing}\n{new_comment}" if existing else new_comment
            self.db_manager.cursor.execute("UPDATE Genome SET comments = ? WHERE accession = ?", (combined, accession))
            self.db_manager.commit()
        except Exception as exc:  # boundary: comment enrichment failure should not block verification.
            self.log(f"Failed to append verify comment for {accession}: {exc}", "WARNING")

    def _validate_gzip(self, path: str) -> bool:
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

    def _compress_if_needed(self, src: str) -> Optional[str]:
        """Compress an uncompressed file if no gz exists; return gz path or None."""
        gz_path = src + ".gz"
        if os.path.exists(gz_path):
            return None
        try:
            with open(src, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            return gz_path
        except (OSError, EOFError, gzip.BadGzipFile) as exc:
            self.log(f"Failed to compress {src}: {exc}", "ERROR")
            return None

    def _find_files(self, folder: str, base_ext: str) -> Tuple[List[str], List[str]]:
        compressed = []
        uncompressed = []
        for fname in os.listdir(folder):
            if fname.endswith(base_ext + ".gz"):
                compressed.append(os.path.join(folder, fname))
            elif fname.endswith(base_ext):
                uncompressed.append(os.path.join(folder, fname))
        return compressed, uncompressed

    def _get_isoforms_cleaned_flag(self, accession: str) -> bool:
        try:
            row = self.db_manager.proteomes.get_default_cleaned_profile(str(accession))
            if row:
                status = str(row[6] or "ready").strip().lower()
                path = self.db_manager.proteomes.resolve_path(row)
                if status == "ready" and path and os.path.exists(path):
                    return True
        except Exception as exc:  # boundary: cleaned-profile lookup is optional for stale-state heuristics.
            self.log(f"Failed to inspect default cleaned proteome profile for {accession}: {exc}", "WARNING")
        try:
            self.db_manager.cursor.execute(
                "SELECT isoforms_cleaned FROM Genome WHERE accession = ?",
                (accession,),
            )
            row = self.db_manager.cursor.fetchone()
            if not row:
                return False
            return bool(row[0])
        except Exception as exc:  # boundary: isoforms-cleaned fallback is optional for stale-state heuristics.
            self.log(f"Failed to inspect isoforms_cleaned flag for {accession}: {exc}", "WARNING")
            return False

    def _strip_ext(self, path: str, base_ext: str) -> str:
        base = os.path.basename(path)
        if base.endswith(".gz"):
            base = base[:-3]
        if base.endswith(base_ext):
            base = base[: -len(base_ext)]
        return base

    def _predict_split_accessions(
        self,
        accession: str,
        nuc_files: List[str],
        prot_files: List[str],
        split_isolated_proteomes: bool,
    ) -> List[str]:
        nuc = sorted(os.path.basename(p) for p in nuc_files)
        prot = sorted(os.path.basename(p) for p in prot_files)
        if len(nuc) + len(prot) <= 1:
            return []

        prot_by_base: Dict[str, List[str]] = {}
        for p in prot:
            base = p.replace(".gz", "").replace(".faa", "")
            prot_by_base.setdefault(base, []).append(p)

        splits: List[str] = []
        for idx, nuc_name in enumerate(nuc[1:], start=1):
            new_acc = f"{accession}.{idx}"
            splits.append(new_acc)
            base = nuc_name.replace(".gz", "").replace(".fna", "")
            prot_by_base.pop(base, None)

        if split_isolated_proteomes:
            for _base in list(prot_by_base.keys()):
                new_acc = f"{accession}.{len(splits) + 1}"
                splits.append(new_acc)
                prot_by_base.pop(_base, None)
        return splits

    def _queue_reacquire(self, accessions: List[str]):
        if not accessions:
            return
        parent_id = self.task_id
        data = {
            "accessions": accessions,
            "protein": True,
            "max_concurrent": 1,
            "force_redownload": True,
        }
        self.queue_subtask(job_type=2, status="P", priority=1, data=data)
        self.log(f"Queued reacquire download for {len(accessions)} accessions", "INFO")

    def _process_assembly_accession(self, acc: str, genome_dir: str, split_isolated: bool) -> Dict[str, Any]:
        folder = self.db_manager.genomes.get_path(acc) or os.path.join(genome_dir, acc)
        folder_exists = os.path.isdir(folder)
        status_before = self.db_manager.genomes.get_status(acc)
        genome_before = self.db_manager.genomes.get(acc)
        col_map = self.db_manager.genomes.get_column_map()
        protein_before = None
        if genome_before and "protein" in col_map:
            try:
                protein_before = genome_before[col_map["protein"]]
            except (TypeError, IndexError):
                protein_before = None
        isoforms_cleaned_before = self._get_isoforms_cleaned_flag(acc)
        if self.repair:
            try:
                self.db_manager.genomes.update_status(acc, 0)
            except Exception as exc:  # boundary: repair pre-mark failure is logged and verification continues.
                self.log(f"{acc}: failed to pre-mark assembly inactive during repair: {exc}", "WARNING")

        nuc_gz: List[str] = []
        nuc_plain: List[str] = []
        prot_gz: List[str] = []
        prot_plain: List[str] = []
        other_files: List[str] = []
        if folder_exists:
            nuc_gz, nuc_plain = self._find_files(folder, ".fna")
            prot_gz, prot_plain = self._find_files(folder, ".faa")
            for fname in os.listdir(folder):
                low = fname.lower()
                if low.endswith((".fna", ".fna.gz", ".faa", ".faa.gz", ".gff", ".gff.gz", ".gff3", ".gff3.gz")):
                    continue
                other_files.append(fname)

        actions: List[str] = []
        reacquire_reasons: List[str] = []
        row: Dict[str, str] = {
            "accession": acc,
            "folder": folder,
            "folder_exists": "yes" if folder_exists else "no",
            "gff_files": "",
            "default_proteome_profile": "",
            "ready_proteome_profiles": "",
            "stale_proteome_profiles": "",
            "proteome_profile_status": "missing",
            "other_files": self._assembly_join(other_files),
            "organise_action": "none",
            "organise_new_accessions": "",
            "discover_marked": "no",
            "reacquire_queued": "no",
        }

        nuc_bases = {self._strip_ext(p, ".fna") for p in (nuc_gz + nuc_plain)}
        prot_bases = {self._strip_ext(p, ".faa") for p in (prot_gz + prot_plain)}
        multi_nuc = len(nuc_bases) > 1
        multi_prot = len(prot_bases) > 1
        row["multi_nuc_basenames"] = self._assembly_join(list(nuc_bases)) if multi_nuc else ""
        row["multi_prot_basenames"] = self._assembly_join(list(prot_bases)) if multi_prot else ""
        if multi_nuc or multi_prot:
            self.log(f"{acc}: multiple fasta files found (nuc {len(nuc_bases)}, prot {len(prot_bases)})", "WARNING")
            row["organise_new_accessions"] = self._assembly_join(
                self._predict_split_accessions(
                    acc,
                    nuc_gz + nuc_plain,
                    prot_gz + prot_plain,
                    split_isolated,
                )
            )
            if self.organise and self.repair and not self.organise_check_only:
                self.queue_subtask(
                    job_type=19,
                    status="P",
                    priority=1,
                    data={
                        "accession": acc,
                        "folder": folder,
                        "split_isolated_proteomes": split_isolated,
                        "check_only": False,
                    },
                )
                row["organise_action"] = "queued"
                actions.append("organise_queued")
                self._mark_assembly_checks_skipped(row)
                row["nuc_gz_total"] = str(len(nuc_gz))
                row["nuc_plain_total"] = str(len(nuc_plain))
                row["prot_gz_total"] = str(len(prot_gz))
                row["prot_plain_total"] = str(len(prot_plain))
                row["nuc_basenames"] = self._assembly_join(list(nuc_bases))
                row["prot_basenames"] = self._assembly_join(list(prot_bases))
                finalized = self._finalize_assembly_row(
                    row,
                    actions,
                    reacquire_reasons,
                    status_before=status_before,
                    protein_before=protein_before,
                    isoforms_cleaned_before=isoforms_cleaned_before,
                )
                finalized.update({"to_reacquire": [], "forced_clean_target": False})
                return finalized
            if self.organise_check_only:
                row["organise_action"] = "check_only"
                actions.append("organise_check_only")
                self._mark_assembly_checks_skipped(row)
                row["nuc_gz_total"] = str(len(nuc_gz))
                row["nuc_plain_total"] = str(len(nuc_plain))
                row["prot_gz_total"] = str(len(prot_gz))
                row["prot_plain_total"] = str(len(prot_plain))
                row["nuc_basenames"] = self._assembly_join(list(nuc_bases))
                row["prot_basenames"] = self._assembly_join(list(prot_bases))
                finalized = self._finalize_assembly_row(
                    row,
                    actions,
                    reacquire_reasons,
                    status_before=status_before,
                    protein_before=protein_before,
                    isoforms_cleaned_before=isoforms_cleaned_before,
                )
                finalized.update({"to_reacquire": [], "forced_clean_target": False})
                return finalized

        nuc_ok = False
        prot_ok = False
        bad_nuc: List[str] = []
        bad_prot: List[str] = []

        if nuc_gz:
            good = [p for p in nuc_gz if self._validate_gzip(p)]
            bad_nuc = [p for p in nuc_gz if p not in good]
            for b in bad_nuc:
                if self.repair:
                    try:
                        os.remove(b)
                    except OSError as exc:
                        self.log(f"{acc}: failed to remove bad nucleotide gzip {b}: {exc}", "WARNING")
                    actions.append(f"removed_bad_nuc_gz:{os.path.basename(b)}")
                    self._append_comment(acc, f"Removed bad nucleotide gzip: {b}")
            nuc_gz = good
            nuc_ok = bool(good)
        elif nuc_plain:
            nuc_ok = True
            if self.tidy and self.repair:
                gz_created = self._compress_if_needed(nuc_plain[0])
                if gz_created and self._validate_gzip(gz_created):
                    nuc_gz.append(gz_created)
                    actions.append(f"compressed_nuc:{os.path.basename(gz_created)}")
                    self._append_comment(acc, f"Compressed nucleotide to {gz_created}")

        if self.tidy and self.repair and nuc_gz and nuc_plain:
            nuc_gz_bases = {self._strip_ext(p, ".fna") for p in nuc_gz}
            kept = []
            for p in nuc_plain:
                if self._strip_ext(p, ".fna") in nuc_gz_bases:
                    try:
                        os.remove(p)
                    except OSError as exc:
                        self.log(f"{acc}: failed to remove uncompressed nucleotide {p}: {exc}", "WARNING")
                    actions.append(f"removed_uncompressed_nuc:{os.path.basename(p)}")
                    self._append_comment(acc, f"Removed uncompressed nucleotide {p}")
                else:
                    kept.append(p)
            nuc_plain = kept

        if prot_gz:
            goodp = [p for p in prot_gz if self._validate_gzip(p)]
            bad_prot = [p for p in prot_gz if p not in goodp]
            for b in bad_prot:
                if self.repair:
                    try:
                        os.remove(b)
                    except OSError as exc:
                        self.log(f"{acc}: failed to remove bad protein gzip {b}: {exc}", "WARNING")
                    actions.append(f"removed_bad_prot_gz:{os.path.basename(b)}")
                    self._append_comment(acc, f"Removed bad protein gzip: {b}")
            prot_gz = goodp
            prot_ok = bool(goodp)
        elif prot_plain:
            prot_ok = True
            if self.tidy and self.repair:
                gz_created = self._compress_if_needed(prot_plain[0])
                if gz_created and self._validate_gzip(gz_created):
                    prot_gz.append(gz_created)
                    actions.append(f"compressed_prot:{os.path.basename(gz_created)}")
                    self._append_comment(acc, f"Compressed protein to {gz_created}")

        if self.tidy and self.repair and prot_gz and prot_plain:
            prot_gz_bases = {self._strip_ext(p, ".faa") for p in prot_gz}
            kept = []
            for p in prot_plain:
                if self._strip_ext(p, ".faa") in prot_gz_bases:
                    try:
                        os.remove(p)
                    except OSError as exc:
                        self.log(f"{acc}: failed to remove uncompressed protein {p}: {exc}", "WARNING")
                    actions.append(f"removed_uncompressed_prot:{os.path.basename(p)}")
                    self._append_comment(acc, f"Removed uncompressed protein {p}")
                else:
                    kept.append(p)
            prot_plain = kept

        gff_files = self._find_gff_files(folder) if folder_exists else []
        if self.tidy and self.repair and gff_files:
            gff_gz = [p for p in gff_files if p.endswith(".gz")]
            gff_plain = [p for p in gff_files if not p.endswith(".gz")]
            for p in list(gff_plain):
                gz_created = self._compress_if_needed(p)
                if gz_created and self._validate_gzip(gz_created):
                    gff_gz.append(gz_created)
                    actions.append(f"compressed_gff:{os.path.basename(gz_created)}")
                    self._append_comment(acc, f"Compressed GFF to {gz_created}")
            gff_gz_bases = {self._strip_gff_ext(p) for p in gff_gz}
            kept = []
            for p in gff_plain:
                if self._strip_gff_ext(p) in gff_gz_bases:
                    try:
                        os.remove(p)
                    except OSError as exc:
                        self.log(f"{acc}: failed to remove uncompressed GFF {p}: {exc}", "WARNING")
                    actions.append(f"removed_uncompressed_gff:{os.path.basename(p)}")
                    self._append_comment(acc, f"Removed uncompressed GFF {p}")
                else:
                    kept.append(p)
            gff_files = sorted(gff_gz + kept)

        forced_clean_target = False

        row["nuc_gz_total"] = str(len(nuc_gz) + len(bad_nuc))
        row["nuc_plain_total"] = str(len(nuc_plain))
        row["prot_gz_total"] = str(len(prot_gz) + len(bad_prot))
        row["prot_plain_total"] = str(len(prot_plain))
        row["nuc_basenames"] = self._assembly_join(list({self._strip_ext(p, ".fna") for p in (nuc_gz + nuc_plain)}))
        row["prot_basenames"] = self._assembly_join(list({self._strip_ext(p, ".faa") for p in (prot_gz + prot_plain)}))
        row["nuc_gz_valid"] = str(len(nuc_gz))
        row["nuc_gz_invalid"] = str(len(bad_nuc))
        row["prot_gz_valid"] = str(len(prot_gz))
        row["prot_gz_invalid"] = str(len(bad_prot))
        row["nuc_ok"] = "yes" if nuc_ok else "no"
        row["nuc_missing_or_invalid"] = "no" if nuc_ok else "yes"
        row["prot_ok"] = "yes" if prot_ok else "no"
        row["prot_missing_or_invalid"] = "no" if prot_ok else "yes"
        row["gff_files"] = self._assembly_join([os.path.basename(p) for p in gff_files])

        artifact_statuses = self._reconcile_owner_artifacts(
            "genome",
            acc,
            {
                "genome_fna": sorted(nuc_gz + nuc_plain)[0] if (nuc_gz or nuc_plain) else None,
                "genome_faa": sorted(prot_gz + prot_plain)[0] if (prot_gz or prot_plain) else None,
                "genome_gff": gff_files[0] if gff_files else None,
                "genome_ncbi_md5checksums": os.path.join(folder, "md5checksums.txt")
                if folder_exists and os.path.exists(os.path.join(folder, "md5checksums.txt"))
                else None,
                "genome_phyloodb_md5checksums": os.path.join(folder, "phyloodb_md5checksums.txt")
                if folder_exists and os.path.exists(os.path.join(folder, "phyloodb_md5checksums.txt"))
                else None,
            },
            metadata={"accession": acc},
            repair=self.repair,
        )
        profile_sync = _sync_profile_state(self, acc, folder) if self.repair and folder_exists else _inspect_profile_state(self, acc, folder)
        row["default_proteome_profile"] = str(profile_sync.get("default_profile") or "")
        row["ready_proteome_profiles"] = self._assembly_join([str(name) for name in (profile_sync.get("ready_profiles") or [])])
        row["stale_proteome_profiles"] = self._assembly_join([str(name) for name in (profile_sync.get("stale_profiles") or [])])
        row["proteome_profile_status"] = str(profile_sync.get("status") or "missing")

        for artifact_type, artifact_status in artifact_statuses.items():
            if artifact_status != "ready":
                actions.append(f"{artifact_type}:{artifact_status}")

        to_reacquire: List[str] = []
        if self.discover and self.repair and nuc_ok:
            try:
                protein_flag = bool(profile_sync.get("has_protein")) if self.discover_protein else False
                self.db_manager.genomes.update_status(acc, 1, (datetime.now(), folder, protein_flag))
                self._append_comment(acc, "Marked as downloaded via discover")
            except Exception as exc:  # boundary: discover status update failure is logged after repair inspection.
                self.log(f"{acc}: failed to mark discovered assembly downloaded: {exc}", "WARNING")
            row["discover_marked"] = "yes"
            actions.append("discover_marked")
            finalized = self._finalize_assembly_row(
                row,
                actions,
                reacquire_reasons,
                status_before=status_before,
                protein_before=protein_before,
                isoforms_cleaned_before=isoforms_cleaned_before,
                artifact_statuses=artifact_statuses,
                profile_state=profile_sync,
            )
            finalized.update({"to_reacquire": to_reacquire, "forced_clean_target": forced_clean_target})
            return finalized

        if not nuc_ok:
            to_reacquire.append(acc)
            if self.repair:
                self._append_comment(acc, "Nucleotide file missing or invalid; flagged for reacquire")
            reacquire_reasons.append("nuc_missing_or_invalid")
            if self.repair:
                try:
                    if prot_ok:
                        self.db_manager.genomes.update_status(acc, 1, (datetime.now(), folder, True))
                    else:
                        self.db_manager.genomes.update_status(acc, 0)
                except Exception as exc:  # boundary: repair status update failure is logged; reacquire can still be queued.
                    self.log(f"{acc}: failed to update status for invalid nucleotide repair: {exc}", "WARNING")
            if self.repair and self.reacquire:
                row["reacquire_queued"] = "yes"
                actions.append("reacquire_queued")
            finalized = self._finalize_assembly_row(
                row,
                actions,
                reacquire_reasons,
                status_before=status_before,
                protein_before=protein_before,
                isoforms_cleaned_before=isoforms_cleaned_before,
                artifact_statuses=artifact_statuses,
                profile_state=profile_sync,
            )
            finalized.update({"to_reacquire": to_reacquire, "forced_clean_target": forced_clean_target})
            return finalized

        if self.repair:
            try:
                self.db_manager.genomes.update_status(acc, 1, (datetime.now(), folder, bool(prot_ok)))
            except Exception as exc:  # boundary: repair status update failure is logged after file checks.
                self.log(f"{acc}: failed to refresh verified assembly status: {exc}", "WARNING")

        if not prot_ok:
            if self.repair:
                self._append_comment(acc, "Protein file missing or invalid")
            reacquire_reasons.append("prot_missing_or_invalid")
            if self.repair and self.reacquire:
                to_reacquire.append(acc)
                row["reacquire_queued"] = "yes"
                actions.append("reacquire_queued")

        finalized = self._finalize_assembly_row(
            row,
            actions,
            reacquire_reasons,
            status_before=status_before,
            protein_before=protein_before,
            isoforms_cleaned_before=isoforms_cleaned_before,
            artifact_statuses=artifact_statuses,
            profile_state=profile_sync,
        )
        finalized.update({"to_reacquire": to_reacquire, "forced_clean_target": forced_clean_target})
        return finalized

    def run(self):
        if self.stage >= 1:
            clean_targets = list(self.data.get("_verify_clean_targets", []) or [])
            # Re-check DB state on resume so stale checkpoint payloads do not re-clean already-cleaned assemblies.
            clean_targets = [acc for acc in clean_targets if not self._get_isoforms_cleaned_flag(acc)]
            if not clean_targets:
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
                self.log(f"Queued prepare-proteome subtask for {len(clean_targets)} verified assemblies (including continuity fixes).", "INFO")
                return True

            outcome = self.manage_subtasks(
                stage=1,
                queue_fn=queue_clean_subtask,
                done_fn=None,
                wait_seconds=0,
                retry_key=None,
                max_retries=0,
                incomplete_message_fn=lambda: ("Verify prepare-proteome subtask did not complete.", ""),
                retry_incomplete=False,
            )
            if outcome == "ERROR":
                return "ERROR"
            if outcome is False:
                return False
            return True

        try:
            genome_dir = self.db_manager.storage.require_root_base("genomes")
        except ValueError as exc:
            return self.handle_exception(str(exc), {})

        try:
            selected = self.prepare_selectors(
                target_key="_verify_selector_accessions",
                taxid=self.data.get("taxid"),
                additional=self.data.get("accessions"),
                root=self.data.get("root"),
                allow_all=self.all or bool(self.data.get("root")),
                downloaded_only=self.downloaded_only,
                released_after=self.data.get("after"),
                released_before=self.data.get("before"),
                level=self.data.get("level"),
                primary_only=self.data.get("primary_only"),
                require_candidates=not bool(self.discover),
            )
        except ValueError as exc:
            return self.handle_exception(str(exc), {})

        to_reacquire: List[str] = []
        report_rows: List[Dict[str, str]] = []
        summary_rows: List[Dict[str, str]] = []
        report_headers = [
            "accession",
            "folder",
            "folder_exists",
            "nuc_gz_total",
            "nuc_plain_total",
            "nuc_basenames",
            "nuc_gz_valid",
            "nuc_gz_invalid",
            "nuc_ok",
            "nuc_missing_or_invalid",
            "prot_gz_total",
            "prot_plain_total",
            "prot_basenames",
            "prot_gz_valid",
            "prot_gz_invalid",
            "prot_ok",
            "prot_missing_or_invalid",
            "multi_nuc_basenames",
            "multi_prot_basenames",
            "gff_files",
            "default_proteome_profile",
            "ready_proteome_profiles",
            "stale_proteome_profiles",
            "proteome_profile_status",
            "other_files",
            "organise_action",
            "organise_new_accessions",
            "discover_marked",
            "reacquire_queued",
            "reacquire_reason",
            "actions",
        ]
        summary_headers = [
            "accession",
            "folder",
            "folder_exists",
            "status_before",
            "status_after",
            "status_changed",
            "assembly_usable",
            "proteome_usable",
            "protein_before",
            "protein_after",
            "isoforms_cleaned_before",
            "isoforms_cleaned_after",
            "assembly_stale",
            "proteome_stale",
            "fna_artifact_status",
            "faa_artifact_status",
            "gff_artifact_status",
            "default_proteome_profile",
            "ready_proteome_profiles",
            "stale_proteome_profiles",
            "proteome_profile_status",
            "deactivated",
            "critical_reasons",
            "artifact_warnings",
            "reacquire_reason",
            "actions",
        ]
        split_isolated = bool(self.data.get("split_isolated_proteomes", False))
        forced_clean_targets = set()
        worker_results, worker_errors = _run_verify_workers(
            self,
            selected,
            lambda acc: self._run_assembly_worker(acc, genome_dir, split_isolated),
            "assembly",
        )
        for result in worker_results:
            if not result:
                continue
            report_rows.append(result["report_row"])
            if result.get("summary_row"):
                summary_rows.append(result["summary_row"])
            to_reacquire.extend(result.get("to_reacquire", []))
            if result.get("forced_clean_target"):
                forced_clean_targets.add(result["report_row"]["accession"])

        if worker_errors:
            to_reacquire = []

        if self.repair and self.reacquire and to_reacquire:
            self._queue_reacquire(to_reacquire)

        clean_targets_to_queue = []
        if not worker_errors and self.repair and ((self.clean_isoforms and not self.skip_clean_isoforms) or forced_clean_targets):
            clean_candidates = [
                row["accession"]
                for row in report_rows
                if row.get("folder_exists") == "yes" and row.get("prot_ok") == "yes"
            ]
            not_cleaned = [acc for acc in clean_candidates if not self._get_isoforms_cleaned_flag(acc)]
            clean_targets = sorted(set(not_cleaned).union(forced_clean_targets))
            if clean_targets:
                clean_targets_to_queue = clean_targets

        try:
            _write_tsv_report(self, "verify_assembly.tsv", report_headers, report_rows)
            _write_tsv_report(self, "verify_assembly_summary.tsv", summary_headers, summary_rows)
        except (OSError, UnicodeError) as exc:
            self.log(f"Failed to write report: {exc}", "ERROR")

        if clean_targets_to_queue:
            self.data["_verify_clean_targets"] = clean_targets_to_queue
            def queue_clean_subtask():
                clean_max_concurrent = self.clean_max_concurrent
                clean_threads_per_job = self.clean_threads_per_job
                parent_threads_budget = self.REQUIRED_THREADS
                if len(clean_targets_to_queue) > 10:
                    clean_max_concurrent = 8
                    clean_threads_per_job = 1
                    parent_threads_budget = 8
                self.queue_subtask(
                    job_type=31,
                    status="P",
                    priority=1,
                    data={
                        "accessions": clean_targets_to_queue,
                        "downloaded_only": True,
                        "skip_gff": self.clean_skip_gff,
                        "skip_cdhit": self.clean_skip_cdhit,
                        "gff_priority": self.clean_gff_priority,
                        "cdhit_identity": self.clean_cdhit_identity,
                        "profile_name": "clean_default",
                        "input_profile": "raw",
                        "set_default": True,
                        "max_concurrent": clean_max_concurrent,
                        "threads_per_job": clean_threads_per_job,
                        "parent_threads_budget": parent_threads_budget,
                    },
                )
                self.log(f"Queued prepare-proteome subtask for {len(clean_targets_to_queue)} verified assemblies (including continuity fixes).", "INFO")
                return True

            outcome = self.manage_subtasks(
                stage=1,
                queue_fn=queue_clean_subtask,
                done_fn=None,
                wait_seconds=0,
                retry_key=None,
                max_retries=0,
                incomplete_message_fn=lambda: ("Verify prepare-proteome subtask did not complete.", ""),
                retry_incomplete=False,
            )
            if outcome == "ERROR":
                return "ERROR"
            if outcome is False:
                return False

        if worker_errors:
            return self.handle_exception(_parallel_error_message("verify-assembly worker failures", worker_errors))
        return True


class _ArtifactVerificationMixin:
    def _artifact_path_exists(self, artifact_row) -> tuple[bool, Optional[str], Optional[int]]:
        resolved = self.db_manager.artifacts.resolve_path(artifact_row)
        if not resolved:
            return False, None, None
        is_dir = bool(artifact_row[9]) if len(artifact_row) > 9 else False
        exists = os.path.isdir(resolved) if is_dir else os.path.exists(resolved)
        if not exists:
            return False, resolved, None
        if is_dir:
            return True, resolved, None
        try:
            size_bytes = os.path.getsize(resolved)
        except OSError:
            size_bytes = None
        return True, resolved, size_bytes

    def _verify_artifact_rows(
        self,
        artifact_rows: List[tuple],
        *,
        repair: bool = False,
        stale_missing: bool = True,
        restore_found: bool = True,
    ) -> dict[str, Any]:
        checked = 0
        stale = 0
        restored = 0
        ready = 0
        staled_updates = 0
        ready_size_updates = 0
        missing_rows: list[tuple[int, str, str, str]] = []
        for row in artifact_rows or []:
            checked += 1
            artifact_id = int(row[0])
            owner_type = str(row[1])
            owner_id = str(row[2])
            artifact_type = str(row[3])
            current_status = str(row[5] or "ready")
            exists, resolved, size_bytes = self._artifact_path_exists(row)
            target_status = "ready" if exists else "stale"
            if exists:
                ready += 1
                if repair and restore_found and current_status != target_status:
                    self.db_manager.artifacts.set_status(artifact_id, "ready", size_bytes=size_bytes)
                    restored += 1
            else:
                stale += 1
                missing_rows.append((artifact_id, owner_type, owner_id, artifact_type))
                if repair and stale_missing and current_status != target_status:
                    self.db_manager.artifacts.set_status(artifact_id, "stale")
                    staled_updates += 1
            if repair and exists and current_status == "ready" and size_bytes is not None:
                self.db_manager.artifacts.set_status(artifact_id, "ready", size_bytes=size_bytes)
                ready_size_updates += 1
        return {
            "checked": checked,
            "ready": ready,
            "stale": stale,
            "restored": restored,
            "staled_updates": staled_updates,
            "ready_size_updates": ready_size_updates,
            "missing_rows": missing_rows,
        }

class VerifyBuscoTask(_ArtifactVerificationMixin, Task):
    """Verify BUSCO results: discover existing outputs, re-ingest if anomalous/missing, optionally queue runs."""

    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads=required_threads)
        self.all = bool(self.data.get("all", False))
        self.discover = bool(self.data.get("discover", False))
        self.queue_missing = bool(self.data.get("queue_missing", False))
        self.reingest = bool(self.data.get("reingest", False))
        self.reingest_all = bool(self.data.get("reingest_all", False))
        self.repair = bool(self.data.get("repair", False))
        self.verify_artifacts = True
        self.stale_missing = bool(self.data.get("stale_missing", True))
        self.restore_found = bool(self.data.get("restore_found", True))
        self.reassign_primary = bool(self.data.get("reassign_primary", True))
        self.library_id = self.data.get("library_id")
        self.library_name = self.data.get("library_name")
        self.run_id = self.data.get("run_id")
        self.downloaded_only = not bool(self.discover)
        self.report = self.data.get("report")

    def _run_busco_worker(self, accession: str, genome_root: str) -> Dict[str, Any]:
        worker = VerifyBuscoTask(
            self.db_manager.get_path(),
            self.task_id,
            data=json.dumps(self.data),
            required_threads=1,
        )
        try:
            return worker._process_busco_accession(accession, genome_root)
        finally:
            worker.db_manager.close()

    def _canonical_busco_path(self, path: Optional[str]) -> Optional[str]:
        if not path:
            return None
        try:
            return os.path.realpath(os.path.abspath(str(path)))
        except (TypeError, ValueError, OSError):
            return os.path.abspath(str(path))

    def _path_within(self, path: Optional[str], root: Optional[str]) -> bool:
        canonical_path = self._canonical_busco_path(path)
        canonical_root = self._canonical_busco_path(root)
        if not canonical_path or not canonical_root:
            return False
        try:
            return os.path.commonpath([canonical_path, canonical_root]) == canonical_root
        except ValueError:
            return False

    def _busco_run_folders(self, genome_path: str, lineage_name: str) -> List[str]:
        candidates = []
        legacy = os.path.join(genome_path, f"{lineage_name}_results")
        if os.path.isdir(legacy):
            candidates.append(legacy)
        candidates.extend(
            sorted(glob.glob(os.path.join(genome_path, f"{lineage_name}_results__*")), reverse=True)
        )
        found: List[str] = []
        seen: set[str] = set()
        for base in candidates:
            runs = sorted(glob.glob(os.path.join(base, "run_*")), reverse=True)
            for run_dir in runs:
                canonical = self._canonical_busco_path(run_dir)
                if (
                    canonical
                    and os.path.isdir(canonical)
                    and self._path_within(canonical, genome_path)
                    and canonical not in seen
                ):
                    seen.add(canonical)
                    found.append(canonical)
        return found

    def _all_busco_run_folders(self, genome_path: str) -> Dict[str, List[str]]:
        found: Dict[str, List[str]] = {}
        seen: Dict[str, set[str]] = defaultdict(set)
        pattern_roots = sorted(
            set(
                glob.glob(os.path.join(genome_path, "*_results"))
                + glob.glob(os.path.join(genome_path, "*_results__*"))
            )
        )
        for base in pattern_roots:
            if not os.path.isdir(base):
                continue
            base_name = os.path.basename(base)
            if "_results" not in base_name:
                continue
            lineage_name = base_name.split("_results", 1)[0]
            run_dirs = sorted(glob.glob(os.path.join(base, "run_*")), reverse=True)
            for run_dir in run_dirs:
                canonical = self._canonical_busco_path(run_dir)
                if (
                    canonical
                    and os.path.isdir(canonical)
                    and self._path_within(canonical, genome_path)
                    and canonical not in seen[lineage_name]
                ):
                    seen[lineage_name].add(canonical)
                    found.setdefault(lineage_name, []).append(canonical)
        return found

    def _scoped_busco_targets(
        self,
        accession: str,
        genome_path: str,
    ) -> List[Tuple[int, str, List[str]]]:
        if self.library_id:
            parent_id = self.db_manager.libraries.get_parent_id(self.library_id)
            busco_lib_id = int(parent_id) if parent_id else int(self.library_id)
            busco_lineage_name = self.db_manager.libraries.get_name(busco_lib_id)
            if not busco_lineage_name:
                return []
            return [(busco_lib_id, str(busco_lineage_name), self._busco_run_folders(genome_path, str(busco_lineage_name)))]

        targets: Dict[int, Tuple[str, List[str]]] = {}
        discovered = self._all_busco_run_folders(genome_path)
        existing_run_rows = self.db_manager.busco.get_runs_for_accessions([accession], library_id=None)
        if self.run_id is not None:
            existing_run_rows = [row for row in existing_run_rows if int(row[0]) == int(self.run_id)]

        for row in existing_run_rows:
            lib_id = int(row[2])
            lineage_name = str(row[4] or row[3] or "")
            targets[lib_id] = (lineage_name, discovered.get(lineage_name, []))

        for lineage_name, run_dirs in discovered.items():
            lib_id = self.db_manager.libraries.get_id(lineage_name, include_inactive=True)
            if lib_id is None:
                self.log(f"{accession}: discovered BUSCO lineage folder {lineage_name} but no matching library exists in the DB.", "WARNING")
                continue
            targets[int(lib_id)] = (lineage_name, run_dirs)

        return [(lib_id, lineage_name, run_dirs) for lib_id, (lineage_name, run_dirs) in sorted(targets.items())]

    def _busco_full_table(self, run_dir: str) -> Optional[str]:
        if not run_dir or not os.path.isdir(run_dir):
            return None
        tables = glob.glob(os.path.join(run_dir, "full_table*.tsv"))
        if not tables:
            return None
        return tables[0]

    def _busco_summary_json(self, run_dir: str) -> Optional[str]:
        if not run_dir or not os.path.isdir(run_dir):
            return None
        result_dir = os.path.dirname(run_dir)
        candidates = []
        candidates.extend(sorted(glob.glob(os.path.join(result_dir, "short_summary*.json"))))
        candidates.extend(sorted(glob.glob(os.path.join(run_dir, "short_summary*.json"))))
        if not candidates:
            return None
        return candidates[0]

    def _busco_seq_path(self, run_dir: str, family_id: str, status: int) -> Optional[str]:
        if not run_dir:
            return None
        if status == 1:
            sub = "single_copy_busco_sequences"
        elif status == 2:
            sub = "multi_copy_busco_sequences"
        elif status == 3:
            sub = "fragmented_busco_sequences"
        else:
            return None
        for ext in (".faa", ".fna", ".fa", ".fasta"):
            cand = os.path.join(run_dir, "busco_sequences", sub, f"{family_id}{ext}")
            if os.path.isfile(cand):
                return cand
        return None

    def _infer_pipeline_from_run_dir(self, run_dir: str) -> str:
        base = os.path.basename(run_dir or "")
        # Expected naming now: run_<pipeline>_<lineage>...
        if base.startswith("run_"):
            token = base[4:].split("_", 1)[0].strip().lower()
            if token in {"miniprot", "metaeuk", "augustus"}:
                return token
        return "miniprot"

    def _infer_input_mode_from_run_dir(self, run_dir: str, *, pipeline: Optional[str] = None) -> str:
        seq_root = os.path.join(run_dir, "busco_sequences")
        pipeline_name = str(pipeline or "").strip().lower()
        if not os.path.isdir(seq_root):
            return "genome" if pipeline_name in {"miniprot", "metaeuk", "augustus"} else "protein"
        saw_protein = False
        for root, _dirs, files in os.walk(seq_root):
            for fn in files:
                low = fn.lower()
                if low.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")):
                    return "genome"
                if low.endswith((".faa", ".pep", ".faa.gz", ".pep.gz")):
                    saw_protein = True
        if saw_protein and pipeline_name != "miniprot":
            return "protein"
        if pipeline_name in {"miniprot", "metaeuk", "augustus"}:
            return "genome"
        return "protein"

    def _sequence_kind_for_path(self, location: Optional[str]) -> Optional[str]:
        if not location:
            return None
        low = str(location).lower()
        if low.endswith((".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            return "nucl"
        if low.endswith((".faa", ".faa.gz", ".pep", ".pep.gz")):
            return "prot"
        return None

    def _input_mode_for_path(self, location: Optional[str]) -> Optional[str]:
        kind = self._sequence_kind_for_path(location)
        if kind == "prot":
            return "protein"
        if kind == "nucl":
            return "genome"
        return None

    def _summary_mode_to_input_mode(self, mode: Optional[str]) -> Optional[str]:
        if not mode:
            return None
        low = str(mode).strip().lower()
        if low in {"protein", "proteins"}:
            return "protein"
        if low in {"genome", "genomes", "transcriptome", "transcriptomes", "nucl", "nucleotide", "nucleotides"}:
            return "genome"
        return None

    def _summary_flags_to_pipeline(self, params: dict) -> Optional[str]:
        true_flags = []
        for flag, name in (
            ("use_metaeuk", "metaeuk"),
            ("use_miniprot", "miniprot"),
            ("use_augustus", "augustus"),
        ):
            raw = params.get(flag)
            if isinstance(raw, str):
                value = raw.strip().lower() == "true"
            else:
                value = bool(raw)
            if value:
                true_flags.append(name)
        if len(true_flags) == 1:
            return true_flags[0]
        return None

    def _extract_summary_counts(self, results: Any) -> Optional[Dict[str, int]]:
        if not isinstance(results, dict):
            return None
        mapping = {
            "Single copy BUSCOs": "Single copy BUSCOs",
            "Multi copy BUSCOs": "Multi copy BUSCOs",
            "Fragmented BUSCOs": "Fragmented BUSCOs",
            "Missing BUSCOs": "Missing BUSCOs",
        }
        counts: Dict[str, int] = {}
        for source_key, target_key in mapping.items():
            if source_key not in results:
                return None
            try:
                counts[target_key] = int(results[source_key])
            except (TypeError, ValueError):
                try:
                    counts[target_key] = int(float(results[source_key]))
                except (TypeError, ValueError):
                    return None
        return counts

    def _run_counts_from_row(self, run_row: Optional[tuple]) -> Optional[Dict[str, int]]:
        if not run_row:
            return None
        try:
            return {
                "Single copy BUSCOs": int(run_row[9] or 0),
                "Multi copy BUSCOs": int(run_row[10] or 0),
                "Fragmented BUSCOs": int(run_row[11] or 0),
                "Missing BUSCOs": int(run_row[12] or 0),
            }
        except (TypeError, ValueError, IndexError):
            return None

    def _counts_equal(self, left: Optional[Dict[str, int]], right: Optional[Dict[str, int]]) -> bool:
        if left is None or right is None:
            return False
        return all(int(left.get(key, -1)) == int(right.get(key, -2)) for key in (
            "Single copy BUSCOs",
            "Multi copy BUSCOs",
            "Fragmented BUSCOs",
            "Missing BUSCOs",
        ))

    def _parse_busco_table_counts(self, accession: str, run_dir: str, library_id: int) -> tuple[Optional[Dict[str, int]], list[str]]:
        table = self._busco_full_table(run_dir)
        notes: list[str] = []
        if not table:
            return None, ["full_table_missing"]
        counts = {"Single copy BUSCOs": 0, "Multi copy BUSCOs": 0, "Fragmented BUSCOs": 0, "Missing BUSCOs": 0}
        seen_family_status: Dict[str, int] = {}
        threshold = 1e5
        try:
            with open(table, "r", encoding="utf-8") as handle:
                reader = csv.reader(handle, delimiter="\t")
                for row in reader:
                    if not row or row[0].startswith("#") or len(row) < 2:
                        continue
                    family_id = row[0]
                    status_str = row[1]
                    status_map = {"Complete": 1, "Duplicated": 2, "Fragmented": 3, "Missing": 4}
                    status = status_map.get(status_str, 0)
                    if family_id not in seen_family_status:
                        counts_key = {
                            1: "Single copy BUSCOs",
                            2: "Multi copy BUSCOs",
                            3: "Fragmented BUSCOs",
                            4: "Missing BUSCOs",
                        }.get(status)
                        if counts_key:
                            counts[counts_key] += 1
                        seen_family_status[family_id] = status
                    if len(row) >= 7 and row[6]:
                        try:
                            score = float(row[6])
                            if score > threshold:
                                notes.append("full_table_anomalous_score")
                        except (TypeError, ValueError):
                            self.log(f"Ignoring non-numeric BUSCO full_table score in {full_table}: {row}", "DEBUG")
        except (OSError, UnicodeError, csv.Error) as exc:
            return None, [f"full_table_unreadable:{type(exc).__name__}"]
        return counts, notes

    def _owner_artifact_path(self, owner_type: str, owner_id: str, artifact_type: str) -> Optional[str]:
        rows = self.db_manager.artifacts.find(owner_type=owner_type, owner_id=owner_id, artifact_type=artifact_type, status="ready")
        for row in rows or []:
            resolved = self.db_manager.artifacts.resolve_path(row)
            if resolved:
                return resolved
        rows = self.db_manager.artifacts.find(owner_type=owner_type, owner_id=owner_id, artifact_type=artifact_type)
        for row in rows or []:
            resolved = self.db_manager.artifacts.resolve_path(row)
            if resolved:
                return resolved
        return None

    def _resolve_profile_row_from_staged_busco_input(
        self,
        accession: str,
        summary_input_path: Optional[str],
    ) -> Optional[tuple]:
        profile_name = staged_busco_input_profile_name(summary_input_path)
        if not profile_name:
            return None
        if profile_name == DEFAULT_CLEAN_PROFILE:
            default_cleaned = self.db_manager.proteomes.get_default_cleaned_profile(str(accession))
            if default_cleaned:
                return default_cleaned
        return self.db_manager.proteomes.get_profile(str(accession), profile_name)

    def _resolve_proteome_profile_row_for_run(
        self,
        accession: str,
        *,
        input_mode: Optional[str],
        summary_input_path: Optional[str] = None,
        matched_run_row: Optional[tuple] = None,
    ) -> tuple[Optional[tuple], list[str]]:
        notes: list[str] = []
        if str(input_mode or "").strip().lower() != "protein":
            return None, notes
        row = None
        matched_from_staged_input = False
        if summary_input_path:
            row = self.db_manager.proteomes.find_profile_by_path(str(accession), str(summary_input_path))
        if row is None and summary_input_path:
            row = self._resolve_profile_row_from_staged_busco_input(str(accession), str(summary_input_path))
            matched_from_staged_input = row is not None
        if row is None and matched_run_row and len(matched_run_row) > 20 and matched_run_row[20]:
            matched_profile_name = str(matched_run_row[20] or "")
            matched_profile_row = self.db_manager.proteomes.get_profile(str(accession), matched_profile_name)
            matched_profile_status = str(matched_profile_row[6] or "").strip().lower() if matched_profile_row and len(matched_profile_row) > 6 else ""
            if matched_profile_row and matched_profile_status == "ready":
                row = matched_profile_row
            else:
                default_cleaned_row = self.db_manager.proteomes.get_default_cleaned_profile(str(accession))
                if default_cleaned_row:
                    default_cleaned_name = str(default_cleaned_row[2] or "")
                    if matched_profile_name == DEFAULT_CLEAN_PROFILE:
                        notes.append(f"legacy_clean_default_rebound:{default_cleaned_name}")
                        row = default_cleaned_row
                    elif matched_profile_row and matched_profile_status == "stale":
                        notes.append(f"stale_profile_rebound:{matched_profile_name}->{default_cleaned_name}")
                        row = default_cleaned_row
                    elif matched_profile_row is None:
                        row = None
                else:
                    row = matched_profile_row
        if row is None:
            row = self.db_manager.proteomes.get_default_profile(str(accession))
        if row is not None:
            row_status = str(row[6] or "").strip().lower() if len(row) > 6 else ""
            if row_status == "stale":
                default_cleaned_row = self.db_manager.proteomes.get_default_cleaned_profile(str(accession))
                if default_cleaned_row and int(default_cleaned_row[0]) != int(row[0]):
                    notes.append(f"stale_profile_rebound:{str(row[2] or '')}->{str(default_cleaned_row[2] or '')}")
                    row = default_cleaned_row
        if summary_input_path:
            if row is None:
                notes.append("summary_input_profile_unresolved")
            else:
                expected_path = self.db_manager.proteomes.resolve_path(row)
                if (
                    not matched_from_staged_input
                    and self._canonical_busco_path(summary_input_path) != self._canonical_busco_path(expected_path)
                ):
                    notes.append(f"summary_input_profile_mismatch:{str(row[2])}")
        return row, notes

    def _inspect_summary_input(self, accession: str, run_dir: str) -> Dict[str, Any]:
        summary_path = self._busco_summary_json(run_dir)
        result: Dict[str, Any] = {
            "summary_path": summary_path,
            "summary_input_path": None,
            "summary_mode_raw": None,
            "summary_input_mode": None,
            "summary_pipeline": None,
            "summary_counts": None,
            "path_input_mode": None,
            "verified_input_mode": None,
            "notes": [],
        }
        if not summary_path or not os.path.isfile(summary_path):
            return result
        try:
            payload = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            result["notes"].append(f"summary_json_unreadable:{type(exc).__name__}")
            return result
        params = payload.get("parameters") if isinstance(payload, dict) else None
        if not isinstance(params, dict):
            result["notes"].append("summary_json_missing_parameters")
            return result
        json_results = payload.get("results") if isinstance(payload, dict) else None
        input_path = params.get("in")
        mode_raw = params.get("mode")
        summary_mode = self._summary_mode_to_input_mode(mode_raw)
        summary_pipeline = self._summary_flags_to_pipeline(params)
        summary_counts = self._extract_summary_counts(json_results)
        path_mode = self._input_mode_for_path(input_path)
        inferred_pipeline = self._infer_pipeline_from_run_dir(run_dir)
        result["summary_input_path"] = input_path
        result["summary_mode_raw"] = mode_raw
        result["summary_input_mode"] = summary_mode
        result["summary_pipeline"] = summary_pipeline
        result["summary_counts"] = summary_counts
        result["path_input_mode"] = path_mode
        if summary_mode and path_mode and summary_mode != path_mode:
            result["notes"].append(f"summary_mode_path_mismatch:{summary_mode}:{path_mode}")
        verified_input_mode = summary_mode or path_mode or self._infer_input_mode_from_run_dir(
            run_dir,
            pipeline=summary_pipeline or inferred_pipeline,
        )
        if verified_input_mode:
            result["verified_input_mode"] = verified_input_mode
        expected_artifact_type = "genome_fna" if verified_input_mode == "genome" else None
        if expected_artifact_type and input_path:
            expected_path = self._owner_artifact_path("genome", accession, expected_artifact_type)
            if expected_path:
                input_canonical = self._canonical_busco_path(input_path)
                expected_canonical = self._canonical_busco_path(expected_path)
                if input_canonical != expected_canonical:
                    result["notes"].append(f"summary_input_artifact_mismatch:{expected_artifact_type}")
            else:
                result["notes"].append(f"summary_expected_artifact_missing:{expected_artifact_type}")
        return result

    def _register_busco_run_artifacts(
        self,
        run_id: int,
        accession: str,
        library_id: int,
        run_dir: str,
        family_locations: List[Tuple[str, int, str, Optional[str]]],
    ) -> None:
        result_dir = os.path.dirname(run_dir)
        self.db_manager.busco.register_run_artifact(run_id, "busco_result_root", result_dir, is_dir=True, format="directory")
        self.db_manager.busco.register_run_artifact(run_id, "busco_run_dir", run_dir, is_dir=True, format="directory")
        sequences_dir = os.path.join(run_dir, "busco_sequences")
        if os.path.isdir(sequences_dir):
            self.db_manager.busco.register_run_artifact(run_id, "busco_sequences_dir", sequences_dir, is_dir=True, format="directory")
        full_table = self._busco_full_table(run_dir)
        if full_table and os.path.isfile(full_table):
            self.db_manager.busco.register_run_artifact(run_id, "busco_full_table_tsv", full_table, format="tsv")
        summary_json = self._busco_summary_json(run_dir)
        if summary_json and os.path.isfile(summary_json):
            self.db_manager.busco.register_run_artifact(run_id, "busco_summary_json", summary_json, format="json")
        for family_id, _lib_id, _acc, location in family_locations:
            if not location or not os.path.exists(location):
                continue
            self.db_manager.busco.register_family_artifact(
                run_id=run_id,
                family_id=family_id,
                library_id=library_id,
                accession=accession,
                path=location,
                sequence_kind=self._sequence_kind_for_path(location),
                role=family_id,
                format="fasta",
                metadata={"source": "verify-busco"},
            )

    def _match_existing_run(
        self,
        run_rows: List[tuple],
        *,
        result_dir: str,
        pipeline: str,
        input_mode: str,
    ) -> Optional[Tuple[int, bool]]:
        target_dir = self._canonical_busco_path(result_dir)
        exact_matches: List[Tuple[Tuple[int, int, int, str, int], int]] = []
        canonical_matches: List[Tuple[Tuple[int, int, int, str, int], int]] = []
        for row in run_rows or []:
            row_id = int(row[0])
            row_pipeline = str(row[6] or "").lower()
            row_mode = str(row[5] or "").lower()
            row_result_dir = self._canonical_busco_path(row[7]) if row[7] else None
            if row_result_dir != target_dir:
                continue
            stored_dir = os.path.abspath(str(row[7])) if row[7] else ""
            score = (
                1 if stored_dir == target_dir else 0,
                1 if int(row[14] or 0) else 0,
                self.db_manager.busco.count_run_family_rows(row_id),
                str(row[13] or ""),
                row_id,
            )
            if row_pipeline == str(pipeline).lower() and row_mode == str(input_mode).lower():
                exact_matches.append((score, row_id))
            else:
                canonical_matches.append((score, row_id))
        if exact_matches:
            return sorted(exact_matches, reverse=True)[0][1], False
        if canonical_matches:
            return sorted(canonical_matches, reverse=True)[0][1], True
        return None

    def _discovered_run_metadata(self, accession: str, run_dirs: List[str]) -> Dict[str, Dict[str, str]]:
        metadata: Dict[str, Dict[str, str]] = {}
        for run_dir in run_dirs or []:
            canonical_run_dir = self._canonical_busco_path(run_dir)
            if not canonical_run_dir:
                continue
            canonical_result_dir = self._canonical_busco_path(os.path.dirname(canonical_run_dir))
            if not canonical_result_dir:
                continue
            summary_input = self._inspect_summary_input(accession, canonical_run_dir)
            inferred_pipeline = self._infer_pipeline_from_run_dir(canonical_run_dir)
            metadata[canonical_result_dir] = {
                "run_dir": canonical_run_dir,
                "result_dir": canonical_result_dir,
                "pipeline": str(summary_input.get("summary_pipeline") or inferred_pipeline),
                "input_mode": str(
                    summary_input.get("verified_input_mode")
                    or self._infer_input_mode_from_run_dir(canonical_run_dir, pipeline=inferred_pipeline)
                ),
            }
        return metadata

    def _tidy_duplicate_runs(
        self,
        accession: str,
        library_id: int,
        run_rows: List[tuple],
        discovered_by_result_dir: Dict[str, Dict[str, str]],
    ) -> List[Dict[str, Any]]:
        duplicate_rows: List[Dict[str, Any]] = []
        delete_run_ids: set[int] = set()
        rows_by_key: Dict[Tuple[Optional[str], str, str], List[tuple]] = defaultdict(list)
        rows_by_result_dir: Dict[Optional[str], List[tuple]] = defaultdict(list)

        for row in run_rows or []:
            canonical_result_dir = self._canonical_busco_path(row[7]) if row[7] else None
            pipeline = str(row[6] or "").lower()
            input_mode = str(row[5] or "").lower()
            rows_by_key[(canonical_result_dir, pipeline, input_mode)].append(row)
            rows_by_result_dir[canonical_result_dir].append(row)

        def _row_score(row: tuple) -> Tuple[int, int, int, str, int]:
            canonical_result_dir = self._canonical_busco_path(row[7]) if row[7] else None
            stored_dir = os.path.abspath(str(row[7])) if row[7] else ""
            return (
                1 if stored_dir == canonical_result_dir else 0,
                1 if int(row[14] or 0) else 0,
                self.db_manager.busco.count_run_family_rows(int(row[0])),
                str(row[13] or ""),
                int(row[0]),
            )

        for key, rows in rows_by_key.items():
            if len(rows) <= 1:
                continue
            keep = sorted(rows, key=_row_score, reverse=True)[0]
            for row in rows:
                run_id = int(row[0])
                if run_id == int(keep[0]):
                    continue
                delete_run_ids.add(run_id)
                duplicate_rows.append(
                    {
                        "accession": accession,
                        "run_id": run_id,
                        "library_name": str(row[3] or ""),
                        "change_type": "duplicate_run_removed",
                        "status_before": str(row[8] or ""),
                        "status_after": "deleted",
                        "artifact_updates": "",
                        "primary_changes": "",
                        "note": f"duplicate_of_run:{int(keep[0])}",
                    }
                )

        for canonical_result_dir, rows in rows_by_result_dir.items():
            discovered = discovered_by_result_dir.get(canonical_result_dir)
            if not discovered:
                continue
            actual_pipeline = str(discovered["pipeline"]).lower()
            actual_input_mode = str(discovered["input_mode"]).lower()
            has_correct_row = any(
                str(row[6] or "").lower() == actual_pipeline
                and str(row[5] or "").lower() == actual_input_mode
                and int(row[0]) not in delete_run_ids
                for row in rows
            )
            if not has_correct_row:
                continue
            for row in rows:
                run_id = int(row[0])
                if run_id in delete_run_ids:
                    continue
                row_pipeline = str(row[6] or "").lower()
                row_input_mode = str(row[5] or "").lower()
                if row_pipeline == actual_pipeline and row_input_mode == actual_input_mode:
                    continue
                delete_run_ids.add(run_id)
                duplicate_rows.append(
                    {
                        "accession": accession,
                        "run_id": run_id,
                        "library_name": str(row[3] or ""),
                        "change_type": "duplicate_run_removed",
                        "status_before": str(row[8] or ""),
                        "status_after": "deleted",
                        "artifact_updates": "",
                        "primary_changes": "",
                        "note": f"mismatched_identity:expected={actual_pipeline}/{actual_input_mode}",
                    }
                )

        if not self.repair:
            return duplicate_rows

        for run_id in sorted(delete_run_ids):
            self.db_manager.busco.delete_run(run_id)
            self.log(
                f"{accession}: removed duplicate BUSCO run {run_id} during verification cleanup.",
                "WARNING",
            )
        return duplicate_rows

    def _artifact_rows_for_run(self, run_id: int) -> List[tuple]:
        return self.db_manager.artifacts.find(owner_type="busco_run", owner_id=int(run_id))

    def _stale_run_artifacts(self, run_id: int, *, repair: bool = False) -> dict[str, Any]:
        artifact_rows = self._artifact_rows_for_run(run_id)
        staled_updates = 0
        for row in artifact_rows:
            if repair and str(row[5] or "ready") != "stale":
                self.db_manager.artifacts.set_status(int(row[0]), "stale")
                staled_updates += 1
        return {
            "checked": len(artifact_rows),
            "ready": 0,
            "stale": len(artifact_rows),
            "restored": 0,
            "staled_updates": staled_updates,
            "ready_size_updates": 0,
            "missing_rows": [],
        }

    def _run_dir_for_row(self, run_row: tuple) -> Optional[str]:
        run_id = int(run_row[0])
        result_dir = run_row[7]
        if result_dir and os.path.isdir(str(result_dir)):
            basename = os.path.basename(str(result_dir))
            if basename.startswith("run_"):
                return str(result_dir)
            run_dirs = sorted(glob.glob(os.path.join(str(result_dir), "run_*")))
            if run_dirs:
                return run_dirs[0]
        artifacts = self.db_manager.artifacts.find(owner_type="busco_run", owner_id=run_id, artifact_type="busco_run_dir")
        for artifact in artifacts:
            resolved = self.db_manager.artifacts.resolve_path(artifact)
            if resolved and os.path.isdir(resolved):
                return resolved
        return None

    def _verify_run_artifacts(
        self,
        run_row: tuple,
        accession: str,
        library_id: int,
        run_dirs_by_result_dir: dict[str, str],
        *,
        genome_path: Optional[str] = None,
    ) -> dict[str, Any]:
        run_id = int(run_row[0])
        run_dir = self._run_dir_for_row(run_row)
        if not run_dir and run_row[7]:
            run_dir = run_dirs_by_result_dir.get(str(self._canonical_busco_path(run_row[7])))
        if run_dir and genome_path and not self._path_within(run_dir, genome_path):
            if self.repair and self.db_manager.busco.get_run_status(run_id) != "stale":
                self.db_manager.busco.set_run_status(run_id, "stale")
                self.log(
                    f"{accession}: BUSCO run {run_id} is outside the current genome binding and was marked stale.",
                    "WARNING",
                )
            summary = self._stale_run_artifacts(run_id, repair=self.repair)
            summary["run_dir"] = run_dir
            summary["usable"] = False
            summary["outside_current_binding"] = True
            return summary
        if run_dir and self.verify_artifacts and self.repair:
            locations = self.db_manager.busco.get_run_family_locations(run_id)
            self._register_busco_run_artifacts(run_id, accession, library_id, run_dir, locations)
        artifact_rows = self._artifact_rows_for_run(run_id)
        summary = self._verify_artifact_rows(
            artifact_rows,
            repair=self.repair,
            stale_missing=self.stale_missing,
            restore_found=self.restore_found,
        )
        summary["run_dir"] = run_dir
        usable = self._run_is_usable_for_purpose(run_id, "default")
        target_status = "completed" if usable else "stale"
        if self.repair and self.db_manager.busco.get_run_status(run_id) != target_status:
            self.db_manager.busco.set_run_status(run_id, target_status)
        summary["usable"] = usable
        return summary

    def _run_has_sequence_kind(self, run_id: int, sequence_kind: str) -> bool:
        artifact_rows = self.db_manager.artifacts.find(
            owner_type="busco_run",
            owner_id=int(run_id),
            artifact_type="busco_family_sequence",
            sequence_kind=sequence_kind,
        )
        for row in artifact_rows:
            exists, _resolved, _size = self._artifact_path_exists(row)
            if exists:
                return True
        for _family_id, _library_id, _accession, location in self.db_manager.busco.get_run_family_locations(run_id):
            kind = self._sequence_kind_for_path(location)
            if kind == sequence_kind and location and os.path.exists(location):
                return True
        return False

    def _run_is_usable_for_purpose(self, run_id: int, purpose: str) -> bool:
        run_row = self.db_manager.busco.get_run(run_id)
        if not run_row:
            return False
        result_dir = run_row[6]
        has_core_dir = bool(result_dir and os.path.isdir(str(result_dir)))
        if not has_core_dir:
            for artifact in self.db_manager.artifacts.find(owner_type="busco_run", owner_id=int(run_id)):
                if str(artifact[3]) in {"busco_result_root", "busco_run_dir"}:
                    exists, _resolved, _size = self._artifact_path_exists(artifact)
                    if exists:
                        has_core_dir = True
                        break
        has_table = False
        if self.db_manager.busco.count_run_family_rows(run_id) > 0:
            has_table = True
        else:
            for artifact in self.db_manager.artifacts.find(
                owner_type="busco_run",
                owner_id=int(run_id),
                artifact_type="busco_full_table_tsv",
            ):
                exists, _resolved, _size = self._artifact_path_exists(artifact)
                if exists:
                    has_table = True
                    break
        if not (has_core_dir and has_table):
            return False
        if purpose == "export_protein":
            return self._run_has_sequence_kind(run_id, "prot")
        if purpose == "export_nucleotide":
            return self._run_has_sequence_kind(run_id, "nucl")
        return True

    def _repair_primary_assignments(self, accession: str, library_id: int) -> list[str]:
        changes: list[str] = []
        run_rows = self.db_manager.busco.get_runs_for_accessions([accession], library_id=library_id)
        default_profile = self.db_manager.proteomes.get_default_profile_name(accession)
        if not run_rows:
            legacy_summary_row = self.db_manager.cursor.execute(
                "SELECT 1 FROM BUSCO_Results WHERE accession = ? AND library_id = ? LIMIT 1",
                (accession, int(library_id)),
            ).fetchone()
            legacy_family_row = self.db_manager.cursor.execute(
                "SELECT 1 FROM BUSCO_Family_Data WHERE accession = ? AND library_id = ? LIMIT 1",
                (accession, int(library_id)),
            ).fetchone()
            has_orphaned_legacy_records = bool(legacy_summary_row or legacy_family_row)
            for purpose in ("default", "export_protein", "export_nucleotide"):
                if self.db_manager.busco.is_manual_primary_override(accession, library_id, purpose=purpose):
                    changes.append(f"preserved_manual:{purpose}")
                    continue
                current = self.db_manager.busco.get_primary_run(accession, library_id, purpose=purpose)
                if current and self.repair:
                    self.db_manager.busco.clear_primary_run(accession, library_id, purpose=purpose)
                if current:
                    changes.append(f"cleared:{purpose}")
            if has_orphaned_legacy_records:
                if self.repair:
                    self.db_manager.busco.delete_records(accession, library_id)
                    changes.append("cleared:legacy_records")
                else:
                    changes.append("orphaned:legacy_records")
            return changes

        for purpose in ("default", "export_protein", "export_nucleotide"):
            if self.db_manager.busco.is_manual_primary_override(accession, library_id, purpose=purpose):
                changes.append(f"preserved_manual:{purpose}")
                continue
            current = self.db_manager.busco.get_primary_run(accession, library_id, purpose=purpose)
            current_run_id = int(current[0]) if current and current[0] is not None else None
            current_ok = bool(current_run_id is not None and self._run_is_usable_for_purpose(current_run_id, purpose))
            if current_ok:
                continue
            best = None
            if default_profile:
                profile_runs = [
                    row
                    for row in self.db_manager.busco.get_runs_for_primary_choice(
                        accession,
                        library_id,
                        proteome_profile=default_profile,
                    )
                    if self._run_is_usable_for_purpose(int(row[0]), purpose)
                ]
                if profile_runs:
                    best = sorted(profile_runs, key=lambda row: (str(row[13] or ""), int(row[0])), reverse=True)[0]
            if best is None:
                best = self.db_manager.busco.choose_best_run(
                    accession,
                    library_id,
                    purpose=purpose,
                    preferred_proteome_profile=default_profile,
                )
            if best:
                best_run_id = int(best[0])
                if self.repair:
                    self.db_manager.busco.set_primary_run(
                        accession=accession,
                        library_id=library_id,
                        run_id=best_run_id,
                        purpose=purpose,
                        policy="verify_repaired",
                        updated_by="verify-busco",
                    )
                changes.append(f"reassigned:{purpose}:{best_run_id}")
            elif current_run_id is not None:
                if self.repair:
                    self.db_manager.busco.clear_primary_run(accession, library_id, purpose=purpose)
                    if purpose == "default":
                        self.db_manager.busco.delete_records(accession, library_id)
                changes.append(f"cleared:{purpose}")
        return changes

    def _parse_busco_table(self, accession: str, run_dir: str, library_id: int):
        """Parse BUSCO full_table and produce rows; returns (family_data, family_locations, summary_counts)."""
        table = self._busco_full_table(run_dir)
        if not table:
            self.log(
                f"{accession}: skipping BUSCO reingest for {run_dir} because full_table.tsv is missing.",
                "WARNING",
            )
            return None
        family_data = []
        family_locations = []
        counts = {"Single copy BUSCOs": 0, "Multi copy BUSCOs": 0, "Fragmented BUSCOs": 0, "Missing BUSCOs": 0}
        seen_family_status: Dict[str, int] = {}
        threshold = 1e5  # guard against mis-parsed coordinates
        try:
            with open(table, "r") as handle:
                reader = csv.reader(handle, delimiter="\t")
                for row in reader:
                    if not row or row[0].startswith("#"):
                        continue
                    if len(row) < 2:
                        continue
                    family_id = row[0]
                    status_str = row[1]
                    status_map = {"Complete": 1, "Duplicated": 2, "Fragmented": 3, "Missing": 4}
                    status = status_map.get(status_str, 0)
                    # Only count each family once toward summary counts to avoid double-counting duplicated copies
                    if family_id not in seen_family_status:
                        counts_key = {
                            1: "Single copy BUSCOs",
                            2: "Multi copy BUSCOs",
                            3: "Fragmented BUSCOs",
                            4: "Missing BUSCOs",
                        }.get(status)
                        if counts_key:
                            counts[counts_key] += 1
                        seen_family_status[family_id] = status
                    sequence = row[2] if len(row) >= 3 and row[2] else None
                    score = None
                    length = None
                    if len(row) >= 7 and row[6]:
                        try:
                            score = float(row[6])
                        except (TypeError, ValueError):
                            score = None
                    if len(row) >= 8 and row[7]:
                        try:
                            length = int(row[7])
                        except (TypeError, ValueError):
                            length = None
                    if score is None and len(row) >= 4 and row[3]:
                        try:
                            score = float(row[3])
                        except (TypeError, ValueError):
                            score = None
                    if length is None and len(row) >= 5 and row[4]:
                        try:
                            length = int(row[4])
                        except (TypeError, ValueError):
                            length = None
                    if score is not None and score > threshold:
                        self.log(f"Discarding anomalous BUSCO score {score} for {family_id} ({accession}); likely a coordinate.", "WARNING")
                        score = None
                    location = self._busco_seq_path(run_dir, family_id, status)
                    family_data.append((family_id, library_id, accession, status, sequence, score, length))
                    family_locations.append((family_id, library_id, accession, location))
        except (OSError, UnicodeError, csv.Error, ValueError) as exc:
            raise RuntimeError(f"Failed to parse BUSCO full_table for {accession} in {run_dir}: {exc}") from exc
        return family_data, family_locations, counts

    def _process_busco_accession(self, acc: str, genome_root: str) -> Dict[str, Any]:
        genome_path = self.db_manager.genomes.resolve_path(acc) or os.path.join(genome_root, acc)
        profile_sync: Dict[str, Any] = {}
        if self.repair and os.path.isdir(genome_path):
            profile_sync = _sync_profile_state(self, acc, genome_path)
        scoped_targets = self._scoped_busco_targets(acc, genome_path)
        result = {
            "reingested": 0,
            "queued": 0,
            "missing": 0,
            "skipped": 0,
            "anomalies_fixed": 0,
            "stale_artifacts": 0,
            "restored_artifacts": 0,
            "primary_repairs": 0,
            "report_rows": [],
            "summary_rows": [],
        }
        removed_staged_inputs = [str(path) for path in (profile_sync.get("removed_staged_inputs") or [])]
        if removed_staged_inputs:
            result["summary_rows"].append(
                {
                    "accession": acc,
                    "run_id": "",
                    "library_name": "",
                    "change_type": "staged_input_removed",
                    "status_before": "",
                    "status_after": "",
                    "artifact_updates": "",
                    "primary_changes": "",
                    "note": f"removed_staged_busco_inputs:{len(removed_staged_inputs)}",
                }
            )
        if not scoped_targets:
            result["missing"] += 1
            result["report_rows"].append({"accession": acc, "status": "missing", "lineage": "", "note": "no_known_busco_scope"})
            return result

        for busco_lib_id, busco_lineage_name, run_dirs in scoped_targets:
            validation_notes: List[str] = []
            if not run_dirs:
                result["missing"] += 1
                self.log(f"{acc}: BUSCO run folder not found under {genome_path} for {busco_lineage_name}.", "WARNING")
                if self.queue_missing:
                    self.queue_subtask(
                        job_type=4,
                        status="P",
                        priority=1,
                        data={
                            "accession": acc,
                            "lineage": busco_lineage_name,
                            "format": "protein",
                            "force": False,
                        },
                    )
                    result["queued"] += 1
                    result["summary_rows"].append(
                        {
                            "accession": acc,
                            "run_id": "",
                            "library_name": busco_lineage_name,
                            "change_type": "queued_missing_run",
                            "status_before": "",
                            "status_after": "",
                            "artifact_updates": "",
                            "primary_changes": "",
                            "note": "queued_missing_busco_run",
                        }
                    )

            existing_rows = self.db_manager.busco.get_family_results_for_library(
                library_id=busco_lib_id,
                accessions=[acc],
            ) or []
            has_rows = bool(existing_rows)
            has_anomalies = any(r[5] is not None and r[5] > 1e5 for r in existing_rows)
            existing_run_rows = self.db_manager.busco.get_runs_for_accessions([acc], library_id=busco_lib_id)
            if self.run_id is not None:
                existing_run_rows = [row for row in existing_run_rows if int(row[0]) == int(self.run_id)]
            existing_run_rows_by_id = {int(row[0]): row for row in existing_run_rows}
            run_dirs_by_result_dir = {
                str(self._canonical_busco_path(os.path.dirname(run_dir))): run_dir
                for run_dir in run_dirs
            }
            discovered_by_result_dir = self._discovered_run_metadata(acc, run_dirs)
            would_reingest_reasons: set[str] = set()
            if not has_rows:
                would_reingest_reasons.add("missing_legacy_family_rows")
            if has_anomalies:
                would_reingest_reasons.add("anomalous_legacy_family_rows")
            if len(existing_run_rows) > len({(self._canonical_busco_path(row[7]) if row[7] else None, str(row[6] or "").lower(), str(row[5] or "").lower()) for row in existing_run_rows}):
                would_reingest_reasons.add("duplicate_run_candidates")

            run_prechecks: List[Dict[str, Any]] = []
            for run_dir in run_dirs or []:
                inferred_pipeline = self._infer_pipeline_from_run_dir(run_dir)
                summary_input = self._inspect_summary_input(acc, run_dir)
                summary_notes = list(summary_input.get("notes") or [])
                validation_notes.extend(summary_notes)
                summary_pipeline = str(summary_input.get("summary_pipeline") or "") or None
                pipeline = summary_pipeline or inferred_pipeline
                inferred_input_mode = self._infer_input_mode_from_run_dir(run_dir, pipeline=pipeline)
                input_mode = str(summary_input.get("verified_input_mode") or inferred_input_mode)
                if summary_input.get("summary_input_mode") and summary_input.get("summary_input_mode") != inferred_input_mode:
                    validation_notes.append(f"summary_verified_input_mode:{summary_input.get('summary_input_mode')}")
                    would_reingest_reasons.add("summary_input_mode_db_mismatch")
                if summary_pipeline and summary_pipeline != inferred_pipeline:
                    validation_notes.append(f"summary_verified_pipeline:{summary_pipeline}")
                    would_reingest_reasons.add("summary_pipeline_dir_mismatch")
                result_dir = self._canonical_busco_path(os.path.dirname(run_dir))
                match = self._match_existing_run(
                    existing_run_rows,
                    result_dir=result_dir,
                    pipeline=pipeline,
                    input_mode=input_mode,
                )
                run_id = None
                identity_mismatch = False
                if match is not None:
                    run_id, identity_mismatch = match
                matched_run_row = existing_run_rows_by_id.get(int(run_id)) if run_id is not None else None
                if run_id is None:
                    would_reingest_reasons.add("missing_run_row")
                if identity_mismatch:
                    would_reingest_reasons.add("run_identity_mismatch")
                profile_row, profile_notes = self._resolve_proteome_profile_row_for_run(
                    acc,
                    input_mode=input_mode,
                    summary_input_path=summary_input.get("summary_input_path"),
                    matched_run_row=matched_run_row,
                )
                validation_notes.extend(profile_notes)
                if any(
                    note.startswith("summary_input_artifact_mismatch:")
                    or note.startswith("summary_expected_artifact_missing:")
                    or note.startswith("summary_mode_path_mismatch:")
                    for note in summary_notes
                ):
                    would_reingest_reasons.add("summary_input_artifact_mismatch")
                if any(
                    note.startswith("summary_input_profile_mismatch:")
                    or note == "summary_input_profile_unresolved"
                    for note in profile_notes
                ):
                    would_reingest_reasons.add("summary_input_profile_mismatch")
                if summary_pipeline and matched_run_row and summary_pipeline != str(matched_run_row[6] or "").lower():
                    would_reingest_reasons.add("summary_pipeline_db_mismatch")
                if summary_input.get("summary_input_mode") and matched_run_row and str(summary_input.get("summary_input_mode")) != str(matched_run_row[5] or "").lower():
                    would_reingest_reasons.add("summary_input_mode_db_mismatch")
                if matched_run_row and profile_row and str(matched_run_row[20] or "") != str(profile_row[2] or ""):
                    would_reingest_reasons.add("run_profile_db_mismatch")
                matched_counts = self._run_counts_from_row(matched_run_row)
                summary_counts = summary_input.get("summary_counts")
                if summary_counts and not self._counts_equal(summary_counts, matched_counts):
                    would_reingest_reasons.add("summary_db_counts_mismatch")
                    validation_notes.append("summary_db_counts_mismatch")
                if matched_run_row and self.db_manager.busco.count_run_family_rows(int(matched_run_row[0])) == 0:
                    would_reingest_reasons.add("canonical_run_family_rows_missing")
                run_prechecks.append(
                    {
                        "run_dir": run_dir,
                        "result_dir": result_dir,
                        "pipeline": pipeline,
                        "inferred_pipeline": inferred_pipeline,
                        "input_mode": input_mode,
                        "inferred_input_mode": inferred_input_mode,
                        "summary_input": summary_input,
                        "matched_run_id": int(run_id) if run_id is not None else None,
                        "identity_mismatch": identity_mismatch,
                        "matched_run_row": matched_run_row,
                        "proteome_profile_row": profile_row,
                    }
                )

            if self.reingest_all:
                reingest_scope = "all"
                would_reingest = bool(run_dirs)
            else:
                reingest_scope = "conditional"
                would_reingest = bool(self.reingest and would_reingest_reasons)
            need_reingest = bool(self.repair and would_reingest)

            if not need_reingest and not self.verify_artifacts:
                result["skipped"] += 1
                result["report_rows"].append({"accession": acc, "status": "ok", "lineage": busco_lineage_name, "note": "kept"})
                continue

            if need_reingest and run_dirs:
                full_table_reason_notes: List[str] = []
                for item in run_prechecks:
                    run_dir = str(item["run_dir"])
                    summary_counts = item["summary_input"].get("summary_counts")
                    needs_full_table_check = bool(self.reingest_all or summary_counts or has_anomalies or item["matched_run_id"] is None)
                    if not needs_full_table_check:
                        continue
                    full_counts, full_notes = self._parse_busco_table_counts(acc, run_dir, busco_lib_id)
                    if full_notes:
                        full_table_reason_notes.extend(full_notes)
                    if full_counts is not None:
                        matched_counts = self._run_counts_from_row(item["matched_run_row"])
                        if summary_counts and not self._counts_equal(summary_counts, full_counts):
                            would_reingest_reasons.add("full_table_summary_mismatch")
                        if matched_counts and not self._counts_equal(matched_counts, full_counts):
                            would_reingest_reasons.add("db_full_table_mismatch")
                validation_notes.extend(full_table_reason_notes)

            run_summaries: List[Tuple[int, Dict[str, int], str, str]] = []
            reingest_items = [item for item in run_prechecks if self._busco_full_table(str(item["run_dir"]))]
            skipped_reingest_runs = max(len(run_prechecks) - len(reingest_items), 0)
            if need_reingest and not reingest_items:
                validation_notes.append("reingest_skipped_no_full_table")
                need_reingest = False
            if need_reingest and not self.db_manager.busco.delete_records(acc, busco_lib_id):
                self.log(f"{acc}: failed to delete existing BUSCO records for {busco_lineage_name}; continuing with re-ingest.", "WARNING")
            if need_reingest:
                for item in reingest_items:
                    run_dir = str(item["run_dir"])
                    parsed = self._parse_busco_table(acc, run_dir, busco_lib_id)
                    if parsed is None:
                        continue
                    family_data, family_locations, counts = parsed
                    run_id = item["matched_run_id"]
                    profile_row = item.get("proteome_profile_row")
                    profile_id = int(profile_row[0]) if profile_row and profile_row[0] is not None else None
                    if run_id is None:
                        run_id = self.db_manager.busco.create_run(
                            accession=acc,
                            library_id=busco_lib_id,
                            lineage_name=busco_lineage_name,
                            input_mode=item["input_mode"],
                            pipeline=item["pipeline"],
                            pipeline_params_effective={},
                            pipeline_params_source={"migration": "verify-busco"},
                            busco_cli_args=[],
                            result_dir=item["result_dir"],
                            proteome_profile_id=profile_id,
                            status="completed",
                        )
                    if run_id is not None:
                        self.db_manager.busco.update_run(
                            run_id,
                            status="completed",
                            result_dir=item["result_dir"],
                            lineage_name=busco_lineage_name,
                            input_mode=item["input_mode"],
                            pipeline=item["pipeline"],
                            proteome_profile_id=profile_id,
                            counts=counts,
                            completed=True,
                        )
                        self.db_manager.busco.add_run_family_data(run_id, family_data)
                        self.db_manager.busco.add_run_family_locations(run_id, list(set(family_locations)))
                        self._register_busco_run_artifacts(int(run_id), acc, busco_lib_id, run_dir, list(set(family_locations)))
                        run_summaries.append((int(run_id), counts, item["pipeline"], item["input_mode"]))
                        if item["identity_mismatch"]:
                            result["summary_rows"].append(
                                {
                                    "accession": acc,
                                    "run_id": int(run_id),
                                    "library_name": busco_lineage_name,
                                    "change_type": "run_identity_repaired",
                                    "status_before": "",
                                    "status_after": "",
                                    "artifact_updates": "",
                                    "primary_changes": "",
                                    "note": f"updated_identity:{item['pipeline']}/{item['input_mode']}"
                                    + (
                                        f";summary={os.path.basename(str(item['summary_input'].get('summary_path') or ''))}"
                                        if item["summary_input"].get("summary_path")
                                        else ""
                                    ),
                                }
                            )
                if run_summaries:
                    best = sorted(
                        run_summaries,
                        key=lambda r: (
                            int((r[1] or {}).get("Single copy BUSCOs", 0) or 0),
                            int((r[1] or {}).get("Single copy BUSCOs", 0) or 0) + int((r[1] or {}).get("Multi copy BUSCOs", 0) or 0),
                            -int((r[1] or {}).get("Multi copy BUSCOs", 0) or 0),
                            r[0],
                        ),
                        reverse=True,
                    )[0]
                    best_run_id, best_counts, best_pipeline, _best_mode = best
                    self.db_manager.busco.refresh_auto_primary_runs_for_accession(
                        acc,
                        busco_lib_id,
                        updated_by="verify-busco",
                        policy="auto_best",
                    )
                    self.db_manager.busco.add_results(acc, busco_lib_id, best_counts, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                    result["summary_rows"].append(
                        {
                            "accession": acc,
                            "run_id": best_run_id,
                            "library_name": busco_lineage_name,
                            "change_type": "primary_selected",
                            "status_before": "",
                            "status_after": "",
                            "artifact_updates": "",
                            "primary_changes": "set:default;set:export_protein" + (";set:export_nucleotide" if best_pipeline in {"augustus", "metaeuk"} else ""),
                            "note": "best_run_selected_during_verify",
                        }
                    )
            if self.repair and not need_reingest:
                for item in run_prechecks:
                    run_id = item["matched_run_id"]
                    if run_id is None:
                        continue
                    matched = item["matched_run_row"]
                    before_pipeline = str(matched[6] or "") if matched else ""
                    before_mode = str(matched[5] or "") if matched else ""
                    before_result_dir = str(matched[7] or "") if matched else ""
                    before_profile = str(matched[20] or "") if matched and len(matched) > 20 else ""
                    after_pipeline = str(item["pipeline"])
                    after_mode = str(item["input_mode"])
                    after_result_dir = str(item["result_dir"])
                    profile_row = item.get("proteome_profile_row")
                    after_profile_id = int(profile_row[0]) if profile_row and profile_row[0] is not None else None
                    after_profile = str(profile_row[2] or "") if profile_row and profile_row[2] is not None else ""
                    if (
                        before_pipeline != after_pipeline
                        or before_mode != after_mode
                        or before_profile != after_profile
                        or self._canonical_busco_path(before_result_dir) != self._canonical_busco_path(after_result_dir)
                    ):
                        self.db_manager.busco.update_run(
                            int(run_id),
                            result_dir=after_result_dir,
                            lineage_name=busco_lineage_name,
                            input_mode=after_mode,
                            pipeline=after_pipeline,
                            proteome_profile_id=after_profile_id,
                        )
                        result["summary_rows"].append(
                            {
                                "accession": acc,
                                "run_id": int(run_id),
                                "library_name": busco_lineage_name,
                                "change_type": "run_identity_repaired",
                                "status_before": "",
                                "status_after": "",
                                "artifact_updates": "",
                                "primary_changes": "",
                                "note": f"updated_identity:{after_pipeline}/{after_mode}/{after_profile or '-'}",
                            }
                        )

            current_run_rows = self.db_manager.busco.get_runs_for_accessions([acc], library_id=busco_lib_id)
            if self.run_id is not None:
                current_run_rows = [row for row in current_run_rows if int(row[0]) == int(self.run_id)]
            if self.repair:
                cleanup_rows = self._tidy_duplicate_runs(acc, busco_lib_id, current_run_rows, discovered_by_result_dir)
                if cleanup_rows:
                    result["summary_rows"].extend(cleanup_rows)
                    current_run_rows = self.db_manager.busco.get_runs_for_accessions([acc], library_id=busco_lib_id)
                    if self.run_id is not None:
                        current_run_rows = [row for row in current_run_rows if int(row[0]) == int(self.run_id)]
            for run_row in current_run_rows:
                run_id = int(run_row[0])
                run_status_before = self.db_manager.busco.get_run_status(run_id)
                artifact_summary = self._verify_run_artifacts(
                    run_row,
                    acc,
                    busco_lib_id,
                    run_dirs_by_result_dir,
                    genome_path=genome_path,
                )
                result["stale_artifacts"] += int(artifact_summary.get("stale", 0) or 0)
                result["restored_artifacts"] += int(artifact_summary.get("restored", 0) or 0)
                run_status_after = self.db_manager.busco.get_run_status(run_id)
                artifact_updates = []
                if int(artifact_summary.get("staled_updates", 0) or 0):
                    artifact_updates.append(f"staled:{int(artifact_summary.get('staled_updates', 0) or 0)}")
                if int(artifact_summary.get("restored", 0) or 0):
                    artifact_updates.append(f"restored:{int(artifact_summary.get('restored', 0) or 0)}")
                if run_status_before != run_status_after or artifact_updates:
                    result["summary_rows"].append(
                        {
                            "accession": acc,
                            "run_id": run_id,
                            "library_name": busco_lineage_name,
                            "change_type": "run_reconciled",
                            "status_before": str(run_status_before or ""),
                            "status_after": str(run_status_after or ""),
                            "artifact_updates": ";".join(artifact_updates),
                            "primary_changes": "",
                            "note": "run_status_or_artifacts_updated",
                        }
                    )

            primary_changes: List[str] = []
            if self.reassign_primary and (self.repair or self.verify_artifacts):
                primary_changes = self._repair_primary_assignments(acc, busco_lib_id)
                if primary_changes:
                    result["primary_repairs"] += len(primary_changes)
                    self.log(f"{acc}: BUSCO primary updates for {busco_lineage_name} -> {','.join(primary_changes)}", "WARNING")
                    result["summary_rows"].append(
                        {
                            "accession": acc,
                            "run_id": "",
                            "library_name": busco_lineage_name,
                            "change_type": "primary_repaired",
                            "status_before": "",
                            "status_after": "",
                            "artifact_updates": "",
                            "primary_changes": ";".join(primary_changes),
                            "note": "primary_selection_updated",
                        }
                    )

            if has_anomalies and need_reingest:
                result["anomalies_fixed"] += 1
            if need_reingest and run_dirs:
                result["reingested"] += 1
                self.log(
                    f"{acc}: reingested BUSCO results from {len(run_summaries)} run folder(s) for {busco_lineage_name}."
                    + (
                        f" Skipped {skipped_reingest_runs} run folder(s) without full_table.tsv."
                        if skipped_reingest_runs
                        else ""
                    ),
                    "INFO",
                )
                reasons = ";".join(sorted(dict.fromkeys(would_reingest_reasons)))
                note = f"reingest_scope={reingest_scope};reingest_reason={reasons or 'requested'};anomalies={has_anomalies}"
                if skipped_reingest_runs:
                    note += f";reingest_skipped_missing_full_table={skipped_reingest_runs}"
                if validation_notes:
                    note += f";{';'.join(dict.fromkeys(validation_notes))}"
                if primary_changes:
                    note += f";primary={','.join(primary_changes)}"
                result["report_rows"].append({"accession": acc, "status": "reingested", "lineage": busco_lineage_name, "note": note})
                result["summary_rows"].append(
                    {
                        "accession": acc,
                        "run_id": "",
                        "library_name": busco_lineage_name,
                        "change_type": "reingested",
                        "status_before": "",
                        "status_after": "",
                        "artifact_updates": "",
                        "primary_changes": ";".join(primary_changes),
                        "note": note,
                    }
                )
            elif would_reingest and run_dirs:
                result["skipped"] += 1
                reasons = ";".join(sorted(dict.fromkeys(would_reingest_reasons)))
                note = f"would_reingest=yes;reingest_scope={reingest_scope};reingest_reason={reasons or 'requested'};anomalies={has_anomalies}"
                if validation_notes:
                    note += f";{';'.join(dict.fromkeys(validation_notes))}"
                result["report_rows"].append({"accession": acc, "status": "ok", "lineage": busco_lineage_name, "note": note})
            elif run_dirs:
                result["skipped"] += 1
                note = "would_reingest=no;reingest_scope=none;reingest_reason="
                if validation_notes:
                    note += f";{';'.join(dict.fromkeys(validation_notes))}"
                if primary_changes:
                    note += f";primary={','.join(primary_changes)}"
                result["report_rows"].append({"accession": acc, "status": "ok", "lineage": busco_lineage_name, "note": note})
            else:
                note = "no_run_folder"
                if primary_changes:
                    note += f";primary={','.join(primary_changes)}"
                result["report_rows"].append({"accession": acc, "status": "missing", "lineage": busco_lineage_name, "note": note})
        return result

    def run(self):
        if self.run_id:
            run_row = self.db_manager.busco.get_run(int(self.run_id))
            if not run_row:
                return self.handle_exception("BUSCO run_id not found for verification.", {"run_id": self.run_id})
            self.library_id = int(run_row[2])
            self.data["accessions"] = [str(run_row[1])]
        if not self.library_id and self.library_name:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
        if not self.library_id and not (self.all or self.data.get("accessions") or self.data.get("taxid") is not None or self.data.get("clade")):
            return self.handle_exception("library_id (or resolvable library_name) is required for BUSCO verification.", {"library_id": self.library_id, "library_name": self.library_name})

        try:
            selected = self.prepare_selectors(
                taxid=self.data.get("taxid"),
                additional=self.data.get("accessions"),
                root=self.data.get("root"),
                allow_all=self.all or bool(self.data.get("root")),
                downloaded_only=self.downloaded_only,
                released_after=self.data.get("after"),
                released_before=self.data.get("before"),
                level=self.data.get("level"),
                primary_only=self.data.get("primary_only"),
                require_candidates=not bool(self.discover),
            )
        except ValueError as exc:
            return self.handle_exception(str(exc), {})

        try:
            genome_root = self.db_manager.storage.require_root_base("genomes")
        except ValueError as exc:
            return self.handle_exception(str(exc), {})

        reingested = 0
        queued = 0
        missing = 0
        skipped = 0
        anomalies_fixed = 0
        stale_artifacts = 0
        restored_artifacts = 0
        primary_repairs = 0
        report_rows: List[Dict[str, Any]] = []
        summary_rows: List[Dict[str, Any]] = []

        worker_results, worker_errors = _run_verify_workers(
            self,
            selected,
            lambda acc: self._run_busco_worker(acc, genome_root),
            "BUSCO accession",
        )
        for worker_result in worker_results:
            if not worker_result:
                continue
            reingested += int(worker_result.get("reingested", 0) or 0)
            queued += int(worker_result.get("queued", 0) or 0)
            missing += int(worker_result.get("missing", 0) or 0)
            skipped += int(worker_result.get("skipped", 0) or 0)
            anomalies_fixed += int(worker_result.get("anomalies_fixed", 0) or 0)
            stale_artifacts += int(worker_result.get("stale_artifacts", 0) or 0)
            restored_artifacts += int(worker_result.get("restored_artifacts", 0) or 0)
            primary_repairs += int(worker_result.get("primary_repairs", 0) or 0)
            report_rows.extend(worker_result.get("report_rows", []))
            summary_rows.extend(worker_result.get("summary_rows", []))

        self.log(
            f"Verify BUSCO summary: reingested={reingested} skipped={skipped} queued={queued} "
            f"missing={missing} anomalies_fixed={anomalies_fixed} stale_artifacts={stale_artifacts} "
            f"restored_artifacts={restored_artifacts} primary_repairs={primary_repairs}",
            "INFO",
        )
        try:
            _write_tsv_report(self, "verify_busco.tsv", ["accession", "status", "lineage", "note"], report_rows)
            _write_tsv_report(
                self,
                "verify_busco_summary.tsv",
                ["accession", "run_id", "library_name", "change_type", "status_before", "status_after", "artifact_updates", "primary_changes", "note"],
                summary_rows,
            )
        except (OSError, UnicodeError) as exc:
            self.log(f"Failed to write BUSCO verify report: {exc}", "ERROR")
        if worker_errors:
            return self.handle_exception(_parallel_error_message("verify-busco worker failures", worker_errors))
        return True


class VerifyLibrariesTask(_ArtifactVerificationMixin, Task):
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads=required_threads)
        self.stage = checkpoint if checkpoint is not None else 0
        self.repair = bool(self.data.get("repair", False))
        self.all = bool(self.data.get("all", False))
        self.root = self.data.get("root")
        self.library_id = self.data.get("library_id")
        self.library_name = self.data.get("library_name")
        self.ref_accessions = normalize_accessions(self.data.get("ref_accessions") or [])

    def _target_library_ids(self) -> List[int]:
        explicit_ids = [int(item) for item in (self.data.get("_verify_target_library_ids") or [])]
        if explicit_ids:
            return explicit_ids
        root_id = None
        if self.root:
            row = self.db_manager.storage.resolve_root_token(self.root, kind="libraries")
            root_id = int(row[0])

        def _filter_by_root(ids: List[int]) -> List[int]:
            if root_id is None:
                return ids
            if not ids:
                return []
            placeholders = ",".join("?" for _ in ids)
            rows = self.db_manager.cursor.execute(
                f"SELECT library_id FROM Libraries WHERE storage_root_id = ? AND library_id IN ({placeholders})",
                (int(root_id), *[int(lib) for lib in ids]),
            ).fetchall()
            return [int(row[0]) for row in rows if row and row[0] is not None]

        if self.all:
            rows = self.db_manager.libraries.get(include_inactive=True) or []
            ids = [int(row[0]) for row in rows if row and row[0] is not None]
            return _filter_by_root(ids)
        if self.library_id:
            return _filter_by_root([int(self.library_id)])
        if self.library_name:
            lib_id = self.db_manager.libraries.get_id(self.library_name, include_inactive=True)
            return _filter_by_root([int(lib_id)]) if lib_id else []
        rows = self.db_manager.libraries.get_by_reference_accessions(self.ref_accessions)
        ids = [int(row[0]) for row in rows if row and row[0] is not None]
        return _filter_by_root(ids)

    def _write_core_set_json(self, library_id: int, path: str) -> bool:
        descriptions = self.db_manager.libraries.get_busco_descriptions(int(library_id)) or []
        family_ids = sorted({str(row[0]) for row in descriptions if row and row[0]})
        if not family_ids:
            return False
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(family_ids, indent=2), encoding="utf-8")
        return True

    def _run_library_worker(self, library_id: int) -> Dict[str, Any]:
        worker = VerifyLibrariesTask(
            self.db_manager.get_path(),
            self.task_id,
            0,
            data=json.dumps(self.data),
            required_threads=1,
        )
        try:
            return worker._process_library_id(library_id)
        finally:
            worker.db_manager.close()

    def _process_library_id(self, library_id: int) -> Dict[str, Any]:
        lib_rows = self.db_manager.libraries.get(library_id, include_inactive=True)
        if not lib_rows:
            return {"row": None, "summary_row": None}
        lib = lib_rows[0]
        name = str(lib[1])
        root = self.db_manager.libraries.resolve_path(library_id)
        parent_id = self.db_manager.libraries.get_parent_id(library_id)
        root_exists = bool(root and os.path.isdir(root))
        dataset_cfg = os.path.join(root, "dataset.cfg") if root else None
        core_set_json = os.path.join(root, "cleaned_busco_families.json") if root else None
        ref_accessions = self.db_manager.libraries.get_reference_assemblies(library_id) or []
        is_lineage = bool((root and "lineages" in Path(root).parts) or (dataset_cfg and os.path.exists(dataset_cfg)))
        status_before = self.db_manager.libraries.get_status(library_id)
        parent_ready = True
        parent_note = ""
        if not is_lineage and parent_id:
            parent_ready = self.db_manager.libraries.get_status(parent_id) == "ready"
            if not parent_ready:
                parent_note = f"parent_not_ready:{parent_id}"
        elif not is_lineage and not parent_id:
            parent_ready = False
            parent_note = "missing_parent"
        missing_refs = [acc for acc in ref_accessions if not self.db_manager.genomes.get(acc)]
        root_artifact_before = bool(self.db_manager.artifacts.find(owner_type="library", owner_id=library_id, artifact_type="library_root"))
        required_type = "library_dataset_cfg" if is_lineage else "library_core_set_json"
        required_artifact_before = bool(self.db_manager.artifacts.find(owner_type="library", owner_id=library_id, artifact_type=required_type))
        root_artifact_staled = 0
        required_artifact_staled = 0
        core_set_rewritten = False
        if root_exists:
            self.db_manager.artifacts.register(
                owner_type="library",
                owner_id=library_id,
                artifact_type="library_root",
                path=root,
                is_dir=True,
                format="directory",
                metadata={"library_id": library_id, "library_name": name},
            )
        else:
            for artifact in self.db_manager.artifacts.find(owner_type="library", owner_id=library_id, artifact_type="library_root"):
                if self.repair:
                    self.db_manager.artifacts.set_status(int(artifact[0]), "stale")
                    root_artifact_staled += 1
        required_path = dataset_cfg if is_lineage else core_set_json
        if (not is_lineage) and root_exists and self.repair and core_set_json and not os.path.exists(core_set_json):
            if self._write_core_set_json(library_id, core_set_json):
                self.log(f"{name}: rewrote cleaned_busco_families.json from BUSCO_descriptions.", "WARNING")
                core_set_rewritten = True
        if required_path and os.path.exists(required_path):
            self.db_manager.artifacts.register(
                owner_type="library",
                owner_id=library_id,
                artifact_type=required_type,
                path=required_path,
                format="json" if required_type == "library_core_set_json" else "cfg",
                metadata={"library_id": library_id, "library_name": name},
            )
        else:
            for artifact in self.db_manager.artifacts.find(owner_type="library", owner_id=library_id, artifact_type=required_type):
                if self.repair:
                    self.db_manager.artifacts.set_status(int(artifact[0]), "stale")
                    required_artifact_staled += 1
        usable = bool(root_exists and required_path and os.path.exists(required_path) and parent_ready and not missing_refs)
        target_status = "ready" if usable else "stale"
        if self.repair:
            self.db_manager.libraries.set_status(library_id, target_status)
        status_after = self.db_manager.libraries.get_status(library_id)
        root_artifact_after = bool(self.db_manager.artifacts.find(owner_type="library", owner_id=library_id, artifact_type="library_root"))
        required_artifact_after = bool(self.db_manager.artifacts.find(owner_type="library", owner_id=library_id, artifact_type=required_type))
        row = {
            "library_id": library_id,
            "library_name": name,
            "status": target_status,
            "root_exists": "yes" if root_exists else "no",
            "required_artifact": required_type,
            "required_present": "yes" if required_path and os.path.exists(required_path) else "no",
            "parent_ok": "yes" if parent_ready else "no",
            "missing_ref_accessions": ",".join(missing_refs),
            "note": parent_note,
        }
        summary_row = None
        if (
            status_before != status_after
            or core_set_rewritten
            or root_artifact_staled
            or required_artifact_staled
            or (not root_artifact_before and root_artifact_after)
            or (not required_artifact_before and required_artifact_after)
        ):
            summary_row = {
                "library_id": library_id,
                "library_name": name,
                "status_before": str(status_before or ""),
                "status_after": str(status_after or ""),
                "root_artifact_backfilled": "yes" if (not root_artifact_before and root_artifact_after) else "no",
                "required_artifact_backfilled": "yes" if (not required_artifact_before and required_artifact_after) else "no",
                "root_artifact_staled": str(root_artifact_staled),
                "required_artifact_staled": str(required_artifact_staled),
                "core_set_rewritten": "yes" if core_set_rewritten else "no",
                "note": parent_note or ("missing_refs:" + ",".join(missing_refs) if missing_refs else ""),
            }
        if not usable:
            self.log(f"{name}: marked stale for library verification.", "WARNING")
        return {"row": row, "summary_row": summary_row}

    def run(self):
        library_ids = self._target_library_ids()
        rows = []
        summary_rows = []
        worker_results, worker_errors = _run_verify_workers(
            self,
            library_ids,
            self._run_library_worker,
            "library",
        )
        for worker_result in worker_results:
            if not worker_result or not worker_result.get("row"):
                continue
            rows.append(worker_result["row"])
            if worker_result.get("summary_row"):
                summary_rows.append(worker_result["summary_row"])
        _write_tsv_report(
            self,
            "verify_libraries.tsv",
            ["library_id", "library_name", "status", "root_exists", "required_artifact", "required_present", "parent_ok", "missing_ref_accessions", "note"],
            rows,
        )
        _write_tsv_report(
            self,
            "verify_libraries_summary.tsv",
            ["library_id", "library_name", "status_before", "status_after", "root_artifact_backfilled", "required_artifact_backfilled", "root_artifact_staled", "required_artifact_staled", "core_set_rewritten", "note"],
            summary_rows,
        )
        if worker_errors:
            return self.handle_exception(_parallel_error_message("verify-libraries worker failures", worker_errors))
        return True


class VerifyOrthofinderTask(_ArtifactVerificationMixin, Task):
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads=required_threads)
        self.stage = checkpoint if checkpoint is not None else 0
        self.repair = bool(self.data.get("repair", False))
        self.all = bool(self.data.get("all", False))
        self.root = self.data.get("root")
        self.library_id = self.data.get("library_id")
        self.library_name = self.data.get("library_name")
        self.ref_accessions = normalize_accessions(self.data.get("ref_accessions") or [])

    def _target_runs(self) -> List[tuple]:
        explicit_run_ids = [int(item) for item in (self.data.get("_verify_target_run_ids") or [])]
        if explicit_run_ids:
            all_runs = self.db_manager.orthofinder.get_many(include_inactive=True) or []
            by_id = {int(row[0]): row for row in all_runs if row and row[0] is not None}
            return [by_id[run_id] for run_id in explicit_run_ids if run_id in by_id]
        root_id = None
        if self.root:
            row = self.db_manager.storage.resolve_root_token(self.root, kind="orthofinder")
            root_id = int(row[0])

        def _filter_by_root(runs: List[tuple]) -> List[tuple]:
            if root_id is None:
                return runs
            filtered: List[tuple] = []
            for run in runs:
                detected_root_id, _ = self.db_manager.storage.detect_root_for_path(run[2], kind="orthofinder")
                if detected_root_id is not None and int(detected_root_id) == int(root_id):
                    filtered.append(run)
            return filtered

        if self.library_name and not self.library_id:
            self.library_id = self.db_manager.libraries.get_id(self.library_name, include_inactive=True)
        if self.all:
            return _filter_by_root(self.db_manager.orthofinder.get_many(include_inactive=True))
        if self.library_id:
            return _filter_by_root(self.db_manager.orthofinder.get_many(library_id=self.library_id, include_inactive=True))
        if self.ref_accessions:
            return _filter_by_root(self.db_manager.orthofinder.get_by_reference_accessions(self.ref_accessions))
        return []

    def _run_orthofinder_worker(self, run_id: int) -> Dict[str, Any]:
        worker = VerifyOrthofinderTask(
            self.db_manager.get_path(),
            self.task_id,
            0,
            data=json.dumps(self.data),
            required_threads=1,
        )
        try:
            return worker._process_orthofinder_run(run_id)
        finally:
            worker.db_manager.close()

    def _process_orthofinder_run(self, run_id: int) -> Dict[str, Any]:
        run = self.db_manager.orthofinder.get(run_id)
        if not run:
            return {"row": None, "summary_row": None}
        orthofinder_id = int(run[0])
        library_id = int(run[1]) if run[1] is not None else None
        location = run[2]
        status_before = self.db_manager.orthofinder.get_status(orthofinder_id)
        location_exists = bool(location and os.path.isdir(str(location)))
        artifact_before = bool(self.db_manager.artifacts.find(owner_type="orthofinder_run", owner_id=orthofinder_id, artifact_type="orthofinder_results_dir"))
        artifact_staled = 0
        accession_list_rewritten = False
        if location_exists:
            self.db_manager.artifacts.register(
                owner_type="orthofinder_run",
                owner_id=orthofinder_id,
                artifact_type="orthofinder_results_dir",
                path=location,
                is_dir=True,
                format="directory",
                metadata={"library_id": library_id},
            )
        else:
            for artifact in self.db_manager.artifacts.find(owner_type="orthofinder_run", owner_id=orthofinder_id, artifact_type="orthofinder_results_dir"):
                if self.repair:
                    self.db_manager.artifacts.set_status(int(artifact[0]), "stale")
                    artifact_staled += 1
        db_accessions = self.db_manager.orthofinder.get_accessions(orthofinder_id)
        accession_file = os.path.join(str(location), "accession_list.txt") if location else None
        file_accessions: List[str] = []
        if accession_file and os.path.isfile(accession_file):
            file_accessions = [line.strip() for line in Path(accession_file).read_text(encoding="utf-8").splitlines() if line.strip()]
        elif self.repair and location_exists and db_accessions:
            Path(accession_file).write_text("\n".join(db_accessions) + "\n", encoding="utf-8")
            file_accessions = list(db_accessions)
            accession_list_rewritten = True
        usable = bool(location_exists and accession_file and os.path.isfile(accession_file))
        target_status = "ready" if usable else "stale"
        if self.repair:
            self.db_manager.orthofinder.set_status(orthofinder_id, target_status)
        status_after = self.db_manager.orthofinder.get_status(orthofinder_id)
        artifact_after = bool(self.db_manager.artifacts.find(owner_type="orthofinder_run", owner_id=orthofinder_id, artifact_type="orthofinder_results_dir"))
        note = "" if sorted(db_accessions) == sorted(file_accessions) else "accession_list_mismatch"
        row = {
            "orthofinder_id": orthofinder_id,
            "library_id": library_id or "",
            "status": target_status,
            "results_dir_exists": "yes" if location_exists else "no",
            "accession_list_exists": "yes" if accession_file and os.path.isfile(accession_file) else "no",
            "db_accessions": len(db_accessions),
            "file_accessions": len(file_accessions),
            "note": note,
        }
        summary_row = None
        if status_before != status_after or accession_list_rewritten or artifact_staled or (not artifact_before and artifact_after):
            summary_row = {
                "orthofinder_id": orthofinder_id,
                "library_id": library_id or "",
                "status_before": str(status_before or ""),
                "status_after": str(status_after or ""),
                "artifact_backfilled": "yes" if (not artifact_before and artifact_after) else "no",
                "artifact_staled": str(artifact_staled),
                "accession_list_rewritten": "yes" if accession_list_rewritten else "no",
                "note": note,
            }
        if note:
            self.log(f"OrthoFinder {orthofinder_id}: accession list differs from DB truth.", "WARNING")
        return {"row": row, "summary_row": summary_row}

    def run(self):
        runs = self._target_runs()
        rows = []
        summary_rows = []
        worker_results, worker_errors = _run_verify_workers(
            self,
            [int(run[0]) for run in runs],
            self._run_orthofinder_worker,
            "orthofinder run",
        )
        for worker_result in worker_results:
            if not worker_result or not worker_result.get("row"):
                continue
            rows.append(worker_result["row"])
            if worker_result.get("summary_row"):
                summary_rows.append(worker_result["summary_row"])
        _write_tsv_report(
            self,
            "verify_orthofinder.tsv",
            ["orthofinder_id", "library_id", "status", "results_dir_exists", "accession_list_exists", "db_accessions", "file_accessions", "note"],
            rows,
        )
        _write_tsv_report(
            self,
            "verify_orthofinder_summary.tsv",
            ["orthofinder_id", "library_id", "status_before", "status_after", "artifact_backfilled", "artifact_staled", "accession_list_rewritten", "note"],
            summary_rows,
        )
        if worker_errors:
            return self.handle_exception(_parallel_error_message("verify-orthofinder worker failures", worker_errors))
        return True


class VerifyTask(Task):
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads=required_threads)
        self.stage = checkpoint if checkpoint is not None else 0

    def _child_payload(self, key: str) -> Dict[str, Any]:
        root_token = self.data.get("root")
        root_kind = None
        if root_token:
            try:
                root_row = self.db_manager.storage.resolve_root_token(root_token)
                root_kind = str(root_row[1])
            except Exception as exc:  # boundary: optional root kind hint failure falls back to generic behavior.
                self.log(f"Failed to resolve root kind for {root_token}: {exc}", "WARNING")
                root_kind = None
        shared = {
            "accessions": self.data.get("accessions") or [],
            "taxid": self.data.get("taxid"),
            "clade": self.data.get("clade"),
            "all": bool(self.data.get("all", False)),
            "downloaded_only": bool(self.data.get("downloaded_only", False)),
            "primary_only": bool(self.data.get("primary_only", False)),
            "after": self.data.get("after"),
            "before": self.data.get("before"),
            "level": self.data.get("level"),
            "filters": self.data.get("filters"),
            "ranks": self.data.get("ranks"),
            "quantities": self.data.get("quantities"),
            "library_id": self.data.get("library_id"),
            "library_name": self.data.get("library_name"),
            "run_id": self.data.get("run_id"),
            "ref_accessions": self.data.get("ref_accessions") or [],
            "repair": bool(self.data.get("repair", False)),
            "reingest": bool(self.data.get("reingest", False)),
            "reingest_all": bool(self.data.get("reingest_all", False)),
            "report_root": self.data.get("report_root"),
        }
        if root_token and root_kind == "genomes" and key in {"verify-assembly", "verify-busco"}:
            shared["root"] = root_token
        elif root_token and root_kind == "libraries" and key == "verify-libraries":
            shared["root"] = root_token
        elif root_token and root_kind == "orthofinder" and key == "verify-orthofinder":
            shared["root"] = root_token
        if key == "verify-assembly":
            shared["discover"] = False
        if key == "verify-busco":
            shared["discover"] = bool(self.data.get("discover", False))
        return shared

    def run(self):
        root_token = self.data.get("root")
        if root_token:
            try:
                root_row = self.db_manager.storage.resolve_root_token(root_token)
            except Exception as exc:  # boundary: required storage root lookup failure becomes this task error.
                return self.handle_exception(str(exc), {"root": root_token})
            root_kind = str(root_row[1])
            if root_kind not in {"genomes", "libraries", "orthofinder"}:
                return self.handle_exception(
                    f"verify does not support scoping by root kind '{root_kind}'. Use a genomes, libraries, or orthofinder root.",
                    {"root": root_token, "kind": root_kind},
                )
        task_map = []
        if bool(self.data.get("include_assembly", True)):
            task_map.append((18, "verify-assembly"))
        if bool(self.data.get("include_libraries", True)):
            task_map.append((28, "verify-libraries"))
        if bool(self.data.get("include_busco", True)):
            task_map.append((20, "verify-busco"))
        if bool(self.data.get("include_orthofinder", True)):
            task_map.append((29, "verify-orthofinder"))

        def queue_verify_subtasks():
            for job_type, key in task_map:
                self.queue_subtask(job_type=job_type, status="P", priority=1, data=self._child_payload(key))
            return True

        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=queue_verify_subtasks,
            done_fn=lambda: self._subtasks_state() == "complete",
            wait_seconds=0,
            retry_key=None,
            max_retries=0,
            incomplete_message_fn=lambda: ("Verify subtasks did not complete.", ""),
            retry_incomplete=False,
        )
        if outcome in ("ERROR", False):
            return outcome
        rows = []
        for child in self.db_manager.tasks.get_subtasks(self.task_id) or []:
            rows.append(
                {
                    "task_id": child[0],
                    "job_type": child[1],
                    "status": child[2],
                    "priority": child[3],
                }
            )
        _write_tsv_report(self, "verify.tsv", ["task_id", "job_type", "status", "priority"], rows)
        return True


VerifyDownloadsTask = VerifyAssemblyTask
