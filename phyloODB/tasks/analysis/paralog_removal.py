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
from pathlib import Path
from typing import Optional, Dict, Any
import concurrent.futures
import threading
import hashlib
import uuid

from ..task import Task
from ..reporting import resolve_report_run_dir
from .blastdb import build_proteome_blastdb_prefix, resolve_proteome_profile_input
from ...database import DBManager
from ...proteome_profile_utils import resolve_profile_selector
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
)

class ParalogRemovalTask(Task):
    '''A class that handles the paralog removal task'''
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=32):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.ref_accessions = list(dict.fromkeys(normalize_accessions(self.data.get("ref_accessions", []))))
        median_raw = list(self.data.get("accessions", []) or [])
        target_raw = list(self.data.get("targets", []) or [])
        self.accessions = list(dict.fromkeys(normalize_accessions(median_raw)))
        self.targets = list(dict.fromkeys(normalize_accessions(target_raw)))
        if not self.accessions and self.targets:
            self.accessions = list(self.targets)
        if not self.targets and self.accessions:
            self.targets = list(self.accessions)
        accessions_set = set(self.accessions)
        missing_targets = [acc for acc in self.targets if acc not in accessions_set]
        if missing_targets:
            raise ValueError(
                "Targets must be included in the paralog-removal accessions set. "
                f"Missing: {', '.join(missing_targets)}"
            )
        self.library_id = self.data.get("library_id")
        self.library_name = self.data.get("library_name")
        self.mode = str(self.data.get("mode") or "median").strip().lower()
        self.percentile = self.data.get("percentile")
        self.bitscore_threshold = self.data.get("bitscore_threshold")
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        self.report_dir = self.data.get("report_dir")
        self.run_label = self.data.get("run_label")
        self.run_id = self.data.get("run_id")
        self.force = self.data.get("force", False)
        self.stage = checkpoint if checkpoint is not None else 0
        self.genomes_checked = self.data.get("genomes_checked", False)
        max_concurrent_raw = self.data.get("max_concurrent")
        if max_concurrent_raw in (None, ""):
            self.max_concurrent = required_threads
        else:
            try:
                self.max_concurrent = max(1, int(max_concurrent_raw))
            except (TypeError, ValueError):
                self.max_concurrent = required_threads
        self.reuse_existing = bool(self.data.get("reuse_existing", False))
        self.rebuild_proteome_dbs = bool(self.data.get("rebuild_proteome_dbs", False))
        self.avoid_unclean_buscos = bool(self.data.get("avoid_unclean_buscos", True))
        self.include_duplicated = bool(self.data.get("include_duplicated", False))
        try:
            self.paralog_report_hit_limit = max(1, int(self.data.get("paralog_report_hit_limit", 5) or 5))
        except (TypeError, ValueError):
            self.paralog_report_hit_limit = 5
        self.blast_threads = 1
        self.blastp_path: Optional[str] = None
        self._report_rows_lock = threading.Lock()

    def _get_missing_downloads(self):
        rows = self.db_manager.genomes.get_downloaded()
        present = [r[0] for r in rows] if rows else []
        missing = [acc for acc in self.ref_accessions if acc not in present] + [acc for acc in self.accessions if acc not in present]
        return missing
    
    def _get_missing_blastDBs(self):
        rows = self.db_manager.filtering.get_blast_dbs()
        self.log(rows, "DEBUG")
        self.log(f"Checking for missing BLAST DBs among reference accessions: {self.ref_accessions}", "DEBUG")
        self.log(f"Existing BLAST DB entries: {rows}", "DEBUG")
        existing_by_accession: dict[str, list[str]] = {}
        for row in rows or []:
            accession = str(row[2]) if row and len(row) > 2 and row[2] is not None else None
            location = str(row[3]) if row and len(row) > 3 and row[3] is not None else None
            if accession and location:
                existing_by_accession.setdefault(accession, []).append(location)
        missing: list[str] = []
        for acc in self.ref_accessions:
            expected_location = self._expected_blastdb_location(acc)
            if expected_location is None:
                missing.append(acc)
                continue
            db_files = [f"{expected_location}.{ext}" for ext in ["phr", "pin", "psq"]]
            if expected_location not in existing_by_accession.get(acc, []) or not all(os.path.exists(path) for path in db_files):
                missing.append(acc)
        present = sorted(existing_by_accession)
        self.log(f"Present BLAST DBs: {present}", "DEBUG")
        self.log(f"Missing BLAST DBs for accessions: {missing}", "DEBUG")
        return missing

    def _get_blastdb_targets(self):
        if self.rebuild_proteome_dbs:
            targets = list(self.ref_accessions)
            self.log(f"Rebuilding proteome BLAST DBs for reference accessions: {targets}", "DEBUG")
            return targets
        return self._get_missing_blastDBs()

    def _blastdb_stage_done(self):
        missing = self._get_missing_blastDBs()
        if missing:
            return False
        if not self.rebuild_proteome_dbs:
            return True

        gen = int(self.data.get("_stage_1_gen", 0) or 0)
        targets = self._get_blastdb_targets()
        if not targets:
            return True
        if gen <= 0:
            return False

        subtasks = self.db_manager.tasks.get_subtasks(self.task_id) or []
        phase_children = []
        for task in subtasks:
            try:
                payload = json.loads(task[6]) if task[6] else {}
            except (TypeError, json.JSONDecodeError):
                payload = {}
            if payload.get("__stage") == 1 and payload.get("__gen") == gen:
                phase_children.append(task)

        if len(phase_children) < len(targets):
            return False
        return all(task[2] == "C" for task in phase_children)
    
    def get_accession_to_blastdb_map(self):
        rows = self.db_manager.filtering.get_blast_dbs()
        accession_to_blastdb = {}
        for row in rows:
            accession = row[2]
            location = row[3]
            expected_location = self._expected_blastdb_location(accession)
            if expected_location and str(location) == str(expected_location):
                accession_to_blastdb[accession] = location
            elif accession not in self.ref_accessions and accession not in accession_to_blastdb:
                accession_to_blastdb[accession] = location
        return accession_to_blastdb
    
    def _queue_blastDB_creation(self):
        targets = self._get_blastdb_targets()
        if targets:
            for acc in targets:
                expected_location = self._expected_blastdb_location(acc)
                self.queue_subtask(
                    job_type=13, status="P", priority=1,
                    data={
                        "accession": acc,
                        "library_id": self.library_id,
                        "proteome_profile": self.proteome_profile,
                        "prefer_proteome_profile": self.prefer_proteome_profile,
                        "output_path": expected_location,
                        "force": True,
                    }
)
            return True
        return False

    def _resolve_reference_profile(self, accession: str) -> tuple[Optional[str], Optional[str]]:
        try:
            selected_profile, _row, profile_path = resolve_proteome_profile_input(
                self.db_manager,
                str(accession),
                proteome_profile=self.proteome_profile,
                prefer_proteome_profile=self.prefer_proteome_profile,
            )
            return str(selected_profile), str(profile_path)
        except (FileNotFoundError, OSError, ValueError) as exc:
            self.log(
                f"Failed to resolve proteome profile for reference accession {accession}: {exc}",
                "WARNING",
            )
            return None, None

    def _expected_blastdb_location(self, accession: str) -> Optional[str]:
        selected_profile, _profile_path = self._resolve_reference_profile(str(accession))
        if selected_profile is None:
            return None
        return build_proteome_blastdb_prefix(self, str(accession), self.library_id, selected_profile)

    def _prepare_paralog_cache_locations(self) -> Dict[str, tuple[str | None, str | None]]:
        locations: Dict[str, tuple[str | None, str | None]] = {}
        fallback_root = self.db_manager.storage.get_root_base(
            "cache",
            fallback=os.path.join(Path(self.db_manager.get_path()).resolve().parent, "cache"),
        ) or "."
        for acc in self.targets:
            genome_path = self.db_manager.genomes.get_path(acc)
            base_dir = None
            if genome_path:
                if os.path.isdir(genome_path):
                    base_dir = genome_path
                else:
                    base_dir = os.path.dirname(genome_path)
            if not base_dir:
                base_dir = os.path.join(fallback_root, "paralog-filtering-cache", acc)
            cache_dir = os.path.join(base_dir, "paralog_cache")
            filename = f"{acc}_paralog_{self.library_id}.json.gz"
            locations[acc] = (cache_dir, os.path.join(cache_dir, filename))
        return locations

    def _build_family_ref_signature(self, family_to_ref_expected: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, str]]:
        signature: Dict[str, Dict[str, str]] = {}
        for family_id, refs in family_to_ref_expected.items():
            normalized_refs = {ref: refs[ref] for ref in sorted(refs)}
            signature[str(family_id)] = normalized_refs
        return signature

    def _get_allowed_family_ids(self) -> Optional[set[str]]:
        """Return the custom-library family set, or None for default BUSCO libraries."""

        if not self.library_id:
            return None
        parent_library_id = self.db_manager.libraries.get_parent_id(self.library_id)
        if not parent_library_id:
            return None
        rows = self.db_manager.cursor.execute(
            "SELECT family_id FROM BUSCO_descriptions WHERE library_id = ?",
            (int(self.library_id),),
        ).fetchall() or []
        return {str(row[0]) for row in rows if row and row[0] is not None}

    @staticmethod
    def _restrict_rows_to_family_set(rows, allowed_family_ids: Optional[set[str]]):
        if not rows or allowed_family_ids is None:
            return rows
        return [row for row in rows if str(row[0]) in allowed_family_ids]

    @staticmethod
    def _dedupe_single_copy_rows(rows):
        if not rows:
            return rows
        deduped: dict[tuple[str, str], tuple] = {}
        for row in rows:
            family_id, _library_id, accession, _status, _sequence, score, _length = row
            key = (str(accession), str(family_id))
            current = deduped.get(key)
            if current is None:
                deduped[key] = row
                continue
            current_score = current[5]
            try:
                row_score = float(score) if score is not None else float("-inf")
            except (TypeError, ValueError):
                row_score = float("-inf")
            try:
                best_score = float(current_score) if current_score is not None else float("-inf")
            except (TypeError, ValueError):
                best_score = float("-inf")
            if row_score > best_score:
                deduped[key] = row
        return list(deduped.values())

    def _load_accession_cache(
        self,
        cache_path: str | None,
        busco_lib_id: int,
        config_signature: str,
    ) -> Optional[Dict[str, Dict[str, Any]]]:
        if not cache_path or not os.path.isfile(cache_path):
            return None
        try:
            with gzip.open(cache_path, "rt") as handle:
                payload = json.load(handle)
        except (OSError, EOFError, gzip.BadGzipFile, json.JSONDecodeError, TypeError) as exc:
            self.log(f"Stage 2: failed to read paralog cache at {cache_path}: {exc}", "DEBUG")
            return None
        if payload.get("version") != 2:
            return None
        if payload.get("target_library_id") != self.library_id:
            return None
        if payload.get("busco_library_id") != busco_lib_id:
            return None
        if payload.get("config_signature") != config_signature:
            return None
        families = {}
        for family_id, info in (payload.get("families") or {}).items():
            if not isinstance(info, dict):
                continue
            clean_flag = info.get("clean")
            if isinstance(clean_flag, bool):
                families[str(family_id)] = {
                    "clean": clean_flag,
                    "selection_signature": info.get("selection_signature"),
                }
        return families

    def _write_accession_cache(
        self,
        cache_dir: str | None,
        cache_path: str | None,
        busco_lib_id: int,
        config_signature: str,
        family_results: Dict[str, Dict[str, Any]],
    ) -> None:
        if not cache_dir or not cache_path:
            return
        try:
            os.makedirs(cache_dir, exist_ok=True)
        except OSError as exc:
            self.log(f"Stage 2: failed to create paralog cache directory {cache_dir}: {exc}", "DEBUG")
            return
        payload = {
            "version": 2,
            "target_library_id": self.library_id,
            "busco_library_id": busco_lib_id,
            "config_signature": config_signature,
            "families": {
                str(fid): {
                    "clean": bool(info.get("clean")),
                    "selection_signature": info.get("selection_signature"),
                }
                for fid, info in family_results.items()
                if isinstance(info, dict)
            },
            "generated_at": datetime.now().isoformat(timespec="seconds"),
        }
        tmp_path = f"{cache_path}.tmp"
        try:
            with gzip.open(tmp_path, "wt") as handle:
                json.dump(payload, handle)
            os.replace(tmp_path, cache_path)
        except (OSError, EOFError, TypeError, ValueError) as exc:
            self.log(f"Stage 2: failed to write paralog cache at {cache_path}: {exc}", "DEBUG")
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except OSError as cleanup_exc:
                self.log(f"Stage 2: failed to remove incomplete paralog cache {tmp_path}: {cleanup_exc}", "DEBUG")

    def get_median(self, scores):
        sorted_scores = sorted(scores)
        n = len(sorted_scores)
        if n == 0:
            return 0
        if n % 2 == 1:
            return sorted_scores[n // 2]
        else:
            mid1 = sorted_scores[n // 2 - 1]
            mid2 = sorted_scores[n // 2]
            return (mid1 + mid2) / 2

    def _quantile(self, scores, quantile: float) -> float:
        ordered = sorted(float(score) for score in scores)
        if not ordered:
            return 0.0
        if len(ordered) == 1:
            return ordered[0]
        q = min(1.0, max(0.0, float(quantile)))
        idx = (len(ordered) - 1) * q
        lower = math.floor(idx)
        upper = math.ceil(idx)
        if lower == upper:
            return ordered[lower]
        fraction = idx - lower
        return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction

    def _selection_params_payload(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "percentile": self.percentile,
            "bitscore_threshold": self.bitscore_threshold,
            "proteome_profile": self.proteome_profile,
            "prefer_proteome_profile": self.prefer_proteome_profile,
            "avoid_unclean_buscos": self.avoid_unclean_buscos,
            "include_duplicated": self.include_duplicated,
            "reuse_existing": self.reuse_existing,
            "rebuild_proteome_dbs": self.rebuild_proteome_dbs,
        }

    def _config_signature(self, busco_lib_id: int) -> str:
        payload = {
            "target_library_id": self.library_id,
            "busco_library_id": busco_lib_id,
            "accessions": sorted(self.accessions),
            "targets": sorted(self.targets),
            "ref_accessions": sorted(self.ref_accessions),
            "selection": self._selection_params_payload(),
        }
        return hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()

    def _selection_threshold_value(self, scores: list[float]) -> float:
        if self.mode == "lower-quartile":
            return self._quantile(scores, 0.25)
        if self.mode == "upper-quartile":
            return self._quantile(scores, 0.75)
        return self._quantile(scores, 0.5)

    def _build_family_selection(
        self,
        *,
        family_scores: dict[str, list[float]],
        ref_rows,
    ) -> tuple[dict[str, dict[str, str]], dict[str, float | None], dict[str, int], dict[str, str]]:
        family_to_ref_expected: dict[str, dict[str, str]] = {}
        family_thresholds: dict[str, float | None] = {}
        family_candidate_counts: dict[str, int] = {}
        family_signatures: dict[str, str] = {}
        ref_rows_by_family: dict[str, list[tuple[str, str, float]]] = {}
        for family_id, _lib_id, acc, _status, seq_id, bitscore, _length in ref_rows:
            if bitscore is None:
                continue
            ref_rows_by_family.setdefault(str(family_id), []).append((str(acc), str(seq_id), float(bitscore)))
        all_family_ids = set(family_scores.keys()) | set(ref_rows_by_family.keys())
        for family_id in sorted(all_family_ids):
            entries = ref_rows_by_family.get(family_id, [])
            selected: dict[str, str] = {}
            threshold_value: float | None = None
            if self.mode in {"median", "lower-quartile", "upper-quartile"}:
                scores = family_scores.get(family_id) or []
                if scores:
                    threshold_value = self._selection_threshold_value(scores)
                    for ref_acc, seq_id, bitscore in entries:
                        if bitscore >= float(threshold_value):
                            selected[ref_acc] = seq_id
            elif self.mode == "percent":
                ranked = sorted(entries, key=lambda item: item[2], reverse=True)
                keep = max(1, int(math.ceil((float(self.percentile or 0.0) / 100.0) * len(ranked)))) if ranked else 0
                for ref_acc, seq_id, _bitscore in ranked[:keep]:
                    selected[ref_acc] = seq_id
                threshold_value = ranked[keep - 1][2] if keep and ranked else None
            elif self.mode == "bitscore":
                threshold_value = float(self.bitscore_threshold or 0.0)
                for ref_acc, seq_id, bitscore in entries:
                    if bitscore >= threshold_value:
                        selected[ref_acc] = seq_id
            family_candidate_counts[family_id] = len(entries)
            family_thresholds[family_id] = threshold_value
            family_to_ref_expected[family_id] = selected
            signature_payload = {
                "family_id": family_id,
                "mode": self.mode,
                "threshold": threshold_value,
                "refs": sorted((acc, seq_id) for acc, seq_id in selected.items()),
                "candidate_count": len(entries),
            }
            family_signatures[family_id] = hashlib.sha1(
                json.dumps(signature_payload, sort_keys=True).encode("utf-8")
            ).hexdigest()
        return family_to_ref_expected, family_thresholds, family_candidate_counts, family_signatures

    def _generate_run_id(self) -> str:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"paralog_{self.library_id}_{stamp}_{uuid.uuid4().hex[:8]}"

    def _write_query_batch(self, query_records: list[dict[str, Any]]) -> str:
        tmp = tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False)
        try:
            for record in query_records:
                tmp.write(f">{record['qseqid']}\n{record['sequence']}\n")
            tmp.flush()
            return tmp.name
        finally:
            tmp.close()

    def _group_hits_by_query(self, blast_output: str) -> dict[str, list[dict[str, Any]]]:
        grouped: dict[str, list[dict[str, Any]]] = {}
        for hit in self._parse_blast_output(blast_output):
            grouped.setdefault(str(hit["qseqid"]), []).append(hit)
        hit_limit = max(1, int(getattr(self, "paralog_report_hit_limit", 5) or 5))
        for hits in grouped.values():
            hits.sort(key=lambda item: (-float(item["bitscore"]), float(item["evalue"])))
            del hits[hit_limit:]
        return grouped

    def _write_tsv_report(self, path: str | Path, rows: list[dict[str, Any]]) -> None:
        headers: list[str] = []
        for row in rows:
            for key in row.keys():
                if key not in headers:
                    headers.append(key)
        if not headers:
            headers = ["note"]
            rows = [{"note": "no rows"}]
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def _blast_hits_report_path(self, report_dir: str, accession: str) -> Path:
        safe_accession = re.sub(r"[^A-Za-z0-9._-]+", "_", str(accession)).strip("._-") or "unknown_accession"
        return Path(report_dir) / "blast_hits" / f"{safe_accession}.tsv"

    def _write_accession_blast_report(
        self,
        report_dir: str,
        accession: str,
        blast_rows: list[dict[str, Any]],
    ) -> None:
        if not report_dir:
            return
        blast_path = self._blast_hits_report_path(report_dir, accession)
        blast_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_tsv_report(blast_path, blast_rows)

    def _write_paralog_reports(
        self,
        report_dir: str,
        summary_rows: list[dict[str, Any]],
        family_rows: list[dict[str, Any]],
        decision_rows: list[dict[str, Any]],
    ) -> None:
        if not report_dir:
            return
        os.makedirs(report_dir, exist_ok=True)

        self._write_tsv_report(Path(report_dir) / "paralog_filtering_summary.tsv", summary_rows)
        self._write_tsv_report(Path(report_dir) / "paralog_filtering_family_selection.tsv", family_rows)
        self._write_tsv_report(Path(report_dir) / "paralog_filtering_decisions.tsv", decision_rows)
        self.log(f"Wrote paralog-filtering report to {report_dir}", "INFO")


    def _run_blastp(self, query_faa, db_path, out_file=None):
        blastp_path = self.blastp_path
        if not blastp_path:
            self.error("BLASTP_PATH was not initialized before running BLASTP.")
            return False

        command = [
            blastp_path,
            "-query", query_faa,
            "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-max_target_seqs", str(self.paralog_report_hit_limit),
        ]
        threads = max(1, int(self.blast_threads)) if getattr(self, "blast_threads", None) else 1
        if threads > 1:
            command.extend(["-num_threads", str(threads)])

        # If an out_file is requested, instruct BLAST to also write it; otherwise capture stdout.
        if out_file:
            command.extend(["-out", out_file])

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode != 0:
            self.error(f"blastp failed: {result.stderr}")
            return False

        # Prefer stdout if BLAST wrote to it; otherwise, if out_file was provided, read and return its contents.
        if result.stdout:
            return result.stdout

        if out_file and os.path.exists(out_file):
            try:
                with open(out_file, "r") as f:
                    return f.read()
            except OSError as e:
                self.error(f"Failed to read BLAST output file {out_file}: {e}")
                return False

        # Nothing returned but command succeeded: return empty string
        return ""

    def _parse_blast_output(self, blast_output):
        hits = []
        for line in blast_output.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 12:
                continue  # Malformed line
            hit = {
                "qseqid": parts[0],
                "sseqid": parts[1],
                "pident": float(parts[2]),
                "length": int(parts[3]),
                "mismatch": int(parts[4]),
                "gapopen": int(parts[5]),
                "qstart": int(parts[6]),
                "qend": int(parts[7]),
                "sstart": int(parts[8]),
                "send": int(parts[9]),
                "evalue": float(parts[10]),
                "bitscore": float(parts[11]),
            }
            hits.append(hit)
        return hits

    def _assert_top_hit_is_sequenced(self, hits, sequence_id, min_score=None, min_identity=None):
        # Ensure the top hit within report is that of sequence_id and optionally above min_score
        if not hits:
            return False
        top_hit = hits[0]
        if top_hit["sseqid"] != sequence_id:
            return False
        if min_identity is not None and top_hit["pident"] < min_identity:
            return False
        if min_score is not None and top_hit["bitscore"] < min_score:
            return False
        return True

    def _iter_query_records(self, fasta_path):
        header = None
        seq_lines = []
        with open(fasta_path) as handle:
            for raw in handle:
                line = raw.rstrip("\n")
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq_lines)
                    header = line[1:].strip()
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            if header is not None:
                yield header, "".join(seq_lines)

    def _pick_duplicate_query_record(self, records, query_id, used):
        query_id = str(query_id or "")
        token = query_id.split()[0] if query_id else ""
        for idx, (header, seq) in enumerate(records):
            if idx in used:
                continue
            first = header.split()[0] if header else ""
            if query_id and (header == query_id or first == query_id or (token and token == first) or query_id in header):
                used.add(idx)
                return header, seq
        for idx, rec in enumerate(records):
            if idx in used:
                continue
            used.add(idx)
            return rec
        return None, None

    def _write_temp_query_record(self, header, sequence):
        tmp = tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False)
        try:
            tmp.write(f">{header}\n{sequence}\n")
            tmp.flush()
            return tmp.name
        finally:
            tmp.close()

    def run(self):
        
        if self.library_name and not self.library_id:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
            if not self.library_id:
                return self.handle_exception(f"Library '{self.library_name}' not found in database.", {"library_name": self.library_name})        
        if not self.ref_accessions:
            if not self.library_id:
                return self.handle_exception(
                    "No reference accessions provided and library_id is unknown.",
                    {"library_name": self.library_name},
                )
            fallback_refs = self.db_manager.libraries.get_reference_assemblies(self.library_id)
            if fallback_refs:
                normalized = list(dict.fromkeys(normalize_accessions(fallback_refs)))
                self.ref_accessions = normalized
                self.data["ref_accessions"] = normalized
                self.log(
                    f"Loaded {len(self.ref_accessions)} reference accessions from library {self.library_id}.",
                    "DEBUG",
                )
            else:
                return self.handle_exception(
                    "No reference accessions provided and none recorded for this library.",
                    {"library_id": self.library_id, "library_name": self.library_name},
                )

        # Check all accesions have been downloaded
        if not self.genomes_checked:
            missing = self._get_missing_downloads()
            if missing:
                return self.handle_exception("Some ref accessions have not been downloaded yet.", {"missing_accessions": missing})
            else:
                self.checkpoint(self.stage, {'genomes_checked': True,})
                self.genomes_checked = True

        # Check there are BUSCO results for all accessions
        parent_library_id = self.db_manager.libraries.get_parent_id(self.library_id)
        busco_lib_id = parent_library_id if parent_library_id else self.library_id
        allowed_family_ids = self._get_allowed_family_ids()
        busco_selector_kwargs = {
            "run_ids": self.data.get("busco_run_ids"),
            "pipeline": self.data.get("busco_pipeline"),
            "input_mode": self.data.get("busco_input_mode"),
            "preferred_pipeline": self.data.get("prefer_busco_pipeline"),
            "preferred_input_mode": self.data.get("prefer_busco_input_mode"),
            "proteome_profile": self.proteome_profile,
            "preferred_proteome_profile": self.prefer_proteome_profile,
            "selection": self.data.get("busco_run_selection") or "primary",
            "purpose": "default",
        }
        target_run_map = self.db_manager.busco.get_effective_run_ids_for_accessions(
            busco_lib_id,
            accessions=self.accessions,
            **busco_selector_kwargs,
        )
        ref_run_map = self.db_manager.busco.get_effective_run_ids_for_accessions(
            busco_lib_id,
            accessions=self.ref_accessions,
            **busco_selector_kwargs,
        )
        hard_busco_selector_requested = any(
            (
                busco_selector_kwargs.get("run_ids"),
                busco_selector_kwargs.get("pipeline"),
                busco_selector_kwargs.get("input_mode"),
            )
        )
        if hard_busco_selector_requested:
            missing_targets = [acc for acc in self.accessions if acc not in target_run_map]
            if missing_targets:
                return self.handle_exception(
                    "No matching BUSCO run found for some paralog-removal targets under the requested BUSCO selector context.",
                    {
                        "missing_accessions": missing_targets,
                        "library_id": busco_lib_id,
                        "pipeline": busco_selector_kwargs.get("pipeline"),
                        "input_mode": busco_selector_kwargs.get("input_mode"),
                        "selection": busco_selector_kwargs.get("selection"),
                        "hint": "Use --busco-run-selection latest if you want the latest matching run rather than only matching primaries.",
                    },
                )
            missing_refs = [acc for acc in self.ref_accessions if acc not in ref_run_map]
            if missing_refs:
                return self.handle_exception(
                    "No matching BUSCO run found for some paralog-removal reference accessions under the requested BUSCO selector context.",
                    {
                        "missing_accessions": missing_refs,
                        "library_id": busco_lib_id,
                        "pipeline": busco_selector_kwargs.get("pipeline"),
                        "input_mode": busco_selector_kwargs.get("input_mode"),
                        "selection": busco_selector_kwargs.get("selection"),
                        "hint": "Use --busco-run-selection latest if you want the latest matching run rather than only matching primaries.",
                    },
                )

        # Stage 1: Ensure proteome BLAST DBs are present
        
        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=self._queue_blastDB_creation,
            done_fn=self._blastdb_stage_done,
            wait_seconds=0,
            retry_key="blastdb_retries",   # <‑ track retry count for this phase
            max_retries=1,                 # allow one additional attempt
            incomplete_message_fn=lambda: (
                f"Proteome BLAST DBs could not be prepared for ({self._get_blastdb_targets() or self._get_missing_blastDBs()}).",
                ""
            ),
            retry_incomplete=True,         # retry when nothing is running but still incomplete
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        self.log("All proteome BLAST DBs are present. Proceeding to paralog removal.", "DEBUG")

        # Stage 2: Paralog removal process
        self.blastp_path = self.db_manager.env.get("BLASTP_PATH")
        if not self.blastp_path:
            return self.handle_exception(
                "BLASTP_PATH is not set in environment variables.",
                {"variable": "BLASTP_PATH"},
            )

        COORD_THRESHOLD = 1e5
        all_busco_rows = self.db_manager.busco.get_family_results_for_library(
            library_id=busco_lib_id,
            accessions=self.accessions,
            status=[1],
            **busco_selector_kwargs,
        )
        all_busco_rows = self._restrict_rows_to_family_set(all_busco_rows, allowed_family_ids)
        all_busco_rows = self._dedupe_single_copy_rows(all_busco_rows)
        if not all_busco_rows:
            return self.handle_exception(
                "No single-copy BUSCO results found for provided paralog-removal accessions.",
                {"library_id": parent_library_id, "accessions": self.accessions},
            )

        unclean_pairs: set[tuple[str, str]] = set()
        if self.avoid_unclean_buscos:
            filter_accessions = list(dict.fromkeys(self.accessions + self.ref_accessions))
            unclean_paralog = self.db_manager.filtering.get_paralog_unclean_families(
                target_library_id=self.library_id,
                busco_library_id=busco_lib_id,
                accessions=filter_accessions,
            )
            supported_decisions = ["support", "weak"]
            accessions_by_run: Dict[str, list[str]] = {}
            active_payload = self.db_manager.env.get(f"ACTIVE_DECONT_RUN_{self.library_id}")
            active_run_id = active_payload.get("run_id") if isinstance(active_payload, dict) else None
            if active_run_id:
                accessions_by_run[str(active_run_id)] = filter_accessions
            else:
                latest_runs = self.db_manager.filtering.get_latest_decontamination_summary(
                    target_library_id=self.library_id,
                    accessions=filter_accessions,
                )
                for acc, (run_id, _decision, _date) in latest_runs.items():
                    if run_id:
                        accessions_by_run.setdefault(str(run_id), []).append(acc)
            unclean_decontam = self.db_manager.filtering.get_decontamination_unclean_families(
                target_library_id=self.library_id,
                busco_library_id=busco_lib_id,
                accessions_by_run=accessions_by_run,
                supported_decisions=supported_decisions,
            ) if accessions_by_run else set()
            unclean_pairs = set(unclean_paralog).union(unclean_decontam)
            if unclean_pairs:
                all_busco_rows = [
                    row for row in all_busco_rows
                    if (str(row[2]), str(row[0])) not in unclean_pairs
                ]
            if not all_busco_rows:
                return self.handle_exception(
                    "No single-copy BUSCO results remain after filtering unclean BUSCOs.",
                    {
                        "library_id": parent_library_id,
                        "accessions": self.accessions,
                        "hint": "Rerun with --use-unclean-buscos to disable filtering.",
                    },
                )

        duplicate_busco_rows = []
        if self.include_duplicated:
            duplicate_busco_rows = self.db_manager.busco.get_family_results_for_library(
                library_id=busco_lib_id,
                accessions=self.targets,
                status=[2],
                **busco_selector_kwargs,
            ) or []
            duplicate_busco_rows = self._restrict_rows_to_family_set(duplicate_busco_rows, allowed_family_ids) or []

        family_scores: dict[str, list[float]] = {}
        acc_to_rows: dict[str, list[tuple[str, str, float, Any]]] = {}
        duplicate_acc_to_rows: dict[str, list[tuple[str, str, Any, Any, Any]]] = {}
        skipped_anomalous = 0
        for family_id, _lib_id, acc, _status, seq_id, bitscore, length in all_busco_rows:
            if bitscore is None:
                continue
            if float(bitscore) > COORD_THRESHOLD:
                skipped_anomalous += 1
                continue
            family_key = str(family_id)
            family_scores.setdefault(family_key, []).append(float(bitscore))
            acc_to_rows.setdefault(str(acc), []).append((family_key, str(seq_id), float(bitscore), length))
        for family_id, _lib_id, acc, status, seq_id, bitscore, length in duplicate_busco_rows:
            duplicate_acc_to_rows.setdefault(str(acc), []).append((str(family_id), str(seq_id), bitscore, length, status))
        if not family_scores:
            return self.handle_exception(
                "No BUSCO family thresholds could be computed for paralog filtering.",
                {"library_id": self.library_id, "mode": self.mode},
            )
        if skipped_anomalous:
            self.log(
                f"Stage 2: skipped {skipped_anomalous} anomalous BUSCO scores above {COORD_THRESHOLD}.",
                "WARNING",
            )

        ref_rows = self.db_manager.busco.get_family_results_for_library(
            library_id=busco_lib_id,
            accessions=self.ref_accessions,
            status=[1],
            **busco_selector_kwargs,
        ) or []
        ref_rows = self._restrict_rows_to_family_set(ref_rows, allowed_family_ids) or []
        ref_rows = self._dedupe_single_copy_rows(ref_rows) or []
        if self.avoid_unclean_buscos and unclean_pairs and ref_rows:
            ref_rows = [
                row for row in ref_rows
                if (str(row[2]), str(row[0])) not in unclean_pairs
            ]
        if not ref_rows:
            return self.handle_exception(
                "No single-copy BUSCO results found for provided reference accessions.",
                {"ref_accessions": self.ref_accessions},
            )

        family_to_ref_expected, family_thresholds, family_candidate_counts, family_signatures = self._build_family_selection(
            family_scores=family_scores,
            ref_rows=ref_rows,
        )
        if not any(family_to_ref_expected.values()):
            return self.handle_exception(
                "No reference proteomes were selected for any BUSCO family.",
                {"mode": self.mode, "ref_accessions": self.ref_accessions},
            )

        report_dir = str(
            resolve_report_run_dir(
                self,
                namespace="paralog-filtering-reports",
                explicit_dir=self.report_dir,
                run_label=self.run_label or self.library_name,
                cache_attr="_paralog_report_dir",
            )
        )
        if not self.run_id:
            self.run_id = self._generate_run_id()
        config_signature = self._config_signature(busco_lib_id)
        self.db_manager.filtering.add_paralog_filtering_run(
            self.run_id,
            self.library_id,
            busco_lib_id,
            json.dumps(self.targets),
            json.dumps(self.accessions),
            json.dumps(self.ref_accessions),
            self.mode,
            json.dumps(self._selection_params_payload(), sort_keys=True),
            config_signature=config_signature,
            run_label=self.run_label,
            report_dir=report_dir,
            datetime=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

        family_selection_rows = []
        for family_id in sorted(family_signatures):
            selected_refs = family_to_ref_expected.get(family_id, {})
            selected_count = len(selected_refs)
            sparse_flag = ""
            if selected_count == 0:
                sparse_flag = "selected_refs=0"
            elif selected_count == 1:
                sparse_flag = "selected_refs=1"
            family_selection_rows.append(
                {
                    "family_id": family_id,
                    "selection_mode": self.mode,
                    "selection_threshold": family_thresholds.get(family_id),
                    "candidate_ref_count": family_candidate_counts.get(family_id, 0),
                    "selected_ref_count": selected_count,
                    "selected_ref_accessions": ",".join(sorted(selected_refs.keys())),
                    "selected_ref_sequences": ",".join(f"{acc}:{seq}" for acc, seq in sorted(selected_refs.items())),
                    "selection_signature": family_signatures.get(family_id),
                    "sparse_flag": sparse_flag,
                }
            )

        cache_locations = self._prepare_paralog_cache_locations()
        cache_lock = threading.Lock()
        per_accession_cache: Dict[str, Dict[str, Dict[str, Any]]] = {acc: {} for acc in self.targets}
        history_lookup: Dict[tuple[str, str, str], bool] = {}
        if self.reuse_existing:
            for row in self.db_manager.filtering.get_paralog_results_history(target_library_id=self.library_id):
                if len(row) < 11:
                    continue
                family_id, stored_busco_lib, _target_lib, acc, clean, _run_id, _sel_count, _sel_thr, _reused, _reason, selection_signature = row
                if int(stored_busco_lib) != int(busco_lib_id) or acc not in per_accession_cache or not selection_signature:
                    continue
                history_lookup[(str(acc), str(family_id), str(selection_signature))] = bool(clean)
            for acc in self.targets:
                _, cache_path = cache_locations.get(acc, (None, None))
                cached_file = self._load_accession_cache(cache_path, busco_lib_id, config_signature)
                if cached_file:
                    per_accession_cache[acc] = cached_file

        accession_to_blastdb = self.get_accession_to_blastdb_map()
        db_file_path = self.db_manager.get_path()
        summary_rows: list[dict[str, Any]] = []
        decision_rows: list[dict[str, Any]] = []
        report_lock = threading.Lock()

        def process_accession(acc: str):
            thread_db = DBManager(db_file_path)
            thread_db.connect()
            try:
                thread_db.set_busco_run_context(
                    pipeline=self.data.get("busco_pipeline"),
                    input_mode=self.data.get("busco_input_mode"),
                    prefer_pipeline=self.data.get("prefer_busco_pipeline"),
                    prefer_input_mode=self.data.get("prefer_busco_input_mode"),
                    proteome_profile=self.proteome_profile,
                    prefer_proteome_profile=self.prefer_proteome_profile,
                    run_ids=self.data.get("busco_run_ids"),
                    selection=self.data.get("busco_run_selection") or "primary",
                )
                effective_busco_run_id = thread_db.busco.get_effective_run_id_for_accession(
                    acc,
                    busco_lib_id,
                    purpose="default",
                )
                local_cache = dict(per_accession_cache.get(acc, {}))
                cache_dir, cache_path = cache_locations.get(acc, (None, None))
                rows = acc_to_rows.get(acc, [])
                duplicate_rows = duplicate_acc_to_rows.get(acc, [])
                now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                clean_count = 0
                dirty_count = 0
                reused_count = 0
                zero_ref_count = 0
                one_ref_count = 0
                local_decisions: list[dict[str, Any]] = []
                local_hits: list[dict[str, Any]] = []
                query_jobs: list[dict[str, Any]] = []
                paralog_rows: list[dict[str, Any]] = []
                paralog_copy_rows: list[dict[str, Any]] = []

                def record_family_result(
                    *,
                    family_id: str,
                    clean_flag: bool,
                    selected_ref_count: int,
                    reason_code: str,
                    selection_signature: str | None,
                    selection_threshold: float | None,
                    reused: bool,
                    query_kind: str = "single_copy",
                    query_id: str | None = None,
                    query_header: str | None = None,
                    query_status: Any = None,
                ) -> None:
                    nonlocal clean_count, dirty_count, reused_count, zero_ref_count, one_ref_count
                    if query_kind == "single_copy":
                        paralog_rows.append(
                            {
                                "family_id": family_id,
                                "busco_library_id": busco_lib_id,
                                "target_library_id": self.library_id,
                                "accession": acc,
                                "run_id": self.run_id,
                                "busco_run_id": effective_busco_run_id,
                                "clean": clean_flag,
                                "selected_ref_count": selected_ref_count,
                                "selection_threshold": selection_threshold,
                                "reused": reused,
                                "reason_code": reason_code,
                                "selection_signature": selection_signature,
                                "date": now,
                            }
                        )
                        local_cache[str(family_id)] = {
                            "clean": bool(clean_flag),
                            "selection_signature": selection_signature,
                        }
                    else:
                        paralog_copy_rows.append(
                            {
                                "family_id": family_id,
                                "library_id": busco_lib_id,
                                "target_library_id": self.library_id,
                                "accession": acc,
                                "run_id": self.run_id,
                                "busco_run_id": effective_busco_run_id,
                                "query_id": query_id,
                                "query_header": query_header,
                                "query_status": query_status,
                                "clean": clean_flag,
                                "selected_ref_count": selected_ref_count,
                                "reused": reused,
                                "reason_code": reason_code,
                                "selection_signature": selection_signature,
                                "date": now,
                            }
                        )
                    if query_kind == "single_copy":
                        if clean_flag:
                            clean_count += 1
                        else:
                            dirty_count += 1
                        if reused:
                            reused_count += 1
                        if selected_ref_count == 0:
                            zero_ref_count += 1
                        elif selected_ref_count == 1:
                            one_ref_count += 1
                    local_decisions.append(
                        {
                            "accession": acc,
                            "family_id": family_id,
                            "query_kind": query_kind,
                            "query_id": query_id or "",
                            "query_header": query_header or "",
                            "query_status": query_status if query_status is not None else "",
                            "clean": int(bool(clean_flag)),
                            "reason_code": reason_code,
                            "selected_ref_count": selected_ref_count,
                            "selection_threshold": selection_threshold,
                            "selection_signature": selection_signature or "",
                            "reused": int(bool(reused)),
                            "run_id": self.run_id,
                        }
                    )

                for family_id, seq_id, _bitscore, _length in rows:
                    selection_signature = family_signatures.get(family_id)
                    selected_ref_count = len(family_to_ref_expected.get(family_id, {}))
                    selection_threshold = family_thresholds.get(family_id)
                    cached_entry = local_cache.get(family_id)
                    if (
                        self.reuse_existing
                        and isinstance(cached_entry, dict)
                        and cached_entry.get("selection_signature") == selection_signature
                        and "clean" in cached_entry
                    ):
                        record_family_result(
                            family_id=family_id,
                            clean_flag=bool(cached_entry.get("clean")),
                            selected_ref_count=selected_ref_count,
                            reason_code="reused_cache",
                            selection_signature=selection_signature,
                            selection_threshold=selection_threshold,
                            reused=True,
                        )
                        continue
                    if self.reuse_existing and selection_signature:
                        prior = history_lookup.get((acc, family_id, selection_signature))
                        if prior is not None:
                            local_cache[family_id] = {"clean": bool(prior), "selection_signature": selection_signature}
                            record_family_result(
                                family_id=family_id,
                                clean_flag=bool(prior),
                                selected_ref_count=selected_ref_count,
                                reason_code="reused_db",
                                selection_signature=selection_signature,
                                selection_threshold=selection_threshold,
                                reused=True,
                            )
                            continue
                    ref_map = family_to_ref_expected.get(family_id, {})
                    if not ref_map:
                        record_family_result(
                            family_id=family_id,
                            clean_flag=False,
                            selected_ref_count=0,
                            reason_code="no_selected_refs",
                            selection_signature=selection_signature,
                            selection_threshold=selection_threshold,
                            reused=False,
                        )
                        continue
                    query_faa = thread_db.busco.get_family_location(
                        family_id,
                        busco_lib_id,
                        acc,
                        sequence_kind="prot",
                    )
                    if not query_faa or not os.path.isfile(query_faa):
                        record_family_result(
                            family_id=family_id,
                            clean_flag=False,
                            selected_ref_count=selected_ref_count,
                            reason_code="missing_query_fasta",
                            selection_signature=selection_signature,
                            selection_threshold=selection_threshold,
                            reused=False,
                        )
                        continue
                    records = list(self._iter_query_records(query_faa))
                    if not records:
                        record_family_result(
                            family_id=family_id,
                            clean_flag=False,
                            selected_ref_count=selected_ref_count,
                            reason_code="empty_query_fasta",
                            selection_signature=selection_signature,
                            selection_threshold=selection_threshold,
                            reused=False,
                        )
                        continue
                    header, sequence = records[0]
                    query_jobs.append(
                        {
                            "qseqid": f"sc|{family_id}|{uuid.uuid4().hex}",
                            "family_id": family_id,
                            "query_kind": "single_copy",
                            "query_id": seq_id,
                            "query_header": header,
                            "query_status": 1,
                            "sequence": sequence,
                            "ref_map": ref_map,
                            "selected_ref_count": selected_ref_count,
                            "selection_threshold": selection_threshold,
                            "selection_signature": selection_signature,
                        }
                    )

                if self.include_duplicated and duplicate_rows:
                    duplicate_groups: dict[str, list[tuple[str, str, Any, Any, Any]]] = {}
                    for entry in duplicate_rows:
                        duplicate_groups.setdefault(str(entry[0]), []).append(entry)
                    for family_id, dup_entries in duplicate_groups.items():
                        selection_signature = family_signatures.get(family_id)
                        ref_map = family_to_ref_expected.get(family_id, {})
                        selected_ref_count = len(ref_map)
                        selection_threshold = family_thresholds.get(family_id)
                        if not ref_map:
                            for _family_id, seq_id, _bitscore, _length, status in dup_entries:
                                record_family_result(
                                    family_id=family_id,
                                    clean_flag=False,
                                    selected_ref_count=0,
                                    reason_code="no_selected_refs",
                                    selection_signature=selection_signature,
                                    selection_threshold=selection_threshold,
                                    reused=False,
                                    query_kind="duplicate_copy",
                                    query_id=seq_id,
                                    query_header=seq_id,
                                    query_status=status,
                                )
                            continue
                        query_faa = thread_db.busco.get_family_location(
                            family_id,
                            busco_lib_id,
                            acc,
                            sequence_kind="prot",
                        )
                        if not query_faa or not os.path.isfile(query_faa):
                            for _family_id, seq_id, _bitscore, _length, status in dup_entries:
                                record_family_result(
                                    family_id=family_id,
                                    clean_flag=False,
                                    selected_ref_count=selected_ref_count,
                                    reason_code="missing_query_fasta",
                                    selection_signature=selection_signature,
                                    selection_threshold=selection_threshold,
                                    reused=False,
                                    query_kind="duplicate_copy",
                                    query_id=seq_id,
                                    query_header=seq_id,
                                    query_status=status,
                                )
                            continue
                        records = list(self._iter_query_records(query_faa))
                        used_records = set()
                        for _family_id, seq_id, _bitscore, _length, status in dup_entries:
                            query_header, query_sequence = self._pick_duplicate_query_record(records, seq_id, used_records)
                            if not query_header or not query_sequence:
                                record_family_result(
                                    family_id=family_id,
                                    clean_flag=False,
                                    selected_ref_count=selected_ref_count,
                                    reason_code="missing_duplicate_query",
                                    selection_signature=selection_signature,
                                    selection_threshold=selection_threshold,
                                    reused=False,
                                    query_kind="duplicate_copy",
                                    query_id=seq_id,
                                    query_header=seq_id,
                                    query_status=status,
                                )
                                continue
                            query_jobs.append(
                                {
                                    "qseqid": f"dup|{family_id}|{seq_id}|{uuid.uuid4().hex}",
                                    "family_id": family_id,
                                    "query_kind": "duplicate_copy",
                                    "query_id": seq_id,
                                    "query_header": query_header,
                                    "query_status": status,
                                    "sequence": query_sequence,
                                    "ref_map": ref_map,
                                    "selected_ref_count": selected_ref_count,
                                    "selection_threshold": selection_threshold,
                                    "selection_signature": selection_signature,
                                }
                            )

                ref_batches: dict[str, list[dict[str, Any]]] = {}
                query_results: dict[str, dict[str, Any]] = {}
                for item in query_jobs:
                    query_results[item["qseqid"]] = {
                        "checked": 0,
                        "failed": False,
                        "blast_failed": False,
                    }
                    for ref_acc in item["ref_map"].keys():
                        blast_db = accession_to_blastdb.get(ref_acc)
                        if not blast_db:
                            continue
                        ref_batches.setdefault(ref_acc, []).append(item)

                for ref_acc, items in ref_batches.items():
                    blast_db = accession_to_blastdb.get(ref_acc)
                    if not blast_db:
                        continue
                    for start in range(0, len(items), 250):
                        chunk = items[start:start + 250]
                        batch_records = [{"qseqid": item["qseqid"], "sequence": item["sequence"]} for item in chunk]
                        query_path = self._write_query_batch(batch_records)
                        try:
                            out = self._run_blastp(query_path, blast_db)
                        finally:
                            try:
                                os.remove(query_path)
                            except OSError as exc:
                                self.log(f"Failed to remove temporary paralog query file {query_path}: {exc}", "WARNING")
                        if out is False:
                            for item in chunk:
                                query_results[item["qseqid"]]["blast_failed"] = True
                                query_results[item["qseqid"]]["failed"] = True
                            continue
                        grouped_hits = self._group_hits_by_query(out)
                        for item in chunk:
                            hits = grouped_hits.get(item["qseqid"], [])
                            expected_seq = item["ref_map"].get(ref_acc)
                            query_results[item["qseqid"]]["checked"] += 1
                            matched = self._assert_top_hit_is_sequenced(hits, expected_seq)
                            if not matched:
                                query_results[item["qseqid"]]["failed"] = True
                            for rank_idx, hit in enumerate(hits, start=1):
                                local_hits.append(
                                    {
                                        "accession": acc,
                                        "family_id": item["family_id"],
                                        "query_kind": item["query_kind"],
                                        "query_id": item["query_id"],
                                        "query_header": item["query_header"],
                                        "reference_accession": ref_acc,
                                        "expected_sequence_id": expected_seq,
                                        "hit_sequence_id": hit["sseqid"],
                                        "rank": rank_idx,
                                        "bitscore": hit["bitscore"],
                                        "identity": hit["pident"],
                                        "evalue": hit["evalue"],
                                        "matched_expected": int(hit["sseqid"] == expected_seq),
                                        "run_id": self.run_id,
                                    }
                                )

                for item in query_jobs:
                    outcome = query_results.get(item["qseqid"], {})
                    checked = int(outcome.get("checked", 0) or 0)
                    if checked <= 0:
                        clean_flag = False
                        reason_code = "no_checked_refs"
                    elif outcome.get("blast_failed"):
                        clean_flag = False
                        reason_code = "blast_failed"
                    elif outcome.get("failed"):
                        clean_flag = False
                        reason_code = "top_hit_mismatch"
                    else:
                        clean_flag = True
                        reason_code = "clean"
                    record_family_result(
                        family_id=item["family_id"],
                        clean_flag=clean_flag,
                        selected_ref_count=item["selected_ref_count"],
                        reason_code=reason_code,
                        selection_signature=item["selection_signature"],
                        selection_threshold=item["selection_threshold"],
                        reused=False,
                        query_kind=item["query_kind"],
                        query_id=item["query_id"],
                        query_header=item["query_header"],
                        query_status=item["query_status"],
                    )

                summary_row = {
                    "accession": acc,
                    "run_id": self.run_id,
                    "selection_mode": self.mode,
                    "percentile": self.percentile if self.percentile is not None else "",
                    "bitscore_threshold": self.bitscore_threshold if self.bitscore_threshold is not None else "",
                    "families_tested": len(rows),
                    "duplicate_queries_tested": len([item for item in query_jobs if item["query_kind"] == "duplicate_copy"]),
                    "duplicate_clean_count": len([row for row in local_decisions if row["query_kind"] == "duplicate_copy" and int(row["clean"]) == 1]),
                    "duplicate_dirty_count": len([row for row in local_decisions if row["query_kind"] == "duplicate_copy" and int(row["clean"]) == 0]),
                    "clean_count": clean_count,
                    "dirty_count": dirty_count,
                    "reused_count": reused_count,
                    "families_zero_selected_refs": zero_ref_count,
                    "families_one_selected_ref": one_ref_count,
                }
                self.log(f"Paralog filtering for {acc}: {clean_count}/{len(rows)} single-copy families clean.", "DEBUG")
                if paralog_rows and not thread_db.filtering.add_paralog_filtering_results(paralog_rows):
                    raise RuntimeError(f"Failed to persist paralog-filtering rows for {acc}")
                if paralog_copy_rows and not thread_db.filtering.add_paralog_filtering_copy_results(paralog_copy_rows):
                    raise RuntimeError(f"Failed to persist paralog-filtering copy rows for {acc}")
                with cache_lock:
                    per_accession_cache[acc] = local_cache
                    self._write_accession_cache(
                        cache_dir,
                        cache_path,
                        busco_lib_id,
                        config_signature,
                        local_cache,
                    )
                self._write_accession_blast_report(report_dir, acc, local_hits)
                with report_lock:
                    summary_rows.append(summary_row)
                    decision_rows.extend(local_decisions)
                return acc, 0
            finally:
                try:
                    thread_db.close()
                except Exception as exc:  # boundary: worker cleanup failure is logged after primary work completes/fails.
                    self.log(f"Failed to close worker database connection for paralog filtering: {exc}", "WARNING")

        effective_max = self.max_concurrent if self.max_concurrent and self.max_concurrent > 0 else self.REQUIRED_THREADS
        max_workers = max(1, min(len(self.targets), effective_max, self.REQUIRED_THREADS))
        self.blast_threads = max(1, math.ceil(self.REQUIRED_THREADS / max_workers)) if max_workers else self.REQUIRED_THREADS
        self.log(
            f"Stage 2: running paralog filtering with mode={self.mode}, max_workers={max_workers}, "
            f"blast_threads={self.blast_threads}, report_hit_limit={self.paralog_report_hit_limit}.",
            "INFO",
        )
        results = []
        progress_total = len(self.targets)
        progress_done = 0
        progress_next = 10
        progress_lock = threading.Lock()

        def _log_progress():
            nonlocal progress_done, progress_next
            if progress_total <= 0:
                return
            with progress_lock:
                progress_done += 1
                pct = int((progress_done / progress_total) * 100)
                while progress_next <= pct and progress_next <= 100:
                    self.log(f"Paralog filtering progress: {progress_next}% ({progress_done}/{progress_total})", "INFO")
                    progress_next += 10

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="PR_") as executor:
            futures = {executor.submit(process_accession, acc): acc for acc in self.targets}
            for future in concurrent.futures.as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:  # boundary: one accession failure is aggregated while independent targets continue.
                    acc = futures[future]
                    self.error(f"Paralog filtering failed for {acc}: {exc}")
                    results.append((acc, 2))
                _log_progress()

        failed = [acc for acc, code in results if code != 0]
        self._write_paralog_reports(
            report_dir,
            summary_rows,
            family_selection_rows,
            decision_rows,
        )
        payload = {
            "run_id": self.run_id,
            "run_label": self.run_label,
            "selection_mode": self.mode,
            "report_dir": report_dir,
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        try:
            self.db_manager.env.set(f"ACTIVE_PARALOG_RUN_{self.library_id}", payload)
        except Exception as exc:  # boundary: active-run convenience pointer is optional; results have already been persisted.
            self.log(f"Failed to persist active paralog run pointer for library {self.library_id}: {exc}", "WARNING")
        if failed:
            preview = ", ".join(failed[:10])
            suffix = "" if len(failed) <= 10 else ", ..."
            self.log(
                f"Paralog filtering completed with issues for {len(failed)}/{len(self.targets)} targets: {preview}{suffix}",
                "WARNING",
            )
        else:
            self.log(f"Paralog filtering completed for {len(self.targets)} targets.", "INFO")
        return True
