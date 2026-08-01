import tempfile
import time
import json
import logging
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
from typing import Optional, Dict, Any
import concurrent.futures
import threading
import hashlib
import uuid
from pathlib import Path

from ..task import Task
from ..reporting import sanitize_report_label
from ...logging_utils import DEFAULT_FORMAT
from ...database import DBManager
from ...proteome_profile_utils import resolve_profile_selector
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
    RANK_HIERARCHY,
)


class _TaskLogMirrorFormatter(logging.Formatter):
    """Formatter that tolerates task extras missing from some records."""

    def format(self, record: logging.LogRecord) -> str:
        for key in ("task_id", "task_type", "task_name", "stage", "task_ref"):
            if not hasattr(record, key):
                setattr(record, key, None)
        if getattr(record, "task_ref", None) in (None, ""):
            task_name = getattr(record, "task_name", None) or "Task"
            task_id = getattr(record, "task_id", None)
            stage = getattr(record, "stage", None)
            if task_id in (None, ""):
                record.task_ref = str(task_name)
            else:
                stage_value = 0 if stage in (None, "") else stage
                record.task_ref = f"{task_name}:{task_id}.{stage_value}"
        return super().format(record)


class _TaskIdFilter(logging.Filter):
    def __init__(self, task_id: int | None):
        super().__init__()
        self.task_id = None if task_id is None else str(task_id)

    def filter(self, record: logging.LogRecord) -> bool:
        if self.task_id is None:
            return False
        return str(getattr(record, "task_id", "")) == self.task_id

class ExportLibraryTask(Task):
    '''Export BUSCO family FASTAs for a library/taxid selection.

    BUSCOs are pulled from the parent lineage if using a custom library, filtered to the
    custom library's cleaned set, and further pruned using paralog removal results plus
    a chosen decontamination run (explicit run id/label or ACTIVE_DECONT_RUN_<library_id>)
    when those results are available. Missing filtering does not block export unless the
    caller explicitly sets require_paralog_filtering/require_decontamination.
    '''
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=4):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        # Accessions and taxids to use.
        self.accession = self.data.get("accession")
        self.taxid = self.data.get("taxid")
        # If its a rule it will have the quanitity and rank - i.e. get me the top 2 species within each family within the taxid provided
        self.rule_quantity = self.data.get("quantity", None)
        self.rule_rank = self.data.get("rank", None)
        if isinstance(self.rule_rank, str) and not self.rule_rank.strip():
            self.rule_rank = None
        self.clade = self.data.get("clade")
        raw_required = self.data.get("require") or []
        if isinstance(raw_required, str):
            raw_required = [raw_required]
        self.require_groups = [str(token).strip() for token in raw_required if str(token).strip()]
        self.resolved_require_clauses: list[dict[str, Any]] = []

        self.out_dir = str(self.data.get("out_dir")) if self.data.get("out_dir") else None

        self.library_name = self.data.get("library_name") or self.data.get("lineage") or "metazoa_odb10"
        self.sequence_type = str(self.data.get("sequence_type") or "protein").strip().lower()
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        raw_profiles = self.data.get("proteome_profiles") or []
        self.proteome_profiles = [str(token).strip() for token in raw_profiles if str(token).strip()]
        self.protein_only = self.data.get("protein_only") # Only accept accessions from rules/accessions list that have protein files
        self.library_id = self.data.get("library_id", None) # Library to use for BUSCO selection/final list for export
        self.rerun = self.data.get("rerun", False) # Redo analysis if the results are already in the result folder.
        self.stage = checkpoint if checkpoint is not None else 0
        self.decont_run_id = self.data.get("decont_run_id")
        self.decont_run_label = self.data.get("decont_run_label")
        self.decontamination_run = self.data.get("decontamination_run")
        self.paralog_run_id = self.data.get("paralog_run_id") or self.data.get("use_paralog_run")
        self.disable_paralog_filter = bool(self.data.get("disable_paralog_filter", False))
        self.disable_decont_filter = bool(self.data.get("disable_decont_filter", False))
        self.require_paralog_filtering = bool(self.data.get("require_paralog_filtering", False))
        self.require_decontamination = bool(self.data.get("require_decontamination", False))
        self.min_occupancy = float(self.data.get("min_occupancy", 0.5))
        self.min_taxa_occupancy = float(self.data.get("min_taxa_occupancy", 0.3))
        self.min_completeness = self.data.get("min_completeness")
        self.min_single_copy_complete = self.data.get("min_single_copy_complete")
        self.write_lineage_csv = bool(self.data.get("write_lineage_csv", True))
        self.write_busco_report = bool(self.data.get("write_busco_report", True))
        self.write_busco_family_matrix = bool(self.data.get("write_busco_family_matrix", True))
        self.lineage_csv_path = self.data.get("lineage_csv_path")
        self.busco_report_path = self.data.get("busco_report_path")
        self.busco_family_matrix_path = self.data.get("busco_family_matrix_path")
        self.report_dir = self.out_dir
        self.busco_report_extended = bool(self.data.get("busco_report_extended", False))
        self.rescue_duplicates = bool(self.data.get("rescue_duplicates", False))
        self.family_ids = self._normalize_family_ids(self.data.get("family_ids"))
        self.include_duplicated = bool(self.data.get("include_duplicated", False))
        self.include_paralog_filtering_in_score = self.data.get("include_paralog_filtering_in_score")
        self.include_decontamination_in_score = self.data.get("include_decontamination_in_score")
        self.allow_ambiguous_contaminants = self.data.get("allow_ambiguous_contaminants")
        self.strict_decontamination = self.data.get("strict_decontamination")
        if self.disable_decont_filter and self.include_decontamination_in_score is None:
            self.include_decontamination_in_score = False
        self.retain_headers = bool(self.data.get("retain_headers", False))
        self.header_template = self.data.get("header")
        self.header_rank = self.data.get("header_rank")
        self._export_log_handler: Optional[logging.Handler] = None
        self._export_log_path: Optional[str] = None
        # self.log(f"BUSCOTask init accession={self.accession} lineage={self.lineage} format={self.format} stage={self.stage}")

    def _default_export_dir_name(self) -> str:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        suffix = sanitize_report_label(self.library_name or self.library_id or self.taxid or "export")
        folder = f"task_{self.task_id}_{stamp}"
        if suffix:
            folder = f"{folder}_{suffix}"
        return folder

    def _ensure_export_paths(self) -> bool:
        if not self.out_dir:
            try:
                exports_root = self.db_manager.storage.get_root_base("exports")
            except Exception as exc:  # boundary: missing exports root is reported by caller as path setup failure.
                self.log(f"Failed to resolve exports storage root: {exc}", "WARNING")
                exports_root = None
            if not exports_root:
                return False
            self.out_dir = os.path.join(str(exports_root), self._default_export_dir_name())
            self.data["out_dir"] = self.out_dir
        self.report_dir = self.out_dir
        return True

    def _enable_export_log_copy(self) -> bool:
        if self._export_log_handler is not None or not self.out_dir:
            return False
        log_path = os.path.join(self.out_dir, "export_task.log")
        try:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
        except OSError as exc:
            self.log(f"Failed to create export log directory for {log_path}: {exc}", "WARNING")
        try:
            base_logger = self.logger.logger if hasattr(self.logger, "logger") else self.logger
            handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
            handler.setLevel(logging.INFO)
            handler.setFormatter(_TaskLogMirrorFormatter(DEFAULT_FORMAT))
            handler.addFilter(_TaskIdFilter(self.task_id))
            base_logger.addHandler(handler)
            self._export_log_handler = handler
            self._export_log_path = log_path
            return True
        except Exception as exc:  # boundary: task log mirroring is optional.
            self.log(f"Failed to enable export log mirror {log_path}: {exc}", "WARNING")
            return False

    def _disable_export_log_copy(self) -> None:
        handler = self._export_log_handler
        if handler is None:
            return
        try:
            base_logger = self.logger.logger if hasattr(self.logger, "logger") else self.logger
            base_logger.removeHandler(handler)
        except Exception as exc:  # boundary: task log mirror cleanup must not fail export finalization.
            self.log(f"Failed to detach export log mirror: {exc}", "WARNING")
        try:
            handler.close()
        except Exception as exc:  # boundary: task log mirror cleanup must not fail export finalization.
            self.log(f"Failed to close export log mirror: {exc}", "WARNING")
        self._export_log_handler = None
        self._export_log_path = None

    def _taxa_occupancy_report_path(self) -> str:
        if self.busco_report_path:
            base, ext = os.path.splitext(self.busco_report_path)
            return f"{base}_taxa_occupancy.tsv" if ext else f"{self.busco_report_path}_taxa_occupancy.tsv"
        return os.path.join(self.report_dir, "taxa_occupancy.tsv")

    def _export_sequence_kind(self) -> str:
        return "nucl" if self.sequence_type == "nucleotide" else "prot"

    def _location_matches_sequence_kind(self, path: Optional[str], sequence_kind: str) -> bool:
        if not path:
            return False
        lowered = str(path).lower()
        if sequence_kind == "nucl":
            return lowered.endswith((".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz"))
        return lowered.endswith((".faa", ".faa.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz"))

    def _parse_export_fasta(self, path: str) -> list[tuple[str, str]]:
        records: list[tuple[str, str]] = []
        header: Optional[str] = None
        seq_lines: list[str] = []
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        records.append((header, "".join(seq_lines)))
                    header = line
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_lines)))
        return records

    def _normalize_family_ids(self, raw_value: Any) -> list[str]:
        if raw_value is None:
            return []
        if isinstance(raw_value, str):
            tokens = [raw_value]
        else:
            tokens = [str(item).strip() for item in raw_value if str(item).strip()]
        if len(tokens) == 1:
            candidate = Path(tokens[0]).expanduser()
            if candidate.is_file():
                loaded: list[str] = []
                for line in candidate.read_text(encoding="utf-8").splitlines():
                    loaded.extend(part.strip() for part in line.split(",") if part.strip())
                tokens = loaded
        normalized: list[str] = []
        seen: set[str] = set()
        for token in tokens:
            for part in str(token).split(","):
                family_id = part.strip()
                if family_id and family_id not in seen:
                    seen.add(family_id)
                    normalized.append(family_id)
        return normalized

    def _resolve_export_sequence_source(self, row: dict[str, Any], sequence_kind: str) -> Optional[str]:
        candidates = [
            row.get("artifact_path"),
            row.get("artifact_location"),
        ]
        legacy_location = row.get("legacy_location")
        if legacy_location and self._location_matches_sequence_kind(str(legacy_location), sequence_kind):
            candidates.append(legacy_location)
        for candidate in candidates:
            if candidate and os.path.exists(str(candidate)):
                return str(candidate)
        return None

    def _resolve_records_for_export_row(
        self,
        row: dict[str, Any],
        *,
        sequence_kind: str,
    ) -> list[tuple[str, str]]:
        source_path = self._resolve_export_sequence_source(row, sequence_kind)
        if not source_path:
            return []
        records = self._parse_export_fasta(source_path)
        if not records:
            return []
        if int(row.get("status") or 0) == 1 and len(records) == 1:
            return records

        sequence_id = str(row.get("sequence_id") or "").strip()
        if not sequence_id:
            return records if len(records) == 1 else []
        wanted_tokens = {
            sequence_id,
            sequence_id.split()[0],
        }
        matched = []
        for header, sequence in records:
            raw_header = header[1:] if header.startswith(">") else header
            first_token = raw_header.split()[0] if raw_header else ""
            if raw_header in wanted_tokens or first_token in wanted_tokens:
                matched.append((header, sequence))
        return matched

    def _load_strict_export_family_entries(
        self,
        selected_busco_runs: dict[str, dict[str, Any]],
        *,
        library_id: int,
    ) -> dict[str, list[tuple[str, str, str]]]:
        selected_run_map: dict[str, int] = {}
        missing_runs: list[str] = []
        for accession, meta in (selected_busco_runs or {}).items():
            run_id = meta.get("run_id")
            if run_id is None:
                missing_runs.append(str(accession))
                continue
            selected_run_map[str(accession)] = int(run_id)
        if missing_runs:
            raise ValueError(
                "selected accessions lack strict BUSCO run metadata: " + ",".join(sorted(missing_runs)[:25])
            )
        sequence_kind = self._export_sequence_kind()
        missing_family_rows = [
            accession
            for accession, run_id in selected_run_map.items()
            if self.db_manager.busco.count_run_family_rows(int(run_id)) <= 0
        ]
        if missing_family_rows:
            raise ValueError(
                "selected BUSCO runs are missing canonical family rows for: "
                + ",".join(sorted(missing_family_rows)[:25])
            )
        family_filter = set(self.family_ids)
        export_statuses = [1, 2] if self.include_duplicated else [1]
        row_map = self.db_manager.busco.get_export_family_rows(
            run_ids=list(selected_run_map.values()),
            sequence_kind=sequence_kind,
            status=export_statuses,
            accessions=list(selected_run_map.keys()),
        )
        rows_by_accession: dict[str, list[dict[str, Any]]] = {}
        for row in row_map:
            accession = str(row.get("accession"))
            run_id = int(row.get("run_id"))
            if selected_run_map.get(accession) != run_id:
                continue
            family_id = str(row.get("family_id"))
            if family_filter and family_id not in family_filter:
                continue
            rows_by_accession.setdefault(accession, []).append(row)
        missing_exportable = [
            accession for accession in selected_run_map if not rows_by_accession.get(accession)
        ]
        if missing_exportable and not family_filter:
            raise ValueError(
                "selected BUSCO runs are missing exportable BUSCO rows for: "
                + ",".join(sorted(missing_exportable)[:25])
            )
        family_entries: dict[str, list[tuple[str, str, str]]] = {}
        invalid_sources: list[str] = []
        for accession, rows in rows_by_accession.items():
            for row in rows:
                family_id = str(row.get("family_id"))
                run_id = int(row.get("run_id"))
                resolved_records = self._resolve_records_for_export_row(
                    row,
                    sequence_kind=sequence_kind,
                )
                if not resolved_records:
                    invalid_sources.append(f"{accession}:{family_id}:run{run_id}:missing_source")
                    continue
                if int(row.get("status") or 0) == 1 and len(resolved_records) != 1:
                    invalid_sources.append(
                        f"{accession}:{family_id}:run{run_id}:record_count={len(resolved_records)}"
                    )
                    continue
                for header, sequence in resolved_records:
                    if not sequence:
                        invalid_sources.append(f"{accession}:{family_id}:run{run_id}:empty_sequence")
                        continue
                    family_entries.setdefault(family_id, []).append((accession, header, sequence))
        if invalid_sources:
            preview = ",".join(invalid_sources[:25])
            raise ValueError(f"invalid BUSCO sequence sources: {preview}")
        return family_entries

    def _active_decont_env_key(self, library_id: int):
        return f"ACTIVE_DECONT_RUN_{library_id}"

    def _load_active_decont_run(self, library_id: int):
        try:
            return self.db_manager.env.get(self._active_decont_env_key(library_id))
        except Exception as exc:  # boundary: active decontamination pointer is optional filter metadata.
            self.log(f"Failed to load active decontamination run for library {library_id}: {exc}", "WARNING")
            return None

    def _list_decont_runs(self, library_id: int, run_label: str | None = None):
        rows = []
        try:
            sql = """
                SELECT run_id, run_label, targets_json, date
                FROM Decontamination_Runs
                WHERE target_library_id = ?
            """
            params = [library_id]
            if run_label:
                sql += " AND run_label = ?"
                params.append(run_label)
            self.db_manager.cursor.execute(sql, tuple(params))
            rows = self.db_manager.cursor.fetchall() or []
        except Exception as exc:  # boundary: decontamination run listing is optional filter metadata.
            self.log(f"Failed to list decontamination runs for library {library_id}: {exc}", "WARNING")
            rows = []
        return rows

    def _resolve_decont_run(self, library_id: int, selected: list[str], parent_id: Optional[int] = None):
        """Pick a decont run deterministically without requiring manual run ids."""
        selected_set = set(selected)
        # 1) Explicit run id
        if self.decont_run_id:
            if self.db_manager.filtering.get_decontamination_run(self.decont_run_id):
                self.log(f"Using explicit decont run_id={self.decont_run_id}", "DEBUG")
                return self.decont_run_id
            self.log(f"Explicit decont run_id {self.decont_run_id} not found; falling back.", "WARNING")
        # 2) Explicit run label
        if self.decont_run_label:
            rows = self._list_decont_runs(library_id, self.decont_run_label)
            if rows:
                chosen = sorted(rows, key=lambda r: r[3] or "", reverse=True)[0]
                self.log(f"Using decont run with label '{self.decont_run_label}': {chosen[0]}", "DEBUG")
                return chosen[0]
            if parent_id:
                rows = self._list_decont_runs(parent_id, self.decont_run_label)
                if rows:
                    chosen = sorted(rows, key=lambda r: r[3] or "", reverse=True)[0]
                    self.log(
                        f"Using parent decont run with label '{self.decont_run_label}': {chosen[0]}",
                        "INFO",
                    )
                    return chosen[0]
            self.log(f"No decont runs found with label '{self.decont_run_label}'.", "WARNING")
        # 3) Active env var
        active = self._load_active_decont_run(library_id)
        if active and active.get("run_id"):
            self.log(f"Using active decont run from env {self._active_decont_env_key(library_id)}: {active}", "DEBUG")
            return active.get("run_id")
        if parent_id:
            active_parent = self._load_active_decont_run(parent_id)
            if active_parent and active_parent.get("run_id"):
                self.log(
                    f"Using active parent decont run from env {self._active_decont_env_key(parent_id)}: {active_parent}",
                    "INFO",
                )
                return active_parent.get("run_id")
        # 4) Exact target match (newest)
        candidates = []
        for run_id, run_label, targets_json, dt in self._list_decont_runs(library_id):
            try:
                targets = set(json.loads(targets_json or "[]"))
            except (TypeError, json.JSONDecodeError):
                continue
            if targets == selected_set:
                candidates.append((run_id, dt or "", run_label))
        if candidates:
            run_id, dt, lbl = sorted(candidates, key=lambda r: r[1], reverse=True)[0]
            self.log(
                f"Using decont run by target match (run_id={run_id}, label={lbl}, date={dt}). "
                f"Override with decont_run_label or decont_run_id if desired.",
                "INFO",
            )
            return run_id
        if parent_id:
            parent_candidates = []
            for run_id, run_label, targets_json, dt in self._list_decont_runs(parent_id):
                try:
                    targets = set(json.loads(targets_json or "[]"))
                except (TypeError, json.JSONDecodeError):
                    continue
                if targets == selected_set:
                    parent_candidates.append((run_id, dt or "", run_label))
            if parent_candidates:
                run_id, dt, lbl = sorted(parent_candidates, key=lambda r: r[1], reverse=True)[0]
                self.log(
                    f"Using parent decont run by target match (run_id={run_id}, label={lbl}, date={dt}). "
                    f"Override with decont_run_label or decont_run_id if desired.",
                    "INFO",
                )
                return run_id
        # 5) Latest run for library
        latest = sorted(self._list_decont_runs(library_id), key=lambda r: r[3] or "", reverse=True)
        if latest:
            run_id, run_label, _targets_json, dt = latest[0]
            self.log(
                f"Using latest decont run for library (run_id={run_id}, label={run_label}, date={dt}).",
                "INFO",
            )
            return run_id
        if parent_id:
            latest_parent = sorted(self._list_decont_runs(parent_id), key=lambda r: r[3] or "", reverse=True)
            if latest_parent:
                run_id, run_label, _targets_json, dt = latest_parent[0]
                self.log(
                    f"Using latest decont run for parent library (run_id={run_id}, label={run_label}, date={dt}).",
                    "INFO",
                )
                return run_id
        # 6) Fallback: infer run id from latest decontamination summaries
        try:
            summaries = self.db_manager.filtering.get_latest_decontamination_summary_with_fallback(
                target_library_id=library_id,
                parent_library_id=parent_id,
                accessions=list(selected_set),
                run_id=None,
            )
        except Exception as exc:  # boundary: fallback decontamination summaries are optional filter metadata.
            self.log(f"Failed to infer decontamination run from summaries: {exc}", "WARNING")
            summaries = {}
        if summaries:
            counts: dict[str, int] = {}
            latest_date: dict[str, str] = {}
            for _acc, (run_id, _decision, date) in summaries.items():
                rid = str(run_id)
                counts[rid] = counts.get(rid, 0) + 1
                date_val = str(date or "")
                if date_val and date_val > latest_date.get(rid, ""):
                    latest_date[rid] = date_val
            if counts:
                best = sorted(counts.items(), key=lambda kv: (kv[1], latest_date.get(kv[0], "")), reverse=True)[0][0]
                self.log(
                    f"Inferred decont run from summaries (run_id={best}, coverage={counts.get(best)}).",
                    "INFO",
                )
                return best
        self.log(
            "No matching decontamination run found; proceeding without decontamination filtering.",
            "WARNING",
        )
        return None

    def _fmt_busco_value(self, value: Optional[float]) -> str:
        if value is None:
            return "NA"
        try:
            return f"{float(value):.2f}"
        except (TypeError, ValueError):
            return str(value)

    def _fmt_header_bitscore(self, value: Optional[float]) -> str:
        if value is None:
            return ""
        try:
            return f"{float(value):g}"
        except (TypeError, ValueError):
            return str(value)

    def _fetch_accession_taxids(self, accessions: list[str]) -> dict[str, int]:
        mapping: dict[str, int] = {}
        if not accessions:
            return mapping
        for chunk in [accessions[i : i + 900] for i in range(0, len(accessions), 900)]:
            placeholders = ",".join("?" for _ in chunk)
            self.db_manager.cursor.execute(
                f"SELECT accession, taxid FROM Genome WHERE accession IN ({placeholders})",
                tuple(chunk),
            )
            for accession, taxid in self.db_manager.cursor.fetchall() or []:
                if accession is None or taxid is None:
                    continue
                mapping[str(accession)] = int(taxid)
        return mapping

    def _split_require_clause(self, token: str) -> list[str]:
        raw = str(token or "").strip()
        if not raw:
            return []
        if raw.startswith("(") and raw.endswith(")"):
            raw = raw[1:-1].strip()
        if "(" in raw or ")" in raw:
            raise ValueError(f"Invalid require clause '{token}': unmatched parentheses.")
        parts = [part.strip() for part in raw.split("|") if part.strip()]
        if not parts:
            raise ValueError(f"Invalid require clause '{token}'.")
        return parts

    def _resolve_required_accession(self, candidate: str) -> Optional[str]:
        raw = str(candidate or "").strip()
        if not raw:
            return None
        lowered = raw.lower()
        if lowered.startswith("acc:"):
            raw = raw[4:].strip()
        elif lowered.startswith("accession:"):
            raw = raw[len("accession:"):].strip()
        else:
            return None
        if not raw:
            return None
        normalized = normalize_accessions([raw])
        accession = normalized[0] if normalized else raw
        self.db_manager.cursor.execute(
            "SELECT accession FROM Genome WHERE accession = ? LIMIT 1",
            (accession,),
        )
        row = self.db_manager.cursor.fetchone()
        if row is None or row[0] is None:
            return None
        return str(row[0])

    def _resolve_untyped_required_accession(self, candidate: str) -> Optional[str]:
        raw = str(candidate or "").strip()
        if not raw:
            return None
        normalized = normalize_accessions([raw])
        accession = normalized[0] if normalized else raw
        self.db_manager.cursor.execute(
            "SELECT accession FROM Genome WHERE accession = ? LIMIT 1",
            (accession,),
        )
        row = self.db_manager.cursor.fetchone()
        if row is None or row[0] is None:
            return None
        return str(row[0])

    def _format_required_alternative_label(self, alternative: dict[str, Any]) -> str:
        if str(alternative.get("kind") or "").lower() == "accession":
            accession = str(alternative.get("accession") or alternative.get("token") or "")
            return f"{accession}[acc]"
        return f"{alternative['name']}[{alternative['taxid']}]"

    def _clause_satisfied(
        self,
        clause: dict[str, Any],
        *,
        present_taxids: set[int],
        present_accessions: set[str],
    ) -> bool:
        for alternative in clause.get("alternatives") or []:
            if str(alternative.get("kind") or "").lower() == "accession":
                accession = str(alternative.get("accession") or "")
                if accession and accession in present_accessions:
                    return True
                continue
            taxid = alternative.get("taxid")
            if taxid is not None and int(taxid) in present_taxids:
                return True
        return False

    def _missing_required_clause_labels(
        self,
        clauses: list[dict[str, Any]],
        *,
        present_taxids: set[int],
        present_accessions: set[str],
    ) -> list[str]:
        missing: list[str] = []
        for clause in clauses:
            if self._clause_satisfied(
                clause,
                present_taxids=present_taxids,
                present_accessions=present_accessions,
            ):
                continue
            alternatives = clause.get("alternatives") or []
            if len(alternatives) == 1:
                missing.append(self._format_required_alternative_label(alternatives[0]))
            else:
                missing.append(
                    "(" + "|".join(self._format_required_alternative_label(alt) for alt in alternatives) + ")"
                )
        return missing

    def _resolve_required_clauses(self) -> list[dict[str, Any]]:
        clauses: list[dict[str, Any]] = []
        for token in self.require_groups:
            candidates = self._split_require_clause(str(token))
            alternatives: list[dict[str, Any]] = []
            seen_keys: set[tuple[str, str]] = set()
            unresolved: list[str] = []
            for candidate in candidates:
                accession = self._resolve_required_accession(candidate)
                if accession is not None:
                    dedupe_key = ("accession", accession)
                    if dedupe_key in seen_keys:
                        continue
                    seen_keys.add(dedupe_key)
                    alternatives.append(
                        {
                            "kind": "accession",
                            "token": candidate,
                            "accession": accession,
                            "name": accession,
                        }
                    )
                    continue
                taxid: Optional[int]
                if candidate.isdigit():
                    taxid = int(candidate)
                else:
                    taxid = resolve_clade_to_taxid(self.db_manager, candidate)
                if taxid is None:
                    untyped_accession = self._resolve_untyped_required_accession(candidate)
                    if untyped_accession is not None:
                        unresolved.append(f"{candidate} (use acc:{untyped_accession})")
                    else:
                        unresolved.append(candidate)
                    continue
                self.db_manager.cursor.execute(
                    "SELECT name FROM Taxonomy WHERE taxid = ? LIMIT 1",
                    (taxid,),
                )
                row = self.db_manager.cursor.fetchone()
                if row is None:
                    unresolved.append(candidate)
                    continue
                dedupe_key = ("taxid", str(int(taxid)))
                if dedupe_key in seen_keys:
                    continue
                seen_keys.add(dedupe_key)
                alternatives.append(
                    {
                        "kind": "taxid",
                        "token": candidate,
                        "taxid": int(taxid),
                        "name": str(row[0] or taxid),
                    }
                )
            if unresolved:
                raise ValueError(
                    f"Could not resolve required group(s) in '{token}': {', '.join(unresolved)}"
                )
            if not alternatives:
                raise ValueError(f"No resolved required groups for clause '{token}'.")
            clauses.append(
                {
                    "raw": str(token),
                    "alternatives": alternatives,
                }
            )
        return clauses

    def _build_required_membership(
        self,
        accessions: list[str],
        required_taxids: set[int],
        acc_to_tax: dict[str, int],
    ) -> dict[str, set[int]]:
        membership: dict[str, set[int]] = {}
        if not accessions or not required_taxids:
            return membership
        lineage_cache: dict[int, set[int]] = {}
        for accession in accessions:
            taxid = acc_to_tax.get(accession)
            if taxid is None:
                continue
            lineage_taxids = lineage_cache.get(taxid)
            if lineage_taxids is None:
                lineage_rows = self.db_manager.genomes.get_lineage_root_to_leaf(int(taxid)) or []
                lineage_taxids = {int(row[0]) for row in lineage_rows if row and row[0] is not None}
                lineage_cache[taxid] = lineage_taxids
            hit = lineage_taxids.intersection(required_taxids)
            if hit:
                membership[accession] = hit
        return membership

    def _resolve_rank_with_fallback(
        self,
        lineage_rows: list[tuple[int, str, str, Optional[int]]],
        rank_token: str,
    ) -> Optional[str]:
        if not lineage_rows:
            return None
        requested = str(rank_token or "").strip().lower()
        if not requested:
            return None
        exact = None
        for tid, name, rank, _parent in lineage_rows:
            if (rank or "").lower() == requested:
                exact = name or str(tid)
        if exact is not None:
            return exact
        order = {rank: idx for idx, rank in enumerate(RANK_HIERARCHY)}
        if requested not in order:
            return None
        req_idx = order[requested]
        best_idx = -1
        best_name = None
        for tid, name, rank, _parent in lineage_rows:
            token = (rank or "").lower()
            idx = order.get(token)
            if idx is None:
                continue
            if idx <= req_idx and idx > best_idx:
                best_idx = idx
                best_name = name or str(tid)
        return best_name

    def _build_header_taxon_info(
        self,
        lineage_rows: list[tuple[int, str, str, Optional[int]]],
        taxid_val: Optional[int],
        *,
        rank_token: Optional[str] = None,
    ) -> dict[str, str]:
        lineage_map = {
            str(rank).lower(): str(name)
            for tid, name, rank, _parent in lineage_rows
            if rank and name
        }
        taxon_name = next(
            (str(name) for tid, name, _rank, _parent in lineage_rows if tid == taxid_val and name),
            "",
        )
        taxon_rank = next(
            (str(rank).lower() for tid, _name, rank, _parent in lineage_rows if tid == taxid_val and rank),
            "",
        )
        species_full = lineage_map.get("species") or ""
        genus_name = lineage_map.get("genus") or ""
        species_token = ""
        if species_full:
            if genus_name and species_full.lower().startswith(f"{genus_name.lower()} "):
                species_token = species_full[len(genus_name):].strip()
            else:
                parts = species_full.split(maxsplit=1)
                species_token = parts[1] if len(parts) > 1 else parts[0]

        if taxon_rank == "species":
            taxon = species_full or taxon_name
            if not genus_name and taxon:
                genus_name = taxon.split(maxsplit=1)[0]
            if not species_token and taxon:
                parts = taxon.split(maxsplit=1)
                species_token = parts[1] if len(parts) > 1 else parts[0]
        elif taxon_rank == "genus":
            genus_name = taxon_name or genus_name
            taxon = f"{genus_name}_sp" if genus_name else ""
            species_token = "sp" if genus_name else ""
        else:
            fallback_taxon = taxon_name or ""
            taxon = f"{fallback_taxon}_sp" if fallback_taxon else ""
            genus_name = ""
            species_token = ""

        rank_val = self._resolve_rank_with_fallback(lineage_rows, rank_token or "") if rank_token else ""
        return {
            "taxon_name": taxon_name,
            "taxon_rank": taxon_rank,
            "taxon": taxon or "",
            "kingdom": str(lineage_map.get("kingdom") or ""),
            "phylum": str(lineage_map.get("phylum") or ""),
            "class": str(lineage_map.get("class") or ""),
            "order": str(lineage_map.get("order") or ""),
            "family": str(lineage_map.get("family") or ""),
            "genus": genus_name or "",
            "species": species_token or "",
            "rank": str(rank_val or ""),
        }

    def _render_export_header(
        self,
        *,
        accession: str,
        family_id: str,
        source_header: str,
        sequence: str,
        taxon_info: dict[str, str],
        acc_to_tax: dict[str, int],
        gene_cache: dict[str, str],
        bitscore_cache: dict[tuple[str, str], Optional[float]],
        header_template: str,
    ) -> str:
        def _clean(val: Optional[str]) -> str:
            if val is None:
                return ""
            return re.sub(r"\s+", "_", str(val).strip())

        sequence_id = str(source_header or "").lstrip(">").split(maxsplit=1)[0]
        values = {
            "ACCESSION": accession,
            "BITSCORE": self._fmt_header_bitscore(bitscore_cache.get((accession, family_id))),
            "TAXON": _clean(taxon_info.get("taxon") or accession),
            "KINGDOM": _clean(taxon_info.get("kingdom") or ""),
            "PHYLUM": _clean(taxon_info.get("phylum") or ""),
            "CLASS": _clean(taxon_info.get("class") or ""),
            "ORDER": _clean(taxon_info.get("order") or ""),
            "FAMILY": _clean(taxon_info.get("family") or ""),
            "GENUS": _clean(taxon_info.get("genus") or ""),
            "SPECIES": _clean(taxon_info.get("species") or ""),
            "RANK": _clean(taxon_info.get("rank") or ""),
            "BUSCO": _clean(family_id),
            "SEQUENCE": _clean(sequence_id),
            "LENGTH": str(len(sequence)),
            "GENE": _clean(gene_cache.get(family_id) or ""),
            "TAXID": str(acc_to_tax.get(accession, "")),
        }
        token_re = re.compile("|".join(sorted(values.keys(), key=len, reverse=True)))
        rendered = token_re.sub(lambda m: values.get(m.group(0), ""), header_template)
        return f">{rendered}"

    def _header_needs_copy_suffix(self, *, header_template: str) -> bool:
        template = str(header_template or "")
        return "GENE" not in template and "SEQUENCE" not in template

    def _effective_header_rank(self, header_template: str) -> Optional[str]:
        raw_rank = str(self.header_rank).strip().lower() if self.header_rank else ""
        if raw_rank:
            return raw_rank
        if "RANK" in str(header_template or ""):
            return "phylum"
        return None

    def _fetch_assembly_info(self, accessions: list[str]) -> dict[str, tuple[Optional[int], str]]:
        info: dict[str, tuple[Optional[int], str]] = {}
        if not accessions:
            return info
        for chunk in [accessions[i : i + 900] for i in range(0, len(accessions), 900)]:
            placeholders = ",".join("?" for _ in chunk)
            self.db_manager.cursor.execute(
                f"""
                SELECT taxid, name, assembly_accession
                FROM TaxonomyAssemblySummary
                WHERE assembly_accession IN ({placeholders})
                """,
                tuple(chunk),
            )
            for taxid, name, accession in self.db_manager.cursor.fetchall() or []:
                info[str(accession)] = (int(taxid) if taxid is not None else None, str(name or ""))
        return info

    def _fetch_busco_scores(
        self,
        accessions: list[str],
        library_id: int,
        *,
        run_refs: Optional[list[tuple[str, int]]] = None,
        decont_run_id: Optional[str],
        paralog_run_id: Optional[str] = None,
        include_decontam_override: Optional[bool] = None,
    ) -> dict[str, tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]]:
        scores: dict[str, tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]] = {}
        if not accessions:
            return scores
        include_decontam = (
            self.include_decontamination_in_score
            if include_decontam_override is None
            else include_decontam_override
        )
        if run_refs:
            ref_rows = self.db_manager.busco.get_display_results_for_runs(
                library_id=library_id,
                run_refs=run_refs,
                include_paralog=self.include_paralog_filtering_in_score,
                paralog_run_id=paralog_run_id,
                include_decontam=include_decontam,
                decont_run_id=decont_run_id,
                allow_ambiguous_contaminants=self.allow_ambiguous_contaminants,
                strict_decontamination=self.strict_decontamination,
                rescue_duplicates=self.rescue_duplicates,
            ) or {}
            for (acc, _run_id), row in ref_rows.items():
                scores[str(acc)] = (
                    row[0],
                    row[1],
                    row[2],
                    row[3],
                    row[4],
                    row[5],
                    row[6],
                )
            return scores

        for chunk in [accessions[i : i + 900] for i in range(0, len(accessions), 900)]:
            rows = self.db_manager.busco.get_results_adjusted(
                library_id=library_id,
                accessions=chunk,
                include_paralog=self.include_paralog_filtering_in_score,
                paralog_run_id=paralog_run_id,
                include_decontam=include_decontam,
                decont_run_id=decont_run_id,
                allow_ambiguous_contaminants=self.allow_ambiguous_contaminants,
                strict_decontamination=self.strict_decontamination,
                rescue_duplicates=self.rescue_duplicates,
            )
            for row in rows or []:
                acc = str(row[0])
                scores[acc] = (
                    row[3],
                    row[4],
                    row[5],
                    row[6],
                    row[7],
                    row[8],
                    row[9],
                )
        return scores

    def _fetch_family_bitscores(
        self,
        accessions: list[str],
        family_ids: list[str],
        library_id: int,
    ) -> dict[tuple[str, str], Optional[float]]:
        scores: dict[tuple[str, str], Optional[float]] = {}
        if not accessions or not family_ids:
            return scores
        for acc_chunk in [accessions[i : i + 900] for i in range(0, len(accessions), 900)]:
            acc_placeholders = ",".join("?" for _ in acc_chunk)
            for fam_chunk in [family_ids[i : i + 900] for i in range(0, len(family_ids), 900)]:
                fam_placeholders = ",".join("?" for _ in fam_chunk)
                self.db_manager.cursor.execute(
                    f"""
                    SELECT accession, family_id, score
                    FROM BUSCO_Family_Data
                    WHERE library_id = ?
                      AND accession IN ({acc_placeholders})
                      AND family_id IN ({fam_placeholders})
                    """,
                    (library_id, *acc_chunk, *fam_chunk),
                )
                for accession, family_id, score in self.db_manager.cursor.fetchall() or []:
                    scores[(str(accession), str(family_id))] = float(score) if score is not None else None
        return scores

    def _fetch_decontamination_breakdown(
        self,
        accessions: list[str],
        library_id: int,
        *,
        run_id: Optional[str],
    ) -> dict[str, tuple[Optional[float], Optional[float], Optional[float]]]:
        if not accessions:
            return {}
        return self.db_manager.filtering.get_decontamination_decision_percentages(
            library_id=library_id,
            accessions=accessions,
            run_id=run_id,
        )

    def _fetch_decontamination_summaries(
        self,
        accessions: list[str],
        library_id: int,
        *,
        run_id: Optional[str],
    ) -> dict[str, tuple[Optional[str], Optional[str], Optional[str]]]:
        if not accessions:
            return {}
        parent_id = self.db_manager.libraries.get_parent_id(library_id)
        return self.db_manager.filtering.get_latest_decontamination_summary_with_fallback(
            target_library_id=library_id,
            parent_library_id=parent_id,
            accessions=accessions,
            run_id=run_id,
        )

    def _write_lineage_csv(self, accessions: list[str]) -> bool:
        if not accessions:
            return False
        output_path = self.lineage_csv_path or os.path.join(self.report_dir, "lineage.csv")
        try:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
        except OSError as exc:
            self.log(f"Failed to create lineage CSV directory for {output_path}: {exc}", "WARNING")

        placeholders = ",".join("?" for _ in accessions)
        try:
            self.db_manager.cursor.execute(
                f"SELECT accession, taxid FROM Genome WHERE accession IN ({placeholders})",
                tuple(accessions),
            )
            acc_tax_rows = self.db_manager.cursor.fetchall() or []
        except Exception as exc:  # boundary: optional lineage CSV database lookup failure disables that report.
            self.log(f"Failed to resolve taxids for lineage CSV: {exc}", "WARNING")
            return False

        acc_to_tax: dict[str, int] = {}
        for acc, tax in acc_tax_rows:
            if acc is None or tax is None:
                continue
            acc_to_tax[str(acc)] = int(tax)

        rank_columns = [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
            "superkingdom",
        ]
        records: list[dict[str, str]] = []
        lineage_cache: dict[int, list[tuple[int, str, str, Optional[int]]]] = {}
        for acc in accessions:
            taxid = acc_to_tax.get(acc)
            lineage_rows = []
            if taxid is not None:
                lineage_rows = lineage_cache.get(taxid)
                if lineage_rows is None:
                    lineage_rows = self.db_manager.genomes.get_lineage_root_to_leaf(taxid) or []
                    lineage_cache[taxid] = lineage_rows
            lineage_map = {str(rank).lower(): name for (tid, name, rank, _parent) in lineage_rows if rank and name}
            taxon_name = next((name for (tid, name, _rank, _parent) in lineage_rows if tid == taxid), "") if lineage_rows else ""
            taxon_rank = next((rank for (tid, _name, rank, _parent) in lineage_rows if tid == taxid and rank), "") if lineage_rows else ""

            row: dict[str, str] = {
                "accession": acc,
                "taxid": str(taxid or ""),
                "taxon_name": taxon_name or "",
                "taxon_rank": str(taxon_rank or ""),
            }
            for col in rank_columns:
                row[col] = lineage_map.get(col, "")
            records.append(row)

        if not records:
            return False

        records.sort(key=lambda r: r["accession"])
        fieldnames = ["accession", "taxid", "taxon_name", "taxon_rank", *rank_columns]
        try:
            with open(output_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(records)
            self.log(f"Wrote lineage CSV to {output_path}", "INFO")
            return True
        except (OSError, csv.Error) as exc:
            self.log(f"Failed to write lineage CSV: {exc}", "WARNING")
            return False

    def _write_busco_report(
        self,
        accessions: list[str],
        *,
        busco_library_id: int,
        decont_run_id: Optional[str],
        custom_library: bool,
    ) -> bool:
        if not accessions:
            return False
        report_path = self.busco_report_path or os.path.join(self.report_dir, "busco_report.tsv")
        try:
            os.makedirs(os.path.dirname(report_path), exist_ok=True)
        except OSError as exc:
            self.log(f"Failed to create BUSCO report directory for {report_path}: {exc}", "WARNING")

        info_map = self._fetch_assembly_info(accessions)
        contam_map = None
        busco_override = None
        run_refs = []
        if hasattr(self, "selected_busco_runs"):
            run_refs = [
                (acc, int(meta["run_id"]))
                for acc, meta in getattr(self, "selected_busco_runs", {}).items()
                if meta.get("run_id") is not None
            ]
        if self.disable_decont_filter and self.include_decontamination_in_score is False:
            busco_override = False
            contam_map = self._fetch_busco_scores(
                accessions,
                busco_library_id,
                run_refs=run_refs,
                decont_run_id=decont_run_id,
                paralog_run_id=self.paralog_run_id,
                include_decontam_override=True,
            )
        busco_map = self._fetch_busco_scores(
            accessions,
            busco_library_id,
            run_refs=run_refs,
            decont_run_id=decont_run_id,
            paralog_run_id=self.paralog_run_id,
            include_decontam_override=busco_override,
        )
        use_decont_scores = self.include_decontamination_in_score is not False
        decont_summary = self._fetch_decontamination_summaries(accessions, busco_library_id, run_id=decont_run_id)
        decont_breakdown = {}
        if self.busco_report_extended:
            decont_breakdown = self._fetch_decontamination_breakdown(
                accessions,
                busco_library_id,
                run_id=decont_run_id,
            )

        headers = ["accession", "species", "busco_run_id", "pipeline", "format", "proteome_profile"]
        if custom_library:
            headers.extend(
                [
                    "complete",
                    "single_copy_complete",
                    "duplicated",
                    "fragmented",
                    "hidden_paralog",
                    "contaminated",
                ]
            )
        else:
            headers.extend(
                [
                    "complete",
                    "single_copy_complete",
                    "duplicated",
                    "fragmented",
                    "contaminated",
                ]
            )
        if self.busco_report_extended:
            headers.extend(["support", "weak", "unknown", "decontamination_run"])
        headers.append("missing")
        if custom_library:
            headers.append("contaminated_assembly")

        try:
            with open(report_path, "w", encoding="utf-8") as handle:
                handle.write("\t".join(headers) + "\n")
                for acc in accessions:
                    tax_info = info_map.get(acc, (None, ""))
                    run_meta = (getattr(self, "selected_busco_runs", {}) or {}).get(acc, {})
                    run_id = run_meta.get("run_id")
                    pipeline = run_meta.get("pipeline")
                    input_mode = str(run_meta.get("input_mode") or "").strip().lower()
                    format_label = "proteome" if input_mode == "protein" else ("genome" if input_mode == "genome" else "NA")
                    profile_name = run_meta.get("proteome_profile")
                    row = [
                        acc,
                        tax_info[1],
                        str(run_id) if run_id is not None else "NA",
                        str(pipeline) if pipeline else "NA",
                        format_label,
                        str(profile_name) if profile_name else ("unset" if input_mode == "protein" else "NA"),
                    ]
                    bvals = busco_map.get(acc)
                    if bvals:
                        complete, single_copy, duplicated, fragmented, missing, hidden_paralog, contaminated = bvals
                        if contam_map:
                            cbvals = contam_map.get(acc)
                            if cbvals:
                                contaminated = cbvals[6]
                        contam_out = "SKIPPED" if not use_decont_scores else self._fmt_busco_value(contaminated)
                        if custom_library:
                            row.extend(
                                self._fmt_busco_value(val)
                                for val in (
                                    complete,
                                    single_copy,
                                    duplicated,
                                    fragmented,
                                    hidden_paralog,
                                )
                            )
                            row.append(contam_out)
                        else:
                            row.extend(
                                self._fmt_busco_value(val)
                                for val in (
                                    complete,
                                    single_copy,
                                    duplicated,
                                    fragmented,
                                )
                            )
                            row.append(contam_out)
                        if self.busco_report_extended:
                            support, weak, unknown = decont_breakdown.get(acc, (None, None, None))
                            row.extend(self._fmt_busco_value(val) for val in (support, weak, unknown))
                            decont_id = decont_summary.get(acc, (None, None, None))[0]
                            row.append(str(decont_id) if decont_id else "NA")
                        row.append(self._fmt_busco_value(missing))
                        if custom_library:
                            decision = decont_summary.get(acc, (None, None, None))[1]
                            row.append(str(decision) if decision else "NA")
                    else:
                        if custom_library:
                            row.extend(["NA", "NA", "NA", "NA", "NA", "NA"])
                        else:
                            row.extend(["NA", "NA", "NA", "NA", "NA"])
                        if self.busco_report_extended:
                            row.extend(["NA", "NA", "NA", "NA"])
                        row.append("NA")
                        if custom_library:
                            row.append("NA")
                    handle.write("\t".join(row) + "\n")
            self.log(f"Wrote BUSCO report to {report_path}", "INFO")
            return True
        except (OSError, csv.Error) as exc:
            self.log(f"Failed to write BUSCO report: {exc}", "WARNING")
            return False

    def _write_busco_family_matrix(
        self,
        accessions: list[str],
        *,
        library_id: int,
        busco_library_id: int,
        decont_run_id: Optional[str],
        custom_library: bool,
    ) -> bool:
        if not accessions:
            return False
        report_path = self.busco_family_matrix_path
        if not report_path:
            if self.busco_report_path:
                base, ext = os.path.splitext(self.busco_report_path)
                report_path = f"{base}_family_matrix.tsv" if ext else f"{self.busco_report_path}_family_matrix.tsv"
            else:
                report_path = os.path.join(self.report_dir, "busco_family_matrix.tsv")
        try:
            os.makedirs(os.path.dirname(report_path), exist_ok=True)
        except OSError as exc:
            self.log(f"Failed to create BUSCO matrix report directory for {report_path}: {exc}", "WARNING")

        info_map = self._fetch_assembly_info(accessions)

        # Families to include (subset for custom libraries, full for default)
        families: list[str] = []
        try:
            self.db_manager.cursor.execute(
                "SELECT family_id FROM BUSCO_descriptions WHERE library_id = ? ORDER BY family_id",
                (library_id,),
            )
            families = [str(r[0]) for r in self.db_manager.cursor.fetchall() or []]
        except Exception as exc:  # boundary: optional family list lookup failure yields no matrix.
            self.log(f"Failed to load BUSCO families for matrix: {exc}", "WARNING")
            families = []
        if not families:
            self.log("No BUSCO families found for family matrix report.", "WARNING")
            return False

        # Base BUSCO family statuses
        status_map: dict[str, dict[str, int]] = {acc: {} for acc in accessions}
        try:
            for chunk in [accessions[i : i + 900] for i in range(0, len(accessions), 900)]:
                placeholders = ",".join("?" for _ in chunk)
                self.db_manager.cursor.execute(
                    f"""
                    SELECT bfd.accession, bfd.family_id, bfd.status
                    FROM BUSCO_Family_Data bfd
                    JOIN BUSCO_descriptions bd
                      ON bd.family_id = bfd.family_id AND bd.library_id = ?
                    WHERE bfd.library_id = ?
                      AND bfd.accession IN ({placeholders})
                    """,
                    (library_id, busco_library_id, *chunk),
                )
                for acc, fam, status in self.db_manager.cursor.fetchall() or []:
                    status_map.setdefault(str(acc), {})[str(fam)] = int(status or 0)
        except Exception as exc:  # boundary: optional BUSCO status matrix annotation.
            self.log(f"Failed to load BUSCO family statuses for matrix: {exc}", "WARNING")

        # Paralog filtering (hidden paralogs)
        hidden_flags: dict[str, set[str]] = {}
        if custom_library:
            try:
                rows = self.db_manager.filtering.get_paralog_results(target_library_id=library_id)
                for fam, lib_id, target_lib, acc, clean in rows:
                    if target_lib != library_id or lib_id != busco_library_id or acc not in accessions:
                        continue
                    if clean is None:
                        continue
                    if int(clean) == 0:
                        hidden_flags.setdefault(str(acc), set()).add(str(fam))
            except Exception as exc:  # boundary: optional hidden-paralog matrix annotation.
                self.log(f"Failed to load hidden paralog flags for matrix: {exc}", "WARNING")

        # Decontamination votes per family
        decont_flags: dict[str, set[str]] = {}
        summaries: dict[str, tuple[str, Optional[str], Optional[str], int]] = {}
        try:
            primary = self.db_manager.filtering.latest_decont_summary(
                target_library_id=library_id,
                accessions=accessions,
                run_id=decont_run_id,
            )
            summaries = {acc: (rid, decision, date, library_id) for acc, (rid, decision, date) in primary.items()}
            parent_id = self.db_manager.libraries.get_parent_id(library_id)
            if parent_id:
                fallback = self.db_manager.filtering.latest_decont_summary(
                    target_library_id=parent_id,
                    accessions=accessions,
                    run_id=decont_run_id,
                )
                for acc, (rid, decision, date) in fallback.items():
                    if acc not in summaries:
                        summaries[acc] = (rid, decision, date, parent_id)
        except Exception as exc:  # boundary: optional decontamination summary lookup failure disables matrix annotation.
            self.log(f"Failed to load decontamination summaries for matrix: {exc}", "WARNING")
            summaries = {}

        allow_ambiguous = bool(self.allow_ambiguous_contaminants)
        strict_decontam = bool(self.strict_decontamination)
        supported_decisions = {"support"}
        if not strict_decontam:
            supported_decisions.add("weak")
        if allow_ambiguous:
            supported_decisions.add("unknown")

        if summaries:
            accessions_by_run: dict[tuple[str, int], list[str]] = {}
            for acc, (run_id, _decision, _date, source_lib) in summaries.items():
                if run_id is None:
                    continue
                accessions_by_run.setdefault((str(run_id), int(source_lib)), []).append(acc)

            for (run_id, target_lib), acc_list in accessions_by_run.items():
                for chunk in [acc_list[i : i + 900] for i in range(0, len(acc_list), 900)]:
                    placeholders = ",".join("?" for _ in chunk)
                    self.db_manager.cursor.execute(
                        f"""
                        SELECT v.accession, v.family_id, COALESCE(v.decision, 'unknown') AS decision
                        FROM Decontamination_Busco_Votes v
                        JOIN BUSCO_descriptions bd
                          ON bd.family_id = v.family_id AND bd.library_id = ?
                        WHERE v.target_library_id = ?
                          AND v.busco_library_id = ?
                          AND v.run_id = ?
                          AND v.accession IN ({placeholders})
                        """,
                        (library_id, target_lib, busco_library_id, run_id, *chunk),
                    )
                    for acc, fam, decision in self.db_manager.cursor.fetchall() or []:
                        if str(decision).lower() not in supported_decisions:
                            decont_flags.setdefault(str(acc), set()).add(str(fam))
        # Write matrix
        try:
            with open(report_path, "w", encoding="utf-8") as handle:
                headers = ["accession", "species", *families]
                handle.write("\t".join(headers) + "\n")
                status_labels = {1: "single_copy", 2: "duplicated", 3: "fragmented", 4: "missing"}
                for acc in accessions:
                    tax_info = info_map.get(acc, (None, ""))
                    row = [acc, tax_info[1]]
                    acc_status = status_map.get(acc, {})
                    acc_hidden = hidden_flags.get(acc, set())
                    acc_decont = decont_flags.get(acc, set())
                    for fam in families:
                        base = status_labels.get(acc_status.get(fam, 4), "missing")
                        flags = []
                        if fam in acc_hidden:
                            flags.append("hidden_paralog")
                        if fam in acc_decont:
                            flags.append("contaminated")
                        if flags:
                            base = base + "|" + "|".join(flags)
                        row.append(base)
                    handle.write("\t".join(row) + "\n")
            self.log(f"Wrote BUSCO family matrix to {report_path}", "INFO")
            return True
        except (OSError, csv.Error) as exc:
            self.log(f"Failed to write BUSCO family matrix: {exc}", "WARNING")
            return False

    def _write_taxa_occupancy_report(
        self,
        selected: list[str],
        *,
        allowed_accessions: set[str],
        accession_counts: dict[str, int],
        total_families_seen: int,
    ) -> bool:
        if not selected:
            return False
        report_path = self._taxa_occupancy_report_path()
        try:
            os.makedirs(os.path.dirname(report_path), exist_ok=True)
        except OSError as exc:
            self.log(f"Failed to create taxa occupancy report directory for {report_path}: {exc}", "WARNING")

        info_map = self._fetch_assembly_info(selected)
        threshold_count = (
            self.min_taxa_occupancy * total_families_seen
            if total_families_seen > 0 and self.min_taxa_occupancy > 0
            else 0.0
        )
        rows = []
        for accession in sorted(selected):
            taxid, species = info_map.get(accession, (None, ""))
            present = int(accession_counts.get(accession, 0))
            absent = max(0, int(total_families_seen) - present)
            pct_present = (present / total_families_seen) if total_families_seen > 0 else 0.0
            retained = accession in allowed_accessions
            rows.append(
                {
                    "accession": accession,
                    "species": str(species or ""),
                    "taxid": str(taxid or ""),
                    "families_evaluated": str(int(total_families_seen)),
                    "buscos_present": str(present),
                    "buscos_absent": str(absent),
                    "pct_present": f"{pct_present:.6f}",
                    "min_taxa_occupancy": f"{self.min_taxa_occupancy:.6f}",
                    "min_taxa_occupancy_threshold_count": f"{threshold_count:.6f}",
                    "retained_after_taxa_occupancy": "yes" if retained else "no",
                    "removed_by_min_taxa_occupancy": "no" if retained else "yes",
                }
            )

        headers = [
            "accession",
            "species",
            "taxid",
            "families_evaluated",
            "buscos_present",
            "buscos_absent",
            "pct_present",
            "min_taxa_occupancy",
            "min_taxa_occupancy_threshold_count",
            "retained_after_taxa_occupancy",
            "removed_by_min_taxa_occupancy",
        ]
        try:
            with open(report_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            self.log(f"Wrote taxa occupancy report to {report_path}", "INFO")
            return True
        except (OSError, csv.Error) as exc:
            self.log(f"Failed to write taxa occupancy report: {exc}", "WARNING")
            return False

    def _write_export_parameters(
        self,
        *,
        selected: list[str],
        busco_library_id: int,
        busco_run_selection: dict[str, dict[str, Any]],
        decont_selection_run: Optional[str],
        resolved_decont_run: Optional[str],
        report_decont_run: Optional[str],
        custom_library: bool,
    ) -> bool:
        if not self.out_dir:
            return False
        params_path = None
        if self.busco_report_path:
            base, ext = os.path.splitext(self.busco_report_path)
            params_path = f"{base}_params.txt" if ext else f"{self.busco_report_path}_params.txt"
        if not params_path:
            params_path = os.path.join(self.report_dir, "export_parameters.txt")
        try:
            os.makedirs(os.path.dirname(params_path), exist_ok=True)
        except OSError as exc:
            self.log(f"Failed to create export parameters directory for {params_path}: {exc}", "WARNING")

        params = {
            "accession": self.accession,
            "accessions_count": len((self.data.get("accessions") or [])),
            "taxid": self.taxid,
            "clade": self.clade,
            "quantity": self.rule_quantity,
            "rank": self.rule_rank,
            "library_id": self.library_id,
            "library_name": self.library_name,
            "sequence_type": self.sequence_type,
            "busco_library_id": busco_library_id,
            "requested_proteome_profile": self.proteome_profile,
            "requested_proteome_profiles": self.proteome_profiles,
            "preferred_proteome_profile": self.prefer_proteome_profile,
            "busco_run_ids": self.data.get("busco_run_ids"),
            "busco_pipeline": self.data.get("busco_pipeline"),
            "busco_input_mode": self.data.get("busco_input_mode"),
            "prefer_busco_pipeline": self.data.get("prefer_busco_pipeline"),
            "prefer_busco_input_mode": self.data.get("prefer_busco_input_mode"),
            "busco_run_selection_policy": self.data.get("busco_run_selection") or "primary",
            "custom_library": custom_library,
            "out_dir": str(self.out_dir),
            "protein_only": self.protein_only,
            "min_completeness": self.min_completeness,
            "min_single_copy_complete": self.min_single_copy_complete,
            "min_occupancy": self.min_occupancy,
            "min_taxa_occupancy": self.min_taxa_occupancy,
            "require": self.require_groups,
            "require_resolved": self.resolved_require_clauses,
            "disable_paralog_filter": self.disable_paralog_filter,
            "disable_decont_filter": self.disable_decont_filter,
            "require_paralog_filtering": self.require_paralog_filtering,
            "require_decontamination": self.require_decontamination,
            "include_paralog_filtering_in_score": self.include_paralog_filtering_in_score,
            "paralog_run_id": self.paralog_run_id,
            "include_decontamination_in_score": self.include_decontamination_in_score,
            "allow_ambiguous_contaminants": self.allow_ambiguous_contaminants,
            "strict_decontamination": self.strict_decontamination,
            "decont_run_id": self.decont_run_id,
            "decont_run_label": self.decont_run_label,
            "decontamination_run": self.decontamination_run,
            "decont_selection_run": decont_selection_run,
            "resolved_decont_run": resolved_decont_run,
            "report_decont_run": report_decont_run,
            "write_lineage_csv": self.write_lineage_csv,
            "lineage_csv_path": str(self.lineage_csv_path) if self.lineage_csv_path else None,
            "write_busco_report": self.write_busco_report,
            "busco_report_path": str(self.busco_report_path) if self.busco_report_path else None,
            "busco_report_extended": self.busco_report_extended,
            "rescue_duplicates": self.rescue_duplicates,
            "write_busco_family_matrix": self.write_busco_family_matrix,
            "busco_family_matrix_path": str(self.busco_family_matrix_path) if self.busco_family_matrix_path else None,
            "taxa_occupancy_report_path": self._taxa_occupancy_report_path() if self.report_dir else None,
            "export_task_log_path": str(self._export_log_path) if self._export_log_path else None,
            "retain_headers": self.retain_headers,
            "header_template": self.header_template,
            "header_rank": self.header_rank,
            "family_ids": self.family_ids,
            "include_duplicated": self.include_duplicated,
            "rerun": self.rerun,
            "selected_accessions": len(selected),
            "selected_busco_runs": busco_run_selection,
        }
        try:
            with open(params_path, "w", encoding="utf-8") as handle:
                for key in sorted(params.keys()):
                    handle.write(f"{key}\t{params[key]}\n")
            self.log(f"Wrote export parameters to {params_path}", "INFO")
            return True
        except (OSError, TypeError, ValueError) as exc:
            self.log(f"Failed to write export parameters: {exc}", "WARNING")
            return False

    def run(self):
        # Basic validation
        if not self._ensure_export_paths():
            return self.handle_exception(
                "Output directory is not specified and no exports root is configured.",
                {"out_dir": self.out_dir},
            )
        if self._enable_export_log_copy() and self._export_log_path:
            self.log(f"Writing export task log copy to {self._export_log_path}", "INFO")

        # Resolve library info from DB
        if not self.library_id and self.library_name:
            self.library_id = self.db_manager.libraries.get_id(self.library_name)
        if not self.library_id:
            return self.handle_exception("library_id (or a resolvable library name) is required for BUSCO context.", {"library_id": self.library_id, "library_name": self.library_name})

        libs = self.db_manager.libraries.get(self.library_id)
        if not libs:
            return self.handle_exception("Library not found in database.", {"library_id": self.library_id})
        # Libraries schema: (library_id, library_name, odb_version, taxid, size, location, parent_id)
        lib_row = libs[0]
        lib_name = lib_row[1]
        lib_location = lib_row[5]
        lib_parent_id = lib_row[6]
        # If custom library (has parent), BUSCO should be run against parent lineage
        busco_lib_id = lib_parent_id if lib_parent_id else self.library_id
        parent_name = None
        if lib_parent_id:
            parent = self.db_manager.libraries.get(lib_parent_id)
            if not parent:
                return self.handle_exception("Parent library not found in database.", {"library_id": self.library_id, "parent_id": lib_parent_id})
            parent_name = parent[0][1]
        busco_lineage_name = parent_name if parent_name else lib_name
        self.log(
            f"Export context: library_id={self.library_id} library_name={lib_name} "
            f"parent_id={lib_parent_id} busco_lineage={busco_lineage_name}",
            "INFO",
        )

        additional_accessions = [self.accession] if self.accession else None
        if self.taxid is None and self.clade:
            try:
                self.taxid = resolve_clade_to_taxid(self.db_manager, self.clade)
            except (LookupError, ValueError):
                self.taxid = None
            if self.taxid is None:
                return self.handle_exception(
                    "Unknown clade; could not resolve to taxid for selection.",
                    {"clade": self.clade},
                )

        if self.rule_rank and self.taxid is None and not (additional_accessions or self.selector_accessions()):
            return self.handle_exception(
                "Rule-based selection by rank requires a taxid to partition within.",
                {"quantity": self.rule_quantity, "rank": self.rule_rank},
            )

        decont_selection_run = None
        report_decont_run_id = None
        if self.decont_run_id or self.decont_run_label or self.decontamination_run:
            decont_selection_run = self.decont_run_id or self.decontamination_run
            if decont_selection_run is None and self.decont_run_label:
                rows = self._list_decont_runs(self.library_id, self.decont_run_label)
                if rows:
                    decont_selection_run = sorted(rows, key=lambda r: r[3] or "", reverse=True)[0][0]
            report_decont_run_id = decont_selection_run

        explicit_selector_accessions = list(dict.fromkeys(normalize_accessions(self.data.get("accessions") or [])))
        if self.accession:
            explicit_selector_accessions.append(str(self.accession))
        explicit_selector_accessions = list(dict.fromkeys(normalize_accessions(explicit_selector_accessions)))

        if explicit_selector_accessions and self.taxid is None and not self.clade and self.rule_quantity is None and self.rule_rank is None and not self.data.get("ranks") and not self.data.get("quantities"):
            selected = explicit_selector_accessions
        else:
            try:
                selected = self.prepare_selectors(
                    taxid=self.taxid,
                    additional=additional_accessions or explicit_selector_accessions,
                    rule_quantity=self.rule_quantity,
                    rule_rank=self.rule_rank,
                    busco_library_id=self.library_id,
                    downloaded_only=False,
                    protein_only=self.protein_only,
                    status_min=1,
                    require_candidates=True,
                    use_rule_selection=True,
                    use_busco=True,
                    min_completeness=self.min_completeness,
                    min_single_copy_complete=self.min_single_copy_complete,
                    include_paralog_filtering_in_score=self.include_paralog_filtering_in_score,
                    paralog_run_id=self.paralog_run_id,
                    include_decontamination_in_score=self.include_decontamination_in_score,
                    decontamination_run_id=decont_selection_run,
                    allow_ambiguous_contaminants=self.allow_ambiguous_contaminants,
                    strict_decontamination=self.strict_decontamination,
                    rescue_duplicates=self.rescue_duplicates,
                    paralog_filtered=bool(self.require_paralog_filtering and not self.disable_paralog_filter and lib_parent_id is not None),
                    decontaminated=bool(self.require_decontamination and not self.disable_decont_filter),
                    decontamination_run=decont_selection_run,
                    ignore_contaminated_assemblies=False if self.disable_decont_filter else None,
                )
            except ValueError:
                return self.handle_exception(
                    "No downloaded accessions matched the provided accession(s)/taxid.",
                    {
                        "accession": self.accession,
                        "taxid": self.taxid,
                        "selectors": self.selector_accessions(),
                    },
                )

        if not selected:
            return self.handle_exception(
                "Rule-based selection yielded no accessions.",
                {"quantity": self.rule_quantity, "rank": self.rule_rank, "taxid": self.taxid},
            )
        self.log(
            f"Selector results: selected={len(selected)} taxid={self.taxid} rank={self.rule_rank} quantity={self.rule_quantity}",
            "INFO",
        )
        if len(selected) <= 25:
            self.log(f"Selected accessions: {','.join(selected)}", "DEBUG")

        required_clauses: list[dict[str, Any]] = []
        required_taxids: set[int] = set()
        if self.require_groups:
            try:
                required_clauses = self._resolve_required_clauses()
            except ValueError as exc:
                return self.handle_exception(
                    "Failed to resolve --require groups.",
                    {"require": self.require_groups, "error": str(exc)},
                )
            self.resolved_require_clauses = list(required_clauses)
            required_taxids = {
                int(alt["taxid"])
                for clause in required_clauses
                for alt in (clause.get("alternatives") or [])
                if str(alt.get("kind") or "").lower() != "accession"
            }
            clause_labels: list[str] = []
            for clause in required_clauses:
                alternatives = clause.get("alternatives") or []
                if len(alternatives) == 1:
                    alt = alternatives[0]
                    clause_labels.append(self._format_required_alternative_label(alt))
                else:
                    clause_labels.append(
                        "(" + "|".join(self._format_required_alternative_label(alt) for alt in alternatives) + ")"
                    )
            self.log(
                "Required-group filter enabled: "
                + ",".join(clause_labels),
                "INFO",
            )

        decont_run_id = None if self.disable_decont_filter else self._resolve_decont_run(
            self.library_id,
            selected,
            parent_id=lib_parent_id,
        )
        if self.disable_decont_filter:
            self.log("Decontamination filtering disabled by request.", "INFO")
        elif self.require_decontamination and decont_run_id is None:
            return self.handle_exception(
                "Decontamination filtering required but no usable decontamination run was found.",
                {"library_id": self.library_id, "selected": len(selected)},
            )

        # Require BUSCO results for selected accessions (parent lineage if custom)
        present = set(self.db_manager.busco.get_processed_accessions(busco_lib_id))
        missing = [acc for acc in selected if acc not in present]
        if missing:
            return self.handle_exception(
                "Missing BUSCO results for selected accessions. Run BUSCO first.",
                {"missing": missing[:25], "missing_count": len(missing), "busco_lineage": busco_lineage_name},
            )

        missing_genomes = []
        missing_requested_profile = []
        selected_busco_runs: dict[str, dict[str, Any]] = {}
        requested_profile = self.proteome_profile or (self.proteome_profiles[0] if self.proteome_profiles else None)
        explicit_profile_requested = bool(str(requested_profile or "").strip())
        run_ids_filter = self.data.get("busco_run_ids")
        pipeline_filter = self.data.get("busco_pipeline")
        input_mode_filter = self.data.get("busco_input_mode")
        preferred_pipeline_filter = self.data.get("prefer_busco_pipeline")
        preferred_input_mode_filter = self.data.get("prefer_busco_input_mode")
        selection_policy = self.data.get("busco_run_selection") or "primary"
        for acc in selected:
            gpath = self.db_manager.genomes.resolve_path(acc)
            if not gpath:
                missing_genomes.append(acc)
                continue
            primary_purpose = "export_nucleotide" if self.sequence_type == "nucleotide" else "export_protein"
            run_filters = {
                "run_ids": run_ids_filter,
                "pipeline": pipeline_filter,
                "input_mode": input_mode_filter,
                "preferred_pipeline": preferred_pipeline_filter,
                "preferred_input_mode": preferred_input_mode_filter,
                "proteome_profile": requested_profile,
                "preferred_proteome_profile": self.prefer_proteome_profile,
                "selection": selection_policy,
            }
            duplicate_candidates = self.db_manager.busco.get_runs_for_primary_choice(
                acc,
                busco_lib_id,
                run_ids=run_ids_filter,
                pipeline=pipeline_filter,
                input_mode=input_mode_filter,
                proteome_profile=requested_profile,
            ) or []
            if len(duplicate_candidates) > 1:
                dup_ids = [str(row[0]) for row in duplicate_candidates[:20] if row and row[0] is not None]
                suffix = "" if len(duplicate_candidates) <= 20 else ",..."
                self.log(
                    f"Multiple BUSCO runs matched export filters for {acc}; selecting one by policy '{selection_policy}'. "
                    f"Candidates ({len(duplicate_candidates)}): {','.join(dup_ids)}{suffix}",
                    "WARNING",
                )
            run_id = self.db_manager.busco.get_effective_run_id_for_accession(
                acc,
                busco_lib_id,
                purpose=primary_purpose,
                **run_filters,
            )
            resolved_purpose = primary_purpose
            if run_id is None and self.sequence_type != "nucleotide":
                run_id = self.db_manager.busco.get_effective_run_id_for_accession(
                    acc,
                    busco_lib_id,
                    purpose="default",
                    **run_filters,
                )
                if run_id is not None:
                    resolved_purpose = "default"

            chosen_run = self.db_manager.busco.get_run(int(run_id)) if run_id is not None else None
            if explicit_profile_requested and run_id is None:
                missing_requested_profile.append(acc)
                continue

            if run_id is not None:
                selected_busco_runs[acc] = {
                    "run_id": run_id,
                    "purpose": resolved_purpose,
                    "pipeline": str(chosen_run[5]) if chosen_run and len(chosen_run) > 5 else None,
                    "input_mode": str(chosen_run[4]) if chosen_run and len(chosen_run) > 4 else None,
                    "proteome_profile": str(chosen_run[15]) if chosen_run and len(chosen_run) > 15 and chosen_run[15] is not None else None,
                    "result_dir": str(chosen_run[7]) if chosen_run and len(chosen_run) > 7 and chosen_run[7] is not None else None,
                    "selection_source": "run_query",
                }
            else:
                selected_busco_runs[acc] = {
                    "run_id": None,
                    "purpose": primary_purpose,
                    "pipeline": None,
                    "input_mode": None,
                    "proteome_profile": None,
                    "result_dir": None,
                    "selection_source": "none",
                }
        if missing_genomes:
            self.log(
                f"Missing genome paths for {len(missing_genomes)} accessions (showing up to 10): "
                f"{','.join(missing_genomes[:10])}",
                "WARNING",
            )
        if missing_requested_profile:
            return self.handle_exception(
                "No BUSCO run matched the requested proteome profile for one or more selected accessions.",
                {
                    "proteome_profile": requested_profile,
                    "missing_accessions": missing_requested_profile[:25],
                    "missing_count": len(missing_requested_profile),
                    "busco_library_id": busco_lib_id,
                },
            )

        selected_run_ids = sorted(
            {
                int(meta.get("run_id"))
                for meta in selected_busco_runs.values()
                if meta.get("run_id") is not None
            }
        )
        selected_profiles = sorted(
            {
                str(meta.get("proteome_profile"))
                for meta in selected_busco_runs.values()
                if meta.get("proteome_profile")
            }
        )
        if selected_run_ids:
            preview_runs = ",".join(str(run_id) for run_id in selected_run_ids[:20])
            run_suffix = "" if len(selected_run_ids) <= 20 else ",..."
            profile_text = ",".join(selected_profiles) if selected_profiles else "none"
            self.log(
                f"Selected BUSCO runs for export ({len(selected_run_ids)} unique run_id): {preview_runs}{run_suffix}; profiles={profile_text}",
                "DEBUG",
            )
        self.selected_busco_runs = dict(selected_busco_runs)
        if selected_run_ids:
            self.log(
                f"Export will emit sequences from strict run-backed metadata ({len(selected_run_ids)} chosen BUSCO run(s)).",
                "INFO",
            )

        # Prepare output directory
        try:
            os.makedirs(self.out_dir, exist_ok=True)
        except OSError as e:
            return self.handle_exception("Failed to create output directory.", {"out_dir": self.out_dir, "error": str(e)})
        self.busco_families_dir = os.path.join(self.out_dir, "busco_families")
        try:
            os.makedirs(self.busco_families_dir, exist_ok=True)
        except OSError as e:
            return self.handle_exception(
                "Failed to create busco_families directory.",
                {"busco_families_dir": self.busco_families_dir, "error": str(e)},
            )
        try:
            self.db_manager.artifacts.register(
                owner_type="export_run",
                owner_id=self.task_id,
                artifact_type="export_root",
                path=self.out_dir,
                is_dir=True,
                format="directory",
                metadata={"library_id": self.library_id, "sequence_type": self.sequence_type},
            )
            self.db_manager.artifacts.register(
                owner_type="export_run",
                owner_id=self.task_id,
                artifact_type="export_busco_families_dir",
                path=self.busco_families_dir,
                is_dir=True,
                format="directory",
                sequence_kind="nucl" if self.sequence_type == "nucleotide" else "prot",
                metadata={"library_id": self.library_id},
            )
        except Exception as exc:  # boundary: optional artifact catalog metadata; export files remain usable.
            self.log(f"Failed to register export artifacts: {exc}", "WARNING")

        try:
            family_kept = self._load_strict_export_family_entries(
                selected_busco_runs,
                library_id=busco_lib_id,
            )
        except Exception as exc:  # boundary: strict BUSCO export source resolution failure becomes this task error.
            return self.handle_exception(
                "Failed to resolve strict BUSCO export sequences.",
                {"error": str(exc), "sequence_type": self.sequence_type},
            )
        self.log(
            f"Resolved strict run-backed BUSCO sequences for {len(family_kept)} families.",
            "INFO",
        )

        # If using a custom library (has parent), filter exported families to the library's cleaned list
        if lib_parent_id is not None:
            # Load cleaned BUSCO families list from library folder
            families_to_keep = set()
            try:
                # Prefer a canonical filename if present
                canonical = os.path.join(lib_location, "cleaned_busco_families.json") if lib_location else None
                if canonical and os.path.isfile(canonical):
                    with open(canonical) as f:
                        families_to_keep = set(json.load(f) or [])
                else:
                    # Fallback: find latest timestamped file
                    pattern = os.path.join(lib_location or "", "cleaned_busco_families_*.json")
                    files = sorted(glob.glob(pattern))
                    if files:
                        with open(files[-1]) as f:
                            families_to_keep = set(json.load(f) or [])
            except (OSError, UnicodeError, json.JSONDecodeError, TypeError) as e:
                return self.handle_exception("Failed to load cleaned BUSCO families for custom library.", {"library": lib_name, "location": lib_location, "error": str(e)})

            if families_to_keep:
                family_kept = {
                    fam: entries
                    for fam, entries in family_kept.items()
                    if fam in families_to_keep
                }

        # Apply paralog removal and decontamination filtering to per-family FASTAs
        paralog_clean: dict[str, set[str]] = {}
        resolved_paralog_run_id = None
        if self.disable_paralog_filter:
            self.log("Paralog filtering disabled by request.", "INFO")
        elif lib_parent_id is not None:
            try:
                resolved_paralog_run_id = self.db_manager.filtering.resolve_paralog_run_id(
                    target_library_id=int(self.library_id),
                    run_id=self.paralog_run_id,
                )
            except Exception as exc:  # boundary: optional paralog run resolution failure falls back to explicit id.
                self.log(f"Failed to resolve paralog run {self.paralog_run_id}: {exc}", "WARNING")
                resolved_paralog_run_id = self.paralog_run_id
            if self.require_paralog_filtering and not resolved_paralog_run_id:
                return self.handle_exception(
                    "Paralog filtering required but no usable paralog-filtering run was found.",
                    {"library_id": self.library_id, "selected": len(selected)},
                )
            try:
                rows = self.db_manager.filtering.get_paralog_results(
                    target_library_id=self.library_id,
                    run_id=resolved_paralog_run_id,
                )
                present_paralog_accessions: set[str] = set()
                for fam_id, busco_lib, target_lib, acc, clean in rows:
                    if target_lib != self.library_id or busco_lib != busco_lib_id:
                        continue
                    present_paralog_accessions.add(str(acc))
                    if clean:
                        paralog_clean.setdefault(str(fam_id), set()).add(acc)
                if self.require_paralog_filtering:
                    missing_paralog = sorted(acc for acc in selected if acc not in present_paralog_accessions)
                    if missing_paralog:
                        return self.handle_exception(
                            "Paralog filtering required but no usable paralog-filtering results were found for one or more selected accessions.",
                            {
                                "library_id": self.library_id,
                                "missing_accessions": missing_paralog[:25],
                                "missing_count": len(missing_paralog),
                                "paralog_run_id": resolved_paralog_run_id,
                            },
                        )
            except Exception as exc:  # boundary: optional paralog filtering is skipped unless explicitly required above.
                self.log(f"Failed to load paralog filtering results; skipping paralog-based filtering: {exc}", "WARNING")
        if paralog_clean:
            self.log(
                f"Paralog filter loaded for {len(paralog_clean)} families (sample: {','.join(list(paralog_clean.keys())[:5])}).",
                "DEBUG",
            )

        decont_allowed: dict[str, set[str]] = {}
        if decont_run_id:
            try:
                votes = self.db_manager.filtering.get_decontamination_votes(target_library_id=self.library_id, run_id=decont_run_id)
                for (
                    fam_id,
                    busco_lib,
                    target_lib,
                    acc,
                    run_id_val,
                    expected_taxid,
                    best_taxid,
                    runner_taxid,
                    rank,
                    best_bitscore,
                    delta_bitscore,
                    decision,
                    top_hits_json,
                    *_vote_extra,
                ) in votes:
                    if decision in ("support", "weak") and target_lib == self.library_id and busco_lib == busco_lib_id:
                        decont_allowed.setdefault(str(fam_id), set()).add(acc)
            except Exception as exc:  # boundary: optional decontamination filtering is skipped unless explicitly required above.
                self.log(f"Failed to load decontamination votes for run {decont_run_id}; skipping decontamination filtering: {exc}", "WARNING")
        if decont_allowed:
            self.log(
                f"Decont filter loaded for {len(decont_allowed)} families (sample: {','.join(list(decont_allowed.keys())[:5])}).",
                "DEBUG",
            )

        rescued_duplicates: dict[tuple[str, str], dict[str, Any]] = {}
        if self.rescue_duplicates:
            supported_decisions = ["support"]
            if not self.strict_decontamination:
                supported_decisions.append("weak")
            if self.allow_ambiguous_contaminants:
                supported_decisions.append("unknown")
            accessions_by_run = {}
            if decont_run_id and not self.disable_decont_filter:
                primary_decont = self.db_manager.filtering.get_decontamination_accessions(
                    target_library_id=self.library_id,
                    run_id=decont_run_id,
                    accessions=selected,
                )
                if primary_decont:
                    accessions_by_run[(str(decont_run_id), int(self.library_id))] = sorted(primary_decont)
                if lib_parent_id:
                    fallback_decont = self.db_manager.filtering.get_decontamination_accessions(
                        target_library_id=lib_parent_id,
                        run_id=decont_run_id,
                        accessions=selected,
                    )
                    fallback_only = sorted(set(fallback_decont) - set(primary_decont))
                    if fallback_only:
                        accessions_by_run[(str(decont_run_id), int(lib_parent_id))] = fallback_only
            try:
                rescued_duplicates = self.db_manager.filtering.get_rescued_duplicate_copies(
                    target_library_id=self.library_id,
                    busco_library_id=busco_lib_id,
                    accessions=selected,
                    include_paralog=bool(not self.disable_paralog_filter and lib_parent_id is not None),
                    include_decontam=bool(decont_run_id and not self.disable_decont_filter),
                    accessions_by_run=accessions_by_run,
                    supported_decisions=supported_decisions,
                )
            except Exception as exc:  # boundary: duplicate rescue is optional export enrichment.
                self.log(f"Failed to load rescued duplicate BUSCOs; continuing without duplicate rescue: {exc}", "WARNING")
                rescued_duplicates = {}

        # First pass: apply filters and collect occupancy info
        strict_family_entries = dict(family_kept)
        family_kept = {}
        accession_counts: dict[str, int] = {}
        trimmed_sequences = 0
        total_families_seen = 0
        family_stats: dict[str, dict[str, float]] = {}
        def _header_matches_rescue(header: str, rescue_info: dict[str, Any]) -> bool:
            raw = header[1:] if header.startswith(">") else header
            first = raw.split()[0] if raw else ""
            query_id = str(rescue_info.get("query_id") or "")
            query_header = str(rescue_info.get("query_header") or "")
            token = query_id.split()[0] if query_id else ""
            if query_header and raw == query_header:
                return True
            if query_header and query_header in raw:
                return True
            if query_id and (raw == query_id or first == query_id or (token and token == first) or query_id in raw):
                return True
            return False

        total_input_sequences = sum(len(entries) for entries in strict_family_entries.values())
        self.log(
            f"Strict export sequence sourcing: families={len(strict_family_entries)} sequences={total_input_sequences}",
            "INFO",
        )
        for fam, source_entries in strict_family_entries.items():
            allow_paralog = paralog_clean.get(fam)
            allow_decont = decont_allowed.get(fam) if decont_allowed else None
            keepers: list[tuple[str, str, str]] = []
            family_present_accessions: set[str] = set()
            total_seq = len(source_entries)
            trimmed_pd = 0
            seen_entries: set[tuple[str, str]] = set()
            try:
                for acc, h, seq_str in source_entries:
                    rescue_info = rescued_duplicates.get((str(acc), str(fam)))
                    if rescue_info:
                        ok_paralog = _header_matches_rescue(h, rescue_info)
                        ok_decont = ok_paralog
                    else:
                        ok_paralog = (allow_paralog is None) or (acc in allow_paralog)
                        ok_decont = (allow_decont is None) or (acc in allow_decont)
                        if allow_paralog is None and allow_decont is None:
                            ok_paralog = True
                            ok_decont = True
                    if ok_paralog and ok_decont:
                        key = (acc, seq_str)
                        if key not in seen_entries:
                            keepers.append((acc, h.rstrip(), seq_str))
                            seen_entries.add(key)
                            if acc not in family_present_accessions:
                                accession_counts[acc] = accession_counts.get(acc, 0) + 1
                                family_present_accessions.add(acc)
                    else:
                        trimmed_sequences += 1
                        trimmed_pd += 1
                if keepers:
                    family_kept[fam] = keepers
                    total_families_seen += 1
                    family_stats[fam] = {
                        "total_seq": total_seq,
                        "kept_after_pd": len(keepers),
                        "trimmed_pd": trimmed_pd,
                    }
            except Exception as exc:  # boundary: required family filtering failure becomes this task error.
                return self.handle_exception("Failed to filter exported BUSCO families.", {"family": fam, "error": str(exc)})

        # Taxa occupancy: drop accessions that occur in too few families
        allowed_accessions = set(selected)
        if total_families_seen > 0 and self.min_taxa_occupancy > 0:
            threshold_count = self.min_taxa_occupancy * total_families_seen
            allowed_accessions = {acc for acc in allowed_accessions if accession_counts.get(acc, 0) >= threshold_count}
        # If nothing remains, keep selected set and warn instead of expanding to non-selected headers
        if not allowed_accessions:
            self.log(
                "No selected accessions passed taxa occupancy; retaining selected set to avoid expanding by header parsing.",
                "WARNING",
            )
            allowed_accessions = set(selected)
        removed_accessions = sorted(set(selected) - set(allowed_accessions))
        self.log(
            f"Export occupancy: selected={len(selected)} total_families_seen={total_families_seen} "
            f"accessions_with_sequences={len(accession_counts)} allowed_after_taxa={len(allowed_accessions)} "
            f"min_taxa_occupancy={self.min_taxa_occupancy}",
            "INFO",
        )
        if removed_accessions:
            preview = ",".join(removed_accessions[:25])
            suffix = "" if len(removed_accessions) <= 25 else ",..."
            self.log(
                f"Taxa removed by min_taxa_occupancy ({len(removed_accessions)}): {preview}{suffix}",
                "INFO",
            )

        report_accessions = sorted(allowed_accessions)
        self._write_taxa_occupancy_report(
            selected,
            allowed_accessions=allowed_accessions,
            accession_counts=accession_counts,
            total_families_seen=total_families_seen,
        )
        acc_to_tax = self._fetch_accession_taxids(report_accessions)
        required_membership = self._build_required_membership(
            report_accessions,
            required_taxids,
            acc_to_tax,
        )
        if self.write_lineage_csv:
            self._write_lineage_csv(report_accessions)
        if self.write_busco_report:
            self._write_busco_report(
                report_accessions,
                busco_library_id=self.library_id,
                decont_run_id=report_decont_run_id,
                custom_library=bool(lib_parent_id is not None),
            )
        if self.write_busco_family_matrix:
            self._write_busco_family_matrix(
                report_accessions,
                library_id=self.library_id,
                busco_library_id=busco_lib_id,
                decont_run_id=report_decont_run_id,
                custom_library=bool(lib_parent_id is not None),
            )
        self._write_export_parameters(
            selected=report_accessions,
            busco_library_id=busco_lib_id,
            busco_run_selection=selected_busco_runs,
            decont_selection_run=decont_selection_run,
            resolved_decont_run=decont_run_id,
            report_decont_run=report_decont_run_id,
            custom_library=bool(lib_parent_id is not None),
        )

        # Second pass: enforce family occupancy and rewrite files
        removed_files = 0
        final_written = 0
        report_rows = []
        header_template = self.header_template or "ACCESSION_TAXON_BUSCO"
        allowed_separators = set(".|_-:[]")
        if not self.retain_headers:
            tokens = {
                "ACCESSION",
                "TAXON",
                "KINGDOM",
                "PHYLUM",
                "CLASS",
                "ORDER",
                "FAMILY",
                "GENUS",
                "SPECIES",
                "RANK",
                "BUSCO",
                "SEQUENCE",
                "LENGTH",
                "GENE",
                "TAXID",
                "BITSCORE",
            }
            idx = 0
            while idx < len(header_template):
                ch = header_template[idx]
                if ch.isalpha():
                    matched = None
                    for token in sorted(tokens, key=len, reverse=True):
                        if header_template.startswith(token, idx):
                            matched = token
                            break
                    if not matched:
                        return self.handle_exception(
                            "Invalid header template token. Use ACCESSION/TAXON/KINGDOM/PHYLUM/CLASS/ORDER/FAMILY/GENUS/SPECIES/RANK/BUSCO/SEQUENCE/LENGTH/GENE/TAXID/BITSCORE.",
                            {"header": header_template},
                        )
                    idx += len(matched)
                else:
                    if ch not in allowed_separators:
                        return self.handle_exception(
                            "Invalid header template separator. Allowed: . | _ - : [ ]",
                            {"header": header_template, "char": ch},
                        )
                    idx += 1

        taxon_cache: dict[str, dict[str, str]] = {}
        if not self.retain_headers:
            lineage_cache: dict[int, list[tuple[int, str, str, Optional[int]]]] = {}
            rank_token = self._effective_header_rank(header_template)
            for acc in report_accessions:
                taxid_val = acc_to_tax.get(acc)
                lineage_rows = []
                if taxid_val is not None:
                    lineage_rows = lineage_cache.get(taxid_val)
                    if lineage_rows is None:
                        lineage_rows = self.db_manager.genomes.get_lineage_root_to_leaf(taxid_val) or []
                        lineage_cache[taxid_val] = lineage_rows
                taxon_cache[acc] = self._build_header_taxon_info(
                    lineage_rows,
                    taxid_val,
                    rank_token=rank_token,
                )

        gene_cache: dict[str, str] = {}
        bitscore_cache: dict[tuple[str, str], Optional[float]] = {}
        if not self.retain_headers:
            fam_ids = list(family_kept.keys())
            rows = self.db_manager.libraries.get_busco_descriptions(busco_lib_id, fam_ids)
            for fam_id, _lib_id, desc, _link in rows or []:
                gene_cache[str(fam_id)] = str(desc or "")
            if lib_parent_id is not None:
                missing = [fam for fam in fam_ids if fam not in gene_cache]
                if missing:
                    rows = self.db_manager.libraries.get_busco_descriptions(self.library_id, missing)
                    for fam_id, _lib_id, desc, _link in rows or []:
                        gene_cache[str(fam_id)] = str(desc or "")
            bitscore_cache = self._fetch_family_bitscores(report_accessions, fam_ids, busco_lib_id)
        for fam, entries in list(family_kept.items()):
            path = os.path.join(self.busco_families_dir, f"{fam}.fasta")
            filtered = [(acc, h, s) for acc, h, s in entries if acc in allowed_accessions]
            if required_clauses:
                present_required: set[int] = set()
                present_accessions = {acc for acc, _header, _seq in filtered}
                for acc, _header, _seq in filtered:
                    present_required.update(required_membership.get(acc, set()))
                missing_clause_labels = self._missing_required_clause_labels(
                    required_clauses,
                    present_taxids=present_required,
                    present_accessions=present_accessions,
                )
                if missing_clause_labels:
                    if os.path.exists(path):
                        try:
                            os.remove(path)
                        except OSError as exc:
                            self.log(f"Failed to remove filtered BUSCO file {path}: {exc}", "WARNING")
                    removed_files += 1
                    report_rows.append(
                        {
                            "family": fam,
                            "total_seq": family_stats.get(fam, {}).get("total_seq", 0),
                            "kept_after_pd": family_stats.get(fam, {}).get("kept_after_pd", 0),
                            "trimmed_pd": family_stats.get(fam, {}).get("trimmed_pd", 0),
                            "removed_due_taxa": len(entries) - len(filtered),
                            "final_kept": len(filtered),
                            "occupancy": len(filtered) / max(1, len(allowed_accessions)),
                            "reason": f"missing_required_clauses:{';'.join(missing_clause_labels)}",
                        }
                    )
                    continue
            # If filtering by allowed_accessions drops everything but we had keepers, fallback to original keepers
            if not filtered and entries:
                filtered = list(entries)
            if not filtered:
                if os.path.exists(path):
                    try:
                        os.remove(path)
                    except OSError as exc:
                        self.log(f"Failed to remove low-occupancy BUSCO file {path}: {exc}", "WARNING")
                removed_files += 1
                report_rows.append(
                    {
                        "family": fam,
                        "total_seq": family_stats.get(fam, {}).get("total_seq", 0),
                        "kept_after_pd": family_stats.get(fam, {}).get("kept_after_pd", 0),
                        "trimmed_pd": family_stats.get(fam, {}).get("trimmed_pd", 0),
                        "removed_due_taxa": len(entries),
                        "final_kept": 0,
                        "occupancy": 0.0,
                        "reason": "no_sequences_after_taxa",
                    }
                )
                continue
            occ = len(filtered) / max(1, len(allowed_accessions))
            if self.min_occupancy and occ < self.min_occupancy:
                if os.path.exists(path):
                    try:
                        os.remove(path)
                    except OSError as exc:
                        self.log(f"Failed to remove filtered BUSCO file {path}: {exc}", "WARNING")
                removed_files += 1
                report_rows.append(
                    {
                        "family": fam,
                        "total_seq": family_stats.get(fam, {}).get("total_seq", 0),
                        "kept_after_pd": family_stats.get(fam, {}).get("kept_after_pd", 0),
                        "trimmed_pd": family_stats.get(fam, {}).get("trimmed_pd", 0),
                        "removed_due_taxa": len(entries) - len(filtered),
                        "final_kept": len(filtered),
                        "occupancy": occ,
                        "reason": "below_min_occupancy",
                    }
                )
                continue
            try:
                copy_totals: dict[str, int] = {}
                for acc, _header, _seq in filtered:
                    copy_totals[acc] = copy_totals.get(acc, 0) + 1
                copy_seen: dict[str, int] = {}
                with open(path, "w") as handle:
                    for acc, header, seq in filtered:
                        if self.retain_headers:
                            out_header = header
                        else:
                            taxon_info = taxon_cache.get(acc, {})
                            out_header = self._render_export_header(
                                accession=acc,
                                family_id=fam,
                                source_header=header,
                                sequence=seq,
                                taxon_info=taxon_info,
                                acc_to_tax=acc_to_tax,
                                gene_cache=gene_cache,
                                bitscore_cache=bitscore_cache,
                                header_template=header_template,
                            )
                            if (
                                self.include_duplicated
                                and copy_totals.get(acc, 0) > 1
                                and self._header_needs_copy_suffix(header_template=header_template)
                            ):
                                copy_seen[acc] = copy_seen.get(acc, 0) + 1
                                out_header = f"{out_header}_{copy_seen[acc]}"
                        handle.write(f"{out_header}\n{seq}\n")
                final_written += 1
                report_rows.append(
                    {
                        "family": fam,
                        "total_seq": family_stats.get(fam, {}).get("total_seq", 0),
                        "kept_after_pd": family_stats.get(fam, {}).get("kept_after_pd", 0),
                        "trimmed_pd": family_stats.get(fam, {}).get("trimmed_pd", 0),
                        "removed_due_taxa": len(entries) - len(filtered),
                        "final_kept": len(filtered),
                        "occupancy": occ,
                        "reason": "kept",
                    }
                )
            except (OSError, UnicodeError) as exc:
                return self.handle_exception("Failed to write filtered BUSCO family.", {"file": path, "error": str(exc)})

        if removed_files or trimmed_sequences:
            self.log(
                f"Filtered export using paralog/decontamination results: removed {removed_files} family files after occupancy checks, "
                f"trimmed {trimmed_sequences} sequences.",
                "INFO",
            )
        self.log(
            f"Export filter summary: selected_accessions={len(selected)} allowed_accessions={len(allowed_accessions)} "
            f"families_seen={total_families_seen} families_written={final_written} removed_files={removed_files}",
            "INFO",
        )
        # Write report if requested
        try:
            report_path = os.path.join(self.report_dir, "export_filter_report.tsv")
            with open(report_path, "w") as rep:
                rep.write("family\ttotal_seq\tkept_after_pd\ttrimmed_pd\tremoved_due_taxa\tfinal_kept\toccupancy\treason\n")
                for row in report_rows:
                    rep.write(
                        "\t".join(
                            [
                                str(row.get("family", "")),
                                str(row.get("total_seq", 0)),
                                str(row.get("kept_after_pd", 0)),
                                str(row.get("trimmed_pd", 0)),
                                str(row.get("removed_due_taxa", 0)),
                                str(row.get("final_kept", 0)),
                                f"{row.get('occupancy', 0):.3f}",
                                str(row.get("reason", "")),
                            ]
                        )
                        + "\n"
                    )
                # Even if empty, ensure file exists with header
            self.log(f"Wrote export filter report to {report_path}", "DEBUG")
        except (OSError, csv.Error) as exc:
            self.log(f"Failed to write export filter report: {exc}", "WARNING")

        try:
            fam_count = len([f for f in os.listdir(self.busco_families_dir) if f.endswith(".fasta")])
        except OSError:
            fam_count = 0
        self.log(
            f"Export completed: {fam_count} BUSCO family files written to {self.busco_families_dir}",
            "INFO",
        )
        return True

    def start(self):
        try:
            return super().start()
        finally:
            self._disable_export_log_copy()
