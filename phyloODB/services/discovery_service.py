from __future__ import annotations

import csv
import glob
import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ..database import DBManager
from ..db.errors import PhyloODBDatabaseError
from ..proteome_profile_utils import DEFAULT_CLEAN_PROFILE, staged_busco_input_profile_name
from ..proteome_state import summarize_proteome_state
from ..tasks.utilities.ncbi_helper import NCBIHelper


NCBI_ASSEMBLY_RE = re.compile(r"^(GCA|GCF)_\d+\.\d+$", re.IGNORECASE)


@dataclass
class DiscoveryEvent:
    root_id: int
    root_label: Optional[str]
    scan_scope: str
    scan_path: str
    accession: str
    action: str
    detail: str


@dataclass
class DiscoveryReport:
    dry_run: bool
    overwrite: bool
    attempt_knowledge_update: bool
    events: List[DiscoveryEvent] = field(default_factory=list)

    def add(
        self,
        *,
        root_id: int,
        root_label: Optional[str],
        scan_scope: str,
        scan_path: str,
        accession: str,
        action: str,
        detail: str,
    ) -> None:
        self.events.append(
            DiscoveryEvent(
                root_id=int(root_id),
                root_label=(str(root_label) if root_label is not None else None),
                scan_scope=str(scan_scope),
                scan_path=os.path.abspath(str(scan_path)),
                accession=accession,
                action=action,
                detail=detail,
            )
        )

    def as_lines(self) -> List[str]:
        lines = ["root_id\troot_label\tscan_scope\tscan_path\taccession\taction\tdetail"]
        for event in self.events:
            lines.append(
                f"{event.root_id}\t{event.root_label or ''}\t{event.scan_scope}\t{event.scan_path}\t{event.accession}\t{event.action}\t{event.detail}"
            )
        return lines


class DiscoveryService:
    def __init__(self, manager: DBManager):
        self.db = manager

    def discover(
        self,
        *,
        root: Optional[str] = None,
        path: Optional[str] = None,
        overwrite: bool = False,
        dry_run: bool = False,
        attempt_knowledge_update: bool = False,
    ) -> DiscoveryReport:
        if root and path:
            raise ValueError("Use only one of root or path for discovery scope.")
        report = DiscoveryReport(
            dry_run=dry_run,
            overwrite=overwrite,
            attempt_knowledge_update=attempt_knowledge_update,
        )
        for root_row, scan_scope, scan_path in self._resolve_scopes(root=root, path=path):
            self._discover_scope(
                root_row=root_row,
                scan_scope=scan_scope,
                scan_path=scan_path,
                overwrite=overwrite,
                dry_run=dry_run,
                attempt_knowledge_update=attempt_knowledge_update,
                report=report,
            )
        return report

    def _resolve_scopes(self, *, root: Optional[str], path: Optional[str]) -> List[tuple]:
        if root:
            root_row = self.db.storage.resolve_root_token(root, kind="genomes")
            base = os.path.abspath(str(root_row[3]))
            if not os.path.isdir(base):
                raise ValueError(f"Discovery root path does not exist or is not a directory: {base}")
            return [(root_row, "root", base)]
        if path:
            abs_path = os.path.abspath(str(path))
            if not os.path.isdir(abs_path):
                raise ValueError(f"Discovery path does not exist or is not a directory: {abs_path}")
            root_id, _rel = self.db.storage.detect_root_for_path(abs_path, kind="genomes")
            if root_id is None:
                raise ValueError(
                    f"Discovery path is not inside a registered genomes root: {abs_path}. "
                    "Register the genomes root first with 'storage add-root --kind genomes --base-path ...'."
                )
            root_row = self.db.storage.get_root(int(root_id))
            if not root_row:
                raise ValueError(f"Registered genomes root {root_id} could not be loaded.")
            base = os.path.abspath(str(root_row[3]))
            scan_scope = "root" if abs_path == base else "subtree"
            return [(root_row, scan_scope, abs_path)]
        rows = self.db.storage.list_roots(kind="genomes") or []
        if not rows:
            raise ValueError("No genomes roots are registered. Add one with 'storage add-root --kind genomes --base-path ...'.")
        scopes: List[tuple] = []
        for row in rows:
            base = os.path.abspath(str(row[3]))
            if os.path.isdir(base):
                scopes.append((row, "root", base))
        return scopes

    def _discover_scope(
        self,
        *,
        root_row,
        scan_scope: str,
        scan_path: str,
        overwrite: bool,
        dry_run: bool,
        attempt_knowledge_update: bool,
        report: DiscoveryReport,
    ) -> None:
        root_id = int(root_row[0])
        root_label = str(root_row[2]) if root_row[2] is not None else None
        if scan_scope == "subtree":
            accession = os.path.basename(os.path.abspath(scan_path)).strip()
            if accession:
                direct_known = self.db.genomes.get(accession)
                direct_taxid = self._read_taxid_file(scan_path)
                if direct_known is not None or direct_taxid is not None or self._looks_like_ncbi_assembly(accession):
                    try:
                        self._discover_accession(
                            root_id=root_id,
                            root_label=root_label,
                            scan_scope=scan_scope,
                            scan_path=scan_path,
                            accession=accession,
                            folder=scan_path,
                            overwrite=overwrite,
                            dry_run=dry_run,
                            attempt_knowledge_update=attempt_knowledge_update,
                            report=report,
                        )
                    except PhyloODBDatabaseError:
                        raise
                    except Exception as exc:  # boundary: isolate one discovered storage root
                        report.add(
                            root_id=root_id,
                            root_label=root_label,
                            scan_scope=scan_scope,
                            scan_path=scan_path,
                            accession=accession,
                            action="error",
                            detail=str(exc),
                        )
                    return
        for entry in sorted(os.listdir(scan_path)):
            folder = os.path.join(scan_path, entry)
            if not os.path.isdir(folder):
                continue
            accession = entry.strip()
            if not accession:
                continue
            try:
                self._discover_accession(
                    root_id=root_id,
                    root_label=root_label,
                    scan_scope=scan_scope,
                    scan_path=scan_path,
                    accession=accession,
                    folder=folder,
                    overwrite=overwrite,
                    dry_run=dry_run,
                    attempt_knowledge_update=attempt_knowledge_update,
                    report=report,
                )
            except PhyloODBDatabaseError:
                raise
            except Exception as exc:  # boundary: isolate one explicitly requested scan
                report.add(
                    root_id=root_id,
                    root_label=root_label,
                    scan_scope=scan_scope,
                    scan_path=scan_path,
                    accession=accession,
                    action="error",
                    detail=str(exc),
                )

    def _discover_accession(
        self,
        *,
        root_id: int,
        root_label: Optional[str],
        scan_scope: str,
        scan_path: str,
        accession: str,
        folder: str,
        overwrite: bool,
        dry_run: bool,
        attempt_knowledge_update: bool,
        report: DiscoveryReport,
    ) -> None:
        known = self.db.genomes.get(accession)
        known_path = self.db.genomes.resolve_path(accession) if known else None
        folder = os.path.abspath(folder)

        if known is None:
            taxid = self._read_taxid_file(folder)
            if taxid is not None:
                report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="new-custom", detail=f"taxid={taxid} path={folder}")
                if not dry_run:
                    self._insert_custom_accession(accession, folder, taxid)
            elif attempt_knowledge_update and self._looks_like_ncbi_assembly(accession):
                report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="knowledge-update", detail=f"attempt path={folder}")
                if not dry_run:
                    success = self._attempt_knowledge_update(accession, folder)
                    if not success:
                        report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="warning", detail="knowledge update failed; accession skipped")
                        return
            else:
                report.add(
                    root_id=root_id,
                    root_label=root_label,
                    scan_scope=scan_scope,
                    scan_path=scan_path,
                    accession=accession,
                    action="warning",
                    detail="unknown accession skipped; add taxid file or use --attempt-knowledge-update",
                )
                return
        elif known_path and os.path.abspath(known_path) != folder and not overwrite:
            report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="warning", detail=f"known accession at different path: current={known_path} discovered={folder}")
            return
        elif overwrite:
            report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="overwrite", detail=f"rebind path to {folder}")
            if not dry_run:
                self.db.storage.bind_genome_location(accession, folder, kind="genomes")
        elif known is not None and known_path is None:
            report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="rebind", detail=f"bind path to {folder}")
            if not dry_run:
                self.db.storage.bind_genome_location(accession, folder, kind="genomes")

        if known is not None and known_path == folder and not overwrite:
            report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="known", detail=f"scan existing binding path={folder}")

        if not dry_run and (known is None or overwrite):
            self.db.storage.bind_genome_location(accession, folder, kind="genomes")
        if not dry_run:
            self._apply_discovered_genome_state(accession, folder)

        discovered_runs = self._discover_busco_runs(
            accession,
            folder,
            root_id=root_id,
            root_label=root_label,
            scan_scope=scan_scope,
            scan_path=scan_path,
            overwrite=overwrite,
            dry_run=dry_run,
            report=report,
        )
        if not discovered_runs:
            report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="scan", detail="no BUSCO runs discovered")

    def _looks_like_ncbi_assembly(self, accession: str) -> bool:
        return bool(NCBI_ASSEMBLY_RE.match(str(accession or "")))

    def _read_taxid_file(self, folder: str) -> Optional[int]:
        path = os.path.join(folder, "taxid")
        if not os.path.isfile(path):
            path = os.path.join(folder, "taxid.txt")
            if not os.path.isfile(path):
                return None
        try:
            token = Path(path).read_text(encoding="utf-8").strip().split()[0]
            value = int(token)
        except (OSError, UnicodeError, IndexError, TypeError, ValueError):
            return None
        return value if value > 0 else None

    def _apply_discovered_genome_state(self, accession: str, folder: str) -> None:
        state = summarize_proteome_state(folder)
        sync = self.db.proteomes.sync_profiles_from_filesystem(str(accession), str(folder), set_default=True)
        has_protein = bool(sync.get("has_protein")) or bool(state.protein_flag)
        isoforms_cleaned = bool(self.db.proteomes.get_default_cleaned_profile_name(str(accession)))
        self.db.genomes.update_status(accession, 1, (datetime.now(), folder, has_protein))
        self.db.genomes.set_protein(accession, has_protein)
        self.db.genomes.set_isoforms_cleaned(accession, isoforms_cleaned)

    def _insert_custom_accession(self, accession: str, folder: str, taxid: int) -> None:
        state = summarize_proteome_state(folder)
        self.db.genomes.insert_assembly(
            {
                "accession": accession,
                "origin": "local",
                "release_date": datetime.now().strftime("%Y-%m-%d"),
            }
        )
        self.db.genomes.insert(
            {
                "accession": accession,
                "taxid": int(taxid),
                "protein": state.protein_flag,
                "isoforms_cleaned": state.isoforms_cleaned_flag,
                "location": folder,
                "status": 1,
                "dl_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }
        )

    def _attempt_knowledge_update(self, accession: str, folder: str) -> bool:
        email = self.db.env.get("EMAIL")
        if not email:
            return False
        api_key = self.db.env.get("NCBI_API_KEY")
        helper = NCBIHelper(email=email, db_manager=self.db, api_key=api_key)
        assembly_dataset, tax_info_dataset, taxonomy_update_dict, genome_dataset = helper.fetch_assemblies_v2([accession], accessions_only=True)
        if not assembly_dataset or not genome_dataset:
            return False
        if taxonomy_update_dict:
            self.db.genomes.insert_taxonomy_information(taxonomy_update_dict)
        self.db.genomes.upsert_assembly(assembly_dataset[0])
        genome_record = dict(genome_dataset[0])
        state = summarize_proteome_state(folder)
        genome_record["location"] = folder
        genome_record["status"] = 1
        genome_record["protein"] = state.protein_flag
        genome_record["isoforms_cleaned"] = state.isoforms_cleaned_flag
        self.db.genomes.upsert(genome_record)
        return True

    def _summary_json_for_profile_inference(self, result_dir: str, run_dir: str) -> Optional[str]:
        candidates = []
        candidates.extend(sorted(glob.glob(os.path.join(result_dir, "short_summary*.json"))))
        candidates.extend(sorted(glob.glob(os.path.join(run_dir, "short_summary*.json"))))
        return candidates[0] if candidates else None

    def _infer_proteome_profile_id(
        self,
        accession: str,
        input_mode: str,
        result_dir: str,
        run_dir: str,
    ) -> tuple[Optional[int], Optional[str]]:
        if str(input_mode or "").strip().lower() != "protein":
            return None, None
        summary_json = self._summary_json_for_profile_inference(result_dir, run_dir)
        input_path = None
        if summary_json and os.path.isfile(summary_json):
            try:
                payload = json.loads(Path(summary_json).read_text(encoding="utf-8"))
            except (OSError, UnicodeError, json.JSONDecodeError):
                payload = None
            params = payload.get("parameters") if isinstance(payload, dict) else None
            if isinstance(params, dict):
                input_path = params.get("in")
        row = self.db.proteomes.find_profile_by_path(str(accession), input_path)
        if row:
            return int(row[0]), None
        profile_name = staged_busco_input_profile_name(input_path)
        if profile_name:
            if profile_name == DEFAULT_CLEAN_PROFILE:
                row = self.db.proteomes.get_default_cleaned_profile(str(accession))
            else:
                row = self.db.proteomes.get_profile(str(accession), profile_name)
            if row:
                return int(row[0]), None
            return None, (
                f"summary input profile '{profile_name}' inferred from '{input_path}' is not registered "
                f"for accession '{accession}'"
            )
        default_row = self.db.proteomes.get_default_profile(str(accession))
        return (int(default_row[0]), None) if default_row and default_row[0] is not None else (None, None)

    def _discover_busco_runs(
        self,
        accession: str,
        folder: str,
        *,
        root_id: int,
        root_label: Optional[str],
        scan_scope: str,
        scan_path: str,
        overwrite: bool,
        dry_run: bool,
        report: DiscoveryReport,
    ) -> int:
        total = 0
        purged: set[int] = set()
        for discovered_result_dir in self._result_dirs(folder):
            lineage_name = self._lineage_name_from_result_dir(discovered_result_dir)
            if not lineage_name:
                report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="warning", detail=f"could not infer lineage from {discovered_result_dir}")
                continue
            library_id = self.db.libraries.get_id(lineage_name)
            if not library_id:
                report.add(root_id=root_id, root_label=root_label, scan_scope=scan_scope, scan_path=scan_path, accession=accession, action="warning", detail=f"unknown BUSCO lineage '{lineage_name}' for {discovered_result_dir}")
                continue
            result_dir = self._canonical_dir_path(discovered_result_dir) or os.path.abspath(discovered_result_dir)
            run_dirs = self._run_dirs(result_dir)
            if overwrite and not dry_run and int(library_id) not in purged:
                self._purge_busco_state(accession, library_id)
                purged.add(int(library_id))
            for run_dir in run_dirs:
                total += 1
                existing_run_id = None if overwrite else self._find_existing_run(accession, int(library_id), result_dir, run_dir)
                action = "busco-run-update" if existing_run_id else "busco-run"
                report.add(
                    root_id=root_id,
                    root_label=root_label,
                    scan_scope=scan_scope,
                    scan_path=scan_path,
                    accession=accession,
                    action=action,
                    detail=f"lineage={lineage_name} run={os.path.basename(run_dir)}" + (f" run_id={existing_run_id}" if existing_run_id else ""),
                )
                if dry_run:
                    continue
                ingested, warning = self._ingest_run(
                    accession,
                    library_id,
                    lineage_name,
                    result_dir,
                    run_dir,
                    existing_run_id=existing_run_id,
                )
                if not ingested:
                    report.add(
                        root_id=root_id,
                        root_label=root_label,
                        scan_scope=scan_scope,
                        scan_path=scan_path,
                        accession=accession,
                        action="warning",
                        detail=(
                            f"skipped BUSCO run ingestion for lineage={lineage_name} "
                            f"run={os.path.basename(run_dir)}: {warning or 'unknown reason'}"
                        ),
                    )
        return total

    def _canonical_dir_path(self, path: str) -> Optional[str]:
        resolved = os.path.realpath(os.path.abspath(path))
        return resolved if os.path.isdir(resolved) else None

    def _result_dirs(self, folder: str) -> List[str]:
        candidates = sorted(glob.glob(os.path.join(folder, "*_results*")))
        chosen: Dict[str, str] = {}
        for candidate in candidates:
            if not os.path.isdir(candidate):
                continue
            canonical = self._canonical_dir_path(candidate)
            if not canonical:
                continue
            current = chosen.get(canonical)
            if current is None:
                chosen[canonical] = candidate
                continue
            current_name = os.path.basename(current)
            candidate_name = os.path.basename(candidate)
            current_score = (
                1 if "_results__" in current_name else 0,
                0 if os.path.islink(current) else 1,
                len(current_name),
            )
            candidate_score = (
                1 if "_results__" in candidate_name else 0,
                0 if os.path.islink(candidate) else 1,
                len(candidate_name),
            )
            if candidate_score > current_score:
                chosen[canonical] = candidate
        return sorted(chosen.values())

    def _run_dirs(self, result_dir: str) -> List[str]:
        candidates = sorted(glob.glob(os.path.join(result_dir, "run_*")))
        chosen: Dict[str, str] = {}
        for candidate in candidates:
            if not os.path.isdir(candidate):
                continue
            canonical = self._canonical_dir_path(candidate)
            if not canonical:
                continue
            current = chosen.get(canonical)
            if current is None:
                chosen[canonical] = candidate
                continue
            current_score = (0 if os.path.islink(current) else 1, len(os.path.basename(current)))
            candidate_score = (0 if os.path.islink(candidate) else 1, len(os.path.basename(candidate)))
            if candidate_score > current_score:
                chosen[canonical] = candidate
        return [self._canonical_dir_path(path) or os.path.abspath(path) for path in sorted(chosen.values())]

    def _lineage_name_from_result_dir(self, result_dir: str) -> Optional[str]:
        name = os.path.basename(result_dir)
        if "_results" not in name:
            return None
        return name.split("_results", 1)[0] or None

    def _purge_busco_state(self, accession: str, library_id: int) -> None:
        self.db.busco.delete_records(accession, library_id)
        self.db.cursor.execute("DELETE FROM BUSCO_Primary WHERE accession = ? AND library_id = ?", (accession, int(library_id)))
        self.db.cursor.execute("DELETE FROM BUSCO_Runs WHERE accession = ? AND library_id = ?", (accession, int(library_id)))
        self.db.commit()
        self.db.busco.invalidate_adjusted_results_for_busco_scope(
            int(library_id),
            accessions=[str(accession)],
            reason="busco_state_purged",
        )

    def _infer_pipeline(self, run_dir: str) -> str:
        base = os.path.basename(run_dir or "")
        if base.startswith("run_"):
            token = base[4:].split("_", 1)[0].strip().lower()
            if token in {"miniprot", "metaeuk", "augustus"}:
                return token
        return "miniprot"

    def _infer_input_mode(self, run_dir: str, *, pipeline: Optional[str] = None) -> str:
        pipeline_name = str(pipeline or "").strip().lower()
        for path in glob.glob(os.path.join(run_dir, "busco_sequences", "**", "*"), recursive=True):
            low = path.lower()
            if low.endswith((".fna", ".fna.gz")):
                return "genome"
            if low.endswith((".faa", ".faa.gz")):
                if pipeline_name != "miniprot":
                    return "protein"
        if pipeline_name in {"miniprot", "metaeuk", "augustus"}:
            return "genome"
        return "protein"

    def _full_table(self, run_dir: str) -> Optional[str]:
        rows = sorted(glob.glob(os.path.join(run_dir, "full_table*.tsv")))
        return rows[0] if rows else None

    def _summary_json_for_result_dir(self, result_dir: str) -> Optional[str]:
        rows = sorted(glob.glob(os.path.join(result_dir, "short_summary*.json")))
        return rows[0] if rows else None

    def _sequence_kind(self, path: Optional[str]) -> Optional[str]:
        if not path:
            return None
        low = str(path).lower()
        if low.endswith((".faa", ".faa.gz", ".pep", ".pep.gz")):
            return "prot"
        if low.endswith((".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            return "nucl"
        return None

    def _input_mode_for_path(self, path: Optional[str]) -> Optional[str]:
        kind = self._sequence_kind(path)
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
            value = raw.strip().lower() == "true" if isinstance(raw, str) else bool(raw)
            if value:
                true_flags.append(name)
        if len(true_flags) == 1:
            return true_flags[0]
        return None

    def _inspect_summary_input(self, result_dir: str, run_dir: str) -> Dict[str, Any]:
        summary_path = self._summary_json_for_profile_inference(result_dir, run_dir)
        result: Dict[str, Any] = {
            "summary_path": summary_path,
            "summary_input_path": None,
            "summary_mode_raw": None,
            "summary_input_mode": None,
            "summary_pipeline": None,
            "verified_input_mode": None,
        }
        if not summary_path or not os.path.isfile(summary_path):
            return result
        try:
            payload = json.loads(Path(summary_path).read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            return result
        params = payload.get("parameters") if isinstance(payload, dict) else None
        if not isinstance(params, dict):
            return result
        input_path = params.get("in")
        summary_mode_raw = params.get("mode")
        summary_input_mode = self._summary_mode_to_input_mode(summary_mode_raw)
        path_input_mode = self._input_mode_for_path(input_path)
        result["summary_input_path"] = input_path
        result["summary_mode_raw"] = summary_mode_raw
        result["summary_input_mode"] = summary_input_mode
        result["summary_pipeline"] = self._summary_flags_to_pipeline(params)
        result["verified_input_mode"] = summary_input_mode or path_input_mode
        return result

    def _parse_full_table(self, accession: str, run_dir: str, library_id: int) -> tuple[list[tuple], list[tuple], dict[str, int]]:
        table = self._full_table(run_dir)
        if not table:
            return [], [], {}
        family_data: list[tuple] = []
        family_locations: list[tuple] = []
        counts = {"Single copy BUSCOs": 0, "Multi copy BUSCOs": 0, "Fragmented BUSCOs": 0, "Missing BUSCOs": 0}
        seen: set[str] = set()
        with open(table, "r", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if not row or row[0].startswith("#") or len(row) < 2:
                    continue
                family_id = str(row[0])
                status = {"Complete": 1, "Duplicated": 2, "Fragmented": 3, "Missing": 4}.get(str(row[1]), 0)
                if family_id not in seen:
                    seen.add(family_id)
                    key = {1: "Single copy BUSCOs", 2: "Multi copy BUSCOs", 3: "Fragmented BUSCOs", 4: "Missing BUSCOs"}.get(status)
                    if key:
                        counts[key] += 1
                sequence = row[2] if len(row) > 2 and row[2] else None
                score = None
                length = None
                if len(row) > 6 and row[6]:
                    try:
                        score = float(row[6])
                    except (TypeError, ValueError):
                        score = None
                if len(row) > 7 and row[7]:
                    try:
                        length = int(row[7])
                    except (TypeError, ValueError):
                        length = None
                location = self._find_family_sequence(run_dir, family_id, status)
                family_data.append((family_id, library_id, accession, status, sequence, score, length))
                if location:
                    family_locations.append((family_id, library_id, accession, location))
        return family_data, family_locations, counts

    def _find_family_sequence(self, run_dir: str, family_id: str, status: int) -> Optional[str]:
        sub = {1: "single_copy_busco_sequences", 2: "multi_copy_busco_sequences", 3: "fragmented_busco_sequences"}.get(status)
        if not sub:
            return None
        for ext in (".faa", ".fna", ".fa", ".fasta"):
            candidate = os.path.join(run_dir, "busco_sequences", sub, f"{family_id}{ext}")
            if os.path.isfile(candidate):
                return candidate
        return None

    def _register_run_artifacts(self, run_id: int, accession: str, library_id: int, result_dir: str, run_dir: str, family_locations: List[Tuple[Any, ...]]) -> None:
        self.db.busco.register_run_artifact(run_id, "busco_result_root", result_dir, is_dir=True, format="directory")
        self.db.busco.register_run_artifact(run_id, "busco_run_dir", run_dir, is_dir=True, format="directory")
        sequences_dir = os.path.join(run_dir, "busco_sequences")
        if os.path.isdir(sequences_dir):
            self.db.busco.register_run_artifact(run_id, "busco_sequences_dir", sequences_dir, is_dir=True, format="directory")
        full_table = self._full_table(run_dir)
        if full_table:
            self.db.busco.register_run_artifact(run_id, "busco_full_table_tsv", full_table, format="tsv")
        summary_json = self._summary_json_for_result_dir(result_dir)
        if summary_json:
            self.db.busco.register_run_artifact(run_id, "busco_summary_json", summary_json, format="json")
        for family_id, _lib, _acc, location in family_locations:
            if not location or not os.path.exists(location):
                continue
            self.db.busco.register_family_artifact(
                run_id=run_id,
                family_id=str(family_id),
                library_id=library_id,
                accession=accession,
                path=location,
                sequence_kind=self._sequence_kind(location),
                role=str(family_id),
                format="fasta",
                metadata={"source": "discover"},
            )

    def _find_existing_run(self, accession: str, library_id: int, result_dir: str, run_dir: str) -> Optional[int]:
        run_dir = os.path.abspath(run_dir)
        result_dir = os.path.abspath(result_dir)
        row = self.db.busco.core.fetchone(
            """
            SELECT r.run_id
            FROM BUSCO_Runs r
            JOIN Artifacts a
              ON a.owner_type = 'busco_run'
             AND CAST(a.owner_id AS INTEGER) = r.run_id
             AND a.artifact_type = 'busco_run_dir'
            WHERE r.accession = ? AND r.library_id = ?
              AND a.absolute_path = ?
            ORDER BY r.run_id DESC
            LIMIT 1
            """,
            (str(accession), int(library_id), run_dir),
        )
        if row and row[0] is not None:
            return int(row[0])
        row = self.db.busco.core.fetchone(
            """
            SELECT run_id
            FROM BUSCO_Runs
            WHERE accession = ? AND library_id = ? AND result_dir = ?
            ORDER BY run_id DESC
            LIMIT 1
            """,
            (str(accession), int(library_id), result_dir),
        )
        return int(row[0]) if row and row[0] is not None else None

    def _reset_run_state(self, run_id: int) -> None:
        self.db.busco.invalidate_adjusted_results_for_busco_run(int(run_id), reason="busco_run_state_reset")
        self.db.busco.core.execute("DELETE FROM BUSCO_Run_Family_Artifacts WHERE run_id = ?", (int(run_id),))
        self.db.busco.core.execute("DELETE FROM BUSCO_Run_Family_Locations WHERE run_id = ?", (int(run_id),))
        self.db.busco.core.execute("DELETE FROM BUSCO_Run_Family_Data WHERE run_id = ?", (int(run_id),))
        self.db.busco.core.execute(
            "DELETE FROM Artifacts WHERE owner_type = 'busco_run' AND owner_id = ?",
            (str(int(run_id)),),
        )
        self.db.commit()

    def _ensure_busco_descriptions(self, rows: List[Tuple[Any, ...]]) -> None:
        pairs = {(str(row[0]), int(row[1])) for row in rows if row and row[0] and row[1] is not None}
        if not pairs:
            return
        self.db.cursor.executemany(
            """
            INSERT OR IGNORE INTO BUSCO_descriptions (family_id, library_id, description, link)
            VALUES (?, ?, ?, ?)
            """,
            [(family_id, library_id, None, "discover") for family_id, library_id in pairs],
        )
        self.db.commit()
        for _family_id, library_id in pairs:
            self.db.busco.invalidate_adjusted_results_for_library(
                int(library_id),
                reason="library_busco_descriptions_updated",
            )

    def _ingest_run(
        self,
        accession: str,
        library_id: int,
        lineage_name: str,
        result_dir: str,
        run_dir: str,
        *,
        existing_run_id: Optional[int] = None,
    ) -> tuple[bool, Optional[str]]:
        family_data, family_locations, counts = self._parse_full_table(accession, run_dir, library_id)
        self._ensure_busco_descriptions(family_data)
        self._ensure_busco_descriptions(family_locations)
        summary_input = self._inspect_summary_input(result_dir, run_dir)
        pipeline = str(summary_input.get("summary_pipeline") or self._infer_pipeline(run_dir))
        input_mode = str(
            summary_input.get("verified_input_mode")
            or self._infer_input_mode(run_dir, pipeline=pipeline)
        )
        proteome_profile_id, profile_warning = self._infer_proteome_profile_id(accession, input_mode, result_dir, run_dir)
        if profile_warning:
            return False, profile_warning
        run_id = existing_run_id
        if run_id is None:
            run_id = self.db.busco.create_run(
                accession=accession,
                library_id=library_id,
                lineage_name=lineage_name,
                input_mode=input_mode,
                pipeline=pipeline,
                pipeline_params_effective={},
                pipeline_params_source={"source": "discover"},
                busco_cli_args=[],
                result_dir=result_dir,
                proteome_profile_id=proteome_profile_id,
                status="completed",
            )
            if run_id is None:
                return False, "failed to create BUSCO run record"
        else:
            self._reset_run_state(int(run_id))
        self.db.busco.update_run(
            run_id,
            status="completed",
            result_dir=result_dir,
            counts=counts,
            completed=True,
            proteome_profile_id=proteome_profile_id,
        )
        if family_data:
            self.db.busco.add_run_family_data(run_id, family_data)
        if family_locations:
            self.db.busco.add_run_family_locations(run_id, list(set(family_locations)))
        self._register_run_artifacts(int(run_id), accession, library_id, result_dir, run_dir, list(set(family_locations)))
        self.db.busco.refresh_auto_primary_runs_for_accession(
            accession,
            library_id,
            updated_by="discover",
            policy="auto_best",
        )
        self.db.busco.add_results(accession, library_id, counts, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        return True, None
