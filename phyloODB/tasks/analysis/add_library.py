import tempfile
from ..utilities import OrthoBuscoAnalyzer
import time
import json
from abc import ABC, abstractmethod
from datetime import datetime
import os
import shutil
import csv
import glob
import re
from typing import Optional, Dict

from ..task import Task
from ...accession_utils import canonicalize_accession
from ...proteome_profile_utils import resolve_profile_selector
from ...selector_utils import (
    normalize_accessions,
    resolve_clade_to_taxid,
)
from .trees import (
    DEFAULT_IQTREE_TASK_THREADS,
    DEFAULT_MAFFT_TASK_THREADS,
    expected_iqtree_tree_dir,
    expected_mafft_output_path,
    valid_iqtree_tree,
    valid_mafft_alignment,
)
from .blastdb import resolve_proteome_profile_input

class AddLibraryTask(Task):
    '''A class that handles the add library task'''
    _BUSCO_STATUS_LABELS = {
        1: "single_copy",
        2: "duplicated",
        3: "fragmented",
        4: "missing",
    }
    _BUSCO_STATUS_MARKERS = {
        1: "BUSCO_SC",
        2: "BUSCO_DUP",
        3: "BUSCO_FG",
        4: "BUSCO_MS",
    }
    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data=data, required_threads=required_threads)
        self.name = self.data.get("name")
        self.library_id = self.data.get("library_id", None)
        self.coverage = self.data.get("coverage", None)
        self.coverage_taxid = self.data.get("coverage_taxid", None)
        self.coverage_label = self.data.get("coverage_label", self.coverage)
        self.accessions = list(dict.fromkeys(normalize_accessions(self.data.get("accessions", []))))
        self.location = self.data.get("location", None)
        self.parent_library_name = self.data.get("parent_library_name")
        self.parent_library_id = self.data.get("parent_library_id", None)
        self.orthofinder_id = self.data.get("orthofinder_id", None)
        self.orthofinder_location = self.data.get("orthofinder_location", None)
        self.rerun_busco = self.data.get("rerun_busco", False)
        self.rerun_orthofinder = self.data.get("rerun_orthofinder", False)
        self.rerun_gene_trees = bool(self.data.get("rerun_gene_trees", False))
        self.skip_paralog_analysis = bool(self.data.get("skip_paralog_analysis", False))
        self.gene_tree_source = str(self.data.get("gene_tree_source") or "iqtree").strip().lower() or "iqtree"
        self.orthofinder_threads = self.data.get("orthofinder_threads", None)
        self.orthofinder_mcl_inflation = self.data.get("orthofinder_mcl_inflation", None)
        self.clean_refs = bool(self.data.get("clean_refs", False))
        self.clean_refs_strict = bool(self.data.get("clean_refs_strict", False))
        self.set_cleaned_primary = bool(self.data.get("set_cleaned_primary", True))
        self.proteome_profile = resolve_profile_selector(
            proteome_profile=self.data.get("proteome_profile"),
            isoforms_cleaned=self.data.get("isoforms_cleaned"),
            raw_proteome=self.data.get("raw_proteome"),
        )
        self.prefer_proteome_profile = str(self.data.get("prefer_proteome_profile") or "").strip() or None
        self.min_species_in_trees = int(self.data.get("min_species_in_trees", 4) or 4)
        self.mafft_threads = int(self.data.get("mafft_threads") or DEFAULT_MAFFT_TASK_THREADS)
        self.iqtree_threads = int(self.data.get("iqtree_threads") or DEFAULT_IQTREE_TASK_THREADS)
        self.mafft_flags = self.data.get("mafft_flags")
        self.iqtree_flags = self.data.get("iqtree_flags")
        self.annotate_og_trees = bool(self.data.get("annotate_og_trees", False))
        self.debug_path = self.data.get("debug_path")
        self.force = bool(self.data.get("force", False))
        self.data["force"] = self.force
        self.stage = checkpoint if checkpoint is not None else 0
        self.log(
            f"AddLibraryTask init name={self.name} coverage={self.coverage_label} parent={self.parent_library_name} stage={self.stage}",
            "DEBUG",
        )

    def _core_set_strategy(self) -> str:
        return "skip_paralog_analysis" if self.skip_paralog_analysis else "paralog_aware"

    def _effective_gene_tree_source(self) -> str:
        if self.skip_paralog_analysis:
            return "none"
        if self.gene_tree_source == "fasttree":
            return "fasttree"
        return "iqtree"

    def _accepted_family_rule(self) -> str:
        if self.skip_paralog_analysis:
            return "exact_1to1_plus_min_species"
        return "exact_1to1_plus_paralog_filter"

    def _rerun_gene_trees_effective(self) -> bool:
        return bool(self.rerun_gene_trees and self._effective_gene_tree_source() == "iqtree")

    def _annotate_og_trees_effective(self) -> bool:
        return bool(self.annotate_og_trees and self._effective_gene_tree_source() != "none")

    # Subtask functions Queue / Done_check / Verify Completion message

    # Phase 1: Ensure assembly metadata present; queue UpdateAssemblyInformation for missing
    def _get_missing_accessions(self):
        rows = self.db_manager.genomes.get_many(self.accessions)
        if not rows:
            present = []
        else:
            present = [r[0] for r in rows]
        missing = [acc for acc in self.accessions if acc not in present]
        return missing

    def queue_metadata(self):
        # Check which accessions are already present in the database.
        # Queue a subtask for missing accession.
        missing = self._get_missing_accessions()
        if missing:
            self.log(f"Queueing metadata subtask for {len(missing)}", "DEBUG")
            payload = {"accessions": missing}
            if self.debug_path:
                payload["debug_path"] = self.debug_path
            self.queue_subtask(job_type=1, status="P", priority=1, data=payload)
            return True
        return False    

    def metadata_done(self):
        # Check which accessions are already present in the database.
        missing = self._get_missing_accessions()
        # present = [acc for acc in self.accessions if acc not in missing]
        # self.log(f"Metadata check complete: {len(missing)} missing assemblies.", "DEBUG")
        # self.log(f"{missing}", "DEBUG")
        # self.log(f"Present assemblies: {len(present)}", "DEBUG")
        # self.log(f"Total assemblies checked: {len(self.accessions)}", "DEBUG")
        self.log(f"Checking if we have downloaded all the missing assemblies... {len(missing)} remain...", "DEBUG")
        for acc in missing:
            self.log(f"Missing accession: {acc}", "DEBUG")
        return len(missing) == 0

    def metadata_incomplete_message(self):
        # No error, and not ongoing but not complete, this function explains why in the error stack
        missing = self._get_missing_accessions()
        summary = f"Missing metadata for {len(missing)} accession(s): {', '.join(missing)}, it is possible these accessions were not found in the database."
        stack = ""
        return (summary, stack)
    
    # Phase 2: Download assemblies using the same orchestration pattern
    def _get_missing_downloads(self):
        rows = self.db_manager.genomes.get_downloaded()
        present = [r[0] for r in rows] if rows else []
        missing = [acc for acc in self.accessions if acc not in present]
        return missing

    def queue_downloads(self):
        missing = self._get_missing_downloads()
        if missing:
            self.queue_subtask(job_type=2, status="P", priority=1, data={"accessions": missing, "protein": True, "max_concurrent": 3})
            return True
        return False

    def downloads_done(self):
        missing = self._get_missing_downloads()
        self.log(f"Download check: {len(missing)} remaining.", "DEBUG")
        return len(missing) == 0

    def downloads_incomplete_message(self):
        missing = self._get_missing_downloads() 
        if not missing:
            return None
        summary = f"Missing BUSCO results for {len(missing)} accession(s): {', '.join(missing)}"
        stack = ""
        return (summary, stack)

    # Phase 3: BUSCO analysis per accession
    # Helpers specific to BUSCO phase

    def _busco_missing(self):
        """Get accessions that are missing BUSCO results."""
        missing = []
        pinned = self._pin_proteome_profile_inputs()
        for accession in self.accessions:
            profile_name = str(pinned[str(accession)]["profile_name"])
            run_id = self.db_manager.busco.get_effective_run_id_for_accession(
                str(accession),
                int(self.parent_library_id),
                pipeline="busco",
                input_mode="protein",
                proteome_profile=profile_name,
                purpose="default",
            )
            if run_id is None:
                missing.append(accession)
        return missing

    def _orthofinder_missing(self):
        """Get true false value as to whether an orthofinder run exists for the list of accessions provided."""
        profile_inputs = self._orthofinder_profile_inputs()
        # Return a tuple (id, location) when both are present; otherwise None so the caller knows it's incomplete.
        rows = self.db_manager.orthofinder.assert_results_exist(
            self.accessions,
            mcl_inflation=self.orthofinder_mcl_inflation,
            profile_inputs=profile_inputs,
        )
        if not rows or len(rows) < 2 or rows[0] is None or rows[1] is None:
            return None
        return rows[0], rows[1]

    def _orthofinder_profile_inputs(self):
        """Resolve the exact profile identity expected for every OrthoFinder input."""
        pinned = self._pin_proteome_profile_inputs()
        return {
            accession: (
                int(values["profile_id"]),
                str(values["checksum"]) if values.get("checksum") is not None else None,
            )
            for accession, values in pinned.items()
        }

    def _pin_proteome_profile_inputs(self):
        """Resolve once and persist the exact input profile used by all child tasks."""
        existing = self.data.get("_resolved_proteome_profile_inputs")
        if isinstance(existing, dict) and set(existing) == set(self.accessions):
            return existing
        pinned = {}
        for accession in self.accessions:
            profile_name, row, path = resolve_proteome_profile_input(
                self.db_manager,
                str(accession),
                proteome_profile=self.proteome_profile,
                prefer_proteome_profile=self.prefer_proteome_profile,
            )
            pinned[str(accession)] = {
                "profile_name": str(profile_name),
                "profile_id": int(row[0]),
                "checksum": str(row[8]) if len(row) > 8 and row[8] is not None else None,
                "path": str(path),
            }
        self.checkpoint(
            int(self.stage or 0),
            {"_resolved_proteome_profile_inputs": pinned},
        )
        self.log(
            "Pinned proteome profiles for add-library: "
            + ", ".join(f"{acc}={values['profile_name']}" for acc, values in sorted(pinned.items())),
            "INFO",
        )
        return pinned

    def queue_busco(self):
        # lineage = self.data.get("busco_lineage", "metazoa_odb10")
        format = "protein"
        # first = True
        missing = self._busco_missing()
        # self.log(f"BUSCO missing: {missing}")
        of_results = self._orthofinder_missing()
        self.log(f"Queueing BUSCO subtask for {len(missing)} missing BUSCO results. Orthofinder results exist: {str(of_results)}", "DEBUG")
        queued = False
        if not of_results:
            orthofinder_root = self.db_manager.storage.get_root_base("orthofinder")
            if not orthofinder_root:
                return False
            orthofinder_payload = {
                "force": self.rerun_orthofinder,
                "accessions": self.accessions,
                "library_id": self.library_id,
                "out_dir": orthofinder_root,
                "proteome_profile_inputs": self._pin_proteome_profile_inputs(),
            }
            if self.orthofinder_mcl_inflation not in (None, ""):
                orthofinder_payload["mcl_inflation"] = float(self.orthofinder_mcl_inflation)
            if self.proteome_profile:
                orthofinder_payload["proteome_profile"] = self.proteome_profile
            if self.prefer_proteome_profile:
                orthofinder_payload["prefer_proteome_profile"] = self.prefer_proteome_profile
            if self.orthofinder_threads not in (None, ""):
                orthofinder_payload["required_threads"] = int(self.orthofinder_threads)
            self.queue_subtask(
                job_type=5, status='P', priority=1,
                data=orthofinder_payload,
            )
            self.log(f"Queued OrthoFinder subtask for accessions: {', '.join(self.accessions)}", "DEBUG")
            queued = True
        for acc in missing:
            pinned_input = self._pin_proteome_profile_inputs()[str(acc)]
            resolved_profile = pinned_input["profile_name"]
            busco_payload = {
                "accession": acc,
                "lineage": self.parent_library_name,
                "format": format,
                "force": self.rerun_busco,
                "proteome_profile": resolved_profile,
                "expected_proteome_profile_id": int(pinned_input["profile_id"]),
            }
            if pinned_input.get("checksum") is not None:
                busco_payload["expected_proteome_checksum"] = str(pinned_input["checksum"])
            if self.prefer_proteome_profile:
                busco_payload["prefer_proteome_profile"] = self.prefer_proteome_profile
            self.queue_subtask(
                job_type=4, status="P", priority=1,
                data=busco_payload,
            )
            # if first:
            #     #Give first task a headstart to download lineage if not available...
            #     time.sleep(10)
            #     first = False
            queued = True
        if not missing and not of_results and not queued:
            # Nothing missing and orthofinder already satisfied; nothing to queue.
            return False
        return queued

    def busco_done(self):
        # Check that all accessions have BUSCO results
        missing = self._busco_missing()
        orthofinder_results = self._orthofinder_missing()

        # Safely assign defaults then override if results present to avoid UnboundLocalError
        if orthofinder_results:
            orthofinder_id, orthofinder_location = orthofinder_results
            # Consider orthofinder done only if both id and location are present/non-empty
            orthofinder_done = bool(orthofinder_id and orthofinder_location)
        else:
            orthofinder_id = None
            orthofinder_location = None
            orthofinder_done = False

        self.log(f"BUSCO check: {len(missing)} remaining. Orthofinder results available: {str(orthofinder_id)}", "DEBUG")
        return (len(missing) == 0) and orthofinder_done

    def busco_incomplete_message(self):
        missing = self._busco_missing()
        of_results = self._orthofinder_missing()
        summary = f"BUSCO phase incomplete. Missing BUSCO results for {len(missing)} accession(s): {', '.join(missing)}"
        if not of_results:
            summary += f" Orthofinder results are also missing."
        stack = ""
        return (summary, stack)

    def handle_exception(self, exc, context = None):
        # remove the library from the database if it was added during this task
        self.db_manager.libraries.purge(self.library_id)
        return super().handle_exception(exc, context)

    def _purge_existing_library_state(self, *, library_id: int, location: Optional[str]) -> bool:
        preserve_orthofinder = not bool(self.rerun_orthofinder)
        if not self.db_manager.libraries.purge(int(library_id), preserve_orthofinder=preserve_orthofinder):
            return False
        if location and os.path.isdir(location):
            try:
                shutil.rmtree(location)
            except OSError as exc:
                self.log(f"Failed to remove existing library directory before rebuild: {exc}", "ERROR")
                return False
        return True

    def _clean_refs_mode(self) -> Optional[str]:
        if bool(self.clean_refs):
            return "out_only"
        if bool(self.clean_refs_strict):
            return "strict_any"
        return None

    def _rewrite_family_status(self, source_status: int, *, has_in: bool, has_out: bool, mode: str) -> int:
        try:
            status_val = int(source_status)
        except (TypeError, ValueError):
            return source_status
        if status_val != 1:
            return status_val
        if mode == "out_only" and has_out:
            return 2
        if mode == "strict_any" and (has_in or has_out):
            return 2
        return status_val

    def _family_status_map_from_rows(self, family_rows):
        status_map: Dict[str, int] = {}
        for family_id, _library_id, _accession, status, *_rest in family_rows or []:
            fam = str(family_id)
            try:
                status_val = int(status)
            except (TypeError, ValueError):
                continue
            previous = status_map.get(fam)
            if previous is None or status_val > previous:
                status_map[fam] = status_val
        return status_map

    def _status_label(self, status: object) -> str:
        try:
            return self._BUSCO_STATUS_LABELS.get(int(status), str(status))
        except (TypeError, ValueError):
            return str(status or "NA")

    def _status_marker(self, status: object) -> str:
        try:
            return self._BUSCO_STATUS_MARKERS.get(int(status), "BUSCO_NA")
        except (TypeError, ValueError):
            return "BUSCO_NA"

    def _source_gene_name_map(self, family_rows, accession: str) -> Dict[str, str]:
        accession_token = canonicalize_accession(accession)
        genes: Dict[str, list[str]] = {}
        for family_id, _library_id, row_accession, _status, sequence, *_rest in family_rows or []:
            if canonicalize_accession(row_accession) != accession_token:
                continue
            family_token = str(family_id)
            gene_name = str(sequence or "").strip()
            if not gene_name:
                continue
            genes.setdefault(family_token, [])
            if gene_name not in genes[family_token]:
                genes[family_token].append(gene_name)
        return {
            family_id: ",".join(values) if values else "NA"
            for family_id, values in genes.items()
        }

    def _resolve_paralog_evidence_file(self, orthogroup: str, paralog_class: str, compare_results: Dict[str, object]) -> str:
        if not orthogroup:
            return "NA"
        files = compare_results.get("files", {}) or {}
        species_tsv = str(files.get("species_paralog_tsv") or "")
        paralog_dir = os.path.join(os.path.dirname(species_tsv), "paralogs") if species_tsv else os.path.join(self.location or "", "core_set_analysis", "paralogs")
        if paralog_class == "In-paralog":
            path = os.path.join(paralog_dir, f"{orthogroup}_inparalogs.txt")
            return path if os.path.exists(path) else "NA"
        if paralog_class == "Out-paralog":
            path = os.path.join(paralog_dir, f"{orthogroup}_outparalogs.txt")
            return path if os.path.exists(path) else "NA"
        return "NA"

    def _resolve_gene_tree_file(self, orthogroup: str) -> str:
        if not orthogroup or not self.orthofinder_location:
            return "NA"
        resolved_dir = self._effective_tree_dir()
        if not resolved_dir:
            return "NA"
        candidates = [
            os.path.join(resolved_dir, f"{orthogroup}_tree.txt"),
            os.path.join(resolved_dir, f"{orthogroup}.txt"),
            os.path.join(resolved_dir, f"{orthogroup}.nex"),
        ]
        for path in candidates:
            if os.path.exists(path):
                return path
        consolidated = os.path.join(resolved_dir, "Resolved_Gene_Trees.txt")
        return consolidated if os.path.exists(consolidated) else "NA"

    def _write_cleaned_reference_report(
        self,
        accession: str,
        source_status_map: Dict[str, int],
        source_gene_map: Dict[str, str],
        rewritten_status_map: Dict[str, int],
        compare_results: Dict[str, object],
    ) -> bool:
        if not self.location:
            return False
        report_dir = os.path.join(self.location, "cleaned_reference_proteomes")
        os.makedirs(report_dir, exist_ok=True)
        report_path = os.path.join(report_dir, f"{accession}.tsv")

        family_to_orthogroup = compare_results.get("family_to_orthogroup_exact", {}) or {}
        family_to_orthogroups_all = compare_results.get("family_to_orthogroups_all", {}) or {}
        family_species_status = compare_results.get("family_species_paralog_status", {}) or {}
        unmapped_summary = compare_results.get("unmapped_busco_summary", {}) or {}
        per_species_unmapped = unmapped_summary.get("per_species_unmapped_counts", {}) or {}
        accession_token = canonicalize_accession(accession)
        families = sorted(set(source_status_map) | set(rewritten_status_map))
        try:
            with open(report_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["# accession", accession])
                writer.writerow(["# unmapped_buscos", int(per_species_unmapped.get(accession_token, 0) or 0)])
                writer.writerow(
                    [
                        "accession",
                        "family_name",
                        "orthogroup",
                        "orthogroup_mapping",
                        "original_status",
                        "new_status",
                        "gene_name",
                        "paralog",
                        "paralog_file",
                        "gene_tree_file",
                    ]
                )
                for family_id in families:
                    family_key = str(family_id)
                    exact_orthogroup = str(family_to_orthogroup.get(family_key) or "")
                    all_orthogroups = [str(item) for item in (family_to_orthogroups_all.get(family_key) or []) if str(item)]
                    if exact_orthogroup:
                        orthogroup = exact_orthogroup
                        mapping_status = "Exact_1to1"
                    elif len(all_orthogroups) == 1:
                        orthogroup = all_orthogroups[0]
                        mapping_status = "Mapped_nonexact"
                    elif len(all_orthogroups) > 1:
                        orthogroup = ",".join(all_orthogroups)
                        mapping_status = "Mapped_multiple_OGs"
                    else:
                        orthogroup = "NA"
                        mapping_status = "Unmapped"
                    species_status = (family_species_status.get(str(family_id), {}) or {}).get(accession_token, {})
                    has_in = bool(species_status.get("has_inparalogs"))
                    has_out = bool(species_status.get("has_outparalogs"))
                    if has_out:
                        paralog_class = "Out-paralog"
                    elif has_in:
                        paralog_class = "In-paralog"
                    elif orthogroup != "NA":
                        paralog_class = "Ortholog"
                    else:
                        paralog_class = "NA"
                    writer.writerow(
                        [
                            accession,
                            family_id,
                            orthogroup,
                            mapping_status,
                            self._status_label(source_status_map.get(str(family_id), "NA")),
                            self._status_label(rewritten_status_map.get(str(family_id), "NA")),
                            str(source_gene_map.get(str(family_id)) or "NA"),
                            paralog_class,
                            self._resolve_paralog_evidence_file(orthogroup if orthogroup != "NA" else "", paralog_class, compare_results),
                            str(self._resolve_gene_tree_file(orthogroup if orthogroup != "NA" else "")),
                        ]
                    )
            return True
        except OSError as exc:
            self.log(f"Failed to write cleaned reference report for {accession}: {exc}", "WARNING")
            return False

    def _resolve_source_busco_run_id(self, accession: str) -> Optional[int]:
        return self.db_manager.busco.get_effective_run_id_for_accession(
            str(accession),
            int(self.parent_library_id),
            proteome_profile=self.proteome_profile,
            preferred_proteome_profile=self.prefer_proteome_profile,
            purpose="default",
        )

    def _resolve_original_source_busco_run_id(self, accession: str) -> Optional[int]:
        runs = self.db_manager.busco.get_runs_for_accessions(
            [str(accession)],
            library_id=int(self.parent_library_id),
            purpose="default",
        )
        if not runs:
            return None

        profile_wanted = self.db_manager.busco._normalize_proteome_profile(self.proteome_profile)
        preferred_profile = self.db_manager.busco._normalize_proteome_profile(self.prefer_proteome_profile)

        candidates = []
        for row in runs:
            row_mode = str(row[5] or "").strip().lower()
            row_pipeline = str(row[6] or "").strip().lower()
            row_profile = self.db_manager.busco._normalize_proteome_profile(row[20] if len(row) > 20 else None)
            if row_mode != "protein":
                continue
            if row_pipeline == "orthofinder":
                continue
            profile_score = 0
            if profile_wanted and row_profile == profile_wanted:
                profile_score = 2
            elif preferred_profile and row_profile == preferred_profile:
                profile_score = 1
            candidates.append((profile_score, row))

        if not candidates:
            return self._resolve_source_busco_run_id(accession)

        candidates.sort(
            key=lambda item: (
                item[0],
                str(item[1][13] or ""),
                int(item[1][0] or 0),
            ),
            reverse=True,
        )
        best = candidates[0][1]
        try:
            return int(best[0])
        except (TypeError, ValueError):
            return None

    def _prepare_exact_busco_orthogroups(self, analyser: OrthoBuscoAnalyzer) -> tuple[list[tuple[str, str]], list[tuple[str, str]], Dict[str, object]]:
        out_prefix = os.path.join(analyser.working_dir, f"{analyser.identifier}_busco_to_orthogroup")
        mapping_tsv = out_prefix + "_map.tsv"
        one_to_one_tsv = out_prefix + "_1to1.tsv"
        og_to_busco_tsv = out_prefix + "_og_to_busco_families.tsv"
        busco_to_og_tsv = out_prefix + "_busco_family_to_ogs.tsv"
        exact_mapping_tsv = out_prefix + "_exact_map.tsv"
        seq_index, family_to_seq_counts = analyser._index_busco_sequences(use_hash=True)
        if seq_index is None or family_to_seq_counts is None:
            raise ValueError("Failed to index BUSCO sequences for add-library analysis.")
        (
            og_to_busco_families,
            busco_family_to_ogs,
            og_family_seq_counts,
            og_total_seq_count,
            og_family_species,
            mapped_busco_ids,
        ) = analyser._scan_orthogroups(
            seq_index=seq_index,
            use_hash=True,
            mapping_tsv=mapping_tsv,
        )
        if og_to_busco_families is None:
            raise ValueError("Failed to map orthogroups to BUSCO families.")
        analyser.log(f"Mapping completed. Outputs: {mapping_tsv}, {one_to_one_tsv}")
        unmapped_summary = analyser._summarize_unmapped_buscos(mapped_busco_ids, out_prefix)
        with open(og_to_busco_tsv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Orthogroup", "BUSCO_family"])
            for orthogroup in sorted(og_to_busco_families.keys()):
                for family_id in sorted(og_to_busco_families[orthogroup]):
                    writer.writerow([orthogroup, family_id])
        with open(busco_to_og_tsv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["BUSCO_family", "Orthogroup"])
            for family_id in sorted(busco_family_to_ogs.keys()):
                for orthogroup in sorted(busco_family_to_ogs[family_id]):
                    writer.writerow([family_id, orthogroup])
        analyser.log(f"Mapping pair TSVs written: {og_to_busco_tsv}, {busco_to_og_tsv}")
        exact_pairs, _relaxed = analyser._write_1to1_report(
            one_to_one_tsv=one_to_one_tsv,
            og_to_busco_families=og_to_busco_families,
            busco_family_to_ogs=busco_family_to_ogs,
            og_family_seq_counts=og_family_seq_counts,
            og_total_seq_count=og_total_seq_count,
            family_to_seq_counts=family_to_seq_counts,
            require_no_extra_sequences=False,
        )
        analyser.log("Summary: Orthogroup -> BUSCO families map (complete list)")
        for orthogroup in sorted(og_to_busco_families.keys()):
            family_ids = sorted(og_to_busco_families[orthogroup])
            analyser.log(f"  {orthogroup} -> {', '.join(family_ids) if family_ids else '(none)'}")

        multi_busco_ogs = [orthogroup for orthogroup, family_ids in og_to_busco_families.items() if len(family_ids) > 1]
        analyser.log(f"Orthogroups containing multiple BUSCO families: {len(multi_busco_ogs)}")
        for orthogroup in sorted(multi_busco_ogs):
            analyser.log(f"  MULTI-BUSCO OG {orthogroup}: {', '.join(sorted(og_to_busco_families[orthogroup]))}")

        multi_og_buscos = [family_id for family_id, orthogroups in busco_family_to_ogs.items() if len(orthogroups) > 1]
        analyser.log(f"BUSCO families mapped to multiple Orthogroups: {len(multi_og_buscos)}")
        for family_id in sorted(multi_og_buscos):
            analyser.log(f"  MULTI-OG BUSCO {family_id}: {', '.join(sorted(busco_family_to_ogs[family_id]))}")
        occupancy_tsv = out_prefix + "_1to1_occupancy.tsv"
        exact_pairs_filtered: list[tuple[str, str]] = []
        below_threshold_families: list[tuple[str, str, int]] = []
        with open(occupancy_tsv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["BUSCO_family", "Orthogroup", "Species_count", "Status"])
            for family_id, orthogroup in exact_pairs:
                species_count = len(og_family_species.get(orthogroup, {}).get(family_id, set()))
                accepted = species_count >= self.min_species_in_trees
                writer.writerow([family_id, orthogroup, species_count, "accepted" if accepted else "rejected"])
                if accepted:
                    exact_pairs_filtered.append((family_id, orthogroup))
                else:
                    below_threshold_families.append((family_id, orthogroup, species_count))
        analyser.log(f"Exact 1-1 BUSCO families total: {len(exact_pairs)}")
        analyser.log(f"Exact 1-1 families with >= min_species ({self.min_species_in_trees}): {len(exact_pairs_filtered)}")
        analyser.log(f"Exact 1-1 families rejected due to low occupancy: {len(below_threshold_families)}")
        analyser.log(f"1-1 occupancy TSV written: {occupancy_tsv}")
        if below_threshold_families:
            analyser.log("BUSCO families rejected due to low occupancy:")
            for family_id, orthogroup, species_count in sorted(below_threshold_families, key=lambda item: (item[2], item[0])):
                analyser.log(f"  {family_id} in {orthogroup}: species={species_count}")
        analyser._write_exact_mapping_tsv(mapping_tsv, exact_mapping_tsv, exact_pairs)
        analyser.log(f"Exact BUSCO mapping written: {exact_mapping_tsv}")
        compare_results: Dict[str, object] = {
            "files": {
                "mapping_tsv": mapping_tsv,
                "one_to_one_tsv": one_to_one_tsv,
                "occupancy_tsv": occupancy_tsv,
                "og_to_busco_tsv": og_to_busco_tsv,
                "busco_to_og_tsv": busco_to_og_tsv,
                "exact_mapping_tsv": exact_mapping_tsv,
                "unmapped_busco_summary_tsv": unmapped_summary.get("summary_tsv", ""),
            },
            "family_to_orthogroup_exact": {str(family_id): str(orthogroup) for family_id, orthogroup in exact_pairs_filtered},
            "family_to_orthogroups_all": {
                str(family_id): sorted(str(orthogroup) for orthogroup in orthogroups)
                for family_id, orthogroups in busco_family_to_ogs.items()
            },
            "unmapped_busco_summary": unmapped_summary,
        }
        return exact_pairs, exact_pairs_filtered, compare_results

    def _tree_workspace_paths(self) -> dict[str, str]:
        if not self.orthofinder_location:
            raise ValueError("OrthoFinder results location is missing.")
        location = str(self.orthofinder_location)
        return {
            "raw_dir": os.path.join(location, "Orthogroup_Sequences"),
            # Keep replacement MAFFT outputs separate from OrthoFinder's own alignments.
            "align_dir": os.path.join(location, "MAFFT_Alignments"),
            "fasttree_dir": os.path.join(location, "Resolved_Gene_Trees"),
            "iqtree_dir": os.path.join(location, "IQ-TREE_Orthogroup_trees"),
            "metadata_dir": os.path.join(location, "IQTREE_Metadata"),
        }

    def _existing_iqtree_tree_path(self, orthogroup: str) -> Optional[str]:
        if not orthogroup or not self.orthofinder_location:
            return None
        resolved_dir = os.path.join(str(self.orthofinder_location), "IQ-TREE_Orthogroup_trees")
        candidates = [
            os.path.join(resolved_dir, f"{orthogroup}_tree.txt"),
            os.path.join(resolved_dir, f"{orthogroup}.txt"),
            os.path.join(resolved_dir, f"{orthogroup}.nex"),
        ]
        for path in candidates:
            if valid_iqtree_tree(path):
                return path
        return None

    def _replacement_tree_rows(self, exact_pairs_filtered: list[tuple[str, str]]) -> list[dict[str, str]]:
        paths = self._tree_workspace_paths()
        rows: list[dict[str, str]] = []
        for family_id, orthogroup in exact_pairs_filtered:
            raw_fasta = os.path.join(paths["raw_dir"], f"{orthogroup}.fa")
            if not os.path.exists(raw_fasta):
                alt = os.path.join(paths["raw_dir"], f"{orthogroup}.faa")
                if os.path.exists(alt):
                    raw_fasta = alt
            if not os.path.exists(raw_fasta):
                raise FileNotFoundError(f"Orthogroup FASTA not found for {orthogroup}.")
            alignment_path = expected_mafft_output_path(
                input_fasta=raw_fasta,
                out_dir=paths["align_dir"],
                output_name=f"{orthogroup}.fa",
            )
            tree_dir, prefix = expected_iqtree_tree_dir(
                input_alignment=alignment_path,
                out_dir=paths["metadata_dir"],
                prefix=orthogroup,
            )
            rows.append(
                {
                    "family_id": str(family_id),
                    "orthogroup": str(orthogroup),
                    "raw_fasta": raw_fasta,
                    "alignment_path": alignment_path,
                    "tree_dir": tree_dir,
                    "prefix": prefix,
                    "canonical_tree": os.path.join(paths["iqtree_dir"], f"{orthogroup}_tree.txt"),
                    "fasttree_tree_path": os.path.join(paths["fasttree_dir"], f"{orthogroup}_tree.txt"),
                }
            )
        return rows

    def _prepare_replacement_tree_workspace(self) -> None:
        paths = self._tree_workspace_paths()
        if self._effective_gene_tree_source() != "iqtree":
            return
        os.makedirs(paths["align_dir"], exist_ok=True)
        os.makedirs(paths["iqtree_dir"], exist_ok=True)
        os.makedirs(paths["metadata_dir"], exist_ok=True)

    def _existing_fasttree_tree_path(self, orthogroup: str) -> Optional[str]:
        if not orthogroup or not self.orthofinder_location:
            return None
        fasttree_dir = os.path.join(str(self.orthofinder_location), "Resolved_Gene_Trees")
        candidates = [
            os.path.join(fasttree_dir, f"{orthogroup}_tree.txt"),
            os.path.join(fasttree_dir, f"{orthogroup}.txt"),
            os.path.join(fasttree_dir, f"{orthogroup}.tree"),
        ]
        for path in candidates:
            if os.path.exists(path):
                return path
        return None

    def _iqtree_tree_satisfies_orthogroup(self, orthogroup: str) -> bool:
        return bool(not self._rerun_gene_trees_effective() and self._existing_iqtree_tree_path(orthogroup))

    def _queue_replacement_mafft_subtasks(self, exact_pairs_filtered: list[tuple[str, str]]) -> bool:
        if self._effective_gene_tree_source() != "iqtree":
            return False
        queued = False
        for row in self._replacement_tree_rows(exact_pairs_filtered):
            if self._iqtree_tree_satisfies_orthogroup(row["orthogroup"]):
                continue
            if valid_mafft_alignment(row["alignment_path"], row["raw_fasta"]):
                continue
            self.queue_subtask(
                job_type=32,
                status="P",
                priority=1,
                data={
                    "input_fasta": row["raw_fasta"],
                    "out_dir": os.path.dirname(row["alignment_path"]),
                    "output_name": os.path.basename(row["alignment_path"]),
                    "mafft_flags": self.mafft_flags,
                    "required_threads": int(self.mafft_threads),
                },
            )
            queued = True
        return queued

    def _replacement_mafft_done(self, exact_pairs_filtered: list[tuple[str, str]]) -> bool:
        if self._effective_gene_tree_source() != "iqtree":
            return True
        rows = self._replacement_tree_rows(exact_pairs_filtered)
        return bool(rows) and all(
            self._iqtree_tree_satisfies_orthogroup(row["orthogroup"])
            or valid_mafft_alignment(row["alignment_path"], row["raw_fasta"])
            for row in rows
        )

    def _replacement_mafft_incomplete_message(self, exact_pairs_filtered: list[tuple[str, str]]) -> str:
        rows = self._replacement_tree_rows(exact_pairs_filtered)
        failures = []
        for row in rows:
            if self._iqtree_tree_satisfies_orthogroup(row["orthogroup"]):
                continue
            if valid_mafft_alignment(row["alignment_path"], row["raw_fasta"]):
                continue
            state = "invalid" if os.path.exists(row["alignment_path"]) else "missing"
            failures.append(f"{row['orthogroup']} ({state}: {row['alignment_path']})")
        if not rows:
            return "MAFFT phase incomplete: no replacement orthogroups were available."
        preview = "; ".join(failures[:10])
        if len(failures) > 10:
            preview += f"; ... and {len(failures) - 10} more"
        return f"MAFFT phase incomplete: {len(failures)} alignment output(s) missing or invalid: {preview}"

    def _queue_replacement_iqtree_subtasks(self, exact_pairs_filtered: list[tuple[str, str]]) -> bool:
        if self._effective_gene_tree_source() != "iqtree":
            return False
        queued = False
        for row in self._replacement_tree_rows(exact_pairs_filtered):
            if self._iqtree_tree_satisfies_orthogroup(row["orthogroup"]):
                continue
            best_tree = self._find_iqtree_tree(row["tree_dir"], row["prefix"])
            if best_tree:
                if self._rerun_gene_trees_effective():
                    pass
                else:
                    continue
            self.queue_subtask(
                job_type=33,
                status="P",
                priority=1,
                data={
                    "input_alignment": row["alignment_path"],
                    "out_dir": os.path.dirname(row["tree_dir"]),
                    "prefix": row["orthogroup"],
                    "iqtree_flags": self.iqtree_flags,
                    "force_restart": bool(
                        self._rerun_gene_trees_effective()
                        or getattr(self, "_phase_meta", {}).get("gen", 1) > 1
                    ),
                    "required_threads": int(self.iqtree_threads),
                },
            )
            queued = True
        return queued

    def _replacement_iqtree_done(self, exact_pairs_filtered: list[tuple[str, str]]) -> bool:
        if self._effective_gene_tree_source() != "iqtree":
            return True
        rows = self._replacement_tree_rows(exact_pairs_filtered)
        if not rows:
            return False
        for row in rows:
            if self._iqtree_tree_satisfies_orthogroup(row["orthogroup"]):
                continue
            best_tree = self._find_iqtree_tree(row["tree_dir"], row["prefix"])
            if not best_tree:
                return False
        return True

    def _replacement_iqtree_incomplete_message(self, exact_pairs_filtered: list[tuple[str, str]]) -> str:
        rows = self._replacement_tree_rows(exact_pairs_filtered)
        failures = []
        for row in rows:
            if self._iqtree_tree_satisfies_orthogroup(row["orthogroup"]):
                continue
            if self._find_iqtree_tree(row["tree_dir"], row["prefix"]):
                continue
            candidates = [
                os.path.join(row["tree_dir"], f"{row['prefix']}.treefile"),
                os.path.join(row["tree_dir"], f"{row['prefix']}.contree"),
            ]
            state = "invalid" if any(os.path.exists(path) for path in candidates) else "missing"
            failures.append(f"{row['orthogroup']} ({state}: {row['tree_dir']})")
        if not rows:
            return "IQ-TREE phase incomplete: no replacement orthogroups were available."
        preview = "; ".join(failures[:10])
        if len(failures) > 10:
            preview += f"; ... and {len(failures) - 10} more"
        return f"IQ-TREE phase incomplete: {len(failures)} tree output(s) missing or invalid: {preview}"

    def _find_iqtree_tree(self, tree_dir: str, prefix: str) -> Optional[str]:
        candidates = [
            os.path.join(tree_dir, f"{prefix}.treefile"),
            os.path.join(tree_dir, f"{prefix}.contree"),
        ]
        for path in candidates:
            if valid_iqtree_tree(path):
                return path
        return None

    def _find_fasttree_tree(self, orthogroup: str) -> Optional[str]:
        return self._existing_fasttree_tree_path(orthogroup)

    def _effective_tree_dir(self) -> str:
        if not self.orthofinder_location:
            return ""
        paths = self._tree_workspace_paths()
        if self._effective_gene_tree_source() == "iqtree":
            return paths["iqtree_dir"]
        if self._effective_gene_tree_source() == "none":
            return ""
        return paths["fasttree_dir"]

    def _annotation_output_dir(self) -> str:
        return os.path.join(str(self.location), "annotated-og-trees")

    def _annotation_input_dir(self) -> str:
        return self._effective_tree_dir()

    def _queue_annotation_subtask(self) -> bool:
        if not self._annotate_og_trees_effective():
            return False
        input_dir = self._annotation_input_dir()
        manifest_path = os.path.join(str(self.location), "orthogroup_tree_manifest.tsv")
        if not input_dir or not os.path.isdir(input_dir):
            raise FileNotFoundError("Gene tree directory is missing for annotation.")
        if not os.path.exists(manifest_path):
            raise FileNotFoundError("Orthogroup tree manifest is missing for annotation.")
        if self._annotation_done():
            return False
        self.queue_subtask(
            job_type=35,
            status="P",
            priority=1,
            data={
                "input_dir": input_dir,
                "out_dir": self._annotation_output_dir(),
                "manifest_tsv": manifest_path,
                "orthofinder_location": self.orthofinder_location,
                "target_library_id": self.library_id,
                "busco_library_id": self.parent_library_id,
            },
        )
        self.log("Queued annotate-orthogroup-tree subtask for replacement gene trees.", "INFO")
        return True

    def _annotation_done(self) -> bool:
        if not self._annotate_og_trees_effective():
            return True
        input_dir = self._annotation_input_dir()
        output_dir = self._annotation_output_dir()
        if not input_dir or not os.path.isdir(input_dir) or not os.path.isdir(output_dir):
            return False
        tree_bases = sorted(
            (
                base[:-9] if base.endswith("_tree.txt") else os.path.splitext(base)[0]
                for base in os.listdir(input_dir)
                if base != "Resolved_Gene_Trees.txt"
                and (base.endswith("_tree.txt") or base.endswith(".treefile"))
            )
        )
        if not tree_bases:
            return False
        return all(os.path.exists(os.path.join(output_dir, f"{orthogroup}.nex")) for orthogroup in tree_bases)

    def _materialize_replacement_orthogroup_trees(self, exact_pairs_filtered: list[tuple[str, str]], compare_results: Dict[str, object]) -> str:
        self._prepare_replacement_tree_workspace()
        rows = self._replacement_tree_rows(exact_pairs_filtered)
        paths = self._tree_workspace_paths()
        manifest_path = os.path.join(self.location, "orthogroup_tree_manifest.tsv")
        source_run_ids = sorted(
            {
                int(run_id)
                for run_id in (self._resolve_original_source_busco_run_id(accession) for accession in self.accessions)
                if run_id is not None
            }
        )
        analyser_working_dir = os.path.join(self.location, "core_set_analysis", "paralogs")
        with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "family_id",
                    "orthogroup",
                    "raw_fasta",
                    "alignment_path",
                    "tree_dir",
                    "tree_path",
                    "fasttree_tree_path",
                    "paralog_in_file",
                    "paralog_out_file",
                    "source_run_ids",
                ]
            )
            for row in rows:
                if self._effective_gene_tree_source() == "fasttree":
                    source_tree = self._find_fasttree_tree(row["orthogroup"])
                    if not source_tree:
                        raise FileNotFoundError(f"FastTree output missing for {row['orthogroup']}.")
                    tree_path_value = source_tree
                else:
                    existing_iqtree = None if self._rerun_gene_trees_effective() else self._existing_iqtree_tree_path(row["orthogroup"])
                    best_tree = self._find_iqtree_tree(row["tree_dir"], row["prefix"])
                    source_tree = existing_iqtree or best_tree
                    if not source_tree:
                        raise FileNotFoundError(f"IQ-TREE output missing for {row['orthogroup']}.")
                    tree_path_value = row["canonical_tree"]
                    if os.path.abspath(source_tree) != os.path.abspath(row["canonical_tree"]):
                        shutil.copyfile(source_tree, row["canonical_tree"])
                writer.writerow(
                    [
                        row["family_id"],
                        row["orthogroup"],
                        row["raw_fasta"],
                        row["alignment_path"] if self._effective_gene_tree_source() == "iqtree" else "",
                        row["tree_dir"] if self._effective_gene_tree_source() == "iqtree" else "",
                        tree_path_value,
                        row["fasttree_tree_path"] if self._effective_gene_tree_source() == "fasttree" else "",
                        os.path.join(analyser_working_dir, f"{row['orthogroup']}_inparalogs.txt"),
                        os.path.join(analyser_working_dir, f"{row['orthogroup']}_outparalogs.txt"),
                        ",".join(str(run_id) for run_id in source_run_ids),
                    ]
                )
        if self._effective_gene_tree_source() == "iqtree":
            consolidated = os.path.join(paths["iqtree_dir"], "Resolved_Gene_Trees.txt")
            with open(consolidated, "w", encoding="utf-8") as out_handle:
                for row in rows:
                    if not os.path.exists(row["canonical_tree"]):
                        continue
                    with open(row["canonical_tree"], "r", encoding="utf-8") as tree_handle:
                        out_handle.write(f"{row['orthogroup']}:{tree_handle.read().strip()}\n")
        compare_results.setdefault("files", {})
        compare_results["files"]["orthogroup_tree_manifest"] = manifest_path
        compare_results["files"]["effective_gene_tree_dir"] = self._effective_tree_dir()
        if self._effective_gene_tree_source() == "fasttree" and os.path.isdir(paths["fasttree_dir"]):
            compare_results["files"]["fasttree_gene_tree_dir"] = paths["fasttree_dir"]
        if self._effective_gene_tree_source() == "iqtree" and os.path.isdir(paths["iqtree_dir"]):
            compare_results["files"]["iqtree_gene_tree_dir"] = paths["iqtree_dir"]
        return manifest_path

    def _build_replacement_orthogroup_trees(self, exact_pairs_filtered: list[tuple[str, str]], compare_results: Dict[str, object]) -> str:
        if not self.orthofinder_location:
            raise ValueError("OrthoFinder results location is missing.")
        return self._materialize_replacement_orthogroup_trees(exact_pairs_filtered, compare_results)

    def _finalize_replacement_tree_analysis(
        self,
        analyser: OrthoBuscoAnalyzer,
        exact_pairs: list[tuple[str, str]],
        exact_pairs_filtered: list[tuple[str, str]],
        compare_results: Dict[str, object],
    ) -> list[str]:
        matched_exact_ogs = {orthogroup for _family_id, orthogroup in exact_pairs}
        matched_exact_ogs_filtered = {orthogroup for _family_id, orthogroup in exact_pairs_filtered}
        paralog_results = analyser._identify_paralogs(
            matched_exact_ogs=matched_exact_ogs,
            min_species=self.min_species_in_trees,
            write_tsv=True,
        )
        if paralog_results is None:
            raise ValueError("Failed to compute paralog summaries from replacement gene trees.")
        og_has_out = set(paralog_results.get("og_has_outparalogs", []) or [])
        good_ogs = {orthogroup for orthogroup in matched_exact_ogs_filtered if orthogroup not in og_has_out}
        cleaned_busco_families = sorted(
            str(family_id)
            for family_id, orthogroup in exact_pairs_filtered
            if orthogroup in good_ogs
        )
        final_list_path = os.path.join(analyser.working_dir, f"{analyser.identifier}_good_busco_families.txt")
        with open(final_list_path, "w", encoding="utf-8") as handle:
            for family_id in cleaned_busco_families:
                handle.write(f"{family_id}\n")
        analyser.log(f"Wrote final list of good BUSCO families ({len(cleaned_busco_families)}): {final_list_path}")
        og_species_in = paralog_results.get("og_species_inparalogs", {}) or {}
        og_species_out = paralog_results.get("og_species_outparalogs", {}) or {}
        family_species_status: Dict[str, Dict[str, Dict[str, bool]]] = {}
        out_prefix = os.path.join(analyser.working_dir, f"{analyser.identifier}_busco_to_orthogroup")
        mapping_tsv = str((compare_results.get("files", {}) or {}).get("mapping_tsv") or "")
        species_tsv = out_prefix + "_species_paralog_status.tsv"
        with open(species_tsv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["BUSCO_family", "Orthogroup", "Accession", "Has_In_Paralogs", "Has_Out_Paralogs"])
            for family_id, orthogroup in exact_pairs_filtered:
                family_species_status.setdefault(str(family_id), {})
                species_rows: Dict[str, Dict[str, bool]] = {}
                for accession_token in sorted(set(og_species_in.get(orthogroup, set())) | set(og_species_out.get(orthogroup, set()))):
                    has_in = accession_token in set(og_species_in.get(orthogroup, set()) or set())
                    has_out = accession_token in set(og_species_out.get(orthogroup, set()) or set())
                    species_rows[str(accession_token)] = {
                        "has_inparalogs": has_in,
                        "has_outparalogs": has_out,
                    }
                    family_species_status[str(family_id)][accession_token] = {
                        "has_inparalogs": has_in,
                        "has_outparalogs": has_out,
                    }
                for accession_token in sorted(species_rows):
                    row = species_rows[accession_token]
                    writer.writerow([
                        family_id,
                        orthogroup,
                        accession_token,
                        "YES" if row.get("has_inparalogs") else "NO",
                        "YES" if row.get("has_outparalogs") else "NO",
                    ])
        augmented_mapping_tsv = out_prefix + "_map_with_paralog_class.tsv"
        paralog_classification = paralog_results.get("og_paralog_classification", {}) or {}
        clean_pairs = {
            (str(family_id), str(orthogroup))
            for family_id, orthogroup in exact_pairs_filtered
            if paralog_classification.get(orthogroup) in ("No Paralogs", "Only In-Paralogs")
        }
        if mapping_tsv and os.path.exists(mapping_tsv):
            with open(mapping_tsv, "r", newline="", encoding="utf-8") as in_f, open(augmented_mapping_tsv, "w", newline="", encoding="utf-8") as out_f:
                reader = csv.DictReader(in_f, delimiter="\t")
                base_fields = list(reader.fieldnames) if reader.fieldnames else []
                fieldnames = base_fields + ["OG_Paralog_Classification", "Clean_1to1_NoParalogs"]
                writer = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                for row in reader:
                    orthogroup = str(row.get("Orthogroup") or "")
                    family_id = str(row.get("BUSCO_family") or "")
                    row["OG_Paralog_Classification"] = paralog_classification.get(orthogroup, "NA")
                    row["Clean_1to1_NoParalogs"] = "YES" if (family_id, orthogroup) in clean_pairs else "NO"
                    writer.writerow(row)
        compare_results.setdefault("files", {})
        compare_results["files"]["species_paralog_tsv"] = species_tsv
        compare_results["files"]["augmented_mapping_tsv"] = augmented_mapping_tsv
        compare_results["files"]["good_families_txt"] = final_list_path
        compare_results["family_species_paralog_status"] = family_species_status
        analyser.log(f"Species-specific paralog mapping written: {species_tsv}")
        analyser.log(f"Augmented mapping written: {augmented_mapping_tsv}")
        analyser.last_compare_results = compare_results
        return cleaned_busco_families

    def _create_cleaned_reference_runs(self, compare_results: Dict[str, object]) -> Dict[str, int]:
        mode = self._clean_refs_mode()
        if not mode:
            return {}

        family_species_status = compare_results.get("family_species_paralog_status", {}) or {}
        created_runs: Dict[str, int] = {}

        for accession in self.accessions:
            source_run_id = self._resolve_source_busco_run_id(accession)
            if source_run_id is None:
                self.log(f"Skipping cleaned BUSCO run for {accession}: no source BUSCO run selected.", "WARNING")
                continue
            run_row = self.db_manager.busco.get_run(int(source_run_id))
            if not run_row:
                self.log(f"Skipping cleaned BUSCO run for {accession}: source BUSCO run {source_run_id} missing.", "WARNING")
                continue
            source_family_rows = self.db_manager.busco.get_run_family_data(int(source_run_id))
            if not source_family_rows:
                self.log(f"Skipping cleaned BUSCO run for {accession}: source BUSCO run {source_run_id} has no family rows.", "WARNING")
                continue
            source_location_rows = self.db_manager.busco.get_run_family_locations(int(source_run_id))
            source_status_map = self._family_status_map_from_rows(source_family_rows)
            source_gene_map = self._source_gene_name_map(source_family_rows, accession)

            rewritten_families: set[str] = set()
            rewritten_rows = []
            accession_token = canonicalize_accession(accession)
            for family_id, library_id, row_accession, status, sequence, score, length in source_family_rows:
                fam = str(family_id)
                row_acc = canonicalize_accession(row_accession)
                row_status = status
                if row_acc == accession_token:
                    species_status = (family_species_status.get(fam, {}) or {}).get(accession_token, {})
                    row_status = self._rewrite_family_status(
                        int(status),
                        has_in=bool(species_status.get("has_inparalogs")),
                        has_out=bool(species_status.get("has_outparalogs")),
                        mode=mode,
                    )
                    if int(row_status) != int(status):
                        rewritten_families.add(fam)
                rewritten_rows.append((family_id, library_id, row_accession, row_status, sequence, score, length))

            if not rewritten_rows:
                continue
            rewritten_status_map = self._family_status_map_from_rows(rewritten_rows)

            counts = {
                "Single copy BUSCOs": int(run_row[8] or 0),
                "Multi copy BUSCOs": int(run_row[9] or 0),
                "Fragmented BUSCOs": int(run_row[10] or 0),
                "Missing BUSCOs": int(run_row[11] or 0),
            }
            converted = sum(1 for fam in rewritten_families if int(source_status_map.get(fam, 0) or 0) == 1)
            if converted:
                counts["Single copy BUSCOs"] = max(0, counts["Single copy BUSCOs"] - converted)
                counts["Multi copy BUSCOs"] += converted

            new_run_id = self.db_manager.busco.create_run(
                accession=str(accession),
                library_id=int(self.parent_library_id),
                lineage_name=str(run_row[3] or self.parent_library_name),
                input_mode="protein",
                pipeline="orthofinder",
                pipeline_params_effective={
                    "derived_from_run_id": int(source_run_id),
                    "orthofinder_id": int(self.orthofinder_id) if self.orthofinder_id is not None else None,
                    "derived_library_id": int(self.library_id) if self.library_id is not None else None,
                    "derived_library_name": self.name,
                    "cleaning_mode": mode,
                    "set_cleaned_primary": bool(self.set_cleaned_primary),
                },
                pipeline_params_source={
                    "source_pipeline": str(run_row[5] or ""),
                    "source_input_mode": str(run_row[4] or ""),
                },
                result_dir=str(run_row[6] or ""),
                proteome_profile_id=int(run_row[14]) if run_row[14] is not None else None,
                status="completed",
            )
            if new_run_id is None:
                self.log(f"Failed to create cleaned BUSCO run for {accession}.", "WARNING")
                continue

            if not self.db_manager.busco.add_run_family_data(int(new_run_id), rewritten_rows):
                self.db_manager.busco.delete_run(int(new_run_id))
                self.log(f"Failed to clone BUSCO family rows for cleaned run {new_run_id}.", "WARNING")
                continue
            if source_location_rows and not self.db_manager.busco.add_run_family_locations(int(new_run_id), source_location_rows):
                self.db_manager.busco.delete_run(int(new_run_id))
                self.log(f"Failed to clone BUSCO family locations for cleaned run {new_run_id}.", "WARNING")
                continue
            if not self.db_manager.busco.update_run(int(new_run_id), counts=counts, completed=True):
                self.db_manager.busco.delete_run(int(new_run_id))
                self.log(f"Failed to finalize cleaned BUSCO run {new_run_id}.", "WARNING")
                continue

            if self.set_cleaned_primary:
                self.db_manager.busco.set_primary_run(
                    accession=str(accession),
                    library_id=int(self.parent_library_id),
                    run_id=int(new_run_id),
                    purpose="default",
                    policy="auto_clean_refs",
                    updated_by="add-library",
                )
            self._write_cleaned_reference_report(
                accession=str(accession),
                source_status_map=source_status_map,
                source_gene_map=source_gene_map,
                rewritten_status_map=rewritten_status_map,
                compare_results=compare_results,
            )
            created_runs[str(accession)] = int(new_run_id)
            self.log(
                f"Created cleaned BUSCO run {new_run_id} for {accession} using mode={mode}; rewritten families={len(rewritten_families)}.",
                "INFO",
            )

        return created_runs

    def _write_library_build_metadata(
        self,
        *,
        accepted_family_count: int,
        cleaned_busco_families_path: str,
        compare_results: Optional[Dict[str, object]] = None,
    ) -> str:
        if not self.location:
            raise ValueError("Library location is required to write build metadata.")
        metadata_path = os.path.join(str(self.location), "library_build_metadata.json")
        compare_results = compare_results or {}
        files = compare_results.get("files", {}) or {}
        payload = {
            "library_id": self.library_id,
            "library_name": self.name,
            "parent_library_id": self.parent_library_id,
            "parent_library_name": self.parent_library_name,
            "core_set_strategy": self._core_set_strategy(),
            "gene_tree_source": self._effective_gene_tree_source(),
            "accepted_family_rule": self._accepted_family_rule(),
            "min_species_in_trees": int(self.min_species_in_trees),
            "rerun_busco": bool(self.rerun_busco),
            "rerun_orthofinder": bool(self.rerun_orthofinder),
            "rerun_gene_trees_requested": bool(self.rerun_gene_trees),
            "rerun_gene_trees_effective": bool(self._rerun_gene_trees_effective()),
            "annotate_og_trees_requested": bool(self.annotate_og_trees),
            "annotate_og_trees_effective": bool(self._annotate_og_trees_effective()),
            "clean_refs_mode_effective": self._clean_refs_mode(),
            "accepted_family_count": int(accepted_family_count),
            "output_paths": {
                "cleaned_busco_families_json": cleaned_busco_families_path,
                "mapping_tsv": files.get("mapping_tsv"),
                "one_to_one_tsv": files.get("one_to_one_tsv"),
                "occupancy_tsv": files.get("occupancy_tsv"),
                "og_to_busco_tsv": files.get("og_to_busco_tsv"),
                "busco_to_og_tsv": files.get("busco_to_og_tsv"),
                "exact_mapping_tsv": files.get("exact_mapping_tsv"),
                "species_paralog_tsv": files.get("species_paralog_tsv"),
                "augmented_mapping_tsv": files.get("augmented_mapping_tsv"),
                "orthogroup_tree_manifest": files.get("orthogroup_tree_manifest"),
                "effective_gene_tree_dir": files.get("effective_gene_tree_dir"),
                "fasttree_gene_tree_dir": files.get("fasttree_gene_tree_dir"),
                "iqtree_gene_tree_dir": files.get("iqtree_gene_tree_dir"),
            },
        }
        with open(metadata_path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=4, sort_keys=True)
        return metadata_path

    def run(self):

        # This is the process to create a new library - part 1 of the pipeline

        # Create the library in the database if it doesn't already exist
        if not self.name:
            return self.handle_exception("Library name is not specified.", {"name": self.name})
        if not self.parent_library_name:
            return self.handle_exception("Parent library name is not specified.", {"parent_library_name": self.parent_library_name})
        if not self.accessions or not isinstance(self.accessions, list):
            return self.handle_exception("Accessions list is not specified or invalid.", {"accessions": self.accessions})
        if len(self.accessions) < 2:
            return self.handle_exception("At least two accessions are required to create a library.", {"accessions": self.accessions})
        
        if not self.coverage_taxid:
            if self.coverage is None:
                return self.handle_exception("Coverage taxid is not specified.", {"coverage_taxid": self.coverage_taxid})

            inferred_taxid: Optional[int]
            if isinstance(self.coverage, int):
                inferred_taxid = self.coverage
            else:
                # Attempt numeric coercion first (e.g., user passed "33208") before resolving by name
                candidate = None
                try:
                    if isinstance(self.coverage, str) and self.coverage.strip():
                        candidate = int(self.coverage.strip())
                except ValueError:
                    candidate = None
                if candidate is not None:
                    inferred_taxid = candidate
                else:
                    inferred_taxid = resolve_clade_to_taxid(self.db_manager, str(self.coverage))

            if not inferred_taxid:
                return self.handle_exception(
                    f"Could not infer coverage taxid from '{self.coverage}'.",
                    {"coverage": self.coverage},
                )
            self.coverage_taxid = inferred_taxid

        try:
            coverage_taxid = int(self.coverage_taxid)
        except (TypeError, ValueError):
            return self.handle_exception("Coverage taxid is not specified or invalid.", {"coverage_taxid": self.coverage_taxid})

        if coverage_taxid <= 0:
            return self.handle_exception("Coverage taxid is not specified or invalid.", {"coverage_taxid": coverage_taxid})

        # Normalise stored values so downstream code (and checkpoints) have an integer taxid
        self.coverage_taxid = coverage_taxid
        self.coverage = coverage_taxid

        # Rather than keep this as a concurrent subtask it will now be a prerequisite to the library being added
        if not self.parent_library_id:
            self.parent_library_id = self.db_manager.libraries.get_id(self.parent_library_name)
            if not self.parent_library_id:
                return self.handle_exception(f"Parent library '{self.parent_library_name}' not found in database.", {"parent_library_name": self.parent_library_name})

        libraries_dir = self.db_manager.storage.get_root_base("libraries")
        if not libraries_dir:
            return self.handle_exception("Libraries directory is not configured.", {})
        libraries_dir = str(libraries_dir)

        if self.library_id and not self.location:
            existing = self.db_manager.libraries.get(self.library_id) or []
            if existing and existing[0][5]:
                self.location = existing[0][5]

        if self.stage < 1:
            existing_record = None
            if self.library_id:
                records = self.db_manager.libraries.get(self.library_id) or []
                existing_record = records[0] if records else None
            else:
                existing_id = self.db_manager.libraries.get_id(self.name)
                if existing_id:
                    if not self.force:
                        return self.handle_exception(
                            f"Library '{self.name}' already exists. Re-run with --force to rebuild.",
                            {"library_name": self.name},
                        )
                    self.library_id = existing_id
                    records = self.db_manager.libraries.get(existing_id) or []
                    existing_record = records[0] if records else None

            # Determine target location for this library
            if existing_record and existing_record[5]:
                location = existing_record[5]
            else:
                location = os.path.join(libraries_dir, self.name)

            if self.force and self.library_id:
                if not self._purge_existing_library_state(library_id=int(self.library_id), location=location):
                    return self.handle_exception(
                        "Failed to purge existing library data before rebuild.",
                        {
                            "library_id": self.library_id,
                            "location": location,
                            "preserve_orthofinder": not bool(self.rerun_orthofinder),
                        },
                    )

            try:
                os.makedirs(location, exist_ok=True)
            except OSError as exc:
                return self.handle_exception("Failed to create library directory.", {"location": location, "error": str(exc)})

            new_id = self.db_manager.libraries.add(
                library_name=self.name,
                taxid=self.coverage_taxid,
                size=len(self.accessions),
                location=location,
                parent_id=self.parent_library_id,
                ref_accessions=None,
            )
            if not new_id:
                return self.handle_exception(
                    "Failed to create or update library record.",
                    {"name": self.name, "parent_library_id": self.parent_library_id},
                )
            self.library_id = new_id
            self.location = location
            if self.force:
                self.log(f"Rebuilding library '{self.name}' (ID {self.library_id}) with force.", "WARNING")
            else:
                self.log(f"Ensured library '{self.name}' exists at {self.location} (ID {self.library_id}).", "DEBUG")

            # Persist identifiers for downstream stages/resume safety
            self.checkpoint(
                0,
                {
                    'parent_library_id': self.parent_library_id,
                    'library_id': self.library_id,
                    'location': self.location,
                    'coverage_taxid': self.coverage_taxid,
                    'coverage': self.coverage,
                    'coverage_label': self.coverage_label,
                    'force': self.force,
                },
            )

        # The user has specified the reference genomes they wish to use for analysis.
        # We must first add the metadata for these genomes to the database.

        outcome = self.manage_subtasks(
            stage=1,
            queue_fn=self.queue_metadata,
            done_fn=self.metadata_done,
            wait_seconds=int(self.data.get("metadata_wait_seconds", 5)),
            # Use internal retry counter; treat provided value as the max
            retry_key=None,
            max_retries=int(self.data.get("assembly_metadata_retries", 3)),
            incomplete_message_fn=self.metadata_incomplete_message,
            retry_incomplete=False,  # do not retry incomplete outcomes; immediately record message
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False
        # outcome == "CONTINUE" -> proceed to next phase

        # If data has been successfully added to the database we can now attach reference accessions
        # to the library (only those present so far), then proceed to download the assemblies.
        if not self.db_manager.libraries.add_reference_assemblies(self.library_id, self.accessions):
            return self.handle_exception("Failed to attach reference accessions to library.", {"library_id": self.library_id, "accessions": self.accessions})

        self.log("Metadata phase complete.", "INFO")

        outcome = self.manage_subtasks(
            stage=2,
            queue_fn=self.queue_downloads,
            done_fn=self.downloads_done,
            wait_seconds=int(self.data.get("download_wait_seconds", 0)),
            # Use internal retry counter; treat provided value as the max
            retry_key=None,
            max_retries=int(self.data.get("download_retries", 2)),
            incomplete_message_fn=self.downloads_incomplete_message,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        # Phase 3+: BUSCO 
        
        self.log("Download phase complete. Starting BUSCO analysis.", "INFO")
        try:
            self._pin_proteome_profile_inputs()
        except (FileNotFoundError, ValueError) as exc:
            return self.handle_exception(
                "Failed to pin proteome profiles for add-library.",
                {"error": str(exc), "accessions": self.accessions},
            )

        outcome = self.manage_subtasks(
            stage=3,
            queue_fn=self.queue_busco,
            done_fn=self.busco_done,
            wait_seconds=int(self.data.get("busco_wait_seconds", 0)),
            retry_key=None,
            max_retries=int(self.data.get("busco_retries", 0)),
            incomplete_message_fn=self.busco_incomplete_message,
            retry_incomplete=False,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        cleaned_busco_families = self.data.get("_prepared_cleaned_busco_families")
        compare_results = self.data.get("_prepared_compare_results")
        if self.stage < 6 or cleaned_busco_families is None:
            self.orthofinder_id, self.orthofinder_location = self.db_manager.orthofinder.assert_results_exist(
                self.accessions,
                mcl_inflation=self.orthofinder_mcl_inflation,
                profile_inputs=self._orthofinder_profile_inputs(),
            )
            if not self.orthofinder_id or not self.orthofinder_location:
                return self.handle_exception("OrthoFinder results not found after BUSCO phase despite checks indicating they should be present.", {"accessions": self.accessions})
            self.checkpoint(self.stage, {"orthofinder_id": self.orthofinder_id, "orthofinder_location": self.orthofinder_location})

            core_set_dir = f"{self.location}/core_set_analysis"
            if self.stage < 4 and os.path.exists(core_set_dir):
                shutil.rmtree(core_set_dir)
                self.log(f"Removing existing core-set analysis directory before rebuild: {core_set_dir}", "WARNING")

            analyser = OrthoBuscoAnalyzer(
                identifier=self.name,
                working_dir=core_set_dir,
                orthofinder_run_folder=self.orthofinder_location,
                gene_tree_dir=self._effective_tree_dir(),
            )

            busco_results_locations = []
            for accession in self.accessions:
                run_dir = self.db_manager.busco.get_primary_result_dir(
                    accession,
                    self.parent_library_id,
                    purpose="default",
                )
                if run_dir:
                    busco_results_locations.append(run_dir)

            for loc in busco_results_locations:
                if not os.path.exists(loc):
                    return self.handle_exception("BUSCO results location does not exist.", {"location": loc, "accession": loc.split("/")[-2], "library": self.parent_library_name})

            busco_fastas_dir = os.path.join(core_set_dir, "busco_fastas")
            if self.stage < 4 or not os.path.isdir(busco_fastas_dir):
                if not analyser.transform_busco_results(busco_results_locations, out_dir=None, update=True, force=True):
                    return self.handle_exception("Failed to transform BUSCO results for analysis.", {"accessions": self.accessions, "busco_results_locations": busco_results_locations})
            else:
                analyser.busco_dir = busco_fastas_dir

            try:
                exact_pairs, exact_pairs_filtered, compare_results = self._prepare_exact_busco_orthogroups(analyser)
            except Exception as exc:  # boundary: required BUSCO/OrthoFinder comparison failure becomes this task error.
                return self.handle_exception("Failed to compare BUSCO and OrthoFinder results.", {"error": str(exc), "orthofinder_id": self.orthofinder_id})
            if self._core_set_strategy() == "skip_paralog_analysis":
                cleaned_busco_families = sorted(str(family_id) for family_id, _orthogroup in exact_pairs_filtered)
            else:
                self._prepare_replacement_tree_workspace()

                outcome = self.manage_subtasks(
                    stage=4,
                    queue_fn=lambda: self._queue_replacement_mafft_subtasks(exact_pairs_filtered),
                    done_fn=lambda: self._replacement_mafft_done(exact_pairs_filtered),
                    wait_seconds=0,
                    max_retries=int(self.data.get("mafft_retries", 1)),
                    incomplete_message_fn=lambda: self._replacement_mafft_incomplete_message(exact_pairs_filtered),
                    retry_incomplete=True,
                )
                if outcome == "ERROR":
                    return "ERROR"
                if outcome is False:
                    return False

                outcome = self.manage_subtasks(
                    stage=5,
                    queue_fn=lambda: self._queue_replacement_iqtree_subtasks(exact_pairs_filtered),
                    done_fn=lambda: self._replacement_iqtree_done(exact_pairs_filtered),
                    wait_seconds=0,
                    max_retries=int(self.data.get("iqtree_retries", 1)),
                    incomplete_message_fn=lambda: self._replacement_iqtree_incomplete_message(exact_pairs_filtered),
                    retry_incomplete=True,
                )
                if outcome == "ERROR":
                    return "ERROR"
                if outcome is False:
                    return False

                try:
                    self._materialize_replacement_orthogroup_trees(exact_pairs_filtered, compare_results)
                    cleaned_busco_families = self._finalize_replacement_tree_analysis(
                        analyser,
                        exact_pairs,
                        exact_pairs_filtered,
                        compare_results,
                    )
                except Exception as exc:  # boundary: required replacement-tree analysis failure becomes this task error.
                    return self.handle_exception("Failed to compare BUSCO and OrthoFinder results.", {"error": str(exc), "orthofinder_id": self.orthofinder_id})

            if not cleaned_busco_families:
                return self.handle_exception("Failed to compare BUSCO and OrthoFinder results.", {"accessions": self.accessions, "orthofinder_id": self.orthofinder_id})
            
            self.log(
                f"BUSCO/OrthoFinder comparison identified {len(cleaned_busco_families)} core families.",
                "INFO",
            )

            cleaned_compare_results = getattr(analyser, "last_compare_results", {}) or compare_results or {}
            cleaned_run_ids = {}
            if self._core_set_strategy() != "skip_paralog_analysis":
                cleaned_run_ids = self._create_cleaned_reference_runs(cleaned_compare_results)
            if cleaned_run_ids:
                self.log(
                    f"Created {len(cleaned_run_ids)} OrthoFinder-derived BUSCO runs for reference taxa.",
                    "INFO",
                )
            self.checkpoint(
                self.stage,
                {
                    "_prepared_cleaned_busco_families": cleaned_busco_families,
                    "_prepared_compare_results": compare_results,
                },
            )
        outcome = self.manage_subtasks(
            stage=6,
            queue_fn=self._queue_annotation_subtask,
            done_fn=self._annotation_done,
            wait_seconds=0,
        )
        if outcome == "ERROR":
            return "ERROR"
        if outcome is False:
            return False

        # Save the cleaned BUSCO families to a file in the library location for reference
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        busco_families_file = f"{self.location}/cleaned_busco_families_{timestamp}.json"
        try:
            with open(busco_families_file, "w") as f:
                json.dump(cleaned_busco_families, f, indent=4)
            # Also write/refresh a canonical filename for easy lookup
            canonical_path = os.path.join(self.location, "cleaned_busco_families.json")
            try:
                with open(canonical_path, "w") as f:
                    json.dump(cleaned_busco_families, f, indent=4)
            except OSError as e:
                self.log(f"Failed to write canonical cleaned BUSCO list: {e}", "WARNING")
        except (OSError, TypeError, ValueError) as e:
            return self.handle_exception("Failed to save cleaned BUSCO families to file.", {"file": busco_families_file, "error": str(e)})
        self.log(f"Cleaned BUSCO families saved to {busco_families_file}", "INFO")
        metadata_path = None
        try:
            metadata_path = self._write_library_build_metadata(
                accepted_family_count=len(cleaned_busco_families),
                cleaned_busco_families_path=os.path.join(self.location, "cleaned_busco_families.json"),
                compare_results=compare_results if isinstance(compare_results, dict) else {},
            )
        except (OSError, TypeError, ValueError) as e:
            return self.handle_exception("Failed to write library build metadata.", {"error": str(e), "library_id": self.library_id})
        try:
            artifact_metadata = {
                "library_id": self.library_id,
                "library_name": self.name,
                "core_set_strategy": self._core_set_strategy(),
                "gene_tree_source": self._effective_gene_tree_source(),
            }
            self.db_manager.artifacts.register(
                owner_type="library",
                owner_id=self.library_id,
                artifact_type="library_root",
                path=self.location,
                is_dir=True,
                format="directory",
                metadata=artifact_metadata,
            )
            self.db_manager.artifacts.register(
                owner_type="library",
                owner_id=self.library_id,
                artifact_type="library_core_set_json",
                path=os.path.join(self.location, "cleaned_busco_families.json"),
                format="json",
                metadata=artifact_metadata,
            )
            if metadata_path:
                self.db_manager.artifacts.register(
                    owner_type="library",
                    owner_id=self.library_id,
                    artifact_type="library_build_metadata_json",
                    path=metadata_path,
                    format="json",
                    metadata=artifact_metadata,
                )
            compare_files = compare_results.get("files", {}) if isinstance(compare_results, dict) else {}
            manifest_path = str(compare_files.get("orthogroup_tree_manifest") or "").strip()
            if manifest_path:
                self.db_manager.artifacts.register(
                    owner_type="library",
                    owner_id=self.library_id,
                    artifact_type="library_orthogroup_tree_manifest_tsv",
                    path=manifest_path,
                    format="tsv",
                    metadata=artifact_metadata,
                )
        except Exception as e:  # boundary: optional artifact catalog metadata; library files and DB records are already present.
            self.log(f"Failed to register custom library artifacts: {e}", "WARNING")

        # Persist the subset of BUSCOs for this derived library into BUSCO_descriptions
        try:
            if isinstance(cleaned_busco_families, dict):
                family_list = list(cleaned_busco_families.keys())
            else:
                family_list = list(cleaned_busco_families)
            family_list = [str(f) for f in family_list]
            parent_desc = self.db_manager.libraries.get_busco_descriptions(self.parent_library_id, family_list) or []
            desc_map = {str(fid): (desc, link) for fid, _lib, desc, link in parent_desc}
            subset_rows = []
            for fam in family_list:
                desc, link = desc_map.get(fam, (None, None))
                subset_rows.append((fam, self.library_id, desc, link))
            if subset_rows:
                if not self.db_manager.libraries.add_busco_descriptions(subset_rows):
                    return self.handle_exception("Failed to record BUSCO subset for custom library.", {"library_id": self.library_id})
                self.db_manager.libraries.update_size(self.library_id, len(family_list))
                self.log(f"Recorded {len(family_list)} BUSCO families for library {self.library_id}.", "INFO")
        except Exception as e:  # boundary: required BUSCO subset persistence failure becomes this task error.
            return self.handle_exception("Failed to persist BUSCO subset for custom library.", {"error": str(e), "library_id": self.library_id})

        # All phases complete
        self.log(
            f"Library build complete for '{self.name}': references={len(self.accessions)}, "
            f"families={len(cleaned_busco_families)}, location={self.location}.",
            "INFO",
        )
        return True
