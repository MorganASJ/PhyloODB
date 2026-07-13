import hashlib
import json
from datetime import datetime

from .internal_decontamination import InternalDecontaminationTask
from ..reporting import resolve_report_base_path


class ExternalDecontaminationApplyTask(InternalDecontaminationTask):
    """Apply external decontamination results to a previous internal run and write under a new run_id."""

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=8):
        super().__init__(db_path, task_id, checkpoint, data, required_threads=required_threads)
        self.source_run_id = self.data.get("source_run_id") or self.data.get("run_id")
        self.new_run_id = self.data.get("new_run_id")
        self.report_path = self.data.get("report_path")
        self.run_label = self.data.get("run_label")
        if not self.report_path:
            self.report_path = str(
                resolve_report_base_path(
                    self,
                    namespace="external-decontamination-reports",
                    default_stem="external_decontamination",
                    run_label=self.run_label or self.library_name or self.library_id,
                    cache_attr="_external_decontamination_report_dir",
                )
            )

    def _hydrate_from_params(self, params: dict):
        if not isinstance(params, dict):
            return
        for key, attr in {
            "rank": "rank",
            "hit_window": "hit_window",
            "p_value_threshold": "p_value_threshold",
            "off_clade_fraction": "off_clade_fraction",
            "min_buscos": "min_buscos",
            "min_identity": "min_identity",
            "min_coverage": "min_coverage",
            "min_alignment_length": "min_alignment_length",
            "min_bitscore": "min_bitscore",
            "max_evalue": "max_evalue",
        }.items():
            if key in params and params[key] is not None:
                setattr(self, attr, params[key])
        if params.get("config_path"):
            self.config_path = params.get("config_path")
        if params.get("config_signature"):
            self.config_signature = params.get("config_signature")
        if params.get("run_label") and not self.run_label:
            self.run_label = params.get("run_label")
        if params.get("external_blast_output_dir") and not self.external_blast_output_dir:
            self.external_blast_output_dir = params.get("external_blast_output_dir")
        if params.get("external_reuse_blast_results") and not self.external_reuse_blast_results:
            self.external_reuse_blast_results = params.get("external_reuse_blast_results")
        if params.get("external_blast_db_path") and not self.external_blast_db_path:
            self.external_blast_db_path = params.get("external_blast_db_path")
        if params.get("external_blast_db_type") and not self.external_blast_db_type:
            self.external_blast_db_type = params.get("external_blast_db_type")
        if params.get("external_blast_program") and not self.external_blast_program:
            self.external_blast_program = params.get("external_blast_program")
        if params.get("external_max_target_seqs") and not self.external_max_target_seqs:
            self.external_max_target_seqs = params.get("external_max_target_seqs")

    def _generate_new_run_id(self, params: dict) -> str:
        fingerprint = {
            "source_run_id": self.source_run_id,
            "external_reuse_blast_results": self.external_reuse_blast_results,
            "external_blast_output_dir": self.external_blast_output_dir,
            "run_label": self.run_label,
            "params": params,
        }
        digest = hashlib.sha1(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()
        return f"idc_ext_{digest}"

    def run(self):
        if not self.source_run_id:
            return self.handle_exception("source_run_id is required for external apply.", {})

        run_row = self.db_manager.filtering.get_decontamination_run(self.source_run_id)
        if not run_row:
            return self.handle_exception("Source decontamination run not found.", {"run_id": self.source_run_id})

        (
            _run_id,
            target_library_id,
            busco_library_id,
            targets_json,
            _refs_json,
            params_json,
            _config_signature,
            _run_label,
            _date,
        ) = run_row

        try:
            params = json.loads(params_json) if params_json else {}
        except (TypeError, json.JSONDecodeError):
            params = {}

        self.library_id = target_library_id
        self.busco_lib_id = busco_library_id
        try:
            self.accessions = json.loads(targets_json) if targets_json else []
        except (TypeError, json.JSONDecodeError):
            self.accessions = []

        self._hydrate_from_params(params)
        if self.config_path:
            self._load_config()

        if not (self.external_reuse_blast_results or self.external_blast_output_dir):
            return self.handle_exception(
                "external_reuse_blast_results or external_blast_output_dir is required to apply external results.",
                {"run_id": self.source_run_id},
            )

        if not self.new_run_id:
            self.new_run_id = self._generate_new_run_id(params)

        # Copy run metadata
        run_label = self.run_label or params.get("run_label") or f"{self.source_run_id}_external"
        new_params = dict(params or {})
        new_params.update(
            {
                "source_run_id": self.source_run_id,
                "external_apply_only": True,
                "external_blast_output_dir": self.external_blast_output_dir,
                "external_reuse_blast_results": self.external_reuse_blast_results,
                "run_label": run_label,
            }
        )
        self.db_manager.filtering.add_decontamination_run(
            run_id=self.new_run_id,
            target_library_id=self.library_id,
            busco_library_id=self.busco_lib_id,
            targets_json=json.dumps(self.accessions),
            refs_json=json.dumps(self.accessions),
            params_json=json.dumps(new_params),
            config_signature=self.config_signature,
            run_label=run_label,
        )

        # Copy votes and summaries into new run_id
        votes = self.db_manager.filtering.get_decontamination_votes(run_id=self.source_run_id)
        for row in votes or []:
            (
                family_id,
                busco_library_id,
                target_library_id,
                accession,
                _run_id,
                expected_taxid,
                best_taxid,
                runner_taxid,
                rank,
                best_bitscore,
                delta_bitscore,
                decision,
                top_hits_json,
                busco_run_id,
            ) = row
            self.db_manager.filtering.add_decontamination_vote(
                family_id,
                busco_library_id,
                target_library_id,
                accession,
                self.new_run_id,
                expected_taxid,
                best_taxid,
                runner_taxid,
                rank,
                best_bitscore,
                delta_bitscore,
                decision,
                top_hits_json,
                busco_run_id=busco_run_id,
            )

        summaries = self.db_manager.filtering.get_decontamination_summary(run_id=self.source_run_id)
        for row in summaries or []:
            (
                accession,
                target_library_id,
                busco_library_id,
                _run_id,
                expected_taxid,
                majority_taxid,
                rank,
                buscos_tested,
                buscos_supporting,
                buscos_outside,
                off_clade_fraction,
                decision,
                params_json,
                busco_run_id,
            ) = row
            self.db_manager.filtering.add_decontamination_summary(
                accession,
                target_library_id,
                busco_library_id,
                self.new_run_id,
                expected_taxid,
                majority_taxid,
                rank,
                buscos_tested,
                buscos_supporting,
                buscos_outside,
                off_clade_fraction,
                decision,
                params_json=params_json,
                busco_run_id=busco_run_id,
                datetime=datetime.utcnow(),
            )

        # Apply external results to new run
        self.run_id = self.new_run_id
        self._apply_external_results()

        # Reports
        if self.report_path:
            acc_taxa = {}
            acc_rank_taxa = {}
            for acc in self.accessions:
                genome = self.db_manager.genomes.get(acc)
                taxid = genome[1] if genome else None
                acc_taxa[acc] = taxid
                acc_rank_taxa[acc] = self._taxon_at_rank(taxid)
            summaries = self.db_manager.filtering.get_decontamination_summary(run_id=self.run_id)
            self._write_external_report(acc_taxa, acc_rank_taxa)
            self._write_external_summary(summaries)

        self.log(
            f"External apply completed: source_run_id={self.source_run_id} new_run_id={self.new_run_id}",
            "INFO",
        )
        return True
