import csv
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path
from typing import Dict, List, Sequence

from ..task import Task

class CreateTaxonomyDB(Task):
    """
    Task to insert NCBI taxonomy data into the central database.
    If no taxdump is provided, downloads the latest one from NCBI.
    """

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads)
        self.stage = checkpoint if checkpoint is not None else 0
        self.path_to_taxdump = self.data.get("path_to_taxdump", None)
        self.retain_taxdump = self.data.get("retain_taxdump", False)
        self.working_dir = self.data.get("working_dir", tempfile.gettempdir())

    def run(self):
        self.log("Starting taxonomy database setup.", "INFO")
        cleanup = False

        # Download latest taxdump if not provided
        if not self.path_to_taxdump:
            self.log("No taxdump path provided; downloading the latest NCBI taxdump.", "INFO")
            try:
                url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
                self.path_to_taxdump = os.path.join(self.working_dir, "taxdump.tar.gz")
                result = urllib.request.urlretrieve(url, self.path_to_taxdump)
                if result:
                    self.log(f"Downloaded NCBI taxdump to {self.path_to_taxdump}.", "DEBUG")
                    cleanup = not self.retain_taxdump
                    self.checkpoint(stage=1, checkpoint_data={"path_to_taxdump": self.path_to_taxdump, "retain_taxdump": self.retain_taxdump, "working_dir": self.working_dir})
                else:
                    return self.handle_exception("Download returned no result.", {'url': url, 'result': result})
            except Exception as e:  # boundary: required network download failure becomes this task error.
                self.log(f"Failed to download taxdump: {e}", level="error")
                return self.handle_exception(e, self.stage)

        # Insert taxonomy data
        try:
            self.log("Inserting taxonomy data into the database.", "INFO")
            self.db_manager.genomes.insert_taxdump(self.path_to_taxdump)
            self.log("Taxonomy database creation complete.", "INFO")
        except Exception as e:  # boundary: required taxonomy ingestion failure becomes this task error.
            self.log(f"Failed to insert taxonomy data: {e}", level="error")
            return self.handle_exception(e, self.stage)

        # Clean up downloaded file if required
        if cleanup:
            try:
                os.remove(self.path_to_taxdump)
                self.log("Temporary taxdump deleted.", "DEBUG")
            except OSError as e:
                self.log(f"Warning: failed to delete temporary taxdump: {e}", level="warning")

        self.log("Taxonomy database setup finished successfully.", "INFO")
        return True

class GenerateLineageCsvTask(Task):
    """Task that resolves selector inputs to accessions and exports lineage details to CSV."""

    _RANK_COLUMNS: List[str] = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
        "superkingdom",
    ]

    def __init__(self, db_path, task_id, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads)
        self.stage = 0
        self.accessions = self.selector_accessions()
        self.taxid = self.data.get("taxid")
        self.rule_quantity = self.data.get("quantity")
        self.rule_rank = self.data.get("rank")
        self.library_id = self.data.get("library_id")
        self.lineage = self.data.get("lineage")
        self.output = self.data.get("output")
        self.protein_only = self.data.get("protein_only")

    def run(self):
        if not self.output:
            return self.handle_exception("Output path is required for lineage CSV export.", {"output": self.output})

        output_path = Path(self.output)
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            return self.handle_exception("Failed to create parent directory for output.", {"output": self.output, "error": str(exc)})

        busco_lib_id = self.library_id
        if busco_lib_id is None and self.lineage:
            busco_lib_id = self.db_manager.libraries.get_id(self.lineage)

        downloaded_only_flag = self.data.get("downloaded_only")
        released_after = self.data.get("after")
        released_before = self.data.get("before")
        level_filter = self.data.get("level")
        primary_only = self.data.get("primary_only")
        use_busco = self.data.get("use_busco")
        min_complete = self.data.get("min_completeness")
        min_sc = self.data.get("min_single_copy_complete")

        try:
            selected = self.prepare_selectors(
                taxid=self.taxid,
                rule_quantity=self.rule_quantity,
                rule_rank=self.rule_rank,
                busco_library_id=busco_lib_id,
                downloaded_only=downloaded_only_flag,
                released_after=released_after,
                released_before=released_before,
                level=level_filter,
                protein_only=self.protein_only,
                primary_only=primary_only,
                use_busco=use_busco,
                min_completeness=min_complete,
                min_single_copy_complete=min_sc,
            )
        except ValueError as exc:
            return self.handle_exception(
                str(exc),
                {
                    "quantity": self.rule_quantity,
                    "rank": self.rule_rank,
                    "taxid": self.taxid,
                    "after": released_after,
                    "before": released_before,
                    "level": level_filter,
                    "primary_only": primary_only,
                    "use_busco": use_busco,
                    "min_completeness": min_complete,
                    "min_single_copy_complete": min_sc,
                },
            )

        placeholders = ",".join("?" for _ in selected)
        try:
            self.db_manager.cursor.execute(
                f"SELECT accession, taxid FROM Genome WHERE accession IN ({placeholders})",
                tuple(selected),
            )
            acc_tax_rows = self.db_manager.cursor.fetchall() or []
        except Exception as exc:  # boundary: required selector taxid query failure becomes this task error.
            return self.handle_exception("Failed to resolve taxids for selected accessions.", {"error": str(exc)})

        acc_to_tax: Dict[str, int] = {}
        for acc, tax in acc_tax_rows:
            if acc is None or tax is None:
                continue
            acc_to_tax[str(acc)] = int(tax)

        records: List[Dict[str, str]] = []
        for acc in selected:
            taxid = acc_to_tax.get(acc)
            lineage_rows = self.db_manager.genomes.get_lineage_root_to_leaf(taxid) if taxid is not None else []
            lineage_map = {str(rank).lower(): name for (tid, name, rank, _parent) in lineage_rows if rank and name}
            taxon_name = next((name for (tid, name, rank, _parent) in lineage_rows if tid == taxid), "") if lineage_rows else ""
            taxon_rank = next((rank for (tid, _name, rank, _parent) in lineage_rows if tid == taxid and rank), "") if lineage_rows else ""

            row: Dict[str, str] = {
                "accession": acc,
                "taxid": str(taxid or ""),
                "taxon_name": taxon_name or "",
                "taxon_rank": str(taxon_rank or ""),
            }
            for col in self._RANK_COLUMNS:
                row[col] = lineage_map.get(col, "")
            records.append(row)

        if not records:
            return self.handle_exception("No lineage information found for selected accessions.", {"selected": selected})

        records.sort(key=lambda r: r["accession"])
        fieldnames = ["accession", "taxid", "taxon_name", "taxon_rank", *self._RANK_COLUMNS]

        try:
            with output_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(records)
        except (OSError, csv.Error) as exc:
            return self.handle_exception("Failed to write lineage CSV.", {"output": self.output, "error": str(exc)})

        self.log(f"Generated lineage CSV with {len(records)} accessions at {self.output}")
        return True


class FinalizeGenomeMoveTask(Task):
    """Copy/rebind genomes, verify the destination, then delete the source."""

    def __init__(self, db_path, task_id, checkpoint, data, required_threads=1):
        super().__init__(db_path, task_id, data, required_threads)
        self.stage = checkpoint if checkpoint is not None else 0
        self.rows = list(self.data.get("rows") or [])
        self.verify = bool(self.data.get("verify", True))
        self.tidy = bool(self.data.get("tidy", True))

    @staticmethod
    def _remove_path(path: str) -> None:
        if not path:
            return
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        elif os.path.exists(path):
            os.remove(path)

    def _copy_and_rebind(self) -> None:
        for row in self.rows:
            accession = str(row["accession"])
            source = os.path.abspath(str(row["source_path"]))
            dest = os.path.abspath(str(row["destination_path"]))
            action = str(row.get("action") or "move-files")
            if action not in {"move-files", "rebind-only"}:
                raise ValueError(f"{accession}: unsupported move action '{action}'.")
            operation_id = row.get("filesystem_operation_id")
            if operation_id is None:
                operation_id = self.db_manager.storage.create_filesystem_operation(
                    operation_type="move-genomes",
                    source_path=source,
                    destination_path=dest,
                    payload={"kind": "genomes", "row": row},
                )
                row["filesystem_operation_id"] = int(operation_id)
                self.data["rows"] = self.rows
                self.db_manager.tasks.update_data(self.task_id, data=self.data)
            if action == "move-files":
                if not os.path.exists(source):
                    raise ValueError(f"{accession}: source missing: {source}")
                if os.path.abspath(source) != os.path.abspath(dest):
                    if os.path.exists(dest):
                        raise ValueError(f"{accession}: destination already exists: {dest}")
                    os.makedirs(os.path.dirname(dest), exist_ok=True)
                    shutil.copytree(source, dest, symlinks=True)
                self.db_manager.storage.update_filesystem_operation(
                    int(operation_id), status="prepared"
                )
            else:
                if not os.path.exists(dest):
                    raise ValueError(f"{accession}: rebind-only destination missing: {dest}")
            with self.db_manager.transaction(operation=f"commit genome move {accession}"):
                self.db_manager.genomes.set_binding(accession, dest, kind="genomes")
                self.db_manager.storage.update_filesystem_operation(
                    int(operation_id), status="db_committed"
                )

    def _rollback_to_source(self) -> None:
        for row in self.rows:
            accession = str(row["accession"])
            source = os.path.abspath(str(row["source_path"]))
            dest = os.path.abspath(str(row["destination_path"]))
            action = str(row.get("action") or "move-files")
            if os.path.exists(source):
                try:
                    self.db_manager.genomes.set_binding(accession, source, kind="genomes")
                except Exception as exc:  # boundary: rollback attempts all rows and logs failures.
                    self.log(f"Rollback failed to restore binding for {accession} to {source}: {exc}", "WARNING")
            if action == "move-files" and os.path.abspath(source) != os.path.abspath(dest):
                try:
                    if os.path.exists(dest):
                        self._remove_path(dest)
                except OSError as exc:
                    self.log(f"Rollback failed to remove staged destination {dest}: {exc}", "WARNING")
            operation_id = row.get("filesystem_operation_id")
            if operation_id is not None:
                try:
                    self.db_manager.storage.update_filesystem_operation(
                        int(operation_id), status="rolled_back"
                    )
                except Exception as exc:  # boundary: rollback attempts all rows and logs journal failures.
                    self.log(f"Rollback failed to mark filesystem operation {operation_id} rolled_back: {exc}", "WARNING")

    def _delete_sources(self) -> list[str]:
        leftovers: list[str] = []
        for row in self.rows:
            source = os.path.abspath(str(row["source_path"]))
            dest = os.path.abspath(str(row["destination_path"]))
            action = str(row.get("action") or "move-files")
            if action != "move-files" or os.path.abspath(source) == os.path.abspath(dest):
                continue
            try:
                if os.path.exists(source):
                    self._remove_path(source)
                operation_id = row.get("filesystem_operation_id")
                if operation_id is not None:
                    self.db_manager.storage.update_filesystem_operation(
                        int(operation_id), status="finalized"
                    )
            except Exception as exc:  # boundary: source cleanup failure is recoverable through storage recover.
                leftovers.append(source)
                operation_id = row.get("filesystem_operation_id")
                if operation_id is not None:
                    self.db_manager.storage.update_filesystem_operation(
                        int(operation_id),
                        status="failed",
                        error_message=f"Could not delete source path: {source}: {exc}",
                    )
        return leftovers

    def _queue_verify_assembly_child(self) -> bool:
        accessions = [str(row["accession"]) for row in self.rows]
        self.queue_subtask(
            job_type=18,
            status="P",
            priority=1,
            data={
                "accessions": accessions,
                "repair": True,
                "tidy": bool(self.tidy),
            },
        )
        return True

    def _queue_verify_busco_child(self) -> bool:
        accessions = [str(row["accession"]) for row in self.rows]
        self.queue_subtask(
            job_type=20,
            status="P",
            priority=1,
            data={
                "accessions": accessions,
                "repair": True,
                "reingest": True,
            },
        )
        return True

    def _child_ids(self) -> Sequence[int]:
        subtasks = self.db_manager.tasks.get_subtasks(self.task_id) or []
        return [int(row[0]) for row in subtasks if row and row[0] is not None]

    def run(self):
        if not self.rows:
            return self.handle_exception("No genome move rows were provided.", {})

        if self.stage < 1:
            try:
                self._copy_and_rebind()
            except Exception as exc:  # boundary: required staging failure triggers rollback and task error.
                self._rollback_to_source()
                return self.handle_exception(exc, context="Genome move staging failed")

        if self.verify:
            outcome = self.manage_subtasks(
                stage=1,
                queue_fn=self._queue_verify_assembly_child,
                done_fn=lambda: self._subtasks_state() == "complete",
                wait_seconds=0,
                retry_key=None,
            )
            if outcome is False:
                return False
            if outcome == "ERROR":
                self._rollback_to_source()
                summary, _stacks = self.aggregate_subtask_errors(
                    "Genome move verification failed.",
                    subtask_ids=list(self._child_ids()),
                )
                return self.handle_exception(
                    f"{summary} Rollback restored original genome bindings and removed copied destinations.",
                    context="Genome move verification failed",
                )
            outcome = self.manage_subtasks(
                stage=2,
                queue_fn=self._queue_verify_busco_child,
                done_fn=lambda: self._subtasks_state() == "complete",
                wait_seconds=0,
                retry_key=None,
            )
            if outcome is False:
                return False
            if outcome == "ERROR":
                self._rollback_to_source()
                summary, _stacks = self.aggregate_subtask_errors(
                    "Genome move BUSCO verification failed.",
                    subtask_ids=list(self._child_ids()),
                )
                return self.handle_exception(
                    f"{summary} Rollback restored original genome bindings and removed copied destinations.",
                    context="Genome move BUSCO verification failed",
                )

        leftovers = self._delete_sources()
        if leftovers:
            self.log(
                "Genome move completed but some source paths could not be deleted: " + ", ".join(leftovers[:10]),
                level="WARNING",
            )
        return True
