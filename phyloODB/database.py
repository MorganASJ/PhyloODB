from __future__ import annotations

import json
import os
import sqlite3
import stat
import threading
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path
from typing import Any, Optional, Sequence

from .db import (
    ArtifactRepository,
    BuscoRepository,
    EnvRepository,
    FilteringRepository,
    GenomeRepository,
    LibraryRepository,
    OrthoFinderRepository,
    ProteomeRepository,
    SelectorPresetRepository,
    StorageRepository,
    TaskRepository,
)
from .db.core import DatabaseCore, sqlite_busy_timeout_ms
from .db.errors import MigrationError, PhyloODBDatabaseError, SchemaCompatibilityError
from .db.schema import (
    CURRENT_SCHEMA_VERSION,
    ensure_assembly_accession_alias_schema,
    ensure_taxonomy_schema,
    ensure_busco_run_schema,
    ensure_environment_variable_schema,
    ensure_proteome_schema,
    ensure_selector_preset_schema,
    ensure_task_queue_schema,
    ensure_storage_schema,
    setup_database as setup_database_schema,
)

sqlite3.register_adapter(date, lambda d: d.isoformat())
sqlite3.register_adapter(datetime, lambda d: d.strftime("%Y-%m-%d %H:%M:%S"))
sqlite3.register_converter("DATE", lambda s: datetime.strptime(s.decode("utf-8"), "%Y-%m-%d").date())
sqlite3.register_converter("DATETIME", lambda s: datetime.strptime(s.decode("utf-8"), "%Y-%m-%d %H:%M:%S"))


class DBManager:
    """Application DB entrypoint composed from domain repositories."""

    def __init__(self, db_path, *, read_only: bool = False):
        self.db_path = db_path
        self.read_only = bool(read_only)
        self._paralog_filtering_has_target = None
        self._paralog_filtering_has_runs = None
        self._env_updated_at = None
        self.conn = None
        self.cursor = None
        self._cursor_lock = threading.RLock()
        self._transaction_depth = 0
        self._savepoint_counter = 0
        self._syncing_storage_env = False
        self._busco_run_context = {
            "pipeline": None,
            "input_mode": None,
            "prefer_pipeline": None,
            "prefer_input_mode": None,
            "proteome_profile": None,
            "prefer_proteome_profile": None,
            "selection": "primary",
        }
        self._busco_compat_audit = {
            "counts": defaultdict(int),
            "samples": defaultdict(list),
        }
        self._init_repositories()

    def _init_repositories(self) -> None:
        self.tasks = TaskRepository(self)
        self.env = EnvRepository(self)
        self.genomes = GenomeRepository(self)
        self.libraries = LibraryRepository(self)
        self.busco = BuscoRepository(self)
        self.orthofinder = OrthoFinderRepository(self)
        self.filtering = FilteringRepository(self)
        self.artifacts = ArtifactRepository(self)
        self.proteomes = ProteomeRepository(self)
        self.selector_presets = SelectorPresetRepository(self)
        self.storage = StorageRepository(self)

    @staticmethod
    def _json_dump(value: Any) -> str:
        return json.dumps(value)

    @staticmethod
    def _json_load(value: Optional[str], default: Any = None) -> Any:
        if value is None:
            return default
        try:
            return json.loads(value)
        except (json.JSONDecodeError, TypeError):
            return default

    def get_path(self):
        return self.db_path

    def connect(self):
        # Bootstrap a persisted shared-project umask before SQLite can create
        # WAL/SHM sidecars. Creation handles the not-yet-existing DB explicitly.
        if not self.read_only and Path(self.db_path).is_file():
            bootstrap_conn = None
            try:
                bootstrap_conn = sqlite3.connect(
                    f"{Path(self.db_path).resolve().as_uri()}?mode=ro",
                    uri=True,
                    timeout=1.0,
                )
                rows = bootstrap_conn.execute(
                    "SELECT var_name, var_value FROM Environment_Variables "
                    "WHERE var_name IN ('PROJECT_PERMISSION_MODE', 'SHARED_GROUP')"
                ).fetchall()
                values = {str(name): json.loads(value) for name, value in rows}
                from .permissions import apply_shared_umask, policy_from_values
                apply_shared_umask(policy_from_values(values))
            except (sqlite3.Error, ValueError):
                pass
            finally:
                if bootstrap_conn is not None:
                    bootstrap_conn.close()
        connect_target = self.db_path
        connect_kwargs = {}
        if self.read_only:
            if not Path(self.db_path).exists():
                raise FileNotFoundError(f"Database does not exist: {self.db_path}")
            database_path = Path(self.db_path).resolve()
            parent_mode = database_path.parent.stat().st_mode
            immutable_fallback = not bool(parent_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
            immutable_option = "&immutable=1" if immutable_fallback else ""
            connect_target = f"{database_path.as_uri()}?mode=ro{immutable_option}"
            connect_kwargs["uri"] = True
        try:
            self.conn = sqlite3.connect(
                connect_target,
                detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
                timeout=sqlite_busy_timeout_ms() / 1000.0,
                check_same_thread=False,
                **connect_kwargs,
            )
            self.cursor = self.conn.cursor()
            if self.read_only:
                self.cursor.execute(f"PRAGMA busy_timeout = {sqlite_busy_timeout_ms()}")
                self.cursor.execute("PRAGMA mmap_size = 0")
                self.cursor.execute("PRAGMA query_only = ON")
            else:
                self.cursor.execute(f"PRAGMA busy_timeout = {sqlite_busy_timeout_ms()}")
                self.cursor.execute("PRAGMA journal_mode = WAL")
                self.cursor.execute("PRAGMA mmap_size = 0")
            self.cursor.execute("PRAGMA foreign_keys = ON")
        except Exception as exc:  # boundary: connection setup failures must close partial handles and add path context.
            self.close()
            raise PhyloODBDatabaseError(f"Failed to open database '{self.db_path}': {exc}") from exc
        return True

    def get_schema_version(self) -> int:
        row = self.cursor.execute("PRAGMA user_version").fetchone()
        return int(row[0] if row else 0)

    def _has_legacy_core_schema(self) -> bool:
        names = {
            str(row[0])
            for row in self.cursor.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table' AND name IN ('Genome', 'Tasks', 'Environment_Variables')"
            ).fetchall()
        }
        return {"Genome", "Tasks", "Environment_Variables"}.issubset(names)

    def validate_schema(self) -> int:
        version = self.get_schema_version()
        if version == CURRENT_SCHEMA_VERSION:
            return version
        if version > CURRENT_SCHEMA_VERSION:
            raise SchemaCompatibilityError(
                f"Database schema version {version} is newer than supported version "
                f"{CURRENT_SCHEMA_VERSION}; upgrade PhyloODB."
            )
        if version == 0 and not self._has_legacy_core_schema():
            raise SchemaCompatibilityError(
                f"Database '{self.db_path}' is uninitialized; run "
                f"'phyloODB {self.db_path} create'."
            )
        raise SchemaCompatibilityError(
            f"Database schema version {version} requires migration; run "
            f"'phyloODB {self.db_path} migrate'."
        )

    def migrate_database(self) -> list[str]:
        version = self.get_schema_version()
        if version > CURRENT_SCHEMA_VERSION:
            raise SchemaCompatibilityError(
                f"Database schema version {version} is newer than supported version "
                f"{CURRENT_SCHEMA_VERSION}; upgrade PhyloODB."
            )
        if version == CURRENT_SCHEMA_VERSION:
            return []
        if version == 0 and not self._has_legacy_core_schema():
            raise SchemaCompatibilityError(
                f"Database '{self.db_path}' is uninitialized; run "
                f"'phyloODB {self.db_path} create'."
            )

        applied: list[str] = []
        try:
            with self.transaction(operation=f"schema migration {version} -> {CURRENT_SCHEMA_VERSION}"):
                ensure_taxonomy_schema(self)
                ensure_assembly_accession_alias_schema(self)
                ensure_environment_variable_schema(self)
                self.ensure_task_queue_schema()
                self.ensure_busco_run_schema()
                self.ensure_storage_schema()
                self.ensure_proteome_schema()
                ensure_selector_preset_schema(self)
                self.storage.ensure_default_roots_from_env()
                self.storage.normalize_root_states(sync_env=True, promote_if_none=False)
                self.storage.backfill_table_locations(table="Genome", key_col="accession", kind="genomes")
                self.storage.backfill_table_locations(table="Libraries", key_col="library_id", kind="libraries")
                self.cursor.execute(f"PRAGMA user_version = {CURRENT_SCHEMA_VERSION}")
            applied.append(f"{version} -> {CURRENT_SCHEMA_VERSION}")
        except Exception as exc:  # boundary: migration failures retain original cause and identify version step.
            raise MigrationError(
                f"Migration from schema version {version} to {CURRENT_SCHEMA_VERSION} failed: {exc}"
            ) from exc
        return applied

    def transaction(self, *, operation: str = "database write"):
        """Open an atomic transaction, nesting through SQLite savepoints."""

        return DatabaseCore(self).transaction(operation=operation)

    def commit(self) -> None:
        DatabaseCore(self).commit()

    def rollback(self) -> None:
        DatabaseCore(self).rollback()

    def close(self):
        if self.conn:
            self.conn.close()
            self.conn = None
            self.cursor = None
        return True

    def _table_exists(self, name: str) -> bool:
        self.cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name = ?",
            (name,),
        )
        return self.cursor.fetchone() is not None

    def _column_exists(self, table: str, column: str) -> bool:
        self.cursor.execute(f"PRAGMA table_info({table})")
        cols = [row[1] for row in self.cursor.fetchall() or []]
        return column in cols

    def _env_has_updated_at(self) -> bool:
        cached = getattr(self, "_env_updated_at", None)
        if cached is None:
            cached = self._column_exists("Environment_Variables", "updated_at")
            self._env_updated_at = cached
        return bool(cached)

    def set_busco_run_context(
        self,
        *,
        pipeline: Optional[str] = None,
        input_mode: Optional[str] = None,
        prefer_pipeline: Optional[str] = None,
        prefer_input_mode: Optional[str] = None,
        proteome_profile: Optional[str] = None,
        prefer_proteome_profile: Optional[str] = None,
        run_ids: Optional[Sequence[Any]] = None,
        selection: Optional[str] = None,
    ) -> None:
        def _normalize_mode(value: Optional[str]) -> Optional[str]:
            if value is None:
                return None
            token = str(value).strip().lower()
            if token in {"nucl", "nucleotide"}:
                return "genome"
            if token == "prot":
                return "protein"
            return token or None

        if pipeline is not None:
            self._busco_run_context["pipeline"] = str(pipeline).strip().lower() or None
        if input_mode is not None:
            self._busco_run_context["input_mode"] = _normalize_mode(input_mode)
        if prefer_pipeline is not None:
            self._busco_run_context["prefer_pipeline"] = str(prefer_pipeline).strip().lower() or None
        if prefer_input_mode is not None:
            self._busco_run_context["prefer_input_mode"] = _normalize_mode(prefer_input_mode)
        if proteome_profile is not None:
            self._busco_run_context["proteome_profile"] = str(proteome_profile).strip() or None
        if prefer_proteome_profile is not None:
            self._busco_run_context["prefer_proteome_profile"] = str(prefer_proteome_profile).strip() or None
        if run_ids is not None:
            normalized_run_ids = []
            for run_id in run_ids or []:
                text = str(run_id).strip()
                if text.startswith("@") and len(text) > 1:
                    value = self.get_environment_variable(text[1:])
                    if isinstance(value, str):
                        tokens = [part.strip() for part in value.split(",") if part.strip()]
                    elif isinstance(value, (list, tuple, set)):
                        tokens = [str(part).strip() for part in value if str(part).strip()]
                    else:
                        tokens = []
                    for token in tokens:
                        if token.isdigit():
                            normalized_run_ids.append(int(token))
                    continue
                if text.isdigit():
                    normalized_run_ids.append(int(text))
            self._busco_run_context["run_ids"] = list(dict.fromkeys(normalized_run_ids))
        if selection is not None:
            self._busco_run_context["selection"] = str(selection).strip().lower() or "primary"

    def _get_busco_context(self) -> dict:
        return dict(self._busco_run_context)

    def record_busco_compat_event(
        self,
        event: str,
        *,
        count: int = 1,
        accession: Optional[str] = None,
        library_id: Optional[int] = None,
        run_id: Optional[int] = None,
        detail: Optional[str] = None,
    ) -> None:
        token = str(event or "").strip()
        if not token:
            return
        try:
            self._busco_compat_audit["counts"][token] += max(int(count), 0)
        except (TypeError, ValueError):
            self._busco_compat_audit["counts"][token] += 1
        sample_parts: list[str] = []
        if accession is not None:
            sample_parts.append(f"acc={accession}")
        if library_id is not None:
            sample_parts.append(f"lib={int(library_id)}")
        if run_id is not None:
            sample_parts.append(f"run={int(run_id)}")
        if detail:
            sample_parts.append(str(detail))
        sample = " ".join(sample_parts).strip()
        if sample:
            samples = self._busco_compat_audit["samples"][token]
            if sample not in samples and len(samples) < 10:
                samples.append(sample)

    def get_busco_compat_audit(self, *, reset: bool = False) -> dict[str, dict[str, Any]]:
        snapshot = {
            "counts": dict(self._busco_compat_audit["counts"]),
            "samples": {key: list(values) for key, values in self._busco_compat_audit["samples"].items()},
        }
        if reset:
            self.reset_busco_compat_audit()
        return snapshot

    def reset_busco_compat_audit(self) -> None:
        self._busco_compat_audit = {
            "counts": defaultdict(int),
            "samples": defaultdict(list),
        }

    def ensure_task_queue_schema(self) -> None:
        ensure_task_queue_schema(self)

    def ensure_busco_run_schema(self) -> None:
        ensure_busco_run_schema(self)

    def ensure_storage_schema(self) -> None:
        ensure_storage_schema(self)

    def ensure_proteome_schema(self) -> None:
        ensure_proteome_schema(self)

    def setup_database(self):
        setup_database_schema(self)
        self.storage.normalize_root_states(sync_env=True, promote_if_none=True)
        self.cursor.execute(f"PRAGMA user_version = {CURRENT_SCHEMA_VERSION}")
        self.conn.commit()

    def reset(self):
        if self.conn:
            self.close()
        if os.path.exists(self.db_path):
            try:
                os.remove(self.db_path)
            except OSError as exc:
                raise RuntimeError(f"Failed to remove database '{self.db_path}': {exc}") from exc
            for ext in ("-wal", "-shm"):
                sidecar = f"{self.db_path}{ext}"
                try:
                    if os.path.exists(sidecar):
                        os.remove(sidecar)
                except OSError as exc:
                    raise RuntimeError(f"Failed to remove database sidecar '{sidecar}': {exc}") from exc
        self.connect()
        self.setup_database()
        return True

    # Compatibility facade methods
    def get_environment_variable(self, var_name):
        return self.env.get(var_name)

    def get_environment_variables(self, var_names=None):
        return self.env.get_many(var_names)

    def get_environment_variable_records(self, var_names=None):
        return self.env.get_records(var_names)

    def set_environment_variable(self, var_name, var_value, *, kind=None):
        return self.env.set(var_name, var_value, kind=kind)

    def set_environment_variables(self, env_vars, *, kind=None, kinds=None):
        return self.env.set_many(env_vars, kind=kind, kinds=kinds)

    def get_libraries(self, library_id=None):
        return self.libraries.get(library_id)

    def get_library_id(self, library_name):
        return self.libraries.get_id(library_name)

    def get_libraries_view(self):
        return self.libraries.get_view()

    def assert_library_has_parent(self, library_id):
        return self.libraries.assert_has_parent(library_id)

    def get_accessions_by_taxid(self, taxid: int, include_descendants: bool = True, status_min: int | None = 1, protein_only: bool = False):
        return self.genomes.get_accessions_by_taxid(
            taxid,
            include_descendants=include_descendants,
            status_min=status_min,
            protein_only=protein_only,
        )

    def get_genome(self, accession):
        return self.genomes.get(accession)

    def get_genomes(self, genome_ids=None, status=None):
        return self.genomes.get_many(genome_ids, status=status)

    def get_lineage_root_to_leaf(self, taxid):
        return self.genomes.get_lineage_root_to_leaf(taxid)

    def get_taxid_for_genus(self, genus):
        return self.genomes.get_taxid_for_genus(genus)

    def get_taxid_for_name(self, name):
        return self.genomes.get_taxid_for_name(name)

    def get_taxid_for_species(self, genus, species=None):
        return self.genomes.get_taxid_for_species(genus, species)

    def get_taxid_for_genus_species(self, name):
        return self.genomes.get_taxid_for_genus_species(name)

    def get_busco_processed_accessions(self, library=None):
        return self.busco.get_processed_accessions(library)

    def get_busco_processed_accessions_any(self):
        return self.busco.get_processed_accessions_any()

    def get_busco_results_percentages(self, *args, **kwargs):
        return self.busco.get_results_percentages(*args, **kwargs)

    def get_busco_results_adjusted(self, *args, **kwargs):
        return self.busco.get_results_adjusted(*args, **kwargs)

    def get_busco_runs_for_accessions(self, *args, **kwargs):
        return self.busco.get_runs_for_accessions(*args, **kwargs)

    def get_busco_family_results_for_library(self, *args, **kwargs):
        return self.busco.get_family_results_for_library(*args, **kwargs)

    def add_busco_results(self, accession, library_id, results, datetime=None):
        return self.busco.add_results(accession, library_id, results, datetime=datetime)

    def add_busco_family_multiple_data(self, rows):
        return self.busco.add_legacy_family_data(rows)

    def add_busco_family_multiple_locations(self, rows):
        return self.busco.add_legacy_family_locations(rows)

    def add_busco_run_family_multiple_data(self, run_id, rows):
        return self.busco.add_run_family_data(run_id, rows)

    def add_busco_run_family_multiple_locations(self, run_id, rows):
        return self.busco.add_run_family_locations(run_id, rows)

    def get_decontamination_decision_percentages(self, *args, **kwargs):
        return self.filtering.get_decontamination_decision_percentages(*args, **kwargs)

    def get_latest_decontamination_summary_with_fallback(self, *args, **kwargs):
        return self.filtering.get_latest_decontamination_summary_with_fallback(*args, **kwargs)

    def get_decontamination_accessions(self, *args, **kwargs):
        return self.filtering.get_decontamination_accessions(*args, **kwargs)

    def get_decontamination_summary(self, *args, **kwargs):
        return self.filtering.get_decontamination_summary(*args, **kwargs)

    def get_decontamination_votes(self, *args, **kwargs):
        return self.filtering.get_decontamination_votes(*args, **kwargs)

    def get_paralog_filtering_accessions(self, *args, **kwargs):
        return self.filtering.get_paralog_filtering_accessions(*args, **kwargs)

    def get_paralog_results(self, *args, **kwargs):
        return self.filtering.get_paralog_results(*args, **kwargs)

    def add_paralog_filtering_copy_result(self, *args, **kwargs):
        return self.filtering.add_paralog_filtering_copy_result(*args, **kwargs)

    def add_paralog_filtering_results(self, *args, **kwargs):
        return self.filtering.add_paralog_filtering_result(*args, **kwargs)

    def add_decontamination_copy_vote(self, *args, **kwargs):
        return self.filtering.add_decontamination_copy_vote(*args, **kwargs)

    def add_decontamination_vote(self, *args, **kwargs):
        return self.filtering.add_decontamination_vote(*args, **kwargs)

    def add_decontamination_summary(self, *args, **kwargs):
        return self.filtering.add_decontamination_summary(*args, **kwargs)

    def add_proteome_blastdb(self, *args, **kwargs):
        return self.filtering.add_proteome_blastdb(*args, **kwargs)

    def _get_library_size(self, library_id: int) -> int:
        return self.busco._get_library_size(library_id)

    def _paralog_hidden_counts(self, *args, **kwargs):
        return self.filtering.paralog_hidden_counts(*args, **kwargs)

    def queue_task(self, *args, **kwargs):
        return self.tasks.queue(*args, **kwargs)

    def get_task_by_id(self, task_id):
        return self.tasks.get(task_id)

    def get_task_status(self, task_id):
        return self.tasks.get_status(task_id)

    def get_tasks(self, task_id=None):
        return self.tasks.get_many(task_id)

    def get_task_error_info(self, task_id):
        return self.tasks.get_error_info(task_id)

    def get_task_errors_from_subtasks(self, task_id):
        return self.tasks.get_errors_from_subtasks(task_id)

    def get_task_blocks(self, task_ids=None, *, unsatisfied_only: bool = True):
        return self.tasks.get_blocks(task_ids, unsatisfied_only=unsatisfied_only)

    def get_task_dependencies(self, task_ids=None):
        return self.tasks.get_dependencies(task_ids)

    def get_task_time_constraints(self, task_ids=None):
        return self.tasks.get_time_constraints(task_ids)

    def cancel_task_and_descendants(self, task_id: int, reason: str = "Canceled by user"):
        return self.tasks.cancel_and_descendants(task_id, reason)

    def kill_task_and_descendants(self, task_id: int, reason: str = "Killed by user"):
        return self.tasks.kill_and_descendants(task_id, reason)

    def reset_tasks(self, full: bool = False):
        return self.tasks.reset(full=full)

    def clear_tasks(self, keep_running: bool = True):
        return self.tasks.clear(keep_running=keep_running)

    def insert_assembly_information(self, assembly_data):
        return self.genomes.insert_assembly(assembly_data)

    def upsert_assembly_information(self, assembly_data):
        return self.genomes.upsert_assembly(assembly_data)

    def insert_taxonomy_information(self, taxonomy_data_dict):
        return self.genomes.insert_taxonomy_information(taxonomy_data_dict)

    def insert_genome_information(self, genome_data):
        return self.genomes.insert(genome_data)

    def upsert_genome_information(self, genome_data):
        return self.genomes.upsert(genome_data)

    def _normalize_library_status(self, status: Any) -> str:
        if status is None:
            return "ready"
        if isinstance(status, bool):
            return "ready" if status else "stale"
        if isinstance(status, (int, float)):
            token = int(status)
            if token == 1:
                return "ready"
            if token == 0:
                return "stale"
        token = str(status).strip().lower()
        if token in {"", "1", "true", "ready"}:
            return "ready"
        if token in {"0", "false", "stale"}:
            return "stale"
        raise ValueError(f"Unsupported library status '{status}'. Use ready/stale or 1/0.")

    def add_library(
        self,
        library_name,
        taxid,
        size,
        location=None,
        parent_id=None,
        ref_accessions=None,
        odb_version=None,
        status=None,
    ):
        resolved_location = location
        if resolved_location is None:
            base = self.storage.require_root_base("libraries")
            resolved_location = os.path.join(str(base), str(library_name))
        normalized_status = self._normalize_library_status(status)
        return self.libraries.add(
            library_name,
            taxid,
            size,
            resolved_location,
            parent_id=parent_id,
            ref_accessions=ref_accessions,
            odb_version=odb_version,
            status=normalized_status,
        )

    def add_busco_descriptions(self, rows):
        return self.libraries.add_busco_descriptions(rows)


__all__ = ["DBManager"]
