"""Service layer responsible for queueing and running tasks via the registry."""
from __future__ import annotations

import traceback
from typing import Any, Dict, Optional, Sequence

from ..database import DBManager
from ..errors import TaskExecutionError
from ..logging_utils import configure_logging_from_db
from ..permissions import preflight_task
from ..registry import registry
from ..selector_utils import expand_accession_variables
from ..schemas import TaskPayload
from ..scheduling import BarrierConstraint, Constraint, DependencyConstraint, ScheduledConstraint, TimeConstraint
from ..thread_defaults import refresh_runtime_thread_defaults, resolve_task_required_threads


class TaskService:
    """Facade for queueing tasks with schema validation.
    """

    def __init__(self, db_path: str, *, db_manager: Optional[DBManager] = None):
        self.db_path = db_path
        self.db_manager = db_manager or DBManager(db_path)

    def _ensure_connection(self) -> None:
        if not getattr(self.db_manager, "conn", None):
            self.db_manager.connect()

    def _update_last_env(self, task_id: int, task_key: str) -> None:
        key = task_key.upper().replace("-", "_")
        self.db_manager.env.set_many(
            {
                "LAST": task_id,
                f"LAST_{key}": task_id,
            }
        )

    def _apply_constraints(self, task_id: int, constraints: Sequence[ScheduledConstraint]) -> None:
        for wrapped in constraints:
            constraint = wrapped.constraint
            block_set = wrapped.block_set
            block_group = wrapped.block_group
            if isinstance(constraint, DependencyConstraint):
                if constraint.depends_on_task_id == task_id:
                    raise ValueError("Task cannot depend on itself.")
                if self.db_manager.tasks.has_dependency_path(constraint.depends_on_task_id, task_id):
                    raise ValueError("Dependency cycle detected.")
                block_id = self.db_manager.tasks.add_block(
                    task_id,
                    "dependency",
                    constraint.condition,
                    constraint.message,
                    block_set=block_set,
                    block_group=block_group,
                )
                self.db_manager.tasks.add_dependency(
                    task_id,
                    constraint.depends_on_task_id,
                    constraint.required_state,
                    block_id=block_id,
                    allow_failed=constraint.allow_failed,
                )
                continue
            if isinstance(constraint, TimeConstraint):
                block_id = self.db_manager.tasks.add_block(
                    task_id,
                    "time",
                    constraint.condition,
                    constraint.message,
                    block_set=block_set,
                    block_group=block_group,
                )
                self.db_manager.tasks.add_time_constraint(
                    task_id,
                    constraint.mode,
                    constraint.not_before.isoformat(),
                    block_id=block_id,
                )
                continue
            if isinstance(constraint, BarrierConstraint):
                self.db_manager.tasks.add_block(
                    task_id,
                    "barrier",
                    constraint.condition,
                    constraint.message,
                    block_set=block_set,
                    block_group=block_group,
                )
                continue
            raise ValueError(f"Unsupported constraint type: {constraint}")

    def _expand_accession_payload(self, data_dict: Dict[str, Any]) -> Dict[str, Any]:
        """Expand accession variables in task payloads (e.g. @VAR)."""
        if not data_dict:
            return data_dict
        if "accessions" in data_dict:
            raw = data_dict.get("accessions") or []
            if isinstance(raw, str):
                raw = [raw]
            data_dict["accessions"] = expand_accession_variables(self.db_manager, raw, allow_bare=True)
        if "exclude_accessions" in data_dict:
            raw = data_dict.get("exclude_accessions") or []
            if isinstance(raw, str):
                raw = [raw]
            data_dict["exclude_accessions"] = expand_accession_variables(self.db_manager, raw, allow_bare=True)
        return data_dict

    @staticmethod
    def _coerce_bool(value: Any, default: bool = False) -> bool:
        if value is None:
            return bool(default)
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        token = str(value).strip().lower()
        if token in {"1", "true", "yes", "y", "on"}:
            return True
        if token in {"0", "false", "no", "n", "off", "none", "null", ""}:
            return False
        return bool(default)

    def _env_bool(self, key: str, default: bool = False) -> bool:
        return self._coerce_bool(self.db_manager.get_environment_variable(key), default)

    def _env_float(self, key: str, default: Optional[float] = None) -> Optional[float]:
        value = self.db_manager.get_environment_variable(key)
        if value is None or str(value).strip() == "":
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    def _env_int(self, key: str, default: int) -> int:
        value = self.db_manager.get_environment_variable(key)
        if value is None or str(value).strip() == "":
            return int(default)
        try:
            return int(value)
        except (TypeError, ValueError):
            return int(default)

    def _env_str(self, key: str, default: str) -> str:
        value = self.db_manager.get_environment_variable(key)
        token = str(value).strip() if value is not None else ""
        return token or default

    def _apply_proteome_env_defaults(
        self,
        task_key: str,
        payload: Dict[str, Any],
        explicit_fields: set[str],
        model_fields: set[str],
    ) -> Dict[str, Any]:
        """Fill omitted proteome-preparation options from project environment defaults."""

        automatic_tasks = {"download", "import-local-assembly", "batch-import-local-assembly", "verify-assembly"}
        if task_key in automatic_tasks:
            clean_isoforms = self._env_bool("DEFAULT_PROTEOME_CLEAN_ISOFORMS", True)
            defaults = {
                "clean_isoforms": clean_isoforms,
                "skip_clean_isoforms": not clean_isoforms,
                "clean_skip_gff": not self._env_bool("DEFAULT_PROTEOME_USE_GFF", True),
                "clean_skip_cdhit": not self._env_bool("DEFAULT_PROTEOME_USE_CDHIT", False),
                "clean_gff_priority": self._env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False),
                "clean_cdhit_identity": self._env_float("DEFAULT_PROTEOME_CDHIT_IDENTITY", None),
                "clean_max_concurrent": self._env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1),
                "clean_threads_per_job": self._env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1),
            }
        elif task_key == "prepare-proteome":
            defaults = {
                "skip_gff": not self._env_bool("DEFAULT_PROTEOME_USE_GFF", True),
                "skip_cdhit": not self._env_bool("DEFAULT_PROTEOME_USE_CDHIT", True),
                "gff_priority": self._env_bool("DEFAULT_PROTEOME_GFF_PRIORITY", False),
                "cdhit_identity": self._env_float("DEFAULT_PROTEOME_CDHIT_IDENTITY", None),
                "max_concurrent": self._env_int("DEFAULT_PROTEOME_MAX_CONCURRENT", 1),
                "threads_per_job": self._env_int("DEFAULT_PROTEOME_THREADS_PER_JOB", 1),
                "set_default": self._env_bool("DEFAULT_PROTEOME_SET_DEFAULT", True),
                "input_profile": self._env_str("DEFAULT_PROTEOME_INPUT_PROFILE", "raw"),
            }
        else:
            return payload

        for field, value in defaults.items():
            if field not in model_fields or field in explicit_fields or value is None:
                continue
            payload[field] = value
        return payload

    def _payload_dict_with_env_defaults(self, spec: Any, payload: Dict[str, Any] | TaskPayload) -> Dict[str, Any]:
        if isinstance(payload, TaskPayload):
            data = payload.model_dump(mode="python", exclude_none=True)
            explicit_fields = set(getattr(payload, "model_fields_set", set()))
        else:
            data = dict(payload)
            explicit_fields = set(data.keys())
        return self._apply_proteome_env_defaults(
            spec.key,
            data,
            explicit_fields,
            set(spec.payload_model.model_fields.keys()),
        )

    def queue(
        self,
        task_key: str,
        *,
        payload: Dict[str, Any] | TaskPayload,
        priority: int | None = None,
        parent_id: Optional[int] = None,
        constraints: Optional[Sequence[ScheduledConstraint]] = None,
    ) -> int:
        """Validate payload with the registry schema and queue the task."""

        self._ensure_connection()
        spec = registry.get_by_key(task_key)
        payload_data = self._payload_dict_with_env_defaults(spec, payload)
        payload_model = spec.payload_model(**payload_data)
        data_dict = payload_model.as_task_data()
        data_dict = self._expand_accession_payload(data_dict)
        constraints = [c if isinstance(c, ScheduledConstraint) else ScheduledConstraint(c) for c in (constraints or [])]
        with self.db_manager.transaction(operation=f"queue {spec.key} task"):
            task_id = self.db_manager.tasks.queue(
                job_type=spec.job_type,
                status="B" if constraints else "P",
                priority=3 if priority is None else priority,
                parent_id=parent_id,
                checkpoint=0,
                data=data_dict,
            )
            if constraints:
                self._apply_constraints(task_id, constraints)
            self._update_last_env(task_id, spec.key)
        return task_id

    def run_immediately(self, task_key: str, *, payload: Dict[str, Any] | TaskPayload, threads: Optional[int] = None) -> int:
        """Instantiate and execute the task synchronously (CLI run mode)."""

        self._ensure_connection()
        spec = registry.get_by_key(task_key)
        payload_data = self._payload_dict_with_env_defaults(spec, payload)
        model = spec.payload_model(**payload_data)
        data_dict = model.as_task_data()
        data_dict = self._expand_accession_payload(data_dict)
        model = spec.payload_model(**data_dict)
        data_json = spec.serialise_payload(model)
        preflight_task(self.db_manager, spec, data_dict)
        with self.db_manager.transaction(operation=f"prepare immediate {spec.key} task"):
            task_id = self.db_manager.tasks.queue(
                job_type=spec.job_type,
                status="H",
                priority=3,
                parent_id=None,
                checkpoint=0,
                data=data_dict,
            )
            self._update_last_env(task_id, spec.key)
        configure_logging_from_db(self.db_manager, force=True)
        runtime = refresh_runtime_thread_defaults(
            self.db_manager,
            explicit_max_threads=threads,
        )
        required_threads = (
            int(threads)
            if threads is not None
            else resolve_task_required_threads(self.db_manager, spec, runtime.max_threads)
        )
        try:
            task = spec.build_task(
                self.db_path,
                task_id=task_id,
                data_json=data_json,
                required_threads=required_threads,
                checkpoint=0,
            )
        except Exception as exc:  # boundary: registry task constructors may be third-party code
            failure = TaskExecutionError(f"Failed to instantiate task '{spec.key}': {exc}")
            with self.db_manager.transaction(operation="record task instantiation failure"):
                self.db_manager.tasks.set_error(task_id, f"Failed to instantiate task: {exc}", traceback.format_exc())
                self.db_manager.tasks.update_status(task_id, "E")
            raise failure from exc
        task.start()
        return task_id

    def describe(self, task_key: str) -> Dict[str, Any]:
        spec = registry.get_by_key(task_key)
        return {
            "key": spec.key,
            "display_name": spec.display_name or spec.key,
            "description": spec.description,
            "job_type": spec.job_type,
            "aliases": spec.aliases,
            "required_threads": spec.daemon.required_threads,
        }
