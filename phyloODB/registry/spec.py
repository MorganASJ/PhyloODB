"""Core definitions for the task registry architecture.

The registry centralises all metadata about runnable tasks, including:
- The numeric job_type used by the database / daemon
- A human-readable key (for CLI routing)
- The concrete Task subclass to execute in the daemon
- The Pydantic payload schema responsible for validation & serialisation
- Optional daemon configuration (threads, scheduling hints)

Later implementation steps will populate this registry for each existing task
and migrate the CLI / daemon to resolve entries from here instead of the large
if/elif blocks currently present in `phyloODB.py` and `task_daemon.py`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from importlib import import_module
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    Iterable,
    Mapping,
    MutableMapping,
    Optional,
    Protocol,
    Type,
    TypeVar,
)

from pydantic import BaseModel

from ..schemas import DaemonConfig, TaskPayload

if TYPE_CHECKING:  # typing helpers without runtime import cycle
    from ..task import Task as BaseTask
else:
    class BaseTask:  # pragma: no cover - runtime placeholder
        def __init__(self, *args: Any, **kwargs: Any) -> None: ...

TaskT = TypeVar("TaskT", bound="BaseTask")


class TaskFactory(Protocol[TaskT]):
    def __call__(
        self,
        db_path: str,
        task_id: Optional[int],
        data_json: str,
        *,
        required_threads: int,
        checkpoint: Optional[int] = None,
    ) -> TaskT:
        ...
PayloadT = TypeVar("PayloadT", bound=TaskPayload)

_BUILTIN_IO: dict[str, tuple[tuple[str, ...], bool, tuple[str, ...]]] = {
    "update-assembly": ((), True, ()),
    "download-assemblies": (("genomes",), True, ()),
    "add-library": (("libraries", "reports"), True, ()),
    "busco-run": (("genomes",), True, ()),
    "orthofinder-run": (("orthofinder",), True, ("out_dir",)),
    "download-busco-library": (("libraries",), True, ()),
    "import-local-assembly": (("genomes",), True, ()),
    "create-taxonomy": (("misc",), True, ("working_dir",)),
    "export-library": (("exports", "reports"), True, ("out_dir", "output_dir")),
    "generate-lineage-csv": (("reports",), False, ("output", "output_path")),
    "finalize-genome-move": (("genomes",), False, ()),
    "batch-import-local-assembly": (("genomes",), True, ()),
    "create-proteome-blast-db": (("genomes", "cache"), True, ("out_dir",)),
    "paralog-removal": (("cache", "reports"), True, ("report_dir",)),
    "BatchBuscoTask": (("genomes",), True, ()),
    "construct-busco-blast-db": (("cache",), True, ("out_dir",)),
    "decontamination": (("reports", "cache"), True, ("output_dir",)),
    "internal-decontamination": (("reports", "cache"), True, ("output_dir",)),
    "external-decontamination-check": (("reports", "cache"), True, ("output_dir",)),
    "external-decontamination-apply": (("reports",), True, ("output_dir",)),
    "split-records": (("misc",), True, ("output", "output_path")),
    "import-custom-library": (("libraries",), True, ("location",)),
    "prepare-proteome": (("genomes",), True, ()),
    "mafft-run": (("misc",), True, ("out_dir",)),
    "iqtree-run": (("misc",), True, ("out_dir",)),
    "build-busco-trees": (("misc",), True, ("out_dir",)),
    "annotate-orthogroup-tree": (("misc",), True, ("output", "output_path")),
}


@dataclass(slots=True)
class TaskSpec:
    """Definition for a registered task."""

    job_type: int
    key: str
    task_cls: Type[TaskT]
    payload_model: Type[PayloadT]
    description: str
    task_builder: Optional[TaskFactory[TaskT]] = None
    daemon: DaemonConfig = field(default_factory=DaemonConfig)
    aliases: tuple[str, ...] = ()
    display_name: Optional[str] = None
    queue_defaults: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)
    write_root_kinds: tuple[str, ...] = ()
    uses_scratch: bool = False
    output_path_fields: tuple[str, ...] = ()
    checkpoint_schema: Optional[Type[BaseModel]] = None
    cli_handler: Optional[str] = None
    requires_checkpoint: bool = False

    def __post_init__(self) -> None:
        defaults = _BUILTIN_IO.get(self.key)
        if defaults is None:
            return
        roots, scratch, fields = defaults
        if not self.write_root_kinds:
            self.write_root_kinds = roots
        if not self.uses_scratch:
            self.uses_scratch = scratch
        if not self.output_path_fields:
            self.output_path_fields = fields

    def serialise_payload(self, payload: Mapping[str, Any] | TaskPayload) -> str:
        """Return a JSON string representation for the database."""
        if isinstance(payload, TaskPayload):
            model = payload
        else:
            model = self.payload_model(**payload)
        return model.model_dump_json(exclude_none=True)

    def instantiate_payload(self, data: Mapping[str, Any] | str | None) -> PayloadT:
        """Build a payload model from raw data (dict or JSON string)."""
        if data is None:
            return self.payload_model()
        if isinstance(data, str):
            return self.payload_model.model_validate_json(data)
        return self.payload_model(**data)

    def build_task(
        self,
        db_path: str,
        *,
        task_id: Optional[int],
        data_json: str,
        required_threads: int | None = None,
        checkpoint: Optional[int] = None,
    ) -> TaskT:
        builder = self.task_builder
        threads = required_threads or self.daemon.required_threads
        if builder:
            return builder(
                db_path,
                task_id,
                data_json,
                required_threads=threads,
                checkpoint=checkpoint,
            )
        if self.requires_checkpoint:
            cp = checkpoint if checkpoint is not None else 0
            return self.task_cls(db_path, task_id, cp, data_json, threads)
        return self.task_cls(db_path, task_id, data_json, threads)

    def resolve_cli_handler(self) -> Optional[Callable[..., Any]]:
        """Resolve the dotted cli handler path, if provided."""
        if not self.cli_handler:
            return None
        module_name, _, attr = self.cli_handler.rpartition(".")
        if not module_name:
            raise ValueError(f"cli_handler must be a dotted path, got {self.cli_handler}")
        module = import_module(module_name)
        return getattr(module, attr)


class TaskRegistry:
    """Runtime registry that stores all TaskSpec entries."""

    def __init__(self):
        self._by_key: MutableMapping[str, TaskSpec] = {}
        self._by_job_type: MutableMapping[int, TaskSpec] = {}

    def register(self, spec: TaskSpec, *, override: bool = False) -> None:
        if not override and spec.key in self._by_key:
            raise KeyError(f"Task key '{spec.key}' already registered")
        if not override and spec.job_type in self._by_job_type:
            raise KeyError(f"Task job_type '{spec.job_type}' already registered")
        self._by_key[spec.key] = spec
        self._by_job_type[spec.job_type] = spec

    def bulk_register(self, specs: Iterable[TaskSpec], *, override: bool = False) -> None:
        for spec in specs:
            self.register(spec, override=override)

    def get_by_key(self, key: str) -> TaskSpec:
        try:
            return self._by_key[key]
        except KeyError as exc:
            raise KeyError(f"Unknown task key '{key}'") from exc

    def get_by_job_type(self, job_type: int) -> TaskSpec:
        try:
            return self._by_job_type[job_type]
        except KeyError as exc:
            raise KeyError(f"Unknown job_type '{job_type}'") from exc

    def find_by_alias(self, token: str) -> Optional[TaskSpec]:
        token_lower = token.lower()
        for spec in self._by_key.values():
            if token_lower == spec.key.lower() or token_lower in spec.aliases:
                return spec
        return None

    def specs(self) -> Iterable[TaskSpec]:
        return self._by_key.values()

    def metadata_view(self) -> Dict[str, Dict[str, Any]]:
        return {
            spec.key: {
                "job_type": spec.job_type,
                "display_name": spec.display_name or spec.key,
                "description": spec.description,
                "aliases": spec.aliases,
                "daemon": spec.daemon.model_dump(),
            }
            for spec in self._by_key.values()
        }


registry = TaskRegistry()

__all__ = ["TaskSpec", "TaskRegistry", "registry"]
