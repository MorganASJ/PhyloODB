"""Central task registry exported for CLI and daemon wiring."""

from .spec import TaskRegistry, TaskSpec, registry
from .specs import builtin_task_specs


def register_builtin_specs() -> None:
    for spec in builtin_task_specs():
        registry.register(spec, override=True)


# Register built-in specs on import so CLI/daemon can rely on populated registry.
register_builtin_specs()

__all__ = ["TaskRegistry", "TaskSpec", "registry", "register_builtin_specs"]
