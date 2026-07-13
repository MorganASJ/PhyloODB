"""Discovery helpers for populating the task registry."""
from __future__ import annotations

import importlib
from importlib.metadata import entry_points
from typing import Iterable, Sequence

from .spec import TaskRegistry, TaskSpec


ENTRY_POINT_GROUP = "phyloodb.tasks"


def load_entry_point_specs(registry: TaskRegistry) -> None:
    """Load TaskSpec providers from Python entry points.

    Third-party packages will expose a callable that returns an iterable of
    TaskSpec instances. The callable is executed lazily so import side-effects
    remain under control.
    """

    eps_collection = entry_points()
    if hasattr(eps_collection, "select"):
        eps: Sequence = eps_collection.select(group=ENTRY_POINT_GROUP)
    else:  # pragma: no cover - Python <3.10 compatibility
        eps = eps_collection.get(ENTRY_POINT_GROUP, [])
    for ep in eps:
        factory = ep.load()
        specs = factory() if callable(factory) else factory
        if isinstance(specs, TaskSpec):
            registry.register(specs, override=True)
        elif isinstance(specs, Iterable):
            registry.bulk_register(specs, override=True)
        else:
            raise TypeError(f"Entry point {ep.name} returned unsupported object {type(specs)!r}")


def load_module_specs(registry: TaskRegistry, dotted_paths: Iterable[str]) -> None:
    """Import dotted callables that return TaskSpec iterables."""

    for dotted in dotted_paths:
        module_name, _, attr = dotted.rpartition(".")
        if not module_name:
            raise ValueError(f"Expected module path for task spec provider, got {dotted}")
        module = importlib.import_module(module_name)
        factory = getattr(module, attr)
        specs = factory()
        if isinstance(specs, TaskSpec):
            registry.register(specs, override=True)
        else:
            registry.bulk_register(specs, override=True)
