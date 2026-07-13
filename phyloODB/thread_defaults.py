"""Thread detection and task-specific thread default helpers."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Any, Callable, Mapping, Optional

from .registry import registry

DEFAULT_THREADS_PREFIX = "DEFAULT_THREADS_"
SET_MAX_THREADS_ON_START = "SET_MAX_THREADS_ON_START"


@dataclass(frozen=True)
class ThreadRuntime:
    detected_threads: int
    max_threads: int
    refreshed: bool


def _coerce_positive_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed > 0 else None


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


def detect_available_threads(
    *,
    env: Mapping[str, str] | None = None,
    cpu_count: Optional[int] = None,
    affinity_count: Optional[int] = None,
) -> int:
    """Return the best conservative estimate of available worker threads."""

    env = os.environ if env is None else env
    candidates: list[int] = []
    if affinity_count is None and hasattr(os, "sched_getaffinity"):
        try:
            affinity_count = len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
        except Exception:
            affinity_count = None
    if affinity_count and affinity_count > 0:
        candidates.append(int(affinity_count))

    for key in ("SLURM_CPUS_PER_TASK", "NSLOTS", "PBS_NP"):
        value = _coerce_positive_int(env.get(key))
        if value:
            candidates.append(value)

    slurm_nodes = str(env.get("SLURM_JOB_CPUS_PER_NODE") or "").strip()
    if slurm_nodes:
        match = re.match(r"^(\d+)", slurm_nodes)
        if match:
            value = _coerce_positive_int(match.group(1))
            if value:
                candidates.append(value)

    if cpu_count is None:
        cpu_count = os.cpu_count()
    if cpu_count and cpu_count > 0:
        candidates.append(int(cpu_count))
    return max(min(candidates), 1) if candidates else 1


def task_thread_env_name(task_key: str) -> str:
    """Build DEFAULT_THREADS_* variable name from a registry key."""

    token = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", str(task_key or "task"))
    token = re.sub(r"[^A-Za-z0-9]+", "_", token).strip("_").upper()
    return f"{DEFAULT_THREADS_PREFIX}{token or 'TASK'}"


def computed_task_thread_default(spec: Any, daemon_max_threads: int) -> int:
    """Ask the task class for its computed default thread count."""

    daemon_max_threads = max(int(daemon_max_threads or 1), 1)
    method = getattr(spec.task_cls, "default_thread_count", None)
    if callable(method):
        return max(int(method(spec.daemon.required_threads or 1, daemon_max_threads)), 1)
    return min(max(int(spec.daemon.required_threads or 1), 1), daemon_max_threads)


def computed_task_thread_defaults(daemon_max_threads: int) -> dict[str, int]:
    return {
        task_thread_env_name(spec.key): computed_task_thread_default(spec, daemon_max_threads)
        for spec in registry.specs()
    }


def _log(logger: Callable[..., Any] | None, message: str, level: str) -> None:
    if logger is None:
        return
    try:
        logger(message, level=level)
    except TypeError:
        try:
            logger(message, level)
        except TypeError:
            logger(message)


def validate_thread_cap(requested: int, detected: int, *, source: str) -> int:
    requested = max(int(requested or 1), 1)
    detected = max(int(detected or 1), 1)
    if requested > detected:
        raise ValueError(
            f"{source} requests {requested} threads, but only {detected} available thread(s) were detected."
        )
    return requested


def refresh_runtime_thread_defaults(
    manager: Any,
    *,
    explicit_max_threads: Optional[int] = None,
    logger: Callable[..., Any] | None = None,
) -> ThreadRuntime:
    """Refresh daemon and task thread defaults when configured to do so."""

    detected = detect_available_threads()
    _log(logger, f"Detected {detected} available thread(s).", "INFO")
    env = manager.get_environment_variables(["SET_MAX_THREADS_ON_START", "DAEMON_MAX_THREADS"]) or {}
    should_refresh = _coerce_bool(env.get("SET_MAX_THREADS_ON_START"), True)
    if explicit_max_threads is not None:
        effective_max = validate_thread_cap(int(explicit_max_threads), detected, source="Explicit daemon thread limit")
        if should_refresh:
            values = {
                "DAEMON_MAX_THREADS": effective_max,
                **computed_task_thread_defaults(effective_max),
            }
            manager.set_environment_variables(values, kind="env")
            _log(logger, f"Refreshed DAEMON_MAX_THREADS and task thread defaults for {effective_max} thread(s).", "INFO")
            return ThreadRuntime(detected_threads=detected, max_threads=effective_max, refreshed=True)
        return ThreadRuntime(detected_threads=detected, max_threads=effective_max, refreshed=False)

    if should_refresh:
        effective_max = detected
        values = {
            "DAEMON_MAX_THREADS": effective_max,
            **computed_task_thread_defaults(effective_max),
        }
        manager.set_environment_variables(values, kind="env")
        _log(logger, f"Refreshed DAEMON_MAX_THREADS and task thread defaults for {effective_max} thread(s).", "INFO")
        return ThreadRuntime(detected_threads=detected, max_threads=effective_max, refreshed=True)

    stored = _coerce_positive_int(env.get("DAEMON_MAX_THREADS")) or detected
    effective_max = validate_thread_cap(stored, detected, source="DAEMON_MAX_THREADS")
    return ThreadRuntime(detected_threads=detected, max_threads=effective_max, refreshed=False)


def resolve_task_required_threads(
    manager: Any,
    spec: Any,
    daemon_max_threads: int,
    *,
    logger: Callable[..., Any] | None = None,
) -> int:
    """Resolve a task's required threads from DEFAULT_THREADS_* or registry fallback."""

    key = task_thread_env_name(spec.key)
    value = manager.get_environment_variable(key)
    parsed = _coerce_positive_int(value)
    if parsed is None:
        fallback = min(max(int(spec.daemon.required_threads or 1), 1), max(int(daemon_max_threads or 1), 1))
        _log(
            logger,
            f"{key} is not set to a positive integer; falling back to registry default {fallback}.",
            "WARNING",
        )
        return fallback
    return min(parsed, max(int(daemon_max_threads or 1), 1))
