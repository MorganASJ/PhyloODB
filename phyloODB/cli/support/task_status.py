"""Helpers for turning task database rows into shell-friendly status results."""
from __future__ import annotations

from typing import Any, Dict, Optional, Sequence

from ...registry import registry

TASK_EXIT_SUCCESS = 0
TASK_EXIT_FAILED = 1
TASK_EXIT_INCOMPLETE = 2
TASK_EXIT_NOT_FOUND = 3
TASK_EXIT_CHECK_ERROR = 4
TASK_EXIT_PARTIAL_RESERVED = 10

TERMINAL_TASK_STATUSES = {"C", "E"}

STATUS_LABELS = {
    "P": "pending",
    "R": "running",
    "S": "suspended",
    "B": "blocked",
    "C": "complete",
    "E": "failed",
}


def task_exit_code(status: Any) -> int:
    """Map a raw task status to the public shell exit code contract."""

    raw_status = str(status or "").upper()
    if raw_status == "C":
        return TASK_EXIT_SUCCESS
    if raw_status == "E":
        return TASK_EXIT_FAILED
    return TASK_EXIT_INCOMPLETE


def normalize_task_state(status: Any) -> str:
    """Return a stable, script-readable state label for a raw task status."""

    raw_status = str(status or "").upper()
    return STATUS_LABELS.get(raw_status, "unknown")


def task_key_for_job_type(job_type: Any) -> Optional[str]:
    """Resolve a registry task key for a database job type when available."""

    try:
        return registry.get_by_job_type(int(job_type)).key
    except (KeyError, TypeError, ValueError):
        return None


def _json_time(value: Any) -> Optional[str]:
    if value in (None, ""):
        return None
    if hasattr(value, "isoformat"):
        return value.isoformat(sep=" ") if callable(getattr(value, "isoformat", None)) else str(value)
    return str(value)


def task_status_payload(row: Sequence[Any]) -> Dict[str, Any]:
    """Build the JSON-serialisable status payload for one Tasks row."""

    status = str(row[2] or "").upper()
    job_type = row[1]
    return {
        "task_id": int(row[0]),
        "task_key": task_key_for_job_type(job_type),
        "job_type": int(job_type) if job_type is not None else None,
        "status": status,
        "state": normalize_task_state(status),
        "normalized_state": normalize_task_state(status),
        "exit_code": task_exit_code(status),
        "queue_time": _json_time(row[7]) if len(row) > 7 else None,
        "start_time": _json_time(row[8]) if len(row) > 8 else None,
        "end_time": _json_time(row[9]) if len(row) > 9 else None,
        "status_updated_at": _json_time(row[13]) if len(row) > 13 else None,
        "error_message": row[10] if len(row) > 10 else None,
    }


__all__ = [
    "TASK_EXIT_CHECK_ERROR",
    "TASK_EXIT_FAILED",
    "TASK_EXIT_INCOMPLETE",
    "TASK_EXIT_NOT_FOUND",
    "TASK_EXIT_PARTIAL_RESERVED",
    "TASK_EXIT_SUCCESS",
    "TERMINAL_TASK_STATUSES",
    "normalize_task_state",
    "task_exit_code",
    "task_key_for_job_type",
    "task_status_payload",
]
