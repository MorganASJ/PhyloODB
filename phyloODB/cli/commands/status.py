"""CLI handler for shell-friendly task status checks."""
from __future__ import annotations

import argparse
import json
import sys
import time
from typing import Any, Optional

from ...scheduling import resolve_task_selector
from ..support.common import _connect_manager
from ..support.task_status import (
    TASK_EXIT_CHECK_ERROR,
    TASK_EXIT_NOT_FOUND,
    TERMINAL_TASK_STATUSES,
    task_status_payload,
)


def _print_status_error(args: argparse.Namespace, message: str, code: int) -> int:
    if not getattr(args, "quiet", False):
        print(f"Error: {message}", file=sys.stderr)
    return code


def _resolve_status_task_id(manager: Any, selector: str) -> int:
    return resolve_task_selector(
        selector,
        env_lookup=lambda key: manager.get_environment_variable(key),
    )


def _render_status(payload: dict[str, Any], *, output_json: bool) -> None:
    if output_json:
        print(json.dumps(payload, sort_keys=True))
        return
    task_name = payload.get("task_key") or f"job_type={payload.get('job_type')}"
    print(
        f"Task {payload['task_id']} ({task_name}) "
        f"is {payload['state']} [{payload['status']}]; exit_code={payload['exit_code']}"
    )
    if payload.get("error_message"):
        print(f"Error: {payload['error_message']}")


def _load_status_payload(args: argparse.Namespace) -> tuple[Optional[dict[str, Any]], int]:
    try:
        manager = _connect_manager(args.database, read_only=True)
    except Exception as exc:
        return None, _print_status_error(args, f"failed to open database: {exc}", TASK_EXIT_CHECK_ERROR)
    try:
        try:
            task_id = _resolve_status_task_id(manager, args.task)
        except (TypeError, ValueError) as exc:
            return None, _print_status_error(args, str(exc), TASK_EXIT_NOT_FOUND)
        row = manager.tasks.get(task_id)
        if row is None:
            return None, _print_status_error(args, f"Task {task_id} not found.", TASK_EXIT_NOT_FOUND)
        return task_status_payload(row), 0
    except Exception as exc:
        return None, _print_status_error(args, f"failed to check task status: {exc}", TASK_EXIT_CHECK_ERROR)
    finally:
        manager.close()


def _handle_status(args: argparse.Namespace) -> int:
    if args.task is None:
        return _print_status_error(args, "Provide a task id or selector.", TASK_EXIT_NOT_FOUND)
    if args.interval <= 0:
        return _print_status_error(args, "--interval must be greater than zero.", TASK_EXIT_CHECK_ERROR)
    if args.timeout is not None and args.timeout < 0:
        return _print_status_error(args, "--timeout must be zero or greater.", TASK_EXIT_CHECK_ERROR)

    deadline = time.monotonic() + args.timeout if args.timeout is not None else None
    while True:
        payload, load_code = _load_status_payload(args)
        if payload is None:
            return load_code
        status = str(payload.get("status") or "").upper()
        if not args.wait or status in TERMINAL_TASK_STATUSES:
            if not args.quiet:
                _render_status(payload, output_json=bool(args.output_json))
            return int(payload["exit_code"])
        if deadline is not None and time.monotonic() >= deadline:
            if not args.quiet:
                _render_status(payload, output_json=bool(args.output_json))
            return int(payload["exit_code"])
        sleep_for = args.interval
        if deadline is not None:
            sleep_for = min(sleep_for, max(0.0, deadline - time.monotonic()))
        if sleep_for <= 0:
            continue
        time.sleep(sleep_for)


def register_status_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "status",
        help="Check a queued task status and return a shell-friendly exit code.",
    )
    parser.add_argument("task", nargs="?", help="Task id or selector such as LAST or LAST_DOWNLOAD.")
    parser.add_argument("--json", dest="output_json", action="store_true", help="Emit machine-readable JSON.")
    parser.add_argument("--quiet", action="store_true", help="Print nothing; use only the process exit code.")
    parser.add_argument("--wait", action="store_true", help="Poll until the task reaches a terminal state.")
    parser.add_argument("--interval", type=float, default=2.0, help="Polling interval in seconds for --wait.")
    parser.add_argument("--timeout", type=float, default=None, help="Maximum seconds to wait before returning incomplete.")
    parser.set_defaults(handler=_handle_status)
    return parser


__all__ = ["_handle_status", "register_status_parser"]
