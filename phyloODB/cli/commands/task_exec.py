"""CLI registration and handlers for queue and run task execution."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from ...scheduling import (
    ScheduledConstraint,
    coerce_timezone,
    parse_schedule_expression,
    parse_schedule_groups,
    resolve_task_selector,
)
from ...services.task_service import TaskService
from ...task_daemon import TaskDaemon
from ...thread_defaults import refresh_runtime_thread_defaults, resolve_task_required_threads
from ..support.argparse_utils import _validate_date
from ..support.common import (
    _connect_manager,
    _format_selector_help,
    _load_selector_defaults,
    _print_error,
    _resolve_task_spec,
)
from ..support.task_parser import (
    _apply_selector_enrichment,
    _build_task_parser,
    _extract_payload,
    _load_json_payload,
)
from ..support.task_status import task_exit_code


def _merge_hidden_categories(existing: Any, *categories: str) -> str:
    merged: list[str] = []

    def add(value: str) -> None:
        normalized = str(value or "").strip().upper()
        if normalized and normalized not in merged:
            merged.append(normalized)

    if isinstance(existing, str):
        for part in existing.replace(";", ",").split(","):
            add(part)
    elif isinstance(existing, (list, tuple, set)):
        for part in existing:
            add(str(part))
    elif existing not in (None, ""):
        add(str(existing))

    for category in categories:
        add(category)
    return ",".join(merged)

def _handle_queue(args: argparse.Namespace) -> int:
    """Queue a task for daemon execution, including optional schedule constraints."""

    if args.task is None:
        return _print_error("Provide a task to queue. Use 'list tasks' to discover options.")

    task_args = list(args.task_args or [])
    # ``queue`` uses argparse.REMAINDER so task-specific options can follow the
    # task name naturally. Keep output-control flags ergonomic too:
    # both ``queue --print-id download`` and ``queue download --print-id`` work.
    if "--print-id" in task_args:
        args.print_id = True
        task_args = [token for token in task_args if token != "--print-id"]
    if "--output-json" in task_args:
        args.output_json = True
        task_args = [token for token in task_args if token != "--output-json"]

    try:
        spec = _resolve_task_spec(args.task)
    except KeyError as exc:
        return _print_error(str(exc))

    output_json = False
    payload_text = args.json_payload
    if payload_text == "__OUTPUT__":
        output_json = True
        payload_text = None
    if getattr(args, "output_json", False):
        output_json = True

    manager = _connect_manager(args.database)
    try:
        service = TaskService(args.database, db_manager=manager)
        if payload_text or args.payload_file:
            payload = _load_json_payload(payload_text, args.payload_file)
        else:
            base_prog = Path(sys.argv[0]).name or "phyloODB"
            selector_defaults = _load_selector_defaults(args.database)
            task_parser, fields = _build_task_parser(
                base_prog,
                args.database,
                "queue",
                spec,
                selector_defaults=selector_defaults,
            )
            parsed = task_parser.parse_args(task_args)
            payload = _extract_payload(parsed, fields)
            payload = _apply_selector_enrichment(manager, spec, payload, parsed)

        def resolve_existing(selector: str) -> int:
            task_id = resolve_task_selector(
                selector,
                env_lookup=lambda key: manager.get_environment_variable(key),
            )
            if manager.get_task_by_id(task_id) is None:
                raise ValueError(f"Task {task_id} not found.")
            return task_id

        parent_id = args.parent
        if args.as_subtask_of:
            if parent_id is not None:
                return _print_error("Use --parent or --as-subtask-of, not both.")
            try:
                parent_id = resolve_existing(args.as_subtask_of)
            except ValueError as exc:
                return _print_error(str(exc))

        constraints = []
        if args.schedule:
            tzinfo = coerce_timezone("local")
            for idx, token in enumerate(args.schedule):
                try:
                    if "|" in token or "&" in token:
                        groups = parse_schedule_groups(
                            token,
                            resolver=resolve_existing,
                            tzinfo=tzinfo,
                            allow_failed=False,
                        )
                        set_key = f"s{idx}"
                        for g_idx, group in enumerate(groups):
                            group_key = f"{set_key}.g{g_idx}"
                            for constraint in group:
                                constraints.append(
                                    ScheduledConstraint(
                                        constraint=constraint,
                                        block_set=set_key,
                                        block_group=group_key,
                                    )
                                )
                    else:
                        constraint = parse_schedule_expression(
                            token,
                            resolver=resolve_existing,
                            tzinfo=tzinfo,
                            allow_failed=False,
                        )
                        constraints.append(ScheduledConstraint(constraint=constraint))
                except ValueError as exc:
                    return _print_error(str(exc))

        task_id = service.queue(
            spec.key,
            payload=payload,
            priority=args.priority,
            parent_id=parent_id,
            constraints=constraints,
        )
    finally:
        manager.close()

    if args.print_id:
        print(task_id)
        return 0
    if output_json:
        status = "BLOCKED" if constraints else "QUEUED"
        response = {
            "task_id": task_id,
            "task_name": spec.key,
            "status": status,
        }
        if constraints:
            response["block_reason"] = [c.constraint.condition for c in constraints]
        print(json.dumps(response))
        return 0
    print(f"Queued task {task_id} ({spec.key}).")
    return 0


def _handle_run(args: argparse.Namespace) -> int:
    """Queue a root task, then follow its chain in a temporary foreground daemon."""

    if args.task is None:
        return _print_error("Provide a task to run. Use 'list tasks' to discover options.")

    try:
        spec = _resolve_task_spec(args.task)
    except KeyError as exc:
        return _print_error(str(exc))

    manager = _connect_manager(args.database)
    try:
        service = TaskService(args.database, db_manager=manager)
        task_args = list(args.task_args or [])
        quiet = bool(getattr(args, "quiet", False))
        show_scheduler = bool(getattr(args, "show_scheduler", False))
        if "--quiet" in task_args:
            task_args = [token for token in task_args if token != "--quiet"]
            quiet = True
        if "--show-scheduler" in task_args:
            task_args = [token for token in task_args if token != "--show-scheduler"]
            show_scheduler = True
        if args.json_payload:
            payload = _load_json_payload(args.json_payload, args.payload_file)
            parsed_threads = args.threads
        else:
            base_prog = Path(sys.argv[0]).name or "phyloODB"
            selector_defaults = _load_selector_defaults(args.database)
            task_parser, fields = _build_task_parser(
                base_prog,
                args.database,
                "run",
                spec,
                include_threads=True,
                selector_defaults=selector_defaults,
            )
            parsed = task_parser.parse_args(task_args)
            payload = _extract_payload(parsed, fields)
            payload = _apply_selector_enrichment(manager, spec, payload, parsed)
            parsed_threads = getattr(parsed, "threads", None)
            if parsed_threads is None:
                parsed_threads = getattr(parsed, "required_threads", None)
        try:
            runtime = refresh_runtime_thread_defaults(
                manager,
                explicit_max_threads=parsed_threads,
            )
        except ValueError as exc:
            return _print_error(str(exc))
        daemon_threads = runtime.max_threads
        task_threads = (
            int(parsed_threads)
            if parsed_threads is not None
            else resolve_task_required_threads(manager, spec, daemon_threads)
        )
        payload["required_threads"] = int(task_threads)
        task_id = service.queue(
            spec.key,
            payload=payload,
            priority=1,
        )
        log_overrides = {
            "LOG_TO_CONSOLE": not quiet,
        }
        if not quiet and not show_scheduler:
            existing_hidden = manager.get_environment_variable("LOG_HIDE_CATEGORIES_CONSOLE")
            log_overrides["LOG_HIDE_CATEGORIES_CONSOLE"] = _merge_hidden_categories(
                existing_hidden,
                "SCHEDULER",
            )
        daemon = TaskDaemon(
            args.database,
            data=None,
            max_threads=parsed_threads,
            root_task_id=task_id,
            scope_mode="task-chain",
            log_overrides=log_overrides,
        )
        print(f"Running task {task_id} ({spec.key}) now.")
        daemon.start()
        root_task = manager.tasks.get(task_id)
    finally:
        manager.close()
    status = (root_task[2] or "").upper() if root_task else ""
    return task_exit_code(status)



def register_task_exec_parsers(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
) -> tuple[argparse.ArgumentParser, argparse.ArgumentParser]:
    """Register the top-level ``queue`` and ``run`` commands."""

    queue_parser = subparsers.add_parser("queue", help="Queue a task for the daemon.")
    queue_options_group = queue_parser.add_argument_group("Queue options")
    queue_selector_group = queue_parser.add_argument_group("Selector options")
    queue_parser.add_argument("task", nargs="?", help="Task key or alias to queue.")
    queue_options_group.add_argument("--priority", type=int, default=3, help="Task priority (lower runs sooner).")
    queue_options_group.add_argument("--parent", type=int, default=None, help="Optional parent task id.")
    queue_options_group.add_argument(
        "--json",
        dest="json_payload",
        nargs="?",
        const="__OUTPUT__",
        help="Raw JSON payload for the task, or omit value to emit JSON output.",
    )
    queue_options_group.add_argument("--payload-file", help="Path to JSON payload file.")
    queue_options_group.add_argument("--config-path", dest="payload_file", help="Alias of --payload-file for JSON task payloads.")
    queue_options_group.add_argument(
        "--schedule",
        action="append",
        help=(
            "Release condition (repeatable). Forms: started|finished|succeeded|failed:<task>, "
            "delay:<Ns|Nm|Nh>, at:HH:MM. Combine alternatives with '|' and requirements with '&'. "
            "Task selectors include numeric ids and LAST."
        ),
    )
    queue_options_group.add_argument("--as-subtask-of", dest="as_subtask_of", help="Parent selector (task id/LAST/...)")
    queue_options_group.add_argument("--print-id", action="store_true", help="Print only the queued task id.")
    queue_options_group.add_argument("--output-json", action="store_true", help="Emit JSON output for the queued task.")
    queue_selector_group.add_argument("--busco-complete-min", type=float, help="Selector: minimum BUSCO complete proportion (0-1).")
    queue_selector_group.add_argument("--busco-single-min", type=float, help="Selector: minimum BUSCO single-copy proportion (0-1).")
    queue_selector_group.add_argument("-c", "--clade", help="Resolve a scientific name to a taxid (compatible tasks only).")
    queue_selector_group.add_argument("-rt", "--root", help="Restrict derived selectors to a specific storage root id or exact label.")
    queue_selector_group.add_argument("-i", "--taxid", type=int, help="Derive accessions or set the task taxid using a numeric identifier.")
    queue_selector_group.add_argument("-d", "--downloaded-only", action="store_true", help="Restrict derived accessions to entries already downloaded.")
    queue_selector_group.add_argument("--not-downloaded", action="store_true", help="Restrict derived accessions to entries not yet downloaded.")
    queue_selector_group.add_argument("-af", "--after", type=lambda value: _validate_date(value, "--after"), help="Limit derived accessions to assemblies released on/after YYYY-MM-DD.")
    queue_selector_group.add_argument("-bf", "--before", type=lambda value: _validate_date(value, "--before"), help="Limit derived accessions to assemblies released on/before YYYY-MM-DD.")
    queue_selector_group.add_argument("--level", choices=["complete genome", "chromosome", "scaffold", "contig"], help="Limit derived accessions to the given assembly level.")
    queue_selector_group.add_argument(
        "--primary-only",
        action="store_true",
        help=_format_selector_help(
            "Restrict selectors to primary assemblies.",
            selector_defaults or {},
            "primary_only",
            fallback=False,
        ),
    )
    queue_parser.add_argument("task_args", nargs=argparse.REMAINDER)
    queue_parser.set_defaults(handler=_handle_queue)

    run_parser = subparsers.add_parser("run", help="Run a task immediately.")
    run_options_group = run_parser.add_argument_group("Run options")
    run_selector_group = run_parser.add_argument_group("Selector options")
    run_parser.add_argument("task", nargs="?", help="Task key or alias to run.")
    run_options_group.add_argument("--json", dest="json_payload", help="Raw JSON payload for the task.")
    run_options_group.add_argument("--payload-file", help="Path to JSON payload file.")
    run_options_group.add_argument("-t", "--threads", type=int, default=None, help="Override thread count when using --json.")
    run_console_group = run_options_group.add_mutually_exclusive_group()
    run_console_group.add_argument("--quiet", action="store_true", help="Suppress console log streaming for this run.")
    run_console_group.add_argument(
        "--show-scheduler",
        action="store_true",
        help="Include scheduler lifecycle logs in the run console output.",
    )
    run_selector_group.add_argument("--busco-complete-min", type=float, help="Selector: minimum BUSCO complete proportion (0-1).")
    run_selector_group.add_argument("--busco-single-min", type=float, help="Selector: minimum BUSCO single-copy proportion (0-1).")
    run_selector_group.add_argument("-c", "--clade", help="Resolve a scientific name to a taxid (compatible tasks only).")
    run_selector_group.add_argument("-rt", "--root", help="Restrict derived selectors to a specific storage root id or exact label.")
    run_selector_group.add_argument("-i", "--taxid", type=int, help="Derive accessions or set the task taxid using a numeric identifier.")
    run_selector_group.add_argument("-d", "--downloaded-only", action="store_true", help="Restrict derived accessions to entries already downloaded.")
    run_selector_group.add_argument("--not-downloaded", action="store_true", help="Restrict derived accessions to entries not yet downloaded.")
    run_selector_group.add_argument("-af", "--after", type=lambda value: _validate_date(value, "--after"), help="Limit derived accessions to assemblies released on/after YYYY-MM-DD.")
    run_selector_group.add_argument("-bf", "--before", type=lambda value: _validate_date(value, "--before"), help="Limit derived accessions to assemblies released on/before YYYY-MM-DD.")
    run_selector_group.add_argument("--level", choices=["complete genome", "chromosome", "scaffold", "contig"], help="Limit derived accessions to the given assembly level.")
    run_selector_group.add_argument(
        "--primary-only",
        action="store_true",
        help=_format_selector_help(
            "Restrict selectors to primary assemblies.",
            selector_defaults or {},
            "primary_only",
            fallback=False,
        ),
    )
    run_parser.add_argument("task_args", nargs=argparse.REMAINDER)
    run_parser.set_defaults(handler=_handle_run)
    return queue_parser, run_parser


__all__ = ["_handle_queue", "_handle_run", "register_task_exec_parsers"]
