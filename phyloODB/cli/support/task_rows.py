"""Renderer-neutral row collectors for task queue and error lists."""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from typing import Any, Dict, Mapping, Sequence

from ...database import DBManager
from ...registry import registry

QUEUE_HEADERS = (
    "Task ID",
    "Task Name",
    "Priority",
    "C",
    "Status",
    "Why",
    "Queue Time",
    "Start Time",
    "End Time",
)

STATUS_STYLE = {
    "P": "red",
    "R": "bold green",
    "C": "grey50",
    "S": "bold blue",
    "E": "bold yellow",
    "B": "dim magenta",
}

ACTIVE_STATUSES = {"R", "P", "S", "B"}
STATUS_ORDER = {"R": 0, "P": 1, "S": 2, "B": 3, "E": 4, "C": 5}
QUEUE_SORT_ALIASES = {
    "changed": "latest",
    "newest": "new",
    "oldest": "old",
    "active": "running",
}


def _job_display_name(job_type: int) -> str:
    if job_type == 0:
        return "TaskDaemon"
    if job_type == -1:
        return "Barrier"
    try:
        spec = registry.get_by_job_type(job_type)
    except KeyError:
        return f"job_type[{job_type}]"
    return spec.display_name or spec.key


def _fmt_time(value: Any) -> str:
    if value in (None, ""):
        return "N/A"
    if hasattr(value, "strftime"):
        return value.strftime("%H:%M:%S")
    text = str(value).strip()
    if not text:
        return "N/A"
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00")).strftime("%H:%M:%S")
    except ValueError:
        pass
    for pattern in (
        "%Y-%m-%d %H:%M:%S.%f",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%H:%M:%S",
        "%H:%M",
    ):
        try:
            return datetime.strptime(text, pattern).strftime("%H:%M:%S")
        except ValueError:
            continue
    return text


def _task_status_updated(task: Sequence[Any]) -> str:
    if len(task) > 13 and task[13]:
        return str(task[13])
    return str(task[9] or task[8] or task[7] or "")


def _sort_time(value: Any) -> float:
    text = str(value or "").strip()
    if not text:
        return 0.0
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        for pattern in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%H:%M:%S", "%H:%M"):
            try:
                parsed = datetime.strptime(text, pattern)
                break
            except ValueError:
                continue
        else:
            return 0.0
    if parsed.tzinfo is not None:
        parsed = parsed.astimezone(timezone.utc).replace(tzinfo=None)
    return parsed.timestamp()


def _parse_not_before(value: Any) -> datetime:
    text = str(value).strip().replace("Z", "+00:00")
    parsed = datetime.fromisoformat(text)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _block_summary_map(db: DBManager, tasks: Sequence[Any]) -> Dict[int, str]:
    task_ids = [task[0] for task in tasks]
    if not task_ids:
        return {}
    blocks = db.get_task_blocks(task_ids, unsatisfied_only=True)
    deps = db.get_task_dependencies(task_ids)
    times = db.get_task_time_constraints(task_ids)
    dep_by_block = {row[4]: row for row in deps if row[4] is not None}
    time_by_block = {row[4]: row for row in times if row[4] is not None}
    now = datetime.now(timezone.utc)
    summary: Dict[int, list[str]] = {}
    for block in sorted(blocks, key=lambda row: row[0]):
        block_id, task_id, block_type, condition, _message, *_ = block
        code = ""
        if block_type == "dependency" and block_id in dep_by_block:
            depends_on = dep_by_block[block_id]
            state = (depends_on[3] or "").lower()
            prefix = {"started": "s", "finished": "f", "succeeded": "ok", "failed": "x"}.get(state, "d")
            code = f"{prefix}:{depends_on[2]}"
        elif block_type == "time" and block_id in time_by_block:
            time_row = time_by_block[block_id]
            if time_row[2] == "delay":
                try:
                    code = f"d:{max(0, int((_parse_not_before(time_row[3]) - now).total_seconds()))}"
                except (TypeError, ValueError):
                    code = "d"
            else:
                code = f"t:{condition.split(':', 1)[1]}" if condition and condition.startswith("at:") else "t"
        elif block_type == "barrier" and condition and condition.startswith("queued-drained"):
            code = "q"
        if code:
            summary.setdefault(task_id, []).append(code)
    return {
        task_id: "+".join(codes) if len(codes) <= 2 else f"{codes[0]}+"
        for task_id, codes in summary.items()
    }


def collect_queue_rows(
    db: DBManager,
    options: Mapping[str, Any],
) -> tuple[list[tuple[str, ...]], list[str]]:
    statuses = {str(status).upper() for status in options.get("status", [])}
    sort_profile = QUEUE_SORT_ALIASES.get(str(options.get("sort") or "latest").lower(), str(options.get("sort") or "latest").lower())
    simple = bool(options.get("simple"))
    hide_complete = bool(options.get("hide_complete"))
    hide_done = bool(options.get("hide_done"))
    try:
        all_tasks = db.get_tasks()
    except sqlite3.OperationalError as exc:
        if "locked" in str(exc).lower():
            return [("Database is locked; waiting...", "", "", "", "", "", "", "", "")], ["yellow"]
        raise
    tasks = [task for task in all_tasks if not statuses or (task[2] or "").upper() in statuses]
    if not tasks:
        return [("No tasks found.", "", "", "", "", "", "", "", "")], ["bold yellow"]

    status_by_id = {int(task[0]): (task[2] or "").upper() for task in all_tasks}
    parent_by_id = {int(task[0]): task[4] for task in all_tasks}
    all_by_id = {int(task[0]): task for task in all_tasks}
    ancestor_cache: Dict[int, bool] = {}
    visiting: set[int] = set()

    def ancestors_complete(task_id: int) -> bool:
        if task_id in ancestor_cache:
            return ancestor_cache[task_id]
        if task_id in visiting:
            return False
        visiting.add(task_id)
        parent = parent_by_id.get(task_id)
        result = parent is None or (status_by_id.get(parent) == "C" and ancestors_complete(parent))
        visiting.remove(task_id)
        ancestor_cache[task_id] = result
        return result

    all_children: Dict[int, list[Any]] = {}
    for task in all_tasks:
        parent = task[4]
        if parent is not None:
            all_children.setdefault(int(parent), []).append(task)

    tree_cache: Dict[int, Dict[str, Any]] = {}
    tree_visiting: set[int] = set()

    def tree_metrics(task_id: int) -> Dict[str, Any]:
        if task_id in tree_cache:
            return tree_cache[task_id]
        task = all_by_id.get(task_id)
        if task is None or task_id in tree_visiting:
            return {"updated": "", "has_error": False, "has_active": False, "status": ""}
        tree_visiting.add(task_id)
        status = status_by_id.get(task_id, "")
        updated = _task_status_updated(task)
        has_error = status == "E"
        has_active = status in ACTIVE_STATUSES
        for child in all_children.get(task_id, []):
            child_metrics = tree_metrics(int(child[0]))
            if child_metrics["updated"] > updated:
                updated = child_metrics["updated"]
            has_error = has_error or bool(child_metrics["has_error"])
            has_active = has_active or bool(child_metrics["has_active"])
        tree_visiting.remove(task_id)
        tree_cache[task_id] = {
            "updated": updated,
            "has_error": has_error,
            "has_active": has_active,
            "status": status,
        }
        return tree_cache[task_id]

    def sort_key(task: Sequence[Any]):
        task_id = int(task[0])
        metrics = tree_metrics(task_id)
        status = status_by_id.get(task_id, "")
        updated = _sort_time(metrics["updated"])
        if sort_profile == "new":
            return (-task_id,)
        if sort_profile == "old":
            return (task_id,)
        if sort_profile == "errors":
            return (0 if metrics["has_error"] else 1, -updated, -task_id)
        if sort_profile == "running":
            return (0 if metrics["has_active"] else 1, -updated, -task_id)
        if sort_profile == "status":
            return (STATUS_ORDER.get(status, 99), -updated, -task_id)
        return (-updated, -task_id)

    children: Dict[int, list[Any]] = {}
    for task in tasks:
        ancestors_complete(task[0])
        if task[4] is not None:
            children.setdefault(task[4], []).append(task)
    all_ids = {task[0] for task in tasks}
    roots = sorted(
        [task for task in tasks if task[4] is None or task[4] not in all_ids],
        key=sort_key,
    )
    blocks = _block_summary_map(db, tasks)
    rows: list[tuple[str, ...]] = []
    styles: list[str] = []
    displayed: set[int] = set()

    def add(task: Sequence[Any], indent: int = 0) -> None:
        task_id = task[0]
        if task_id in displayed:
            return
        displayed.add(task_id)
        status = (task[2] or "").upper()
        hidden = (
            (hide_done and status in {"C", "E"})
            or (hide_complete and status == "C")
            or (simple and status == "C" and ancestor_cache.get(task_id, False))
        )
        if not hidden:
            parent = task[4]
            task_id_text = str(task_id) if parent is None else f"{task_id} ({parent})"
            rows.append(
                (
                    f"{'  ' * indent}{task_id_text}",
                    _job_display_name(task[1]),
                    str(task[3]),
                    str(task[5] if task[5] is not None else ""),
                    status or "?",
                    blocks.get(task_id, ""),
                    _fmt_time(task[7]),
                    _fmt_time(task[8]),
                    _fmt_time(task[9]),
                )
            )
            styles.append(STATUS_STYLE.get(status, ""))
        for child in sorted(children.get(task_id, []), key=sort_key):
            add(child, indent + 1)

    for root in roots:
        add(root)
    if not rows:
        return [("No tasks to display (filtered).", "", "", "", "", "", "", "", "")], ["bold yellow"]
    return rows, styles


def collect_error_rows(
    db: DBManager,
    options: Mapping[str, Any],
) -> tuple[tuple[str, ...], list[tuple[str, ...]]]:
    include_stack = bool(options.get("include_stack"))
    headers = (
        ("Task ID", "Task Name", "Error Message", "Error Stack")
        if include_stack
        else ("Task ID", "Task Name", "Error Message")
    )
    sql = """
        SELECT task_id, job_type, error_message, error_stack
        FROM Tasks
        WHERE (error_message IS NOT NULL AND TRIM(error_message) != '')
           OR (error_stack IS NOT NULL AND TRIM(error_stack) != '')
        ORDER BY task_id DESC
    """
    params: tuple[Any, ...] = ()
    if options.get("limit") is not None:
        sql += " LIMIT ?"
        params = (int(options["limit"]),)
    db.cursor.execute(sql, params)
    raw_rows = db.cursor.fetchall()
    if not raw_rows:
        empty = ("No task errors recorded.", "", "", "") if include_stack else ("No task errors recorded.", "", "")
        return headers, [empty]
    rows = []
    for task_id, job_type, message, stack in raw_rows:
        base = (str(task_id), _job_display_name(job_type), str(message or ""))
        rows.append((*base, str(stack or "")) if include_stack else base)
    return headers, rows
