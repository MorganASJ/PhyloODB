"""Scheduling helpers shared by CLI and workflow runners."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone, time as dt_time
from typing import Callable, Optional, Union
import re

try:  # Python 3.9+
    from zoneinfo import ZoneInfo  # type: ignore
except ImportError:  # pragma: no cover - fallback
    ZoneInfo = None  # type: ignore


SCHEDULE_STATES = {"started", "finished", "succeeded", "failed"}
TIME_OF_DAY_RE = re.compile(r"^(?P<hour>[0-2]\d):(?P<minute>[0-5]\d)(?::(?P<second>[0-5]\d))?$")
DURATION_RE = re.compile(r"^(?P<value>\d+)(?P<unit>[smhSMH])$")


@dataclass(frozen=True)
class DependencyConstraint:
    depends_on_task_id: int
    required_state: str
    allow_failed: bool
    condition: str
    message: str


@dataclass(frozen=True)
class TimeConstraint:
    mode: str
    not_before: datetime
    condition: str
    message: str


@dataclass(frozen=True)
class BarrierConstraint:
    condition: str
    message: str


Constraint = Union[DependencyConstraint, TimeConstraint, BarrierConstraint]


@dataclass(frozen=True)
class ScheduledConstraint:
    constraint: Constraint
    block_set: Optional[str] = None
    block_group: Optional[str] = None


def _local_tz() -> timezone:
    return datetime.now().astimezone().tzinfo or timezone.utc


def coerce_timezone(token: Optional[str]) -> timezone:
    if token is None or str(token).lower() == "local":
        return _local_tz()
    if ZoneInfo is None:
        raise ValueError("Timezone support requires Python 3.9+ (zoneinfo).")
    try:
        return ZoneInfo(str(token))
    except Exception as exc:  # boundary: ZoneInfo may raise implementation-specific lookup errors.
        raise ValueError(f"Unknown timezone '{token}'.") from exc


def parse_duration(value: str) -> timedelta:
    match = DURATION_RE.match(value.strip())
    if not match:
        raise ValueError(f"Invalid duration '{value}'. Expected <N>s|m|h.")
    amount = int(match.group("value"))
    unit = match.group("unit").lower()
    if unit == "s":
        return timedelta(seconds=amount)
    if unit == "m":
        return timedelta(minutes=amount)
    if unit == "h":
        return timedelta(hours=amount)
    raise ValueError(f"Unsupported duration unit in '{value}'.")


def parse_time_of_day(value: str) -> dt_time:
    match = TIME_OF_DAY_RE.match(value.strip())
    if not match:
        raise ValueError(f"Invalid time-of-day '{value}'. Expected HH:MM or HH:MM:SS.")
    parts = match.groupdict()
    return dt_time(
        hour=int(parts["hour"]),
        minute=int(parts["minute"]),
        second=int(parts["second"] or 0),
    )


def parse_timestamp(value: str, tzinfo: Optional[timezone]) -> datetime:
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(f"Invalid timestamp '{value}'. Expected ISO-8601.") from exc
    if parsed.tzinfo is None:
        tz = tzinfo or _local_tz()
        parsed = parsed.replace(tzinfo=tz)
    return parsed


def next_time_of_day(target: dt_time, tzinfo: timezone) -> datetime:
    now = datetime.now(tzinfo)
    candidate = datetime.combine(now.date(), target, tzinfo)
    if candidate <= now:
        candidate = candidate + timedelta(days=1)
    return candidate


def normalize_env_key(value: str) -> str:
    return value.strip().upper().replace("-", "_")


def resolve_task_selector(
    selector: str,
    *,
    env_lookup: Callable[[str], Optional[int]],
    workflow_ids: Optional[dict[str, int]] = None,
    workflow_last: Optional[int] = None,
) -> int:
    token = selector.strip()
    if workflow_ids and token in workflow_ids:
        return workflow_ids[token]
    lowered = token.lower()
    if lowered == "last":
        if workflow_last is not None:
            return workflow_last
        env_value = env_lookup("LAST")
        if env_value is None:
            raise ValueError("LAST is not set.")
        return int(env_value)
    if token.isdigit():
        return int(token)
    env_key = normalize_env_key(token)
    env_value = env_lookup(env_key)
    if env_value is None:
        raise ValueError(f"Selector '{selector}' is not set.")
    return int(env_value)


def build_dependency_constraint(
    required_state: str,
    selector: str,
    *,
    resolver: Callable[[str], int],
    allow_failed: bool = False,
) -> DependencyConstraint:
    state = required_state.lower()
    if state not in SCHEDULE_STATES:
        raise ValueError(f"Unsupported dependency state '{required_state}'.")
    task_id = resolver(selector)
    condition = f"{state}:{task_id}"
    verb = {
        "started": "start",
        "finished": "finish",
        "succeeded": "succeed",
        "failed": "fail",
    }[state]
    message = f"Waiting for task {task_id} to {verb}."
    return DependencyConstraint(
        depends_on_task_id=task_id,
        required_state=state,
        allow_failed=allow_failed,
        condition=condition,
        message=message,
    )


def build_time_constraint_from_delay(value: str, tzinfo: timezone) -> TimeConstraint:
    duration = parse_duration(value)
    not_before = datetime.now(tzinfo) + duration
    return TimeConstraint(
        mode="delay",
        not_before=not_before,
        condition=f"delay:{value}",
        message=f"Waiting {value}.",
    )


def build_time_constraint_from_at(value: str, tzinfo: timezone) -> TimeConstraint:
    text = value.strip()
    try:
        tod = parse_time_of_day(text)
        not_before = next_time_of_day(tod, tzinfo)
        mode = "at_time"
    except ValueError:
        not_before = parse_timestamp(text, tzinfo)
        mode = "at_timestamp"
    return TimeConstraint(
        mode=mode,
        not_before=not_before,
        condition=f"at:{text}",
        message=f"Waiting until {text}.",
    )


def build_barrier_constraint(condition: str, message: str) -> BarrierConstraint:
    return BarrierConstraint(condition=condition, message=message)


def parse_schedule_expression(
    expression: str,
    *,
    resolver: Callable[[str], int],
    tzinfo: timezone,
    allow_failed: bool = False,
) -> Constraint:
    text = expression.strip()
    if not text:
        raise ValueError("Empty schedule expression.")
    if "&" in text or "|" in text:
        raise ValueError("Schedule expression contains operators; use grouped parsing.")
    lowered = text.lower()
    if lowered in {"queued-drained", "queue-drained", "drained"}:
        return build_barrier_constraint("queued-drained:global", "Waiting for queue to drain.")
    if lowered.startswith("delay:"):
        return build_time_constraint_from_delay(text.split(":", 1)[1], tzinfo)
    if lowered.startswith("at:"):
        return build_time_constraint_from_at(text.split(":", 1)[1], tzinfo)
    if ":" not in text:
        raise ValueError(f"Invalid schedule expression '{expression}'.")
    state, selector = text.split(":", 1)
    return build_dependency_constraint(state, selector, resolver=resolver, allow_failed=allow_failed)


def parse_schedule_groups(
    expression: str,
    *,
    resolver: Callable[[str], int],
    tzinfo: timezone,
    allow_failed: bool = False,
) -> list[list[Constraint]]:
    text = expression.strip()
    if not text:
        raise ValueError("Empty schedule expression.")
    or_terms = [term.strip() for term in text.split("|")]
    groups: list[list[Constraint]] = []
    for term in or_terms:
        if not term:
            raise ValueError("Invalid schedule expression near '|'.")
        and_terms = [t.strip() for t in term.split("&")]
        if not all(and_terms):
            raise ValueError("Invalid schedule expression near '&'.")
        constraints = [
            parse_schedule_expression(
                token,
                resolver=resolver,
                tzinfo=tzinfo,
                allow_failed=allow_failed,
            )
            for token in and_terms
        ]
        groups.append(constraints)
    return groups
