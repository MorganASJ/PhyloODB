from __future__ import annotations

import logging
import os
from logging.handlers import RotatingFileHandler

_LOGGER_NAME = "phyloodb"
_INITIALIZED = False
DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] [%(task_ref)s] %(message)s"
DEFAULT_LOG_FILENAME = "phyloodb.log"

# ANSI color codes for levels and task types
LEVEL_COLORS = {
    "DEBUG": "\033[36m",    # Cyan
    "INFO": "\033[32m",     # Green
    "WARNING": "\033[33m",  # Yellow
    "ERROR": "\033[31m",    # Red
    "CRITICAL": "\033[41m", # Red background
}
TASK_TYPE_COLORS = {
    "analysis": "\033[35m",   # Magenta
    "import": "\033[34m",     # Blue
    "export": "\033[36m",     # Cyan
    # Add more task types as needed
}
RESET = "\033[0m"


def _normalize_category(value) -> str:
    text = str(value or "").strip().upper()
    return text or "TASK"


def _parse_hidden_categories(value) -> set[str]:
    if value in (None, "", [], (), set(), {}):
        return set()
    if isinstance(value, str):
        raw_items = [part.strip() for part in value.replace(";", ",").split(",")]
    elif isinstance(value, (list, tuple, set)):
        raw_items = [str(part).strip() for part in value]
    else:
        raw_items = [str(value).strip()]
    return {_normalize_category(item) for item in raw_items if str(item).strip()}


def _task_ref(task_name, task_id, stage) -> str:
    label = str(task_name or "Task").strip() or "Task"
    if task_id in (None, ""):
        return label
    stage_value = 0 if stage in (None, "") else stage
    return f"{label}:{task_id}.{stage_value}"


class _CategoryFilter(logging.Filter):
    def __init__(self, hidden_categories: set[str] | None = None):
        super().__init__()
        self.hidden_categories = set(hidden_categories or set())

    def filter(self, record: logging.LogRecord) -> bool:
        category = _normalize_category(getattr(record, "log_category", "TASK"))
        return category not in self.hidden_categories

class _SafeExtraFormatter(logging.Formatter):
    """Formatter that tolerates missing extra fields used in the format string."""

    def format(self, record: logging.LogRecord) -> str:
        for key in ("task_id", "task_type", "task_name", "stage", "log_category"):
            if not hasattr(record, key):
                setattr(record, key, None)
        record.log_category = _normalize_category(getattr(record, "log_category", "TASK"))
        record.task_ref = _task_ref(
            getattr(record, "task_name", None),
            getattr(record, "task_id", None),
            getattr(record, "stage", None),
        )
        return super().format(record)

class _ColorFormatter(_SafeExtraFormatter):
    def __init__(self, fmt: str):
        super().__init__(fmt)

    def format(self, record: logging.LogRecord) -> str:
        # Color levelname
        color = LEVEL_COLORS.get(record.levelname, "")
        record.levelname = f"{color}{record.levelname}{RESET}"

        # Color task_type if present
        if hasattr(record, "task_type") and record.task_type:
            tcolor = TASK_TYPE_COLORS.get(str(record.task_type).lower(), "")
            record.task_type = f"{tcolor}{record.task_type}{RESET}"
        return super().format(record)


def init_logging(
    level: str = "INFO",
    to_console: bool = True,
    log_file: str | None = None,
    max_bytes: int = 10 * 1024 * 1024,
    backups: int = 3,
    fmt: str | None = None,
    use_color: bool = True,
    hide_console_categories: set[str] | None = None,
    hide_file_categories: set[str] | None = None,
    force: bool = False,
) -> None:
    """Initialize the package logger once.

    - level: logging level name
    - to_console: whether to log to stderr
    - log_file: optional rotating file path
    - max_bytes/backups: rotation controls
    """
    global _INITIALIZED
    logger = logging.getLogger(_LOGGER_NAME)
    if _INITIALIZED and not force:
        return
    if force:
        for handler in list(logger.handlers):
            logger.removeHandler(handler)

    logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    format_str = fmt or DEFAULT_FORMAT
    formatter_cls = _ColorFormatter if use_color else _SafeExtraFormatter
    formatter = formatter_cls(format_str)

    if to_console:
        ch = logging.StreamHandler()
        ch.setLevel(getattr(logging, level.upper(), logging.INFO))
        ch.setFormatter(formatter)
        if hide_console_categories:
            ch.addFilter(_CategoryFilter(hide_console_categories))
        logger.addHandler(ch)

    if log_file:
        try:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
        except OSError:
            # If dirname is empty or cannot be created, let handler raise later if invalid
            logger.debug("Could not create log directory for %s; handler creation will report invalid paths.", log_file)
        fh = RotatingFileHandler(log_file, maxBytes=max_bytes, backupCount=backups)
        fh.setLevel(getattr(logging, level.upper(), logging.INFO))
        fh.setFormatter(formatter)
        if hide_file_categories:
            fh.addFilter(_CategoryFilter(hide_file_categories))
        logger.addHandler(fh)

    _INITIALIZED = True


def _coerce_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def resolve_log_file_from_db(dbm, overrides: dict | None = None) -> str | None:
    """Resolve the effective log file path, preferring explicit overrides then logs roots."""

    overrides = overrides or {}
    override_log_file = overrides.get("LOG_FILE")
    if override_log_file:
        return os.path.abspath(os.path.expanduser(str(override_log_file)))
    try:
        root_base = dbm.storage.get_root_base("logs", ensure_from_env=False)
    except Exception:  # boundary: logging path resolution falls back when storage metadata is unavailable.
        root_base = None
    if root_base:
        return os.path.abspath(os.path.join(str(root_base), DEFAULT_LOG_FILENAME))
    try:
        cfg = dbm.get_environment_variables(["LOG_DIR", "LOG_FILE"]) or {}
    except Exception:  # boundary: logging path resolution falls back when environment settings are unavailable.
        cfg = {}
    log_dir = cfg.get("LOG_DIR")
    if log_dir:
        return os.path.abspath(os.path.join(os.path.expanduser(str(log_dir)), DEFAULT_LOG_FILENAME))
    log_file = cfg.get("LOG_FILE")
    if log_file:
        return os.path.abspath(os.path.expanduser(str(log_file)))
    return None


def configure_logging_from_db(dbm, overrides: dict | None = None, *, force: bool = False) -> None:
    """Configure logging using Environment_Variables if present.

    Expected keys (all optional): LOG_LEVEL, LOG_TO_CONSOLE, LOG_DIR, LOG_FILE, LOG_MAX_BYTES, LOG_BACKUPS
    Fallbacks: INFO, True, None, 10MB, 3
    """
    try:
        cfg = dbm.get_environment_variables([
            'LOG_LEVEL',
            'LOG_TO_CONSOLE',
            'LOG_DIR',
            'LOG_FILE',
            'LOG_MAX_BYTES',
            'LOG_BACKUPS',
            'LOG_HIDE_CATEGORIES_CONSOLE',
            'LOG_HIDE_CATEGORIES_FILE',
        ])
    except Exception:  # boundary: logging configuration must fall back safely if DB settings are unavailable.
        cfg = {}

    merged = dict(cfg)
    if overrides:
        merged.update({k: v for k, v in overrides.items() if v is not None})

    level = merged.get('LOG_LEVEL', 'INFO')
    to_console = _coerce_bool(merged.get('LOG_TO_CONSOLE', True))
    log_file = resolve_log_file_from_db(dbm, overrides=overrides)
    max_bytes = merged.get('LOG_MAX_BYTES', 10 * 1024 * 1024)
    backups = merged.get('LOG_BACKUPS', 3)
    fmt = merged.get('LOG_FORMAT', DEFAULT_FORMAT)
    use_color = _coerce_bool(merged.get('LOG_USE_COLOR', True))
    hide_console_categories = _parse_hidden_categories(merged.get('LOG_HIDE_CATEGORIES_CONSOLE'))
    hide_file_categories = _parse_hidden_categories(merged.get('LOG_HIDE_CATEGORIES_FILE'))

    init_logging(
        level,
        to_console,
        log_file,
        int(max_bytes),
        int(backups),
        fmt=fmt,
        use_color=use_color,
        hide_console_categories=hide_console_categories,
        hide_file_categories=hide_file_categories,
        force=force,
    )


def get_task_logger(task_id=None, task_type: str | None = None, task_name: str | None = None):
    """Return a LoggerAdapter injecting task context into records."""
    base = logging.getLogger(_LOGGER_NAME)

    class _Adapter(logging.LoggerAdapter):
        def process(self, msg, kwargs):
            extra = kwargs.get('extra', {})
            extra.setdefault('task_id', task_id)
            extra.setdefault('task_type', task_type)
            extra.setdefault('task_name', task_name)
            # stage can be passed per-call via extra if available
            kwargs['extra'] = extra
            return msg, kwargs

    return _Adapter(base, {})
