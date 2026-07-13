"""Shared non-command-specific helpers for the PhyloODB CLI."""
from __future__ import annotations

import os
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from pydantic import BaseModel

from ...database import DBManager
from ...registry import TaskSpec, registry
from ...selector_utils import expand_accession_variables


# ---------------------------------------------------------------------------
# Shared constants
# Purpose: Keep database-backed CLI defaults and storage-kind constants in one
# importable place for parsers, handlers, and entrypoint bootstrap.
# ---------------------------------------------------------------------------

SELECTOR_DEFAULT_ENV_KEYS = {
    "downloaded_only": "SELECTOR_DEFAULT_DOWNLOADED_ONLY",
    "primary_only": "SELECTOR_DEFAULT_PRIMARY_ONLY",
    "use_busco": "SELECTOR_DEFAULT_USE_BUSCO",
    "status_min": "SELECTOR_DEFAULT_STATUS_MIN",
    "protein_only": "SELECTOR_DEFAULT_PROTEIN_ONLY",
}

LIST_COLOR_ENV_KEYS = [
    "LIST_USE_COLOR",
    "LIST_GROUP_COLORS",
    "LIST_BUSCO_GRADIENT",
    "LIST_BUSCO_GRADIENT_POS",
    "LIST_BUSCO_GRADIENT_NEG",
    "LIST_BUSCO_STEEP_MAX",
    "LIST_BUSCO_POS_STOPS",
    "LIST_BUSCO_NEG_STOPS",
    "LIST_BUSCO_STEEP_STOPS",
]

STORAGE_ROOT_KINDS = ["genomes", "libraries", "orthofinder", "exports", "reports", "logs", "cache", "misc"]
STRICT_WORKING_ROOT_KINDS = {"genomes", "libraries", "orthofinder", "exports", "logs"}


# ---------------------------------------------------------------------------
# Bootstrap and environment helpers
# Purpose: Load database-backed defaults before parser construction and keep
# low-level environment coercion out of the command modules.
# ---------------------------------------------------------------------------

def _infer_db_path(args: Sequence[str]) -> Optional[str]:
    """Infer the database path from raw argv before the parser is built."""

    if not args:
        return None
    candidate = args[0]
    if candidate.startswith("-"):
        return None
    return candidate


def _load_selector_defaults(db_path: Optional[str]) -> Dict[str, Any]:
    """Load selector default values from the database environment table."""

    if not db_path:
        return {}
    if not os.path.exists(db_path):
        return {}
    manager = DBManager(db_path, read_only=True)
    try:
        manager.connect()
        if not manager._table_exists("Environment_Variables"):
            return {}
        env_keys = list(SELECTOR_DEFAULT_ENV_KEYS.values())
        return manager.get_environment_variables(env_keys) or {}
    finally:
        manager.close()


def _load_list_color_defaults(manager: DBManager) -> Dict[str, Any]:
    """Load list colour preferences from database-backed environment values."""

    return manager.get_environment_variables(LIST_COLOR_ENV_KEYS) or {}


def _coerce_bool(value: Any, default: bool = False) -> bool:
    """Coerce loosely typed CLI or environment values into a boolean."""

    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _format_selector_help(help_text: str, defaults: Mapping[str, Any], field: str, fallback: Any = None) -> str:
    """Append selector default information to help text when available."""

    env_key = SELECTOR_DEFAULT_ENV_KEYS.get(field)
    if not env_key:
        return help_text
    if defaults and env_key in defaults:
        val = defaults.get(env_key)
        return f"{help_text} (default: {val}, {env_key}={val})"
    if fallback is None:
        return f"{help_text}"
    return f"{help_text} (default: {fallback})"


# ---------------------------------------------------------------------------
# Registry and database convenience helpers
# Purpose: Centralise the small DB/task-registry lookups used by multiple
# command modules so those modules stay focused on command flow.
# ---------------------------------------------------------------------------

def _normalize_action_alias(argv: Sequence[str]) -> List[str]:
    """Expand short action aliases before argparse sees the command."""

    tokens = list(argv)
    if len(tokens) >= 2:
        alias = tokens[1]
        mapping = {"-q": "queue", "-r": "run"}
        replacement = mapping.get(alias)
        if replacement:
            tokens[1] = replacement
    if len(tokens) >= 3 and tokens[1] == "watch" and tokens[2] in {"queue", "errors"}:
        tokens[1] = "list"
        tokens.insert(3, "--watch")
    idx = 2
    while idx < len(tokens):
        token = tokens[idx]
        if token in {"-q", "-r"}:
            values: List[str] = []
            probe = idx + 1
            while probe < len(tokens):
                candidate = str(tokens[probe])
                if candidate.startswith("-"):
                    break
                values.append(candidate)
                probe += 1
            plural = len(values) > 1 or any("," in value for value in values)
            if token == "-q":
                tokens[idx] = "--quantities" if plural else "--quantity"
            else:
                tokens[idx] = "--ranks" if plural else "--rank"
        idx += 1
    if (
        len(tokens) >= 3
        and tokens[1] == "set"
        and tokens[2] not in {"var", "env", "busco-primary", "proteome-profile"}
        and not str(tokens[2]).startswith("-")
    ):
        tokens.insert(2, "var")
        tokens.insert(3, "--legacy-syntax")
    return tokens


def _connect_manager(db_path: str, *, read_only: bool = False) -> DBManager:
    """Open and return a connected database manager for the target database."""

    manager = DBManager(db_path, read_only=read_only)
    manager.connect()
    return manager


def _apply_busco_context_from_args(manager: DBManager, args: Any) -> None:
    """Apply BUSCO run selection context from parsed CLI args to the DB manager."""

    manager.set_busco_run_context(
        pipeline=getattr(args, "busco_pipeline", None),
        input_mode=getattr(args, "busco_input_mode", None) or getattr(args, "format", None),
        prefer_pipeline=getattr(args, "prefer_busco_pipeline", None),
        prefer_input_mode=getattr(args, "prefer_busco_input_mode", None) or getattr(args, "prefer_format", None),
        proteome_profile=getattr(args, "proteome_profile", None),
        prefer_proteome_profile=getattr(args, "prefer_proteome_profile", None),
        run_ids=getattr(args, "busco_run_ids", None) or getattr(args, "run_ids", None),
        selection=getattr(args, "busco_run_selection", None),
    )


def _expand_accessions(manager: DBManager, accessions: Sequence[Any]) -> List[str]:
    """Expand literal accessions and stored accession variables into a flat list."""

    return expand_accession_variables(manager, accessions or [], allow_bare=True)


def _resolve_task_spec(token: str) -> TaskSpec:
    """Resolve a task key or alias to its registered task specification."""

    spec = registry.find_by_alias(token)
    if spec:
        return spec
    return registry.get_by_key(token)


def _resolve_library_id(
    manager: DBManager,
    value: Optional[Union[int, str]],
    *,
    option: str,
) -> Optional[int]:
    """Resolve a BUSCO library selector given either an id or a library name."""

    if value is None:
        return None
    if isinstance(value, int):
        return value
    token = str(value).strip()
    if not token:
        return None
    if token.isdigit():
        return int(token)
    lib_id = manager.get_library_id(token)
    if lib_id is None:
        raise ValueError(f"{option} expects a library id or name; unknown library '{token}'.")
    return lib_id


def _resolve_library_selector(
    manager: DBManager,
    *,
    library_id: Optional[int],
    library_name: Optional[str],
    legacy: Optional[Union[int, str]] = None,
) -> Optional[int]:
    """Resolve the preferred library selector arguments into a single library id."""

    if library_id is not None and library_name is not None:
        resolved = _resolve_library_id(manager, library_name, option="--library-name")
        if int(library_id) != int(resolved):
            raise ValueError("--library-id and --library-name refer to different libraries.")
        return int(library_id)
    if library_id is not None:
        return int(library_id)
    if library_name is not None:
        return _resolve_library_id(manager, library_name, option="--library-name")
    if legacy is not None:
        return _resolve_library_id(manager, legacy, option="--busco-library")
    return None


def _resolve_library_selector_alias(
    manager: DBManager,
    *,
    library_id: Optional[int],
    library_name: Optional[str],
    legacy: Optional[Union[int, str]] = None,
    library: Optional[Union[int, str]] = None,
) -> Optional[int]:
    """Resolve library selectors while honouring the older ``--library`` alias."""

    resolved = _resolve_library_selector(
        manager,
        library_id=library_id,
        library_name=library_name,
        legacy=legacy,
    )
    if library is not None:
        alias = _resolve_library_id(manager, library, option="--library")
        if resolved is not None and alias is not None and int(resolved) != int(alias):
            raise ValueError("--library conflicts with --library-id/--library-name.")
        resolved = alias if alias is not None else resolved
    return resolved


def _print_error(message: str) -> int:
    """Write a standard CLI error message and return the failing exit code."""

    print(f"Error: {message}", file=sys.stderr)
    return 1


def _iter_payload_fields(model: type[BaseModel]) -> Iterable[Tuple[str, Mapping[str, Any]]]:
    """Yield task payload field metadata for the ``info`` command."""

    for name, field in model.model_fields.items():
        yield name, {
            "annotation": field.annotation,
            "required": field.is_required(),
            "default": None if field.is_required() else field.get_default(call_default_factory=True),
            "description": field.description,
        }


__all__ = [
    "LIST_COLOR_ENV_KEYS",
    "SELECTOR_DEFAULT_ENV_KEYS",
    "STORAGE_ROOT_KINDS",
    "STRICT_WORKING_ROOT_KINDS",
    "_apply_busco_context_from_args",
    "_coerce_bool",
    "_connect_manager",
    "_expand_accessions",
    "_format_selector_help",
    "_infer_db_path",
    "_iter_payload_fields",
    "_load_list_color_defaults",
    "_load_selector_defaults",
    "_normalize_action_alias",
    "_print_error",
    "_resolve_library_id",
    "_resolve_library_selector",
    "_resolve_library_selector_alias",
    "_resolve_task_spec",
]
