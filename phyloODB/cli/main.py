"""Main command-line entry point for PhyloODB."""
from __future__ import annotations

import argparse
import logging
import sys
from typing import Mapping, Optional, Sequence

from .commands.admin import init_default_environment, register_admin_parsers
from .commands.listing import (
    _handle_assemblies,
    _handle_count,
    _handle_list_assemblies,
    _handle_list_busco_runs,
    _handle_list_buscos,
    _handle_watch,
    register_assemblies_parser,
    register_count_parser,
    register_list_parser,
    register_watch_parser,
)
from .commands.purge import register_purge_parser
from .commands.status import register_status_parser
from .commands.storage import _handle_list_roots, register_storage_parser
from .commands.selectors import register_selector_parser
from .commands.task_exec import register_task_exec_parsers
from .commands.taxonomic_tree import register_tree_parser
from .support.common import _infer_db_path, _load_selector_defaults, _normalize_action_alias
from ..database import DBManager
from ..db.errors import PhyloODBDatabaseError
from ..errors import PhyloODBError


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Compatibility exports
# Purpose: Preserve the small surface still imported by tests and adjacent
# modules while the command implementation is split into dedicated files.
# ---------------------------------------------------------------------------

def _handle_list(args: argparse.Namespace) -> int:
    """Dispatch ``list`` subcommands using the split command modules."""

    if args.choice == "roots":
        return _handle_list_roots(args)
    if args.choice == "assemblies":
        return _handle_list_assemblies(args)
    if args.choice in {"busco", "results"}:
        args.busco = True
        args.has_busco_results = True
        return _handle_list_assemblies(args)
    if args.choice == "busco-runs":
        return _handle_list_busco_runs(args)
    if args.choice == "buscos":
        return _handle_list_buscos(args)

    from .commands import listing as listing_commands

    return listing_commands._handle_list(args)


def _dispatch_watch(args: argparse.Namespace) -> int:
    from .commands import listing as listing_commands

    return listing_commands._handle_watch(args)


# ---------------------------------------------------------------------------
# Parser construction
# Purpose: Build the root parser, then hand each top-level command module its
# own parser registration work.
# ---------------------------------------------------------------------------

def build_parser(selector_defaults: Optional[Mapping[str, object]] = None) -> argparse.ArgumentParser:
    """Build the root CLI parser and register top-level command modules."""

    parser = argparse.ArgumentParser(
        prog="phyloODB",
        description="CLI for working with the PhyloODB database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("database", help="Path to the PhyloODB SQLite database.")
    subparsers = parser.add_subparsers(dest="action", required=True)

    register_list_parser(subparsers, selector_defaults=selector_defaults, handler=_handle_list)
    register_watch_parser(subparsers, handler=_dispatch_watch)
    register_storage_parser(subparsers, selector_defaults=selector_defaults)
    register_tree_parser(subparsers, selector_defaults=selector_defaults)
    register_count_parser(subparsers, selector_defaults=selector_defaults)
    register_assemblies_parser(subparsers, selector_defaults=selector_defaults)
    register_task_exec_parsers(subparsers, selector_defaults=selector_defaults)
    register_selector_parser(subparsers)
    register_status_parser(subparsers)
    register_admin_parsers(subparsers)
    register_purge_parser(subparsers)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# Purpose: Normalise argv, pre-load selector defaults, and dispatch the chosen
# command handler through the parser built above.
# ---------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    """Parse CLI arguments and dispatch to the selected command handler."""

    raw_args = list(sys.argv[1:] if argv is None else argv)
    normalized_args = _normalize_action_alias(raw_args)
    db_path = _infer_db_path(normalized_args)
    try:
        selector_defaults = _load_selector_defaults(db_path)
    except (PhyloODBError, ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except Exception:  # boundary: CLI startup must return a stable exit status
        logger.exception("Unexpected failure while loading selector defaults")
        print("Error: unexpected failure while loading selector defaults; see the log for details.", file=sys.stderr)
        return 1
    parser = build_parser(selector_defaults=selector_defaults)
    try:
        args = parser.parse_args(normalized_args)
    except SystemExit as exc:
        return exc.code

    try:
        if "--help" not in normalized_args and "-h" not in normalized_args and args.action not in {"create", "reset", "migrate", "status"}:
            manager = DBManager(args.database, read_only=True)
            try:
                manager.connect()
                manager.validate_schema()
            finally:
                manager.close()
        return args.handler(args)
    except SystemExit:
        raise
    except (PhyloODBError, PhyloODBDatabaseError, ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except Exception:  # boundary: CLI must return a stable nonzero exit status
        logger.exception("Unexpected CLI failure")
        print("Error: unexpected internal failure; see the log for details.", file=sys.stderr)
        return 1


__all__ = [
    "_handle_assemblies",
    "_handle_count",
    "_handle_list",
    "_handle_list_assemblies",
    "_handle_list_busco_runs",
    "_handle_list_buscos",
    "_handle_list_roots",
    "_handle_watch",
    "build_parser",
    "init_default_environment",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
