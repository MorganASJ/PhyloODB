"""CLI registration and handlers for storage administration commands."""
from __future__ import annotations

import argparse
import sys
from typing import Any, Mapping, Optional, Sequence

from ...services.storage_admin_service import StorageAdminService
from ..support.argparse_utils import AppendCommaSeparated
from ..support.common import STORAGE_ROOT_KINDS, _print_error
from ..support.output import _format_tsv, _render_list_output
from ..support.selectors import _add_selector_arguments, _selector_request_from_args

def _handle_list_roots(args: argparse.Namespace) -> int:
    """List registered storage roots."""

    service = StorageAdminService(args.database)
    try:
        kind = getattr(args, "kind", None)
        if kind and kind not in STORAGE_ROOT_KINDS:
            return _print_error(f"--kind={kind} is invalid for list roots.")
        rows = service.list_roots(kind=kind)
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    if not rows:
        print("No storage roots registered.")
        return 0
    rendered = []
    for row in rows:
        rendered.append(
            (
                str(row[0]),
                str(row[1]),
                str(row[2] or ""),
                str(row[3] or ""),
                "yes" if bool(row[4]) else "no",
                "yes" if bool(row[5]) else "no",
            )
        )
    return _render_list_output(args, ("root_id", "kind", "label", "base_path", "writable", "active"), rendered, default_tidy=True)


def _print_plan_issues(issues: Sequence[str]) -> None:
    """Render preflight/storage plan issues."""

    if not issues:
        return
    print("Issues:")
    for issue in issues:
        print(f"- {issue}")


def _suspension_warning(kind: str) -> str:
    return (
        f"No active {kind} root remains. "
        f"Program operations that create new {kind} data are suspended until a root is activated."
    )


def _handle_storage_roots(args: argparse.Namespace) -> int:
    """Alias handler for ``storage roots``."""

    return _handle_list_roots(args)


def _handle_storage_add_root(args: argparse.Namespace) -> int:
    """Create or reuse a storage root."""

    service = StorageAdminService(args.database)
    try:
        result = service.add_root(
            kind=args.kind,
            base_path=args.base_path,
            label=args.label,
        )
        root_id = int(result["root_id"])
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    print(f"Storage root {root_id} ready.")
    if bool(result.get("created_inactive")):
        print(
            f"Created new inactive {args.kind} root. "
            f"It will not be used for new writes until you run 'storage activate-root {root_id}'."
        )
    return 0


def _handle_storage_rename_root(args: argparse.Namespace) -> int:
    """Rename a storage root label without changing its path or state."""

    service = StorageAdminService(args.database)
    try:
        result = service.rename_root(root_id=args.root_id, label=args.label)
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    old_label = result["old_label"] or "(unlabelled)"
    print(f"Renamed storage root {result['root_id']} from '{old_label}' to '{result['new_label']}'.")
    return 0


def _handle_storage_activate_root(args: argparse.Namespace) -> int:
    """Promote a root to be the active working root for its kind."""

    service = StorageAdminService(args.database)
    try:
        result = service.activate_root(root_id=args.root_id)
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    print(f"Activated storage root {result['root_id']} for kind '{result['kind']}'.")
    return 0


def _handle_storage_deactivate_root(args: argparse.Namespace) -> int:
    """Deactivate a root without promoting a replacement."""

    service = StorageAdminService(args.database)
    try:
        result = service.deactivate_root(root_id=args.root_id)
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    print(f"Deactivated storage root {result['root_id']} for kind '{result['kind']}'.")
    if bool(result.get("suspended")):
        print(_suspension_warning(str(result["kind"])))
    return 0


def _handle_storage_rebind_root(args: argparse.Namespace) -> int:
    """Preview or apply a whole-root rebind."""

    service = StorageAdminService(args.database)
    try:
        plan = service.plan_rebind_root(root_id=args.root_id, new_base_path=args.base_path)
        rows = [
            (
                str(plan.root_id),
                plan.kind,
                plan.label or "",
                plan.old_base_path,
                plan.new_base_path,
                str(plan.genome_count),
                str(plan.library_count),
                str(plan.artifact_count),
            )
        ]
        print(_format_tsv(("root_id", "kind", "label", "old_base", "new_base", "genomes", "libraries", "artifact_rows"), rows))
        if not getattr(args, "apply", False):
            print("Dry run only. Re-run with --apply to update the storage root.")
            return 0
        result = service.apply_rebind_root(
            root_id=args.root_id,
            new_base_path=args.base_path,
            verify=bool(getattr(args, "verify", True)),
        )
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    print(f"Rebound storage root {plan.root_id} to {plan.new_base_path}.")
    verify_task_ids = result.get("verify_task_ids") or []
    if verify_task_ids:
        print("Queued verify tasks: " + ", ".join(str(task_id) for task_id in verify_task_ids))
    return 0


def _handle_storage_flush_cache(args: argparse.Namespace) -> int:
    """Preview or apply a cache flush across cache roots."""

    service = StorageAdminService(args.database)
    try:
        plan = service.plan_flush_cache(root_id=getattr(args, "root_id", None))
        if not plan.root_ids:
            print("No cache roots matched.")
            return 0
        rows = [
            (
                ",".join(str(root_id) for root_id in plan.root_ids),
                ",".join(plan.root_paths),
                str(plan.artifact_count),
                str(plan.blastdb_count),
                str(plan.filesystem_entries),
            )
        ]
        print(_format_tsv(("root_ids", "base_paths", "artifact_rows", "blastdb_rows", "filesystem_entries"), rows))
        if not getattr(args, "apply", False):
            print("Dry run only. Re-run with --apply to flush cache.")
            return 0
        result = service.flush_cache(root_id=getattr(args, "root_id", None))
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    print(
        f"Flushed cache: deleted_paths={result['deleted_paths']} "
        f"deleted_artifacts={result['deleted_artifacts']} deleted_blastdbs={result['deleted_blastdbs']}"
    )
    return 0


def _handle_storage_move_genomes(args: argparse.Namespace) -> int:
    """Preview or apply genome storage moves."""

    service = StorageAdminService(args.database)
    try:
        service._ensure_connection()
        manager = service.db_manager
        selectors = _selector_request_from_args(args, profile="assembly_with_exclusions", manager=manager)
        plan = service.plan_move_genomes(
            request=selectors,
            target_root_id=args.to_root,
            rebind_only=bool(getattr(args, "rebind_only", False)),
        )
        rows = [
            (
                str(row["accession"]),
                str(row["source_path"] or ""),
                str(row["destination_path"] or ""),
                "" if row["source_root_id"] is None else str(row["source_root_id"]),
                str(row["destination_root_id"]),
                str(row["action"]),
            )
            for row in plan["rows"]
        ]
        if rows:
            print(_format_tsv(("accession", "source_path", "destination_path", "source_root_id", "destination_root_id", "action"), rows))
        else:
            print("No genomes matched the provided selectors.")
        _print_plan_issues(plan.get("issues") or [])
        if not getattr(args, "apply", False):
            print("Dry run only. Re-run with --apply to perform the move.")
            return 0
        if plan.get("issues"):
            return _print_error("Cannot apply move while preflight issues remain.")
        result = service.apply_move_genomes(
            plan,
            rebind_only=bool(getattr(args, "rebind_only", False)),
            verify=bool(getattr(args, "verify", True)),
        )
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    queued_task_id = result.get("queued_task_id")
    if queued_task_id is not None:
        print(f"Queued genome move task {queued_task_id}.")
        rebind_only = bool(getattr(args, "rebind_only", False))
        if bool(getattr(args, "verify", True)) and rebind_only:
            print(
                "The move task will rebind matched genomes to the destination root without copying files, "
                "run verify-assembly --repair --tidy first, then verify-busco --repair --reingest, and keep "
                "the source directory in place."
            )
        elif bool(getattr(args, "verify", True)):
            print(
                "The move task will copy data to the destination, run verify-assembly --repair --tidy first, "
                "then verify-busco --repair --reingest, suspend until they complete, and delete the original "
                "source only after successful verification."
            )
        elif rebind_only:
            print("The move task will rebind matched genomes to the destination root without copying files.")
        else:
            print("The move task will apply the move without verification.")
    return 0


def _handle_storage_move_libraries(args: argparse.Namespace) -> int:
    """Preview or apply library storage moves."""

    service = StorageAdminService(args.database)
    try:
        library_ids = service.resolve_library_ids(
            library_id=getattr(args, "library_id", None),
            library_name=getattr(args, "library_name", None),
            ref_accessions=getattr(args, "ref_accessions", None),
            all=bool(getattr(args, "all", False)),
        )
        plan = service.plan_move_libraries(
            library_ids=library_ids,
            target_root_id=args.to_root,
            rebind_only=bool(getattr(args, "rebind_only", False)),
        )
        rows = [
            (
                str(row["library_id"]),
                str(row["library_name"]),
                str(row["source_path"] or ""),
                str(row["destination_path"] or ""),
                "" if row["source_root_id"] is None else str(row["source_root_id"]),
                str(row["destination_root_id"]),
                str(row["action"]),
            )
            for row in plan["rows"]
        ]
        if rows:
            print(_format_tsv(("library_id", "library_name", "source_path", "destination_path", "source_root_id", "destination_root_id", "action"), rows))
        else:
            print("No libraries matched the provided selectors.")
        _print_plan_issues(plan.get("issues") or [])
        if not getattr(args, "apply", False):
            print("Dry run only. Re-run with --apply to perform the move.")
            return 0
        if plan.get("issues"):
            return _print_error("Cannot apply move while preflight issues remain.")
        result = service.apply_move_libraries(
            plan,
            rebind_only=bool(getattr(args, "rebind_only", False)),
            verify=bool(getattr(args, "verify", True)),
        )
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        try:
            service.db_manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace command result.
            print(f"Warning: failed to close storage service: {exc}", file=sys.stderr)
    moved = result.get("applied_library_ids") or []
    print(f"Moved {len(moved)} librar{'y' if len(moved) == 1 else 'ies'}.")
    verify_task_ids = result.get("verify_task_ids") or []
    if verify_task_ids:
        print("Queued verify tasks: " + ", ".join(str(task_id) for task_id in verify_task_ids))
    return 0


def _handle_storage_recover(args: argparse.Namespace) -> int:
    """Preview or retry pending journaled filesystem operations."""

    service = StorageAdminService(args.database)
    try:
        results = service.recover_filesystem_operations(
            operation_id=getattr(args, "operation_id", None),
            apply=bool(getattr(args, "apply", False)),
        )
    finally:
        service.db_manager.close()
    rows = [
        (
            str(row["operation_id"]),
            str(row["operation_type"]),
            str(row["status"]),
            str(row["action"]),
            "yes" if row["recoverable"] else "no",
            str(row["source_path"] or ""),
            str(row["staging_path"] or ""),
            str(row["destination_path"] or ""),
            str(row["error"] or ""),
        )
        for row in results
    ]
    if not rows:
        print("No pending filesystem operations.")
        return 0
    print(
        _format_tsv(
            (
                "operation_id",
                "type",
                "status",
                "action",
                "recoverable",
                "source",
                "staging",
                "destination",
                "error",
            ),
            rows,
        )
    )
    if not getattr(args, "apply", False):
        print("Dry run only. Re-run with --apply to retry safe finalization.")
    return 0


def _handle_storage(args: argparse.Namespace) -> int:
    """Dispatch storage admin subcommands."""

    if args.storage_action == "roots":
        return _handle_storage_roots(args)
    if args.storage_action == "add-root":
        return _handle_storage_add_root(args)
    if args.storage_action == "rename-root":
        return _handle_storage_rename_root(args)
    if args.storage_action == "activate-root":
        return _handle_storage_activate_root(args)
    if args.storage_action == "deactivate-root":
        return _handle_storage_deactivate_root(args)
    if args.storage_action == "rebind-root":
        return _handle_storage_rebind_root(args)
    if args.storage_action == "flush-cache":
        return _handle_storage_flush_cache(args)
    if args.storage_action == "move-genomes":
        return _handle_storage_move_genomes(args)
    if args.storage_action == "move-libraries":
        return _handle_storage_move_libraries(args)
    if args.storage_action == "recover":
        return _handle_storage_recover(args)
    return _print_error(f"Unknown storage action '{args.storage_action}'.")



def register_storage_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
) -> argparse.ArgumentParser:
    """Register the top-level ``storage`` command and its subcommands."""

    storage_parser = subparsers.add_parser("storage", help="Manage storage roots and move bound data.")
    storage_sub = storage_parser.add_subparsers(dest="storage_action", required=True)

    storage_roots = storage_sub.add_parser("roots", help="Alias of 'list roots'.")
    storage_roots.add_argument("--kind", choices=STORAGE_ROOT_KINDS, help="Filter listed roots by logical kind.")
    storage_roots.set_defaults(handler=_handle_storage)

    storage_add_root = storage_sub.add_parser("add-root", help="Create or reuse a storage root.")
    storage_add_root.add_argument("--kind", choices=STORAGE_ROOT_KINDS, required=True, help="Logical kind for the root.")
    storage_add_root.add_argument("--base-path", required=True, help="Base directory for the storage root.")
    storage_add_root.add_argument("--label", default=None, help="Optional human-readable label.")
    storage_add_root.set_defaults(handler=_handle_storage)

    storage_rename = storage_sub.add_parser("rename-root", help="Change a storage root's human-readable label.")
    storage_rename.add_argument("root_id", help="Storage root id or exact current label.")
    storage_rename.add_argument("--label", required=True, help="New unique, non-empty label.")
    storage_rename.set_defaults(handler=_handle_storage)

    storage_activate = storage_sub.add_parser("activate-root", help="Make a root the sole active write target for its kind.")
    storage_activate.add_argument("root_id", help="Storage root id or exact label to activate.")
    storage_activate.set_defaults(handler=_handle_storage)

    storage_deactivate = storage_sub.add_parser("deactivate-root", help="Deactivate a root. Strict kinds may be left without an active write target.")
    storage_deactivate.add_argument("root_id", help="Storage root id or exact label to deactivate.")
    storage_deactivate.set_defaults(handler=_handle_storage)

    storage_rebind = storage_sub.add_parser("rebind-root", help="Preview or apply a whole-root rebind.")
    storage_rebind.add_argument("root_id", help="Storage root id or exact label to rebind.")
    storage_rebind.add_argument("--base-path", required=True, help="New base path for the storage root.")
    storage_rebind.add_argument("--verify", dest="verify", action="store_true", default=True, help="Queue targeted verify tasks after apply.")
    storage_rebind.add_argument("--no-verify", dest="verify", action="store_false", help="Skip verify queueing after apply.")
    storage_rebind.add_argument("--apply", action="store_true", help="Apply the root rebind. Default is dry-run.")
    storage_rebind.set_defaults(handler=_handle_storage)

    storage_flush_cache = storage_sub.add_parser("flush-cache", help="Preview or clear files and DB rows under cache roots.")
    storage_flush_cache.add_argument("--root-id", help="Optional cache storage root id or exact label to flush.")
    storage_flush_cache.add_argument("--apply", action="store_true", help="Apply the cache flush. Default is dry-run.")
    storage_flush_cache.set_defaults(handler=_handle_storage)

    storage_move_genomes = storage_sub.add_parser("move-genomes", help="Preview or apply genome moves to another root.")
    storage_move_genomes_selector = storage_move_genomes.add_argument_group("Selector options")
    _add_selector_arguments(
        storage_move_genomes_selector,
        profile="assembly_with_exclusions",
        selector_defaults=selector_defaults,
        context_label="move genomes",
    )
    storage_move_genomes.add_argument("--to-root", required=True, help="Destination genomes storage root id or exact label.")
    storage_move_genomes.add_argument("--rebind-only", action="store_true", help="Only update bindings; assume files are already at the destination.")
    storage_move_genomes.add_argument("--verify", dest="verify", action="store_true", default=True, help="After apply, queue a move task that runs verify-assembly --repair --tidy first, then verify-busco --repair --reingest, before deleting the source.")
    storage_move_genomes.add_argument("--no-verify", dest="verify", action="store_false", help="Skip verify queueing after apply.")
    storage_move_genomes.add_argument("--apply", action="store_true", help="Apply the move. Default is dry-run.")
    storage_move_genomes.set_defaults(handler=_handle_storage)

    storage_move_libraries = storage_sub.add_parser("move-libraries", help="Preview or apply library moves to another root.")
    storage_move_libraries.add_argument("--library-id", type=int, help="Move a specific library id.")
    storage_move_libraries.add_argument("--library-name", help="Move a specific library by name.")
    storage_move_libraries.add_argument("--ref-accessions", action=AppendCommaSeparated, help="Resolve libraries by their reference accession set.")
    storage_move_libraries.add_argument("--all", action="store_true", help="Select all libraries.")
    storage_move_libraries.add_argument("--to-root", required=True, help="Destination libraries storage root id or exact label.")
    storage_move_libraries.add_argument("--rebind-only", action="store_true", help="Only update bindings; assume files are already at the destination.")
    storage_move_libraries.add_argument("--verify", dest="verify", action="store_true", default=True, help="Queue verify-libraries --repair after apply.")
    storage_move_libraries.add_argument("--no-verify", dest="verify", action="store_false", help="Skip verify queueing after apply.")
    storage_move_libraries.add_argument("--apply", action="store_true", help="Apply the move. Default is dry-run.")
    storage_move_libraries.set_defaults(handler=_handle_storage)

    storage_recover = storage_sub.add_parser("recover", help="Inspect or retry pending filesystem operations.")
    storage_recover.add_argument("--operation-id", type=int, help="Limit recovery to one operation id.")
    storage_recover.add_argument("--apply", action="store_true", help="Retry safe finalization. Default is dry-run.")
    storage_recover.set_defaults(handler=_handle_storage)
    return storage_parser


__all__ = [
    "_handle_list_roots",
    "_handle_storage",
    "_handle_storage_activate_root",
    "_handle_storage_add_root",
    "_handle_storage_deactivate_root",
    "_handle_storage_flush_cache",
    "_handle_storage_move_genomes",
    "_handle_storage_move_libraries",
    "_handle_storage_rename_root",
    "_handle_storage_rebind_root",
    "_handle_storage_recover",
    "register_storage_parser",
]
