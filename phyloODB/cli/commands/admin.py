"""CLI registration and handlers for database administration commands."""
from __future__ import annotations

import argparse
import contextlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ...database import DBManager
from ...db.errors import StorageOperationError
from ...permissions import (
    PermissionPolicy,
    apply_shared_umask,
    comprehensive_directory_probe,
    resolve_scratch_dir,
    resolve_group,
    sqlite_wal_probe,
    validate_config_updates,
)
from ...selector_utils import resolve_selector_accessions
from ...services.discovery_service import DiscoveryService
from ...services.task_service import TaskService
from ...thread_defaults import computed_task_thread_defaults, detect_available_threads
from ...variable_definitions import environment_variable_definitions, validate_variables_json_document
from ...variable_kinds import SYSTEM_VARIABLE_EXACT, normalize_variable_kind
from ..support.argparse_utils import AppendCommaSeparated
from ..support.common import (
    _connect_manager,
    _iter_payload_fields,
    _print_error,
    _resolve_library_selector,
    _resolve_task_spec,
)
from ..support.output import _format_table, _render_list_output
from ..support.selectors import _add_basic_list_output_options, _add_selector_arguments, _selector_request_from_args

def _handle_set(args: argparse.Namespace) -> int:
    """Persist a database-backed environment variable from CLI input."""
    subject = getattr(args, "set_subject", "var")
    if subject in {"var", "env"}:
        json_path = getattr(args, "json_file", None)
        if json_path:
            if getattr(args, "assignment", None):
                return _print_error("Use either --json PATH or VARIABLE VALUE, not both.")
            if getattr(args, "kind", None):
                return _print_error("Do not combine --json with --kind; JSON sections define each variable kind.")
            try:
                raw = Path(json_path).read_text(encoding="utf-8")
                payload = json.loads(raw)
                values, kinds = validate_variables_json_document(payload)
            except OSError as exc:
                return _print_error(f"Could not read variable JSON '{json_path}': {exc}")
            except json.JSONDecodeError as exc:
                return _print_error(f"Invalid JSON in '{json_path}': {exc}")
            except ValueError as exc:
                return _print_error(str(exc))
            if not values:
                return _print_error("Variable JSON did not contain any importable variables.")
            manager = _connect_manager(args.database)
            try:
                values = validate_config_updates(manager, values)
                manager.set_environment_variables(values, kinds=kinds)
            except (ValueError, StorageOperationError) as exc:
                return _print_error(str(exc))
            finally:
                manager.close()
            print(f"Imported {len(values)} variable(s) from {json_path}.")
            return 0

        assignment = [str(part) for part in (getattr(args, "assignment", None) or []) if str(part).strip()]
        if not assignment:
            return _print_error("Provide VARIABLE VALUE or VARIABLE=VALUE.")
        if len(assignment) == 1 and "=" in assignment[0]:
            key, value = assignment[0].split("=", 1)
        elif len(assignment) >= 2:
            key = assignment[0]
            value = " ".join(assignment[1:])
        else:
            return _print_error("Invalid set syntax. Use VARIABLE VALUE or VARIABLE=VALUE.")
        key = str(key or "").strip().upper()
        value = str(value or "").strip()
        if not key:
            return _print_error("Provide a variable name.")
        if value == "":
            return _print_error("Provide a variable value.")
        try:
            value_json = json.loads(value)
        except json.JSONDecodeError:
            value_json = value
        kind = None
        raw_kind = getattr(args, "kind", None)
        if raw_kind:
            kind = normalize_variable_kind(raw_kind)
            if kind is None:
                return _print_error(f"Unknown variable kind '{raw_kind}'.")

        manager = _connect_manager(args.database)
        try:
            normalized = validate_config_updates(manager, {key: value_json})
            manager.set_environment_variable(key, normalized[key], kind=kind)
        except (ValueError, StorageOperationError) as exc:
            return _print_error(str(exc))
        finally:
            manager.close()
        if getattr(args, "legacy_syntax", False):
            print("Warning: legacy 'set VARIABLE=VALUE' syntax is deprecated. Use 'set var NAME VALUE'.")
        print(f"Set {key}.")
        return 0
    if subject == "busco-primary":
        return _handle_set_busco_primary(args)
    if subject == "proteome-profile":
        return _handle_set_proteome_profile(args)
    return _print_error(f"Unknown set subject '{subject}'.")


def _expand_run_id_tokens(manager: DBManager, raw_tokens: Sequence[Any]) -> list[int]:
    resolved: list[int] = []
    for token in raw_tokens or []:
        if token is None:
            continue
        text = str(token).strip()
        if not text:
            continue
        if text.startswith("@") and len(text) > 1:
            value = manager.get_environment_variable(text[1:])
            if value is None:
                raise ValueError(f"Unknown BUSCO run variable '{text[1:]}'.")
            if isinstance(value, str):
                items = [part.strip() for part in value.split(",") if part.strip()]
            elif isinstance(value, (list, tuple, set)):
                items = [str(part).strip() for part in value if str(part).strip()]
            else:
                raise ValueError(f"BUSCO run variable '{text[1:]}' must contain a run-id list.")
            for item in items:
                if not item.isdigit():
                    raise ValueError(f"BUSCO run variable '{text[1:]}' contains non-numeric run id '{item}'.")
                resolved.append(int(item))
            continue
        if not text.isdigit():
            raise ValueError(f"Invalid BUSCO run id '{text}'.")
        resolved.append(int(text))
    return list(dict.fromkeys(resolved))


def _run_score_summary(row: Optional[Sequence[Any]]) -> str:
    if not row:
        return "NA"
    try:
        if len(row) >= 20:
            input_mode = row[5]
            pipeline = row[6]
            sc = int(row[9] or 0)
        else:
            input_mode = row[4] if len(row) > 4 else None
            pipeline = row[5] if len(row) > 5 else None
            sc = int(row[8] or 0) if len(row) > 8 else 0
        pipe = str(pipeline or "").strip().lower()
        if pipe == "augustus":
            pipe_code = "Au"
        elif pipe == "miniprot":
            pipe_code = "Mi"
        elif pipe == "metaeuk":
            pipe_code = "Me"
        else:
            pipe_code = str(pipeline or "?")[:2] or "?"
        mode = str(input_mode or "").strip().lower()
        mode_code = "G" if mode == "genome" else "P"
        return f"{pipe_code},{mode_code},SC={sc}"
    except (TypeError, ValueError, IndexError):
        return f"run={row[0]}" if row and row[0] is not None else "NA"


def _format_busco_primary_note(statuses: Dict[str, str]) -> str:
    purpose_labels = {
        "default": "D",
        "export_protein": "EP",
        "export_nucleotide": "EN",
    }
    grouped: Dict[str, list[str]] = {}
    for purpose, status in statuses.items():
        grouped.setdefault(str(status), []).append(purpose_labels.get(purpose, str(purpose)))
    if not grouped:
        return "unchanged"
    if set(grouped.keys()) == {"unchanged"}:
        return "unchanged"
    order = ("refreshed", "manual_override", "cleared", "unsupported", "no_match", "unchanged")
    parts: list[str] = []
    for key in order:
        labels = grouped.get(key)
        if not labels:
            continue
        if key == "unchanged" and len(grouped) == 1:
            parts.append("unchanged")
        elif key == "unchanged":
            continue
        else:
            parts.append(f"{key}:{','.join(labels)}")
    return ";".join(parts) or "unchanged"


def _busco_primary_selector_scope_requested(args: argparse.Namespace) -> bool:
    fields = (
        "accessions",
        "root",
        "clade",
        "taxid",
        "all",
        "exclude_accessions",
        "exclude_clades",
        "exclude_taxids",
        "filter",
        "downloaded_only",
        "not_downloaded",
        "local_only",
        "not_local",
        "primary_only",
        "after",
        "before",
        "level",
        "has_busco_results",
        "missing_busco_results",
        "library_id",
        "library_name",
        "busco_library",
        "proteome_profile",
        "prefer_proteome_profile",
        "isoforms_cleaned",
        "raw_proteome",
        "busco_run_selection",
        "quantity",
        "rank",
        "ranks",
        "quantities",
        "sample_strategy",
        "sample_seed",
        "allow_duplicate_species",
        "use_busco",
    )
    for field in fields:
        value = getattr(args, field, None)
        if value in (None, False, "", [], (), set(), {}):
            continue
        return True
    return False


def _resolve_busco_primary_refresh_targets(
    manager: DBManager,
    args: argparse.Namespace,
    *,
    busco_library: Optional[int],
) -> list[tuple[str, int, str]]:
    scope_requested = _busco_primary_selector_scope_requested(args)
    selected: Optional[list[str]] = None
    if scope_requested:
        selectors = _selector_request_from_args(
            args,
            profile="assembly_with_exclusions",
            busco_library_id=busco_library,
            manager=manager,
        )
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=True,
            require_candidates=False,
            use_rule_selection=True,
        )
        if not selected:
            return []

    query = """
        SELECT DISTINCT r.accession, r.library_id, l.library_name
        FROM BUSCO_Runs r
        JOIN Libraries l ON l.library_id = r.library_id
        WHERE COALESCE(r.status, 'completed') = 'completed'
    """
    params: list[Any] = []
    if busco_library is not None:
        query += " AND r.library_id = ?"
        params.append(int(busco_library))
    if selected is not None:
        placeholders = ",".join("?" for _ in selected)
        query += f" AND r.accession IN ({placeholders})"
        params.extend([str(acc) for acc in selected])
    query += " ORDER BY r.accession, r.library_id"
    manager.cursor.execute(query, tuple(params))
    return [(str(accession), int(library_id), str(library_name or "")) for accession, library_id, library_name in (manager.cursor.fetchall() or [])]


def _handle_refresh_busco_primary(
    manager: DBManager,
    args: argparse.Namespace,
    *,
    busco_library: Optional[int],
) -> int:
    conflicting = []
    if getattr(args, "run_id", None) is not None:
        conflicting.append("--run-id")
    if getattr(args, "run_ids", None):
        conflicting.append("--run-ids")
    if getattr(args, "busco_pipeline", None):
        conflicting.append("--busco-pipeline")
    if getattr(args, "prefer_busco_pipeline", None):
        conflicting.append("--prefer-busco-pipeline")
    if getattr(args, "format", None) or getattr(args, "busco_input_mode", None):
        conflicting.append("--format")
    if getattr(args, "prefer_format", None) or getattr(args, "prefer_busco_input_mode", None):
        conflicting.append("--prefer-format")
    if getattr(args, "orthofinder_target_library", None):
        conflicting.append("--orthofinder-target-library")
    if conflicting:
        joined = ", ".join(conflicting)
        return _print_error(f"--refresh cannot be combined with manual run-pinning options: {joined}.")

    targets = _resolve_busco_primary_refresh_targets(manager, args, busco_library=busco_library)
    if not targets:
        print("No BUSCO runs matched the requested selectors.")
        return 0

    purposes = ("default", "export_protein", "export_nucleotide")
    output_rows: list[tuple[str, ...]] = []
    any_changes = False
    for accession, library_id, library_name in targets:
        current_runs = {purpose: manager.busco.get_primary_run(accession, library_id, purpose=purpose) for purpose in purposes}
        current_assignments = {purpose: manager.busco.get_primary_assignment(accession, library_id, purpose=purpose) for purpose in purposes}

        proposed_runs: dict[str, Optional[Sequence[Any]]] = {}
        note_statuses: Dict[str, str] = {}
        pair_changed = False

        for purpose in purposes:
            current_assignment = current_assignments[purpose]
            current_run_id = int(current_assignment[0]) if current_assignment and current_assignment[0] is not None else None
            if manager.busco.is_manual_primary_override(accession, library_id, purpose=purpose):
                proposed_runs[purpose] = current_runs[purpose]
                note_statuses[purpose] = "manual_override"
                continue

            best = manager.busco.choose_best_run(
                accession,
                library_id,
                purpose=purpose,
                preferred_proteome_profile=manager.proteomes.get_default_profile_name(accession),
            )
            if best is None:
                proposed_runs[purpose] = None
                if current_run_id is not None:
                    note_statuses[purpose] = "cleared"
                    pair_changed = True
                else:
                    note_statuses[purpose] = "unchanged"
                continue

            proposed_runs[purpose] = best
            proposed_run_id = int(best[0])
            if current_run_id != proposed_run_id:
                note_statuses[purpose] = "refreshed"
                pair_changed = True
            else:
                note_statuses[purpose] = "unchanged"

        if not bool(getattr(args, "dry", False)):
            manager.busco.refresh_auto_primary_runs_for_accession(
                accession,
                library_id,
                updated_by="set-busco-primary-refresh",
                policy="auto_best",
            )
        if pair_changed:
            any_changes = True
        output_rows.append(
            (
                accession,
                library_name,
                _run_score_summary(current_runs["default"]),
                _run_score_summary(proposed_runs["default"]),
                _run_score_summary(current_runs["export_protein"]),
                _run_score_summary(proposed_runs["export_protein"]),
                _run_score_summary(current_runs["export_nucleotide"]),
                _run_score_summary(proposed_runs["export_nucleotide"]),
                _format_busco_primary_note(note_statuses),
            )
        )

    headers = (
        "accession",
        "library_name",
        "current_default",
        "new_default",
        "current_export_protein",
        "new_export_protein",
        "current_export_nucleotide",
        "new_export_nucleotide",
        "note",
    )
    _render_list_output(args, headers, output_rows, default_tidy=True)
    if getattr(args, "dry", False):
        print("Dry run only. Re-run without --dry to apply BUSCO primary refresh.")
    elif not any_changes:
        print("No BUSCO primary changes were required.")
    return 0


def _matches_orthofinder_target_library(
    pipeline: Any,
    pipeline_params_json: Any,
    target_library: str,
) -> bool:
    """Return whether a BUSCO run came from the requested derived-library build."""

    if str(pipeline or "").strip().lower() != "orthofinder":
        return False
    try:
        pipeline_params = json.loads(pipeline_params_json) if pipeline_params_json else {}
    except (TypeError, ValueError, json.JSONDecodeError):
        pipeline_params = {}
    return str(pipeline_params.get("derived_library_name") or "").strip() == target_library


def _handle_set_busco_primary(args: argparse.Namespace) -> int:
    manager = _connect_manager(args.database, read_only=bool(getattr(args, "dry", False)))
    try:
        required_pipeline = str(getattr(args, "busco_pipeline", "") or "").strip().lower() or None
        orthofinder_target_library = str(
            getattr(args, "orthofinder_target_library", "") or ""
        ).strip() or None
        if orthofinder_target_library:
            if required_pipeline not in {None, "orthofinder"}:
                return _print_error(
                    "--orthofinder-target-library can only be combined with "
                    "--busco-pipeline orthofinder."
                )
            required_pipeline = "orthofinder"
        preferred_pipeline = str(getattr(args, "prefer_busco_pipeline", "") or "").strip().lower() or None
        required_mode = str(getattr(args, "format", None) or getattr(args, "busco_input_mode", None) or "").strip().lower() or None
        preferred_mode = str(getattr(args, "prefer_format", None) or getattr(args, "prefer_busco_input_mode", None) or "").strip().lower() or None
        required_profile = str(getattr(args, "proteome_profile", "") or "").strip() or None
        preferred_profile = str(getattr(args, "prefer_proteome_profile", "") or "").strip() or None
        if required_mode in {"nucl", "nucleotide"}:
            required_mode = "genome"
        if preferred_mode in {"nucl", "nucleotide"}:
            preferred_mode = "genome"

        explicit_run_ids = _expand_run_id_tokens(
            manager,
            ([getattr(args, "run_id", None)] if getattr(args, "run_id", None) is not None else [])
            + list(getattr(args, "run_ids", None) or []),
        )
        busco_library = _resolve_library_selector(
            manager,
            library_id=getattr(args, "library_id", None),
            library_name=getattr(args, "library_name", None),
            legacy=getattr(args, "busco_library", None),
        )
        if bool(getattr(args, "refresh", False)):
            return _handle_refresh_busco_primary(manager, args, busco_library=busco_library)
        if not explicit_run_ids and required_pipeline is None and required_mode is None:
            return _print_error(
                "Without --run-id/--run-ids, provide at least --format, "
                "--busco-pipeline, or --orthofinder-target-library."
            )
        selectors = _selector_request_from_args(
            args,
            profile="assembly_with_exclusions",
            busco_library_id=busco_library,
            manager=manager,
        )
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=True,
            require_candidates=False,
            use_rule_selection=True,
        )

        rows: list[tuple] = []
        query = """
            SELECT r.run_id, r.accession, r.library_id, l.library_name,
                   r.pipeline, r.pipeline_params_effective_json
            FROM BUSCO_Runs r
            JOIN Libraries l ON l.library_id = r.library_id
            WHERE COALESCE(r.status, 'completed') = 'completed'
        """
        params: list[Any] = []
        if explicit_run_ids:
            placeholders = ",".join("?" for _ in explicit_run_ids)
            query += f" AND r.run_id IN ({placeholders})"
            params.extend(explicit_run_ids)
        if busco_library is not None:
            query += " AND r.library_id = ?"
            params.append(int(busco_library))
        if selected:
            placeholders = ",".join("?" for _ in selected)
            query += f" AND r.accession IN ({placeholders})"
            params.extend([str(acc) for acc in selected])
        elif not explicit_run_ids:
            return _print_error("No accessions matched selectors.")
        query += " ORDER BY r.accession, r.library_id, r.run_id"
        manager.cursor.execute(query, tuple(params))
        rows = manager.cursor.fetchall() or []
        if orthofinder_target_library:
            rows = [
                row
                for row in rows
                if _matches_orthofinder_target_library(
                    row[4], row[5], orthofinder_target_library
                )
            ]
        if not rows:
            print("No BUSCO runs matched the requested selectors.")
            return 0

        candidate_run_ids = [int(row[0]) for row in rows]

        grouped: Dict[tuple[str, int], Dict[str, Any]] = {}
        for run_id, accession, library_id, library_name in rows:
            key = (str(accession), int(library_id))
            grouped.setdefault(key, {"library_name": str(library_name or "")})

        output_rows: list[tuple[str, ...]] = []
        any_changes = False
        for (accession, library_id), meta in sorted(grouped.items(), key=lambda item: (item[0][0], item[0][1])):
            chosen = manager.busco.choose_best_run(
                accession,
                library_id,
                purpose="default",
                run_ids=candidate_run_ids,
                pipeline=required_pipeline,
                input_mode=required_mode,
                preferred_pipeline=preferred_pipeline,
                preferred_input_mode=preferred_mode,
                proteome_profile=required_profile,
                preferred_proteome_profile=preferred_profile,
            )
            current_default = manager.busco.get_primary_run(accession, library_id, purpose="default")
            current_protein = manager.busco.get_primary_run(accession, library_id, purpose="export_protein")
            current_nucl = manager.busco.get_primary_run(accession, library_id, purpose="export_nucleotide")
            note_statuses: Dict[str, str] = {}
            if not chosen:
                note_statuses["default"] = "no_match"
                output_rows.append(
                    (
                        accession,
                        meta["library_name"],
                        _run_score_summary(current_default),
                        "NA",
                        _run_score_summary(current_protein),
                        "NA",
                        _run_score_summary(current_nucl),
                        "NA",
                        "no_match",
                    )
                )
                continue
            chosen_run_id = int(chosen[0])
            set_default = True
            set_protein = manager.busco.run_supports_purpose(chosen_run_id, purpose="export_protein")
            set_nucl = manager.busco.run_supports_purpose(chosen_run_id, purpose="export_nucleotide")
            if not set_protein:
                note_statuses["export_protein"] = "unsupported"
            if not set_nucl:
                note_statuses["export_nucleotide"] = "unsupported"
            if not bool(getattr(args, "dry", False)):
                manager.busco.set_primary_run(
                    accession=accession,
                    library_id=library_id,
                    run_id=chosen_run_id,
                    purpose="default",
                    policy="manual_override",
                    updated_by="set-busco-primary",
                )
                if set_protein:
                    manager.busco.set_primary_run(
                        accession=accession,
                        library_id=library_id,
                        run_id=chosen_run_id,
                        purpose="export_protein",
                        policy="manual_override",
                        updated_by="set-busco-primary",
                    )
                if set_nucl:
                    manager.busco.set_primary_run(
                        accession=accession,
                        library_id=library_id,
                        run_id=chosen_run_id,
                        purpose="export_nucleotide",
                        policy="manual_override",
                        updated_by="set-busco-primary",
                    )
            if current_default is None or int(current_default[0]) != chosen_run_id or (set_protein and (current_protein is None or int(current_protein[0]) != chosen_run_id)) or (set_nucl and (current_nucl is None or int(current_nucl[0]) != chosen_run_id)):
                any_changes = True
            output_rows.append(
                (
                    accession,
                    meta["library_name"],
                    _run_score_summary(current_default),
                    _run_score_summary(chosen if set_default else None),
                    _run_score_summary(current_protein),
                    _run_score_summary(chosen if set_protein else None),
                    _run_score_summary(current_nucl),
                    _run_score_summary(chosen if set_nucl else None),
                    _format_busco_primary_note(note_statuses) if note_statuses else "ok",
                )
            )

        headers = (
            "accession",
            "library_name",
            "current_default",
            "new_default",
            "current_export_protein",
            "new_export_protein",
            "current_export_nucleotide",
            "new_export_nucleotide",
            "note",
        )
        _render_list_output(args, headers, output_rows, default_tidy=True)
        if getattr(args, "dry", False):
            print("Dry run only. Re-run without --dry to apply BUSCO primary overrides.")
        elif not any_changes:
            print("No BUSCO primary changes were required.")
        return 0
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        manager.close()


def _handle_set_proteome_profile(args: argparse.Namespace) -> int:
    manager = _connect_manager(args.database, read_only=bool(getattr(args, "dry", False)))
    try:
        profile_name = str(getattr(args, "profile_name", "") or "").strip()
        if not profile_name:
            return _print_error("Provide --profile-name.")
        selectors = _selector_request_from_args(args, profile="assembly_with_exclusions", manager=manager)
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=True,
            require_candidates=False,
            use_rule_selection=True,
        )
        if not selected:
            return _print_error("No accessions matched selectors.")
        rows: list[tuple[str, str, str]] = []
        for accession in selected:
            profile = manager.proteomes.get_profile(str(accession), profile_name)
            if not profile:
                rows.append((str(accession), profile_name, "missing"))
                continue
            if not bool(getattr(args, "dry", False)):
                manager.proteomes.set_default_profile(str(accession), profile_id=int(profile[0]))
                library_rows = manager.cursor.execute(
                    "SELECT DISTINCT library_id FROM BUSCO_Runs WHERE accession = ?",
                    (str(accession),),
                ).fetchall() or []
                for (library_id,) in library_rows:
                    if library_id is None:
                        continue
                    manager.busco.refresh_auto_primary_runs_for_accession(
                        str(accession),
                        int(library_id),
                        updated_by="set-proteome-profile",
                        policy="auto_best",
                    )
            rows.append((str(accession), profile_name, "default"))
        headers = ("accession", "profile_name", "status")
        _render_list_output(args, headers, rows, default_tidy=True)
        if getattr(args, "dry", False):
            print("Dry run only. Re-run without --dry to apply proteome profile defaults.")
        return 0
    finally:
        manager.close()

def _handle_info(args: argparse.Namespace) -> int:
    """Show database stats, action help, or task specification details."""

    if args.subject is None:
        manager = _connect_manager(args.database)
        try:
            info = _database_overview(manager)
        finally:
            manager.close()
        print(info)
        return 0

    token = args.subject.lower()
    if token in {"list", "watch", "storage", "discover", "queue", "run", "set", "info", "clear", "reset", "create", "count", "assemblies", "purge"}:
        print(_action_help(token))
        return 0

    try:
        spec = _resolve_task_spec(token)
    except KeyError as exc:
        return _print_error(str(exc))

    print(_describe_task_spec(spec))
    return 0


# ---------------------------------------------------------------------------
# Database setup and maintenance
# Purpose: Group lifecycle helpers for database creation, reset, discovery,
# maintenance commands, and operator-facing diagnostics.
# ---------------------------------------------------------------------------

def init_default_environment(
    mgr: DBManager,
    db_path: str,
    working_dir: Optional[str] = None,
    *,
    cache_dir: Optional[str] = None,
    scratch_dir: Optional[str] = None,
    permission_policy: Optional[PermissionPolicy] = None,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> None:
    """Initialise default environment variables and ensure expected paths exist."""

    base_dir = working_dir or os.path.dirname(os.path.abspath(db_path)) or os.getcwd()
    genomes_dir = os.path.join(base_dir, "genomes")
    libraries_dir = os.path.join(base_dir, "libraries")
    orthofinder_dir = os.path.join(base_dir, "orthofinder")
    exports_dir = os.path.join(base_dir, "exports")
    reports_dir = os.path.join(base_dir, "reports")
    misc_dir = os.path.join(base_dir, "misc")
    cache_dir = cache_dir or os.path.join(base_dir, "cache")
    logs_dir = os.path.join(base_dir, "logs")
    for d in (genomes_dir, libraries_dir, orthofinder_dir, exports_dir, reports_dir, misc_dir, cache_dir, logs_dir):
        os.makedirs(d, exist_ok=True)

    # Defaults should be portable across fresh installs; tasks that require
    # external tools validate their configuration when invoked.
    detected_threads = detect_available_threads()
    env_defaults = {
        "GENOME_DIR": genomes_dir,
        "DAEMON_MAX_THREADS": detected_threads,
        "DAEMON_PROCESS_POLLING_TIME": 2,
        "BLOCKED_TASK_QUEUE_POLLING_TIME": 2,
        "SET_MAX_THREADS_ON_START": True,
        "LIBRARIES_DIR": libraries_dir,
        "EXPORTS_DIR": exports_dir,
        "REPORTS_DIR": reports_dir,
        "MISC_DIR": misc_dir,
        "CACHE_DIR": cache_dir,
        "SCRATCH_DIR": scratch_dir,
        "PROJECT_PERMISSION_MODE": (permission_policy or PermissionPolicy()).mode,
        "SHARED_GROUP": (permission_policy or PermissionPolicy()).group,
        "LOG_DIR": logs_dir,
        "BUSCO_BINARIES_PATH": "busco",
        "ORTHOFINDER_BINARIES_PATH": "orthofinder",
        "BUSCO_LINEAGE_DIR": os.path.join(libraries_dir, "lineages"),
        "ORTHOFINDER_OUTPUT_DIR": orthofinder_dir,
        "MAKEBLASTDB_PATH": "makeblastdb",
        "LOG_LEVEL": "DEBUG",
        "LOG_TO_CONSOLE": False,
        "LOG_FILE": os.path.join(logs_dir, "phyloodb.log"),
        "LOG_MAX_BYTES": 5 * 1024 * 1024,
        "LOG_BACKUPS": 5,
        "BLASTP_PATH": "blastp",
        "SELECTOR_SCORE_ORDER": ["busco", "refseq", "level", "n50", "date", "accession"],
        "SELECTOR_BUSCO_BUCKETS": [[0, -5],[10, -4],[20, -3],[30, -2],[40, -1],[50, 1],[60, 2],[70, 3],[80, 4],[90, 5],],
        "SELECTOR_DEFAULT_DOWNLOADED_ONLY": False,
        "SELECTOR_DEFAULT_PRIMARY_ONLY": False,
        "SELECTOR_DEFAULT_USE_BUSCO": False,
        "SELECTOR_DEFAULT_STATUS_MIN": None,
        "SELECTOR_DEFAULT_PROTEIN_ONLY": False,
        "DEFAULT_PROTEOME_CLEAN_ISOFORMS": True,
        "DEFAULT_PROTEOME_USE_GFF": True,
        "DEFAULT_PROTEOME_USE_CDHIT": False,
        "DEFAULT_PROTEOME_GFF_PRIORITY": False,
        "DEFAULT_PROTEOME_CDHIT_IDENTITY": 0.96,
        "DEFAULT_PROTEOME_MAX_CONCURRENT": 1,
        "DEFAULT_PROTEOME_THREADS_PER_JOB": 1,
        "DEFAULT_PROTEOME_SET_DEFAULT": True,
        "DEFAULT_PROTEOME_INPUT_PROFILE": "raw",
        "DEFAULT_BUSCO_FORMAT": "protein",
        "DEFAULT_BUSCO_PIPELINE": "miniprot",
        "BUSCO_AUGUSTUS_EVALUE": None,
        "BUSCO_AUGUSTUS_LIMIT": 3,
        "BUSCO_AUGUSTUS_LONG": False,
        "BUSCO_AUGUSTUS_SPECIES": None,
        "BUSCO_AUGUSTUS_PARAMETERS": None,
        "BUSCO_METAEUK_PARAMETERS": None,
        "BUSCO_METAEUK_RERUN_PARAMETERS": None,
        "BUSCO_MINIPROT_PARAMETERS": None,
        "BUSCO_MINIPROT_KEEP_REF_FILE": False,
    }
    definition_defaults = environment_variable_definitions()
    for name in sorted(SYSTEM_VARIABLE_EXACT):
        if name not in env_defaults:
            env_defaults[name] = definition_defaults.get(name, {}).get("default")
    if email:
        env_defaults["EMAIL"] = str(email).strip()
    if api_key:
        env_defaults["NCBI_API_KEY"] = str(api_key).strip()
    env_defaults.update(computed_task_thread_defaults(detected_threads))
    mgr.set_environment_variables(env_defaults, kind="env")


def _handle_clear(args: argparse.Namespace) -> int:
    """Clear queued tasks, optionally including running task records."""

    manager = _connect_manager(args.database)
    try:
        if getattr(args, "full", False):
            manager.reset_tasks(full=True)
        else:
            manager.clear_tasks(keep_running=True)
    finally:
        manager.close()
    if getattr(args, "full", False):
        print("Cleared all tasks.")
    else:
        print("Cleared queued tasks (running tasks preserved).")
    return 0


def _handle_create(args: argparse.Namespace) -> int:
    """Create a new database, seed defaults, and import taxonomy immediately."""

    db_path = os.path.abspath(args.database)
    if os.path.exists(db_path) and not args.force:
        return _print_error(f"Database already exists at {db_path}. Use --force to overwrite.")
    shared = bool(getattr(args, "shared", False))
    group_name = getattr(args, "group", None)
    if group_name and not shared:
        return _print_error("--group requires --shared.")
    policy = PermissionPolicy("shared" if shared else "private", group_name if shared else None)
    try:
        if shared:
            resolve_group(group_name)
        apply_shared_umask(policy)
        db_dir = os.path.dirname(db_path) or os.getcwd()
        base_dir = os.path.abspath(getattr(args, "working_dir", None) or db_dir)
        cache_dir = os.path.abspath(getattr(args, "cache_dir", None) or os.path.join(base_dir, "cache"))
        scratch_value = getattr(args, "scratch_dir", None)
        scratch_dir = resolve_scratch_dir(scratch_value)
        durable_dirs = [
            os.path.join(base_dir, name)
            for name in ("genomes", "libraries", "orthofinder", "exports", "reports", "misc", "logs")
        ] + [cache_dir]
        preexisting = {
            path: Path(path).exists()
            for path in [db_dir, *durable_dirs, scratch_dir]
        }
        checks = [sqlite_wal_probe(db_dir, policy=policy)]
        checks.extend(comprehensive_directory_probe(path, policy=policy, create=True) for path in durable_dirs)
        # Scratch is job-local and deliberately does not need the project's group/setgid policy.
        checks.append(comprehensive_directory_probe(scratch_dir, create=True))
        failures = [check for check in checks if not check.ok]
        if failures:
            for path, existed in reversed(list(preexisting.items())):
                if not existed:
                    with contextlib.suppress(OSError):
                        Path(path).rmdir()
            details = "\n".join(f"- {check.path}: {check.message}" for check in failures)
            return _print_error(f"Creation preflight failed; the existing database was not changed:\n{details}")
    except (OSError, StorageOperationError, ValueError) as exc:
        return _print_error(str(exc))
    if os.path.exists(db_path) and args.force:
        try:
            os.remove(db_path)
            for ext in ("-wal", "-shm"):
                sidecar = f"{db_path}{ext}"
                if os.path.exists(sidecar):
                    os.remove(sidecar)
        except OSError as exc:
            return _print_error(f"Failed to remove existing database: {exc}")
    manager = DBManager(db_path)
    try:
        manager.connect()
        manager.setup_database()
        init_default_environment(
            manager,
            db_path,
            getattr(args, "working_dir", None),
            cache_dir=cache_dir,
            scratch_dir=scratch_value,
            permission_policy=policy,
            email=getattr(args, "email", None),
            api_key=getattr(args, "api_key", None),
        )

        payload: Dict[str, Any] = {}
        working_dir = getattr(args, "working_dir", None)
        if working_dir:
            payload["working_dir"] = Path(working_dir)
        taxdump_path = getattr(args, "taxdump", None)
        if taxdump_path:
            payload["path_to_taxdump"] = Path(taxdump_path)
        retain_flag = getattr(args, "retain_taxdump", None)
        if retain_flag is not None:
            payload["retain_taxdump"] = bool(retain_flag)

        service = TaskService(db_path, db_manager=manager)
        service.run_immediately("create-taxonomy", payload=payload or {})
        if policy.shared:
            group = resolve_group(policy.group)
            os.chown(db_path, -1, group.gr_gid)
            os.chmod(db_path, 0o660)
    finally:
        manager.close()
    print(f"Created database at {db_path}.")
    if not getattr(args, "email", None) or not getattr(args, "api_key", None):
        print(
            "Warning: EMAIL and/or NCBI_API_KEY were not set during create. "
            "Add them with 'phyloODB <db> set env EMAIL ...' and "
            "'phyloODB <db> set env NCBI_API_KEY ...' before NCBI add/download workflows."
        )
    return 0


def _handle_reset(args: argparse.Namespace) -> int:
    """Reset the database schema and rebuild the initial taxonomy state."""

    if not args.confirm and not args.force:
        print("Reset will erase the entire database. Re-run with --confirm to proceed.")
        return 1

    manager = DBManager(args.database)
    try:
        manager.connect()
        saved = manager.get_environment_variables(
            ["PROJECT_PERMISSION_MODE", "SHARED_GROUP", "CACHE_DIR", "SCRATCH_DIR"]
        ) or {}
        policy = PermissionPolicy(
            str(saved.get("PROJECT_PERMISSION_MODE") or "private"),
            saved.get("SHARED_GROUP"),
        )
        apply_shared_umask(policy)
        manager.reset()
        init_default_environment(
            manager,
            args.database,
            getattr(args, "working_dir", None),
            cache_dir=saved.get("CACHE_DIR"),
            scratch_dir=saved.get("SCRATCH_DIR"),
            permission_policy=policy,
        )

        payload: Dict[str, Any] = {}
        working_dir = getattr(args, "working_dir", None)
        if working_dir:
            payload["working_dir"] = Path(working_dir)
        taxdump_path = getattr(args, "taxdump", None)
        if taxdump_path:
            payload["path_to_taxdump"] = Path(taxdump_path)
        retain_flag = getattr(args, "retain_taxdump", None)
        if retain_flag is not None:
            payload["retain_taxdump"] = bool(retain_flag)

        service = TaskService(args.database, db_manager=manager)
        service.run_immediately("create-taxonomy", payload=payload or {})

    finally:
        try:
            manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace reset result.
            print(f"Warning: failed to close database after reset: {exc}", file=sys.stderr)
    return 0


def _handle_migrate(args: argparse.Namespace) -> int:
    """Explicitly upgrade a recognized legacy database to the current schema."""

    manager = DBManager(args.database)
    try:
        manager.connect()
        starting_version = manager.get_schema_version()
        print(f"Database schema version: {starting_version}")
        applied = manager.migrate_database()
        if not applied:
            print("Database schema is already current.")
        else:
            for step in applied:
                print(f"Applied migration {step}.")
            print(f"Database schema version: {manager.get_schema_version()}")
        return 0
    finally:
        manager.close()


def _handle_discover(args: argparse.Namespace) -> int:
    """Discover assemblies and BUSCO runs from registered genomes roots."""

    manager = _connect_manager(args.database)
    try:
        service = DiscoveryService(manager)
        report = service.discover(
            root=getattr(args, "root", None),
            path=args.path,
            overwrite=bool(args.overwrite),
            dry_run=bool(args.dry_run),
            attempt_knowledge_update=bool(args.attempt_knowledge_update),
        )
        print("\n".join(report.as_lines()))
        return 0
    finally:
        manager.close()


def _handle_kill(args: argparse.Namespace) -> int:
    """Kill a task tree and mark it as errored."""

    manager = _connect_manager(args.database)
    try:
        ok = manager.kill_task_and_descendants(args.task_id, args.reason)
    finally:
        try:
            manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace kill result.
            print(f"Warning: failed to close database after kill: {exc}", file=sys.stderr)
    if not ok:
        print(f"Failed to kill task {args.task_id}")
        return 1
    print(f"Killed task {args.task_id} and descendants.")
    return 0


def _handle_cancel(args: argparse.Namespace) -> int:
    """Cancel a pending or suspended task tree without marking it as errored."""

    manager = _connect_manager(args.database)
    try:
        ok = manager.cancel_task_and_descendants(args.task_id, args.reason)
    finally:
        try:
            manager.close()
        except Exception as exc:  # boundary: CLI cleanup failure should not replace cancel result.
            print(f"Warning: failed to close database after cancel: {exc}", file=sys.stderr)
    if not ok:
        print(f"Failed to cancel task {args.task_id} (ensure it is pending/suspended).")
        return 1
    print(f"Canceled task {args.task_id} and eligible subtasks.")
    return 0

def _database_overview(manager: DBManager) -> str:
    """Return a compact database overview table for the ``info`` command."""

    stats: List[Tuple[str, str]] = []
    manager.cursor.execute("SELECT COUNT(*) FROM Genome")
    stats.append(("Genomes", str(manager.cursor.fetchone()[0])))
    manager.cursor.execute("SELECT COUNT(*) FROM Genome WHERE status > 0")
    stats.append(("Downloaded genomes", str(manager.cursor.fetchone()[0])))
    manager.cursor.execute("SELECT COUNT(*) FROM Libraries")
    stats.append(("Libraries", str(manager.cursor.fetchone()[0])))
    manager.cursor.execute("SELECT COUNT(*) FROM Tasks")
    stats.append(("Tasks", str(manager.cursor.fetchone()[0])))
    manager.cursor.execute("SELECT MIN(queue_time) FROM Tasks")
    created = manager.cursor.fetchone()[0]
    if created:
        stats.insert(0, ("Queue initialised", str(created)))
    return _format_table(("Metric", "Value"), stats)


def _action_help(action: str) -> str:
    """Return one-line help text for a top-level CLI action."""

    descriptions = {
        "list": "List tasks, queue/errors, assemblies, libraries, roots, or variables.",
        "watch": "Watch the live queue or task errors.",
        "storage": "Manage storage roots and preview/apply genome or library moves.",
        "discover": "Discover assemblies and BUSCO runs from registered genomes roots.",
        "count": "Count items matching selector filters (e.g., assemblies).",
        "assemblies": "Print assemblies matching selector filters and ranking rules.",
        "queue": "Queue a task for background execution via the daemon.",
        "run": "Execute a task immediately in the foreground.",
        "set": "Update an environment variable stored in the database.",
        "info": "Display metadata about a task, action, or the database.",
        "clear": "Remove queued tasks while preserving running tasks.",
        "create": "Create a new database file and initialize defaults.",
        "reset": "Drop and recreate the database schema (destructive).",
        "kill": "Kill a task and its subtasks (mark errored).",
        "purge": "Purge selected data (dry-run by default).",
    }
    return f"{action}: {descriptions.get(action, 'No description available.')}"


def _describe_task_spec(spec: TaskSpec) -> str:
    """Render a detailed task-spec description for the ``info`` command."""

    lines = [
        f"Key:          {spec.key}",
        f"Display name: {spec.display_name or spec.key}",
        f"Job type:     {spec.job_type}",
        f"Description:  {spec.description}",
    ]
    if spec.aliases:
        lines.append(f"Aliases:      {', '.join(spec.aliases)}")
    lines.append(f"Threads:      {spec.daemon.required_threads}")
    lines.append("Payload fields:")
    for name, info in _iter_payload_fields(spec.payload_model):
        annotation = info["annotation"]
        readable = getattr(annotation, "__name__", repr(annotation))
        requirement = "required" if info["required"] else f"default={info['default']!r}"
        description = info["description"] or ""
        lines.append(f"  - {name}: {readable} ({requirement}) {description}")
    return "\n".join(lines)


def _load_json_payload(text: Optional[str], file_path: Optional[str]) -> Dict[str, Any]:
    """Load a JSON object from inline text or a file path."""

    if text and file_path:
        raise ValueError("Use either --json or --payload-file, not both.")
    if file_path:
        raw = Path(file_path).read_text(encoding="utf-8")
    elif text:
        raw = text
    else:
        return {}
    payload = json.loads(raw)
    if not isinstance(payload, dict):
        raise ValueError("Payload JSON must decode to an object.")
    return payload



def register_admin_parsers(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    """Register top-level setup and maintenance commands."""

    discover_parser = subparsers.add_parser("discover", help="Discover assemblies and BUSCO runs from registered genomes roots.")
    discover_scope = discover_parser.add_mutually_exclusive_group()
    discover_scope.add_argument("--root", help="Registered genomes root id or exact label to scan.")
    discover_scope.add_argument(
        "--path",
        help="Path to scan within a registered genomes root. Errors if the path is not inside a registered genomes root.",
    )
    overwrite_group = discover_parser.add_mutually_exclusive_group()
    overwrite_group.add_argument(
        "--overwrite",
        action="store_true",
        help="Rebind known accessions to the discovered folder and rebuild BUSCO state for discovered runs.",
    )
    overwrite_group.add_argument(
        "--no-overwrite",
        dest="overwrite",
        action="store_false",
        help="Do not modify known accessions found at a different path.",
    )
    discover_parser.set_defaults(overwrite=False)
    discover_parser.add_argument("--dry-run", action="store_true", help="Show what would change without writing to the database.")
    discover_parser.add_argument(
        "--attempt-knowledge-update",
        action="store_true",
        help="For unknown NCBI-style accessions, fetch assembly knowledge before ingesting.",
    )
    discover_parser.set_defaults(handler=_handle_discover)

    set_parser = subparsers.add_parser("set", help="Set database variables or BUSCO primary selections.")
    set_subparsers = set_parser.add_subparsers(dest="set_subject", required=True)

    set_var = set_subparsers.add_parser("var", aliases=["env"], help="Set a database-backed variable.")
    set_var.add_argument("--legacy-syntax", action="store_true", help=argparse.SUPPRESS)
    set_var.add_argument("--json", dest="json_file", help="Import variables from a kinded JSON file.")
    set_var.add_argument(
        "--kind",
        choices=[
            "assemblies",
            "assembly",
            "accessions",
            "accession",
            "busco-runs",
            "busco-run",
            "runs",
            "env",
            "environment",
            "environmental",
            "enviornmental",
        ],
        help="Store an explicit variable kind instead of inferring from the value.",
    )
    set_var.add_argument("assignment", nargs="*", help="VARIABLE VALUE or VARIABLE=VALUE.")
    set_var.set_defaults(handler=_handle_set)

    set_busco_primary = set_subparsers.add_parser("busco-primary", help="Set BUSCO primary run overrides or refresh automatic BUSCO primaries.")
    set_busco_primary_selectors = set_busco_primary.add_argument_group("Selector options")
    set_busco_primary_output = set_busco_primary.add_argument_group("Output options")
    _add_selector_arguments(
        set_busco_primary_selectors,
        profile="assembly_with_exclusions",
        context_label="BUSCO primary selection",
        skip_fields={"busco_run_ids"},
    )
    set_busco_primary_selectors.add_argument("--run-id", type=int, help="Explicit BUSCO run id to set.")
    set_busco_primary_selectors.add_argument("--run-ids", action=AppendCommaSeparated, nargs="+", help="Explicit BUSCO run ids to consider or apply.")
    set_busco_primary_selectors.add_argument(
        "--refresh",
        action="store_true",
        help="Recompute automatic BUSCO primary assignments for matched accessions/libraries instead of setting a manual override.",
    )
    set_busco_primary_selectors.add_argument(
        "--orthofinder-target-library",
        help=(
            "Select OrthoFinder-derived BUSCO runs created for this target library. "
            "This implies --busco-pipeline orthofinder."
        ),
    )
    set_busco_primary_output.add_argument("--dry", action="store_true", help="Preview primary changes without writing.")
    _add_basic_list_output_options(set_busco_primary_output)
    set_busco_primary.set_defaults(handler=_handle_set)

    set_proteome_profile = set_subparsers.add_parser("proteome-profile", help="Set the default proteome profile for matched accessions.")
    set_proteome_profile_selectors = set_proteome_profile.add_argument_group("Selector options")
    set_proteome_profile_output = set_proteome_profile.add_argument_group("Output options")
    _add_selector_arguments(
        set_proteome_profile_selectors,
        profile="assembly_with_exclusions",
        context_label="proteome profile selection",
        skip_fields={"has_busco_results", "missing_busco_results", "library_id", "library_name", "busco_pipeline", "prefer_busco_pipeline", "busco_input_mode", "prefer_busco_input_mode", "busco_export_format", "busco_run_ids", "busco_run_selection", "busco_complete_min", "busco_single_min", "include_paralog_filtering_in_score", "include_decontamination_in_score", "use_decontamination_run", "use_paralog_run", "allow_ambiguous_contaminants", "strict_decontamination", "rescue_duplicates", "paralog_filtered", "not_paralog_filtered", "min_hidden_paralogs", "max_hidden_paralogs", "decontaminated", "not_decontaminated", "contaminated", "decontamination_run", "ignore_contaminated_assemblies", "proteome_profile", "prefer_proteome_profile", "isoforms_cleaned", "raw_proteome"},
    )
    set_proteome_profile_selectors.add_argument("--profile-name", required=True, help="Proteome profile to mark as the default for each matched accession.")
    set_proteome_profile_output.add_argument("--dry", action="store_true", help="Preview changes without writing.")
    _add_basic_list_output_options(set_proteome_profile_output)
    set_proteome_profile.set_defaults(handler=_handle_set)

    info_parser = subparsers.add_parser("info", help="Describe tasks, actions, or the database.")
    info_parser.add_argument("subject", nargs="?", help="Task key, action name, or omit for stats.")
    info_parser.set_defaults(handler=_handle_info)

    create_parser = subparsers.add_parser("create", help="Create a new database and initialize defaults.")
    create_parser.add_argument("--force", action="store_true", help="Overwrite if the database already exists.")
    create_parser.add_argument("--taxdump", default=None, help="Optional path to new_taxdump.tar.gz.")
    create_parser.add_argument("--retain-taxdump", dest="retain_taxdump", action="store_true", help="Retain downloaded taxdump after import.")
    create_parser.add_argument("--working-dir", default=None, help="Working directory for taxonomy import/default paths.")
    create_parser.add_argument("--cache-dir", default=None, help="Durable cache directory (default: WORKING_DIR/cache).")
    create_parser.add_argument("--scratch-dir", default=None, help="Disposable job scratch directory (default: runtime TMPDIR).")
    create_parser.add_argument("--shared", action="store_true", help="Create a POSIX group-shared project.")
    create_parser.add_argument("--group", default=None, help="Unix group for --shared projects.")
    create_parser.add_argument("--email", default=None, help="Optional email to store as EMAIL for NCBI-backed workflows.")
    create_parser.add_argument("--api-key", dest="api_key", default=None, help="Optional NCBI API key to store as NCBI_API_KEY.")
    create_parser.set_defaults(handler=_handle_create)

    clear_parser = subparsers.add_parser("clear", help="Clear all queued tasks.")
    clear_parser.add_argument("--full", action="store_true", help="Remove all tasks, including running ones.")
    clear_parser.set_defaults(handler=_handle_clear)

    reset_parser = subparsers.add_parser("reset", help="Reset the database (destructive).")
    reset_parser.add_argument("--confirm", action="store_true", help="Acknowledge that reset is destructive.")
    reset_parser.add_argument("--force", action="store_true", help="Alias for --confirm.")
    reset_parser.add_argument("--taxdump", default=None, help="Optional path to new_taxdump.tar.gz.")
    reset_parser.add_argument("--retain-taxdump", dest="retain_taxdump", action="store_true", help="Retain downloaded taxdump after import.")
    reset_parser.add_argument("--working-dir", default=None, help="Working directory for taxonomy import/default paths.")
    reset_parser.set_defaults(handler=_handle_reset)

    migrate_parser = subparsers.add_parser("migrate", help="Upgrade an existing database schema explicitly.")
    migrate_parser.set_defaults(handler=_handle_migrate)

    kill_parser = subparsers.add_parser("kill", help="Kill a task and its subtasks (marks errored).")
    kill_parser.add_argument("task_id", type=int, help="Task ID to kill.")
    kill_parser.add_argument("--reason", default="Killed by user", help="Reason to record.")
    kill_parser.set_defaults(handler=_handle_kill)

    cancel_parser = subparsers.add_parser("cancel", help="Cancel a pending/suspended task and its subtasks.")
    cancel_parser.add_argument("task_id", type=int, help="Task ID to cancel.")
    cancel_parser.add_argument("--reason", default="Canceled by user", help="Reason to record.")
    cancel_parser.set_defaults(handler=_handle_cancel)


__all__ = [
    "_handle_cancel",
    "_handle_clear",
    "_handle_create",
    "_handle_discover",
    "_handle_info",
    "_handle_kill",
    "_handle_migrate",
    "_handle_reset",
    "_handle_set",
    "init_default_environment",
    "register_admin_parsers",
]
