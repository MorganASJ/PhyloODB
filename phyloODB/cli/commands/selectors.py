"""CLI commands for named selector presets."""
from __future__ import annotations

import argparse
import json
from typing import Any, Mapping, Sequence

from ...selector_utils import (
    SelectorRequest,
    prune_selector_mapping,
    resolve_selector_accessions,
    validate_filter_expressions,
)
from ..support.common import _connect_manager, _print_error
from ..support.selectors import _add_basic_list_output_options, _add_selector_arguments
from .listing import _render_list_output, _store_accession_variable


def _selector_summary(selector: Mapping[str, Any]) -> str:
    parts: list[str] = []
    for key in sorted(selector.keys()):
        value = selector[key]
        if isinstance(value, list):
            rendered = ",".join(str(item) for item in value)
        else:
            rendered = str(value)
        parts.append(f"{key}={rendered}")
    return " ".join(parts)


def _resolve_preset(manager, name: str) -> list[str]:
    preset = manager.selector_presets.get(name)
    if preset is None:
        raise ValueError(f"Unknown selector preset '{name}'.")
    selectors = SelectorRequest.from_mapping(preset.get("selector") or {})
    return resolve_selector_accessions(
        manager,
        selectors,
        allow_all=True,
        require_candidates=False,
        use_rule_selection=True,
    )


def _validate_selector_request(manager, request: SelectorRequest) -> None:
    if request.downloaded_only and request.not_downloaded:
        raise ValueError("Use only one of --downloaded-only or --not-downloaded.")
    if request.local_only and request.not_local:
        raise ValueError("Use only one of --local-only or --not-local.")
    if request.missing_busco_results and (
        request.busco_complete_min is not None or request.busco_single_min is not None
    ):
        raise ValueError("BUSCO threshold selectors cannot be combined with --missing-busco-results.")
    validate_filter_expressions(manager, request.filters)


def _handle_selector(args: argparse.Namespace) -> int:
    manager = _connect_manager(args.database, read_only=getattr(args, "selector_action", "") in {"list", "show", "preview"})
    try:
        action = getattr(args, "selector_action", None)
        if action == "save":
            request = SelectorRequest.from_namespace(args)
            _validate_selector_request(manager, request)
            selector = prune_selector_mapping(request.as_mapping())
            if not selector:
                return _print_error("selector save requires at least one selector flag.")
            name = manager.selector_presets.save(args.name, selector, description=getattr(args, "description", None))
            print(f"Saved selector preset '{name}'.")
            return 0

        if action == "list":
            presets = manager.selector_presets.list()
            if not presets:
                print("No selector presets saved.")
                return 0
            rows = [
                (
                    preset["preset_name"],
                    preset.get("description") or "",
                    preset.get("updated_at") or "",
                    _selector_summary(preset.get("selector") or {}),
                )
                for preset in presets
            ]
            return _render_list_output(args, ("Preset", "Description", "Updated", "Selector"), rows, default_tidy=True)

        if action == "show":
            preset = manager.selector_presets.get(args.name)
            if preset is None:
                return _print_error(f"Unknown selector preset '{args.name}'.")
            if getattr(args, "json", False):
                print(json.dumps(preset, sort_keys=True))
                return 0
            rows = [
                ("name", preset["preset_name"]),
                ("description", preset.get("description") or ""),
                ("created_at", preset.get("created_at") or ""),
                ("updated_at", preset.get("updated_at") or ""),
            ]
            for key, value in sorted((preset.get("selector") or {}).items()):
                rows.append((key, json.dumps(value) if isinstance(value, (list, dict)) else value))
            return _render_list_output(args, ("Field", "Value"), rows, default_tidy=True)

        if action in {"preview", "resolve"}:
            selected = _resolve_preset(manager, args.name)
            if action == "resolve" and getattr(args, "store_variable", None):
                _store_accession_variable(manager, args.store_variable, selected)
            for accession in selected:
                print(accession)
            return 0

        if action == "delete":
            if not manager.selector_presets.delete(args.name):
                return _print_error(f"Unknown selector preset '{args.name}'.")
            print(f"Deleted selector preset '{args.name}'.")
            return 0
    except ValueError as exc:
        return _print_error(str(exc))
    finally:
        manager.close()
    return _print_error("Unknown selector command.")


def register_selector_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> argparse.ArgumentParser:
    parser = subparsers.add_parser("selector", help="Save, inspect, and resolve named selector presets.")
    selector_subparsers = parser.add_subparsers(dest="selector_action", required=True)

    save = selector_subparsers.add_parser("save", help="Save a named selector preset.")
    save.add_argument("name", help="Preset name.")
    save.add_argument("--description", help="Optional preset description.")
    _add_selector_arguments(
        save.add_argument_group("Selector options"),
        profile="assembly_with_exclusions",
        context_label="selector preset",
        suppress_defaults=True,
        skip_fields={"preset"},
    )
    save.set_defaults(handler=_handle_selector)

    list_parser = selector_subparsers.add_parser("list", help="List selector presets.")
    _add_basic_list_output_options(list_parser.add_argument_group("Output options"))
    list_parser.set_defaults(handler=_handle_selector)

    show = selector_subparsers.add_parser("show", help="Show a selector preset.")
    show.add_argument("name", help="Preset name.")
    show.add_argument("--json", action="store_true", help="Emit preset data as JSON.")
    _add_basic_list_output_options(show.add_argument_group("Output options"))
    show.set_defaults(handler=_handle_selector)

    preview = selector_subparsers.add_parser("preview", help="Resolve a preset without storing it.")
    preview.add_argument("name", help="Preset name.")
    preview.set_defaults(handler=_handle_selector)

    resolve = selector_subparsers.add_parser("resolve", help="Resolve a preset and optionally store a panel.")
    resolve.add_argument("name", help="Preset name.")
    resolve.add_argument("-S", "--store", "--save-set", dest="store_variable", help="Store resolved accessions in a named set.")
    resolve.set_defaults(handler=_handle_selector)

    delete = selector_subparsers.add_parser("delete", help="Delete a selector preset.")
    delete.add_argument("name", help="Preset name.")
    delete.set_defaults(handler=_handle_selector)

    return parser


__all__ = ["register_selector_parser"]
