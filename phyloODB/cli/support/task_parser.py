"""Shared task-parser construction helpers for queue, run, and workflows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Literal, Mapping, Optional, Tuple, Union, get_args, get_origin

from ...database import DBManager
from ...registry import TaskSpec
from ...selector_utils import active_selector_overrides, apply_selector_enrichment, merge_selector_preset
from .argparse_utils import AppendCommaSeparated, _with_default_help
from .common import SELECTOR_DEFAULT_ENV_KEYS, _format_selector_help
from .selectors import COMMON_SHORT_ALIASES, _add_selector_arguments

VERIFY_HIDDEN_FIELDS = {
    "allow_duplicate_species",
    "allow_ambiguous_contaminants",
    "strict_decontamination",
    "include_paralog_filtering_in_score",
    "include_decontamination_in_score",
    "use_paralog_run",
    "use_decontamination_run",
    "paralog_filtered",
    "not_paralog_filtered",
    "min_hidden_paralogs",
    "max_hidden_paralogs",
    "decontaminated",
    "not_decontaminated",
    "contaminated",
    "decontamination_run",
    "ignore_contaminated_assemblies",
    "busco_pipeline",
    "busco_input_mode",
    "busco_run_selection",
}

VERIFY_LIBRARY_ONLY_FIELDS = {
    "accessions",
    "taxid",
    "clade",
    "downloaded_only",
    "after",
    "before",
    "level",
    "primary_only",
    "filters",
    "ranks",
    "quantities",
}


def _payload_option_strings(name: str, long_option: Optional[str] = None) -> Tuple[str, ...]:
    """Return shared short aliases for common payload-derived CLI options."""

    option = long_option or f"--{name.replace('_', '-')}"
    return (*COMMON_SHORT_ALIASES.get(name, ()), option)

def _add_payload_arguments(
    parser: argparse.ArgumentParser,
    spec: TaskSpec,
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
    selector_group: Optional[argparse._ArgumentGroup] = None,
    task_group: Optional[argparse._ArgumentGroup] = None,
    runtime_default_overrides: Optional[Mapping[str, Any]] = None,
    hidden_fields: Optional[set[str]] = None,
) -> Tuple[List[str], List[str]]:
    """Populate a parser from a task payload model and return tracked fields."""

    payload_fields: List[str] = []
    boolean_fields: List[str] = []
    field_names = set(spec.payload_model.model_fields.keys())
    selector_defaults = selector_defaults or {}
    runtime_default_overrides = runtime_default_overrides or {}
    hidden_fields = hidden_fields or set()
    selector_fields = {
        "accessions",
        "root",
        "exclude_accessions",
        "exclude_clades",
        "exclude_taxids",
        "clade",
        "taxid",
        "downloaded_only",
        "after",
        "before",
        "level",
        "primary_only",
        "quantity",
        "rank",
        "allow_duplicate_species",
        "library_id",
        "library_name",
        "busco_library_id",
        "busco_pipeline",
        "prefer_busco_pipeline",
        "busco_input_mode",
        "prefer_busco_input_mode",
        "proteome_profile",
        "prefer_proteome_profile",
        "isoforms_cleaned",
        "raw_proteome",
        "busco_run_selection",
        "ref_accessions",
        "busco_complete_min",
        "busco_single_min",
        "has_busco_results",
        "use_busco",
        "protein_only",
        "status_min",
        "include_paralog_filtering_in_score",
        "include_decontamination_in_score",
        "allow_ambiguous_contaminants",
        "strict_decontamination",
        "rescue_duplicates",
        "paralog_run_id",
        "decontamination_run_id",
        "use_paralog_run",
        "decontamination_run",
        "use_decontamination_run",
        "paralog_filtered",
        "not_paralog_filtered",
        "min_hidden_paralogs",
        "max_hidden_paralogs",
        "decontaminated",
        "not_decontaminated",
        "contaminated",
        "ignore_contaminated_assemblies",
    }
    for name, field in spec.payload_model.model_fields.items():
        if name in hidden_fields:
            continue
        option = f"--{name.replace('_', '-')}"
        annotation = field.annotation
        origin = get_origin(annotation)
        args = get_args(annotation)
        required = field.is_required()
        help_text = field.description or ""
        if name in SELECTOR_DEFAULT_ENV_KEYS:
            fallback = False if name in {"downloaded_only", "primary_only", "use_busco", "protein_only"} else None
            help_text = _format_selector_help(help_text, selector_defaults, name, fallback=fallback)
        field_default = None
        if not required:
            field_default = field.get_default(call_default_factory=True)
            if name in runtime_default_overrides:
                field_default = runtime_default_overrides[name]
            help_text = _with_default_help(help_text, field_default)

        kwargs: Dict[str, Any] = {"dest": name}
        if not required:
            kwargs["default"] = argparse.SUPPRESS
        target = parser
        if selector_group is not None and name in selector_fields:
            target = selector_group
        elif task_group is not None:
            target = task_group

        if origin is Union and type(None) in args:
            sub_args = tuple(arg for arg in args if arg is not type(None))  # noqa: E721
            if len(sub_args) == 1:
                annotation = sub_args[0]
                origin = get_origin(annotation)
                args = get_args(annotation)
            else:
                annotation = Union[sub_args]  # type: ignore

        if origin in {list, List}:
            element_type = args[0] if args else str
            kwargs.update(
                {
                    "action": AppendCommaSeparated,
                    "nargs": "+",
                    "metavar": getattr(element_type, "__name__", str(element_type)),
                    "help": help_text,
                }
            )
            target.add_argument(*_payload_option_strings(name, option), **kwargs)
            if name == "accessions" and "accession" not in field_names:
                target.add_argument(
                    "--accession",
                    action=AppendCommaSeparated,
                    dest=name,
                    nargs="+",
                    metavar="ACC",
                    help=argparse.SUPPRESS,
                )
            # Provide a singular alias for other accession lists (e.g. exclude_accessions)
            if name.endswith("accessions") and name != "accessions":
                alias = name[:-1] if name.endswith("s") else name  # basic singularisation
                target.add_argument(
                    f"--{alias.replace('_', '-')}",
                    action=AppendCommaSeparated,
                    dest=name,
                    nargs="+",
                    metavar="ACC",
                    help=argparse.SUPPRESS,
                )
            if name.endswith("clades"):
                alias = name[:-1] if name.endswith("s") else name
                target.add_argument(
                    f"--{alias.replace('_', '-')}",
                    action=AppendCommaSeparated,
                    dest=name,
                    nargs="+",
                    metavar="CLADE",
                    help=argparse.SUPPRESS,
                )
            if name.endswith("taxids"):
                alias = name[:-1] if name.endswith("s") else name
                target.add_argument(
                    f"--{alias.replace('_', '-')}",
                    action=AppendCommaSeparated,
                    dest=name,
                    nargs="+",
                    metavar="TAXID",
                    help=argparse.SUPPRESS,
                )
            if name == "busco_run_ids":
                target.add_argument(
                    "--run-ids",
                    action=AppendCommaSeparated,
                    dest=name,
                    nargs="+",
                    metavar="RUN_ID",
                    help="Limit BUSCO-aware selection to specific BUSCO run ids (comma-separated or @VARIABLE).",
                )
            payload_fields.append(name)
            continue

        if origin is Literal:
            kwargs["choices"] = list(args)

        base_type = annotation if origin is None else origin

        if base_type in (int, float, str, Path):
            kwargs["type"] = base_type
            kwargs["help"] = help_text
            option_strings = _payload_option_strings(name, option)
            if name == "accession" and (
                spec.key == "busco-run" or any(alias == "busco" for alias in (spec.aliases or ()))
            ):
                option_strings = ("-a", "--accession")
            target.add_argument(*option_strings, **kwargs)
            if name == "busco_pipeline":
                target.add_argument(
                    "--require-busco-pipeline",
                    type=base_type,
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=argparse.SUPPRESS,
                )
            if name == "busco_input_mode":
                target.add_argument(
                    "--format",
                    type=base_type,
                    choices=["protein", "genome", "nucleotide"],
                    default=argparse.SUPPRESS,
                    dest=name,
                    help="Require BUSCO runs by input format (protein or genome).",
                )
                target.add_argument(
                    "--require-format",
                    type=base_type,
                    choices=["protein", "genome", "nucleotide"],
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=argparse.SUPPRESS,
                )
            if name == "prefer_busco_input_mode":
                target.add_argument(
                    "--prefer-format",
                    type=base_type,
                    choices=["protein", "genome", "nucleotide"],
                    default=argparse.SUPPRESS,
                    dest=name,
                    help="Prefer BUSCO runs by input format while allowing fallback to other matching runs.",
                )
            payload_fields.append(name)
            continue

        if base_type is bool or (isinstance(base_type, type) and issubclass(base_type, bool)):
            if name == "force_redownload":
                target.add_argument(
                    "--redownload",
                    "--force",
                    action="store_true",
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=help_text,
                )
            elif name in {"skip_gff", "skip_cdhit"}:
                feature = name.removeprefix("skip_")
                toggle = target.add_mutually_exclusive_group()
                toggle.add_argument(
                    f"--{feature}",
                    action="store_false",
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=f"Enable {feature.upper()} proteome filtering.",
                )
                toggle.add_argument(
                    f"--skip-{feature}",
                    action="store_true",
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=f"Disable {feature.upper()} proteome filtering.",
                )
            elif name in {"include_paralog_filtering_in_score", "include_decontamination_in_score"}:
                target.add_argument(
                    *_payload_option_strings(name, option),
                    action="store_true",
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=help_text,
                )
                if name == "include_paralog_filtering_in_score":
                    target.add_argument(
                        "--ignore-paralog-filtering",
                        action="store_false",
                        dest=name,
                        default=argparse.SUPPRESS,
                        help=_with_default_help(
                            "Ignore paralog filtering when calculating BUSCO selector scores.",
                            field_default,
                        ),
                    )
                if name == "include_decontamination_in_score":
                    target.add_argument(
                        "--ignore-decontamination",
                        action="store_false",
                        dest=name,
                        default=argparse.SUPPRESS,
                        help=_with_default_help(
                            "Ignore decontamination when calculating BUSCO selector scores.",
                            field_default,
                        ),
                    )
            else:
                target.add_argument(
                    *_payload_option_strings(name, option),
                    action=argparse.BooleanOptionalAction,
                    default=argparse.SUPPRESS,
                    dest=name,
                    help=help_text,
                )
            if name == "use_busco":
                target.add_argument(
                    "--with-busco",
                    action="store_true",
                    dest=name,
                    default=argparse.SUPPRESS,
                    help=argparse.SUPPRESS,
                )
                target.add_argument(
                    "--without-busco",
                    action="store_false",
                    dest=name,
                    default=argparse.SUPPRESS,
                    help=argparse.SUPPRESS,
                )
            if name == "avoid_unclean_buscos":
                target.add_argument(
                    "--use-unclean-buscos",
                    action="store_false",
                    dest=name,
                    default=argparse.SUPPRESS,
                    help=_with_default_help(
                        "Use all BUSCOs when computing medians, even if previously flagged unclean.",
                        field_default,
                    ),
                )
            if name == "ignore_contaminated_assemblies":
                target.add_argument(
                    "--include-contaminated-assemblies",
                    action="store_false",
                    dest=name,
                    default=argparse.SUPPRESS,
                    help=_with_default_help(
                        "Include assemblies flagged as contaminated by decontamination.",
                        field_default,
                    ),
                )
            boolean_fields.append(name)
            payload_fields.append(name)
            continue

        kwargs["type"] = str
        kwargs["help"] = help_text
        target.add_argument(*_payload_option_strings(name, option), **kwargs)
        payload_fields.append(name)

    return payload_fields, boolean_fields


def _extract_payload(values: argparse.Namespace, fields: Iterable[str]) -> Dict[str, Any]:
    """Extract only task payload fields from a parsed argparse namespace."""

    namespace = vars(values)
    data: Dict[str, Any] = {}
    for field in fields:
        if field in namespace:
            value = namespace[field]
            if value is None:
                continue
            data[field] = value
    return data


def _build_task_parser(
    base_prog: str,
    db_path: str,
    action: str,
    spec: TaskSpec,
    *,
    include_threads: bool = False,
    selector_defaults: Optional[Mapping[str, Any]] = None,
) -> Tuple[argparse.ArgumentParser, List[str]]:
    """Build a task-specific parser for ``queue`` or ``run`` subcommands."""

    parser = argparse.ArgumentParser(
        prog=f"{base_prog} {db_path} {action} {spec.key}",
        description=spec.description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    selector_group = parser.add_argument_group("Selector options")
    task_group = parser.add_argument_group("Task options")
    hidden_fields: set[str] = set()
    hidden_fields.update({"selector_requested_accessions", "selector_skipped_accessions"})
    selector_profile = str((spec.metadata or {}).get("selector_profile", "")).strip().lower()
    hidden_fields.add("required_threads")
    if spec.key.startswith("verify"):
        hidden_fields.update(VERIFY_HIDDEN_FIELDS)
    if selector_profile in {"library_only", "orthofinder_only"}:
        hidden_fields.update(VERIFY_LIBRARY_ONLY_FIELDS)
    busco_selector_skip_fields: set[str] = set()
    if spec.key in {"busco-run", "BatchBuscoTask"} or any(alias in {"busco", "batch-busco"} for alias in (spec.aliases or ())):
        busco_selector_skip_fields.update(
            {
                "busco_pipeline",
                "prefer_busco_pipeline",
                "busco_input_mode",
                "prefer_busco_input_mode",
                "busco_export_format",
                "busco_run_ids",
                "busco_run_selection",
            }
        )
        hidden_fields.update(busco_selector_skip_fields)

    payload_fields, _ = _add_payload_arguments(
        parser,
        spec,
        selector_defaults=selector_defaults,
        selector_group=selector_group,
        task_group=task_group,
        runtime_default_overrides={"required_threads": spec.daemon.required_threads},
        hidden_fields=hidden_fields,
    )
    if spec.key == "add-library":
        task_group.add_argument(
            "--fast-tree",
            action="store_const",
            const="fasttree",
            dest="gene_tree_source",
            default=argparse.SUPPRESS,
            help="Shortcut for --gene-tree-source fasttree.",
        )
    if "accessions" in spec.payload_model.model_fields and selector_profile not in {"library_only", "orthofinder_only"}:
        _add_selector_arguments(
            selector_group,
            profile="task_dynamic",
            selector_defaults=selector_defaults,
            context_label="derived task selectors",
            skip_fields=set(payload_fields) | busco_selector_skip_fields,
            suppress_defaults=True,
            include_library_scope=False,
        )
    if include_threads:
        task_group.add_argument(
            *_payload_option_strings("threads", "--threads"),
            type=int,
            default=argparse.SUPPRESS,
            help=_with_default_help(
                "Override the required thread count for foreground execution.",
                spec.daemon.required_threads,
            ),
        )
    if not include_threads:
        task_group.add_argument(
            *_payload_option_strings("threads", "--threads"),
            type=int,
            dest="required_threads",
            default=argparse.SUPPRESS,
            help=_with_default_help(
                "Override the daemon required thread count for this task.",
                spec.daemon.required_threads,
            ),
        )
        if "required_threads" not in payload_fields:
            payload_fields.append("required_threads")
    return parser, payload_fields


# ---------------------------------------------------------------------------
# Selector enrichment and list rendering helpers
# Purpose: Prepare resolved selector payloads and format assembly-centric CLI
# output in plain text or rich-colour tables.
# ---------------------------------------------------------------------------

def _apply_selector_enrichment(
    manager: DBManager,
    spec: TaskSpec,
    payload: Dict[str, Any],
    parsed: argparse.Namespace,
) -> Dict[str, Any]:
    """Populate derived selector fields before a task payload is queued or run."""

    source = vars(parsed)
    preset_name = source.get("preset_name")
    if preset_name:
        source = merge_selector_preset(manager, preset_name, active_selector_overrides(source)).as_mapping()
    return apply_selector_enrichment(manager, spec, payload, source)


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


__all__ = [
    "VERIFY_HIDDEN_FIELDS",
    "VERIFY_LIBRARY_ONLY_FIELDS",
    "_add_payload_arguments",
    "_apply_selector_enrichment",
    "_build_task_parser",
    "_extract_payload",
    "_load_json_payload",
]
