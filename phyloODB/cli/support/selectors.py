"""Shared selector argument registration for CLI command modules."""
from __future__ import annotations

import argparse
from typing import Any, Dict, Mapping, Optional, Tuple

from ...selector_utils import SelectorRequest, merge_selector_preset, prune_selector_mapping
from .argparse_utils import AppendCommaSeparated, _validate_date, _with_default_help
from .common import _format_selector_help


# ---------------------------------------------------------------------------
# Selector profiles
# Purpose: Keep selector flag groupings in one place so command modules can
# opt into the exact selector surface they need.
# ---------------------------------------------------------------------------

SELECTOR_PROFILE_GROUPS: Dict[str, Tuple[str, ...]] = {
    "assembly": ("base", "busco", "rules"),
    "assembly_with_exclusions": ("base", "exclusions", "busco", "rules"),
    "busco_listing": ("base", "busco"),
    "task_dynamic": ("base", "busco", "rules"),
    "view_assemblies": ("base", "busco"),
}

COMMON_SHORT_ALIASES: Dict[str, Tuple[str, ...]] = {
    "accessions": ("-a",),
    "lineage": ("-l",),
    "format": ("-f",),
    "root": ("-rt",),
    "clade": ("-c",),
    "taxid": ("-i",),
    "downloaded_only": ("-d",),
    "after": ("-af",),
    "before": ("-bf",),
    "library_id": ("-li",),
    "library_name": ("-l",),
    "pipeline": ("-pl",),
    "busco_pipeline": ("-pl",),
    "quantity": ("-q",),
    "rank": ("-r",),
    "threads": ("-t",),
    "max_concurrent": ("-mc",),
}

LIST_OUTPUT_SHORT_ALIASES: Dict[str, Tuple[str, ...]] = {
    "tidy": ("-y",),
    "pretty": ("-p",),
    "busco": ("-b",),
    "meta": ("-m",),
    "store": ("-S",),
    "append": ("-A",),
}


def _selector_default(suppress_defaults: bool, value: Any) -> Any:
    """Return argparse defaults for shared selector helpers."""

    return argparse.SUPPRESS if suppress_defaults else value


def _option_strings(long_option: str, *aliases: str) -> Tuple[str, ...]:
    """Return argparse option strings with short aliases before the long flag."""

    return (*aliases, long_option)


# ---------------------------------------------------------------------------
# Selector parser helpers
# Purpose: Register selector flags consistently and convert parsed CLI values
# back into the shared selector request model.
# ---------------------------------------------------------------------------

def _add_selector_arguments(
    group: argparse._ArgumentGroup,
    *,
    profile: str,
    selector_defaults: Optional[Mapping[str, Any]] = None,
    context_label: str,
    skip_fields: Optional[set[str]] = None,
    suppress_defaults: bool = False,
    include_library_scope: bool = True,
) -> None:
    """Register a consistent selector option profile on an argparse group."""

    selector_defaults = selector_defaults or {}
    skip_fields = skip_fields or set()
    groups = SELECTOR_PROFILE_GROUPS.get(profile, ())

    def _skip(name: str) -> bool:
        return name in skip_fields

    if "base" in groups:
        if not _skip("preset"):
            group.add_argument(
                "--preset",
                dest="preset_name",
                default=_selector_default(suppress_defaults, None),
                help=f"Load a named selector preset for {context_label}.",
            )
        if not _skip("accessions"):
            group.add_argument(
                *_option_strings("--accessions", *COMMON_SHORT_ALIASES["accessions"]),
                action=AppendCommaSeparated,
                help=f"Accessions selector for {context_label}.",
            )
        if not _skip("root"):
            group.add_argument(
                *_option_strings("--root", *COMMON_SHORT_ALIASES["root"]),
                default=_selector_default(suppress_defaults, None),
                help=f"Limit {context_label} to entries bound under a specific storage root id or exact label.",
            )
        if not _skip("clade"):
            group.add_argument(
                *_option_strings("--clade", *COMMON_SHORT_ALIASES["clade"]),
                default=_selector_default(suppress_defaults, None),
                help=f"Resolve a scientific name to a taxid for {context_label}.",
            )
        if not _skip("taxid"):
            group.add_argument(
                *_option_strings("--taxid", *COMMON_SHORT_ALIASES["taxid"]),
                type=int,
                default=_selector_default(suppress_defaults, None),
                help=f"Restrict {context_label} to a taxid (and descendants).",
            )
        if profile == "assembly_with_exclusions" and not _skip("allow_all"):
            group.add_argument("--all", action="store_true", default=_selector_default(suppress_defaults, False), help=f"Preselect all assemblies before applying filters for {context_label}.")
        if "exclusions" in groups:
            if not _skip("exclude_accessions"):
                group.add_argument("--exclude-accessions", action=AppendCommaSeparated, help=f"Accessions to exclude from {context_label}.")
                group.add_argument("--exclude-accession", action=AppendCommaSeparated, dest="exclude_accessions", help=argparse.SUPPRESS)
            if not _skip("exclude_clades"):
                group.add_argument("--exclude-clades", action=AppendCommaSeparated, help=f"Clades to exclude from {context_label}.")
                group.add_argument("--exclude-clade", action=AppendCommaSeparated, dest="exclude_clades", help=argparse.SUPPRESS)
            if not _skip("exclude_taxids"):
                group.add_argument("--exclude-taxids", action=AppendCommaSeparated, help=f"Taxids to exclude from {context_label}.")
                group.add_argument("--exclude-taxid", action=AppendCommaSeparated, dest="exclude_taxids", help=argparse.SUPPRESS)
        if not _skip("filters"):
            group.add_argument(
                "--filter",
                action="append",
                default=_selector_default(suppress_defaults, None),
                help=(
                    "Filter by metadata/BUSCO values. Operators: =, !=, <, <=, >, >=, "
                    "~ or contains, !~ or not contains, in, not in, exists, missing. "
                    "Use ',' for AND and '|' for OR; repeat for additional filters."
                ),
            )
        if not _skip("downloaded_only"):
            group.add_argument(
                *_option_strings("--downloaded-only", *COMMON_SHORT_ALIASES["downloaded_only"]),
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=_format_selector_help(
                    f"Limit {context_label} to downloaded genomes.",
                    selector_defaults,
                    "downloaded_only",
                    fallback=False,
                ),
            )
        if not _skip("not_downloaded"):
            group.add_argument(
                "--not-downloaded",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=_with_default_help(
                    f"Limit {context_label} to genomes not yet downloaded.",
                    False,
                ),
            )
        if not _skip("local_only"):
            group.add_argument(
                "--local-only",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=_with_default_help(
                    f"Limit {context_label} to assemblies with origin 'local'.",
                    False,
                ),
            )
        if not _skip("not_local"):
            group.add_argument(
                "--not-local",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=_with_default_help(
                    f"Exclude assemblies with origin 'local' from {context_label}.",
                    False,
                ),
            )
        if not _skip("primary_only"):
            group.add_argument(
                "--primary-only",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=_format_selector_help(
                    f"Limit {context_label} to primary assemblies only.",
                    selector_defaults,
                    "primary_only",
                    fallback=False,
                ),
            )
        if not _skip("after"):
            group.add_argument(
                *_option_strings("--after", *COMMON_SHORT_ALIASES["after"]),
                type=lambda value: _validate_date(value, "--after"),
                default=_selector_default(suppress_defaults, None),
                help=f"Assemblies released on/after YYYY-MM-DD for {context_label}.",
            )
        if not _skip("before"):
            group.add_argument(
                *_option_strings("--before", *COMMON_SHORT_ALIASES["before"]),
                type=lambda value: _validate_date(value, "--before"),
                default=_selector_default(suppress_defaults, None),
                help=f"Assemblies released on/before YYYY-MM-DD for {context_label}.",
            )
        if not _skip("level"):
            group.add_argument(
                "--level",
                choices=["complete genome", "chromosome", "scaffold", "contig"],
                default=_selector_default(suppress_defaults, None),
                help=f"Assembly level filter for {context_label}.",
            )

    if "busco" in groups:
        if not _skip("has_busco_results"):
            group.add_argument(
                "--has-busco-results",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=f"Filter {context_label} to accessions with BUSCO results.",
            )
        if not _skip("missing_busco_results"):
            group.add_argument(
                "--missing-busco-results",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help=f"Filter {context_label} to accessions missing BUSCO results (use --library-name/--library-id to scope).",
            )
        if include_library_scope:
            if not _skip("library_id"):
                group.add_argument(
                    *_option_strings("--library-id", *COMMON_SHORT_ALIASES["library_id"]),
                    type=int,
                    default=_selector_default(suppress_defaults, None),
                    help="Limit BUSCO filters/output to a specific library id.",
                )
            if not _skip("library_name"):
                group.add_argument(
                    *_option_strings("--library-name", *COMMON_SHORT_ALIASES["library_name"]),
                    default=_selector_default(suppress_defaults, None),
                    help="Limit BUSCO filters/output to a specific library name.",
                )
            if not _skip("busco_library"):
                group.add_argument("--busco-library", type=str, default=_selector_default(suppress_defaults, None), help=argparse.SUPPRESS)
        if not _skip("busco_pipeline"):
            group.add_argument(
                *_option_strings("--busco-pipeline", *COMMON_SHORT_ALIASES["busco_pipeline"]),
                "--require-busco-pipeline",
                default=_selector_default(suppress_defaults, None),
                dest="busco_pipeline",
                help="Require BUSCO runs from a specific pipeline (e.g. miniprot, metaeuk, augustus).",
            )
        if not _skip("prefer_busco_pipeline"):
            group.add_argument(
                "--prefer-busco-pipeline",
                default=_selector_default(suppress_defaults, None),
                dest="prefer_busco_pipeline",
                help="Prefer BUSCO runs from a specific pipeline while allowing fallback to other matching runs.",
            )
        if not _skip("busco_input_mode"):
            group.add_argument(
                "--format",
                "--require-format",
                default=_selector_default(suppress_defaults, None),
                dest="format",
                choices=["protein", "genome", "nucleotide"],
                help="Require BUSCO runs by input format (protein or genome).",
            )
            group.add_argument(
                "--busco-input-mode",
                default=_selector_default(suppress_defaults, None),
                dest="busco_input_mode",
                help=argparse.SUPPRESS,
            )
        if not _skip("prefer_busco_input_mode"):
            group.add_argument(
                "--prefer-format",
                default=_selector_default(suppress_defaults, None),
                dest="prefer_format",
                choices=["protein", "genome", "nucleotide"],
                help="Prefer BUSCO runs by input format while allowing fallback to other matching runs.",
            )
        if not _skip("proteome_profile"):
            group.add_argument(
                "--proteome-profile",
                default=_selector_default(suppress_defaults, None),
                dest="proteome_profile",
                help="Require BUSCO/proteome-aware selection to use a specific proteome profile.",
            )
        if not _skip("prefer_proteome_profile"):
            group.add_argument(
                "--prefer-proteome-profile",
                default=_selector_default(suppress_defaults, None),
                dest="prefer_proteome_profile",
                help="Prefer BUSCO runs from a specific proteome profile while allowing fallback to other matching runs.",
            )
        if not _skip("isoforms_cleaned"):
            group.add_argument(
                "--isoforms-cleaned",
                action=argparse.BooleanOptionalAction,
                default=_selector_default(suppress_defaults, None),
                dest="isoforms_cleaned",
                help="Shortcut for --proteome-profile clean_default. Use --no-isoforms-cleaned to request raw proteomes.",
            )
        if not _skip("raw_proteome"):
            group.add_argument(
                "--raw-proteome",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                dest="raw_proteome",
                help="Shortcut for --proteome-profile raw.",
            )
        if not _skip("busco_export_format"):
            group.add_argument(
                "--export-format",
                default=_selector_default(suppress_defaults, None),
                dest="busco_export_format",
                choices=["protein", "nucleotide"],
                help="Require the selected BUSCO run to support protein or nucleotide export.",
            )
        if not _skip("busco_run_ids"):
            group.add_argument(
                "--run-ids",
                action=AppendCommaSeparated,
                default=_selector_default(suppress_defaults, None),
                dest="busco_run_ids",
                help="Limit BUSCO-aware selection to specific BUSCO run ids (comma-separated or @VARIABLE).",
            )
        if not _skip("busco_run_selection"):
            group.add_argument(
                "--busco-run-selection",
                default=_selector_default(suppress_defaults, None),
                help="BUSCO run-selection policy for downstream analyses (primary or latest).",
            )
        if not _skip("busco_complete_min"):
            group.add_argument("--busco-complete-min", type=float, default=_selector_default(suppress_defaults, None), help="Minimum BUSCO complete proportion (0-1 or 0-100).")
        if not _skip("busco_single_min"):
            group.add_argument("--busco-single-min", type=float, default=_selector_default(suppress_defaults, None), help="Minimum BUSCO single-copy proportion (0-1 or 0-100).")
        if not _skip("include_paralog_filtering_in_score"):
            group.add_argument(
                "--include-paralog-filtering-in-score",
                action="store_true",
                default=_selector_default(suppress_defaults, None),
                help="Include paralog filtering when calculating BUSCO scores (custom libraries only).",
            )
            group.add_argument(
                "--ignore-paralog-filtering",
                action="store_false",
                dest="include_paralog_filtering_in_score",
                default=_selector_default(suppress_defaults, None),
                help="Ignore paralog filtering when calculating BUSCO scores (custom libraries only).",
            )
        if not _skip("include_decontamination_in_score"):
            group.add_argument(
                "--include-decontamination-in-score",
                action="store_true",
                default=_selector_default(suppress_defaults, None),
                help="Include decontamination when calculating BUSCO scores.",
            )
            group.add_argument(
                "--ignore-decontamination",
                action="store_false",
                dest="include_decontamination_in_score",
                default=_selector_default(suppress_defaults, None),
                help="Ignore decontamination when calculating BUSCO scores (custom libraries only).",
            )
        if not _skip("use_decontamination_run"):
            group.add_argument("--use-decontamination-run", default=_selector_default(suppress_defaults, None), help="Use a specific decontamination run id for scoring.")
        if not _skip("use_paralog_run"):
            group.add_argument("--use-paralog-run", default=_selector_default(suppress_defaults, None), help="Use a specific paralog filtering run id for scoring.")
        if not _skip("allow_ambiguous_contaminants"):
            group.add_argument(
                "--allow-ambiguous-contaminants",
                action=argparse.BooleanOptionalAction,
                default=_selector_default(suppress_defaults, None),
                help="Treat decontamination 'unknown' BUSCOs as supported.",
            )
        if not _skip("strict_decontamination"):
            group.add_argument(
                "--strict-decontamination",
                action=argparse.BooleanOptionalAction,
                default=_selector_default(suppress_defaults, None),
                help="Treat decontamination 'weak' BUSCOs as contaminated.",
            )
        if not _skip("rescue_duplicates"):
            group.add_argument(
                "--rescue-duplicates",
                action="store_true",
                default=_selector_default(suppress_defaults, False),
                help="Treat duplicated BUSCO families with exactly one copy passing the active filters as effective single-copy.",
            )
        if not _skip("paralog_filtered") or not _skip("not_paralog_filtered"):
            paralog_group = group.add_mutually_exclusive_group()
            if not _skip("paralog_filtered"):
                paralog_group.add_argument("--paralog-filtered", action="store_true", default=_selector_default(suppress_defaults, False), help=f"Only include paralog-filtered accessions for {context_label}.")
            if not _skip("not_paralog_filtered"):
                paralog_group.add_argument("--not-paralog-filtered", action="store_true", default=_selector_default(suppress_defaults, False), help=f"Only include accessions without paralog filtering for {context_label}.")
        if not _skip("min_hidden_paralogs"):
            group.add_argument("--min-hidden-paralogs", type=float, default=_selector_default(suppress_defaults, None), help="Minimum hidden paralog proportion (0-1).")
        if not _skip("max_hidden_paralogs"):
            group.add_argument("--max-hidden-paralogs", type=float, default=_selector_default(suppress_defaults, None), help="Maximum hidden paralog proportion (0-1).")
        if not _skip("decontaminated") or not _skip("not_decontaminated") or not _skip("contaminated"):
            decont_group = group.add_mutually_exclusive_group()
            if not _skip("decontaminated"):
                decont_group.add_argument("--decontaminated", action="store_true", default=_selector_default(suppress_defaults, False), help=f"Only include decontaminated accessions for {context_label}.")
            if not _skip("not_decontaminated"):
                decont_group.add_argument("--not-decontaminated", action="store_true", default=_selector_default(suppress_defaults, False), help=f"Only include accessions without decontamination for {context_label}.")
            if not _skip("contaminated"):
                decont_group.add_argument("--contaminated", action="store_true", default=_selector_default(suppress_defaults, False), help=f"Only include accessions marked contaminated for {context_label}.")
        if not _skip("decontamination_run"):
            group.add_argument("--decontamination-run", default=_selector_default(suppress_defaults, None), help=f"Only include accessions in a specific decontamination run for {context_label}.")
        if not _skip("ignore_contaminated_assemblies"):
            contam_group = group.add_mutually_exclusive_group()
            if not suppress_defaults:
                contam_group.set_defaults(ignore_contaminated_assemblies=None)
            contam_group.add_argument(
                "--ignore-contaminated-assemblies",
                action="store_true",
                dest="ignore_contaminated_assemblies",
                default=_selector_default(suppress_defaults, None),
                help="Exclude assemblies flagged as contaminated by the latest decontamination run.",
            )
            contam_group.add_argument(
                "--include-contaminated-assemblies",
                action="store_false",
                dest="ignore_contaminated_assemblies",
                default=_selector_default(suppress_defaults, None),
                help="Include assemblies flagged as contaminated by decontamination.",
            )

    if "rules" in groups:
        if not _skip("quantity"):
            group.add_argument(
                *_option_strings("--quantity", *COMMON_SHORT_ALIASES["quantity"]),
                type=int,
                default=_selector_default(suppress_defaults, None),
                help="Top-N selector (overall, or per rank when --rank is set on explicit accessions or taxid/clade-derived candidates).",
            )
        if not _skip("rank"):
            group.add_argument(
                *_option_strings("--rank", *COMMON_SHORT_ALIASES["rank"]),
                default=_selector_default(suppress_defaults, None),
                help="Taxonomic rank for rule-based groups within the resolved candidate set (e.g., genus, family).",
            )
        if not _skip("ranks"):
            group.add_argument(
                *_option_strings("--ranks", *COMMON_SHORT_ALIASES.get("ranks", ())),
                action=AppendCommaSeparated,
                default=_selector_default(suppress_defaults, None),
                help="Multi-stage rank list for subsampling within the resolved candidate set (e.g., phylum,genus).",
            )
        if not _skip("quantities"):
            group.add_argument(
                *_option_strings("--quantities", *COMMON_SHORT_ALIASES.get("quantities", ())),
                action=AppendCommaSeparated,
                default=_selector_default(suppress_defaults, None),
                help="Quantities per rank in --ranks for the resolved candidate set (use 'all' or '*' to keep all).",
            )
        if not _skip("sample_strategy"):
            group.add_argument(
                "--sample-strategy",
                choices=["rank", "random"],
                default=_selector_default(suppress_defaults, "rank"),
                help="Sampling strategy for rank-based selection.",
            )
        if not _skip("sample_seed"):
            group.add_argument("--sample-seed", type=int, default=_selector_default(suppress_defaults, None), help="Random seed for sampling.")
        if not _skip("allow_duplicate_species"):
            group.add_argument("--allow-duplicate-species", action="store_true", default=_selector_default(suppress_defaults, False), help="Allow duplicate species in rule-based selection.")


def _selector_request_from_args(
    args: argparse.Namespace,
    *,
    profile: str,
    busco_library_id: Optional[int] = None,
    manager: Any = None,
) -> SelectorRequest:
    """Build the canonical selector request for a CLI handler/profile."""

    request = SelectorRequest.from_namespace(args)
    overrides: Dict[str, Any] = {}
    if "exclusions" not in SELECTOR_PROFILE_GROUPS.get(profile, ()):
        overrides.update(
            exclude_accessions=[],
            exclude_taxids=[],
            exclude_clades=[],
        )
    if "rules" not in SELECTOR_PROFILE_GROUPS.get(profile, ()):
        overrides.update(
            quantity=None,
            rank=None,
            ranks=[],
            quantities=[],
            sample_strategy=None,
            sample_seed=None,
            allow_duplicate_species=None,
        )
    if "busco" not in SELECTOR_PROFILE_GROUPS.get(profile, ()):
        overrides.update(
            has_busco_results=None,
            missing_busco_results=None,
            busco_complete_min=None,
            busco_single_min=None,
        )
    if profile != "assembly_with_exclusions":
        overrides["allow_all"] = False
    elif getattr(args, "all", False):
        overrides["allow_all"] = True
    if busco_library_id is not None:
        overrides["busco_library_id"] = busco_library_id
    request = request.with_overrides(**overrides) if overrides else request
    preset_name = getattr(args, "preset_name", None)
    if preset_name:
        if manager is None:
            raise ValueError("--preset requires a database-backed selector context.")
        request = merge_selector_preset(manager, preset_name, prune_selector_mapping(request.as_mapping()))
    return request


def _add_basic_list_output_options(group: argparse._ArgumentGroup) -> None:
    """Register shared formatting options for simple list-style outputs."""

    group.add_argument(
        *_option_strings("--tidy", *LIST_OUTPUT_SHORT_ALIASES["tidy"]),
        action="store_true",
        help="Align columns instead of TSV output.",
    )
    group.add_argument(
        *_option_strings("--pretty", *LIST_OUTPUT_SHORT_ALIASES["pretty"]),
        dest="list_color",
        action="store_true",
        help="Render pretty colored list output (TTY only).",
    )
    group.add_argument(
        "--no-pretty",
        dest="list_color",
        action="store_false",
        help="Disable pretty colored list output.",
    )
    group.add_argument(
        "--no-pager",
        "--no-pagination",
        dest="no_pager",
        action="store_true",
        help="Print the complete pretty table without interactive pagination.",
    )
    group.add_argument(
        "--no-header",
        action="store_true",
        help="Omit the column header row from list output.",
    )
    group.add_argument(
        "--colour",
        dest="list_color",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    group.add_argument(
        "--no-colour",
        dest="list_color",
        action="store_false",
        help=argparse.SUPPRESS,
    )


def _add_assembly_list_output_options(group: argparse._ArgumentGroup) -> None:
    """Register assembly-specific list output options."""

    group.add_argument(
        *_option_strings("--busco", *LIST_OUTPUT_SHORT_ALIASES["busco"]),
        action="store_true",
        help="Include BUSCO summary columns in list assemblies.",
    )
    group.add_argument(
        "--all-runs",
        action="store_true",
        help="With --busco, expand one assembly row into one row per BUSCO run.",
    )
    group.add_argument("--group-by-rank", action="store_true", help="Group list assemblies output by rank header.")
    _add_basic_list_output_options(group)
    group.add_argument("--output-path", help="Write list assemblies output to a file.")
    group.add_argument(
        *LIST_OUTPUT_SHORT_ALIASES["store"],
        "--store",
        "--save-set",
        dest="store_variable",
        help="Store resolved accessions in a named set (use with @NAME later).",
    )
    group.add_argument(
        *LIST_OUTPUT_SHORT_ALIASES["append"],
        "--append-to",
        "--append-set",
        dest="append_to_variable",
        help="Union resolved accessions into a named set.",
    )
    group.add_argument(
        "--intersection",
        action=AppendCommaSeparated,
        help="Intersect resolved accessions with an explicit accession or @VARIABLE set.",
    )
    group.add_argument("--one-line", action="store_true", help="Emit a comma-separated accession list.")
    group.add_argument(
        "--one-line-quotes",
        action="store_true",
        help="Emit a comma-separated accession list wrapped in quotes.",
    )
    group.add_argument(
        "--extended-decontamination-headers",
        action="store_true",
        help="Include support/weak/unknown decontamination percentages in list output.",
    )


__all__ = [
    "SELECTOR_PROFILE_GROUPS",
    "_add_assembly_list_output_options",
    "_add_basic_list_output_options",
    "_add_selector_arguments",
    "_selector_request_from_args",
]
