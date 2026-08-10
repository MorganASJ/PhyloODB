"""CLI registration and handlers for list, count, and assemblies commands."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ...database import DBManager
from ...proteome_profile_utils import RAW_PROFILE
from ...registry import registry
from ...selector_utils import (
    collect_accession_candidates,
    _evaluate_filter_condition,
    _parse_filter_expression,
    expand_accession_variables,
    expand_busco_run_id_variables,
    filter_accessions_by_busco_selectors,
    normalize_rank_list,
    RANK_HIERARCHY,
    resolve_clade_to_taxid,
    resolve_selector_accessions,
)
from ...variable_kinds import (
    infer_variable_kind,
    normalize_variable_kind,
)
from ...variable_definitions import build_variables_json_document
from ..support.output import TableColumn, TableData
from ..support.task_rows import collect_error_rows, collect_queue_rows, QUEUE_HEADERS
from .storage import _handle_list_roots
from ..support.argparse_utils import AppendCommaSeparated
from ..support.common import (
    STORAGE_ROOT_KINDS,
    _apply_busco_context_from_args,
    _coerce_bool,
    _connect_manager,
    _expand_accessions,
    _load_list_color_defaults,
    _print_error,
    _resolve_library_selector,
)
from ..support.output import _coerce_float, _color_from_gradient, _fmt_busco_value, _parse_color_tokens, _parse_float_tokens, _render_grouped_rows, _render_grouped_rows_rich, _render_list_output
from ..support.selectors import (
    LIST_OUTPUT_SHORT_ALIASES,
    _add_assembly_list_output_options,
    _add_basic_list_output_options,
    _add_selector_arguments,
    _selector_request_from_args,
)

def _add_row_sort_option(group: argparse._ArgumentGroup, *, fields: str) -> None:
    group.add_argument(
        "-s",
        "--sort",
        help=(
            "Sort rows by FIELD[:asc|desc][,FIELD[:asc|desc]...]. "
            f"Available fields: {fields}."
        ),
    )


def _normalize_sort_field(value: str) -> str:
    return str(value or "").strip().lower().replace("_", ".").replace("-", ".")


def _sort_value(value: Any) -> tuple[int, int, Any]:
    text = str(value or "").strip()
    if text == "" or text.upper() == "NA":
        return (1, 2, "")
    try:
        return (0, 0, float(text.replace(",", "")))
    except ValueError:
        return (0, 1, text.lower())


class _DescendingSortValue:
    def __init__(self, value: Any):
        self.value = value

    def __lt__(self, other: "_DescendingSortValue") -> bool:
        return self.value > other.value


def _parse_row_sort(sort_spec: Any, headers: Sequence[str], aliases: Mapping[str, str] | None = None):
    if not sort_spec:
        return []
    aliases = aliases or {}
    header_map = {_normalize_sort_field(header): idx for idx, header in enumerate(headers)}
    alias_map: Dict[str, tuple[int, Optional[str]]] = {}
    for alias, target in aliases.items():
        target_text = str(target or "").strip()
        target_direction: Optional[str] = None
        if ":" in target_text:
            target_field, target_direction_raw = target_text.rsplit(":", 1)
            target_direction_candidate = target_direction_raw.strip().lower()
            if target_direction_candidate in {"asc", "desc"}:
                target_text = target_field
                target_direction = target_direction_candidate
        target_key = _normalize_sort_field(target_text)
        if target_key in header_map:
            alias_map[_normalize_sort_field(alias)] = (header_map[target_key], target_direction)
    parsed = []
    for raw_token in str(sort_spec).split(","):
        token = raw_token.strip()
        if not token:
            continue
        explicit_direction = False
        if ":" in token:
            field, direction = token.rsplit(":", 1)
            direction = direction.strip().lower()
            explicit_direction = True
        else:
            field, direction = token, "asc"
        if direction not in {"asc", "desc"}:
            raise ValueError(f"Invalid sort direction '{direction}' for '{field}'. Use asc or desc.")
        key = _normalize_sort_field(field)
        if key in alias_map:
            idx, alias_direction = alias_map[key]
            if alias_direction and not explicit_direction:
                direction = alias_direction
        elif key in header_map:
            idx = header_map[key]
        else:
            available = ", ".join(headers)
            raise ValueError(f"Unknown sort field '{field}'. Available fields: {available}.")
        parsed.append((idx, direction == "desc"))
    return parsed


def _sort_rows(rows: Sequence[Sequence[Any]], headers: Sequence[str], sort_spec: Any, aliases: Mapping[str, str] | None = None) -> List[Sequence[Any]]:
    parsed = _parse_row_sort(sort_spec, headers, aliases)
    sorted_rows = list(rows)
    for idx, descending in reversed(parsed):
        def key(row: Sequence[Any]):
            missing, type_rank, value = _sort_value(row[idx] if idx < len(row) else "")
            return (missing, type_rank, _DescendingSortValue(value) if descending else value)

        sorted_rows.sort(key=key)
    return sorted_rows


def _apply_row_sort_or_error(
    rows: Sequence[Sequence[Any]],
    headers: Sequence[str],
    args: argparse.Namespace,
    aliases: Mapping[str, str] | None = None,
):
    try:
        return _sort_rows(rows, headers, getattr(args, "sort", None), aliases)
    except ValueError as exc:
        return _print_error(str(exc))


def _list_requires_write(args: argparse.Namespace) -> bool:
    """Return whether the requested listing action persists any state."""

    return bool(
        getattr(args, "store_variable", None)
        or getattr(args, "append_to_variable", None)
        or getattr(args, "store_results", None)
    )


def _classify_variable_kind(value: Any, known_accessions: set[str], name: Any = "") -> str:
    """Classify stored variables as assembly selector sets or environment config."""

    return infer_variable_kind(name, value, known_accessions)


def _chunked(values: Sequence[str], size: int) -> Iterable[List[str]]:
    """Yield a sequence in fixed-size chunks for batch database queries."""

    if size <= 0:
        raise ValueError("chunk size must be positive")
    for idx in range(0, len(values), size):
        yield list(values[idx:idx + size])


def _normalize_variable_target(raw_name: Any, *, option: str, uppercase: bool = False) -> str:
    """Normalize CLI variable target names and reject selector-only syntax."""

    name = str(raw_name or "").strip()
    if name.startswith("@"):
        name = name[1:].strip()
    if uppercase:
        name = name.upper()
    if not name:
        raise ValueError(f"{option} requires a non-empty variable name.")
    if "@" in name:
        raise ValueError(f"{option} variable names cannot contain '@'; use @NAME only when referencing an existing variable.")
    return name


def _store_accession_variable(manager: DBManager, raw_name: Any, values: Sequence[str], *, append: bool = False) -> str:
    """Write an accession variable, optionally unioning into the existing stored order."""

    option = "--append-to" if append else "--store"
    var_name = _normalize_variable_target(raw_name, option=option)
    stored = list(dict.fromkeys(values))
    if append:
        existing_value = manager.get_environment_variable(var_name)
        existing_items: List[str] = []
        if existing_value is None:
            existing_items = []
        elif isinstance(existing_value, (list, tuple, set)):
            existing_items = [str(item) for item in existing_value if item is not None]
        elif isinstance(existing_value, str):
            existing_items = [part.strip() for part in existing_value.split(",") if part.strip()]
        else:
            raise ValueError(f"Accession variable '{var_name}' must be a list or comma-separated string.")
        stored = list(dict.fromkeys([*existing_items, *stored]))
    if not manager.set_environment_variable(var_name, stored, kind="assemblies"):
        raise ValueError(f"Failed to store accessions in variable '{var_name}'.")
    return var_name


def _apply_accession_set_ops(manager: DBManager, args: argparse.Namespace, selected: Sequence[str]) -> List[str]:
    """Apply post-resolution set operations to a resolved accession list."""

    result = list(dict.fromkeys(selected))
    intersection_tokens = getattr(args, "intersection", None) or []
    if intersection_tokens:
        rhs = set(expand_accession_variables(manager, intersection_tokens, allow_bare=False))
        result = [acc for acc in result if acc in rhs]
    return result

METADATA_FIELDS: Dict[str, Dict[str, str]] = {
    "release_date": {"label": "release_date", "sql": "a.release_date", "source": "Assembly", "description": "Assembly release date."},
    "origin": {"label": "origin", "sql": "a.origin", "source": "Assembly", "description": "Assembly origin (for example local, refseq, genbank)."},
    "level": {"label": "level", "sql": "g.assembly_level", "source": "Genome", "description": "Assembly level."},
    "n50": {"label": "n50", "sql": "a.contig_n50", "source": "Assembly", "description": "Contig N50."},
    "comments": {"label": "comments", "sql": "a.comments", "source": "Assembly", "description": "Assembly comments."},
    "submitter": {"label": "submitter", "sql": "a.submitter", "source": "Assembly", "description": "Submitter name."},
    "assembly_method": {"label": "assembly_method", "sql": "a.assembly_method", "source": "Assembly", "description": "Assembly method."},
    "assembly_type": {"label": "assembly_type", "sql": "a.assembly_type", "source": "Assembly", "description": "Assembly type."},
    "assembly_status": {"label": "assembly_status", "sql": "a.assembly_status", "source": "Assembly", "description": "Assembly status."},
    "warnings": {"label": "warnings", "sql": "a.warnings", "source": "Assembly", "description": "Assembly warnings."},
    "bioproject_accession": {"label": "bioproject_accession", "sql": "a.bioproject_accession", "source": "Assembly", "description": "BioProject accession."},
    "biosample_accession": {"label": "biosample_accession", "sql": "a.biosample_accession", "source": "Assembly", "description": "BioSample accession."},
    "diploid_role": {"label": "diploid_role", "sql": "a.diploid_role", "source": "Assembly", "description": "Diploid role."},
    "refseq_category": {"label": "refseq_category", "sql": "a.refseq_category", "source": "Assembly", "description": "RefSeq category."},
    "sequencing_tech": {"label": "sequencing_tech", "sql": "a.sequencing_tech", "source": "Assembly", "description": "Sequencing technology."},
    "contig_l50": {"label": "contig_l50", "sql": "a.contig_l50", "source": "Assembly", "description": "Contig L50."},
    "gc_count": {"label": "gc_count", "sql": "a.gc_count", "source": "Assembly", "description": "GC count."},
    "gc_percent": {"label": "gc_percent", "sql": "a.gc_percent", "source": "Assembly", "description": "GC percent."},
    "genome_coverage": {"label": "genome_coverage", "sql": "a.genome_coverage", "source": "Assembly", "description": "Genome coverage."},
    "number_of_component_sequences": {"label": "number_of_component_sequences", "sql": "a.number_of_component_sequences", "source": "Assembly", "description": "Component sequence count."},
    "number_of_contigs": {"label": "number_of_contigs", "sql": "a.number_of_contigs", "source": "Assembly", "description": "Contig count."},
    "number_of_organelles": {"label": "number_of_organelles", "sql": "a.number_of_organelles", "source": "Assembly", "description": "Organelle count."},
    "number_of_scaffolds": {"label": "number_of_scaffolds", "sql": "a.number_of_scaffolds", "source": "Assembly", "description": "Scaffold count."},
    "scaffold_l50": {"label": "scaffold_l50", "sql": "a.scaffold_l50", "source": "Assembly", "description": "Scaffold L50."},
    "scaffold_n50": {"label": "scaffold_n50", "sql": "a.scaffold_n50", "source": "Assembly", "description": "Scaffold N50."},
    "total_number_of_chromosomes": {"label": "total_number_of_chromosomes", "sql": "a.total_number_of_chromosomes", "source": "Assembly", "description": "Chromosome count."},
    "total_sequence_length": {"label": "total_sequence_length", "sql": "a.total_sequence_length", "source": "Assembly", "description": "Total sequence length."},
    "total_ungapped_length": {"label": "total_ungapped_length", "sql": "a.total_ungapped_length", "source": "Assembly", "description": "Total ungapped length."},
    "assembly_name": {"label": "assembly_name", "sql": "g.assembly_name", "source": "Genome", "description": "Assembly name."},
    "assembly_properties": {"label": "assembly_properties", "sql": "g.assembly_properties", "source": "Genome", "description": "Assembly properties."},
    "protein": {"label": "protein", "sql": "g.protein", "source": "Genome", "description": "Protein flag."},
    "isoforms_cleaned": {"label": "isoforms_cleaned", "sql": "g.isoforms_cleaned", "source": "Genome", "description": "Isoform-cleaned flag (0/1)."},
    "proteome_profile": {"label": "proteome_profile", "sql": "", "source": "BUSCORun", "description": "Proteome profile used by the listed BUSCO run when BUSCO rows are shown."},
    "orthofinder_target_library": {"label": "orthofinder_target_library", "sql": "", "source": "BUSCORun", "description": "Library constructed using the listed OrthoFinder-derived BUSCO run."},
    "default_proteome_profile": {"label": "default_proteome_profile", "sql": "", "source": "ProteomeProfile", "description": "Default proteome profile for the accession."},
    "genome_comments": {"label": "genome_comments", "sql": "g.comments", "source": "Genome", "description": "Genome comments."},
    "download_date": {"label": "download_date", "sql": "g.dl_date", "source": "Genome", "description": "Download date."},
    "location": {"label": "location", "sql": "g.location", "source": "Genome", "description": "Local storage path."},
    "status": {"label": "status", "sql": "g.status", "source": "Genome", "description": "Download status."},
}

METADATA_ALIASES = {
    "n50": "n50",
    "contig_n50": "n50",
    "level": "level",
    "assembly_level": "level",
    "release": "release_date",
    "release_date": "release_date",
    "origin": "origin",
    "source": "origin",
    "comments": "comments",
    "submitter": "submitter",
    "isoforms_cleaned": "isoforms_cleaned",
    "cleaned_isoforms": "isoforms_cleaned",
    "cleaned_isoform": "isoforms_cleaned",
    "isoform_cleaned": "isoforms_cleaned",
    "proteome_profile": "proteome_profile",
    "profile": "proteome_profile",
    "orthofinder_target_library": "orthofinder_target_library",
    "of_target_library": "orthofinder_target_library",
    "default_proteome_profile": "default_proteome_profile",
    "default_profile": "default_proteome_profile",
}

BUSCO_METADATA_TERMS: Dict[str, Dict[str, str]] = {
    "busco.complete": {
        "label": "busco.complete",
        "aliases": "complete",
    },
    "busco.single_copy_complete": {
        "label": "busco.single_copy_complete",
        "aliases": "quality, busco.single_copy, busco.sc, single_copy, single_copy_complete, sc",
    },
    "busco.duplicated": {
        "label": "busco.duplicated",
        "aliases": "duplicated",
    },
    "busco.fragmented": {
        "label": "busco.fragmented",
        "aliases": "fragmented",
    },
    "busco.missing": {
        "label": "busco.missing",
        "aliases": "missing",
    },
}

ASSEMBLY_SORT_ALIASES = {
    "taxon": "species",
    "name": "species",
    "busco.complete": "complete",
    "busco.single_copy": "single_copy",
    "busco.single.copy": "single_copy",
    "busco.single_copy_complete": "single_copy",
    "busco.single.copy.complete": "single_copy",
    "busco.sc": "single_copy",
    "quality": "single_copy:desc",
    "busco.quality": "single_copy:desc",
    "busco.duplicated": "duplicated",
    "busco.fragmented": "fragmented",
    "busco.missing": "missing",
    "paralog.hidden": "hidden_paralog",
    "decontamination.contaminated": "contaminated",
    "release": "release_date",
    "date": "release_date",
    "latest": "release_date:desc",
    "newest": "release_date:desc",
    "oldest": "release_date:asc",
    "level": "assembly_level",
    "status": "status",
}


def _metadata_aliases_by_target() -> Dict[str, List[str]]:
    aliases_by_target: Dict[str, List[str]] = {key: [] for key in METADATA_FIELDS}
    for alias, target in sorted(METADATA_ALIASES.items()):
        if alias != target and target in aliases_by_target:
            aliases_by_target[target].append(alias)
    return aliases_by_target


def _split_aliases(value: str) -> List[str]:
    return [alias.strip() for alias in value.split(",") if alias.strip()]


def _format_term_line(term: str, aliases: Sequence[str]) -> str:
    unique_aliases = list(dict.fromkeys(alias for alias in aliases if alias and alias != term))
    if not unique_aliases:
        return f"  {term}"
    return f"  {term} (aliases: {', '.join(unique_aliases)})"


def _metadata_term_sections() -> List[Tuple[str, List[Tuple[str, List[str]]]]]:
    aliases_by_target = _metadata_aliases_by_target()
    metadata_terms = [
        (METADATA_FIELDS[key]["label"], aliases_by_target.get(key, []))
        for key in sorted(METADATA_FIELDS)
    ]
    busco_terms = [
        (field["label"], _split_aliases(field.get("aliases", "")))
        for _, field in sorted(BUSCO_METADATA_TERMS.items())
    ]
    return [("Metadata", metadata_terms), ("BUSCO", busco_terms)]


def _metadata_sort_profiles() -> List[Tuple[str, str]]:
    profiles = []
    for alias, target in sorted(ASSEMBLY_SORT_ALIASES.items()):
        if ":" in target:
            display_target = target
            if target.startswith("single_copy:"):
                display_target = target.replace("single_copy", "busco.single_copy_complete", 1)
            profiles.append((alias, display_target))
    return profiles


def _print_metadata_terms() -> int:
    print("Accepted list assemblies --filter / --sort terms.")
    print("Use canonical terms when possible; aliases in parentheses are accepted shortcuts.")
    print()

    for heading, terms in _metadata_term_sections():
        print(f"{heading}:")
        for term, aliases in terms:
            print(_format_term_line(term, aliases))
        print()

    print("Sort profiles:")
    for alias, target in _metadata_sort_profiles():
        print(f"  {alias} -> {target}")
    return 0


# ---------------------------------------------------------------------------
# Assembly metadata and BUSCO lookup helpers
# Purpose: Keep data-fetching helpers close to the list/count handlers that use
# them so assembly reporting logic is easier to follow.
# ---------------------------------------------------------------------------

def _normalize_meta_token(token: str) -> str:
    """Normalise user-facing metadata tokens into canonical field keys."""

    key = str(token or "").strip().lower()
    key = key.replace(" ", "_")
    key = key.replace("-", "_")
    if key in METADATA_ALIASES:
        return METADATA_ALIASES[key]
    return key


def _split_meta_tokens(value: Any) -> List[str]:
    """Split metadata option values into individual field tokens."""

    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        items = value
    else:
        items = [value]
    tokens: List[str] = []
    for entry in items:
        if entry is None:
            continue
        if isinstance(entry, str):
            parts = [part.strip() for part in entry.split(",")]
            tokens.extend([part for part in parts if part])
        else:
            tokens.append(str(entry))
    return tokens


def _resolve_meta_fields(tokens: Sequence[Any]) -> List[str]:
    """Validate requested metadata fields and preserve user order."""

    requested = _split_meta_tokens(tokens)
    fields: List[str] = []
    seen: set[str] = set()
    for token in requested:
        key = _normalize_meta_token(token)
        if key in METADATA_ALIASES:
            key = METADATA_ALIASES[key]
        if key not in METADATA_FIELDS:
            raise ValueError(f"Unknown metadata field '{token}'.")
        if key in seen:
            continue
        seen.add(key)
        fields.append(key)
    return fields


def _fetch_metadata(
    manager: DBManager,
    accessions: Sequence[str],
    fields: Sequence[str],
) -> Dict[str, Dict[str, str]]:
    """Fetch selected assembly metadata columns for a set of accessions."""

    if not accessions or not fields:
        return {}
    data: Dict[str, Dict[str, str]] = {}
    selected_fields = [METADATA_FIELDS[field] for field in fields if field in METADATA_FIELDS]
    sql_fields = [field for field in selected_fields if field.get("sql")]
    custom_default_profile_field = any(field["label"] in {"proteome_profile", "default_proteome_profile"} for field in selected_fields)

    if sql_fields:
        cols = ", ".join(f"{field['sql']} AS {field['label']}" for field in sql_fields)
        for chunk in _chunked(list(accessions), 900):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"""
                SELECT g.accession, {cols}
                FROM Genome g
                LEFT JOIN Assembly a ON a.accession = g.accession
                WHERE g.accession IN ({placeholders})
                """,
                tuple(chunk),
            )
            for row in manager.cursor.fetchall() or []:
                acc = str(row[0])
                values = data.setdefault(acc, {})
                for idx, field in enumerate(sql_fields, start=1):
                    val = row[idx]
                    values[field["label"]] = "" if val is None else str(val)

    if custom_default_profile_field:
        profile_info = _fetch_proteome_profile_display_info(manager, accessions)
        for acc in accessions:
            values = data.setdefault(str(acc), {})
            default_profile = ""
            for _profile_name, details in (profile_info.get(str(acc), {}) or {}).items():
                if details.get("default"):
                    default_profile = str(details.get("display") or "")
                    break
            values["proteome_profile"] = default_profile
            values["default_proteome_profile"] = default_profile
    return data


def _format_cdhit_label(identity: Optional[Any]) -> str:
    if identity is None:
        return "cdhit"
    try:
        value = float(identity)
    except (TypeError, ValueError):
        return f"cdhit{identity}"
    pct = value * 100.0 if value <= 1.0 else value
    if abs(pct - round(pct)) < 1e-6:
        return f"cdhit{int(round(pct))}"
    return f"cdhit{pct:.2f}".rstrip("0").rstrip(".")


def _profile_display_base(profile_name: str, prep_row: Optional[Sequence[Any]]) -> str:
    profile = str(profile_name or "").strip()
    if not profile:
        return ""
    if profile == RAW_PROFILE:
        return RAW_PROFILE
    if not prep_row:
        return profile
    used_gff = bool(prep_row[5]) if len(prep_row) > 5 else False
    skip_cdhit = bool(prep_row[8]) if len(prep_row) > 8 else False
    cdhit_identity = prep_row[10] if len(prep_row) > 10 else None
    tokens: List[str] = []
    if used_gff:
        tokens.append("gff")
    if not skip_cdhit:
        tokens.append(_format_cdhit_label(cdhit_identity))
    return ",".join(tokens) if tokens else profile


def _fetch_proteome_profile_display_info(
    manager: DBManager,
    accessions: Sequence[str],
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    info: Dict[str, Dict[str, Dict[str, Any]]] = {}
    if not accessions:
        return info
    for row in manager.proteomes.list_profiles(accessions=list(accessions)) or []:
        profile_id = int(row[0])
        accession = str(row[1])
        profile_name = str(row[2] or "")
        status = str(row[6] or "")
        is_default = bool(row[9])
        if not profile_name or status != "ready":
            continue
        prep = manager.proteomes.latest_preparation_for_output(profile_id)
        base = _profile_display_base(profile_name, prep)
        info.setdefault(accession, {})[profile_name] = {
            "base": base,
            "default": is_default,
            "display": f"{base}*" if is_default and base else base,
        }
    return info


def _fetch_assembly_info(manager: DBManager, accessions: Sequence[str]) -> Dict[str, Tuple[Optional[int], str]]:
    """Fetch taxid and species name for the supplied assembly accessions."""

    info: Dict[str, Tuple[Optional[int], str]] = {}
    if not accessions:
        return info
    for chunk in _chunked(list(accessions), 900):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"""
            SELECT taxid, name, assembly_accession
            FROM TaxonomyAssemblySummary
            WHERE assembly_accession IN ({placeholders})
            """,
            tuple(chunk),
        )
        for taxid, name, accession in manager.cursor.fetchall() or []:
            info[str(accession)] = (int(taxid) if taxid is not None else None, str(name or ""))
    return info


def _fetch_busco_scores(
    manager: DBManager,
    accessions: Sequence[str],
    library_id: int,
    *,
    include_paralog_filtering_in_score: Optional[bool] = None,
    include_decontamination_in_score: Optional[bool] = None,
    paralog_run_id: Optional[str] = None,
    decontamination_run_id: Optional[str] = None,
    allow_ambiguous_contaminants: Optional[bool] = None,
    strict_decontamination: Optional[bool] = None,
    rescue_duplicates: Optional[bool] = None,
) -> Dict[str, Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]]:
    """Fetch adjusted BUSCO summary values keyed by assembly accession."""

    scores: Dict[str, Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]] = {}
    if not accessions:
        return scores
    for chunk in _chunked(list(accessions), 900):
        rows = manager.get_busco_results_adjusted(
            library_id=library_id,
            accessions=chunk,
            include_paralog=include_paralog_filtering_in_score,
            paralog_run_id=paralog_run_id,
            include_decontam=include_decontamination_in_score,
            decont_run_id=decontamination_run_id,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )
        for row in rows or []:
            acc = str(row[0])
            scores[acc] = (
                row[3],
                row[4],
                row[5],
                row[6],
                row[7],
                row[8],
                row[9],
            )
    return scores


def _busco_run_tag(pipeline: Optional[str], input_mode: Optional[str], supports_nt: bool) -> str:
    """Build a compact BUSCO run tag for assembly list output."""

    mode = (input_mode or "").strip().lower()
    pipe = (pipeline or "").strip().lower()
    if mode == "protein":
        if pipe == "orthofinder":
            return "orthofinder-faa"
        return "proteome-faa"
    if mode == "genome":
        suffix = "fna" if supports_nt else "faa"
    else:
        suffix = "faa"
    return f"{pipe}-{suffix}"


def _busco_mode_code(pipeline: Optional[str], input_mode: Optional[str], *, is_primary: bool = False) -> str:
    """Return the short pipeline code shown in ``list assemblies --all-runs``."""

    mode = (input_mode or "").strip().lower()
    pipe = (pipeline or "").strip().lower()
    suffix = "*" if is_primary else ""
    if mode == "protein":
        if pipe == "orthofinder":
            return f"O{suffix}"
        return f"P{suffix}"
    if pipe == "augustus":
        return f"Au{suffix}"
    if pipe == "miniprot":
        return f"Mi{suffix}"
    if pipe == "metaeuk":
        return f"Me{suffix}"
    return suffix


def _busco_pipeline_alias(pipeline: Optional[str], input_mode: Optional[str]) -> str:
    mode = (input_mode or "").strip().lower()
    pipe = (pipeline or "").strip().lower()
    if mode == "protein":
        if pipe == "orthofinder":
            return "orthofinder"
        return "proteome"
    return pipe


def _busco_pipeline_display(pipeline: Optional[str], input_mode: Optional[str]) -> str:
    alias = _busco_pipeline_alias(pipeline, input_mode)
    if alias == "proteome":
        return "Proteome"
    if alias == "orthofinder":
        return "OrthoFinder"
    return str(pipeline or "")


def _busco_format_code(
    input_mode: Optional[str],
    *,
    pipeline: Optional[str] = None,
    assembly_has_protein: Optional[bool] = None,
) -> str:
    """Return the short BUSCO input-format code used in list output."""

    mode = (input_mode or "").strip().lower()
    if mode == "genome":
        return "G"
    if mode == "protein":
        pipe = (pipeline or "").strip().lower()
        # Legacy migrated miniprot rows can be tagged as protein even when
        # the assembly is genome-only; present those as genome mode.
        if pipe == "miniprot" and assembly_has_protein is False:
            return "G"
        return "P"
    return "P"


def _fetch_genome_protein_flags(manager: DBManager, accessions: Sequence[str]) -> Dict[str, Optional[bool]]:
    """Fetch per-accession protein availability for BUSCO run labelling."""

    result: Dict[str, Optional[bool]] = {}
    if not accessions:
        return result
    for chunk in _chunked(list(accessions), 900):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"SELECT accession, protein FROM Genome WHERE accession IN ({placeholders})",
            tuple(chunk),
        )
        for acc, protein in manager.cursor.fetchall() or []:
            if protein is None:
                result[str(acc)] = None
            else:
                try:
                    result[str(acc)] = bool(int(protein))
                except (TypeError, ValueError):
                    result[str(acc)] = bool(protein)
    return result


def _fetch_busco_runs(
    manager: DBManager,
    accessions: Sequence[str],
    library_id: int,
    *,
    purpose: str = "default",
    run_ids: Optional[Sequence[int]] = None,
    pipeline: Optional[str] = None,
    input_mode: Optional[str] = None,
    preferred_pipeline: Optional[str] = None,
    preferred_input_mode: Optional[str] = None,
    proteome_profile: Optional[str] = None,
    preferred_proteome_profile: Optional[str] = None,
    export_format: Optional[str] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """Fetch BUSCO run history for the supplied accessions."""

    result: Dict[str, List[Dict[str, Any]]] = {}
    allowed_run_ids = {int(run_id) for run_id in (run_ids or []) if run_id is not None}
    requested_pipeline = str(pipeline).strip().lower() if pipeline else None
    requested_input_mode = str(input_mode).strip().lower() if input_mode else None
    preferred_pipeline_norm = str(preferred_pipeline).strip().lower() if preferred_pipeline else None
    preferred_input_mode_norm = str(preferred_input_mode).strip().lower() if preferred_input_mode else None
    requested_profile = str(proteome_profile).strip() if proteome_profile else None
    preferred_profile_norm = str(preferred_proteome_profile).strip() if preferred_proteome_profile else None
    export_purpose = None
    if export_format:
        token = str(export_format).strip().lower()
        export_purpose = "export_nucleotide" if token == "nucleotide" else "export_protein"
    rows = manager.get_busco_runs_for_accessions(accessions, library_id=library_id, purpose=purpose)
    for row in rows or []:
        (
            run_id,
            accession,
            _lib_id,
            _lib_name,
            _lineage_name,
            row_input_mode,
            row_pipeline,
            _result_dir,
            status,
            _sc,
            _dup,
            _frag,
            _miss,
            completed_at,
            is_primary,
            complete,
            single_copy,
            duplicated,
            fragmented,
            missing,
            proteome_profile_name,
            is_default_profile,
        ) = row
        if allowed_run_ids and int(run_id) not in allowed_run_ids:
            continue
        if export_purpose and not manager.busco.run_supports_purpose(int(run_id), purpose=export_purpose):
            continue
        pipe_norm = str(row_pipeline or "").strip().lower()
        supports_nt = pipe_norm in {"augustus", "metaeuk"}
        display_pipeline = _busco_pipeline_display(row_pipeline, row_input_mode)
        alias_pipeline = _busco_pipeline_alias(row_pipeline, row_input_mode)
        entry = {
            "run_id": run_id,
            "pipeline": display_pipeline,
            "pipeline_raw": row_pipeline,
            "pipeline_alias": alias_pipeline,
            "input_mode": row_input_mode,
            "seq_support": "aa+nt" if supports_nt else "aa-only",
            "run_tag": _busco_run_tag(row_pipeline, row_input_mode, supports_nt),
            "is_primary": bool(is_primary),
            "status": status,
            "completed_at": completed_at,
            "complete": complete,
            "single_copy": single_copy,
            "duplicated": duplicated,
            "fragmented": fragmented,
            "missing": missing,
            "proteome_profile": proteome_profile_name,
            "is_default_profile": bool(is_default_profile),
        }
        if requested_pipeline and str(entry.get("pipeline_alias") or "").strip().lower() != requested_pipeline:
            continue
        if requested_input_mode and str(entry.get("input_mode") or "").strip().lower() != requested_input_mode:
            continue
        if requested_profile and not manager.proteomes.profile_matches_selector(
            str(accession),
            str(entry.get("proteome_profile") or "").strip() or None,
            requested_profile,
        ):
            continue
        result.setdefault(str(accession), []).append(entry)
    for accession, entries in result.items():
        result[accession] = sorted(
            entries,
            key=lambda entry: (
                1 if preferred_pipeline_norm and str(entry.get("pipeline_alias") or "").strip().lower() == preferred_pipeline_norm else 0,
                1 if preferred_input_mode_norm and str(entry.get("input_mode") or "").strip().lower() == preferred_input_mode_norm else 0,
                1 if preferred_profile_norm and str(entry.get("proteome_profile") or "").strip() == preferred_profile_norm else 0,
                1 if entry.get("is_primary") else 0,
                str(entry.get("completed_at") or ""),
                int(entry.get("run_id") or 0),
            ),
            reverse=True,
        )
    return result


def _fetch_orthofinder_target_libraries(
    manager: DBManager,
    run_ids: Sequence[int],
) -> Dict[int, str]:
    """Return target-library provenance for OrthoFinder-derived BUSCO runs."""

    result: Dict[int, str] = {}
    unique_run_ids = sorted({int(run_id) for run_id in run_ids if run_id is not None})
    for chunk in _chunked(unique_run_ids, 900):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"""
            SELECT run_id, pipeline_params_effective_json
            FROM BUSCO_Runs
            WHERE pipeline = 'orthofinder'
              AND run_id IN ({placeholders})
            """,
            tuple(chunk),
        )
        for run_id, params_json in manager.cursor.fetchall() or []:
            try:
                params = json.loads(params_json) if params_json else {}
            except (TypeError, ValueError, json.JSONDecodeError):
                params = {}
            target_library = str(params.get("derived_library_name") or "").strip()
            if target_library:
                result[int(run_id)] = target_library
    return result


def _fetch_decontamination_breakdown(
    manager: DBManager,
    accessions: Sequence[str],
    library_id: int,
    *,
    run_id: Optional[str] = None,
) -> Dict[str, Tuple[Optional[float], Optional[float], Optional[float]]]:
    """Fetch decontamination support percentages for list output."""

    if not accessions:
        return {}
    return manager.get_decontamination_decision_percentages(
        library_id=library_id,
        accessions=accessions,
        run_id=run_id,
    )


def _fetch_decontamination_summaries(
    manager: DBManager,
    accessions: Sequence[str],
    library_id: int,
    *,
    run_id: Optional[str] = None,
) -> Dict[str, Tuple[Optional[str], Optional[str], Optional[str]]]:
    """Fetch decontamination summary labels, with parent-library fallback."""

    if not accessions:
        return {}
    parent_id = manager.assert_library_has_parent(library_id)
    return manager.get_latest_decontamination_summary_with_fallback(
        target_library_id=library_id,
        parent_library_id=parent_id,
        accessions=accessions,
        run_id=run_id,
    )


def _fmt_busco_value(value: Optional[float]) -> str:
    """Format BUSCO values consistently for terminal output."""

    if value is None:
        return "NA"
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return str(value)


def _infer_busco_library_from_run_ids(manager: DBManager, run_ids: Sequence[int]) -> Optional[int]:
    """Infer a single BUSCO library id from explicit run ids."""

    run_vals = [int(run_id) for run_id in run_ids or [] if run_id is not None]
    if not run_vals:
        return None
    placeholders = ",".join("?" for _ in run_vals)
    rows = manager.cursor.execute(
        f"SELECT DISTINCT library_id FROM BUSCO_Runs WHERE run_id IN ({placeholders})",
        tuple(run_vals),
    ).fetchall() or []
    library_ids = [int(row[0]) for row in rows if row and row[0] is not None]
    if not library_ids:
        return None
    if len(library_ids) > 1:
        raise ValueError(
            "Selected BUSCO run ids span multiple libraries. Use --library-id/--library-name to disambiguate."
        )
    return library_ids[0]


def _split_busco_run_metadata_filters(
    filters: Any,
) -> Tuple[List[str], List[List[List[Dict[str, Any]]]]]:
    """Separate run-provenance filters from accession-level selector filters."""

    ordinary: List[str] = []
    target_filters: List[List[List[Dict[str, Any]]]] = []
    for expression in filters or []:
        groups = _parse_filter_expression(str(expression))
        fields = {
            str(condition.get("field") or "").strip().lower()
            for group in groups
            for condition in group
        }
        target_names = {"orthofinder_target_library", "of_target_library"}
        if fields & target_names:
            if not fields <= target_names:
                raise ValueError(
                    "orthofinder_target_library cannot currently be combined with other fields "
                    "inside the same --filter expression; pass them as separate --filter options."
                )
            target_filters.append(groups)
        else:
            ordinary.append(str(expression))
    return ordinary, target_filters


def _matches_busco_run_metadata_filters(
    value: str,
    filters: Sequence[List[List[Dict[str, Any]]]],
) -> bool:
    """Evaluate repeated run-provenance filters using normal filter semantics."""

    return all(
        any(all(_evaluate_filter_condition(value, condition) for condition in group) for group in groups)
        for groups in filters
    )


def _handle_list_busco_runs(args: argparse.Namespace) -> int:
    """List BUSCO run records matching the active selectors."""

    manager = _connect_manager(args.database, read_only=not _list_requires_write(args))
    try:
        try:
            ordinary_filters, target_library_filters = _split_busco_run_metadata_filters(
                getattr(args, "filters", None) or []
            )
        except ValueError as exc:
            return _print_error(str(exc))
        selector_args = argparse.Namespace(**vars(args))
        selector_args.filters = ordinary_filters or None
        busco_library = _resolve_library_selector(
            manager,
            library_id=getattr(args, "library_id", None),
            library_name=getattr(args, "library_name", None),
            legacy=getattr(args, "busco_library", None),
        )
        selectors = _selector_request_from_args(
            selector_args,
            profile="busco_listing",
            busco_library_id=busco_library,
            manager=manager,
        )
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=True,
            require_candidates=False,
            use_rule_selection=False,
        )
        explicit_run_ids = expand_busco_run_id_variables(manager, getattr(args, "busco_run_ids", None) or [])
        params: list[Any] = []
        where = ["1=1"]
        if busco_library is not None:
            where.append("r.library_id = ?")
            params.append(int(busco_library))
        if selected:
            placeholders = ",".join("?" for _ in selected)
            where.append(f"r.accession IN ({placeholders})")
            params.extend([str(a) for a in selected])
        if explicit_run_ids:
            placeholders = ",".join("?" for _ in explicit_run_ids)
            where.append(f"r.run_id IN ({placeholders})")
            params.extend(explicit_run_ids)
        manager.cursor.execute(
            f"""
            SELECT r.run_id, r.accession, l.library_name, r.pipeline, r.input_mode,
                   r.status, r.completed_at, p.policy, pp.profile_name, COALESCE(pp.is_default, 0),
                   CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * (COALESCE(r.no_sc_complete,0)+COALESCE(r.no_duplicated_complete,0))/l.size, 2) END AS complete,
                   CASE WHEN l.size IS NULL OR l.size = 0 THEN NULL ELSE ROUND(100.0 * COALESCE(r.no_sc_complete,0)/l.size, 2) END AS single_copy,
                   CASE WHEN p.run_id IS NULL THEN 0 ELSE 1 END AS is_primary
            FROM BUSCO_Runs r
            JOIN Libraries l ON l.library_id = r.library_id
            LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
            LEFT JOIN BUSCO_Primary p
              ON p.run_id = r.run_id AND p.purpose = 'default'
            WHERE {' AND '.join(where)}
            ORDER BY r.accession, is_primary DESC, r.completed_at DESC, r.run_id DESC
            """,
            tuple(params),
        )
        rows = manager.cursor.fetchall() or []
        orthofinder_run_ids = [
            int(row[0])
            for row in rows
            if row and str(row[3] or "").strip().lower() == "orthofinder"
        ]
        target_libraries = (
            _fetch_orthofinder_target_libraries(manager, orthofinder_run_ids)
            if orthofinder_run_ids and (target_library_filters or not getattr(args, "ids_only", False))
            else {}
        )
        if target_library_filters:
            rows = [
                row
                for row in rows
                if _matches_busco_run_metadata_filters(
                    target_libraries.get(int(row[0]), ""),
                    target_library_filters,
                )
            ]
        run_accessions = [str(row[1]) for row in rows if row and row[1] is not None]
        assembly_info = _fetch_assembly_info(manager, run_accessions)
        profile_display_info = _fetch_proteome_profile_display_info(manager, run_accessions)
        out = []
        preferred_pipeline = str(getattr(args, "prefer_busco_pipeline", "") or "").strip().lower() or None
        preferred_input_mode = str(getattr(args, "prefer_busco_input_mode", None) or getattr(args, "prefer_format", None) or "").strip().lower() or None
        preferred_profile = str(getattr(args, "prefer_proteome_profile", "") or "").strip() or None
        required_input_mode = str(getattr(args, "busco_input_mode", None) or getattr(args, "format", None) or "").strip().lower() or None
        required_profile = str(getattr(args, "proteome_profile", "") or "").strip() or None
        if required_input_mode in {"nucleotide", "nucl"}:
            required_input_mode = "genome"
        if preferred_input_mode in {"nucleotide", "nucl"}:
            preferred_input_mode = "genome"
        filtered_rows = []
        for row in rows:
            run_id, acc, lib_name, pipeline, input_mode, status, completed_at, primary_policy, profile_name, is_default_profile, complete, single_copy, is_primary = row
            if required_input_mode and str(input_mode or "").strip().lower() != required_input_mode:
                continue
            if getattr(args, "busco_pipeline", None) and _busco_pipeline_alias(pipeline, input_mode) != str(getattr(args, "busco_pipeline")).strip().lower():
                continue
            if required_profile and not manager.proteomes.profile_matches_selector(
                str(acc),
                str(profile_name or "").strip() or None,
                required_profile,
            ):
                continue
            filtered_rows.append(row)
        rows = sorted(
            filtered_rows,
            key=lambda row: (
                1 if preferred_pipeline and _busco_pipeline_alias(row[3], row[4]) == preferred_pipeline else 0,
                1 if preferred_input_mode and str(row[4] or "").strip().lower() == preferred_input_mode else 0,
                1 if preferred_profile and str(row[8] or "").strip() == preferred_profile else 0,
                1 if row[12] else 0,
                str(row[6] or ""),
                int(row[0] or 0),
            ),
            reverse=True,
        )
        if not rows:
            print("No BUSCO runs found.")
            return 0
        run_ids_only = [str(int(row[0])) for row in rows if row and row[0] is not None]
        if getattr(args, "store_results", None):
            try:
                var_name = _normalize_variable_target(args.store_results, option="--store-results", uppercase=True)
            except ValueError as exc:
                return _print_error(str(exc))
            if not manager.set_environment_variable(var_name, run_ids_only, kind="busco_runs"):
                return _print_error(f"Failed to store BUSCO run ids in variable '{var_name}'.")
        if getattr(args, "ids_only", False):
            print(",".join(run_ids_only))
            return 0
        out = []
        for run_id, acc, lib_name, pipeline, input_mode, status, completed_at, primary_policy, profile_name, is_default_profile, complete, single_copy, is_primary in rows:
            species = assembly_info.get(str(acc), (None, ""))[1] or ""
            supports_nt = str(pipeline or "").lower() in {"augustus", "metaeuk"}
            run_tag = _busco_run_tag(pipeline, input_mode, supports_nt)
            input_mode_norm = str(input_mode or "").strip().lower()
            profile_display = ""
            details = profile_display_info.get(str(acc), {}).get(str(profile_name or ""))
            if details:
                profile_display = str(details.get("display") or "")
            elif profile_name:
                profile_display = str(profile_name)
                if is_default_profile:
                    profile_display += "*"
            elif input_mode_norm == "protein":
                profile_display = "unset"
            out.append(
                (
                    str(run_id),
                    str(acc),
                    str(species),
                    str(lib_name or ""),
                    _busco_pipeline_display(pipeline, input_mode),
                    target_libraries.get(int(run_id), ""),
                    str(input_mode or ""),
                    profile_display,
                    "1" if details and details.get("default") else ("1" if profile_display not in {"", "unset"} and is_default_profile else ""),
                    run_tag,
                    "1" if is_primary else "0",
                    str(primary_policy or ""),
                    str(status or ""),
                    _fmt_busco_value(complete),
                    _fmt_busco_value(single_copy),
                    str(completed_at or ""),
                )
            )
        headers = (
            "run_id",
            "accession",
            "species",
            "library",
            "pipeline",
            "orthofinder_target_library",
            "input_mode",
            "proteome_profile",
            "profile_default",
            "run_tag",
            "is_primary",
            "primary_policy",
            "status",
            "complete",
            "single_copy",
            "completed_at",
        )
        sorted_out = _apply_row_sort_or_error(
            out,
            headers,
            args,
            aliases={
                "id": "run_id",
                "lib": "library",
                "of_target_library": "orthofinder_target_library",
                "busco.complete": "complete",
                "busco.single_copy": "single_copy",
                "busco.single_copy_complete": "single_copy",
                "busco.sc": "single_copy",
                "quality": "single_copy:desc",
                "busco.quality": "single_copy:desc",
                "date": "completed_at",
            },
        )
        if isinstance(sorted_out, int):
            return sorted_out
        out = sorted_out
        return _render_list_output(args, headers, out, default_tidy=False)
    finally:
        manager.close()


def _handle_list_buscos(args: argparse.Namespace) -> int:
    """List BUSCO family rows for matching runs and assemblies."""

    manager = _connect_manager(args.database, read_only=True)
    try:
        busco_library = _resolve_library_selector(
            manager,
            library_id=getattr(args, "library_id", None),
            library_name=getattr(args, "library_name", None),
            legacy=getattr(args, "busco_library", None),
        )
        family_ids = []
        for tok in getattr(args, "family_id", []) or []:
            if tok is None:
                continue
            family_ids.append(str(tok))
        selectors = _selector_request_from_args(
            args,
            profile="busco_listing",
            busco_library_id=busco_library,
            manager=manager,
        )
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=True,
            require_candidates=False,
            use_rule_selection=False,
        )
        explicit_run_ids = expand_busco_run_id_variables(
            manager,
            getattr(args, "busco_run_ids", None)
            or ([getattr(args, "run_id", None)] if getattr(args, "run_id", None) is not None else []),
        )
        params: list[Any] = []
        where = ["1=1"]
        if explicit_run_ids:
            placeholders = ",".join("?" for _ in explicit_run_ids)
            where.append(f"d.run_id IN ({placeholders})")
            params.extend(explicit_run_ids)
        if busco_library is not None:
            where.append("d.library_id = ?")
            params.append(int(busco_library))
        if selected:
            placeholders = ",".join("?" for _ in selected)
            where.append(f"d.accession IN ({placeholders})")
            params.extend([str(a) for a in selected])
        if family_ids:
            placeholders = ",".join("?" for _ in family_ids)
            where.append(f"d.family_id IN ({placeholders})")
            params.extend(family_ids)
        manager.cursor.execute(
            f"""
            SELECT d.run_id, d.accession, d.family_id, d.status, d.sequence, d.score, d.length,
                   r.pipeline, r.input_mode
            FROM BUSCO_Run_Family_Data d
            JOIN BUSCO_Runs r ON r.run_id = d.run_id
            WHERE {' AND '.join(where)}
            ORDER BY d.run_id DESC, d.accession, d.family_id
            """,
            tuple(params),
        )
        rows = manager.cursor.fetchall() or []
        if not rows:
            print("No BUSCO family rows found.")
            return 0

        paralog_votes: Dict[tuple[int, str, str, int], str] = {}
        if manager._table_exists("Paralog_Filtering"):
            has_busco_link = manager._column_exists("Paralog_Filtering", "busco_run_id")
            run_ids = sorted({int(row[0]) for row in rows if row[0] is not None})
            accessions = sorted({str(row[1]) for row in rows if row[1]})
            family_tokens = sorted({str(row[2]) for row in rows if row[2]})
            library_ids = sorted({int(busco_library)}) if busco_library is not None else sorted({int(manager.busco.get_run(int(row[0]))[2]) for row in rows if row[0] is not None})
            if run_ids and accessions and family_tokens and library_ids:
                run_placeholders = ",".join("?" for _ in run_ids)
                acc_placeholders = ",".join("?" for _ in accessions)
                fam_placeholders = ",".join("?" for _ in family_tokens)
                lib_placeholders = ",".join("?" for _ in library_ids)
                select_busco_run = "pf.busco_run_id," if has_busco_link else "NULL AS busco_run_id,"
                manager.cursor.execute(
                    f"""
                    SELECT
                        pf.rowid,
                        pf.family_id,
                        pf.library_id,
                        pf.accession,
                        pf.clean,
                        pf.date,
                        {select_busco_run}
                        pf.run_id
                    FROM Paralog_Filtering pf
                    WHERE pf.accession IN ({acc_placeholders})
                      AND pf.family_id IN ({fam_placeholders})
                      AND pf.library_id IN ({lib_placeholders})
                    """,
                    tuple([*accessions, *family_tokens, *library_ids]),
                )
                candidate_rows = manager.cursor.fetchall() or []
                requested_keys = {
                    (int(run_id_v), str(acc), str(fam), int(busco_library if busco_library is not None else manager.busco.get_run(int(run_id_v))[2]))
                    for run_id_v, acc, fam, *_rest in rows
                }
                best_rank: Dict[tuple[int, str, str, int], tuple[str, int]] = {}
                for rowid_v, family_id_v, library_id_v, accession_v, clean_v, date_v, busco_run_id_v, _pf_run_id in candidate_rows:
                    family_key = str(family_id_v)
                    accession_key = str(accession_v)
                    try:
                        library_key = int(library_id_v)
                    except (TypeError, ValueError):
                        continue
                    decision = "CLEAN" if bool(clean_v) else "FAILED"
                    rank = (str(date_v or ""), int(rowid_v or 0))
                    if has_busco_link and busco_run_id_v is not None:
                        try:
                            key = (int(busco_run_id_v), accession_key, family_key, library_key)
                        except (TypeError, ValueError):
                            key = None
                        if key in requested_keys and rank >= best_rank.get(key, ("", -1)):
                            best_rank[key] = rank
                            paralog_votes[key] = decision
                        continue
                    for key in requested_keys:
                        run_key, acc_key, fam_key, lib_key = key
                        if accession_key != acc_key or family_key != fam_key or library_key != lib_key:
                            continue
                        if rank >= best_rank.get(key, ("", -1)):
                            best_rank[key] = rank
                            paralog_votes[key] = decision

        out = []
        status_map = {1: "Complete", 2: "Duplicated", 3: "Fragmented", 4: "Missing"}
        for row in rows:
            run_id_v, acc, fam, status, seq, score, length, pipeline, input_mode = row
            run_row = manager.busco.get_run(int(run_id_v))
            busco_lib_id = int(run_row[2]) if run_row and run_row[2] is not None else (int(busco_library) if busco_library is not None else -1)
            out.append(
                (
                    str(run_id_v),
                    str(acc),
                    str(fam),
                    status_map.get(int(status or 0), str(status or "")),
                    str(seq or ""),
                    "" if score is None else str(score),
                    "" if length is None else str(length),
                    str(pipeline or ""),
                    str(input_mode or ""),
                    paralog_votes.get((int(run_id_v), str(acc), str(fam), busco_lib_id), ""),
                )
            )
        headers = ("run_id", "accession", "family_id", "status", "sequence", "score", "length", "pipeline", "input_mode", "paralog_vote")
        return _render_list_output(args, headers, out, default_tidy=False)
    finally:
        manager.close()


def _handle_list_libraries(args: argparse.Namespace) -> int:
    """List libraries with optional library, parent, and reference filters."""

    manager = _connect_manager(args.database, read_only=True)
    try:
        join_params: list[Any] = []
        where_params: list[Any] = []
        where: list[str] = ["1=1"]
        if getattr(args, "library_id", None) is not None:
            where.append("l.library_id = ?")
            where_params.append(int(args.library_id))
        if getattr(args, "library_name", None):
            where.append("LOWER(l.library_name) = LOWER(?)")
            where_params.append(str(args.library_name))
        if getattr(args, "parent_id", None) is not None:
            where.append("l.parent_id = ?")
            where_params.append(int(args.parent_id))
        if getattr(args, "parent_name", None):
            where.append("LOWER(p.library_name) = LOWER(?)")
            where_params.append(str(args.parent_name))
        if getattr(args, "status", None):
            where.append("LOWER(COALESCE(l.status, 'ready')) = ?")
            where_params.append(str(args.status).strip().lower())

        ref_accessions = expand_accession_variables(
            manager,
            getattr(args, "ref_accessions", None) or [],
            allow_bare=True,
        )
        join_sql = ""
        having_sql = ""
        if ref_accessions:
            placeholders = ",".join("?" for _ in ref_accessions)
            join_sql = (
                "JOIN Reference_Assemblies raf ON raf.library_id = l.library_id "
                f"AND raf.accession IN ({placeholders})"
            )
            join_params.extend([str(acc) for acc in ref_accessions])
            having_sql = f"HAVING COUNT(DISTINCT raf.accession) = {len(set(ref_accessions))}"

        rows = manager.cursor.execute(
            f"""
            SELECT l.library_id,
                   l.library_name,
                   l.taxid,
                   p.library_name AS parent_name,
                   COALESCE(l.status, 'ready') AS status,
                   (
                       SELECT GROUP_CONCAT(ordered.accession)
                       FROM (
                           SELECT accession
                           FROM Reference_Assemblies
                           WHERE library_id = l.library_id
                           ORDER BY accession
                       ) AS ordered
                   ) AS accessions,
                   (
                       SELECT COUNT(*)
                       FROM Reference_Assemblies ra_count
                       WHERE ra_count.library_id = l.library_id
                   ) AS accession_count
            FROM Libraries l
            LEFT JOIN Libraries p ON p.library_id = l.parent_id
            {join_sql}
            WHERE {' AND '.join(where)}
            GROUP BY l.library_id, l.library_name, l.taxid, p.library_name, COALESCE(l.status, 'ready')
            {having_sql}
            ORDER BY l.library_id ASC
            """,
            tuple([*join_params, *where_params]),
        ).fetchall() or []
        if not rows:
            print("No libraries found.")
            return 0
        taxid_name_cache: dict[int, str] = {}

        def _coverage_name(taxid: Any) -> str:
            if taxid is None:
                return ""
            try:
                taxid_int = int(taxid)
            except (TypeError, ValueError):
                return ""
            if taxid_int in taxid_name_cache:
                return taxid_name_cache[taxid_int]
            rows_local = manager.get_lineage_root_to_leaf(taxid_int) or []
            name = ""
            if rows_local:
                leaf = rows_local[-1]
                if len(leaf) > 1 and leaf[1] is not None:
                    name = str(leaf[1])
            taxid_name_cache[taxid_int] = name
            return name

        extended = bool(getattr(args, "extended_library_info", False))
        out = []
        for library_id, library_name, taxid, parent_name, status, accessions, accession_count in rows:
            base = [
                str(library_id or ""),
                str(library_name or ""),
                _coverage_name(taxid),
                str(parent_name or ""),
                str(status or ""),
                str(accession_count or 0),
            ]
            if extended:
                base.extend(
                    [
                        "" if taxid is None else str(taxid),
                        str(accessions or ""),
                    ]
                )
            out.append(tuple(base))
        headers = [
            "library_id",
            "library_name",
            "coverage_clade",
            "parent_name",
            "status",
            "ref_count",
        ]
        if extended:
            headers.extend(["coverage_taxid", "ref_accessions"])
        sorted_out = _apply_row_sort_or_error(
            out,
            tuple(headers),
            args,
            aliases={
                "id": "library_id",
                "name": "library_name",
                "parent": "parent_name",
                "taxid": "coverage_taxid",
                "accessions": "ref_accessions",
                "size": "ref_count",
            },
        )
        if isinstance(sorted_out, int):
            return sorted_out
        out = sorted_out
        return _render_list_output(args, tuple(headers), out, default_tidy=True)
    finally:
        manager.close()


def _handle_list_proteome_profiles(args: argparse.Namespace) -> int:
    """List registered proteome profiles and preparation provenance."""

    manager = _connect_manager(args.database, read_only=True)
    try:
        selectors = _selector_request_from_args(args, profile="assembly_with_exclusions", manager=manager)
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=True,
            require_candidates=False,
            use_rule_selection=False,
        )
        rows = manager.proteomes.list_profiles(accessions=selected or None)
        if not rows:
            print("No proteome profiles found.")
            return 0
        out = []
        for row in rows:
            profile_id, accession, profile_name, kind, parent_profile_id, _artifact_id, status, sequence_count, checksum, is_default, created_at, _updated_at, _metadata_json = row
            prep = manager.proteomes.latest_preparation_for_output(int(profile_id))
            prep_type = prep[4] if prep else ""
            used_gff = "1" if prep and prep[5] else "0"
            cdhit_identity = "" if not prep or prep[10] is None else str(prep[10])
            out.append(
                (
                    str(profile_id),
                    str(accession),
                    str(profile_name),
                    str(kind),
                    "" if parent_profile_id is None else str(parent_profile_id),
                    "1" if is_default else "0",
                    str(status or ""),
                    "" if sequence_count is None else str(sequence_count),
                    prep_type or "",
                    used_gff,
                    cdhit_identity,
                    str(created_at or ""),
                    str(checksum or ""),
                )
            )
        headers = (
            "profile_id",
            "accession",
            "profile_name",
            "kind",
            "parent_profile_id",
            "is_default",
            "status",
            "sequence_count",
            "preparation_type",
            "used_gff",
            "cdhit_identity",
            "created_at",
            "checksum",
        )
        sorted_out = _apply_row_sort_or_error(
            out,
            headers,
            args,
            aliases={
                "id": "profile_id",
                "profile": "profile_name",
                "default": "is_default",
                "count": "sequence_count",
                "created": "created_at",
                "cdhit": "cdhit_identity",
            },
        )
        if isinstance(sorted_out, int):
            return sorted_out
        out = sorted_out
        return _render_list_output(args, headers, out, default_tidy=False)
    finally:
        manager.close()


def _handle_list(args: argparse.Namespace) -> int:
    """Dispatch ``list`` subcommands to the appropriate listing workflow."""

    if args.choice == "roots":
        return _handle_list_roots(args)
    if args.choice == "libraries":
        return _handle_list_libraries(args)
    if args.choice == "queue":
        return _handle_list_queue(args)
    if args.choice == "errors":
        return _handle_list_errors(args)
    if args.choice == "assemblies":
        return _handle_list_assemblies(args)
    if args.choice in {"busco", "results"}:
        args.busco = True
        args.has_busco_results = True
        return _handle_list_assemblies(args)
    if args.choice == "busco-runs":
        return _handle_list_busco_runs(args)
    if args.choice == "proteome-profiles":
        return _handle_list_proteome_profiles(args)
    if args.choice == "buscos":
        return _handle_list_buscos(args)
    if args.choice == "tasks":
        rows: List[Sequence[str]] = []
        for spec_entry in sorted(registry.specs(), key=lambda item: item.key):
            alias_suffix = f" (aliases: {', '.join(spec_entry.aliases)})" if spec_entry.aliases else ""
            rows.append(
                (
                    spec_entry.key,
                    str(spec_entry.job_type),
                    f"{spec_entry.description}{alias_suffix}",
                )
            )
        if not rows:
            print("No tasks registered.")
            return 0
        return _render_list_output(args, ("Task", "Job", "Description"), rows, default_tidy=True)

    if args.choice == "variables":
        manager = _connect_manager(args.database, read_only=True)
        try:
            records = manager.get_environment_variable_records()
        finally:
            manager.close()
        kind_filter_raw = str(getattr(args, "kind", "") or "").strip().lower()
        kind_filter = None
        if kind_filter_raw:
            kind_filter = normalize_variable_kind(kind_filter_raw)
            if kind_filter is None:
                return _print_error(f"--kind={kind_filter_raw} is invalid for list variables.")
        if getattr(args, "json", False):
            document = build_variables_json_document(records, source=args.database, kind_filter=kind_filter)
            print(json.dumps(document, indent=2, sort_keys=True))
            return 0
        rows = []
        for name, record in sorted(records.items()):
            value = record["value"]
            kind = normalize_variable_kind(record.get("kind")) or "env"
            if kind_filter and kind != kind_filter:
                continue
            display_name = f"@{name}" if kind == "assemblies" else name
            if kind_filter:
                rows.append((display_name, json.dumps(value)))
            else:
                rows.append((display_name, kind, json.dumps(value)))
        if not rows:
            print("No variables set in the database.")
            return 0
        headers = ("Variable", "Value") if kind_filter else ("Variable", "Kind", "Value")
        return _render_list_output(args, headers, rows, default_tidy=True)

    if args.choice == "ranks":
        manager = _connect_manager(args.database, read_only=True)
        try:
            manager.cursor.execute(
                "SELECT DISTINCT rank FROM Taxonomy WHERE rank IS NOT NULL AND TRIM(rank) != '' ORDER BY rank"
            )
            raw = [str(row[0]) for row in (manager.cursor.fetchall() or [])]
        finally:
            manager.close()
        if not raw:
            print("No taxonomic ranks found.")
            return 0
        hierarchy = {rank: idx for idx, rank in enumerate(RANK_HIERARCHY)}
        ordered = sorted(raw, key=lambda r: (hierarchy.get(r.lower(), 999), r.lower()))
        rows = [(rank,) for rank in ordered]
        return _render_list_output(args, ("Rank",), rows, default_tidy=True)

    if args.choice == "metadata":
        if not METADATA_FIELDS and not BUSCO_METADATA_TERMS:
            print("No metadata fields registered.")
            return 0
        return _print_metadata_terms()

    return _print_error(f"Unknown list choice '{args.choice}'.")


def _handle_list_queue(args: argparse.Namespace) -> int:
    options: Dict[str, Any] = {}
    queue_visibility_explicit = False
    if getattr(args, "status", None):
        options["status"] = args.status
    if getattr(args, "simple", False):
        options["simple"] = True
        queue_visibility_explicit = True
    if getattr(args, "hide_complete", False):
        options["hide_complete"] = True
        queue_visibility_explicit = True
    if getattr(args, "hide_done", False):
        options["hide_done"] = True
        queue_visibility_explicit = True
    if getattr(args, "sort", None):
        options["sort"] = args.sort
    if not queue_visibility_explicit and not bool(getattr(args, "show_all", False)):
        options["hide_done"] = True
    refresh = float(getattr(args, "refresh", 0) or 0)
    if refresh and not getattr(args, "watch", False):
        return _print_error("--refresh requires --watch.")
    if getattr(args, "watch", False) and refresh <= 0:
        refresh = 2.0

    columns = (
        TableColumn("Task ID", priority=100, min_width=7, max_width=18, no_wrap=True),
        TableColumn("Task Name", priority=90, min_width=10, max_width=32),
        TableColumn("Priority", priority=35, min_width=3, max_width=8, no_wrap=True),
        TableColumn("C", priority=30, min_width=1, max_width=3, no_wrap=True),
        TableColumn("Status", priority=95, min_width=6, max_width=8, no_wrap=True),
        TableColumn("Why", priority=45, min_width=3, max_width=10),
        TableColumn("Queue Time", priority=20, min_width=8, max_width=10, no_wrap=True),
        TableColumn("Start Time", priority=15, min_width=8, max_width=10, no_wrap=True),
        TableColumn("End Time", priority=10, min_width=8, max_width=10, no_wrap=True),
    )

    def provider() -> TableData:
        watched_manager = _connect_manager(args.database, read_only=True)
        try:
            provider_rows, provider_styles = collect_queue_rows(watched_manager, options)
            return TableData(QUEUE_HEADERS, provider_rows, provider_styles, columns)
        finally:
            watched_manager.close()

    manager = _connect_manager(args.database, read_only=True)
    try:
        color_defaults = _load_list_color_defaults(manager)
        rows, row_styles = collect_queue_rows(manager, options)
    finally:
        manager.close()
    return _render_list_output(
        args,
        QUEUE_HEADERS,
        rows,
        default_tidy=True,
        color_defaults=color_defaults,
        row_styles=row_styles,
        columns=columns,
        watch_provider=provider if getattr(args, "watch", False) else None,
        refresh=refresh or 2.0,
    )


def _handle_list_errors(args: argparse.Namespace) -> int:
    options: Dict[str, Any] = {}
    if getattr(args, "stack", False):
        options["include_stack"] = True
    if getattr(args, "limit", None) is not None:
        options["limit"] = args.limit
    refresh = float(getattr(args, "refresh", 0) or 0)
    if refresh and not getattr(args, "watch", False):
        return _print_error("--refresh requires --watch.")
    if getattr(args, "watch", False) and refresh <= 0:
        refresh = 2.0

    def provider() -> TableData:
        watched_manager = _connect_manager(args.database, read_only=True)
        try:
            provider_headers, provider_rows = collect_error_rows(watched_manager, options)
            return TableData(provider_headers, provider_rows, columns=_error_columns(provider_headers))
        finally:
            watched_manager.close()

    manager = _connect_manager(args.database, read_only=True)
    try:
        color_defaults = _load_list_color_defaults(manager)
        headers, rows = collect_error_rows(manager, options)
    finally:
        manager.close()
    return _render_list_output(
        args,
        headers,
        rows,
        default_tidy=True,
        color_defaults=color_defaults,
        columns=_error_columns(headers),
        watch_provider=provider if getattr(args, "watch", False) else None,
        refresh=refresh or 2.0,
    )


def _error_columns(headers: Sequence[str]) -> tuple[TableColumn, ...]:
    columns = [
        TableColumn("Task ID", priority=100, min_width=7, max_width=10, no_wrap=True, style="white"),
        TableColumn("Task Name", priority=90, min_width=10, max_width=30, style="white"),
        TableColumn("Error Message", priority=95, min_width=18, max_width=60, style="red"),
    ]
    if len(headers) > 3:
        columns.append(TableColumn("Error Stack", priority=10, min_width=20, max_width=80, style="yellow"))
    return tuple(columns)


def _handle_watch(args: argparse.Namespace) -> int:
    args.watch = True
    if args.watch_choice == "queue":
        return _handle_list_queue(args)
    if args.watch_choice == "errors":
        return _handle_list_errors(args)
    return _print_error(f"Unknown watch target '{args.watch_choice}'.")


# ---------------------------------------------------------------------------
# Assembly selector helpers and list/count workflows
# Purpose: Resolve assembly-focused selectors, enrich list output, and keep the
# higher-volume reporting commands grouped in one place.
# ---------------------------------------------------------------------------

def _rank_group_name(
    manager: DBManager,
    taxid: Optional[int],
    rank_token: str,
    cache: Dict[int, Optional[str]],
) -> Optional[str]:
    """Resolve the display name for a requested rank within a lineage."""

    if taxid is None:
        return None
    if taxid in cache:
        return cache[taxid]
    rows = manager.get_lineage_root_to_leaf(int(taxid)) or []
    if not rows:
        cache[taxid] = None
        return None

    order = {rank: idx for idx, rank in enumerate(RANK_HIERARCHY)}
    requested = (rank_token or "").lower()

    # Prefer an exact match when possible.
    exact = None
    for tid, name, rank, _parent in rows:
        if (rank or "").lower() == requested:
            exact = name or str(tid)
    if exact is not None:
        cache[taxid] = exact
        return exact

    # Fallback: pick the nearest available rank ABOVE (less specific) or AT the requested rank.
    if requested in order:
        req_idx = order[requested]
        best_idx = -1
        best_name = None
        for tid, name, rank, _parent in rows:
            token = (rank or "").lower()
            idx = order.get(token)
            if idx is None:
                continue
            if idx <= req_idx and idx > best_idx:
                best_idx = idx
                best_name = name or str(tid)
        cache[taxid] = best_name
        return best_name

    cache[taxid] = None
    return None


def _handle_list_assemblies(args: argparse.Namespace) -> int:
    """Render assembly listings with optional metadata, grouping, and BUSCO data."""

    manager = _connect_manager(args.database, read_only=not _list_requires_write(args))
    try:
        if getattr(args, "store_variable", None) and getattr(args, "append_to_variable", None):
            return _print_error("--store and --append-to cannot be used together.")
        _apply_busco_context_from_args(manager, args)
        one_line = bool(getattr(args, "one_line", False) or getattr(args, "one_line_quotes", False))
        if one_line:
            if not (args.accessions or args.clade or args.taxid is not None or getattr(args, "all", False) or getattr(args, "preset_name", None)):
                return _print_error("--one-line/--one-line-quotes requires --accessions, --clade, --taxid, --preset, or --all.")
            # --one-line overrides incompatible output flags.
            args.group_by_rank = False
            args.busco = False
            args.tidy = False
            args.output_path = None
            if hasattr(args, "extended_decontamination_headers"):
                args.extended_decontamination_headers = False

        try:
            busco_library = _resolve_library_selector(
                manager,
                library_id=args.library_id,
                library_name=args.library_name,
                legacy=args.busco_library,
            )
        except ValueError as exc:
            return _print_error(str(exc))
        try:
            explicit_run_ids = expand_busco_run_id_variables(manager, getattr(args, "busco_run_ids", None) or [])
        except ValueError as exc:
            return _print_error(str(exc))
        if args.busco and busco_library is None and explicit_run_ids:
            try:
                busco_library = _infer_busco_library_from_run_ids(manager, explicit_run_ids)
            except ValueError as exc:
                return _print_error(str(exc))

        if args.busco and busco_library is None:
            return _print_error("Use --busco with --library-id or --library-name.")
        if getattr(args, "all_runs", False) and not args.busco:
            return _print_error("--all-runs requires --busco.")
        if getattr(args, "extended_decontamination_headers", False) and not args.busco:
            return _print_error("--extended-decontamination-headers requires --busco.")

        if getattr(args, "group_by_rank", False) and not (args.rank or getattr(args, "ranks", None)):
            return _print_error("--group-by-rank requires --rank or --ranks.")

        selectors = _selector_request_from_args(
            args,
            profile="assembly_with_exclusions",
            busco_library_id=busco_library,
            manager=manager,
        )

        try:
            selected = resolve_selector_accessions(
                manager,
                selectors,
                allow_all=True,
                require_candidates=False,
                use_rule_selection=True,
            )
        except ValueError as exc:
            return _print_error(str(exc))

        try:
            selected = _apply_accession_set_ops(manager, args, selected)
        except ValueError as exc:
            return _print_error(str(exc))

        try:
            if args.store_variable:
                _store_accession_variable(manager, args.store_variable, selected)
            elif getattr(args, "append_to_variable", None):
                _store_accession_variable(manager, args.append_to_variable, selected, append=True)
        except ValueError as exc:
            return _print_error(str(exc))

        if getattr(args, "one_line_quotes", False):
            print(",".join(f"\"{acc}\"" for acc in selected))
            return 0
        if getattr(args, "one_line", False):
            print(",".join(selected))
            return 0

        if not selected:
            if args.output_path:
                out_path = Path(args.output_path)
                out_path.parent.mkdir(parents=True, exist_ok=True)
                out_path.write_text("", encoding="utf-8")
            else:
                print("No assemblies matched selectors.", file=sys.stderr)
            return 0

        info_map = _fetch_assembly_info(manager, selected)
        meta_fields: List[str] = []
        if getattr(args, "meta", None) is not None:
            if args.meta == "__DEFAULT__":
                meta_fields = ["release_date", "level", "n50", "comments"]
            else:
                try:
                    meta_fields = _resolve_meta_fields([args.meta])
                except ValueError as exc:
                    return _print_error(str(exc))
        meta_values = _fetch_metadata(manager, selected, meta_fields) if meta_fields else {}
        parent_library_id = manager.assert_library_has_parent(busco_library) if busco_library else None
        custom_library = bool(parent_library_id)
        show_hidden = bool(busco_library)
        show_contam = True
        busco_run_map = (
            _fetch_busco_runs(
                manager,
                selected,
                int(parent_library_id) if custom_library and parent_library_id is not None else busco_library,
                purpose=(
                    "export_nucleotide"
                    if getattr(args, "busco_export_format", None) == "nucleotide"
                    else "export_protein"
                    if getattr(args, "busco_export_format", None) == "protein"
                    else "default"
                ),
                run_ids=explicit_run_ids or None,
                pipeline=getattr(args, "busco_pipeline", None),
                input_mode=getattr(args, "busco_input_mode", None) or getattr(args, "format", None),
                preferred_pipeline=getattr(args, "prefer_busco_pipeline", None),
                preferred_input_mode=getattr(args, "prefer_busco_input_mode", None) or getattr(args, "prefer_format", None),
                proteome_profile=getattr(args, "proteome_profile", None),
                preferred_proteome_profile=getattr(args, "prefer_proteome_profile", None),
                export_format=getattr(args, "busco_export_format", None),
            )
            if args.busco
            else {}
        )
        protein_flag_map = (
            _fetch_genome_protein_flags(manager, selected)
            if args.busco
            else {}
        )
        if args.busco:
            score_accessions = list(selected)
            if getattr(args, "all_runs", False):
                score_accessions = [acc for acc in selected if not busco_run_map.get(acc)]
            busco_map = (
                _fetch_busco_scores(
                    manager,
                    score_accessions,
                    busco_library,
                    include_paralog_filtering_in_score=args.include_paralog_filtering_in_score,
                    include_decontamination_in_score=args.include_decontamination_in_score,
                    paralog_run_id=args.use_paralog_run,
                    decontamination_run_id=args.use_decontamination_run,
                    allow_ambiguous_contaminants=args.allow_ambiguous_contaminants,
                    strict_decontamination=args.strict_decontamination,
                    rescue_duplicates=getattr(args, "rescue_duplicates", False),
                )
                if score_accessions
                else {}
            )
        else:
            busco_map = {}
        busco_run_display_map: Dict[Tuple[str, int], Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], bool, bool, Optional[float], Optional[float], Optional[float], Optional[str], Optional[str], Optional[str]]] = {}
        if args.busco and busco_run_map:
            run_refs = []
            for acc in selected:
                run_items = busco_run_map.get(acc) or []
                if getattr(args, "all_runs", False):
                    run_refs.extend(
                        (acc, int(run_item["run_id"]))
                        for run_item in run_items
                        if run_item.get("run_id") is not None
                    )
                elif run_items and run_items[0].get("run_id") is not None:
                    run_refs.append((acc, int(run_items[0]["run_id"])))
            if run_refs:
                try:
                    busco_run_display_map = manager.busco.get_display_results_for_runs(
                        library_id=int(busco_library),
                        run_refs=run_refs,
                        include_paralog=getattr(args, "include_paralog_filtering_in_score", None),
                        paralog_run_id=getattr(args, "use_paralog_run", None),
                        include_decontam=getattr(args, "include_decontamination_in_score", None),
                        decont_run_id=getattr(args, "use_decontamination_run", None),
                        allow_ambiguous_contaminants=getattr(args, "allow_ambiguous_contaminants", None),
                        strict_decontamination=getattr(args, "strict_decontamination", None),
                        rescue_duplicates=getattr(args, "rescue_duplicates", False),
                    )
                except Exception as exc:  # boundary: run-specific display enrichment failure becomes a CLI validation error.
                    return _print_error(f"Failed to compute run-specific BUSCO display values: {exc}")
        decont_summary = (
            _fetch_decontamination_summaries(
                manager,
                selected,
                busco_library,
                run_id=args.use_decontamination_run or args.decontamination_run,
            )
            if args.busco
            else {}
        )
        decont_breakdown = (
            _fetch_decontamination_breakdown(
                manager,
                selected,
                busco_library,
                run_id=args.use_decontamination_run or args.decontamination_run,
            )
            if args.busco and getattr(args, "extended_decontamination_headers", False)
            else {}
        )
        rank_labels: List[str] = []
        if getattr(args, "ranks", None):
            rank_labels = normalize_rank_list(args.ranks)
        elif args.rank:
            rank_labels = [str(args.rank).strip()]
        rank_tokens = [label.lower() for label in rank_labels]
        rank_cache: Dict[str, Dict[int, Optional[str]]] = {token: {} for token in rank_tokens}
        rank_value_by_acc: Dict[str, List[str]] = {}
        if rank_labels:
            for acc in selected:
                taxid_val = info_map.get(acc, (None, ""))[0]
                values: List[str] = []
                for token in rank_tokens:
                    cache = rank_cache.get(token, {})
                    values.append(_rank_group_name(manager, taxid_val, token, cache) or "")
                rank_value_by_acc[acc] = values
        profile_display_info = _fetch_proteome_profile_display_info(manager, selected) if (meta_fields and ("proteome_profile" in meta_fields or "default_proteome_profile" in meta_fields)) or args.busco else {}
        proteome_profile_meta_index: Optional[int] = None
        if meta_fields and "proteome_profile" in meta_fields:
            proteome_profile_meta_index = 1 + (2 if args.busco else 0) + 1 + len(rank_labels) + meta_fields.index("proteome_profile")
        orthofinder_target_meta_index: Optional[int] = None
        orthofinder_target_libraries: Dict[int, str] = {}
        if args.busco and meta_fields and "orthofinder_target_library" in meta_fields:
            orthofinder_target_meta_index = 1 + 2 + 1 + len(rank_labels) + meta_fields.index("orthofinder_target_library")
            displayed_orthofinder_run_ids: List[int] = []
            for acc in selected:
                run_items = busco_run_map.get(acc) or []
                if not getattr(args, "all_runs", False):
                    run_items = run_items[:1]
                displayed_orthofinder_run_ids.extend(
                    int(run_item["run_id"])
                    for run_item in run_items
                    if run_item.get("run_id") is not None
                    and str(run_item.get("pipeline_raw") or "").strip().lower() == "orthofinder"
                )
            if displayed_orthofinder_run_ids:
                orthofinder_target_libraries = _fetch_orthofinder_target_libraries(
                    manager,
                    displayed_orthofinder_run_ids,
                )

        row_map: Dict[str, List[List[str]]] = {}
        for acc in selected:
            tax_info = info_map.get(acc, (None, ""))
            base = [acc]
            if args.busco:
                base.extend(["", ""])
            base.append(tax_info[1])
            if rank_labels:
                base.extend(rank_value_by_acc.get(acc, [""] * len(rank_labels)))
            if meta_fields:
                values = meta_values.get(acc, {})
                for field in meta_fields:
                    label = METADATA_FIELDS[field]["label"]
                    base.append(values.get(label, ""))

            busco_rows_for_acc = busco_run_map.get(acc) if getattr(args, "all_runs", False) else None
            if busco_rows_for_acc:
                expanded_rows = []
                for run_item in busco_rows_for_acc:
                    row = list(base)
                    if args.busco:
                        if proteome_profile_meta_index is not None:
                            input_mode_norm = str(run_item.get("input_mode") or "").strip().lower()
                            run_profile = str(run_item.get("proteome_profile") or "")
                            details = profile_display_info.get(acc, {}).get(run_profile)
                            row[proteome_profile_meta_index] = (
                                str(details.get("display") or "")
                                if details
                                else run_profile or ("unset" if input_mode_norm == "protein" else "")
                            )
                        if orthofinder_target_meta_index is not None:
                            run_id = run_item.get("run_id")
                            row[orthofinder_target_meta_index] = (
                                orthofinder_target_libraries.get(int(run_id), "")
                                if run_id is not None
                                else ""
                            )
                        row[1] = _busco_mode_code(
                            run_item.get("pipeline_raw"),
                            run_item.get("input_mode"),
                            is_primary=bool(run_item.get("is_primary")),
                        )
                        row[2] = _busco_format_code(
                            run_item.get("input_mode"),
                            pipeline=run_item.get("pipeline_raw"),
                            assembly_has_protein=protein_flag_map.get(acc),
                        )
                    run_id = run_item.get("run_id")
                    display_entry = busco_run_display_map.get((acc, int(run_id))) if run_id is not None else None
                    bvals = display_entry[:7] if display_entry is not None else None
                    if bvals is None:
                        bvals = (
                            run_item.get("complete"),
                            run_item.get("single_copy"),
                            run_item.get("duplicated"),
                            run_item.get("fragmented"),
                            run_item.get("missing"),
                            None,
                            None,
                        )
                    if show_hidden:
                        complete, single_copy, duplicated, fragmented, missing, hidden_paralog, contaminated = bvals
                        row.extend(
                            _fmt_busco_value(val)
                            for val in (
                                complete,
                                single_copy,
                                duplicated,
                                fragmented,
                                hidden_paralog,
                                contaminated,
                            )
                        )
                    else:
                        complete, single_copy, duplicated, fragmented, missing, _hidden, contaminated = bvals
                        row.extend(
                            _fmt_busco_value(val)
                            for val in (
                                complete,
                                single_copy,
                                duplicated,
                                fragmented,
                                contaminated,
                            )
                        )
                    if getattr(args, "extended_decontamination_headers", False):
                        if display_entry is not None:
                            support, weak, unknown = display_entry[9], display_entry[10], display_entry[11]
                            paralog_run = display_entry[12]
                            decont_display_run = display_entry[13]
                        else:
                            support, weak, unknown = decont_breakdown.get(acc, (None, None, None))
                            paralog_run = None
                            decont_display_run = None
                        row.extend(_fmt_busco_value(val) for val in (support, weak, unknown))
                        row.append(str(paralog_run) if paralog_run else "NA")
                        row.append(str(decont_display_run) if decont_display_run else "NA")
                    row.append(_fmt_busco_value(missing))
                    if custom_library:
                        decision = display_entry[14] if display_entry is not None else decont_summary.get(acc, (None, None, None))[1]
                        row.append(str(decision) if decision else "NA")
                    expanded_rows.append(row)
                row_map[acc] = expanded_rows
                continue

            row = list(base)
            if args.busco:
                primary_run = (busco_run_map.get(acc) or [None])[0]
                if primary_run:
                    if proteome_profile_meta_index is not None:
                        input_mode_norm = str(primary_run.get("input_mode") or "").strip().lower()
                        run_profile = str(primary_run.get("proteome_profile") or "")
                        details = profile_display_info.get(acc, {}).get(run_profile)
                        row[proteome_profile_meta_index] = (
                            str(details.get("display") or "")
                            if details
                            else run_profile or ("unset" if input_mode_norm == "protein" else "")
                        )
                    if orthofinder_target_meta_index is not None:
                        run_id = primary_run.get("run_id")
                        row[orthofinder_target_meta_index] = (
                            orthofinder_target_libraries.get(int(run_id), "")
                            if run_id is not None
                            else ""
                        )
                    row[1] = _busco_mode_code(
                        primary_run.get("pipeline_raw"),
                        primary_run.get("input_mode"),
                        is_primary=bool(primary_run.get("is_primary")),
                    )
                    row[2] = _busco_format_code(
                        primary_run.get("input_mode"),
                        pipeline=primary_run.get("pipeline_raw"),
                        assembly_has_protein=protein_flag_map.get(acc),
                    )
                display_entry = None
                if primary_run and primary_run.get("run_id") is not None:
                    display_entry = busco_run_display_map.get((acc, int(primary_run["run_id"])))
                bvals = display_entry[:7] if display_entry is not None else busco_map.get(acc)
                if bvals:
                    complete, single_copy, duplicated, fragmented, missing, hidden_paralog, contaminated = bvals
                    if show_hidden:
                        row.extend(
                            _fmt_busco_value(val)
                            for val in (
                                complete,
                                single_copy,
                                duplicated,
                                fragmented,
                                hidden_paralog,
                                contaminated,
                            )
                        )
                    else:
                        row.extend(
                            _fmt_busco_value(val)
                            for val in (
                                complete,
                                single_copy,
                                duplicated,
                                fragmented,
                                contaminated,
                            )
                        )
                    if getattr(args, "extended_decontamination_headers", False):
                        if display_entry is not None:
                            support, weak, unknown = display_entry[9], display_entry[10], display_entry[11]
                            paralog_run = display_entry[12]
                            decont_display_run = display_entry[13]
                        else:
                            support, weak, unknown = decont_breakdown.get(acc, (None, None, None))
                            paralog_run = None
                            decont_display_run = decont_summary.get(acc, (None, None, None))[0]
                        row.extend(_fmt_busco_value(val) for val in (support, weak, unknown))
                        row.append(str(paralog_run) if paralog_run else "NA")
                        row.append(str(decont_display_run) if decont_display_run else "NA")
                    row.append(_fmt_busco_value(missing))
                    if custom_library:
                        decision = display_entry[14] if display_entry is not None else decont_summary.get(acc, (None, None, None))[1]
                        row.append(str(decision) if decision else "NA")
                else:
                    if show_hidden:
                        row.extend(["", "", "", "", "", ""])
                    else:
                        row.extend(["", "", "", "", ""])
                    if getattr(args, "extended_decontamination_headers", False):
                        row.extend(["", "", "", "", ""])
                    row.append("")
                    if custom_library:
                        row.append("")
            row_map[acc] = [row]

        headers = ["accession"]
        if args.busco:
            headers.extend(["M", "F"])
        headers.append("species")
        if rank_labels:
            headers.extend(rank_labels)
        if meta_fields:
            headers.extend(METADATA_FIELDS[field]["label"] for field in meta_fields)
        if args.busco:
            if show_hidden:
                headers.extend(
                    [
                        "complete",
                        "single_copy",
                        "duplicated",
                        "fragmented",
                        "hidden_paralog",
                        "contaminated",
                    ]
                )
            else:
                headers.extend(
                    [
                        "complete",
                        "single_copy",
                        "duplicated",
                        "fragmented",
                        "contaminated",
                    ]
                )
            if getattr(args, "extended_decontamination_headers", False):
                headers.extend(["support", "weak", "unknown", "paralog_run", "decontamination_run"])
            headers.append("missing")
            if custom_library:
                headers.append("contaminated_assembly")

        groups: List[Tuple[Optional[str], Sequence[Sequence[str]]]] = []
        group_rank_label = str(args.rank).strip() if args.rank else (rank_labels[0] if rank_labels else None)
        group_rank_index = 0
        if group_rank_label and rank_labels:
            try:
                group_rank_index = rank_labels.index(group_rank_label)
            except ValueError:
                group_rank_index = 0
        if getattr(args, "group_by_rank", False) and group_rank_label and rank_labels:
            if len(rank_labels) == 1:
                grouped: Dict[str, List[Sequence[str]]] = {}
                order: List[str] = []
                for acc in selected:
                    values = rank_value_by_acc.get(acc, [])
                    group_name = (values[group_rank_index] if values else "") or "unranked"
                    if group_name not in grouped:
                        grouped[group_name] = []
                        order.append(group_name)
                    grouped[group_name].extend(row_map[acc])
                for group_name in order:
                    count = len(grouped[group_name])
                    groups.append((f"{group_name} [{count}]", grouped[group_name]))
            else:
                def _group_hierarchy(accessions: Sequence[str], depth: int) -> None:
                    grouped_accs: Dict[str, List[str]] = {}
                    order: List[str] = []
                    for acc in accessions:
                        values = rank_value_by_acc.get(acc, [])
                        group_name = (values[depth] if len(values) > depth else "") or "unranked"
                        if group_name not in grouped_accs:
                            grouped_accs[group_name] = []
                            order.append(group_name)
                        grouped_accs[group_name].append(acc)

                    for group_name in order:
                        group_accs = grouped_accs[group_name]
                        count = len(group_accs)
                        header = f"{group_name} [{count}]"
                        indent = "  " * depth
                        if depth >= len(rank_labels) - 1:
                            row_indent = "  " * (depth + 1)
                            indented_rows = []
                            for a in group_accs:
                                for row in row_map[a]:
                                    if row:
                                        indented_rows.append([row_indent + row[0]] + list(row[1:]))
                                    else:
                                        indented_rows.append(row)
                            groups.append((indent + header, indented_rows))
                        else:
                            groups.append((indent + header, []))
                            _group_hierarchy(group_accs, depth + 1)

                _group_hierarchy(selected, 0)
        else:
            flat_rows: List[List[str]] = []
            for acc in selected:
                flat_rows.extend(row_map[acc])
            groups.append((None, flat_rows))

        if getattr(args, "sort", None):
            try:
                groups = [
                    (group_name, _sort_rows(group_rows, headers, args.sort, ASSEMBLY_SORT_ALIASES))
                    for group_name, group_rows in groups
                ]
            except ValueError as exc:
                return _print_error(str(exc))

        color_defaults = _load_list_color_defaults(manager)
        list_color = args.list_color if args.list_color is not None else _coerce_bool(color_defaults.get("LIST_USE_COLOR", False))
        use_rich = bool(list_color) and args.output_path is None and sys.stdout.isatty()

        if use_rich:
            group_colors = _parse_color_tokens(
                color_defaults.get("LIST_GROUP_COLORS"),
                ["bright_cyan", "bright_yellow", "cyan", "blue", "magenta"],
            )
            base_grad = _parse_color_tokens(
                color_defaults.get("LIST_BUSCO_GRADIENT"),
                ["#d73027", "#fee08b", "#1a9850"],
            )
            pos_grad = _parse_color_tokens(
                color_defaults.get("LIST_BUSCO_GRADIENT_POS"),
                base_grad,
            )
            neg_grad = _parse_color_tokens(
                color_defaults.get("LIST_BUSCO_GRADIENT_NEG"),
                list(reversed(base_grad)),
            )
            pos_stops = _parse_float_tokens(
                color_defaults.get("LIST_BUSCO_POS_STOPS"),
                [30.0, 65.0, 100.0],
            )
            neg_stops = _parse_float_tokens(
                color_defaults.get("LIST_BUSCO_NEG_STOPS"),
                [0.0, 35.0, 70.0],
            )
            steep_stops = _parse_float_tokens(
                color_defaults.get("LIST_BUSCO_STEEP_STOPS"),
                [0.0, 10.0, 20.0],
            )

            busco_pos_keys = {"complete", "single_copy"}
            busco_neg_keys = {"duplicated", "fragmented", "missing"}
            busco_steep_keys = {"hidden_paralog", "contaminated"}
            busco_bold_keys = {"complete", "single_copy", "duplicated"}
            busco_pos_indices = [i for i, h in enumerate(headers) if h in busco_pos_keys]
            busco_neg_indices = [i for i, h in enumerate(headers) if h in busco_neg_keys]
            busco_steep_indices = [i for i, h in enumerate(headers) if h in busco_steep_keys]
            busco_bold_indices = [i for i, h in enumerate(headers) if h in busco_bold_keys]
            busco_steep_max = _coerce_float(color_defaults.get("LIST_BUSCO_STEEP_MAX", 20.0), 20.0)
            rank_start = headers.index(rank_labels[0]) if rank_labels else 1
            rank_color_map = {rank_start + i: i for i in range(len(rank_labels))}

            _render_grouped_rows_rich(
                groups,
                headers=headers,
                group_colors=group_colors,
                busco_pos_gradient=pos_grad,
                busco_neg_gradient=neg_grad,
                busco_pos_stops=pos_stops,
                busco_neg_stops=neg_stops,
                busco_steep_stops=steep_stops,
                busco_pos_indices=busco_pos_indices,
                busco_neg_indices=busco_neg_indices,
                busco_steep_indices=busco_steep_indices,
                busco_steep_max=busco_steep_max,
                rank_color_map=rank_color_map,
                busco_bold_indices=busco_bold_indices,
                paginate=not bool(getattr(args, "no_pager", False)),
                show_header=not bool(getattr(args, "no_header", False)),
            )
        else:
            output = _render_grouped_rows(
                groups,
                headers=headers,
                tidy=bool(args.tidy),
                show_header=not bool(getattr(args, "no_header", False)),
            )
            if args.output_path:
                out_path = Path(args.output_path)
                out_path.parent.mkdir(parents=True, exist_ok=True)
                suffix = "" if output.endswith("\n") or not output else "\n"
                out_path.write_text(output + suffix, encoding="utf-8")
            else:
                print(output)
        return 0
    finally:
        manager.close()


def _handle_count(args: argparse.Namespace) -> int:
    """Count assemblies after applying the same selector logic used by listing."""

    if args.subject != "assemblies":
        return _print_error(f"Unknown count subject '{args.subject}'.")

    manager = _connect_manager(args.database, read_only=True)
    try:
        _apply_busco_context_from_args(manager, args)
        if getattr(args, "group_by_rank", False) and not (args.rank or getattr(args, "ranks", None)):
            return _print_error("--group-by-rank requires --rank or --ranks.")

        busco_library = _resolve_library_selector(
            manager,
            library_id=args.library_id,
            library_name=args.library_name,
            legacy=args.busco_library,
        )

        selectors = _selector_request_from_args(
            args,
            profile="assembly_with_exclusions",
            busco_library_id=busco_library,
            manager=manager,
        )

        try:
            candidates = resolve_selector_accessions(
                manager,
                selectors,
                allow_all=True,
                require_candidates=False,
                use_rule_selection=True,
            )
        except ValueError as exc:
            return _print_error(str(exc))

        if getattr(args, "group_by_rank", False) and (args.rank or getattr(args, "ranks", None)):
            info_map = _fetch_assembly_info(manager, candidates)
            rank_labels = normalize_rank_list(getattr(args, "ranks", None) or [])
            rank_label = str(args.rank).strip() if args.rank else (rank_labels[0] if rank_labels else "rank")
            rank_token = rank_label.lower()
            cache: Dict[int, Optional[str]] = {}
            grouped: Dict[str, int] = {}
            order: List[str] = []
            for acc in candidates:
                taxid_val = info_map.get(acc, (None, ""))[0]
                group_name = _rank_group_name(manager, taxid_val, rank_token, cache) or "unranked"
                if group_name not in grouped:
                    grouped[group_name] = 0
                    order.append(group_name)
                grouped[group_name] += 1
            for group_name in order:
                print(f"{rank_label}: {group_name} [{grouped[group_name]}]")
            return 0

        print(str(len(candidates)))
        return 0
    finally:
        manager.close()


def _handle_assemblies(args: argparse.Namespace) -> int:
    """Print resolved assembly accessions without the richer list formatting layer."""

    manager = _connect_manager(args.database, read_only=not _list_requires_write(args))
    try:
        if getattr(args, "store_variable", None) and getattr(args, "append_to_variable", None):
            return _print_error("--store and --append-to cannot be used together.")
        if not (args.accessions or args.taxid is not None or args.clade or getattr(args, "preset_name", None)):
            return _print_error("Provide --accessions, --taxid, --clade, or --preset to select assemblies.")

        busco_library = _resolve_library_selector(
            manager,
            library_id=args.library_id,
            library_name=args.library_name,
            legacy=args.busco_library,
        )

        selectors = _selector_request_from_args(
            args,
            profile="assembly",
            busco_library_id=busco_library,
            manager=manager,
        )

        try:
            selected = resolve_selector_accessions(
                manager,
                selectors,
                allow_all=False,
                require_candidates=True,
                use_rule_selection=True,
            )
        except ValueError as exc:
            return _print_error(str(exc))

        try:
            selected = _apply_accession_set_ops(manager, args, selected)
        except ValueError as exc:
            return _print_error(str(exc))

        try:
            if args.store_variable:
                _store_accession_variable(manager, args.store_variable, selected)
            elif getattr(args, "append_to_variable", None):
                _store_accession_variable(manager, args.append_to_variable, selected, append=True)
        except ValueError as exc:
            return _print_error(str(exc))

        for acc in selected:
            print(acc)
        return 0
    finally:
        manager.close()



def _positive_refresh(value: str) -> float:
    refresh = float(value)
    if refresh <= 0:
        raise argparse.ArgumentTypeError("refresh interval must be positive")
    return refresh


def _add_queue_monitor_options(group: argparse._ArgumentGroup, *, include_watch: bool) -> None:
    if include_watch:
        group.add_argument("-w", "--watch", action="store_true", help="Refresh the queue continuously until interrupted.")
    group.add_argument("--refresh", type=_positive_refresh, help="Refresh interval in seconds for --watch.")
    group.add_argument("--all", dest="show_all", action="store_true", help="Show completed and errored tasks too.")
    group.add_argument("--status", action="append", help="Filter to specific task statuses (repeatable).")
    group.add_argument("--simple", action="store_true", help="Hide completed tasks once their parent chain is completed.")
    group.add_argument("--hide-complete", action="store_true", help="Hide all completed tasks.")
    group.add_argument("--hide-done", action="store_true", help="Hide completed and errored tasks.")
    group.add_argument(
        "-s",
        "--sort",
        choices=("latest", "changed", "new", "newest", "old", "oldest", "errors", "running", "active", "status"),
        default="latest",
        help=(
            "Queue sort profile: latest/changed = newest status change anywhere in each parent/subtask tree; "
            "new/newest = newest task id; old/oldest = oldest task id; errors = blocks with errors first; "
            "running/active = blocks with active tasks first; status = group by task status. "
            "Children are sorted recursively inside each block."
        ),
    )


def _add_error_monitor_options(group: argparse._ArgumentGroup, *, include_watch: bool) -> None:
    if include_watch:
        group.add_argument("-w", "--watch", action="store_true", help="Refresh the errors list continuously until interrupted.")
    group.add_argument("--refresh", type=_positive_refresh, help="Refresh interval in seconds for --watch.")
    group.add_argument("--stack", action="store_true", help="Include stored traceback text.")
    group.add_argument("--limit", type=int, help="Limit the number of error rows.")


def register_list_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
    handler=None,
) -> argparse.ArgumentParser:
    """Register the ``list`` command and its subject-specific subparsers."""

    list_handler = handler or _handle_list
    list_parser = subparsers.add_parser("list", help="List tasks, queue/errors, assemblies, libraries, roots, or variables.")
    list_subparsers = list_parser.add_subparsers(dest="choice", required=True)

    list_tasks = list_subparsers.add_parser("tasks", help="List registered tasks.")
    _add_basic_list_output_options(list_tasks.add_argument_group("Output options"))
    list_tasks.set_defaults(handler=list_handler)

    list_queue = list_subparsers.add_parser("queue", help="List queued, running, blocked, and completed tasks.")
    list_queue_output_group = list_queue.add_argument_group("Output options")
    _add_queue_monitor_options(list_queue_output_group, include_watch=True)
    _add_basic_list_output_options(list_queue_output_group)
    list_queue.set_defaults(handler=list_handler)

    list_errors = list_subparsers.add_parser("errors", help="List task failures.")
    list_errors_output_group = list_errors.add_argument_group("Output options")
    _add_error_monitor_options(list_errors_output_group, include_watch=True)
    _add_basic_list_output_options(list_errors_output_group)
    list_errors.set_defaults(handler=list_handler)

    list_variables = list_subparsers.add_parser("variables", help="List stored environment, accession, and BUSCO run-id variables.")
    list_variables.add_argument("--kind", choices=["assemblies", "assembly", "accessions", "accession", "busco-runs", "busco-run", "runs", "env", "environment", "environmental", "enviornmental"], help="Filter variables by kind.")
    list_variables.add_argument("--json", action="store_true", help="Emit variables as kinded JSON.")
    _add_basic_list_output_options(list_variables.add_argument_group("Output options"))
    list_variables.set_defaults(handler=list_handler)

    list_ranks = list_subparsers.add_parser("ranks", help="List known taxonomy ranks.")
    _add_basic_list_output_options(list_ranks.add_argument_group("Output options"))
    list_ranks.set_defaults(handler=list_handler)

    list_metadata = list_subparsers.add_parser("metadata", help="List supported metadata field keys.")
    _add_basic_list_output_options(list_metadata.add_argument_group("Output options"))
    list_metadata.set_defaults(handler=list_handler)

    list_roots = list_subparsers.add_parser("roots", help="List registered storage roots.")
    list_roots.add_argument("--kind", choices=STORAGE_ROOT_KINDS, help="Filter listed roots by logical kind.")
    _add_basic_list_output_options(list_roots.add_argument_group("Output options"))
    list_roots.set_defaults(handler=list_handler)

    list_libraries = list_subparsers.add_parser("libraries", help="List libraries and their reference accessions.")
    list_libraries_selector_group = list_libraries.add_argument_group("Selector options")
    list_libraries_output_group = list_libraries.add_argument_group("Output options")
    list_libraries_selector_group.add_argument("-li", "--library-id", type=int, help="Limit output to one library id.")
    list_libraries_selector_group.add_argument("-l", "--library-name", help="Limit output to one library name.")
    list_libraries_selector_group.add_argument("--parent-id", type=int, help="Limit output to libraries with the given parent library id.")
    list_libraries_selector_group.add_argument("--parent-name", help="Limit output to libraries with the given parent library name.")
    list_libraries_selector_group.add_argument("--status", choices=["ready", "stale"], help="Limit output to libraries with the given status.")
    list_libraries_selector_group.add_argument("--ref-accessions", action=AppendCommaSeparated, help="Limit output to libraries containing all specified reference accessions.")
    list_libraries_output_group.add_argument("--extended-library-info", action="store_true", help="Include raw coverage taxid and full reference accession list.")
    _add_basic_list_output_options(list_libraries_output_group)
    _add_row_sort_option(
        list_libraries_output_group,
        fields="library_id, library_name, coverage_clade, parent_name, status, ref_count",
    )
    list_libraries.set_defaults(handler=list_handler)

    list_assemblies = list_subparsers.add_parser("assemblies", help="List assemblies matching selector filters.")
    list_assemblies_selector_group = list_assemblies.add_argument_group("Selector options")
    list_assemblies_output_group = list_assemblies.add_argument_group("Output options")
    _add_selector_arguments(
        list_assemblies_selector_group,
        profile="assembly_with_exclusions",
        selector_defaults=selector_defaults,
        context_label="list assemblies",
    )
    list_assemblies_selector_group.add_argument(
        "-m",
        "--meta",
        "--metadata",
        dest="meta",
        nargs="?",
        const="__DEFAULT__",
        help="Include metadata columns. Omit value for defaults (release_date, level, n50, comments).",
    )
    _add_assembly_list_output_options(list_assemblies_output_group)
    _add_row_sort_option(
        list_assemblies_output_group,
        fields="accession, species/taxon, rank columns, metadata columns, latest, quality, BUSCO columns when shown",
    )
    list_assemblies.set_defaults(handler=list_handler)

    list_busco = list_subparsers.add_parser(
        "busco",
        aliases=["results"],
        help="Alias of 'list assemblies --busco --has-busco-results'.",
    )
    list_busco_selector_group = list_busco.add_argument_group("Selector options")
    list_busco_output_group = list_busco.add_argument_group("Output options")
    _add_selector_arguments(
        list_busco_selector_group,
        profile="assembly_with_exclusions",
        selector_defaults=selector_defaults,
        context_label="list assemblies",
    )
    list_busco_selector_group.add_argument(
        "-m",
        "--meta",
        "--metadata",
        dest="meta",
        nargs="?",
        const="__DEFAULT__",
        help="Include metadata columns. Omit value for defaults (release_date, level, n50, comments).",
    )
    _add_assembly_list_output_options(list_busco_output_group)
    _add_row_sort_option(
        list_busco_output_group,
        fields="accession, species/taxon, rank columns, metadata columns, latest, quality, BUSCO columns when shown",
    )
    list_busco.set_defaults(handler=list_handler)

    list_busco_runs = list_subparsers.add_parser("busco-runs", help="List BUSCO run records.")
    list_busco_runs_selector_group = list_busco_runs.add_argument_group("Selector options")
    list_busco_runs_output_group = list_busco_runs.add_argument_group("Output options")
    _add_selector_arguments(
        list_busco_runs_selector_group,
        profile="busco_listing",
        selector_defaults=selector_defaults,
        context_label="BUSCO run listing",
    )
    list_busco_runs_output_group.add_argument("--ids-only", action="store_true", help="Emit only the filtered BUSCO run ids as a comma-separated list.")
    list_busco_runs_output_group.add_argument("--store-results", help="Store filtered BUSCO run ids in a database variable.")
    _add_basic_list_output_options(list_busco_runs_output_group)
    _add_row_sort_option(
        list_busco_runs_output_group,
        fields="run_id, accession, species, library, pipeline, orthofinder_target_library, input_mode, status, complete, single_copy, busco.single_copy_complete, quality, completed_at",
    )
    list_busco_runs.set_defaults(handler=list_handler)

    list_proteome_profiles = list_subparsers.add_parser("proteome-profiles", help="List registered proteome profiles.")
    list_proteome_profiles_selector_group = list_proteome_profiles.add_argument_group("Selector options")
    list_proteome_profiles_output_group = list_proteome_profiles.add_argument_group("Output options")
    _add_selector_arguments(
        list_proteome_profiles_selector_group,
        profile="assembly_with_exclusions",
        selector_defaults=selector_defaults,
        context_label="proteome profile listing",
        skip_fields={"has_busco_results", "missing_busco_results", "library_id", "library_name", "busco_pipeline", "prefer_busco_pipeline", "busco_input_mode", "prefer_busco_input_mode", "busco_export_format", "busco_run_ids", "busco_run_selection", "busco_complete_min", "busco_single_min", "include_paralog_filtering_in_score", "include_decontamination_in_score", "use_decontamination_run", "use_paralog_run", "allow_ambiguous_contaminants", "strict_decontamination", "rescue_duplicates", "paralog_filtered", "not_paralog_filtered", "min_hidden_paralogs", "max_hidden_paralogs", "decontaminated", "not_decontaminated", "contaminated", "decontamination_run", "ignore_contaminated_assemblies", "proteome_profile", "prefer_proteome_profile", "isoforms_cleaned", "raw_proteome"},
    )
    _add_basic_list_output_options(list_proteome_profiles_output_group)
    _add_row_sort_option(
        list_proteome_profiles_output_group,
        fields="profile_id, accession, profile_name, kind, is_default, status, sequence_count, preparation_type, created_at",
    )
    list_proteome_profiles.set_defaults(handler=list_handler)

    list_buscos = list_subparsers.add_parser("buscos", help="List BUSCO family rows.")
    list_buscos_selector_group = list_buscos.add_argument_group("Selector options")
    list_buscos_output_group = list_buscos.add_argument_group("Output options")
    _add_selector_arguments(
        list_buscos_selector_group,
        profile="busco_listing",
        selector_defaults=selector_defaults,
        context_label="BUSCO family listing",
    )
    list_buscos_selector_group.add_argument("--family-id", action=AppendCommaSeparated, help="Limit BUSCO family listing to one or more BUSCO family ids.")
    _add_basic_list_output_options(list_buscos_output_group)
    list_buscos.set_defaults(handler=list_handler)

    return list_parser


def register_watch_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    handler=None,
) -> argparse.ArgumentParser:
    """Register the ``watch`` command as a thin live alias over queue/error listings."""

    watch_handler = handler or _handle_watch
    watch_parser = subparsers.add_parser("watch", help="Watch the live queue or task errors.")
    watch_subparsers = watch_parser.add_subparsers(dest="watch_choice", required=True)

    watch_queue = watch_subparsers.add_parser("queue", help="Watch the live task queue.")
    watch_queue_group = watch_queue.add_argument_group("Options")
    _add_queue_monitor_options(watch_queue_group, include_watch=False)
    watch_queue.set_defaults(handler=watch_handler)

    watch_errors = watch_subparsers.add_parser("errors", help="Watch live task errors.")
    watch_errors_group = watch_errors.add_argument_group("Options")
    _add_error_monitor_options(watch_errors_group, include_watch=False)
    watch_errors.set_defaults(handler=watch_handler)

    return watch_parser


def register_count_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
) -> argparse.ArgumentParser:
    """Register the top-level ``count`` command."""

    count_parser = subparsers.add_parser("count", help="Count items matching selector filters.")
    count_selector_group = count_parser.add_argument_group("Selector options")
    count_output_group = count_parser.add_argument_group("Output options")
    count_parser.add_argument("subject", choices=["assemblies"], help="What to count.")
    _add_selector_arguments(
        count_selector_group,
        profile="assembly_with_exclusions",
        selector_defaults=selector_defaults,
        context_label="assemblies count",
    )
    count_output_group.add_argument("--group-by-rank", action="store_true", help="Group count output by rank header.")
    count_parser.set_defaults(handler=_handle_count)
    return count_parser


def register_assemblies_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
) -> argparse.ArgumentParser:
    """Register the top-level ``assemblies`` selector command."""

    assemblies_parser = subparsers.add_parser("assemblies", help="List assemblies matching selector filters.")
    assemblies_selector_group = assemblies_parser.add_argument_group("Selector options")
    assemblies_output_group = assemblies_parser.add_argument_group("Output options")
    _add_selector_arguments(
        assemblies_selector_group,
        profile="assembly",
        selector_defaults=selector_defaults,
        context_label="assemblies output",
    )
    assemblies_output_group.add_argument(
        *LIST_OUTPUT_SHORT_ALIASES["store"],
        "--store",
        "--save-set",
        dest="store_variable",
        help="Store resolved accessions in a named set (use with @NAME later).",
    )
    assemblies_output_group.add_argument(
        *LIST_OUTPUT_SHORT_ALIASES["append"],
        "--append-to",
        "--append-set",
        dest="append_to_variable",
        help="Union resolved accessions into a named set.",
    )
    assemblies_output_group.add_argument(
        "--intersection",
        action=AppendCommaSeparated,
        help="Intersect resolved accessions with an explicit accession or @VARIABLE set.",
    )
    assemblies_parser.set_defaults(handler=_handle_assemblies)
    return assemblies_parser


__all__ = [
    "METADATA_FIELDS",
    "_ACCESSIONISH_TOKEN_RE",
    "_classify_variable_kind",
    "_handle_assemblies",
    "_handle_count",
    "_handle_list",
    "_handle_list_assemblies",
    "_handle_list_busco_runs",
    "_handle_list_buscos",
    "_handle_watch",
    "register_assemblies_parser",
    "register_count_parser",
    "register_list_parser",
    "register_watch_parser",
]
