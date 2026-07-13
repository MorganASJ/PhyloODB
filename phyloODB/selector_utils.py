"""Shared helpers for resolving accession selectors across CLI and tasks."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import re
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, TYPE_CHECKING

from .database import DBManager
from .proteome_profile_utils import resolve_profile_selector

if TYPE_CHECKING:  # pragma: no cover - typing helper only
    from .registry import spec as spec_module  # type: ignore
    TaskSpec = spec_module.TaskSpec
else:  # pragma: no cover - simplify typing at runtime
    TaskSpec = object


def _ensure_sequence(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return list(value)
    if isinstance(value, (tuple, set)):
        return list(value)
    return [value]


def _normalize_busco_input_mode(value: Any) -> Optional[str]:
    if value is None:
        return None
    token = str(value).strip().lower()
    if not token:
        return None
    if token in {"nucl", "nucleotide"}:
        return "genome"
    if token == "prot":
        return "protein"
    return token


def expand_busco_run_id_variables(manager, tokens: Sequence[Any]) -> List[int]:
    resolved: List[int] = []
    for token in tokens or []:
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


def _coerce_optional_flag(value: Any) -> Optional[bool]:
    coerced = _coerce_bool(value)
    if coerced is None:
        return None
    return True if coerced else None


@dataclass(slots=True)
class SelectorRequest:
    """Canonical selector input shared by CLI, tasks, and views."""

    accessions: List[Any] = field(default_factory=list)
    root: Optional[str] = None
    clade: List[str] = field(default_factory=list)
    taxid: Any = None
    allow_all: bool = False
    exclude_accessions: List[Any] = field(default_factory=list)
    exclude_taxids: List[Any] = field(default_factory=list)
    exclude_clades: List[Any] = field(default_factory=list)
    downloaded_only: Optional[bool] = None
    not_downloaded: Optional[bool] = None
    local_only: Optional[bool] = None
    not_local: Optional[bool] = None
    after: Optional[str] = None
    before: Optional[str] = None
    level: Optional[str] = None
    primary_only: Optional[bool] = None
    protein_only: Optional[bool] = None
    status_min: Any = None
    filters: Any = None
    has_busco_results: Optional[bool] = None
    missing_busco_results: Optional[bool] = None
    busco_library_id: Any = None
    library_name: Optional[str] = None
    busco_pipeline: Optional[str] = None
    busco_input_mode: Optional[str] = None
    prefer_busco_pipeline: Optional[str] = None
    prefer_busco_input_mode: Optional[str] = None
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    busco_export_format: Optional[str] = None
    busco_run_ids: List[str] = field(default_factory=list)
    busco_run_selection: Optional[str] = None
    busco_complete_min: Any = None
    busco_single_min: Any = None
    include_paralog_filtering_in_score: Optional[bool] = None
    include_decontamination_in_score: Optional[bool] = None
    paralog_run_id: Optional[str] = None
    decontamination_run_id: Optional[str] = None
    allow_ambiguous_contaminants: Optional[bool] = None
    strict_decontamination: Optional[bool] = None
    rescue_duplicates: Optional[bool] = None
    paralog_filtered: Optional[bool] = None
    not_paralog_filtered: Optional[bool] = None
    min_hidden_paralogs: Any = None
    max_hidden_paralogs: Any = None
    decontaminated: Optional[bool] = None
    not_decontaminated: Optional[bool] = None
    contaminated: Optional[bool] = None
    decontamination_run: Optional[str] = None
    ignore_contaminated_assemblies: Optional[bool] = None
    quantity: Any = None
    rank: Optional[str] = None
    ranks: List[str] = field(default_factory=list)
    quantities: List[Any] = field(default_factory=list)
    sample_strategy: Optional[str] = None
    sample_seed: Any = None
    allow_duplicate_species: Optional[bool] = None
    use_busco: Optional[bool] = None

    @classmethod
    def from_mapping(cls, source: Optional[Mapping[str, Any]]) -> "SelectorRequest":
        """Build a request from a mapping, collapsing known selector aliases."""

        data = dict(source or {})
        filters = data.get("filters")
        if filters is None:
            filters = data.get("filter")
        after = data.get("after")
        if after is None:
            after = data.get("released_after")
        before = data.get("before")
        if before is None:
            before = data.get("released_before")
        busco_library_id = data.get("busco_library_id")
        if busco_library_id is None:
            busco_library_id = data.get("library_id")
        busco_pipeline = data.get("busco_pipeline")
        if busco_pipeline is None:
            busco_pipeline = data.get("require_busco_pipeline")
        busco_input_mode = data.get("busco_input_mode")
        if busco_input_mode is None:
            busco_input_mode = data.get("format")
        if busco_input_mode is None:
            busco_input_mode = data.get("require_format")
        if busco_input_mode is None:
            busco_input_mode = data.get("input_mode")
        prefer_busco_pipeline = data.get("prefer_busco_pipeline")
        prefer_busco_input_mode = data.get("prefer_busco_input_mode")
        if prefer_busco_input_mode is None:
            prefer_busco_input_mode = data.get("prefer_format")
        proteome_profile = resolve_profile_selector(
            proteome_profile=data.get("proteome_profile"),
            isoforms_cleaned=_coerce_bool(data.get("isoforms_cleaned")),
            raw_proteome=_coerce_bool(data.get("raw_proteome")),
        )
        prefer_proteome_profile = data.get("prefer_proteome_profile")
        busco_run_ids = data.get("busco_run_ids")
        if busco_run_ids is None:
            busco_run_ids = data.get("run_ids")
        if busco_run_ids is None:
            busco_run_ids = data.get("run_id")
        busco_export_format = data.get("busco_export_format")
        if busco_export_format is None:
            busco_export_format = data.get("export_format")
        busco_run_selection = data.get("busco_run_selection")
        quantity = data.get("quantity")
        if quantity is None:
            quantity = data.get("rule_quantity")
        rank = data.get("rank")
        if rank is None:
            rank = data.get("rule_rank")
        decontamination_run_id = data.get("decontamination_run_id")
        if decontamination_run_id is None:
            decontamination_run_id = data.get("use_decontamination_run")
        paralog_run_id = data.get("paralog_run_id")
        if paralog_run_id is None:
            paralog_run_id = data.get("use_paralog_run")
        busco_complete_min = data.get("busco_complete_min")
        if busco_complete_min is None:
            busco_complete_min = data.get("min_completeness")
        busco_single_min = data.get("busco_single_min")
        if busco_single_min is None:
            busco_single_min = data.get("min_single_copy_complete")

        return cls(
            accessions=_ensure_sequence(data.get("accessions")),
            root=str(data.get("root")).strip() if data.get("root") not in (None, "") else None,
            clade=_split_tokens(data.get("clade")),
            taxid=data.get("taxid"),
            allow_all=bool(data.get("allow_all", data.get("all", False))),
            exclude_accessions=_ensure_sequence(data.get("exclude_accessions")),
            exclude_taxids=_ensure_sequence(data.get("exclude_taxids")),
            exclude_clades=_ensure_sequence(data.get("exclude_clades")),
            downloaded_only=_coerce_optional_flag(data.get("downloaded_only")),
            not_downloaded=_coerce_optional_flag(data.get("not_downloaded")),
            local_only=_coerce_optional_flag(data.get("local_only")),
            not_local=_coerce_optional_flag(data.get("not_local")),
            after=after,
            before=before,
            level=data.get("level"),
            primary_only=_coerce_optional_flag(data.get("primary_only")),
            protein_only=_coerce_optional_flag(data.get("protein_only")),
            status_min=data.get("status_min"),
            filters=filters,
            has_busco_results=_coerce_optional_flag(data.get("has_busco_results")),
            missing_busco_results=_coerce_optional_flag(data.get("missing_busco_results")),
            busco_library_id=busco_library_id,
            library_name=data.get("library_name"),
            busco_pipeline=str(busco_pipeline).strip().lower() if busco_pipeline not in (None, "") else None,
            busco_input_mode=_normalize_busco_input_mode(busco_input_mode),
            prefer_busco_pipeline=str(prefer_busco_pipeline).strip().lower() if prefer_busco_pipeline not in (None, "") else None,
            prefer_busco_input_mode=_normalize_busco_input_mode(prefer_busco_input_mode),
            proteome_profile=str(proteome_profile).strip() if proteome_profile not in (None, "") else None,
            prefer_proteome_profile=str(prefer_proteome_profile).strip() if prefer_proteome_profile not in (None, "") else None,
            busco_export_format=str(busco_export_format).strip().lower() if busco_export_format not in (None, "") else None,
            busco_run_ids=[str(token).strip() for token in _ensure_sequence(busco_run_ids) if str(token).strip()],
            busco_run_selection=busco_run_selection,
            busco_complete_min=busco_complete_min,
            busco_single_min=busco_single_min,
            include_paralog_filtering_in_score=_coerce_bool(data.get("include_paralog_filtering_in_score")),
            include_decontamination_in_score=_coerce_bool(data.get("include_decontamination_in_score")),
            paralog_run_id=paralog_run_id,
            decontamination_run_id=decontamination_run_id,
            allow_ambiguous_contaminants=_coerce_bool(data.get("allow_ambiguous_contaminants")),
            strict_decontamination=_coerce_bool(data.get("strict_decontamination")),
            rescue_duplicates=_coerce_optional_flag(data.get("rescue_duplicates")),
            paralog_filtered=_coerce_optional_flag(data.get("paralog_filtered")),
            not_paralog_filtered=_coerce_optional_flag(data.get("not_paralog_filtered")),
            min_hidden_paralogs=data.get("min_hidden_paralogs"),
            max_hidden_paralogs=data.get("max_hidden_paralogs"),
            decontaminated=_coerce_optional_flag(data.get("decontaminated")),
            not_decontaminated=_coerce_optional_flag(data.get("not_decontaminated")),
            contaminated=_coerce_optional_flag(data.get("contaminated")),
            decontamination_run=data.get("decontamination_run"),
            ignore_contaminated_assemblies=_coerce_bool(data.get("ignore_contaminated_assemblies")),
            quantity=quantity,
            rank=str(rank).strip() if rank is not None else None,
            ranks=_split_tokens(data.get("ranks")),
            quantities=_ensure_sequence(data.get("quantities")),
            sample_strategy=data.get("sample_strategy"),
            sample_seed=data.get("sample_seed"),
            allow_duplicate_species=_coerce_optional_flag(data.get("allow_duplicate_species")),
            use_busco=_coerce_bool(data.get("use_busco")),
        )

    @classmethod
    def from_namespace(cls, namespace: Any) -> "SelectorRequest":
        """Build a request from an argparse-style namespace."""

        return cls.from_mapping(vars(namespace))

    def with_overrides(self, **overrides: Any) -> "SelectorRequest":
        """Return a copy of the request with normalised override values applied."""

        merged = self.as_mapping()
        merged.update(overrides)
        return type(self).from_mapping(merged)

    def as_mapping(self) -> Dict[str, Any]:
        """Return the canonical selector mapping used by shared resolver helpers."""

        return {
            "accessions": list(self.accessions),
            "root": self.root,
            "clade": list(self.clade),
            "taxid": self.taxid,
            "allow_all": self.allow_all,
            "exclude_accessions": list(self.exclude_accessions),
            "exclude_taxids": list(self.exclude_taxids),
            "exclude_clades": list(self.exclude_clades),
            "downloaded_only": self.downloaded_only,
            "not_downloaded": self.not_downloaded,
            "local_only": self.local_only,
            "not_local": self.not_local,
            "after": self.after,
            "before": self.before,
            "level": self.level,
            "primary_only": self.primary_only,
            "protein_only": self.protein_only,
            "status_min": self.status_min,
            "filters": self.filters,
            "has_busco_results": self.has_busco_results,
            "missing_busco_results": self.missing_busco_results,
            "busco_library_id": self.busco_library_id,
            "library_name": self.library_name,
            "busco_pipeline": self.busco_pipeline,
            "busco_input_mode": self.busco_input_mode,
            "prefer_busco_pipeline": self.prefer_busco_pipeline,
            "prefer_busco_input_mode": self.prefer_busco_input_mode,
            "proteome_profile": self.proteome_profile,
            "prefer_proteome_profile": self.prefer_proteome_profile,
            "busco_export_format": self.busco_export_format,
            "busco_run_ids": list(self.busco_run_ids),
            "busco_run_selection": self.busco_run_selection,
            "busco_complete_min": self.busco_complete_min,
            "busco_single_min": self.busco_single_min,
            "include_paralog_filtering_in_score": self.include_paralog_filtering_in_score,
            "include_decontamination_in_score": self.include_decontamination_in_score,
            "paralog_run_id": self.paralog_run_id,
            "decontamination_run_id": self.decontamination_run_id,
            "allow_ambiguous_contaminants": self.allow_ambiguous_contaminants,
            "strict_decontamination": self.strict_decontamination,
            "rescue_duplicates": self.rescue_duplicates,
            "paralog_filtered": self.paralog_filtered,
            "not_paralog_filtered": self.not_paralog_filtered,
            "min_hidden_paralogs": self.min_hidden_paralogs,
            "max_hidden_paralogs": self.max_hidden_paralogs,
            "decontaminated": self.decontaminated,
            "not_decontaminated": self.not_decontaminated,
            "contaminated": self.contaminated,
            "decontamination_run": self.decontamination_run,
            "ignore_contaminated_assemblies": self.ignore_contaminated_assemblies,
            "quantity": self.quantity,
            "rank": self.rank,
            "ranks": list(self.ranks),
            "quantities": list(self.quantities),
            "sample_strategy": self.sample_strategy,
            "sample_seed": self.sample_seed,
            "allow_duplicate_species": self.allow_duplicate_species,
            "use_busco": self.use_busco,
        }

    def has_busco_filters(self) -> bool:
        """Return True when BUSCO or decontamination-related selectors are active."""

        return any(
            [
                self.has_busco_results,
                self.missing_busco_results,
                self.busco_library_id is not None,
                self.library_name,
                self.busco_pipeline,
                self.busco_input_mode,
                self.prefer_busco_pipeline,
                self.prefer_busco_input_mode,
                self.proteome_profile,
                self.prefer_proteome_profile,
                self.busco_export_format,
                bool(self.busco_run_ids),
                self.busco_complete_min is not None,
                self.busco_single_min is not None,
                self.include_paralog_filtering_in_score is not None,
                self.include_decontamination_in_score is not None,
                self.decontamination_run_id is not None,
                self.allow_ambiguous_contaminants is not None,
                self.strict_decontamination is not None,
                self.rescue_duplicates,
                self.paralog_filtered,
                self.not_paralog_filtered,
                self.min_hidden_paralogs is not None,
                self.max_hidden_paralogs is not None,
                self.decontaminated,
                self.not_decontaminated,
                self.contaminated,
                self.decontamination_run is not None,
                self.ignore_contaminated_assemblies is not None,
            ]
        )

    def has_rule_selectors(self) -> bool:
        """Return True when rank- or quantity-based selection is active."""

        return any(
            [
                self.quantity is not None,
                bool(self.rank),
                bool(self.ranks),
                bool(self.quantities),
                self.sample_strategy is not None,
                self.sample_seed is not None,
                self.allow_duplicate_species,
            ]
        )

    def has_candidate_selectors(self) -> bool:
        """Return True when base accession-selection inputs are present."""

        return any(
            [
                bool(self.accessions),
                self.root is not None,
                bool(self.clade),
                self.taxid is not None,
                self.allow_all,
                bool(self.exclude_accessions),
                bool(self.exclude_taxids),
                bool(self.exclude_clades),
                self.downloaded_only,
                self.not_downloaded,
                self.local_only,
                self.not_local,
                self.after,
                self.before,
                self.level,
                self.primary_only,
                self.protein_only,
                self.status_min is not None,
                bool(self.filters),
            ]
        )

    def has_any_selector_input(self) -> bool:
        """Return True when any selector input is active."""

        return self.has_candidate_selectors() or self.has_busco_filters() or self.has_rule_selectors()


def normalize_accessions(accessions: Sequence[Any]) -> List[str]:
    """Return a flat list of accession strings, discarding None values."""
    cleaned: List[str] = []
    for entry in accessions:
        if entry is None:
            continue
        if isinstance(entry, (list, tuple)):
            if entry:
                cleaned.append(str(entry[0]))
        else:
            cleaned.append(str(entry))
    return cleaned


def _coerce_bool(value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    token = str(value).strip().lower()
    if token in {"1", "true", "yes", "y", "on"}:
        return True
    if token in {"0", "false", "no", "n", "off"}:
        return False
    return None


FILTER_FIELD_ALIASES = {
    "level": "assembly_level",
    "assembly_level": "assembly_level",
    "n50": "contig_n50",
    "contig_n50": "contig_n50",
    "origin": "origin",
    "source": "origin",
    "release": "release_date",
    "release_date": "release_date",
    "download_date": "dl_date",
    "dl_date": "dl_date",
    "proteome_profile": "proteome_profile",
    "profile": "proteome_profile",
    "default_proteome_profile": "default_proteome_profile",
    "default_profile": "default_proteome_profile",
}

BUSCO_FIELD_ALIASES = {
    "complete": "complete",
    "quality": "single_copy_complete",
    "single_copy": "single_copy_complete",
    "single_copy_complete": "single_copy_complete",
    "duplicated": "duplicated",
    "fragmented": "fragmented",
    "missing": "missing",
    "sc": "single_copy_complete",
    "proteome_profile": "proteome_profile",
    "profile": "proteome_profile",
    "default_proteome_profile": "default_proteome_profile",
    "default_profile": "default_proteome_profile",
}

FILTER_OP_PATTERN = re.compile(
    r"^\s*([A-Za-z0-9_.]+)\s*(>=|<=|!=|==|=|>|<|contains|notcontains|not\\s+contains|isnull|notnull|~|!~)\s*(.*)\s*$",
    re.IGNORECASE,
)


def _split_outside_quotes(text: str, delimiter: str) -> List[str]:
    if text is None:
        return []
    buf: List[str] = []
    parts: List[str] = []
    in_single = False
    in_double = False
    for ch in str(text):
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        if ch == delimiter and not in_single and not in_double:
            token = "".join(buf).strip()
            if token:
                parts.append(token)
            buf = []
        else:
            buf.append(ch)
    token = "".join(buf).strip()
    if token:
        parts.append(token)
    return parts


def _normalize_filter_tokens(filters: Any) -> List[str]:
    if not filters:
        return []
    tokens: List[str] = []
    if isinstance(filters, (list, tuple, set)):
        items = filters
    else:
        items = [filters]
    for entry in items:
        if entry is None:
            continue
        if isinstance(entry, (list, tuple, set)):
            tokens.extend(_normalize_filter_tokens(entry))
        else:
            tokens.append(str(entry))
    return tokens


def _strip_quotes(value: str) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    if len(s) >= 2 and ((s[0] == s[-1]) and s[0] in {"'", '"'}):
        return s[1:-1]
    return s


def _parse_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    s = str(value).strip()
    if not s:
        return None
    try:
        return float(s)
    except (TypeError, ValueError, OverflowError):
        return None


def _parse_filter_condition(token: str) -> Dict[str, Any]:
    match = FILTER_OP_PATTERN.match(token or "")
    if not match:
        raise ValueError(f"Invalid filter expression '{token}'.")
    raw_field = match.group(1).strip()
    op_raw = match.group(2).strip().lower()
    value_raw = match.group(3).strip()

    if op_raw in {"contains", "~"}:
        op = "contains"
    elif op_raw in {"notcontains", "not contains", "!~"}:
        op = "not_contains"
    elif op_raw in {"=", "=="}:
        op = "eq"
    elif op_raw == "!=":
        op = "ne"
    elif op_raw in {">", ">=", "<", "<="}:
        op = op_raw
    elif op_raw in {"isnull", "notnull"}:
        op = op_raw
    else:
        raise ValueError(f"Unsupported filter operator '{op_raw}'.")

    prefix = None
    field_name = raw_field
    if "." in raw_field:
        prefix, field_name = raw_field.split(".", 1)
        prefix = prefix.strip().lower()
        field_name = field_name.strip()
        if prefix not in {"busco", "genome", "assembly"}:
            raise ValueError(f"Unknown filter prefix '{prefix}' in '{token}'.")
    elif field_name.strip().lower() == "quality":
        prefix = "busco"

    if prefix == "busco":
        field_key = BUSCO_FIELD_ALIASES.get(field_name.lower(), field_name.lower())
    else:
        if prefix is None:
            field_key = FILTER_FIELD_ALIASES.get(field_name.lower(), field_name.lower())
        else:
            field_key = field_name.lower()

    if op in {"isnull", "notnull"}:
        value_clean = None
    else:
        if not value_raw:
            raise ValueError(f"Filter '{token}' is missing a value.")
        value_clean = _strip_quotes(value_raw)

    value_num = None
    if op in {">", ">=", "<", "<="}:
        value_num = _parse_float(value_clean)
        if value_num is None:
            raise ValueError(f"Filter '{token}' requires a numeric value.")
    elif op in {"eq", "ne"}:
        value_num = _parse_float(value_clean)

    return {
        "raw": token,
        "prefix": prefix,
        "field": field_key,
        "op": op,
        "value": value_clean,
        "value_num": value_num,
        "value_str": str(value_clean).lower() if value_clean is not None else "",
    }


def _parse_filter_expression(expr: str) -> List[List[Dict[str, Any]]]:
    if not expr:
        return []
    or_groups = _split_outside_quotes(expr, "|")
    groups: List[List[Dict[str, Any]]] = []
    for group in or_groups:
        and_parts = _split_outside_quotes(group, ",")
        if not and_parts:
            continue
        groups.append([_parse_filter_condition(part) for part in and_parts])
    return groups


def validate_filter_expressions(
    manager: DBManager,
    filters: Any,
) -> None:
    """Validate selector filter syntax and field names without resolving rows."""

    tokens = _normalize_filter_tokens(filters)
    if not tokens:
        return
    parsed_filters: List[List[List[Dict[str, Any]]]] = []
    for token in tokens:
        groups = _parse_filter_expression(token)
        if groups:
            parsed_filters.append(groups)
    if not parsed_filters:
        return

    genome_cols = set(_get_table_columns(manager, "Genome"))
    assembly_cols = set(_get_table_columns(manager, "Assembly"))
    custom_busco_fields = {"proteome_profile", "default_proteome_profile"}
    for groups in parsed_filters:
        for group in groups:
            for cond in group:
                prefix = cond["prefix"]
                field = cond["field"]
                if prefix == "busco":
                    if field not in set(BUSCO_FIELD_ALIASES.values()) | custom_busco_fields:
                        raise ValueError(f"Unknown BUSCO filter field '{field}'.")
                    continue
                if field in custom_busco_fields:
                    continue
                if prefix == "genome":
                    if field not in genome_cols:
                        raise ValueError(f"Unknown genome filter field '{field}'.")
                    continue
                if prefix == "assembly":
                    if field not in assembly_cols:
                        raise ValueError(f"Unknown assembly filter field '{field}'.")
                    continue
                if field not in genome_cols and field not in assembly_cols:
                    raise ValueError(f"Unknown filter field '{field}'.")


def _get_table_columns(manager: DBManager, table: str) -> List[str]:
    cache_key = f"_col_cache_{table}"
    cached = getattr(manager, cache_key, None)
    if cached is not None:
        return cached
    manager.cursor.execute(f"PRAGMA table_info({table})")
    cols = [row[1].lower() for row in (manager.cursor.fetchall() or [])]
    setattr(manager, cache_key, cols)
    return cols


def _fetch_table_rows(
    manager: DBManager,
    *,
    table: str,
    accessions: Sequence[str],
    columns: Sequence[str],
) -> Dict[str, Dict[str, Any]]:
    if not accessions or not columns:
        return {}
    cols = ["accession"] + sorted(set(columns) - {"accession"})
    placeholders = ",".join("?" for _ in accessions)
    sql = f"SELECT {', '.join(cols)} FROM {table} WHERE accession IN ({placeholders})"
    manager.cursor.execute(sql, tuple(accessions))
    rows = manager.cursor.fetchall() or []
    result: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        acc = str(row[0])
        data: Dict[str, Any] = {}
        for idx, col in enumerate(cols[1:], start=1):
            data[col] = row[idx]
        result[acc] = data
    return result


def _compute_best_busco_map(
    manager: DBManager,
    accessions: Sequence[str],
    *,
    library_id: Optional[int] = None,
    library_name: Optional[str] = None,
) -> Dict[str, Dict[str, Optional[float]]]:
    if not accessions:
        return {}
    rows = manager.get_busco_results_percentages(
        accessions=list(accessions),
        library_id=library_id,
        library_name=library_name,
    )
    if not rows:
        return {}

    if library_id is not None or library_name is not None:
        best_map: Dict[str, Dict[str, Optional[float]]] = {}
        for acc, _species, _lib_name, complete, sc_complete, duplicated, fragmented, missing in rows:
            best_map[str(acc)] = {
                "complete": complete,
                "single_copy_complete": sc_complete,
                "duplicated": duplicated,
                "fragmented": fragmented,
                "missing": missing,
            }
        return best_map

    libraries = manager.get_libraries() or []
    library_meta: Dict[str, Dict[str, Any]] = {}
    for row in libraries:
        lib_id = row[0]
        lib_name = row[1]
        taxid = row[3]
        library_meta[str(lib_name)] = {"library_id": lib_id, "taxid": taxid}

    acc_taxid: Dict[str, Optional[int]] = {}
    placeholders = ",".join("?" for _ in accessions)
    manager.cursor.execute(
        f"SELECT accession, taxid FROM Genome WHERE accession IN ({placeholders})",
        tuple(accessions),
    )
    for acc, taxid in manager.cursor.fetchall() or []:
        acc_taxid[str(acc)] = int(taxid) if taxid is not None else None

    busco_rows: Dict[str, List[tuple]] = {}
    for row in rows:
        acc = str(row[0])
        busco_rows.setdefault(acc, []).append(row[2:])

    lineage_cache: Dict[int, Dict[int, int]] = {}

    def _lineage_depth_map(tid: int) -> Dict[int, int]:
        if tid in lineage_cache:
            return lineage_cache[tid]
        depth_map: Dict[int, int] = {}
        rows_lineage = manager.get_lineage_root_to_leaf(tid) or []
        for depth, row in enumerate(rows_lineage):
            try:
                depth_map[int(row[0])] = depth
            except (TypeError, ValueError, IndexError):
                continue
        lineage_cache[tid] = depth_map
        return depth_map

    best_map: Dict[str, Dict[str, Optional[float]]] = {}
    for acc, rows_acc in busco_rows.items():
        if not rows_acc:
            continue
        acc_tid = acc_taxid.get(acc)
        preferred = None
        preferred_depth = -1
        if acc_tid is not None:
            depth_map = _lineage_depth_map(int(acc_tid))
            for lib_name, complete, sc_complete, duplicated, fragmented, missing in rows_acc:
                lib_taxid = library_meta.get(str(lib_name), {}).get("taxid")
                if lib_taxid is None:
                    continue
                if lib_taxid in depth_map and depth_map[lib_taxid] > preferred_depth:
                    preferred = (complete, sc_complete, duplicated, fragmented, missing)
                    preferred_depth = depth_map[lib_taxid]
        if preferred is None:
            best = None
            for _lib_name, complete, sc_complete, duplicated, fragmented, missing in rows_acc:
                cand = (complete, sc_complete, duplicated, fragmented, missing)
                if best is None:
                    best = cand
                    continue
                sc_best = best[1] if best[1] is not None else -1.0
                sc_cand = cand[1] if cand[1] is not None else -1.0
                c_best = best[0] if best[0] is not None else -1.0
                c_cand = cand[0] if cand[0] is not None else -1.0
                if sc_cand > sc_best or (sc_cand == sc_best and c_cand > c_best):
                    best = cand
            preferred = best
        if preferred is None:
            continue
        best_map[acc] = {
            "complete": preferred[0],
            "single_copy_complete": preferred[1],
            "duplicated": preferred[2],
            "fragmented": preferred[3],
            "missing": preferred[4],
        }
    return best_map


def _compute_busco_profile_map(
    manager: DBManager,
    accessions: Sequence[str],
    *,
    library_id: Optional[int] = None,
    library_name: Optional[str] = None,
) -> Dict[str, Dict[str, Optional[str]]]:
    if not accessions:
        return {}

    resolved_library_id = int(library_id) if library_id is not None else None
    if resolved_library_id is None and library_name:
        resolved = manager.get_library_id(str(library_name))
        if resolved is not None:
            resolved_library_id = int(resolved)

    result: Dict[str, Dict[str, Optional[str]]] = {str(acc): {} for acc in accessions}

    placeholders = ",".join("?" for _ in accessions)
    manager.cursor.execute(
        f"""
        SELECT accession, profile_name
        FROM Proteome_Profiles
        WHERE accession IN ({placeholders}) AND COALESCE(is_default, 0) = 1
        ORDER BY accession ASC, proteome_profile_id DESC
        """,
        tuple(str(acc) for acc in accessions),
    )
    for accession, profile_name in manager.cursor.fetchall() or []:
        result.setdefault(str(accession), {})["default_proteome_profile"] = str(profile_name) if profile_name else None

    if resolved_library_id is None:
        return result

    run_map = manager.busco._resolve_busco_runs_for_query(
        int(resolved_library_id),
        accessions=[str(acc) for acc in accessions],
    )

    run_ids = sorted({int(run_id) for run_id in run_map.values() if run_id is not None})
    if not run_ids:
        return result

    run_placeholders = ",".join("?" for _ in run_ids)
    manager.cursor.execute(
        f"""
        SELECT r.accession, pp.profile_name
        FROM BUSCO_Runs r
        LEFT JOIN Proteome_Profiles pp ON pp.proteome_profile_id = r.proteome_profile_id
        WHERE r.run_id IN ({run_placeholders})
        """,
        tuple(run_ids),
    )
    for accession, profile_name in manager.cursor.fetchall() or []:
        result.setdefault(str(accession), {})["proteome_profile"] = str(profile_name) if profile_name else None

    return result


def _evaluate_filter_condition(value: Any, cond: Dict[str, Any]) -> bool:
    op = cond["op"]
    if op == "isnull":
        return value is None or value == ""
    if op == "notnull":
        return value is not None and value != ""
    if op == "contains":
        if value is None:
            return False
        return cond["value_str"] in str(value).lower()
    if op == "not_contains":
        if value is None:
            return True
        return cond["value_str"] not in str(value).lower()

    if op in {">", ">=", "<", "<="}:
        lhs = _parse_float(value)
        rhs = cond.get("value_num")
        if lhs is None or rhs is None:
            return False
        if op == ">":
            return lhs > rhs
        if op == ">=":
            return lhs >= rhs
        if op == "<":
            return lhs < rhs
        return lhs <= rhs

    if op in {"eq", "ne"}:
        lhs_num = _parse_float(value)
        rhs_num = cond.get("value_num")
        if lhs_num is not None and rhs_num is not None:
            result = lhs_num == rhs_num
        else:
            lhs_str = "" if value is None else str(value).lower()
            result = lhs_str == cond["value_str"]
        return result if op == "eq" else not result

    return False


def filter_accessions_by_expressions(
    manager: DBManager,
    accessions: Sequence[Any],
    filters: Any,
    *,
    busco_library_id: Optional[int] = None,
    busco_library_name: Optional[str] = None,
) -> List[str]:
    tokens = _normalize_filter_tokens(filters)
    if not tokens:
        return normalize_accessions(accessions)

    parsed_filters: List[List[List[Dict[str, Any]]]] = []
    for token in tokens:
        groups = _parse_filter_expression(token)
        if groups:
            parsed_filters.append(groups)

    if not parsed_filters:
        return normalize_accessions(accessions)

    pool = normalize_accessions(accessions)
    if not pool:
        return []

    genome_cols = set(_get_table_columns(manager, "Genome"))
    assembly_cols = set(_get_table_columns(manager, "Assembly"))
    custom_busco_fields = {"proteome_profile", "default_proteome_profile"}
    required_genome: set[str] = set()
    required_assembly: set[str] = set()
    needs_busco = False

    for groups in parsed_filters:
        for group in groups:
            for cond in group:
                prefix = cond["prefix"]
                field = cond["field"]
                if prefix == "busco":
                    if field not in set(BUSCO_FIELD_ALIASES.values()) | custom_busco_fields:
                        raise ValueError(f"Unknown BUSCO filter field '{field}'.")
                    needs_busco = True
                    continue
                if field in custom_busco_fields:
                    needs_busco = True
                    continue
                if prefix == "genome":
                    if field not in genome_cols:
                        raise ValueError(f"Unknown genome filter field '{field}'.")
                    required_genome.add(field)
                    continue
                if prefix == "assembly":
                    if field not in assembly_cols:
                        raise ValueError(f"Unknown assembly filter field '{field}'.")
                    required_assembly.add(field)
                    continue
                if field in genome_cols:
                    required_genome.add(field)
                if field in assembly_cols:
                    required_assembly.add(field)
                if field not in genome_cols and field not in assembly_cols:
                    raise ValueError(f"Unknown filter field '{field}'.")

    genome_rows = _fetch_table_rows(
        manager,
        table="Genome",
        accessions=pool,
        columns=sorted(required_genome),
    )
    assembly_rows = _fetch_table_rows(
        manager,
        table="Assembly",
        accessions=pool,
        columns=sorted(required_assembly),
    )

    record_map: Dict[str, Dict[str, Any]] = {acc: {} for acc in pool}
    for acc, row in genome_rows.items():
        record = record_map.setdefault(acc, {})
        for col, value in row.items():
            record[f"genome.{col}"] = value
            if col not in record or record[col] is None:
                record[col] = value
    for acc, row in assembly_rows.items():
        record = record_map.setdefault(acc, {})
        for col, value in row.items():
            record[f"assembly.{col}"] = value
            if value is not None:
                record[col] = value
            elif col not in record:
                record[col] = value

    busco_map: Dict[str, Dict[str, Optional[float]]] = {}
    busco_profile_map: Dict[str, Dict[str, Optional[str]]] = {}
    if needs_busco:
        busco_map = _compute_best_busco_map(
            manager,
            pool,
            library_id=busco_library_id,
            library_name=busco_library_name,
        )
        busco_profile_map = _compute_busco_profile_map(
            manager,
            pool,
            library_id=busco_library_id,
            library_name=busco_library_name,
        )

    for groups in parsed_filters:
        if not pool:
            break
        filtered: List[str] = []
        for acc in pool:
            record = record_map.get(acc, {})
            def _matches_group(group: List[Dict[str, Any]]) -> bool:
                for cond in group:
                    if cond["prefix"] == "busco" and cond["field"] in custom_busco_fields:
                        value = busco_profile_map.get(acc, {}).get(cond["field"])
                    elif cond["prefix"] == "busco":
                        value = busco_map.get(acc, {}).get(cond["field"])
                    elif cond["field"] in custom_busco_fields:
                        value = busco_profile_map.get(acc, {}).get(cond["field"])
                    elif cond["prefix"] == "genome":
                        value = record.get(f"genome.{cond['field']}")
                    elif cond["prefix"] == "assembly":
                        value = record.get(f"assembly.{cond['field']}")
                    else:
                        value = record.get(cond["field"])
                    if cond["field"] in custom_busco_fields and cond["op"] in {"eq", "ne"} and cond.get("value") not in (None, ""):
                        matches = manager.proteomes.profile_matches_selector(
                            acc,
                            str(value).strip() or None,
                            str(cond["value"]).strip(),
                        )
                        if cond["op"] == "ne":
                            matches = not matches
                        if not matches:
                            return False
                        continue
                    if not _evaluate_filter_condition(value, cond):
                        return False
                return True

            if any(_matches_group(group) for group in groups):
                filtered.append(acc)
        pool = filtered

    return list(dict.fromkeys(pool))


def _split_tokens(value: Any) -> List[str]:
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


RANK_HIERARCHY = [
    "domain",
    "superkingdom",
    "kingdom",
    "subkingdom",
    "clade",
    "superphylum",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "cohort",
    "subcohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species",
    "subspecies",
]


def _normalize_rank_token(token: str) -> str:
    return str(token or "").strip().lower()


def normalize_rank_list(ranks: Sequence[Any]) -> List[str]:
    tokens = _split_tokens(ranks)
    normalized: List[str] = []
    seen: set[str] = set()
    for token in tokens:
        norm = _normalize_rank_token(token)
        if not norm:
            continue
        if norm not in RANK_HIERARCHY:
            raise ValueError(f"Unknown rank '{token}'.")
        if norm in seen:
            continue
        seen.add(norm)
        normalized.append(norm)
    if not normalized:
        return []
    order = {rank: idx for idx, rank in enumerate(RANK_HIERARCHY)}
    return sorted(normalized, key=lambda r: order.get(r, 0))


def _parse_quantity_tokens(values: Any) -> List[Optional[int]]:
    tokens = _split_tokens(values)
    parsed: List[Optional[int]] = []
    for token in tokens:
        key = str(token or "").strip().lower()
        if not key:
            continue
        if key in {"*", "all"}:
            parsed.append(None)
            continue
        try:
            qty = int(key)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid quantity '{token}'.") from exc
        if qty <= 0:
            raise ValueError("Quantities must be positive when provided.")
        parsed.append(qty)
    return parsed


def expand_accession_variables(
    manager: DBManager,
    accessions: Sequence[Any],
    *,
    allow_bare: bool = True,
) -> List[str]:
    """Expand accession tokens that reference environment variables (e.g. @VAR)."""
    normalized = normalize_accessions(accessions or [])
    if not normalized:
        return []

    explicit_vars: List[str] = []
    bare_candidates: List[str] = []
    for token in normalized:
        if token.startswith("@") and len(token) > 1:
            explicit_vars.append(token[1:])
        elif allow_bare:
            bare_candidates.append(token)

    env_lookup: Dict[str, Any] = {}
    if explicit_vars:
        env_lookup.update(manager.get_environment_variables(explicit_vars))
    if allow_bare and bare_candidates:
        env_lookup.update(manager.get_environment_variables(bare_candidates))

    expanded: List[str] = []
    for token in normalized:
        var_name: Optional[str] = None
        if token.startswith("@") and len(token) > 1:
            var_name = token[1:]
        elif allow_bare and token in env_lookup:
            var_name = token

        if var_name is None:
            expanded.append(token)
            continue

        value = env_lookup.get(var_name)
        if value is None:
            raise ValueError(f"Unknown accession variable '{var_name}'.")
        if isinstance(value, (list, tuple, set)):
            expanded.extend([str(item) for item in value if item is not None])
        elif isinstance(value, str):
            expanded.extend([part for part in value.split(",") if part.strip()])
        else:
            raise ValueError(f"Accession variable '{var_name}' must be a list or comma-separated string.")

    return normalize_accessions(expanded)


def resolve_clade_to_taxid(manager: DBManager, clade: str) -> Optional[int]:
    """Translate a clade/scientific name into a taxid using the taxonomy table."""
    token = str(clade or "").strip()
    if not token:
        return None
    manager.cursor.execute(
        "SELECT taxid FROM Taxonomy WHERE name = ? LIMIT 1",
        (token,),
    )
    row = manager.cursor.fetchone()
    if row:
        return int(row[0])
    manager.cursor.execute(
        "SELECT taxid FROM Taxonomy WHERE name = ? COLLATE NOCASE LIMIT 1",
        (token,),
    )
    row = manager.cursor.fetchone()
    return int(row[0]) if row else None


def _coerce_selector_request(selectors: SelectorRequest | Mapping[str, Any] | None) -> SelectorRequest:
    if isinstance(selectors, SelectorRequest):
        return selectors
    return SelectorRequest.from_mapping(selectors)


def prune_selector_mapping(mapping: Mapping[str, Any]) -> Dict[str, Any]:
    """Return a compact selector mapping suitable for persistence."""

    pruned: Dict[str, Any] = {}
    for key, value in dict(mapping or {}).items():
        if value is None:
            continue
        if value is False:
            continue
        if value == [] or value == {} or value == "":
            continue
        if key == "sample_strategy" and value == "rank":
            continue
        pruned[key] = value
    return pruned


def active_selector_overrides(mapping: Mapping[str, Any]) -> Dict[str, Any]:
    """Return user-supplied selector values that should override a preset."""

    overrides: Dict[str, Any] = {}
    for key, value in dict(mapping or {}).items():
        if key in {"preset", "preset_name", "description"}:
            continue
        if value is None:
            continue
        if value == [] or value == {} or value == "":
            continue
        if value is False:
            # False is meaningful for argparse parsers using SUPPRESS defaults
            # (for example --no-isoforms-cleaned), but not for store_true flags
            # that are always present. Keep only values explicitly present in
            # a sparse namespace/mapping; dense parsers are handled by callers.
            overrides[key] = value
            continue
        overrides[key] = value
    return overrides


def merge_selector_preset(
    manager: DBManager,
    preset_name: Any,
    overrides: Mapping[str, Any] | SelectorRequest | None = None,
) -> SelectorRequest:
    """Load a named selector preset and merge active overrides on top."""

    preset = manager.selector_presets.get(preset_name)
    if preset is None:
        raise ValueError(f"Unknown selector preset '{preset_name}'.")
    base = dict(preset.get("selector") or {})
    if isinstance(overrides, SelectorRequest):
        override_mapping = overrides.as_mapping()
    else:
        override_mapping = dict(overrides or {})
    base.update(active_selector_overrides(override_mapping))
    return SelectorRequest.from_mapping(base)


def _load_selector_env_defaults(manager: DBManager, scope: Optional[str] = None) -> Dict[str, Any]:
    defaults = {
        "downloaded_only": None,
        "not_downloaded": None,
        "primary_only": None,
        "use_busco": None,
        "status_min": None,
        "protein_only": None,
    }
    env_keys = {
        "downloaded_only": "SELECTOR_DEFAULT_DOWNLOADED_ONLY",
        "primary_only": "SELECTOR_DEFAULT_PRIMARY_ONLY",
        "use_busco": "SELECTOR_DEFAULT_USE_BUSCO",
        "status_min": "SELECTOR_DEFAULT_STATUS_MIN",
        "protein_only": "SELECTOR_DEFAULT_PROTEIN_ONLY",
    }

    scoped_env_keys: Dict[str, str] = {}
    if scope:
        scope_key = str(scope).strip().upper()
        for field, key in env_keys.items():
            scoped_env_keys[field] = f"SELECTOR_DEFAULT_{scope_key}_{key.split('SELECTOR_DEFAULT_', 1)[-1]}"

    keys = list(env_keys.values()) + list(scoped_env_keys.values())
    env = manager.get_environment_variables(keys) if keys else {}

    for field, key in env_keys.items():
        if key in env:
            defaults[field] = env.get(key)
    for field, key in scoped_env_keys.items():
        if key in env:
            defaults[field] = env.get(key)

    # Coerce known types
    defaults["downloaded_only"] = _coerce_bool(defaults.get("downloaded_only"))
    defaults["not_downloaded"] = _coerce_bool(defaults.get("not_downloaded"))
    defaults["primary_only"] = _coerce_bool(defaults.get("primary_only"))
    defaults["use_busco"] = _coerce_bool(defaults.get("use_busco"))
    defaults["protein_only"] = _coerce_bool(defaults.get("protein_only"))
    try:
        defaults["status_min"] = int(defaults["status_min"]) if defaults.get("status_min") is not None else None
    except (TypeError, ValueError):
        defaults["status_min"] = None

    return defaults


def _merge_selector_env_defaults(
    manager: DBManager,
    request: SelectorRequest,
    *,
    scope: Optional[str] = None,
) -> SelectorRequest:
    defaults = _load_selector_env_defaults(manager, scope=scope)
    overrides = {
        field: value
        for field, value in defaults.items()
        if getattr(request, field) is None and value is not None
    }
    if not overrides:
        return request
    return request.with_overrides(**overrides)


def _filter_accessions(
    manager: DBManager,
    accessions: Sequence[Any],
    *,
    root_id: Optional[int] = None,
    downloaded_only: bool,
    not_downloaded: bool = False,
    local_only: bool = False,
    not_local: bool = False,
    released_after: Optional[str],
    released_before: Optional[str],
    level: Optional[str],
    protein_only: bool = False,
    status_min: Optional[int] = None,
    primary_only: bool = False,
) -> List[str]:
    pool = normalize_accessions(accessions)
    if not pool:
        return []
    placeholders = ",".join("?" for _ in pool)
    params: List[Any] = list(pool)
    conditions = [f"g.accession IN ({placeholders})"]
    joins = " LEFT JOIN Hidden_Genomes h ON h.accession = g.accession"
    if root_id is not None:
        conditions.append("g.storage_root_id = ?")
        params.append(int(root_id))
    if downloaded_only and not_downloaded:
        raise ValueError("Use only one of downloaded_only or not_downloaded.")
    if local_only and not_local:
        raise ValueError("Use only one of local_only or not_local.")
    if status_min is not None and not downloaded_only and not not_downloaded:
        conditions.append("g.status >= ?")
        params.append(status_min)
    if downloaded_only:
        conditions.append("g.status > 0 AND g.location IS NOT NULL")
    if not_downloaded:
        conditions.append("(g.status <= 0 OR g.location IS NULL)")
    conditions.append("h.accession IS NULL")
    join_needed = bool(released_after or released_before or level or primary_only or local_only or not_local)
    if join_needed:
        joins += " LEFT JOIN Assembly a ON a.accession = g.accession"
        if released_after:
            conditions.append("a.release_date IS NOT NULL AND a.release_date >= ?")
            params.append(released_after)
        if released_before:
            conditions.append("a.release_date IS NOT NULL AND a.release_date <= ?")
            params.append(released_before)
        if local_only:
            conditions.append("LOWER(COALESCE(a.origin, '')) = 'local'")
        if not_local:
            conditions.append("LOWER(COALESCE(a.origin, '')) != 'local'")
        if primary_only:
            conditions.append(
                "COALESCE(LOWER(a.refseq_category), '') NOT LIKE '%alternate%' "
                "AND COALESCE(LOWER(a.assembly_type), '') NOT LIKE '%alternate%' "
                "AND COALESCE(LOWER(a.diploid_role), '') NOT LIKE '%alternate%'"
            )
    if level:
        conditions.append("LOWER(g.assembly_level) = ?")
        params.append(level.lower())
    if protein_only:
        conditions.append("g.protein = 1")
    sql = "SELECT DISTINCT g.accession FROM Genome g" + joins + " WHERE " + " AND ".join(conditions)
    manager.cursor.execute(sql, params)
    rows = manager.cursor.fetchall() or []
    return [str(row[0]) for row in rows]


def _select_accessions_by_filters(
    manager: DBManager,
    *,
    root_id: Optional[int] = None,
    downloaded_only: bool,
    not_downloaded: bool = False,
    local_only: bool = False,
    not_local: bool = False,
    released_after: Optional[str],
    released_before: Optional[str],
    level: Optional[str],
    protein_only: bool = False,
    status_min: Optional[int] = None,
    primary_only: bool = False,
) -> List[str]:
    params: List[Any] = []
    conditions = ["1=1"]
    joins = " LEFT JOIN Hidden_Genomes h ON h.accession = g.accession"
    if root_id is not None:
        conditions.append("g.storage_root_id = ?")
        params.append(int(root_id))
    if downloaded_only and not_downloaded:
        raise ValueError("Use only one of downloaded_only or not_downloaded.")
    if local_only and not_local:
        raise ValueError("Use only one of local_only or not_local.")
    if status_min is not None and not downloaded_only and not not_downloaded:
        conditions.append("g.status >= ?")
        params.append(status_min)
    if downloaded_only:
        conditions.append("g.status > 0 AND g.location IS NOT NULL")
    if not_downloaded:
        conditions.append("(g.status <= 0 OR g.location IS NULL)")
    conditions.append("h.accession IS NULL")
    join_needed = bool(released_after or released_before or level or primary_only or local_only or not_local)
    if join_needed:
        joins += " LEFT JOIN Assembly a ON a.accession = g.accession"
        if released_after:
            conditions.append("a.release_date IS NOT NULL AND a.release_date >= ?")
            params.append(released_after)
        if released_before:
            conditions.append("a.release_date IS NOT NULL AND a.release_date <= ?")
            params.append(released_before)
        if local_only:
            conditions.append("LOWER(COALESCE(a.origin, '')) = 'local'")
        if not_local:
            conditions.append("LOWER(COALESCE(a.origin, '')) != 'local'")
        if primary_only:
            conditions.append(
                "COALESCE(LOWER(a.refseq_category), '') NOT LIKE '%alternate%' "
                "AND COALESCE(LOWER(a.assembly_type), '') NOT LIKE '%alternate%' "
                "AND COALESCE(LOWER(a.diploid_role), '') NOT LIKE '%alternate%'"
            )
    if level:
        conditions.append("LOWER(g.assembly_level) = ?")
        params.append(level.lower())
    if protein_only:
        conditions.append("g.protein = 1")
    sql = "SELECT DISTINCT g.accession FROM Genome g" + joins + " WHERE " + " AND ".join(conditions)
    manager.cursor.execute(sql, params)
    rows = manager.cursor.fetchall() or []
    return [str(row[0]) for row in rows]


def _accessions_for_taxid(
    manager: DBManager,
    taxid: int,
    *,
    root_id: Optional[int] = None,
    downloaded_only: bool,
    not_downloaded: bool = False,
    local_only: bool = False,
    not_local: bool = False,
    released_after: Optional[str],
    released_before: Optional[str],
    level: Optional[str],
    protein_only: bool = False,
    status_min: Optional[int] = None,
    primary_only: bool = False,
) -> List[str]:
    base = manager.get_accessions_by_taxid(
        taxid,
        include_descendants=True,
        status_min=status_min,
        protein_only=protein_only,
    )
    if not base:
        return []
    seeds = [row[0] if isinstance(row, (list, tuple)) else row for row in base]
    return _filter_accessions(
        manager,
        seeds,
        root_id=root_id,
        downloaded_only=downloaded_only,
        not_downloaded=not_downloaded,
        local_only=local_only,
        not_local=not_local,
        released_after=released_after,
        released_before=released_before,
        level=level,
        protein_only=protein_only,
        status_min=status_min,
        primary_only=primary_only,
    )


def _ensure_list(payload: MutableMapping[str, Any], key: str) -> List[Any]:
    value = payload.get(key)
    if value is None:
        payload[key] = []
        return payload[key]  # type: ignore[return-value]
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        converted = list(value)
        payload[key] = converted
        return converted
    new_list = [value]
    payload[key] = new_list
    return new_list


def apply_selector_enrichment(
    manager: DBManager,
    spec: TaskSpec,
    payload: Mapping[str, Any],
    selector_source: SelectorRequest | Mapping[str, Any],
) -> MutableMapping[str, Any]:
    """Combine selectors with payload fields, returning a new enriched payload."""
    enriched: MutableMapping[str, Any] = dict(payload)
    request = _coerce_selector_request(selector_source)

    selector_config = {}
    if hasattr(spec, "metadata") and isinstance(spec.metadata, Mapping):
        selector_config = dict(spec.metadata.get("selectors", {}))
    auto_expand_taxid = selector_config.get("auto_expand_taxid_accessions", True)
    status_min = selector_config.get("status_min")
    try:
        status_min = int(status_min) if status_min is not None else None
    except (TypeError, ValueError):
        status_min = None

    downloaded_only = bool(request.downloaded_only)
    not_downloaded = bool(request.not_downloaded)
    local_only = bool(request.local_only)
    not_local = bool(request.not_local)
    root_id: Optional[int] = None
    if request.root is not None:
        row = manager.storage.resolve_root_token(request.root, kind="genomes")
        root_id = int(row[0])
    require_busco = bool(request.has_busco_results if request.has_busco_results is not None else selector_config.get("require_busco_results", False))
    missing_busco_results = bool(request.missing_busco_results)
    busco_complete_min = request.busco_complete_min
    busco_single_min = request.busco_single_min
    released_after = request.after
    released_before = request.before
    level = request.level
    filters = request.filters
    busco_library_id = request.busco_library_id
    busco_library_name = request.library_name
    if downloaded_only and not_downloaded:
        raise ValueError("Use only one of downloaded_only or not_downloaded.")
    if local_only and not_local:
        raise ValueError("Use only one of local_only or not_local.")
    filters_requested = bool(
        downloaded_only
        or not_downloaded
        or local_only
        or not_local
        or root_id is not None
        or released_after
        or released_before
        or level
        or status_min is not None
        or filters
    )
    protein_only_flag = bool(enriched.get("protein_only", False))
    primary_only_flag = bool(enriched.get("primary_only", False) or request.primary_only)

    # Expand/normalise existing accession selectors early so filtering and
    # queue-time validation operate on concrete accessions instead of @VAR tokens.
    def _expand_list_field(field_name: str) -> None:
        if field_name not in enriched:
            return
        enriched[field_name] = expand_accession_variables(
            manager,
            enriched.get(field_name) or [],
            allow_bare=True,
        )

    _expand_list_field("accessions")
    _expand_list_field("exclude_accessions")
    _expand_list_field("targets")
    _expand_list_field("ref_accessions")
    _expand_list_field("refs")

    # If selector thresholds are provided, map them to the payload keys expected by tasks
    try:
        if busco_complete_min is not None:
            enriched["min_completeness"] = float(busco_complete_min)
    except (TypeError, ValueError, OverflowError):
        pass
    try:
        if busco_single_min is not None:
            enriched["min_single_copy_complete"] = float(busco_single_min)
    except (TypeError, ValueError, OverflowError):
        pass

    payload_fields = spec.payload_model.model_fields
    # Tasks without accessions are treated as direct-input tasks; skip selector enrichment.
    if "accessions" not in payload_fields:
        return enriched

    def _extend_accessions(items: Iterable[str]) -> None:
        if "accessions" not in payload_fields:
            raise ValueError("This task does not accept accessions selectors.")
        store = _ensure_list(enriched, "accessions")
        store.extend(items)

    if request.accessions:
        _extend_accessions(
            expand_accession_variables(
                manager,
                request.accessions,
                allow_bare=True,
            )
        )

    clades = list(request.clade)

    for clade in clades:
        taxid = resolve_clade_to_taxid(manager, clade)
        if taxid is None:
            raise ValueError(f"Unknown clade '{clade}'.")
        if "taxid" in payload_fields and "taxid" not in enriched:
            enriched["taxid"] = taxid
        else:
            if not auto_expand_taxid:
                raise ValueError(
                    "This task does not support deriving accessions from clades when a taxid is already provided."
                )
            derived = _accessions_for_taxid(
                manager,
                taxid,
                downloaded_only=downloaded_only,
                root_id=root_id,
                not_downloaded=not_downloaded,
                local_only=local_only,
                not_local=not_local,
                released_after=released_after,
                released_before=released_before,
                level=level,
                protein_only=protein_only_flag,
                primary_only=primary_only_flag,
                status_min=status_min,
            )
            if filters:
                derived = filter_accessions_by_expressions(
                    manager,
                    derived,
                    filters,
                    busco_library_id=busco_library_id,
                    busco_library_name=busco_library_name,
                )
            if not derived:
                raise ValueError(f"No accessions resolved for clade '{clade}'.")
            _extend_accessions(derived)

    explicit_taxid = request.taxid
    if explicit_taxid is not None:
        taxid_val = int(explicit_taxid)
        if "taxid" in payload_fields:
            existing_taxid = enriched.get("taxid")
            if existing_taxid is not None:
                try:
                    existing_val = int(existing_taxid)
                except (TypeError, ValueError):
                    existing_val = None
                if existing_val is not None and existing_val != taxid_val:
                    raise ValueError(
                        f"Conflicting taxid selectors: payload/clade resolved to '{existing_val}', but --taxid is '{taxid_val}'."
                    )
            enriched["taxid"] = taxid_val
            filters_requested = True
        else:
            derived = _accessions_for_taxid(
                manager,
                taxid_val,
                downloaded_only=downloaded_only,
                root_id=root_id,
                not_downloaded=not_downloaded,
                local_only=local_only,
                not_local=not_local,
                released_after=released_after,
                released_before=released_before,
                level=level,
                protein_only=protein_only_flag,
                primary_only=primary_only_flag,
                status_min=status_min,
            )
            if filters:
                derived = filter_accessions_by_expressions(
                    manager,
                    derived,
                    filters,
                    busco_library_id=busco_library_id,
                    busco_library_name=busco_library_name,
                )
            if not derived:
                raise ValueError(f"No accessions resolved for taxid '{taxid_val}'.")
            _extend_accessions(derived)
            filters_requested = True

    if (
        auto_expand_taxid
        and "taxid" in enriched
        and "accessions" in payload_fields
        and "accessions" not in enriched
    ):
        derived = _accessions_for_taxid(
            manager,
            int(enriched["taxid"]),
            downloaded_only=downloaded_only,
            root_id=root_id,
            not_downloaded=not_downloaded,
            local_only=local_only,
            not_local=not_local,
            released_after=released_after,
            released_before=released_before,
            level=level,
            protein_only=protein_only_flag,
            primary_only=primary_only_flag,
            status_min=status_min,
        )
        if filters and derived:
            derived = filter_accessions_by_expressions(
                manager,
                derived,
                filters,
                busco_library_id=busco_library_id,
                busco_library_name=busco_library_name,
            )
        if derived:
            _extend_accessions(derived)
            filters_requested = True

    if require_busco and missing_busco_results:
        raise ValueError("Use only one of has_busco_results or missing_busco_results.")

    if "accessions" in enriched:
        unique = list(dict.fromkeys(normalize_accessions(enriched["accessions"])))
        requested_accessions = list(unique)
        if require_busco or missing_busco_results:
            # Determine library to scope BUSCO results: prefer explicit library_id/name, else any BUSCO result.
            busco_lib_id = None
            lib_id = enriched.get("library_id")
            if not lib_id and enriched.get("library_name"):
                lib_id = manager.get_library_id(enriched["library_name"])
            if lib_id:
                parent = manager.assert_library_has_parent(lib_id)
                busco_lib_id = parent if parent else lib_id
            if busco_lib_id is None:
                allowed = set(manager.get_busco_processed_accessions_any())
            else:
                allowed = set(manager.get_busco_processed_accessions(busco_lib_id))
            if require_busco:
                unique = [acc for acc in unique if acc in allowed]
            else:
                unique = [acc for acc in unique if acc not in allowed]
        if require_busco and not unique:
            raise ValueError("No accessions remain after requiring BUSCO results for the selected library.")
        if missing_busco_results and not unique:
            raise ValueError("No accessions remain after filtering to missing BUSCO results.")
        if filters_requested:
            filtered = _filter_accessions(
                manager,
                unique,
                downloaded_only=downloaded_only,
                root_id=root_id,
                not_downloaded=not_downloaded,
                local_only=local_only,
                not_local=not_local,
                released_after=released_after,
                released_before=released_before,
                level=level,
                protein_only=protein_only_flag,
                status_min=status_min,
                primary_only=primary_only_flag,
            )
            if filters:
                filtered = filter_accessions_by_expressions(
                    manager,
                    filtered,
                    filters,
                    busco_library_id=busco_library_id,
                    busco_library_name=busco_library_name,
            )
            if not filtered:
                raise ValueError("No accessions matched the provided filters.")
            unique = filtered
            enriched["accessions"] = filtered
        else:
            if unique:
                enriched["accessions"] = unique
            else:
                enriched.pop("accessions", None)
        skipped_accessions = [acc for acc in requested_accessions if acc not in set(unique)]
        if skipped_accessions:
            enriched["selector_requested_accessions"] = requested_accessions
            enriched["selector_skipped_accessions"] = skipped_accessions

    return enriched


def resolve_selector_candidates(
    manager: DBManager,
    selectors: SelectorRequest | Mapping[str, Any],
    *,
    allow_all: bool = False,
    require_candidates: bool = True,
    scope: Optional[str] = None,
    allow_env_defaults: bool = True,
    allow_bare_variables: bool = True,
) -> List[str]:
    """Resolve selectors into a candidate list (filters + exclusions, no BUSCO or rule selection)."""
    request = _coerce_selector_request(selectors)
    if allow_env_defaults:
        request = _merge_selector_env_defaults(manager, request, scope=scope)

    downloaded_only = bool(request.downloaded_only)
    not_downloaded = bool(request.not_downloaded)
    local_only = bool(request.local_only)
    not_local = bool(request.not_local)
    released_after = request.after
    released_before = request.before
    level = request.level
    primary_only = bool(request.primary_only)
    protein_only = bool(request.protein_only)
    status_min = request.status_min
    try:
        status_min = int(status_min) if status_min is not None else None
    except (TypeError, ValueError):
        status_min = None
    filters = request.filters
    busco_library_id = request.busco_library_id
    busco_library_name = request.library_name
    busco_run_ids = expand_busco_run_id_variables(manager, request.busco_run_ids)
    root_id: Optional[int] = None
    if request.root is not None:
        row = manager.storage.resolve_root_token(request.root, kind="genomes")
        root_id = int(row[0])

    accessions = expand_accession_variables(
        manager,
        request.accessions,
        allow_bare=allow_bare_variables,
    )
    exclude_accessions = expand_accession_variables(
        manager,
        request.exclude_accessions,
        allow_bare=allow_bare_variables,
    )

    taxid_values: List[int] = []
    if request.taxid is not None:
        try:
            taxid_values.append(int(request.taxid))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid taxid '{request.taxid}'.") from exc

    for clade in request.clade:
        taxid = resolve_clade_to_taxid(manager, clade)
        if taxid is None:
            raise ValueError(f"Unknown clade '{clade}'.")
        taxid_values.append(int(taxid))

    # Base pool
    pool: List[str] = list(accessions)
    for taxid in taxid_values:
        rows = manager.get_accessions_by_taxid(
            int(taxid),
            include_descendants=True,
            status_min=status_min,
            protein_only=protein_only,
        )
        pool.extend(row[0] if isinstance(row, (list, tuple)) else row for row in rows or [])

    busco_scope = None
    if busco_library_id is not None:
        busco_scope = manager.assert_library_has_parent(busco_library_id) or busco_library_id
    busco_filters_requested = request.has_busco_filters()
    if not pool and busco_filters_requested:
        if busco_run_ids:
            pool = manager.busco.get_accessions_for_run_ids(busco_run_ids, library_id=busco_scope)
        elif busco_scope is not None:
            pool = manager.get_busco_processed_accessions(busco_scope)
        else:
            pool = manager.get_busco_processed_accessions_any()

    if downloaded_only and not_downloaded:
        raise ValueError("Use only one of --downloaded-only or --not-downloaded.")
    if local_only and not_local:
        raise ValueError("Use only one of --local-only or --not-local.")
    filters_requested = bool(
        downloaded_only
        or not_downloaded
        or local_only
        or not_local
        or root_id is not None
        or released_after
        or released_before
        or level
        or primary_only
        or status_min is not None
        or protein_only
    )

    if not pool and (allow_all or request.allow_all):
        pool = _select_accessions_by_filters(
            manager,
            root_id=root_id,
            downloaded_only=downloaded_only,
            not_downloaded=not_downloaded,
            local_only=local_only,
            not_local=not_local,
            released_after=released_after,
            released_before=released_before,
            level=level,
            protein_only=protein_only,
            status_min=status_min,
            primary_only=primary_only,
        )
    elif pool:
        pool = _filter_accessions(
            manager,
            pool,
            root_id=root_id,
            downloaded_only=downloaded_only,
            not_downloaded=not_downloaded,
            local_only=local_only,
            not_local=not_local,
            released_after=released_after,
            released_before=released_before,
            level=level,
            protein_only=protein_only,
            status_min=status_min,
            primary_only=primary_only,
        )

    if filters:
        pool = filter_accessions_by_expressions(
            manager,
            pool,
            filters,
            busco_library_id=busco_library_id,
            busco_library_name=busco_library_name,
        )

    # Exclusions
    exclusion_set = set(exclude_accessions)
    exclude_taxids = _split_tokens(request.exclude_taxids)
    exclude_clades = _split_tokens(request.exclude_clades)
    for token in exclude_taxids:
        try:
            taxid_values = [int(token)]
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid exclude taxid '{token}'.") from exc
        for taxid in taxid_values:
            rows = manager.get_accessions_by_taxid(
                int(taxid),
                include_descendants=True,
                status_min=status_min,
                protein_only=protein_only,
            )
            derived = [row[0] if isinstance(row, (list, tuple)) else row for row in rows or []]
            if derived:
                filtered = _filter_accessions(
                    manager,
                    derived,
                    root_id=root_id,
                    downloaded_only=downloaded_only,
                    not_downloaded=not_downloaded,
                    local_only=local_only,
                    not_local=not_local,
                    released_after=released_after,
                    released_before=released_before,
                    level=level,
                    protein_only=protein_only,
                    status_min=status_min,
                    primary_only=primary_only,
                )
                if filters:
                    filtered = filter_accessions_by_expressions(
                        manager,
                        filtered,
                        filters,
                        busco_library_id=busco_library_id,
                        busco_library_name=busco_library_name,
                    )
                exclusion_set.update(filtered)
    for clade in exclude_clades:
        taxid = resolve_clade_to_taxid(manager, clade)
        if taxid is None:
            raise ValueError(f"Unknown exclude clade '{clade}'.")
        rows = manager.get_accessions_by_taxid(
            int(taxid),
            include_descendants=True,
            status_min=status_min,
            protein_only=protein_only,
        )
        derived = [row[0] if isinstance(row, (list, tuple)) else row for row in rows or []]
        if derived:
            filtered = _filter_accessions(
                manager,
                derived,
                root_id=root_id,
                downloaded_only=downloaded_only,
                not_downloaded=not_downloaded,
                local_only=local_only,
                not_local=not_local,
                released_after=released_after,
                released_before=released_before,
                level=level,
                protein_only=protein_only,
                status_min=status_min,
                primary_only=primary_only,
            )
            if filters:
                filtered = filter_accessions_by_expressions(
                    manager,
                    filtered,
                    filters,
                    busco_library_id=busco_library_id,
                    busco_library_name=busco_library_name,
                )
            exclusion_set.update(filtered)

    if exclusion_set:
        pool = [acc for acc in pool if acc not in exclusion_set]

    pool = list(dict.fromkeys(normalize_accessions(pool)))
    if not pool and require_candidates:
        raise ValueError("No accessions matched the provided selectors.")
    return pool


def resolve_selector_accessions(
    manager: DBManager,
    selectors: SelectorRequest | Mapping[str, Any],
    *,
    allow_all: bool = False,
    require_candidates: bool = True,
    use_rule_selection: bool = True,
    scope: Optional[str] = None,
    allow_env_defaults: bool = True,
    allow_bare_variables: bool = True,
) -> List[str]:
    """Resolve selectors into a final accession list (BUSCO filters + rule selection)."""
    request = _coerce_selector_request(selectors)
    if allow_env_defaults:
        request = _merge_selector_env_defaults(manager, request, scope=scope)

    candidates = resolve_selector_candidates(
        manager,
        request,
        allow_all=allow_all or request.allow_all,
        require_candidates=require_candidates,
        scope=scope,
        allow_env_defaults=False,
        allow_bare_variables=allow_bare_variables,
    )

    has_busco_results = bool(request.has_busco_results)
    missing_busco_results = bool(request.missing_busco_results)
    busco_library_id = request.busco_library_id
    busco_complete_min = request.busco_complete_min
    busco_single_min = request.busco_single_min
    ranks_raw = request.ranks
    quantities_raw = request.quantities
    sample_strategy = request.sample_strategy or "rank"
    sample_seed = request.sample_seed

    include_paralog_filtering_in_score = request.include_paralog_filtering_in_score
    include_decontamination_in_score = request.include_decontamination_in_score
    paralog_run_id = request.paralog_run_id
    decontamination_run_id = request.decontamination_run_id
    allow_ambiguous_contaminants = request.allow_ambiguous_contaminants
    strict_decontamination = request.strict_decontamination
    rescue_duplicates = request.rescue_duplicates
    busco_pipeline = request.busco_pipeline
    busco_input_mode = request.busco_input_mode
    prefer_busco_pipeline = request.prefer_busco_pipeline
    prefer_busco_input_mode = request.prefer_busco_input_mode
    busco_export_format = request.busco_export_format
    busco_run_ids = request.busco_run_ids
    busco_run_selection = request.busco_run_selection
    paralog_filtered = bool(request.paralog_filtered)
    not_paralog_filtered = bool(request.not_paralog_filtered)
    min_hidden_paralogs = request.min_hidden_paralogs
    max_hidden_paralogs = request.max_hidden_paralogs
    decontaminated = bool(request.decontaminated)
    not_decontaminated = bool(request.not_decontaminated)
    contaminated = bool(request.contaminated)
    decontamination_run = request.decontamination_run
    ignore_contaminated_assemblies = request.ignore_contaminated_assemblies

    candidates = filter_accessions_by_busco_selectors(
        manager,
        candidates,
        has_busco_results=has_busco_results,
        missing_busco_results=missing_busco_results,
        busco_library_id=busco_library_id,
        busco_complete_min=busco_complete_min,
        busco_single_min=busco_single_min,
        include_paralog_filtering_in_score=include_paralog_filtering_in_score,
        include_decontamination_in_score=include_decontamination_in_score,
        paralog_run_id=paralog_run_id,
        decontamination_run_id=decontamination_run_id,
        allow_ambiguous_contaminants=allow_ambiguous_contaminants,
        strict_decontamination=strict_decontamination,
        rescue_duplicates=rescue_duplicates,
        busco_pipeline=busco_pipeline,
        busco_input_mode=busco_input_mode,
        prefer_busco_pipeline=prefer_busco_pipeline,
        prefer_busco_input_mode=prefer_busco_input_mode,
        busco_export_format=busco_export_format,
        busco_run_ids=busco_run_ids,
        busco_run_selection=busco_run_selection,
        paralog_filtered=paralog_filtered,
        not_paralog_filtered=not_paralog_filtered,
        min_hidden_paralogs=min_hidden_paralogs,
        max_hidden_paralogs=max_hidden_paralogs,
        decontaminated=decontaminated,
        not_decontaminated=not_decontaminated,
        contaminated=contaminated,
        decontamination_run=decontamination_run,
        ignore_contaminated_assemblies=ignore_contaminated_assemblies,
    )

    use_busco = request.use_busco
    if use_busco is None:
        use_busco = False
    use_busco = bool(use_busco)

    rule_quantity = request.quantity
    rule_rank = request.rank

    explicit_accessions = expand_accession_variables(
        manager,
        request.accessions,
        allow_bare=allow_bare_variables,
    )
    has_explicit_accessions = bool(explicit_accessions)

    taxid = None if has_explicit_accessions else request.taxid

    if rule_rank is not None and rule_quantity is not None and taxid is None:
        clade_tokens = list(request.clade)
        if len(clade_tokens) == 1 and not has_explicit_accessions:
            resolved = resolve_clade_to_taxid(manager, clade_tokens[0])
            if resolved is None:
                raise ValueError(f"Unknown clade '{clade_tokens[0]}'.")
            taxid = int(resolved)

    def _as_pct(value: Optional[float]) -> Optional[float]:
        if value is None:
            return None
        try:
            val = float(value)
        except (TypeError, ValueError, OverflowError):
            return None
        return val * 100.0 if val <= 1.0 else val

    min_complete_pct = _as_pct(busco_complete_min)
    min_single_pct = _as_pct(busco_single_min)

    if use_rule_selection:
        if quantities_raw or ranks_raw:
            if rule_rank is not None or rule_quantity is not None:
                raise ValueError("Use either --rank/--quantity or --ranks/--quantities, not both.")
            if not quantities_raw:
                # Ranks provided without quantities are treated as display-only (no subsampling).
                return candidates
            rank_tokens = _split_tokens(ranks_raw)
            qty_tokens = _parse_quantity_tokens(quantities_raw)
            if not rank_tokens:
                raise ValueError("--ranks requires at least one rank token.")
            if not qty_tokens:
                raise ValueError("--quantities requires at least one quantity token.")
            if len(rank_tokens) != len(qty_tokens):
                raise ValueError("--ranks and --quantities must have the same number of entries.")
            normalized_ranks: List[str] = []
            seen_ranks: set[str] = set()
            for token in rank_tokens:
                norm = _normalize_rank_token(token)
                if not norm:
                    continue
                if norm not in RANK_HIERARCHY:
                    raise ValueError(f"Unknown rank '{token}'.")
                if norm in seen_ranks:
                    raise ValueError(f"Duplicate rank '{token}' in --ranks.")
                seen_ranks.add(norm)
                normalized_ranks.append(norm)
            if len(normalized_ranks) != len(qty_tokens):
                raise ValueError("--ranks and --quantities must have the same number of entries.")

            if taxid is None and not has_explicit_accessions:
                clade_tokens = list(request.clade)
                if len(clade_tokens) == 1:
                    resolved = resolve_clade_to_taxid(manager, clade_tokens[0])
                    if resolved is None:
                        raise ValueError(f"Unknown clade '{clade_tokens[0]}'.")
                    taxid = int(resolved)

            candidates = apply_multistage_selection(
                manager,
                candidates,
                taxid=taxid,
                ranks=normalized_ranks,
                quantities=qty_tokens,
                sample_strategy=sample_strategy,
                sample_seed=sample_seed,
                busco_library_id=busco_library_id,
                use_busco=use_busco,
                min_completeness=min_complete_pct,
                min_single_copy_complete=min_single_pct,
                allow_duplicate_species=bool(request.allow_duplicate_species),
                include_paralog_filtering_in_score=include_paralog_filtering_in_score,
                include_decontamination_in_score=include_decontamination_in_score,
                decontamination_run_id=decontamination_run_id,
                allow_ambiguous_contaminants=allow_ambiguous_contaminants,
                strict_decontamination=strict_decontamination,
                rescue_duplicates=rescue_duplicates,
            )
        else:
            candidates = apply_rule_selection(
                manager,
                candidates,
                taxid=taxid,
                rule_quantity=rule_quantity,
                rule_rank=rule_rank,
                busco_library_id=busco_library_id,
                use_busco=use_busco,
                min_completeness=min_complete_pct,
                min_single_copy_complete=min_single_pct,
                allow_duplicate_species=bool(request.allow_duplicate_species),
                include_paralog_filtering_in_score=include_paralog_filtering_in_score,
                include_decontamination_in_score=include_decontamination_in_score,
                decontamination_run_id=decontamination_run_id,
                allow_ambiguous_contaminants=allow_ambiguous_contaminants,
                strict_decontamination=strict_decontamination,
                rescue_duplicates=rescue_duplicates,
            )

    selected = list(dict.fromkeys(normalize_accessions(candidates)))
    if not selected and require_candidates:
        raise ValueError("No accessions matched the provided selectors.")
    return selected


def collect_accession_candidates(
    manager: DBManager,
    explicit: Sequence[str],
    *,
    taxid: Optional[int],
    downloaded_only: bool = False,
    not_downloaded: bool = False,
    local_only: bool = False,
    not_local: bool = False,
    released_after: Optional[str] = None,
    released_before: Optional[str] = None,
    level: Optional[str] = None,
    protein_only: bool = False,
    status_min: Optional[int] = None,
    primary_only: bool = False,
    filters: Any = None,
    busco_library_id: Optional[int] = None,
    busco_library_name: Optional[str] = None,
) -> List[str]:
    """Return a sorted, unique list of candidate accessions for downstream processing."""

    pool = set(normalize_accessions(explicit))
    if taxid is not None:
        derived = manager.get_accessions_by_taxid(
            int(taxid),
            include_descendants=True,
            status_min=status_min,
            protein_only=protein_only,
        )
        pool.update(row[0] for row in derived or [])

    if not pool:
        return []

    filtered = _filter_accessions(
        manager,
        list(pool),
        downloaded_only=downloaded_only,
        not_downloaded=not_downloaded,
        local_only=local_only,
        not_local=not_local,
        released_after=released_after,
        released_before=released_before,
        level=level,
        protein_only=protein_only,
        status_min=status_min,
        primary_only=primary_only,
    )
    if filters:
        filtered = filter_accessions_by_expressions(
            manager,
            filtered,
            filters,
            busco_library_id=busco_library_id,
            busco_library_name=busco_library_name,
        )
    # Preserve deterministic ordering for reproducibility
    return list(dict.fromkeys(sorted(filtered)))


def filter_accessions_by_busco_selectors(
    manager: DBManager,
    candidates: Sequence[str],
    *,
    has_busco_results: bool = False,
    missing_busco_results: bool = False,
    busco_library_id: Optional[int] = None,
    busco_complete_min: Optional[float] = None,
    busco_single_min: Optional[float] = None,
    include_paralog_filtering_in_score: Optional[bool] = None,
    include_decontamination_in_score: Optional[bool] = None,
    paralog_run_id: Optional[str] = None,
    decontamination_run_id: Optional[str] = None,
    allow_ambiguous_contaminants: Optional[bool] = None,
    strict_decontamination: Optional[bool] = None,
    rescue_duplicates: Optional[bool] = None,
    busco_pipeline: Optional[str] = None,
    busco_input_mode: Optional[str] = None,
    prefer_busco_pipeline: Optional[str] = None,
    prefer_busco_input_mode: Optional[str] = None,
    proteome_profile: Optional[str] = None,
    prefer_proteome_profile: Optional[str] = None,
    busco_export_format: Optional[str] = None,
    busco_run_ids: Optional[Sequence[Any]] = None,
    busco_run_selection: Optional[str] = None,
    paralog_filtered: bool = False,
    not_paralog_filtered: bool = False,
    min_hidden_paralogs: Optional[float] = None,
    max_hidden_paralogs: Optional[float] = None,
    decontaminated: bool = False,
    not_decontaminated: bool = False,
    contaminated: bool = False,
    decontamination_run: Optional[str] = None,
    ignore_contaminated_assemblies: Optional[bool] = None,
) -> List[str]:
    """Apply BUSCO/paralog/decontamination selector filters to a candidate list."""

    base = list(dict.fromkeys(normalize_accessions(candidates)))
    if not base:
        return []

    def _as_pct(value: Optional[float]) -> Optional[float]:
        if value is None:
            return None
        try:
            val = float(value)
        except (TypeError, ValueError) as exc:
            from .errors import SelectorError

            raise SelectorError("BUSCO selector thresholds must be numeric.") from exc
        if val < 0:
            raise ValueError("BUSCO selector thresholds must be non-negative.")
        return val * 100.0 if val <= 1.0 else val

    complete_thresh = _as_pct(busco_complete_min)
    single_thresh = _as_pct(busco_single_min)

    parent_id = None
    if busco_library_id is not None:
        parent_id = manager.assert_library_has_parent(busco_library_id)

    paralog_required = any(
        [
            include_paralog_filtering_in_score is not None,
            paralog_run_id is not None,
            paralog_filtered,
            not_paralog_filtered,
            min_hidden_paralogs is not None,
            max_hidden_paralogs is not None,
        ]
    )
    decont_required = any(
        [
            include_decontamination_in_score is not None,
            decontamination_run_id is not None,
            allow_ambiguous_contaminants is not None,
            strict_decontamination is not None,
            decontaminated,
            not_decontaminated,
            contaminated,
            decontamination_run is not None,
            ignore_contaminated_assemblies is not None,
        ]
    )

    if paralog_required:
        if busco_library_id is None:
            raise ValueError("Paralog selectors require --library-id or --library-name.")
    if decont_required and busco_library_id is None:
        raise ValueError("Decontamination selectors require --library-id or --library-name.")

    if has_busco_results and missing_busco_results:
        raise ValueError("Use only one of --has-busco-results or --missing-busco-results.")
    if missing_busco_results and (busco_complete_min is not None or busco_single_min is not None):
        raise ValueError("BUSCO threshold selectors cannot be combined with --missing-busco-results.")

    run_query_requested = any(
        (
            busco_run_ids,
            busco_pipeline is not None,
            busco_input_mode is not None,
            busco_run_selection is not None,
            busco_export_format is not None,
        )
    )
    if busco_library_id is not None and base and run_query_requested:
        busco_scope = parent_id if parent_id is not None else busco_library_id
        run_map = manager.busco._resolve_busco_runs_for_query(
            int(busco_scope),
            accessions=base,
            run_ids=expand_busco_run_id_variables(manager, busco_run_ids or []),
            pipeline=busco_pipeline,
            input_mode=busco_input_mode,
            preferred_pipeline=prefer_busco_pipeline,
            preferred_input_mode=prefer_busco_input_mode,
            proteome_profile=proteome_profile,
            preferred_proteome_profile=prefer_proteome_profile,
            selection=busco_run_selection or "primary",
            purpose="default",
        )
        if busco_export_format:
            purpose = "export_nucleotide" if str(busco_export_format).strip().lower() == "nucleotide" else "export_protein"
            run_map = {
                acc: run_id
                for acc, run_id in run_map.items()
                if manager.busco.run_supports_purpose(int(run_id), purpose=purpose)
            }
        allowed_accessions = set(run_map.keys())
        base = [acc for acc in base if acc in allowed_accessions]

    if has_busco_results or missing_busco_results:
        busco_scope = busco_library_id
        if parent_id:
            busco_scope = parent_id
        if busco_scope is None:
            allowed = set(manager.get_busco_processed_accessions_any())
        else:
            allowed = set(manager.get_busco_processed_accessions(busco_scope))
        if has_busco_results:
            base = [acc for acc in base if acc in allowed]
        else:
            base = [acc for acc in base if acc not in allowed]

    if complete_thresh is not None or single_thresh is not None:
        allowed: set[str] = set()
        if busco_library_id is None:
            rows = manager.get_busco_results_percentages(accessions=base)
            for row in rows or []:
                acc = str(row[0])
                complete_pct = row[3]
                single_pct = row[4]
                ok = True
                if complete_thresh is not None:
                    ok = ok and (complete_pct is not None and complete_pct >= complete_thresh)
                if single_thresh is not None:
                    ok = ok and (single_pct is not None and single_pct >= single_thresh)
                if ok:
                    allowed.add(acc)
        elif (
            include_paralog_filtering_in_score is not None
            or parent_id is not None
            or paralog_run_id is not None
            or include_decontamination_in_score is not False
            or decontamination_run_id is not None
            or allow_ambiguous_contaminants is not None
            or strict_decontamination is not None
        ):
            rows = manager.get_busco_results_adjusted(
                library_id=busco_library_id,
                accessions=base,
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
                complete_pct = row[3]
                single_pct = row[4]
                ok = True
                if complete_thresh is not None:
                    ok = ok and (complete_pct is not None and complete_pct >= complete_thresh)
                if single_thresh is not None:
                    ok = ok and (single_pct is not None and single_pct >= single_thresh)
                if ok:
                    allowed.add(acc)
        else:
            rows = manager.get_busco_results_percentages(accessions=base, library_id=busco_library_id)
            for row in rows or []:
                acc = str(row[0])
                complete_pct = row[3]
                single_pct = row[4]
                ok = True
                if complete_thresh is not None:
                    ok = ok and (complete_pct is not None and complete_pct >= complete_thresh)
                if single_thresh is not None:
                    ok = ok and (single_pct is not None and single_pct >= single_thresh)
                if ok:
                    allowed.add(acc)
        base = [acc for acc in base if acc in allowed]

    if not base:
        return []

    if busco_library_id is not None:
        if paralog_filtered and not_paralog_filtered:
            raise ValueError("Use only one of --paralog-filtered or --not-paralog-filtered.")
        if decontaminated and not_decontaminated:
            raise ValueError("Use only one of --decontaminated or --not-decontaminated.")
        if contaminated and (not_decontaminated or decontaminated):
            raise ValueError("Use only one of --contaminated, --decontaminated, or --not-decontaminated.")
        if decontamination_run and not_decontaminated:
            raise ValueError("--decontamination-run cannot be combined with --not-decontaminated.")

        if parent_id is not None and (paralog_filtered or not_paralog_filtered):
            have_paralog = manager.get_paralog_filtering_accessions(
                target_library_id=busco_library_id,
                busco_library_id=parent_id,
                accessions=base,
                run_id=paralog_run_id,
            )
            if paralog_filtered:
                base = [acc for acc in base if acc in have_paralog]
            elif not_paralog_filtered:
                base = [acc for acc in base if acc not in have_paralog]

        if parent_id is not None and (min_hidden_paralogs is not None or max_hidden_paralogs is not None):
            size = manager._get_library_size(busco_library_id)
            if size <= 0:
                base = []
            else:
                hidden_counts = manager._paralog_hidden_counts(
                    target_library_id=busco_library_id,
                    busco_library_id=parent_id,
                    accessions=base,
                    run_id=paralog_run_id,
                )
                min_hidden_pct = _as_pct(min_hidden_paralogs)
                max_hidden_pct = _as_pct(max_hidden_paralogs)
                filtered: List[str] = []
                for acc in base:
                    counts = hidden_counts.get(acc)
                    if counts is None:
                        continue
                    hidden_val = counts[0]
                    pct = (float(hidden_val) / size) * 100.0
                    if min_hidden_pct is not None and pct < min_hidden_pct:
                        continue
                    if max_hidden_pct is not None and pct > max_hidden_pct:
                        continue
                    filtered.append(acc)
                base = filtered

        if decontaminated or not_decontaminated or contaminated or decontamination_run:
            have_decont = manager.get_decontamination_accessions(
                target_library_id=busco_library_id,
                run_id=decontamination_run,
                accessions=base,
            )
            if parent_id:
                have_decont.update(
                    manager.get_decontamination_accessions(
                        target_library_id=parent_id,
                        run_id=decontamination_run,
                        accessions=base,
                    )
                )
            if decontaminated or decontamination_run:
                base = [acc for acc in base if acc in have_decont]
            elif not_decontaminated:
                base = [acc for acc in base if acc not in have_decont]

        if contaminated:
            run_id = decontamination_run or decontamination_run_id
            summaries = manager.get_latest_decontamination_summary_with_fallback(
                target_library_id=busco_library_id,
                parent_library_id=parent_id,
                accessions=base,
                run_id=run_id,
            )
            contaminated_only = {
                acc for acc, (_run, decision, _date) in summaries.items()
                if str(decision or "").strip().lower() == "contaminated"
            }
            base = [acc for acc in base if acc in contaminated_only]

        if contaminated:
            ignore_contaminated_assemblies = False
        elif ignore_contaminated_assemblies is None:
            ignore_contaminated_assemblies = False
        if ignore_contaminated_assemblies:
            run_id = decontamination_run or decontamination_run_id
            summaries = manager.get_latest_decontamination_summary_with_fallback(
                target_library_id=busco_library_id,
                parent_library_id=parent_id,
                accessions=base,
                run_id=run_id,
            )
            contaminated = {
                acc for acc, (_run, decision, _date) in summaries.items()
                if str(decision or "").strip().lower() == "contaminated"
            }
            if contaminated:
                base = [acc for acc in base if acc not in contaminated]

    return base


def _build_scoring_context(
    manager: DBManager,
    base: Sequence[str],
    *,
    busco_library_id: Optional[int],
    use_busco: bool,
    min_completeness: Optional[float],
    min_single_copy_complete: Optional[float],
    include_paralog_filtering_in_score: Optional[bool],
    include_decontamination_in_score: Optional[bool],
    paralog_run_id: Optional[str],
    decontamination_run_id: Optional[str],
    allow_ambiguous_contaminants: Optional[bool],
    strict_decontamination: Optional[bool],
    rescue_duplicates: Optional[bool],
) -> tuple[List[str], callable, callable, Dict[str, Dict[str, object]]]:
    base_list = list(dict.fromkeys(normalize_accessions(base)))
    if not base_list:
        return [], lambda _a: (), lambda _a: "", {}

    def _parse_order(val: Any) -> List[str]:
        if val is None:
            return []
        if isinstance(val, (list, tuple)):
            return [str(x) for x in val if x]
        try:
            import json

            parsed = json.loads(val)
            if isinstance(parsed, (list, tuple)):
                return [str(x) for x in parsed if x]
        except (TypeError, json.JSONDecodeError):
            pass
        return [token.strip() for token in str(val).split(",") if token.strip()]

    def _parse_buckets(val: Any) -> List[tuple]:
        if val is None:
            return []
        try:
            import json

            parsed = json.loads(val)
            if isinstance(parsed, list):
                cleaned = []
                for entry in parsed:
                    if isinstance(entry, (list, tuple)) and len(entry) == 2:
                        try:
                            cleaned.append((float(entry[0]), float(entry[1])))
                        except (TypeError, ValueError):
                            continue
                return cleaned
        except (TypeError, json.JSONDecodeError):
            pass
        pairs = []
        for token in str(val).split(","):
            if ":" not in token:
                continue
            left, right = token.split(":", 1)
            try:
                pairs.append((float(left), float(right)))
            except (TypeError, ValueError):
                continue
        return pairs

    env_cfg = manager.get_environment_variables(["SELECTOR_SCORE_ORDER", "SELECTOR_BUSCO_BUCKETS"])

    score_order_env = _parse_order(env_cfg.get("SELECTOR_SCORE_ORDER"))
    default_order = ["busco", "refseq", "level", "n50", "date", "accession"]
    score_order = score_order_env or default_order

    buckets_env = _parse_buckets(env_cfg.get("SELECTOR_BUSCO_BUCKETS"))
    default_buckets = [
        (0, -5),
        (10, -4),
        (20, -3),
        (30, -2),
        (40, -1),
        (50, 1),
        (60, 2),
        (70, 3),
        (80, 4),
        (90, 5),
    ]
    busco_buckets = buckets_env or default_buckets

    placeholders = ",".join("?" for _ in base_list)
    manager.cursor.execute(
        f"""
        SELECT g.accession, g.taxid, g.assembly_level, a.release_date, a.contig_n50, a.origin, t.name
        FROM Genome g
        LEFT JOIN Assembly a ON a.accession = g.accession
        LEFT JOIN Taxonomy t ON t.taxid = g.taxid
        WHERE g.accession IN ({placeholders})
        """,
        tuple(base_list),
    )
    acc_meta_rows = manager.cursor.fetchall() or []
    acc_meta: Dict[str, Dict[str, object]] = {}
    for acc, tid, level, release_date, contig_n50, origin, species in acc_meta_rows:
        acc_meta[str(acc)] = {
            "taxid": int(tid) if tid is not None else None,
            "assembly_level": (level or "").lower(),
            "release_date": release_date,
            "contig_n50": contig_n50,
            "origin": (origin or "").lower(),
            "species": str(species) if species else None,
        }

    libraries: Dict[str, Dict[str, object]] = {}
    manager.cursor.execute("SELECT library_id, library_name, taxid FROM Libraries")
    for lid, lname, ltaxid in manager.cursor.fetchall() or []:
        libraries[str(lname)] = {"library_id": int(lid), "taxid": int(ltaxid) if ltaxid is not None else None}

    parent_id = None
    if busco_library_id is not None:
        parent_id = manager.assert_library_has_parent(busco_library_id)
    if (
        include_paralog_filtering_in_score is not None
        or paralog_run_id is not None
        or include_decontamination_in_score is not None
        or decontamination_run_id is not None
        or allow_ambiguous_contaminants is not None
        or strict_decontamination is not None
    ):
        if busco_library_id is None:
            raise ValueError("Paralog/decontamination scoring requires --library-id or --library-name.")

    busco_rows: Dict[str, list] = {}
    use_adjusted = False
    if busco_library_id is not None:
        use_adjusted = (
            include_paralog_filtering_in_score is not None
            or
            bool(parent_id)
            or paralog_run_id is not None
            or include_decontamination_in_score is not False
            or decontamination_run_id is not None
            or allow_ambiguous_contaminants is not None
            or strict_decontamination is not None
        )
    if busco_library_id is not None and use_adjusted:
        rows = manager.get_busco_results_adjusted(
            library_id=busco_library_id,
            accessions=base_list,
            include_paralog=include_paralog_filtering_in_score,
            paralog_run_id=paralog_run_id,
            include_decontam=include_decontamination_in_score,
            decont_run_id=decontamination_run_id,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )
        for row in rows or []:
            acc = row[0]
            lib_name = row[2]
            pct_c = row[3]
            pct_sc = row[4]
            busco_rows.setdefault(str(acc), []).append((str(lib_name), pct_c, pct_sc))
    else:
        rows = manager.get_busco_results_percentages(accessions=base_list, library_id=busco_library_id)
        for row in rows or []:
            acc = row[0]
            lib_name = row[2]
            pct_c = row[3]
            pct_sc = row[4]
            busco_rows.setdefault(str(acc), []).append((str(lib_name), pct_c, pct_sc))

    lineage_cache: Dict[int, Dict[int, int]] = {}

    def _lineage_depth_map(tid: int) -> Dict[int, int]:
        if tid in lineage_cache:
            return lineage_cache[tid]
        depth_map: Dict[int, int] = {}
        rows = manager.get_lineage_root_to_leaf(tid) or []
        for depth, row in enumerate(rows):
            try:
                depth_map[int(row[0])] = depth
            except (TypeError, ValueError, IndexError):
                continue
        lineage_cache[tid] = depth_map
        return depth_map

    def _best_busco_scores(acc: str) -> tuple[Optional[float], Optional[float]]:
        rows = busco_rows.get(acc) or []
        meta = acc_meta.get(acc, {})
        acc_tid = meta.get("taxid")

        preferred: tuple[Optional[float], Optional[float]] | None = None
        preferred_depth = -1
        if acc_tid is not None:
            depth_map = _lineage_depth_map(int(acc_tid))
            for lib_name, pct_c, pct_sc in rows:
                lib_info = libraries.get(lib_name, {})
                lib_taxid = lib_info.get("taxid")
                lib_id = lib_info.get("library_id")
                if busco_library_id is not None and lib_id == busco_library_id:
                    return (
                        float(pct_c or 0.0) if pct_c is not None else None,
                        float(pct_sc or 0.0) if pct_sc is not None else None,
                    )
                if lib_taxid is None:
                    continue
                if lib_taxid in depth_map and depth_map[lib_taxid] > preferred_depth:
                    preferred = (
                        float(pct_c or 0.0) if pct_c is not None else None,
                        float(pct_sc or 0.0) if pct_sc is not None else None,
                    )
                    preferred_depth = depth_map[lib_taxid]
        if preferred is not None:
            return preferred
        best: tuple[Optional[float], Optional[float]] | None = None
        for _lib, pct_c, pct_sc in rows:
            cand = (
                float(pct_c or 0.0) if pct_c is not None else None,
                float(pct_sc or 0.0) if pct_sc is not None else None,
            )
            if best is None:
                best = cand
                continue
            sc_best = best[1] if best[1] is not None else -1.0
            sc_cand = cand[1] if cand[1] is not None else -1.0
            c_best = best[0] if best[0] is not None else -1.0
            c_cand = cand[0] if cand[0] is not None else -1.0
            if sc_cand > sc_best or (sc_cand == sc_best and c_cand > c_best):
                best = cand
        return best or (None, None)

    def _busco_bucket(sc_value: Optional[float]) -> int:
        if sc_value is None:
            return -5
        try:
            v = float(sc_value)
        except (TypeError, ValueError, OverflowError):
            return -5
        bucket_score = -5
        for threshold, score_val in sorted(busco_buckets, key=lambda x: x[0]):
            if v >= threshold:
                bucket_score = int(score_val)
        return bucket_score

    busco_bucket_map: Dict[str, Optional[int]] = {}
    if use_busco:
        for acc in base_list:
            comp, sc_comp = _best_busco_scores(acc)
            if min_completeness is not None and (comp is None or comp < float(min_completeness)):
                busco_bucket_map[acc] = None
                continue
            if min_single_copy_complete is not None and (sc_comp is None or sc_comp < float(min_single_copy_complete)):
                busco_bucket_map[acc] = None
                continue
            busco_bucket_map[acc] = _busco_bucket(sc_comp if sc_comp is not None else comp)

    level_priority = {
        "complete genome": 4,
        "chromosome": 3,
        "scaffold": 2,
        "contig": 1,
    }
    origin_priority = {
        "refseq": 1,
    }

    def score(accession: str) -> tuple:
        meta = acc_meta.get(accession, {})
        busco_score = busco_bucket_map.get(accession, 0) if use_busco else 0
        lvl = meta.get("assembly_level") or ""
        level_score = level_priority.get(str(lvl).lower(), 0)
        origin = meta.get("origin") or ""
        refcat_score = origin_priority.get(str(origin).lower(), 0)
        try:
            n50 = float(meta.get("contig_n50") or 0)
        except (TypeError, ValueError, OverflowError):
            n50 = 0.0
        rd = meta.get("release_date")
        try:
            rd_val = datetime.strptime(str(rd).split(" ")[0], "%Y-%m-%d")
        except (TypeError, ValueError):
            rd_val = datetime.min
        components = {
            "busco": busco_score,
            "refseq": refcat_score,
            "level": level_score,
            "n50": n50,
            "date": rd_val,
            "accession": accession,
        }
        return tuple(components[k] for k in score_order if k in components)

    eligible = [acc for acc in base_list if not (use_busco and busco_bucket_map.get(acc) is None)]

    def species_key(accession: str) -> str:
        meta = acc_meta.get(accession, {})
        label = meta.get("species")
        if label:
            return str(label).strip().lower()
        return f"__acc__:{accession}"

    return eligible, score, species_key, acc_meta


def _select_top_ranked(
    entries: Sequence[str],
    *,
    score,
    species_key,
    quantity: Optional[int],
    allow_duplicate_species: bool,
    strategy: str = "rank",
    rng=None,
) -> List[str]:
    if strategy == "random":
        ordered = list(entries)
        if rng is None:
            raise ValueError("Random selection requires an RNG instance.")
        rng.shuffle(ordered)
    else:
        ordered = sorted(entries, key=score, reverse=True)
    if quantity is None:
        return ordered
    if allow_duplicate_species:
        return ordered[:quantity]
    selected: List[str] = []
    seen: set[str] = set()
    for accession in ordered:
        key = species_key(accession)
        if key in seen:
            continue
        seen.add(key)
        selected.append(accession)
        if len(selected) >= quantity:
            break
    return selected


def _build_accession_rank_map(
    manager: DBManager,
    eligible: Sequence[str],
    acc_meta: Mapping[str, Dict[str, object]],
) -> Dict[str, Dict[str, tuple]]:
    lineage_cache: Dict[int, Dict[str, tuple]] = {}

    def _rank_map(tid: int) -> Dict[str, tuple]:
        cached = lineage_cache.get(tid)
        if cached is not None:
            return cached
        mapping: Dict[str, tuple] = {}
        rows = manager.get_lineage_root_to_leaf(tid) or []
        for tid_val, name, rank, _parent in rows:
            token = (rank or "").lower()
            if token:
                mapping[token] = (tid_val, name or str(tid_val))
        lineage_cache[tid] = mapping
        return mapping

    acc_rank_map: Dict[str, Dict[str, tuple]] = {}
    for accession in eligible:
        tid = acc_meta.get(accession, {}).get("taxid")
        if tid is None:
            continue
        acc_rank_map[accession] = _rank_map(int(tid))
    return acc_rank_map


def apply_rule_selection(
    manager: DBManager,
    candidates: Sequence[str],
    *,
    taxid: Optional[int],
    rule_quantity: Optional[int],
    rule_rank: Optional[str],
    busco_library_id: Optional[int] = None,
    use_busco: bool = True,
    min_completeness: Optional[float] = None,
    min_single_copy_complete: Optional[float] = None,
    allow_duplicate_species: bool = False,
    include_paralog_filtering_in_score: Optional[bool] = None,
    include_decontamination_in_score: Optional[bool] = None,
    paralog_run_id: Optional[str] = None,
    decontamination_run_id: Optional[str] = None,
    allow_ambiguous_contaminants: Optional[bool] = None,
    strict_decontamination: Optional[bool] = None,
    rescue_duplicates: Optional[bool] = None,
) -> List[str]:
    """Apply rule-based selection to a candidate list, returning the selected accessions.

    If rule_rank is provided, selection is top-N per rank across the candidate list.
    When taxid is also provided, it is used only to validate that the requested rank
    sits below the supplied anchor taxon.
    If rule_rank is omitted, selection is top-N across the full candidate list.
    Ranking considers (in order): BUSCO (optional, decile bucketed on single-copy completeness),
    RefSeq origin, assembly level, contig N50, release date, then accession as a deterministic
    tiebreak."""

    base = list(dict.fromkeys(normalize_accessions(candidates)))
    if not base:
        return []

    if rule_quantity is None:
        return base

    try:
        top_n = int(rule_quantity)
    except (TypeError, ValueError) as exc:
        raise ValueError("Rule quantity must be an integer.") from exc
    if top_n <= 0:
        raise ValueError("Rule quantity must be positive.")

    rank_token = None
    if rule_rank is not None:
        rank_token = str(rule_rank).strip().lower()
        if not rank_token:
            raise ValueError("rule_rank must be a non-empty string when quantity is provided.")
        if taxid is not None:
            lineage = manager.get_lineage_root_to_leaf(int(taxid)) or []
            ancestor_ranks = {(row[2] or "").lower() for row in lineage}
            if rank_token in ancestor_ranks:
                raise ValueError("rule_rank is at or above the provided taxid's rank; cannot partition in such a way.")

    eligible, score, species_key, acc_meta = _build_scoring_context(
        manager,
        base,
        busco_library_id=busco_library_id,
        use_busco=use_busco,
        min_completeness=min_completeness,
        min_single_copy_complete=min_single_copy_complete,
        include_paralog_filtering_in_score=include_paralog_filtering_in_score,
        include_decontamination_in_score=include_decontamination_in_score,
        paralog_run_id=paralog_run_id,
        decontamination_run_id=decontamination_run_id,
        allow_ambiguous_contaminants=allow_ambiguous_contaminants,
        strict_decontamination=strict_decontamination,
        rescue_duplicates=rescue_duplicates,
    )
    if not eligible:
        raise ValueError("Rule-based selection yielded no accessions.")

    if rank_token:
        acc_rank_map = _build_accession_rank_map(manager, eligible, acc_meta)
        grouped: Dict[str, List[str]] = {}
        group_order: List[str] = []
        for acc in eligible:
            info = acc_rank_map.get(acc, {}).get(rank_token)
            if not info:
                continue
            group_name = str(info[1])
            if group_name not in grouped:
                grouped[group_name] = []
                group_order.append(group_name)
            grouped[group_name].append(acc)

        if not grouped:
            raise ValueError("Rule-based selection yielded no accessions.")

        selected: List[str] = []
        for group_name in group_order:
            entries = grouped[group_name]
            selected.extend(
                _select_top_ranked(
                    entries,
                    score=score,
                    species_key=species_key,
                    quantity=top_n,
                    allow_duplicate_species=allow_duplicate_species,
                )
            )
        if not selected:
            raise ValueError("Rule-based selection yielded no accessions.")
        return list(dict.fromkeys(selected))

    selected = _select_top_ranked(
        eligible,
        score=score,
        species_key=species_key,
        quantity=top_n,
        allow_duplicate_species=allow_duplicate_species,
    )
    if not selected:
        raise ValueError("Rule-based selection yielded no accessions.")
    return list(dict.fromkeys(selected))


def apply_multistage_selection(
    manager: DBManager,
    candidates: Sequence[str],
    *,
    taxid: Optional[int],
    ranks: Sequence[str],
    quantities: Sequence[Optional[int]],
    sample_strategy: str,
    sample_seed: Optional[int],
    busco_library_id: Optional[int] = None,
    use_busco: bool = True,
    min_completeness: Optional[float] = None,
    min_single_copy_complete: Optional[float] = None,
    allow_duplicate_species: bool = False,
    include_paralog_filtering_in_score: Optional[bool] = None,
    include_decontamination_in_score: Optional[bool] = None,
    paralog_run_id: Optional[str] = None,
    decontamination_run_id: Optional[str] = None,
    allow_ambiguous_contaminants: Optional[bool] = None,
    strict_decontamination: Optional[bool] = None,
    rescue_duplicates: Optional[bool] = None,
) -> List[str]:
    base = list(dict.fromkeys(normalize_accessions(candidates)))
    if not base:
        return []

    if not ranks:
        return base
    if len(ranks) != len(quantities):
        raise ValueError("ranks and quantities must have the same length.")

    strategy = (sample_strategy or "rank").strip().lower()
    if strategy not in {"rank", "random"}:
        raise ValueError("sample_strategy must be 'rank' or 'random'.")

    eligible, score, species_key, acc_meta = _build_scoring_context(
        manager,
        base,
        busco_library_id=busco_library_id,
        use_busco=use_busco,
        min_completeness=min_completeness,
        min_single_copy_complete=min_single_copy_complete,
        include_paralog_filtering_in_score=include_paralog_filtering_in_score,
        include_decontamination_in_score=include_decontamination_in_score,
        paralog_run_id=paralog_run_id,
        decontamination_run_id=decontamination_run_id,
        allow_ambiguous_contaminants=allow_ambiguous_contaminants,
        strict_decontamination=strict_decontamination,
        rescue_duplicates=rescue_duplicates,
    )
    if not eligible:
        return []

    order = {rank: idx for idx, rank in enumerate(RANK_HIERARCHY)}
    ordered_ranks = sorted(list(dict.fromkeys(ranks)), key=lambda r: order.get(r, 0))
    rank_to_qty = {rank: qty for rank, qty in zip(ranks, quantities)}
    ordered_quantities = [rank_to_qty.get(rank) for rank in ordered_ranks]

    if taxid is not None:
        lineage = manager.get_lineage_root_to_leaf(int(taxid)) or []
        ancestor_ranks = {(row[2] or "").lower() for row in lineage}
        for rank in ordered_ranks:
            if rank in ancestor_ranks:
                raise ValueError("Rule-based selection by rank requires a taxid below the requested ranks.")

    acc_rank_map = _build_accession_rank_map(manager, eligible, acc_meta)

    rng = None
    if strategy == "random":
        import random

        rng = random.Random(sample_seed)

    lowest_rank = ordered_ranks[-1]
    lowest_qty = ordered_quantities[-1]
    grouped: Dict[str, List[str]] = {}
    group_order: List[str] = []
    for acc in eligible:
        rank_map = acc_rank_map.get(acc)
        if not rank_map:
            continue
        info = rank_map.get(lowest_rank)
        if not info:
            continue
        group_id = str(info[0])
        if group_id not in grouped:
            grouped[group_id] = []
            group_order.append(group_id)
        grouped[group_id].append(acc)

    units: List[Dict[str, object]] = []
    for group_id in group_order:
        members = _select_top_ranked(
            grouped[group_id],
            score=score,
            species_key=species_key,
            quantity=lowest_qty,
            allow_duplicate_species=allow_duplicate_species,
            strategy=strategy,
            rng=rng,
        )
        if not members:
            continue
        unit_score = max((score(acc) for acc in members), default=())
        units.append(
            {
                "rank": lowest_rank,
                "group_id": group_id,
                "members": members,
                "score": unit_score,
            }
        )

    if not units:
        raise ValueError("Rule-based selection yielded no accessions.")

    higher_ranks = list(reversed(ordered_ranks[:-1]))
    higher_quantities = list(reversed(ordered_quantities[:-1]))

    for rank_token, qty in zip(higher_ranks, higher_quantities):
        if qty is None:
            continue
        grouped_units: Dict[str, List[Dict[str, object]]] = {}
        order_units: List[str] = []
        for unit in units:
            members = unit.get("members") or []
            if not members:
                continue
            acc = str(members[0])
            rank_map = acc_rank_map.get(acc, {})
            info = rank_map.get(rank_token)
            if not info:
                continue
            parent_id = str(info[0])
            if parent_id not in grouped_units:
                grouped_units[parent_id] = []
                order_units.append(parent_id)
            grouped_units[parent_id].append(unit)

        next_units: List[Dict[str, object]] = []
        for parent_id in order_units:
            entries = grouped_units[parent_id]
            if strategy == "rank":
                entries = sorted(entries, key=lambda u: u.get("score") or (), reverse=True)
            else:
                entries = list(entries)
                rng.shuffle(entries)  # type: ignore[union-attr]
            next_units.extend(entries[:qty])
        units = next_units
        if not units:
            break

    if not units:
        raise ValueError("Rule-based selection yielded no accessions.")

    selected: List[str] = []
    seen_acc: set[str] = set()
    for unit in units:
        for acc in unit.get("members") or []:
            if acc in seen_acc:
                continue
            seen_acc.add(acc)
            selected.append(acc)

    if not selected:
        raise ValueError("Rule-based selection yielded no accessions.")
    return selected


__all__ = [
    "SelectorRequest",
    "apply_selector_enrichment",
    "active_selector_overrides",
    "merge_selector_preset",
    "normalize_accessions",
    "normalize_rank_list",
    "prune_selector_mapping",
    "expand_accession_variables",
    "resolve_clade_to_taxid",
    "collect_accession_candidates",
    "resolve_selector_candidates",
    "resolve_selector_accessions",
    "validate_filter_expressions",
    "filter_accessions_by_expressions",
    "filter_accessions_by_busco_selectors",
    "apply_rule_selection",
    "apply_multistage_selection",
]
