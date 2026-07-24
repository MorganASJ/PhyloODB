from abc import ABC, abstractmethod
from datetime import datetime
import json
import re
import time
import traceback
from typing import Any, List, Optional, Sequence

from ..database import DBManager
from ..errors import BatchFailure, BatchItemError, TaskExecutionError
from ..logging_utils import configure_logging_from_db, get_task_logger

_TASK_ACRONYMS = {
    "busco": "BUSCO",
    "iqtree": "IQ-TREE",
    "mafft": "MAFFT",
    "ncbi": "NCBI",
    "blast": "BLAST",
    "cdhit": "CD-HIT",
    "orthofinder": "OrthoFinder",
}


def _humanize_task_name(raw: str) -> str:
    token = str(raw or "").strip()
    if not token:
        return "Task"
    hyphenated = re.sub(r"([a-z0-9])([A-Z])", r"\1-\2", token).replace("_", "-")
    parts = [part for part in hyphenated.split("-") if part]
    if len(parts) == 2 and parts[1].lower() == "run":
        base = _TASK_ACRONYMS.get(parts[0].lower(), parts[0].capitalize())
        return f"Run {base}"
    words = []
    for part in parts:
        lower = part.lower()
        words.append(_TASK_ACRONYMS.get(lower, part.capitalize()))
    return " ".join(words) or "Task"


class Task(ABC):
    """A class that includes the basic structure and methods of a task"""
    @classmethod
    def default_thread_count(cls, registry_required_threads: int, daemon_max_threads: int) -> int:
        """Return the default required threads for this task under a daemon cap."""

        registry_required_threads = max(int(registry_required_threads or 1), 1)
        daemon_max_threads = max(int(daemon_max_threads or 1), 1)
        return min(registry_required_threads, daemon_max_threads)

    def __init__(self, db_path, task_id, data, required_threads=1):
        self.db_manager = DBManager(db_path)
        self.db_manager.connect()
        self.db_manager.validate_schema()
        self.task_id = task_id
        self.REQUIRED_THREADS = required_threads
        self.default_log_category = "TASK"
        self.status = "P"  # P, R, C, S, B, E
        if data:
            self.data = json.loads(data)
        else:
            self.data = {}
        self.task_display_name = self._resolve_task_display_name()
        # Initialize a task-scoped logger adapter (task_id may be None for daemon before queueing)
        self.logger = get_task_logger(task_id=self.task_id, task_name=self.task_display_name)
        self.stage = 0
        self._batch_failures: list[BatchFailure] = []

    @staticmethod
    def _coerce_bool_value(value: Any, default: bool = False) -> bool:
        if value is None:
            return bool(default)
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        token = str(value).strip().lower()
        if token in {"1", "true", "yes", "y", "on"}:
            return True
        if token in {"0", "false", "no", "n", "off", "none", "null", ""}:
            return False
        return bool(default)

    def env_bool(self, key: str, default: bool = False) -> bool:
        return self._coerce_bool_value(self.db_manager.get_environment_variable(key), default)

    def env_float(self, key: str, default: Optional[float] = None) -> Optional[float]:
        value = self.db_manager.get_environment_variable(key)
        if value is None or str(value).strip() == "":
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    def env_int(self, key: str, default: int) -> int:
        value = self.db_manager.get_environment_variable(key)
        if value is None or str(value).strip() == "":
            return int(default)
        try:
            return int(value)
        except (TypeError, ValueError):
            return int(default)

    def env_str(self, key: str, default: str) -> str:
        value = self.db_manager.get_environment_variable(key)
        token = str(value).strip() if value is not None else ""
        return token or default

    def payload_bool(self, key: str, default: bool = False) -> bool:
        if key in self.data:
            return self._coerce_bool_value(self.data.get(key), default)
        return bool(default)

    def _resolve_task_display_name(self) -> str:
        if self.task_id is None:
            return self.__class__.__name__
        row = self.db_manager.tasks.get(self.task_id)
        if row:
            job_type = row[1]
            if job_type == 0:
                return "Scheduler"
            try:
                from ..registry import registry

                spec = registry.get_by_job_type(job_type)
            except KeyError:
                spec = None
            if spec is not None:
                if spec.display_name:
                    return str(spec.display_name)
                key = str(spec.key or "").strip()
                if key:
                    return _humanize_task_name(key)
        return self.__class__.__name__

    def update_status(self, new_status):
        self.status = new_status
        self.db_manager.tasks.update_status(self.task_id, new_status)

    def selector_accessions(self):
        """Return a deduplicated list of accessions resolved by CLI selectors."""
        from ..selector_utils import expand_accession_variables

        raw = self.data.get("accessions") or []
        return list(dict.fromkeys(expand_accession_variables(self.db_manager, raw)))

    def resolve_assembly_accession(self, accession):
        """Resolve a stored RefSeq/GenBank alias and log canonicalization."""

        requested = str(accession or "").strip()
        resolved = self.db_manager.genomes.resolve_accession(requested)
        if requested and resolved and requested != resolved:
            self.log(
                f"Resolved requested accession {requested} to canonical assembly {resolved}.",
                "INFO",
            )
        return resolved

    def _build_selector_request(
        self,
        *,
        root: Optional[str] = None,
        taxid: Optional[int] = None,
        additional: Optional[Sequence[str]] = None,
        allow_all: bool = False,
        downloaded_only: Optional[bool] = None,
        not_downloaded: Optional[bool] = None,
        released_after: Optional[str] = None,
        released_before: Optional[str] = None,
        level: Optional[str] = None,
        protein_only: Optional[bool] = None,
        status_min: Optional[int] = None,
        primary_only: Optional[bool] = None,
        rule_quantity: Optional[int] = None,
        rule_rank: Optional[str] = None,
        busco_library_id: Optional[int] = None,
        use_busco: Optional[bool] = None,
        min_completeness: Optional[float] = None,
        min_single_copy_complete: Optional[float] = None,
        allow_duplicate_species: Optional[bool] = None,
        include_paralog_filtering_in_score: Optional[bool] = None,
        include_decontamination_in_score: Optional[bool] = None,
        paralog_run_id: Optional[str] = None,
        decontamination_run_id: Optional[str] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
        paralog_filtered: Optional[bool] = None,
        not_paralog_filtered: Optional[bool] = None,
        min_hidden_paralogs: Optional[float] = None,
        max_hidden_paralogs: Optional[float] = None,
        decontaminated: Optional[bool] = None,
        not_decontaminated: Optional[bool] = None,
        contaminated: Optional[bool] = None,
        decontamination_run: Optional[str] = None,
        ignore_contaminated_assemblies: Optional[bool] = None,
    ):
        """Build the canonical selector request for task-side resolution."""

        from ..selector_utils import SelectorRequest, expand_accession_variables

        request = SelectorRequest.from_mapping(self.data)
        accessions = list(self.selector_accessions())
        if additional:
            accessions.extend(expand_accession_variables(self.db_manager, additional))

        overrides = {
            "accessions": accessions,
            "root": root if root is not None else request.root,
            "taxid": taxid if taxid is not None else request.taxid,
            "allow_all": allow_all or request.allow_all,
        }
        if downloaded_only is not None:
            overrides["downloaded_only"] = downloaded_only
        if not_downloaded is not None:
            overrides["not_downloaded"] = not_downloaded
        if released_after is not None:
            overrides["after"] = released_after
        if released_before is not None:
            overrides["before"] = released_before
        if level is not None:
            overrides["level"] = level
        if protein_only is not None:
            overrides["protein_only"] = protein_only
        if status_min is not None:
            overrides["status_min"] = status_min
        if primary_only is not None:
            overrides["primary_only"] = primary_only
        if rule_quantity is not None:
            overrides["quantity"] = rule_quantity
        if rule_rank is not None:
            overrides["rank"] = rule_rank
        if busco_library_id is not None:
            overrides["busco_library_id"] = busco_library_id
        if use_busco is not None:
            overrides["use_busco"] = use_busco
        if min_completeness is not None:
            overrides["busco_complete_min"] = min_completeness
        if min_single_copy_complete is not None:
            overrides["busco_single_min"] = min_single_copy_complete
        if allow_duplicate_species is not None:
            overrides["allow_duplicate_species"] = allow_duplicate_species
        if include_paralog_filtering_in_score is not None:
            overrides["include_paralog_filtering_in_score"] = include_paralog_filtering_in_score
        if paralog_run_id is not None:
            overrides["paralog_run_id"] = paralog_run_id
        if include_decontamination_in_score is not None:
            overrides["include_decontamination_in_score"] = include_decontamination_in_score
        if decontamination_run_id is not None:
            overrides["decontamination_run_id"] = decontamination_run_id
        if allow_ambiguous_contaminants is not None:
            overrides["allow_ambiguous_contaminants"] = allow_ambiguous_contaminants
        if strict_decontamination is not None:
            overrides["strict_decontamination"] = strict_decontamination
        if rescue_duplicates is not None:
            overrides["rescue_duplicates"] = rescue_duplicates
        if paralog_filtered is not None:
            overrides["paralog_filtered"] = paralog_filtered
        if not_paralog_filtered is not None:
            overrides["not_paralog_filtered"] = not_paralog_filtered
        if min_hidden_paralogs is not None:
            overrides["min_hidden_paralogs"] = min_hidden_paralogs
        if max_hidden_paralogs is not None:
            overrides["max_hidden_paralogs"] = max_hidden_paralogs
        if decontaminated is not None:
            overrides["decontaminated"] = decontaminated
        if not_decontaminated is not None:
            overrides["not_decontaminated"] = not_decontaminated
        if contaminated is not None:
            overrides["contaminated"] = contaminated
        if decontamination_run is not None:
            overrides["decontamination_run"] = decontamination_run
        if ignore_contaminated_assemblies is not None:
            overrides["ignore_contaminated_assemblies"] = ignore_contaminated_assemblies
        return request.with_overrides(**overrides)

    def selector_candidates(
        self,
        *,
        root: Optional[str] = None,
        taxid: Optional[int] = None,
        additional: Optional[Sequence[str]] = None,
        allow_all: bool = False,
        downloaded_only: bool = False,
        not_downloaded: bool = False,
        released_after: Optional[str] = None,
        released_before: Optional[str] = None,
        level: Optional[str] = None,
        protein_only: bool = False,
        status_min: Optional[int] = None,
        primary_only: bool = False,
    ) -> List[str]:
        """Combine selector accessions with optional taxid-derived entries and filters."""
        from ..selector_utils import resolve_selector_candidates

        selectors = self._build_selector_request(
            root=root,
            taxid=taxid,
            additional=additional,
            allow_all=allow_all,
            downloaded_only=downloaded_only,
            not_downloaded=not_downloaded,
            released_after=released_after,
            released_before=released_before,
            level=level,
            protein_only=protein_only,
            status_min=status_min,
            primary_only=primary_only,
        )
        return resolve_selector_candidates(
            self.db_manager,
            selectors,
            allow_all=allow_all or bool(self.data.get("all", False)),
            require_candidates=False,
        )

    def apply_selector_rules(
        self,
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
        """Apply rule-based selection (quantity/rank) to a candidate list."""
        from ..selector_utils import apply_rule_selection

        return apply_rule_selection(
            self.db_manager,
            candidates,
            taxid=taxid,
            rule_quantity=rule_quantity,
            rule_rank=rule_rank,
            busco_library_id=busco_library_id,
            use_busco=use_busco,
            min_completeness=min_completeness,
            min_single_copy_complete=min_single_copy_complete,
            allow_duplicate_species=allow_duplicate_species,
            include_paralog_filtering_in_score=include_paralog_filtering_in_score,
            include_decontamination_in_score=include_decontamination_in_score,
            paralog_run_id=paralog_run_id,
            decontamination_run_id=decontamination_run_id,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
        )

    def resolve_selector_accessions(
        self,
        *,
        root: Optional[str] = None,
        taxid: Optional[int] = None,
        additional: Optional[Sequence[str]] = None,
        allow_all: bool = False,
        rule_quantity: Optional[int] = None,
        rule_rank: Optional[str] = None,
        busco_library_id: Optional[int] = None,
        downloaded_only: bool = False,
        not_downloaded: bool = False,
        released_after: Optional[str] = None,
        released_before: Optional[str] = None,
        level: Optional[str] = None,
        protein_only: bool = False,
        status_min: Optional[int] = None,
        primary_only: bool = False,
        require_candidates: bool = True,
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
        """Resolve selectors into a final accession list, applying filters and rules."""
        from ..selector_utils import resolve_selector_accessions

        selectors = self._build_selector_request(
            root=root,
            taxid=taxid,
            additional=additional,
            allow_all=allow_all,
            downloaded_only=downloaded_only,
            not_downloaded=not_downloaded,
            released_after=released_after,
            released_before=released_before,
            level=level,
            protein_only=protein_only,
            status_min=status_min,
            primary_only=primary_only,
            rule_quantity=rule_quantity,
            rule_rank=rule_rank,
            busco_library_id=busco_library_id,
            use_busco=use_busco,
            min_completeness=min_completeness,
            min_single_copy_complete=min_single_copy_complete,
            allow_duplicate_species=allow_duplicate_species,
            include_paralog_filtering_in_score=include_paralog_filtering_in_score,
            include_decontamination_in_score=include_decontamination_in_score,
            paralog_run_id=paralog_run_id,
            decontamination_run_id=decontamination_run_id,
            allow_ambiguous_contaminants=allow_ambiguous_contaminants,
            strict_decontamination=strict_decontamination,
            rescue_duplicates=rescue_duplicates,
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

        selected = resolve_selector_accessions(
            self.db_manager,
            selectors,
            allow_all=allow_all or bool(self.data.get("all", False)),
            require_candidates=require_candidates,
            use_rule_selection=True,
        )
        if not selected and require_candidates:
            raise ValueError("Rule-based selection yielded no accessions.")
        return list(dict.fromkeys(selected))

    def prepare_selectors(
        self,
        *,
        target_key: str = "_selector_accessions",
        root: Optional[str] = None,
        taxid: Optional[int] = None,
        additional: Optional[Sequence[str]] = None,
        rule_quantity: Optional[int] = None,
        rule_rank: Optional[str] = None,
        busco_library_id: Optional[int] = None,
        downloaded_only: Optional[bool] = None,
        not_downloaded: Optional[bool] = None,
        released_after: Optional[str] = None,
        released_before: Optional[str] = None,
        level: Optional[str] = None,
        protein_only: Optional[bool] = None,
        status_min: Optional[int] = None,
        require_candidates: bool = True,
        use_rule_selection: bool = True,
        checkpoint_stage: Optional[int] = None,
        extra_data: Optional[dict] = None,
        persist: bool = True,
        allow_all: bool = False,
        primary_only: Optional[bool] = None,
        use_busco: Optional[bool] = None,
        min_completeness: Optional[float] = None,
        min_single_copy_complete: Optional[float] = None,
        allow_duplicate_species: Optional[bool] = None,
        include_paralog_filtering_in_score: Optional[bool] = None,
        include_decontamination_in_score: Optional[bool] = None,
        paralog_run_id: Optional[str] = None,
        decontamination_run_id: Optional[str] = None,
        allow_ambiguous_contaminants: Optional[bool] = None,
        strict_decontamination: Optional[bool] = None,
        rescue_duplicates: Optional[bool] = None,
        paralog_filtered: Optional[bool] = None,
        not_paralog_filtered: Optional[bool] = None,
        min_hidden_paralogs: Optional[float] = None,
        max_hidden_paralogs: Optional[float] = None,
        decontaminated: Optional[bool] = None,
        not_decontaminated: Optional[bool] = None,
        contaminated: Optional[bool] = None,
        decontamination_run: Optional[str] = None,
        ignore_contaminated_assemblies: Optional[bool] = None,
    ) -> List[str]:
        """Resolve and persist selector-derived accessions for reuse on resume.

        Returns the cached selection when available; otherwise resolves using the
        shared selector helpers and stores the result under ``target_key``. When
        ``checkpoint_stage`` is provided the stage and payload are persisted via
        :meth:`checkpoint`; otherwise task data is updated in-place (when
        ``persist`` is True).
        """

        from ..selector_utils import normalize_accessions

        if target_key in self.data:
            cached_raw = self.data.get(target_key) or []
            cached = list(dict.fromkeys(normalize_accessions(cached_raw)))
            if cached or not require_candidates:
                return cached

        if allow_duplicate_species is None:
            allow_duplicate_species = bool(self.data.get("allow_duplicate_species", False))
        if include_paralog_filtering_in_score is None:
            include_paralog_filtering_in_score = self.data.get("include_paralog_filtering_in_score")
        if include_decontamination_in_score is None:
            include_decontamination_in_score = self.data.get("include_decontamination_in_score")
        if paralog_run_id is None:
            paralog_run_id = self.data.get("use_paralog_run")
        if decontamination_run_id is None:
            decontamination_run_id = self.data.get("use_decontamination_run")
        if allow_ambiguous_contaminants is None:
            allow_ambiguous_contaminants = self.data.get("allow_ambiguous_contaminants")
        if strict_decontamination is None:
            strict_decontamination = self.data.get("strict_decontamination")
        if rescue_duplicates is None:
            rescue_duplicates = self.data.get("rescue_duplicates")
        if paralog_filtered is None:
            paralog_filtered = bool(self.data.get("paralog_filtered", False))
        if not_paralog_filtered is None:
            not_paralog_filtered = bool(self.data.get("not_paralog_filtered", False))
        if min_hidden_paralogs is None:
            min_hidden_paralogs = self.data.get("min_hidden_paralogs")
        if max_hidden_paralogs is None:
            max_hidden_paralogs = self.data.get("max_hidden_paralogs")
        if decontaminated is None:
            decontaminated = bool(self.data.get("decontaminated", False))
        if not_decontaminated is None:
            not_decontaminated = bool(self.data.get("not_decontaminated", False))
        if contaminated is None:
            contaminated = bool(self.data.get("contaminated", False))
        if decontamination_run is None:
            decontamination_run = self.data.get("decontamination_run")
        if ignore_contaminated_assemblies is None:
            ignore_contaminated_assemblies = self.data.get("ignore_contaminated_assemblies")

        simple_candidate_only = (
            use_rule_selection
            and rule_quantity is None
            and rule_rank is None
            and not self.data.get("ranks")
            and not self.data.get("quantities")
            and not self.data.get("has_busco_results", False)
            and not self.data.get("missing_busco_results", False)
            and busco_library_id is None
            and use_busco in (None, False)
            and min_completeness is None
            and min_single_copy_complete is None
            and not allow_duplicate_species
            and include_paralog_filtering_in_score is None
            and include_decontamination_in_score is None
            and paralog_run_id is None
            and decontamination_run_id is None
            and allow_ambiguous_contaminants is None
            and strict_decontamination is None
            and rescue_duplicates is None
            and not paralog_filtered
            and not not_paralog_filtered
            and min_hidden_paralogs is None
            and max_hidden_paralogs is None
            and not decontaminated
            and not not_decontaminated
            and not contaminated
            and decontamination_run is None
            and ignore_contaminated_assemblies is None
        )

        if use_rule_selection:
            if simple_candidate_only:
                selected = self.selector_candidates(
                    root=root,
                    taxid=taxid,
                    additional=additional,
                    allow_all=allow_all,
                    downloaded_only=downloaded_only,
                    not_downloaded=not_downloaded or False,
                    released_after=released_after,
                    released_before=released_before,
                    level=level,
                    protein_only=protein_only,
                    status_min=status_min,
                    primary_only=primary_only,
                )
                if not selected and require_candidates:
                    raise ValueError("No accessions matched the provided selectors.")
            else:
                selected = self.resolve_selector_accessions(
                    root=root,
                    taxid=taxid,
                    additional=additional,
                    allow_all=allow_all,
                    rule_quantity=rule_quantity,
                    rule_rank=rule_rank,
                    busco_library_id=busco_library_id,
                    downloaded_only=downloaded_only,
                    not_downloaded=not_downloaded or False,
                    released_after=released_after,
                    released_before=released_before,
                    level=level,
                    protein_only=protein_only,
                    status_min=status_min,
                    primary_only=primary_only,
                    use_busco=use_busco,
                    min_completeness=min_completeness,
                    min_single_copy_complete=min_single_copy_complete,
                    allow_duplicate_species=allow_duplicate_species,
                    include_paralog_filtering_in_score=include_paralog_filtering_in_score,
                    include_decontamination_in_score=include_decontamination_in_score,
                    paralog_run_id=paralog_run_id,
                    decontamination_run_id=decontamination_run_id,
                    allow_ambiguous_contaminants=allow_ambiguous_contaminants,
                    strict_decontamination=strict_decontamination,
                    rescue_duplicates=rescue_duplicates,
                    paralog_filtered=paralog_filtered,
                    not_paralog_filtered=not_paralog_filtered,
                    min_hidden_paralogs=min_hidden_paralogs,
                    max_hidden_paralogs=max_hidden_paralogs,
                    decontaminated=decontaminated,
                    not_decontaminated=not_decontaminated,
                    contaminated=contaminated,
                    decontamination_run=decontamination_run,
                    ignore_contaminated_assemblies=ignore_contaminated_assemblies,
                    require_candidates=require_candidates,
                )
        else:
            selected = self.selector_candidates(
                root=root,
                taxid=taxid,
                additional=additional,
                allow_all=allow_all,
                downloaded_only=downloaded_only,
                not_downloaded=not_downloaded or False,
                released_after=released_after,
                released_before=released_before,
                level=level,
                protein_only=protein_only,
                status_min=status_min,
                primary_only=primary_only,
            )
            if not selected and require_candidates:
                raise ValueError("No accessions matched the provided selectors.")

        final = list(dict.fromkeys(normalize_accessions(selected)))
        if not final and require_candidates:
            raise ValueError("No accessions matched the provided selectors.")

        payload = dict(extra_data or {})
        payload[target_key] = final

        if checkpoint_stage is not None:
            self.checkpoint(checkpoint_stage, payload)
        else:
            self.data.update(payload)
            if persist:
                self.db_manager.tasks.update_data(self.task_id, data=self.data)

        return final

    # Lightweight logging helpers
    def log(self, msg: str, level: str = "INFO", category: str | None = None, **kwargs):
        extra = {"stage": getattr(self, "stage", None)}
        extra.setdefault("log_category", category or getattr(self, "default_log_category", "TASK"))
        extra.update(kwargs.pop("extra", {}))
        lvl = level.upper()
        if lvl == "DEBUG":
            self.logger.debug(msg, extra=extra)
        elif lvl in ("WARN", "WARNING"):
            self.logger.warning(msg, extra=extra)
        elif lvl == "ERROR":
            self.logger.error(msg, extra=extra)
        else:
            self.logger.info(msg, extra=extra)

    def error(self, msg: str, exc: bool = False, category: str | None = None, **kwargs):
        extra = {"stage": getattr(self, "stage", None)}
        extra.setdefault("log_category", category or getattr(self, "default_log_category", "TASK"))
        extra.update(kwargs.pop("extra", {}))
        self.logger.error(msg, exc_info=exc, extra=extra)

    @abstractmethod
    def run(self):
        """Implement the main logic of the task. Return False to skip marking complete."""
        raise NotImplementedError

    def start(self):
        """Wrapper around run() with lifecycle and error handling."""
        try:
            if self.db_manager.cursor is None:
                self.db_manager.connect()
            configure_logging_from_db(self.db_manager, force=False)
            self.db_manager.set_busco_run_context(
                pipeline=self.data.get("busco_pipeline"),
                input_mode=self.data.get("busco_input_mode"),
                prefer_pipeline=self.data.get("prefer_busco_pipeline"),
                prefer_input_mode=self.data.get("prefer_busco_input_mode"),
                run_ids=self.data.get("busco_run_ids"),
                selection=self.data.get("busco_run_selection") or "primary",
            )
            self.update_on_start()
            result = self.run()
            # Only mark complete on explicit True. False/None => suspended or deferred.
            if result is True:
                self.update_on_complete()
            elif result == "ERROR":
                # Ensure task is marked errored if not already
                # The error will have already been set we just need to set the status to error
                if self.status != "E":
                    self.update_on_error()
        except Exception as e:  # boundary: convert every unhandled task failure into task state
            # In this case the error has not already been caught!
            self.handle_exception(e, context="Unhandled exception in task")
        finally:
            try:
                self.db_manager.close()
            except Exception as close_exc:  # boundary: closing must not hide the task outcome
                self.logger.warning(
                    "Failed to close task database connection: %s",
                    close_exc,
                    exc_info=True,
                    extra={"stage": getattr(self, "stage", None), "log_category": "DATABASE"},
                )

    def handle_exception(self, exc: Exception, context: str | None = None):
        """Log and persist error info to DB, then mark task as errored."""
        # Normalize and build message
        raw_exc = exc if isinstance(exc, BaseException) else Exception(str(exc))
        msg = self._format_error_message(raw_exc, context)
        stack = self._format_exception_stack(raw_exc)
        # Log once
        self.error(msg, exc=bool(stack))
        # Persist to DB and mark E
        try:
            with self.db_manager.transaction(operation="persist task failure"):
                self.db_manager.tasks.set_error(self.task_id, msg, stack)
                self.update_on_error()
        except Exception as persistence_exc:  # boundary: task failure persistence failure is reported with original failure context
            self.logger.critical(
                "Failed to persist task error state for task %s: %s; original failure: %s",
                self.task_id,
                persistence_exc,
                msg,
                exc_info=True,
                extra={"stage": getattr(self, "stage", None), "log_category": "DATABASE"},
            )
            failure = TaskExecutionError(
                f"Task {self.task_id} failed and its error state could not be persisted: {persistence_exc}"
            )
            failure.original_error = raw_exc
            raise failure from persistence_exc
        return "ERROR"

    def collect_batch_failure(self, item, operation: str, exc: BaseException) -> BatchItemError:
        """Record one independent item failure while allowing the batch to continue."""
        error = exc if isinstance(exc, BatchItemError) else BatchItemError(item, operation, exc)
        stack = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        self._batch_failures.append(
            BatchFailure(item=str(item), operation=str(operation), message=str(exc), stack=stack)
        )
        self.error(str(error), exc=bool(exc.__traceback__), category="TASK")
        return error

    def fail_if_batch_failures(self, summary: str = "Batch processing failed"):
        """Persist all collected item failures and return the task error sentinel."""
        if not self._batch_failures:
            return None
        preview = "; ".join(
            f"{failure.item}: {failure.message}" for failure in self._batch_failures[:5]
        )
        if len(self._batch_failures) > 5:
            preview += f"; ... and {len(self._batch_failures) - 5} more"
        stacks = "\n\n".join(
            (
                f"[{failure.operation}] {failure.item}: {failure.message}\n"
                f"{failure.stack}"
            ).rstrip()
            for failure in self._batch_failures
        )
        error = TaskExecutionError(
            f"{summary}: {len(self._batch_failures)} item(s) failed. {preview}"
        )
        error.batch_stack = stacks
        return self.handle_exception(error, context={"failed_items": len(self._batch_failures)})

    def _format_error_message(self, exc: BaseException, context: object | None) -> str:
        base = str(exc).strip() or exc.__class__.__name__
        if context is None:
            return base
        if isinstance(context, str):
            token = context.strip()
            return f"{token}: {base}" if token else base
        try:
            context_blob = json.dumps(context, sort_keys=True, default=str)
        except (TypeError, ValueError, RecursionError):
            context_blob = str(context)
        context_blob = context_blob.strip()
        if not context_blob:
            return base
        return f"{base} | context={context_blob}"

    def _format_exception_stack(self, exc: BaseException) -> str:
        batch_stack = getattr(exc, "batch_stack", None)
        if batch_stack:
            return str(batch_stack)
        tb = getattr(exc, "__traceback__", None)
        if tb is not None:
            return "".join(traceback.format_exception(type(exc), exc, tb))
        return ""

    def activate(self):
        """Activate the task"""
        if self.db_manager.cursor is None:
            self.db_manager.connect()
        self.db_manager.set_busco_run_context(
            pipeline=self.data.get("busco_pipeline"),
            input_mode=self.data.get("busco_input_mode"),
            prefer_pipeline=self.data.get("prefer_busco_pipeline"),
            prefer_input_mode=self.data.get("prefer_busco_input_mode"),
            run_ids=self.data.get("busco_run_ids"),
            selection=self.data.get("busco_run_selection") or "primary",
        )
        self.update_on_start()

    def complete(self):
        self.update_on_complete()
        self.db_manager.close()

    def checkpoint(self, stage, checkpoint_data=None):
        """Create a checkpoint for the task, modifying the data so that it can be resumed"""
        # Ensure stage is recorded correctly
        self.stage = stage
        if checkpoint_data:
            self.data.update(checkpoint_data)

        self.db_manager.tasks.update_data(self.task_id, checkpoint=stage, data=self.data)

    def update_start_time(self, start_time):
        """Update the start time of the task"""
        self.db_manager.tasks.update_start_time(self.task_id, start_time)
    
    def update_end_time(self, end_time):
        """Update the end time of the task"""
        self.db_manager.tasks.update_end_time(self.task_id, end_time)

    def update_on_start(self):
        """Update the task status to running"""
        self.update_status("R")
        self.update_start_time(datetime.now())
        self.log("Task started.", category="SCHEDULER")

    def update_on_complete(self):
        """Update the task status to completed"""
        self.update_status("C")
        self.update_end_time(datetime.now())
        self.log("Task completed.", category="SCHEDULER")

    def update_on_error(self):
        """Mark the task as errored and set end time."""
        self.update_status("E")
        self.update_end_time(datetime.now())
        self.log("Task errored.", level="ERROR", category="SCHEDULER")

    def _subtasks_state(self):
        """Return current state of subtasks: 'none' | 'complete' | 'pending' | 'error'."""
        subtasks = self.db_manager.tasks.get_subtasks(self.task_id)
        if not subtasks:
            return "none"
        statuses = [t[2] for t in subtasks]
        if any(s == "E" for s in statuses):
            return "error"
        if all(s == "C" for s in statuses):
            return "complete"
        return "pending"

    def wait_for_subtasks(self, timeout=None):
        """Non-blocking wait: suspend this task until its subtasks complete.

        Returns True if subtasks are complete (or none exist) and the task can continue now.
        Returns False if subtasks are still pending (task will be marked 'S').
        Returns "ERROR" if any subtask errored.
        """

        # This method is called following the queuing of a subtask.
        # We can wait the timeout period for subtasks to complete.
        # Wait at least the provided timeout (or 0 for immediate suspend in tests)
        min_wait = 0  # immediate suspend if timeout not given
        wait_for = max(min_wait, int(timeout) if timeout is not None else min_wait)
        time.sleep(wait_for)

        # Now we check if the subtasks have finished within this expected period
        state = self._subtasks_state()
        if state in ("none", "complete"):
            return True
        if state == "error":
            # Subtasks errored; let the parent handle on resume (or immediately if desired)
            return "ERROR"
        # Pending: suspend and let the daemon reactivate later
        self.update_status("S")
        self.log("Suspended awaiting subtasks...", level="DEBUG", category="SCHEDULER")
        return False

    def queue_subtask(self, job_type, status, priority, data=None):
        """Queue a subtask in the database"""
        # Auto-inject phase metadata so parents can attribute status to the latest attempt only
        payload = {}
        if isinstance(data, dict):
            payload.update(data)
        elif data is not None:
            # If non-dict, still carry as-is under 'payload'
            payload = {"payload": data}
        # Attach phase metadata if manage_subtasks set it
        phase_meta = getattr(self, "_phase_meta", None)
        if phase_meta:
            payload.setdefault("__stage", phase_meta.get("stage"))
            payload.setdefault("__gen", phase_meta.get("gen"))
        self.db_manager.tasks.queue(
            job_type=job_type,
            status=status,
            priority=priority,
            parent_id=self.task_id,
            data=payload or None
        )
        self.log(f"Subtask queued job_type={job_type} status={status}", category="SCHEDULER")

    # Generic helper: handle subtask error with bounded retries
    def retry_subtasks_on_error(self, retry_key: str, max_retries: int, queue_retry_fn):
        """If any subtask errored, increment retry counter and requeue via queue_retry_fn.

        Returns:
        - False: a retry was queued and task should suspend
        - "ERROR": retries exhausted; caller should error the task
        - None: no error condition detected
        """
        subtasks = self.db_manager.tasks.get_subtasks(self.task_id) or []
        statuses = [t[2] for t in subtasks]
        if any(s == "E" for s in statuses):
            retries = int(self.data.get(retry_key, 0))
            if retries < max_retries:
                # Queue retry and suspend
                try:
                    queue_retry_fn()
                except Exception as e:  # boundary: caller-supplied retry queue callback
                    self.error(f"Failed to queue retry: {e}", category="SCHEDULER")
                    return "ERROR"
                self.data[retry_key] = retries + 1
                # Persist data but keep same stage
                self.checkpoint(stage=getattr(self, "stage", None) or 0, checkpoint_data={})
                self.update_status("S")
                self.log(
                    f"Retry {self.data[retry_key]}/{max_retries} queued for errored subtask(s).",
                    level="WARNING",
                    category="SCHEDULER",
                )
                return False
            else:
                return "ERROR"
        return None

    # Helper: aggregate all child subtask errors and persist to this task
    def aggregate_subtask_errors(self, default_summary: str = "Subtask failed after retries", *,
                                 subtask_ids: list[int] | None = None) -> tuple[str, str]:
        """Aggregate errors from child tasks; optionally limit to specific subtask_ids."""
        try:
            if subtask_ids is None:
                errs = self.db_manager.tasks.get_errors_from_subtasks(self.task_id) or []
            else:
                errs = []
                for sid in subtask_ids:
                    info = self.db_manager.tasks.get_error_info(sid)
                    if info:
                        errs.append((sid, info[0], info[1]))
            snippets = []
            for entry in errs:
                if len(entry) == 3:
                    sid, message, _stack = entry
                else:
                    sid, message, _stack = None, entry[0], entry[1]
                text = str(message or "").strip()
                if not text:
                    continue
                snippets.append(f"{sid}: {text}" if sid is not None else text)
            if snippets:
                preview = "; ".join(snippets[:3])
                if len(snippets) > 3:
                    preview += f"; ... and {len(snippets) - 3} more"
                summary = f"{default_summary}. Child errors: {preview}"
            else:
                summary = default_summary
            stacks = ""
        except (IndexError, TypeError, ValueError) as exc:
            self.log(f"Could not format child task errors: {exc}", level="WARNING", category="TASK")
            summary = default_summary
            stacks = ""
        self.db_manager.tasks.set_error(self.task_id, summary, stacks)
        return summary, stacks

    def manage_subtasks(
        self,
        *,
        stage: int,
        queue_fn,
        done_fn,
        wait_seconds: int = 0,
    retry_key: str | None = None,
        max_retries: int = 0,
        incomplete_message_fn=None,
    retry_incomplete: bool = False,
    ):
        """One-stop orchestration for a subtask phase with retry and incomplete detection.

        Simplified flow with per-phase generation tagging so only the latest attempt's
        child tasks are considered when deciding status and aggregating errors.
        """
        # Helpers
        phase_gen_key = f"_stage_{stage}_gen"
        # If caller didn't pass a retry_key, use an internal per-stage key to avoid colliding with user config
        if retry_key is None:
            retry_key = f"_retries_stage_{stage}"
        def _set_phase_meta(gen: int):
            self._phase_meta = {"stage": stage, "gen": gen}
        def _clear_phase_meta():
            if hasattr(self, "_phase_meta"):
                delattr(self, "_phase_meta")
        def _current_gen() -> int:
            return int(self.data.get(phase_gen_key, 0))
        def _filter_phase_children(gen: int):
            tasks = self.db_manager.tasks.get_subtasks(self.task_id) or []
            filtered = []
            for t in tasks:
                # t[6] is data JSON
                try:
                    d = json.loads(t[6]) if t[6] else {}
                except (TypeError, json.JSONDecodeError):
                    d = {}
                if d.get("__stage") == stage and d.get("__gen") == gen:
                    filtered.append(t)
            return filtered
        def _state(subtasks):
            if not subtasks:
                return {"has_active": False, "has_error": False, "all_complete": False}
            statuses = [t[2] for t in subtasks]
            return {
                "has_active": any(s in ("P", "R") for s in statuses),
                "has_error": any(s == "E" for s in statuses),
                "all_complete": all(s == "C" for s in statuses),
            }
        def _done_for_phase(phase_tasks) -> bool:
            if done_fn is None:
                return _state(phase_tasks)["all_complete"]
            return bool(done_fn())
        def _queue_retry_or_error(reason: str, *, use_incomplete_msg: bool = False, phase_tasks=None, allow_retry: bool = True):
            # Retry if allowed
            retries = int(self.data.get(retry_key, 0)) if retry_key else 0
            if allow_retry and retry_key and retries < max_retries:
                try:
                    # bump generation
                    new_gen = _current_gen() + 1
                    self.data[phase_gen_key] = new_gen
                    _set_phase_meta(new_gen)
                    queue_fn()
                except Exception as e:  # boundary: caller-supplied retry queue callback
                    self.error(f"Failed to queue retry: {e}", category="SCHEDULER")
                    # Fall through to error recording below
                else:
                    self.data[retry_key] = retries + 1
                    self.checkpoint(stage=stage, checkpoint_data={phase_gen_key: self.data[phase_gen_key]})
                    self.update_status("S")
                    self.log(
                        f"Retry {self.data[retry_key]}/{max_retries} queued ({reason}).",
                        level="WARNING",
                        category="SCHEDULER",
                    )
                    _clear_phase_meta()
                    return False
            # Exhausted or not configured: persist error
            if use_incomplete_msg and incomplete_message_fn:
                try:
                    msg = incomplete_message_fn()
                except Exception as exc:  # boundary: optional diagnostic callback
                    self.log(
                        f"Failed to build incomplete-phase diagnostic: {exc}",
                        level="WARNING",
                        category="SCHEDULER",
                    )
                    msg = None
                if msg:
                    if isinstance(msg, tuple):
                        self.db_manager.tasks.set_error(self.task_id, msg[0], msg[1])
                    else:
                        self.db_manager.tasks.set_error(self.task_id, str(msg), "")
                else:
                    # Fallback aggregate
                    self.aggregate_subtask_errors("Phase incomplete after max retries", subtask_ids=[t[0] for t in (phase_tasks or [])])
            else:
                # Prefer last failed subtask's error for summary, include all stacks
                ids = [t[0] for t in (phase_tasks or []) if t[2] == "E"]
                if ids:
                    # Last failed = highest task_id among failures
                    last_id = max(ids)
                    last_err = self.db_manager.tasks.get_error_info(last_id) or (None, None)
                    all_summary, all_stacks = self.aggregate_subtask_errors(
                        f"All subtasks failed after {self.data.get(retry_key, 0)} retries.",
                        subtask_ids=ids,
                    )
                    summary = f"All subtasks failed after {self.data.get(retry_key, 0)} retries. Last error: {last_err[0] or ''}"
                    self.db_manager.tasks.set_error(self.task_id, summary, all_stacks)
                else:
                    self.aggregate_subtask_errors("Subtask failure after retries", subtask_ids=[t[0] for t in (phase_tasks or [])])
            _clear_phase_meta()
            return "ERROR"

        # Initialize retry counter if applicable
        if retry_key is not None:
            self.data.setdefault(retry_key, 0)

        queued = False

        # First-time entry into this phase
        if getattr(self, "stage", 0) < stage:
            try:
                # Start generation 1 for this stage
                new_gen = _current_gen() + 1
                self.data[phase_gen_key] = new_gen
                _set_phase_meta(new_gen)
                queued = queue_fn()
            except Exception as e:  # boundary: caller-supplied subtask queue callback
                self.handle_exception(e, context="Failed to queue subtasks")
                _clear_phase_meta()
                return "ERROR"
            # checkpoint to this stage so resumes land here; persist gen
            self.checkpoint(stage=stage, checkpoint_data={phase_gen_key: self.data[phase_gen_key]})
            # brief wait for fast subtasks
            decision = not queued or self.wait_for_subtasks(wait_seconds)
            if decision is True:
                # If work finished quickly: check completion; if not satisfied, treat as incomplete now
                try:
                    phase_tasks = _filter_phase_children(_current_gen())
                    if _done_for_phase(phase_tasks):
                        _clear_phase_meta()
                        return "CONTINUE"
                except Exception as e:  # boundary: caller-supplied completion predicate
                    self.handle_exception(e, context="Error evaluating done_fn")
                    _clear_phase_meta()
                    return "ERROR"
                # Not satisfied: handle as incomplete outcome immediately (retry or error)
                phase_tasks = _filter_phase_children(_current_gen())
                outcome = _queue_retry_or_error("incomplete outcome", use_incomplete_msg=True, phase_tasks=phase_tasks, allow_retry=retry_incomplete)
                _clear_phase_meta()
                return outcome
            if decision == "ERROR":
                # Retry or error using current generation's children
                phase_tasks = _filter_phase_children(_current_gen())
                return _queue_retry_or_error("error", use_incomplete_msg=False, phase_tasks=phase_tasks, allow_retry=True)
            # decision is False -> suspended awaiting children
            _clear_phase_meta()
            return False

        # Resumed or reentered this phase: consider only current generation
        gen = _current_gen() or 1
        phase_tasks = _filter_phase_children(gen)
        st = _state(phase_tasks)

        # If success criteria satisfied, move on
        try:
            if _done_for_phase(phase_tasks):
                return "CONTINUE"
        except Exception as e:  # boundary: caller-supplied completion predicate
            return self.handle_exception(e, context="Error evaluating done_fn")
            

        # If work is still in progress, stay suspended
        if st["has_active"]:
            self.update_status("S")
            self.log("Suspended awaiting active subtasks...", level="DEBUG", category="SCHEDULER")
            return False

        # Handle error/retry when no active children
        if st["has_error"]:
            return _queue_retry_or_error("error", use_incomplete_msg=False, phase_tasks=phase_tasks, allow_retry=True)

        # No active and no explicit error, but not done -> incomplete outcome
        if not st["has_active"] and not st["has_error"]:
            return _queue_retry_or_error("incomplete outcome", use_incomplete_msg=True, phase_tasks=phase_tasks, allow_retry=retry_incomplete)

        # Fallback: suspend
        self.update_status("S")
        self.log("Suspended awaiting subtasks (pending)...", level="DEBUG", category="SCHEDULER")
        return False
