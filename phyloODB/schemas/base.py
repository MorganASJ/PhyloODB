"""Shared Pydantic models and utilities for task payload validation."""
from __future__ import annotations

from typing import Any, Dict, Optional, List

from pydantic import BaseModel, ConfigDict, Field


class TaskPayload(BaseModel):
    """Base class for task payload schemas.

    All concrete task schemas should inherit from this base so we can control
    JSON serialisation and future common helpers (e.g. checkpoint merging).
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    selector_requested_accessions: Optional[List[str]] = Field(
        default=None,
        description="Internal: accessions resolved before task-specific selector filters were applied.",
    )
    selector_skipped_accessions: Optional[List[str]] = Field(
        default=None,
        description="Internal: resolved accessions removed by task-specific selector filters.",
    )
    allow_duplicate_species: bool = False
    allow_ambiguous_contaminants: Optional[bool] = Field(
        default=None,
        description="Treat decontamination 'unknown' BUSCOs as supported.",
    )
    strict_decontamination: Optional[bool] = Field(
        default=None,
        description="Treat decontamination 'weak' BUSCOs as contaminated.",
    )
    include_paralog_filtering_in_score: Optional[bool] = Field(
        default=None,
        description="Include paralog filtering results when calculating BUSCO selector scores.",
    )
    include_decontamination_in_score: Optional[bool] = Field(
        default=None,
        description="Include decontamination results when calculating BUSCO selector scores.",
    )
    use_decontamination_run: Optional[str] = Field(
        default=None,
        description="Use a specific decontamination run id when calculating BUSCO scores.",
    )
    use_paralog_run: Optional[str] = Field(
        default=None,
        description="Use a specific paralog filtering run id when calculating BUSCO scores.",
    )
    paralog_filtered: Optional[bool] = None
    not_paralog_filtered: Optional[bool] = None
    min_hidden_paralogs: Optional[float] = Field(
        default=None,
        ge=0,
        le=1,
        description="Minimum hidden paralog proportion (0-1).",
    )
    max_hidden_paralogs: Optional[float] = Field(
        default=None,
        ge=0,
        le=1,
        description="Maximum hidden paralog proportion (0-1).",
    )
    decontaminated: Optional[bool] = None
    not_decontaminated: Optional[bool] = None
    contaminated: Optional[bool] = None
    decontamination_run: Optional[str] = Field(
        default=None,
        description="Filter selectors to accessions included in a decontamination run id.",
    )
    ignore_contaminated_assemblies: Optional[bool] = Field(
        default=None,
        description="Exclude assemblies whose latest decontamination decision is CONTAMINATED.",
    )
    required_threads: Optional[int] = Field(
        default=None,
        ge=1,
        description="Override the daemon required thread count for this task.",
    )
    busco_pipeline: Optional[str] = Field(
        default=None,
        description="Require BUSCO runs from a specific pipeline (e.g. miniprot, metaeuk, augustus).",
    )
    busco_input_mode: Optional[str] = Field(
        default=None,
        description="Require BUSCO runs by input format (protein or genome).",
    )
    prefer_busco_pipeline: Optional[str] = Field(
        default=None,
        description="Prefer BUSCO runs from a specific pipeline while allowing fallback to other matching runs.",
    )
    prefer_busco_input_mode: Optional[str] = Field(
        default=None,
        description="Prefer BUSCO runs by input format (protein or genome) while allowing fallback to other matching runs.",
    )
    busco_export_format: Optional[str] = Field(
        default=None,
        description="Require the selected BUSCO run to support protein or nucleotide export.",
    )
    busco_run_ids: Optional[List[str]] = Field(
        default=None,
        description="Limit BUSCO-aware selection to specific BUSCO run ids.",
    )
    busco_run_selection: Optional[str] = Field(
        default=None,
        description="BUSCO run-selection policy for downstream analyses (primary or latest).",
    )

    def as_task_data(self) -> Dict[str, Any]:
        """Return data ready to persist into the database."""
        return self.model_dump(mode="json", exclude_none=True)

    def with_overrides(self, **updates: Any) -> "TaskPayload":
        """Return a clone of the payload with patched fields.

        This is useful for staging retries or checkpoint restores while keeping
        schema validation in place.
        """

        data = self.model_dump()
        data.update(updates)
        return self.model_copy(update=data)


class DaemonConfig(BaseModel):
    """Lightweight schema describing daemon runtime requirements.

    This will remain optional but allows specs to communicate constraints like
    required CPU threads or mutually exclusive scheduling rules.
    """

    model_config = ConfigDict(extra="forbid")

    required_threads: int = 1
    description: Optional[str] = None
    tags: Optional[tuple[str, ...]] = None
