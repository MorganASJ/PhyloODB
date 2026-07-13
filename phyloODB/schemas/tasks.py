"""Pydantic payload schemas for built-in PhyloODB tasks."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Literal, Optional
from datetime import datetime

from pydantic import BaseModel, Field, PositiveInt, field_validator, model_validator

from .base import TaskPayload


class UpdateAssemblyInformationPayload(TaskPayload):
    taxid: Optional[int] = Field(default=None, ge=1)
    accessions: Optional[List[str]] = Field(default=None, min_length=1)
    force_update: bool = False
    after: Optional[str] = Field(default=None, description="Include assemblies released on/after YYYY-MM-DD.")
    before: Optional[str] = Field(default=None, description="Include assemblies released on/before YYYY-MM-DD.")
    level: Optional[Literal["complete genome", "chromosome", "scaffold", "contig"]] = Field(
        default=None, description="Filter assemblies by assembly level."
    )
    exclude_accessions: Optional[List[str]] = Field(default=None, description="Comma-separated accessions to exclude.")
    exclude_taxids: Optional[List[int]] = Field(default=None, description="Taxids to exclude (descendants included).")
    exclude_clades: Optional[List[str]] = Field(default=None, description="Clades to exclude (descendants included).")
    primary_only: bool = False
    debug_path: Optional[Path] = Field(default=None, description="Write fetched metadata JSON to this directory.")

    @model_validator(mode="after")
    def _check_source(self) -> "UpdateAssemblyInformationPayload":
        if bool(self.taxid) == bool(self.accessions):
            raise ValueError("Provide exactly one of 'taxid' or 'accessions'.")
        return self

    @field_validator("after", "before")
    @classmethod
    def _validate_dates(cls, value: Optional[str]) -> Optional[str]:
        if value is None:
            return value
        try:
            datetime.strptime(value, "%Y-%m-%d")
        except ValueError as exc:
            raise ValueError("Dates must be in YYYY-MM-DD format.") from exc
        return value


class DownloadAssembliesPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    taxid: Optional[int] = Field(default=None, ge=1)
    protein: bool = False
    max_concurrent: PositiveInt = Field(default=1)
    force_redownload: bool = False
    download_retries: int = Field(default=0, ge=0, description="Number of retry attempts per accession on failure.")
    quantity: Optional[int] = Field(default=None, ge=1)
    rank: Optional[str] = None
    exclude_accessions: Optional[List[str]] = Field(default=None, description="Comma-separated accessions to exclude.")
    exclude_taxids: Optional[List[int]] = Field(default=None, description="Taxids to exclude (descendants included).")
    exclude_clades: Optional[List[str]] = Field(default=None, description="Clades to exclude (descendants included).")
    use_busco: bool = Field(default=False, description="Use BUSCO results when ranking (default off for downloads).")
    min_completeness: Optional[float] = Field(default=None, ge=0, le=100, description="Minimum BUSCO completeness (percent).")
    min_single_copy_complete: Optional[float] = Field(default=None, ge=0, le=100, description="Minimum BUSCO single-copy completeness (percent).")
    primary_only: bool = False
    clean_isoforms: bool = True
    skip_clean_isoforms: bool = False
    clean_skip_gff: bool = False
    clean_skip_cdhit: bool = False
    clean_gff_priority: bool = False
    clean_cdhit_identity: Optional[float] = Field(default=None, gt=0.0, le=1.0)
    clean_max_concurrent: PositiveInt = Field(default=1, le=8)
    clean_threads_per_job: PositiveInt = Field(default=1)

    @model_validator(mode="after")
    def _ensure_sources(self) -> "DownloadAssembliesPayload":
        if not self.accessions and self.taxid is None:
            raise ValueError("Provide at least one accession or a taxid.")
        if self.rank and self.taxid is None and not self.accessions:
            raise ValueError("Rule-based selection by rank requires a taxid.")
        return self


class BatchImportLocalAssemblyTaskPayload(TaskPayload):
    assembly_dir: Path
    accessions_for_import: Optional[List[str]] = Field(default=None)
    clean_isoforms: bool = True
    skip_clean_isoforms: bool = False
    clean_skip_gff: bool = False
    clean_skip_cdhit: bool = False
    clean_gff_priority: bool = False
    clean_cdhit_identity: Optional[float] = Field(default=None, gt=0.0, le=1.0)
    clean_max_concurrent: PositiveInt = Field(default=1, le=8)
    clean_threads_per_job: PositiveInt = Field(default=1)

class ImportLocalAssemblyPayload(TaskPayload):
    fna: Optional[Path] = None
    faa: Optional[Path] = None
    gff: Optional[Path] = None
    others: List[Path] = Field(default_factory=list)
    accession: Optional[str] = None
    metadata: Dict[str, object] = Field(default_factory=dict)
    taxid: Optional[int] = Field(default=None, ge=1)
    taxon_name: Optional[str] = None
    genus: Optional[str] = None
    species: Optional[str] = None
    location: Optional[Path] = None
    copy_to_genome_dir: bool = True
    clean_isoforms: bool = True
    skip_clean_isoforms: bool = False
    clean_skip_gff: bool = False
    clean_skip_cdhit: bool = False
    clean_gff_priority: bool = False
    clean_cdhit_identity: Optional[float] = Field(default=None, gt=0.0, le=1.0)
    clean_max_concurrent: PositiveInt = Field(default=1, le=8)
    clean_threads_per_job: PositiveInt = Field(default=1)

    @field_validator("others", mode="before")
    @classmethod
    def _parse_others(cls, value):
        if value is None:
            return []
        if isinstance(value, (str, Path)):
            return [value]
        if isinstance(value, (list, tuple, set)):
            return list(value)
        raise TypeError("Invalid type for 'others'.")

    @model_validator(mode="after")
    def _check_files(self) -> "ImportLocalAssemblyPayload":
        if not self.fna and not self.faa:
            raise ValueError("Provide at least one of 'fna' or 'faa'.")
        id_fields = [bool(self.taxid), bool(self.taxon_name), bool(self.genus and self.species)]
        if not any(id_fields):
            raise ValueError("Provide either taxid, taxon_name, or both genus and species.")
        if self.location and self.copy_to_genome_dir:
            raise ValueError("'location' and 'copy_to_genome_dir' cannot both be set.")
        return self


class DownloadBuscoLibraryPayload(TaskPayload):
    lineage: str
    libraries_dir: Optional[Path] = None
    busco_path: Optional[Path] = None
    parent_library_name: Optional[str] = None
    coverage: Optional[int] = Field(default=None, ge=1)
    size: Optional[int] = Field(default=None, ge=1)
    debug_skip_dl: bool = False


class BuscoTaskPayload(TaskPayload):
    lineage: str
    library: Optional[str] = None
    format: Literal["auto", "protein", "genome", "nucleotide"] = "auto"
    pipeline: Optional[Literal["auto", "miniprot", "metaeuk", "augustus"]] = "auto"
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    augustus_evalue: Optional[float] = Field(default=None, ge=0.0)
    augustus_limit: Optional[int] = Field(default=None, ge=1)
    augustus_long: Optional[bool] = None
    augustus_species: Optional[str] = None
    augustus_parameters: Optional[str] = None
    metaeuk_parameters: Optional[str] = None
    metaeuk_rerun_parameters: Optional[str] = None
    miniprot_parameters: Optional[str] = None
    accession: str
    output_path: Optional[Path] = None
    force: bool = False
    keep_miniprot_ref_file: bool = Field(
        default=False,
        description="Retain miniprot ref.mpi file in BUSCO output.",
    )
    busco_lib_wait_seconds: int = Field(default=0, ge=0)
    busco_lib_retries: int = Field(default=0, ge=0)


class OrthoFinderPayload(TaskPayload):
    input_dir: Optional[Path] = None
    out_dir: Optional[Path] = Field(default=None, description="Optional export directory. Defaults under the active exports root.")
    force: bool = False
    accessions: List[str] = Field(default_factory=list)
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    check_for_previous_run_folders: bool = Field(default=True)
    mcl_inflation: Optional[float] = Field(
        default=None,
        gt=0.0,
        description="Optional OrthoFinder MCL inflation value passed through as -I.",
    )

    @model_validator(mode="after")
    def _validate_inputs(self) -> "OrthoFinderPayload":
        if not self.out_dir:
            raise ValueError("out_dir is required")
        if not self.input_dir and not self.accessions:
            raise ValueError("input_dir is required when accessions are not provided")
        return self


class MafftRunPayload(TaskPayload):
    input_fasta: Path
    out_dir: Path
    output_name: Optional[str] = None
    mafft_flags: Optional[str] = None


class IQTreeRunPayload(TaskPayload):
    input_alignment: Path
    out_dir: Path
    prefix: Optional[str] = None
    iqtree_flags: Optional[str] = None


class AddLibraryPayload(TaskPayload):
    name: str
    coverage: Optional[str] = None
    coverage_taxid: Optional[int] = None
    accessions: List[str] = Field(min_length=2)
    parent_library_name: str
    parent_library_id: Optional[int] = Field(default=None, ge=1)
    library_id: Optional[int] = Field(default=None, ge=1)
    location: Optional[Path] = None
    metadata_wait_seconds: int = Field(default=5, ge=0)
    assembly_metadata_retries: int = Field(default=3, ge=0)
    download_wait_seconds: int = Field(default=0, ge=0)
    download_retries: int = Field(default=2, ge=0)
    busco_retries: int = Field(default=0, ge=0)
    busco_wait_seconds: int = Field(default=0, ge=0)
    min_species_in_trees: PositiveInt = Field(
        default=4,
        description="Minimum species occupancy required for 1-to-1 BUSCO/orthogroup matches.",
    )
    rerun_busco: bool = Field(
        default=False,
        description="Force fresh BUSCO runs for the reference accessions instead of reusing existing matching runs.",
    )
    rerun_orthofinder: bool = Field(
        default=False,
        description="Force a fresh OrthoFinder run for the reference accession set instead of reusing an existing matching run.",
    )
    rerun_gene_trees: bool = Field(
        default=False,
        description="Force fresh IQ-TREE orthogroup trees even when reusable trees are already present in IQ-TREE_Orthogroup_trees.",
    )
    skip_paralog_analysis: bool = Field(
        default=False,
        description="Accept exact 1-to-1 BUSCO/orthogroup families directly after occupancy filtering and skip paralog analysis.",
    )
    gene_tree_source: Literal["iqtree", "fasttree"] = Field(
        default="iqtree",
        description="Gene tree source for paralog-aware add-library analysis: iqtree builds IQ-TREE_Orthogroup_trees, fasttree reuses OrthoFinder Resolved_Gene_Trees.",
    )
    orthofinder_threads: Optional[PositiveInt] = Field(
        default=None,
        description="Override the OrthoFinder child task thread allocation.",
    )
    orthofinder_mcl_inflation: Optional[float] = Field(
        default=None,
        gt=0.0,
        description="Optional OrthoFinder MCL inflation value for the child orthofinder-run.",
    )
    clean_refs: bool = Field(
        default=False,
        description="Create OrthoFinder-derived BUSCO runs for reference taxa and mark accession-specific out-paralog BUSCOs as duplicated.",
    )
    clean_refs_strict: bool = Field(
        default=True,
        description="Create OrthoFinder-derived BUSCO runs for reference taxa and mark accession-specific in- or out-paralog BUSCOs as duplicated.",
    )
    set_cleaned_primary: bool = Field(
        default=True,
        description="Set the derived OrthoFinder BUSCO run as the default primary BUSCO run for each cleaned reference accession.",
    )
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    mafft_threads: Optional[PositiveInt] = Field(
        default=None,
        description="Override per-gene MAFFT thread count for replacement orthogroup alignments.",
    )
    iqtree_threads: Optional[PositiveInt] = Field(
        default=None,
        description="Override per-gene IQ-TREE thread count for replacement orthogroup trees.",
    )
    mafft_flags: Optional[str] = None
    iqtree_flags: Optional[str] = None
    annotate_og_trees: bool = Field(
        default=False,
        description="Run annotate-orthogroup-tree after replacement gene trees are built.",
    )
    force: bool = Field(
        default=False,
        description="Allow rebuilding an existing derived library by purging its stored library state and output directory. This does not by itself rerun BUSCO or OrthoFinder; use --rerun-busco and/or --rerun-orthofinder for that.",
    )
    debug_path: Optional[Path] = Field(default=None, description="Write fetched metadata JSON to this directory.")

    def as_task_data(self) -> Dict[str, object]:
        data = super().as_task_data()
        # Keep help/defaults strict-by-default, but do not activate cleaning unless
        # the user explicitly requested one of the clean flags.
        strict_requested = "clean_refs_strict" in self.model_fields_set
        clean_requested = "clean_refs" in self.model_fields_set and bool(self.clean_refs)
        if not clean_requested and not strict_requested:
            data.pop("clean_refs_strict", None)
        if not clean_requested and not strict_requested:
            data.pop("set_cleaned_primary", None)
        return data
    

class ImportCustomLibraryPayload(TaskPayload):
    library_name: str
    coverage: Optional[str] = None
    coverage_taxid: Optional[int] = None
    parent_library_name: Optional[str] = None
    parent_library_id: Optional[int] = Field(default=None, ge=1)
    busco_ids: List[str] = Field(
        min_length=1,
        description="BUSCO ids (comma-separated or a file path when a single value is supplied).",
    )
    ref_accessions: Optional[List[str]] = Field(
        default=None,
        description="Optional reference accessions to attach to the library.",
    )
    force: bool = False

    @model_validator(mode="after")
    def _validate_inputs(self) -> "ImportCustomLibraryPayload":
        if not self.library_name or not str(self.library_name).strip():
            raise ValueError("library_name is required.")
        if not self.busco_ids:
            raise ValueError("busco_ids is required.")
        if not (self.parent_library_id or (self.parent_library_name and str(self.parent_library_name).strip())):
            raise ValueError("Provide parent_library_id or parent_library_name.")
        if self.coverage_taxid is None and (self.coverage is None or not str(self.coverage).strip()):
            raise ValueError("Provide coverage or coverage_taxid.")
        return self


class ExportLibraryPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    accession: Optional[str] = None
    taxid: Optional[int] = Field(default=None, ge=1)
    clade: Optional[str] = None
    require: List[str] = Field(
        default_factory=list,
        description="Required clause list for export headers/families. Comma-separated clauses are ANDed; each clause may use '|' for OR (e.g. x,y,(w|u)).",
    )
    quantity: Optional[int] = Field(default=None, ge=1)
    rank: Optional[str] = None
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    sequence_type: Literal["protein", "nucleotide"] = "protein"
    proteome_profile: Optional[str] = None
    proteome_profiles: List[str] = Field(default_factory=list)
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    protein_only: bool = False
    out_dir: Optional[Path] = Field(
        default=None,
        description="Optional export directory. Defaults under the active exports root.",
    )
    rerun: bool = False
    disable_paralog_filter: bool = Field(default=False, description="Skip paralog-filter results when exporting.")
    disable_decont_filter: bool = Field(default=False, description="Skip decontamination run filtering when exporting.")
    require_paralog_filtering: bool = Field(
        default=False,
        description="Require paralog-filtering results for selected accessions. By default export uses those results when present but does not require them.",
    )
    require_decontamination: bool = Field(
        default=False,
        description="Require decontamination results for selected accessions. By default export uses those results when present but does not require them.",
    )
    write_lineage_csv: bool = Field(default=True, description="Write lineage CSV for selected accessions.")
    write_busco_report: bool = Field(default=True, description="Write BUSCO/decontamination report for selected accessions.")
    write_busco_family_matrix: bool = Field(
        default=True,
        description="Write per-family BUSCO status matrix for selected accessions.",
    )
    lineage_csv_path: Optional[Path] = Field(default=None, description="Optional path for lineage CSV output.")
    busco_report_path: Optional[Path] = Field(default=None, description="Optional path for BUSCO report output.")
    busco_family_matrix_path: Optional[Path] = Field(
        default=None,
        description="Optional path for BUSCO family matrix output.",
    )
    busco_report_extended: bool = Field(default=False, description="Include extended decontamination breakdown columns in BUSCO report.")
    rescue_duplicates: bool = Field(
        default=False,
        description="Treat duplicated BUSCO families with exactly one copy passing the active filters as effective single-copy for this export.",
    )
    retain_headers: bool = Field(default=False, description="Preserve original BUSCO headers in exported FASTAs.")
    header: Optional[str] = Field(
        default=None,
        description="Custom header template for exported FASTAs. Supports ACCESSION/TAXON/GENUS/SPECIES/RANK/FAMILY/LENGTH/GENE/TAXID/BITSCORE; TAXON falls back to <taxon>_sp for non-species assemblies.",
    )
    header_rank: Optional[str] = Field(default=None, description="Rank token used when HEADER includes RANK.")
    family_ids: List[str] = Field(
        default_factory=list,
        description="Optional BUSCO family ids to export. Accepts repeated values, a comma-separated value, or a file path when a single value is supplied.",
    )
    include_duplicated: bool = Field(
        default=False,
        description="Include duplicated BUSCO copies in exported family FASTAs instead of restricting export to single-copy hits only.",
    )
    min_completeness: Optional[float] = Field(default=None, ge=0, le=1, description="Minimum BUSCO completeness (0-1) for selector filtering.")
    min_single_copy_complete: Optional[float] = Field(default=None, ge=0, le=1, description="Minimum BUSCO single-copy completeness (0-1) for selector filtering.")
    min_occupancy: float = Field(default=0.5, ge=0, le=1, description="Minimum fraction of selected accessions that must be present in a family to retain the family.")
    min_taxa_occupancy: float = Field(default=0.3, ge=0, le=1, description="Minimum fraction of families an accession must appear in to be retained.")
    min_completeness: Optional[float] = Field(default=None, ge=0, le=1, description="Minimum BUSCO completeness (0-1) for selector filtering.")
    min_single_copy_complete: Optional[float] = Field(default=None, ge=0, le=1, description="Minimum BUSCO single-copy completeness (0-1) for selector filtering.")

    @model_validator(mode="after")
    def _check_inputs(self) -> "ExportLibraryPayload":
        candidates = list(self.accessions)
        if self.accession:
            candidates.append(self.accession)
        if not candidates and self.taxid is None and not (self.clade and str(self.clade).strip()):
            raise ValueError("Provide accessions, a taxid, or a clade for export.")
        if self.quantity is not None and self.quantity <= 0:
            raise ValueError("quantity must be positive when provided.")
        if self.rank and self.taxid is None and not candidates and not (self.clade and str(self.clade).strip()):
            raise ValueError("Rule-based selection by rank requires a taxid or clade.")
        return self


class BuildBuscoTreesPayload(ExportLibraryPayload):
    mafft_threads: PositiveInt = Field(
        default=2,
        description="Thread count for queued MAFFT child tasks.",
    )
    iqtree_threads: PositiveInt = Field(
        default=4,
        description="Thread count for queued IQ-TREE child tasks.",
    )
    mafft_flags: Optional[str] = None
    iqtree_flags: Optional[str] = None


class AnnotateOrthogroupTreePayload(TaskPayload):
    input_tree: Optional[Path] = None
    input_dir: Optional[Path] = None
    out_dir: Optional[Path] = None
    manifest_tsv: Optional[Path] = None
    mapping_tsv: Optional[Path] = None
    species_paralog_tsv: Optional[Path] = None
    orthofinder_location: Optional[Path] = None
    source_run_ids: List[int] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check_tree_inputs(self) -> "AnnotateOrthogroupTreePayload":
        if bool(self.input_tree) == bool(self.input_dir):
            raise ValueError("Provide exactly one of input_tree or input_dir.")
        return self


class GenerateLineageCsvPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    taxid: Optional[int] = Field(default=None, ge=1)
    quantity: Optional[int] = Field(default=None, ge=1)
    rank: Optional[str] = None
    library_id: Optional[int] = Field(default=None, ge=1)
    lineage: Optional[str] = None
    protein_only: bool = False
    output: Path
    exclude_accessions: Optional[List[str]] = Field(default=None, description="Comma-separated accessions to exclude.")
    exclude_taxids: Optional[List[int]] = Field(default=None, description="Taxids to exclude (descendants included).")
    exclude_clades: Optional[List[str]] = Field(default=None, description="Clades to exclude (descendants included).")
    primary_only: bool = False
    use_busco: bool = Field(default=True, description="Use BUSCO results when ranking.")
    min_completeness: Optional[float] = Field(default=None, ge=0, le=100, description="Minimum BUSCO completeness (percent).")
    min_single_copy_complete: Optional[float] = Field(default=None, ge=0, le=100, description="Minimum BUSCO single-copy completeness (percent).")

    @model_validator(mode="after")
    def _validate_sources(self) -> "GenerateLineageCsvPayload":
        candidates = list(self.accessions)
        if not candidates and self.taxid is None:
            raise ValueError("Provide accessions or a taxid for lineage export.")
        if self.quantity is not None and self.quantity <= 0:
            raise ValueError("quantity must be positive when provided.")
        if self.rank and self.taxid is None and not candidates:
            raise ValueError("Rule-based selection by rank requires a taxid.")
        return self


class GenomeMoveRowPayload(BaseModel):
    accession: str
    source_path: str
    destination_path: str
    source_root_id: Optional[int] = None
    destination_root_id: int
    action: Literal["move-files", "rebind-only"] = "move-files"


class FinalizeGenomeMovePayload(TaskPayload):
    rows: List[GenomeMoveRowPayload] = Field(min_length=1)
    verify: bool = True
    tidy: bool = True


class CreateTaxonomyPayload(TaskPayload):
    path_to_taxdump: Optional[Path] = None
    retain_taxdump: bool = False
    working_dir: Optional[Path] = None


class ExampleTaskPayload(TaskPayload):
    note: Optional[str] = None

class CreateProteomeBlastDBPayload(TaskPayload):
    accession: str
    output_path: Optional[Path] = None
    location: Optional[Path] = None
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    force: bool = False

    @model_validator(mode="after")
    def _validate_library(self) -> "CreateProteomeBlastDBPayload":
        if not self.library_id and not self.library_name:
            raise ValueError("Either library_id or library_name must be provided.")
        return self
    
class ParalogRemovalPayload(TaskPayload):
    mode: Literal["median", "percent", "bitscore", "lower-quartile", "upper-quartile"] = "median"
    percentile: Optional[float] = Field(default=None, gt=0.0, le=100.0)
    bitscore_threshold: Optional[float] = Field(default=None, ge=0.0)
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    ref_accessions: List[str] = Field(default_factory=list)
    accessions: List[str] = Field(
        default_factory=list,
        description="Accessions used to compute BUSCO family medians and reference expectations.",
    )
    targets: List[str] = Field(
        default_factory=list,
        description="Subset of accessions to write/update paralog-filtering results for.",
    )
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    report_dir: Optional[Path] = None
    run_label: Optional[str] = None
    force: bool = False
    max_concurrent: Optional[PositiveInt] = None
    reuse_existing: bool = Field(default=False)
    rebuild_proteome_dbs: bool = Field(
        default=False,
        description="Force rebuilding reference proteome BLAST databases before paralog removal.",
    )
    avoid_unclean_buscos: bool = Field(
        default=True,
        description="Ignore BUSCOs flagged unclean by prior paralog filtering or decontamination when computing medians.",
    )
    include_duplicated: bool = Field(
        default=False,
        description="Evaluate duplicated BUSCO copies and store per-copy paralog evidence without changing default family-level outputs.",
    )

    @model_validator(mode="after")
    def _validate_library(self) -> "ParalogRemovalPayload":
        if not self.library_id and not self.library_name:
            raise ValueError("Either library_id or library_name must be provided.")
        if not self.accessions and not self.targets:
            raise ValueError("Provide accessions, targets, or selector-derived accessions for paralog removal.")
        if self.mode == "percent" and self.percentile is None:
            raise ValueError("--percentile is required when --mode percent is selected.")
        if self.mode != "percent" and self.percentile is not None:
            raise ValueError("--percentile is only valid with --mode percent.")
        if self.mode == "bitscore" and self.bitscore_threshold is None:
            raise ValueError("--bitscore-threshold is required when --mode bitscore is selected.")
        if self.mode != "bitscore" and self.bitscore_threshold is not None:
            raise ValueError("--bitscore-threshold is only valid with --mode bitscore.")
        if not self.accessions and self.targets:
            self.accessions = list(self.targets)
        if not self.targets and self.accessions:
            self.targets = list(self.accessions)
        if self.targets:
            accessions_set = {str(acc) for acc in self.accessions}
            missing = [str(acc) for acc in self.targets if str(acc) not in accessions_set]
            if missing:
                raise ValueError(
                    "Paralog-removal targets must be a subset of accessions used for median calculation. "
                    f"Missing from accessions: {', '.join(missing)}"
                )
        return self
    
class BatchBuscoTaskPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    lineage: str
    format: Literal["auto", "protein", "genome", "nucleotide"] = "auto"
    pipeline: Optional[Literal["auto", "miniprot", "metaeuk", "augustus"]] = "auto"
    proteome_profile: Optional[str] = None
    prefer_proteome_profile: Optional[str] = None
    isoforms_cleaned: Optional[bool] = None
    raw_proteome: bool = False
    augustus_evalue: Optional[float] = Field(default=None, ge=0.0)
    augustus_limit: Optional[int] = Field(default=None, ge=1)
    augustus_long: Optional[bool] = None
    augustus_species: Optional[str] = None
    augustus_parameters: Optional[str] = None
    metaeuk_parameters: Optional[str] = None
    metaeuk_rerun_parameters: Optional[str] = None
    miniprot_parameters: Optional[str] = None
    output_dir: Optional[Path] = None
    busco_path: Optional[Path] = None
    force: bool = False
    keep_miniprot_ref_file: bool = Field(
        default=False,
        description="Retain miniprot ref.mpi file in BUSCO output.",
    )
    max_concurrent: PositiveInt = Field(default=1)
    busco_lib_wait_seconds: int = Field(default=0, ge=0)
    busco_lib_retries: int = Field(default=0, ge=0)


class DecontaminationPayload(TaskPayload):
    targets: List[str] = Field(default_factory=list, description="Target accessions to decontaminate.")
    accessions: List[str] = Field(default_factory=list, description="Alias for targets (legacy).")
    refs: List[str] = Field(default_factory=list, description="Explicit reference accessions or clade names.")
    ref_accessions: List[str] = Field(default_factory=list, description="Alias for refs (legacy).")
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    ref_taxid: Optional[int] = Field(default=None, ge=1)
    ref_clade: Optional[str] = None
    ref_rule_rank: Optional[str] = None
    ref_rule_quantity: Optional[int] = Field(default=None, ge=1)
    ref_select_clade: Optional[str] = None
    ref_select_rank: Optional[str] = None
    ref_select_top: Optional[int] = Field(default=None, ge=1)
    rank: str = Field(default="phylum")
    off_clade_fraction: float = Field(default=0.1, ge=0.0, le=1.0)
    min_buscos: PositiveInt = Field(default=20)
    min_identity: float = Field(default=0.0, ge=0.0, le=100.0)
    min_coverage: float = Field(default=0.0, ge=0.0, le=100.0)
    min_delta_bitscore: float = Field(default=0.0, ge=0.0)
    min_bitscore: float = Field(default=0.0, ge=0.0)
    max_evalue: Optional[float] = Field(default=None, ge=0.0)
    min_hits: PositiveInt = Field(default=1)
    hit_window: PositiveInt = Field(
        default=1,
        description="Number of top hits considered when checking min_hits.",
    )
    config_path: Optional[Path] = None
    run_label: Optional[str] = None
    use_paralog_filtered_refs: bool = False
    include_duplicated: bool = Field(
        default=False,
        description="Evaluate duplicated BUSCO copies and store per-copy decontamination votes without changing default assembly classification.",
    )
    allow_same_species: bool = False
    allow_sparse_references: bool = False
    report_path: Optional[Path] = None
    force: bool = False
    max_concurrent: Optional[PositiveInt] = None
    busco_wait_seconds: int = Field(default=0, ge=0)
    busco_retries: int = Field(default=0, ge=0)


class InternalDecontaminationPayload(TaskPayload):
    targets: List[str] = Field(default_factory=list, description="Target accessions to screen internally.")
    accessions: List[str] = Field(default_factory=list, description="Alias for targets (legacy).")
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    rank: str = Field(default="phylum")
    hit_window: PositiveInt = Field(
        default=8,
        description="Number of top hits to inspect after removing the self-hit.",
    )
    p_value_threshold: float = Field(default=0.05, ge=0.0, le=1.0)
    off_clade_fraction: float = Field(
        default=0.05,
        ge=0.0,
        le=1.0,
        description="Assembly-level contamination threshold: fraction of BUSCOs called outside.",
    )
    min_buscos: PositiveInt = Field(default=20)
    min_identity: float = Field(default=0.0, ge=0.0, le=100.0)
    min_coverage: float = Field(default=0.0, ge=0.0, le=100.0)
    min_alignment_length: Optional[PositiveInt] = Field(default=None)
    min_bitscore: float = Field(default=0.0, ge=0.0)
    max_evalue: Optional[float] = Field(default=None, ge=0.0)
    max_target_seqs: Optional[PositiveInt] = Field(default=None)
    blast_program: Optional[Literal["blastp", "blastn"]] = None
    blast_db_type: Optional[Literal["prot", "nucl"]] = None
    save_blast_output: Optional[Path] = Field(
        default=None,
        description="Directory to save raw BLAST outfmt 6 outputs for reuse.",
    )
    reuse_blast_results: Optional[Path] = Field(
        default=None,
        description="Directory containing previously saved raw BLAST outfmt 6 outputs.",
    )
    external_blast_db_path: Optional[Path] = Field(
        default=None,
        description="Path to external BLAST database (e.g., local NR) for follow-up checks.",
    )
    external_blast_db_type: Optional[Literal["prot", "nucl"]] = None
    external_blast_program: Optional[Literal["blastp", "blastn", "diamond"]] = None
    external_blast_output_dir: Optional[Path] = Field(
        default=None,
        description="Directory to save external BLAST outputs keyed by accession/family.",
    )
    external_reuse_blast_results: Optional[Path] = Field(
        default=None,
        description="Directory containing external BLAST outputs to reuse.",
    )
    external_max_target_seqs: Optional[PositiveInt] = Field(
        default=None,
        description="Max target sequences for external BLAST (defaults to hit_window + 1).",
    )
    config_path: Optional[Path] = None
    run_label: Optional[str] = None
    use_paralog_filtered_buscos: bool = False
    report_path: Optional[Path] = None
    force: bool = False
    max_concurrent: Optional[PositiveInt] = None
    busco_wait_seconds: int = Field(default=0, ge=0)
    busco_retries: int = Field(default=0, ge=0)

    @model_validator(mode="after")
    def _validate_inputs(self) -> "InternalDecontaminationPayload":
        if not (self.targets or self.accessions) and not self.config_path:
            raise ValueError("At least one accession (or a config_path with members) is required for internal decontamination.")
        if not self.library_id and not self.library_name:
            raise ValueError("Either library_id or library_name must be provided.")
        return self


class ExternalDecontaminationCheckPayload(TaskPayload):
    run_id: str
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    blast_db_path: Optional[Path] = Field(
        default=None,
        description="Path to external BLAST database (e.g., local NR).",
    )
    blast_db_type: Optional[Literal["prot", "nucl"]] = None
    blast_program: Optional[Literal["blastp", "blastn", "diamond"]] = None
    output_dir: Optional[Path] = Field(
        default=None,
        description="Directory to write external BLAST outputs.",
    )
    reuse_blast_results: Optional[Path] = Field(
        default=None,
        description="Directory containing previously saved external BLAST outputs.",
    )
    max_target_seqs: Optional[PositiveInt] = Field(default=None)
    hit_window: Optional[PositiveInt] = Field(default=None)
    force: bool = False
    max_concurrent: Optional[PositiveInt] = None

    @model_validator(mode="after")
    def _validate_inputs(self) -> "ExternalDecontaminationCheckPayload":
        if not self.run_id:
            raise ValueError("run_id is required for external decontamination checks.")
        return self


class ExternalDecontaminationApplyPayload(TaskPayload):
    source_run_id: str
    new_run_id: Optional[str] = None
    run_label: Optional[str] = None
    report_path: Optional[Path] = None
    external_blast_output_dir: Optional[Path] = None
    external_reuse_blast_results: Optional[Path] = None

    @model_validator(mode="after")
    def _validate_inputs(self) -> "ExternalDecontaminationApplyPayload":
        if not self.source_run_id:
            raise ValueError("source_run_id is required for external apply.")
        if not self.external_blast_output_dir and not self.external_reuse_blast_results:
            # Allow fallback to run params; do not hard fail here.
            return self
        return self


class ConstructBuscoBlastDBPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    busco_library_id: Optional[int] = Field(default=None, ge=1)
    target_library_id: Optional[int] = Field(default=None, ge=1, description="Library to use for paralog filtering if enabled.")
    family_ids: List[str] = Field(default_factory=list)
    output_path: Optional[Path] = None
    use_paralog_filtered: bool = False
    force: bool = False
    id_mode: Optional[Literal["legacy", "internal"]] = Field(default=None)
    id_map_path: Optional[Path] = None
    db_type: Optional[Literal["prot", "nucl"]] = None

    @model_validator(mode="after")
    def _validate(self) -> "ConstructBuscoBlastDBPayload":
        if not self.accessions:
            raise ValueError("At least one accession is required to build a BUSCO BLAST DB.")
        return self


class VerifyAssemblyPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    root: Optional[str] = None
    taxid: Optional[int] = Field(default=None, ge=1)
    clade: Optional[str] = None
    all: bool = False
    downloaded_only: bool = False
    reacquire: bool = False
    discover: bool = False
    discover_protein: bool = False
    tidy: bool = False
    organise: bool = False
    organise_check_only: bool = False
    split_isolated_proteomes: bool = False
    report: Optional[Path] = None
    report_root: Optional[Path] = None
    primary_only: bool = False
    after: Optional[str] = None
    before: Optional[str] = None
    level: Optional[str] = None
    filters: Optional[List[str]] = None
    ranks: Optional[List[str]] = None
    quantities: Optional[List[str]] = None
    repair: bool = False
    clean_isoforms: bool = False
    skip_clean_isoforms: bool = False
    clean_skip_gff: bool = False
    clean_skip_cdhit: bool = False
    clean_gff_priority: bool = False
    clean_revert_from_archive: bool = False
    clean_cdhit_identity: Optional[float] = Field(default=None, gt=0.0, le=1.0)
    clean_max_concurrent: PositiveInt = Field(default=1, le=8)
    clean_threads_per_job: PositiveInt = Field(default=1)


VerifyDownloadsPayload = VerifyAssemblyPayload
class PrepareProteomePayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    taxid: Optional[int] = Field(default=None, ge=1)
    downloaded_only: bool = True
    after: Optional[str] = None
    before: Optional[str] = None
    level: Optional[str] = None
    primary_only: bool = False
    profile_name: str = "clean_default"
    input_profile: str = "raw"
    skip_gff: bool = False
    skip_cdhit: bool = False
    gff_priority: bool = False
    max_concurrent: PositiveInt = Field(default=1, le=8)
    threads_per_job: PositiveInt = Field(default=1)
    cdhit_identity: Optional[float] = Field(default=None, gt=0.0, le=1.0)
    replace_existing: bool = False
    set_default: bool = True


class VerifyBuscoPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    root: Optional[str] = None
    taxid: Optional[int] = Field(default=None, ge=1)
    clade: Optional[str] = None
    all: bool = False
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    run_id: Optional[int] = Field(default=None, ge=1)
    downloaded_only: bool = False
    discover: bool = False
    queue_missing: bool = False
    reingest: bool = False
    reingest_all: bool = False
    repair: bool = False
    stale_missing: bool = True
    restore_found: bool = True
    reassign_primary: bool = True
    primary_only: bool = False
    after: Optional[str] = None
    before: Optional[str] = None
    level: Optional[str] = None
    report: Optional[Path] = None
    report_root: Optional[Path] = None
    filters: Optional[List[str]] = None
    ranks: Optional[List[str]] = None
    quantities: Optional[List[str]] = None

    @model_validator(mode="after")
    def _validate_inputs(self) -> "VerifyBuscoPayload":
        if self.reingest and self.reingest_all:
            raise ValueError("Choose either reingest or reingest_all, not both.")
        if (self.reingest or self.reingest_all) and not self.repair:
            raise ValueError("reingest and reingest_all require repair=True.")
        if self.run_id is not None and (self.library_id is None and self.library_name is None):
            return self
        if self.all:
            return self
        if self.accessions or self.taxid is not None or self.clade:
            return self
        if not self.library_id and not self.library_name:
            raise ValueError("Either library_id or library_name must be provided for BUSCO verification.")
        return self


class VerifyLibrariesPayload(TaskPayload):
    root: Optional[str] = None
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    ref_accessions: List[str] = Field(default_factory=list)
    all: bool = False
    repair: bool = False
    report: Optional[Path] = None
    report_root: Optional[Path] = None

    @model_validator(mode="after")
    def _validate_inputs(self) -> "VerifyLibrariesPayload":
        if self.all:
            return self
        if self.library_id or self.library_name or self.ref_accessions:
            return self
        raise ValueError("Provide --all, a library identifier, or ref_accessions for library verification.")


class VerifyOrthofinderPayload(TaskPayload):
    root: Optional[str] = None
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    ref_accessions: List[str] = Field(default_factory=list)
    all: bool = False
    repair: bool = False
    report: Optional[Path] = None
    report_root: Optional[Path] = None

    @model_validator(mode="after")
    def _validate_inputs(self) -> "VerifyOrthofinderPayload":
        if self.all:
            return self
        if self.library_id or self.library_name or self.ref_accessions:
            return self
        raise ValueError("Provide --all, a library identifier, or ref_accessions for OrthoFinder verification.")


class VerifyPayload(TaskPayload):
    accessions: List[str] = Field(default_factory=list)
    root: Optional[str] = None
    taxid: Optional[int] = Field(default=None, ge=1)
    clade: Optional[str] = None
    all: bool = False
    downloaded_only: bool = False
    primary_only: bool = False
    after: Optional[str] = None
    before: Optional[str] = None
    level: Optional[str] = None
    filters: Optional[List[str]] = None
    ranks: Optional[List[str]] = None
    quantities: Optional[List[str]] = None
    library_id: Optional[int] = Field(default=None, ge=1)
    library_name: Optional[str] = None
    run_id: Optional[int] = Field(default=None, ge=1)
    ref_accessions: List[str] = Field(default_factory=list)
    repair: bool = False
    include_assembly: bool = True
    include_libraries: bool = True
    include_busco: bool = True
    include_orthofinder: bool = True
    report: Optional[Path] = None
    report_root: Optional[Path] = None


class SplitRecordsPayload(TaskPayload):
    accession: str
    folder: Optional[Path] = None
    split_isolated_proteomes: bool = False
    check_only: bool = False
    report: Optional[Path] = None
    @model_validator(mode="after")
    def _validate_inputs(self) -> "DecontaminationPayload":
        # Allow accessions to come from JSON config when config_path is provided
        if not (self.targets or self.accessions) and not self.config_path:
            raise ValueError("At least one accession (or a config_path with members) is required for decontamination.")
        if not self.library_id and not self.library_name:
            raise ValueError("Either library_id or library_name must be provided.")
        return self
