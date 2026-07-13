"""Definitions and JSON helpers for database-backed PhyloODB variables."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .variable_kinds import (
    SYSTEM_VARIABLE_EXACT,
    VARIABLE_KIND_ASSEMBLIES,
    VARIABLE_KIND_BUSCO_RUNS,
    VARIABLE_KIND_ENV,
    normalize_variable_kind,
)

VARIABLE_JSON_FORMAT = "phyloodb.variables.v1"

VARIABLE_SECTION_BY_KIND = {
    VARIABLE_KIND_ENV: "environment",
    VARIABLE_KIND_ASSEMBLIES: "assemblies",
    VARIABLE_KIND_BUSCO_RUNS: "busco_runs",
}

VARIABLE_KIND_BY_SECTION = {section: kind for kind, section in VARIABLE_SECTION_BY_KIND.items()}

VARIABLE_JSON_SECTIONS = tuple(VARIABLE_KIND_BY_SECTION)

_DYNAMIC_ENV_DEFINITIONS: dict[str, dict[str, Any]] = {
    "DEFAULT_THREADS_<TASK>": {
        "type": "integer-or-null",
        "default": None,
        "description": "Default required thread count for a registered task; null falls back to the registry default.",
    },
    "LAST_<TASK>": {
        "type": "integer",
        "default": None,
        "description": "Automatically maintained pointer to the most recent task id for a task family.",
    },
    "ACTIVE_DECONT_RUN_<library_id>": {
        "type": "object",
        "default": None,
        "description": "Automatically maintained pointer to the active external decontamination run for a library.",
    },
    "ACTIVE_INTERNAL_DECONT_RUN_<library_id>": {
        "type": "object",
        "default": None,
        "description": "Automatically maintained pointer to the active internal decontamination run for a library.",
    },
    "ACTIVE_PARALOG_RUN_<library_id>": {
        "type": "object",
        "default": None,
        "description": "Automatically maintained pointer to the active paralog-removal run for a library.",
    },
}

_ENV_DEFINITION_OVERRIDES: dict[str, dict[str, Any]] = {
    "BLASTN_PATH": {"type": "path-or-command", "default": "blastn", "description": "BLASTN executable used by nucleotide decontamination checks."},
    "BLASTP_PATH": {"type": "path-or-command", "default": "blastp", "description": "BLASTP executable used by protein searches and decontamination checks."},
    "BLOCKED_TASK_QUEUE_POLLING_TIME": {"type": "number", "default": 2, "description": "Seconds between blocked task queue polling checks."},
    "BUSCO_AUGUSTUS_EVALUE": {"type": "number-or-null", "default": None, "description": "Optional Augustus e-value override passed to BUSCO."},
    "BUSCO_AUGUSTUS_LIMIT": {"type": "integer", "default": 3, "description": "Augustus retraining limit passed to BUSCO."},
    "BUSCO_AUGUSTUS_LONG": {"type": "boolean", "default": False, "description": "Whether BUSCO Augustus mode should use the long option."},
    "BUSCO_AUGUSTUS_PARAMETERS": {"type": "string-or-null", "default": None, "description": "Additional raw Augustus parameters for BUSCO."},
    "BUSCO_AUGUSTUS_SPECIES": {"type": "string-or-null", "default": None, "description": "Optional Augustus species model passed to BUSCO."},
    "BUSCO_BINARIES_PATH": {"type": "path-or-command", "default": "busco", "description": "BUSCO executable or command name."},
    "BUSCO_LINEAGE_DIR": {"type": "path", "default": "/path/to/libraries/lineages", "description": "Directory containing BUSCO lineage datasets."},
    "BUSCO_METAEUK_PARAMETERS": {"type": "string-or-null", "default": None, "description": "Additional raw MetaEuk parameters for BUSCO."},
    "BUSCO_METAEUK_RERUN_PARAMETERS": {"type": "string-or-null", "default": None, "description": "Additional raw MetaEuk rerun parameters for BUSCO."},
    "BUSCO_MINIPROT_KEEP_REF_FILE": {"type": "boolean", "default": False, "description": "Whether BUSCO Miniprot reference files should be retained."},
    "BUSCO_MINIPROT_PARAMETERS": {"type": "string-or-null", "default": None, "description": "Additional raw Miniprot parameters for BUSCO."},
    "CACHE_DIR": {"type": "path", "default": "/path/to/cache", "description": "Shared cache directory for reusable intermediate files."},
    "DAEMON_MAX_THREADS": {"type": "integer", "default": 1, "description": "Maximum total threads the task daemon may allocate."},
    "DAEMON_PROCESS_POLLING_TIME": {"type": "number", "default": 2, "description": "Seconds between daemon task polling checks."},
    "DEFAULT_BUSCO_FORMAT": {"type": "string", "default": "protein", "description": "Default BUSCO input format when one is not supplied."},
    "DEFAULT_BUSCO_PIPELINE": {"type": "string", "default": "miniprot", "description": "Default BUSCO pipeline when one is not supplied."},
    "DEFAULT_PROTEOME_CDHIT_IDENTITY": {"type": "number", "default": 0.96, "description": "Default CD-HIT identity for automatic proteome preparation."},
    "DEFAULT_PROTEOME_CLEAN_ISOFORMS": {"type": "boolean", "default": True, "description": "Whether automatic proteome preparation should clean isoforms."},
    "DEFAULT_PROTEOME_GFF_PRIORITY": {"type": "boolean", "default": False, "description": "Whether GFF-supported isoform choices should be prioritized when preparing proteomes."},
    "DEFAULT_PROTEOME_INPUT_PROFILE": {"type": "string", "default": "raw", "description": "Default source proteome profile for automatic proteome preparation."},
    "DEFAULT_PROTEOME_MAX_CONCURRENT": {"type": "integer", "default": 1, "description": "Default maximum concurrent proteome preparation jobs."},
    "DEFAULT_PROTEOME_SET_DEFAULT": {"type": "boolean", "default": True, "description": "Whether prepared proteomes become the accession default profile."},
    "DEFAULT_PROTEOME_THREADS_PER_JOB": {"type": "integer", "default": 1, "description": "Default threads per proteome preparation job."},
    "DEFAULT_PROTEOME_USE_CDHIT": {"type": "boolean", "default": False, "description": "Whether automatic proteome preparation uses CD-HIT by default."},
    "DEFAULT_PROTEOME_USE_GFF": {"type": "boolean", "default": True, "description": "Whether automatic proteome preparation uses associated GFF files by default."},
    "DIAMOND_PATH": {"type": "path-or-command", "default": "diamond", "description": "DIAMOND executable used by translated decontamination checks."},
    "EMAIL": {"type": "string-or-null", "default": None, "description": "Email address used for NCBI Entrez requests."},
    "EXPORTS_DIR": {"type": "path", "default": "/path/to/exports", "description": "Default export output directory."},
    "GENOME_DIR": {"type": "path", "default": "/path/to/genomes", "description": "Directory used to store downloaded genome/proteome files."},
    "IQTREE_FLAGS": {"type": "string-or-null", "default": None, "description": "Additional IQ-TREE flags for tree-building tasks."},
    "IQTREE_PATH": {"type": "path-or-command", "default": "iqtree2", "description": "IQ-TREE executable used by tree-building tasks."},
    "LAST": {"type": "integer", "default": None, "description": "Automatically maintained pointer to the most recent task id."},
    "LIBRARIES_DIR": {"type": "path", "default": "/path/to/libraries", "description": "Directory used to store PhyloODB libraries."},
    "LIST_BUSCO_GRADIENT": {"type": "boolean-or-null", "default": None, "description": "Legacy/default toggle for BUSCO color gradients in list output."},
    "LIST_BUSCO_GRADIENT_NEG": {"type": "boolean-or-null", "default": None, "description": "Toggle for negative BUSCO color gradients in list output."},
    "LIST_BUSCO_GRADIENT_POS": {"type": "boolean-or-null", "default": None, "description": "Toggle for positive BUSCO color gradients in list output."},
    "LIST_BUSCO_NEG_STOPS": {"type": "array-or-null", "default": None, "description": "Color stops for negative BUSCO gradients."},
    "LIST_BUSCO_POS_STOPS": {"type": "array-or-null", "default": None, "description": "Color stops for positive BUSCO gradients."},
    "LIST_BUSCO_STEEP_MAX": {"type": "number", "default": 20.0, "description": "Maximum delta used for steep BUSCO gradient scaling."},
    "LIST_BUSCO_STEEP_STOPS": {"type": "array-or-null", "default": None, "description": "Color stops for steep BUSCO gradients."},
    "LIST_GROUP_COLORS": {"type": "object-or-null", "default": None, "description": "Group color overrides for list output."},
    "LIST_USE_COLOR": {"type": "boolean-or-null", "default": None, "description": "Whether list output should use color when supported."},
    "LOG_BACKUPS": {"type": "integer", "default": 5, "description": "Number of rotated log files to keep."},
    "LOG_DIR": {"type": "path", "default": "/path/to/logs", "description": "Directory used for PhyloODB log files."},
    "LOG_FILE": {"type": "path", "default": "/path/to/logs/phyloodb.log", "description": "Main PhyloODB log file."},
    "LOG_FORMAT": {"type": "string-or-null", "default": None, "description": "Optional logging format override."},
    "LOG_HIDE_CATEGORIES_CONSOLE": {"type": "array-or-string-or-null", "default": None, "description": "Log categories hidden from console output."},
    "LOG_HIDE_CATEGORIES_FILE": {"type": "array-or-string-or-null", "default": None, "description": "Log categories hidden from file output."},
    "LOG_LEVEL": {"type": "string", "default": "DEBUG", "description": "Default logging level."},
    "LOG_MAX_BYTES": {"type": "integer", "default": 5242880, "description": "Maximum size of a log file before rotation."},
    "LOG_TO_CONSOLE": {"type": "boolean", "default": False, "description": "Whether normal logging should also be emitted to the console."},
    "LOG_USE_COLOR": {"type": "boolean", "default": True, "description": "Whether logging output should use color when supported."},
    "MAFFT_FLAGS": {"type": "string-or-null", "default": None, "description": "Additional MAFFT flags for alignment tasks."},
    "MAFFT_PATH": {"type": "path-or-command", "default": "mafft", "description": "MAFFT executable used by alignment tasks."},
    "MAKEBLASTDB_PATH": {"type": "path-or-command", "default": "makeblastdb", "description": "makeblastdb executable used when building BLAST databases."},
    "MISC_DIR": {"type": "path", "default": "/path/to/misc", "description": "Directory for miscellaneous derived artifacts."},
    "NCBI_API_KEY": {"type": "string-or-null", "default": None, "description": "Optional NCBI API key for Entrez requests."},
    "ORTHOFINDER_BINARIES_PATH": {"type": "path-or-command", "default": "orthofinder", "description": "OrthoFinder executable or command name."},
    "ORTHOFINDER_OUTPUT_DIR": {"type": "path", "default": "/path/to/orthofinder", "description": "Directory used for OrthoFinder outputs."},
    "REPORTS_DIR": {"type": "path", "default": "/path/to/reports", "description": "Directory used for task reports."},
    "SELECTOR_BUSCO_BUCKETS": {"type": "array", "default": [[0, -5], [10, -4], [20, -3], [30, -2], [40, -1], [50, 1], [60, 2], [70, 3], [80, 4], [90, 5]], "description": "BUSCO score buckets used by default selector scoring."},
    "SELECTOR_DEFAULT_DOWNLOADED_ONLY": {"type": "boolean", "default": False, "description": "Default selector setting for requiring downloaded assemblies."},
    "SELECTOR_DEFAULT_PRIMARY_ONLY": {"type": "boolean", "default": False, "description": "Default selector setting for primary assemblies only."},
    "SELECTOR_DEFAULT_PROTEIN_ONLY": {"type": "boolean", "default": False, "description": "Default selector setting for requiring protein data."},
    "SELECTOR_DEFAULT_STATUS_MIN": {"type": "integer-or-null", "default": None, "description": "Default minimum assembly status used by selectors."},
    "SELECTOR_DEFAULT_USE_BUSCO": {"type": "boolean", "default": False, "description": "Default selector setting for using BUSCO scores."},
    "SELECTOR_SCORE_ORDER": {"type": "array", "default": ["busco", "refseq", "level", "n50", "date", "accession"], "description": "Default selector ranking criteria."},
    "SET_MAX_THREADS_ON_START": {"type": "boolean", "default": True, "description": "Refresh detected daemon and task thread defaults at run or daemon startup."},
}


def environment_variable_definitions() -> dict[str, dict[str, Any]]:
    """Return definitions for all known static and dynamic environment variables."""

    definitions: dict[str, dict[str, Any]] = {}
    for name in sorted(SYSTEM_VARIABLE_EXACT):
        definitions[name] = _ENV_DEFINITION_OVERRIDES.get(
            name,
            {
                "type": "json",
                "default": None,
                "description": f"PhyloODB environment variable {name}.",
            },
        )
    try:
        from .thread_defaults import computed_task_thread_defaults, detect_available_threads

        for name, default in computed_task_thread_defaults(detect_available_threads()).items():
            definitions[name] = {
                "type": "integer-or-null",
                "default": default,
                "description": "Default required thread count for this registered task; null falls back to the registry default.",
            }
    except Exception:
        pass
    definitions.update(_DYNAMIC_ENV_DEFINITIONS)
    return definitions


def build_variables_json_document(
    records: Mapping[str, Mapping[str, Any]],
    *,
    source: str | None = None,
    kind_filter: str | None = None,
) -> dict[str, Any]:
    """Build the public kinded variable JSON document from stored records."""

    normalized_filter = normalize_variable_kind(kind_filter) if kind_filter else None
    document: dict[str, Any] = {
        "_metadata": {
            "format": VARIABLE_JSON_FORMAT,
        },
        "_definitions": {
            "environment": environment_variable_definitions(),
            "assemblies": {
                "_example": {
                    "type": "array[string]",
                    "description": "Named accession lists reusable as @NAME selectors.",
                }
            },
            "busco_runs": {
                "_example": {
                    "type": "array[integer]",
                    "description": "Named BUSCO run-id lists reusable by run-id aware commands.",
                }
            },
        },
    }
    if source:
        document["_metadata"]["source"] = Path(source).name

    included_sections = []
    for kind, section in VARIABLE_SECTION_BY_KIND.items():
        if normalized_filter and kind != normalized_filter:
            continue
        document[section] = {}
        included_sections.append(section)

    for name, record in sorted(records.items()):
        kind = normalize_variable_kind(record.get("kind")) or VARIABLE_KIND_ENV
        if normalized_filter and kind != normalized_filter:
            continue
        section = VARIABLE_SECTION_BY_KIND.get(kind, "environment")
        if section not in document:
            continue
        document[section][str(name)] = record.get("value")

    if normalized_filter:
        document["_definitions"] = {
            section: document["_definitions"][section]
            for section in included_sections
            if section in document["_definitions"]
        }
    return document


def validate_variables_json_document(payload: Any) -> tuple[dict[str, Any], dict[str, str]]:
    """Validate an import document and return values plus explicit per-key kinds."""

    if not isinstance(payload, dict):
        raise ValueError("Variable JSON must be a top-level object.")
    unknown = sorted(key for key in payload if not str(key).startswith("_") and key not in VARIABLE_JSON_SECTIONS)
    if unknown:
        raise ValueError(f"Unknown variable JSON section(s): {', '.join(unknown)}.")

    values: dict[str, Any] = {}
    kinds: dict[str, str] = {}
    for section, kind in VARIABLE_KIND_BY_SECTION.items():
        raw_section = payload.get(section)
        if raw_section is None:
            continue
        if not isinstance(raw_section, dict):
            raise ValueError(f"Variable JSON section '{section}' must be an object.")
        for raw_name, value in raw_section.items():
            name = str(raw_name or "").strip()
            if section == "assemblies" and name.startswith("@"):
                name = name[1:].strip()
            name = name.upper()
            if not name:
                raise ValueError(f"Variable JSON section '{section}' contains an empty variable name.")
            if name in kinds:
                raise ValueError(f"Variable '{name}' appears in more than one JSON variable section.")
            if section == "assemblies":
                if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
                    raise ValueError(f"Assembly variable '{raw_name}' must be an array of strings.")
            elif section == "busco_runs":
                if not isinstance(value, list) or not all(isinstance(item, int) or (isinstance(item, str) and item.strip().isdigit()) for item in value):
                    raise ValueError(f"BUSCO run variable '{raw_name}' must be an array of integers or numeric strings.")
                value = [int(item) for item in value]
            values[name] = value
            kinds[name] = kind
    return values, kinds


def example_variables_json_document() -> dict[str, Any]:
    """Return an importable example document for docs/variables.example.json."""

    try:
        from .thread_defaults import detect_available_threads, computed_task_thread_defaults

        detected_threads = detect_available_threads()
        computed_threads = computed_task_thread_defaults(detected_threads)
    except Exception:
        detected_threads = 1
        computed_threads = {}
    env_values = {
        name: definition.get("default")
        for name, definition in environment_variable_definitions().items()
        if "<" not in name and ">" not in name
    }
    env_values["DAEMON_MAX_THREADS"] = detected_threads
    env_values.update(computed_threads)
    return {
        "_metadata": {
            "format": VARIABLE_JSON_FORMAT,
            "source": "docs/variables.example.json",
        },
        "_definitions": {
            "environment": environment_variable_definitions(),
            "assemblies": {
                "_example": {
                    "type": "array[string]",
                    "description": "Named accession lists reusable as @NAME selectors.",
                }
            },
            "busco_runs": {
                "_example": {
                    "type": "array[integer]",
                    "description": "Named BUSCO run-id lists reusable by run-id aware commands.",
                }
            },
        },
        "environment": env_values,
        "assemblies": {
            "EXAMPLE_ASSEMBLIES": ["GCF_000001405.40"],
        },
        "busco_runs": {
            "EXAMPLE_BUSCO_RUNS": [101],
        },
    }
