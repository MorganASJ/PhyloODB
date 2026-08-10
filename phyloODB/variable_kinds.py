"""Shared helpers for Environment_Variables typing."""
from __future__ import annotations

import re
from typing import Any, Mapping

from .accession_utils import canonicalize_accession


VARIABLE_KIND_ENV = "env"
VARIABLE_KIND_ASSEMBLIES = "assemblies"
VARIABLE_KIND_BUSCO_RUNS = "busco_runs"
VARIABLE_KINDS = {VARIABLE_KIND_ENV, VARIABLE_KIND_ASSEMBLIES, VARIABLE_KIND_BUSCO_RUNS}

VARIABLE_KIND_ALIASES = {
    "env": VARIABLE_KIND_ENV,
    "environment": VARIABLE_KIND_ENV,
    "environmental": VARIABLE_KIND_ENV,
    "enviornmental": VARIABLE_KIND_ENV,
    "assemblies": VARIABLE_KIND_ASSEMBLIES,
    "assembly": VARIABLE_KIND_ASSEMBLIES,
    "accessions": VARIABLE_KIND_ASSEMBLIES,
    "accession": VARIABLE_KIND_ASSEMBLIES,
    "busco-runs": VARIABLE_KIND_BUSCO_RUNS,
    "busco-run": VARIABLE_KIND_BUSCO_RUNS,
    "busco_runs": VARIABLE_KIND_BUSCO_RUNS,
    "busco_run": VARIABLE_KIND_BUSCO_RUNS,
    "runs": VARIABLE_KIND_BUSCO_RUNS,
}

SYSTEM_VARIABLE_EXACT = {
    "BLASTN_PATH",
    "BLASTP_PATH",
    "BLOCKED_TASK_QUEUE_POLLING_TIME",
    "BUSCO_AUGUSTUS_EVALUE",
    "BUSCO_AUGUSTUS_LIMIT",
    "BUSCO_AUGUSTUS_LONG",
    "BUSCO_AUGUSTUS_PARAMETERS",
    "BUSCO_AUGUSTUS_SPECIES",
    "BUSCO_BINARIES_PATH",
    "BUSCO_LINEAGE_DIR",
    "BUSCO_METAEUK_PARAMETERS",
    "BUSCO_METAEUK_RERUN_PARAMETERS",
    "BUSCO_MINIPROT_KEEP_REF_FILE",
    "BUSCO_MINIPROT_PARAMETERS",
    "CACHE_DIR",
    "PROJECT_PERMISSION_MODE",
    "SCRATCH_DIR",
    "SHARED_GROUP",
    "DAEMON_MAX_THREADS",
    "DAEMON_PROCESS_POLLING_TIME",
    "DEFAULT_BUSCO_FORMAT",
    "DEFAULT_BUSCO_PIPELINE",
    "DEFAULT_PROTEOME_CDHIT_IDENTITY",
    "DEFAULT_PROTEOME_CLEAN_ISOFORMS",
    "DEFAULT_PROTEOME_GFF_PRIORITY",
    "DEFAULT_PROTEOME_INPUT_PROFILE",
    "DEFAULT_PROTEOME_MAX_CONCURRENT",
    "DEFAULT_PROTEOME_SET_DEFAULT",
    "DEFAULT_PROTEOME_THREADS_PER_JOB",
    "DEFAULT_PROTEOME_USE_CDHIT",
    "DEFAULT_PROTEOME_USE_GFF",
    "EMAIL",
    "EXPORTS_DIR",
    "GENOME_DIR",
    "DIAMOND_PATH",
    "IQTREE_FLAGS",
    "IQTREE_PATH",
    "LAST",
    "LIBRARIES_DIR",
    "LIST_BUSCO_GRADIENT",
    "LIST_BUSCO_GRADIENT_NEG",
    "LIST_BUSCO_GRADIENT_POS",
    "LIST_BUSCO_NEG_STOPS",
    "LIST_BUSCO_POS_STOPS",
    "LIST_BUSCO_STEEP_MAX",
    "LIST_BUSCO_STEEP_STOPS",
    "LIST_GROUP_COLORS",
    "LIST_USE_COLOR",
    "LOG_BACKUPS",
    "LOG_DIR",
    "LOG_FILE",
    "LOG_FORMAT",
    "LOG_HIDE_CATEGORIES_CONSOLE",
    "LOG_HIDE_CATEGORIES_FILE",
    "LOG_LEVEL",
    "LOG_MAX_BYTES",
    "LOG_TO_CONSOLE",
    "LOG_USE_COLOR",
    "MAFFT_FLAGS",
    "MAFFT_PATH",
    "MAKEBLASTDB_PATH",
    "MISC_DIR",
    "NCBI_API_KEY",
    "ORTHOFINDER_BINARIES_PATH",
    "ORTHOFINDER_OUTPUT_DIR",
    "REPORTS_DIR",
    "SELECTOR_BUSCO_BUCKETS",
    "SELECTOR_DEFAULT_DOWNLOADED_ONLY",
    "SELECTOR_DEFAULT_PRIMARY_ONLY",
    "SELECTOR_DEFAULT_PROTEIN_ONLY",
    "SELECTOR_DEFAULT_STATUS_MIN",
    "SELECTOR_DEFAULT_USE_BUSCO",
    "SELECTOR_SCORE_ORDER",
    "SET_MAX_THREADS_ON_START",
}

SYSTEM_VARIABLE_PREFIXES = (
    "ACTIVE_DECONT_RUN_",
    "ACTIVE_INTERNAL_DECONT_RUN_",
    "ACTIVE_PARALOG_RUN_",
    "DEFAULT_THREADS_",
    "LAST_",
    "LIST_",
    "SELECTOR_DEFAULT_",
)

_ACCESSIONISH_TOKEN_RE = re.compile(r"^[A-Za-z0-9_.@-]+$")


def normalize_variable_kind(raw_kind: Any | None) -> str | None:
    if raw_kind is None:
        return None
    token = str(raw_kind or "").strip().lower()
    if not token:
        return None
    return VARIABLE_KIND_ALIASES.get(token)


def is_system_variable(name: Any) -> bool:
    token = str(name or "").strip()
    if not token:
        return True
    if token in SYSTEM_VARIABLE_EXACT:
        return True
    return any(token.startswith(prefix) for prefix in SYSTEM_VARIABLE_PREFIXES)


def _looks_like_accession_token(token: Any, known_accessions: set[str]) -> bool:
    raw = str(token or "").strip()
    if not raw:
        return False
    if raw.startswith("@") and len(raw) > 1:
        return True
    canonical = canonicalize_accession(raw)
    if canonical in known_accessions:
        return True
    if not _ACCESSIONISH_TOKEN_RE.fullmatch(raw):
        return False
    return bool(re.search(r"\d", raw) or "_" in raw)


def infer_variable_kind(name: Any, value: Any, known_accessions: set[str] | None = None) -> str:
    if is_system_variable(name):
        return VARIABLE_KIND_ENV
    known_accessions = known_accessions or set()
    if isinstance(value, list) and value and all(str(token).strip().isdigit() for token in value):
        return VARIABLE_KIND_BUSCO_RUNS
    if isinstance(value, list) and value and all(_looks_like_accession_token(token, known_accessions) for token in value):
        return VARIABLE_KIND_ASSEMBLIES
    return VARIABLE_KIND_ENV


def infer_variable_kinds(values: Mapping[str, Any], known_accessions: set[str] | None = None) -> dict[str, str]:
    return {str(name): infer_variable_kind(name, value, known_accessions) for name, value in values.items()}
