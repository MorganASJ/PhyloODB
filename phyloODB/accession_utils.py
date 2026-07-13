"""Helpers for canonical accession handling without DB-layer imports."""
from __future__ import annotations

import re
from typing import Any, Sequence


_NCBI_ACCESSION_RE = re.compile(r"^(GC[AF])_?(\d+)\.(\d+)$", re.IGNORECASE)


def canonicalize_accession(value: Any) -> str:
    """Canonicalize common accession spellings without altering custom ids."""

    if value is None:
        return ""
    token = str(value).strip()
    if not token:
        return ""
    match = _NCBI_ACCESSION_RE.fullmatch(token)
    if not match:
        return token
    prefix, digits, version = match.groups()
    return f"{prefix.upper()}_{digits}.{version}"


def canonicalize_accessions(accessions: Sequence[Any]) -> list[str]:
    """Canonicalize a sequence of accessions while preserving order."""

    cleaned: list[str] = []
    for accession in accessions or []:
        if accession is None:
            continue
        if isinstance(accession, (list, tuple)):
            if not accession:
                continue
            token = accession[0]
        else:
            token = accession
        canonical = canonicalize_accession(token)
        if canonical:
            cleaned.append(canonical)
    return cleaned
