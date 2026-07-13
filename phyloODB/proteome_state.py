from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

from .proteome_profile_utils import is_staged_busco_input_path


@dataclass(frozen=True)
class ProteomeState:
    genome_path: str
    active_faa: Optional[str]
    archive_faa: Optional[str]

    @property
    def has_active_proteome(self) -> bool:
        return bool(self.active_faa)

    @property
    def has_archive_proteome(self) -> bool:
        return bool(self.archive_faa)

    @property
    def protein_flag(self) -> int:
        return 1 if self.has_active_proteome else 0

    @property
    def isoforms_cleaned_flag(self) -> int:
        return 1 if self.has_active_proteome and self.has_archive_proteome else 0


def summarize_proteome_state(genome_path: str) -> ProteomeState:
    active_faa = None
    archive_faa = None
    if not genome_path or not os.path.isdir(genome_path):
        return ProteomeState(str(genome_path or ""), None, None)

    for fname in sorted(os.listdir(genome_path)):
        low = fname.lower()
        path = os.path.join(genome_path, fname)
        if (
            active_faa is None
            and low.endswith((".faa", ".faa.gz"))
            and ".archive" not in low
            and not is_staged_busco_input_path(fname)
        ):
            active_faa = path
        elif archive_faa is None and low.endswith((".faa.archive", ".faa.archive.gz")):
            archive_faa = path
    return ProteomeState(str(genome_path), active_faa, archive_faa)
