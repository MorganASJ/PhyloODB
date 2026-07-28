from __future__ import annotations

import hashlib
import json
import os
from typing import Any, Iterable, Optional, Sequence

from .base import BaseRepository, transactional
from ...proteome_profile_utils import DEFAULT_CLEAN_PROFILE, RAW_PROFILE, is_staged_busco_input_path
from ...proteome_state import summarize_proteome_state


def _sha256_file(path: str) -> Optional[str]:
    try:
        hasher = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                if not chunk:
                    break
                hasher.update(chunk)
        return hasher.hexdigest()
    except OSError:
        return None


def _count_fasta_headers(path: str) -> Optional[int]:
    import gzip

    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        count = 0
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    count += 1
        return count
    except OSError:
        return None


class ProteomeRepository(BaseRepository):
    def get_default_cleaned_profile(self, accession: str):
        accession = str(accession)
        row = self.core.fetchone(
            """
            SELECT proteome_profile_id, accession, profile_name, kind, parent_profile_id,
                   artifact_id, status, sequence_count, checksum, is_default,
                   created_at, updated_at, metadata_json
            FROM Proteome_Profiles
            WHERE accession = ?
              AND kind != 'raw'
              AND COALESCE(status, 'ready') = 'ready'
            ORDER BY is_default DESC,
                     CASE WHEN profile_name = ? THEN 0 ELSE 1 END,
                     proteome_profile_id DESC
            LIMIT 1
            """,
            (accession, DEFAULT_CLEAN_PROFILE),
        )
        return row

    def get_default_cleaned_profile_name(self, accession: str) -> Optional[str]:
        row = self.get_default_cleaned_profile(str(accession))
        return str(row[2]) if row and row[2] is not None else None

    def resolve_selector_profile_name(self, accession: str, profile_name: Optional[str]) -> Optional[str]:
        requested = str(profile_name or "").strip() or None
        if requested is None:
            return None
        if requested == DEFAULT_CLEAN_PROFILE:
            actual = self.get_default_cleaned_profile_name(str(accession))
            if actual:
                return actual
        return requested

    def profile_matches_selector(self, accession: str, actual_profile_name: Optional[str], requested_profile_name: Optional[str]) -> bool:
        requested = str(requested_profile_name or "").strip() or None
        actual = str(actual_profile_name or "").strip() or None
        if requested is None:
            return True
        if actual is None:
            return False
        resolved = self.resolve_selector_profile_name(str(accession), requested)
        return bool(resolved and actual == resolved)

    @staticmethod
    def _canonical_path(path: Optional[str]) -> Optional[str]:
        if not path:
            return None
        try:
            return os.path.realpath(os.path.abspath(str(path)))
        except OSError:
            return os.path.abspath(str(path))

    @staticmethod
    def _profile_name_from_filename(path: str) -> Optional[str]:
        base = os.path.basename(str(path))
        low = base.lower()
        for suffix in (".faa.gz", ".faa"):
            if low.endswith(suffix):
                token = base[: -len(suffix)].strip()
                return token or None
        return None

    def get(self, profile_id: int):
        return self.core.fetchone(
            """
            SELECT proteome_profile_id, accession, profile_name, kind, parent_profile_id,
                   artifact_id, status, sequence_count, checksum, is_default,
                   created_at, updated_at, metadata_json
            FROM Proteome_Profiles
            WHERE proteome_profile_id = ?
            """,
            (int(profile_id),),
        )

    def get_profile(self, accession: str, profile_name: str):
        return self.core.fetchone(
            """
            SELECT proteome_profile_id, accession, profile_name, kind, parent_profile_id,
                   artifact_id, status, sequence_count, checksum, is_default,
                   created_at, updated_at, metadata_json
            FROM Proteome_Profiles
            WHERE accession = ? AND profile_name = ?
            ORDER BY proteome_profile_id DESC
            LIMIT 1
            """,
            (str(accession), str(profile_name)),
        )

    def get_default_profile(self, accession: str):
        row = self.core.fetchone(
            """
            SELECT proteome_profile_id, accession, profile_name, kind, parent_profile_id,
                   artifact_id, status, sequence_count, checksum, is_default,
                   created_at, updated_at, metadata_json
            FROM Proteome_Profiles
            WHERE accession = ? AND is_default = 1
            ORDER BY proteome_profile_id DESC
            LIMIT 1
            """,
            (str(accession),),
        )
        if row:
            return row
        return self.get_default_cleaned_profile(str(accession)) or self.get_profile(str(accession), RAW_PROFILE)

    def get_default_profile_name(self, accession: str) -> Optional[str]:
        row = self.get_default_profile(str(accession))
        return str(row[2]) if row and row[2] is not None else None

    def list_profiles(
        self,
        *,
        accessions: Optional[Sequence[str]] = None,
        profile_names: Optional[Sequence[str]] = None,
    ):
        clauses = ["1=1"]
        params: list[Any] = []
        if accessions:
            vals = [str(item) for item in accessions if item is not None]
            if vals:
                placeholders = ",".join("?" for _ in vals)
                clauses.append(f"p.accession IN ({placeholders})")
                params.extend(vals)
        if profile_names:
            vals = [str(item) for item in profile_names if item is not None]
            if vals:
                placeholders = ",".join("?" for _ in vals)
                clauses.append(f"p.profile_name IN ({placeholders})")
                params.extend(vals)
        return self.core.fetchall(
            f"""
            SELECT p.proteome_profile_id, p.accession, p.profile_name, p.kind, p.parent_profile_id,
                   p.artifact_id, p.status, p.sequence_count, p.checksum, p.is_default,
                   p.created_at, p.updated_at, p.metadata_json
            FROM Proteome_Profiles p
            WHERE {' AND '.join(clauses)}
            ORDER BY p.accession, p.is_default DESC, p.profile_name, p.proteome_profile_id
            """,
            tuple(params),
        ) or []

    def resolve_path(self, profile_row_or_id) -> Optional[str]:
        row = profile_row_or_id
        if isinstance(profile_row_or_id, int):
            row = self.get(profile_row_or_id)
        if not row:
            return None
        artifact_id = row[5] if len(row) > 5 else None
        if artifact_id is None:
            return None
        return self.manager.artifacts.resolve_path(int(artifact_id))

    @transactional("set proteome profile status")
    def set_status(self, profile_id: int, status: str) -> bool:
        row = self.get(int(profile_id))
        if not row:
            return False
        self.core.execute(
            "UPDATE Proteome_Profiles SET status = ?, updated_at = datetime('now') WHERE proteome_profile_id = ?",
            (str(status), int(profile_id)),
        )
        artifact_id = row[5] if len(row) > 5 else None
        if artifact_id is not None:
            self.manager.artifacts.set_status(int(artifact_id), str(status))
        self.conn.commit()
        return True

    def find_profile_by_path(self, accession: str, path: Optional[str]):
        canonical_target = self._canonical_path(path)
        if not canonical_target:
            return None
        for row in self.list_profiles(accessions=[str(accession)]):
            resolved = self.resolve_path(row)
            if self._canonical_path(resolved) == canonical_target:
                return row
        return None

    @transactional("set default proteome profile")
    def set_default_profile(self, accession: str, *, profile_name: Optional[str] = None, profile_id: Optional[int] = None) -> bool:
        if profile_name is None and profile_id is None:
            raise ValueError("Provide either profile_name or profile_id to set the default proteome profile.")
        target = None
        if profile_id is not None:
            target = self.get(int(profile_id))
            if not target:
                raise ValueError(f"Unknown proteome profile id {int(profile_id)}.")
            if str(target[1]) != str(accession):
                raise ValueError("Proteome profile id does not belong to the requested accession.")
        else:
            target = self.get_profile(str(accession), str(profile_name))
            if not target:
                raise ValueError(f"Unknown proteome profile '{profile_name}' for accession '{accession}'.")
        self.core.execute("UPDATE Proteome_Profiles SET is_default = 0, updated_at = datetime('now') WHERE accession = ?", (str(accession),))
        self.core.execute(
            "UPDATE Proteome_Profiles SET is_default = 1, updated_at = datetime('now') WHERE proteome_profile_id = ?",
            (int(target[0]),),
        )
        self.conn.commit()
        return True

    def _ensure_profile_row(
        self,
        *,
        accession: str,
        profile_name: str,
        kind: str,
        parent_profile_id: Optional[int],
        status: str,
        sequence_count: Optional[int],
        checksum: Optional[str],
        metadata: Optional[dict],
        is_default: bool,
    ) -> int:
        existing = self.get_profile(str(accession), str(profile_name))
        metadata_json = json.dumps(metadata or {})
        if existing:
            self.core.execute(
                """
                UPDATE Proteome_Profiles
                SET kind = ?,
                    parent_profile_id = ?,
                    status = ?,
                    sequence_count = COALESCE(?, sequence_count),
                    checksum = COALESCE(?, checksum),
                    metadata_json = ?,
                    is_default = CASE WHEN ? THEN 1 ELSE is_default END,
                    updated_at = datetime('now')
                WHERE proteome_profile_id = ?
                """,
                (
                    str(kind),
                    int(parent_profile_id) if parent_profile_id is not None else None,
                    str(status),
                    sequence_count,
                    checksum,
                    metadata_json,
                    1 if is_default else 0,
                    int(existing[0]),
                ),
            )
            self.conn.commit()
            return int(existing[0])
        self.core.execute(
            """
            INSERT INTO Proteome_Profiles (
                accession, profile_name, kind, parent_profile_id, status,
                sequence_count, checksum, is_default, metadata_json, created_at, updated_at
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'), datetime('now'))
            """,
            (
                str(accession),
                str(profile_name),
                str(kind),
                int(parent_profile_id) if parent_profile_id is not None else None,
                str(status),
                sequence_count,
                checksum,
                1 if is_default else 0,
                metadata_json,
            ),
        )
        self.conn.commit()
        return int(self.cursor.lastrowid)

    @transactional("register proteome profile")
    def register_profile(
        self,
        *,
        accession: str,
        profile_name: str,
        path: str,
        kind: str,
        parent_profile_id: Optional[int] = None,
        is_default: bool = False,
        status: str = "ready",
        metadata: Optional[dict] = None,
        format: str = "fasta",
    ) -> int:
        abs_path = os.path.abspath(str(path))
        sequence_count = _count_fasta_headers(abs_path)
        checksum = _sha256_file(abs_path)
        profile_id = self._ensure_profile_row(
            accession=str(accession),
            profile_name=str(profile_name),
            kind=str(kind),
            parent_profile_id=parent_profile_id,
            status=status,
            sequence_count=sequence_count,
            checksum=checksum,
            metadata=metadata,
            is_default=is_default,
        )
        artifact_id = self.manager.artifacts.register(
            owner_type="proteome_profile",
            owner_id=profile_id,
            artifact_type="proteome_faa",
            path=abs_path,
            format=format,
            sequence_kind="prot",
            metadata={"accession": accession, "profile_name": profile_name, **(metadata or {})},
        )
        self.core.execute(
            """
            UPDATE Proteome_Profiles
            SET artifact_id = ?, sequence_count = COALESCE(?, sequence_count), checksum = COALESCE(?, checksum),
                updated_at = datetime('now')
            WHERE proteome_profile_id = ?
            """,
            (
                int(artifact_id) if artifact_id is not None else None,
                sequence_count,
                checksum,
                int(profile_id),
            ),
        )
        if is_default:
            self.core.execute("UPDATE Proteome_Profiles SET is_default = 0 WHERE accession = ? AND proteome_profile_id != ?", (str(accession), int(profile_id)))
        self.conn.commit()
        return int(profile_id)

    def ensure_raw_profile(self, accession: str, *, path: Optional[str] = None, is_default: bool = False) -> Optional[int]:
        existing = self.get_profile(str(accession), RAW_PROFILE)
        if existing and existing[5] is not None:
            existing_path = self.resolve_path(existing)
            if existing_path and os.path.isfile(existing_path):
                if is_default:
                    self.set_default_profile(str(accession), profile_id=int(existing[0]))
                return int(existing[0])
        genome_path = self.manager.genomes.resolve_path(str(accession))
        candidate = str(path or "").strip() or None
        if candidate is None and genome_path and os.path.isdir(genome_path):
            for fname in sorted(os.listdir(genome_path)):
                low = fname.lower()
                if (
                    low.endswith((".faa", ".faa.gz"))
                    and ".archive" not in low
                    and not is_staged_busco_input_path(fname)
                ):
                    candidate = os.path.join(genome_path, fname)
                    break
        if not candidate or not os.path.exists(candidate):
            return None
        return self.register_profile(
            accession=str(accession),
            profile_name=RAW_PROFILE,
            path=str(candidate),
            kind="raw",
            parent_profile_id=None,
            is_default=is_default or self.get_default_profile(str(accession)) is None,
            metadata={"source": "ensure_raw_profile"},
        )

    def sync_profiles_from_filesystem(
        self,
        accession: str,
        genome_path: str,
        *,
        set_default: bool = True,
    ) -> dict[str, Any]:
        accession = str(accession)
        genome_path = str(genome_path or "")
        existing_rows = self.list_profiles(accessions=[accession])
        existing_by_name = {str(row[2]): row for row in existing_rows if row and row[2] is not None}
        state = summarize_proteome_state(genome_path)
        seen_names: set[str] = set()
        ready_names: set[str] = set()
        profile_ids: dict[str, int] = {}

        raw_path = state.active_faa if state.active_faa and os.path.exists(state.active_faa) else None
        raw_id: Optional[int] = None
        if raw_path:
            raw_id = self.register_profile(
                accession=accession,
                profile_name=RAW_PROFILE,
                path=str(raw_path),
                kind="raw",
                parent_profile_id=None,
                is_default=False,
                metadata={"source": "filesystem_sync", "kind": "raw"},
            )
            profile_ids[RAW_PROFILE] = int(raw_id)
            seen_names.add(RAW_PROFILE)
            ready_names.add(RAW_PROFILE)

        profiles_dir = os.path.join(genome_path, "proteome_profiles")
        if os.path.isdir(profiles_dir):
            for fname in sorted(os.listdir(profiles_dir)):
                path = os.path.join(profiles_dir, fname)
                if not os.path.isfile(path):
                    continue
                profile_name = self._profile_name_from_filename(path)
                if not profile_name:
                    continue
                parent_profile_id = None if profile_name == RAW_PROFILE else raw_id
                kind = "raw" if profile_name == RAW_PROFILE else "derived"
                profile_id = self.register_profile(
                    accession=accession,
                    profile_name=profile_name,
                    path=path,
                    kind=kind,
                    parent_profile_id=parent_profile_id,
                    is_default=False,
                    metadata={"source": "filesystem_sync", "kind": kind},
                )
                profile_ids[str(profile_name)] = int(profile_id)
                seen_names.add(str(profile_name))
                ready_names.add(str(profile_name))

        stale_names = set(existing_by_name) - seen_names
        for profile_name in sorted(stale_names):
            row = existing_by_name.get(profile_name)
            if row:
                self.set_status(int(row[0]), "stale")

        default_profile_name: Optional[str] = None
        if set_default:
            current_default = self.get_default_profile_name(accession)
            nonraw_ready = sorted(name for name in ready_names if name != RAW_PROFILE)
            if current_default in ready_names:
                default_profile_name = current_default
            elif DEFAULT_CLEAN_PROFILE in ready_names:
                default_profile_name = DEFAULT_CLEAN_PROFILE
            elif nonraw_ready:
                default_profile_name = nonraw_ready[0]
            elif RAW_PROFILE in ready_names:
                default_profile_name = RAW_PROFILE
            elif ready_names:
                default_profile_name = sorted(ready_names)[0]
            if default_profile_name:
                self.set_default_profile(accession, profile_name=default_profile_name)

        for profile_name in ready_names:
            row = self.get_profile(accession, profile_name)
            if row and str(row[6] or "") != "ready":
                self.set_status(int(row[0]), "ready")

        if default_profile_name is None:
            default_row = self.get_default_profile(accession)
            default_profile_name = str(default_row[2]) if default_row and default_row[2] is not None else None

        return {
            "profiles": dict(profile_ids),
            "ready_profiles": sorted(ready_names),
            "stale_profiles": sorted(stale_names),
            "default_profile": default_profile_name,
            "has_protein": bool(ready_names),
            "raw_profile_id": int(raw_id) if raw_id is not None else None,
        }

    @transactional("record proteome preparation")
    def record_preparation(
        self,
        *,
        accession: str,
        input_profile_id: int,
        output_profile_id: int,
        preparation_type: str,
        used_gff: bool = False,
        gff_artifact_id: Optional[int] = None,
        skip_gff: bool = False,
        skip_cdhit: bool = False,
        gff_priority: bool = False,
        cdhit_identity: Optional[float] = None,
        cdhit_threads: Optional[int] = None,
        input_count: Optional[int] = None,
        output_count: Optional[int] = None,
        gff_removed: Optional[int] = None,
        cdhit_removed: Optional[int] = None,
        total_removed: Optional[int] = None,
        status: str = "completed",
        params: Optional[dict] = None,
    ) -> int:
        self.core.execute(
            """
            INSERT INTO Proteome_Preparations (
                accession, input_profile_id, output_profile_id, preparation_type,
                used_gff, gff_artifact_id, skip_gff, skip_cdhit, gff_priority,
                cdhit_identity, cdhit_threads, input_count, output_count,
                gff_removed, cdhit_removed, total_removed, status, params_json, created_at
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))
            """,
            (
                str(accession),
                int(input_profile_id),
                int(output_profile_id),
                str(preparation_type),
                1 if used_gff else 0,
                int(gff_artifact_id) if gff_artifact_id is not None else None,
                1 if skip_gff else 0,
                1 if skip_cdhit else 0,
                1 if gff_priority else 0,
                float(cdhit_identity) if cdhit_identity is not None else None,
                int(cdhit_threads) if cdhit_threads is not None else None,
                int(input_count) if input_count is not None else None,
                int(output_count) if output_count is not None else None,
                int(gff_removed) if gff_removed is not None else None,
                int(cdhit_removed) if cdhit_removed is not None else None,
                int(total_removed) if total_removed is not None else None,
                str(status),
                json.dumps(params or {}),
            ),
        )
        self.conn.commit()
        return int(self.cursor.lastrowid)

    def latest_preparation_for_output(self, output_profile_id: int):
        return self.core.fetchone(
            """
            SELECT preparation_id, accession, input_profile_id, output_profile_id, preparation_type,
                   used_gff, gff_artifact_id, skip_gff, skip_cdhit, gff_priority,
                   cdhit_identity, cdhit_threads, input_count, output_count,
                   gff_removed, cdhit_removed, total_removed, status, params_json, created_at
            FROM Proteome_Preparations
            WHERE output_profile_id = ?
            ORDER BY preparation_id DESC
            LIMIT 1
            """,
            (int(output_profile_id),),
        )
