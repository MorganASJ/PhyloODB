from __future__ import annotations

import os
import sqlite3
from typing import Any, Optional

from .base import BaseRepository, transactional


class ArtifactRepository(BaseRepository):
    NATURAL_ROOT_KINDS = {
        "busco_run": "genomes",
        "orthofinder_run": "orthofinder",
        "export_run": "exports",
        "library": "libraries",
        "genome": "genomes",
        "proteome_profile": "genomes",
    }

    @transactional("register artifact")
    def register(
        self,
        *,
        owner_type: str,
        owner_id: Any,
        artifact_type: str,
        path: Optional[str] = None,
        storage_root_id: Optional[int] = None,
        relative_path: Optional[str] = None,
        role: Optional[str] = None,
        status: str = "ready",
        is_dir: bool = False,
        format: Optional[str] = None,
        sequence_kind: Optional[str] = None,
        checksum: Optional[str] = None,
        size_bytes: Optional[int] = None,
        metadata: Optional[dict] = None,
    ) -> Optional[int]:
        root_id = storage_root_id
        rel_path = relative_path
        abs_path = os.path.abspath(str(path)) if path else None
        if abs_path:
            root_id, rel_path = self._resolve_binding_for_path(
                owner_type=str(owner_type),
                owner_id=owner_id,
                abs_path=abs_path,
                storage_root_id=root_id,
                relative_path=rel_path,
            )
        if not abs_path and root_id is not None and rel_path is not None:
            abs_path = self.manager.storage.resolve_path(storage_root_id=int(root_id), relative_path=str(rel_path))
        role_value = role or ""
        sequence_kind_value = sequence_kind or ""
        rel_path_value = rel_path if rel_path is not None else ""
        metadata_json = self.manager._json_dump(metadata or {})

        exact = self._find_artifact_by_unique_key(
            owner_type=str(owner_type),
            owner_id=owner_id,
            artifact_type=str(artifact_type),
            role=role_value,
            sequence_kind=sequence_kind_value,
            relative_path=rel_path_value,
        )
        existing = exact or self._find_existing_artifact(
            owner_type=str(owner_type),
            owner_id=owner_id,
            artifact_type=str(artifact_type),
            role=role_value,
            sequence_kind=sequence_kind_value,
        )
        target_row = existing
        try:
            if target_row:
                self._update_artifact_row(
                    artifact_id=int(target_row[0]),
                    status=status,
                    storage_root_id=root_id,
                    relative_path=rel_path_value,
                    absolute_path=abs_path,
                    is_dir=is_dir,
                    format=format,
                    sequence_kind=sequence_kind_value,
                    checksum=checksum,
                    size_bytes=size_bytes,
                    metadata_json=metadata_json,
                )
            else:
                self._insert_artifact_row(
                    owner_type=str(owner_type),
                    owner_id=owner_id,
                    artifact_type=str(artifact_type),
                    role=role_value,
                    status=status,
                    storage_root_id=root_id,
                    relative_path=rel_path_value,
                    absolute_path=abs_path,
                    is_dir=is_dir,
                    format=format,
                    sequence_kind=sequence_kind_value,
                    checksum=checksum,
                    size_bytes=size_bytes,
                    metadata_json=metadata_json,
                )
        except sqlite3.IntegrityError:
            conflict = self._find_artifact_by_unique_key(
                owner_type=str(owner_type),
                owner_id=owner_id,
                artifact_type=str(artifact_type),
                role=role_value,
                sequence_kind=sequence_kind_value,
                relative_path=rel_path_value,
            )
            if conflict is None:
                raise
            if target_row and int(target_row[0]) != int(conflict[0]):
                self.core.execute("DELETE FROM Artifacts WHERE artifact_id = ?", (int(target_row[0]),))
            self._update_artifact_row(
                artifact_id=int(conflict[0]),
                status=status,
                storage_root_id=root_id,
                relative_path=rel_path_value,
                absolute_path=abs_path,
                is_dir=is_dir,
                format=format,
                sequence_kind=sequence_kind_value,
                checksum=checksum,
                size_bytes=size_bytes,
                metadata_json=metadata_json,
            )
        self.conn.commit()
        row = self._find_artifact_by_unique_key(
            owner_type=str(owner_type),
            owner_id=owner_id,
            artifact_type=str(artifact_type),
            role=role_value,
            sequence_kind=sequence_kind_value,
            relative_path=rel_path_value,
        ) or self._find_existing_artifact(
            owner_type=str(owner_type),
            owner_id=owner_id,
            artifact_type=str(artifact_type),
            role=role_value,
            sequence_kind=sequence_kind_value,
        )
        return int(row[0]) if row and row[0] is not None else None

    def _find_existing_artifact(
        self,
        *,
        owner_type: str,
        owner_id: Any,
        artifact_type: str,
        role: str,
        sequence_kind: str,
    ):
        return self.core.fetchone(
            """
            SELECT * FROM Artifacts
            WHERE owner_type = ? AND owner_id = ? AND artifact_type = ?
              AND COALESCE(role, '') = COALESCE(?, '')
              AND COALESCE(sequence_kind, '') = COALESCE(?, '')
            ORDER BY artifact_id DESC LIMIT 1
            """,
            (str(owner_type), str(owner_id), str(artifact_type), role, sequence_kind),
        )

    def _find_artifact_by_unique_key(
        self,
        *,
        owner_type: str,
        owner_id: Any,
        artifact_type: str,
        role: str,
        sequence_kind: str,
        relative_path: str,
    ):
        return self.core.fetchone(
            """
            SELECT * FROM Artifacts
            WHERE owner_type = ? AND owner_id = ? AND artifact_type = ?
              AND COALESCE(role, '') = COALESCE(?, '')
              AND COALESCE(sequence_kind, '') = COALESCE(?, '')
              AND COALESCE(relative_path, '') = COALESCE(?, '')
            ORDER BY artifact_id DESC LIMIT 1
            """,
            (str(owner_type), str(owner_id), str(artifact_type), role, sequence_kind, relative_path),
        )

    def _update_artifact_row(
        self,
        *,
        artifact_id: int,
        status: str,
        storage_root_id: Optional[int],
        relative_path: str,
        absolute_path: Optional[str],
        is_dir: bool,
        format: Optional[str],
        sequence_kind: str,
        checksum: Optional[str],
        size_bytes: Optional[int],
        metadata_json: str,
    ) -> None:
        self.core.execute(
            """
            UPDATE Artifacts
            SET status = ?,
                storage_root_id = ?,
                relative_path = ?,
                absolute_path = ?,
                is_dir = ?,
                format = ?,
                sequence_kind = ?,
                checksum = COALESCE(?, checksum),
                size_bytes = COALESCE(?, size_bytes),
                metadata_json = ?,
                updated_at = datetime('now')
            WHERE artifact_id = ?
            """,
            (
                status,
                int(storage_root_id) if storage_root_id is not None else None,
                relative_path,
                absolute_path,
                1 if is_dir else 0,
                format,
                sequence_kind,
                checksum,
                size_bytes,
                metadata_json,
                int(artifact_id),
            ),
        )

    def _insert_artifact_row(
        self,
        *,
        owner_type: str,
        owner_id: Any,
        artifact_type: str,
        role: str,
        status: str,
        storage_root_id: Optional[int],
        relative_path: str,
        absolute_path: Optional[str],
        is_dir: bool,
        format: Optional[str],
        sequence_kind: str,
        checksum: Optional[str],
        size_bytes: Optional[int],
        metadata_json: str,
    ) -> None:
        self.core.execute(
            """
            INSERT INTO Artifacts (
                owner_type, owner_id, artifact_type, role, status,
                storage_root_id, relative_path, absolute_path, is_dir,
                format, sequence_kind, checksum, size_bytes, metadata_json
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                owner_type,
                str(owner_id),
                artifact_type,
                role,
                status,
                int(storage_root_id) if storage_root_id is not None else None,
                relative_path,
                absolute_path,
                1 if is_dir else 0,
                format,
                sequence_kind,
                checksum,
                size_bytes,
                metadata_json,
            ),
        )

    def _relative_to_root(self, root_id: int, abs_path: str) -> str:
        root = self.manager.storage.get_root(int(root_id))
        if not root or not root[3]:
            raise ValueError(f"Unknown storage root {int(root_id)}.")
        base_path = os.path.abspath(str(root[3]))
        try:
            common = os.path.commonpath([abs_path, base_path])
        except ValueError as exc:
            raise ValueError(f"Path '{abs_path}' is not under storage root {int(root_id)} at '{base_path}'.") from exc
        if common != base_path:
            raise ValueError(f"Path '{abs_path}' is not under storage root {int(root_id)} at '{base_path}'.")
        rel = os.path.relpath(abs_path, base_path)
        return "" if rel == "." else rel

    def _owner_bound_root(self, owner_type: str, owner_id: Any) -> tuple[Optional[int], Optional[str]]:
        owner_type = str(owner_type)
        if owner_type == "genome":
            row = self.core.fetchone(
                "SELECT storage_root_id, location FROM Genome WHERE accession = ?",
                (str(owner_id),),
            )
            return (int(row[0]), str(row[1])) if row and row[0] is not None else (None, str(row[1]) if row and row[1] else None)
        if owner_type == "library":
            row = self.core.fetchone(
                "SELECT storage_root_id, location FROM Libraries WHERE library_id = ?",
                (int(owner_id),),
            )
            return (int(row[0]), str(row[1])) if row and row[0] is not None else (None, str(row[1]) if row and row[1] else None)
        if owner_type == "busco_run":
            row = self.core.fetchone(
                """
                SELECT g.storage_root_id, g.location
                FROM BUSCO_Runs br
                LEFT JOIN Genome g ON g.accession = br.accession
                WHERE br.run_id = ?
                """,
                (int(owner_id),),
            )
            return (int(row[0]), str(row[1])) if row and row[0] is not None else (None, str(row[1]) if row and row[1] else None)
        if owner_type == "proteome_profile":
            row = self.core.fetchone(
                """
                SELECT g.storage_root_id, g.location
                FROM Proteome_Profiles pp
                LEFT JOIN Genome g ON g.accession = pp.accession
                WHERE pp.proteome_profile_id = ?
                """,
                (int(owner_id),),
            )
            return (int(row[0]), str(row[1])) if row and row[0] is not None else (None, str(row[1]) if row and row[1] else None)
        return (None, None)

    def _fallback_root_id_for_owner(self, owner_type: str) -> Optional[int]:
        owner_type = str(owner_type)
        if owner_type == "blast_db":
            root_id = self.manager.storage.get_default_root_id("cache")
            if root_id is not None:
                return int(root_id)
            self.manager.storage.ensure_cache_root()
            root_id = self.manager.storage.get_default_root_id("cache")
            return int(root_id) if root_id is not None else None
        kind = self.NATURAL_ROOT_KINDS.get(owner_type)
        if kind is None:
            return self.manager.storage.get_default_root_id("misc")
        return self.manager.storage.get_default_root_id(kind)

    def _resolve_binding_for_path(
        self,
        *,
        owner_type: str,
        owner_id: Any,
        abs_path: str,
        storage_root_id: Optional[int],
        relative_path: Optional[str],
    ) -> tuple[Optional[int], Optional[str]]:
        if storage_root_id is not None:
            rel = relative_path if relative_path is not None else self._relative_to_root(int(storage_root_id), abs_path)
            return (int(storage_root_id), rel)

        bound_root_id, _bound_location = self._owner_bound_root(owner_type, owner_id)
        if bound_root_id is not None:
            return (int(bound_root_id), self._relative_to_root(int(bound_root_id), abs_path))

        detected_root_id, detected_rel = self.manager.storage.detect_root_for_path(abs_path, kind=None)
        if detected_root_id is not None:
            return (int(detected_root_id), detected_rel)

        fallback_root_id = self._fallback_root_id_for_owner(owner_type)
        if fallback_root_id is None:
            natural_kind = self.NATURAL_ROOT_KINDS.get(str(owner_type))
            if natural_kind:
                raise ValueError(
                    f"No active {natural_kind} storage root is configured for {owner_type} artifacts."
                )
            raise ValueError("No misc storage root is configured for generic artifacts.")
        return (int(fallback_root_id), self._relative_to_root(int(fallback_root_id), abs_path))

    def find(
        self,
        *,
        owner_type: Optional[str] = None,
        owner_id: Optional[Any] = None,
        owner_ids: Optional[list[Any]] = None,
        artifact_type: Optional[str] = None,
        role: Optional[str] = None,
        sequence_kind: Optional[str] = None,
        status: Optional[str] = None,
    ):
        clauses = []
        params = []
        if owner_type is not None:
            clauses.append("owner_type = ?")
            params.append(str(owner_type))
        if owner_id is not None:
            clauses.append("owner_id = ?")
            params.append(str(owner_id))
        elif owner_ids:
            owner_vals = [str(value) for value in owner_ids]
            placeholders = ",".join("?" for _ in owner_vals)
            clauses.append(f"owner_id IN ({placeholders})")
            params.extend(owner_vals)
        if artifact_type is not None:
            clauses.append("artifact_type = ?")
            params.append(str(artifact_type))
        if role is not None:
            clauses.append("COALESCE(role, '') = COALESCE(?, '')")
            params.append(role)
        if sequence_kind is not None:
            clauses.append("COALESCE(sequence_kind, '') = COALESCE(?, '')")
            params.append(sequence_kind)
        if status is not None:
            clauses.append("status = ?")
            params.append(str(status))
        sql = "SELECT * FROM Artifacts"
        if clauses:
            sql += " WHERE " + " AND ".join(clauses)
        sql += " ORDER BY artifact_id ASC"
        return self.core.fetchall(sql, tuple(params))

    def get(self, artifact_id: int):
        return self.core.fetchone("SELECT * FROM Artifacts WHERE artifact_id = ?", (int(artifact_id),))

    @transactional("set artifact status")
    def set_status(
        self,
        artifact_id: int,
        status: str,
        *,
        size_bytes: Optional[int] = None,
        checksum: Optional[str] = None,
    ) -> bool:
        self.core.execute(
            """
            UPDATE Artifacts
            SET status = ?,
                size_bytes = COALESCE(?, size_bytes),
                checksum = COALESCE(?, checksum),
                updated_at = datetime('now')
            WHERE artifact_id = ?
            """,
            (str(status), size_bytes, checksum, int(artifact_id)),
        )
        self.conn.commit()
        return True

    @transactional("move artifact binding")
    def move_binding(
        self,
        artifact_id: int,
        *,
        new_root_id: int,
        new_relative_path: str,
        update_absolute_path: bool = True,
    ) -> bool:
        absolute_path = None
        if update_absolute_path:
            absolute_path = self.manager.storage.resolve_path(
                storage_root_id=int(new_root_id),
                relative_path=str(new_relative_path),
            )
        self.core.execute(
            """
            UPDATE Artifacts
            SET storage_root_id = ?, relative_path = ?, absolute_path = COALESCE(?, absolute_path), updated_at = datetime('now')
            WHERE artifact_id = ?
            """,
            (int(new_root_id), str(new_relative_path), absolute_path, int(artifact_id)),
        )
        self.conn.commit()
        return True

    def resolve_path(self, artifact_row_or_id) -> Optional[str]:
        row = artifact_row_or_id
        if isinstance(artifact_row_or_id, int):
            row = self.get(artifact_row_or_id)
        if not row:
            return None
        storage_root_id = row[6] if len(row) > 6 else None
        relative_path = row[7] if len(row) > 7 else None
        absolute_path = row[8] if len(row) > 8 else None
        return self.manager.storage.resolve_path(
            storage_root_id=storage_root_id,
            relative_path=relative_path,
            fallback_path=absolute_path,
        )
