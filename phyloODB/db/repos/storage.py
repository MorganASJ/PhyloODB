from __future__ import annotations

import json
import os
from typing import Any, Optional, Tuple

from .base import BaseRepository, transactional_methods


@transactional_methods(
    "ensure_default_roots_from_env",
    "normalize_root_states",
    "ensure_root",
    "ensure_cache_root",
    "rename_root",
    "activate_root",
    "deactivate_root",
    "delete_root",
    "rebind_root",
    "bind_genome_location",
    "bind_library_location",
    "move_binding",
    "backfill_table_locations",
)
class StorageRepository(BaseRepository):
    STRICT_ACTIVE_KINDS = {"genomes", "libraries", "orthofinder", "exports", "logs"}
    DEFAULT_ENV_ROOTS = {
        "genomes": "GENOME_DIR",
        "libraries": "LIBRARIES_DIR",
        "orthofinder": "ORTHOFINDER_OUTPUT_DIR",
        "exports": "EXPORTS_DIR",
        "cache": "CACHE_DIR",
        "reports": "REPORTS_DIR",
        "logs": "LOG_DIR",
        "misc": "MISC_DIR",
    }

    def _default_cache_base_path(self) -> str:
        db_path = os.path.abspath(str(self.manager.get_path()) or "phyloodb")
        return os.path.join(os.path.dirname(db_path), "cache")

    def ensure_default_roots_from_env(self) -> None:
        if not getattr(self.manager, "cursor", None):
            return
        env = self.manager.get_environment_variables([*self.DEFAULT_ENV_ROOTS.values(), "LOG_FILE"]) or {}
        for kind, env_key in self.DEFAULT_ENV_ROOTS.items():
            path = env.get(env_key)
            derived_from_log_file = False
            if not path and kind == "logs":
                log_file = env.get("LOG_FILE")
                if log_file:
                    path = os.path.dirname(os.path.abspath(str(log_file)))
                    derived_from_log_file = True
            if path:
                abs_path = os.path.abspath(str(path))
                existing = self.get_root_by_base_path(kind=kind, base_path=abs_path)
                roots = self.list_roots(kind=kind)
                if kind in self.STRICT_ACTIVE_KINDS:
                    if not roots:
                        self.ensure_root(kind=kind, base_path=abs_path, label=env_key, writable=1, is_active=1)
                    elif existing is None:
                        self.ensure_root(kind=kind, base_path=abs_path, label=env_key, writable=0, is_active=0)
                else:
                    self.ensure_root(kind=kind, base_path=abs_path, label=env_key, writable=1, is_active=1)
                if derived_from_log_file:
                    if bool(getattr(self.manager, "_syncing_storage_env", False)):
                        continue
                    self.manager._syncing_storage_env = True
                    try:
                        self.manager.set_environment_variable(env_key, abs_path)
                    finally:
                        self.manager._syncing_storage_env = False
        self.normalize_root_states(promote_if_none=False)

    def _strict_kind(self, kind: str) -> bool:
        return str(kind) in self.STRICT_ACTIVE_KINDS

    @staticmethod
    def _normalize_label(label: Optional[str]) -> Optional[str]:
        if label is None:
            return None
        value = str(label).strip()
        return value or None

    def _assert_label_unique(self, label: Optional[str], *, exclude_root_id: Optional[int] = None) -> None:
        normalized = self._normalize_label(label)
        if normalized is None:
            return
        row = self.core.fetchone(
            "SELECT storage_root_id FROM StorageRoots WHERE label = ? ORDER BY storage_root_id ASC LIMIT 1",
            (normalized,),
        )
        if not row:
            return
        existing_root_id = int(row[0])
        if exclude_root_id is not None and existing_root_id == int(exclude_root_id):
            return
        raise ValueError(f"Storage root label '{normalized}' is already in use by root {existing_root_id}. Labels must be unique.")

    def _assert_base_path_non_overlapping(self, base_path: str, *, exclude_root_id: Optional[int] = None) -> None:
        candidate = os.path.abspath(str(base_path))
        for row in self.list_roots():
            root_id = int(row[0])
            if exclude_root_id is not None and root_id == int(exclude_root_id):
                continue
            existing = os.path.abspath(str(row[3]))
            try:
                common = os.path.commonpath([candidate, existing])
            except ValueError:
                continue
            if common != candidate and common != existing:
                continue
            raise ValueError(
                f"Storage root base path '{candidate}' overlaps existing root {root_id} "
                f"({row[1]}) at '{existing}'. Root base paths must be unique and non-overlapping."
            )

    def _sync_strict_siblings(self, kind: str, active_root_id: int) -> None:
        self.core.execute(
            """
            UPDATE StorageRoots
            SET is_active = CASE WHEN storage_root_id = ? THEN 1 ELSE 0 END,
                writable = CASE WHEN storage_root_id = ? THEN 1 ELSE 0 END,
                updated_at = datetime('now')
            WHERE logical_kind = ?
            """,
            (int(active_root_id), int(active_root_id), str(kind)),
        )

    def _sync_env_for_kind(self, kind: str) -> None:
        env_key = self.DEFAULT_ENV_ROOTS.get(str(kind))
        if not env_key:
            return
        root_id = self.get_default_root_id(str(kind))
        value = None
        if root_id is not None:
            row = self.get_root(int(root_id))
            if row and row[3]:
                value = os.path.abspath(str(row[3]))
        if bool(getattr(self.manager, "_syncing_storage_env", False)):
            return
        self.manager._syncing_storage_env = True
        try:
            self.manager.set_environment_variable(env_key, value)
        finally:
            self.manager._syncing_storage_env = False

    def normalize_root_states(
        self,
        *,
        kind: Optional[str] = None,
        sync_env: bool = False,
        promote_if_none: bool = True,
    ) -> None:
        kinds = [str(kind)] if kind else sorted(self.STRICT_ACTIVE_KINDS)
        changed = False
        for current_kind in kinds:
            if current_kind not in self.STRICT_ACTIVE_KINDS:
                continue
            rows = self.list_roots(kind=current_kind) or []
            if not rows:
                continue
            active_rows = [row for row in rows if bool(row[5])]
            chosen_root_id: Optional[int] = None
            if len(active_rows) == 1:
                chosen_root_id = int(active_rows[0][0])
            elif len(active_rows) > 1:
                chosen_root_id = min(int(row[0]) for row in active_rows)
            elif promote_if_none:
                chosen_root_id = min(int(row[0]) for row in rows)
            for row in rows:
                root_id = int(row[0])
                should_active = 1 if chosen_root_id is not None and root_id == chosen_root_id else 0
                should_writable = should_active
                if int(bool(row[4])) != should_writable or int(bool(row[5])) != should_active:
                    self.core.execute(
                        """
                        UPDATE StorageRoots
                        SET writable = ?, is_active = ?, updated_at = datetime('now')
                        WHERE storage_root_id = ?
                        """,
                        (should_writable, should_active, root_id),
                    )
                    changed = True
            if sync_env:
                self._sync_env_for_kind(current_kind)
        if changed:
            self.conn.commit()

    def ensure_root(
        self,
        *,
        kind: str,
        base_path: str,
        label: Optional[str] = None,
        writable: int = 1,
        is_active: int = 1,
        metadata: Optional[dict] = None,
    ) -> Optional[int]:
        if not base_path:
            return None
        kind = str(kind)
        label = self._normalize_label(label)
        base_path = os.path.abspath(str(base_path))
        row = self.core.fetchone(
            "SELECT storage_root_id FROM StorageRoots WHERE logical_kind = ? AND base_path = ?",
            (kind, base_path),
        )
        existing_kind_roots = self.list_roots(kind=kind) or []
        if self._strict_kind(kind):
            requested_active = 1 if bool(is_active) else 0
            if row is None and existing_kind_roots:
                requested_active = 0
            is_active = requested_active
            writable = is_active
        meta_json = self.manager._json_dump(metadata or {})
        if row:
            root_id = int(row[0])
            self._assert_label_unique(label, exclude_root_id=root_id)
            self._assert_base_path_non_overlapping(base_path, exclude_root_id=root_id)
            self.core.execute(
                """
                UPDATE StorageRoots
                SET label = COALESCE(?, label), writable = ?, is_active = ?, metadata_json = ?, updated_at = datetime('now')
                WHERE storage_root_id = ?
                """,
                (label, int(bool(writable)), int(bool(is_active)), meta_json, root_id),
            )
            if self._strict_kind(kind) and bool(is_active):
                self._sync_strict_siblings(kind, root_id)
            self.conn.commit()
            return root_id
        self._assert_label_unique(label)
        self._assert_base_path_non_overlapping(base_path)
        self.core.execute(
            """
            INSERT INTO StorageRoots (logical_kind, label, base_path, writable, is_active, metadata_json)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (kind, label, base_path, int(bool(writable)), int(bool(is_active)), meta_json),
        )
        root_id = int(self.cursor.lastrowid)
        if self._strict_kind(kind) and bool(is_active):
            self._sync_strict_siblings(kind, root_id)
        self.conn.commit()
        return root_id

    def get_root(self, root_id: int):
        return self.core.fetchone("SELECT * FROM StorageRoots WHERE storage_root_id = ?", (int(root_id),))

    def get_default_root_id(self, kind: str) -> Optional[int]:
        row = self.core.fetchone(
            "SELECT storage_root_id FROM StorageRoots WHERE logical_kind = ? AND is_active = 1 ORDER BY storage_root_id ASC LIMIT 1",
            (str(kind),),
        )
        return int(row[0]) if row and row[0] is not None else None

    def get_root_base(
        self,
        kind: str,
        *,
        fallback: Optional[str] = None,
        ensure_from_env: bool = True,
    ) -> Optional[str]:
        kind = str(kind)
        root_id = self.get_default_root_id(kind)
        if root_id is not None:
            row = self.get_root(int(root_id))
            if row and row[3]:
                return os.path.abspath(str(row[3]))
        env_key = self.DEFAULT_ENV_ROOTS.get(str(kind))
        roots = self.list_roots(kind=kind) or []
        if self._strict_kind(kind) and roots:
            if fallback:
                return os.path.abspath(str(fallback))
            return None
        if env_key:
            env_path = self.manager.env.get(env_key)
            if env_path:
                abs_path = os.path.abspath(str(env_path))
                if ensure_from_env:
                    self.ensure_root(kind=kind, base_path=abs_path, label=env_key, writable=1, is_active=1)
                return abs_path
        if fallback:
            return os.path.abspath(str(fallback))
        return None

    def require_root_base(self, kind: str, *, message: Optional[str] = None) -> str:
        base = self.get_root_base(kind, ensure_from_env=True)
        if base:
            return str(base)
        env_key = self.DEFAULT_ENV_ROOTS.get(str(kind))
        detail = f" Configure a {kind} storage root"
        if env_key:
            detail += f" or set {env_key}."
        else:
            detail += "."
        raise ValueError(message or f"No active {kind} storage root is configured.{detail}")

    def ensure_cache_root(self) -> str:
        base = self.get_root_base("cache", ensure_from_env=True)
        if base:
            return str(base)
        cache_base = self._default_cache_base_path()
        root_id = self.ensure_root(kind="cache", base_path=cache_base, label="CACHE_DIR", writable=1, is_active=1)
        if root_id is not None:
            self.manager.set_environment_variable("CACHE_DIR", cache_base)
        return cache_base

    def get_path_for_kind(
        self,
        kind: str,
        *parts: str,
        fallback: Optional[str] = None,
        ensure_from_env: bool = True,
    ) -> Optional[str]:
        base = self.get_root_base(str(kind), fallback=fallback, ensure_from_env=ensure_from_env)
        if not base:
            return None
        if not parts:
            return base
        return os.path.abspath(os.path.join(base, *[str(part) for part in parts]))

    def list_roots(self, kind: Optional[str] = None):
        if kind:
            return self.core.fetchall("SELECT * FROM StorageRoots WHERE logical_kind = ? ORDER BY storage_root_id ASC", (str(kind),))
        return self.core.fetchall("SELECT * FROM StorageRoots ORDER BY logical_kind ASC, storage_root_id ASC")

    def get_root_by_base_path(self, *, kind: str, base_path: str):
        return self.core.fetchone(
            "SELECT * FROM StorageRoots WHERE logical_kind = ? AND base_path = ? ORDER BY storage_root_id ASC LIMIT 1",
            (str(kind), os.path.abspath(str(base_path))),
        )

    def resolve_root_token(self, token: Any, *, kind: Optional[str] = None):
        raw = str(token or "").strip()
        if not raw:
            raise ValueError("Storage root identifier is required.")
        if raw.isdigit():
            row = self.get_root(int(raw))
            if not row:
                raise ValueError(f"Unknown storage root id {int(raw)}.")
            if kind and str(row[1]) != str(kind):
                raise ValueError(f"Storage root {int(raw)} is kind '{row[1]}', expected '{kind}'.")
            return row
        rows = self.list_roots(kind=kind) or []
        matches = [row for row in rows if str(row[2] or "") == raw]
        if not matches:
            scope = f" for kind '{kind}'" if kind else ""
            raise ValueError(f"Unknown storage root label '{raw}'{scope}.")
        if len(matches) > 1:
            ids = ", ".join(str(int(row[0])) for row in matches)
            raise ValueError(f"Storage root label '{raw}' is ambiguous; matches root ids: {ids}.")
        return matches[0]

    def rename_root(self, root_id: int, label: str) -> str:
        row = self.get_root(int(root_id))
        if not row:
            raise ValueError(f"Unknown storage root id {root_id}.")
        normalized = self._normalize_label(label)
        if normalized is None:
            raise ValueError("A non-empty storage root label is required.")
        self._assert_label_unique(normalized, exclude_root_id=int(root_id))
        self.core.execute(
            "UPDATE StorageRoots SET label = ?, updated_at = datetime('now') WHERE storage_root_id = ?",
            (normalized, int(root_id)),
        )
        self.conn.commit()
        return normalized

    def activate_root(self, root_id: int) -> int:
        row = self.get_root(int(root_id))
        if not row:
            raise ValueError(f"Unknown storage root id {root_id}.")
        kind = str(row[1])
        if self._strict_kind(kind):
            self._sync_strict_siblings(kind, int(root_id))
            self.conn.commit()
            self._sync_env_for_kind(kind)
            return int(root_id)
        self.core.execute(
            "UPDATE StorageRoots SET is_active = 1, writable = 1, updated_at = datetime('now') WHERE storage_root_id = ?",
            (int(root_id),),
        )
        self.conn.commit()
        self._sync_env_for_kind(kind)
        return int(root_id)

    def deactivate_root(self, root_id: int) -> int:
        row = self.get_root(int(root_id))
        if not row:
            raise ValueError(f"Unknown storage root id {root_id}.")
        kind = str(row[1])
        writable = 0 if self._strict_kind(kind) else int(bool(row[4]))
        self.core.execute(
            "UPDATE StorageRoots SET is_active = 0, writable = ?, updated_at = datetime('now') WHERE storage_root_id = ?",
            (writable, int(root_id)),
        )
        self.conn.commit()
        self._sync_env_for_kind(kind)
        return int(root_id)

    def delete_root(self, root_id: int) -> bool:
        row = self.get_root(int(root_id))
        kind = str(row[1]) if row else None
        self.core.execute("DELETE FROM StorageRoots WHERE storage_root_id = ?", (int(root_id),))
        deleted = int(self.cursor.rowcount or 0)
        self.conn.commit()
        if deleted and kind:
            self._sync_env_for_kind(kind)
        return bool(deleted)

    def resolve_path(
        self,
        *,
        storage_root_id: Optional[int] = None,
        relative_path: Optional[str] = None,
        fallback_path: Optional[str] = None,
    ) -> Optional[str]:
        if storage_root_id is not None and relative_path is not None:
            row = self.get_root(int(storage_root_id))
            if row and row[3]:
                return os.path.abspath(os.path.join(str(row[3]), str(relative_path)))
        if fallback_path:
            return str(fallback_path)
        return None

    def rebind_root(self, root_id: int, new_base_path: str) -> bool:
        row = self.get_root(int(root_id))
        self._assert_base_path_non_overlapping(new_base_path, exclude_root_id=int(root_id))
        self.core.execute(
            "UPDATE StorageRoots SET base_path = ?, updated_at = datetime('now') WHERE storage_root_id = ?",
            (os.path.abspath(str(new_base_path)), int(root_id)),
        )
        self.conn.commit()
        if row:
            self._sync_env_for_kind(str(row[1]))
        return True

    def detect_root_for_path(self, path: Optional[str], *, kind: Optional[str] = None) -> Tuple[Optional[int], Optional[str]]:
        if not path:
            return (None, None)
        abs_path = os.path.abspath(str(path))
        rows = self.list_roots(kind=kind)
        rows = sorted(rows, key=lambda r: len(str(r[3] or "")), reverse=True)
        for row in rows:
            root_id = int(row[0])
            base_path = os.path.abspath(str(row[3]))
            try:
                common = os.path.commonpath([abs_path, base_path])
            except ValueError:
                continue
            if common != base_path:
                continue
            rel = os.path.relpath(abs_path, base_path)
            if rel == ".":
                rel = ""
            return (root_id, rel)
        return (None, None)

    def _update_binding(self, table: str, key_col: str, key_val: Any, *, root_id: Optional[int], rel_path: Optional[str], location: Optional[str]) -> bool:
        self.core.execute(
            f"UPDATE {table} SET storage_root_id = ?, relative_path = ?, location = ? WHERE {key_col} = ?",
            (root_id, rel_path, location, key_val),
        )
        self.conn.commit()
        return True

    def bind_genome_location(self, accession: str, path: str, *, kind: str = "genomes") -> bool:
        if not path:
            return False
        root_id, rel = self.detect_root_for_path(path, kind=kind)
        if root_id is None:
            default_root = self.get_default_root_id(kind)
            if default_root is None:
                raise ValueError(f"No active {kind} storage root is configured.")
            row = self.get_root(int(default_root)) if default_root is not None else None
            if not row:
                raise ValueError(f"Storage root {int(default_root)} is missing.")
            base_path = os.path.abspath(str(row[3]))
            abs_path = os.path.abspath(str(path))
            try:
                common = os.path.commonpath([abs_path, base_path])
            except ValueError as exc:
                raise ValueError(f"Path '{abs_path}' is not under the configured {kind} root '{base_path}'.") from exc
            if common != base_path:
                raise ValueError(f"Path '{abs_path}' is not under the configured {kind} root '{base_path}'.")
            rel = os.path.relpath(abs_path, base_path)
            root_id = default_root
        return self._update_binding("Genome", "accession", accession, root_id=root_id, rel_path=rel, location=os.path.abspath(str(path)))

    def bind_library_location(self, library_id: int, path: str, *, kind: str = "libraries") -> bool:
        if not path:
            return False
        root_id, rel = self.detect_root_for_path(path, kind=kind)
        if root_id is None:
            default_root = self.get_default_root_id(kind)
            if default_root is None:
                raise ValueError(f"No active {kind} storage root is configured.")
            row = self.get_root(int(default_root)) if default_root is not None else None
            if not row:
                raise ValueError(f"Storage root {int(default_root)} is missing.")
            base_path = os.path.abspath(str(row[3]))
            abs_path = os.path.abspath(str(path))
            try:
                common = os.path.commonpath([abs_path, base_path])
            except ValueError as exc:
                raise ValueError(f"Path '{abs_path}' is not under the configured {kind} root '{base_path}'.") from exc
            if common != base_path:
                raise ValueError(f"Path '{abs_path}' is not under the configured {kind} root '{base_path}'.")
            rel = os.path.relpath(abs_path, base_path)
            root_id = default_root
        return self._update_binding("Libraries", "library_id", int(library_id), root_id=root_id, rel_path=rel, location=os.path.abspath(str(path)))

    def resolve_genome_location(self, accession: str) -> Optional[str]:
        row = self.core.fetchone(
            "SELECT storage_root_id, relative_path, location FROM Genome WHERE accession = ?",
            (accession,),
        )
        if not row:
            return None
        return self.resolve_path(storage_root_id=row[0], relative_path=row[1], fallback_path=row[2])

    def resolve_library_location(self, library_id: int) -> Optional[str]:
        row = self.core.fetchone(
            "SELECT storage_root_id, relative_path, location FROM Libraries WHERE library_id = ?",
            (int(library_id),),
        )
        if not row:
            return None
        return self.resolve_path(storage_root_id=row[0], relative_path=row[1], fallback_path=row[2])

    def move_binding(
        self,
        *,
        table: str,
        key_col: str,
        key_val: Any,
        new_root_id: int,
        new_relative_path: str,
        update_location: bool = True,
    ) -> bool:
        location = None
        if update_location:
            location = self.resolve_path(storage_root_id=int(new_root_id), relative_path=str(new_relative_path))
        self.core.execute(
            f"UPDATE {table} SET storage_root_id = ?, relative_path = ?, location = COALESCE(?, location) WHERE {key_col} = ?",
            (int(new_root_id), str(new_relative_path), location, key_val),
        )
        self.conn.commit()
        return True

    def backfill_table_locations(self, *, table: str, key_col: str, kind: str) -> int:
        rows = self.core.fetchall(
            f"SELECT {key_col}, location, storage_root_id, relative_path FROM {table} WHERE location IS NOT NULL"
        )
        updated = 0
        for key_val, location, root_id, rel_path in rows:
            if root_id is not None and rel_path:
                continue
            detected_root, detected_rel = self.detect_root_for_path(location, kind=kind)
            if detected_root is None:
                default_root = self.get_default_root_id(kind)
                if default_root is None:
                    continue
                detected_root = default_root
                row = self.get_root(int(default_root))
                if not row:
                    continue
                detected_rel = os.path.relpath(os.path.abspath(str(location)), os.path.abspath(str(row[3])))
            self.core.execute(
                f"UPDATE {table} SET storage_root_id = ?, relative_path = ? WHERE {key_col} = ?",
                (int(detected_root), str(detected_rel), key_val),
            )
            updated += 1
        if updated:
            self.conn.commit()
        return updated

    def list_genome_bindings(self, *, root_id: Optional[int] = None):
        sql = "SELECT accession, storage_root_id, relative_path, location FROM Genome"
        params: list[Any] = []
        if root_id is not None:
            sql += " WHERE storage_root_id = ?"
            params.append(int(root_id))
        sql += " ORDER BY accession"
        return self.core.fetchall(sql, tuple(params))

    def list_library_bindings(self, *, root_id: Optional[int] = None):
        sql = "SELECT library_id, library_name, storage_root_id, relative_path, location FROM Libraries"
        params: list[Any] = []
        if root_id is not None:
            sql += " WHERE storage_root_id = ?"
            params.append(int(root_id))
        sql += " ORDER BY library_id"
        return self.core.fetchall(sql, tuple(params))

    def create_filesystem_operation(
        self,
        *,
        operation_type: str,
        source_path: Optional[str] = None,
        staging_path: Optional[str] = None,
        destination_path: Optional[str] = None,
        status: str = "preparing",
        payload: Optional[dict[str, Any]] = None,
    ) -> int:
        with self.core.transaction(operation=f"create filesystem operation {operation_type}"):
            self.core.execute(
                """
                INSERT INTO FilesystemOperations (
                    operation_type, source_path, staging_path, destination_path,
                    status, payload_json
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    str(operation_type),
                    source_path,
                    staging_path,
                    destination_path,
                    str(status),
                    json.dumps(payload or {}, sort_keys=True),
                ),
            )
            return int(self.cursor.lastrowid)

    def update_filesystem_operation(
        self,
        operation_id: int,
        *,
        status: str,
        error_message: Optional[str] = None,
        staging_path: Optional[str] = None,
    ) -> None:
        with self.core.transaction(operation=f"update filesystem operation {operation_id}"):
            self.core.execute(
                """
                UPDATE FilesystemOperations
                SET status = ?, error_message = ?,
                    staging_path = COALESCE(?, staging_path),
                    updated_at = datetime('now')
                WHERE operation_id = ?
                """,
                (str(status), error_message, staging_path, int(operation_id)),
            )

    def list_filesystem_operations(self, *, pending_only: bool = False):
        sql = """
            SELECT operation_id, operation_type, source_path, staging_path,
                   destination_path, status, payload_json, error_message,
                   created_at, updated_at
            FROM FilesystemOperations
        """
        params: tuple[Any, ...] = ()
        if pending_only:
            sql += " WHERE status IN ('preparing', 'prepared', 'db_committed', 'failed')"
        sql += " ORDER BY operation_id"
        return self.core.fetchall(sql, params)

    def get_filesystem_operation(self, operation_id: int):
        return self.core.fetchone(
            """
            SELECT operation_id, operation_type, source_path, staging_path,
                   destination_path, status, payload_json, error_message,
                   created_at, updated_at
            FROM FilesystemOperations WHERE operation_id = ?
            """,
            (int(operation_id),),
        )

    def count_artifacts_for_root(self, root_id: int) -> int:
        row = self.core.fetchone(
            "SELECT COUNT(*) FROM Artifacts WHERE storage_root_id = ?",
            (int(root_id),),
        )
        return int(row[0]) if row and row[0] is not None else 0

    def count_bound_entities(self, root_id: int) -> dict[str, int]:
        return {
            "genomes": len(self.list_genome_bindings(root_id=int(root_id)) or []),
            "libraries": len(self.list_library_bindings(root_id=int(root_id)) or []),
            "artifact_rows": self.count_artifacts_for_root(int(root_id)),
        }
