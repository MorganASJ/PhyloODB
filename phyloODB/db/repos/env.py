from __future__ import annotations

import json
from .base import BaseRepository
from ...variable_kinds import infer_variable_kind, normalize_variable_kind


class EnvRepository(BaseRepository):
    def get(self, key):
        if not key:
            return None
        row = self.core.fetchone(
            "SELECT var_value FROM Environment_Variables WHERE var_name = ?",
            (str(key),),
        )
        if not row:
            return None
        return json.loads(row[0])

    def get_many(self, keys=None):
        if keys:
            keys = list(keys)
            placeholders = ", ".join("?" for _ in keys)
            rows = self.core.fetchall(
                f"SELECT var_name, var_value FROM Environment_Variables WHERE var_name IN ({placeholders})",
                tuple(keys),
            )
        else:
            rows = self.core.fetchall("SELECT var_name, var_value FROM Environment_Variables")
        return {var_name: json.loads(var_value) for var_name, var_value in rows}

    def get_records(self, keys=None):
        has_kind = self.manager._column_exists("Environment_Variables", "var_kind")
        kind_expr = "var_kind" if has_kind else "'env'"
        if keys:
            keys = list(keys)
            placeholders = ", ".join("?" for _ in keys)
            rows = self.core.fetchall(
                f"SELECT var_name, var_value, {kind_expr} FROM Environment_Variables WHERE var_name IN ({placeholders})",
                tuple(keys),
            )
        else:
            rows = self.core.fetchall(f"SELECT var_name, var_value, {kind_expr} FROM Environment_Variables")
        records = {}
        for var_name, var_value, var_kind in rows:
            value = json.loads(var_value)
            records[var_name] = {
                "value": value,
                "kind": normalize_variable_kind(var_kind) or infer_variable_kind(var_name, value),
            }
        return records

    def set(self, key, value, *, kind=None):
        return self.set_many({key: value}, kind=kind)

    def set_many(self, values, *, kind=None, kinds=None):
        explicit_kind = normalize_variable_kind(kind)
        if kind is not None and explicit_kind is None:
            raise ValueError(f"Unknown variable kind '{kind}'.")
        kinds = kinds or {}
        with self.core.transaction(operation="set environment variables"):
            has_updated_at = self.manager._column_exists("Environment_Variables", "updated_at")
            has_kind = self.manager._column_exists("Environment_Variables", "var_kind")
            if has_kind and has_updated_at:
                sql = (
                    """
                    INSERT INTO Environment_Variables (var_name, var_value, var_kind)
                    VALUES (?, ?, ?)
                    ON CONFLICT(var_name) DO UPDATE SET
                        var_value = excluded.var_value,
                        var_kind = excluded.var_kind,
                        updated_at = datetime('now')
                    """
                )
            elif has_kind:
                sql = (
                    """
                    INSERT INTO Environment_Variables (var_name, var_value, var_kind)
                    VALUES (?, ?, ?)
                    ON CONFLICT(var_name) DO UPDATE SET
                        var_value = excluded.var_value,
                        var_kind = excluded.var_kind
                    """
                )
            elif has_updated_at:
                sql = (
                    """
                    INSERT INTO Environment_Variables (var_name, var_value)
                    VALUES (?, ?)
                    ON CONFLICT(var_name) DO UPDATE SET var_value = excluded.var_value, updated_at = datetime('now')
                    """
                )
            else:
                sql = (
                    """
                    INSERT INTO Environment_Variables (var_name, var_value)
                    VALUES (?, ?)
                    ON CONFLICT(var_name) DO UPDATE SET var_value = excluded.var_value
                    """
                )
            for var_name, var_value in values.items():
                if has_kind:
                    per_key_kind = normalize_variable_kind(kinds.get(var_name)) if var_name in kinds else None
                    if var_name in kinds and per_key_kind is None:
                        raise ValueError(f"Unknown variable kind '{kinds.get(var_name)}'.")
                    var_kind = per_key_kind or explicit_kind or infer_variable_kind(var_name, var_value)
                    self.core.execute(sql, (str(var_name), json.dumps(var_value), var_kind))
                else:
                    self.core.execute(sql, (str(var_name), json.dumps(var_value)))
            if not bool(getattr(self.manager, "_syncing_storage_env", False)):
                self.manager.storage.ensure_default_roots_from_env()
            return True

    def all(self):
        return self.get_many()

    def all_records(self):
        return self.get_records()
