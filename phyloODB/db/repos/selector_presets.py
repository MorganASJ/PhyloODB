from __future__ import annotations

import json
import re
from typing import Any, Mapping, Optional

from .base import BaseRepository, transactional


_PRESET_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def _display_timestamp(value: Any) -> Optional[str]:
    if value is None:
        return None
    return str(value)


def normalize_preset_name(name: Any) -> str:
    token = str(name or "").strip()
    if not token:
        raise ValueError("Selector preset name cannot be empty.")
    if token.startswith("@"):
        raise ValueError("Selector preset names do not use '@'; reserve @NAME for stored accession panels.")
    if not _PRESET_NAME_RE.match(token):
        raise ValueError("Selector preset names may contain only letters, numbers, dots, underscores, and hyphens.")
    return token


class SelectorPresetRepository(BaseRepository):
    def get(self, name: Any) -> Optional[dict[str, Any]]:
        preset_name = normalize_preset_name(name)
        row = self.core.fetchone(
            """
            SELECT preset_name, selector_json, description, created_at, updated_at
            FROM Selector_Presets
            WHERE preset_name = ?
            """,
            (preset_name,),
        )
        if not row:
            return None
        selector = json.loads(row[1] or "{}")
        if not isinstance(selector, dict):
            selector = {}
        return {
            "preset_name": row[0],
            "selector": selector,
            "description": row[2],
            "created_at": _display_timestamp(row[3]),
            "updated_at": _display_timestamp(row[4]),
        }

    def list(self) -> list[dict[str, Any]]:
        rows = self.core.fetchall(
            """
            SELECT preset_name, selector_json, description, created_at, updated_at
            FROM Selector_Presets
            ORDER BY LOWER(preset_name), preset_name
            """
        )
        presets: list[dict[str, Any]] = []
        for row in rows:
            selector = json.loads(row[1] or "{}")
            if not isinstance(selector, dict):
                selector = {}
            presets.append(
                {
                    "preset_name": row[0],
                    "selector": selector,
                    "description": row[2],
                    "created_at": _display_timestamp(row[3]),
                    "updated_at": _display_timestamp(row[4]),
                }
            )
        return presets

    @transactional("save selector preset")
    def save(self, name: Any, selector: Mapping[str, Any], *, description: Optional[str] = None) -> str:
        preset_name = normalize_preset_name(name)
        payload = json.dumps(dict(selector), sort_keys=True)
        self.core.execute(
            """
            INSERT INTO Selector_Presets (preset_name, selector_json, description)
            VALUES (?, ?, ?)
            ON CONFLICT(preset_name) DO UPDATE SET
                selector_json = excluded.selector_json,
                description = excluded.description,
                updated_at = datetime('now')
            """,
            (preset_name, payload, description),
        )
        return preset_name

    @transactional("delete selector preset")
    def delete(self, name: Any) -> bool:
        preset_name = normalize_preset_name(name)
        self.core.execute("DELETE FROM Selector_Presets WHERE preset_name = ?", (preset_name,))
        return self.core.cursor.rowcount > 0
