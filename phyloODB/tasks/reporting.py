from __future__ import annotations

import os
import re
from datetime import datetime
from pathlib import Path
from typing import Optional

from .task import Task


def sanitize_report_label(value: object | None) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    return text.strip("._-")


def default_reports_root(task: Task) -> Path:
    base = task.db_manager.storage.get_root_base("reports")
    if not base:
        raise ValueError("No reports storage root is configured. Configure a reports root or set REPORTS_DIR.")
    return Path(base)


def resolve_report_run_dir(
    task: Task,
    *,
    namespace: str,
    explicit_dir: Optional[os.PathLike | str] = None,
    explicit_root: Optional[os.PathLike | str] = None,
    run_label: object | None = None,
    cache_attr: Optional[str] = None,
) -> Path:
    if explicit_dir:
        return Path(explicit_dir)
    cache_key = cache_attr or f"_report_dir_{sanitize_report_label(namespace)}"
    cached = getattr(task, cache_key, None)
    if cached is not None:
        return Path(cached)
    root = Path(explicit_root) if explicit_root else default_reports_root(task) / str(namespace)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = sanitize_report_label(run_label)
    folder = f"task_{task.task_id}_{stamp}"
    if suffix:
        folder = f"{folder}_{suffix}"
    out = root / folder
    setattr(task, cache_key, str(out))
    return out


def resolve_report_base_path(
    task: Task,
    *,
    namespace: str,
    default_stem: str,
    explicit_path: Optional[os.PathLike | str] = None,
    explicit_dir: Optional[os.PathLike | str] = None,
    explicit_root: Optional[os.PathLike | str] = None,
    run_label: object | None = None,
    cache_attr: Optional[str] = None,
) -> Path:
    if explicit_path:
        return Path(explicit_path)
    return resolve_report_run_dir(
        task,
        namespace=namespace,
        explicit_dir=explicit_dir,
        explicit_root=explicit_root,
        run_label=run_label,
        cache_attr=cache_attr,
    ) / default_stem
