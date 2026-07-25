"""Purge subcommands for selective data cleanup.

The purge command is dry-run by default. Use --apply to execute changes.
"""
from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from ...database import DBManager
from ...proteome_state import summarize_proteome_state
from ...selector_utils import expand_busco_run_id_variables, resolve_selector_candidates
from ...variable_kinds import is_system_variable
from ..support.argparse_utils import AppendCommaSeparated, _validate_date
from ..support.common import STRICT_WORKING_ROOT_KINDS


CHUNK_SIZE = 800


@dataclass
class PurgeResult:
    subject: str
    dry_run: bool
    counts: List[Tuple[str, int]] = field(default_factory=list)
    files: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    deleted_files: int = 0
    skipped_files: int = 0
    failed_files: int = 0


# ---------------------------------------------------------------------------
# Database and filesystem purge helpers
# Purpose: Gather purge-specific selectors, deletion plans, and constrained
# filesystem cleanup logic for immediate maintenance commands.
# ---------------------------------------------------------------------------


def _chunked(values: Sequence[Any], size: int = CHUNK_SIZE) -> Iterable[List[Any]]:
    if size <= 0:
        raise ValueError("chunk size must be positive")
    for idx in range(0, len(values), size):
        yield list(values[idx:idx + size])


def _table_exists(manager: DBManager, table: str) -> bool:
    manager.cursor.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name = ?",
        (table,),
    )
    return manager.cursor.fetchone() is not None


def _column_exists(manager: DBManager, table: str, column: str) -> bool:
    if not _table_exists(manager, table):
        return False
    manager.cursor.execute(f"PRAGMA table_info({table})")
    return any(str(row[1]) == column for row in (manager.cursor.fetchall() or []))


def _fetch_ids_in(
    manager: DBManager,
    table: str,
    id_column: str,
    filter_column: str,
    values: Sequence[Any],
) -> List[Any]:
    if not values or not _table_exists(manager, table):
        return []
    ids: List[Any] = []
    for chunk in _chunked(values):
        placeholders = ",".join("?" for _ in chunk)
        sql = f"SELECT DISTINCT {id_column} FROM {table} WHERE {filter_column} IN ({placeholders})"
        manager.cursor.execute(sql, tuple(chunk))
        ids.extend(row[0] for row in (manager.cursor.fetchall() or []))
    return list(dict.fromkeys(ids))


def _count_where_in(
    manager: DBManager,
    table: str,
    column: str,
    values: Sequence[Any],
    *,
    extra_where: str = "",
    extra_params: Sequence[Any] = (),
) -> int:
    if not values or not _table_exists(manager, table):
        return 0
    total = 0
    for chunk in _chunked(values):
        placeholders = ",".join("?" for _ in chunk)
        sql = f"SELECT COUNT(*) FROM {table} WHERE {column} IN ({placeholders})"
        params: List[Any] = list(chunk)
        if extra_where:
            sql += f" AND ({extra_where})"
            params.extend(extra_params)
        manager.cursor.execute(sql, tuple(params))
        row = manager.cursor.fetchone()
        total += int(row[0] if row else 0)
    return total


def _count_where(
    manager: DBManager,
    table: str,
    where_sql: str,
    params: Sequence[Any],
) -> int:
    if not _table_exists(manager, table):
        return 0
    sql = f"SELECT COUNT(*) FROM {table}"
    if where_sql.strip():
        sql += f" WHERE {where_sql}"
    manager.cursor.execute(sql, tuple(params))
    row = manager.cursor.fetchone()
    return int(row[0] if row else 0)


def _delete_where_in(
    manager: DBManager,
    table: str,
    column: str,
    values: Sequence[Any],
    *,
    extra_where: str = "",
    extra_params: Sequence[Any] = (),
) -> int:
    if not values or not _table_exists(manager, table):
        return 0
    deleted = 0
    for chunk in _chunked(values):
        placeholders = ",".join("?" for _ in chunk)
        sql = f"DELETE FROM {table} WHERE {column} IN ({placeholders})"
        params: List[Any] = list(chunk)
        if extra_where:
            sql += f" AND ({extra_where})"
            params.extend(extra_params)
        manager.cursor.execute(sql, tuple(params))
        deleted += int(manager.cursor.rowcount or 0)
    return deleted


def _delete_where(manager: DBManager, table: str, where_sql: str, params: Sequence[Any]) -> int:
    if not _table_exists(manager, table):
        return 0
    sql = f"DELETE FROM {table}"
    if where_sql.strip():
        sql += f" WHERE {where_sql}"
    manager.cursor.execute(sql, tuple(params))
    return int(manager.cursor.rowcount or 0)


def _count_filtered_by_accessions(
    manager: DBManager,
    table: str,
    accessions: Sequence[str],
    *,
    extra_where: str = "",
    extra_params: Sequence[Any] = (),
) -> int:
    if accessions:
        return _count_where_in(
            manager,
            table,
            "accession",
            accessions,
            extra_where=extra_where,
            extra_params=extra_params,
        )
    return _count_where(manager, table, extra_where, extra_params)


def _delete_filtered_by_accessions(
    manager: DBManager,
    table: str,
    accessions: Sequence[str],
    *,
    extra_where: str = "",
    extra_params: Sequence[Any] = (),
) -> int:
    if accessions:
        return _delete_where_in(
            manager,
            table,
            "accession",
            accessions,
            extra_where=extra_where,
            extra_params=extra_params,
        )
    return _delete_where(manager, table, extra_where, extra_params)


def _build_full_accession_where(
    accessions: Sequence[str],
    *,
    extra_where: str = "",
    extra_params: Sequence[Any] = (),
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []
    if accessions:
        placeholders = ",".join("?" for _ in accessions)
        clauses.append(f"accession IN ({placeholders})")
        params.extend([str(acc) for acc in accessions])
    if extra_where:
        clauses.append(f"({extra_where})")
        params.extend(extra_params)
    return " AND ".join(clauses), params


def _fetch_distinct_column_where(
    manager: DBManager,
    table: str,
    column: str,
    where_sql: str,
    params: Sequence[Any],
) -> List[Any]:
    if not _table_exists(manager, table):
        return []
    sql = f"SELECT DISTINCT {column} FROM {table}"
    if where_sql.strip():
        sql += f" WHERE {where_sql}"
    manager.cursor.execute(sql, tuple(params))
    return [row[0] for row in (manager.cursor.fetchall() or []) if row and row[0] is not None]


def _fetch_affected_library_accessions(
    manager: DBManager,
    table: str,
    where_sql: str,
    params: Sequence[Any],
    *,
    library_column: str = "target_library_id",
) -> List[Tuple[int, str]]:
    if not _table_exists(manager, table):
        return []
    sql = f"SELECT DISTINCT {library_column}, accession FROM {table}"
    if where_sql.strip():
        sql += f" WHERE {where_sql}"
    manager.cursor.execute(sql, tuple(params))
    out: List[Tuple[int, str]] = []
    for row in manager.cursor.fetchall() or []:
        if not row or row[0] is None or row[1] is None:
            continue
        out.append((int(row[0]), str(row[1])))
    return out


def _filter_empty_runs(
    manager: DBManager,
    run_ids: Sequence[str],
    *,
    row_tables: Sequence[str],
) -> List[str]:
    remaining: List[str] = []
    for run_id in run_ids:
        has_rows = False
        for table in row_tables:
            if not _table_exists(manager, table):
                continue
            manager.cursor.execute(f"SELECT 1 FROM {table} WHERE run_id = ? LIMIT 1", (str(run_id),))
            if manager.cursor.fetchone() is not None:
                has_rows = True
                break
        if not has_rows:
            remaining.append(str(run_id))
    return remaining


def _expand_file_token(token: str) -> List[Path]:
    raw = str(token or "").strip()
    if not raw:
        return []
    p = Path(raw)
    candidates: List[Path] = []
    if p.exists():
        return [p]
    parent = p.parent if str(p.parent) != "" else Path(".")
    pattern = f"{p.name}.*"
    try:
        for child in parent.glob(pattern):
            candidates.append(child)
    except OSError:
        return []
    return candidates


def _path_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except (OSError, RuntimeError, ValueError):
        return False


def _allowed_roots(manager: DBManager) -> List[Path]:
    roots: List[Path] = []
    db_dir = Path(manager.get_path()).resolve().parent
    roots.append(db_dir)
    for row in manager.storage.list_roots():
        base_path = row[3] if row and len(row) > 3 else None
        if isinstance(base_path, str) and base_path.strip():
            try:
                roots.append(Path(base_path).expanduser().resolve())
            except (OSError, RuntimeError, ValueError):
                continue
    env = manager.get_environment_variables(["BUSCO_LINEAGE_DIR"]) or {}
    lineage_dir = env.get("BUSCO_LINEAGE_DIR")
    if isinstance(lineage_dir, str) and lineage_dir.strip():
        try:
            roots.append(Path(lineage_dir).expanduser().resolve())
        except (OSError, RuntimeError, ValueError):
            pass
    unique: List[Path] = []
    seen: Set[str] = set()
    for root in roots:
        key = str(root)
        if key not in seen:
            seen.add(key)
            unique.append(root)
    return unique


def _collect_locations_for_accessions(manager: DBManager, accessions: Sequence[str]) -> Set[str]:
    files: Set[str] = set()
    if not accessions or not _table_exists(manager, "Genome"):
        return files
    for chunk in _chunked(accessions):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"SELECT location FROM Genome WHERE accession IN ({placeholders})",
            tuple(chunk),
        )
        for (loc,) in (manager.cursor.fetchall() or []):
            if isinstance(loc, str) and loc.strip():
                files.add(loc.strip())
    return files


def _collect_busco_run_ids_for_accessions(manager: DBManager, accessions: Sequence[str]) -> List[int]:
    if not accessions or not _table_exists(manager, "BUSCO_Runs"):
        return []
    run_ids: List[int] = []
    for chunk in _chunked(accessions):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"SELECT run_id FROM BUSCO_Runs WHERE accession IN ({placeholders}) ORDER BY run_id",
            tuple(chunk),
        )
        run_ids.extend(int(row[0]) for row in (manager.cursor.fetchall() or []) if row and row[0] is not None)
    return list(dict.fromkeys(run_ids))


def _collect_busco_artifact_ids_for_run_ids(manager: DBManager, run_ids: Sequence[int]) -> List[int]:
    if not run_ids or not _table_exists(manager, "Artifacts"):
        return []
    artifact_ids: List[int] = []
    owner_tokens = [str(run_id) for run_id in run_ids]
    for chunk in _chunked(owner_tokens):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"""
            SELECT artifact_id
            FROM Artifacts
            WHERE owner_type = 'busco_run'
              AND owner_id IN ({placeholders})
            ORDER BY artifact_id
            """,
            tuple(chunk),
        )
        artifact_ids.extend(int(row[0]) for row in (manager.cursor.fetchall() or []) if row and row[0] is not None)
    return list(dict.fromkeys(artifact_ids))


def _collect_genome_artifact_ids(manager: DBManager, accessions: Sequence[str]) -> List[int]:
    if not accessions or not _table_exists(manager, "Artifacts"):
        return []
    artifact_ids: List[int] = []
    owner_tokens = [str(accession) for accession in accessions]
    for chunk in _chunked(owner_tokens):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"""
            SELECT artifact_id
            FROM Artifacts
            WHERE owner_type = 'genome'
              AND owner_id IN ({placeholders})
            ORDER BY artifact_id
            """,
            tuple(chunk),
        )
        artifact_ids.extend(int(row[0]) for row in (manager.cursor.fetchall() or []) if row and row[0] is not None)
    return list(dict.fromkeys(artifact_ids))


def _collect_busco_locations_for_run_ids(manager: DBManager, run_ids: Sequence[int]) -> Set[str]:
    files: Set[str] = set()
    if not run_ids:
        return files
    for chunk in _chunked(run_ids):
        placeholders = ",".join("?" for _ in chunk)
        if _table_exists(manager, "BUSCO_Runs"):
            manager.cursor.execute(
                f"SELECT result_dir FROM BUSCO_Runs WHERE run_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
            manager.cursor.execute(
                f"SELECT accession, library_id, result_dir FROM BUSCO_Runs WHERE run_id IN ({placeholders})",
                tuple(chunk),
            )
            for accession, library_id, loc in (manager.cursor.fetchall() or []):
                if accession is None or library_id is None or not isinstance(loc, str) or not loc.strip():
                    continue
                genome_path = manager.genomes.resolve_path(str(accession))
                library_name = manager.libraries.get_name(int(library_id))
                if not genome_path or not library_name:
                    continue
                legacy_path = Path(str(genome_path)) / f"{library_name}_results"
                try:
                    result_path = Path(str(loc)).expanduser().resolve()
                    if legacy_path.exists() or legacy_path.is_symlink():
                        if legacy_path.resolve() == result_path:
                            files.add(str(legacy_path))
                except (OSError, RuntimeError, ValueError):
                    continue
        if _table_exists(manager, "BUSCO_Run_Family_Locations"):
            manager.cursor.execute(
                f"SELECT location FROM BUSCO_Run_Family_Locations WHERE run_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
        if _table_exists(manager, "BUSCO_Run_Family_Artifacts"):
            manager.cursor.execute(
                f"SELECT location FROM BUSCO_Run_Family_Artifacts WHERE run_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    artifact_ids = _collect_busco_artifact_ids_for_run_ids(manager, run_ids)
    if artifact_ids and _table_exists(manager, "Artifacts"):
        for chunk in _chunked(artifact_ids):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT absolute_path FROM Artifacts WHERE artifact_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    return files


def _collect_genome_artifact_locations(manager: DBManager, accessions: Sequence[str]) -> Set[str]:
    files: Set[str] = set()
    if not accessions or not _table_exists(manager, "Artifacts"):
        return files
    for chunk in _chunked([str(accession) for accession in accessions]):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"""
            SELECT absolute_path
            FROM Artifacts
            WHERE owner_type = 'genome'
              AND owner_id IN ({placeholders})
            """,
            tuple(chunk),
        )
        for (loc,) in (manager.cursor.fetchall() or []):
            if isinstance(loc, str) and loc.strip():
                files.add(loc.strip())
    return files


def _collect_assembly_purge_files(
    manager: DBManager,
    accessions: Sequence[str],
    *,
    run_ids: Sequence[int],
) -> Set[str]:
    files = set(_collect_locations_for_accessions(manager, accessions))
    files.update(_collect_genome_artifact_locations(manager, accessions))
    files.update(_collect_busco_locations_for_run_ids(manager, run_ids))
    if _table_exists(manager, "BUSCO_Family_Locations"):
        for chunk in _chunked(accessions):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT location FROM BUSCO_Family_Locations WHERE accession IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    return files


def _current_genome_locations(manager: DBManager, accessions: Sequence[str]) -> Dict[str, Optional[str]]:
    locations: Dict[str, Optional[str]] = {}
    if not accessions or not _table_exists(manager, "Genome"):
        return locations
    for chunk in _chunked(accessions):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"SELECT accession, location FROM Genome WHERE accession IN ({placeholders})",
            tuple(chunk),
        )
        for accession, location in (manager.cursor.fetchall() or []):
            locations[str(accession)] = str(location) if isinstance(location, str) and location.strip() else None
    return locations


def _recompute_genome_state(manager: DBManager, accession: str, previous_location: Optional[str]) -> None:
    state = summarize_proteome_state(previous_location or "")
    assignments = ["status = 0", "location = NULL", "dl_date = NULL", "protein = ?", "isoforms_cleaned = ?"]
    params: List[Any] = [state.protein_flag, state.isoforms_cleaned_flag]
    if _column_exists(manager, "Genome", "storage_root_id"):
        assignments.append("storage_root_id = NULL")
    if _column_exists(manager, "Genome", "relative_path"):
        assignments.append("relative_path = NULL")
    sql = "UPDATE Genome SET " + ", ".join(assignments) + " WHERE accession = ?"
    params.append(str(accession))
    manager.cursor.execute(sql, tuple(params))


def _verify_no_accession_references(manager: DBManager, accessions: Sequence[str]) -> List[str]:
    if not accessions:
        return []
    manager.cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name")
    table_names = [str(row[0]) for row in (manager.cursor.fetchall() or [])]
    lingering: List[str] = []
    for table in table_names:
        if not _column_exists(manager, table, "accession"):
            continue
        count = _count_where_in(manager, table, "accession", accessions)
        if count > 0:
            lingering.append(f"{table}:{count}")
    return lingering


def _collect_locations_for_libraries(
    manager: DBManager,
    library_ids: Sequence[int],
    *,
    include_library_locations: bool = False,
) -> Set[str]:
    files: Set[str] = set()
    if not library_ids:
        return files
    if _table_exists(manager, "Proteome_BlastDBs"):
        for chunk in _chunked(library_ids):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT location FROM Proteome_BlastDBs WHERE library_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    if _table_exists(manager, "OrthoFinder_Results"):
        for chunk in _chunked(library_ids):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT location FROM OrthoFinder_Results WHERE library_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    if _table_exists(manager, "BUSCO_Family_Locations"):
        for chunk in _chunked(library_ids):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT location FROM BUSCO_Family_Locations WHERE library_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    if include_library_locations and _table_exists(manager, "Libraries"):
        for chunk in _chunked(library_ids):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT location FROM Libraries WHERE library_id IN ({placeholders})",
                tuple(chunk),
            )
            for (loc,) in (manager.cursor.fetchall() or []):
                if isinstance(loc, str) and loc.strip():
                    files.add(loc.strip())
    return files


def _delete_files(
    manager: DBManager,
    file_tokens: Sequence[str],
) -> Tuple[int, int, int, List[str]]:
    if not file_tokens:
        return (0, 0, 0, [])
    roots = _allowed_roots(manager)
    deleted = 0
    skipped = 0
    failed = 0
    failures: List[str] = []
    seen: Set[str] = set()
    for token in file_tokens:
        for target in _expand_file_token(token):
            try:
                target_abs = target.expanduser().absolute()
                resolved = target_abs.resolve()
            except (OSError, RuntimeError, ValueError):
                skipped += 1
                continue
            key = str(target_abs)
            if key in seen:
                continue
            seen.add(key)
            if not any(_path_within(target_abs, root) or _path_within(resolved, root) for root in roots):
                skipped += 1
                continue
            try:
                if target_abs.is_symlink():
                    target_abs.unlink()
                elif target_abs.is_dir():
                    shutil.rmtree(target_abs)
                elif target_abs.exists():
                    target_abs.unlink()
                deleted += 1
            except OSError as exc:
                failed += 1
                failures.append(f"{target_abs}: {exc}")
    return (deleted, skipped, failed, failures)


def _selector_requested(args: argparse.Namespace) -> bool:
    return any(
        (
            bool(getattr(args, "accessions", None)),
            bool(getattr(args, "clade", None)),
            getattr(args, "taxid", None) is not None,
            bool(getattr(args, "downloaded_only", False)),
            bool(getattr(args, "not_downloaded", False)),
            bool(getattr(args, "after", None)),
            bool(getattr(args, "before", None)),
            bool(getattr(args, "level", None)),
            bool(getattr(args, "primary_only", False)),
            bool(getattr(args, "protein_only", False)),
            getattr(args, "status_min", None) is not None,
        )
    )


def _selector_config(args: argparse.Namespace) -> Dict[str, Any]:
    cfg: Dict[str, Any] = {}
    if getattr(args, "accessions", None):
        cfg["accessions"] = list(args.accessions)
    if getattr(args, "clade", None):
        cfg["clade"] = args.clade
    if getattr(args, "taxid", None) is not None:
        cfg["taxid"] = int(args.taxid)
    if getattr(args, "downloaded_only", False):
        cfg["downloaded_only"] = True
    if getattr(args, "not_downloaded", False):
        cfg["not_downloaded"] = True
    if getattr(args, "after", None):
        cfg["after"] = args.after
    if getattr(args, "before", None):
        cfg["before"] = args.before
    if getattr(args, "level", None):
        cfg["level"] = args.level
    if getattr(args, "primary_only", False):
        cfg["primary_only"] = True
    if getattr(args, "protein_only", False):
        cfg["protein_only"] = True
    if getattr(args, "status_min", None) is not None:
        cfg["status_min"] = int(args.status_min)
    return cfg


def _resolve_accessions(manager: DBManager, args: argparse.Namespace) -> List[str]:
    selector_cfg = _selector_config(args)
    if not selector_cfg and not getattr(args, "all", False):
        raise ValueError("Provide selectors or use --all.")
    return resolve_selector_candidates(
        manager,
        selector_cfg,
        allow_all=bool(getattr(args, "all", False)),
        require_candidates=not bool(getattr(args, "all", False)),
        allow_env_defaults=False,
        allow_bare_variables=True,
    )


def _resolve_library_ids(
    manager: DBManager,
    raw_values: Sequence[str],
) -> List[int]:
    resolved: List[int] = []
    for raw in raw_values:
        token = str(raw).strip()
        if not token:
            continue
        if token.isdigit():
            resolved.append(int(token))
            continue
        lib_id = manager.get_library_id(token)
        if lib_id is None:
            raise ValueError(f"Unknown library '{token}'.")
        resolved.append(int(lib_id))
    return list(dict.fromkeys(resolved))


def _resolve_subject_libraries(manager: DBManager, args: argparse.Namespace) -> List[int]:
    raw_tokens = list(getattr(args, "libraries", None) or [])
    if raw_tokens:
        return _resolve_library_ids(manager, raw_tokens)
    if not getattr(args, "all", False):
        raise ValueError("Provide --library/--libraries or use --all.")
    include_core = bool(getattr(args, "include_core", False))
    if not _table_exists(manager, "Libraries"):
        return []
    if include_core:
        manager.cursor.execute("SELECT library_id FROM Libraries ORDER BY library_id")
    else:
        manager.cursor.execute("SELECT library_id FROM Libraries WHERE parent_id IS NOT NULL ORDER BY library_id")
    return [int(row[0]) for row in (manager.cursor.fetchall() or [])]


def _running_task_count(manager: DBManager) -> int:
    if not _table_exists(manager, "Tasks"):
        return 0
    manager.cursor.execute(
        """
        SELECT COUNT(*)
        FROM Tasks
        WHERE status = 'R'
          AND (job_type IS NULL OR job_type != 0)
        """
    )
    row = manager.cursor.fetchone()
    return int(row[0] if row else 0)


def _plan_and_execute_assemblies(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    accessions = _resolve_accessions(manager, args)
    result = PurgeResult(subject="assemblies", dry_run=not bool(args.apply))
    if not accessions:
        result.notes.append("No assemblies matched.")
        return result

    include_metadata = bool(getattr(args, "include_metadata", False))
    run_ids = _collect_busco_run_ids_for_accessions(manager, accessions)
    busco_artifact_ids = _collect_busco_artifact_ids_for_run_ids(manager, run_ids)
    genome_artifact_ids = _collect_genome_artifact_ids(manager, accessions)
    ortho_ids = _fetch_ids_in(manager, "OrthoFinder_Accessions", "orthofinder_id", "accession", accessions)
    busco_run_ids = _fetch_ids_in(manager, "BUSCO_Runs", "run_id", "accession", accessions)
    table_counts: List[Tuple[str, int]] = [
        ("OrthoFinder_Accessions", _count_where_in(manager, "OrthoFinder_Accessions", "orthofinder_id", ortho_ids)),
        ("OrthoFinder_Results", _count_where_in(manager, "OrthoFinder_Results", "orthofinder_id", ortho_ids)),
        ("BUSCO_Run_Family_Artifacts", _count_where_in(manager, "BUSCO_Run_Family_Artifacts", "run_id", busco_run_ids)),
        ("BUSCO_Run_Family_Locations", _count_where_in(manager, "BUSCO_Run_Family_Locations", "run_id", busco_run_ids)),
        ("BUSCO_Run_Family_Data", _count_where_in(manager, "BUSCO_Run_Family_Data", "run_id", busco_run_ids)),
        ("BUSCO_Runs", _count_where_in(manager, "BUSCO_Runs", "accession", accessions)),
        ("BUSCO_Adjusted_Results", _count_where_in(manager, "BUSCO_Adjusted_Results", "accession", accessions)),
        ("Decontamination_Busco_Votes", _count_where_in(manager, "Decontamination_Busco_Votes", "accession", accessions)),
        ("Decontamination_Busco_Copy_Votes", _count_where_in(manager, "Decontamination_Busco_Copy_Votes", "accession", accessions)),
        ("Decontamination_Summary", _count_where_in(manager, "Decontamination_Summary", "accession", accessions)),
        ("Paralog_Filtering", _count_where_in(manager, "Paralog_Filtering", "accession", accessions)),
        ("Paralog_Filtering_Copy", _count_where_in(manager, "Paralog_Filtering_Copy", "accession", accessions)),
        ("BUSCO_Family_Locations", _count_where_in(manager, "BUSCO_Family_Locations", "accession", accessions)),
        ("BUSCO_Family_Data", _count_where_in(manager, "BUSCO_Family_Data", "accession", accessions)),
        ("BUSCO_Results", _count_where_in(manager, "BUSCO_Results", "accession", accessions)),
        ("BUSCO_Primary", _count_where_in(manager, "BUSCO_Primary", "accession", accessions)),
        ("BUSCO_Adjusted_Results", _count_where_in(manager, "BUSCO_Adjusted_Results", "accession", accessions)),
        ("BUSCO_Run_Family_Data", _count_where_in(manager, "BUSCO_Run_Family_Data", "run_id", run_ids)),
        ("BUSCO_Run_Family_Locations", _count_where_in(manager, "BUSCO_Run_Family_Locations", "run_id", run_ids)),
        ("BUSCO_Run_Family_Artifacts", _count_where_in(manager, "BUSCO_Run_Family_Artifacts", "run_id", run_ids)),
        ("BUSCO_Runs", _count_where_in(manager, "BUSCO_Runs", "run_id", run_ids)),
        ("Artifacts(busco_run)", _count_where_in(manager, "Artifacts", "artifact_id", busco_artifact_ids)),
        ("Artifacts(genome)", _count_where_in(manager, "Artifacts", "artifact_id", genome_artifact_ids)),
        ("Proteome_BlastDBs", _count_where_in(manager, "Proteome_BlastDBs", "accession", accessions)),
        ("Reference_Assemblies", _count_where_in(manager, "Reference_Assemblies", "accession", accessions)),
        ("Hidden_Genomes", _count_where_in(manager, "Hidden_Genomes", "accession", accessions)),
        ("Genome", _count_where_in(manager, "Genome", "accession", accessions)),
    ]
    if include_metadata:
        table_counts.append(("Assembly", _count_where_in(manager, "Assembly", "accession", accessions)))
    result.counts = [(name, count) for name, count in table_counts if count > 0]

    if args.delete_files:
        result.files = sorted(_collect_assembly_purge_files(manager, accessions, run_ids=run_ids))

    if result.dry_run:
        if include_metadata:
            result.notes.append("Mode: purge all knowledge of selected assemblies.")
        else:
            result.notes.append("Mode: purge downloaded/workflow state; metadata retained.")
        result.notes.append(f"Matched assemblies: {len(accessions)}")
        return result

    previous_locations = _current_genome_locations(manager, accessions)
    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        _delete_where_in(manager, "OrthoFinder_Accessions", "orthofinder_id", ortho_ids)
        _delete_where_in(manager, "OrthoFinder_Results", "orthofinder_id", ortho_ids)
        _delete_where_in(manager, "BUSCO_Adjusted_Results", "accession", accessions)
        _delete_where_in(manager, "BUSCO_Runs", "accession", accessions)
        _delete_where_in(manager, "Decontamination_Busco_Votes", "accession", accessions)
        _delete_where_in(manager, "Decontamination_Busco_Copy_Votes", "accession", accessions)
        _delete_where_in(manager, "Decontamination_Summary", "accession", accessions)
        _delete_where_in(manager, "Paralog_Filtering", "accession", accessions)
        _delete_where_in(manager, "Paralog_Filtering_Copy", "accession", accessions)
        _delete_where_in(manager, "BUSCO_Primary", "accession", accessions)
        _delete_where_in(manager, "BUSCO_Adjusted_Results", "accession", accessions)
        _delete_where_in(manager, "BUSCO_Run_Family_Artifacts", "run_id", run_ids)
        _delete_where_in(manager, "Artifacts", "artifact_id", busco_artifact_ids)
        _delete_where_in(manager, "BUSCO_Run_Family_Locations", "run_id", run_ids)
        _delete_where_in(manager, "BUSCO_Run_Family_Data", "run_id", run_ids)
        _delete_where_in(manager, "BUSCO_Runs", "run_id", run_ids)
        _delete_where_in(manager, "BUSCO_Family_Locations", "accession", accessions)
        _delete_where_in(manager, "BUSCO_Family_Data", "accession", accessions)
        _delete_where_in(manager, "BUSCO_Results", "accession", accessions)
        _delete_where_in(manager, "Proteome_BlastDBs", "accession", accessions)
        _delete_where_in(manager, "Reference_Assemblies", "accession", accessions)
        _delete_where_in(manager, "Hidden_Genomes", "accession", accessions)
        _delete_where_in(manager, "Artifacts", "artifact_id", genome_artifact_ids)
        if include_metadata:
            _delete_where_in(manager, "Genome", "accession", accessions)
            _delete_where_in(manager, "Assembly", "accession", accessions)
            lingering = _verify_no_accession_references(manager, accessions)
            if lingering:
                raise ValueError("Remaining accession references after metadata purge: " + ", ".join(lingering))
        else:
            for accession in accessions:
                _recompute_genome_state(manager, accession, previous_locations.get(accession))
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise

    if include_metadata:
        result.notes.append(f"Purged assemblies: {len(accessions)}")
    else:
        result.notes.append(f"Purged downloaded/workflow state for assemblies: {len(accessions)}")
    if args.delete_files and result.files:
        deleted, skipped, failed, failures = _delete_files(manager, result.files)
        result.deleted_files = deleted
        result.skipped_files = skipped
        result.failed_files = failed
        if not include_metadata:
            manager.conn.execute("BEGIN IMMEDIATE")
            try:
                for accession in accessions:
                    _recompute_genome_state(manager, accession, previous_locations.get(accession))
                manager.conn.commit()
            except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
                manager.conn.rollback()
                raise
        result.notes.extend(failures[:5])
        if failures and len(failures) > 5:
            result.notes.append(f"... {len(failures) - 5} more file deletion errors.")
    return result


def _build_decontam_where(
    *,
    run_ids: Sequence[str],
    busco_run_ids: Sequence[int],
    target_library_ids: Sequence[int],
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []
    if run_ids:
        placeholders = ",".join("?" for _ in run_ids)
        clauses.append(f"run_id IN ({placeholders})")
        params.extend([str(run_id) for run_id in run_ids])
    if busco_run_ids:
        placeholders = ",".join("?" for _ in busco_run_ids)
        clauses.append(f"busco_run_id IN ({placeholders})")
        params.extend([int(run_id) for run_id in busco_run_ids])
    if target_library_ids:
        placeholders = ",".join("?" for _ in target_library_ids)
        clauses.append(f"target_library_id IN ({placeholders})")
        params.extend(target_library_ids)
    where_sql = " AND ".join(clauses) if clauses else ""
    return where_sql, params


def _plan_and_execute_decontamination(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="decontamination", dry_run=not bool(args.apply))
    use_selector = _selector_requested(args)
    accessions: List[str] = []
    if use_selector:
        accessions = _resolve_accessions(manager, args)

    run_ids = [str(value).strip() for value in (getattr(args, "run_ids", None) or []) if str(value).strip()]
    busco_run_ids = expand_busco_run_id_variables(manager, getattr(args, "busco_run_ids", None) or [])
    target_library_ids: List[int] = []
    if getattr(args, "libraries", None):
        target_library_ids = _resolve_library_ids(manager, args.libraries)

    if not (accessions or run_ids or busco_run_ids or target_library_ids or args.all):
        raise ValueError("Provide --run-id/--busco-run-id, library filter, selectors, or --all.")

    where_sql, params = _build_decontam_where(
        run_ids=run_ids,
        busco_run_ids=busco_run_ids,
        target_library_ids=target_library_ids,
    )
    if args.all and not where_sql and not accessions:
        where_sql, params = "", []
    if not (args.all or where_sql or accessions):
        raise ValueError("Provide --run-id/--busco-run-id, library filter, selectors, or --all.")

    full_where, full_params = _build_full_accession_where(accessions, extra_where=where_sql, extra_params=params)
    affected = (
        _fetch_affected_library_accessions(manager, "Decontamination_Busco_Votes", full_where, full_params)
        + _fetch_affected_library_accessions(manager, "Decontamination_Summary", full_where, full_params)
    )
    affected = list(dict.fromkeys(affected))

    vote_count = _count_where(manager, "Decontamination_Busco_Votes", full_where, full_params)
    copy_vote_count = _count_where(manager, "Decontamination_Busco_Copy_Votes", full_where, full_params)
    summary_count = _count_where(manager, "Decontamination_Summary", full_where, full_params)
    result.counts.extend(
        [
            ("Decontamination_Busco_Votes", vote_count),
            ("Decontamination_Busco_Copy_Votes", copy_vote_count),
            ("Decontamination_Summary", summary_count),
        ]
    )

    delete_runs = not bool(getattr(args, "keep_runs", False))
    matched_run_ids = list(
        dict.fromkeys(
            [str(value) for value in _fetch_distinct_column_where(manager, "Decontamination_Busco_Votes", "run_id", full_where, full_params)]
            + [str(value) for value in _fetch_distinct_column_where(manager, "Decontamination_Busco_Copy_Votes", "run_id", full_where, full_params)]
            + [str(value) for value in _fetch_distinct_column_where(manager, "Decontamination_Summary", "run_id", full_where, full_params)]
        )
    )
    deletable_run_ids: List[str] = []
    if delete_runs:
        if result.dry_run and matched_run_ids:
            if not full_where:
                deletable_run_ids = list(matched_run_ids)
            else:
                for run_id in matched_run_ids:
                    if _count_where(manager, "Decontamination_Busco_Votes", "run_id = ? AND NOT (" + full_where + ")", [run_id, *full_params]) > 0:
                        continue
                    if _count_where(manager, "Decontamination_Busco_Copy_Votes", "run_id = ? AND NOT (" + full_where + ")", [run_id, *full_params]) > 0:
                        continue
                    if _count_where(manager, "Decontamination_Summary", "run_id = ? AND NOT (" + full_where + ")", [run_id, *full_params]) > 0:
                        continue
                    deletable_run_ids.append(run_id)
                deletable_run_ids = list(dict.fromkeys(deletable_run_ids))
        result.counts.append(("Decontamination_Runs", len(deletable_run_ids)))

    result.counts = [(name, count) for name, count in result.counts if count > 0]
    if result.dry_run:
        return result

    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        _delete_where(manager, "Decontamination_Busco_Votes", full_where, full_params)
        _delete_where(manager, "Decontamination_Busco_Copy_Votes", full_where, full_params)
        _delete_where(manager, "Decontamination_Summary", full_where, full_params)
        if delete_runs:
            deletable_run_ids = _filter_empty_runs(
                manager,
                matched_run_ids,
                row_tables=("Decontamination_Busco_Votes", "Decontamination_Busco_Copy_Votes", "Decontamination_Summary"),
            )
            if deletable_run_ids:
                _delete_where_in(manager, "Decontamination_Runs", "run_id", deletable_run_ids)
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise
    for library_id, accession in affected:
        manager.busco.invalidate_adjusted_results_for_library(
            int(library_id),
            accessions=[str(accession)],
            reason="purge-decontamination",
        )
    return result


def _build_busco_where(
    *,
    library_ids: Sequence[int],
) -> Tuple[str, List[Any]]:
    clauses: List[str] = []
    params: List[Any] = []
    if library_ids:
        placeholders = ",".join("?" for _ in library_ids)
        clauses.append(f"library_id IN ({placeholders})")
        params.extend(library_ids)
    return " AND ".join(clauses), params


def _plan_and_execute_busco(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="busco", dry_run=not bool(args.apply))
    accessions: List[str] = []
    if _selector_requested(args):
        accessions = _resolve_accessions(manager, args)
    library_ids: List[int] = []
    if getattr(args, "libraries", None):
        library_ids = _resolve_library_ids(manager, args.libraries)
    explicit_run_ids = expand_busco_run_id_variables(manager, getattr(args, "run_ids", None) or [])
    run_status = str(getattr(args, "status", "") or "").strip().lower() or None
    protein_only = bool(getattr(args, "protein_only", False))

    if not (accessions or library_ids or explicit_run_ids or run_status or args.all):
        raise ValueError("Provide selectors, library filter, or --all.")

    # Mode-filtered purges must operate on BUSCO_Runs directly so they do not
    # delete shared legacy summary tables for other run modes.
    run_targeted = bool(explicit_run_ids or run_status or protein_only)
    if run_targeted:
        clauses = ["1=1"]
        params: List[Any] = []
        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            clauses.append(f"accession IN ({placeholders})")
            params.extend([str(acc) for acc in accessions])
        if library_ids:
            placeholders = ",".join("?" for _ in library_ids)
            clauses.append(f"library_id IN ({placeholders})")
            params.extend(library_ids)
        if explicit_run_ids:
            placeholders = ",".join("?" for _ in explicit_run_ids)
            clauses.append(f"run_id IN ({placeholders})")
            params.extend(explicit_run_ids)
        if run_status:
            clauses.append("LOWER(COALESCE(status, '')) = ?")
            params.append(run_status)
        if protein_only:
            clauses.append("LOWER(COALESCE(input_mode, '')) = 'protein'")
        run_where = " AND ".join(clauses)

        run_rows = []
        if _table_exists(manager, "BUSCO_Runs"):
            manager.cursor.execute(
                f"SELECT run_id, accession, library_id, result_dir, status FROM BUSCO_Runs WHERE {run_where} ORDER BY run_id",
                tuple(params),
            )
            run_rows = manager.cursor.fetchall() or []
        run_ids = [int(row[0]) for row in run_rows if row and row[0] is not None]
        if not run_ids:
            result.notes.append("No BUSCO runs matched.")
            return result
        affected_pairs = sorted(
            {
                (str(row[1]), int(row[2]))
                for row in run_rows
                if row and row[1] is not None and row[2] is not None
            }
        )

        if args.delete_files:
            result.files = sorted(_collect_busco_locations_for_run_ids(manager, run_ids))

        if _table_exists(manager, "BUSCO_Run_Family_Data"):
            result.counts.append(("BUSCO_Run_Family_Data", _count_where_in(manager, "BUSCO_Run_Family_Data", "run_id", run_ids)))
        if _table_exists(manager, "BUSCO_Run_Family_Locations"):
            result.counts.append(("BUSCO_Run_Family_Locations", _count_where_in(manager, "BUSCO_Run_Family_Locations", "run_id", run_ids)))
        if _table_exists(manager, "BUSCO_Run_Family_Artifacts"):
            result.counts.append(("BUSCO_Run_Family_Artifacts", _count_where_in(manager, "BUSCO_Run_Family_Artifacts", "run_id", run_ids)))
        if _table_exists(manager, "BUSCO_Primary"):
            result.counts.append(("BUSCO_Primary", _count_where_in(manager, "BUSCO_Primary", "run_id", run_ids)))
        result.counts.append(("BUSCO_Runs", len(run_ids)))
        result.counts = [(name, count) for name, count in result.counts if count > 0]
        result.notes.append(
            f"Matched BUSCO run ids: {', '.join(str(run_id) for run_id in run_ids[:20])}"
            + (f" ... and {len(run_ids) - 20} more" if len(run_ids) > 20 else "")
        )
        result.notes.append(
            "Run-targeted BUSCO purge removes matching BUSCO_Runs and run-owned tables/files only; legacy BUSCO summary tables are unchanged."
        )

        if result.dry_run:
            return result

        manager.conn.execute("BEGIN IMMEDIATE")
        try:
            _delete_where_in(manager, "BUSCO_Runs", "run_id", run_ids)
            manager.conn.commit()
        except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
            manager.conn.rollback()
            raise
        for accession, library_id in affected_pairs:
            refreshed = manager.busco.refresh_auto_primary_runs_for_accession(
                accession,
                library_id,
                updated_by="purge-busco",
                policy="auto_best",
            )
            if refreshed.get("default") is None:
                manager.busco.delete_records(accession, library_id)

        if args.delete_files and result.files:
            deleted, skipped, failed, failures = _delete_files(manager, result.files)
            result.deleted_files = deleted
            result.skipped_files = skipped
            result.failed_files = failed
            result.notes.extend(failures[:5])
            if failures and len(failures) > 5:
                result.notes.append(f"... {len(failures) - 5} more file deletion errors.")
        return result

    where_sql, params = _build_busco_where(library_ids=library_ids)
    if args.all and not where_sql and not accessions:
        where_sql, params = "", []
    if not (args.all or where_sql or accessions):
        raise ValueError("Provide selectors, library filter, or --all.")

    matched_run_ids = _fetch_distinct_column_where(manager, "BUSCO_Runs", "run_id", where_sql, params)
    if accessions:
        matched_run_ids = [
            int(run_id)
            for run_id in matched_run_ids
            if run_id is not None
        ]
        if matched_run_ids:
            filtered_run_ids: List[int] = []
            for chunk in _chunked(matched_run_ids):
                placeholders = ",".join("?" for _ in chunk)
                acc_placeholders = ",".join("?" for _ in accessions)
                manager.cursor.execute(
                    f"SELECT run_id FROM BUSCO_Runs WHERE run_id IN ({placeholders}) AND accession IN ({acc_placeholders})",
                    tuple(chunk) + tuple(accessions),
                )
                filtered_run_ids.extend(int(row[0]) for row in (manager.cursor.fetchall() or []) if row and row[0] is not None)
            matched_run_ids = list(dict.fromkeys(filtered_run_ids))
    else:
        matched_run_ids = [int(run_id) for run_id in matched_run_ids if run_id is not None]
    affected_pairs: List[Tuple[str, int]] = []
    if matched_run_ids and _table_exists(manager, "BUSCO_Runs"):
        for chunk in _chunked(matched_run_ids):
            placeholders = ",".join("?" for _ in chunk)
            manager.cursor.execute(
                f"SELECT accession, library_id FROM BUSCO_Runs WHERE run_id IN ({placeholders})",
                tuple(chunk),
            )
            for accession, library_id in (manager.cursor.fetchall() or []):
                if accession is None or library_id is None:
                    continue
                affected_pairs.append((str(accession), int(library_id)))
    affected_pairs = sorted(set(affected_pairs))
    if args.delete_files:
        result.files = sorted(_collect_busco_locations_for_run_ids(manager, matched_run_ids))

    result.counts = [
        (
            "BUSCO_Family_Locations",
            _count_filtered_by_accessions(
                manager,
                "BUSCO_Family_Locations",
                accessions,
                extra_where=where_sql,
                extra_params=params,
            ),
        ),
        (
            "BUSCO_Family_Data",
            _count_filtered_by_accessions(
                manager,
                "BUSCO_Family_Data",
                accessions,
                extra_where=where_sql,
                extra_params=params,
            ),
        ),
        (
            "BUSCO_Results",
            _count_filtered_by_accessions(
                manager,
                "BUSCO_Results",
                accessions,
                extra_where=where_sql,
                extra_params=params,
            ),
        ),
        (
            "BUSCO_Runs",
            _count_filtered_by_accessions(
                manager,
                "BUSCO_Runs",
                accessions,
                extra_where=where_sql,
                extra_params=params,
            ),
        ),
    ]
    if library_ids and not accessions and _table_exists(manager, "BUSCO_descriptions"):
        desc_sql = f"library_id IN ({','.join('?' for _ in library_ids)})"
        desc_count = _count_where(manager, "BUSCO_descriptions", desc_sql, library_ids)
        if desc_count > 0:
            result.counts.append(("BUSCO_descriptions", desc_count))
    result.counts = [(name, count) for name, count in result.counts if count > 0]

    if result.dry_run:
        return result

    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        _delete_filtered_by_accessions(
            manager,
            "BUSCO_Family_Locations",
            accessions,
            extra_where=where_sql,
            extra_params=params,
        )
        _delete_filtered_by_accessions(
            manager,
            "BUSCO_Family_Data",
            accessions,
            extra_where=where_sql,
            extra_params=params,
        )
        _delete_filtered_by_accessions(
            manager,
            "BUSCO_Results",
            accessions,
            extra_where=where_sql,
            extra_params=params,
        )
        _delete_filtered_by_accessions(
            manager,
            "BUSCO_Runs",
            accessions,
            extra_where=where_sql,
            extra_params=params,
        )
        if library_ids and not accessions and _table_exists(manager, "BUSCO_descriptions"):
            where_desc = f"library_id IN ({','.join('?' for _ in library_ids)})"
            _delete_where(manager, "BUSCO_descriptions", where_desc, library_ids)
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise
    for accession, library_id in affected_pairs:
        refreshed = manager.busco.refresh_auto_primary_runs_for_accession(
            accession,
            library_id,
            updated_by="purge-busco",
            policy="auto_best",
        )
        if refreshed.get("default") is None:
            manager.busco.delete_records(accession, library_id)

    if args.delete_files and result.files:
        deleted, skipped, failed, failures = _delete_files(manager, result.files)
        result.deleted_files = deleted
        result.skipped_files = skipped
        result.failed_files = failed
        result.notes.extend(failures[:5])
        if failures and len(failures) > 5:
            result.notes.append(f"... {len(failures) - 5} more file deletion errors.")
    return result


def _plan_and_execute_busco_primary(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="busco-primary", dry_run=not bool(args.apply))
    accessions: List[str] = []
    if _selector_requested(args):
        accessions = _resolve_accessions(manager, args)
    library_ids: List[int] = []
    if getattr(args, "libraries", None):
        library_ids = _resolve_library_ids(manager, args.libraries)

    clauses = ["LOWER(COALESCE(policy, '')) = 'manual_override'"]
    params: List[Any] = []
    if accessions:
        placeholders = ",".join("?" for _ in accessions)
        clauses.append(f"accession IN ({placeholders})")
        params.extend([str(acc) for acc in accessions])
    if library_ids:
        placeholders = ",".join("?" for _ in library_ids)
        clauses.append(f"library_id IN ({placeholders})")
        params.extend(library_ids)
    if not (accessions or library_ids or args.all):
        raise ValueError("Provide selectors, library filter, or --all.")
    where_sql = " AND ".join(clauses)
    count = _count_where(manager, "BUSCO_Primary", where_sql, params)
    result.counts = [("BUSCO_Primary(manual_override)", count)] if count > 0 else []
    if count == 0:
        result.notes.append("No manual BUSCO primary overrides matched.")
        return result

    manager.cursor.execute(
        f"""
        SELECT accession, library_id, purpose, run_id
        FROM BUSCO_Primary
        WHERE {where_sql}
        ORDER BY accession, library_id, purpose
        """,
        tuple(params),
    )
    rows = manager.cursor.fetchall() or []
    for accession, library_id, purpose, run_id in rows[:20]:
        result.notes.append(f"{accession} lib={library_id} purpose={purpose} run_id={run_id}")
    if len(rows) > 20:
        result.notes.append(f"... {len(rows) - 20} more manual BUSCO primary override rows.")
    if result.dry_run:
        return result

    affected = sorted({(str(accession), int(library_id)) for accession, library_id, _purpose, _run_id in rows})
    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        _delete_where(manager, "BUSCO_Primary", where_sql, params)
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise

    reassigned = 0
    for accession, library_id in affected:
        refreshed = manager.busco.refresh_auto_primary_runs_for_accession(
            accession,
            library_id,
            updated_by="purge-busco-primary",
            policy="auto_best",
        )
        if any(value is not None for value in refreshed.values()):
            reassigned += 1
        else:
            result.notes.append(f"{accession} lib={library_id}: no usable automatic replacement run.")
    result.notes.append(f"Unlocked manual BUSCO primary overrides for {len(affected)} accession/library pairs.")
    result.notes.append(f"Reassigned automatic primaries for {reassigned} accession/library pairs.")
    return result


def _plan_and_execute_hidden_paralog(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="hidden-paralog", dry_run=not bool(args.apply))
    accessions: List[str] = []
    if _selector_requested(args):
        accessions = _resolve_accessions(manager, args)
    target_library_ids: List[int] = []
    if getattr(args, "libraries", None):
        target_library_ids = _resolve_library_ids(manager, args.libraries)
    busco_library_ids: List[int] = []
    if getattr(args, "busco_libraries", None):
        busco_library_ids = _resolve_library_ids(manager, args.busco_libraries)
    run_ids = [str(value).strip() for value in (getattr(args, "run_ids", None) or []) if str(value).strip()]
    busco_run_ids = expand_busco_run_id_variables(manager, getattr(args, "busco_run_ids", None) or [])

    supports_target = _column_exists(manager, "Paralog_Filtering", "target_library_id")
    supports_busco_run = _column_exists(manager, "Paralog_Filtering", "busco_run_id")
    if not (accessions or target_library_ids or busco_library_ids or run_ids or busco_run_ids or args.all):
        raise ValueError("Provide selectors/library/run filters or --all.")

    clauses: List[str] = []
    params: List[Any] = []
    if busco_library_ids:
        placeholders = ",".join("?" for _ in busco_library_ids)
        clauses.append(f"library_id IN ({placeholders})")
        params.extend(busco_library_ids)
    if target_library_ids:
        placeholders = ",".join("?" for _ in target_library_ids)
        if supports_target:
            clauses.append(f"target_library_id IN ({placeholders})")
        else:
            clauses.append(f"library_id IN ({placeholders})")
        params.extend(target_library_ids)
    if run_ids:
        placeholders = ",".join("?" for _ in run_ids)
        clauses.append(f"run_id IN ({placeholders})")
        params.extend([str(run_id) for run_id in run_ids])
    if busco_run_ids:
        if not supports_busco_run:
            raise ValueError("This database does not support BUSCO-run-linked paralog purge.")
        placeholders = ",".join("?" for _ in busco_run_ids)
        clauses.append(f"busco_run_id IN ({placeholders})")
        params.extend([int(run_id) for run_id in busco_run_ids])
    where_sql = " AND ".join(clauses)
    if args.all and not where_sql and not accessions:
        where_sql, params = "", []
    if not (args.all or where_sql or accessions):
        raise ValueError("Provide selectors/library/run filters or --all.")

    full_where, full_params = _build_full_accession_where(accessions, extra_where=where_sql, extra_params=params)
    main_count = _count_where(manager, "Paralog_Filtering", full_where, full_params)
    copy_count = _count_where(manager, "Paralog_Filtering_Copy", full_where, full_params)
    if main_count > 0:
        result.counts.append(("Paralog_Filtering", main_count))
    if copy_count > 0:
        result.counts.append(("Paralog_Filtering_Copy", copy_count))

    affected = (
        _fetch_affected_library_accessions(manager, "Paralog_Filtering", full_where, full_params)
        + _fetch_affected_library_accessions(manager, "Paralog_Filtering_Copy", full_where, full_params)
    )
    affected = list(dict.fromkeys(affected))
    matched_run_ids = list(
        dict.fromkeys(
            [str(value) for value in _fetch_distinct_column_where(manager, "Paralog_Filtering", "run_id", full_where, full_params)]
            + [str(value) for value in _fetch_distinct_column_where(manager, "Paralog_Filtering_Copy", "run_id", full_where, full_params)]
        )
    )
    delete_runs = not bool(getattr(args, "keep_runs", False))
    deletable_run_ids: List[str] = []
    if delete_runs:
        if not full_where:
            deletable_run_ids = list(matched_run_ids)
        else:
            for run_id in matched_run_ids:
                if _count_where(manager, "Paralog_Filtering", "run_id = ? AND NOT (" + full_where + ")", [run_id, *full_params]) > 0:
                    continue
                if _count_where(manager, "Paralog_Filtering_Copy", "run_id = ? AND NOT (" + full_where + ")", [run_id, *full_params]) > 0:
                    continue
                deletable_run_ids.append(run_id)
        if deletable_run_ids:
            result.counts.append(("Paralog_Filtering_Runs", len(deletable_run_ids)))
    if result.dry_run:
        return result

    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        _delete_where(manager, "Paralog_Filtering", full_where, full_params)
        _delete_where(manager, "Paralog_Filtering_Copy", full_where, full_params)
        if delete_runs:
            deletable_run_ids = _filter_empty_runs(
                manager,
                matched_run_ids,
                row_tables=("Paralog_Filtering", "Paralog_Filtering_Copy"),
            )
            if deletable_run_ids:
                _delete_where_in(manager, "Paralog_Filtering_Runs", "run_id", deletable_run_ids)
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise
    for library_id, accession in affected:
        manager.busco.invalidate_adjusted_results_for_library(
            int(library_id),
            accessions=[str(accession)],
            reason="purge-hidden-paralog",
        )
    return result


def _expand_library_descendants(manager: DBManager, library_ids: Sequence[int]) -> List[int]:
    if not library_ids:
        return []
    placeholders = ",".join("?" for _ in library_ids)
    manager.cursor.execute(
        f"""
        WITH RECURSIVE descendants(library_id) AS (
            SELECT library_id FROM Libraries WHERE library_id IN ({placeholders})
            UNION ALL
            SELECT l.library_id
            FROM Libraries l
            JOIN descendants d ON l.parent_id = d.library_id
        )
        SELECT DISTINCT library_id FROM descendants ORDER BY library_id
        """,
        tuple(library_ids),
    )
    return [int(row[0]) for row in (manager.cursor.fetchall() or [])]


def _plan_and_execute_libraries(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="libraries", dry_run=not bool(args.apply))
    library_ids = _resolve_subject_libraries(manager, args)
    if not library_ids:
        result.notes.append("No libraries matched.")
        return result

    if getattr(args, "recursive", False):
        library_ids = _expand_library_descendants(manager, library_ids)

    drop_library = bool(getattr(args, "drop_library", False))
    supports_target = _column_exists(manager, "Paralog_Filtering", "target_library_id")

    if drop_library and not getattr(args, "recursive", False):
        placeholders = ",".join("?" for _ in library_ids)
        manager.cursor.execute(
            f"SELECT library_id, library_name, parent_id FROM Libraries WHERE parent_id IN ({placeholders}) ORDER BY library_id",
            tuple(library_ids),
        )
        children = [row for row in (manager.cursor.fetchall() or []) if int(row[0]) not in set(library_ids)]
        if children:
            names = ", ".join(f"{row[1]}({row[0]})" for row in children[:5])
            raise ValueError(
                f"Cannot drop library with children outside selection: {names}. Use --recursive."
            )

    if args.delete_files:
        result.files = sorted(
            _collect_locations_for_libraries(
                manager,
                library_ids,
                include_library_locations=drop_library,
            )
        )

    counts: List[Tuple[str, int]] = [
        ("Reference_Assemblies", _count_where_in(manager, "Reference_Assemblies", "library_id", library_ids)),
        ("Proteome_BlastDBs", _count_where_in(manager, "Proteome_BlastDBs", "library_id", library_ids)),
        ("BUSCO_Family_Locations", _count_where_in(manager, "BUSCO_Family_Locations", "library_id", library_ids)),
        ("BUSCO_Family_Data", _count_where_in(manager, "BUSCO_Family_Data", "library_id", library_ids)),
        ("BUSCO_Results", _count_where_in(manager, "BUSCO_Results", "library_id", library_ids)),
        ("BUSCO_descriptions", _count_where_in(manager, "BUSCO_descriptions", "library_id", library_ids)),
        ("OrthoFinder_Results", _count_where_in(manager, "OrthoFinder_Results", "library_id", library_ids)),
    ]
    ortho_ids = _fetch_ids_in(manager, "OrthoFinder_Results", "orthofinder_id", "library_id", library_ids)
    counts.append(("OrthoFinder_Accessions", _count_where_in(manager, "OrthoFinder_Accessions", "orthofinder_id", ortho_ids)))

    if supports_target:
        counts.append(("Paralog_Filtering(target)", _count_where_in(manager, "Paralog_Filtering", "target_library_id", library_ids)))
    else:
        counts.append(("Paralog_Filtering(target)", _count_where_in(manager, "Paralog_Filtering", "library_id", library_ids)))
    counts.append(("Decontamination_Busco_Votes(target)", _count_where_in(manager, "Decontamination_Busco_Votes", "target_library_id", library_ids)))
    counts.append(("Decontamination_Summary(target)", _count_where_in(manager, "Decontamination_Summary", "target_library_id", library_ids)))
    counts.append(("Decontamination_Runs(target)", _count_where_in(manager, "Decontamination_Runs", "target_library_id", library_ids)))

    if drop_library:
        counts.append(("Paralog_Filtering(busco)", _count_where_in(manager, "Paralog_Filtering", "library_id", library_ids)))
        counts.append(("Decontamination_Busco_Votes(busco)", _count_where_in(manager, "Decontamination_Busco_Votes", "busco_library_id", library_ids)))
        counts.append(("Decontamination_Summary(busco)", _count_where_in(manager, "Decontamination_Summary", "busco_library_id", library_ids)))
        counts.append(("Decontamination_Runs(busco)", _count_where_in(manager, "Decontamination_Runs", "busco_library_id", library_ids)))
        counts.append(("Libraries", _count_where_in(manager, "Libraries", "library_id", library_ids)))

    result.counts = [(name, count) for name, count in counts if count > 0]
    ids = ", ".join(str(v) for v in library_ids)
    result.notes.append(f"Matched libraries count={len(library_ids)} ids=[{ids}]")
    if not drop_library:
        result.notes.append(
            "Library definitions are retained by default; use --drop-library to delete the matched Libraries rows."
        )
        if not result.counts:
            result.notes.append("No dependent library data remains to purge.")
    if result.dry_run:
        return result

    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        if supports_target:
            _delete_where_in(manager, "Paralog_Filtering", "target_library_id", library_ids)
        else:
            _delete_where_in(manager, "Paralog_Filtering", "library_id", library_ids)
        _delete_where_in(manager, "Decontamination_Busco_Votes", "target_library_id", library_ids)
        _delete_where_in(manager, "Decontamination_Summary", "target_library_id", library_ids)
        _delete_where_in(manager, "Decontamination_Runs", "target_library_id", library_ids)

        if drop_library:
            _delete_where_in(manager, "Paralog_Filtering", "library_id", library_ids)
            _delete_where_in(manager, "Decontamination_Busco_Votes", "busco_library_id", library_ids)
            _delete_where_in(manager, "Decontamination_Summary", "busco_library_id", library_ids)
            _delete_where_in(manager, "Decontamination_Runs", "busco_library_id", library_ids)

        _delete_where_in(manager, "Reference_Assemblies", "library_id", library_ids)
        _delete_where_in(manager, "Proteome_BlastDBs", "library_id", library_ids)
        _delete_where_in(manager, "OrthoFinder_Accessions", "orthofinder_id", ortho_ids)
        _delete_where_in(manager, "OrthoFinder_Results", "library_id", library_ids)
        _delete_where_in(manager, "BUSCO_Family_Locations", "library_id", library_ids)
        _delete_where_in(manager, "BUSCO_Family_Data", "library_id", library_ids)
        _delete_where_in(manager, "BUSCO_Results", "library_id", library_ids)
        _delete_where_in(manager, "BUSCO_descriptions", "library_id", library_ids)

        if drop_library:
            _delete_where_in(manager, "Libraries", "library_id", library_ids)

        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise

    if args.delete_files and result.files:
        deleted, skipped, failed, failures = _delete_files(manager, result.files)
        result.deleted_files = deleted
        result.skipped_files = skipped
        result.failed_files = failed
        result.notes.extend(failures[:5])
        if failures and len(failures) > 5:
            result.notes.append(f"... {len(failures) - 5} more file deletion errors.")
    return result


def _is_system_variable(name: str) -> bool:
    return is_system_variable(name)


def _plan_and_execute_variables(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="variables", dry_run=not bool(args.apply))
    if not _table_exists(manager, "Environment_Variables"):
        result.notes.append("Environment_Variables table not present.")
        return result

    manager.cursor.execute("SELECT var_name FROM Environment_Variables ORDER BY var_name")
    names = [str(row[0]) for row in (manager.cursor.fetchall() or [])]
    requested = set(str(v).strip() for v in (getattr(args, "var_names", None) or []) if str(v).strip())
    prefix = str(getattr(args, "prefix", "") or "").strip()
    custom_only = bool(getattr(args, "custom_only", True))

    selected: List[str] = []
    for name in names:
        if requested and name not in requested:
            continue
        if prefix and not name.startswith(prefix):
            continue
        if custom_only and _is_system_variable(name):
            continue
        selected.append(name)

    if args.all and not requested and not prefix:
        selected = names if not custom_only else [name for name in names if not _is_system_variable(name)]

    result.counts = [("Environment_Variables", len(selected))]
    if result.dry_run:
        result.notes.append("Variables to purge: " + (", ".join(selected[:10]) if selected else "(none)"))
        if len(selected) > 10:
            result.notes.append(f"... and {len(selected) - 10} more")
        return result

    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        _delete_where_in(manager, "Environment_Variables", "var_name", selected)
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise
    return result


def _plan_and_execute_roots(manager: DBManager, args: argparse.Namespace) -> PurgeResult:
    result = PurgeResult(subject="roots", dry_run=not bool(args.apply))
    rows = manager.storage.list_roots(kind=getattr(args, "kind", None)) or []
    root_id = getattr(args, "root_id", None)
    if root_id is not None:
        rows = [row for row in rows if int(row[0]) == int(root_id)]
    if getattr(args, "inactive_only", False):
        rows = [row for row in rows if not bool(row[5])]
    if not rows:
        result.notes.append("No storage roots matched.")
        return result

    selected = []
    blocked_active: list[str] = []
    blocked_bound: list[str] = []
    selected_ids = {int(row[0]) for row in rows}
    suspension_kinds: set[str] = set()
    for row in rows:
        current_root_id = int(row[0])
        kind = str(row[1])
        label = str(row[2] or "")
        base_path = str(row[3] or "")
        writable = bool(row[4])
        active = bool(row[5])
        counts = manager.storage.count_bound_entities(current_root_id)
        selected.append((current_root_id, kind, label, base_path, writable, active, counts))
        result.notes.append(
            f"root {current_root_id}: kind={kind} active={'yes' if active else 'no'} "
            f"writable={'yes' if writable else 'no'} base={base_path} "
            f"bound(genomes={counts['genomes']}, libraries={counts['libraries']}, artifact_rows={counts['artifact_rows']})"
        )
        if active and not bool(getattr(args, "force_active", False)):
            blocked_active.append(str(current_root_id))
        if any(int(counts[name]) > 0 for name in ("genomes", "libraries", "artifact_rows")):
            blocked_bound.append(str(current_root_id))
    for kind in sorted({str(row[1]) for row in rows if str(row[1]) in STRICT_WORKING_ROOT_KINDS}):
        selected_kind_rows = [row for row in rows if str(row[1]) == kind]
        if not any(bool(row[5]) for row in selected_kind_rows):
            continue
        remaining_rows = [
            row
            for row in (manager.storage.list_roots(kind=kind) or [])
            if int(row[0]) not in selected_ids and bool(row[5])
        ]
        if not remaining_rows:
            suspension_kinds.add(kind)

    result.counts = [
        ("StorageRoots", len(selected)),
        ("Active roots selected", sum(1 for row in selected if row[5])),
        ("Bound genomes", sum(int(row[6]["genomes"]) for row in selected)),
        ("Bound libraries", sum(int(row[6]["libraries"]) for row in selected)),
        ("Bound artifact rows", sum(int(row[6]["artifact_rows"]) for row in selected)),
    ]
    result.counts = [(name, count) for name, count in result.counts if count > 0]

    if blocked_active:
        result.notes.append(
            "Active roots are protected by default: " + ", ".join(blocked_active) + ". Re-run with --force-active to purge them."
        )
    if blocked_bound:
        result.notes.append(
            "Bound roots cannot be purged: " + ", ".join(blocked_bound) + ". Move or purge bound data first."
        )
    for kind in sorted(suspension_kinds):
        result.notes.append(
            f"No active {kind} root would remain after purge. Program operations that create new {kind} data will be suspended until a root is activated."
        )

    if result.dry_run:
        return result

    if blocked_active:
        raise ValueError("Refusing to purge active root(s) without --force-active.")
    if blocked_bound:
        raise ValueError("Refusing to purge root(s) that still have bound genomes, libraries, or artifact rows.")

    deleted = 0
    manager.conn.execute("BEGIN IMMEDIATE")
    try:
        for current_root_id, _kind, _label, _base_path, _writable, _active, _counts in selected:
            if manager.storage.delete_root(int(current_root_id)):
                deleted += 1
        manager.conn.commit()
    except Exception:  # boundary: purge transaction must rollback and re-raise original failure.
        manager.conn.rollback()
        raise
    result.notes.append(f"Purged storage roots: {deleted}")
    return result


def _add_purge_common(parser: argparse.ArgumentParser) -> None:
    """Register arguments shared by all purge subjects."""

    parser.add_argument("--apply", action="store_true", help="Execute deletions. Default is dry-run.")
    parser.add_argument(
        "--delete-files",
        action="store_true",
        help="Also delete files/dirs referenced by purged records (restricted to known data roots).",
    )
    parser.add_argument("--all", action="store_true", help="Apply subject purge across all matching rows.")
    parser.add_argument(
        "--force-running",
        action="store_true",
        help="Allow purge while tasks are running.",
    )


def _add_purge_selector_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the reduced selector surface used by purge subjects."""

    parser.add_argument("-a", "--accessions", "--accession", action=AppendCommaSeparated, nargs="+", help="Target accessions.")
    parser.add_argument("-c", "--clade", help="Resolve a scientific name to a taxid.")
    parser.add_argument("-i", "--taxid", type=int, help="Target taxid.")
    parser.add_argument("-d", "--downloaded-only", action="store_true", help="Only include downloaded assemblies.")
    parser.add_argument("--not-downloaded", action="store_true", help="Only include assemblies not downloaded.")
    parser.add_argument("-af", "--after", type=lambda value: _validate_date(value, "--after"), help="Release date lower bound (YYYY-MM-DD).")
    parser.add_argument("-bf", "--before", type=lambda value: _validate_date(value, "--before"), help="Release date upper bound (YYYY-MM-DD).")
    parser.add_argument("--level", choices=["complete genome", "chromosome", "scaffold", "contig"], help="Assembly level filter.")
    parser.add_argument("--primary-only", action="store_true", help="Only include primary assemblies.")
    parser.add_argument("--protein-only", action="store_true", help="Only include protein assemblies.")
    parser.add_argument("--status-min", type=int, help="Minimum Genome.status selector.")


def register_purge_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    """Register the top-level ``purge`` command and its subject parsers."""

    purge_parser = subparsers.add_parser("purge", help="Purge selected data (dry-run by default).")
    purge_sub = purge_parser.add_subparsers(dest="purge_subject", required=True)

    assemblies = purge_sub.add_parser("assemblies", help="Purge assembly-linked records.")
    _add_purge_selector_arguments(assemblies)
    assemblies.add_argument(
        "--include-metadata",
        action="store_true",
        help="Also delete Assembly/Genome metadata rows after purging downloaded/workflow state.",
    )
    _add_purge_common(assemblies)

    decontam = purge_sub.add_parser("decontamination", help="Purge decontamination results.")
    _add_purge_selector_arguments(decontam)
    decontam.add_argument(
        "--run-ids",
        "--run-id",
        dest="run_ids",
        action=AppendCommaSeparated,
        nargs="+",
        help="Filter to one or more decontamination run ids.",
    )
    decontam.add_argument(
        "--busco-run-ids",
        "--busco-run-id",
        dest="busco_run_ids",
        action=AppendCommaSeparated,
        nargs="+",
        help="Filter to one or more linked BUSCO run ids (supports @RUNSET variables).",
    )
    decontam.add_argument(
        "--libraries",
        "--library",
        "--library-name",
        action=AppendCommaSeparated,
        nargs="+",
        help="Target library IDs/names (target_library filter).",
    )
    decontam.add_argument(
        "--delete-runs",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    decontam.add_argument(
        "--keep-runs",
        action="store_true",
        help="Keep Decontamination_Runs envelopes even if they become empty after purge.",
    )
    _add_purge_common(decontam)

    busco = purge_sub.add_parser("busco", help="Purge BUSCO result rows.")
    _add_purge_selector_arguments(busco)
    busco.add_argument(
        "--libraries",
        "--library",
        "--library-name",
        action=AppendCommaSeparated,
        nargs="+",
        help="Library IDs/names to scope BUSCO purge.",
    )
    busco.add_argument(
        "--run-ids",
        "--run-id",
        dest="run_ids",
        action=AppendCommaSeparated,
        nargs="+",
        help="Target BUSCO run ids directly (supports @RUNSET variables).",
    )
    busco.add_argument(
        "--status",
        choices=["running", "completed", "failed", "stale"],
        help="Limit BUSCO purge to runs with the given status.",
    )
    _add_purge_common(busco)

    busco_primary = purge_sub.add_parser("busco-primary", help="Purge manual BUSCO primary overrides.")
    _add_purge_selector_arguments(busco_primary)
    busco_primary.add_argument(
        "--libraries",
        "--library",
        "--library-name",
        action=AppendCommaSeparated,
        nargs="+",
        help="Library IDs/names to scope BUSCO primary purge.",
    )
    _add_purge_common(busco_primary)

    paralog = purge_sub.add_parser("hidden-paralog", help="Purge hidden-paralog (Paralog_Filtering) rows.")
    _add_purge_selector_arguments(paralog)
    paralog.add_argument(
        "--run-ids",
        "--run-id",
        dest="run_ids",
        action=AppendCommaSeparated,
        nargs="+",
        help="Filter to one or more paralog-filtering run ids.",
    )
    paralog.add_argument(
        "--busco-run-ids",
        "--busco-run-id",
        dest="busco_run_ids",
        action=AppendCommaSeparated,
        nargs="+",
        help="Filter to one or more linked BUSCO run ids (supports @RUNSET variables).",
    )
    paralog.add_argument(
        "--libraries",
        "--library",
        "--library-name",
        action=AppendCommaSeparated,
        nargs="+",
        help="Target library IDs/names.",
    )
    paralog.add_argument(
        "--busco-libraries",
        "--busco-library",
        action=AppendCommaSeparated,
        nargs="+",
        help="BUSCO library IDs/names.",
    )
    paralog.add_argument(
        "--delete-runs",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    paralog.add_argument(
        "--keep-runs",
        action="store_true",
        help="Keep Paralog_Filtering_Runs envelopes even if they become empty after purge.",
    )
    _add_purge_common(paralog)

    libraries = purge_sub.add_parser("libraries", help="Purge data scoped to one or more libraries.")
    libraries.add_argument(
        "--libraries",
        "--library",
        "--library-name",
        action=AppendCommaSeparated,
        nargs="+",
        help="Library IDs/names to purge.",
    )
    libraries.add_argument(
        "--drop-library",
        action="store_true",
        help="Delete the matched library definitions themselves; without this flag only dependent data is purged.",
    )
    libraries.add_argument(
        "--recursive",
        action="store_true",
        help="Include child libraries recursively (recommended with --drop-library).",
    )
    libraries.add_argument(
        "--include-core",
        action="store_true",
        help="When --all is used, include root/core libraries (default all custom libraries only).",
    )
    _add_purge_common(libraries)

    variables = purge_sub.add_parser("variables", help="Purge environment variables.")
    variables.add_argument(
        "--custom-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Purge custom variables only (default: true).",
    )
    variables.add_argument(
        "--var-names",
        "--var-name",
        action=AppendCommaSeparated,
        nargs="+",
        help="Explicit variable names to purge.",
    )
    variables.add_argument("--prefix", help="Only purge variables with this prefix.")
    _add_purge_common(variables)

    roots = purge_sub.add_parser("roots", help="Purge storage root definitions.")
    roots.add_argument("--root-id", type=int, help="Purge one specific storage root id.")
    roots.add_argument("--kind", help="Restrict to one storage root kind.")
    roots.add_argument("--inactive-only", action="store_true", help="Restrict purge to inactive roots.")
    roots.add_argument("--force-active", action="store_true", help="Allow purging active roots.")
    _add_purge_common(roots)

    purge_parser.set_defaults(handler=handle_purge)


def _render_result(result: PurgeResult) -> None:
    """Print a dry-run or apply summary for a purge operation."""

    mode = "DRY-RUN" if result.dry_run else "APPLY"
    print(f"Purge {mode}: {result.subject}")
    if result.counts:
        for table, count in result.counts:
            print(f"  {table}: {count}")
    else:
        if result.subject == "libraries" and any(
            note.startswith("Matched libraries count=") for note in result.notes
        ):
            print("  No dependent rows matched.")
        else:
            print("  No matching rows.")
    if result.files:
        print(f"  Files matched: {len(result.files)}")
        preview = result.files[:5]
        for path in preview:
            print(f"    {path}")
        if len(result.files) > 5:
            print(f"    ... {len(result.files) - 5} more")
    if not result.dry_run and (result.deleted_files or result.skipped_files or result.failed_files):
        print(
            f"  Files deleted={result.deleted_files} skipped={result.skipped_files} failed={result.failed_files}"
        )
    for note in result.notes:
        print(f"  Note: {note}")


def handle_purge(args: argparse.Namespace) -> int:
    """Execute the selected purge subject immediately against the database."""

    manager = DBManager(args.database)
    manager.connect()
    try:
        if args.apply and not args.force_running:
            running = _running_task_count(manager)
            if running > 0:
                print(
                    f"Error: {running} task(s) are currently running. Re-run with --force-running to purge anyway."
                )
                return 1

        subject = args.purge_subject
        if subject == "assemblies":
            result = _plan_and_execute_assemblies(manager, args)
        elif subject == "decontamination":
            result = _plan_and_execute_decontamination(manager, args)
        elif subject == "busco":
            result = _plan_and_execute_busco(manager, args)
        elif subject == "busco-primary":
            result = _plan_and_execute_busco_primary(manager, args)
        elif subject == "hidden-paralog":
            result = _plan_and_execute_hidden_paralog(manager, args)
        elif subject == "libraries":
            result = _plan_and_execute_libraries(manager, args)
        elif subject == "variables":
            result = _plan_and_execute_variables(manager, args)
        elif subject == "roots":
            result = _plan_and_execute_roots(manager, args)
        else:  # pragma: no cover - guarded by argparse
            print(f"Error: unknown purge subject '{subject}'.")
            return 1

        _render_result(result)
        return 0
    except ValueError as exc:
        print(f"Error: {exc}")
        return 1
    finally:
        manager.close()


__all__ = ["handle_purge", "register_purge_parser"]
