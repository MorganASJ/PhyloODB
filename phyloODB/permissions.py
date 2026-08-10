"""Project permission policy, scratch resolution, and filesystem preflights."""
from __future__ import annotations

import contextlib
import os
import shutil
import sqlite3
import stat
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

try:
    import grp
except ImportError:  # pragma: no cover - native Windows
    grp = None  # type: ignore[assignment]

from .db.errors import StorageOperationError

PERMISSION_MODE_VAR = "PROJECT_PERMISSION_MODE"
SHARED_GROUP_VAR = "SHARED_GROUP"
SCRATCH_DIR_VAR = "SCRATCH_DIR"
POLICY_VARIABLES = {PERMISSION_MODE_VAR, SHARED_GROUP_VAR, SCRATCH_DIR_VAR, "CACHE_DIR"}

_VERIFIED_PATHS: set[tuple[str, int]] = set()


@dataclass(frozen=True)
class PermissionPolicy:
    mode: str = "private"
    group: str | None = None

    @property
    def shared(self) -> bool:
        return self.mode == "shared"


@dataclass(frozen=True)
class PathCheck:
    path: str
    ok: bool
    message: str = "ok"
    mode: int | None = None
    group: str | None = None


def is_native_windows() -> bool:
    return os.name == "nt"


def resolve_group(group_name: str | None) -> Any:
    name = str(group_name or "").strip()
    if not name:
        raise StorageOperationError("Shared projects require a non-empty SHARED_GROUP/--group value.")
    if is_native_windows():
        raise StorageOperationError("Shared POSIX permissions are unsupported on native Windows; use WSL2.")
    if grp is None:
        raise StorageOperationError("Shared POSIX permissions are unsupported on native Windows; use WSL2.")
    try:
        entry = grp.getgrnam(name)
    except KeyError as exc:
        raise StorageOperationError(f"Unix group '{name}' does not exist.") from exc
    groups = set(os.getgroups()) | {os.getgid(), os.getegid()}
    if entry.gr_gid not in groups:
        raise StorageOperationError(
            f"Current user is not an active member of group '{name}'. Log in again or run 'newgrp {name}'."
        )
    return entry


def policy_from_values(values: Mapping[str, Any]) -> PermissionPolicy:
    mode = str(values.get(PERMISSION_MODE_VAR, "private") or "private").strip().lower()
    if mode not in {"private", "shared"}:
        raise StorageOperationError(f"{PERMISSION_MODE_VAR} must be 'private' or 'shared'.")
    group = values.get(SHARED_GROUP_VAR)
    group = str(group).strip() if group not in (None, "") else None
    if mode == "shared":
        resolve_group(group)
    return PermissionPolicy(mode=mode, group=group)


def policy_from_manager(manager: Any) -> PermissionPolicy:
    values = manager.get_environment_variables([PERMISSION_MODE_VAR, SHARED_GROUP_VAR]) or {}
    return policy_from_values(values)


def apply_shared_umask(policy: PermissionPolicy) -> None:
    if policy.shared:
        os.umask(0o007)


def resolve_scratch_dir(value: Any = None) -> str:
    raw = str(value or "").strip()
    if raw:
        return str(Path(raw).expanduser().resolve(strict=False))
    return str(Path(tempfile.gettempdir()).resolve(strict=False))


def scratch_dir_from_manager(manager: Any) -> str:
    return resolve_scratch_dir(manager.get_environment_variable(SCRATCH_DIR_VAR))


def _group_name(gid: int) -> str:
    if grp is None:
        return str(gid)
    try:
        return grp.getgrgid(gid).gr_name
    except KeyError:
        return str(gid)


def remediation(path: Path, group: str) -> str:
    quoted = str(path).replace("'", "'\\''")
    return f"Ask an administrator to run: chgrp '{group}' '{quoted}' && chmod g+rwx,g+s '{quoted}'"


def quick_check_directory(
    path_value: str | os.PathLike[str],
    *,
    policy: PermissionPolicy = PermissionPolicy(),
    require_setgid: bool = True,
) -> PathCheck:
    path = Path(path_value).expanduser().resolve(strict=False)
    try:
        info = path.stat()
        if not stat.S_ISDIR(info.st_mode):
            raise NotADirectoryError(str(path))
        if not os.access(path, os.R_OK | os.W_OK | os.X_OK):
            raise PermissionError("directory is not readable, writable, and searchable")
        usage = shutil.disk_usage(path)
        if usage.free <= 0:
            raise OSError("filesystem reports no free space")
        group_name = _group_name(info.st_gid)
        mode = stat.S_IMODE(info.st_mode)
        if policy.shared:
            wanted = resolve_group(policy.group)
            if info.st_gid != wanted.gr_gid or not bool(mode & stat.S_IWGRP) or (require_setgid and not bool(mode & stat.S_ISGID)):
                raise PermissionError(
                    f"expected group '{policy.group}', group-write, and setgid; found group '{group_name}' mode {mode:04o}. "
                    + remediation(path, str(policy.group))
                )
        return PathCheck(str(path), True, mode=mode, group=group_name)
    except (OSError, RuntimeError, StorageOperationError) as exc:
        return PathCheck(str(path), False, str(exc))


def quick_check_database_file(path_value: str | os.PathLike[str], *, policy: PermissionPolicy) -> PathCheck:
    path = Path(path_value).expanduser().resolve(strict=False)
    try:
        info = path.stat()
        if not stat.S_ISREG(info.st_mode):
            raise OSError("database path is not a regular file")
        if not os.access(path, os.R_OK | os.W_OK):
            raise PermissionError("database is not readable and writable")
        mode = stat.S_IMODE(info.st_mode)
        group_name = _group_name(info.st_gid)
        if policy.shared:
            wanted = resolve_group(policy.group)
            if info.st_gid != wanted.gr_gid or not bool(mode & stat.S_IWGRP):
                quoted = str(path).replace("'", "'\\''")
                raise PermissionError(
                    f"expected group '{policy.group}' and group-write; found group '{group_name}' mode {mode:04o}. "
                    f"Ask an administrator to run: chgrp '{policy.group}' '{quoted}' && chmod g+rw '{quoted}'"
                )
        return PathCheck(str(path), True, mode=mode, group=group_name)
    except (OSError, StorageOperationError) as exc:
        return PathCheck(str(path), False, str(exc))


def comprehensive_directory_probe(
    path_value: str | os.PathLike[str],
    *,
    policy: PermissionPolicy = PermissionPolicy(),
    create: bool = False,
) -> PathCheck:
    path = Path(path_value).expanduser().resolve(strict=False)
    created = False
    created_paths: list[Path] = []
    try:
        if not path.exists():
            if not create:
                raise FileNotFoundError(str(path))
            cursor = path
            while not cursor.exists():
                created_paths.append(cursor)
                cursor = cursor.parent
            path.mkdir(parents=True, mode=0o770 if policy.shared else 0o777, exist_ok=False)
            created = True
        if policy.shared and created:
            group = resolve_group(policy.group)
            for created_path in reversed(created_paths):
                os.chown(created_path, -1, group.gr_gid)
                os.chmod(created_path, 0o2770)
        checked = quick_check_directory(path, policy=policy)
        if not checked.ok:
            return checked
        fd, raw_probe = tempfile.mkstemp(prefix=".phyloodb-write-test-", dir=path)
        probe = Path(raw_probe)
        renamed = probe.with_name(probe.name + ".renamed")
        try:
            os.write(fd, b"phyloodb\n")
            os.fsync(fd)
            os.close(fd)
            fd = -1
            os.replace(probe, renamed)
            renamed.unlink()
        finally:
            if fd >= 0:
                os.close(fd)
            with contextlib.suppress(OSError):
                probe.unlink()
            with contextlib.suppress(OSError):
                renamed.unlink()
        return checked
    except (OSError, RuntimeError, StorageOperationError) as exc:
        if created:
            for created_path in created_paths:
                with contextlib.suppress(OSError):
                    created_path.rmdir()
        return PathCheck(str(path), False, str(exc))


def sqlite_wal_probe(parent_value: str | os.PathLike[str], *, policy: PermissionPolicy = PermissionPolicy()) -> PathCheck:
    parent = Path(parent_value).expanduser().resolve(strict=False)
    base = comprehensive_directory_probe(parent, policy=policy, create=True)
    if not base.ok:
        return base
    fd, raw_path = tempfile.mkstemp(prefix=".phyloodb-wal-test-", suffix=".db", dir=parent)
    os.close(fd)
    db_path = Path(raw_path)
    conn = None
    try:
        conn = sqlite3.connect(db_path)
        mode = str(conn.execute("PRAGMA journal_mode=WAL").fetchone()[0]).lower()
        if mode != "wal":
            raise OSError(f"SQLite refused WAL mode (reported {mode!r})")
        conn.execute("CREATE TABLE probe(value TEXT)")
        conn.execute("INSERT INTO probe VALUES ('ok')")
        conn.commit()
        if conn.execute("SELECT value FROM probe").fetchone() != ("ok",):
            raise OSError("SQLite WAL write could not be read back")
        return PathCheck(str(parent), True, mode=base.mode, group=base.group)
    except (OSError, sqlite3.Error) as exc:
        return PathCheck(str(parent), False, f"SQLite WAL probe failed: {exc}")
    finally:
        if conn is not None:
            with contextlib.suppress(sqlite3.Error):
                conn.close()
        for suffix in ("", "-wal", "-shm"):
            with contextlib.suppress(OSError):
                Path(str(db_path) + suffix).unlink()


def ensure_quick_path(
    path_value: str | os.PathLike[str],
    *,
    policy: PermissionPolicy,
    probe_if_new: bool = False,
) -> None:
    path = Path(path_value).expanduser().resolve(strict=False)
    checked = quick_check_directory(path, policy=policy)
    if not checked.ok:
        raise StorageOperationError(f"Path preflight failed for '{checked.path}': {checked.message}")
    identity = (str(path), path.stat().st_dev)
    if probe_if_new and identity not in _VERIFIED_PATHS:
        checked = comprehensive_directory_probe(path, policy=policy)
        if not checked.ok:
            raise StorageOperationError(f"Path probe failed for '{checked.path}': {checked.message}")
        _VERIFIED_PATHS.add(identity)


def validate_policy_change(manager: Any, updates: Mapping[str, Any]) -> PermissionPolicy:
    current = manager.get_environment_variables([PERMISSION_MODE_VAR, SHARED_GROUP_VAR]) or {}
    proposed = dict(current)
    proposed.update({key: value for key, value in updates.items() if key in {PERMISSION_MODE_VAR, SHARED_GROUP_VAR}})
    policy = policy_from_values(proposed)
    if SHARED_GROUP_VAR in updates and updates[SHARED_GROUP_VAR] not in (None, ""):
        resolve_group(str(updates[SHARED_GROUP_VAR]))
    if policy.shared:
        paths: list[str] = [str(Path(manager.get_path()).resolve().parent)]
        paths.extend(str(row[3]) for row in (manager.storage.list_roots() or []))
        failures = [result for result in (quick_check_directory(path, policy=policy) for path in dict.fromkeys(paths)) if not result.ok]
        db_check = quick_check_database_file(manager.get_path(), policy=policy)
        if not db_check.ok:
            failures.insert(0, db_check)
        if failures:
            details = "\n".join(f"- {item.path}: {item.message}" for item in failures)
            raise StorageOperationError(f"Existing project paths do not satisfy the proposed shared policy:\n{details}")
    return policy


def validate_config_updates(manager: Any, updates: Mapping[str, Any]) -> dict[str, Any]:
    normalized = dict(updates)
    if "CACHE_DIR" in updates:
        raise StorageOperationError("CACHE_DIR is managed by the registered cache root; use 'storage rebind-root CACHE_DIR --base-path PATH --apply'.")
    if SCRATCH_DIR_VAR in updates and updates[SCRATCH_DIR_VAR] not in (None, ""):
        scratch = Path(str(updates[SCRATCH_DIR_VAR])).expanduser()
        if "\x00" in str(scratch):
            raise StorageOperationError("SCRATCH_DIR contains a NUL byte.")
        checked = comprehensive_directory_probe(scratch, create=True)
        if not checked.ok:
            raise StorageOperationError(f"SCRATCH_DIR preflight failed for '{checked.path}': {checked.message}")
        normalized[SCRATCH_DIR_VAR] = checked.path
    elif SCRATCH_DIR_VAR in updates:
        normalized[SCRATCH_DIR_VAR] = None
    if PERMISSION_MODE_VAR in updates:
        normalized[PERMISSION_MODE_VAR] = str(updates[PERMISSION_MODE_VAR]).strip().lower()
    if SHARED_GROUP_VAR in updates:
        normalized[SHARED_GROUP_VAR] = str(updates[SHARED_GROUP_VAR]).strip() if updates[SHARED_GROUP_VAR] not in (None, "") else None
    if {PERMISSION_MODE_VAR, SHARED_GROUP_VAR} & set(updates):
        validate_policy_change(manager, normalized)
    return normalized


def task_paths(manager: Any, spec: Any, payload: Mapping[str, Any]) -> tuple[list[str], list[str]]:
    metadata = dict(getattr(spec, "metadata", {}) or {})
    roots: list[str] = []
    for kind in (getattr(spec, "write_root_kinds", ()) or metadata.get("write_roots", ())):
        base = manager.storage.get_root_base(str(kind), ensure_from_env=True)
        if base:
            roots.append(str(base))
    explicit: list[str] = []
    for field in (getattr(spec, "output_path_fields", ()) or metadata.get("output_path_fields", ())):
        value = payload.get(str(field))
        if value:
            candidate = Path(str(value)).expanduser().resolve(strict=False)
            explicit.append(str(candidate if candidate.exists() and candidate.is_dir() else candidate.parent))
    return list(dict.fromkeys(roots)), list(dict.fromkeys(explicit))


def preflight_task(manager: Any, spec: Any, payload: Mapping[str, Any]) -> None:
    policy = policy_from_manager(manager)
    apply_shared_umask(policy)
    ensure_quick_path(Path(manager.get_path()).resolve().parent, policy=policy)
    db_check = quick_check_database_file(manager.get_path(), policy=policy)
    if not db_check.ok:
        raise StorageOperationError(f"Database preflight failed for '{db_check.path}': {db_check.message}")
    roots, explicit = task_paths(manager, spec, payload)
    for path in roots:
        ensure_quick_path(path, policy=policy)
    for path in explicit:
        ensure_quick_path(path, policy=policy, probe_if_new=True)
    if bool(getattr(spec, "uses_scratch", False) or metadata.get("uses_scratch")):
        ensure_quick_path(scratch_dir_from_manager(manager), policy=PermissionPolicy(), probe_if_new=True)


@contextlib.contextmanager
def task_scratch_directory(manager: Any, *, prefix: str = "phyloodb-"):
    base = scratch_dir_from_manager(manager)
    ensure_quick_path(base, policy=PermissionPolicy(), probe_if_new=True)
    with tempfile.TemporaryDirectory(prefix=prefix, dir=base) as path:
        yield path
