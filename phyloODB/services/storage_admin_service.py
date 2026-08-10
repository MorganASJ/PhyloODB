from __future__ import annotations

import os
import shutil
import json
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

from ..database import DBManager
from ..db.errors import StorageOperationError
from ..permissions import (
    PermissionPolicy,
    comprehensive_directory_probe,
    policy_from_manager,
    quick_check_database_file,
    resolve_scratch_dir,
    sqlite_wal_probe,
)
from ..selector_utils import SelectorRequest, resolve_selector_accessions
from .task_service import TaskService


@dataclass
class RootRebindPlan:
    root_id: int
    kind: str
    label: Optional[str]
    old_base_path: str
    new_base_path: str
    genome_count: int
    library_count: int
    artifact_count: int
    affected_accessions: list[str]
    affected_library_ids: list[int]


@dataclass
class CacheFlushPlan:
    root_ids: list[int]
    root_paths: list[str]
    artifact_count: int
    blastdb_count: int
    filesystem_entries: int


class StorageAdminService:
    def __init__(self, db_path: str, *, db_manager: Optional[DBManager] = None):
        self.db_path = db_path
        self.db_manager = db_manager or DBManager(db_path)

    def _ensure_connection(self) -> None:
        if not getattr(self.db_manager, "conn", None):
            self.db_manager.connect()

    def _task_service(self) -> TaskService:
        return TaskService(self.db_path, db_manager=self.db_manager)

    def resolve_root_token(self, token: Any, *, kind: Optional[str] = None) -> int:
        self._ensure_connection()
        row = self.db_manager.storage.resolve_root_token(token, kind=kind)
        return int(row[0])

    def _root_row(self, root_id: int):
        row = self.db_manager.storage.get_root(int(root_id))
        if not row:
            raise ValueError(f"Unknown storage root id {root_id}.")
        return row

    @staticmethod
    def _root_kind(row) -> str:
        return str(row[1])

    @staticmethod
    def _root_label(row) -> Optional[str]:
        return str(row[2]) if row[2] is not None else None

    @staticmethod
    def _root_base(row) -> str:
        return os.path.abspath(str(row[3]))

    @staticmethod
    def _path_under_base(path: Optional[str], base_path: str) -> bool:
        if not path:
            return False
        abs_path = os.path.abspath(str(path))
        abs_base = os.path.abspath(str(base_path))
        try:
            return os.path.commonpath([abs_path, abs_base]) == abs_base
        except ValueError:
            return False

    def list_roots(self, *, kind: Optional[str] = None):
        self._ensure_connection()
        return self.db_manager.storage.list_roots(kind=kind)

    def check_storage(self) -> list[dict[str, Any]]:
        """Comprehensively validate SQLite WAL, durable roots, and job scratch."""
        self._ensure_connection()
        policy = policy_from_manager(self.db_manager)
        db_parent = Path(self.db_path).expanduser().resolve(strict=False).parent
        db_result = sqlite_wal_probe(db_parent, policy=policy)
        db_file = quick_check_database_file(self.db_path, policy=policy)
        if db_result.ok and not db_file.ok:
            db_result = db_file
        results = [{
            "kind": "database",
            "source": "database parent",
            "path": db_result.path,
            "ok": db_result.ok,
            "mode": db_result.mode,
            "group": db_result.group,
            "message": db_result.message,
        }]
        for row in self.db_manager.storage.list_roots() or []:
            checked = comprehensive_directory_probe(str(row[3]), policy=policy)
            results.append({
                "kind": str(row[1]),
                "source": str(row[2] or f"root {row[0]}"),
                "path": checked.path,
                "ok": checked.ok,
                "mode": checked.mode,
                "group": checked.group,
                "message": checked.message,
            })
        scratch_value = self.db_manager.get_environment_variable("SCRATCH_DIR")
        checked = comprehensive_directory_probe(resolve_scratch_dir(scratch_value))
        results.append({
            "kind": "scratch",
            "source": "SCRATCH_DIR" if scratch_value else "runtime TMPDIR",
            "path": checked.path,
            "ok": checked.ok,
            "mode": checked.mode,
            "group": checked.group,
            "message": checked.message,
        })
        return results

    def _prepare_root_directory(self, base_path: str) -> str:
        raw_path = str(base_path or "").strip()
        if not raw_path:
            raise StorageOperationError("A storage root base path is required.")
        try:
            path = Path(raw_path).expanduser().resolve(strict=False)
            checked = comprehensive_directory_probe(
                path,
                policy=policy_from_manager(self.db_manager),
                create=True,
            )
            if not checked.ok:
                raise OSError(checked.message)
        except (OSError, RuntimeError) as exc:
            raise StorageOperationError(
                f"Storage root '{raw_path}' could not be created or verified as readable and writable: {exc}"
            ) from exc
        return str(path)

    def add_root(
        self,
        *,
        kind: str,
        base_path: str,
        label: Optional[str] = None,
    ) -> dict[str, Any]:
        self._ensure_connection()
        kind = str(kind)
        abs_base = self._prepare_root_directory(base_path)
        existing = self.db_manager.storage.list_roots(kind=kind) or []
        existing_same_path = self.db_manager.storage.get_root_by_base_path(kind=kind, base_path=abs_base)
        if existing_same_path:
            active = bool(existing_same_path[5])
            writable = bool(existing_same_path[4])
        else:
            active = not existing if kind in self.db_manager.storage.STRICT_ACTIVE_KINDS else True
            writable = active if kind in self.db_manager.storage.STRICT_ACTIVE_KINDS else True
        root_id = self.db_manager.storage.ensure_root(
            kind=kind,
            base_path=abs_base,
            label=label,
            is_active=1 if active else 0,
            writable=1 if writable else 0,
        )
        if root_id is None:
            raise ValueError("Failed to create storage root.")
        if active:
            self.db_manager.storage.activate_root(int(root_id))
        row = self.db_manager.storage.get_root(int(root_id))
        return {
            "root_id": int(root_id),
            "kind": kind,
            "created_inactive": existing_same_path is None and bool(existing) and kind in self.db_manager.storage.STRICT_ACTIVE_KINDS and not active,
            "row": row,
        }

    def rename_root(self, *, root_id: Any, label: str) -> dict[str, Any]:
        self._ensure_connection()
        resolved_root_id = self.resolve_root_token(root_id)
        row = self._root_row(resolved_root_id)
        old_label = self._root_label(row)
        new_label = str(label or "").strip()
        if not new_label:
            raise ValueError("A non-empty storage root label is required.")
        self.db_manager.storage.rename_root(resolved_root_id, new_label)
        return {
            "root_id": resolved_root_id,
            "kind": self._root_kind(row),
            "old_label": old_label,
            "new_label": new_label,
        }

    def plan_flush_cache(self, *, root_id: Optional[int] = None) -> CacheFlushPlan:
        self._ensure_connection()
        if root_id is not None:
            root_row = self.db_manager.storage.resolve_root_token(root_id, kind="cache")
            cache_roots = [root_row]
        else:
            cache_roots = self.db_manager.storage.list_roots(kind="cache") or []
        if not cache_roots:
            return CacheFlushPlan(root_ids=[], root_paths=[], artifact_count=0, blastdb_count=0, filesystem_entries=0)

        root_ids = [int(row[0]) for row in cache_roots]
        root_paths = [self._root_base(row) for row in cache_roots]
        artifact_count = 0
        filesystem_entries = 0
        for row in cache_roots:
            artifact_count += int(self.db_manager.storage.count_artifacts_for_root(int(row[0])))
            base_path = self._root_base(row)
            if os.path.isdir(base_path):
                try:
                    filesystem_entries += len(os.listdir(base_path))
                except OSError:
                    pass

        blastdb_count = 0
        for blastdb_id, _library_id, _accession, location, *_rest in self.db_manager.filtering.get_blast_dbs() or []:
            if any(self._path_under_base(location, base_path) for base_path in root_paths):
                blastdb_count += 1

        return CacheFlushPlan(
            root_ids=root_ids,
            root_paths=root_paths,
            artifact_count=int(artifact_count),
            blastdb_count=int(blastdb_count),
            filesystem_entries=int(filesystem_entries),
        )

    def flush_cache(self, *, root_id: Optional[int] = None) -> dict[str, Any]:
        self._ensure_connection()
        plan = self.plan_flush_cache(root_id=root_id)
        if not plan.root_ids:
            return {"plan": plan, "deleted_artifacts": 0, "deleted_blastdbs": 0, "deleted_paths": 0}

        deleted_paths = 0
        for base_path in plan.root_paths:
            if not os.path.isdir(base_path):
                continue
            for name in os.listdir(base_path):
                full = os.path.join(base_path, name)
                try:
                    if os.path.isdir(full) and not os.path.islink(full):
                        shutil.rmtree(full)
                    else:
                        os.unlink(full)
                    deleted_paths += 1
                except FileNotFoundError:
                    continue

        with self.db_manager.transaction(operation="flush cache database records"):
            artifact_ids: list[int] = []
            if plan.root_ids:
                placeholders = ",".join("?" for _ in plan.root_ids)
                artifact_rows = self.db_manager.cursor.execute(
                    f"SELECT artifact_id FROM Artifacts WHERE storage_root_id IN ({placeholders})",
                    tuple(int(value) for value in plan.root_ids),
                ).fetchall()
                artifact_ids = [int(row[0]) for row in artifact_rows if row and row[0] is not None]
            deleted_artifacts = 0
            if artifact_ids:
                placeholders = ",".join("?" for _ in artifact_ids)
                self.db_manager.cursor.execute(f"DELETE FROM Artifacts WHERE artifact_id IN ({placeholders})", tuple(artifact_ids))
                deleted_artifacts = int(self.db_manager.cursor.rowcount or 0)

            blastdb_ids = []
            for row in self.db_manager.filtering.get_blast_dbs() or []:
                blastdb_id = int(row[0])
                location = row[3] if len(row) > 3 else None
                if any(self._path_under_base(location, base_path) for base_path in plan.root_paths):
                    blastdb_ids.append(blastdb_id)
            deleted_blastdbs = self.db_manager.filtering.delete_proteome_blastdb_ids(blastdb_ids)
        return {
            "plan": plan,
            "deleted_artifacts": int(deleted_artifacts),
            "deleted_blastdbs": int(deleted_blastdbs),
            "deleted_paths": int(deleted_paths),
        }

    def activate_root(self, *, root_id: int) -> dict[str, Any]:
        self._ensure_connection()
        resolved_root_id = self.resolve_root_token(root_id)
        row = self._root_row(int(resolved_root_id))
        self.db_manager.storage.activate_root(int(resolved_root_id))
        return {"root_id": int(resolved_root_id), "kind": self._root_kind(row), "base_path": self._root_base(row)}

    def deactivate_root(self, *, root_id: int) -> dict[str, Any]:
        self._ensure_connection()
        resolved_root_id = self.resolve_root_token(root_id)
        row = self._root_row(int(resolved_root_id))
        kind = self._root_kind(row)
        self.db_manager.storage.deactivate_root(int(resolved_root_id))
        active_remaining = self.db_manager.storage.get_default_root_id(kind)
        return {
            "root_id": int(resolved_root_id),
            "kind": kind,
            "base_path": self._root_base(row),
            "active_remaining": active_remaining,
            "suspended": kind in self.db_manager.storage.STRICT_ACTIVE_KINDS and active_remaining is None,
        }

    def plan_rebind_root(self, *, root_id: int, new_base_path: str) -> RootRebindPlan:
        self._ensure_connection()
        resolved_root_id = self.resolve_root_token(root_id)
        row = self._root_row(int(resolved_root_id))
        new_base = os.path.abspath(str(new_base_path))
        genome_rows = self.db_manager.storage.list_genome_bindings(root_id=int(resolved_root_id))
        library_rows = self.db_manager.storage.list_library_bindings(root_id=int(resolved_root_id))
        artifact_count = self.db_manager.storage.count_artifacts_for_root(int(resolved_root_id))
        return RootRebindPlan(
            root_id=int(resolved_root_id),
            kind=self._root_kind(row),
            label=self._root_label(row),
            old_base_path=self._root_base(row),
            new_base_path=new_base,
            genome_count=len(genome_rows),
            library_count=len(library_rows),
            artifact_count=int(artifact_count),
            affected_accessions=[str(r[0]) for r in genome_rows if r and r[0] is not None],
            affected_library_ids=[int(r[0]) for r in library_rows if r and r[0] is not None],
        )

    def apply_rebind_root(self, *, root_id: int, new_base_path: str, verify: bool = True) -> dict[str, Any]:
        self._ensure_connection()
        plan = self.plan_rebind_root(root_id=root_id, new_base_path=new_base_path)
        if not os.path.isdir(plan.new_base_path):
            raise ValueError(f"Destination base path does not exist: {plan.new_base_path}")
        checked = comprehensive_directory_probe(
            plan.new_base_path,
            policy=policy_from_manager(self.db_manager),
        )
        if not checked.ok:
            raise StorageOperationError(
                f"Destination root '{checked.path}' failed preflight: {checked.message}"
            )
        self.db_manager.storage.rebind_root(int(plan.root_id), plan.new_base_path)
        verify_task_ids = self.queue_verify_for_root_plan(plan, verify=verify)
        return {"plan": plan, "verify_task_ids": verify_task_ids}

    def resolve_genome_accessions(self, request: SelectorRequest) -> list[str]:
        self._ensure_connection()
        selected = resolve_selector_accessions(
            self.db_manager,
            request,
            allow_all=True,
            require_candidates=True,
            use_rule_selection=True,
        )
        return list(dict.fromkeys(str(acc) for acc in selected))

    def resolve_library_ids(
        self,
        *,
        library_id: Optional[int] = None,
        library_name: Optional[str] = None,
        ref_accessions: Optional[Sequence[str]] = None,
        all: bool = False,
    ) -> list[int]:
        self._ensure_connection()
        if all:
            rows = self.db_manager.libraries.get(include_inactive=True) or []
            return [int(row[0]) for row in rows if row and row[0] is not None]
        resolved: list[int] = []
        if library_id is not None or library_name is not None:
            if library_id is not None and library_name is not None:
                by_name = self.db_manager.libraries.get_id(str(library_name), include_inactive=True)
                if by_name is None or int(by_name) != int(library_id):
                    raise ValueError("--library-id and --library-name refer to different libraries.")
                resolved.append(int(library_id))
            elif library_id is not None:
                resolved.append(int(library_id))
            else:
                by_name = self.db_manager.libraries.get_id(str(library_name), include_inactive=True)
                if by_name is None:
                    raise ValueError(f"Unknown library '{library_name}'.")
                resolved.append(int(by_name))
        if ref_accessions:
            rows = self.db_manager.libraries.get_by_reference_accessions(ref_accessions) or []
            ids = [int(row[0]) for row in rows if row and row[0] is not None]
            if resolved:
                resolved = [lib_id for lib_id in resolved if lib_id in ids]
                if not resolved:
                    raise ValueError("Provided library selector and ref-accessions did not resolve to the same library.")
            else:
                resolved = ids
        if not resolved:
            raise ValueError("Provide --all, a library identifier, or --ref-accessions.")
        return list(dict.fromkeys(resolved))

    @staticmethod
    def _binding_relative_path(rel_path: Optional[str], fallback_path: str) -> str:
        if rel_path:
            return str(rel_path)
        return os.path.basename(os.path.abspath(str(fallback_path)))

    @staticmethod
    def _dir_size_bytes(path: str) -> int:
        if not os.path.exists(path):
            return 0
        if os.path.isfile(path):
            return int(os.path.getsize(path))
        total = 0
        for root, _dirs, files in os.walk(path):
            for name in files:
                full = os.path.join(root, name)
                try:
                    if not os.path.islink(full):
                        total += int(os.path.getsize(full))
                except OSError:
                    continue
        return total

    @staticmethod
    def _ensure_space(dest_base: str, required_bytes: int) -> tuple[bool, int]:
        usage = shutil.disk_usage(dest_base)
        return usage.free >= required_bytes, int(usage.free)

    def _make_move_rows(
        self,
        *,
        kind: str,
        target_root_id: int,
        rebind_only: bool,
        entities: Sequence[tuple[Any, ...]],
    ) -> list[dict[str, Any]]:
        target_root = self._root_row(int(target_root_id))
        if self._root_kind(target_root) != str(kind):
            raise ValueError(f"Destination root {target_root_id} is kind '{self._root_kind(target_root)}', expected '{kind}'.")
        target_base = self._root_base(target_root)
        rows: list[dict[str, Any]] = []
        for entity in entities:
            if kind == "genomes":
                accession, source_root_id, rel_path, location = entity
                label = str(accession)
                source_path = self.db_manager.storage.resolve_path(
                    storage_root_id=source_root_id,
                    relative_path=rel_path,
                    fallback_path=location,
                )
                relative = self._binding_relative_path(rel_path, source_path or label)
                dest_path = os.path.abspath(os.path.join(target_base, relative))
                if (
                    source_path
                    and source_root_id is not None
                    and int(source_root_id) == int(target_root_id)
                    and os.path.abspath(str(source_path)) == dest_path
                ):
                    continue
                rows.append(
                    {
                        "accession": label,
                        "source_path": source_path,
                        "destination_path": dest_path,
                        "source_root_id": int(source_root_id) if source_root_id is not None else None,
                        "destination_root_id": int(target_root_id),
                        "action": "rebind-only" if rebind_only else "move-files",
                    }
                )
            else:
                library_id, library_name, source_root_id, rel_path, location = entity
                source_path = self.db_manager.storage.resolve_path(
                    storage_root_id=source_root_id,
                    relative_path=rel_path,
                    fallback_path=location,
                )
                relative = self._binding_relative_path(rel_path, source_path or str(library_name))
                dest_path = os.path.abspath(os.path.join(target_base, relative))
                if (
                    source_path
                    and source_root_id is not None
                    and int(source_root_id) == int(target_root_id)
                    and os.path.abspath(str(source_path)) == dest_path
                ):
                    continue
                rows.append(
                    {
                        "library_id": int(library_id),
                        "library_name": str(library_name),
                        "source_path": source_path,
                        "destination_path": dest_path,
                        "source_root_id": int(source_root_id) if source_root_id is not None else None,
                        "destination_root_id": int(target_root_id),
                        "action": "rebind-only" if rebind_only else "move-files",
                    }
                )
        return rows

    def _validate_move_rows(self, rows: Sequence[dict[str, Any]], *, rebind_only: bool) -> list[str]:
        issues: list[str] = []
        if not rows:
            return ["No matching entities were selected."]
        if rebind_only:
            for row in rows:
                dest = row["destination_path"]
                if not os.path.exists(dest):
                    key = row.get("accession") or row.get("library_name")
                    issues.append(f"{key}: destination missing: {dest}")
            return issues
        total_size = 0
        dest_roots: dict[int, str] = {}
        for row in rows:
            source = row["source_path"]
            dest = row["destination_path"]
            dest_root_id = int(row["destination_root_id"])
            key = row.get("accession") or row.get("library_name")
            if not source or not os.path.exists(source):
                issues.append(f"{key}: source missing: {source}")
                continue
            if os.path.abspath(source) == os.path.abspath(dest):
                continue
            if os.path.exists(dest):
                issues.append(f"{key}: destination already exists: {dest}")
                continue
            total_size += self._dir_size_bytes(source)
            dest_roots.setdefault(dest_root_id, os.path.dirname(dest) or dest)
        seen_root_sizes: set[int] = set()
        for row in rows:
            dest_root_id = int(row["destination_root_id"])
            if dest_root_id in seen_root_sizes:
                continue
            seen_root_sizes.add(dest_root_id)
            root = self._root_row(dest_root_id)
            base = self._root_base(root)
            if not os.path.isdir(base):
                issues.append(f"destination root base missing: {base}")
                continue
            ok, free_bytes = self._ensure_space(base, total_size)
            if not ok:
                issues.append(
                    f"destination root {dest_root_id} lacks free space: need {total_size} bytes, have {free_bytes} bytes"
                )
        return issues

    def plan_move_genomes(
        self,
        *,
        request: SelectorRequest,
        target_root_id: int,
        rebind_only: bool = False,
    ) -> dict[str, Any]:
        resolved_target_root_id = self.resolve_root_token(target_root_id, kind="genomes")
        accessions = self.resolve_genome_accessions(request)
        entities = []
        for acc in accessions:
            row = self.db_manager.cursor.execute(
                "SELECT accession, storage_root_id, relative_path, location FROM Genome WHERE accession = ?",
                (acc,),
            ).fetchone()
            if row:
                entities.append(row)
        rows = self._make_move_rows(kind="genomes", target_root_id=int(resolved_target_root_id), rebind_only=rebind_only, entities=entities)
        issues = self._validate_move_rows(rows, rebind_only=rebind_only)
        return {"rows": rows, "issues": issues, "accessions": accessions}

    def plan_move_libraries(
        self,
        *,
        library_ids: Sequence[int],
        target_root_id: int,
        rebind_only: bool = False,
    ) -> dict[str, Any]:
        resolved_target_root_id = self.resolve_root_token(target_root_id, kind="libraries")
        entities = []
        for library_id in library_ids:
            row = self.db_manager.cursor.execute(
                "SELECT library_id, library_name, storage_root_id, relative_path, location FROM Libraries WHERE library_id = ?",
                (int(library_id),),
            ).fetchone()
            if row:
                entities.append(row)
        rows = self._make_move_rows(kind="libraries", target_root_id=int(resolved_target_root_id), rebind_only=rebind_only, entities=entities)
        issues = self._validate_move_rows(rows, rebind_only=rebind_only)
        return {"rows": rows, "issues": issues, "library_ids": [int(lib) for lib in library_ids]}

    def _apply_row_move(self, row: dict[str, Any], *, kind: str, rebind_only: bool) -> None:
        dest_path = str(row["destination_path"])
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        if not rebind_only:
            source_path = str(row["source_path"])
            if os.path.abspath(source_path) != os.path.abspath(dest_path):
                shutil.move(source_path, dest_path)
        if kind == "genomes":
            self.db_manager.genomes.set_binding(str(row["accession"]), dest_path, kind="genomes")
        else:
            self.db_manager.libraries.set_binding(int(row["library_id"]), dest_path, kind="libraries")

    @staticmethod
    def _copy_to_stage(source_path: str, staging_path: str) -> None:
        if os.path.isdir(source_path):
            shutil.copytree(source_path, staging_path)
        else:
            os.makedirs(os.path.dirname(staging_path), exist_ok=True)
            shutil.copy2(source_path, staging_path)

    @staticmethod
    def _remove_source(path: str) -> None:
        if not os.path.exists(path):
            return
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.unlink(path)

    def _apply_journaled_row_move(self, row: dict[str, Any], *, kind: str) -> int:
        source = os.path.abspath(str(row["source_path"]))
        destination = os.path.abspath(str(row["destination_path"]))
        operation_id = self.db_manager.storage.create_filesystem_operation(
            operation_type=f"move-{kind}",
            source_path=source,
            destination_path=destination,
            payload={"kind": kind, "row": row},
        )
        staging = f"{destination}.phyloodb-stage-{operation_id}"
        self.db_manager.storage.update_filesystem_operation(
            operation_id,
            status="preparing",
            staging_path=staging,
        )
        try:
            self._copy_to_stage(source, staging)
            self.db_manager.storage.update_filesystem_operation(operation_id, status="prepared")
            with self.db_manager.transaction(operation=f"commit filesystem move {operation_id}"):
                if kind == "genomes":
                    self.db_manager.genomes.set_binding(str(row["accession"]), destination, kind="genomes")
                else:
                    self.db_manager.libraries.set_binding(int(row["library_id"]), destination, kind="libraries")
                self.db_manager.storage.update_filesystem_operation(operation_id, status="db_committed")
            os.makedirs(os.path.dirname(destination), exist_ok=True)
            os.replace(staging, destination)
            self._remove_source(source)
            self.db_manager.storage.update_filesystem_operation(operation_id, status="finalized")
            return operation_id
        except Exception as exc:  # boundary: journal every filesystem finalization failure
            current = self.db_manager.storage.get_filesystem_operation(operation_id)
            current_status = str(current[5]) if current else "failed"
            if current_status != "db_committed":
                self._remove_source(staging)
            self.db_manager.storage.update_filesystem_operation(
                operation_id,
                status="failed",
                error_message=str(exc),
            )
            raise StorageOperationError(
                f"Filesystem move {operation_id} failed ({source} -> {destination}): {exc}"
            ) from exc

    def recover_filesystem_operations(
        self,
        *,
        operation_id: Optional[int] = None,
        apply: bool = False,
    ) -> list[dict[str, Any]]:
        self._ensure_connection()
        rows = (
            [self.db_manager.storage.get_filesystem_operation(int(operation_id))]
            if operation_id is not None
            else self.db_manager.storage.list_filesystem_operations(pending_only=True)
        )
        results: list[dict[str, Any]] = []
        for record in [row for row in rows if row]:
            op_id, op_type, source, staging, destination, status, payload_json, error, *_ = record
            action = "inspect"
            recoverable = bool(destination and (os.path.exists(str(staging or "")) or os.path.exists(str(destination))))
            if apply and recoverable:
                try:
                    if staging and os.path.exists(str(staging)):
                        if os.path.exists(str(destination)):
                            raise StorageOperationError(f"Destination already exists: {destination}")
                        os.makedirs(os.path.dirname(str(destination)), exist_ok=True)
                        os.replace(str(staging), str(destination))
                    if source and destination and os.path.exists(str(destination)):
                        self._remove_source(str(source))
                    self.db_manager.storage.update_filesystem_operation(int(op_id), status="finalized")
                    status = "finalized"
                    error = None
                    action = "finalized"
                except Exception as exc:  # boundary: recover remaining journal entries independently
                    self.db_manager.storage.update_filesystem_operation(
                        int(op_id), status="failed", error_message=str(exc)
                    )
                    status = "failed"
                    error = str(exc)
                    action = "failed"
            results.append(
                {
                    "operation_id": int(op_id),
                    "operation_type": str(op_type),
                    "status": str(status),
                    "source_path": source,
                    "staging_path": staging,
                    "destination_path": destination,
                    "recoverable": recoverable,
                    "action": action,
                    "error": error,
                    "payload": json.loads(payload_json or "{}"),
                }
            )
        return results

    def queue_verify_for_accessions(self, accessions: Sequence[str], *, tidy: bool = True, verify: bool = True) -> list[int]:
        if not verify or not accessions:
            return []
        task_id = self._task_service().queue(
            "verify-assembly",
            payload={
                "accessions": list(dict.fromkeys([str(acc) for acc in accessions])),
                "repair": True,
                "tidy": bool(tidy),
            },
        )
        return [int(task_id)]

    def queue_finalize_genome_move(self, plan: dict[str, Any], *, verify: bool = True, tidy: bool = True) -> int:
        rows = list(plan.get("rows") or [])
        if not rows:
            raise ValueError("No genome move rows were planned.")
        return int(
            self._task_service().queue(
                "finalize-genome-move",
                payload={
                    "rows": rows,
                    "verify": bool(verify),
                    "tidy": bool(tidy),
                },
            )
        )

    def queue_verify_for_libraries(self, library_ids: Sequence[int], *, verify: bool = True) -> list[int]:
        if not verify:
            return []
        task_ids: list[int] = []
        service = self._task_service()
        for library_id in list(dict.fromkeys(int(lib) for lib in library_ids)):
            task_ids.append(
                int(
                    service.queue(
                        "verify-libraries",
                        payload={
                            "library_id": int(library_id),
                            "repair": True,
                        },
                    )
                )
            )
        return task_ids

    def queue_verify_for_root_plan(self, plan: RootRebindPlan, *, verify: bool = True) -> list[int]:
        if not verify:
            return []
        task_ids: list[int] = []
        if plan.kind == "genomes":
            task_ids.extend(self.queue_verify_for_accessions(plan.affected_accessions, tidy=True, verify=True))
        elif plan.kind == "libraries":
            task_ids.extend(self.queue_verify_for_libraries(plan.affected_library_ids, verify=True))
        return task_ids

    def apply_move_genomes(self, plan: dict[str, Any], *, rebind_only: bool = False, verify: bool = True) -> dict[str, Any]:
        issues = list(plan.get("issues") or [])
        if issues:
            raise ValueError("Cannot apply move with preflight issues.")
        task_id = self.queue_finalize_genome_move(plan, verify=verify, tidy=True)
        return {"queued_task_id": int(task_id)}

    def apply_move_libraries(self, plan: dict[str, Any], *, rebind_only: bool = False, verify: bool = True) -> dict[str, Any]:
        issues = list(plan.get("issues") or [])
        if issues:
            raise ValueError("Cannot apply move with preflight issues.")
        applied: list[int] = []
        operation_ids: list[int] = []
        for row in plan.get("rows", []):
            if rebind_only:
                with self.db_manager.transaction(operation=f"rebind library {row['library_id']}"):
                    self._apply_row_move(row, kind="libraries", rebind_only=True)
            else:
                operation_ids.append(self._apply_journaled_row_move(row, kind="libraries"))
            applied.append(int(row["library_id"]))
        verify_task_ids = self.queue_verify_for_libraries(applied, verify=verify)
        return {
            "applied_library_ids": applied,
            "verify_task_ids": verify_task_ids,
            "filesystem_operation_ids": operation_ids,
        }
