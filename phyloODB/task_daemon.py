from __future__ import annotations

import json
import sqlite3
import threading
import time
import traceback
from datetime import datetime, timezone
from typing import Any, Optional, Set

from .tasks.task import Task
from .logging_utils import configure_logging_from_db, get_task_logger
from .registry import registry
from .thread_defaults import refresh_runtime_thread_defaults, resolve_task_required_threads

class TaskDaemon(Task):
    def __init__(
        self,
        db_path,
        data,
        max_threads: int | None = None,
        root_task_id: int | None = None,
        scope_mode: str | None = None,
        *,
        polling_time: float | None = None,
        env_overrides: dict | None = None,
        log_overrides: dict | None = None,
        stop_on_error: bool = False,
        stop_after: float | None = None,
    ):
        super().__init__(db_path, task_id=None, data=data, required_threads=1)
        self.default_log_category = "SCHEDULER"
        self.task_display_name = "Scheduler"
        self.root_task_id = int(root_task_id) if root_task_id is not None else None
        self.scope_mode = str(scope_mode or "").strip().lower() or None
        self.running_threads = []
        self.running_task_ids = set()  # Track running task IDs
        self.stop_on_error = bool(stop_on_error)
        self.stop_after = float(stop_after) if stop_after is not None else None
        self.draining = False
        self._stop_deadline: Optional[float] = None
        self._known_error_ids: Set[int] = set()
        # Add self to the Tasks table in the database
        self.db_manager.connect()
        # Configure process logging once from DB env vars (safe no-op if already configured)
        configure_logging_from_db(self.db_manager, overrides=log_overrides, force=True)
        self.task_id = self.db_manager.tasks.queue(
            job_type=0,  # Daemon job type
            status="P",  # Pending
            priority=1,  # Low priority for daemon
            parent_id=None,  # No parent
            data={"daemon": True}  # Indicate this is a daemon task
        )
        # Set daemon logger context now that we have an ID
        self.logger = get_task_logger(task_id=self.task_id, task_type="daemon", task_name=self.task_display_name)
        self.stop_event = threading.Event()
        env_overrides = env_overrides or {}
        runtime = refresh_runtime_thread_defaults(
            self.db_manager,
            explicit_max_threads=max_threads,
            logger=self.log,
        )
        self.max_threads = runtime.max_threads
        env_vars = self.db_manager.env.get_many(
            ['DAEMON_PROCESS_POLLING_TIME', 'BLOCKED_TASK_QUEUE_POLLING_TIME']
        )
        combined = dict(env_vars)
        combined.update({k: v for k, v in env_overrides.items() if v is not None})
        self.log(f"Daemon thread limit set to {self.max_threads}.", level="INFO")
        poll_value = polling_time if polling_time is not None else combined.get('DAEMON_PROCESS_POLLING_TIME', 2)
        self.daemon_process_polling_time = max(0.05, self._coerce_float(poll_value, default=2.0))
        self.log(f"Polling interval set to {self.daemon_process_polling_time}s", level="DEBUG")
        blocked_poll_value = combined.get('BLOCKED_TASK_QUEUE_POLLING_TIME', self.daemon_process_polling_time)
        self.blocked_task_queue_polling_time = max(
            0.05,
            self._coerce_float(blocked_poll_value, default=self.daemon_process_polling_time),
        )
        self.log(f"Blocked-task polling interval set to {self.blocked_task_queue_polling_time}s", level="DEBUG")
        self._next_blocked_poll = 0.0

        # Startup recovery and initial rule enforcement
        try:
            # First, adopt orphaned tasks from previous daemons
            self._adopt_orphaned_tasks()
            # Then run recovery on the unified tree
            self._startup_recovery()
        except Exception as e:  # boundary: startup recovery must not prevent daemon startup
            self.error(f"Recovery step failed: {e}")
        finally:
            # Close initial connection before worker thread starts; runtime thread will reopen as needed
            try:
                self.db_manager.close()
            except Exception as exc:  # boundary: closing initial daemon connection must not hide startup outcome
                self.error(f"Failed to close startup database connection: {exc}")

    def get_overseer_task_id(self):
        """Get the task ID of the overseer task"""
        return self.task_id

    def run(self):
        scope_message = ""
        if self.scope_mode == "task-chain" and self.root_task_id is not None:
            scope_message = f" Scope=task-chain root={self.root_task_id}."
        self.log(f"TaskDaemon started. Monitoring task queue...{scope_message}")

        try:
            self.db_manager.connect()
        except Exception as exc:  # boundary: daemon run startup reports connection failure as task result
            self.error(f"Failed to connect to database: {exc}")
            return False
        self._known_error_ids = self._collect_error_ids()
        if self.stop_after is not None:
            self._stop_deadline = time.time() + self.stop_after

        while not self.stop_event.is_set():
            self._check_stop_conditions()
            if self.stop_event.is_set():
                break
            now_monotonic = time.monotonic()
            # Clean up finished threads
            self.running_threads = [t for t in self.running_threads if t.is_alive()]
            self.running_task_ids = {t.task_id for t in self.running_threads if hasattr(t, "task_id")}

            # Calculate currently used threads
            used_threads = sum(getattr(t, "required_threads", 1) for t in self.running_threads)
            avail_threads = max(self.max_threads - used_threads, 0)

            # Enforce blocking/orphan rules each iteration
            if now_monotonic >= self._next_blocked_poll:
                self._call_with_db_retry(
                    self._enforce_blocking_rules,
                    context="Failed enforcing blocking rules",
                )
                self._next_blocked_poll = now_monotonic + self.blocked_task_queue_polling_time

            # Drain: stop scheduling new tasks and exit when none running
            if self.draining and not self.running_threads:
                self.log("Draining complete; exiting daemon.", level="INFO")
                break

            if self._scope_is_complete():
                self.log("Scoped run complete; exiting daemon.", level="INFO")
                break

            # Check if we can run more tasks
            if avail_threads <= 0:
                # self.log(f"Maximum threads in use ({self.max_threads}); waiting...", level="DEBUG")
                self.stop_event.wait(min(self.daemon_process_polling_time, self.blocked_task_queue_polling_time))
                continue

            # Fetch all runnable tasks sorted by priority
            tasks = self._call_with_db_retry(
                self._load_runnable_tasks,
                context="Failed loading runnable tasks",
                default=[],
            )  # status in ('P','S') ordered by priority, task_id
            if not tasks:
                self.stop_event.wait(min(self.daemon_process_polling_time, self.blocked_task_queue_polling_time))
                continue

            # Try to start tasks where resources allow
            for task_data in tasks:
                if self.draining:
                    break
                task_id = task_data[0]
                job_type = task_data[1]
                task_status = task_data[2]
                required_threads = self._extract_required_threads_override(task_data)
                if required_threads is None:
                    required_threads = self._resolve_required_threads(job_type)
                if required_threads is None:
                    self.error(f"Task {task_id} has unknown job_type {job_type}; marking errored.")
                    try:
                        self.db_manager.tasks.set_error(task_id, f"Unknown job_type {job_type}", "")
                    except Exception as exc:  # boundary: isolate one malformed queue row
                        self._log_scheduler_update_failure(task_id, "record unknown job type", exc)
                    continue
                if required_threads <= 0:
                    required_threads = 1
                required_threads = min(required_threads, int(self.max_threads or 1))
                # self.log(f"Considering task {task_id} (type {job_type}, status {task_status})", level="DEBUG")
                # Never run daemon tasks (job_type 0) from within the daemon
                if job_type == 0:
                    continue

                # Skip tasks that are already running
                if task_id in self.running_task_ids:
                    continue

                # If task is suspended, only resume when its subtasks have finished (or errored)
                if task_status == "S":
                    subtasks = self.db_manager.tasks.get_subtasks(task_id)
                    if subtasks:
                        sub_statuses = [t[2] for t in subtasks]
                        all_complete = all(s == "C" for s in sub_statuses)
                        has_error = any(s in ("E", "B") for s in sub_statuses)
                        # A suspended child is still unfinished: it may be
                        # waiting for its own descendants.  Resuming the parent
                        # while that subtree is active lets later phases race
                        # ahead of required nested work.
                        has_active = any(s in ("P", "R", "S") for s in sub_statuses)
                        if not (all_complete or (has_error and not has_active)):
                            # Still pending subtasks or retry in progress; skip for now
                            continue

                # Check if there are enough threads available
                if required_threads > max(self.max_threads - used_threads, 0):
                    continue

                # Instantiate and start the task
                try:
                    task = self._instantiate_task(task_data, required_threads=required_threads)
                except Exception as e:  # boundary: isolate one task instantiation failure
                    self.error(f"Failed to instantiate task {task_id}: {e}")
                    try:
                        self.db_manager.tasks.set_error(task_id, f"Failed to instantiate task: {e}", traceback.format_exc())
                        self.db_manager.tasks.update_status(task_id, "E")
                    except Exception as exc:  # boundary: isolate failure while recording instantiation failure
                        self._log_scheduler_update_failure(task_id, "record task instantiation failure", exc)
                    continue

                # Set logger context for the child task
                task.logger = get_task_logger(
                    task_id=task_id,
                    task_type=str(job_type),
                    task_name=getattr(task, "task_display_name", None),
                )
                thread = threading.Thread(target=task.start, name=f"task-{task_id}")
                thread.task_id = task_id  # Attach task ID for tracking
                thread.required_threads = required_threads  # Attach info for accounting
                if task_status == "S":
                    self.log(f"Resuming suspended task {task_id}", level="WARNING")
                thread.start()
                self.running_threads.append(thread)
                self.running_task_ids.add(task_id)
                used_threads += required_threads
                avail_threads = max(self.max_threads - used_threads, 0)
                if avail_threads <= 0:
                    break
                # Quiet log; task started

            # Wait for polling time before reassessing
            self.stop_event.wait(min(self.daemon_process_polling_time, self.blocked_task_queue_polling_time))
            self._check_stop_conditions()


        self.log("TaskDaemon stopped.")
        return True

    def _reconnect_db(self) -> None:
        try:
            self.db_manager.close()
        except Exception as exc:  # boundary: reconnect cleanup should not hide the retry attempt
            self.error(f"Failed to close database before reconnect: {exc}")
        self.db_manager.connect()

    def _call_with_db_retry(self, func, *args, context: str, default=None, retries: int = 1, **kwargs):
        attempts = max(0, int(retries)) + 1
        for attempt in range(attempts):
            try:
                return func(*args, **kwargs)
            except sqlite3.OperationalError as exc:
                if "database is locked" not in str(exc).lower():
                    raise
                if attempt >= attempts - 1:
                    self.error(f"{context}: {exc}")
                    return default
                self.log(f"{context}: {exc}; reconnecting and retrying.", level="WARNING")
                try:
                    self._reconnect_db()
                except Exception as reconnect_exc:  # boundary: database-lock retry path
                    self.error(f"{context}: failed to reconnect after lock: {reconnect_exc}")
                    return default
                time.sleep(0.1)
            except Exception as exc:  # boundary: scheduler loop helper isolates one polling operation
                self.error(f"{context}: {exc}")
                return default
        return default

    def _log_scheduler_update_failure(self, task_id: int | None, operation: str, exc: BaseException) -> None:
        target = f"task {task_id}" if task_id is not None else "scheduler state"
        self.error(f"Failed to {operation} for {target}: {exc}")

    def stop(self):
        self.log("Stopping TaskDaemon...")
        self.stop_event.set()

    def drain(self):
        """Stop scheduling new tasks and exit once running tasks finish."""
        self.log("Drain requested: will not schedule new tasks and will exit after running tasks complete.", level="WARNING")
        self.draining = True

    def get_active_task_count(self) -> int:
        """Get the number of currently active (running) tasks."""
        return len(self.running_threads)

    # ===== Control APIs =====
    def block_job(self, job_id: int):
        """Block/pause a job and all its descendants. Safe to call from outer layer."""
        try:
            self.db_manager.tasks.update_status(job_id, "B")
            self.log(f"Blocked job {job_id}", level="WARNING")
            self._block_descendants(job_id)
        except Exception as e:  # boundary: control API must report failure without crashing caller
            self.error(f"Failed to block job {job_id}: {e}")

    def cancel_job(self, job_id: int, reason: str | None = None):
        """Cancel a job and all its descendants by marking them errored with a message."""
        summary = reason or "Cancelled by user"
        try:
            self.db_manager.tasks.set_error(job_id, summary, "")
            self.db_manager.tasks.update_status(job_id, "E")
            self.log(f"Cancelled job {job_id}", level="WARNING")
            # Propagate cancel to descendants
            subtasks = self.db_manager.tasks.get_subtasks(job_id) or []
            for t in subtasks:
                tid = t[0]
                status = t[2]
                if status not in ("C", "E"):
                    try:
                        self.db_manager.tasks.set_error(tid, f"Parent {job_id} cancelled", "")
                    except Exception as exc:  # boundary: continue cancelling descendants while logging one failed error write
                        self._log_scheduler_update_failure(tid, "record cancellation reason", exc)
                    self.db_manager.tasks.update_status(tid, "E")
                    self.log(f"Cancelled child job {tid} due to parent {job_id}", level="WARNING")
                # Recurse
                self.cancel_job(tid, reason=f"Ancestor {job_id} cancelled")
        except Exception as e:  # boundary: control API must report failure without crashing caller
            self.error(f"Failed to cancel job {job_id}: {e}")

    def check_job_status(self, job_id: int) -> str | None:
        """Check the status of a job by its ID."""
        try:
            task = self.db_manager.tasks.get(job_id)
            if task:
                return task[2]  # status field
            else:
                self.error(f"Job {job_id} not found in database.")
                return None
        except Exception as e:  # boundary: control API status check returns None on failure
            self.error(f"Failed to check status of job {job_id}: {e}")
            return None

    # ===== Internals =====
    def _block_descendants(self, job_id: int):
        """Set status B for all descendants (recursive), non-terminal."""
        try:
            subtasks = self.db_manager.tasks.get_subtasks(job_id) or []
        except Exception as exc:  # boundary: recursive block control isolates lookup failure
            self._log_scheduler_update_failure(job_id, "load descendants for blocking", exc)
            subtasks = []
        for t in subtasks:
            tid = t[0]
            status = t[2]
            if status not in ("C", "E"):
                try:
                    self.db_manager.tasks.update_status(tid, "B")
                    self.log(f"Blocked child job {tid} due to parent {job_id}", level="WARNING")
                except Exception as exc:  # boundary: continue blocking other descendants
                    self._log_scheduler_update_failure(tid, "block descendant", exc)
            # Recurse
            self._block_descendants(tid)

    def _startup_recovery(self):
        """Normalize states after crash/restart and block inconsistent children.

        - Flip all Running (R) to Suspended (S), so they can be resumed.
        - Block orphaned subtasks (no parent row).
        - Propagate parent-invalid states: if parent is B/E/C, block children.
        - Emit a recovery summary to the default log.
        """
        tasks = self.db_manager.tasks.get_many() or []
        id_to_task = {t[0]: t for t in tasks}
        flipped = []
        orphans = []
        blocked = []

        # 1) Flip R->S
        for t in tasks:
            tid, status = t[0], t[2]
            if status == "R":
                try:
                    self.db_manager.tasks.update_status(tid, "S")
                    flipped.append(tid)
                except Exception as exc:  # boundary: startup recovery continues with other tasks
                    self._log_scheduler_update_failure(tid, "recover running task to suspended", exc)

        # 2) Orphans & parent-invalid blocking
        for t in tasks:
            tid = t[0]
            status = t[2]
            parent_id = t[4]
            if status in ("C", "E"):
                continue
            if parent_id:
                parent = id_to_task.get(parent_id)
                if parent is None:
                    # Orphan: block
                    try:
                        self.db_manager.tasks.update_status(tid, "B")
                        orphans.append(tid)
                    except Exception as exc:  # boundary: startup recovery continues with other tasks
                        self._log_scheduler_update_failure(tid, "block orphaned task during startup recovery", exc)
                else:
                    p_status = parent[2]
                    if p_status in ("B", "E", "C") and status not in ("B", "C", "E"):
                        try:
                            self.db_manager.tasks.update_status(tid, "B")
                            blocked.append(tid)
                        except Exception as exc:  # boundary: startup recovery continues with other tasks
                            self._log_scheduler_update_failure(tid, "block child during startup recovery", exc)

        # Emit per-item warnings
        for tid in flipped:
            self.log(f"Recovered task {tid}: set R→S for resume", level="WARNING")
        for tid in orphans:
            self.log(f"Blocked orphaned task {tid} (no valid parent)", level="WARNING")
        for tid in blocked:
            self.log(f"Blocked task {tid} due to parent state", level="WARNING")

        # Summary
        self.log(
            f"Recovery summary: flipped R→S={len(flipped)}, orphans blocked={len(orphans)}, children blocked (parent B/E/C)={len(blocked)}"
        )

    def _enforce_blocking_rules(self):
        """Ensure parent/child constraints at runtime.

        - Orphaned subtasks are blocked.
        - Parent blocked -> child blocked.
        - Parent E/C -> block pending/suspended children.
        """
        self._release_blocked_tasks()
        tasks = self.db_manager.tasks.get_many() or []
        id_to_task = {t[0]: t for t in tasks}
        for t in tasks:
            tid = t[0]
            status = t[2]
            if status in ("C", "E", "B"):
                continue
            parent_id = t[4]
            if not parent_id:
                continue
            parent = id_to_task.get(parent_id)
            if parent is None:
                # Orphan
                self.db_manager.tasks.update_status(tid, "B")
                self.log(f"Runtime: blocked orphaned task {tid}", level="WARNING")
                continue
            p_status = parent[2]
            if p_status in ("B", "E", "C"):
                self.db_manager.tasks.update_status(tid, "B")
                self.log(f"Runtime: blocked task {tid} due to parent {parent_id} status {p_status}", level="WARNING")

    def _release_blocked_tasks(self) -> None:
        """Evaluate scheduling constraints and release blocked tasks when satisfied."""
        blocked = self.db_manager.tasks.get_by_status(["B"]) or []
        if not blocked:
            return
        task_ids = [t[0] for t in blocked]
        blocks = self.db_manager.tasks.get_blocks(task_ids, unsatisfied_only=False)
        if not blocks:
            return
        block_by_task: dict[int, list[tuple]] = {}
        block_by_id: dict[int, tuple] = {}
        for block in blocks:
            block_by_id[block[0]] = block
            block_by_task.setdefault(block[1], []).append(block)

        deps = self.db_manager.tasks.get_dependencies(task_ids)
        dep_by_task: dict[int, list[tuple]] = {}
        dep_by_block: dict[int, tuple] = {}
        for dep in deps:
            dep_by_task.setdefault(dep[1], []).append(dep)
            if dep[4] is not None:
                dep_by_block[dep[4]] = dep

        times = self.db_manager.tasks.get_time_constraints(task_ids)
        time_by_task: dict[int, list[tuple]] = {}
        time_by_block: dict[int, tuple] = {}
        for row in times:
            time_by_task.setdefault(row[1], []).append(row)
            if row[4] is not None:
                time_by_block[row[4]] = row

        all_tasks = self.db_manager.tasks.get_many() or []
        status_by_id = {t[0]: (t[2] or "").upper() for t in all_tasks}
        now = datetime.now(timezone.utc)

        for task in blocked:
            task_id = task[0]
            is_barrier = bool(task[12]) if len(task) > 12 else False
            task_blocks = block_by_task.get(task_id)
            if not task_blocks:
                continue
            task_failed = False

            # Dependencies
            for dep in dep_by_task.get(task_id, []):
                dep_id, _tid, depends_on_id, required_state, block_id, allow_failed, satisfied_at = dep
                if satisfied_at:
                    continue
                required_state = (required_state or "").lower()
                dep_status = (status_by_id.get(depends_on_id) or "").upper()
                if not dep_status:
                    self._fail_blocked_task(task_id, f"Dependency {depends_on_id} not found.")
                    task_failed = True
                    break

                satisfied = False
                terminal_mismatch = False
                if required_state == "started":
                    satisfied = dep_status in ("R", "S", "C", "E")
                elif required_state == "finished":
                    satisfied = dep_status in ("C", "E")
                elif required_state == "succeeded":
                    if dep_status == "C":
                        satisfied = True
                    elif dep_status == "E":
                        if allow_failed:
                            satisfied = True
                        else:
                            terminal_mismatch = True
                elif required_state == "failed":
                    if dep_status == "E":
                        satisfied = True
                    elif dep_status == "C":
                        terminal_mismatch = True

                if satisfied:
                    if not satisfied_at:
                        try:
                            self.db_manager.tasks.mark_dependency_satisfied(dep_id)
                        except Exception as exc:  # boundary: continue evaluating other blocked tasks
                            self._log_scheduler_update_failure(task_id, "mark dependency satisfied", exc)
                    if block_id:
                        try:
                            self.db_manager.tasks.mark_block_satisfied(block_id)
                        except Exception as exc:  # boundary: continue evaluating other blocked tasks
                            self._log_scheduler_update_failure(task_id, "mark dependency block satisfied", exc)
                else:
                    if terminal_mismatch:
                        block = block_by_id.get(block_id) if block_id else None
                        block_set = block[8] if block else None
                        if not block_set:
                            self._fail_blocked_task(
                                task_id,
                                f"Dependency {depends_on_id} status {dep_status} does not meet requirement '{required_state}'.",
                            )
                            task_failed = True
                            break
                    if block_id:
                        try:
                            self.db_manager.tasks.mark_block_checked(block_id)
                        except Exception as exc:  # boundary: continue evaluating other blocked tasks
                            self._log_scheduler_update_failure(task_id, "mark dependency block checked", exc)

            if task_failed:
                continue

            # Time constraints
            for constraint in time_by_task.get(task_id, []):
                constraint_id, _tid, _mode, not_before, block_id, satisfied_at = constraint
                if satisfied_at:
                    continue
                try:
                    target = self._coerce_not_before(not_before)
                except (TypeError, ValueError):
                    self._fail_blocked_task(task_id, f"Invalid time constraint '{not_before}'.")
                    task_failed = True
                    break
                if now >= target:
                    if not satisfied_at:
                        try:
                            self.db_manager.tasks.mark_time_constraint_satisfied(constraint_id)
                        except Exception as exc:  # boundary: continue evaluating other blocked tasks
                            self._log_scheduler_update_failure(task_id, "mark time constraint satisfied", exc)
                    if block_id:
                        try:
                            self.db_manager.tasks.mark_block_satisfied(block_id)
                        except Exception as exc:  # boundary: continue evaluating other blocked tasks
                            self._log_scheduler_update_failure(task_id, "mark time block satisfied", exc)
                else:
                    if block_id:
                        try:
                            self.db_manager.tasks.mark_block_checked(block_id)
                        except Exception as exc:  # boundary: continue evaluating other blocked tasks
                            self._log_scheduler_update_failure(task_id, "mark time block checked", exc)

            if task_failed:
                continue

            # Barrier blocks
            for block in task_blocks:
                block_id, _tid, block_type, condition, _message, _created, block_satisfied, _checked, _set, _group = block
                if (block_type or "") != "barrier":
                    continue
                if block_satisfied:
                    continue
                if self._barrier_satisfied(condition or "", task):
                    try:
                        self.db_manager.tasks.mark_block_satisfied(block_id)
                    except Exception as exc:  # boundary: continue evaluating other blocked tasks
                        self._log_scheduler_update_failure(task_id, "mark barrier block satisfied", exc)
                else:
                    try:
                        self.db_manager.tasks.mark_block_checked(block_id)
                    except Exception as exc:  # boundary: continue evaluating other blocked tasks
                        self._log_scheduler_update_failure(task_id, "mark barrier block checked", exc)

            ready, impossible = self._task_blocks_state(
                task_blocks,
                dep_by_block,
                time_by_block,
                status_by_id,
                now,
            )
            if impossible:
                self._fail_blocked_task(task_id, "No remaining schedule option can be satisfied.")
                continue
            if not ready:
                continue

            if is_barrier:
                try:
                    self.db_manager.tasks.update_status(task_id, "C")
                    self.db_manager.tasks.update_end_time(task_id, datetime.utcnow())
                    self.log(f"Released barrier task {task_id}: all constraints satisfied.", level="INFO")
                except Exception as exc:  # boundary: continue releasing other blocked tasks
                    self._log_scheduler_update_failure(task_id, "complete released barrier task", exc)
            else:
                try:
                    self.db_manager.tasks.update_status(task_id, "P")
                    self.log(f"Released blocked task {task_id}: all constraints satisfied.", level="INFO")
                except Exception as exc:  # boundary: continue releasing other blocked tasks
                    self._log_scheduler_update_failure(task_id, "release blocked task", exc)

    def _task_blocks_state(
        self,
        task_blocks: list[tuple],
        dep_by_block: dict[int, tuple],
        time_by_block: dict[int, tuple],
        status_by_id: dict[int, str],
        now: datetime,
    ) -> tuple[bool, bool]:
        grouped = [b for b in task_blocks if b[8]]
        ungrouped = [b for b in task_blocks if not b[8]]

        def block_state(block: tuple) -> str:
            block_id, _tid, block_type, condition, _message, _created, satisfied_at, _checked, _set, _group = block
            if satisfied_at:
                return "satisfied"
            if block_type == "dependency":
                dep = dep_by_block.get(block_id)
                if not dep:
                    return "impossible"
                _dep_id, _tid2, depends_on_id, required_state, _block_id, allow_failed, _sat = dep
                dep_status = (status_by_id.get(depends_on_id) or "").upper()
                if not dep_status:
                    return "impossible"
                state = (required_state or "").lower()
                if state == "started":
                    return "satisfied" if dep_status in ("R", "S", "C", "E") else "pending"
                if state == "finished":
                    return "satisfied" if dep_status in ("C", "E") else "pending"
                if state == "succeeded":
                    if dep_status == "C":
                        return "satisfied"
                    if dep_status == "E" and allow_failed:
                        return "satisfied"
                    if dep_status in ("C", "E"):
                        return "impossible"
                    return "pending"
                if state == "failed":
                    if dep_status == "E":
                        return "satisfied"
                    if dep_status in ("C", "E"):
                        return "impossible"
                    return "pending"
                return "pending"
            if block_type == "time":
                constraint = time_by_block.get(block_id)
                if not constraint:
                    return "impossible"
                _cid, _tid2, _mode, not_before, _block_id, _sat = constraint
                try:
                    target = self._coerce_not_before(not_before)
                except (TypeError, ValueError):
                    return "impossible"
                return "satisfied" if now >= target else "pending"
            if block_type == "barrier":
                return "satisfied" if self._barrier_satisfied(condition or "", (block[1],)) else "pending"
            return "pending"

        if ungrouped:
            if any(block_state(b) == "impossible" for b in ungrouped):
                return False, True
            if any(block_state(b) != "satisfied" for b in ungrouped):
                return False, False

        if not grouped:
            return True, False

        sets: dict[str, dict[str, list[tuple]]] = {}
        for block in grouped:
            set_key = block[8]
            group_key = block[9] or set_key
            sets.setdefault(set_key, {}).setdefault(group_key, []).append(block)

        for groups in sets.values():
            group_satisfied = False
            group_possible = False
            for group_blocks in groups.values():
                states = [block_state(b) for b in group_blocks]
                if all(state == "satisfied" for state in states):
                    group_satisfied = True
                    group_possible = True
                    break
                if "impossible" not in states:
                    group_possible = True
            if not group_possible:
                return False, True
            if not group_satisfied:
                return False, False

        return True, False

    def _coerce_not_before(self, value: str) -> datetime:
        text = str(value).strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        parsed = datetime.fromisoformat(text)
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=timezone.utc)
        return parsed.astimezone(timezone.utc)

    def _barrier_satisfied(self, condition: str, task: tuple) -> bool:
        if not condition:
            return True
        parts = condition.split(":")
        if parts[0] != "queued-drained":
            return True
        return self._queue_drained_global(task[0])

    def _queue_drained_global(self, exclude_task_id: int) -> bool:
        self.db_manager.cursor.execute(
            """
            SELECT COUNT(*) FROM Tasks
            WHERE status IN ('P','R','S','B')
              AND task_id != ?
              AND job_type != 0
            """,
            (exclude_task_id,),
        )
        return int(self.db_manager.cursor.fetchone()[0]) == 0

    def _fail_blocked_task(self, task_id: int, message: str) -> None:
        try:
            self.db_manager.tasks.set_error(task_id, message, "")
        except Exception as exc:  # boundary: still try to set terminal status after error-message failure
            self._log_scheduler_update_failure(task_id, "set blocked-task error message", exc)
        try:
            self.db_manager.tasks.update_status(task_id, "E")
            self.db_manager.tasks.update_end_time(task_id, datetime.utcnow())
        except Exception as exc:  # boundary: continue scheduler even if one blocked task cannot be failed
            self._log_scheduler_update_failure(task_id, "mark blocked task failed", exc)

    def _adopt_orphaned_tasks(self):
        """Reclaim only previously running tasks from earlier overseers.

        - Identify previous overseer tasks (job_type 0) other than the current one.
        - For each previous overseer, reparent only children that were Running (R) at startup.
          Leave P/S/B children under the previous overseer; they'll be blocked by parent E/C.
        - If a previous overseer is not completed (status != 'C'), mark it errored (E) with a
          message listing reclaimed child IDs. If it's completed, leave as-is.
        """
        tasks = self.db_manager.tasks.get_many() or []
        # Identify overseers
        overseers = [t for t in tasks if t[1] == 0]
        if not overseers:
            return
        current_overseer_id = self.task_id
        previous = [t for t in overseers if t[0] != current_overseer_id]
        if not previous:
            return

        # Build child lists per previous overseer and reclaim running ones
        rescued_by_prev: dict[int, list[int]] = {}
        prev_ids = {p[0] for p in previous}
        for t in tasks:
            tid = t[0]
            status = t[2]
            parent_id = t[4]
            if not parent_id or parent_id == current_overseer_id:
                continue
            if parent_id in prev_ids and status == "R":
                # Reclaim only tasks that were running when we started up
                self.db_manager.tasks.update_parent(tid, current_overseer_id)
                rescued_by_prev.setdefault(parent_id, []).append(tid)
                self.log(f"Adopted running task {tid} from overseer {parent_id} → {current_overseer_id}")

        # Mark previous overseers accordingly and attach error messages
        for p in previous:
            pid = p[0]
            p_status = p[2]
            children = rescued_by_prev.get(pid, [])
            if children:
                msg = (
                    f"Task daemon {pid} superseded; reclaimed {len(children)} running subtask(s): "
                    + ", ".join(map(str, sorted(children)))
                    + f"; new overseer {current_overseer_id}"
                )
                try:
                    self.db_manager.tasks.set_error(pid, msg, "")
                except Exception as exc:  # boundary: continue adopting other overseer state
                    self._log_scheduler_update_failure(pid, "record superseded daemon message", exc)
            # If previous overseer isn't completed, mark as errored
            if p_status != "C":
                try:
                    self.db_manager.tasks.update_status(pid, "E")
                    self.log(f"Marked previous overseer {pid} as errored", level="WARNING")
                except Exception as exc:  # boundary: continue startup after previous-daemon status failure
                    self._log_scheduler_update_failure(pid, "mark previous overseer errored", exc)
            elif children:
                # Completed overseer but still reclaimed running tasks (edge case)
                self.log(
                    f"Previous overseer {pid} completed; reclaimed {len(children)} running task(s) anyway.",
                    level="WARNING",
                )
    def _instantiate_task(self, task_data, required_threads=None):
        task_id = task_data[0]
        job_type = task_data[1]
        checkpoint = task_data[5]
        data_json = task_data[6]

        try:
            spec = registry.get_by_job_type(job_type)
        except KeyError as exc:
            raise ValueError(f"Unknown task type: {job_type}") from exc

        if required_threads is None:
            required_threads = self._resolve_required_threads(job_type)
        if required_threads is None:
            required_threads = spec.daemon.required_threads or 1

        if data_json is None:
            payload = {}
        elif isinstance(data_json, dict):
            payload = dict(data_json)
        elif isinstance(data_json, str):
            try:
                payload = json.loads(data_json)
            except (TypeError, json.JSONDecodeError):
                payload = {}
        else:
            try:
                payload = dict(data_json)
            except (TypeError, ValueError):
                payload = {}

        if not isinstance(payload, dict):
            payload = {}
        payload["required_threads"] = int(required_threads)
        data_json = json.dumps(payload)

        return spec.build_task(
            self.db_manager.get_path(),
            task_id=task_id,
            data_json=data_json,
            required_threads=required_threads,
            checkpoint=checkpoint,
        )

    # ----- helpers -----

    def _coerce_float(self, value: Any, default: float) -> float:
        if value is None:
            return float(default)
        try:
            return float(value)
        except (TypeError, ValueError):
            return float(default)

    def _extract_required_threads_override(self, task_data) -> Optional[int]:
        if not task_data or len(task_data) < 7:
            return None
        data_json = task_data[6]
        if not data_json:
            return None
        payload = None
        if isinstance(data_json, dict):
            payload = data_json
        elif isinstance(data_json, str):
            try:
                payload = json.loads(data_json)
            except (TypeError, json.JSONDecodeError):
                payload = None
        if not isinstance(payload, dict):
            return None
        override = payload.get("required_threads")
        if override is None:
            return None
        try:
            override_int = int(override)
        except (TypeError, ValueError):
            return None
        if override_int <= 0:
            return None
        return min(override_int, int(self.max_threads or 1))

    def _load_runnable_tasks(self):
        if self.scope_mode == "task-chain" and self.root_task_id is not None:
            return self.db_manager.tasks.get_runnable_in_subtree(self.root_task_id)
        return self.db_manager.tasks.get_all_runnable()

    def _scope_is_complete(self) -> bool:
        if self.scope_mode != "task-chain" or self.root_task_id is None:
            return False
        try:
            root = self.db_manager.tasks.get(self.root_task_id)
        except Exception as exc:  # boundary: scoped daemon run should keep polling after transient status read failure
            self.error(f"Failed to check scoped root task {self.root_task_id}: {exc}")
            return False
        if not root:
            return True
        root_status = (root[2] or "").upper()
        if root_status not in {"C", "E"}:
            return False
        active = self.db_manager.tasks.count_active_in_subtree(self.root_task_id, statuses=["P", "R", "S"])
        return active == 0

    def _resolve_required_threads(self, job_type: int) -> Optional[int]:
        try:
            spec = registry.get_by_job_type(job_type)
        except KeyError:
            return None
        return resolve_task_required_threads(
            self.db_manager,
            spec,
            self.max_threads,
            logger=self.log,
        )

    def _collect_error_ids(self) -> Set[int]:
        try:
            tasks = self.db_manager.tasks.get_many() or []
        except Exception as exc:  # boundary: stop-on-error check must not kill daemon polling
            self.error(f"Failed to collect task error ids: {exc}")
            return set()
        return {int(t[0]) for t in tasks if len(t) > 2 and t[2] == "E"}

    def _check_stop_conditions(self) -> None:
        if self.stop_event.is_set():
            return
        if self._stop_deadline is not None and time.time() >= self._stop_deadline:
            self.log("Timeout reached; requesting daemon stop.", level="WARNING")
            self.stop_event.set()
            return
        if not self.stop_on_error:
            return
        current_errors = self._collect_error_ids()
        new_errors = sorted(current_errors - self._known_error_ids)
        if new_errors:
            self.log(
                f"Detected new errored task(s) {new_errors}; stop-on-error engaged.",
                level="ERROR",
            )
            self.stop_event.set()
        self._known_error_ids = current_errors
    
