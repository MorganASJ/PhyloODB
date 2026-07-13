from __future__ import annotations

import json
from datetime import datetime

from .base import BaseRepository, transactional_methods


@transactional_methods(
    "queue",
    "update_status",
    "update_data",
    "set_error",
    "mark_dependency_satisfied",
    "mark_block_satisfied",
    "mark_block_checked",
    "mark_time_constraint_satisfied",
    "update_end_time",
    "update_start_time",
    "update_parent",
    "clear",
    "reset",
    "kill_and_descendants",
    "cancel_and_descendants",
    "add_block",
    "add_dependency",
    "add_time_constraint",
    "get_next_runnable",
)
class TaskRepository(BaseRepository):
    @staticmethod
    def _normalize_dt(value):
        if isinstance(value, datetime):
            return value.strftime("%Y-%m-%d %H:%M:%S")
        if value is None:
            return None
        return str(value).split(".")[0]

    def queue(
        self,
        job_type,
        status,
        priority,
        parent_id,
        checkpoint=0,
        data=None,
        is_barrier=0,
    ):
        self.core.execute(
            """
            INSERT INTO Tasks (
                job_type, status, priority, parent_id, checkpoint, data, is_barrier, status_updated_at
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, datetime('now'))
            """,
            (
                job_type,
                status,
                priority,
                parent_id,
                checkpoint,
                json.dumps(data) if data is not None else None,
                int(is_barrier or 0),
            ),
        )
        self.conn.commit()
        return self.cursor.lastrowid

    def update_status(self, task_id, status):
        self.core.execute(
            "UPDATE Tasks SET status = ?, status_updated_at = datetime('now') WHERE task_id = ?",
            (status, task_id),
        )
        self.conn.commit()
        return True

    def update_data(self, task_id, checkpoint=None, data=None):
        if checkpoint is not None and data is not None:
            self.core.execute(
                "UPDATE Tasks SET checkpoint = ?, data = ? WHERE task_id = ?",
                (checkpoint, json.dumps(data), task_id),
            )
        elif data is not None:
            self.core.execute(
                "UPDATE Tasks SET data = ? WHERE task_id = ?",
                (json.dumps(data), task_id),
            )
        elif checkpoint is not None:
            self.core.execute(
                "UPDATE Tasks SET checkpoint = ? WHERE task_id = ?",
                (checkpoint, task_id),
            )
        else:
            return False
        self.conn.commit()
        return True

    def set_error(self, task_id, message, stack):
        self.core.execute(
            "UPDATE Tasks SET error_message = ?, error_stack = ? WHERE task_id = ?",
            (message, stack, task_id),
        )
        self.conn.commit()
        return True

    def get_status(self, task_id):
        row = self.core.fetchone("SELECT status FROM Tasks WHERE task_id = ?", (task_id,))
        return row[0] if row else None

    def get(self, task_id):
        return self.core.fetchone("SELECT * FROM Tasks WHERE task_id = ?", (task_id,))

    def get_many(self, task_id=None):
        if task_id is None:
            return self.core.fetchall("SELECT * FROM Tasks")
        return self.core.fetchall("SELECT * FROM Tasks WHERE task_id = ?", (task_id,))

    def get_subtasks(self, task_id):
        return self.core.fetchall("SELECT * FROM Tasks WHERE parent_id = ?", (task_id,))

    def get_all_runnable(self):
        return self.core.fetchall(
            """
            SELECT *
            FROM Tasks
            WHERE status = 'P' OR status = 'S'
            ORDER BY priority ASC, task_id ASC
            """
        )

    def get_subtree_tasks(self, root_task_id: int):
        return self.core.fetchall(
            """
            WITH RECURSIVE subtree(task_id) AS (
                SELECT ?
                UNION ALL
                SELECT t.task_id
                FROM Tasks t
                JOIN subtree s ON t.parent_id = s.task_id
            )
            SELECT t.*
            FROM Tasks t
            JOIN subtree s ON t.task_id = s.task_id
            ORDER BY t.task_id ASC
            """,
            (root_task_id,),
        )

    def get_runnable_in_subtree(self, root_task_id: int):
        return self.core.fetchall(
            """
            WITH RECURSIVE subtree(task_id) AS (
                SELECT ?
                UNION ALL
                SELECT t.task_id
                FROM Tasks t
                JOIN subtree s ON t.parent_id = s.task_id
            )
            SELECT t.*
            FROM Tasks t
            JOIN subtree s ON t.task_id = s.task_id
            WHERE t.status = 'P' OR t.status = 'S'
            ORDER BY t.priority ASC, t.task_id ASC
            """,
            (root_task_id,),
        )

    def count_active_in_subtree(self, root_task_id: int, statuses=None) -> int:
        use_statuses = list(statuses or ["P", "R", "S", "B"])
        placeholders = ",".join("?" for _ in use_statuses)
        row = self.core.fetchone(
            f"""
            WITH RECURSIVE subtree(task_id) AS (
                SELECT ?
                UNION ALL
                SELECT t.task_id
                FROM Tasks t
                JOIN subtree s ON t.parent_id = s.task_id
            )
            SELECT COUNT(*)
            FROM Tasks t
            JOIN subtree s ON t.task_id = s.task_id
            WHERE t.status IN ({placeholders})
            """,
            (root_task_id, *use_statuses),
        )
        return int(row[0]) if row and row[0] is not None else 0

    def get_by_status(self, statuses):
        statuses = list(statuses or [])
        if not statuses:
            return []
        placeholders = ",".join("?" for _ in statuses)
        return self.core.fetchall(
            f"SELECT * FROM Tasks WHERE status IN ({placeholders}) ORDER BY task_id ASC",
            tuple(statuses),
        )

    def get_blocks(self, task_ids, unsatisfied_only=False):
        params = []
        if task_ids:
            task_ids = list(task_ids)
            placeholders = ",".join("?" for _ in task_ids)
            where = f"task_id IN ({placeholders})"
            params.extend(task_ids)
        else:
            where = "1=1"
        if unsatisfied_only:
            where += " AND satisfied_at IS NULL"
        return self.core.fetchall(
            f"""
            SELECT block_id, task_id, block_type, condition, message, created_at, satisfied_at, last_checked_at,
                   block_set, block_group
            FROM TaskBlocks
            WHERE {where}
            """,
            tuple(params),
        )

    def get_dependencies(self, task_ids):
        params = []
        if task_ids:
            task_ids = list(task_ids)
            placeholders = ",".join("?" for _ in task_ids)
            where = f"task_id IN ({placeholders})"
            params.extend(task_ids)
        else:
            where = "1=1"
        return self.core.fetchall(
            f"""
            SELECT dependency_id, task_id, depends_on_task_id, required_state, block_id, allow_failed, satisfied_at
            FROM TaskDependencies
            WHERE {where}
            """,
            tuple(params),
        )

    def get_time_constraints(self, task_ids):
        params = []
        if task_ids:
            task_ids = list(task_ids)
            placeholders = ",".join("?" for _ in task_ids)
            where = f"task_id IN ({placeholders})"
            params.extend(task_ids)
        else:
            where = "1=1"
        return self.core.fetchall(
            f"""
            SELECT constraint_id, task_id, mode, not_before, block_id, satisfied_at
            FROM TaskTimeConstraints
            WHERE {where}
            """,
            tuple(params),
        )

    def mark_dependency_satisfied(self, dep_id):
        self.core.execute(
            "UPDATE TaskDependencies SET satisfied_at = datetime('now') WHERE dependency_id = ?",
            (dep_id,),
        )
        self.conn.commit()
        return True

    def mark_block_satisfied(self, block_id):
        self.core.execute(
            "UPDATE TaskBlocks SET satisfied_at = datetime('now') WHERE block_id = ?",
            (block_id,),
        )
        self.conn.commit()
        return True

    def mark_block_checked(self, block_id):
        self.core.execute(
            "UPDATE TaskBlocks SET last_checked_at = datetime('now') WHERE block_id = ?",
            (block_id,),
        )
        self.conn.commit()
        return True

    def mark_time_constraint_satisfied(self, constraint_id):
        self.core.execute(
            "UPDATE TaskTimeConstraints SET satisfied_at = datetime('now') WHERE constraint_id = ?",
            (constraint_id,),
        )
        self.conn.commit()
        return True

    def update_end_time(self, task_id, end_time):
        self.core.execute(
            "UPDATE Tasks SET end_time = ? WHERE task_id = ?",
            (self._normalize_dt(end_time), task_id),
        )
        self.conn.commit()
        return True

    def update_start_time(self, task_id, start_time):
        self.core.execute(
            "UPDATE Tasks SET start_time = ? WHERE task_id = ?",
            (self._normalize_dt(start_time), task_id),
        )
        self.conn.commit()
        return True

    def update_parent(self, task_id, parent_id):
        self.core.execute(
            "UPDATE Tasks SET parent_id = ? WHERE task_id = ?",
            (parent_id, task_id),
        )
        self.conn.commit()
        return True

    def has_dependency_path(self, start_task_id: int, target_task_id: int) -> bool:
        row = self.core.fetchone(
            """
            WITH RECURSIVE chain(task_id) AS (
                SELECT depends_on_task_id FROM TaskDependencies WHERE task_id = ?
                UNION ALL
                SELECT td.depends_on_task_id
                FROM TaskDependencies td
                JOIN chain c ON td.task_id = c.task_id
            )
            SELECT 1 FROM chain WHERE task_id = ? LIMIT 1
            """,
            (start_task_id, target_task_id),
        )
        return row is not None

    def add_block(
        self,
        task_id: int,
        block_type: str,
        condition: str,
        message: str | None = None,
        *,
        block_set: str | None = None,
        block_group: str | None = None,
    ):
        self.core.execute(
            """
            INSERT INTO TaskBlocks (task_id, block_type, condition, message, block_set, block_group)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (task_id, block_type, condition, message, block_set, block_group),
        )
        self.conn.commit()
        return self.cursor.lastrowid

    def add_dependency(
        self,
        task_id: int,
        depends_on_task_id: int,
        required_state: str,
        *,
        block_id: int | None = None,
        allow_failed: bool = False,
    ):
        self.core.execute(
            """
            INSERT INTO TaskDependencies (task_id, depends_on_task_id, required_state, block_id, allow_failed)
            VALUES (?, ?, ?, ?, ?)
            """,
            (task_id, depends_on_task_id, required_state, block_id, int(bool(allow_failed))),
        )
        self.conn.commit()
        return self.cursor.lastrowid

    def add_time_constraint(
        self,
        task_id: int,
        mode: str,
        not_before: str,
        *,
        block_id: int | None = None,
    ):
        self.core.execute(
            """
            INSERT INTO TaskTimeConstraints (task_id, mode, not_before, block_id)
            VALUES (?, ?, ?, ?)
            """,
            (task_id, mode, not_before, block_id),
        )
        self.conn.commit()
        return self.cursor.lastrowid

    def get_errors_from_subtasks(self, task_id):
        return self.core.fetchall(
            "SELECT task_id, error_message, error_stack FROM Tasks WHERE parent_id = ?",
            (task_id,),
        )

    def get_error_info(self, task_id):
        return self.core.fetchone(
            "SELECT error_message, error_stack FROM Tasks WHERE task_id = ?",
            (task_id,),
        )

    def get_next_runnable(self):
        return self.core.fetchone(
            """
            SELECT *
            FROM Tasks
            WHERE status = 'P' OR status = 'S'
            ORDER BY priority ASC, task_id ASC
            LIMIT 1
            """
        )

    def kill_and_descendants(self, task_id: int, reason: str = "Killed by user"):
        rows = self.core.fetchall(
            """
            WITH RECURSIVE children(id) AS (
                SELECT ? UNION ALL
                SELECT t.task_id FROM Tasks t JOIN children c ON t.parent_id = c.id
            )
            SELECT id FROM children
            """,
            (task_id,),
        )
        ids = [row[0] for row in rows]
        if not ids:
            return False
        placeholders = ",".join("?" for _ in ids)
        self.core.execute(
            f"""
            UPDATE Tasks
            SET status = 'E',
                end_time = datetime('now'),
                status_updated_at = datetime('now'),
                error_message = ?
            WHERE task_id IN ({placeholders})
            """,
            (reason, *ids),
        )
        self.conn.commit()
        return True

    def cancel_and_descendants(self, task_id: int, reason: str = "Canceled by user"):
        rows = self.core.fetchall(
            """
            WITH RECURSIVE children(id) AS (
                SELECT ? UNION ALL
                SELECT t.task_id FROM Tasks t JOIN children c ON t.parent_id = c.id
            )
            SELECT id FROM children
            """,
            (task_id,),
        )
        ids = [row[0] for row in rows]
        if not ids:
            return False
        placeholders = ",".join("?" for _ in ids)
        self.core.execute(
            f"""
            UPDATE Tasks
            SET status = 'E',
                end_time = datetime('now'),
                status_updated_at = datetime('now'),
                error_message = ?
            WHERE task_id IN ({placeholders}) AND status IN ('P','S','B')
            """,
            (reason, *ids),
        )
        self.conn.commit()
        return True

    def reset(self, full: bool = False):
        if full:
            self.core.execute("DELETE FROM TaskDependencies")
            self.core.execute("DELETE FROM TaskTimeConstraints")
            self.core.execute("DELETE FROM TaskBlocks")
            self.core.execute("DELETE FROM Tasks")
            self.conn.commit()
            return True

        keep_cte = """
            WITH RECURSIVE keep(task_id) AS (
                SELECT task_id
                FROM Tasks
                WHERE status IS NULL OR UPPER(status) NOT IN ('E','B','C')
                UNION
                SELECT t.parent_id
                FROM Tasks t
                JOIN keep k ON t.task_id = k.task_id
                WHERE t.parent_id IS NOT NULL
                UNION
                SELECT t.task_id
                FROM Tasks t
                JOIN keep k ON t.parent_id = k.task_id
                UNION
                SELECT td.depends_on_task_id
                FROM TaskDependencies td
                JOIN keep k ON td.task_id = k.task_id
            ),
            deletable AS (
                SELECT task_id
                FROM Tasks
                WHERE UPPER(status) IN ('E','B','C')
                AND task_id NOT IN (SELECT task_id FROM keep)
            )
        """
        self.core.execute(
            keep_cte
            + """
            DELETE FROM TaskDependencies
            WHERE task_id IN (SELECT task_id FROM deletable)
            """
        )
        self.core.execute(
            keep_cte
            + """
            DELETE FROM TaskTimeConstraints
            WHERE task_id IN (SELECT task_id FROM deletable)
            """
        )
        self.core.execute(
            keep_cte
            + """
            DELETE FROM TaskBlocks
            WHERE task_id IN (SELECT task_id FROM deletable)
            """
        )
        self.core.execute(
            keep_cte
            + """
            DELETE FROM Tasks
            WHERE task_id IN (SELECT task_id FROM deletable)
            """
        )
        self.conn.commit()
        return True

    def clear(self, keep_running: bool = True):
        if not keep_running:
            self.core.execute("DELETE FROM TaskDependencies")
            self.core.execute("DELETE FROM TaskTimeConstraints")
            self.core.execute("DELETE FROM TaskBlocks")
            self.core.execute("DELETE FROM Tasks")
            self.conn.commit()
            return True

        rows = self.core.fetchall("SELECT task_id FROM Tasks WHERE status = 'R'")
        running_ids = [row[0] for row in rows]
        if not running_ids:
            return self.clear(keep_running=False)

        placeholders = ",".join("?" for _ in running_ids)
        keep_cte = f"""
            WITH RECURSIVE keep(id) AS (
                SELECT task_id FROM Tasks WHERE task_id IN ({placeholders})
                UNION
                SELECT t.task_id FROM Tasks t JOIN keep k ON t.parent_id = k.id
                UNION
                SELECT t.parent_id FROM Tasks t JOIN keep k ON t.task_id = k.id
                WHERE t.parent_id IS NOT NULL
            ),
            deletable AS (
                SELECT task_id
                FROM Tasks
                WHERE task_id NOT IN (SELECT id FROM keep WHERE id IS NOT NULL)
            ),
            deletable_blocks AS (
                SELECT block_id
                FROM TaskBlocks
                WHERE task_id IN (SELECT task_id FROM deletable)
            )
        """
        self.core.execute(
            keep_cte
            + """
            DELETE FROM TaskDependencies
            WHERE task_id IN (SELECT task_id FROM deletable)
               OR depends_on_task_id IN (SELECT task_id FROM deletable)
               OR block_id IN (SELECT block_id FROM deletable_blocks)
            """,
            running_ids,
        )
        self.core.execute(
            keep_cte
            + """
            DELETE FROM TaskTimeConstraints
            WHERE task_id IN (SELECT task_id FROM deletable)
               OR block_id IN (SELECT block_id FROM deletable_blocks)
            """,
            running_ids,
        )
        self.core.execute(
            keep_cte
            + """
            DELETE FROM TaskBlocks
            WHERE task_id IN (SELECT task_id FROM deletable)
            """,
            running_ids,
        )
        self.core.execute(
            keep_cte
            + """
            DELETE FROM Tasks
            WHERE task_id IN (SELECT task_id FROM deletable)
            """,
            running_ids,
        )
        self.conn.commit()
        return True
