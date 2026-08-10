from __future__ import annotations

from contextlib import contextmanager
import os
import sqlite3
import threading
import time
from typing import Any, Iterable, Iterator

from .errors import RepositoryConflictError, RepositoryReadError, RepositoryWriteError


DEFAULT_SQLITE_BUSY_TIMEOUT_MS = 300_000
_WRITE_LOCK = threading.RLock()


def sqlite_busy_timeout_ms() -> int:
    raw = os.environ.get("PHYOODB_SQLITE_BUSY_TIMEOUT_MS")
    if raw is None:
        return DEFAULT_SQLITE_BUSY_TIMEOUT_MS
    try:
        return max(int(raw), 0)
    except (TypeError, ValueError):
        return DEFAULT_SQLITE_BUSY_TIMEOUT_MS


def _busy_timeout_seconds() -> float:
    return sqlite_busy_timeout_ms() / 1000.0


def _is_locked_error(exc: BaseException) -> bool:
    if not isinstance(exc, sqlite3.OperationalError):
        return False
    message = str(exc).lower()
    return "database is locked" in message or "database table is locked" in message


class DatabaseCore:
    """Small wrapper around an existing DBManager connection."""

    def __init__(self, manager: Any):
        self.manager = manager

    @property
    def conn(self):
        return self.manager.conn

    @property
    def cursor(self):
        return self.manager.cursor

    @property
    def cursor_lock(self):
        """Serialize complete cursor operations for managers shared by threads."""

        lock = getattr(self.manager, "_cursor_lock", None)
        if lock is None:
            lock = threading.RLock()
            self.manager._cursor_lock = lock
        return lock

    def _retry_locked(self, operation):
        deadline = time.monotonic() + _busy_timeout_seconds()
        delay = 0.1
        while True:
            try:
                return operation()
            except sqlite3.OperationalError as exc:
                if not _is_locked_error(exc) or time.monotonic() >= deadline:
                    raise
                time.sleep(min(delay, max(deadline - time.monotonic(), 0.0)))
                delay = min(delay * 1.5, 5.0)

    def execute(self, sql: str, params: Iterable[Any] = ()):
        bound = tuple(params)
        with self.cursor_lock:
            return self._retry_locked(lambda: self.cursor.execute(sql, bound))

    def executemany(self, sql: str, seq_of_params: Iterable[Iterable[Any]]):
        bound = tuple(tuple(params) for params in seq_of_params)
        with self.cursor_lock:
            return self._retry_locked(lambda: self.cursor.executemany(sql, bound))

    def fetchone(self, sql: str, params: Iterable[Any] = ()):
        try:
            with self.cursor_lock:
                self.execute(sql, params)
                return self.cursor.fetchone()
        except sqlite3.Error as exc:
            raise RepositoryReadError(f"Database read failed: {exc}") from exc

    def fetchall(self, sql: str, params: Iterable[Any] = ()):
        try:
            with self.cursor_lock:
                self.execute(sql, params)
                return self.cursor.fetchall() or []
        except sqlite3.Error as exc:
            raise RepositoryReadError(f"Database read failed: {exc}") from exc

    def commit(self) -> None:
        if int(getattr(self.manager, "_transaction_depth", 0) or 0) == 0:
            self._retry_locked(self.conn.commit)

    def rollback(self) -> None:
        if int(getattr(self.manager, "_transaction_depth", 0) or 0) == 0:
            self.conn.rollback()

    @contextmanager
    def transaction(self, *, operation: str = "database write") -> Iterator[None]:
        depth = int(getattr(self.manager, "_transaction_depth", 0) or 0)
        counter = int(getattr(self.manager, "_savepoint_counter", 0) or 0) + 1
        self.manager._savepoint_counter = counter
        savepoint = f"phyloodb_sp_{counter}"
        lock = _WRITE_LOCK if depth == 0 else None
        if lock is not None:
            lock.acquire()
        try:
            if depth == 0:
                if not self.conn.in_transaction:
                    self.execute("BEGIN")
            else:
                self.execute(f"SAVEPOINT {savepoint}")
            self.manager._transaction_depth = depth + 1
            yield
            if depth == 0:
                self._retry_locked(self.conn.commit)
            else:
                self.execute(f"RELEASE SAVEPOINT {savepoint}")
        except sqlite3.IntegrityError as exc:
            if depth == 0:
                self.conn.rollback()
            else:
                self.execute(f"ROLLBACK TO SAVEPOINT {savepoint}")
                self.execute(f"RELEASE SAVEPOINT {savepoint}")
            raise RepositoryConflictError(f"{operation} failed: {exc}") from exc
        except Exception as exc:  # boundary: transaction manager converts write failures and rolls back/savepoints.
            if depth == 0:
                self.conn.rollback()
            else:
                self.execute(f"ROLLBACK TO SAVEPOINT {savepoint}")
                self.execute(f"RELEASE SAVEPOINT {savepoint}")
            if isinstance(exc, RepositoryWriteError):
                raise
            raise RepositoryWriteError(f"{operation} failed: {exc}") from exc
        finally:
            self.manager._transaction_depth = depth
            if depth == 0:
                self.manager._savepoint_counter = 0
            if lock is not None:
                lock.release()
