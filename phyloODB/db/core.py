from __future__ import annotations

from contextlib import contextmanager
import sqlite3
from typing import Any, Iterable, Iterator

from .errors import RepositoryConflictError, RepositoryReadError, RepositoryWriteError


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

    def execute(self, sql: str, params: Iterable[Any] = ()):
        return self.cursor.execute(sql, tuple(params))

    def executemany(self, sql: str, seq_of_params: Iterable[Iterable[Any]]):
        return self.cursor.executemany(sql, seq_of_params)

    def fetchone(self, sql: str, params: Iterable[Any] = ()):
        try:
            self.execute(sql, params)
            return self.cursor.fetchone()
        except sqlite3.Error as exc:
            raise RepositoryReadError(f"Database read failed: {exc}") from exc

    def fetchall(self, sql: str, params: Iterable[Any] = ()):
        try:
            self.execute(sql, params)
            return self.cursor.fetchall() or []
        except sqlite3.Error as exc:
            raise RepositoryReadError(f"Database read failed: {exc}") from exc

    def commit(self) -> None:
        if int(getattr(self.manager, "_transaction_depth", 0) or 0) == 0:
            self.conn.commit()

    def rollback(self) -> None:
        if int(getattr(self.manager, "_transaction_depth", 0) or 0) == 0:
            self.conn.rollback()

    @contextmanager
    def transaction(self, *, operation: str = "database write") -> Iterator[None]:
        depth = int(getattr(self.manager, "_transaction_depth", 0) or 0)
        counter = int(getattr(self.manager, "_savepoint_counter", 0) or 0) + 1
        self.manager._savepoint_counter = counter
        savepoint = f"phyloodb_sp_{counter}"
        try:
            if depth == 0:
                if not self.conn.in_transaction:
                    self.cursor.execute("BEGIN")
            else:
                self.cursor.execute(f"SAVEPOINT {savepoint}")
            self.manager._transaction_depth = depth + 1
            yield
            if depth == 0:
                self.conn.commit()
            else:
                self.cursor.execute(f"RELEASE SAVEPOINT {savepoint}")
        except sqlite3.IntegrityError as exc:
            if depth == 0:
                self.conn.rollback()
            else:
                self.cursor.execute(f"ROLLBACK TO SAVEPOINT {savepoint}")
                self.cursor.execute(f"RELEASE SAVEPOINT {savepoint}")
            raise RepositoryConflictError(f"{operation} failed: {exc}") from exc
        except Exception as exc:  # boundary: transaction manager converts write failures and rolls back/savepoints.
            if depth == 0:
                self.conn.rollback()
            else:
                self.cursor.execute(f"ROLLBACK TO SAVEPOINT {savepoint}")
                self.cursor.execute(f"RELEASE SAVEPOINT {savepoint}")
            if isinstance(exc, RepositoryWriteError):
                raise
            raise RepositoryWriteError(f"{operation} failed: {exc}") from exc
        finally:
            self.manager._transaction_depth = depth
            if depth == 0:
                self.manager._savepoint_counter = 0
