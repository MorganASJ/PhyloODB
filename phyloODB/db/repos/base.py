from __future__ import annotations

from functools import wraps
from typing import Any

from ..core import DatabaseCore


def transactional(operation: str | None = None):
    """Make a repository write atomic and composable with outer transactions."""

    def decorate(func):
        @wraps(func)
        def wrapped(self, *args, **kwargs):
            label = operation or func.__name__.replace("_", " ")
            with self.core.transaction(operation=label):
                return func(self, *args, **kwargs)
        return wrapped
    return decorate


def transactional_methods(*method_names: str):
    """Class decorator for repositories with many small write methods."""

    def decorate(cls):
        for name in method_names:
            method = getattr(cls, name)
            setattr(cls, name, transactional(name.replace("_", " "))(method))
        return cls
    return decorate


class BaseRepository:
    def __init__(self, manager: Any):
        self.manager = manager
        self.core = DatabaseCore(manager)

    @property
    def conn(self):
        # Repositories route commit/rollback through DatabaseCore so an outer
        # DBManager transaction retains ownership of the real connection.
        return self.core

    @property
    def cursor(self):
        return self.core.cursor
