"""Application-level errors that may cross service, task, and CLI boundaries."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


class PhyloODBError(RuntimeError):
    """Base class for expected operational failures."""


class SelectorError(PhyloODBError, ValueError):
    """A selector is invalid or cannot be resolved as requested."""


class TaskExecutionError(PhyloODBError):
    """A required task stage could not be completed."""


class ExternalToolError(TaskExecutionError):
    """An external executable failed or returned invalid output."""


@dataclass(frozen=True)
class BatchFailure:
    item: str
    operation: str
    message: str
    stack: str = ""


class BatchItemError(TaskExecutionError):
    """A required operation failed for one item in a batch."""

    def __init__(self, item: Any, operation: str, cause: BaseException):
        self.item = str(item)
        self.operation = str(operation)
        self.cause = cause
        super().__init__(f"{self.operation} failed for {self.item}: {cause}")

