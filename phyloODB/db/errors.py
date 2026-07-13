"""Typed database errors exposed by PhyloODB's persistence layer."""
from __future__ import annotations

from ..errors import PhyloODBError


class PhyloODBDatabaseError(PhyloODBError):
    """Base class for contextual database failures."""


class SchemaCompatibilityError(PhyloODBDatabaseError):
    """The database schema cannot be used by this application version."""


class MigrationError(PhyloODBDatabaseError):
    """A versioned schema migration failed."""


class RepositoryReadError(PhyloODBDatabaseError):
    """A repository could not read required state."""


class RepositoryWriteError(PhyloODBDatabaseError):
    """A repository write failed and was rolled back."""


class RepositoryConflictError(RepositoryWriteError):
    """A repository write violated a uniqueness or integrity constraint."""


class StorageOperationError(PhyloODBDatabaseError):
    """A combined filesystem/database operation could not be completed."""
