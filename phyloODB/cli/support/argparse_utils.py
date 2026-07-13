"""Shared argparse helpers used across CLI command modules."""
from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
from typing import Any, List


# ---------------------------------------------------------------------------
# Help text formatting
# Purpose: Keep argparse default rendering and help text tweaks consistent
# across command modules and task-specific parsers.
# ---------------------------------------------------------------------------

def _format_help_default(value: Any) -> str:
    """Render parser defaults consistently for help text."""

    if isinstance(value, Path):
        return str(value)
    return str(value)


def _with_default_help(help_text: str, default: Any) -> str:
    """Append a formatted default value to argument help text when needed."""

    if help_text is argparse.SUPPRESS:
        return help_text
    text = (help_text or "").strip()
    if "default:" in text.lower():
        return text
    suffix = f"(default: {_format_help_default(default)})"
    if text:
        return f"{text} {suffix}"
    return suffix


# ---------------------------------------------------------------------------
# Argument coercion
# Purpose: Provide small reusable argparse actions and validators for the CLI
# command modules so common flags behave identically everywhere.
# ---------------------------------------------------------------------------

def _validate_date(token: str, option: str) -> str:
    """Validate a CLI date argument and preserve the original token."""

    try:
        datetime.strptime(token, "%Y-%m-%d")
    except ValueError as exc:  # pragma: no cover - user input validation
        raise argparse.ArgumentTypeError(f"{option} expects YYYY-MM-DD, got '{token}'") from exc
    return token


class AppendCommaSeparated(argparse.Action):
    """Collect repeated CLI values while expanding comma-delimited tokens."""

    def __call__(self, parser, namespace, values, option_string=None):  # pragma: no cover - invoked via CLI
        collected = getattr(namespace, self.dest, None)
        if collected is None:
            collected = []
        if isinstance(values, (list, tuple)):
            raw_values = values
        else:
            raw_values = [values]
        for value in raw_values:
            for token in str(value).split(","):
                token = token.strip()
                if token:
                    collected.append(token)
        setattr(namespace, self.dest, collected)


__all__ = [
    "AppendCommaSeparated",
    "_format_help_default",
    "_validate_date",
    "_with_default_help",
]
