from __future__ import annotations

import os
import re
from typing import Optional


DEFAULT_CLEAN_PROFILE = "clean_default"
RAW_PROFILE = "raw"
_STAGED_BUSCO_INPUT_RE = re.compile(r"\.busco_input_(.+)\.faa(?:\.gz)?$", flags=re.IGNORECASE)


def _format_cdhit_token(identity: Optional[float]) -> str:
    if identity is None:
        return "cdhit"
    try:
        pct = float(identity) * 100.0
    except (TypeError, ValueError):
        return "cdhit"
    token = f"{pct:.2f}".rstrip("0").rstrip(".")
    return f"cdhit{token}"


def derive_profile_name_from_recipe(
    *,
    used_gff: bool,
    used_cdhit: bool,
    cdhit_identity: Optional[float] = None,
    fallback: str = DEFAULT_CLEAN_PROFILE,
) -> str:
    tokens: list[str] = []
    if used_gff:
        tokens.append("gff")
    if used_cdhit:
        tokens.append(_format_cdhit_token(cdhit_identity))
    if tokens:
        return "_".join(tokens)
    return str(fallback or DEFAULT_CLEAN_PROFILE)


def resolve_profile_selector(
    *,
    proteome_profile: Optional[str] = None,
    isoforms_cleaned: Optional[bool] = None,
    raw_proteome: Optional[bool] = None,
    default_clean_profile: str = DEFAULT_CLEAN_PROFILE,
) -> Optional[str]:
    explicit = str(proteome_profile or "").strip() or None
    cleaned_flag = None if isoforms_cleaned is None else bool(isoforms_cleaned)
    raw_flag = None if raw_proteome is None else bool(raw_proteome)

    if cleaned_flag and raw_flag:
        raise ValueError("Choose either --isoforms-cleaned or --raw-proteome, not both.")

    shortcut_profile: Optional[str] = None
    if raw_flag or cleaned_flag is False:
        shortcut_profile = RAW_PROFILE
    elif cleaned_flag:
        shortcut_profile = str(default_clean_profile).strip() or DEFAULT_CLEAN_PROFILE

    if explicit and shortcut_profile and explicit != shortcut_profile:
        raise ValueError(
            f"--proteome-profile={explicit} conflicts with the requested convenience selector '{shortcut_profile}'."
        )
    return explicit or shortcut_profile


def staged_busco_input_profile_name(path: Optional[str]) -> Optional[str]:
    base = os.path.basename(str(path or ""))
    match = _STAGED_BUSCO_INPUT_RE.fullmatch(base)
    if not match:
        return None
    profile_name = str(match.group(1) or "").strip()
    return profile_name or None


def is_staged_busco_input_path(path: Optional[str]) -> bool:
    return staged_busco_input_profile_name(path) is not None
