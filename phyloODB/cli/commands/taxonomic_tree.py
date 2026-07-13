"""Standalone taxonomic tree CLI command."""
from __future__ import annotations

import argparse
import colorsys
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from ...database import DBManager
from ...selector_utils import (
    RANK_HIERARCHY,
    expand_busco_run_id_variables,
    normalize_rank_list,
    resolve_selector_accessions,
)
from ..support.argparse_utils import AppendCommaSeparated
from ..support.common import (
    _apply_busco_context_from_args,
    _connect_manager,
    _print_error,
    _resolve_library_selector,
)
from ..support.selectors import _add_selector_arguments, _selector_request_from_args

_RANK_ORDER = {rank: idx for idx, rank in enumerate(RANK_HIERARCHY)}


@dataclass
class _TreeNode:
    taxid: Optional[int]
    name: str
    rank: str
    children: Dict[int, "_TreeNode"] = field(default_factory=dict)
    tip_label: Optional[str] = None


@dataclass
class _DisplayNode:
    taxid: Optional[int]
    name: str
    rank: str
    tip_label: Optional[str] = None
    children: List["_DisplayNode"] = field(default_factory=list)


@dataclass
class _RenderRow:
    node: _DisplayNode
    depth: int
    ancestor_continuations: List[bool]
    is_last: bool
    is_tip: bool
    label: str
    colour: Optional[str]


def _chunked(values: Sequence[str], size: int = 900) -> Iterable[List[str]]:
    for idx in range(0, len(values), size):
        yield list(values[idx:idx + size])


def _abbreviate_species(name: str) -> str:
    parts = [part for part in str(name or "").strip().split() if part]
    if len(parts) < 2:
        return str(name or "").strip()
    return f"{parts[0][0]}. {' '.join(parts[1:])}"


def _format_tip_label(
    species_name: str,
    accessions: Sequence[str],
    *,
    full_name: bool = False,
    show_accession: bool = False,
) -> str:
    label = str(species_name or "").strip() if full_name else _abbreviate_species(species_name)
    if not show_accession:
        return label
    accession_tokens = [str(accession).strip() for accession in accessions or [] if str(accession).strip()]
    if not accession_tokens:
        return label
    return f"{label} ({','.join(dict.fromkeys(accession_tokens))})"


def _sanitize_newick_label(label: str) -> str:
    return "'" + str(label or "").replace("'", "''") + "'"


def _node_sort_key(node: _TreeNode | _DisplayNode) -> tuple[int, str, str]:
    return (
        _RANK_ORDER.get(str(node.rank or "").lower(), len(_RANK_ORDER)),
        str(node.name or "").lower(),
        str(node.tip_label or "").lower(),
    )


def _fetch_species_taxa(
    manager: DBManager,
    accessions: Sequence[str],
) -> tuple[Dict[int, str], Dict[int, List[str]], List[str]]:
    species_by_taxid: Dict[int, str] = {}
    accessions_by_taxid: Dict[int, List[str]] = {}
    found_accessions: set[str] = set()
    if not accessions:
        return species_by_taxid, accessions_by_taxid, []

    for chunk in _chunked(list(accessions)):
        placeholders = ",".join("?" for _ in chunk)
        manager.cursor.execute(
            f"""
            SELECT assembly_accession, taxid, name
            FROM TaxonomyAssemblySummary
            WHERE assembly_accession IN ({placeholders})
            """,
            tuple(chunk),
        )
        for accession, taxid, name in manager.cursor.fetchall() or []:
            if accession is None or taxid is None:
                continue
            found_accessions.add(str(accession))
            normalized_taxid = int(taxid)
            species_by_taxid.setdefault(normalized_taxid, str(name or ""))
            accessions_by_taxid.setdefault(normalized_taxid, []).append(str(accession))

    missing = [str(accession) for accession in accessions if str(accession) not in found_accessions]
    return species_by_taxid, accessions_by_taxid, missing


def _infer_busco_library_from_run_ids(manager: DBManager, run_ids: Sequence[int]) -> Optional[int]:
    run_vals = [int(run_id) for run_id in run_ids or [] if run_id is not None]
    if not run_vals:
        return None
    placeholders = ",".join("?" for _ in run_vals)
    rows = manager.cursor.execute(
        f"SELECT DISTINCT library_id FROM BUSCO_Runs WHERE run_id IN ({placeholders})",
        tuple(run_vals),
    ).fetchall() or []
    library_ids = [int(row[0]) for row in rows if row and row[0] is not None]
    if not library_ids:
        return None
    if len(library_ids) > 1:
        raise ValueError(
            "Selected BUSCO run ids span multiple libraries. Use --library-id/--library-name to disambiguate."
        )
    return library_ids[0]


def _format_busco_tip_suffix(value: Optional[float]) -> Optional[str]:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if abs(number - round(number)) < 0.005:
        rendered = str(int(round(number)))
    else:
        rendered = f"{number:.2f}".rstrip("0").rstrip(".")
    return f"[{rendered}%]"


def _resolve_busco_tip_annotations(
    manager: DBManager,
    args: argparse.Namespace,
    *,
    accessions_by_taxid: Mapping[int, Sequence[str]],
) -> tuple[Dict[int, str], int]:
    explicit_run_ids = expand_busco_run_id_variables(manager, getattr(args, "busco_run_ids", None) or [])
    busco_library = _resolve_library_selector(
        manager,
        library_id=getattr(args, "library_id", None),
        library_name=getattr(args, "library_name", None),
        legacy=getattr(args, "busco_library", None),
    )
    if busco_library is None and explicit_run_ids:
        busco_library = _infer_busco_library_from_run_ids(manager, explicit_run_ids)
    if busco_library is None:
        raise ValueError("Use --busco with --library-id or --library-name.")

    unique_tip_accessions: Dict[int, str] = {}
    ambiguous_taxids = 0
    for taxid, accession_group in accessions_by_taxid.items():
        unique_accessions = list(dict.fromkeys(str(accession) for accession in accession_group if accession is not None))
        if len(unique_accessions) != 1:
            ambiguous_taxids += 1
            continue
        unique_tip_accessions[int(taxid)] = unique_accessions[0]
    if not unique_tip_accessions:
        return {}, ambiguous_taxids

    selection = str(getattr(args, "busco_run_selection", None) or "primary").strip().lower() or "primary"
    score_library_id = int(busco_library)
    run_library_id = manager.assert_library_has_parent(score_library_id) or score_library_id
    accessions = list(unique_tip_accessions.values())
    run_map = manager.busco.get_effective_run_ids_for_accessions(
        int(run_library_id),
        accessions=accessions,
        run_ids=explicit_run_ids or None,
        pipeline=getattr(args, "busco_pipeline", None),
        input_mode=getattr(args, "busco_input_mode", None) or getattr(args, "format", None),
        preferred_pipeline=getattr(args, "prefer_busco_pipeline", None),
        preferred_input_mode=getattr(args, "prefer_busco_input_mode", None) or getattr(args, "prefer_format", None),
        proteome_profile=getattr(args, "proteome_profile", None),
        preferred_proteome_profile=getattr(args, "prefer_proteome_profile", None),
        selection=selection,
        purpose="default",
    )

    display_results = manager.busco.get_display_results_for_runs(
        library_id=score_library_id,
        run_refs=[(accession, int(run_id)) for accession, run_id in run_map.items()],
        include_paralog=getattr(args, "include_paralog_filtering_in_score", None),
        paralog_run_id=getattr(args, "use_paralog_run", None),
        include_decontam=getattr(args, "include_decontamination_in_score", None),
        decont_run_id=getattr(args, "use_decontamination_run", None),
        allow_ambiguous_contaminants=getattr(args, "allow_ambiguous_contaminants", None),
        strict_decontamination=getattr(args, "strict_decontamination", None),
        rescue_duplicates=getattr(args, "rescue_duplicates", False),
    )

    annotations: Dict[int, str] = {}
    for taxid, accession in unique_tip_accessions.items():
        run_id = run_map.get(accession)
        if run_id is None:
            continue
        row = display_results.get((accession, int(run_id)))
        if not row:
            continue
        suffix = _format_busco_tip_suffix(row[1])
        if suffix:
            annotations[int(taxid)] = suffix
    return annotations, ambiguous_taxids


def _build_taxonomic_tree(manager: DBManager, species_by_taxid: Mapping[int, str]) -> _TreeNode:
    root = _TreeNode(taxid=None, name="", rank="")
    for taxid, species_name in sorted(species_by_taxid.items(), key=lambda item: str(item[1] or "").lower()):
        lineage = manager.get_lineage_root_to_leaf(int(taxid)) or []
        if not lineage:
            continue
        current = root
        for lineage_taxid, lineage_name, lineage_rank, _parent_taxid in lineage:
            child_key = int(lineage_taxid)
            child = current.children.get(child_key)
            if child is None:
                child = _TreeNode(
                    taxid=child_key,
                    name=str(lineage_name or ""),
                    rank=str(lineage_rank or ""),
                )
                current.children[child_key] = child
            current = child
        current.tip_label = _abbreviate_species(species_name or current.name)
    return root


def _should_display_node(node: _TreeNode, annotated_ranks: set[str], *, force: bool = False) -> bool:
    if force or node.tip_label:
        return True
    if str(node.rank or "").lower() in annotated_ranks:
        return True
    return len(node.children) != 1


def _compress_tree(node: _TreeNode, annotated_ranks: set[str], *, force: bool = False) -> List[_DisplayNode]:
    children: List[_DisplayNode] = []
    for child in sorted(node.children.values(), key=_node_sort_key):
        children.extend(_compress_tree(child, annotated_ranks))

    if _should_display_node(node, annotated_ranks, force=force):
        return [
            _DisplayNode(
                taxid=node.taxid,
                name=node.name,
                rank=node.rank,
                tip_label=node.tip_label,
                children=children,
            )
        ]
    return children


def _display_label(node: _DisplayNode, annotated_ranks: set[str]) -> str:
    if node.tip_label:
        return str(node.tip_label)
    return str(node.name or "")


def _collect_rank_nodes(nodes: Sequence[_DisplayNode], target_rank: str) -> List[_DisplayNode]:
    matches: List[_DisplayNode] = []

    def walk(node: _DisplayNode) -> None:
        if str(node.rank or "").lower() == target_rank:
            matches.append(node)
            return
        for child in node.children:
            walk(child)

    for node in nodes:
        walk(node)
    return sorted(matches, key=_node_sort_key)


def _hsv_hex(hue: float, saturation: float, value: float) -> str:
    r, g, b = colorsys.hsv_to_rgb(hue % 1.0, max(0.0, min(1.0, saturation)), max(0.0, min(1.0, value)))
    return f"#{int(round(r * 255)):02x}{int(round(g * 255)):02x}{int(round(b * 255)):02x}"


def _spread_hues(center: float, count: int, span: float) -> List[float]:
    if count <= 1:
        return [center % 1.0]
    start = center - (span / 2.0)
    step = span / float(count - 1)
    return [(start + step * idx) % 1.0 for idx in range(count)]


def _assign_rank_colours(root_nodes: Sequence[_DisplayNode], ranks: Sequence[str]) -> Dict[str, Dict[int, str]]:
    colour_map: Dict[str, Dict[int, str]] = {rank: {} for rank in ranks}
    if not ranks:
        return colour_map
    if len(ranks) > 3:
        raise ValueError("--colour-by-ranks supports at most 3 ranks.")

    top_nodes = _collect_rank_nodes(root_nodes, ranks[0])
    if top_nodes:
        top_hues = [(0.58 + (idx / max(1, len(top_nodes)))) % 1.0 for idx in range(len(top_nodes))]
        frontier: List[tuple[_DisplayNode, float]] = []
        for node, hue in zip(top_nodes, top_hues):
            if node.taxid is not None:
                colour_map[ranks[0]][int(node.taxid)] = _hsv_hex(hue, 0.72, 0.92)
                frontier.append((node, hue))
    else:
        frontier = [(node, (0.58 + (idx / max(1, len(root_nodes)))) % 1.0) for idx, node in enumerate(root_nodes)]

    spans = [0.12, 0.07]
    sat_vals = [(0.68, 0.86), (0.64, 0.80)]
    for rank_index, rank in enumerate(ranks[1:], start=1):
        next_frontier: List[tuple[_DisplayNode, float]] = []
        seen_taxids: set[int] = set()
        for parent, parent_hue in frontier:
            matches = _collect_rank_nodes([parent], rank)
            if not matches:
                continue
            hues = _spread_hues(parent_hue, len(matches), spans[min(rank_index - 1, len(spans) - 1)])
            sat, val = sat_vals[min(rank_index - 1, len(sat_vals) - 1)]
            for node, hue in zip(matches, hues):
                if node.taxid is None or int(node.taxid) in seen_taxids:
                    continue
                seen_taxids.add(int(node.taxid))
                colour_map[rank][int(node.taxid)] = _hsv_hex(hue, sat, val)
                next_frontier.append((node, hue))
        frontier = next_frontier

    return colour_map


def _newick_for_node(node: _TreeNode, *, is_root: bool = False) -> str:
    children = sorted(node.children.values(), key=_node_sort_key)
    if children:
        inner = ",".join(_newick_for_node(child) for child in children)
        label = "" if is_root or not node.name else _sanitize_newick_label(node.name)
        return f"({inner}){label}"
    label = node.tip_label or node.name or str(node.taxid or "")
    return _sanitize_newick_label(label)


def _render_tree(
    root_nodes: Sequence[_DisplayNode],
    *,
    annotated_ranks: Sequence[str],
    colour_map: Mapping[str, Mapping[int, str]],
    no_colour: bool = False,
) -> None:
    from rich.console import Console
    from rich.text import Text

    annotated_rank_set = {rank.lower() for rank in annotated_ranks}
    rows: List[_RenderRow] = []
    tip_spans: Dict[int, tuple[int, int]] = {}
    nodes_by_taxid: Dict[int, _DisplayNode] = {}

    def walk(
        node: _DisplayNode,
        depth: int,
        ancestor_continuations: List[bool],
        is_last: bool,
        inherited_colour: Optional[str],
    ) -> tuple[int, int]:
        own_colour = None
        rank = str(node.rank or "").lower()
        if node.taxid is not None and rank in colour_map:
            own_colour = colour_map[rank].get(int(node.taxid))
        active_colour = own_colour or inherited_colour
        row_index = len(rows)
        row = _RenderRow(
            node=node,
            depth=depth,
            ancestor_continuations=list(ancestor_continuations),
            is_last=is_last,
            is_tip=bool(node.tip_label),
            label=_display_label(node, annotated_rank_set),
            colour=active_colour,
        )
        rows.append(row)
        if node.taxid is not None:
            nodes_by_taxid[int(node.taxid)] = node

        if row.is_tip:
            span = (row_index, row_index)
        else:
            child_spans: List[tuple[int, int]] = []
            child_count = len(node.children)
            for idx, child in enumerate(node.children):
                child_spans.append(
                    walk(
                        child,
                        depth + 1,
                        ancestor_continuations + ([not is_last] if depth > 0 else []),
                        idx == child_count - 1,
                        active_colour,
                    )
                )
            span = (child_spans[0][0], child_spans[-1][1]) if child_spans else (row_index, row_index)

        if node.taxid is not None:
            tip_spans[int(node.taxid)] = span
        return span

    for idx, node in enumerate(root_nodes):
        walk(node, 0, [], idx == len(root_nodes) - 1, None)

    internal_end = 0
    max_tip_len = 0
    for row in rows:
        if row.is_tip:
            max_tip_len = max(max_tip_len, len(row.label))
            continue
        label_x = 0 if row.depth == 0 else row.depth * 2 + 1
        internal_end = max(internal_end, label_x + len(row.label))
    tip_col = max(internal_end + 3, 8)
    tip_end = tip_col + max_tip_len

    rank_annotations: List[tuple[str, int, List[tuple[_DisplayNode, tuple[int, int], Optional[str]]]]] = []
    next_col = tip_end + 3
    for rank in annotated_ranks:
        nodes = _collect_rank_nodes(root_nodes, rank)
        annotated_nodes = [
            (node, tip_spans.get(int(node.taxid)), colour_map.get(rank, {}).get(int(node.taxid)))
            for node in nodes
            if node.taxid is not None and int(node.taxid) in tip_spans
        ]
        annotated_nodes = [(node, span, colour) for node, span, colour in annotated_nodes if span is not None]
        if not annotated_nodes:
            continue
        width = max(len(node.name) for node, _span, _colour in annotated_nodes) + 3
        rank_annotations.append((rank, next_col, annotated_nodes))
        next_col += width + 2

    total_width = max(next_col + 1, tip_end + 1)
    char_rows = [[" "] * total_width for _ in rows]
    style_rows = [[""] * total_width for _ in rows]

    def paint(row_idx: int, col_idx: int, char: str, style: Optional[str]) -> None:
        if row_idx < 0 or row_idx >= len(char_rows) or col_idx < 0 or col_idx >= total_width:
            return
        char_rows[row_idx][col_idx] = char
        if style:
            style_rows[row_idx][col_idx] = style

    for row_idx, row in enumerate(rows):
        for level, keep_open in enumerate(row.ancestor_continuations):
            if keep_open:
                paint(row_idx, level * 2, "│", row.colour)
        if row.depth > 0:
            junction_x = (row.depth - 1) * 2
            paint(row_idx, junction_x, "└" if row.is_last else "├", row.colour)
            if row.is_tip:
                for col in range(junction_x + 1, tip_col):
                    paint(row_idx, col, "─", row.colour)
            else:
                paint(row_idx, junction_x + 1, "─", row.colour)
                label_x = row.depth * 2 + 1
                for offset, char in enumerate(row.label):
                    paint(row_idx, label_x + offset, char, row.colour)
        else:
            for offset, char in enumerate(row.label):
                paint(row_idx, offset, char, row.colour)

        if row.is_tip:
            for offset, char in enumerate(row.label):
                paint(row_idx, tip_col + offset, char, row.colour)

    for _rank, bracket_col, annotations in rank_annotations:
        for node, span, colour in annotations:
            start, end = span
            if start == end:
                paint(start, bracket_col, "┤", colour)
            else:
                paint(start, bracket_col, "┐", colour)
                for row_idx in range(start + 1, end):
                    paint(row_idx, bracket_col, "│", colour)
                paint(end, bracket_col, "┘", colour)
            label = str(node.name or "")
            paint(start, bracket_col + 1, " ", "")
            for offset, char in enumerate(label):
                paint(start, bracket_col + 2 + offset, char, colour)

    console = Console(no_color=no_colour)
    for chars, styles in zip(char_rows, style_rows):
        text = Text()
        current_style = styles[0]
        current_chars = [chars[0]]
        for idx in range(1, len(chars)):
            if styles[idx] == current_style:
                current_chars.append(chars[idx])
                continue
            text.append("".join(current_chars), style=current_style)
            current_style = styles[idx]
            current_chars = [chars[idx]]
        text.append("".join(current_chars), style=current_style)
        text.rstrip()
        console.print(text)


def _iter_nodes(nodes: Sequence[_DisplayNode]) -> Iterable[_DisplayNode]:
    for node in nodes:
        yield node
        yield from _iter_nodes(node.children)


def _apply_tip_labels(
    root_nodes: Sequence[_DisplayNode],
    *,
    species_by_taxid: Mapping[int, str],
    accessions_by_taxid: Mapping[int, Sequence[str]],
    full_name: bool = False,
    show_accession: bool = False,
) -> None:
    for node in _iter_nodes(root_nodes):
        if not node.tip_label or node.taxid is None:
            continue
        node.tip_label = _format_tip_label(
            str(species_by_taxid.get(int(node.taxid), node.name or "")),
            accessions_by_taxid.get(int(node.taxid), []),
            full_name=full_name,
            show_accession=show_accession,
        )


def _write_newick(path: str, root: _TreeNode) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(_newick_for_node(root, is_root=True) + ";\n", encoding="utf-8")


def _handle_tree(args: argparse.Namespace) -> int:
    manager = _connect_manager(args.database)
    try:
        if not (args.accessions or args.taxid is not None or args.clade or getattr(args, "preset_name", None)):
            return _print_error("Provide --accessions, --taxid, --clade, or --preset to select assemblies.")

        _apply_busco_context_from_args(manager, args)
        selectors = _selector_request_from_args(args, profile="view_assemblies", manager=manager)
        selected = resolve_selector_accessions(
            manager,
            selectors,
            allow_all=False,
            require_candidates=True,
            use_rule_selection=False,
        )
        species_by_taxid, accessions_by_taxid, missing_accessions = _fetch_species_taxa(manager, selected)
        if missing_accessions and len(missing_accessions) == len(selected):
            return _print_error("No taxonomy information found for the selected accessions.")

        full_tree = _build_taxonomic_tree(manager, species_by_taxid)
        annotated_ranks = normalize_rank_list(getattr(args, "colour_by_ranks", None) or [])
        compressed_root_nodes: List[_DisplayNode] = []
        for child in sorted(full_tree.children.values(), key=_node_sort_key):
            compressed_root_nodes.extend(_compress_tree(child, set(annotated_ranks), force=True))

        if not compressed_root_nodes:
            return _print_error("No taxonomic tree could be built from the selected accessions.")

        _apply_tip_labels(
            compressed_root_nodes,
            species_by_taxid=species_by_taxid,
            accessions_by_taxid=accessions_by_taxid,
            full_name=bool(getattr(args, "full_name", False)),
            show_accession=bool(getattr(args, "show_accession", False)),
        )

        ambiguous_busco_taxa = 0
        if getattr(args, "busco", False):
            busco_annotations, ambiguous_busco_taxa = _resolve_busco_tip_annotations(
                manager,
                args,
                accessions_by_taxid=accessions_by_taxid,
            )
            for node in _iter_nodes(compressed_root_nodes):
                if node.tip_label and node.taxid is not None:
                    suffix = busco_annotations.get(int(node.taxid))
                    if suffix:
                        node.tip_label = f"{node.tip_label} {suffix}"

        colour_map = _assign_rank_colours(compressed_root_nodes, annotated_ranks)
        _render_tree(
            compressed_root_nodes,
            annotated_ranks=annotated_ranks,
            colour_map=colour_map,
            no_colour=bool(getattr(args, "no_colour", False)),
        )

        if args.output:
            _write_newick(str(args.output), full_tree)
        if missing_accessions:
            print(
                f"Skipped {len(missing_accessions)} accession(s) with no taxonomy information.",
                file=sys.stderr,
            )
        if ambiguous_busco_taxa:
            print(
                f"Skipped BUSCO annotation for {ambiguous_busco_taxa} collapsed taxon/taxa representing multiple selected accessions.",
                file=sys.stderr,
            )
        return 0
    except ValueError as exc:
        return _print_error(str(exc))
    except OSError as exc:
        return _print_error(f"Failed to write Newick output: {exc}")
    finally:
        manager.close()


def register_tree_parser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    *,
    selector_defaults: Optional[Mapping[str, Any]] = None,
) -> argparse.ArgumentParser:
    """Register the top-level ``tree`` command."""

    tree_parser = subparsers.add_parser("tree", help="Render a simple taxonomic tree for selected assemblies.")
    selector_group = tree_parser.add_argument_group("Selector options")
    output_group = tree_parser.add_argument_group("Output options")
    _add_selector_arguments(
        selector_group,
        profile="view_assemblies",
        selector_defaults=selector_defaults,
        context_label="tree selection",
    )
    output_group.add_argument(
        "--colour-by-ranks",
        action=AppendCommaSeparated,
        default=None,
        help="Colour represented clades for the supplied ranks and print a legend.",
    )
    output_group.add_argument(
        "-b",
        "--busco",
        action="store_true",
        help="Append BUSCO single-copy percentages to terminal leaf labels.",
    )
    output_group.add_argument(
        "--show-accession",
        action="store_true",
        help="Append accession ids to terminal leaf labels.",
    )
    output_group.add_argument(
        "--full-name",
        action="store_true",
        help="Show full species names in terminal leaf labels instead of abbreviated names.",
    )
    output_group.add_argument("-o", "--output", help="Write the rendered taxonomy as a Newick file.")
    output_group.add_argument("--no-colour", action="store_true", dest="no_colour", help="Disable colour in the console tree output.")
    tree_parser.set_defaults(handler=_handle_tree)
    return tree_parser


__all__ = ["_handle_tree", "register_tree_parser"]
