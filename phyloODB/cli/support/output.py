"""Shared CLI output helpers for table and list rendering."""
from __future__ import annotations

import argparse
import math
import queue
import select
import sys
import threading
import textwrap
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from .common import _coerce_bool, _connect_manager, _load_list_color_defaults


@dataclass(frozen=True)
class TableColumn:
    """Display hints for one pretty-table column."""

    heading: str
    priority: int = 0
    min_width: int = 4
    max_width: int = 60
    no_wrap: bool = False
    style: Optional[str] = None


@dataclass
class TableData:
    """Renderer-neutral table content."""

    headers: Sequence[str]
    rows: Sequence[Sequence[Any]]
    row_styles: Optional[Sequence[Any]] = None
    columns: Optional[Sequence[TableColumn]] = None

def _format_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    *,
    show_header: bool = True,
) -> str:
    """Render rows as a wrapped plain-text table for terminal output."""

    if not headers:
        return ""
    max_col_width = 60
    widths = []
    for i, header in enumerate(headers):
        column_values = [str(row[i]) for row in rows] if rows else []
        longest = max(len(str(header)), *(len(value) for value in column_values)) if column_values else len(str(header))
        widths.append(min(longest, max_col_width))

    def _wrap_cell(text: str, width: int) -> List[str]:
        return textwrap.fill(text, width=width, subsequent_indent="", replace_whitespace=False).split("\n") if text else [""]

    rendered = []
    if show_header:
        header_lines = [_wrap_cell(str(header), widths[idx]) for idx, header in enumerate(headers)]
        header_height = max(len(lines) for lines in header_lines)
        for line_idx in range(header_height):
            parts = []
            for col_idx, lines in enumerate(header_lines):
                cell_line = lines[line_idx] if line_idx < len(lines) else ""
                parts.append(cell_line.ljust(widths[col_idx]))
            rendered.append("  ".join(parts))
        divider = "  ".join("-" * widths[idx] for idx in range(len(headers)))
        rendered.append(divider)

    for row in rows:
        wrapped_cols = [_wrap_cell(str(row[idx]), widths[idx]) for idx in range(len(headers))]
        height = max(len(col) for col in wrapped_cols)
        for line_idx in range(height):
            parts = []
            for col_idx, lines in enumerate(wrapped_cols):
                cell_line = lines[line_idx] if line_idx < len(lines) else ""
                parts.append(cell_line.ljust(widths[col_idx]))
            rendered.append("  ".join(parts))

    return "\n".join(rendered)


def _format_tsv(headers: Sequence[str], rows: Sequence[Sequence[Any]], *, show_header: bool = True) -> str:
    """Render rows without wrapping so paths remain exact."""

    rendered = ["\t".join(str(header) for header in headers)] if show_header else []
    for row in rows:
        rendered.append("\t".join("" if value is None else str(value) for value in row))
    return "\n".join(rendered)


def _column_specs(data: TableData) -> List[TableColumn]:
    if data.columns:
        return list(data.columns)
    count = len(data.headers)
    return [
        TableColumn(str(header), priority=count - idx, min_width=min(max(len(str(header)), 4), 12))
        for idx, header in enumerate(data.headers)
    ]


def _column_layout(data: TableData, width: int) -> tuple[List[int], int, List[int]]:
    specs = _column_specs(data)
    desired: List[int] = []
    for idx, spec in enumerate(specs):
        values = [str(row[idx] if idx < len(row) else "") for row in data.rows[:200]]
        longest = max([len(spec.heading), *[len(value) for value in values]])
        desired.append(min(max(longest, spec.min_width), spec.max_width))
    visible = list(range(len(specs)))
    # Rich needs separators/padding in addition to cell contents.
    while len(visible) > 1 and sum(desired[idx] + 3 for idx in visible) > max(width - 2, 20):
        removable = min(visible[1:] or visible, key=lambda idx: (specs[idx].priority, -idx))
        visible.remove(removable)
    return visible, len(specs) - len(visible), desired


def _visible_column_indices(data: TableData, width: int) -> tuple[List[int], int]:
    visible, hidden, _widths = _column_layout(data, width)
    return visible, hidden


def _rich_table(
    data: TableData,
    *,
    width: int,
    page: int,
    page_size: int,
    interactive: bool,
    show_header: bool = True,
):
    try:
        from rich.table import Table
        from rich.text import Text
    except ImportError as exc:
        raise PhyloODBError("The 'rich' package is required for colored list output.") from exc

    specs = _column_specs(data)
    visible, hidden, column_widths = _column_layout(data, width)
    pages: List[List[int]] = [[]]
    used_lines = 0
    for row_index, row in enumerate(data.rows):
        row_height = max(
            1,
            max(
                math.ceil(len(str(row[idx] if idx < len(row) else "")) / max(column_widths[idx], 1))
                for idx in visible
            ),
        )
        if pages[-1] and used_lines + row_height > max(page_size, 1):
            pages.append([])
            used_lines = 0
        pages[-1].append(row_index)
        used_lines += min(row_height, max(page_size, 1))
    total_pages = max(1, len(pages))
    page = max(1, min(page, total_pages))
    page_indices = pages[page - 1] if pages else []
    styles = list(data.row_styles or [])

    table = Table(show_header=show_header, header_style="bold magenta", expand=False)
    for idx in visible:
        spec = specs[idx]
        table.add_column(
            spec.heading,
            width=column_widths[idx],
            min_width=spec.min_width,
            max_width=spec.max_width,
            no_wrap=spec.no_wrap,
            overflow="ellipsis" if spec.no_wrap else "fold",
            style=spec.style,
        )
    for absolute_idx in page_indices:
        row = data.rows[absolute_idx]
        style = styles[absolute_idx] if absolute_idx < len(styles) else None
        table.add_row(*["" if row[idx] is None else row[idx] for idx in visible], style=style)
    if interactive:
        hidden_note = f" • {hidden} columns hidden" if hidden else ""
        table.caption = Text(
            f"Page {page}/{total_pages} • [/], ←/→, or PgUp/PgDn • Home/End • q to quit{hidden_note}"
        )
    return table, total_pages


def _start_pager_input(actions_queue: queue.SimpleQueue[str]) -> tuple[threading.Event, threading.Thread]:
    import termios
    import tty

    fd = sys.stdin.fileno()
    original = termios.tcgetattr(fd)
    stop = threading.Event()
    tty.setcbreak(fd)

    def restore() -> None:
        try:
            termios.tcsetattr(fd, termios.TCSADRAIN, original)
        except (OSError, termios.error):
            pass

    def monitor() -> None:
        sequence = ""
        try:
            while not stop.is_set():
                ready, _, _ = select.select([sys.stdin], [], [], 0.1)
                if not ready:
                    continue
                char = sys.stdin.read(1)
                if sequence or char == "\x1b":
                    sequence += char
                    key_actions = {
                        "\x1b[C": "next",
                        "\x1b[D": "previous",
                        "\x1bOC": "next",
                        "\x1bOD": "previous",
                        "\x1b[5~": "previous",
                        "\x1b[6~": "next",
                        "\x1b[H": "first",
                        "\x1b[F": "last",
                        "\x1b[1~": "first",
                        "\x1b[4~": "last",
                    }
                    if sequence in key_actions:
                        actions_queue.put(key_actions[sequence])
                        sequence = ""
                    elif len(sequence) >= 4 and not any(key.startswith(sequence) for key in key_actions):
                        sequence = ""
                    continue
                if char == "[":
                    actions_queue.put("previous")
                    continue
                if char == "]":
                    actions_queue.put("next")
                    continue
                if char.lower() == "q":
                    stop.set()
        finally:
            restore()

    worker = threading.Thread(target=monitor, name="PhyloODBListPager", daemon=True)
    worker.start()
    return stop, worker


def _render_pretty(
    data_or_provider: TableData | Callable[[], TableData],
    *,
    watch: bool,
    refresh: float,
    paginate: bool = True,
    show_header: bool = True,
) -> None:
    try:
        from rich.console import Console
        from rich.live import Live
    except ImportError as exc:
        raise PhyloODBError("The 'rich' package is required for colored list output.") from exc

    console = Console()
    provider = data_or_provider if callable(data_or_provider) else lambda: data_or_provider
    initial = provider()
    height = max(console.size.height - 6, 3) if paginate else 1_000_000
    initial_table, initial_pages = _rich_table(
        initial,
        width=console.size.width,
        page=1,
        page_size=height,
        interactive=watch or len(initial.rows) > height,
        show_header=show_header,
    )
    interactive = watch or (paginate and initial_pages > 1)
    if not interactive:
        console.print(initial_table)
        return
    initial_table, initial_pages = _rich_table(
        initial,
        width=console.size.width,
        page=1,
        page_size=height,
        interactive=True,
        show_header=show_header,
    )
    if not sys.stdin.isatty() or not sys.stdout.isatty():
        raise PhyloODBError("--watch and pretty pagination require an interactive terminal.")

    state: Dict[str, Any] = {"page": 1}
    actions_queue: queue.SimpleQueue[str] = queue.SimpleQueue()
    listener = _start_pager_input(actions_queue)
    last_refresh = time.monotonic()
    current = initial
    try:
        with Live(initial_table, console=console, refresh_per_second=10, transient=False) as live:
            while not listener[0].is_set():
                now = time.monotonic()
                if watch and now - last_refresh >= refresh:
                    current = provider()
                    last_refresh = now
                page_size = max(console.size.height - 6, 3) if paginate else 1_000_000
                try:
                    action = actions_queue.get_nowait()
                except queue.Empty:
                    action = None
                _, total_pages = _rich_table(
                    current,
                    width=console.size.width,
                    page=int(state["page"]),
                    page_size=page_size,
                    interactive=True,
                    show_header=show_header,
                )
                if action == "next":
                    state["page"] = min(total_pages, int(state["page"]) + 1)
                elif action == "previous":
                    state["page"] = max(1, int(state["page"]) - 1)
                elif action == "first":
                    state["page"] = 1
                elif action == "last":
                    state["page"] = total_pages
                state["page"] = min(int(state["page"]), total_pages)
                table, _ = _rich_table(
                    current,
                    width=console.size.width,
                    page=int(state["page"]),
                    page_size=page_size,
                    interactive=True,
                    show_header=show_header,
                )
                live.update(table)
                time.sleep(0.1)
    except KeyboardInterrupt:
        pass
    finally:
        listener[0].set()
        listener[1].join(timeout=0.5)


def _render_list_output(
    args: argparse.Namespace,
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    *,
    default_tidy: bool,
    color_defaults: Optional[Dict[str, Any]] = None,
    row_styles: Optional[Sequence[Any]] = None,
    columns: Optional[Sequence[TableColumn]] = None,
    watch_provider: Optional[Callable[[], TableData]] = None,
    refresh: float = 2.0,
) -> int:
    """Render list output consistently across simple list subcommands."""

    if color_defaults is None:
        manager = _connect_manager(args.database, read_only=True)
        try:
            color_defaults = _load_list_color_defaults(manager)
        finally:
            manager.close()
    list_color = args.list_color if getattr(args, "list_color", None) is not None else _coerce_bool(color_defaults.get("LIST_USE_COLOR", False))
    watch = bool(getattr(args, "watch", False))
    use_rich = bool(list_color or watch) and sys.stdout.isatty()
    if watch and not use_rich:
        raise PhyloODBError("--watch requires pretty output in an interactive terminal.")
    if use_rich:
        show_header = not bool(getattr(args, "no_header", False))
        data = TableData(headers, rows, row_styles=row_styles, columns=columns)
        _render_pretty(
            watch_provider or data,
            watch=watch,
            refresh=refresh,
            paginate=not bool(getattr(args, "no_pager", False)),
            show_header=show_header,
        )
        return 0
    tidy = bool(getattr(args, "tidy", False) or default_tidy)
    show_header = not bool(getattr(args, "no_header", False))
    rendered = (
        _format_table(headers, rows, show_header=show_header)
        if tidy
        else _format_tsv(headers, rows, show_header=show_header)
    )
    print(rendered)
    return 0

def _render_grouped_rows(
    groups: Sequence[Tuple[Optional[str], Sequence[Sequence[str]]]],
    *,
    headers: Sequence[str],
    tidy: bool,
    show_header: bool = True,
) -> str:
    """Render grouped rows as either tidy columns or raw TSV text."""

    rows = [row for _header, group_rows in groups for row in group_rows]
    if not rows:
        return ""
    col_count = len(headers)
    if tidy:
        widths = [len(str(h)) for h in headers]
        for row in rows:
            for idx, value in enumerate(row):
                widths[idx] = max(widths[idx], len(value))

        def render(row: Sequence[str]) -> str:
            return "  ".join(row[idx].ljust(widths[idx]) for idx in range(col_count)).rstrip()
    else:
        def render(row: Sequence[str]) -> str:
            return "\t".join(row)

    lines: List[str] = []
    if headers and show_header:
        lines.append(render([str(h) for h in headers]))
    for header, group_rows in groups:
        if header:
            lines.append(header)
        for row in group_rows:
            lines.append(render(row))
    return "\n".join(lines)


def _parse_color_tokens(raw: Any, fallback: Sequence[str]) -> List[str]:
    """Parse a colour list from stored environment settings."""

    if raw is None:
        return list(fallback)
    if isinstance(raw, (list, tuple)):
        tokens = [str(v).strip() for v in raw if str(v).strip()]
        return tokens or list(fallback)
    tokens = [t.strip() for t in str(raw).split(",") if t.strip()]
    return tokens or list(fallback)


def _parse_rgb(token: str) -> Optional[Tuple[int, int, int]]:
    """Resolve a named or hexadecimal colour token to an RGB tuple."""

    if not token:
        return None
    t = token.strip().lower()
    if t.startswith("#"):
        hex_val = t[1:]
        if len(hex_val) == 3:
            hex_val = "".join(ch * 2 for ch in hex_val)
        if len(hex_val) == 6:
            try:
                r = int(hex_val[0:2], 16)
                g = int(hex_val[2:4], 16)
                b = int(hex_val[4:6], 16)
                return (r, g, b)
            except ValueError:
                return None
        return None
    named = {
        "red": (220, 38, 38),
        "green": (22, 163, 74),
        "yellow": (234, 179, 8),
        "blue": (37, 99, 235),
        "cyan": (14, 116, 144),
        "bright_cyan": (6, 182, 212),
        "magenta": (192, 38, 211),
        "orange": (249, 115, 22),
        "orange1": (255, 135, 0),
        "grey50": (128, 128, 128),
        "gray50": (128, 128, 128),
        "chartreuse3": (102, 189, 99),
    }
    return named.get(t)


def _parse_float_tokens(raw: Any, fallback: Sequence[float]) -> List[float]:
    """Parse comma-delimited float tokens from environment settings."""

    if raw is None:
        return list(fallback)
    if isinstance(raw, (list, tuple)):
        vals = []
        for v in raw:
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                continue
        return vals or list(fallback)
    tokens = [t.strip() for t in str(raw).split(",") if t.strip()]
    vals = []
    for t in tokens:
        try:
            vals.append(float(t))
        except (TypeError, ValueError):
            continue
    return vals or list(fallback)


def _interpolate_rgb(a: Tuple[int, int, int], b: Tuple[int, int, int], t: float) -> Tuple[int, int, int]:
    """Interpolate between two RGB colours for gradient rendering."""

    t = max(0.0, min(1.0, t))
    return (
        int(a[0] + (b[0] - a[0]) * t),
        int(a[1] + (b[1] - a[1]) * t),
        int(a[2] + (b[2] - a[2]) * t),
    )


def _color_from_gradient(
    value: Optional[float],
    gradient: Sequence[str],
    *,
    max_value: float = 100.0,
    stops: Optional[Sequence[float]] = None,
) -> Optional[str]:
    """Map a numeric value onto a configured terminal colour gradient."""

    if value is None:
        return None
    if not gradient:
        return None
    rgb_stops: List[Tuple[int, int, int]] = []
    for token in gradient:
        rgb = _parse_rgb(token)
        if rgb is not None:
            rgb_stops.append(rgb)
    if not rgb_stops:
        return None
    if len(rgb_stops) == 1:
        r, g, b = rgb_stops[0]
        return f"#{r:02x}{g:02x}{b:02x}"
    if max_value <= 0:
        max_value = 100.0
    v = max(0.0, min(max_value, value))
    if stops is not None and len(stops) == len(gradient):
        numeric_stops = list(stops)
        numeric_stops.sort()
        if v <= numeric_stops[0]:
            rgb = rgb_stops[0]
            return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
        if v >= numeric_stops[-1]:
            rgb = rgb_stops[-1]
            return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
        for i in range(len(numeric_stops) - 1):
            lo = numeric_stops[i]
            hi = numeric_stops[i + 1]
            if v <= hi:
                denom = (hi - lo) if hi != lo else 1.0
                local_t = (v - lo) / denom
                rgb = _interpolate_rgb(rgb_stops[i], rgb_stops[i + 1], local_t)
                return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    t = max(0.0, min(1.0, v / max_value))
    segment = 1.0 / (len(rgb_stops) - 1)
    idx = min(int(t / segment), len(rgb_stops) - 2)
    local_t = (t - idx * segment) / segment if segment > 0 else 0.0
    rgb = _interpolate_rgb(rgb_stops[idx], rgb_stops[idx + 1], local_t)
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def _normalize_percent(value: Any, *, allow_fraction: bool = True) -> Optional[float]:
    """Normalise percentage-like values onto a 0-100 scale."""

    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if val < 0:
        val = 0
    if allow_fraction and val <= 1.0:
        val = val * 100.0
    if val > 100:
        val = 100.0
    return val


def _coerce_float(value: Any, default: float) -> float:
    """Coerce a value to float while preserving a caller-provided default."""

    try:
        if value is None:
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _render_grouped_rows_rich(
    groups: Sequence[Tuple[Optional[str], Sequence[Sequence[str]]]],
    *,
    headers: Sequence[str],
    group_colors: Sequence[str],
    busco_pos_gradient: Sequence[str],
    busco_neg_gradient: Sequence[str],
    busco_pos_stops: Sequence[float],
    busco_neg_stops: Sequence[float],
    busco_steep_stops: Sequence[float],
    busco_pos_indices: Sequence[int],
    busco_neg_indices: Sequence[int],
    busco_steep_indices: Sequence[int],
    busco_steep_max: float,
    rank_color_map: Dict[int, int],
    busco_bold_indices: Sequence[int],
    paginate: bool = True,
    show_header: bool = True,
) -> None:
    """Render grouped assembly rows with Rich styling when colour output is enabled."""

    try:
        from rich.text import Text
        from rich.style import Style
    except ImportError:
        raise PhyloODBError("The 'rich' package is required for colored list output. Install it via 'pip install rich'.")

    rendered_rows: List[Sequence[Any]] = []
    for header, group_rows in groups:
        if header:
            leading = len(header) - len(header.lstrip(" "))
            depth = leading // 2
            color = group_colors[min(depth, len(group_colors) - 1)] if group_colors else "cyan"
            style = Style(color=color, bold=True)
            header_row = [Text(header, style=style)] + [Text("") for _ in range(len(headers) - 1)]
            rendered_rows.append(header_row)
        for row in group_rows:
            cells = []
            for idx, value in enumerate(row):
                # Decision highlighting
                if isinstance(value, str) and value.upper() in {"CONTAMINATED", "UNCERTAIN"} and idx < len(headers):
                    if value.upper() == "CONTAMINATED":
                        cells.append(Text(value, style=Style(color="yellow", bgcolor="red", bold=True)))
                        continue
                    if value.upper() == "UNCERTAIN":
                        cells.append(Text(value, style=Style(color="green", bgcolor="#ffbf00", bold=True)))
                        continue
                if idx in rank_color_map:
                    depth = rank_color_map.get(idx, 0)
                    color = group_colors[min(depth, len(group_colors) - 1)] if group_colors else "cyan"
                    cells.append(Text(str(value), style=Style(color=color, bold=False)))
                    continue
                if idx in busco_pos_indices or idx in busco_neg_indices or idx in busco_steep_indices:
                    pct = _normalize_percent(value, allow_fraction=False)
                    if idx in busco_pos_indices:
                        gradient = busco_pos_gradient
                        max_val = 100.0
                        stops = busco_pos_stops
                    elif idx in busco_steep_indices:
                        gradient = busco_neg_gradient
                        max_val = busco_steep_max
                        stops = busco_steep_stops
                    else:
                        gradient = busco_neg_gradient
                        max_val = 100.0
                        stops = busco_neg_stops
                    color = _color_from_gradient(pct, gradient, max_value=max_val, stops=stops) if pct is not None else None
                    bold = False
                    if idx in busco_bold_indices and pct is not None:
                        bold = pct <= 10.0 or pct >= 90.0
                    if color or bold:
                        cells.append(Text(str(value), style=Style(color=color, bold=bold)))
                    else:
                        cells.append(Text(str(value)))
                else:
                    cells.append(Text(str(value)))
            rendered_rows.append(cells)

    _render_pretty(
        TableData(headers, rendered_rows),
        watch=False,
        refresh=2.0,
        paginate=paginate,
        show_header=show_header,
    )

def _fmt_busco_value(value: Optional[float]) -> str:
    """Format BUSCO values consistently for terminal output."""

    if value is None:
        return "NA"
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return str(value)


__all__ = [
    "_coerce_float",
    "_color_from_gradient",
    "_fmt_busco_value",
    "_format_table",
    "_format_tsv",
    "_interpolate_rgb",
    "_normalize_percent",
    "_parse_color_tokens",
    "_parse_float_tokens",
    "_parse_rgb",
    "TableColumn",
    "TableData",
    "_render_grouped_rows",
    "_render_grouped_rows_rich",
    "_render_list_output",
]
from ...errors import PhyloODBError
