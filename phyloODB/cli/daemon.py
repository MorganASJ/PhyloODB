from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..database import DBManager
from ..logging_utils import resolve_log_file_from_db
from ..task_daemon import TaskDaemon
from ..thread_defaults import detect_available_threads, validate_thread_cap

PID_SUFFIX = ".daemon.pid"
DEFAULT_STOP_TIMEOUT = 10.0
DEFAULT_POLL_INTERVAL = 2.0


def build_parser() -> argparse.ArgumentParser:
    description = textwrap.dedent(
        """
        Control the PhyloODB task daemon.

        The daemon monitors the Tasks table, enforcing parent/child constraints and
        launching registered jobs as resources allow. Use this CLI to start or stop the
        overseer process with custom logging and runtime behaviour.
        """
    ).strip()
    parser = argparse.ArgumentParser(
        prog="phyloODB-daemon",
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("database", help="Path to the PhyloODB SQLite database")

    subparsers = parser.add_subparsers(dest="command", required=True)

    start = subparsers.add_parser(
        "start",
        help="Launch the task daemon",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    location = start.add_mutually_exclusive_group()
    location.add_argument(
        "--background",
        dest="background",
        action="store_true",
        help="Run daemon in the background and report the PID",
    )
    location.add_argument(
        "--here",
        dest="background",
        action="store_false",
        help="Run daemon in the foreground, streaming logs to stdout",
    )
    start.set_defaults(background=True)
    start.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="Automatically request a graceful stop after N seconds",
    )
    start.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop the daemon if a new task enters the error state",
    )
    start.add_argument(
        "--max-threads",
        dest="max_threads",
        type=int,
        default=None,
        help="Maximum concurrent worker threads the daemon may consume",
    )
    start.add_argument(
        "--max-concurrent",
        dest="max_threads",
        type=int,
        help="Alias for --max-threads",
    )
    start.add_argument(
        "--threads",
        dest="max_threads",
        type=int,
        help="Alias for --max-threads",
    )
    start.add_argument(
        "--polling",
        type=float,
        default=None,
        help="Override DAEMON_PROCESS_POLLING_TIME (seconds)",
    )
    start.add_argument(
        "--blocked-polling",
        dest="blocked_polling",
        type=float,
        default=None,
        help="Override BLOCKED_TASK_QUEUE_POLLING_TIME (seconds)",
    )
    # Logging controls
    start.add_argument("--logformat", dest="log_format", help="Override logging format string", default=None)
    start.add_argument("--logfile", dest="log_file", help="Write logs to the provided path", default=None)
    start.add_argument("--log-level", dest="log_level", help="Logging level (e.g. INFO, DEBUG)", default=None)
    start.add_argument(
        "--log-max-bytes",
        dest="log_max_bytes",
        type=int,
        default=None,
        help="Maximum size of log file before rotation",
    )
    start.add_argument(
        "--log-backups",
        dest="log_backups",
        type=int,
        default=None,
        help="Number of rotated log files to keep",
    )
    console_group = start.add_mutually_exclusive_group()
    console_group.add_argument(
        "--log-console",
        dest="log_console",
        action="store_true",
        help="Force console logging on",
    )
    console_group.add_argument(
        "--no-log-console",
        dest="log_console",
        action="store_false",
        help="Disable console logging",
    )
    color_group = start.add_mutually_exclusive_group()
    color_group.add_argument(
        "--color",
        dest="log_color",
        action="store_true",
        help="Force ANSI colour output in console logs",
    )
    color_group.add_argument(
        "--no-color",
        dest="log_color",
        action="store_false",
        help="Disable ANSI colour output",
    )
    start.add_argument(
        "--internal",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    stop = subparsers.add_parser(
        "stop",
        help="Stop the running daemon",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stop.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_STOP_TIMEOUT,
        help="Seconds to wait for a graceful stop",
    )
    stop.add_argument(
        "--force",
        action="store_true",
        help="Send SIGKILL if the daemon does not exit in time",
    )
    stop.add_argument(
        "--drain",
        dest="drain",
        action="store_true",
        help="Stop scheduling new tasks and exit after current tasks finish (safe stop).",
    )
    stop.add_argument(
        "--safe",
        dest="drain",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "start":
        return _handle_start(args)
    if args.command == "stop":
        return _handle_stop(args)
    parser.print_help()
    return 1


def _handle_start(args: argparse.Namespace) -> int:
    db_path = _normalize_db_path(args.database)
    if not Path(db_path).exists():
        print(f"Database not found: {db_path}")
        return 1

    existing = _read_pidfile(db_path)
    current_pid = os.getpid()
    if existing and existing != current_pid and _is_pid_alive(existing):
        print(f"Daemon already running (pid {existing}).")
        return 0
    if existing and not _is_pid_alive(existing):
        _cleanup_pidfile(Path(f"{db_path}{PID_SUFFIX}"), expected_pid=existing)

    env_overrides: Dict[str, Any] = {}
    if args.polling is not None:
        env_overrides['DAEMON_PROCESS_POLLING_TIME'] = args.polling
    if args.blocked_polling is not None:
        env_overrides['BLOCKED_TASK_QUEUE_POLLING_TIME'] = args.blocked_polling
    log_overrides = _build_log_overrides(args)
    max_threads = args.max_threads if args.max_threads is not None else None
    if max_threads is not None:
        try:
            validate_thread_cap(max_threads, detect_available_threads(), source="Explicit daemon thread limit")
        except ValueError as exc:
            print(f"Error: {exc}")
            return 1
    timeout = args.timeout
    stop_on_error = bool(args.stop_on_error)

    if args.background and not args.internal:
        return _start_background(
            db_path,
            env_overrides=env_overrides,
            log_overrides=log_overrides,
            max_threads=max_threads,
            timeout=timeout,
            stop_on_error=stop_on_error,
        )

    return _start_foreground(
        db_path,
        env_overrides=env_overrides,
        log_overrides=log_overrides,
        max_threads=max_threads,
        timeout=timeout,
        stop_on_error=stop_on_error,
    )


def _start_foreground(
    db_path: str,
    *,
    env_overrides: Dict[str, Any],
    log_overrides: Dict[str, Any],
    max_threads: Optional[int],
    timeout: Optional[float],
    stop_on_error: bool,
) -> int:
    pid_path = Path(f"{db_path}{PID_SUFFIX}")
    applied_log_overrides = dict(log_overrides)
    applied_log_overrides.setdefault('LOG_TO_CONSOLE', True)
    try:
        daemon = TaskDaemon(
            db_path,
            data=None,
            max_threads=max_threads,
            polling_time=env_overrides.get('DAEMON_PROCESS_POLLING_TIME'),
            env_overrides=env_overrides,
            log_overrides=applied_log_overrides,
            stop_on_error=stop_on_error,
            stop_after=timeout,
        )
    except ValueError as exc:
        print(f"Error: {exc}")
        return 1

    _write_pidfile(pid_path, os.getpid())

    def _shutdown_handler(signum, _frame):  # pragma: no cover - invoked via signal
        print("\nReceived termination signal; stopping daemon…")
        daemon.stop()
    def _drain_handler(signum, _frame):  # pragma: no cover - invoked via signal
        print("\nReceived drain signal; draining daemon…")
        daemon.drain()

    previous_sigterm = None
    previous_sigusr1 = None
    try:
        previous_sigterm = signal.getsignal(signal.SIGTERM)
        signal.signal(signal.SIGTERM, _shutdown_handler)
        if hasattr(signal, "SIGUSR1"):
            previous_sigusr1 = signal.getsignal(signal.SIGUSR1)
            signal.signal(signal.SIGUSR1, _drain_handler)
    except (AttributeError, ValueError):
        previous_sigterm = None

    print("Starting daemon in foreground (Ctrl-C to stop)…")
    try:
        daemon.start()
    except KeyboardInterrupt:
        print("\nStopping daemon…")
        daemon.stop()
    finally:
        if not daemon.stop_event.is_set():
            daemon.stop()
        try:
            if previous_sigterm is not None:
                signal.signal(signal.SIGTERM, previous_sigterm)
            if previous_sigusr1 is not None and hasattr(signal, "SIGUSR1"):
                signal.signal(signal.SIGUSR1, previous_sigusr1)
        except ValueError:
            pass
        _cleanup_pidfile(pid_path, os.getpid())
    return 0


def _start_background(
    db_path: str,
    *,
    env_overrides: Dict[str, Any],
    log_overrides: Dict[str, Any],
    max_threads: Optional[int],
    timeout: Optional[float],
    stop_on_error: bool,
) -> int:
    child_args = [
        sys.executable,
        "-m",
        "phyloODB.cli.daemon",
        db_path,
        "start",
        "--here",
        "--internal",
    ]
    if timeout is not None:
        child_args += ["--timeout", str(timeout)]
    if stop_on_error:
        child_args.append("--stop-on-error")
    if max_threads is not None:
        child_args += ["--max-threads", str(max_threads)]
    if (poll := env_overrides.get('DAEMON_PROCESS_POLLING_TIME')) is not None:
        child_args += ["--polling", str(poll)]
    if (blocked := env_overrides.get('BLOCKED_TASK_QUEUE_POLLING_TIME')) is not None:
        child_args += ["--blocked-polling", str(blocked)]
    for flag, key in [
        ("--logformat", 'LOG_FORMAT'),
        ("--logfile", 'LOG_FILE'),
        ("--log-level", 'LOG_LEVEL'),
        ("--log-max-bytes", 'LOG_MAX_BYTES'),
        ("--log-backups", 'LOG_BACKUPS'),
    ]:
        if key in log_overrides and log_overrides[key] is not None:
            child_args += [flag, str(log_overrides[key])]
    if 'LOG_TO_CONSOLE' in log_overrides:
        child_args.append("--log-console" if log_overrides['LOG_TO_CONSOLE'] else "--no-log-console")
    if 'LOG_USE_COLOR' in log_overrides:
        child_args.append("--color" if log_overrides['LOG_USE_COLOR'] else "--no-color")

    log_path = _determine_log_path(db_path, log_overrides)
    log_handle = None
    if log_path:
        os.makedirs(log_path.parent, exist_ok=True)
        log_handle = open(log_path, "a", buffering=1)

    proc = subprocess.Popen(
        child_args,
        stdout=log_handle or subprocess.DEVNULL,
        stderr=log_handle or subprocess.DEVNULL,
        stdin=subprocess.DEVNULL,
        close_fds=True,
        start_new_session=True,
    )
    if log_handle:
        log_handle.close()
    _write_pidfile(Path(f"{db_path}{PID_SUFFIX}"), proc.pid)
    print(f"Daemon started in background (pid {proc.pid}).")
    if log_path:
        print(f"Logs: {log_path}")
    return 0


def _handle_stop(args: argparse.Namespace) -> int:
    db_path = _normalize_db_path(args.database)
    pid = _read_pidfile(db_path)
    if not pid:
        print("No daemon pidfile found.")
        return 1
    if not _is_pid_alive(pid):
        _cleanup_pidfile(Path(f"{db_path}{PID_SUFFIX}"))
        print("Daemon not running; cleaned stale pidfile.")
        return 0

    sig = signal.SIGUSR1 if args.drain and hasattr(signal, "SIGUSR1") else signal.SIGTERM
    try:
        os.kill(pid, sig)
    except ProcessLookupError:
        _cleanup_pidfile(Path(f"{db_path}{PID_SUFFIX}"))
        print("Daemon already exited.")
        return 0
    except OSError as exc:
        print(f"Failed to signal daemon: {exc}")
        return 1

    deadline = time.time() + (args.timeout or DEFAULT_STOP_TIMEOUT)
    while time.time() < deadline:
        if not _is_pid_alive(pid):
            _cleanup_pidfile(Path(f"{db_path}{PID_SUFFIX}"))
            print("Daemon stopped.")
            return 0
        time.sleep(0.2)

    if args.force and _is_pid_alive(pid):
        try:
            os.kill(pid, signal.SIGKILL)
        except OSError as exc:
            print(f"Failed to force-kill daemon: {exc}")
            return 1
        _cleanup_pidfile(Path(f"{db_path}{PID_SUFFIX}"))
        print("Daemon force-stopped.")
        return 0

    if _is_pid_alive(pid):
        print("Daemon did not exit within timeout. Use --force to terminate.")
        return 1

    _cleanup_pidfile(Path(f"{db_path}{PID_SUFFIX}"))
    print("Daemon stopped.")
    return 0


# ----- helpers -----

def _normalize_db_path(db: str) -> str:
    return os.path.abspath(db)


def _build_log_overrides(args: argparse.Namespace) -> Dict[str, Any]:
    overrides: Dict[str, Any] = {}
    if args.log_format is not None:
        overrides['LOG_FORMAT'] = args.log_format
    if args.log_file is not None:
        overrides['LOG_FILE'] = args.log_file
    if args.log_level is not None:
        overrides['LOG_LEVEL'] = args.log_level
    if args.log_max_bytes is not None:
        overrides['LOG_MAX_BYTES'] = args.log_max_bytes
    if args.log_backups is not None:
        overrides['LOG_BACKUPS'] = args.log_backups
    if args.log_console is not None:
        overrides['LOG_TO_CONSOLE'] = bool(args.log_console)
    if args.log_color is not None:
        overrides['LOG_USE_COLOR'] = bool(args.log_color)
    return overrides


def _determine_log_path(db_path: str, log_overrides: Dict[str, Any]) -> Optional[Path]:
    manager = DBManager(db_path)
    try:
        manager.connect()
        log_file = resolve_log_file_from_db(manager, overrides=log_overrides)
    finally:
        manager.close()
    if log_file:
        return Path(log_file).expanduser().resolve()
    return None


def _write_pidfile(pid_path: Path, pid: int) -> None:
    pid_path.write_text(str(pid), encoding="utf-8")


def _read_pidfile(db_path: str) -> Optional[int]:
    pid_path = Path(f"{db_path}{PID_SUFFIX}")
    if not pid_path.exists():
        return None
    try:
        return int(pid_path.read_text(encoding="utf-8").strip())
    except (ValueError, OSError):
        return None


def _cleanup_pidfile(pid_path: Path, expected_pid: Optional[int] = None) -> None:
    try:
        if pid_path.exists():
            if expected_pid is not None:
                try:
                    contents = pid_path.read_text(encoding="utf-8").strip()
                    if contents and int(contents) != expected_pid:
                        return
                except (ValueError, OSError):
                    pass
            pid_path.unlink()
    except OSError:
        pass


def _is_pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    proc_stat = Path(f"/proc/{pid}/stat")
    try:
        status_fields = proc_stat.read_text(encoding="utf-8").split()
    except OSError:
        return True
    if len(status_fields) > 2 and status_fields[2] == "Z":
        return False
    return True


if __name__ == "__main__":
    raise SystemExit(main())
