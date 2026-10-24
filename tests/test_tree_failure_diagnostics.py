import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from phyloODB.tasks.analysis.trees import (
    ensure_shared_output_permissions,
    run_iqtree_analysis,
    run_mafft_alignment,
)
from phyloODB.tasks.task import Task


class ExternalTask:
    REQUIRED_THREADS = 2

    def __init__(self, executable: str, *, force_restart: bool = False):
        self.db_manager = SimpleNamespace(env=SimpleNamespace(get=lambda _key: executable))
        self.force_restart = force_restart

    def payload_bool(self, _key, _default=False):
        return self.force_restart

    def log(self, *_args, **_kwargs):
        return None


def test_mafft_failure_reports_status_command_and_stdout_fallback(tmp_path: Path):
    source = tmp_path / "OG0001.faa"
    source.write_text(">a\nAAAA\n>b\nAAAA\n", encoding="utf-8")

    failed = SimpleNamespace(returncode=23, stderr="", stdout="fatal alignment diagnostic")
    with patch("phyloODB.tasks.analysis.trees.subprocess.run", return_value=failed):
        with pytest.raises(RuntimeError) as caught:
            run_mafft_alignment(
                ExternalTask("/bin/false"),
                input_fasta=str(source),
                out_dir=str(tmp_path / "alignments"),
            )

    message = str(caught.value)
    assert "MAFFT exited with status 23" in message
    assert "Command: /bin/false" in message
    assert "stdout:\nfatal alignment diagnostic" in message


def test_iqtree_failure_includes_log_when_stderr_is_empty(tmp_path: Path):
    alignment = tmp_path / "OG0002.aln.faa"
    alignment.write_text(">a\nAAAA\n>b\nAAAC\n>c\nAAAG\n>d\nAAAT\n", encoding="utf-8")

    def fail_with_log(command, **_kwargs):
        prefix = Path(command[command.index("--prefix") + 1])
        prefix.with_suffix(".log").write_text("model selection failed", encoding="utf-8")
        return SimpleNamespace(returncode=2, stderr="", stdout="")

    with patch("phyloODB.tasks.analysis.trees.subprocess.run", side_effect=fail_with_log):
        with pytest.raises(RuntimeError) as caught:
            run_iqtree_analysis(
                ExternalTask("/bin/false"),
                input_alignment=str(alignment),
                out_dir=str(tmp_path / "trees"),
                prefix="OG0002",
            )

    message = str(caught.value)
    assert "IQ-TREE exited with status 2" in message
    assert "model selection failed" in message
    assert "OG0002.log" in message


def test_shared_output_permissions_add_group_read_and_write(tmp_path: Path):
    output = tmp_path / "alignment.fa"
    output.write_text(">a\nAAAA\n", encoding="utf-8")
    output.chmod(0o600)

    ensure_shared_output_permissions(str(output))

    assert output.stat().st_mode & 0o060 == 0o060


class ParentTask(Task):
    def run(self):
        return True


class TaskRepo:
    def __init__(self, children, errors):
        self.children = children
        self.errors = errors
        self.recorded = []

    def get_subtasks(self, _parent_id):
        return self.children

    def get_error_info(self, task_id):
        return self.errors.get(task_id)

    def set_error(self, task_id, message, stack):
        self.recorded.append((task_id, message, stack))


def _child(task_id, generation, input_name):
    return (
        task_id,
        32,
        "E",
        1,
        99,
        None,
        json.dumps({"__stage": 4, "__gen": generation, "input_fasta": input_name}),
    )


def test_subtask_phase_aggregates_diagnostics_across_retry_generations():
    repo = TaskRepo(
        [_child(10, 1, "OG0001.faa"), _child(11, 2, "OG0002.faa"), _child(12, 2, "OG0001.faa")],
        {
            10: ("MAFFT failed for OG0001: superseded diagnostic", "old stack"),
            11: ("MAFFT failed for OG0002: out of memory", "stack two"),
            12: ("MAFFT failed for OG0001: bad residue", "stack one"),
        },
    )
    parent = ParentTask.__new__(ParentTask)
    parent.task_id = 99
    parent.stage = 4
    parent.status = "R"
    parent.data = {"_stage_4_gen": 2, "_retries_stage_4": 1}
    parent.db_manager = SimpleNamespace(tasks=repo)
    parent.log = lambda *_args, **_kwargs: None
    parent.update_status = lambda status: setattr(parent, "status", status)

    outcome = parent.manage_subtasks(
        stage=4,
        queue_fn=lambda: False,
        done_fn=lambda: False,
        max_retries=1,
    )

    assert outcome == "ERROR"
    _, summary, stacks = repo.recorded[-1]
    assert "failed after 2 attempt(s)" in summary
    assert "OG0001: bad residue" in summary
    assert "OG0002: out of memory" in summary
    assert "superseded diagnostic" not in summary
    assert "stack one" in stacks
    assert "stack two" in stacks
    assert "old stack" not in stacks


def test_incomplete_phase_uses_output_specific_diagnostic():
    complete_child = (20, 32, "C", 1, 99, None, json.dumps({"__stage": 4, "__gen": 2}))
    repo = TaskRepo([complete_child], {})
    parent = ParentTask.__new__(ParentTask)
    parent.task_id = 99
    parent.stage = 4
    parent.status = "R"
    parent.data = {"_stage_4_gen": 2, "_retries_stage_4": 1}
    parent.db_manager = SimpleNamespace(tasks=repo)
    parent.log = lambda *_args, **_kwargs: None
    parent.update_status = lambda status: setattr(parent, "status", status)

    outcome = parent.manage_subtasks(
        stage=4,
        queue_fn=lambda: False,
        done_fn=lambda: False,
        max_retries=1,
        incomplete_message_fn=lambda: "MAFFT phase incomplete: OG0001 alignment is invalid",
        retry_incomplete=True,
    )

    assert outcome == "ERROR"
    assert repo.recorded[-1][1] == "MAFFT phase incomplete: OG0001 alignment is invalid"
