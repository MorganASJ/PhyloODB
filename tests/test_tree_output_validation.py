from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from phyloODB.tasks.analysis.trees import (
    run_iqtree_analysis,
    valid_iqtree_tree,
    valid_mafft_alignment,
)


def test_valid_mafft_alignment_requires_complete_matching_records(tmp_path: Path):
    source = tmp_path / "source.faa"
    source.write_text(">taxon_a\nAAAA\n>taxon_b\nAAA\n", encoding="utf-8")

    alignment = tmp_path / "alignment.faa"
    alignment.write_text(">taxon_a\nAAAA\n>taxon_b\nAAA-\n", encoding="utf-8")
    assert valid_mafft_alignment(str(alignment), str(source))

    alignment.write_text(">taxon_a\nAAAA\n>taxon_b\nAA\n", encoding="utf-8")
    assert not valid_mafft_alignment(str(alignment), str(source))

    alignment.write_text(">taxon_a\nAAAA\n", encoding="utf-8")
    assert not valid_mafft_alignment(str(alignment), str(source))


def test_valid_iqtree_tree_requires_complete_parseable_newick(tmp_path: Path):
    tree = tmp_path / "OG0001.treefile"
    tree.write_text("(taxon_a:0.1,taxon_b:0.2);\n", encoding="utf-8")
    assert valid_iqtree_tree(str(tree))

    tree.write_text("(taxon_a:0.1,taxon_b:0.2)\n", encoding="utf-8")
    assert not valid_iqtree_tree(str(tree))

    tree.write_text("not a tree;\n", encoding="utf-8")
    assert not valid_iqtree_tree(str(tree))

    nexus = tmp_path / "OG0001.nex"
    nexus.write_text(
        "#NEXUS\nBegin trees;\nTree tree1 = (taxon_a:0.1,taxon_b:0.2);\nEnd;\n",
        encoding="utf-8",
    )
    assert valid_iqtree_tree(str(nexus))


def test_iqtree_clean_retry_adds_redo_but_normal_resume_does_not(tmp_path: Path):
    alignment = tmp_path / "alignment.faa"
    alignment.write_text(">taxon_a\nAAAA\n>taxon_b\nAAA-\n>taxon_c\nAAAC\n>taxon_d\nAAAG\n", encoding="utf-8")

    class FakeTask:
        REQUIRED_THREADS = 2
        db_manager = SimpleNamespace(env=SimpleNamespace(get=lambda _key: "/bin/true"))

        def __init__(self, force_restart: bool):
            self.force_restart = force_restart

        def payload_bool(self, _key, _default=False):
            return self.force_restart

        def log(self, *_args, **_kwargs):
            return None

    commands = []

    def completed_iqtree(command, **_kwargs):
        commands.append(command)
        prefix = Path(command[command.index("--prefix") + 1])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        prefix.with_suffix(".treefile").write_text(
            "(taxon_a:0.1,taxon_b:0.2);\n",
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0, stderr="")

    with patch("phyloODB.tasks.analysis.trees.subprocess.run", side_effect=completed_iqtree):
        run_iqtree_analysis(
            FakeTask(False),
            input_alignment=str(alignment),
            out_dir=str(tmp_path / "normal"),
            prefix="normal",
        )
        run_iqtree_analysis(
            FakeTask(True),
            input_alignment=str(alignment),
            out_dir=str(tmp_path / "retry"),
            prefix="retry",
        )

    assert "-redo" not in commands[0]
    assert "-redo" in commands[1]
