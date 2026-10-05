import json
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from phyloODB.tasks.analysis.trees import (run_iqtree_analysis, tree_omitted, _without_support_flags,
    _tree_row_best, IQTreeTask, _write_tree_outcome, alignment_tree_counts)
from phyloODB.tasks.analysis.add_library import AddLibraryTask


IQTREE = shutil.which("iqtree") or str(Path(sys.executable).with_name("iqtree"))


class ExternalTask:
    REQUIRED_THREADS = 1
    data = {}
    def __init__(self, flags=None):
        self.db_manager = SimpleNamespace(env=SimpleNamespace(get=lambda key: IQTREE if key == 'IQTREE_PATH' else flags))
    def payload_bool(self, *args): return False
    def log(self, *args): pass


def alignment(tmp_path, count):
    path = tmp_path / 'OG.fa'
    sequences = ['ACDEFGHIKLMNPQRSTVWY'*3, 'ACDEFGHIKLMNPQRSTVWA'*3,
                 'ACDEFGHIKLMNPQRSTVWC'*3, 'ACDEFGHIKLMNPQRSTVWD'*3]
    path.write_text(''.join(f'>s{i}\n{seq}\n' for i, seq in enumerate(sequences[:count] + [sequences[0]]*5)))
    return path


@pytest.mark.parametrize('count', [2, 3, 4])
def test_thresholds_and_stale_outputs(tmp_path, count):
    source = alignment(tmp_path, count)
    output = tmp_path/'trees'/'OG'
    output.mkdir(parents=True)
    (output/'OG.treefile').write_text('(stale1,stale2);')
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        (output/'OG.treefile').write_text('(' + ','.join(f's{i}:0.1' for i in range(count+5)) + ');')
        return SimpleNamespace(returncode=0)
    with patch('phyloODB.tasks.analysis.trees.subprocess.run', side_effect=run):
        run_iqtree_analysis(ExternalTask(), input_alignment=str(source), out_dir=str(tmp_path/'trees'), prefix='OG', iqtree_flags='-m LG -B 1000 -alrt 1000')
    outcome = json.loads((output/'OG.outcome.json').read_text())
    assert outcome['input_count'] == count+5
    assert outcome['distinct_count'] == count
    if count < 3:
        assert not commands
        assert tree_omitted(str(output), 'OG', str(source))
        assert not (output/'OG.treefile').exists()
    else:
        assert ('-B' in commands[0]) == (count == 4)
        assert ('-alrt' in commands[0]) == (count == 4)
        assert '-m' in commands[0] and 'LG' in commands[0]
    source.write_text(source.read_text() + '>changed\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
    assert not tree_omitted(str(output), 'OG', str(source))
    assert _tree_row_best({'tree_dir': str(output), 'prefix':'OG', 'alignment_path':str(source)}) is None


def test_support_aliases():
    assert _without_support_flags(['-m', 'LG', '--ufboot=1000', '-b', '100', '-abayes', '-bnni', '-nt', '2']) == ['-m', 'LG', '-nt', '2']


@pytest.mark.skipif(not Path(IQTREE).is_file(), reason="IQ-TREE is not installed")
def test_real_iqtree_duplicate_restoration(tmp_path):
    source = alignment(tmp_path, 3)
    directory, tree, command = run_iqtree_analysis(ExternalTask(), input_alignment=str(source), out_dir=str(tmp_path/'trees'), iqtree_flags='-m LG -B 1000')
    from ete3 import Tree
    assert set(Tree(Path(tree).read_text(), format=1).get_leaf_names()) == {f's{i}' for i in range(8)}
    assert '-B' not in command


def test_all_omitted_families_fail_clearly(tmp_path):
    source = alignment(tmp_path, 2)
    directory = tmp_path/'trees'/'OG'
    _write_tree_outcome(str(directory), 'OG', {**alignment_tree_counts(str(source)), 'outcome':'omitted', 'reason':'insufficient_distinct_sequences'})
    row = {'orthogroup':'OG', 'tree_dir':str(directory), 'prefix':'OG', 'alignment_path':str(source)}
    task = SimpleNamespace(_replacement_tree_rows=lambda pairs:[row], _effective_gene_tree_source=lambda:'iqtree')
    with pytest.raises(ValueError, match='No eligible core families'):
        AddLibraryTask._finalize_replacement_tree_analysis(task, None, [('f','OG')], [('f','OG')], {})


def test_real_failure_cannot_adopt_stale_tree(tmp_path):
    source = alignment(tmp_path, 4)
    directory = tmp_path/'trees'/'OG'
    directory.mkdir(parents=True)
    (directory/'OG.treefile').write_text('(s0,s1,s2,s3);')
    with patch('phyloODB.tasks.analysis.trees.subprocess.run', return_value=SimpleNamespace(returncode=2, stderr='real failure', stdout='')):
        with pytest.raises(RuntimeError, match='real failure'):
            run_iqtree_analysis(ExternalTask(), input_alignment=str(source), out_dir=str(tmp_path/'trees'), prefix='OG')
    row = {'tree_dir':str(directory), 'prefix':'OG', 'alignment_path':str(source)}
    assert _tree_row_best(row) is None
    assert not tree_omitted(str(directory), 'OG', str(source))


def test_omission_completes_parent_phases_without_queueing(tmp_path):
    from phyloODB.tasks.analysis.trees import BuildBuscoTreesTask
    source = alignment(tmp_path, 2)
    directory = tmp_path/'trees'/'OG'
    _write_tree_outcome(str(directory), 'OG', {**alignment_tree_counts(str(source)), 'outcome':'omitted', 'reason':'insufficient_distinct_sequences'})
    row = {'orthogroup':'OG', 'family_id':'f', 'tree_dir':str(directory), 'prefix':'OG', 'alignment_path':str(source)}
    task = SimpleNamespace(_replacement_tree_rows=lambda pairs:[row], _family_rows=lambda:[row],
                           _effective_gene_tree_source=lambda:'iqtree', _rerun_gene_trees_effective=lambda:False,
                           queue_subtask=lambda **kwargs:pytest.fail('omissions must not be queued'))
    assert AddLibraryTask._replacement_iqtree_done(task, [('f','OG')])
    assert not AddLibraryTask._queue_replacement_iqtree_subtasks(task, [('f','OG')])
    assert BuildBuscoTreesTask._iqtree_done(task)
    assert not BuildBuscoTreesTask._queue_iqtree_subtasks(task)


def test_materialization_clears_stale_canonical_and_records_omission(tmp_path):
    import csv
    source = alignment(tmp_path, 2)
    directory = tmp_path/'metadata'/'OG'
    canonical_dir = tmp_path/'canonical'
    canonical_dir.mkdir()
    canonical = canonical_dir/'OG_tree.txt'
    canonical.write_text('(stale1,stale2);')
    annotations = tmp_path/'annotations'
    annotations.mkdir()
    (annotations/'OG.nex').write_text('stale')
    _write_tree_outcome(str(directory), 'OG', {**alignment_tree_counts(str(source)), 'outcome':'omitted', 'reason':'insufficient_distinct_sequences'})
    row = {'orthogroup':'OG', 'family_id':'f', 'tree_dir':str(directory), 'prefix':'OG', 'alignment_path':str(source), 'raw_fasta':str(source), 'canonical_tree':str(canonical)}
    task = SimpleNamespace(_prepare_replacement_tree_workspace=lambda:None, _replacement_tree_rows=lambda pairs:[row],
                           _tree_workspace_paths=lambda:{'iqtree_dir':str(canonical_dir)},
                           _effective_gene_tree_source=lambda:'iqtree', _effective_tree_dir=lambda:str(canonical_dir),
                           _annotation_output_dir=lambda:str(annotations), location=str(tmp_path), accessions=[], log=lambda *args:None)
    manifest = AddLibraryTask._materialize_replacement_orthogroup_trees(task, [('f','OG')], {})
    with open(manifest) as handle: entries = list(csv.DictReader(handle, delimiter='\t'))
    assert entries[0]['outcome'] == 'omitted' and entries[0]['tree_path'] == ''
    assert not canonical.exists() and not (annotations/'OG.nex').exists()
    assert (canonical_dir/'Resolved_Gene_Trees.txt').read_text() == ''


def test_empty_annotation_set_finishes(tmp_path):
    (tmp_path/'orthogroup_tree_manifest.tsv').write_text('orthogroup\toutcome\nOG\tomitted\n')
    task = SimpleNamespace(location=str(tmp_path), _annotate_og_trees_effective=lambda:True,
                           _annotation_input_dir=lambda:str(tmp_path/'missing'), _annotation_output_dir=lambda:str(tmp_path/'missing-output'))
    assert AddLibraryTask._annotation_done(task)
