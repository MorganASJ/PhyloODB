from types import SimpleNamespace
from ete3 import Tree
from phyloODB.tasks.analysis.trees import _OrthogroupTreeAnnotationMixin, _flags_from_env_or_task
from phyloODB.tasks.utilities.orthobusco_analyzer import OrthoBuscoAnalyzer


def test_flag_precedence():
    assert _flags_from_env_or_task(env_value=None, task_value=None, default_flags=['-m', 'MFP']) == ['-m', 'MFP']
    assert _flags_from_env_or_task(env_value='-m LG -B 1000', task_value='-m WAG', default_flags=[]) == ['-m', 'WAG']


def test_annotation_support_and_leaf_names(tmp_path):
    class Renderer(_OrthogroupTreeAnnotationMixin):
        def _load_busco_markers(self, **kwargs): return {}
        def _load_hidden_paralog_leaf_names(self, **kwargs): return set()
        def _resolve_leaf_metadata(self, name, metadata): return {'accession': 'a', 'sequence': name}
        def _taxonomy_labels_for_accession(self, accession): return {'taxon': 'Taxon', 'phylum': 'Phylum'}
    for label in ('', '97.50', '80/99'):
        path = tmp_path / 'tree.nwk'
        path.write_text(f'((a:1,b:2){label}:3,c:4);')
        out = Renderer()._render_nexus(tree_path=str(path), family_id='', source_run_ids=[], in_leaves=set(), out_leaves=set(), leaf_metadata={})
        assert "a_Taxon_Phylum" in out
        assert '97.50_Taxon' not in out
        if label: assert f"'{label}':3.0" in out
        else: assert '):3.0' in out


def test_paralogy_independent_of_support():
    analyzer = SimpleNamespace(_species_key=lambda leaf: leaf[0], log=lambda *args: None)
    results = [OrthoBuscoAnalyzer._paralogs_for_tree(analyzer, 'OG', Tree(nwk, format=1), '') for nwk in ('((a1,a2),b1,c1);', '((a1,a2)95,b1,c1);')]
    assert results[0] == results[1]
