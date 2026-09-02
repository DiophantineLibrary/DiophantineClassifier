"""Registry consistency: DAG well-formed, matcher flags honest, docs in sync."""

from diophantine_classifier import families, validate
from diophantine_classifier.registry import ancestors, depth


def test_validate():
    assert validate()


def test_dag_depths():
    fams = families()
    assert depth("general-polynomial") == 0
    assert depth("linear") > depth("general-polynomial")
    assert "general-polynomial" in ancestors("linear")


# --- ancestry is paths and edges, not a flattened chain (brief 4.1) -------

def test_path_traversal_on_a_diamond_invents_no_edge():
    """child -> {left, right} -> root has two paths, and the two middle
    nodes are never adjacent -- which is exactly what a flattened ancestor
    list joined by arrows would claim."""
    from diophantine_classifier.registry import _paths_to_roots
    diamond = {"child": ("left", "right"), "left": ("root",),
               "right": ("root",), "root": ()}
    paths = _paths_to_roots("child", diamond.__getitem__)
    assert paths == [["child", "left", "root"], ["child", "right", "root"]]
    for path in paths:
        adjacent = set(zip(path, path[1:]))
        assert ("left", "right") not in adjacent
        assert ("right", "left") not in adjacent


def test_path_traversal_is_deterministic():
    """Order does not depend on how the parents happen to be listed."""
    from diophantine_classifier.registry import _paths_to_roots
    one = {"child": ("left", "right"), "left": ("root",),
           "right": ("root",), "root": ()}
    other = {"child": ("right", "left"), "left": ("root",),
             "right": ("root",), "root": ()}
    runs = [_paths_to_roots("child", one.__getitem__) for _ in range(5)]
    assert all(run == runs[0] for run in runs)
    assert _paths_to_roots("child", other.__getitem__) == runs[0]


def test_lineage_paths_use_only_real_edges():
    """Every consecutive pair of every path is a genuine parent edge, and
    every path ends at a root."""
    from diophantine_classifier.registry import lineage_paths
    fams = families()
    for slug in fams:
        paths = lineage_paths(slug)
        assert paths
        for path in paths:
            assert path[0] == slug
            for child, parent in zip(path, path[1:]):
                assert parent in fams[child].parents, (slug, child, parent)
            assert not fams[path[-1]].parents


def test_lineage_paths_reach_every_ancestor():
    """The paths carry the same information as the flat ancestor set."""
    from diophantine_classifier.registry import ancestors as anc
    from diophantine_classifier.registry import lineage_paths
    for slug in families():
        reached = {node for path in lineage_paths(slug) for node in path[1:]}
        assert reached == set(anc(slug)), slug


def test_lineage_graph_edges_are_real():
    from diophantine_classifier.registry import lineage_graph
    fams = families()
    for slug in fams:
        graph = lineage_graph(slug)
        assert slug in graph["nodes"]
        for child, parent in graph["edges"]:
            assert parent in fams[child].parents
            assert child in graph["nodes"]
            assert parent in graph["nodes"]


def test_multi_parent_family_gets_one_path_per_parent():
    """Inert until a family with two parents is registered on this branch."""
    from diophantine_classifier.registry import lineage_paths
    fams = families()
    for slug, fam in fams.items():
        if len(fam.parents) < 2:
            continue
        assert {path[1] for path in lineage_paths(slug)} == set(fam.parents)


def test_ancestors_order_is_deterministic_and_by_depth():
    from diophantine_classifier.registry import ancestors as anc
    for slug in families():
        order = anc(slug)
        assert order == anc(slug)
        depths = [depth(s) for s in order]
        assert depths == sorted(depths, reverse=True)
        assert order == sorted(order, key=lambda s: (-depth(s), s))

# --- registry hygiene (PR review) ------------------------------------------

def test_list_entries_have_balanced_parentheses():
    """YAML flow lists split on commas: '(I, J)' must not become two items."""
    for fam in families().values():
        for entry in fam.methods + fam.aliases + fam.examples:
            assert entry.count("(") == entry.count(")"), (fam.slug, entry)

