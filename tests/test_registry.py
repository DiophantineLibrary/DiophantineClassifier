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


def test_matcher_flags_match_implementation():
    """Every registered slug the matcher code can emit is flagged matcher: true.

    A recognizer whose family has not landed in the registry yet is inert:
    :func:`~diophantine_classifier.classify.classify` drops matches naming an
    unregistered family.  The check tightens to *every* emitted slug once the
    registry is complete.
    """
    import re
    import diophantine_classifier.matchers as m
    source = open(m.__file__.replace(".pyc", ".py")).read()
    emitted = set(re.findall(r'Match\(\s*\n?\s*"([a-z0-9-]+)"', source))
    emitted |= set(re.findall(r'Match\("([a-z0-9-]+)"', source))
    fams = families()
    for slug in sorted(emitted & set(fams)):
        assert fams[slug].matcher, f"{slug}: emitted but matcher flag is false"


def test_examples_classify_to_their_family():
    """Registry examples must classify to the family that lists them (or a
    descendant, e.g. the fermat example listed under generalized-fermat)."""
    from diophantine_classifier import classify
    from diophantine_classifier.registry import ancestors as anc
    fams = families()
    for fam in fams.values():
        if not fam.matcher:
            continue
        for ex in fam.examples:
            cls = classify(ex, params=_infer_params(ex))
            ok = (cls.slug == fam.slug or fam.slug in anc(cls.slug)
                  or cls.slug in anc(fam.slug))
            assert ok, (f"{fam.slug}: example {ex!r} classified as "
                        f"{cls.slug!r}")


def _infer_params(example):
    """Parameters used in registry examples: single letters not serving as
    unknowns in that example's family form."""
    KNOWN = {"x^2 + y^2 + z^2 = n": "n", "x^2 + y^2 + z^2 + w^2 = n": "n",
             "x^4 + y^4 + z^4 + w^4 = n": "n",
             "4/n = 1/x + 1/y + 1/z": "n"}
    return KNOWN.get(example, "")


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


def test_sage_code_templates_run_on_their_examples():
    """Every filled Sage template must execute on the family's own examples
    (this is what a user copies from the equation page)."""
    import sage.all
    from sage.repl.preparse import preparse
    from diophantine_classifier import classify
    ran = 0
    for fam in families().values():
        if "sage" not in fam.code or not fam.matcher:
            continue
        for ex in fam.examples:
            cls = classify(ex, params=_infer_params(ex))
            match = next((m for m in cls.matches if m.slug == fam.slug), None)
            if match is None:
                continue
            code = fam.fill_code(match.data)["sage"]
            assert "{" not in code, (fam.slug, code)
            env = dict(vars(sage.all))
            exec(preparse(code), env)
            ran += 1
    assert ran >= 1

