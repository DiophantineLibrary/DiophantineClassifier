"""Classification corpus: famous equations -> expected most-specific family.

Run with:  sage -python -m pytest tests/ -q
"""

import pytest

from diophantine_classifier import classify

# (equation, params, expected primary slug)
CORPUS = [
    # linear / univariate
    ("3*x + 5*y = 1", "", "linear"),
    ("12*x - 21*y + 30*z = 9", "", "linear"),
    ("x^2 - 5*x + 6 = 0", "", "univariate"),
    # quadratic, two variables
    ("x^2 - 61*y^2 = 5", "", "pell-like"),
    ("2*x^2 + 3*x*y - 5*y^2 + x - 7 = 0", "", "binary-quadratic"),
    ("3*x^2 + 7*y^2 = 19", "", "binary-qf-representation"),
    # quadratic, more variables
    ("x^2 + 3*y^2 = 7*z^2", "", "legendre"),
    ("x^2 - 3*y^2 + 5*z^2 - 7*w^2 = 0", "", "quadratic-form-zero"),
    ("x^2 + x*y + y^2 + z^2 = 14", "", "quadratic-form-representation"),
    ("x^2 + y^2 - z^2 + 3*x - 7 = 0", "", "quadric"),
    # genus one
    ("y^2 + y = x^3 - x^2 - 10*x - 20", "", "elliptic-weierstrass"),
    # higher-genus curves and binary forms
    ("x^3 + 2*y^3 = 11", "", "thue"),
    ("x^4 - 2*y^4 = 1", "", "thue"),
    ("y^2 = x^7 + 3", "", "hyperelliptic"),
    ("y^3 = x^4 + 2", "", "superelliptic"),
    ("x^3*y + y^3*z + z^3*x = 0", "", "general-curve"),     # Klein quartic
    ("x^2*y^2 = x^3 + 1", "", "genus-one-curve"),
    # Fermat-type
    ("x^2 + y^4 = z^3", "", "generalized-fermat"),
    ("2*x^3 + 3*y^3 = 5*z^3", "", "generalized-fermat"),
    ("3*x^3 + 4*y^3 + 5*z^3 = 0", "", "generalized-fermat"),  # Selmer
    ("x^p + y^q = z^r", "", "generalized-fermat"),            # Beal
    # diagonal / surfaces
    ("x^4 + y^4 + z^4 = w^4", "", "equal-sums-like-powers"),   # Elkies
    # polynomial-exponential
    ("x^3 - 4 = y^n", "", "power-values"),
    ("x^p - y^q = 1", "", "catalan"),
    # unit fractions
    ("1/x + 1/y + 1/z = 1", "", "egyptian-fractions"),
]


@pytest.mark.parametrize("equation,params,expected", CORPUS,
                         ids=[c[0] for c in CORPUS])
def test_corpus(equation, params, expected):
    cls = classify(equation, params=params)
    assert cls.slug == expected, (
        f"{equation!r}: got {cls.slug!r} (matches: "
        f"{[m.slug for m in cls.matches]}), expected {expected!r}")


def test_reducible():
    cls = classify("(x^2 - 2)*(y^2 - 3) = 0")
    assert cls.slug == "reducible"
    assert sorted(c.slug for c in cls.components) == ["univariate", "univariate"]


def test_parametric_quadratic_form_classifies():
    """A parametric quadratic form must classify, not raise (matchers._gram)."""
    cls = classify("x^2 + y^2 = D*z^2", params="D")
    assert cls.slug in {"quadratic-form-zero", "general-polynomial"}


def test_gen_fermat_regimes():
    spherical = classify("x^2 + y^4 = z^3")
    assert spherical.data["regime"] == "spherical"
    hyperbolic = classify("x^2 + y^7 = z^3")
    assert hyperbolic.data["regime"] == "hyperbolic"


def test_match_lookup_by_slug():
    cls = classify("3*x + 5*y = 1")
    assert cls.match_for("linear").slug == "linear"
    assert cls.data_for("linear")["b"] == "1"
    assert cls.match_for("no-such-family") is None
    assert cls.data_for("no-such-family") is None


def test_code_is_filled_from_each_match_not_the_primary():
    """An ancestor's template gets the ancestor's own data."""
    cls = classify("x^2 + 3*y^2 = 7*z^2")
    assert cls.slug == "legendre"
    code = cls.code()
    assert "[1, 3, -7]" in code["sage (legendre)"]
    # quadratic-form-zero wants a Gram matrix, which the Legendre match's
    # a/b/c cannot supply; it was matched too, so it fills from its own data
    gram = code["sage (quadratic-form-zero)"]
    assert "{gram}" not in gram
    assert "[[1, 0, 0], [0, 3, 0], [0, 0, -7]]" in gram


def test_repeated_factor_preserves_input():
    text = "(x + y)^2 = 0"
    cls = classify(text, domain="QQ")
    assert cls.parsed.original == text
    assert cls.as_dict()["equation"] == text
    assert cls.parsed.domain == "QQ"
    assert cls.slug == "linear"


def test_repeated_factor_records_the_working_model():
    cls = classify("(x + y)^2 = 0", domain="QQ")
    assert cls.working is not cls.parsed
    assert str(cls.working.poly) == "x + y"
    assert cls.working.domain == "QQ"
    assert cls.reduction is not None and cls.reduction.is_identity


def test_repeated_factor_preserves_unknown_order():
    cls = classify("(y + x)^2 = 0")
    assert cls.parsed.unknowns == ("y", "x")
    assert cls.as_dict()["unknowns"] == ["y", "x"]


def test_components_inherit_the_domain():
    cls = classify("(x^2 - 2)*(y^2 - 3) = 0", domain="QQ")
    assert cls.slug == "reducible"
    assert all(comp.parsed.domain == "QQ" for comp in cls.components)


def test_components_inherit_the_parent_conditions():
    """A factor of a rational equation is still subject to its conditions."""
    cls = classify("(x - 2)*(x - y)/y = 0")
    parent = {str(c) for c in cls.parsed.conditions}
    assert "y != 0" in parent
    for comp in cls.components:
        checkable = {str(c) for c in comp.parsed.conditions}
        if "y" in comp.parsed.unknowns:
            assert "y != 0" in checkable


# --- ancestry is displayed and serialized as paths (brief 4.1 / 6.3) ------

def _lineage_lines(cls):
    return [line.strip()[len("lineage:"):].strip()
            for line in cls.explain().splitlines()
            if line.strip().startswith("lineage:")]


def test_explain_draws_only_real_edges():
    """One arrow line per genuine path.  Rendering the flattened ancestor
    set as a single chain would draw edges the DAG does not have."""
    from diophantine_classifier.registry import families
    fams = families()
    cls = classify("3*x + 5*y = 1")
    lines = _lineage_lines(cls)
    assert lines
    for line in lines:
        chain = [part.strip() for part in line.split("→")]
        assert chain[0] == cls.slug
        for child, parent in zip(chain, chain[1:]):
            assert parent in fams[child].parents, (child, parent)


def test_as_dict_carries_paths_and_graph():
    import json
    cls = classify("3*x + 5*y = 1")
    d = cls.as_dict()
    assert d["lineage_paths"] == [[cls.slug, "general-polynomial"]]
    assert d["lineage_graph"]["edges"] == [[cls.slug, "general-polynomial"]]
    json.dumps(d)  # still the website contract


def test_lineage_output_is_deterministic():
    runs = [classify("3*x + 5*y = 1").as_dict()["lineage_paths"]
            for _ in range(5)]
    assert all(run == runs[0] for run in runs)


def test_conditions_serialize_structurally():
    import json
    d = classify("1/x = y").as_dict()
    condition, = d["conditions"]
    assert condition == {"type": "nonzero", "expression": "x",
                         "variables": ["x"], "source": "denominator"}
    json.dumps(d)


def test_multi_parent_lineage_has_no_false_edge():
    """legendre specializes both quadratic-form-zero and diagonal-form;
    those two are siblings, and no output may put an arrow between them."""
    cls = classify("x^2 + 3*y^2 = 7*z^2")
    assert cls.slug == "legendre"
    paths = cls.as_dict()["lineage_paths"]
    assert {path[1] for path in paths} == {"quadratic-form-zero",
                                           "diagonal-form"}
    for line in _lineage_lines(cls):
        assert "quadratic-form-zero → diagonal-form" not in line
        assert "diagonal-form → quadratic-form-zero" not in line


# --- conditional identities never reach the matchers (brief 3.3) ---------

def test_conditional_identity_does_not_enter_polynomial_matchers():
    cls = classify("x/x = 1")
    assert cls.parsed.is_conditional_identity
    assert not cls.matches
    assert cls.special_kind == "conditional-identity"


def test_conditional_identity_is_not_called_a_polynomial():
    cls = classify("x/x = 1")
    assert cls.slug != "general-polynomial"
    d = cls.as_dict()
    assert d["equation"] == "x/x = 1"
    assert d["conditions"] == [{"type": "nonzero", "expression": "x",
                                "variables": ["x"], "source": "denominator"}]
    assert d["special_kind"] == "conditional-identity"


def test_conditional_identity_explains_itself():
    text = classify("x/x = 1").explain()
    assert "condition: x != 0" in text
    assert "conditional identity" in text


def test_conditional_identity_keeps_the_domain_and_order():
    cls = classify("(y - x)/(y - x) = 1", domain="QQ")
    assert cls.parsed.unknowns == ("y", "x")
    assert cls.as_dict()["domain"] == "QQ"


# --- the repeated-factor path is filtered too (brief 6.1) ----------------

def test_repeated_factor_drops_unregistered_matches_before_ranking():
    from diophantine_classifier import families
    cls = classify("(x^2 + y^2 - z^2)^2 = 0")
    assert cls.parsed.original == "(x^2 + y^2 - z^2)^2 = 0"
    assert all(m.slug in families() for m in cls.matches)


def test_every_classification_path_filters_the_registry():
    from diophantine_classifier import families
    known = families()
    for equation in ["3*x + 5*y = 1", "(x + y)^2 = 0",
                     "(x^2 - 2)*(y - 3) = 0", "(x^2 + y^2 - z^2)^2 = 0"]:
        cls = classify(equation)
        assert all(m.slug in known for m in cls.matches), equation
        for comp in cls.components:
            assert all(m.slug in known for m in comp.matches), equation


# --- components live in the ambient space (brief 6.2) --------------------

def component_with_factor(cls, text):
    return next(c for c in cls.components if str(c.working.poly) == text)


def test_component_retains_source_conditions_and_free_variables():
    cls = classify("(x - 2)*(x - y)/y = 0")
    comp = component_with_factor(cls, "x - 2")
    assert comp.parsed.original == cls.parsed.original
    assert comp.free_variables == ("y",)
    assert "y != 0" in {str(c) for c in comp.parsed.conditions}


def test_component_keeps_the_ambient_unknown_order():
    cls = classify("(x^2 - 2)*(y - 3) = 0")
    for comp in cls.components:
        assert comp.parsed.unknowns == ("x", "y")
        assert comp.embedding.source_variables == ("x", "y")


def test_univariate_factor_in_a_larger_ambient_space_is_not_collapsed():
    cls = classify("(x^2 - 2)*(y - 3) = 0")
    comp = component_with_factor(cls, "x^2 - 2")
    assert comp.free_variables == ("y",)
    assert comp.embedding.active_variables == ("x",)
    assert comp.reduction is None      # a projection does not invert


def test_component_embedding_serializes():
    import json
    cls = classify("(x^2 - 2)*(y - 3) = 0")
    comp = component_with_factor(cls, "x^2 - 2")
    blob = comp.as_dict()["component"]
    assert blob["free_variables"] == ["y"]
    json.dumps(blob)


# --- the effective transform (brief 6.3) ---------------------------------

def test_effective_transform_is_the_match_transform_when_nothing_reduced():
    cls = classify("3*x + 5*y = 1")
    assert cls.effective_transform(cls.primary) is cls.primary.transform


def test_effective_transform_composes_the_reduction():
    cls = classify("(x + y)^2 = 0")
    assert cls.reduction is not None
    effective = cls.effective_transform(cls.primary)
    assert effective.source_variables == ("x", "y")
    assert "reduced" in effective.description
    source = {"x": 2, "y": -2}
    assert effective.pull_back(effective.push_forward(source)) == source


def test_as_dict_exposes_the_effective_map_and_its_provenance():
    cls = classify("(x + y)^2 = 0")
    d = cls.as_dict()
    assert "reduced" in d["transform"]["description"]
    assert d["reduction"]["description"] == d["transform"]["description"]
