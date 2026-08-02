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
    ("x^2 - 61*y^2 = 1", "", "pell"),
    ("x^2 - 2*y^2 = -1", "", "pell"),
    ("x^2 - 61*y^2 = 5", "", "pell-like"),
    ("x^2 - D*y^2 = 1", "D", "pell"),
    ("x^2 + y^2 = 610", "", "sum-of-two-squares"),
    ("2*x^2 + 3*x*y - 5*y^2 + x - 7 = 0", "", "binary-quadratic"),
    ("3*x^2 + 7*y^2 = 19", "", "binary-qf-representation"),
    # quadratic, more variables
    ("x^2 + y^2 = z^2", "", "pythagorean"),
    ("x^2 + 3*y^2 = 7*z^2", "", "legendre"),
    ("x^2 + y^2 + z^2 = n", "n", "sum-of-three-squares"),
    ("x^2 + y^2 + z^2 + w^2 = n", "n", "sum-of-four-squares"),
    ("x^2 - 3*y^2 + 5*z^2 - 7*w^2 = 0", "", "quadratic-form-zero"),
    ("x^2 + x*y + y^2 + z^2 = 14", "", "quadratic-form-representation"),
    ("x^2 + y^2 - z^2 + 3*x - 7 = 0", "", "quadric"),
    # genus one
    ("y^2 = x^3 - 2", "", "mordell"),
    ("y^2 = x^3 + k", "k", "mordell"),
    ("x^3 = y^2 - 2", "", "mordell"),
    ("y^2 + y = x^3 - x^2 - 10*x - 20", "", "elliptic-weierstrass"),
    ("y^2 = x^4 + 3*x + 1", "", "elliptic-quartic"),
    # higher-genus curves and binary forms
    ("x^3 + 2*y^3 = 11", "", "thue"),
    ("x^4 - 2*y^4 = 1", "", "thue"),
    ("x^3 + y^3 = 1729", "", "binary-form-reducible"),
    ("y^2 = x^5 - x + 1", "", "genus-two"),
    ("y^2 = x^7 + 3", "", "hyperelliptic"),
    ("y^3 = x^4 + 2", "", "superelliptic"),
    ("x^3*y + y^3*z + z^3*x = 0", "", "general-curve"),     # Klein quartic
    ("x^2*y^2 = x^3 + 1", "", "genus-one-curve"),
    ("y^2 = x^3", "", "genus-zero-curve"),                  # cuspidal
    # Fermat-type
    ("x^4 + y^4 = z^4", "", "fermat"),
    ("x^n + y^n = z^n", "", "fermat"),
    ("x^2 + y^4 = z^3", "", "generalized-fermat"),
    ("2*x^3 + 3*y^3 = 5*z^3", "", "generalized-fermat"),
    ("3*x^3 + 4*y^3 + 5*z^3 = 0", "", "generalized-fermat"),  # Selmer
    ("x^p + y^q = z^r", "", "generalized-fermat"),            # Beal
    # diagonal / surfaces
    ("x^3 + y^3 + z^3 = 42", "", "sum-of-three-cubes"),
    ("x^2 + y^2 + z^2 = 3*x*y*z", "", "markov-hurwitz"),
    ("x^2 + y^2 + z^2 + w^2 = 4*x*y*z*w", "", "markov-hurwitz"),
    ("x^4 + y^4 + z^4 = w^4", "", "equal-sums-like-powers"),   # Elkies
    ("x^4 + y^4 + z^4 + w^4 = n", "n", "waring"),
    # polynomial-exponential
    ("x^2 + 7 = 2^n", "", "ramanujan-nagell"),
    ("x^2 + 11 = 3^n", "", "ramanujan-nagell"),
    ("x^2 + 2 = y^n", "", "lebesgue-nagell"),
    ("x^3 - 4 = y^n", "", "power-values"),
    ("x^p - y^q = 1", "", "catalan"),
    ("3^m - 2^n = 1", "", "pillai"),
    ("3^m - 2^n = 5", "", "pillai"),
    ("2^a + 3^b = 5^c", "", "s-unit"),
    ("x^3 + 2*y^3 = 5^a * 11^b", "", "thue-mahler"),
    # unit fractions
    ("4/n = 1/x + 1/y + 1/z", "n", "erdos-straus"),
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


def test_lineage_pell():
    cls = classify("x^2 - 61*y^2 = 1")
    assert cls.lineage[0] == "pell-like"
    assert "binary-qf-representation" in cls.lineage
    assert cls.data["D"] == "61"


def test_mordell_data():
    cls = classify("y^2 = x^3 - 2")
    assert cls.data["k"] == "-2"
    assert "elliptic-weierstrass" in cls.lineage


def test_gen_fermat_regimes():
    spherical = classify("x^2 + y^4 = z^3")
    assert spherical.data["regime"] == "spherical"
    hyperbolic = classify("x^2 + y^7 = z^3")
    assert hyperbolic.data["regime"] == "hyperbolic"


def test_params_as_list():
    cls = classify("y^2 = x^3 + k", params=["k"])
    assert cls.slug == "mordell"


def test_as_dict_roundtrip():
    import json
    cls = classify("x^2 - 61*y^2 = 1")
    blob = json.dumps(cls.as_dict())
    assert "pell" in blob


def test_explain_smoke():
    text = classify("x^2 + 7 = 2^n").explain()
    assert "ramanujan-nagell" in text
    assert "Nagell" in text or "1948" in text


def test_conditions_recorded():
    cls = classify("4/n = 1/x + 1/y + 1/z", params="n")
    assert cls.parsed.conditions  # denominators were cleared
