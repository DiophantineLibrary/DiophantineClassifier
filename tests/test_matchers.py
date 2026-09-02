"""Matchers: what a named family requires, and what it must not swallow."""

from diophantine_classifier.matchers import run
from diophantine_classifier.parsing import parse


def slugs(equation, **kwargs):
    return {m.slug for m in run(parse(equation, **kwargs))}


def test_parametric_quadratic_form_does_not_crash():
    """Coefficients live in QQ[params]; a Gram matrix must not force QQ."""
    found = slugs("x^2 + y^2 = D*z^2", params="D")
    assert "quadratic-form-zero" in found


def test_weighted_diagonal_is_not_equal_sums():
    """Equal sums of like powers means the powers are summed, unweighted."""
    found = slugs("2*x^3 + y^3 = z^3 + 7*w^3")
    assert "equal-sums-like-powers" not in found
    assert "diagonal-form" in found


def test_unweighted_equal_sums_still_matches():
    assert "equal-sums-like-powers" in slugs("x^4 + y^4 + z^4 = w^4")


def test_sum_of_three_cubes_still_matches():
    assert "sum-of-three-cubes" in slugs("x^3 + y^3 + z^3 = 42")


# --- the Pell orientation swap is an equivalence ---------------------------

def test_swapped_pell_keeps_the_same_equation():
    """5*x^2 - y^2 = 1 reads as y^2 - 5*x^2 = -1, so N flips with the swap."""
    match, = [m for m in run(parse("5*x^2 - y^2 = 1")) if m.slug == "pell"]
    assert match.data == {"D": "5", "N": "-1"}


def test_unswapped_pell_is_unchanged():
    match, = [m for m in run(parse("x^2 - 61*y^2 = 1")) if m.slug == "pell"]
    assert match.data == {"D": "61", "N": "1"}


def test_swapped_pell_like_keeps_the_same_equation():
    match, = [m for m in run(parse("3*x^2 - y^2 = 6")) if m.slug == "pell-like"]
    assert match.data == {"D": "3", "N": "-6"}


# --- parametric coefficients are part of the coefficient (brief 5.1) -------

def test_parametric_multiplier_blocks_pillai():
    """A*2^n - 3^m = 1 is not Pillai's equation: the multiplier is unknown."""
    assert "pillai" not in slugs("A*2^n - 3^m = 1", params="A")


def test_concrete_multiplier_still_matches_pillai():
    assert "pillai" in slugs("2^n - 3^m = 1")


def test_parametric_leading_coefficient_blocks_ramanujan_nagell():
    assert "ramanujan-nagell" not in slugs("A*x^2 + 7 = 2^n", params="A")


def test_concrete_leading_coefficient_still_matches_ramanujan_nagell():
    assert "ramanujan-nagell" in slugs("x^2 + 7 = 2^n")


def test_parametric_multiplier_blocks_variable_base_pillai():
    """The variable-base branch reads |coeff| = 1 too."""
    assert "pillai" not in slugs("A*x^p - y^q = 2", params="A")
    assert "catalan" not in slugs("A*x^p - y^q = 1", params="A")


def test_parametric_coefficient_blocks_symbolic_fermat():
    assert "fermat" not in slugs("A*x^n + y^n = z^n", params="A")


def test_parametric_coefficient_blocks_lebesgue_nagell():
    assert "lebesgue-nagell" not in slugs("A*x^2 + 3 = y^n", params="A")
    assert "lebesgue-nagell" in slugs("x^2 + 3 = y^n")


def test_parametric_multiplier_blocks_thue_mahler():
    assert "thue-mahler" not in slugs("x^3 + 2*y^3 = A*2^n", params="A")
    assert "thue-mahler" in slugs("x^3 + 2*y^3 = 5*2^n")


def test_parametric_coefficient_survives_in_the_broad_fallback():
    """A blocked classical shape still gets its root family."""
    assert "polynomial-exponential" in slugs("A*x^2 + 7 = 2^n", params="A")


def test_full_coefficient_is_available_on_terms():
    from diophantine_classifier.parsing import parse
    pe = parse("A*x^2 + 7 = 2^n", params="A")
    quadratic, = [t for t in pe.terms if t.powers]
    assert quadratic.has_parametric_coefficient
    assert quadratic.coefficient_string() == "A"


# --- symbolic generalized-Fermat signatures (brief 5.4b) ------------------

def test_repeated_symbolic_signature_is_generalized_fermat():
    """x^p + y^p = z^q: the exponents need not be pairwise distinct."""
    assert "generalized-fermat" in slugs("x^p + y^p = z^q")


def test_distinct_symbolic_signature_is_still_generalized_fermat():
    assert "generalized-fermat" in slugs("x^p + y^q = z^r")


def test_all_equal_symbolic_signature_is_fermat():
    found = slugs("x^n + y^n = z^n")
    assert "fermat" in found
    assert "generalized-fermat" not in found
