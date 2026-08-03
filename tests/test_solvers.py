"""Solver layer: concrete instances actually solved via Sage/PARI."""

import pytest

from diophantine_classifier import classify, solve, SolverUnavailable


def test_pell_61():
    s = solve("x^2 - 61*y^2 = 1")
    assert s.kind == "infinite"
    assert s.solutions[0] == (1766319049, 226153980)


def test_negative_pell():
    s = solve("x^2 - 2*y^2 = -1")
    assert s.solutions[0] == (1, 1)


def test_pell_like():
    s = solve("x^2 - 2*y^2 = 7")
    assert s.solutions
    for x, y in s.solutions:
        assert x ** 2 - 2 * y ** 2 == 7
    assert s.complete    # orbit representatives + automorph action


def test_two_squares():
    s = solve("x^2 + y^2 = 610")
    assert s.solutions == [(9, 23), (13, 21)]
    assert s.complete


def test_three_squares_obstruction():
    s = solve("x^2 + y^2 + z^2 = 7")
    assert s.kind == "empty" and s.complete


def test_four_squares():
    s = solve("x^2 + y^2 + z^2 + w^2 = 7")
    assert s.solutions == [(1, 1, 1, 2)]
    assert s.complete


def test_bqf():
    s = solve("3*x^2 + 7*y^2 = 19")
    assert s.solutions == [(-2, -1), (-2, 1), (2, -1), (2, 1)]
    assert s.complete


def test_legendre():
    s = solve("x^2 + y^2 = 2*z^2")
    ((x, y, z),) = s.solutions
    assert x ** 2 + y ** 2 == 2 * z ** 2
    assert (x, y, z) != (0, 0, 0)


def test_legendre_obstruction():
    s = solve("x^2 + y^2 = 3*z^2")
    assert s.kind == "empty" and s.complete


def test_linear():
    s = solve("3*x + 5*y = 1")
    ((x, y),) = s.solutions
    assert 3 * x + 5 * y == 1


def test_linear_empty():
    s = solve("6*x + 9*y = 5")
    assert s.kind == "empty" and s.complete


def test_univariate():
    s = solve("x^2 - 5*x + 6 = 0")
    assert [t[0] for t in s.solutions] == [2, 3]


def test_mordell_integral_points():
    s = solve("y^2 = x^3 - 2")
    # unknowns are ordered (y, x) by appearance in the input
    assert set(s.solutions) == {(5, 3), (-5, 3)}
    assert s.complete


def test_thue():
    s = solve("x^3 + 2*y^3 = 11")
    assert (3, -2) in s.solutions
    assert s.complete


def test_catalan():
    s = solve("x^p - y^q = 1")
    # unknowns ordered (x, p, y, q)
    assert s.solutions == [(3, 2, 2, 3)]
    assert s.complete


def test_fermat():
    s = solve("x^4 + y^4 = z^4")
    assert s.kind == "empty" and s.complete


def test_ramanujan_nagell():
    s = solve("x^2 + 7 = 2^n")
    assert (11, 7) in s.solutions and len(s.solutions) == 5
    assert s.complete


def test_unavailable_carries_hints():
    with pytest.raises(SolverUnavailable) as err:
        solve("y^2 = x^7 + 3")
    assert "magma" in str(err.value).lower() or "Chabauty" in str(err.value)


def test_solve_accepts_classification():
    cls = classify("x^2 - 61*y^2 = 1")
    s = solve(cls)
    assert s.solutions
