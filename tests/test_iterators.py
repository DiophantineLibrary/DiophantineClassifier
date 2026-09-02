"""Iteration over solution sets: infinite families yield streams."""

from sage.all import gcd

from diophantine_classifier import solve


def test_pell_like_stream():
    s = solve("x^2 - 2*y^2 = 7")
    sols = s.first(12)
    assert len(set(sols)) == 12
    assert all(x ** 2 - 2 * y ** 2 == 7 for x, y in sols)


def test_linear_stream():
    s = solve("3*x + 5*y = 1")
    sols = s.first(9)
    assert len(set(sols)) == 9
    assert all(3 * x + 5 * y == 1 for x, y in sols)


def test_markov_stream():
    s = solve("x^2 + y^2 + z^2 = 3*x*y*z")
    assert s.first(5) == [(1, 1, 1), (1, 1, 2), (1, 2, 5), (1, 5, 13),
                          (2, 5, 29)]
    assert all(x * x + y * y + z * z == 3 * x * y * z
               for x, y, z in s.first(12))


def test_hurwitz_stream():
    s = solve("x^2 + y^2 + z^2 + w^2 = 4*x*y*z*w")
    sols = s.first(3)
    assert sols[0] == (1, 1, 1, 1)
    assert all(a * a + b * b + c * c + d * d == 4 * a * b * c * d
               for a, b, c, d in sols)


def test_egyptian_complete():
    s = solve("1/x + 1/y + 1/z = 1")
    assert s.solutions == [(2, 3, 6), (2, 4, 4), (3, 3, 3)]
    assert s.complete


def test_finite_iteration_matches_list():
    s = solve("x^3 + 2*y^3 = 11")
    assert list(iter(s)) == s.solutions
