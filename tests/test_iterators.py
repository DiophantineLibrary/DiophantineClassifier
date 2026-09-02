"""Iteration over solution sets: infinite families yield streams."""

from sage.all import gcd

from diophantine_classifier import solve


def test_linear_stream():
    s = solve("3*x + 5*y = 1")
    sols = s.first(9)
    assert len(set(sols)) == 9
    assert all(3 * x + 5 * y == 1 for x, y in sols)


def test_egyptian_complete():
    s = solve("1/x + 1/y + 1/z = 1")
    assert s.solutions == [(2, 3, 6), (2, 4, 4), (3, 3, 3)]
    assert s.complete


