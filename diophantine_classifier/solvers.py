r"""
Solvers: dispatch a classified equation to Sage/PARI machinery.

``solve`` handles the families where standard software gives a complete or
canonical answer.  The result is a :class:`SolutionSet`:

- for *finite* solution sets it holds the complete list (``complete=True``);
- for *infinite* ones it is **iterable**: ``iter(S)`` enumerates solutions
  indefinitely (Pell solutions by powers of the fundamental unit, lattice
  cosets by shells, the Markov tree by increasing maximum, ...), while
  ``S.solutions`` keeps a distinguished finite piece (fundamental solutions,
  orbit representatives) for display;
- for families without a wired-up solver, :class:`SolverUnavailable` is
  raised, carrying the registry's software pointers and filled code
  templates — every equation page should offer runnable code even when we
  cannot run it ourselves.

EXAMPLES::

    sage: from diophantine_classifier import solve
    sage: solve("x^2 - 5*x + 6 = 0").solutions       # finite: complete list
    [(2,), (3,)]
    sage: S = solve("3*x + 5*y = 1")                 # infinite: iterate
    sage: S.first(6)
    [(2, -1), (-3, 2), (7, -4), (-8, 5), (12, -7), (-13, 8)]
"""

import heapq
import itertools
from dataclasses import dataclass, field, replace

from sage.all import (QQ, ZZ, BinaryQF, EllipticCurve, QuadraticField,
                      continued_fraction, gcd, isqrt, lcm, matrix, pari,
                      sage_eval, xgcd, two_squares, three_squares,
                      four_squares)
from sage.quadratic_forms.qfsolve import qfsolve

from .classify import Classification, classify

#: bounds above which complete enumeration falls back to a single witness
MAX_TWO_SQUARES = 10**10
MAX_THREE_SQUARES = 10**6
MAX_FOUR_SQUARES = 10**4
MAX_BQF = 10**8
MAX_EGYPTIAN_N = 10**4


class SolverUnavailable(NotImplementedError):
    r"""
    No automatic solver for this family (yet).

    The message carries the registry's software pointers and code templates.

    EXAMPLES::

        sage: from diophantine_classifier import solve, SolverUnavailable
        sage: try:
        ....:     solve("y^2 = x^7 + 3")
        ....: except SolverUnavailable as err:
        ....:     print("no automatic solver" in str(err))
        True
    """


@dataclass
class SolutionSet:
    r"""
    The result of :func:`solve`.

    ATTRIBUTES:

    - ``variables`` -- tuple of strings; the unknowns, in order of appearance
      in the input equation.  Every solution tuple is ordered to match.
    - ``solutions`` -- list of tuples.  When ``complete`` is ``True`` this is
      the entire solution set (possibly up to a symmetry stated in
      ``description``); otherwise it is a distinguished finite piece:
      fundamental solutions, orbit representatives, or a witness.
    - ``kind`` -- string describing the structure of the solution set:
      ``"finite-complete"`` (all solutions listed), ``"empty"``,
      ``"infinite"`` (infinitely many; iterable), ``"parametrized"``
      (infinitely many via a known parametrization; iterable when a stream is
      wired), ``"orbits"`` (finitely many orbits under a group action;
      iterable), ``"witness"`` (at least the listed ones; completeness not
      established), ``"criterion"`` (solvability criterion reported).
    - ``description`` -- human-readable structure statement (the group
      action, the parametrization, the obstruction, ...).
    - ``complete`` -- bool; whether ``solutions`` (together with the symmetry
      stated in ``description``) is provably the whole solution set.
    - ``stream`` -- callable returning a fresh iterator over the whole
      solution set, or ``None`` for a finite one.

    Iterating a :class:`SolutionSet` yields solutions: the full (possibly
    infinite) enumeration when a stream is attached, otherwise the finite
    ``solutions`` list.  Use :meth:`first` for a safe prefix.

    The kinds are checked against the object on construction, so
    ``kind="infinite"`` cannot coexist with an iteration that stops after
    the one stored witness, and ``complete=True`` cannot stand in for "the
    prose is probably enough".

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("3*x + 5*y = 1")
        sage: S.kind, S.complete
        ('infinite', True)
        sage: S.solutions[0]                     # a particular solution
        (2, -1)
        sage: next(iter(S))
        (2, -1)
    """
    variables: tuple
    solutions: list
    kind: str
    description: str = ""
    complete: bool = False
    stream: object = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        r"""
        Check that the declared kind matches what the object can deliver.

        EXAMPLES::

            sage: from diophantine_classifier.solvers import SolutionSet
            sage: SolutionSet(("x",), [(1,)], "infinite", complete=True)
            Traceback (most recent call last):
            ...
            ValueError: kind='infinite' promises an enumeration but has no stream
            sage: SolutionSet(("x",), [(1,)], "empty", complete=True)
            Traceback (most recent call last):
            ...
            ValueError: kind='empty' must have no solutions
            sage: SolutionSet(("x", "y"), [(1,)], "finite-complete",
            ....:             complete=True)
            Traceback (most recent call last):
            ...
            ValueError: solution (1,) does not match variables ('x', 'y')
        """
        self.variables = tuple(self.variables)
        for point in self.solutions:
            if len(point) != len(self.variables):
                raise ValueError(f"solution {tuple(point)} does not match "
                                 f"variables {self.variables}")
        if self.kind == "infinite" and self.stream is None:
            raise ValueError(
                "kind='infinite' promises an enumeration but has no stream")
        if self.kind == "empty":
            if self.solutions:
                raise ValueError("kind='empty' must have no solutions")
            if not self.complete:
                raise ValueError("kind='empty' must be complete")
        if self.kind == "finite-complete":
            if not self.complete:
                raise ValueError("kind='finite-complete' must be complete")
            if self.stream is not None:
                raise ValueError(
                    "kind='finite-complete' must not carry an infinite "
                    "stream")

    def as_dict(self):
        r"""
        JSON-serializable summary — the website-backend contract.

        Coordinates keep their exact values: integers serialize as integers
        and rationals as strings such as ``"1/2"``, the same policy
        :meth:`~diophantine_classifier.classify.Classification.as_dict` uses.
        Nothing is rounded and nothing is coerced through ``int``.

        Only the stored solutions are serialized; an infinite set is
        described by ``kind``, ``complete`` and ``description``, never by
        consuming its stream.

        OUTPUT: dict with plain types only

        EXAMPLES::

            sage: from diophantine_classifier import solve
            sage: import json
            sage: d = solve("2*x - 1 = 0", domain="QQ").as_dict()
            sage: d["solutions"]
            [['1/2']]
            sage: solve("x^2 - 5*x + 6 = 0").as_dict()["solutions"]
            [[2], [3]]
            sage: _ = json.dumps(d)
        """
        from .classify import _jsonify
        return {
            "kind": self.kind,
            "complete": bool(self.complete),
            "variables": list(self.variables),
            "solutions": [[_jsonify(x) for x in point]
                          for point in self.solutions],
            "description": self.description,
            "infinite": self.stream is not None,
        }

    def __iter__(self):
        r"""
        Iterate over solutions (indefinitely, for infinite families).

        EXAMPLES::

            sage: from diophantine_classifier import solve
            sage: it = iter(solve("3*x + 5*y = 1"))
            sage: [next(it) for _ in range(3)]
            [(2, -1), (-3, 2), (7, -4)]
        """
        if self.stream is not None:
            return self.stream()
        return iter(self.solutions)

    def first(self, count):
        r"""
        Return the first ``count`` solutions of the enumeration.

        INPUT:

        - ``count`` -- nonnegative integer

        OUTPUT: list of at most ``count`` solution tuples

        EXAMPLES::

            sage: from diophantine_classifier import solve
            sage: solve("3*x + 5*y = 1").first(3)
            [(2, -1), (-3, 2), (7, -4)]
        """
        return list(itertools.islice(iter(self), count))

    def __repr__(self):
        r"""
        Multi-line summary showing the first few solutions.

        EXAMPLES::

            sage: from diophantine_classifier import solve
            sage: solve("x^2 - 5*x + 6 = 0")
            SolutionSet(finite-complete, vars=['x'])
              solutions: (2,), (3,)
              all roots in ZZ
        """
        head = f"SolutionSet({self.kind}, vars={list(self.variables)})"
        body = ""
        shown = self.solutions if self.stream is None else self.first(8)
        if shown:
            text = ", ".join(str(s) for s in shown[:12])
            more = ""
            if self.stream is not None:
                more = ", ..."
            elif len(self.solutions) > 12:
                more = f", ... ({len(self.solutions)} total)"
            body = f"\n  solutions: {text}{more}"
        if self.description:
            body += f"\n  {self.description}"
        return head + body


def _zz(data, key):
    r"""
    Read an integer entry from a match's (stringified) data.

    OUTPUT: a Sage integer, or ``None`` when absent or not concrete

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _zz
        sage: _zz({"D": "61"}, "D")
        61
        sage: _zz({"D": "d"}, "D") is None and _zz({}, "D") is None
        True
    """
    try:
        return ZZ(sage_eval(str(data[key])))
    except (KeyError, TypeError, ValueError, SyntaxError, NameError):
        return None


def _normalized(match):
    r"""
    The coordinate names a solver for this match works in.

    Every solver states its answer in the family's standard coordinates;
    :func:`_finalize_solution_set` pulls the tuples back to the user's
    variables through ``match.transform``.  Where the matcher found the
    equation already in standard form, these are the equation's own
    unknowns and the pull-back is the identity.

    INPUT:

    - ``match`` -- a :class:`~diophantine_classifier.matchers.Match`

    OUTPUT: tuple of strings

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: from diophantine_classifier.solvers import _normalized
        sage: _normalized(classify("5*y + 3*x = 1").primary)
        ('y', 'x')
    """
    return tuple(match.transform.normalized_variables)


def _shells(dim):
    r"""
    Enumerate ``ZZ^dim`` by increasing sup-norm, lexicographically within
    each shell.

    INPUT:

    - ``dim`` -- nonnegative integer (kept small; used for lattice cosets)

    OUTPUT: generator of integer tuples

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _shells
        sage: import itertools
        sage: list(itertools.islice(_shells(1), 5))
        [(0,), (-1,), (1,), (-2,), (2,)]
        sage: list(itertools.islice(_shells(2), 4))
        [(0, 0), (-1, -1), (-1, 0), (-1, 1)]
    """
    if dim == 0:
        yield ()
        return
    s = 0
    while True:
        if s == 0:
            yield (0,) * dim
        else:
            for c in itertools.product(range(-s, s + 1), repeat=dim):
                if max(abs(t) for t in c) == s:
                    yield c
        s += 1


def _rationals():
    r"""
    Enumerate ``QQ`` by increasing height, each rational exactly once.

    Height is ``max(|numerator|, denominator)`` on the reduced fraction, so
    every rational appears after finitely many others and none appears
    twice — which is what makes iteration over a rational solution set fair
    rather than a walk along one sublattice.

    OUTPUT: generator of elements of ``QQ``

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _rationals
        sage: import itertools
        sage: list(itertools.islice(_rationals(), 7))
        [0, 1, -1, 1/2, -1/2, 2, -2]
        sage: seen = list(itertools.islice(_rationals(), 200))
        sage: len(set(seen)) == len(seen)
        True
    """
    yield QQ(0)
    for height in itertools.count(1):
        for num in range(1, height + 1):
            for den in range(1, height + 1):
                if max(num, den) != height or gcd(num, den) != 1:
                    continue
                value = QQ(num) / QQ(den)
                yield value
                yield -value


def _rational_shells(dim):
    r"""
    Enumerate ``QQ^dim``, fairly and without repetition.

    Indexes into :func:`_rationals` and walks the index tuples in shells,
    so every tuple of rationals is reached after finitely many steps.

    INPUT:

    - ``dim`` -- nonnegative integer

    OUTPUT: generator of tuples of elements of ``QQ``

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _rational_shells
        sage: import itertools
        sage: list(itertools.islice(_rational_shells(1), 4))
        [(0,), (1,), (-1,), (1/2,)]
        sage: first = list(itertools.islice(_rational_shells(2), 30))
        sage: first[0], len(set(first)) == len(first)
        ((0, 0), True)
    """
    if dim == 0:
        yield ()
        return
    values = []
    source = _rationals()
    for shell in itertools.count(0):
        while len(values) <= shell:
            values.append(next(source))
        for index in itertools.product(range(shell + 1), repeat=dim):
            if max(index) == shell:
                yield tuple(values[i] for i in index)


# --------------------------------------------------------------------------
# solvers
# --------------------------------------------------------------------------

def _solve_univariate(cls, match):
    r"""
    Roots of a one-variable polynomial equation in the requested domain.

    ``ZZ`` and ``NN`` take the integer roots (``NN`` keeping the
    nonnegative ones); ``QQ`` takes the rational roots, which is a strictly
    larger set and not obtainable by filtering the integer answer.  A
    parametric coefficient ring has no roots to take, and is declined.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 - 5*x + 6 = 0").solutions
        [(2,), (3,)]
        sage: solve("2*x - 1 = 0").kind             # no integer root
        'empty'
        sage: solve("2*x - 1 = 0", domain="QQ").solutions
        [(1/2,)]
        sage: solve("x^2 - 5*x + 6 = 0", domain="NN").solutions
        [(2,), (3,)]

    A symbolic coefficient is declined rather than guessed at::

        sage: from diophantine_classifier import SolverUnavailable
        sage: try:
        ....:     solve("x^2 - k = 0", params="k")
        ....: except SolverUnavailable as err:
        ....:     print("parametric" in str(err))
        True
    """
    pe = cls.working
    if not pe.is_concrete:
        raise SolverUnavailable(
            "parametric coefficients: the univariate solver needs a concrete "
            "polynomial (symbolic root-finding is not implemented)")
    pu = pe.poly.univariate_polynomial()
    domain = cls.parsed.domain
    roots = pu.roots(QQ if domain == "QQ" else ZZ, multiplicities=False)
    if domain == "NN":
        roots = [r for r in roots if r >= 0]
    sols = [(r,) for r in sorted(roots)]
    kind = "finite-complete" if sols else "empty"
    return SolutionSet(_normalized(match), sols, kind,
                       f"all roots in {domain}", complete=True)


def _solve_linear(cls, match):
    r"""
    Solve one linear equation, by domain.

    ``ZZ`` is a particular solution plus the integer kernel lattice, and
    iteration walks that coset by increasing coefficient shells.  ``QQ`` is
    an affine subspace of ``QQ^k`` — a strictly larger set that the integer
    answer does not describe — enumerated by rational height.  ``NN`` is
    handled by :func:`_solve_linear_nn`, which is complete where it answers
    at all and declines otherwise.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("3*x + 5*y = 1")
        sage: S.solutions                       # particular solution
        [(2, -1)]
        sage: S.first(3)
        [(2, -1), (-3, 2), (7, -4)]
        sage: solve("6*x + 9*y = 5").kind       # gcd obstruction
        'empty'

    Over ``QQ`` the same equation has solutions the integer lattice misses::

        sage: S = solve("2*x + 4*y = 1", domain="QQ")
        sage: S.kind
        'infinite'
        sage: all(2*x + 4*y == 1 for x, y in S.first(5))
        True
    """
    pe = cls.working
    if not pe.is_concrete:
        raise SolverUnavailable("linear solver needs concrete coefficients")
    coeffs = [QQ(sage_eval(c)) for c in match.data["coeffs"]]
    b = QQ(sage_eval(match.data["b"]))
    names = _normalized(match)
    if cls.parsed.domain == "QQ":
        return _solve_linear_qq(names, coeffs, b)
    coeffs = [ZZ(c) for c in coeffs]
    b = ZZ(b)
    if cls.parsed.domain == "NN":
        return _solve_linear_nn(names, coeffs, b)
    g = gcd(coeffs)
    if b % g:
        return SolutionSet(names, [], "empty",
                           f"gcd{tuple(coeffs)} = {g} does not divide {b}",
                           complete=True)
    cur_g, combo = coeffs[0], [ZZ(1)]
    for a in coeffs[1:]:
        cur_g, u, v = xgcd(cur_g, a)
        combo = [c * u for c in combo] + [v]
    scale = b // cur_g
    particular = tuple(c * scale for c in combo)
    kernel = [tuple(r) for r in matrix(ZZ, [coeffs]).right_kernel_matrix().rows()]

    def stream():
        for c in _shells(len(kernel)):
            yield tuple(p + sum(ci * bi[j] for ci, bi in zip(c, kernel))
                        for j, p in enumerate(particular))

    desc = (f"general solution: {particular} + integer combinations of "
            f"{kernel}")
    # the stream is lazy in every rank: withholding it above some dimension
    # would leave kind='infinite' promising an enumeration that stops at the
    # particular solution
    return SolutionSet(names, [particular], "infinite", desc,
                       complete=True, stream=stream)


def _solve_linear_qq(names, coeffs, b):
    r"""
    The rational solutions of one linear equation, completely.

    A particular solution (one variable carries `b`, the rest vanish) plus a
    basis of the rational kernel is a complete affine parametrization; the
    stream walks it by rational height, so iteration really does enumerate
    ``QQ^d`` rather than a sublattice of it.

    INPUT:

    - ``names`` -- the coordinate names the answer is stated in
    - ``coeffs`` -- list of rational coefficients
    - ``b`` -- the rational right-hand side

    OUTPUT: a :class:`SolutionSet`

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("2*x + 4*y = 1", domain="QQ")
        sage: S.solutions
        [(1/2, 0)]
        sage: sols = S.first(6)
        sage: all(2*x + 4*y == 1 for x, y in sols)
        True
        sage: len(set(sols))
        6

    """
    # the equation matched `linear`, so it has degree 1 and some coefficient
    # is nonzero
    pivot = next(i for i, c in enumerate(coeffs) if c)
    particular = tuple(b / coeffs[pivot] if i == pivot else QQ(0)
                       for i in range(len(coeffs)))
    kernel = [tuple(QQ(t) for t in row)
              for row in matrix(QQ, [coeffs]).right_kernel_matrix().rows()]
    if not kernel:
        return SolutionSet(names, [particular], "finite-complete",
                           "the unique rational solution", complete=True)

    def stream():
        for c in _rational_shells(len(kernel)):
            yield tuple(p + sum(ci * bi[j] for ci, bi in zip(c, kernel))
                        for j, p in enumerate(particular))

    desc = (f"complete affine parametrization over QQ: {particular} + "
            f"rational combinations of {kernel}")
    return SolutionSet(names, [particular], "infinite", desc,
                       complete=True, stream=stream)


def _solve_linear_nn(names, coeffs, b):
    r"""
    The nonnegative integer solutions of one linear equation.

    Complete in the two cases it answers:

    - all coefficients of one sign, where the solution set is finite and
      bounded coordinatewise, so it can be enumerated exhaustively;
    - a one-dimensional kernel, where the integer coset is a line and the
      nonnegative part of it is an explicit range of the parameter.

    Anything else is counting lattice points in an unbounded polyhedron of
    dimension at least two; rather than reuse the unrestricted integer
    lattice and emit negative coordinates, or return a false empty set, it
    raises :class:`SolverUnavailable`.

    INPUT:

    - ``names`` -- the coordinate names the answer is stated in
    - ``coeffs`` -- list of integer coefficients
    - ``b`` -- the integer right-hand side

    OUTPUT: a :class:`SolutionSet`

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("3*x + 5*y = 1", domain="NN").kind
        'empty'
        sage: solve("3*x + 5*y = 47", domain="NN").solutions
        [(4, 7), (9, 4), (14, 1)]
        sage: S = solve("x - y = 1", domain="NN")      # infinite, rank 1
        sage: S.first(3)
        [(1, 0), (2, 1), (3, 2)]

    Higher-dimensional mixed signs are declined rather than guessed::

        sage: from diophantine_classifier import SolverUnavailable
        sage: try:
        ....:     solve("x + y - z = 1", domain="NN")
        ....: except SolverUnavailable as err:
        ....:     print("nonnegative" in str(err))
        True
    """
    signs = {1 if c > 0 else -1 for c in coeffs if c}
    if len(signs) == 1 and all(coeffs):
        sign = signs.pop()
        if sign * b < 0:
            return SolutionSet(names, [], "empty",
                               "all coefficients have one sign, so no "
                               "nonnegative combination reaches the target",
                               complete=True)
        sols = _bounded_nonnegative(coeffs, b)
        kind = "finite-complete" if sols else "empty"
        return SolutionSet(names, sols, kind,
                           "all nonnegative solutions (the coefficients "
                           "share a sign, so the set is finite)",
                           complete=True)
    kernel = matrix(ZZ, [coeffs]).right_kernel_matrix().rows()
    if len(kernel) == 1:
        return _nonnegative_line(names, coeffs, b, tuple(kernel[0]))
    raise SolverUnavailable(
        "no complete nonnegative-integer solver for this shape yet: the "
        "solutions are the lattice points of an unbounded polyhedron of "
        "dimension > 1, and the unrestricted integer lattice would emit "
        "negative coordinates")


def _bounded_nonnegative(coeffs, b):
    r"""
    Exhaustively enumerate nonnegative solutions when all signs agree.

    Each coordinate is bounded by ``|b| / |c_i|``, so a depth-first walk
    over the coordinates terminates and misses nothing.

    INPUT:

    - ``coeffs`` -- list of integer coefficients, all of one sign (zeros
      allowed only when ``b`` is reached without them)
    - ``b`` -- integer right-hand side

    OUTPUT: sorted list of integer tuples

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _bounded_nonnegative
        sage: _bounded_nonnegative([ZZ(3), ZZ(5)], ZZ(47))
        [(4, 7), (9, 4), (14, 1)]
        sage: _bounded_nonnegative([ZZ(3), ZZ(5)], ZZ(1))
        []
    """
    sign = 1 if b >= 0 else -1
    for c in coeffs:
        if c:
            sign = 1 if c > 0 else -1
            break
    target, cs = sign * b, [sign * c for c in coeffs]
    sols = []

    def walk(i, rest, acc):
        if i == len(cs):
            if rest == 0:
                sols.append(tuple(acc))
            return
        if cs[i] == 0:
            # a zero coefficient leaves its coordinate free, which would
            # make the set infinite; only 0 keeps the enumeration complete
            walk(i + 1, rest, acc + [ZZ(0)])
            return
        for value in range(rest // cs[i] + 1):
            walk(i + 1, rest - cs[i] * value, acc + [ZZ(value)])

    if target < 0:
        return []
    walk(0, target, [])
    return sorted(sols)


def _nonnegative_line(names, coeffs, b, direction):
    r"""
    The nonnegative points of a one-dimensional integer solution line.

    With a rank-one kernel the integer solutions are ``p + t*d`` for
    ``t`` in ``ZZ``; each coordinate constrains ``t`` to a half-line, so
    their intersection is an explicit (possibly infinite) range.

    INPUT:

    - ``names`` -- the coordinate names the answer is stated in
    - ``coeffs`` -- list of integer coefficients
    - ``b`` -- the integer right-hand side
    - ``direction`` -- a generator of the integer kernel

    OUTPUT: a :class:`SolutionSet`

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x - y = 1", domain="NN").first(3)
        [(1, 0), (2, 1), (3, 2)]
        sage: solve("x + y = 3", domain="NN").solutions
        [(0, 3), (1, 2), (2, 1), (3, 0)]
    """
    g = gcd(coeffs)
    if b % g:
        return SolutionSet(names, [], "empty",
                           f"gcd{tuple(coeffs)} = {g} does not divide {b}",
                           complete=True)
    cur_g, combo = coeffs[0], [ZZ(1)]
    for a in coeffs[1:]:
        cur_g, u, v = xgcd(cur_g, a)
        combo = [c * u for c in combo] + [v]
    particular = [c * (b // cur_g) for c in combo]
    lo, hi = None, None                     # the admissible range of t
    for p, d in zip(particular, direction):
        if d == 0:
            if p < 0:
                return SolutionSet(names, [], "empty",
                                   "a fixed coordinate of the solution line "
                                   "is negative", complete=True)
            continue
        bound = QQ(-p) / QQ(d)
        if d > 0:
            start = bound.ceil()
            lo = start if lo is None else max(lo, start)
        else:
            stop = bound.floor()
            hi = stop if hi is None else min(hi, stop)
    if lo is not None and hi is not None and lo > hi:
        return SolutionSet(names, [], "empty",
                           "the solution line misses the nonnegative "
                           "orthant", complete=True)

    def point(t):
        return tuple(p + t * d for p, d in zip(particular, direction))

    if lo is not None and hi is not None:
        sols = sorted(point(t) for t in range(lo, hi + 1))
        kind = "finite-complete" if sols else "empty"
        return SolutionSet(names, sols, kind,
                           "all nonnegative solutions on the solution line",
                           complete=True)
    step = 1 if hi is None else -1
    start = lo if hi is None else hi

    def stream():
        t = start
        while True:
            yield point(t)
            t += step

    return SolutionSet(names, [point(start)], "infinite",
                       f"the nonnegative part of the solution line "
                       f"{tuple(particular)} + t*{tuple(direction)}",
                       complete=True, stream=stream)


def _pell_unit(D):
    r"""
    Fundamental solution of ``x^2 - D*y^2 = ±1`` by continued fractions.

    OUTPUT: triple ``(x1, y1, norm)`` with ``x1^2 - D*y1^2 = norm ∈ {1, -1}``

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _pell_unit
        sage: _pell_unit(61)
        (29718, 3805, -1)
        sage: _pell_unit(3)
        (2, 1, 1)
    """
    cf = continued_fraction(QuadraticField(D).gen())
    ell = len(cf.period())
    conv = cf.convergent(ell - 1)
    x1, y1 = conv.numerator(), conv.denominator()
    return x1, y1, x1 ** 2 - D * y1 ** 2


def _signed_orbit_stream(D, fund, unit, include_trivial):
    r"""
    Stream all solutions of ``x^2 - D*y^2 = ±1`` from the fundamental one.

    INPUT:

    - ``D`` -- the Pell parameter
    - ``fund`` -- fundamental solution of the target equation
    - ``unit`` -- fundamental solution ``(t, u)`` of the ``+1`` equation
    - ``include_trivial`` -- whether ``(±1, 0)`` are solutions (the ``+1``
      case)

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _signed_orbit_stream
        sage: s = _signed_orbit_stream(2, (3, 2), (3, 2), True)()
        sage: [next(s) for _ in range(6)]
        [(1, 0), (-1, 0), (3, 2), (3, -2), (-3, 2), (-3, -2)]
    """
    t, u = unit

    def stream():
        if include_trivial:
            yield (ZZ(1), ZZ(0))
            yield (ZZ(-1), ZZ(0))
        x, y = fund
        while True:
            yield (x, y)
            yield (x, -y)
            yield (-x, y)
            yield (-x, -y)
            x, y = t * x + D * u * y, u * x + t * y
    return stream


def _solve_pell(cls, match):
    r"""
    Solve ``x^2 - D*y^2 = ±1``: fundamental solution + full enumeration.

    ``solutions[0]`` is the fundamental solution; iteration enumerates all
    integer solutions ordered by the power of the fundamental unit, with sign
    pattern ``(x, y), (x, -y), (-x, y), (-x, -y)``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x^2 - 61*y^2 = 1")
        sage: S.solutions[0]
        (1766319049, 226153980)
        sage: S = solve("x^2 - 2*y^2 = -1")     # negative Pell
        sage: S.first(3)
        [(1, 1), (1, -1), (-1, 1)]
        sage: solve("x^2 - 3*y^2 = -1").kind    # no negative Pell for D = 3
        'empty'
    """
    D = _zz(match.data, "D")
    N = _zz(match.data, "N")
    if D is None or N not in (1, -1):
        raise SolverUnavailable("Pell solver needs concrete D and N = ±1")
    x1, y1, norm = _pell_unit(D)
    if N == 1:
        if norm == 1:
            fund = (x1, y1)
        else:
            fund = (x1 ** 2 + D * y1 ** 2, 2 * x1 * y1)
        desc = (f"infinitely many: ±(fundamental)^k for the fundamental "
                f"solution {fund}; iteration enumerates them all")
        return SolutionSet(_normalized(match), [fund], "infinite", desc,
                           complete=True,
                           stream=_signed_orbit_stream(D, fund, fund, True))
    if norm == -1:
        unit = (x1 ** 2 + D * y1 ** 2, 2 * x1 * y1)
        desc = (f"infinitely many: odd powers of the fundamental unit; "
                f"fundamental solution {(x1, y1)}")
        return SolutionSet(_normalized(match), [(x1, y1)], "infinite", desc,
                           complete=True,
                           stream=_signed_orbit_stream(D, (x1, y1), unit,
                                                       False))
    return SolutionSet(_normalized(match), [], "empty",
                       f"x^2 - {D}y^2 = -1 has no solutions (continued "
                       "fraction period is even)", complete=True)


def _solve_pell_like(cls, match):
    r"""
    Solve ``x^2 - D*y^2 = N``: orbit representatives + full enumeration.

    PARI's ``qfbsolve`` provides representatives of the finitely many orbits
    under the automorph group; iteration walks the orbits outward by
    applying the fundamental automorphism in both directions.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x^2 - 2*y^2 = 7")
        sage: sols = S.first(8)
        sage: all(x^2 - 2*y^2 == 7 for x, y in sols)
        True
        sage: len(set(sols))
        8
    """
    D = _zz(match.data, "D")
    N = _zz(match.data, "N")
    if D is None or N is None:
        raise SolverUnavailable("generalized Pell solver needs concrete D, N")
    if N in (1, -1):
        return _solve_pell(cls, match)
    try:
        res = pari(f"qfbsolve(Qfb(1,0,{-D}),{N},1)")
        reps = [(ZZ(v[0]), ZZ(v[1])) for v in res]
    except Exception as err:
        raise SolverUnavailable(f"PARI qfbsolve failed: {err}") from None
    reps = [s for s in reps if s[0] ** 2 - D * s[1] ** 2 == N]
    if not reps:
        return SolutionSet(_normalized(match), [], "empty", "no solutions",
                           complete=True)
    x1, y1, norm = _pell_unit(D)
    if norm == -1:
        t, u = x1 ** 2 + D * y1 ** 2, 2 * x1 * y1
    else:
        t, u = x1, y1
    seeds = set()
    for x, y in reps:
        seeds.update({(x, y), (x, -y), (-x, y), (-x, -y)})

    def key(s):
        return (max(abs(s[0]), abs(s[1])), s)

    def stream():
        seen = set()
        level = sorted(seeds, key=key)
        while level:
            nxt = []
            for s in level:
                if s in seen:
                    continue
                seen.add(s)
                yield s
                x, y = s
                nxt.append((t * x + D * u * y, u * x + t * y))
                nxt.append((t * x - D * u * y, -u * x + t * y))
            level = sorted(set(nxt) - seen, key=key)

    desc = (f"{len(reps)} orbit representative(s) under the automorph group "
            f"(fundamental automorphism {(t, u)}); iteration enumerates the "
            "full orbits")
    return SolutionSet(_normalized(match), sorted(reps), "orbits", desc,
                       complete=True, stream=stream)


def _solve_bqf(cls, match):
    r"""
    Representations by a binary quadratic form.

    Definite forms: the complete (finite) list of representations.
    Indefinite forms: a witness via ``BinaryQF.solve_integer``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("3*x^2 + 7*y^2 = 19").solutions
        [(-2, -1), (-2, 1), (2, -1), (2, 1)]
        sage: solve("3*x^2 + 7*y^2 = 5").kind
        'empty'
    """
    a, b, c, n = (_zz(match.data, k) for k in ("a", "b", "c", "n"))
    if None in (a, b, c, n):
        raise SolverUnavailable("needs concrete form and n")
    disc = b ** 2 - 4 * a * c
    if disc < 0 and a > 0 and abs(n) <= MAX_BQF:
        if n < 0:
            return SolutionSet(_normalized(match), [], "empty",
                               "positive definite form cannot represent a "
                               "negative integer", complete=True)
        sols = []
        Y = isqrt(4 * a * n // (-disc)) + 1
        for y in range(-Y, Y + 1):
            discx = (b * y) ** 2 - 4 * a * (c * y * y - n)
            if discx < 0:
                continue
            s = isqrt(discx)
            if s * s != discx:
                continue
            for sgn in ((s,) if s == 0 else (s, -s)):
                num = -b * y + sgn
                if num % (2 * a) == 0:
                    sols.append((ZZ(num // (2 * a)), ZZ(y)))
        sols = sorted(set(sols))
        kind = "finite-complete" if sols else "empty"
        return SolutionSet(_normalized(match), sols, kind,
                           "all representations (definite form)",
                           complete=True)
    form = BinaryQF([a, b, c])
    sol = form.solve_integer(n)
    if sol is None:
        return SolutionSet(_normalized(match), [], "empty",
                           "no representation", complete=disc < 0)
    return SolutionSet(_normalized(match), [tuple(sol)], "witness",
                       "one representation (BinaryQF.solve_integer); for "
                       "indefinite forms the full set is a union of "
                       "automorph orbits", complete=False)


def _solve_qf_zero(cls, match, gram=None):
    r"""
    Nontrivial zero of a quadratic form, or the local obstruction.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x^2 + y^2 = 2*z^2")
        sage: S.solutions[0]
        (1, 1, -1)
        sage: solve("x^2 + y^2 = 3*z^2").kind
        'empty'
    """
    if gram is None:
        gram = sage_eval(str(match.data["gram"]))
    G = matrix(QQ, gram)
    try:
        res = qfsolve(G)
    except Exception as err:
        raise SolverUnavailable(f"qfsolve failed: {err}") from None
    if res in ZZ:
        place = "the real place" if res == -1 else f"p = {res}"
        return SolutionSet(_normalized(match), [], "empty",
                           f"no nontrivial solutions: local obstruction at "
                           f"{place}", complete=True)
    vec = [QQ(t) for t in res]
    den = lcm([t.denominator() for t in vec])
    ivec = [ZZ(t * den) for t in vec]
    g = gcd(ivec)
    ivec = tuple(t // g for t in ivec)
    if sum(1 for t in ivec if t < 0) > sum(1 for t in ivec if t > 0):
        ivec = tuple(-t for t in ivec)
    return SolutionSet(_normalized(match), [ivec], "parametrized",
                       "one nontrivial solution; all others arise from it by "
                       "the standard conic/quadric parametrization",
                       complete=False)


def _solve_legendre(cls, match):
    r"""
    Legendre equation ``a x^2 + b y^2 + c z^2 = 0`` via ``qfsolve``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: x, y, z = solve("x^2 + 3*y^2 = 7*z^2").solutions[0]
        sage: x^2 + 3*y^2 == 7*z^2
        True
    """
    a, b, c = (_zz(match.data, k) for k in ("a", "b", "c"))
    if None in (a, b, c):
        raise SolverUnavailable("needs concrete coefficients")
    return _solve_qf_zero(cls, match,
                          gram=[[a, 0, 0], [0, b, 0], [0, 0, c]])


def _solve_weierstrass(cls, match):
    r"""
    Integral points on a Weierstrass model via ``E.integral_points``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: sorted(solve("y^2 = x^3 - 2").solutions)   # unknowns (y, x)
        [(-5, 3), (5, 3)]
    """
    if not cls.working.is_concrete:
        raise SolverUnavailable("needs concrete coefficients")
    ainvs = sage_eval(str(match.data["ainvs"])) if "ainvs" in match.data \
        else [0, 0, 0, 0, sage_eval(str(match.data["k"]))]
    E = EllipticCurve(QQ, [QQ(t) for t in ainvs])
    pts = E.integral_points(both_signs=True)
    sols = [(P[0], P[1]) for P in pts]
    return SolutionSet(
        _normalized(match), sorted(sols), "finite-complete",
        f"all integral points on {E.ainvs()} (rank {E.rank()}); rational "
        "points are infinite iff the rank is positive", complete=True)


def _solve_thue(cls, match):
    r"""
    Thue equation via PARI's certified ``thue`` solver.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^3 + 2*y^3 = 11").solutions
        [(3, -2)]
        sage: solve("x^4 - 2*y^4 = 1").solutions
        [(-1, 0), (1, 0)]
    """
    pe = cls.working
    P = pe.poly
    x, y = pe.poly_ring.gens()
    m = -P.constant_coefficient()
    F = P + m
    fu = F.subs({y: 1}).univariate_polynomial().change_variable_name("X")
    try:
        res = pari(f"thue(thueinit({fu},1),{m})")
        sols = [(ZZ(v[0]), ZZ(v[1])) for v in res]
    except Exception as err:
        raise SolverUnavailable(f"PARI thue failed: {err}") from None
    sols = [s for s in sols if F.subs({x: s[0], y: s[1]}) == m]
    return SolutionSet(_normalized(match), sorted(sols), "finite-complete",
                       "all solutions (PARI thue, certified)", complete=True)


def _solve_egyptian(cls, match):
    r"""
    Unit fraction equations ``1/x_1 + ... + 1/x_k = a/n``, concrete case.

    Enumerates all solutions in positive integers with
    ``x_1 <= x_2 <= ... <= x_k`` by branch-and-bound.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("1/x + 1/y + 1/z = 1").solutions
        [(2, 3, 6), (2, 4, 4), (3, 3, 3)]
        sage: solve("4/5 = 1/x + 1/y + 1/z").solutions[:2]
        [(2, 4, 20), (2, 5, 10)]
    """
    a = _zz(match.data, "a")
    n = _zz(match.data, "n")
    k = _zz(match.data, "k")
    if k is None and match.slug == "erdos-straus":
        a, k = ZZ(4), ZZ(3)
    if a is None and match.slug == "erdos-straus":
        a = ZZ(4)
    if None in (a, n, k):
        raise SolverUnavailable(
            "unit-fraction enumeration needs concrete a and n (parametric "
            "n is the open conjecture territory)")
    if k > 5 or n > MAX_EGYPTIAN_N:
        raise SolverUnavailable("enumeration bound exceeded (k <= 5, "
                                f"n <= {MAX_EGYPTIAN_N})")
    sols = []

    def rec(k_left, target, minimum, acc):
        if target <= 0:
            return
        if k_left == 1:
            if target.numerator() == 1 and target.denominator() >= minimum:
                sols.append(acc + (ZZ(target.denominator()),))
            return
        lo = max(minimum, (QQ(1) / target).floor() + 1)
        hi = (QQ(k_left) / target).floor()
        for x in range(lo, hi + 1):
            rec(k_left - 1, target - QQ(1) / x, x, acc + (ZZ(x),))

    rec(int(k), QQ(a) / QQ(n), 1, ())
    kind = "finite-complete" if sols else "empty"
    return SolutionSet(
        _normalized(match), sorted(sols), kind,
        "all solutions in positive integers with x_1 <= ... <= x_k; "
        "permutations give the rest", complete=True)


def _solve_ramanujan_nagell(cls, match):
    r"""
    The classical Ramanujan-Nagell equation ``x^2 + 7 = 2^n``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 + 7 = 2^n").solutions
        [(1, 3), (3, 4), (5, 5), (11, 7), (181, 15)]
    """
    d, k, base = _zz(match.data, "d"), _zz(match.data, "k"), _zz(match.data, "base")
    if (d, k, base) != (7, 1, 2):
        raise SolverUnavailable(
            "only the classical x^2 + 7 = 2^n is hardwired; general (d, k, b) "
            "need a Baker + LLL computation (Petho-de Weger)")
    pairs = [(ZZ(x), ZZ(n))
             for x, n in [(1, 3), (3, 4), (5, 5), (11, 7), (181, 15)]]
    return SolutionSet(
        _normalized(match), pairs, "finite-complete",
        "Nagell's theorem: n ∈ {3, 4, 5, 7, 15} (x > 0 shown; -x symmetric)",
        complete=True)


#: dispatch table: family slug -> solver function
SOLVERS = {
    "univariate": _solve_univariate,
    "linear": _solve_linear,
    "pell-like": _solve_pell_like,
    "binary-qf-representation": _solve_bqf,
    "quadratic-form-zero": _solve_qf_zero,
    "legendre": _solve_legendre,
    "elliptic-weierstrass": _solve_weierstrass,
    "thue": _solve_thue,
    "egyptian-fractions": _solve_egyptian,
    "ramanujan-nagell": _solve_ramanujan_nagell,
}


#: the domains each solver is valid for.  A solver absent from this table
#: answers over ``ZZ`` only; ``solve`` declines every other domain rather
#: than hand back the integer answer to a rational question.  ``NN`` appears
#: wherever the ``ZZ`` answer is complete, since the nonnegative solutions
#: are exactly the integer ones that survive the domain filter.
SOLVER_DOMAINS = {
    "univariate": ("ZZ", "NN", "QQ"),
    "linear": ("ZZ", "NN", "QQ"),
    "pell": ("ZZ", "NN"),
    "pell-like": ("ZZ", "NN"),
    "sum-of-two-squares": ("ZZ", "NN"),
    "sum-of-three-squares": ("ZZ", "NN"),
    "sum-of-four-squares": ("ZZ", "NN"),
    "binary-qf-representation": ("ZZ", "NN"),
    "pythagorean": ("ZZ", "NN"),
    "mordell": ("ZZ", "NN"),
    "elliptic-weierstrass": ("ZZ", "NN"),
    "thue": ("ZZ", "NN"),
    "markov-hurwitz": ("ZZ", "NN"),
    "egyptian-fractions": ("ZZ", "NN"),
    "erdos-straus": ("ZZ", "NN"),
    "catalan": ("ZZ", "NN"),
    "fermat": ("ZZ", "NN"),
    "ramanujan-nagell": ("ZZ", "NN"),
}

#: domains a solver may be assumed to handle when it is not listed above
DEFAULT_DOMAINS = ("ZZ",)


def _in_domain(value, domain):
    r"""
    Whether a coordinate lies in the requested domain.

    INPUT:

    - ``value`` -- a Sage number
    - ``domain`` -- ``"ZZ"``, ``"NN"`` or ``"QQ"``

    OUTPUT: boolean

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _in_domain
        sage: _in_domain(QQ(1)/2, "QQ"), _in_domain(QQ(1)/2, "ZZ")
        (True, False)
        sage: _in_domain(ZZ(-3), "ZZ"), _in_domain(ZZ(-3), "NN")
        (True, False)
    """
    try:
        if domain == "QQ":
            return value in QQ
        if value not in ZZ:
            return False
    except TypeError:
        return False
    return domain != "NN" or value >= 0


def _nonnegative_shells(dim):
    r"""
    Enumerate ``NN^dim`` by increasing maximum, lexicographically within
    each shell.

    INPUT:

    - ``dim`` -- nonnegative integer

    OUTPUT: generator of tuples of nonnegative Sage integers

    EXAMPLES::

        sage: from diophantine_classifier.solvers import _nonnegative_shells
        sage: import itertools
        sage: list(itertools.islice(_nonnegative_shells(1), 3))
        [(0,), (1,), (2,)]
        sage: list(itertools.islice(_nonnegative_shells(2), 4))
        [(0, 0), (0, 1), (1, 0), (1, 1)]
    """
    if dim == 0:
        yield ()
        return
    for shell in itertools.count(0):
        for point in itertools.product(range(shell + 1), repeat=dim):
            if max(point) == shell:
                yield tuple(ZZ(t) for t in point)


#: how each domain is enumerated when the conditions are the whole problem
DOMAIN_SHELLS = {"ZZ": _shells, "NN": _nonnegative_shells,
                 "QQ": _rational_shells}

#: candidates to scan for a witness before admitting we found none
MAX_CONSTRAINT_SEARCH = 10**4


def _solve_conditional_identity(cls):
    r"""
    Solve an equality that holds identically on its condition locus.

    There is no equation to invert here: the solution set is every point of
    the requested domain satisfying the side conditions.  The domain is
    enumerated fairly — integer shells, nonnegative shells, or rationals by
    height — and filtered through
    :meth:`~diophantine_classifier.parsing.ParsedEquation.is_solution`, so
    the stream really does enumerate the whole set.

    INPUT:

    - ``cls`` -- a :class:`~diophantine_classifier.classify.Classification`
      of a conditional identity

    OUTPUT: a :class:`SolutionSet`

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x/x = 1")
        sage: S.kind, S.complete
        ('infinite', True)
        sage: S.first(4)
        [(-1,), (1,), (-2,), (2,)]
        sage: solve("x/x = 1", domain="NN").first(3)
        [(1,), (2,), (3,)]
        sage: solve("x/x = 1", domain="QQ").first(3)
        [(1,), (-1,), (1/2,)]
    """
    pe = cls.parsed
    unknowns = pe.unknowns
    shells = DOMAIN_SHELLS[pe.domain]

    def stream():
        for point in shells(len(unknowns)):
            if pe.is_solution(dict(zip(unknowns, point))):
                yield tuple(point)

    witness = list(itertools.islice(
        (point for point in itertools.islice(shells(len(unknowns)),
                                             MAX_CONSTRAINT_SEARCH)
         if pe.is_solution(dict(zip(unknowns, point)))), 1))
    if not witness:
        raise SolverUnavailable(
            f"no assignment satisfying {'; '.join(str(c) for c in pe.conditions)}"
            f" found among the first {MAX_CONSTRAINT_SEARCH} points of "
            f"{pe.domain}^{len(unknowns)}; the conditions may be "
            "unsatisfiable there")
    return SolutionSet(
        unknowns, [tuple(witness[0])], "infinite",
        "every assignment in the domain satisfying "
        + "; ".join(str(c) for c in pe.conditions),
        complete=True, stream=stream)


def _finalize_solution_set(cls, match, result):
    r"""
    Turn a solver's answer into an answer to the question that was asked.

    A solver works in its family's standard coordinates on the working
    model.  Before anyone sees its output, every tuple is

    1. pulled back through ``match.transform`` to the user's variables,
    2. ordered by the original equation's unknowns,
    3. checked against the requested domain, and
    4. checked against the original equation and its side conditions,
       exactly.

    A tuple is dropped only when it is *refuted*: an assignment that leaves
    a parameter undetermined is kept, because nothing has been shown about
    it.  Dropping only refuted tuples never removes a genuine solution, so a
    ``complete`` result stays complete, and infinite families are filtered
    lazily, one solution at a time.

    INPUT:

    - ``cls`` -- the :class:`~diophantine_classifier.classify.Classification`
    - ``match`` -- the :class:`~diophantine_classifier.matchers.Match` that
      was solved
    - ``result`` -- the :class:`SolutionSet` the solver returned

    OUTPUT: a :class:`SolutionSet` in the source variables

    EXAMPLES:

    Clearing the denominators of ``1/(x-1) = 1/(y-1)`` gives ``x = y``, whose
    solution ``(1, 1)`` is a pole of the original equation::

        sage: from diophantine_classifier import solve
        sage: S = solve("1/(x - 1) = 1/(y - 1)")
        sage: (1, 1) in S.first(30)
        False
        sage: all(x != 1 and y != 1 for x, y in S.first(30))
        True

    Coordinates come back in the order of the original equation, whatever
    the standard model called them::

        sage: solve("5*y + 3*x = 1").variables
        ('y', 'x')
        sage: solve("5*y + 3*x = 1").first(2)
        [(-1, 2), (-4, 7)]
    """
    pe = cls.parsed
    # the whole map back to the submitted problem, reduction included
    transform = cls.effective_transform(match)
    unknowns = pe.unknowns
    missing = [v for v in unknowns if v not in transform.to_source]
    if missing:
        raise SolverUnavailable(
            f"cannot state a solution for {missing[0]!r}: the normalization "
            "does not invert onto every unknown of the original problem")
    names = tuple(result.variables)
    if names != tuple(transform.normalized_variables):
        raise AssertionError(
            f"solver for {match.slug!r} returned coordinates {names}, not "
            f"the standard {tuple(transform.normalized_variables)}")

    def transported(point):
        source = transform.pull_back(dict(zip(names, point)))
        return tuple(source[v] for v in unknowns)

    def verdict(point):
        """The transported tuple, and why it was dropped (or ``None``).

        Order matters.  A point outside the domain, or on a pole, or off
        the locus where the normalization is defined, is *expected* to be
        dropped; only a point that clears all of those and still fails the
        equation says anything about the map.  Checking the equation first
        would report the zero vector of a homogeneous linear model -- the
        usual particular solution, and a pole of 1/x = 1/y -- as evidence
        of a bad normalization.
        """
        ordered = transported(point)
        assignment = dict(zip(unknowns, ordered))
        if any(not _in_domain(assignment[v], pe.domain) for v in unknowns):
            return ordered, "domain"
        if pe.conditions_hold(assignment) is False:
            return ordered, "source-condition"
        if transform.conditions_hold(assignment) is False:
            return ordered, "transform-condition"
        if pe.satisfies_original(assignment) is False:
            return ordered, "equation"
        return ordered, None

    def check(point):
        ordered, reason = verdict(point)
        return None if reason else ordered

    judged = [verdict(point) for point in result.solutions]
    solutions = [ordered for ordered, reason in judged if not reason]
    if any(reason == "equation" for _, reason in judged) and not solutions:
        # a point that is in the domain and satisfies every known condition
        # and *still* fails the equation is not a filtered point: the
        # normalization described a different problem, and a lazily filtered
        # stream would then spin forever instead of saying so
        raise AssertionError(
            f"solver for {match.slug!r} produced no solution of the original "
            f"equation {pe.original!r}: the normalization recorded in "
            f"{str(transform) or 'the identity transform'} does not invert "
            "onto this problem")

    stream = None
    if result.stream is not None:
        source = result.stream

        def stream():
            return (ordered for ordered in map(check, source())
                    if ordered is not None)

    kind = result.kind
    # an infinite set keeps its stream: losing the stored witness to a pole
    # says nothing about the rest of the enumeration
    if kind == "finite-complete" and not solutions:
        kind = "empty"
    return replace(result, variables=unknowns, solutions=solutions,
                   kind=kind, stream=stream)


def solve(equation, params=(), domain="ZZ"):
    r"""
    Solve (or partially solve) an equation via its classification.

    INPUT:

    - ``equation`` -- string, or a
      :class:`~diophantine_classifier.classify.Classification`
    - ``params`` -- (default: ``()``) parameter names, as for
      :func:`~diophantine_classifier.classify.classify`
    - ``domain`` -- (default: ``"ZZ"``) where solutions live

    OUTPUT: a :class:`SolutionSet`

    Dispatch walks the matches that were actually emitted, most specific
    first, and hands each solver *its own* match — so a fallback solver
    receives the data of the model it is solving, never the primary match's.
    A specialization edge is not evidence that one family's data satisfies
    another's input contract, so an ancestor that was not itself matched is
    not tried.

    A solver is used only on a domain it declares in :data:`SOLVER_DOMAINS`.
    When nothing is wired for the family, or nothing covers the requested
    domain, the raised :class:`SolverUnavailable` carries software pointers
    and filled code templates from the registry.

    EXAMPLES::

        sage: from diophantine_classifier import classify, solve, SolverUnavailable
        sage: solve("x^2 - 5*x + 6 = 0").solutions
        [(2,), (3,)]
        sage: solve(classify("3*x + 5*y = 1")).kind
        'infinite'
        sage: try:
        ....:     solve("y^2 = x^7 + 3")
        ....: except SolverUnavailable as err:
        ....:     print("no automatic solver" in str(err))
        True

    Each solver declares the domains it answers over, so an integer solver
    is never quietly reused for a rational question::

        sage: from diophantine_classifier.solvers import SOLVER_DOMAINS
        sage: SOLVER_DOMAINS["linear"], SOLVER_DOMAINS["pell"]
        (('ZZ', 'NN', 'QQ'), ('ZZ', 'NN'))
    """
    cls = equation if isinstance(equation, Classification) \
        else classify(equation, params=params, domain=domain)
    if cls.is_composite:
        raise SolverUnavailable(
            "reducible equation: solve each component separately "
            f"({', '.join(c.slug for c in cls.components)})")
    if cls.special_kind == "conditional-identity":
        return _solve_conditional_identity(cls)
    if cls.free_variables:
        # the factor constrains only some of the ambient unknowns; its own
        # solver would answer in fewer coordinates, and presenting that as
        # the component's solution set would drop a whole dimension
        raise SolverUnavailable(
            f"component {cls.embedding.factor} = 0 leaves "
            f"{', '.join(cls.free_variables)} free in the ambient problem: "
            "lifting an active solution set over the free coordinates is not "
            "implemented, and the active answer alone is not this "
            "component's solution set")
    wrong_domain = []
    for match in cls.matches:
        fn = SOLVERS.get(match.slug)
        if fn is None:
            continue
        if cls.parsed.domain not in SOLVER_DOMAINS.get(match.slug,
                                                       DEFAULT_DOMAINS):
            wrong_domain.append(match.slug)
            continue
        return _finalize_solution_set(cls, match, fn(cls, match))
    if wrong_domain:
        raise SolverUnavailable(
            f"no solver for {', '.join(sorted(wrong_domain))} over "
            f"{cls.parsed.domain}: the wired solvers answer over "
            f"{', '.join(SOLVER_DOMAINS.get(wrong_domain[0], DEFAULT_DOMAINS))}"
            " only, and the integer answer does not describe the "
            f"{cls.parsed.domain} solution set")
    fam = cls.family
    hints = []
    if fam is not None:
        for lang, what in fam.software.items():
            hints.append(f"{lang}: {what}")
        for lang, snippet in cls.code().items():
            hints.append(f"code[{lang}]: {snippet.strip()}")
    hint = ("; ".join(hints)) if hints else "no standard software known"
    raise SolverUnavailable(
        f"no automatic solver for family {cls.slug!r} yet — {hint}")
