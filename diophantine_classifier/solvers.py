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
    sage: solve("x^3 + 2*y^3 = 11").solutions        # finite: complete list
    [(3, -2)]
    sage: S = solve("x^2 - 2*y^2 = 1")               # infinite: iterate
    sage: S.first(6)
    [(1, 0), (-1, 0), (3, 2), (3, -2), (-3, 2), (-3, -2)]
    sage: S = solve("x^2 + y^2 + z^2 = 3*x*y*z")     # Markov tree
    sage: S.first(4)
    [(1, 1, 1), (1, 1, 2), (1, 2, 5), (1, 5, 13)]
"""

import heapq
import itertools
from dataclasses import dataclass, field

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
        ....:     print("hyperelliptic" in str(err))
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

    Iterating a :class:`SolutionSet` yields solutions: the full (possibly
    infinite) enumeration when a stream is attached, otherwise the finite
    ``solutions`` list.  Use :meth:`first` for a safe prefix.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x^2 - 61*y^2 = 1")
        sage: S.kind, S.complete
        ('infinite', True)
        sage: S.solutions[0]                     # the fundamental solution
        (1766319049, 226153980)
        sage: next(iter(S))
        (1, 0)
    """
    variables: tuple
    solutions: list
    kind: str
    description: str = ""
    complete: bool = False
    stream: object = field(default=None, repr=False, compare=False)

    def __iter__(self):
        r"""
        Iterate over solutions (indefinitely, for infinite families).

        EXAMPLES::

            sage: from diophantine_classifier import solve
            sage: it = iter(solve("x^2 - 2*y^2 = 1"))
            sage: [next(it) for _ in range(3)]
            [(1, 0), (-1, 0), (3, 2)]
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
              all roots in the domain
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


def _ordered(cls, assignment):
    r"""
    Order an assignment dict into a tuple matching the equation's unknowns.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: from diophantine_classifier.solvers import _ordered
        sage: cls = classify("y^2 = x^3 - 2")     # unknowns appear as (y, x)
        sage: _ordered(cls, {"x": 3, "y": 5})
        (5, 3)
    """
    return tuple(assignment[v] for v in cls.parsed.unknowns)


def _roles(cls):
    r"""
    The role mapping recorded by the matcher, if any.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: from diophantine_classifier.solvers import _roles
        sage: _roles(classify("x^p - y^q = 1"))["x"]
        'x'
    """
    return cls.data.get("roles", {})


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


# --------------------------------------------------------------------------
# solvers
# --------------------------------------------------------------------------

def _solve_univariate(cls):
    r"""
    Roots of a one-variable polynomial equation in the given domain.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 - 5*x + 6 = 0").solutions
        [(2,), (3,)]
        sage: solve("x^2 - 5*x + 6 = 0", domain="QQ").solutions
        [(2,), (3,)]
    """
    pu = cls.parsed.poly.univariate_polynomial()
    ring = QQ if cls.parsed.domain == "QQ" else ZZ
    roots = pu.roots(ring, multiplicities=False)
    if cls.parsed.domain == "NN":
        roots = [r for r in roots if r >= 0]
    return SolutionSet(cls.parsed.unknowns, [(r,) for r in sorted(roots)],
                       "finite-complete", "all roots in the domain",
                       complete=True)


def _solve_linear(cls):
    r"""
    Solve a linear equation: particular solution plus solution lattice.

    Iteration enumerates the coset ``particular + lattice`` by increasing
    coefficient shells.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("3*x + 5*y = 1")
        sage: S.solutions                       # particular solution
        [(2, -1)]
        sage: S.first(3)
        [(2, -1), (-3, 2), (7, -4)]
        sage: solve("6*x + 9*y = 5").kind       # gcd obstruction
        'empty'
    """
    if not cls.parsed.is_concrete:
        raise SolverUnavailable("linear solver needs concrete coefficients")
    coeffs = [ZZ(sage_eval(c)) for c in cls.data["coeffs"]]
    b = ZZ(sage_eval(cls.data["b"]))
    g = gcd(coeffs)
    if b % g:
        return SolutionSet(cls.parsed.unknowns, [], "empty",
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
    return SolutionSet(cls.parsed.unknowns, [particular], "infinite", desc,
                       complete=True,
                       stream=stream if len(kernel) <= 8 else None)


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


def _solve_pell(cls):
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
    D = _zz(cls.data, "D")
    N = _zz(cls.data, "N")
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
        return SolutionSet(cls.parsed.unknowns, [fund], "infinite", desc,
                           complete=True,
                           stream=_signed_orbit_stream(D, fund, fund, True))
    if norm == -1:
        unit = (x1 ** 2 + D * y1 ** 2, 2 * x1 * y1)
        desc = (f"infinitely many: odd powers of the fundamental unit; "
                f"fundamental solution {(x1, y1)}")
        return SolutionSet(cls.parsed.unknowns, [(x1, y1)], "infinite", desc,
                           complete=True,
                           stream=_signed_orbit_stream(D, (x1, y1), unit,
                                                       False))
    return SolutionSet(cls.parsed.unknowns, [], "empty",
                       f"x^2 - {D}y^2 = -1 has no solutions (continued "
                       "fraction period is even)", complete=True)


def _solve_pell_like(cls):
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
    D = _zz(cls.data, "D")
    N = _zz(cls.data, "N")
    if D is None or N is None:
        raise SolverUnavailable("generalized Pell solver needs concrete D, N")
    if N in (1, -1):
        return _solve_pell(cls)
    try:
        res = pari(f"qfbsolve(Qfb(1,0,{-D}),{N},1)")
        reps = [(ZZ(v[0]), ZZ(v[1])) for v in res]
    except Exception as err:
        raise SolverUnavailable(f"PARI qfbsolve failed: {err}") from None
    reps = [s for s in reps if s[0] ** 2 - D * s[1] ** 2 == N]
    if not reps:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "no solutions",
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
    return SolutionSet(cls.parsed.unknowns, sorted(reps), "orbits", desc,
                       complete=True, stream=stream)


def _solve_two_squares(cls):
    r"""
    All representations ``n = x^2 + y^2`` with ``0 <= x <= y``.

    For very large ``n`` (beyond ``MAX_TWO_SQUARES``), falls back to a single
    witness from Sage's ``two_squares``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 + y^2 = 610").solutions
        [(9, 23), (13, 21)]
        sage: solve("x^2 + y^2 = 21").kind
        'empty'
    """
    n = _zz(cls.data, "n")
    if n is None:
        raise SolverUnavailable("needs concrete n")
    if n < 0:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "n < 0",
                           complete=True)
    if n > MAX_TWO_SQUARES:
        try:
            x, y = two_squares(n)
        except ValueError:
            return SolutionSet(cls.parsed.unknowns, [], "empty",
                               "not a sum of two squares (a prime p ≡ 3 mod 4 "
                               "divides n to an odd power)", complete=True)
        return SolutionSet(cls.parsed.unknowns, [(x, y)], "witness",
                           "one representation (n too large for full "
                           "enumeration)", complete=False)
    sols = []
    x = 0
    while 2 * x * x <= n:
        y2 = n - x * x
        y = isqrt(y2)
        if y * y == y2:
            sols.append((ZZ(x), ZZ(y)))
        x += 1
    kind = "finite-complete" if sols else "empty"
    desc = ("all representations with 0 <= x <= y; the rest differ by signs "
            "and order" if sols else
            "not a sum of two squares (a prime p ≡ 3 mod 4 divides n to an "
            "odd power)")
    return SolutionSet(cls.parsed.unknowns, sols, kind, desc, complete=True)


def _solve_three_squares(cls):
    r"""
    All representations ``n = x^2 + y^2 + z^2`` with ``x <= y <= z``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 + y^2 + z^2 = 62").solutions
        [(1, 5, 6), (2, 3, 7)]
        sage: S = solve("x^2 + y^2 + z^2 = 7")    # 7 ≡ 7 mod 8
        sage: S.kind, S.complete
        ('empty', True)
    """
    n = _zz(cls.data, "n")
    if n is None:
        raise SolverUnavailable("needs concrete n")
    if n < 0:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "n < 0",
                           complete=True)
    m = n
    while m % 4 == 0:
        m //= 4
    if m % 8 == 7:
        return SolutionSet(cls.parsed.unknowns, [], "empty",
                           f"{n} = 4^a(8b+7): excluded by the Legendre-Gauss "
                           "criterion", complete=True)
    if n > MAX_THREE_SQUARES:
        sol = three_squares(n)
        return SolutionSet(cls.parsed.unknowns, [tuple(sol)], "witness",
                           "one representation (n too large for full "
                           "enumeration)", complete=False)
    sols = []
    x = 0
    while 3 * x * x <= n:
        y = x
        while x * x + 2 * y * y <= n:
            z2 = n - x * x - y * y
            z = isqrt(z2)
            if z * z == z2 and z >= y:
                sols.append((ZZ(x), ZZ(y), ZZ(z)))
            y += 1
        x += 1
    return SolutionSet(cls.parsed.unknowns, sols, "finite-complete",
                       "all representations with 0 <= x <= y <= z; the rest "
                       "differ by signs and order", complete=True)


def _solve_four_squares(cls):
    r"""
    All representations ``n = x^2 + y^2 + z^2 + w^2``, ``x <= y <= z <= w``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 + y^2 + z^2 + w^2 = 7").solutions
        [(1, 1, 1, 2)]
    """
    n = _zz(cls.data, "n")
    if n is None:
        raise SolverUnavailable("needs concrete n")
    if n < 0:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "n < 0",
                           complete=True)
    if n > MAX_FOUR_SQUARES:
        sol = four_squares(n)
        return SolutionSet(cls.parsed.unknowns, [tuple(sol)], "witness",
                           "one representation (Lagrange: always solvable; "
                           "n too large for full enumeration)",
                           complete=False)
    sols = []
    x = 0
    while 4 * x * x <= n:
        y = x
        while x * x + 3 * y * y <= n:
            z = y
            while x * x + y * y + 2 * z * z <= n:
                w2 = n - x * x - y * y - z * z
                w = isqrt(w2)
                if w * w == w2 and w >= z:
                    sols.append((ZZ(x), ZZ(y), ZZ(z), ZZ(w)))
                z += 1
            y += 1
        x += 1
    return SolutionSet(cls.parsed.unknowns, sols, "finite-complete",
                       "all representations with 0 <= x <= y <= z <= w; the "
                       "rest differ by signs and order (Lagrange: always "
                       "solvable)", complete=True)


def _solve_bqf(cls):
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
    a, b, c, n = (_zz(cls.data, k) for k in ("a", "b", "c", "n"))
    if None in (a, b, c, n):
        raise SolverUnavailable("needs concrete form and n")
    disc = b ** 2 - 4 * a * c
    if disc < 0 and a > 0 and abs(n) <= MAX_BQF:
        if n < 0:
            return SolutionSet(cls.parsed.unknowns, [], "empty",
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
        return SolutionSet(cls.parsed.unknowns, sols, kind,
                           "all representations (definite form)",
                           complete=True)
    form = BinaryQF([a, b, c])
    sol = form.solve_integer(n)
    if sol is None:
        return SolutionSet(cls.parsed.unknowns, [], "empty",
                           "no representation", complete=disc < 0)
    return SolutionSet(cls.parsed.unknowns, [tuple(sol)], "witness",
                       "one representation (BinaryQF.solve_integer); for "
                       "indefinite forms the full set is a union of "
                       "automorph orbits", complete=False)


def _solve_qf_zero(cls, gram=None):
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
        gram = sage_eval(str(cls.data["gram"]))
    G = matrix(QQ, gram)
    try:
        res = qfsolve(G)
    except Exception as err:
        raise SolverUnavailable(f"qfsolve failed: {err}") from None
    if res in ZZ:
        place = "the real place" if res == -1 else f"p = {res}"
        return SolutionSet(cls.parsed.unknowns, [], "empty",
                           f"no nontrivial solutions: local obstruction at "
                           f"{place}", complete=True)
    vec = [QQ(t) for t in res]
    den = lcm([t.denominator() for t in vec])
    ivec = [ZZ(t * den) for t in vec]
    g = gcd(ivec)
    ivec = tuple(t // g for t in ivec)
    if sum(1 for t in ivec if t < 0) > sum(1 for t in ivec if t > 0):
        ivec = tuple(-t for t in ivec)
    return SolutionSet(cls.parsed.unknowns, [ivec], "parametrized",
                       "one nontrivial solution; all others arise from it by "
                       "the standard conic/quadric parametrization",
                       complete=False)


def _solve_legendre(cls):
    r"""
    Legendre equation ``a x^2 + b y^2 + c z^2 = 0`` via ``qfsolve``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: x, y, z = solve("x^2 + 3*y^2 = 7*z^2").solutions[0]
        sage: x^2 + 3*y^2 == 7*z^2
        True
    """
    a, b, c = (_zz(cls.data, k) for k in ("a", "b", "c"))
    if None in (a, b, c):
        raise SolverUnavailable("needs concrete coefficients")
    return _solve_qf_zero(cls, gram=[[a, 0, 0], [0, b, 0], [0, 0, c]])


def _solve_pythagorean(cls):
    r"""
    Pythagorean triples: Euclid's parametrization as a stream.

    Iteration yields the primitive triples ``(m^2 - n^2, 2mn, m^2 + n^2)``
    for ``m > n >= 1`` coprime of opposite parity, ordered by ``m``; the full
    solution set consists of their multiples, sign changes, and leg swaps.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x^2 + y^2 = z^2")
        sage: S.first(4)
        [(3, 4, 5), (5, 12, 13), (15, 8, 17), (7, 24, 25)]
    """
    roles = _roles(cls)
    legs = roles.get("legs")
    hyp = roles.get("hypotenuse")

    def make(triple):
        a, b, c = triple
        if legs and hyp:
            return _ordered(cls, {legs[0]: a, legs[1]: b, hyp: c})
        return (a, b, c)

    def stream():
        for m in itertools.count(2):
            for n in range(1, m):
                if gcd(m, n) == 1 and (m - n) % 2 == 1:
                    yield make((m * m - n * n, 2 * m * n, m * m + n * n))

    first = list(itertools.islice(stream(), 4))
    return SolutionSet(
        cls.parsed.unknowns, first, "parametrized",
        "primitive solutions (m^2 - n^2, 2mn, m^2 + n^2), gcd(m, n) = 1, "
        "m ≢ n (mod 2); all solutions are multiples, sign changes and swaps "
        "of these; iteration enumerates the primitive triples",
        complete=True, stream=stream)


def _solve_weierstrass(cls):
    r"""
    Integral points on a Weierstrass model via ``E.integral_points``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: sorted(solve("y^2 = x^3 - 2").solutions)   # unknowns (y, x)
        [(-5, 3), (5, 3)]
    """
    if not cls.parsed.is_concrete:
        raise SolverUnavailable("needs concrete coefficients")
    ainvs = sage_eval(str(cls.data["ainvs"])) if "ainvs" in cls.data \
        else [0, 0, 0, 0, sage_eval(str(cls.data["k"]))]
    E = EllipticCurve(QQ, [QQ(t) for t in ainvs])
    pts = E.integral_points(both_signs=True)
    xname = cls.data.get("x", "x")
    yname = cls.data.get("y", "y")
    sols = [_ordered(cls, {xname: P[0], yname: P[1]}) for P in pts]
    return SolutionSet(
        cls.parsed.unknowns, sorted(sols), "finite-complete",
        f"all integral points on {E.ainvs()} (rank {E.rank()}); rational "
        "points are infinite iff the rank is positive", complete=True)


def _solve_thue(cls):
    r"""
    Thue equation via PARI's certified ``thue`` solver.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^3 + 2*y^3 = 11").solutions
        [(3, -2)]
        sage: solve("x^4 - 2*y^4 = 1").solutions
        [(-1, 0), (1, 0)]
    """
    pe = cls.parsed
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
    return SolutionSet(pe.unknowns, sorted(sols), "finite-complete",
                       "all solutions (PARI thue, certified)", complete=True)


def _solve_markov(cls):
    r"""
    Markov/Hurwitz equations: enumerate the Vieta tree.

    Iteration yields ascending-ordered positive tuples, ordered by largest
    entry; every solution is a permutation (with sign changes when the number
    of negative entries is even... for the classical positive tree, a
    permutation) of an enumerated tuple.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 + y^2 + z^2 = 3*x*y*z").first(5)
        [(1, 1, 1), (1, 1, 2), (1, 2, 5), (1, 5, 13), (2, 5, 29)]
        sage: solve("x^2 + y^2 + z^2 + w^2 = 4*x*y*z*w").first(2)
        [(1, 1, 1, 1), (1, 1, 1, 3)]
    """
    a = _zz(cls.data, "a")
    k = _zz(cls.data, "k")
    if a is None or k is None:
        raise SolverUnavailable("needs concrete Hurwitz parameters")
    from sage.misc.misc_c import prod as _prod
    seeds = set()
    if a == k:
        seeds.add((ZZ(1),) * k)
    if not seeds:
        for cand in itertools.combinations_with_replacement(range(1, 6), k):
            if sum(t * t for t in cand) == a * _prod(cand):
                seeds.add(tuple(ZZ(t) for t in cand))
    if not seeds:
        return SolutionSet(
            cls.parsed.unknowns, [], "witness",
            f"no fundamental solutions with entries <= 5 found for the "
            f"Hurwitz equation with a = {a}, k = {k} (Hurwitz classified "
            "the admissible a)", complete=False)

    def stream():
        heap = [(max(s), s) for s in seeds]
        heapq.heapify(heap)
        seen = set(seeds)
        while heap:
            _, s = heapq.heappop(heap)
            yield s
            others = _prod(s)
            for i in range(k):
                rest = others // s[i]
                child = tuple(sorted(s[:i] + (a * rest - s[i],) + s[i + 1:]))
                if child not in seen and all(t > 0 for t in child):
                    seen.add(child)
                    heapq.heappush(heap, (max(child), child))

    first = list(itertools.islice(stream(), 4))
    return SolutionSet(
        cls.parsed.unknowns, first, "infinite",
        "the Vieta/Markov tree: ascending tuples ordered by largest entry; "
        "all solutions are permutations of these", complete=True,
        stream=stream)


def _solve_egyptian(cls):
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
    a = _zz(cls.data, "a")
    n = _zz(cls.data, "n")
    k = _zz(cls.data, "k")
    if k is None and cls.slug == "erdos-straus":
        a, k = ZZ(4), ZZ(3)
    if a is None and cls.slug == "erdos-straus":
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
        cls.parsed.unknowns, sorted(sols), kind,
        "all solutions in positive integers with x_1 <= ... <= x_k; "
        "permutations give the rest", complete=True)


def _solve_catalan(cls):
    r"""
    Catalan's equation: Mihailescu's theorem.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^p - y^q = 1").solutions    # unknowns (x, p, y, q)
        [(3, 2, 2, 3)]
    """
    roles = _roles(cls)
    assignment = {roles["x"]: 3, roles["p"]: 2, roles["y"]: 2, roles["q"]: 3}
    return SolutionSet(
        cls.parsed.unknowns, [_ordered(cls, assignment)], "finite-complete",
        "Mihailescu's theorem: 3^2 - 2^3 = 1 is the only solution in "
        "integers > 1", complete=True)


def _solve_fermat(cls):
    r"""
    Fermat's equation: no nontrivial solutions (Wiles).

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: S = solve("x^4 + y^4 = z^4")
        sage: S.kind, S.complete
        ('empty', True)
    """
    return SolutionSet(
        cls.parsed.unknowns, [], "empty",
        "no solutions with xyz ≠ 0 for exponent ≥ 3 (Wiles); only the "
        "trivial solutions with a zero coordinate", complete=True)


def _solve_ramanujan_nagell(cls):
    r"""
    The classical Ramanujan-Nagell equation ``x^2 + 7 = 2^n``.

    EXAMPLES::

        sage: from diophantine_classifier import solve
        sage: solve("x^2 + 7 = 2^n").solutions
        [(1, 3), (3, 4), (5, 5), (11, 7), (181, 15)]
    """
    d, k, base = _zz(cls.data, "d"), _zz(cls.data, "k"), _zz(cls.data, "base")
    if (d, k, base) != (7, 1, 2):
        raise SolverUnavailable(
            "only the classical x^2 + 7 = 2^n is hardwired; general (d, k, b) "
            "need a Baker + LLL computation (Petho-de Weger)")
    roles = _roles(cls)
    pairs = [(1, 3), (3, 4), (5, 5), (11, 7), (181, 15)]
    sols = [_ordered(cls, {roles["x"]: x, roles["n"]: n}) for x, n in pairs]
    return SolutionSet(
        cls.parsed.unknowns, sols, "finite-complete",
        "Nagell's theorem: n ∈ {3, 4, 5, 7, 15} (x > 0 shown; -x symmetric)",
        complete=True)


#: dispatch table: family slug -> solver function
SOLVERS = {
    "univariate": _solve_univariate,
    "linear": _solve_linear,
    "pell": _solve_pell,
    "pell-like": _solve_pell_like,
    "sum-of-two-squares": _solve_two_squares,
    "sum-of-three-squares": _solve_three_squares,
    "sum-of-four-squares": _solve_four_squares,
    "binary-qf-representation": _solve_bqf,
    "quadratic-form-zero": _solve_qf_zero,
    "legendre": _solve_legendre,
    "pythagorean": _solve_pythagorean,
    "mordell": _solve_weierstrass,
    "elliptic-weierstrass": _solve_weierstrass,
    "thue": _solve_thue,
    "markov-hurwitz": _solve_markov,
    "egyptian-fractions": _solve_egyptian,
    "erdos-straus": _solve_egyptian,
    "catalan": _solve_catalan,
    "fermat": _solve_fermat,
    "ramanujan-nagell": _solve_ramanujan_nagell,
}


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

    The dispatch tries the primary family first, then walks up the lineage
    (so a genus-two equation with no dedicated solver falls through to the
    hyperelliptic one, etc.).  When nothing is wired, the raised
    :class:`SolverUnavailable` carries software pointers and filled code
    templates from the registry.

    EXAMPLES::

        sage: from diophantine_classifier import classify, solve, SolverUnavailable
        sage: solve("x^2 + y^2 = 610").solutions
        [(9, 23), (13, 21)]
        sage: solve(classify("x^2 - 61*y^2 = 1")).kind
        'infinite'
        sage: try:
        ....:     solve("x^2 + 5 = y^n")
        ....: except SolverUnavailable as err:
        ....:     print("lebesgue-nagell" in str(err))
        True
    """
    cls = equation if isinstance(equation, Classification) \
        else classify(equation, params=params, domain=domain)
    if cls.is_composite:
        raise SolverUnavailable(
            "reducible equation: solve each component separately "
            f"({', '.join(c.slug for c in cls.components)})")
    for slug in [cls.slug] + cls.lineage:
        fn = SOLVERS.get(slug)
        if fn is not None:
            return fn(cls)
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
