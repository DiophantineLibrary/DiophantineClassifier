"""Solvers: dispatch a classified equation to Sage/PARI machinery.

``solve`` handles the families where standard software gives a complete or
canonical answer (linear, Pell, quadratic forms, conics, Weierstrass integral
points, Thue, ...).  For families without a wired-up solver it raises
:class:`SolverUnavailable` carrying the registry's software pointers and code
templates, matching the Library design: every equation page should offer
runnable code even when we cannot run it ourselves.
"""

from dataclasses import dataclass, field

from sage.all import (QQ, ZZ, BinaryQF, EllipticCurve, QuadraticField,
                      continued_fraction, gcd, lcm, matrix, pari, sage_eval,
                      xgcd, two_squares, three_squares, four_squares)
from sage.quadratic_forms.qfsolve import qfsolve

from .classify import Classification, classify


class SolverUnavailable(NotImplementedError):
    pass


@dataclass
class SolutionSet:
    variables: tuple
    solutions: list
    kind: str                  # finite-complete | witness | orbits | parametrized | infinite | empty
    description: str = ""
    complete: bool = False

    def __repr__(self):
        head = f"SolutionSet({self.kind}, vars={list(self.variables)})"
        body = ""
        if self.solutions:
            shown = ", ".join(str(s) for s in self.solutions[:12])
            more = "" if len(self.solutions) <= 12 else \
                f", ... ({len(self.solutions)} total)"
            body = f"\n  solutions: {shown}{more}"
        if self.description:
            body += f"\n  {self.description}"
        return head + body


def _zz(data, key):
    try:
        return ZZ(sage_eval(str(data[key])))
    except (KeyError, TypeError, ValueError, SyntaxError):
        return None


def _ordered(cls, assignment):
    """Tuple of values ordered like the equation's unknowns."""
    return tuple(assignment[v] for v in cls.parsed.unknowns)


def _roles(cls):
    return cls.data.get("roles", {})


# ---------------------------------------------------------------- solvers

def _solve_univariate(cls):
    pu = cls.parsed.poly.univariate_polynomial()
    ring = QQ if cls.parsed.domain == "QQ" else ZZ
    roots = pu.roots(ring, multiplicities=False)
    if cls.parsed.domain == "NN":
        roots = [r for r in roots if r >= 0]
    return SolutionSet(cls.parsed.unknowns, [(r,) for r in sorted(roots)],
                       "finite-complete", "all roots in the domain",
                       complete=True)


def _solve_linear(cls):
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
    particular = [c * scale for c in combo]
    kernel = matrix(ZZ, [coeffs]).right_kernel_matrix().rows()
    desc = (f"general solution: {tuple(particular)} + integer combinations "
            f"of {[tuple(r) for r in kernel]}")
    return SolutionSet(cls.parsed.unknowns, [tuple(particular)],
                       "parametrized", desc, complete=True)


def _solve_pell(cls):
    D = _zz(cls.data, "D")
    N = _zz(cls.data, "N")
    if D is None or N not in (1, -1):
        raise SolverUnavailable("Pell solver needs concrete D and N = ±1")
    cf = continued_fraction(QuadraticField(D).gen())
    ell = len(cf.period())
    conv = cf.convergent(ell - 1)
    x1, y1 = conv.numerator(), conv.denominator()
    norm = x1 ** 2 - D * y1 ** 2
    if N == 1:
        if norm == 1:
            fund = (x1, y1)
        else:
            fund = (x1 ** 2 + D * y1 ** 2, 2 * x1 * y1)
        desc = (f"infinitely many; all ±(x + y√{D}) = ±(fundamental)^k with "
                f"fundamental solution {fund}")
        return SolutionSet(cls.parsed.unknowns, [fund], "infinite", desc,
                           complete=True)
    if norm == -1:
        desc = (f"infinitely many; fundamental solution {(x1, y1)}, "
                "odd powers of the fundamental unit")
        return SolutionSet(cls.parsed.unknowns, [(x1, y1)], "infinite", desc,
                           complete=True)
    return SolutionSet(cls.parsed.unknowns, [], "empty",
                       f"x^2 - {D}y^2 = -1 has no solutions (continued "
                       "fraction period is even)", complete=True)


def _solve_pell_like(cls):
    D = _zz(cls.data, "D")
    N = _zz(cls.data, "N")
    if D is None or N is None:
        raise SolverUnavailable("generalized Pell solver needs concrete D, N")
    if N in (1, -1):
        return _solve_pell(cls)
    try:
        res = pari(f"qfbsolve(Qfb(1,0,{-D}),{N},1)")
        sols = [(ZZ(v[0]), ZZ(v[1])) for v in res]
    except Exception as err:
        raise SolverUnavailable(f"PARI qfbsolve failed: {err}") from None
    sols = [s for s in sols if s[0] ** 2 - D * s[1] ** 2 == N]
    kind = "orbits" if sols else "empty"
    desc = ("orbit representatives; the full solution set is their orbit "
            f"under the automorph group of x^2 - {D}y^2 (powers of the "
            "fundamental Pell solution)") if sols else "no solutions"
    return SolutionSet(cls.parsed.unknowns, sols, kind, desc,
                       complete=not sols)


def _solve_two_squares(cls):
    n = _zz(cls.data, "n")
    if n is None:
        raise SolverUnavailable("needs concrete n")
    if n < 0:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "n < 0",
                           complete=True)
    try:
        x, y = two_squares(n)
    except ValueError:
        return SolutionSet(
            cls.parsed.unknowns, [], "empty",
            f"{n} is not a sum of two squares (a prime p ≡ 3 mod 4 divides "
            "it to an odd power)", complete=True)
    return SolutionSet(cls.parsed.unknowns, [(x, y)], "witness",
                       "one representation; all others via signs, swaps and "
                       "Gaussian-integer factorization", complete=False)


def _solve_three_squares(cls):
    n = _zz(cls.data, "n")
    if n is None:
        raise SolverUnavailable("needs concrete n")
    if n < 0:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "n < 0",
                           complete=True)
    try:
        sol = three_squares(n)
    except ValueError:
        return SolutionSet(cls.parsed.unknowns, [], "empty",
                           f"{n} = 4^a(8b+7): excluded by the "
                           "Legendre-Gauss criterion", complete=True)
    return SolutionSet(cls.parsed.unknowns, [tuple(sol)], "witness",
                       "one representation", complete=False)


def _solve_four_squares(cls):
    n = _zz(cls.data, "n")
    if n is None:
        raise SolverUnavailable("needs concrete n")
    if n < 0:
        return SolutionSet(cls.parsed.unknowns, [], "empty", "n < 0",
                           complete=True)
    sol = four_squares(n)
    return SolutionSet(cls.parsed.unknowns, [tuple(sol)], "witness",
                       "one representation (Lagrange: always solvable)",
                       complete=False)


def _solve_bqf(cls):
    a, b, c, n = (_zz(cls.data, k) for k in ("a", "b", "c", "n"))
    if None in (a, b, c, n):
        raise SolverUnavailable("needs concrete form and n")
    form = BinaryQF([a, b, c])
    sol = form.solve_integer(n)
    if sol is None:
        definite = b ** 2 - 4 * a * c < 0
        return SolutionSet(cls.parsed.unknowns, [], "empty",
                           "no representation", complete=definite)
    return SolutionSet(cls.parsed.unknowns, [tuple(sol)], "witness",
                       "one representation (BinaryQF.solve_integer)",
                       complete=False)


def _solve_qf_zero(cls, gram=None):
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
    return SolutionSet(cls.parsed.unknowns, [ivec], "parametrized",
                       "one nontrivial solution; all others from it by the "
                       "standard conic/quadric parametrization",
                       complete=False)


def _solve_legendre(cls):
    a, b, c = (_zz(cls.data, k) for k in ("a", "b", "c"))
    if None in (a, b, c):
        raise SolverUnavailable("needs concrete coefficients")
    return _solve_qf_zero(cls, gram=[[a, 0, 0], [0, b, 0], [0, 0, c]])


def _solve_pythagorean(cls):
    sols = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25)]
    return SolutionSet(
        cls.parsed.unknowns, sols, "parametrized",
        "all primitive solutions: (m^2 - n^2, 2mn, m^2 + n^2) up to order and "
        "signs, with gcd(m, n) = 1, m ≢ n (mod 2); first few shown",
        complete=True)


def _solve_weierstrass(cls):
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
        f"all integral points on {E.ainvs()} (rank {E.rank()}); "
        "rational points are infinite iff the rank is positive",
        complete=True)


def _solve_thue(cls):
    pe = cls.parsed
    P = pe.poly
    R = pe.poly_ring
    x, y = R.gens()
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


def _solve_catalan(cls):
    roles = _roles(cls)
    assignment = {roles["x"]: 3, roles["p"]: 2, roles["y"]: 2, roles["q"]: 3}
    return SolutionSet(
        cls.parsed.unknowns, [_ordered(cls, assignment)], "finite-complete",
        "Mihailescu's theorem: 3^2 - 2^3 = 1 is the only solution in "
        "integers > 1", complete=True)


def _solve_fermat(cls):
    return SolutionSet(
        cls.parsed.unknowns, [], "empty",
        "no solutions with xyz ≠ 0 for exponent ≥ 3 (Wiles); only the "
        "trivial solutions with a zero coordinate", complete=True)


def _solve_ramanujan_nagell(cls):
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
    "catalan": _solve_catalan,
    "fermat": _solve_fermat,
    "ramanujan-nagell": _solve_ramanujan_nagell,
}


def solve(equation, params=(), domain="ZZ"):
    """Solve (or partially solve) an equation via its classification.

    Accepts an equation string or an existing :class:`Classification`.
    Raises :class:`SolverUnavailable` — with software pointers — for families
    without a wired-up solver.
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
