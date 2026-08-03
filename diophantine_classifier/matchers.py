r"""
Structural matchers: recognize which families a parsed equation belongs to.

Each matcher inspects a :class:`~diophantine_classifier.parsing.ParsedEquation`
and emits :class:`Match` objects (family slug + extracted data).  :func:`run`
collects matches from every applicable matcher; ranking by DAG specificity
happens in :mod:`diophantine_classifier.classify`.

Wave 1 recognizes equations presented in (or very near) a family's standard
form, plus genus-based routing for plane curves.  Deeper normalization
(unimodular reduction, completing the square, Nagell's algorithm) is the
wave-2 transformation layer; see ``docs/DESIGN.md``.

EXAMPLES::

    sage: from diophantine_classifier.parsing import parse
    sage: from diophantine_classifier.matchers import run
    sage: [m.slug for m in run(parse("x^2 - 61*y^2 = 1"))][-2:]
    ['pell-like', 'pell']
"""

from dataclasses import dataclass, field

from sage.all import QQ, ZZ, Curve, EllipticCurve, prod

#: skip genus computations for inputs of absurdly large degree
MAX_GENUS_DEGREE = 20


@dataclass
class Match:
    r"""
    One structural match: a family together with extracted data.

    ATTRIBUTES:

    - ``slug`` -- string; the matched family's registry slug.
    - ``data`` -- dict of family-specific extracted data, with stringified
      values (e.g. ``{"D": "61", "N": "1"}`` for a Pell match, a-invariants
      for a Weierstrass match).  A ``"roles"`` entry, when present, maps
      structural roles to the user's variable names (used by solvers to
      order solution tuples).
    - ``summary`` -- string; one-line human-readable description with the
      parameters filled in.
    - ``transform`` -- string; description of any normalization applied
      before the pattern matched (sign flip, variable swap, ...); empty when
      the equation was already in standard form.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: m = classify("x^2 - 61*y^2 = 1").primary
        sage: m.slug, m.data["D"]
        ('pell', '61')
    """
    slug: str
    data: dict = field(default_factory=dict)
    summary: str = ""
    transform: str = ""

    def __repr__(self):
        r"""
        Terse representation.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("x^2 - 61*y^2 = 1").primary
            Match('pell')
        """
        return f"Match({self.slug!r})"


def _sign_of(c):
    r"""
    Return -1/0/+1 for constant coefficients, ``None`` if parametric.

    INPUT:

    - ``c`` -- a coefficient: an element of ``QQ`` or of ``QQ[params]``

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _sign_of
        sage: from sage.all import QQ, PolynomialRing
        sage: _sign_of(QQ(-3)), _sign_of(QQ(0)), _sign_of(QQ(2))
        (-1, 0, 1)
        sage: R = PolynomialRing(QQ, "k")
        sage: _sign_of(R.gen()) is None
        True
    """
    try:
        q = QQ(c)
    except (TypeError, ValueError):
        return None
    return 0 if q == 0 else (1 if q > 0 else -1)


def _as_int(c):
    r"""
    Coerce a coefficient to a Sage integer, or return ``None``.

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _as_int
        sage: from sage.all import QQ
        sage: _as_int(QQ(7))
        7
        sage: _as_int(QQ(1)/2) is None
        True
    """
    try:
        return ZZ(c)
    except (TypeError, ValueError):
        return None


def _homogeneous_parts(P):
    r"""
    Decompose a polynomial into its homogeneous components.

    INPUT:

    - ``P`` -- multivariate polynomial

    OUTPUT: dict mapping total degree to the homogeneous component

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _homogeneous_parts
        sage: R.<x, y> = QQ[]
        sage: parts = _homogeneous_parts(y^2 - x^3 + 2)
        sage: sorted(parts)
        [0, 2, 3]
        sage: parts[3]
        -x^3
    """
    parts = {}
    for c, m in zip(P.coefficients(), P.monomials()):
        d = m.degree()
        parts[d] = parts.get(d, P.parent().zero()) + c * m
    return parts


def _diagonal_scan(P, R):
    r"""
    Detect a diagonal shape: every nonconstant monomial a pure power.

    INPUT:

    - ``P`` -- multivariate polynomial
    - ``R`` -- its parent ring

    OUTPUT: ``(entries, const)`` with ``entries`` a list of triples
    ``(coefficient, variable name, exponent)``, one per variable used, and
    ``const`` the constant term — or ``None`` if some monomial mixes
    variables or a variable occurs with two different exponents

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _diagonal_scan
        sage: R.<x, y, z> = QQ[]
        sage: entries, const = _diagonal_scan(x^3 + y^3 + z^3 - 42, R)
        sage: [(str(v), e) for _, v, e in entries]
        [('x', 3), ('y', 3), ('z', 3)]
        sage: const
        -42
        sage: _diagonal_scan(x*y + z, R) is None
        True
    """
    entries = []
    seen = set()
    const = P.parent().base_ring().zero()
    gens = R.gens()
    for c, m in zip(P.coefficients(), P.monomials()):
        degs = m.degrees()
        nonzero = [(i, e) for i, e in enumerate(degs) if e]
        if not nonzero:
            const += c
            continue
        if len(nonzero) != 1:
            return None
        i, e = nonzero[0]
        if i in seen:
            return None
        seen.add(i)
        entries.append((c, str(gens[i]), ZZ(e)))
    return entries, const


def _gram(P, R):
    r"""
    Gram matrix of a quadratic form, as a list of lists of rationals.

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _gram
        sage: R.<x, y> = QQ[]
        sage: _gram(x^2 + 3*x*y - y^2, R)
        [[1, 3/2], [3/2, -1]]
    """
    gens = R.gens()
    k = len(gens)
    G = [[QQ(0)] * k for _ in range(k)]
    for i in range(k):
        G[i][i] = QQ(P.monomial_coefficient(gens[i] ** 2))
        for j in range(i + 1, k):
            half = QQ(P.monomial_coefficient(gens[i] * gens[j])) / 2
            G[i][j] = G[j][i] = half
    return G


def _plane_curve_genus(P):
    r"""
    Geometric genus of the plane curve ``P = 0``, or ``None``.

    Returns ``None`` for degrees above :data:`MAX_GENUS_DEGREE` or when the
    curve construction fails (reducible input, etc.).

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _plane_curve_genus
        sage: R.<x, y> = QQ[]
        sage: _plane_curve_genus(y^2 - x^5 + x - 1)
        2
        sage: _plane_curve_genus(y^2 - x^3)     # cuspidal cubic
        0
    """
    if P.total_degree() > MAX_GENUS_DEGREE:
        return None
    try:
        return ZZ(Curve(P).genus())
    except Exception:
        return None


def _is_irreducible(P):
    r"""
    Whether ``P`` is irreducible (up to constants), or ``None`` on failure.

    EXAMPLES::

        sage: from diophantine_classifier.matchers import _is_irreducible
        sage: R.<x, y> = QQ[]
        sage: _is_irreducible(x^2 + y^2 - 1)
        True
        sage: _is_irreducible(x^2 - y^2)
        False
    """
    try:
        fac = P.factor()
    except Exception:
        return None
    nontrivial = [(g, e) for g, e in fac if g.degree() > 0]
    return len(nontrivial) == 1 and nontrivial[0][1] == 1


# --------------------------------------------------------------------------
# polynomial equations
# --------------------------------------------------------------------------

def _match_polynomial(pe):
    r"""
    Matches for a purely polynomial equation.

    Always emits ``general-polynomial`` with basic invariants, then
    dispatches on the number of unknowns.

    INPUT:

    - ``pe`` -- a polynomial :class:`~diophantine_classifier.parsing.ParsedEquation`

    OUTPUT: list of :class:`Match`

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _match_polynomial
        sage: [m.slug for m in _match_polynomial(parse("3*x + 5*y = 1"))]
        ['general-polynomial', 'linear']
        sage: [m.slug for m in _match_polynomial(parse("x^2 - 5*x + 6 = 0"))]
        ['general-polynomial', 'univariate']
    """
    P, R = pe.poly, pe.poly_ring
    k = R.ngens()
    d = P.total_degree()
    out = [Match(
        "general-polynomial",
        data={"variables": k, "degree": ZZ(d),
              "homogeneous": P.is_homogeneous()},
        summary=f"polynomial equation in {k} unknowns of degree {d}",
    )]

    if k == 1:
        out.append(Match("univariate", data={"f": str(P)},
                         summary="single-variable polynomial: solve by factoring"))
        return out

    if d == 1:
        coeffs = [P.monomial_coefficient(g) for g in R.gens()]
        b = -P.constant_coefficient()
        out.append(Match("linear",
                         data={"coeffs": [str(c) for c in coeffs], "b": str(b)},
                         summary="linear Diophantine equation"))
        return out

    if k == 2:
        out.extend(_match_binary(pe, P, R, d))
    else:
        out.extend(_match_multivar(pe, P, R, k, d))
    return out


def _match_binary(pe, P, R, d):
    r"""
    Matches for polynomial equations in two unknowns.

    Degree 2 goes through the binary-quadratic battery (Pell, sums of two
    squares, representation by a form); degree ≥ 3 through binary forms
    (Thue), curve shapes (Weierstrass, quartic, hyperelliptic,
    superelliptic), then genus routing.

    INPUT:

    - ``pe`` -- the parsed equation; ``P`` its polynomial in ring ``R``;
      ``d`` the total degree

    OUTPUT: list of :class:`Match`

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _match_binary
        sage: pe = parse("x^2 - 61*y^2 = 5")
        sage: [m.slug for m in _match_binary(pe, pe.poly, pe.poly_ring, 2)]
        ['binary-quadratic', 'binary-qf-representation', 'pell-like']
        sage: pe = parse("x^3 + 2*y^3 = 11")
        sage: [m.slug for m in _match_binary(pe, pe.poly, pe.poly_ring, 3)]
        ['thue']
    """
    x, y = R.gens()
    out = []
    if d == 2:
        a = P.monomial_coefficient(x ** 2)
        b = P.monomial_coefficient(x * y)
        c = P.monomial_coefficient(y ** 2)
        dd = P.monomial_coefficient(x)
        ee = P.monomial_coefficient(y)
        f0 = P.constant_coefficient()
        disc = b ** 2 - 4 * a * c
        out.append(Match(
            "binary-quadratic",
            data={"a": str(a), "b": str(b), "c": str(c), "d": str(dd),
                  "e": str(ee), "f": str(f0), "disc": str(disc)},
            summary=f"binary quadratic with discriminant {disc}",
        ))
        if dd == 0 and ee == 0 and f0 != 0:
            n = -f0
            transform = ""
            # orient: prefer positive coefficient on the first square
            if _sign_of(a) == -1:
                a, b, c, n = -a, -b, -c, -n
                transform = "multiplied by -1"
            if b == 0 and _sign_of(a) == 1 and _sign_of(c) == -1 \
                    and a != 1 and c == -1:
                # D x^2 - y^2 = N: swap roles so it reads x^2 - D y^2
                a, c = -c, -a
                transform = (transform + "; " if transform else "") \
                    + "swapped variables"
            out.append(Match(
                "binary-qf-representation",
                data={"a": str(a), "b": str(b), "c": str(c), "n": str(n),
                      "disc": str(disc)},
                summary=f"representation of {n} by the form ({a}, {b}, {c})",
                transform=transform,
            ))
            if b == 0 and a == 1:
                if c == 1:
                    out.append(Match("sum-of-two-squares", data={"n": str(n)},
                                     summary=f"x^2 + y^2 = {n}",
                                     transform=transform))
                else:
                    D = -c
                    Dz = _as_int(D)
                    if Dz is not None and Dz > 0:
                        if Dz.is_square():
                            out.append(Match(
                                "binary-form-reducible",
                                data={"factors":
                                      f"(x - {Dz.sqrt()}*y)*(x + {Dz.sqrt()}*y)",
                                      "m": str(n)},
                                summary=f"x^2 - {Dz}y^2 factors (D is a "
                                        "square): solve by divisor "
                                        "enumeration",
                                transform=transform,
                            ))
                        else:
                            Nz = _as_int(n)
                            out.append(Match(
                                "pell-like", data={"D": str(Dz), "N": str(n)},
                                summary=f"x^2 - {Dz}y^2 = {n}",
                                transform=transform,
                            ))
                            if Nz is not None and Nz in (1, -1):
                                out.append(Match(
                                    "pell", data={"D": str(Dz), "N": str(Nz)},
                                    summary=f"Pell equation with D = {Dz}"
                                            + (" (negative Pell)"
                                               if Nz == -1 else ""),
                                    transform=transform,
                                ))
                    elif Dz is None:
                        # parametric D
                        out.append(Match("pell-like",
                                         data={"D": str(D), "N": str(n)},
                                         summary=f"x^2 - ({D})y^2 = {n}",
                                         transform=transform))
                        if n == 1:
                            out.append(Match("pell",
                                             data={"D": str(D), "N": "1"},
                                             summary=f"Pell equation with "
                                                     f"D = {D}",
                                             transform=transform))
        elif dd == 0 and ee == 0 and f0 == 0:
            out.append(Match(
                "binary-form-reducible",
                data={"m": "0"},
                summary="homogeneous quadratic = 0: rational lines exist iff "
                        "the discriminant is a square",
            ))
        return out

    # degree >= 3 in two variables
    parts = _homogeneous_parts(P)
    c0 = P.constant_coefficient()
    nonconst = sorted(deg for deg in parts if deg > 0)
    if len(nonconst) == 1 and nonconst[0] == d:
        # F_d(x, y) = m
        F = parts[d]
        m = -c0
        data = {"form": str(F), "m": str(m), "degree": ZZ(d),
                "fx": str(F.subs({y: 1})), "x": str(x), "y": str(y)}
        if not pe.is_concrete:
            out.append(Match("binary-form", data=data,
                             summary=f"binary form of degree {d} = {m}"))
            return out
        fac = F.factor()
        nontrivial = [(g, e) for g, e in fac if g.degree() > 0]
        irreducible = len(nontrivial) == 1 and nontrivial[0][1] == 1
        if m == 0:
            out.append(Match(
                "binary-form-reducible", data=data,
                summary="homogeneous form = 0: solutions come from rational "
                        "linear factors" + ("" if not irreducible else
                                            " (none here: only (0,0))"),
            ))
        elif irreducible:
            mz = _as_int(m)
            out.append(Match(
                "thue", data=dict(data, m=str(mz if mz is not None else m)),
                summary=f"Thue equation of degree {d}",
            ))
        else:
            data["factors"] = str(fac)
            out.append(Match(
                "binary-form-reducible", data=data,
                summary="reducible binary form: enumerate factorizations of "
                        f"{m} across the factors",
            ))
        return out

    shape = _match_two_var_shapes(pe, P, R, d)
    if shape:
        out.extend(shape)
        return out

    out.extend(_genus_route(pe, P, affine=True))
    return out


def _match_two_var_shapes(pe, P, R, d):
    r"""
    Weierstrass / quartic / hyperelliptic / superelliptic shapes.

    Two passes over both variable orderings: the quadratic-in-``v`` shapes
    (Weierstrass, quartic, hyperelliptic) are tried for *both* orientations
    before any superelliptic shape, so that e.g. ``x^2 = y^3 - k`` is
    recognized as a Mordell equation rather than as ``y^3 = x^2 + k``.

    OUTPUT: list of :class:`Match`, or ``None`` when no shape applies (the
    caller then falls through to genus routing)

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _match_two_var_shapes
        sage: pe = parse("y^2 = x^3 - 2")
        sage: [m.slug for m in
        ....:  _match_two_var_shapes(pe, pe.poly, pe.poly_ring, 3)]
        ['elliptic-weierstrass', 'mordell']
        sage: pe = parse("y^3 = x^4 + 2")
        sage: [m.slug for m in
        ....:  _match_two_var_shapes(pe, pe.poly, pe.poly_ring, 4)]
        ['superelliptic']
        sage: pe = parse("y^2 = x^3")     # cusp: shapes decline, genus routes
        sage: _match_two_var_shapes(pe, pe.poly, pe.poly_ring, 3) is None
        True
    """
    orientations = (tuple(R.gens()), tuple(reversed(R.gens())))
    for v, u in orientations:
        if P.degree(v) != 2:
            continue
        coeffs = [P.coefficient({v: j}) for j in range(3)]
        A = coeffs[2]
        if not A.is_constant():
            continue
        B1, C = coeffs[1], coeffs[0]
        if A == -1:
            A, B1, C = 1, -B1, -C
        if A != 1:
            continue
        # v^2 + B1*v + C = 0  <=>  v^2 + B1*v = g(u), g = -C
        g = -C
        n = g.degree(u)
        if B1.degree(u) <= 1 and n == 3 \
                and g.monomial_coefficient(u ** 3) == 1:
            a1 = B1.monomial_coefficient(u)
            a3 = B1.constant_coefficient()
            a2 = g.monomial_coefficient(u ** 2)
            a4 = g.monomial_coefficient(u)
            a6 = g.constant_coefficient()
            ainvs = [a1, a2, a3, a4, a6]
            if pe.is_concrete:
                try:
                    EllipticCurve(QQ, [QQ(t) for t in ainvs])
                except (ArithmeticError, TypeError, ValueError):
                    continue  # singular: fall through to genus routing
            data = {"ainvs": "[%s]" % ", ".join(str(t) for t in ainvs),
                    "magma_ainvs": "[%s]" % ", ".join(str(t) for t in ainvs),
                    "x": str(u), "y": str(v)}
            out = [Match("elliptic-weierstrass", data=data,
                         summary=f"Weierstrass equation with "
                                 f"a-invariants {data['ainvs']}")]
            if a1 == 0 and a2 == 0 and a3 == 0 and a4 == 0:
                out.append(Match("mordell", data={"k": str(a6),
                                                  "x": str(u), "y": str(v)},
                                 summary=f"Mordell equation with k = {a6}"))
            return out
        if B1 == 0 and n == 4:
            qu = g.univariate_polynomial()
            if pe.is_concrete and qu.discriminant() == 0:
                return None
            qcoeffs = [qu[i] for i in range(5)]
            return [Match(
                "elliptic-quartic",
                data={"q": str(g), "x": str(u), "y": str(v),
                      "qcoeffs": str(list(reversed(qcoeffs)))},
                summary=f"genus-one quartic {v}^2 = {g}",
            )]
        if B1 == 0 and n >= 5:
            qu = g.univariate_polynomial()
            if pe.is_concrete and not qu.is_squarefree():
                return None
            genus = (n - 1) // 2 if n % 2 else (n - 2) // 2
            slug = "genus-two" if genus == 2 else "hyperelliptic"
            return [Match(
                slug,
                data={"f": str(g), "genus": ZZ(genus), "x": str(u),
                      "y": str(v),
                      "magma_coeffs": str([qu[i] for i in range(n + 1)])},
                summary=f"hyperelliptic: {v}^2 = degree-{n} polynomial, "
                        f"genus {genus}",
            )]
    # second pass: superelliptic shapes v^m = f(u), m >= 3
    for v, u in orientations:
        m = P.degree(v)
        if m < 3:
            continue
        coeffs = [P.coefficient({v: j}) for j in range(m + 1)]
        A = coeffs[m]
        if not A.is_constant():
            continue
        if any(c != 0 for c in coeffs[1:m]):
            continue
        C = coeffs[0]
        if A == -1:
            A, C = 1, -C
        if A != 1:
            continue
        g = -C
        n = g.degree(u)
        if n >= 2:
            qu = g.univariate_polynomial()
            if pe.is_concrete and not qu.is_squarefree():
                continue  # e.g. x^3 = y^2: let genus routing handle it
            return [Match(
                "superelliptic",
                data={"m": ZZ(m), "f": str(g), "x": str(u), "y": str(v)},
                summary=f"superelliptic: {v}^{m} = {g}",
            )]
    return None


def _genus_route(pe, P, affine=True):
    r"""
    Route an irreducible plane curve by genus.

    Mirrors the Library's classification plan: genus 0 → parametrize,
    genus 1 → find a point then reduce to Weierstrass form, genus ≥ 2 →
    Faltings finiteness and Chabauty-type methods.

    OUTPUT: list of at most one :class:`Match` (empty for parametric input,
    reducible polynomials, or when the genus computation is skipped)

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _genus_route
        sage: pe = parse("x^2*y^2 = x^3 + 1")
        sage: [m.slug for m in _genus_route(pe, pe.poly)]
        ['genus-one-curve']
        sage: pe = parse("y^2 = x^3")
        sage: [m.slug for m in _genus_route(pe, pe.poly)]
        ['genus-zero-curve']
    """
    if not pe.is_concrete:
        return []
    if _is_irreducible(P) is not True:
        return []
    g = _plane_curve_genus(P)
    if g is None:
        return []
    data = {"genus": g, "model": "affine" if affine else "projective"}
    if g == 0:
        return [Match("genus-zero-curve", data=data,
                      summary="plane curve of genus 0: parametrize")]
    if g == 1:
        return [Match("genus-one-curve", data=data,
                      summary="plane curve of genus 1: find a point, then "
                              "reduce to Weierstrass form")]
    return [Match("general-curve", data=data,
                  summary=f"plane curve of genus {g}: Faltings finiteness; "
                          "Chabauty-type methods apply")]


def _match_multivar(pe, P, R, k, d):
    r"""
    Matches for polynomial equations in three or more unknowns.

    Handles quadratic forms (isotropy, representation, affine quadrics),
    Markov–Hurwitz shapes, diagonal equations (generalized Fermat, sums of
    cubes, Waring, equal sums of like powers), and plane projective curves.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _match_multivar
        sage: pe = parse("x^2 + y^2 = z^2")
        sage: [m.slug for m in
        ....:  _match_multivar(pe, pe.poly, pe.poly_ring, 3, 2)]
        ['quadratic-form-zero', 'legendre', 'pythagorean']
        sage: pe = parse("x^3 + y^3 + z^3 = 42")
        sage: [m.slug for m in
        ....:  _match_multivar(pe, pe.poly, pe.poly_ring, 3, 3)]
        ['sum-of-three-cubes']
    """
    out = []
    gens = R.gens()
    parts = _homogeneous_parts(P)
    c0 = P.constant_coefficient()

    if d == 2:
        has_linear = 1 in parts
        if not has_linear:
            Q2 = parts.get(2, P.parent().zero())
            gram = _gram(Q2, R)
            if c0 == 0:
                out.append(Match("quadratic-form-zero",
                                 data={"gram": str(gram), "k": k},
                                 summary=f"isotropy of a quadratic form in "
                                         f"{k} variables (Hasse-Minkowski)"))
                scan = _diagonal_scan(Q2, R)
                if scan and k == 3:
                    entries, _ = scan
                    if len(entries) == 3:
                        coeffs = [c for c, _, _ in entries]
                        signs = [_sign_of(c) for c in coeffs]
                        if None not in signs:
                            if signs.count(1) == 1:
                                coeffs = [-c for c in coeffs]
                                signs = [-s for s in signs]
                            a, b, c = coeffs
                            out.append(Match(
                                "legendre",
                                data={"a": str(a), "b": str(b), "c": str(c)},
                                summary=f"Legendre equation "
                                        f"({a})x^2 + ({b})y^2 + ({c})z^2 = 0",
                            ))
                            if sorted(coeffs) == [-1, 1, 1]:
                                legs = [entries[i][1] for i in range(3)
                                        if coeffs[i] == 1]
                                hyp = [entries[i][1] for i in range(3)
                                       if coeffs[i] == -1][0]
                                out.append(Match(
                                    "pythagorean",
                                    data={"roles": {"legs": legs,
                                                    "hypotenuse": hyp}},
                                    summary="Pythagorean equation: solutions "
                                            "(m^2-n^2, 2mn, m^2+n^2)",
                                ))
            else:
                n = -c0
                out.append(Match("quadratic-form-representation",
                                 data={"gram": str(gram), "n": str(n),
                                       "k": k},
                                 summary=f"representation of {n} by a "
                                         f"quadratic form in {k} variables"))
                scan = _diagonal_scan(P, R)
                if scan:
                    entries, _ = scan
                    if len(entries) == k and all(c == 1 for c, _, _ in entries):
                        if k == 3:
                            out.append(Match("sum-of-three-squares",
                                             data={"n": str(n)},
                                             summary=f"three squares: "
                                                     f"n = {n}"))
                        elif k == 4:
                            out.append(Match("sum-of-four-squares",
                                             data={"n": str(n)},
                                             summary=f"four squares: "
                                                     f"n = {n}"))
        else:
            out.append(Match("quadric",
                             data={"k": k},
                             summary="affine quadric: complete the square, "
                                     "then Hasse-Minkowski / "
                                     "Grunewald-Segal"))
        return out

    # --- degree >= 3 ---

    # Markov-Hurwitz: sum of squares = a * product of all variables
    prod_mon = prod(gens)
    expected = {g ** 2 for g in gens} | {prod_mon}
    mons = set(P.monomials())
    if k >= 3 and mons == expected and c0 == 0:
        sq = {P.monomial_coefficient(g ** 2) for g in gens}
        pcoeff = P.monomial_coefficient(prod_mon)
        if len(sq) == 1:
            s = sq.pop()
            if _sign_of(s) == -1:
                s, pcoeff = -s, -pcoeff
            a = -pcoeff
            az = _as_int(a)
            if s == 1 and az is not None and az > 0:
                out.append(Match(
                    "markov-hurwitz", data={"a": str(az), "k": k},
                    summary=("Markov equation" if (az == 3 and k == 3) else
                             f"Hurwitz equation x_1^2+...+x_{k}^2 = "
                             f"{az} x_1...x_{k}"),
                ))
                return out

    scan = _diagonal_scan(P, R)
    if scan:
        entries, _ = scan
        if len(entries) == k:
            exps = [e for _, _, e in entries]
            coeffs = [c for c, _, _ in entries]
            signs = [_sign_of(c) for c in coeffs]
            if c0 == 0 and k == 3 and all(e >= 2 for e in exps) \
                    and None not in signs:
                transform = ""
                if abs(sum(signs)) == 3:
                    # all terms on one side: for an odd exponent we may flip
                    # the sign of that variable (x -> -x)
                    flippable = [i for i, e in enumerate(exps) if e % 2 == 1]
                    if flippable:
                        i = flippable[0]
                        coeffs = list(coeffs)
                        signs = list(signs)
                        coeffs[i] = -coeffs[i]
                        signs[i] = -signs[i]
                        transform = (f"substituted {entries[i][1]} -> "
                                     f"-{entries[i][1]} (odd exponent)")
                if abs(sum(signs)) == 1:
                    # orient to a x^p + b y^q - c z^r = 0 with a, b, c > 0
                    if sum(signs) == -1:
                        coeffs = [-c for c in coeffs]
                        signs = [-s for s in signs]
                    triples = list(zip(coeffs, exps, signs))
                    pos = sorted((e, c) for c, e, s in triples if s > 0)
                    neg = [(e, -c) for c, e, s in triples if s < 0]
                    (p, ca), (q, cb) = pos
                    (r, cc) = neg[0]
                    chi = QQ(1) / p + QQ(1) / q + QQ(1) / r
                    regime = ("spherical" if chi > 1 else
                              "euclidean" if chi == 1 else "hyperbolic")
                    data = {"signature": f"({p}, {q}, {r})", "chi": str(chi),
                            "regime": regime, "a": str(ca), "b": str(cb),
                            "c": str(cc)}
                    out.append(Match(
                        "generalized-fermat", data=data, transform=transform,
                        summary=f"generalized Fermat equation of signature "
                                f"({p}, {q}, {r}), chi = {chi} ({regime})",
                    ))
                    if p == q == r >= 3 and abs(ca) == 1 and abs(cb) == 1 \
                            and abs(cc) == 1:
                        out.append(Match(
                            "fermat", data={"n": str(p)}, transform=transform,
                            summary=f"Fermat equation with exponent {p}",
                        ))
                    return out
            if c0 != 0 and len(set(exps)) == 1 and None not in signs:
                e = exps[0]
                n = -c0
                if e == 3 and k == 3 and all(abs(c) == 1 for c in coeffs):
                    out.append(Match(
                        "sum-of-three-cubes", data={"n": str(n)},
                        summary=f"sum of three cubes: n = {n}",
                    ))
                    return out
                if e >= 3 and all(c == 1 for c in coeffs):
                    out.append(Match(
                        "waring", data={"k": ZZ(e), "s": k, "n": str(n)},
                        summary=f"Waring-type: sum of {k} {e}-th powers "
                                f"= {n}",
                    ))
                    return out
                if e >= 3:
                    out.append(Match(
                        "diagonal-form",
                        data={"k": ZZ(e), "s": k,
                              "coeffs": [str(c) for c in coeffs],
                              "n": str(n)},
                        summary=f"diagonal equation of degree {e}",
                    ))
                    return out
            if c0 == 0 and k >= 4 and len(set(exps)) == 1 \
                    and None not in signs and exps[0] >= 3:
                e = exps[0]
                s_count = signs.count(1)
                t_count = signs.count(-1)
                if 0 < s_count and 0 < t_count:
                    out.append(Match(
                        "equal-sums-like-powers",
                        data={"k": ZZ(e), "s": min(s_count, t_count),
                              "t": max(s_count, t_count)},
                        summary=f"equal sums of {e}-th powers "
                                f"({min(s_count, t_count)} vs "
                                f"{max(s_count, t_count)} terms)",
                    ))
                    return out
                out.append(Match(
                    "diagonal-form",
                    data={"k": ZZ(e), "s": k,
                          "coeffs": [str(c) for c in coeffs], "n": "0"},
                    summary=f"diagonal form of degree {e} = 0"
                            + (" (only trivial real solutions)"
                               if e % 2 == 0 else ""),
                ))
                return out

    if k == 3 and P.is_homogeneous():
        if d == 3 and pe.is_concrete and _is_irreducible(P) is True:
            try:
                smooth = Curve(P).is_smooth()
            except Exception:
                smooth = None
            if smooth:
                out.append(Match(
                    "plane-cubic", data={"F": str(P)},
                    summary="smooth plane cubic: genus 1; find a point, "
                            "then reduce to Weierstrass form (Nagell)",
                ))
                return out
        routed = _genus_route(pe, P, affine=False)
        if routed:
            out.extend(routed)
            return out

    return out


# --------------------------------------------------------------------------
# exponential equations
# --------------------------------------------------------------------------

def _match_exponential(pe):
    r"""
    Matches for equations with exponential content.

    Dispatches on the composition of the terms: purely exponential (Pillai,
    S-unit), variable powers (Catalan, symbolic Fermat, Lebesgue–Nagell,
    Schinzel–Tijdeman power values), polynomial + one exponential
    (Ramanujan–Nagell, Thue–Mahler), with ``polynomial-exponential`` as the
    root fallback.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _match_exponential
        sage: [m.slug for m in _match_exponential(parse("x^2 + 7 = 2^n"))]
        ['polynomial-exponential', 'ramanujan-nagell']
        sage: [m.slug for m in _match_exponential(parse("x^p - y^q = 1"))]
        ['polynomial-exponential', 'catalan']
        sage: [m.slug for m in _match_exponential(parse("3^m - 2^n = 5"))]
        ['polynomial-exponential', 'exponential-diophantine', 'pillai']
    """
    terms = pe.terms
    pure_exp = [t for t in terms
                if t.exp_terms and not t.var_powers and not t.powers]
    var_pow = [t for t in terms if t.var_powers]
    poly_terms = [t for t in terms if t.is_polynomial]
    const_terms = [t for t in poly_terms if t.is_constant]
    nonconst_poly = [t for t in poly_terms if not t.is_constant]
    mixed = [t for t in terms
             if (t.exp_terms and t.powers) or (t.var_powers and t.powers)
             or (t.exp_terms and t.var_powers) or len(t.var_powers) > 1]

    out = [Match(
        "polynomial-exponential",
        data={"exp_terms": len(pure_exp) + len(var_pow) + len(mixed)},
        summary="mixed polynomial-exponential equation",
    )]

    c0 = sum((t.coeff for t in const_terms if not t.param_powers), QQ(0))
    param_const = [t for t in const_terms if t.param_powers]

    # ---- purely exponential: sum of c * prod(b^n) terms and constants
    if not var_pow and not nonconst_poly and not mixed and pure_exp:
        bases = sorted({b for t in pure_exp for b, _ in t.exp_terms},
                       key=str)
        out.append(Match(
            "exponential-diophantine",
            data={"terms": len(pure_exp), "bases": [str(b) for b in bases]},
            summary=f"purely exponential equation in bases {bases}",
        ))
        single = all(len(t.exp_terms) == 1 for t in pure_exp)
        if len(pure_exp) == 2 and single and not param_const and c0 != 0:
            t1, t2 = pure_exp
            s1, s2 = _sign_of(t1.coeff), _sign_of(t2.coeff)
            if abs(t1.coeff) == 1 and abs(t2.coeff) == 1 and s1 == -s2 \
                    and t1.exp_terms[0][1] != t2.exp_terms[0][1]:
                if s1 < 0:
                    t1, t2 = t2, t1
                (a, xvar), (b, yvar) = t1.exp_terms[0], t2.exp_terms[0]
                out.append(Match(
                    "pillai",
                    data={"a": str(a), "b": str(b), "c": str(-c0),
                          "roles": {"x": xvar, "y": yvar}},
                    summary=f"Pillai equation {a}^{xvar} - {b}^{yvar} "
                            f"= {-c0}",
                ))
            else:
                out.append(Match(
                    "s-unit", data={"bases": [str(b) for b in bases]},
                    summary="three-term S-unit-type equation",
                ))
        elif len(pure_exp) + (1 if (c0 != 0 or param_const) else 0) == 3:
            out.append(Match(
                "s-unit", data={"bases": [str(b) for b in bases]},
                summary="three-term S-unit-type equation "
                        f"(S generated by {bases})",
            ))
        return out

    # ---- Thue-Mahler: homogeneous binary form = m * prod(p^N)
    if not var_pow and not mixed and len(pure_exp) == 1 and nonconst_poly \
            and not const_terms:
        polyvars = sorted({v for t in nonconst_poly for v, _ in t.powers})
        degs = {sum(e for _, e in t.powers) for t in nonconst_poly}
        if len(polyvars) == 2 and len(degs) == 1:
            d = degs.pop()
            if d >= 3:
                t = pure_exp[0]
                m = -t.coeff
                primes = [b for b, _ in t.exp_terms]
                form = " + ".join(
                    f"{tt.coeff}*{'*'.join(f'{v}^{e}' for v, e in tt.powers)}"
                    for tt in nonconst_poly)
                out.append(Match(
                    "thue-mahler",
                    data={"form": form, "m": str(m), "degree": ZZ(d),
                          "primes": [str(p) for p in primes]},
                    summary=f"Thue-Mahler equation of degree {d} with "
                            f"primes {primes}",
                ))
                return out

    # ---- var^var families
    simple_vv = [t for t in var_pow
                 if len(t.var_powers) == 1 and not t.powers
                 and not t.exp_terms]
    if len(simple_vv) == len(var_pow) and var_pow and not pure_exp \
            and not mixed:
        if len(var_pow) == 2 and not nonconst_poly and not param_const:
            t1, t2 = var_pow
            s1, s2 = _sign_of(t1.coeff), _sign_of(t2.coeff)
            if abs(t1.coeff) == 1 and abs(t2.coeff) == 1 and s1 == -s2 \
                    and c0 != 0:
                if s1 < 0:
                    t1, t2 = t2, t1
                (xb, pv), (yb, qv) = t1.var_powers[0], t2.var_powers[0]
                c = -c0
                if xb != yb and pv != qv and abs(c) == 1:
                    if c == -1:
                        (xb, pv), (yb, qv) = (yb, qv), (xb, pv)
                    out.append(Match(
                        "catalan",
                        data={"roles": {"x": xb, "p": pv, "y": yb,
                                        "q": qv}},
                        summary=f"Catalan equation {xb}^{pv} - {yb}^{qv} "
                                "= 1: only 3^2 - 2^3 = 1 (Mihailescu)",
                    ))
                    return out
                if xb != yb:
                    out.append(Match(
                        "pillai",
                        data={"c": str(c), "variable_bases": True,
                              "roles": {"x": xb, "p": pv, "y": yb,
                                        "q": qv}},
                        summary=f"perfect-power difference {xb}^{pv} - "
                                f"{yb}^{qv} = {c} (Pillai's conjecture)",
                    ))
                    return out
        if len(var_pow) == 3 and not nonconst_poly and c0 == 0 \
                and not param_const:
            coeffs = [t.coeff for t in var_pow]
            signs = [_sign_of(c) for c in coeffs]
            expvars = [t.var_powers[0][1] for t in var_pow]
            bases = [t.var_powers[0][0] for t in var_pow]
            if all(abs(c) == 1 for c in coeffs) and None not in signs \
                    and abs(sum(signs)) == 1 and len(set(bases)) == 3:
                if len(set(expvars)) == 1:
                    out.append(Match(
                        "fermat", data={"n": expvars[0], "symbolic": True},
                        summary=f"Fermat equation with unknown exponent "
                                f"{expvars[0]} (no solutions for "
                                f"{expvars[0]} >= 3, Wiles)",
                    ))
                    return out
                if len(set(expvars)) == 3:
                    out.append(Match(
                        "generalized-fermat",
                        data={"signature": f"({', '.join(expvars)})",
                              "symbolic": True},
                        summary="generalized Fermat with unknown exponents "
                                "(Beal/Fermat-Catalan territory)",
                    ))
                    return out
        if len(var_pow) == 1 and nonconst_poly and not param_const:
            t = var_pow[0]
            base, expvar = t.var_powers[0]
            polyvars = sorted({v for tt in nonconst_poly
                               for v, _ in tt.powers})
            if len(polyvars) == 1 and polyvars[0] != base \
                    and abs(t.coeff) == 1:
                u = polyvars[0]
                # f(u) = c * base^expvar with c = -t.coeff
                fcoeffs = {}
                for tt in poly_terms:
                    e = tt.powers[0][1] if tt.powers else 0
                    fcoeffs[e] = fcoeffs.get(e, QQ(0)) + tt.coeff
                if _sign_of(t.coeff) == 1:
                    fcoeffs = {e: -c for e, c in fcoeffs.items()}
                degf = max(fcoeffs)
                fstr = " + ".join(f"({c})*{u}^{e}" if e else f"({c})"
                                  for e, c in sorted(fcoeffs.items(),
                                                     reverse=True))
                if degf == 2 and fcoeffs.get(2) == 1 and not fcoeffs.get(1):
                    dval = fcoeffs.get(0, QQ(0))
                    out.append(Match(
                        "lebesgue-nagell",
                        data={"d": str(dval),
                              "roles": {"x": u, "y": base, "n": expvar}},
                        summary=f"Lebesgue-Nagell equation {u}^2 + ({dval}) "
                                f"= {base}^{expvar}",
                    ))
                    return out
                if degf >= 2:
                    out.append(Match(
                        "power-values",
                        data={"f": fstr,
                              "roles": {"x": u, "y": base, "n": expvar}},
                        summary=f"power values of a polynomial: {fstr} = "
                                f"{base}^{expvar} (Schinzel-Tijdeman)",
                    ))
                    return out

    # ---- polynomial part + one true exponential: Ramanujan-Nagell type
    if not var_pow and not mixed and len(pure_exp) == 1 and nonconst_poly:
        polyvars = sorted({v for t in nonconst_poly for v, _ in t.powers})
        if len(polyvars) == 1:
            u = polyvars[0]
            t = pure_exp[0]
            k_coeff = -t.coeff
            fcoeffs = {}
            for tt in poly_terms:
                e = tt.powers[0][1] if tt.powers else 0
                fcoeffs[e] = fcoeffs.get(e, QQ(0)) + tt.coeff
            degf = max(fcoeffs)
            if degf == 2:
                if _sign_of(fcoeffs[2]) == -1:
                    fcoeffs = {e: -c for e, c in fcoeffs.items()}
                    k_coeff = -k_coeff
                base_desc = "*".join(f"{b}^{n}" for b, n in t.exp_terms)
                dval = fcoeffs.get(0, QQ(0))
                monic_pure = (fcoeffs.get(2) == 1 and not fcoeffs.get(1))
                data = {"k": str(k_coeff), "d": str(dval),
                        "base": str(t.exp_terms[0][0]),
                        "exp": base_desc,
                        "roles": {"x": u, "n": t.exp_terms[0][1]}}
                classical = (monic_pure and k_coeff == 1 and dval == 7
                             and len(t.exp_terms) == 1
                             and t.exp_terms[0][0] == 2)
                out.append(Match(
                    "ramanujan-nagell", data=data,
                    summary=("the Ramanujan-Nagell equation x^2 + 7 = 2^n: "
                             "n in {3, 4, 5, 7, 15}" if classical else
                             f"generalized Ramanujan-Nagell: quadratic in "
                             f"{u} = ({k_coeff})*{base_desc}"),
                ))
                return out

    return out


# --------------------------------------------------------------------------
# unit fractions
# --------------------------------------------------------------------------

def _match_fractional(pe):
    r"""
    Matches from the unit-fraction structure of the *uncleared* input.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import _match_fractional
        sage: pe = parse("4/n = 1/x + 1/y + 1/z", params="n")
        sage: [m.slug for m in _match_fractional(pe)]
        ['egyptian-fractions', 'erdos-straus']
        sage: _match_fractional(parse("x^2 + y^2 = 610"))
        []
    """
    fs = pe.fractional
    out = []
    if fs is None:
        return out
    k, a, n = fs["k"], fs["a"], fs["n"]
    out.append(Match(
        "egyptian-fractions",
        data={"k": ZZ(k), "a": str(a), "n": str(n),
              "unit_vars": list(fs["unit_vars"])},
        summary=f"sum of {k} unit fractions = {a}/{n}",
    ))
    if k == 3 and a == 4:
        out.append(Match(
            "erdos-straus", data={"n": str(n), "a": "4", "k": "3"},
            summary=f"Erdos-Straus equation 4/{n} = 1/x + 1/y + 1/z",
        ))
    elif k == 3 and a == 5:
        out[-1].summary += " (Sierpinski's 5/n problem)"
    return out


def run(pe):
    r"""
    Collect matches from all applicable matchers, deduplicated by slug.

    INPUT:

    - ``pe`` -- a :class:`~diophantine_classifier.parsing.ParsedEquation`

    OUTPUT: list of :class:`Match`, in emission order (ranking by
    specificity happens in :mod:`~diophantine_classifier.classify`)

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import run
        sage: [m.slug for m in run(parse("y^2 = x^3 - 2"))]
        ['general-polynomial', 'elliptic-weierstrass', 'mordell']
        sage: [m.slug for m in run(parse("2^a + 3^b = 5^c"))]
        ['polynomial-exponential', 'exponential-diophantine', 's-unit']
    """
    matches = []
    matches.extend(_match_fractional(pe))
    if pe.is_polynomial:
        matches.extend(_match_polynomial(pe))
    else:
        matches.extend(_match_exponential(pe))
    seen = set()
    result = []
    for m in matches:
        if m.slug not in seen:
            seen.add(m.slug)
            result.append(m)
    return result
