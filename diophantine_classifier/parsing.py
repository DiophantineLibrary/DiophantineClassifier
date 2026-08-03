r"""Parse a Diophantine equation string into a structured term model.

The model supports polynomial terms and three kinds of exponential content:

* exponential terms: integer (or parameter) base with an unknown exponent,
  e.g. ``2^n``;
* variable powers: unknown base with unknown exponent, e.g. ``y^q`` in
  Catalan's equation;
* ordinary monomials in the unknowns, with coefficients that may involve the
  declared parameters (e.g. ``D`` in ``x^2 - D*y^2 = 1``).

Everything is normalized to ``f = 0`` with denominators cleared (recording
the nonvanishing conditions this introduces), terms merged, and — in the
purely polynomial case — ``f`` realized in a Sage polynomial ring over
``QQ`` or ``QQ[params]``.

Every identifier in the input becomes a variable: parsing supplies an
explicit symbol table, so Sage's global names (the constant ``e``, the
imaginary unit ``I``, interface objects like ``r`` and ``gp``) are never
consulted and cannot leak into equations.

EXAMPLES::

    sage: from diophantine_classifier.parsing import parse
    sage: pe = parse("x^2 - 61*y^2 = 1")
    sage: pe.is_polynomial, pe.unknowns
    (True, ('x', 'y'))
    sage: pe.poly
    x^2 - 61*y^2 - 1

    sage: pe = parse("x^2 + 7 = 2^n")           # variable exponent
    sage: pe.is_polynomial, pe.exp_unknowns
    (False, ('n',))

    sage: pe = parse("y^2 = x^3 + k", params="k")
    sage: pe.params, pe.base_ring
    (('k',), Univariate Polynomial Ring in k over Rational Field)
"""

import re
from dataclasses import dataclass, field

import operator as _op

from sage.all import QQ, ZZ, SR, PolynomialRing, lcm, prod
from sage.symbolic.operators import add_vararg, mul_vararg


class ParseError(ValueError):
    r"""The input string could not be interpreted as an equation.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: parse("x - x = 0")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.ParseError: equation simplifies to 0 = 0
    """


class UnsupportedEquationError(ParseError):
    r"""Parsed fine, but uses features outside the supported equation model.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: parse("sin(x) = 0")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.UnsupportedEquationError: sin(...) looks like a function call...
    """


_IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def _preprocess(s):
    r"""Normalize an input string before symbolic parsing.

    Handles unicode operators, ``**`` vs ``^``, ``==`` vs ``=``, and the
    unambiguous cases of implicit multiplication (digit-letter,
    digit-parenthesis, parenthesis-letter).

    INPUT:

    - ``s`` -- string

    OUTPUT: normalized string

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _preprocess
        sage: _preprocess("x**2 == 61y^2")
        'x^2 = 61*y^2'
        sage: _preprocess("3(x+1) − 2·y")
        '3*(x+1) - 2*y'
    """
    s = s.strip()
    s = s.replace("−", "-").replace("·", "*").replace("×", "*")
    s = s.replace("**", "^")
    s = s.replace("==", "=")
    s = re.sub(r"(\d)\s*\(", r"\1*(", s)
    s = re.sub(r"\)\s*([A-Za-z0-9_(])", r")*\1", s)
    s = re.sub(r"(\d)\s+([A-Za-z_])", r"\1*\2", s)  # "61 y^2" -> "61*y^2"
    s = re.sub(r"(\d)([A-Za-z_])", r"\1*\2", s)     # "61y^2"  -> "61*y^2"
    return s


def _check_identifiers(s):
    r"""Reject function-call syntax and return the set of identifiers.

    Every identifier is treated as a variable — see the module docstring.

    INPUT:

    - ``s`` -- preprocessed equation string

    OUTPUT: set of identifier strings

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _check_identifiers
        sage: sorted(_check_identifiers("x^2 - 61*y^2 = 1"))
        ['x', 'y']
        sage: _check_identifiers("f(x) = 0")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.UnsupportedEquationError: f(...) looks like a function call...
    """
    names = set()
    for m in _IDENT.finditer(s):
        name = m.group(0)
        names.add(name)
        if s[m.end():m.end() + 1] == "(":
            raise UnsupportedEquationError(
                f"{name}(...) looks like a function call; functions are not "
                "supported (for a product, write an explicit '*')"
            )
    return names


@dataclass(frozen=True)
class Term:
    r"""One additive term of a normalized equation.

    Represents ``coeff · Π params^e · Π unknowns^e · Π base^expvar ·
    Π var^expvar``.

    ATTRIBUTES:

    - ``coeff`` -- rational number (element of ``QQ``); the numerical
      coefficient.
    - ``param_powers`` -- sorted tuple of pairs ``(parameter name, positive
      integer exponent)``; the parameter part of the coefficient (e.g.
      ``(("D", 1),)`` for the term ``-D*y^2``).
    - ``powers`` -- sorted tuple of pairs ``(unknown name, positive integer
      exponent)``; the ordinary monomial part.
    - ``exp_terms`` -- sorted tuple of pairs ``(base, exponent unknown
      name)`` with ``base`` an integer of absolute value ≥ 2 or a parameter
      name; exponential factors such as ``2^n``.
    - ``var_powers`` -- sorted tuple of pairs ``(base unknown name, exponent
      unknown name)``; variable-power factors such as ``y^q``.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: t = parse("x^2 + 7 = 2^n").exponential_terms()[0]
        sage: t.coeff, t.exp_terms
        (-1, ((2, 'n'),))
        sage: t.is_polynomial
        False
    """
    coeff: object
    param_powers: tuple = ()
    powers: tuple = ()
    exp_terms: tuple = ()
    var_powers: tuple = ()

    @property
    def signature(self):
        r"""The term's structure with the coefficient stripped (merge key).

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: t = parse("3*x^2 = 1").terms[0]
            sage: t.signature
            ((), (('x', 2),), (), ())
        """
        return (self.param_powers, self.powers, self.exp_terms,
                self.var_powers)

    @property
    def is_polynomial(self):
        r"""Whether the term is free of exponential content.

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: [t.is_polynomial for t in parse("x^2 + 7 = 2^n").terms]
            [True, False, True]
        """
        return not self.exp_terms and not self.var_powers

    @property
    def is_constant(self):
        r"""Whether the term involves no unknowns at all (parameters allowed).

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: [t.is_constant for t in parse("x^2 + 7 = 2^n").terms]
            [False, False, True]
        """
        return not self.powers and self.is_polynomial

    @property
    def is_pure_power(self):
        r"""Whether the term is ``c * v^e`` for a single unknown and ``e ≥ 2``.

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: [t.is_pure_power for t in parse("x^3 + x*y = 1").terms]
            [True, False, False]
        """
        return self.is_polynomial and len(self.powers) == 1 \
            and self.powers[0][1] >= 2


@dataclass
class ParsedEquation:
    r"""A fully parsed and normalized Diophantine equation.

    Produced by :func:`parse`; consumed by the matchers and solvers.

    ATTRIBUTES:

    - ``original`` -- string; the input as given by the user.
    - ``lhs``, ``rhs`` -- strings; the two sides after preprocessing.
    - ``params`` -- tuple of strings; symbols treated as parameters (they
      live in the coefficient ring, not among the unknowns).
    - ``domain`` -- string in ``{"ZZ", "NN", "QQ"}``; where solutions are
      sought.  Threaded through to classification output.
    - ``terms`` -- list of :class:`Term`; the equation as ``sum(terms) = 0``
      after clearing denominators and merging.
    - ``unknowns`` -- tuple of strings; all non-parameter symbols, in order
      of first appearance in the input.  Solution tuples follow this order.
    - ``exp_unknowns`` -- tuple of strings; the unknowns that occur as
      exponents.
    - ``conditions`` -- list of strings; hypotheses introduced by
      normalization (e.g. denominators assumed nonzero).  Never silently
      dropped: the classification reports them.
    - ``fractional`` -- dict or ``None``; the unit-fraction structure of the
      *uncleared* input when it has one (keys ``k``, ``unit_vars``, ``a``,
      ``n``), used by the Egyptian-fraction matcher.
    - ``poly`` -- Sage polynomial or ``None``; for purely polynomial
      equations, the cleared, content-normalized polynomial whose vanishing
      is the equation.
    - ``poly_ring`` -- its parent, a polynomial ring in ``unknowns`` over
      ``base_ring``.
    - ``base_ring`` -- ``QQ``, or ``QQ[params]`` when parameters are present.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: pe = parse("4/n = 1/x + 1/y + 1/z", params="n")
        sage: pe.unknowns
        ('x', 'y', 'z')
        sage: pe.fractional["a"], pe.fractional["k"]
        (4, 3)
        sage: pe.conditions
        ['multiplied by n*x*y*z (valid where n, x, y, z nonzero)']
    """
    original: str
    lhs: str
    rhs: str
    params: tuple
    domain: str
    terms: list
    unknowns: tuple
    exp_unknowns: tuple
    conditions: list = field(default_factory=list)
    fractional: dict = None
    poly = None
    poly_ring = None
    base_ring = None

    @property
    def is_polynomial(self):
        r"""Whether every term is polynomial (no variable exponents).

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: parse("x^2 = 2").is_polynomial
            True
            sage: parse("x^2 = 2^n").is_polynomial
            False
        """
        return all(t.is_polynomial for t in self.terms)

    @property
    def is_concrete(self):
        r"""Whether there are no parameters (coefficients are rationals).

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: parse("y^2 = x^3 - 2").is_concrete
            True
            sage: parse("y^2 = x^3 + k", params="k").is_concrete
            False
        """
        return not self.params

    @property
    def nunknowns(self):
        r"""The number of unknowns.

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: parse("x^2 + y^2 = z^2").nunknowns
            3
        """
        return len(self.unknowns)

    def poly_part_terms(self):
        r"""The polynomial terms (including constants).

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: len(parse("x^2 + 7 = 2^n").poly_part_terms())
            2
        """
        return [t for t in self.terms if t.is_polynomial]

    def exponential_terms(self):
        r"""The terms with exponential content.

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: len(parse("x^2 + 7 = 2^n").exponential_terms())
            1
        """
        return [t for t in self.terms if not t.is_polynomial]

    def __repr__(self):
        r"""Terse representation.

        EXAMPLES::

            sage: from diophantine_classifier.parsing import parse
            sage: parse("x^2 - 61*y^2 = 1")
            ParsedEquation('x^2 - 61*y^2 = 1')
        """
        return f"ParsedEquation({self.original!r})"


def _split_equation(s):
    r"""Split a preprocessed string at its (single) equality sign.

    A string without ``=`` is interpreted as ``s = 0``.

    OUTPUT: pair of side strings

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _split_equation
        sage: _split_equation("x^2 = 2")
        ('x^2', '2')
        sage: _split_equation("x^2 - 2")
        ('x^2 - 2', '0')
        sage: _split_equation("x = y = z")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.UnsupportedEquationError: multiple '=' signs...
    """
    parts = s.split("=")
    if len(parts) == 1:
        return s, "0"
    if len(parts) == 2:
        lhs, rhs = parts[0].strip(), parts[1].strip()
        if not lhs or not rhs:
            raise ParseError("empty side of equation")
        return lhs, rhs
    raise UnsupportedEquationError(
        "multiple '=' signs: systems of equations are not yet supported"
    )


def _term_list(f):
    r"""The additive terms of a symbolic expression.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _term_list
        sage: x, y = SR.var("x"), SR.var("y")
        sage: _term_list(x + 2*y)
        [x, 2*y]
        sage: _term_list(x), _term_list(SR(0))
        ([x], [])
    """
    if f.operator() is add_vararg:
        return list(f.operands())
    if f.is_zero():
        return []
    return [f]


def _factor_list(t):
    r"""The multiplicative factors of a symbolic term.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _factor_list
        sage: x, y = SR.var("x"), SR.var("y")
        sage: _factor_list(2*x*y^2)
        [x, y^2, 2]
        sage: _factor_list(x)
        [x]
    """
    if t.operator() is mul_vararg:
        return list(t.operands())
    return [t]


def _as_QQ(ex):
    r"""Coerce a symbolic expression to ``QQ``, or return ``None``.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _as_QQ
        sage: _as_QQ(SR(3)/2)
        3/2
        sage: _as_QQ(SR.var("x")) is None
        True
    """
    try:
        return QQ(ex)
    except (TypeError, ValueError):
        return None


def _as_ZZ(ex):
    r"""Coerce a symbolic expression to ``ZZ``, or return ``None``.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _as_ZZ
        sage: _as_ZZ(SR(-7))
        -7
        sage: _as_ZZ(SR(1)/2) is None
        True
    """
    try:
        return ZZ(ex)
    except (TypeError, ValueError):
        return None


def _fractional_structure(f0, params):
    r"""Detect a 'sum of unit fractions = a/n' structure in the uncleared
    input.

    INPUT:

    - ``f0`` -- the expanded symbolic expression ``lhs - rhs`` *before*
      clearing denominators
    - ``params`` -- tuple of parameter names

    OUTPUT: dict with keys ``k`` (number of unit fractions), ``unit_vars``
    (their variables), ``a`` and ``n`` (the target ``a/n``, with ``n`` an
    integer or a parameter name) — or ``None`` if the shape doesn't match

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: parse("1/x + 1/y + 1/z = 1").fractional["n"]
        1
        sage: parse("4/n = 1/x + 1/y + 1/z", params="n").fractional["n"]
        'n'
        sage: parse("x^2 = 2").fractional is None
        True
    """
    units = []
    consts = []  # (QQ value, param name or None)
    for t in _term_list(f0):
        coeff = QQ(1)
        varinv = None
        paraminv = None
        ok = True
        for fac in _factor_list(t):
            q = _as_QQ(fac)
            if q is not None:
                coeff *= q
                continue
            if fac.operator() is _op.pow:
                base, expo = fac.operands()
                if _as_ZZ(expo) == -1 and base.is_symbol():
                    name = str(base)
                    if name in params:
                        if paraminv is None:
                            paraminv = name
                        else:
                            ok = False
                    elif varinv is None:
                        varinv = name
                    else:
                        ok = False
                    continue
            ok = False
        if not ok:
            return None
        if varinv is not None:
            if paraminv is not None or abs(coeff) != 1:
                return None
            units.append((1 if coeff > 0 else -1, varinv))
        else:
            consts.append((coeff, paraminv))
    if len(units) < 2 or len(consts) != 1:
        return None
    signs = {s for s, _ in units}
    if len(signs) != 1:
        return None
    sign = signs.pop()
    cval, cparam = consts[0]
    if cval == 0 or (cval > 0) == (sign > 0):
        return None  # target must sit on the other side of the equation
    a = abs(cval)
    if cparam is None:
        if a.denominator() == 1:
            target_n = 1
            target_a = ZZ(a)
        else:
            target_a, target_n = a.numerator(), a.denominator()
    else:
        if a.denominator() != 1:
            return None
        target_a, target_n = ZZ(a), cparam
    return {
        "k": len(units),
        "unit_vars": tuple(v for _, v in units),
        "a": target_a,
        "n": target_n,
    }


def _canonicalize_exp(base, expo):
    r"""Handle exponents like ``3*n`` by absorbing constants into the base.

    INPUT:

    - ``base``, ``expo`` -- the operands of a symbolic power with
      non-integer exponent

    OUTPUT: pair ``(base, expo)`` with ``expo`` a plain symbol

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: pe = parse("x + 1 = 2^(3*n)")   # exponent 3n: base becomes 8
        sage: [t.exp_terms for t in pe.exponential_terms()]
        [((8, 'n'),)]
        sage: parse("x = 2^(n*m)")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.UnsupportedEquationError: unsupported exponent...
    """
    if expo.is_symbol():
        return base, expo
    if expo.operator() is mul_vararg:
        ops = expo.operands()
        if len(ops) == 2:
            c0, c1 = _as_ZZ(ops[0]), _as_ZZ(ops[1])
            if c0 is not None and c0 > 0 and ops[1].is_symbol():
                return base ** c0, ops[1]
            if c1 is not None and c1 > 0 and ops[0].is_symbol():
                return base ** c1, ops[0]
    raise UnsupportedEquationError(
        f"unsupported exponent {expo} (only 'c*n' with c a positive integer "
        "is folded into the base)"
    )


def _walk_terms(f, params):
    r"""Decompose a cleared symbolic expression into :class:`Term` objects.

    Terms with identical structure are merged; terms with zero coefficient
    are dropped.

    INPUT:

    - ``f`` -- expanded symbolic expression with no denominators
    - ``params`` -- tuple of parameter names

    OUTPUT: list of :class:`Term`

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _walk_terms
        sage: x = SR.var("x")
        sage: sorted((t.coeff, t.powers) for t in _walk_terms(x^2 + 3*x - x, ()))
        [(1, (('x', 2),)), (2, (('x', 1),))]

    TESTS:

    Zero terms are dropped::

        sage: _walk_terms(x - x, ())
        []
    """
    terms = []
    for t in _term_list(f):
        coeff = QQ(1)
        param_powers = {}
        powers = {}
        exp_terms = []
        var_powers = []
        for fac in _factor_list(t):
            q = _as_QQ(fac)
            if q is not None:
                coeff *= q
                continue
            if fac.is_symbol():
                name = str(fac)
                target = param_powers if name in params else powers
                target[name] = target.get(name, 0) + 1
                continue
            if fac.operator() is _op.pow:
                base, expo = fac.operands()
                e = _as_ZZ(expo)
                if e is not None:
                    if e <= 0:
                        raise UnsupportedEquationError(
                            f"unexpected exponent {e} in {fac} after "
                            "clearing denominators"
                        )
                    if base.is_symbol():
                        name = str(base)
                        target = param_powers if name in params else powers
                        target[name] = target.get(name, 0) + e
                        continue
                    qb = _as_QQ(base)
                    if qb is not None:
                        coeff *= qb ** e
                        continue
                    raise UnsupportedEquationError(f"unsupported factor {fac}")
                # symbolic exponent
                base, expo = _canonicalize_exp(base, expo)
                expname = str(expo)
                if base.is_symbol():
                    bname = str(base)
                    if bname in params:
                        exp_terms.append((bname, expname))
                    else:
                        var_powers.append((bname, expname))
                    continue
                zb = _as_ZZ(base)
                if zb is None:
                    raise UnsupportedEquationError(
                        f"unsupported exponential base in {fac} (integer or "
                        "single-parameter bases only)"
                    )
                if abs(zb) <= 1:
                    raise UnsupportedEquationError(
                        f"exponential base {zb} in {fac}: bases 0, 1, -1 are "
                        "degenerate; simplify the equation first"
                    )
                exp_terms.append((zb, expname))
                continue
            raise UnsupportedEquationError(
                f"unsupported factor {fac} in term {t}")
        if coeff == 0:
            continue
        terms.append(Term(
            coeff=coeff,
            param_powers=tuple(sorted(param_powers.items())),
            powers=tuple(sorted(powers.items())),
            exp_terms=tuple(sorted(exp_terms, key=str)),
            var_powers=tuple(sorted(var_powers)),
        ))
    # merge terms with identical structure
    merged = {}
    order = []
    for t in terms:
        key = t.signature
        if key in merged:
            merged[key] = Term(merged[key].coeff + t.coeff, *key)
        else:
            merged[key] = t
            order.append(key)
    return [merged[k] for k in order if merged[k].coeff != 0]


def _appearance_order(s, names):
    r"""Sort ``names`` by first appearance in the string ``s``.

    Determines the ordering of the unknowns, hence of solution tuples.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import _appearance_order
        sage: _appearance_order("y^2 = x^3 - 2", {"x", "y"})
        ('y', 'x')
    """
    firsts = {}
    for m in _IDENT.finditer(s):
        name = m.group(0)
        if name in names and name not in firsts:
            firsts[name] = m.start()
    return tuple(sorted(names, key=lambda n: firsts.get(n, len(s))))


def parse(equation, params=(), domain="ZZ"):
    r"""Parse an equation string into a :class:`ParsedEquation`.

    INPUT:

    - ``equation`` -- string; one equation.  ``=`` and ``==`` both work; a
      missing right-hand side means ``= 0``; ``^`` and ``**`` both denote
      powers; implicit multiplication is resolved where unambiguous
      (``61y^2`` means ``61*y^2``).  Functions and systems are not (yet)
      supported.
    - ``params`` -- (default: ``()``) names of symbols to treat as
      parameters rather than unknowns; either an iterable of strings or a
      comma/space-separated string
    - ``domain`` -- (default: ``"ZZ"``) one of ``"ZZ"``, ``"NN"``, ``"QQ"``:
      where solutions are sought (recorded, and used by solvers)

    OUTPUT: a :class:`ParsedEquation`

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: parse("x^2 - 61*y^2 = 1").poly
        x^2 - 61*y^2 - 1
        sage: parse("x^2 - D*y^2 = 1", params="D").base_ring
        Univariate Polynomial Ring in D over Rational Field
        sage: parse("x^2 + 7 = 2^n").exp_unknowns
        ('n',)

    Denominators are cleared, with the introduced hypothesis recorded::

        sage: pe = parse("1/x + 1/y = 1/2")
        sage: pe.poly
        -x*y + 2*x + 2*y
        sage: pe.conditions
        ['multiplied by 2*x*y (valid where x, y nonzero)']

    Sage's global names cannot leak in — ``e``, ``r``, ``I`` are ordinary
    variables here::

        sage: parse("e^2 + r^2 = I^2").unknowns
        ('e', 'r', 'I')

    TESTS::

        sage: parse("x = y = z")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.UnsupportedEquationError: multiple '=' signs...
        sage: parse("x^2 = 1", domain="RR")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.ParseError: unknown domain 'RR'
        sage: parse("3 = 5", params="")
        Traceback (most recent call last):
        ...
        diophantine_classifier.parsing.ParseError: no unknowns found...
    """
    if isinstance(params, str):
        params = tuple(p.strip() for p in re.split(r"[,\s]+", params)
                       if p.strip())
    else:
        params = tuple(str(p) for p in params)
    if domain not in ("ZZ", "NN", "QQ"):
        raise ParseError(f"unknown domain {domain!r}")
    s = _preprocess(str(equation))
    idents = _check_identifiers(s)
    lhs, rhs = _split_equation(s)
    syms = {name: SR.var(name) for name in idents}
    try:
        from sage.calculus.calculus import symbolic_expression_from_string
        L = symbolic_expression_from_string(lhs, syms)
        R = symbolic_expression_from_string(rhs, syms)
        f0 = (L - R).expand()
    except (TypeError, ValueError, SyntaxError, RuntimeError) as err:
        raise ParseError(f"could not parse {equation!r}: {err}") from None
    if f0.is_zero():
        raise ParseError("equation simplifies to 0 = 0")

    fractional = _fractional_structure(f0, params)

    conditions = []
    f = f0
    den = f.denominator()
    if not den.is_constant():
        f = f.simplify_rational()
        den = f.denominator()
        f = f.numerator().expand()
        denvars = sorted(str(v) for v in den.variables())
        conditions.append(
            "multiplied by %s (valid where %s nonzero)"
            % (den, ", ".join(denvars))
        )

    terms = _walk_terms(f, params)
    if not terms:
        raise ParseError("equation simplifies to 0 = 0")

    names = set()
    exp_names = set()
    for t in terms:
        names.update(n for n, _ in t.powers)
        for base, expname in t.exp_terms:
            exp_names.add(expname)
        for bname, expname in t.var_powers:
            names.add(bname)
            exp_names.add(expname)
    exp_unknowns = tuple(sorted(n for n in exp_names if n not in params))
    names.update(exp_unknowns)
    if not names:
        raise ParseError("no unknowns found (all symbols are parameters)")
    unknowns = _appearance_order(s, names)

    pe = ParsedEquation(
        original=str(equation), lhs=lhs, rhs=rhs, params=params,
        domain=domain, terms=terms, unknowns=unknowns,
        exp_unknowns=exp_unknowns, conditions=conditions,
        fractional=fractional,
    )

    if pe.is_polynomial:
        base = QQ if not params else PolynomialRing(QQ, list(params))
        ring = PolynomialRing(base, len(unknowns), list(unknowns))
        gens = dict(zip(unknowns, ring.gens()))
        if params:
            pgens = dict(zip(params, base.gens()))
        poly = ring.zero()
        denominator = lcm([t.coeff.denominator() for t in terms])
        for t in terms:
            c = base(t.coeff * denominator)
            if params and t.param_powers:
                c *= prod(pgens[p] ** e for p, e in t.param_powers)
            mon = prod([gens[v] ** e for v, e in t.powers], ring.one())
            poly += ring(c) * mon
        if not params:
            content = poly.content()
            if content:
                poly = poly / content
        pe.poly = poly
        pe.poly_ring = ring
        pe.base_ring = base
    return pe
