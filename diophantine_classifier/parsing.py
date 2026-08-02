"""Parse a Diophantine equation string into a structured term model.

The model supports polynomial terms and three kinds of exponential content:

* ``ExpTerm``: integer (or parameter) base with an unknown exponent, e.g. ``2^n``;
* var-power pairs: unknown base with unknown exponent, e.g. ``y^q`` in
  Catalan's equation;
* ordinary monomials in the unknowns, with coefficients that may involve the
  declared parameters (e.g. ``D`` in ``x^2 - D*y^2 = 1``).

Everything is normalized to ``f = 0`` with denominators cleared (recording the
nonvanishing conditions this introduces), terms merged, and — in the purely
polynomial case — ``f`` realized in a Sage polynomial ring over ``QQ`` or
``QQ[params]``.
"""

import re
from dataclasses import dataclass, field

import operator as _op

from sage.all import QQ, ZZ, SR, PolynomialRing, lcm, prod
from sage.symbolic.operators import add_vararg, mul_vararg


class ParseError(ValueError):
    """The input string could not be interpreted as an equation."""


class UnsupportedEquationError(ParseError):
    """Parsed fine, but uses features outside the supported equation model."""


_IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def _preprocess(s):
    s = s.strip()
    s = s.replace("−", "-").replace("·", "*").replace("×", "*")
    s = s.replace("**", "^")
    s = s.replace("==", "=")
    # implicit multiplication we can resolve unambiguously
    s = re.sub(r"(\d)\s*\(", r"\1*(", s)
    s = re.sub(r"\)\s*([A-Za-z0-9_(])", r")*\1", s)
    s = re.sub(r"(\d)\s+([A-Za-z_])", r"\1*\2", s)  # "61 y^2" -> "61*y^2"
    s = re.sub(r"(\d)([A-Za-z_])", r"\1*\2", s)     # "61y^2"  -> "61*y^2"
    return s


def _check_identifiers(s):
    """Reject function-call syntax and return the set of identifiers.

    Every identifier is treated as a variable: parsing supplies an explicit
    symbol table, so Sage constants (``e``, ``pi``, ``I``) and interface
    objects (``r``, ``gp``) in the global namespace are never consulted.
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
    """coeff * prod(params^e) * prod(unknowns^e) * prod(base^expvar) * prod(var^expvar)"""
    coeff: object                # element of QQ
    param_powers: tuple = ()     # ((param name, positive int), ...)
    powers: tuple = ()           # ((unknown name, positive int), ...)
    exp_terms: tuple = ()        # ((base: Integer or param name, exponent unknown name), ...)
    var_powers: tuple = ()       # ((base unknown name, exponent unknown name), ...)

    @property
    def signature(self):
        return (self.param_powers, self.powers, self.exp_terms, self.var_powers)

    @property
    def is_polynomial(self):
        return not self.exp_terms and not self.var_powers

    @property
    def is_constant(self):
        return not self.powers and self.is_polynomial

    @property
    def is_pure_power(self):
        """A single unknown raised to a power >= 2, times a coefficient."""
        return self.is_polynomial and len(self.powers) == 1 and self.powers[0][1] >= 2


@dataclass
class ParsedEquation:
    original: str
    lhs: str
    rhs: str
    params: tuple
    domain: str
    terms: list
    unknowns: tuple            # all non-parameter symbols, in order of appearance
    exp_unknowns: tuple        # unknowns appearing as exponents
    conditions: list = field(default_factory=list)
    fractional: dict = None    # unit-fraction structure of the original, if any
    poly = None                # Sage polynomial (cleared, primitive) when polynomial
    poly_ring = None
    base_ring = None

    @property
    def is_polynomial(self):
        return all(t.is_polynomial for t in self.terms)

    @property
    def is_concrete(self):
        """No parameters: coefficients are honest rationals."""
        return not self.params

    @property
    def nunknowns(self):
        return len(self.unknowns)

    def poly_part_terms(self):
        return [t for t in self.terms if t.is_polynomial]

    def exponential_terms(self):
        return [t for t in self.terms if not t.is_polynomial]

    def __repr__(self):
        return f"ParsedEquation({self.original!r})"


def _split_equation(s):
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
    if f.operator() is add_vararg:
        return list(f.operands())
    if f.is_zero():
        return []
    return [f]


def _factor_list(t):
    if t.operator() is mul_vararg:
        return list(t.operands())
    return [t]


def _as_QQ(ex):
    try:
        return QQ(ex)
    except (TypeError, ValueError):
        return None


def _as_ZZ(ex):
    try:
        return ZZ(ex)
    except (TypeError, ValueError):
        return None


def _fractional_structure(f0, params):
    """Detect a 'sum of unit fractions = a/n' structure in the uncleared input.

    Returns a dict with keys ``unit`` (list of (sign, varname)), ``target``
    ((a, n) with n an int or parameter name) or None if the shape doesn't match.
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
    """Handle exponents like 3*n by absorbing constants into the base."""
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
                            f"unexpected exponent {e} in {fac} after clearing "
                            "denominators"
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
            raise UnsupportedEquationError(f"unsupported factor {fac} in term {t}")
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
    firsts = {}
    for m in _IDENT.finditer(s):
        name = m.group(0)
        if name in names and name not in firsts:
            firsts[name] = m.start()
    return tuple(sorted(names, key=lambda n: firsts.get(n, len(s))))


def parse(equation, params=(), domain="ZZ"):
    """Parse an equation string into a :class:`ParsedEquation`.

    ``params`` names symbols to treat as parameters (coefficients) rather than
    unknowns; ``domain`` records where solutions are sought (``"ZZ"``, ``"NN"``
    or ``"QQ"``) and is threaded through to classification output.
    """
    if isinstance(params, str):
        params = tuple(p.strip() for p in re.split(r"[,\s]+", params) if p.strip())
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
            "multiplied by %s (valid where %s nonzero)" % (den, ", ".join(denvars))
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
        original=str(equation), lhs=lhs, rhs=rhs, params=params, domain=domain,
        terms=terms, unknowns=unknowns, exp_unknowns=exp_unknowns,
        conditions=conditions, fractional=fractional,
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
            if params:
                c *= prod(pgens[p] ** e for p, e in t.param_powers) if t.param_powers else base.one()
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
