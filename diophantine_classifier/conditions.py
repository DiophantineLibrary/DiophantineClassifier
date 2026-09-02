r"""
Side conditions attached to a parsed equation.

Clearing denominators turns one problem into another: ``1/(x-1) = 1/(y-1)``
becomes ``y = x`` only *where* ``x - 1`` and ``y - 1`` are nonzero, and
``(1, 1)`` solves the cleared equation without solving the original one.  The
hypothesis introduced by that step is recorded as a structured, executable
predicate rather than as prose, so that solvers can enforce it and the
website backend can serialize it.

Conditions are collected from the *source syntax*, before Sage simplifies
anything: ``x/x = y`` normalizes to ``y = 1``, but ``x != 0`` still belongs to
the problem (and ``x`` is still one of its unknowns).

Evaluation is exact and three-valued: an assignment that leaves a symbol
undetermined yields ``None`` rather than a guess.

EXAMPLES::

    sage: from diophantine_classifier.conditions import NonzeroCondition
    sage: c = NonzeroCondition("x - 1", ("x",))
    sage: str(c)
    'x - 1 != 0'
    sage: c.evaluate({"x": 2}), c.evaluate({"x": 1}), c.evaluate({}) is None
    (True, False, True)
    sage: c.as_dict()
    {'expression': 'x - 1', 'source': 'denominator', 'type': 'nonzero', 'variables': ['x']}
"""

from dataclasses import dataclass, field

from sage.all import SR
from sage.calculus.calculus import symbolic_expression_from_string


def _expression(text, variables):
    r"""
    Build the Sage expression for a condition, with a private symbol table.

    Identifiers are looked up in a table built from ``variables`` only, so
    Sage's global names (``e``, ``I``, ``pi``, ...) cannot leak in.

    INPUT:

    - ``text`` -- string; the condition's expression source
    - ``variables`` -- iterable of strings; the identifiers occurring in it

    OUTPUT: a symbolic expression

    EXAMPLES::

        sage: from diophantine_classifier.conditions import _expression
        sage: _expression("x - 1", ("x",))
        x - 1
        sage: _expression("e", ("e",)).variables()      # not Euler's number
        (e,)
    """
    table = {name: SR.var(name) for name in variables}
    return symbolic_expression_from_string(text, table)


@dataclass(frozen=True)
class NonzeroCondition:
    r"""
    The hypothesis that an expression does not vanish.

    ATTRIBUTES:

    - ``expression`` -- string; the expression that must not vanish, as it
      appeared (up to canonicalization) in the source equation, e.g.
      ``"x - 1"``.
    - ``variables`` -- tuple of strings; the identifiers it involves, in
      sorted order.
    - ``source`` -- string; why the condition exists (``"denominator"`` for
      cleared denominators).

    EXAMPLES::

        sage: from diophantine_classifier.conditions import NonzeroCondition
        sage: c = NonzeroCondition("x*y", ("x", "y"))
        sage: c
        NonzeroCondition('x*y != 0')
        sage: c.evaluate({"x": 3, "y": 0})
        False
    """
    expression: str
    variables: tuple = ()
    source: str = "denominator"
    _cache: dict = field(default_factory=dict, repr=False, compare=False,
                         hash=False)

    def __str__(self):
        r"""
        Render the condition as it is displayed to users.

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: str(NonzeroCondition("x - 1", ("x",)))
            'x - 1 != 0'
        """
        return f"{self.expression} != 0"

    def __repr__(self):
        r"""
        Terse representation.

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: NonzeroCondition("x", ("x",))
            NonzeroCondition('x != 0')
        """
        return f"NonzeroCondition({str(self)!r})"

    def expression_sr(self):
        r"""
        The condition's expression as a Sage expression.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: NonzeroCondition("x - 1", ("x",)).expression_sr()
            x - 1
        """
        if "sr" not in self._cache:
            self._cache["sr"] = _expression(self.expression, self.variables)
        return self._cache["sr"]

    def evaluate(self, assignment):
        r"""
        Decide the condition on an assignment, exactly and three-valued.

        INPUT:

        - ``assignment`` -- dict mapping variable names to values

        OUTPUT: ``True`` (holds), ``False`` (violated), or ``None`` when the
        assignment does not determine the value

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: c = NonzeroCondition("x*y - 2", ("x", "y"))
            sage: c.evaluate({"x": 2, "y": 3})
            True
            sage: c.evaluate({"x": 1, "y": 2})
            False
            sage: c.evaluate({"x": 1}) is None
            True
        """
        if any(name not in assignment for name in self.variables):
            return None
        value = self.expression_sr().subs(
            {SR.var(name): assignment[name] for name in self.variables})
        try:
            return not bool(value == 0)
        except TypeError:
            return None

    def substitute(self, renaming):
        r"""
        Rewrite the condition through a renaming of variables.

        INPUT:

        - ``renaming`` -- dict mapping variable names to expression strings

        OUTPUT: a new :class:`NonzeroCondition`

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: NonzeroCondition("x - 1", ("x",)).substitute({"x": "y"})
            NonzeroCondition('y - 1 != 0')
        """
        names = set(self.variables) | {
            str(v) for text in renaming.values()
            for v in _expression(text, _identifiers(text)).variables()}
        table = {name: SR.var(name) for name in names}
        image = self.expression_sr().subs(
            {SR.var(k): _expression(v, _identifiers(v))
             for k, v in renaming.items() if k in self.variables})
        del table
        return NonzeroCondition(str(image),
                                tuple(sorted(str(v) for v in image.variables())),
                                self.source)

    def as_dict(self):
        r"""
        JSON-serializable form (the website-backend contract).

        OUTPUT: dict with plain types only

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: import json
            sage: d = NonzeroCondition("x - 1", ("x",)).as_dict()
            sage: d["type"], d["expression"]
            ('nonzero', 'x - 1')
            sage: _ = json.dumps(d)
        """
        return {"type": "nonzero", "expression": self.expression,
                "variables": list(self.variables), "source": self.source}


def _identifiers(text):
    r"""
    The identifiers occurring in an expression string.

    EXAMPLES::

        sage: from diophantine_classifier.conditions import _identifiers
        sage: sorted(_identifiers("x*y - 2"))
        ['x', 'y']
    """
    import re
    return set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text))


def hold(conditions, assignment):
    r"""
    Decide a collection of conditions on one assignment.

    INPUT:

    - ``conditions`` -- iterable of conditions
    - ``assignment`` -- dict mapping variable names to values

    OUTPUT: ``True`` if all hold, ``False`` if one is violated, ``None`` if
    none is violated but some are undetermined

    EXAMPLES::

        sage: from diophantine_classifier.conditions import NonzeroCondition, hold
        sage: cs = [NonzeroCondition("x", ("x",)), NonzeroCondition("y", ("y",))]
        sage: hold(cs, {"x": 1, "y": 2})
        True
        sage: hold(cs, {"x": 0, "y": 2})
        False
        sage: hold(cs, {"x": 1}) is None
        True
        sage: hold([], {})
        True
    """
    verdict = True
    for condition in conditions:
        value = condition.evaluate(assignment)
        if value is False:
            return False
        if value is None:
            verdict = None
    return verdict
