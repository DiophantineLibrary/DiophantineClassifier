r"""
Coordinate maps between a user's equation and the form a matcher recognized.

A matcher rarely finds an equation in exactly the standard form of its
family.  ``5*x^2 - y^2 = 1`` is a Pell equation, but only after reading the
user's ``y`` as the Pell ``x``; ``x^3 + y^3 = z^3`` becomes a generalized
Fermat equation of signature `(3, 3, 3)` only after substituting
``z -> -z``.  Recording those normalizations as prose is enough to explain
what happened and useless for anything else: a solver that works in the
standard coordinates has to carry its answer *back*, and a string cannot
transport a solution.

:class:`CoordinateTransform` is that map, in both directions and executable.
It also carries the structural :attr:`~CoordinateTransform.roles` a family
assigns to variables, so solvers can build an assignment in canonical names
and pull it back without reaching into a match's ``data`` dict.

Operations that change the *equation* but not the coordinates — multiplying
through by `-1`, say — are not coordinate changes and live in
:attr:`~CoordinateTransform.operations` instead.

EXAMPLES::

    sage: from diophantine_classifier.transforms import rename, identity
    sage: swap = rename([("x", "b"), ("y", "a")], description="swapped variables")
    sage: swap.push_forward({"a": 2, "b": 3})       # into standard coordinates
    {'x': 3, 'y': 2}
    sage: swap.pull_back({"x": 3, "y": 2})          # and back again
    {'a': 2, 'b': 3}
    sage: identity(("x", "y")).is_identity
    True
"""

from dataclasses import dataclass, field

from sage.all import QQ, SR, ZZ
from sage.calculus.calculus import symbolic_expression_from_string


def _identifiers(text):
    r"""
    The identifiers occurring in an expression string.

    EXAMPLES::

        sage: from diophantine_classifier.transforms import _identifiers
        sage: sorted(_identifiers("-x + y_1"))
        ['x', 'y_1']
    """
    import re
    return set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text))


def _substitute(text, mapping):
    r"""
    Rewrite an expression by replacing its identifiers with expressions.

    INPUT:

    - ``text`` -- string; an expression
    - ``mapping`` -- dict mapping identifier names to expression strings

    OUTPUT: string; the rewritten expression

    EXAMPLES::

        sage: from diophantine_classifier.transforms import _substitute
        sage: _substitute("x", {"x": "-v"})
        '-v'
        sage: _substitute("-z", {"z": "-w"})
        'w'
        sage: _substitute("x + 1", {"y": "u"})            # nothing to do
        'x + 1'
    """
    names = _identifiers(text)
    replace = {name: str(mapping[name]) for name in names if name in mapping}
    if not replace:
        return text
    table = {name: SR.var(name) for name in names}
    expr = symbolic_expression_from_string(text, table)
    subs = {}
    for name, repl in replace.items():
        rnames = _identifiers(repl)
        subs[SR.var(name)] = symbolic_expression_from_string(
            repl, {other: SR.var(other) for other in rnames})
    return str(expr.subs(subs))


def _evaluate(text, assignment):
    r"""
    Evaluate one coordinate expression on an assignment, exactly.

    Identifiers resolve against ``assignment`` only, so Sage's global names
    (``e``, ``I``, ``pi``, ...) cannot leak into a coordinate map.  Integer
    and rational results come back as Sage numbers rather than as symbolic
    expressions, so transported solutions compare and print like the ones a
    solver produced directly.

    INPUT:

    - ``text`` -- string; an expression in the other coordinate system
    - ``assignment`` -- dict mapping variable names to values

    OUTPUT: the value of the expression

    EXAMPLES::

        sage: from diophantine_classifier.transforms import _evaluate
        sage: _evaluate("-x", {"x": 5})
        -5
        sage: _evaluate("x", {"x": QQ(1)/2})
        1/2
        sage: _evaluate("x + y", {"x": 1})
        Traceback (most recent call last):
        ...
        KeyError: "coordinate expression 'x + y' needs a value for 'y'"
    """
    names = _identifiers(text)
    missing = sorted(name for name in names if name not in assignment)
    if missing:
        raise KeyError(f"coordinate expression {text!r} needs a value for "
                       f"{missing[0]!r}")
    table = {name: SR.var(name) for name in names}
    value = symbolic_expression_from_string(text, table)
    if names:
        value = value.subs({SR.var(name): assignment[name]
                            for name in names})
    for ring in (ZZ, QQ):
        try:
            return ring(value)
        except (TypeError, ValueError):
            continue
    return value


def _rewrite_role(value, mapping):
    r"""
    Rewrite a role's variable (or list of variables) through a mapping.

    INPUT:

    - ``value`` -- a variable name, or a list/tuple of them
    - ``mapping`` -- dict from those names to expressions

    OUTPUT: the rewritten name, or list of names

    EXAMPLES::

        sage: from diophantine_classifier.transforms import _rewrite_role
        sage: _rewrite_role("v", {"v": "b"})
        'b'
        sage: _rewrite_role(["u", "v"], {"u": "a", "v": "b"})
        ['a', 'b']
    """
    if isinstance(value, (list, tuple)):
        return [_substitute(str(item), mapping) for item in value]
    return _substitute(str(value), mapping)


@dataclass(frozen=True)
class CoordinateTransform:
    r"""
    An invertible change of variables between source and standard form.

    Both directions are stored explicitly, each keyed by *its own* system's
    variables and written in terms of the other's: ``to_normalized["x"]`` is
    the formula for the standard coordinate ``x`` in the user's variables,
    and ``to_source["a"]`` is the formula for the user's ``a`` in standard
    coordinates.  Expression maps cover the swaps and sign changes the
    matchers perform today and leave room for affine or birational maps
    later.

    ATTRIBUTES:

    - ``description`` -- string; what the normalization did, for display.
    - ``source_variables`` -- tuple of strings; the unknowns as the user
      wrote them.
    - ``normalized_variables`` -- tuple of strings; the family's standard
      coordinate names.
    - ``to_normalized`` -- dict mapping each standard coordinate to an
      expression in the source variables.
    - ``to_source`` -- dict mapping each source variable to an expression in
      the standard coordinates.
    - ``conditions`` -- tuple of
      :class:`~diophantine_classifier.conditions.NonzeroCondition`; where
      the map is defined, when that is not everywhere.
    - ``roles`` -- dict mapping a structural role name to the source
      variable (or list of source variables) playing it, e.g.
      ``{"legs": ["x", "y"], "hypotenuse": "z"}``.  Solver-relevant, hence
      here rather than buried in a match's ``data``.
    - ``operations`` -- tuple of strings; normalizations that changed the
      equation but no coordinate, such as multiplying through by `-1`.

    EXAMPLES::

        sage: from diophantine_classifier.transforms import CoordinateTransform
        sage: flip = CoordinateTransform(
        ....:     description="substituted z -> -z (odd exponent)",
        ....:     source_variables=("z",), normalized_variables=("z",),
        ....:     to_normalized={"z": "-z"}, to_source={"z": "-z"})
        sage: flip.push_forward({"z": 4})
        {'z': -4}
        sage: flip.pull_back(flip.push_forward({"z": 4}))
        {'z': 4}
    """
    description: str = ""
    source_variables: tuple = ()
    normalized_variables: tuple = ()
    to_normalized: dict = field(default_factory=dict)
    to_source: dict = field(default_factory=dict)
    conditions: tuple = ()
    roles: dict = field(default_factory=dict)
    operations: tuple = ()

    @property
    def is_identity(self):
        r"""
        Whether the map leaves every coordinate alone.

        Equation-level :attr:`operations` do not count: they change no
        coordinate, so a sign-flipped equation still has an identity map.

        EXAMPLES::

            sage: from diophantine_classifier.transforms import identity, rename
            sage: identity(("x", "y"), operations=("multiplied by -1",)).is_identity
            True
            sage: rename([("x", "y"), ("y", "x")]).is_identity
            False
        """
        return all(name == expr for name, expr in self.to_normalized.items())

    def push_forward(self, assignment):
        r"""
        Carry an assignment in the user's variables to standard coordinates.

        INPUT:

        - ``assignment`` -- dict mapping source variable names to values

        OUTPUT: dict mapping standard coordinate names to values

        EXAMPLES::

            sage: from diophantine_classifier.transforms import rename
            sage: t = rename([("x", "b"), ("y", "a")])
            sage: t.push_forward({"a": 2, "b": 3})
            {'x': 3, 'y': 2}
        """
        return {name: _evaluate(expr, assignment)
                for name, expr in self.to_normalized.items()}

    def pull_back(self, assignment):
        r"""
        Carry an assignment in standard coordinates back to the user's.

        This is what turns a solver's answer into a solution of the problem
        that was submitted.

        INPUT:

        - ``assignment`` -- dict mapping standard coordinate names to values

        OUTPUT: dict mapping source variable names to values

        EXAMPLES::

            sage: from diophantine_classifier.transforms import rename
            sage: t = rename([("x", "b"), ("y", "a")])
            sage: t.pull_back({"x": 3, "y": 2})
            {'a': 2, 'b': 3}
        """
        return {name: _evaluate(expr, assignment)
                for name, expr in self.to_source.items()}

    def conditions_hold(self, assignment):
        r"""
        Decide this map's own side conditions on a *source* assignment.

        :attr:`conditions` are stored in the transform's source coordinates
        (see :meth:`then`), so they are checked against an assignment in the
        user's variables — the same one the parser's conditions are checked
        against.

        INPUT:

        - ``assignment`` -- dict mapping source variable names to values

        OUTPUT: ``True``, ``False``, or ``None`` when undetermined

        EXAMPLES::

            sage: from diophantine_classifier.conditions import NonzeroCondition
            sage: from diophantine_classifier.transforms import rename
            sage: t = rename([("x", "u"), ("y", "v")],
            ....:            conditions=(NonzeroCondition("u", ("u",)),))
            sage: t.conditions_hold({"u": 1, "v": 0})
            True
            sage: t.conditions_hold({"u": 0, "v": 1})
            False
        """
        from .conditions import hold
        return hold(self.conditions, assignment)

    def then(self, after):
        r"""
        Compose two maps: ``self`` from `A` to `B`, ``after`` from `B` to
        `C`, giving the map from `A` to `C`.

        The composite satisfies, for every assignment ``a`` in `A` and ``c``
        in `C`::

            self.then(after).push_forward(a)
                == after.push_forward(self.push_forward(a))
            self.then(after).pull_back(c)
                == self.pull_back(after.pull_back(c))

        Conditions live in each map's **source** coordinates, so ``after``'s
        are pulled back through ``self`` and the composite's conditions are
        all checkable on an assignment in `A`.  Equation-level
        :attr:`operations` are concatenated in the order they happened, and
        roles naming `B`-variables are rewritten to their `A`-expressions.

        INPUT:

        - ``after`` -- a :class:`CoordinateTransform` whose source
          coordinates are this one's normalized coordinates

        OUTPUT: a :class:`CoordinateTransform`

        EXAMPLES::

            sage: from diophantine_classifier.transforms import negate, rename
            sage: first = negate(("u", "v"), ("v",))       # A -> B
            sage: second = rename([("x", "v"), ("y", "u")])  # B -> C
            sage: both = first.then(second)
            sage: both.to_normalized == {"x": "-v", "y": "u"}
            True
            sage: both.push_forward({"u": 2, "v": 3})
            {'x': -3, 'y': 2}
            sage: both.pull_back(both.push_forward({"u": 2, "v": 3}))
            {'u': 2, 'v': 3}

        Composing with an identity changes nothing::

            sage: from diophantine_classifier.transforms import identity
            sage: identity(("u", "v")).then(second).to_normalized
            {'x': 'v', 'y': 'u'}
        """
        parts = [text for text in (self.description, after.description)
                 if text]
        return CoordinateTransform(
            description="; ".join(parts),
            source_variables=self.source_variables,
            normalized_variables=after.normalized_variables,
            # a C-coordinate is written in B; rewrite each B-name in A
            to_normalized={name: _substitute(expr, self.to_normalized)
                           for name, expr in after.to_normalized.items()},
            # an A-variable is written in B; rewrite each B-name in C
            to_source={name: _substitute(expr, after.to_source)
                       for name, expr in self.to_source.items()},
            conditions=(tuple(self.conditions)
                        + tuple(c.substitute(self.to_normalized)
                                for c in after.conditions)),
            roles={**{key: _rewrite_role(value, self.to_normalized)
                      for key, value in after.roles.items()},
                   **dict(self.roles)},
            operations=tuple(self.operations) + tuple(after.operations),
        )

    def as_dict(self):
        r"""
        JSON-serializable form (the website-backend contract).

        OUTPUT: dict with plain types only

        EXAMPLES::

            sage: from diophantine_classifier.transforms import rename
            sage: import json
            sage: d = rename([("x", "y"), ("y", "x")],
            ....:            description="swapped variables").as_dict()
            sage: d["to_normalized"], d["to_source"]
            ({'x': 'y', 'y': 'x'}, {'x': 'y', 'y': 'x'})
            sage: _ = json.dumps(d)
        """
        out = {
            "description": self.description,
            "identity": self.is_identity,
            "source_variables": list(self.source_variables),
            "normalized_variables": list(self.normalized_variables),
            "to_normalized": {str(k): str(v)
                              for k, v in self.to_normalized.items()},
            "to_source": {str(k): str(v) for k, v in self.to_source.items()},
        }
        if self.conditions:
            out["conditions"] = [c.as_dict() for c in self.conditions]
        if self.roles:
            out["roles"] = {str(k): (list(v) if isinstance(v, (list, tuple))
                                     else str(v))
                            for k, v in self.roles.items()}
        if self.operations:
            out["operations"] = list(self.operations)
        return out

    def __str__(self):
        r"""
        One-line rendering: the coordinate change and the equation-level
        operations, or the empty string when neither happened.

        EXAMPLES::

            sage: from diophantine_classifier.transforms import identity, rename
            sage: str(rename([("x", "y"), ("y", "x")],
            ....:            description="swapped variables",
            ....:            operations=("multiplied by -1",)))
            'multiplied by -1; swapped variables'
            sage: str(identity(("x",)))
            ''
        """
        parts = list(self.operations)
        if self.description:
            parts.append(self.description)
        return "; ".join(parts)

    def __bool__(self):
        r"""
        Whether anything happened at all (used for "was there a transform?").

        EXAMPLES::

            sage: from diophantine_classifier.transforms import identity, rename
            sage: bool(identity(("x", "y")))
            False
            sage: bool(rename([("x", "y"), ("y", "x")]))
            True
        """
        return bool(self.operations) or not self.is_identity


def identity(variables, description="", operations=(), roles=None,
             normalized_variables=None, conditions=()):
    r"""
    The transform that renames nothing.

    The default for a match found already in standard form — including one
    reached by an equation-level operation, which changes no coordinate.

    INPUT:

    - ``variables`` -- iterable of source variable names
    - ``description`` -- (default: ``""``) display text
    - ``operations`` -- (default: ``()``) equation-level normalizations
    - ``roles`` -- (default: ``None``) structural role map
    - ``normalized_variables`` -- (default: ``None``) standard coordinate
      names, when they should be recorded as something other than
      ``variables``
    - ``conditions`` -- (default: ``()``) where the map is defined

    OUTPUT: a :class:`CoordinateTransform`

    EXAMPLES::

        sage: from diophantine_classifier.transforms import identity
        sage: t = identity(("x", "y"), operations=("multiplied by -1",))
        sage: t.push_forward({"x": 1, "y": 2})
        {'x': 1, 'y': 2}
        sage: str(t)
        'multiplied by -1'
    """
    variables = tuple(str(v) for v in variables)
    same = {v: v for v in variables}
    return CoordinateTransform(
        description=description,
        source_variables=variables,
        normalized_variables=tuple(normalized_variables)
        if normalized_variables is not None else variables,
        to_normalized=dict(same), to_source=dict(same),
        conditions=tuple(conditions),
        roles=dict(roles or {}), operations=tuple(operations),
    )


def rename(pairs, description="", operations=(), roles=None, conditions=()):
    r"""
    A transform that only renames or permutes coordinates.

    Covers every variable-role reversal the matchers perform: reading the
    user's second variable as the standard ``x``, calling the base of a
    curve ``x`` and the squared variable ``y``, and so on.

    INPUT:

    - ``pairs`` -- iterable of ``(normalized name, source name)``
    - ``description`` -- (default: ``""``) display text
    - ``operations`` -- (default: ``()``) equation-level normalizations
    - ``roles`` -- (default: ``None``) structural role map; defaults to the
      pairing itself
    - ``conditions`` -- (default: ``()``) where the map is defined

    OUTPUT: a :class:`CoordinateTransform`

    EXAMPLES::

        sage: from diophantine_classifier.transforms import rename
        sage: t = rename([("x", "u"), ("y", "v")])
        sage: t.push_forward({"u": 7, "v": 8})
        {'x': 7, 'y': 8}
        sage: t.pull_back({"x": 7, "y": 8})
        {'u': 7, 'v': 8}
        sage: t.roles
        {'x': 'u', 'y': 'v'}
    """
    pairs = [(str(n), str(s)) for n, s in pairs]
    return CoordinateTransform(
        description=description,
        source_variables=tuple(s for _, s in pairs),
        normalized_variables=tuple(n for n, _ in pairs),
        to_normalized={n: s for n, s in pairs},
        to_source={s: n for n, s in pairs},
        conditions=tuple(conditions),
        roles=dict(roles) if roles is not None else {n: s for n, s in pairs},
        operations=tuple(operations),
    )


def negate(variables, flipped, description="", operations=(), roles=None,
           conditions=()):
    r"""
    A transform substituting ``v -> -v`` for the named variables.

    Used to orient diagonal equations: at an odd exponent, flipping the sign
    of a variable moves its term to the other side of the equation, which is
    how ``x^3 + y^3 + z^3 = 0`` is read as a generalized Fermat equation.
    The map is its own inverse.

    INPUT:

    - ``variables`` -- iterable of all source variable names
    - ``flipped`` -- iterable of the names to negate
    - ``description`` -- (default: ``""``) display text
    - ``operations`` -- (default: ``()``) equation-level normalizations
    - ``roles`` -- (default: ``None``) structural role map
    - ``conditions`` -- (default: ``()``) where the map is defined

    OUTPUT: a :class:`CoordinateTransform`

    EXAMPLES::

        sage: from diophantine_classifier.transforms import negate
        sage: t = negate(("x", "y", "z"), ("z",))
        sage: t.push_forward({"x": 1, "y": 2, "z": 3})
        {'x': 1, 'y': 2, 'z': -3}
        sage: t.pull_back(t.push_forward({"x": 1, "y": 2, "z": 3}))
        {'x': 1, 'y': 2, 'z': 3}
    """
    variables = tuple(str(v) for v in variables)
    flipped = {str(v) for v in flipped}
    mapping = {v: (f"-{v}" if v in flipped else v) for v in variables}
    return CoordinateTransform(
        description=description,
        source_variables=variables, normalized_variables=variables,
        to_normalized=dict(mapping), to_source=dict(mapping),
        conditions=tuple(conditions),
        roles=dict(roles or {}), operations=tuple(operations),
    )
