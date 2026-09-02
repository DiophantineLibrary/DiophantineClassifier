r"""
Classification pipeline: parse, split reducible equations, run matchers,
rank by specificity in the family DAG, and package the result.

EXAMPLES::

    sage: from diophantine_classifier import classify
    sage: classify("3*x + 5*y = 1")
    Classification('3*x + 5*y = 1' -> linear)
    sage: classify("x^2 - 5*x + 6 = 0").slug
    'univariate'
    sage: cls = classify("(x^2 - 2)*(y^2 - 3) = 0")   # reducible: components
    sage: cls.slug, [c.slug for c in cls.components]
    ('reducible', ['univariate', 'univariate'])
"""

from dataclasses import dataclass, field

from sage.all import ZZ

from . import matchers
from .parsing import ParsedEquation, parse
from .registry import (ancestors, depth, families, lineage_graph,
                       lineage_paths)
from .transforms import identity


def _jsonify(x):
    r"""
    Recursively convert match data to plain JSON-serializable types.

    Sage integers become Python ints; anything else non-primitive becomes a
    string.  Keeps :meth:`Classification.as_dict` honest as the website
    backend contract.

    EXAMPLES::

        sage: from diophantine_classifier.classify import _jsonify
        sage: _jsonify({"D": ZZ(61), "roles": {"x": "x"}})
        {'D': 61, 'roles': {'x': 'x'}}
        sage: _jsonify([QQ(1)/2, True])
        ['1/2', True]
    """
    if isinstance(x, (str, bool, int, type(None))):
        return x
    if isinstance(x, dict):
        return {str(k): _jsonify(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonify(v) for v in x]
    try:
        return int(ZZ(x))
    except (TypeError, ValueError):
        return str(x)


@dataclass
class Classification:
    r"""
    The result of :func:`classify`.

    ATTRIBUTES:

    - ``parsed`` -- the problem *as submitted*: the
      :class:`~diophantine_classifier.parsing.ParsedEquation` built from the
      user's string.  Display, serialization and solution validation all go
      through this one.
    - ``model`` -- the working
      :class:`~diophantine_classifier.parsing.ParsedEquation` the matchers
      actually saw, when normalization produced a different equation (a
      repeated factor reduced to its base).  ``None`` when it is just
      ``parsed``; read it through :attr:`working`.
    - ``matches`` -- list of :class:`~diophantine_classifier.matchers.Match`,
      sorted by decreasing specificity (DAG depth, ties by emission order).
      The first entry is the primary match.
    - ``components`` -- list of :class:`Classification`; nonempty exactly
      when the equation's polynomial factors, in which case the solution set
      is the union over the components and there is no primary match.
    - ``reduction`` -- the
      :class:`~diophantine_classifier.transforms.CoordinateTransform` from
      the source problem to :attr:`working`, or ``None`` when they coincide.
      Reducing `g^e = 0` to `g = 0` changes no coordinate, so that one is
      the identity — but it is recorded, because a solver has to know which
      problem its answer belongs to.
    - ``note`` -- string; pipeline remark (e.g. that a repeated factor was
      reduced).

    Derived properties: :attr:`primary`, :attr:`slug`, :attr:`family`,
    :attr:`lineage`, :attr:`lineage_paths`, :attr:`data`, :attr:`working`.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: cls = classify("3*x + 5*y = 1")
        sage: cls.slug
        'linear'
        sage: cls.lineage
        ['general-polynomial']
        sage: cls.data["b"]
        '1'

    A reduction keeps the problem that was asked::

        sage: cls = classify("(x + y)^2 = 0", domain="QQ")
        sage: cls.slug, cls.parsed.original, cls.parsed.domain
        ('linear', '(x + y)^2 = 0', 'QQ')
        sage: cls.working.poly
        x + y
    """
    parsed: ParsedEquation
    model: ParsedEquation = None
    matches: list = field(default_factory=list)
    components: list = field(default_factory=list)
    reduction: object = None
    embedding: object = None
    special_kind: str = ""
    note: str = ""

    # ------------------------------------------------------------------ api

    @property
    def free_variables(self):
        r"""
        Ambient unknowns this component leaves unconstrained.

        Empty unless the classification is a component of a reducible
        equation whose factor mentions fewer variables than the problem.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1").free_variables
            ()
            sage: cls = classify("(x - 2)*(x - y)/y = 0")
            sage: next(c for c in cls.components
            ....:      if str(c.working.poly) == "x - 2").free_variables
            ('y',)
        """
        return () if self.embedding is None else self.embedding.free_variables

    def effective_transform(self, match):
        r"""
        The whole map a solver's answer follows back to the submitted
        problem: source → working model → the family's standard coordinates.

        When a normalization replaced the problem by an equivalent one, that
        reduction is composed with the match's own transform, so a caller
        never has to know there were two steps — or in which order to undo
        them.

        INPUT:

        - ``match`` -- a :class:`~diophantine_classifier.matchers.Match` of
          this classification

        OUTPUT: a
        :class:`~diophantine_classifier.transforms.CoordinateTransform`

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: cls = classify("3*x + 5*y = 1")
            sage: cls.effective_transform(cls.primary) is cls.primary.transform
            True
            sage: cls = classify("(x + y)^2 = 0")
            sage: cls.effective_transform(cls.primary).source_variables
            ('x', 'y')
        """
        if self.reduction is None:
            return match.transform
        return self.reduction.then(match.transform)

    @property
    def working(self):
        r"""
        The equation the matchers and solvers work on.

        Equal to :attr:`parsed` unless a normalization replaced the problem
        by an equivalent smaller one.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: cls = classify("3*x + 5*y = 1")
            sage: cls.working is cls.parsed
            True
            sage: classify("(x + y)^2 = 0").working.original
            'x + y = 0'
        """
        return self.parsed if self.model is None else self.model

    @property
    def is_composite(self):
        r"""
        Whether the equation split into factor components.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("(x^2 - 2)*(y^2 - 3) = 0").is_composite
            True
            sage: classify("x^2 - 2*y^2 = 1").is_composite
            False
        """
        return bool(self.components)

    @property
    def primary(self):
        r"""
        The most specific match, or ``None`` for composite equations.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1").primary
            Match('linear')
        """
        return self.matches[0] if self.matches else None

    @property
    def slug(self):
        r"""
        The primary family's slug (``"reducible"`` for composites).

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("x^2 - 5*x + 6 = 0").slug
            'univariate'
        """
        if self.is_composite:
            return "reducible"
        return self.primary.slug if self.primary else "unclassified"

    @property
    def family(self):
        r"""
        The primary :class:`~diophantine_classifier.registry.Family`.

        ``None`` for composite equations.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1").family
            Family('linear')
        """
        return families().get(self.slug)

    @property
    def lineage(self):
        r"""
        Ancestors of the primary family, most specific first.

        A flattened *set*, kept for membership tests and for walking
        candidate ancestor solvers.  It is not a chain — see
        :attr:`lineage_paths` for ancestry that can be displayed.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1").lineage[0]
            'general-polynomial'
        """
        if self.primary is None:
            return []
        return ancestors(self.primary.slug)

    @property
    def lineage_paths(self):
        r"""
        The genuine specialization paths from the primary family to a root.

        One list per path, each starting at the primary family; consecutive
        entries are real DAG edges, so a path may be rendered as an arrow
        chain.  A family with two parents has two paths.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1").lineage_paths
            [['linear', 'general-polynomial']]
        """
        if self.primary is None:
            return []
        return lineage_paths(self.primary.slug)

    @property
    def data(self):
        r"""
        The primary match's extracted data.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1").data["b"]
            '1'
        """
        return self.primary.data if self.primary else {}

    def match_for(self, slug):
        r"""
        The emitted match for a family, or ``None``.

        INPUT:

        - ``slug`` -- string; a family slug

        OUTPUT: a :class:`~diophantine_classifier.matchers.Match` or ``None``

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: cls = classify("3*x + 5*y = 1")
            sage: cls.match_for("linear")
            Match('linear')
            sage: cls.match_for("pell") is None
            True
        """
        for match in self.matches:
            if match.slug == slug:
                return match
        return None

    def data_for(self, slug):
        r"""
        The data of the emitted match for a family, or ``None``.

        A family in the lineage that was never actually matched has no data
        here: a DAG edge is not evidence that the primary match's data meets
        the ancestor's input contract.

        INPUT:

        - ``slug`` -- string; a family slug

        OUTPUT: dict or ``None``

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: cls = classify("3*x + 5*y = 1")
            sage: cls.data_for("linear")["b"]
            '1'
            sage: cls.data_for("general-polynomial")["degree"]
            1
            sage: cls.data_for("pell") is None
            True
        """
        match = self.match_for(slug)
        return None if match is None else match.data

    def code(self):
        r"""
        Filled code templates, each from the data of its own match.

        Only families that were actually matched appear.  Filling an
        ancestor's template with the primary match's data used to produce
        confident nonsense — a Pell ``D`` handed to a quadratic-form
        template — because a specialization edge says nothing about whose
        data satisfies whose input contract.

        OUTPUT: dict mapping ``"<system> (<slug>)"`` to a code string

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: cls = classify("3*x + 5*y = 1")
            sage: sorted(cls.code())               # only matched families
            ['sage (linear)']
            sage: "matrix(ZZ, [['3', '5']])" in cls.code()["sage (linear)"]
            True
        """
        out = {}
        for match in self.matches:
            fam = families().get(match.slug)
            if fam and fam.code:
                for lang, snippet in fam.fill_code(match.data).items():
                    out.setdefault(f"{lang} ({match.slug})", snippet)
        return out

    def as_dict(self):
        r"""
        JSON-serializable summary — the website-backend contract.

        OUTPUT: dict with plain types only (tested); keys include
        ``equation``, ``family``, ``status``, ``data``, ``lineage``,
        ``references`` (list of ``{key, why, formatted}`` dicts), ``code``,
        and ``components`` for composites

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: d = classify("3*x + 5*y = 1").as_dict()
            sage: d["family"], d["priority"]
            ('linear', 1)
            sage: d["references"][0]["key"]
            'NivenZuckermanMontgomery1991'
            sage: import json
            sage: _ = json.dumps(d)          # round-trips
        """
        d = {
            "equation": self.parsed.original,
            "domain": self.parsed.domain,
            "params": list(self.parsed.params),
            "unknowns": list(self.parsed.unknowns),
            "conditions": [c.as_dict() for c in self.parsed.conditions],
            "family": self.slug,
        }
        if self.note:
            d["note"] = self.note
        if self.special_kind:
            d["special_kind"] = self.special_kind
        if self.embedding is not None:
            d["component"] = self.embedding.as_dict()
        if self.reduction is not None:
            # provenance: the two halves of the effective map, separately
            d["reduction"] = self.reduction.as_dict()
        if self.is_composite:
            d["components"] = [c.as_dict() for c in self.components]
            return d
        fam = self.family
        if fam:
            d.update({
                "name": fam.name,
                "priority": fam.priority,
                "status": fam.status,
                "software": dict(fam.software),
                "references": [
                    {"key": key, "why": why, "formatted": formatted}
                    for key, why, formatted in fam.formatted_references()],
            })
            if fam.lmfdb:
                d["lmfdb"] = fam.lmfdb
            if fam.constraints:
                d["constraints"] = fam.constraints
        if self.primary:
            d["summary"] = self.primary.summary
            d["data"] = _jsonify(self.primary.data)
            # the complete map back to the submitted equation -- reduction
            # and match transform composed -- not just the matcher's half
            d["transform"] = self.effective_transform(self.primary).as_dict()
        d["lineage"] = self.lineage
        if self.primary is not None:
            # genuine ancestry: paths whose consecutive entries are real
            # edges, plus the subgraph for consumers that draw the DAG
            d["lineage_paths"] = self.lineage_paths
            d["lineage_graph"] = lineage_graph(self.primary.slug)
        d["all_matches"] = [m.slug for m in self.matches]
        d["code"] = self.code()
        return d

    # ------------------------------------------------------------- display

    def explain(self):
        r"""
        Multi-line human-readable report (equation-homepage prototype).

        OUTPUT: string

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: text = classify("3*x + 5*y = 1").explain()
            sage: "family: linear" in text
            True
            sage: "reference [NivenZuckermanMontgomery1991]" in text
            True
        """
        lines = [self.parsed.original]
        if self.parsed.params:
            lines.append(f"  parameters: {', '.join(self.parsed.params)}")
        if self.note:
            lines.append(f"  note: {self.note}")
        if self.parsed.cleared_by:
            lines.append("  cleared denominators: multiplied through by "
                         f"{self.parsed.cleared_by}")
        for cond in self.parsed.conditions:
            lines.append(f"  condition: {cond}")
        if self.free_variables:
            lines.append("  free in this component: "
                         + ", ".join(self.free_variables)
                         + " (unconstrained by the factor, still subject to "
                           "the conditions above)")
        if self.special_kind == "conditional-identity":
            lines.append("  conditional identity: every assignment "
                         "satisfying the conditions above solves the "
                         "equality; there is no equation left to match")
            return "\n".join(lines)
        if self.is_composite:
            lines.append("  reducible: the solution set is the union over "
                         "the factors")
            for i, comp in enumerate(self.components, 1):
                sub = comp.explain().splitlines()
                lines.append(f"  component {i}: {sub[0]}")
                lines.extend("  " + s for s in sub[1:])
            return "\n".join(lines)
        fam = self.family
        if fam is None:
            lines.append("  family: unclassified")
            return "\n".join(lines)
        lines.append(f"  family: {self.slug} — {fam.name}   "
                     f"[P{fam.priority}, {fam.status}]")
        if self.primary.summary:
            lines.append(f"  summary: {self.primary.summary}")
        if self.primary.transform:
            lines.append(f"  transform: {self.primary.transform}")
        if self.primary.data:
            lines.append("  data: " + ", ".join(
                f"{k}={v}" for k, v in self.primary.data.items()))
        # one line per genuine path: a flat ancestor list rendered as a
        # single chain would draw edges the DAG does not have
        paths = [path for path in self.lineage_paths if len(path) > 1]
        for path in paths:
            lines.append("  lineage: " + " → ".join(path))
        others = [m.slug for m in self.matches[1:]
                  if m.slug not in self.lineage]
        if others:
            lines.append("  also matches: " + ", ".join(others))
        if fam.notes:
            lines.append(f"  notes: {fam.notes}")
        if fam.methods:
            lines.append("  methods: " + "; ".join(fam.methods))
        for lang, what in fam.software.items():
            lines.append(f"  software[{lang}]: {what}")
        code = self.code()
        for lang, snippet in code.items():
            snippet = snippet.strip()
            if "\n" in snippet:
                lines.append(f"  code[{lang}]:")
                lines.extend("    " + s for s in snippet.splitlines())
            else:
                lines.append(f"  code[{lang}]: {snippet}")
        for key, why, formatted in fam.formatted_references():
            lines.append(f"  reference [{key}]: {formatted}")
            if why:
                lines.append(f"      relevance: {why}")
        return "\n".join(lines)

    def __repr__(self):
        r"""
        Terse representation.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("3*x + 5*y = 1")
            Classification('3*x + 5*y = 1' -> linear)
        """
        if self.is_composite:
            inner = ", ".join(c.slug for c in self.components)
            return (f"Classification({self.parsed.original!r} -> reducible: "
                    f"[{inner}])")
        return f"Classification({self.parsed.original!r} -> {self.slug})"


def _rank(match_list):
    r"""
    Sort matches by decreasing DAG depth, ties by emission order.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import run
        sage: from diophantine_classifier.classify import _rank
        sage: [m.slug for m in _rank(run(parse("3*x + 5*y = 1")))]
        ['linear', 'general-polynomial']
    """
    indexed = list(enumerate(match_list))
    indexed.sort(key=lambda pair: (-depth(pair[1].slug), pair[0]))
    return [m for _, m in indexed]


def _registered_matches(pe):
    r"""
    The matches for an equation, ranked, minus families the registry does
    not (yet) know.

    Every classification path goes through this: a recognizer whose family
    has not landed is inert, and :func:`~diophantine_classifier.registry.depth`
    would raise on its slug.

    INPUT:

    - ``pe`` -- a :class:`~diophantine_classifier.parsing.ParsedEquation`

    OUTPUT: list of :class:`~diophantine_classifier.matchers.Match`

    EXAMPLES::

        sage: from diophantine_classifier.classify import _registered_matches
        sage: from diophantine_classifier.parsing import parse
        sage: [m.slug for m in _registered_matches(parse("3*x + 5*y = 1"))]
        ['linear', 'general-polynomial']
    """
    known = families()
    return _rank([m for m in matchers.run(pe) if m.slug in known])


@dataclass(frozen=True)
class ComponentEmbedding:
    r"""
    How one factor of a reducible equation sits in the ambient problem.

    A factor of a polynomial in several variables cuts out a subset of the
    *same* space: the variables it does not mention are free on that
    component, and every side condition of the original problem still
    applies.  For ``(x - 2)*(x - y)/y = 0`` the component ``x - 2 = 0`` is
    "``x = 2``, with ``y`` arbitrary and nonzero" — not the one-variable
    point ``x = 2``.

    This is deliberately *not* a
    :class:`~diophantine_classifier.transforms.CoordinateTransform`:
    projecting an ambient cylinder onto its active factor is not a change of
    coordinates, and pretending it is invertible is how the free variables
    get lost.

    ATTRIBUTES:

    - ``source_variables`` -- the ambient unknowns, in the original order.
    - ``active_variables`` -- those the factor actually constrains.
    - ``free_variables`` -- the rest: unconstrained by this component, and
      still ranging over the requested domain.
    - ``factor`` -- string; the factor itself.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: cls = classify("(x - 2)*(x - y)/y = 0")
        sage: comp = next(c for c in cls.components
        ....:             if str(c.working.poly) == "x - 2")
        sage: comp.embedding.free_variables
        ('y',)
        sage: comp.embedding.active_variables
        ('x',)
    """
    source_variables: tuple = ()
    active_variables: tuple = ()
    free_variables: tuple = ()
    factor: str = ""

    def as_dict(self):
        r"""
        JSON-serializable form.

        OUTPUT: dict with plain types only

        EXAMPLES::

            sage: from diophantine_classifier.classify import ComponentEmbedding
            sage: ComponentEmbedding(("x", "y"), ("x",), ("y",),
            ....:                    "x - 2").as_dict()["free_variables"]
            ['y']
        """
        return {"factor": self.factor,
                "source_variables": list(self.source_variables),
                "active_variables": list(self.active_variables),
                "free_variables": list(self.free_variables)}


def _component(pe, g, reduction=None, note=""):
    r"""
    Classify one factor of a reducible equation *inside the ambient
    problem*.

    The result keeps ``pe`` as its source problem — the user's text, domain,
    parameters, unknown order and **all** its side conditions — and carries
    the factor as the working model.  Re-parsing the factor alone would drop
    the ambient variables and any condition mentioning them, which is a
    different (smaller) solution set.

    INPUT:

    - ``pe`` -- the reducible :class:`~diophantine_classifier.parsing.ParsedEquation`
    - ``g`` -- one of its polynomial factors
    - ``reduction`` -- (default: ``None``) the map from the source problem to
      the factor model, when the two are equivalent
    - ``note`` -- (default: ``""``) pipeline remark

    OUTPUT: a :class:`Classification`

    EXAMPLES::

        sage: from diophantine_classifier.classify import _component
        sage: from diophantine_classifier.parsing import parse
        sage: pe = parse("(x - 2)*(x - y)/y = 0")
        sage: comp = _component(pe, pe.poly.factor()[1][0])   # x - 2
        sage: comp.parsed.original, str(comp.working.poly)
        ('(x - 2)*(x - y)/y = 0', 'x - 2')

    The factor says nothing about ``y``, so ``y`` stays free — and the
    condition on it, which re-parsing the factor alone would have dropped,
    stays with the problem::

        sage: comp.free_variables
        ('y',)
        sage: sorted(str(c) for c in comp.parsed.conditions)
        ['y != 0']
    """
    model = parse(f"{g} = 0", params=pe.params, domain=pe.domain)
    active = tuple(v for v in pe.unknowns if v in set(model.unknowns))
    free = tuple(v for v in pe.unknowns if v not in set(model.unknowns))
    return Classification(
        parsed=pe, model=model, matches=_registered_matches(model),
        # a reduction is a claim of equivalence; with a free ambient
        # variable the model is a projection, and projections do not invert
        reduction=None if free else reduction,
        embedding=ComponentEmbedding(source_variables=pe.unknowns,
                                     active_variables=active,
                                     free_variables=free, factor=str(g)),
        note=note,
    )


def classify(equation, params=(), domain="ZZ"):
    r"""
    Classify a Diophantine equation into the most specific known family.

    INPUT:

    - ``equation`` -- string, or an already-parsed
      :class:`~diophantine_classifier.parsing.ParsedEquation`
    - ``params`` -- (default: ``()``) names of symbols to treat as
      parameters; iterable of strings or a comma/space-separated string
    - ``domain`` -- (default: ``"ZZ"``) one of ``"ZZ"``, ``"NN"``, ``"QQ"``

    OUTPUT: a :class:`Classification`

    Reducible polynomial equations (in more than one variable) split into
    components: the solution set is the union over the factors.  A repeated
    factor is reduced to the underlying one.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: classify("3*x + 5*y = 1").slug
        'linear'
        sage: classify("x^2 - 5*x + 6 = 0").slug
        'univariate'
        sage: classify("x^2 + y^3 + z^5 = 7").slug
        'general-polynomial'

    TESTS:

    Repeated factors are reduced, and the reduction keeps the problem that
    was submitted -- its text, its domain, its unknowns and its conditions::

        sage: cls = classify("(x + y)^2 = 0", domain="QQ")
        sage: cls.slug
        'linear'
        sage: cls.parsed.original, cls.parsed.domain
        ('(x + y)^2 = 0', 'QQ')
        sage: "reduced to the underlying factor" in cls.note
        True
    """
    if isinstance(equation, ParsedEquation):
        pe = equation
    else:
        pe = parse(equation, params=params, domain=domain)

    if pe.is_conditional_identity:
        # there is no polynomial here: the matchers would take a degree or a
        # factorization of an equation that does not exist
        return Classification(
            parsed=pe, special_kind="conditional-identity",
            note="the equality holds identically wherever its side "
                 "conditions do; the solution set is cut out by the "
                 "conditions alone",
        )

    # reducible polynomial = 0 splits into components (univariate equations
    # stay whole: factoring *is* their solution method)
    if pe.is_polynomial and pe.is_concrete and pe.poly_ring.ngens() > 1:
        try:
            fac = pe.poly.factor()
        except Exception:
            fac = None
        if fac is not None:
            nontrivial = [(g, e) for g, e in fac if g.total_degree() > 0]
            if len(nontrivial) >= 2:
                comps = [_component(pe, g, note=f"component {g} = 0 of the "
                                                "factored equation")
                         for g, _ in nontrivial]
                return Classification(
                    parsed=pe, components=comps,
                    note="the polynomial factors; the solution set is the "
                         "union of the components",
                )
            if len(nontrivial) == 1 and nontrivial[0][1] > 1:
                g, exponent = nontrivial[0]
                # g^e = 0 and g = 0 have the same solutions in the same
                # coordinates, so the reduction is the identity -- recorded
                # so that a solver knows which problem it answered
                return _component(
                    pe, g,
                    reduction=identity(
                        pe.unknowns,
                        description=f"reduced ({g})^{exponent} = 0 to "
                                    f"{g} = 0"),
                    note=f"input is ({g})^{exponent} = 0; "
                         "reduced to the underlying factor",
                )

    return Classification(parsed=pe, matches=_registered_matches(pe))
