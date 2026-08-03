r"""Classification pipeline: parse, split reducible equations, run matchers,
rank by specificity in the family DAG, and package the result.

EXAMPLES::

    sage: from diophantine_classifier import classify
    sage: classify("x^2 - 61*y^2 = 1")
    Classification('x^2 - 61*y^2 = 1' -> pell)
    sage: classify("y^2 = x^3 + k", params="k").slug
    'mordell'
    sage: cls = classify("(x^2 - 2)*(y^2 - 3) = 0")   # reducible: components
    sage: cls.slug, [c.slug for c in cls.components]
    ('reducible', ['univariate', 'univariate'])
"""

from dataclasses import dataclass, field

from sage.all import ZZ

from . import matchers
from .parsing import ParsedEquation, parse
from .registry import ancestors, depth, families


def _jsonify(x):
    r"""Recursively convert match data to plain JSON-serializable types.

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
    r"""The result of :func:`classify`.

    ATTRIBUTES:

    - ``parsed`` -- the underlying
      :class:`~diophantine_classifier.parsing.ParsedEquation`.
    - ``matches`` -- list of :class:`~diophantine_classifier.matchers.Match`,
      sorted by decreasing specificity (DAG depth, ties by emission order).
      The first entry is the primary match.
    - ``components`` -- list of :class:`Classification`; nonempty exactly
      when the equation's polynomial factors, in which case the solution set
      is the union over the components and there is no primary match.
    - ``note`` -- string; pipeline remark (e.g. that a repeated factor was
      reduced).

    Derived properties: :attr:`primary`, :attr:`slug`, :attr:`family`,
    :attr:`lineage`, :attr:`data`.

    EXAMPLES::

        sage: from diophantine_classifier import classify
        sage: cls = classify("x^2 - 61*y^2 = 1")
        sage: cls.slug
        'pell'
        sage: cls.lineage[:2]
        ['pell-like', 'binary-qf-representation']
        sage: cls.data["D"]
        '61'
    """
    parsed: ParsedEquation
    matches: list = field(default_factory=list)
    components: list = field(default_factory=list)
    note: str = ""

    # ------------------------------------------------------------------ api

    @property
    def is_composite(self):
        r"""Whether the equation split into factor components.

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
        r"""The most specific match, or ``None`` for composite equations.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("x^2 + 7 = 2^n").primary
            Match('ramanujan-nagell')
        """
        return self.matches[0] if self.matches else None

    @property
    def slug(self):
        r"""The primary family's slug (``"reducible"`` for composites).

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("x^3 + 2*y^3 = 11").slug
            'thue'
        """
        if self.is_composite:
            return "reducible"
        return self.primary.slug if self.primary else "unclassified"

    @property
    def family(self):
        r"""The primary :class:`~diophantine_classifier.registry.Family`.

        ``None`` for composite equations.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("x^2 - 61*y^2 = 1").family
            Family('pell')
        """
        return families().get(self.slug)

    @property
    def lineage(self):
        r"""Ancestors of the primary family, most specific first.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("y^2 = x^3 - 2").lineage[0]
            'elliptic-weierstrass'
        """
        if self.primary is None:
            return []
        return ancestors(self.primary.slug)

    @property
    def data(self):
        r"""The primary match's extracted data.

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("y^2 = x^3 - 2").data["k"]
            '-2'
        """
        return self.primary.data if self.primary else {}

    def code(self):
        r"""Filled code templates for the primary family and its ancestors.

        OUTPUT: dict mapping ``"<system> (<slug>)"`` to a code string

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: classify("x^2 - 61*y^2 = 1").code()["pari (pell)"]
            'quadunit(4*61)'
        """
        out = {}
        for slug in ([self.slug] + self.lineage) if self.primary else []:
            fam = families().get(slug)
            if fam and fam.code:
                filled = fam.fill_code(self.data)
                for lang, snippet in filled.items():
                    out.setdefault(f"{lang} ({slug})", snippet)
        return out

    def as_dict(self):
        r"""JSON-serializable summary — the website-backend contract.

        OUTPUT: dict with plain types only (tested); keys include
        ``equation``, ``family``, ``status``, ``data``, ``lineage``,
        ``references`` (list of ``{key, why, formatted}`` dicts), ``code``,
        and ``components`` for composites

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: d = classify("x^2 - 61*y^2 = 1").as_dict()
            sage: d["family"], d["priority"]
            ('pell', 1)
            sage: d["references"][0]["key"]
            'Lenstra2002'
            sage: import json
            sage: _ = json.dumps(d)          # round-trips
        """
        d = {
            "equation": self.parsed.original,
            "domain": self.parsed.domain,
            "params": list(self.parsed.params),
            "unknowns": list(self.parsed.unknowns),
            "conditions": list(self.parsed.conditions),
            "family": self.slug,
        }
        if self.note:
            d["note"] = self.note
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
        if self.primary:
            d["summary"] = self.primary.summary
            d["data"] = _jsonify(self.primary.data)
            if self.primary.transform:
                d["transform"] = self.primary.transform
        d["lineage"] = self.lineage
        d["all_matches"] = [m.slug for m in self.matches]
        d["code"] = self.code()
        return d

    # ------------------------------------------------------------- display

    def explain(self):
        r"""Multi-line human-readable report (equation-homepage prototype).

        OUTPUT: string

        EXAMPLES::

            sage: from diophantine_classifier import classify
            sage: text = classify("x^2 + 7 = 2^n").explain()
            sage: "family: ramanujan-nagell" in text
            True
            sage: "reference [Nagell1961]" in text
            True
        """
        lines = [self.parsed.original]
        if self.parsed.params:
            lines.append(f"  parameters: {', '.join(self.parsed.params)}")
        if self.note:
            lines.append(f"  note: {self.note}")
        for cond in self.parsed.conditions:
            lines.append(f"  condition: {cond}")
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
        showable = {k: v for k, v in self.primary.data.items()
                    if k not in ("roles",)}
        if showable:
            lines.append("  data: " + ", ".join(f"{k}={v}"
                                                for k, v in showable.items()))
        if self.lineage:
            lines.append("  lineage: " + " → ".join(self.lineage))
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
        r"""Terse representation.

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
    r"""Sort matches by decreasing DAG depth, ties by emission order.

    EXAMPLES::

        sage: from diophantine_classifier.parsing import parse
        sage: from diophantine_classifier.matchers import run
        sage: from diophantine_classifier.classify import _rank
        sage: [m.slug for m in _rank(run(parse("x^2 - 61*y^2 = 1")))][:2]
        ['pell', 'pell-like']
    """
    indexed = list(enumerate(match_list))
    indexed.sort(key=lambda pair: (-depth(pair[1].slug), pair[0]))
    return [m for _, m in indexed]


def classify(equation, params=(), domain="ZZ"):
    r"""Classify a Diophantine equation into the most specific known family.

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
        sage: classify("x^2 - 61*y^2 = 1").slug
        'pell'
        sage: classify("y^2 = x^3 + k", params="k").slug
        'mordell'
        sage: classify("x^2 + 7 = 2^n").slug
        'ramanujan-nagell'
        sage: classify("4/n = 1/x + 1/y + 1/z", params="n").slug
        'erdos-straus'
        sage: classify("3*x^3 + 4*y^3 + 5*z^3 = 0").slug     # Selmer
        'generalized-fermat'

    TESTS:

    Repeated factors are reduced::

        sage: cls = classify("(x + y)^2 = 0")
        sage: cls.slug
        'linear'
        sage: "reduced to the underlying factor" in cls.note
        True
    """
    if isinstance(equation, ParsedEquation):
        pe = equation
    else:
        pe = parse(equation, params=params, domain=domain)

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
                comps = [classify(f"{g} = 0", domain=domain)
                         for g, _ in nontrivial]
                return Classification(
                    parsed=pe, components=comps,
                    note="the polynomial factors; the solution set is the "
                         "union of the components",
                )
            if len(nontrivial) == 1 and nontrivial[0][1] > 1:
                g = nontrivial[0][0]
                sub = classify(f"{g} = 0", domain=domain)
                sub.note = (f"input is ({g})^{nontrivial[0][1]} = 0; "
                            "reduced to the underlying factor")
                return sub

    matches = _rank(matchers.run(pe))
    return Classification(parsed=pe, matches=matches)
