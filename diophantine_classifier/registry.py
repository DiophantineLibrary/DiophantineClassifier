r"""
Family registry: load ``data/families/*.yaml`` and expose the DAG.

Each family lives in its own YAML file named ``<slug>.yaml`` (one file per
family keeps merge conflicts local when contributors add families).  The
registry is the single source of truth about families — names, status,
software pointers, code templates, and annotated references into
``data/references.bib``.  Matchers refer to families by slug; the tests
enforce that the two stay consistent.

EXAMPLES::

    sage: from diophantine_classifier import families, family
    sage: fam = family("pell")
    sage: fam.name
    'Pell equation'
    sage: fam.parents
    ('pell-like',)
    sage: len(families()) > 60
    True
"""

import os
from dataclasses import dataclass, field
from functools import lru_cache

import yaml

from .references import bibliography, format_reference

#: directory holding one YAML file per family
DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "families")

#: allowed values of :attr:`Family.status`
STATUSES = ("solved", "algorithmic", "effective", "ineffective", "partial",
            "open", "undecidable")


@dataclass(frozen=True)
class Family:
    r"""
    One family of Diophantine equations, as recorded in the registry.

    Instances are loaded from ``data/families/<slug>.yaml`` and are immutable.

    ATTRIBUTES:

    - ``slug`` -- string; unique kebab-case identifier (equals the YAML
      filename stem), e.g. ``"pell"``.  Stable: the future website will use
      it in URLs.
    - ``name`` -- string; display name, e.g. ``"Pell equation"``.
    - ``priority`` -- integer in {1, 2, 3}; rollout tier for the Diophantine
      Library (1 = launch set, 2 = second wave, 3 = aspirational).
    - ``status`` -- string from :data:`STATUSES`; solvability status of the
      family as a whole: ``solved`` (complete description of solutions),
      ``algorithmic`` (terminating implemented algorithm), ``effective``
      (effective finiteness, instance-by-instance work), ``ineffective``
      (finiteness without algorithm), ``partial``, ``open``, ``undecidable``.
    - ``klass`` -- string; coarse group used for organizing displays
      (``linear``/``quadratic``/``genus1``/``curve``/``fermat``/``expdioph``/
      ``normform``/``surface``/``special``/``boundary``).
    - ``form`` -- string; the family's standard form, with unknowns as
      displayed and remaining letters understood as parameters.
    - ``parents`` -- tuple of slugs; immediate generalizations.  These edges
      form the specialization DAG used to rank matches (deeper = more
      specific).
    - ``matcher`` -- bool; whether :mod:`~diophantine_classifier.matchers`
      can emit this slug (tests enforce consistency).
    - ``aliases`` -- tuple of strings; alternative names.
    - ``methods`` -- tuple of strings; the standard solution/analysis methods.
    - ``software`` -- dict mapping system name (``sage``/``pari``/``magma``/
      ``other``) to a short description of the relevant functionality.
    - ``code`` -- dict mapping system name to a code template with
      ``{placeholder}`` slots filled from a match's ``data`` by
      :meth:`fill_code`.
    - ``examples`` -- tuple of strings; equations parsable by
      :func:`~diophantine_classifier.parsing.parse` that belong to the family
      (tests classify each and check the result).
    - ``references`` -- tuple of dicts ``{"key": ..., "why": ...}``: BibTeX
      keys into ``data/references.bib`` together with an annotation of what
      the reference contributes *to this family specifically*.
    - ``notes`` -- string; free-form mathematical notes shown to users.
    - ``lmfdb`` -- string; slug of a related LMFDB collection (e.g. ``ec.q``),
      empty if none.
    - ``finiteness`` -- string; optional note on the shape of the solution
      set (``finite``/``infinite``/``mixed``), empty if unset.

    EXAMPLES::

        sage: from diophantine_classifier import family
        sage: fam = family("thue")
        sage: fam.priority, fam.status
        (1, 'algorithmic')
        sage: sorted(ref["key"] for ref in fam.references)[-1]
        'TzanakisDeWeger1989'
        sage: all(ref["why"] for ref in fam.references)
        True
    """
    slug: str
    name: str
    priority: int
    status: str
    klass: str
    form: str
    parents: tuple = ()
    matcher: bool = False
    aliases: tuple = ()
    methods: tuple = ()
    software: dict = field(default_factory=dict)
    code: dict = field(default_factory=dict)
    examples: tuple = ()
    references: tuple = ()
    notes: str = ""
    lmfdb: str = ""
    finiteness: str = ""

    def __repr__(self):
        r"""
        Terse representation.

        EXAMPLES::

            sage: from diophantine_classifier import family
            sage: family("pell")
            Family('pell')
        """
        return f"Family({self.slug!r})"

    def fill_code(self, data):
        r"""
        Return code templates with ``{placeholders}`` filled from ``data``.

        Templates whose placeholders are not all available are returned with
        the placeholders left intact (still useful as a recipe).

        INPUT:

        - ``data`` -- dict of match data (values are stringified)

        OUTPUT: dict mapping system name to code string

        EXAMPLES::

            sage: from diophantine_classifier import family
            sage: family("pell").fill_code({"D": 61})["pari"]
            'quadunit(4*61)'
            sage: family("pell").fill_code({})["pari"]   # unfilled
            'quadunit(4*{D})'
        """
        out = {}
        for lang, template in self.code.items():
            try:
                out[lang] = template.format(**{k: str(v) for k, v in data.items()})
            except (KeyError, IndexError):
                out[lang] = template
        return out

    def formatted_references(self):
        r"""
        Return references as triples ``(key, why, formatted)``.

        OUTPUT: list of 3-tuples of strings; ``formatted`` is the plain-text
        rendering of the BibTeX entry

        EXAMPLES::

            sage: from diophantine_classifier import family
            sage: key, why, fmt = family("catalan").formatted_references()[0]
            sage: key
            'Mihailescu2004'
            sage: "Catalan" in fmt
            True
        """
        bib = bibliography()
        out = []
        for ref in self.references:
            key = ref["key"]
            formatted = format_reference(key, bib) if key in bib \
                else f"[unknown key {key}]"
            out.append((key, ref.get("why", ""), formatted))
        return out


def _tuple(x):
    r"""
    Coerce a YAML scalar-or-list field to a tuple.

    EXAMPLES::

        sage: from diophantine_classifier.registry import _tuple
        sage: _tuple(None), _tuple("a"), _tuple(["a", "b"])
        ((), ('a',), ('a', 'b'))
    """
    if x is None:
        return ()
    if isinstance(x, (list, tuple)):
        return tuple(x)
    return (x,)


@lru_cache(maxsize=None)
def families():
    r"""
    Return the registry as a dict slug -> :class:`Family`.

    Families are loaded from ``data/families/*.yaml`` in alphabetical order;
    a file whose ``slug`` field does not match its filename is an error.

    OUTPUT: dict mapping slug to :class:`Family` (cached after first call)

    EXAMPLES::

        sage: from diophantine_classifier import families
        sage: fams = families()
        sage: "ramanujan-nagell" in fams
        True
        sage: all(fams[s].slug == s for s in fams)
        True
    """
    result = {}
    for fname in sorted(os.listdir(DATA_DIR)):
        if not fname.endswith(".yaml"):
            continue
        path = os.path.join(DATA_DIR, fname)
        with open(path) as fobj:
            entry = yaml.safe_load(fobj)
        stem = fname[:-len(".yaml")]
        if entry.get("slug") != stem:
            raise ValueError(
                f"{fname}: slug {entry.get('slug')!r} does not match filename")
        refs = tuple(
            ref if isinstance(ref, dict) else {"key": str(ref), "why": ""}
            for ref in (entry.get("references") or ()))
        fam = Family(
            slug=entry["slug"],
            name=entry["name"],
            priority=int(entry["priority"]),
            status=entry["status"],
            klass=entry.get("class", ""),
            form=entry.get("form", ""),
            parents=_tuple(entry.get("parents")),
            matcher=bool(entry.get("matcher", False)),
            aliases=_tuple(entry.get("aliases")),
            methods=_tuple(entry.get("methods")),
            software=entry.get("software") or {},
            code=entry.get("code") or {},
            examples=_tuple(entry.get("examples")),
            references=refs,
            notes=(entry.get("notes") or "").strip(),
            lmfdb=entry.get("lmfdb", "") or "",
            finiteness=entry.get("finiteness", "") or "",
        )
        result[fam.slug] = fam
    return result


def family(slug):
    r"""
    Return the :class:`Family` with the given slug.

    INPUT:

    - ``slug`` -- string

    EXAMPLES::

        sage: from diophantine_classifier import family
        sage: family("mordell").name
        'Mordell equation'
        sage: family("no-such-family")
        Traceback (most recent call last):
        ...
        KeyError: 'no-such-family'
    """
    return families()[slug]


@lru_cache(maxsize=None)
def depth(slug):
    r"""
    Length of the longest specialization chain from a root to ``slug``.

    Used as the specificity score when ranking matches: an equation matching
    several families is filed under the deepest one.

    INPUT:

    - ``slug`` -- string; a family slug

    OUTPUT: nonnegative integer

    EXAMPLES::

        sage: from diophantine_classifier.registry import depth
        sage: depth("general-polynomial")
        0
        sage: depth("pell") > depth("pell-like")
        True
    """
    fam = families()[slug]
    if not fam.parents:
        return 0
    return 1 + max(depth(p) for p in fam.parents)


def ancestors(slug):
    r"""
    All strict ancestors of ``slug`` in the DAG, most specific first.

    INPUT:

    - ``slug`` -- string; a family slug

    OUTPUT: list of slugs, sorted by decreasing :func:`depth`

    EXAMPLES::

        sage: from diophantine_classifier import ancestors
        sage: ancestors("pell")[:2]
        ['pell-like', 'binary-qf-representation']
        sage: ancestors("general-polynomial")
        []
    """
    seen = []
    frontier = list(families()[slug].parents)
    while frontier:
        nxt = frontier.pop(0)
        if nxt not in seen:
            seen.append(nxt)
            frontier.extend(families()[nxt].parents)
    return sorted(set(seen), key=lambda s: -depth(s))


def validate():
    r"""
    Sanity-check the registry; raises on inconsistency.

    Checks: parents exist, priorities and statuses come from the fixed
    vocabularies, the DAG is acyclic, and every reference key resolves in
    ``references.bib`` with a nonempty ``why`` annotation.

    OUTPUT: ``True`` (or an exception)

    EXAMPLES::

        sage: from diophantine_classifier import validate
        sage: validate()
        True
    """
    fams = families()
    bib = bibliography()
    for fam in fams.values():
        for parent in fam.parents:
            if parent not in fams:
                raise ValueError(f"{fam.slug}: unknown parent {parent!r}")
        if fam.priority not in (1, 2, 3):
            raise ValueError(f"{fam.slug}: bad priority {fam.priority}")
        if fam.status not in STATUSES:
            raise ValueError(f"{fam.slug}: bad status {fam.status!r}")
        for ref in fam.references:
            if ref["key"] not in bib:
                raise ValueError(
                    f"{fam.slug}: reference key {ref['key']!r} not in "
                    "references.bib")
            if not ref.get("why", "").strip():
                raise ValueError(
                    f"{fam.slug}: reference {ref['key']} lacks a 'why' "
                    "annotation")
        depth(fam.slug)  # raises on a cycle
    return True
