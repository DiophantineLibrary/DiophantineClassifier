"""Family registry: loads data/families.yaml and exposes the specialization DAG.

The registry is the single source of truth about families (names, status,
references, software, code templates).  Matchers refer to families by slug;
tests enforce that the two stay consistent.
"""

import os
from dataclasses import dataclass, field
from functools import lru_cache

import yaml

DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "families.yaml")

STATUSES = ("solved", "algorithmic", "effective", "ineffective", "partial", "open", "undecidable")


@dataclass(frozen=True)
class Family:
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
        return f"Family({self.slug!r})"

    def fill_code(self, data):
        """Return code templates with {placeholders} filled from a match's data.

        Templates whose placeholders are not all available are returned with the
        placeholders left intact (still useful as a recipe).
        """
        out = {}
        for lang, template in self.code.items():
            try:
                out[lang] = template.format(**{k: str(v) for k, v in data.items()})
            except (KeyError, IndexError):
                out[lang] = template
        return out


def _tuple(x):
    if x is None:
        return ()
    if isinstance(x, (list, tuple)):
        return tuple(x)
    return (x,)


@lru_cache(maxsize=None)
def families():
    """Return the registry as an ordered dict slug -> Family."""
    with open(DATA_PATH) as fobj:
        raw = yaml.safe_load(fobj)
    result = {}
    for entry in raw:
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
            references=_tuple(entry.get("references")),
            notes=(entry.get("notes") or "").strip(),
            lmfdb=entry.get("lmfdb", "") or "",
            finiteness=entry.get("finiteness", "") or "",
        )
        if fam.slug in result:
            raise ValueError(f"duplicate family slug {fam.slug!r}")
        result[fam.slug] = fam
    return result


def family(slug):
    return families()[slug]


@lru_cache(maxsize=None)
def depth(slug):
    """Length of the longest specialization chain from a root down to ``slug``.

    Used as the specificity score when ranking matches (deeper = more specific).
    """
    fam = families()[slug]
    if not fam.parents:
        return 0
    return 1 + max(depth(p) for p in fam.parents)


def ancestors(slug):
    """All strict ancestors of ``slug`` in the DAG, most specific first."""
    seen = []
    frontier = list(families()[slug].parents)
    while frontier:
        nxt = frontier.pop(0)
        if nxt not in seen:
            seen.append(nxt)
            frontier.extend(families()[nxt].parents)
    return sorted(set(seen), key=lambda s: -depth(s))


def validate():
    """Sanity-check the registry; raises on inconsistency.  Used by the tests."""
    fams = families()
    for fam in fams.values():
        for parent in fam.parents:
            if parent not in fams:
                raise ValueError(f"{fam.slug}: unknown parent {parent!r}")
        if fam.priority not in (1, 2, 3):
            raise ValueError(f"{fam.slug}: bad priority {fam.priority}")
        if fam.status not in STATUSES:
            raise ValueError(f"{fam.slug}: bad status {fam.status!r}")
        depth(fam.slug)  # raises RecursionError on a cycle
    return True
