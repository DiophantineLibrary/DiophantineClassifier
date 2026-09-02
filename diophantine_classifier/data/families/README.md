# Family registry files

One YAML file per family, named `<slug>.yaml` (the loader errors if the
`slug` field disagrees with the filename). One file per family keeps merge
conflicts local when several people add families at once.

Prose documentation with full mathematical context lives in
`docs/FAMILIES.md` — keep the two in sync. References live in
`../references.bib`; entries here point to them by key.

## Schema

| field        | type            | meaning |
|--------------|-----------------|---------|
| `slug`       | string          | unique kebab-case id; equals the filename stem; stable (future site URLs) |
| `name`       | string          | display name |
| `priority`   | 1, 2, 3         | Library rollout tier (1 = launch set) |
| `status`     | enum            | `solved`, `algorithmic`, `effective`, `ineffective`, `partial`, `open`, `undecidable` |
| `class`      | string          | coarse group (`linear`/`quadratic`/`genus1`/`curve`/`fermat`/`expdioph`/`normform`/`surface`/`special`/`boundary`) |
| `form`       | string          | standard form; unknowns as displayed, other letters are parameters |
| `constraints`| string          | side conditions under which the family's statements hold: domain of the unknowns, sign/size conditions on parameters (e.g. `x, y, p, q >= 2`); empty means the bare equation over the domain |
| `aliases`    | list of strings | alternative names (optional) |
| `parents`    | list of slugs   | immediate generalizations; DAG edges used for match ranking |
| `matcher`    | bool            | whether `matchers.py` can emit this slug (tests enforce) |
| `methods`    | list of strings | standard solution/analysis methods |
| `software`   | map             | system (`sage`/`pari`/`magma`/`other`) → short description |
| `code`       | map             | system → code template with `{placeholder}` slots filled from match data |
| `examples`   | list of strings | equations parsable by `parse()` that belong here (tests classify them) |
| `references` | list of maps    | `{key: <bibtex key>, why: <what this reference contributes to THIS family>}` |
| `notes`      | string          | free-form mathematical notes shown to users |
| `lmfdb`      | string          | related LMFDB collection slug (optional) |
| `finiteness` | string          | optional note on the solution set shape |

## The two exponential parents

`polynomial-exponential` is the umbrella for every equation mixing
polynomial and exponential terms, *including* unknown bases raised to
unknown exponents (`x^p - y^q = 1`).  `exponential-diophantine` is its
purely exponential subclass: the bases are fixed integers and only the
exponents are unknown (`2^a + 3^b = 5^c`).  A family whose form has an
unknown base therefore descends from `polynomial-exponential`, not from
`exponential-diophantine`, even though both carry `class: expdioph`.

## Adding a family

1. Create `<slug>.yaml` here (and an entry in `docs/FAMILIES.md`).
2. Add references to `../references.bib` (with `doi` and a legally-free
   `url` when available) and cite them by key with a `why` annotation.
3. If the standard form is structurally recognizable, add a matcher in
   `matchers.py` and set `matcher: true`; add a corpus row in
   `tests/test_classify.py`.
4. Run `make test doctest` and `sage -python tools/check_references.py`.
