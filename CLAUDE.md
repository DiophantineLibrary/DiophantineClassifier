# DiophantineClassifier — notes for Claude sessions

Classification engine for the Diophantine Library
(github.com/DiophantineLibrary). Runs under **Sage's Python only**.

## Commands

- Unit tests: `sage -python -m pytest tests -q` (or `make test`)
- Doctests: `make doctest` (= `PYTHONPATH=. sage -t diophantine_classifier/`)
- Docstring coverage: `make coverage` (must stay 100%)
- References pipeline: `make references` (regenerates references/REPORT.md)
- Smoke: `make smoke`
- CI: `.github/workflows/ci.yml` runs unit tests, doctests, coverage and the
  references pipeline in the `sagemath/sagemath` container on every PR; that
  image has no `make`, so the steps spell out the commands
- CLI: `sage -python -m diophantine_classifier.cli "x^2 - 61*y^2 = 1" --solve`  *(later in the series)*

## Architecture (see docs/DESIGN.md)

The architecture of the finished classifier; entries marked *(later in the series)*
are specified here but land in a later PR of the series.

- `parsing.py` — string → `ParsedEquation` (term model: polynomial + `2^n` +
  `y^q` terms; params live in the coefficient ring; denominator clearing is
  recorded in `conditions`).
- `conditions.py` — the side conditions themselves: structured, exactly
  evaluated, three-valued; collected from the source syntax before Sage can
  cancel a denominator away.
- `data/families/<slug>.yaml` — one registry file per family (slug ==
  filename; DAG `parents`, priority, status, software, code templates,
  annotated references). `registry.py` loads the directory, and exposes
  `lineage_paths` / `lineage_graph` (genuine edges) beside flat `ancestors`.
- `data/references.bib` — bibliography; `references.py` parses/formats it;
  `tools/check_references.py` validates it (and local PDFs in
  `references/pdf/<key>.pdf`).
- `matchers.py` — shape recognizers emitting `Match(slug, data, transform)`.
- `transforms.py` — `CoordinateTransform`: the executable, invertible map
  from the user's variables to a family's standard coordinates, plus the
  structural roles. Solvers work normalized and `pull_back()`.
- `classify.py` — factor-split, run matchers, rank by DAG depth (most
  specific family wins), `explain()` / `as_dict()`.
- `solvers.py` *(later in the series)* — per-family solvers (Sage/PARI); `SolutionSet` is iterable
  (streams for infinite families); `SolverUnavailable` carries code
  templates for the rest.
- `docs/FAMILIES.md` — the human-readable, referenced enumeration.

## Invariants (tests enforce most of these)

- Every function (private helpers included) has a Sage-convention docstring
  with INPUT/OUTPUT (where nontrivial) and EXAMPLES that pass `sage -t`;
  doctests import what they need explicitly.
- Every slug emitted by `matchers.py` exists in the registry with
  `matcher: true`.
- Registry `parents` form a DAG; priorities in {1,2,3}; statuses from
  `registry.STATUSES`; every reference key resolves in `references.bib` with
  a nonempty `why`; every family cites at least one reference.
- `docs/FAMILIES.md` and the per-family YAML describe the same families —
  update both when adding one.
- `as_dict()` must stay JSON-serializable (website backend contract).
- Adding a family = `data/families/<slug>.yaml` + FAMILIES.md entry + bib
  entries + matcher (+ solver if standard software is definitive) + corpus
  row in `tests/test_classify.py`.

## Conventions

- Solutions/tuples are ordered by the unknowns' order of appearance in the
  input equation; infinite families expose their enumeration via
  `iter(solution_set)` / `.first(n)`.
- Matchers never mutate the parsed equation; normalizations live in
  `Match.transform` as an executable map, never only as prose.
- Bibliography `url` fields must point to legally free copies only.
- No web dependencies; no Magma requirement (templates only).
- Do not push or open PRs without being asked.
