# DiophantineClassifier — notes for Claude sessions

Classification engine for the Diophantine Library
(github.com/DiophantineLibrary). Runs under **Sage's Python only**.

## Commands

- Tests: `sage -python -m pytest tests -q` (or `make test`)
- Smoke: `make smoke`
- CLI: `sage -python -m diophantine_classifier.cli "x^2 - 61*y^2 = 1" --solve`

## Architecture (see docs/DESIGN.md)

- `parsing.py` — string → `ParsedEquation` (term model: polynomial + `2^n` +
  `y^q` terms; params live in the coefficient ring; denominator clearing is
  recorded in `conditions`).
- `data/families.yaml` — the family registry (slugs, DAG `parents`, priority,
  status, software, code templates). `registry.py` loads it.
- `matchers.py` — shape recognizers emitting `Match(slug, data)`.
- `classify.py` — factor-split, run matchers, rank by DAG depth (most
  specific family wins), `explain()` / `as_dict()`.
- `solvers.py` — per-family solvers (Sage/PARI); `SolverUnavailable` carries
  code templates for the rest.
- `docs/FAMILIES.md` — the human-readable, referenced enumeration.

## Invariants (tests enforce most of these)

- Every slug emitted by `matchers.py` exists in the YAML with `matcher: true`.
- YAML `parents` form a DAG; priorities in {1,2,3}; statuses from the fixed
  vocabulary (`registry.STATUSES`).
- `docs/FAMILIES.md` and `families.yaml` describe the same families — update
  both when adding one.
- `as_dict()` must stay JSON-serializable (website backend contract).
- Adding a family = YAML entry + FAMILIES.md entry + matcher (+ solver if
  standard software is definitive) + corpus row in `tests/test_classify.py`.

## Conventions

- Solutions/tuples are ordered by the unknowns' order of appearance in the
  input equation.
- Matchers never mutate the parsed equation; normalizations are described in
  `Match.transform`.
- No web dependencies; no Magma requirement (templates only).
- Do not push or open PRs without being asked.
