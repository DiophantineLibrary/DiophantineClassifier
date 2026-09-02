# DiophantineClassifier

Classify Diophantine equations into named families — the classification engine
for the [Diophantine Library](https://github.com/DiophantineLibrary), an
online resource (in planning) to help people **solve, study and classify
Diophantine equations**.

Given an equation, the classifier reports the most specific known *family*
containing it — Pell, Thue, Mordell, generalized Fermat, Ramanujan–Nagell, … —
together with the family's status (solved / algorithmic / effective / open),
the standard methods and references, runnable code for Sage/PARI/Magma, and,
where standard software gives a complete answer, the solutions themselves.

Built on [SageMath](https://www.sagemath.org); some solver suggestions target
Magma or PARI. Designed to serve both as the website backend and as a
standalone library.

> **Status of this branch.** This is one PR of a stacked series that builds the
> classifier layer by layer; the description above is where the series lands.
> What runs *here* is the parser, the registry, the bibliography, the matchers
> and the classification pipeline — the examples below all work on this branch,
> and they use the three families registered so far. The solvers and the
> remaining 61 families arrive in the later PRs; the closing PR restores the
> full README.

## Quick start

The package must run under Sage's Python:

```bash
sage -pip install -e .        # or: use sage -python from the repo root
```

```python
sage: from diophantine_classifier import classify, parse, families

sage: classify("3*x + 5*y = 1")
Classification('3*x + 5*y = 1' -> linear)
sage: classify("x^2 - 5*x + 6 = 0").slug
'univariate'
sage: classify("x^2 + y^3 + z^5 = 7").slug     # last-resort answer
'general-polynomial'

sage: cls = classify("(x^2 - 2)*(y^2 - 3) = 0")    # reducible: components
sage: cls.slug, [c.slug for c in cls.components]
('reducible', ['univariate', 'univariate'])

sage: sorted(families())                       # the registry, so far
['general-polynomial', 'linear', 'univariate']
```

## The family registry

The mathematical content lives in synchronized places:

- **[docs/FAMILIES.md](docs/FAMILIES.md)** — a prioritized, referenced
  enumeration of families of Diophantine equations: standard forms,
  solvability status, methods, software, and the specialization DAG, from
  linear equations through the undecidability boundary of Hilbert's tenth
  problem. Three families are enumerated on this branch; each of the
  remaining ~60 arrives in its own PR, with its prose entry, its registry
  file and its references reviewed together.
- **[diophantine_classifier/data/families/](diophantine_classifier/data/families)**
  — the machine-readable registry driving the classifier: one YAML file per
  family (merge-friendly), with DAG edges (`parents`), priority tier (P1 =
  launch set for the Library), status, software pointers, fillable code
  templates, and annotated references. Schema in that directory's README.
- **[diophantine_classifier/data/references.bib](diophantine_classifier/data/references.bib)**
  — the bibliography. Families cite entries by BibTeX key, each use carrying
  a `why` annotation (what the reference contributes *to that family*);
  entries carry DOIs and links to legally free copies where available.
  `sage -python tools/check_references.py` validates the whole apparatus and
  checks locally downloaded PDFs against their entries (see
  [references/README.md](references/README.md) and the generated
  [references/REPORT.md](references/REPORT.md)).

An equation can belong to many families (`x² − 61y² = 1` is a Pell equation,
hence a binary quadratic form representation, hence a conic…); the classifier
reports the most specific match and the full lineage.

## What runs on this branch

- **Parsing**: polynomial equations over ℤ with optional named parameters
  (`y^2 = x^3 + k`), variable exponents (`2^n`, `y^q`), and unit fractions
  (`4/n = 1/x + 1/y + 1/z`), with denominator-clearing tracked as conditions.
- **The registry**: `data/families/*.yaml` loaded and validated (parents
  resolve, the DAG is acyclic, priorities and statuses come from fixed
  vocabularies, every reference key resolves with a `why` annotation), with
  `depth`/`ancestors` giving the specificity order used to rank matches.
- **The bibliography and its pipeline**: BibTeX parsing and display
  formatting, plus `tools/check_references.py` and its monotone verification
  ledger.
- **Structural matchers** for ~45 families, and a genus-based geometry
  fallback routing irreducible plane curves by genus.
- **Classification**: reducible equations split into components, matches are
  ranked by depth in the family DAG so the most specific family wins, and
  `explain()` / `as_dict()` give the human report and the JSON contract for
  the website backend.

## Coming in the rest of the series

- **Solvers** for the families where standard software is definitive, with
  iterable solution sets for infinite families, and the command-line
  interface.
- **The remaining families**, one PR each.

See [docs/DESIGN.md](docs/DESIGN.md) for the target architecture and the
roadmap beyond this series (wave 2: transformations — completing the square,
GL₂(ℤ) reduction of binary forms, Nagell's algorithm after point search;
wave 3: Igusa-invariant identification in genus 2, systems, positivity
domains; wave 4: L-function hashing for higher genus, database layer).

## Tests

```bash
make test          # unit tests:  sage -python -m pytest tests -q
make doctest       # doctests:    sage -t diophantine_classifier/
make coverage      # docstring coverage: sage --coverage (100%)
make references    # bibliography + local-PDF validation pipeline
make check         # test + doctest
```

GitHub Actions runs those same checks on every pull request, inside the
official `sagemath/sagemath` container ([`.github/workflows/ci.yml`](.github/workflows/ci.yml));
that image ships no `make`, so the workflow spells out the commands and names
the target each one mirrors.

Every function carries a Sage-convention docstring (INPUT/OUTPUT/EXAMPLES)
whose examples run under `sage -t`.

## License

GPL-3.0 (see LICENSE), matching the Sage ecosystem.
