# DiophantineClassifier

Classify Diophantine equations into named families — the classification engine
for the [Diophantine Library](https://github.com/DiophantineLibrary), an
online resource (in planning) to help people **solve, study and classify
Diophantine equations** ([vision talk, LuCaNT 2 (2025)](https://app.icerm.brown.edu/assets/515/9255/9255_5454_Roe_070720250930_Slides.pdf)).

Given an equation, the classifier reports the most specific known *family*
containing it — Pell, Thue, Mordell, generalized Fermat, Ramanujan–Nagell, … —
together with the family's status (solved / algorithmic / effective / open),
the standard methods and references, runnable code for Sage/PARI/Magma, and,
where standard software gives a complete answer, the solutions themselves.

Built on [SageMath](https://www.sagemath.org); some solver suggestions target
Magma or PARI. Designed to serve both as the website backend and as a
standalone library.

## Quick start

The package must run under Sage's Python:

```bash
sage -pip install -e .        # or: use sage -python from the repo root
```

```python
sage: from diophantine_classifier import classify, solve

sage: print(classify("x^2 - 61*y^2 = 1").explain())
x^2 - 61*y^2 = 1
  family: pell — Pell equation   [P1, solved]
  summary: Pell equation with D = 61
  data: D=61, N=1
  lineage: pell-like → binary-qf-representation → binary-quadratic → quadric → general-polynomial
  ...

sage: solve("x^2 - 61*y^2 = 1").solutions[0]
(1766319049, 226153980)

sage: classify("y^2 = x^3 + k", params="k")        # parametric families
Classification('y^2 = x^3 + k' -> mordell)

sage: classify("x^2 + 7 = 2^n")                    # variable exponents
Classification('x^2 + 7 = 2^n' -> ramanujan-nagell)

sage: solve("x^3 + 2*y^3 = 11").solutions          # PARI thue, certified
[(3, -2)]
```

Command line:

```bash
sage -python -m diophantine_classifier.cli "x^2 + 7 = 2^n" --solve
dioclassify "4/n = 1/x + 1/y + 1/z" --params n --json   # after install
```

## The family registry

The mathematical content lives in two synchronized places:

- **[docs/FAMILIES.md](docs/FAMILIES.md)** — a prioritized, referenced
  enumeration of ~60 families of Diophantine equations: standard forms,
  solvability status, methods, software, and the specialization DAG, from
  linear equations through the undecidability boundary of Hilbert's tenth
  problem.
- **[diophantine_classifier/data/families.yaml](diophantine_classifier/data/families.yaml)**
  — the machine-readable registry driving the classifier: one entry per
  family with DAG edges (`parents`), priority tier (P1 = launch set for the
  Library), status, software pointers, and fillable code templates.

An equation can belong to many families (`x² − 61y² = 1` is a Pell equation,
hence a binary quadratic form representation, hence a conic…); the classifier
reports the most specific match and the full lineage.

## What works today (wave 1)

- **Parsing**: polynomial equations over ℤ with optional named parameters
  (`y^2 = x^3 + k`), variable exponents (`2^n`, `y^q`), and unit fractions
  (`4/n = 1/x + 1/y + 1/z`), with denominator-clearing tracked as conditions.
- **Structural matchers** for ~45 families: linear, univariate, Pell/Pell-like,
  binary and n-ary quadratic forms (isotropy and representation), sums of 2/3/4
  squares, Legendre/Pythagorean, Weierstrass/Mordell/quartic genus-one models,
  Thue and reducible binary forms, Thue–Mahler, hyperelliptic/genus-2/
  superelliptic, generalized Fermat (with the χ trichotomy and odd-exponent
  sign normalization), sums of three cubes, Markov–Hurwitz, Waring/diagonal,
  equal sums of like powers, Catalan/Pillai, Ramanujan–Nagell, Lebesgue–Nagell,
  Schinzel–Tijdeman power values, S-unit-type, Egyptian fractions/Erdős–Straus.
- **Geometry fallback**: reducible equations split into components; irreducible
  plane curves are routed by genus (0 → parametrize, 1 → needs a point +
  Nagell, ≥ 2 → Faltings/Chabauty), computed via Sage.
- **Solvers** for the families where standard software is definitive:
  linear (Bezout + lattice), Pell (continued fractions), generalized Pell
  (PARI `qfbsolve`), quadratic form isotropy (`qfsolve`, with local
  obstructions reported), 2/3/4 squares, `BinaryQF.solve_integer`,
  integral points on Weierstrass models (`E.integral_points`), Thue
  (PARI `thue`, certified), plus literature-complete answers (Catalan,
  Fermat, classical Ramanujan–Nagell).
- **Output for the website**: `Classification.as_dict()` is JSON-ready;
  `explain()` is the human-readable equation-homepage prototype.

See [docs/DESIGN.md](docs/DESIGN.md) for the architecture and the roadmap
(wave 2: transformations — completing the square, GL₂(ℤ) reduction of binary
forms, Nagell's algorithm after point search; wave 3: Igusa-invariant
identification in genus 2, systems, positivity domains; wave 4: L-function
hashing for higher genus, database layer).

## Tests

```bash
make test          # = sage -python -m pytest tests -q
```

The corpus in `tests/test_classify.py` doubles as a showcase: fifty famous
equations (Pell, Selmer's cubic, the Klein quartic, Elkies' quartic, 33 and 42
as sums of three cubes, …) with their expected families.

## License

GPL-3.0 (see LICENSE), matching the Sage ecosystem.
