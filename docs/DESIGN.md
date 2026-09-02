# Design

> **Where this branch stands.** This document describes the architecture of
> the finished classifier, which the stacked PR series assembles one layer at
> a time. Sections about modules that have not landed in this branch are
> design specifications for the later PRs, not descriptions of code you can
> run here; `CLAUDE.md` marks which layers are present.

Goal (from the LuCaNT 2 talk): *given an input Diophantine equation, find a
change of coordinates transforming it into the form of a stored family*, then
serve everything the family knows — status, methods, software, solutions.
For non-experts, a name to Google and a transformation into standard form is
most of the value.

## Pipeline

```
input string
   │  parsing.parse            (SR-based; params, variable exponents,
   │                            denominator clearing with conditions)
   ▼
ParsedEquation  ── polynomial? ──► Sage polynomial ring over QQ or QQ[params]
   │
   │  classify.classify
   ├─ reducible over QQ?  → split into components (union of solution sets)
   ├─ matchers.run        → structural Matches (slug + extracted data)
   │     unit fractions │ polynomial branch │ exponential branch
   │     + genus routing fallback for irreducible plane curves
   ▼
ranked by DAG specificity (registry.depth); most specific = primary family
   │
   ▼
Classification ── explain() (equation-homepage prototype)
               ── as_dict() (JSON for the website backend)
               ── solvers.solve() (Sage/PARI where definitive)
               ── code() (filled Sage/PARI/Magma templates otherwise)
```

## The equation model

`parsing.Term` represents `c · Πparams^e · Πunknowns^e · Πb^n · Πv^n`:
rational coefficient, parameter powers (part of the coefficient), polynomial
powers, integer-base exponentials (`2^n`), and variable-base exponentials
(`y^q`). This covers every wave-1 family, including mixed
polynomial–exponential shapes (Ramanujan–Nagell, Lebesgue–Nagell) and
all-exponential ones (Pillai, S-unit).

Everything is normalized to `f = 0`. Multiplying through by denominators is
recorded in `conditions` (the transformation must transport solutions, so
introduced hypotheses are tracked, never silently dropped). The original
pre-clearing shape is kept for the unit-fraction matcher.

Parameters (`params="k"`) move symbols into the coefficient base ring
`QQ[params]`, so `y^2 = x^3 + k` classifies as the Mordell *family* rather
than a single equation.

## Matching and ranking

Matchers are shape recognizers: they only fire on (near-)standard forms, plus
cheap normalizations — global sign, variable swap in `D x² − y² = N`,
odd-exponent sign flips for diagonal equations, orientation of `v² = g(u)`
shapes over both variable orderings. Every match carries the extracted data
(D and N for Pell, a-invariants for Weierstrass, signature and χ for
generalized Fermat, …) and a `CoordinateTransform` (`transforms.py`): the
executable, invertible map from the user's variables to the family's standard
coordinates, together with the structural roles the family assigns them. A
prose note cannot transport a solution, and every wave-1 normalization is
already a coordinate change someone has to undo — `5x² − y² = 1` is a Pell
equation only after reading the user's `y` as the standard `x`. Operations
that change the equation but no coordinate (multiplying through by −1) are
recorded separately.

An equation matches many families; the registry's `parents` edges form a
specialization DAG and the *deepest* match wins. Ancestry is reported as
genuine paths (`lineage_paths`, `lineage_graph`) so the site can offer every
applicable page (Pell ⊂ Pell-like ⊂ binary quadratic forms ⊂ conics) without
drawing an edge the DAG does not have: a family with two parents has two
paths, and the flat `ancestors()` set is never rendered as one chain.

`classify` keeps the problem as submitted apart from the model the matchers
saw. Reducing `g^e = 0` to `g = 0` produces a different equation, so the
`Classification` holds both (`parsed` and `working`) plus the reduction
between them; display, serialization and solution validation all use the
original.

The genus fallback mirrors the talk's routing: conics and genus 0 are easy;
genus 1 needs an initial point (the hard step) after which Weierstrass
machinery is standard; genus 2 will get Igusa-invariant identification;
beyond that, an L-function-based hash is the plan (wave 4).

## Transformation ladder (the hard part, staged)

- **T0 (done)**: expand, clear denominators (tracked), merge, sign/orientation
  normalizations, variable permutations in fixed shapes.
- **T1 (wave 2)**: affine-unimodular changes of variables — completing the
  square to reduce any binary quadratic to Pell-like canonical form
  (solution-set bijections recorded as explicit substitutions).
- **T2 (wave 2)**: GL₂(ℤ) reduction of binary forms (Julia / Cremona–Stoll)
  so disguised Thue equations are recognized; minimization/reduction of
  genus-one models (Cremona–Stoll).
- **T3 (wave 3)**: birational transformations requiring search — point-finding
  on genus-1 models then Nagell's algorithm to Weierstrass form; two-covers
  for quartics. These change the model, so pullback maps must be stored with
  the classification.
- **T4 (wave 4)**: identification up to isomorphism/twist for genus ≥ 2 —
  Igusa (and Shioda) invariants in genus 2–3, L-function hash lookups against
  the Library's stored curves in higher genus.

Every transformation must come with a solution-transport recipe (a
substitution both ways, plus side conditions); "the classification is a map
of solution sets, not just a label."

## Solvers and code templates

`solvers.solve` dispatches over the matches that were actually emitted, most
specific first, handing each solver the match it is solving: a specialization
edge is not evidence that one family's data satisfies another's input
contract. Wired solvers use Sage/PARI (continued fractions, `qfbsolve`,
`qfsolve`, `two_squares`…, `BinaryQF.solve_integer`, `E.integral_points`,
PARI `thue` certified). Each declares the domains it answers over
(`SOLVER_DOMAINS`), so a rational question is declined rather than answered
with the integer solution set; `linear` and `univariate` implement ℤ, ℕ and ℚ
separately.

A solver works in its family's standard coordinates and never states the
answer itself: one wrapper pulls every tuple back through the match's
transform, orders it by the original equation's unknowns, and checks the
domain, the side conditions and the original equation — exactly, and lazily
for infinite families. Nothing reaches a user that does not solve the problem
they submitted.

Families without a definitive solver raise `SolverUnavailable`
carrying the registry's software pointers and *filled* code templates (e.g.
Magma `IntegralQuarticPoints`, Chabauty pipelines) — matching the site plan
that every equation homepage offers runnable code even where we don't run it.

Magma is optional everywhere: templates are emitted, nothing shells out in
wave 1. (Local Magma exists on this dev machine; a `magma` runner can be
added behind a feature check without changing the API.)

## Solution sets

`solve` returns a `SolutionSet` whose semantics are explicit and checked on
construction: `complete` says whether the listed solutions (with any symmetry
stated in the description) are provably all of them, an `infinite` set must
carry a stream rather than promise an enumeration it cannot deliver, and
`as_dict()` serializes coordinates exactly — a rational root stays `"1/2"`.
**Iteration** enumerates the full solution set even when it is infinite — Pell solutions by powers of the
fundamental unit, generalized-Pell automorph orbits walked in both
directions, lattice cosets by sup-norm shells, primitive Pythagorean triples
by Euclid's parametrization, the Markov tree in order of largest entry.
Finite representation problems (definite forms, unit fractions) are
enumerated completely rather than witnessed.

## References

The registry cites `data/references.bib` by key; every use carries a `why`
annotation (what the reference contributes to that family).  Bibliography
policy: `doi` when known-correct, `url` only for legally free copies.
`tools/check_references.py` is the checking/annotation pipeline: structural
validation, DOI/URL resolution (`--online`), and verification of locally
downloaded PDFs (`references/pdf/<key>.pdf`) against their entries via text
extraction.  Output: `references/REPORT.md` with per-entry status and TODO
lists.

## Website backend contract

- `Classification.as_dict()` returns plain JSON types only (tested), so the
  future site (LMFDB-style, psycodict/Postgres) can store and render results
  directly; references arrive formatted with their annotations.
- The per-family registry files are the site's families table in embryo:
  slugs are stable identifiers, `priority` is the rollout order, `lmfdb`
  fields link to existing LMFDB collections (elliptic curves, genus 2
  curves).  One file per family keeps community contributions
  merge-conflict-free.
- No web dependencies in this package, ever; the site imports us, not vice
  versa.

## Wave roadmap

1. **(this commit)** registry + FAMILIES.md; parser; ~45 matchers; genus
   routing; solvers; CLI; tests.
2. Transformations T1–T2; systems of equations (congruent number,
   simultaneous Pell); positivity/domain constraints (ℕ, Frobenius);
   S-unit solver wiring over ℚ (bases → S, via `solve_S_unit_equation`);
   Thue–Mahler runner (Gherga–Siksek code detection).
3. T3; genus-2 Igusa lookup (needs the Library's curve tables); Nagell–
   Ljunggren and recurrence matchers; Magma runner behind a feature flag;
   LMFDB label lookups (Mordell → Cremona label via the elliptic curve db).
4. T4 (L-function hash), database layer, the web front end (separate repo,
   LMFDB-based template per the talk).
