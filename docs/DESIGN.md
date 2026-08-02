# Design

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
generalized Fermat, …) and a human-readable `transform` note when a
normalization was applied.

An equation matches many families; the registry's `parents` edges form a
specialization DAG and the *deepest* match wins. The lineage is reported so
the site can offer every applicable page (Pell ⊂ Pell-like ⊂ binary quadratic
forms ⊂ conics).

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

`solvers.solve` dispatches on the primary family, then up the lineage; wired
solvers use Sage/PARI (continued fractions, `qfbsolve`, `qfsolve`,
`two_squares`…, `BinaryQF.solve_integer`, `E.integral_points`, PARI `thue`
certified). Families without a definitive solver raise `SolverUnavailable`
carrying the registry's software pointers and *filled* code templates (e.g.
Magma `IntegralQuarticPoints`, Chabauty pipelines) — matching the site plan
that every equation homepage offers runnable code even where we don't run it.

Magma is optional everywhere: templates are emitted, nothing shells out in
wave 1. (Local Magma exists on this dev machine; a `magma` runner can be
added behind a feature check without changing the API.)

## Website backend contract

- `Classification.as_dict()` returns plain JSON types only (tested), so the
  future site (LMFDB-style, psycodict/Postgres) can store and render results
  directly.
- The registry YAML is the site's families table in embryo: slugs are stable
  identifiers, `priority` is the rollout order, `lmfdb` fields link to
  existing LMFDB collections (elliptic curves, genus 2 curves).
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
