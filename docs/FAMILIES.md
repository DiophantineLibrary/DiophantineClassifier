# Families of Diophantine equations

A prioritized enumeration of families for the Diophantine Library, and the source of
truth behind the machine-readable registry
([`diophantine_classifier/data/families/`](../diophantine_classifier/data/families),
one YAML file per family) that drives the classifier. References given as
short strings here appear with full bibliographic data — and per-family
annotations of *why* each reference matters — in
[`references.bib`](../diophantine_classifier/data/references.bib) and the
per-family files.

**What counts as a family.** A parametrized collection of equations that is treated
uniformly in the literature: it has a standard form, a body of theory describing its
solution sets (finiteness, structure, effectivity), and ideally an algorithm or at
least a named method. Families form a DAG under specialization (Pell ⊂ Pell-like ⊂
binary quadratic forms ⊂ conics), and the classifier reports the *most specific*
family containing a given equation, together with its ancestors.

**Fields.** Each entry records: the standard form (unknowns are Latin letters near
the end of the alphabet or as displayed; parameters are the remaining letters), the
solvability status, the main methods, software that solves or partially solves the
family, and one or two key references. The YAML registry mirrors these fields.

**Status vocabulary.**

- `solved` — complete description of all solutions (possibly a parametrization).
- `algorithmic` — a terminating algorithm exists and is implemented in standard software.
- `effective` — finiteness with effectively computable bounds (usually Baker's method);
  practical algorithms exist for moderate parameters but may need per-instance work.
- `ineffective` — finiteness known (Thue–Siegel–Roth, subspace theorem, Faltings) with
  no general algorithm; individual instances handled by descent/Chabauty/modularity.
- `partial` — solved for sub-ranges of parameters or conditionally.
- `open` — status of the solution set is a named open problem.

**Priorities.** `P1` = launch set for the Diophantine Library (well-understood theory,
existing software, high query volume); `P2` = second wave (needs more infrastructure
or number-field input); `P3` = aspirational / research-level exhibits.

Solutions are sought in ℤ (or ℕ where stated) throughout; rational-point versions are
noted where they differ. Number-field and S-integer generalizations exist for nearly
every family and are deferred to per-family notes, matching the site plan (first ℤ
and ℚ, eventually number fields and S-integers).

---

## 1. Linear equations and lattice problems

### `linear` — Linear Diophantine equation — P1, solved
**Form.** a₁x₁ + ⋯ + a_k x_k = b.
**Status.** Solvable iff gcd(a₁,…,a_k) ∣ b; solution set is a coset of a rank-(k−1)
lattice, computed by the extended Euclidean algorithm / Hermite normal form.
**Software.** Sage: `xgcd`, `matrix(ZZ,...).hermite_form()`, `solve_diophantine`;
PARI: `matsolvemod`; every CAS.
**References.** Classical (Bachet 1621, Euler). See Niven–Zuckerman–Montgomery ch. 5.


---

## 2. Quadratic equations

### `quadric` — General quadratic Diophantine equation — P2, algorithmic
**Form.** Q(x₁,…,x_k) + L(x₁,…,x_k) + c = 0 (arbitrary quadratic, k ≥ 3).
**Status.** Decidable in general — the deepest case of the quadratic theory
(Grunewald–Segal give a general algorithm for orbit problems covering all
quadratics); in practice: complete the square and reduce to
representation/isotropy plus congruences.
**Software.** case-by-case reduction in Sage/PARI (`qfsolve` after completion);
sympy `diophantine` handles some shapes.
**References.** Siegel 1972; Grunewald–Segal, "How to solve a quadratic equation
in integers" (Math. Proc. Camb. Phil. Soc. 89, 1981).


---

## 3. Cubic equations and genus one

### `cubic-surface` — Cubic surfaces / del Pezzo — P3, research
**Form.** F(x, y, z, w) = 0 cubic (e.g. diagonal ax³+by³+cz³+dw³ = 0).
**Status.** Rational points conjecturally dense once one exists (unirationality);
Brauer–Manin obstructions; decidability unknown. Research-level exhibits
(e.g. x³+y³+z³ = n sits here as an affine slice family).
**References.** Colliot-Thélène–Kanevsky–Sansuc; Poonen, *Rational Points on Varieties*.


---

## 4. Curves of higher genus and binary forms

### `general-curve` — Integral/rational points on a general curve — P1 (as fallback), ineffective
**Form.** C(x, y) = 0 irreducible, genus g.
**Status.** g = 0: reducible to conics/parametrization (integral points via
Pell-type analysis of the parametrization); g = 1: effective for integral points,
rank machinery for rational; g ≥ 2: rational points finite (Faltings 1983,
ineffective — the central open effectivity problem), integral points effective
only in special shapes (hyper/superelliptic); Chabauty–Coleman when rank < g;
in higher genus the Library plans isomorphism identification via L-function hashes.
**Software.** genus: Sage `Curve(...).genus()`; points: `Curve.rational_points`,
PARI `hyperellratpoints`, Magma `Chabauty`, `PointSearch`.
**References.** Siegel 1929; Faltings 1983; Bombieri–Gubler, *Heights in Diophantine
Geometry*; Stoll, "Rational points on curves" (survey, 2011).


---

## 5. Fermat-type equations

### `generalized-fermat` — Generalized Fermat equation — P1, partial (the frontier)
**Form.** a·x^p + b·y^q = c·z^r with gcd(x,y,z) = 1 (primitive), signature (p,q,r).
**Status.** Governed by χ = 1/p + 1/q + 1/r:
χ > 1 (spherical): solutions form finitely many polynomial parametrizations
(Beukers 1998; explicit lists by Edwards 2004 for e.g. (2,3,5));
χ = 1 (euclidean): reduces to elliptic curves, completely understood classically;
χ < 1 (hyperbolic): finitely many primitive solutions for each signature
(Darmon–Granville 1995, via Faltings — ineffective); resolved signatures include
(n,n,n) (Wiles), (2,3,7) (Poonen–Schaefer–Stoll 2007), (2,3,8), (2,3,9), (2,4,n),
(3,3,n) partial, and a growing list via the modular method; ten known sporadic
primitive solutions with min ≥ 2 exponents; Fermat–Catalan/Beal conjectures open
($1M Beal prize for the coefficient-free version with all exponents ≥ 3).
**Software.** No general solver — per-signature Frey-curve computations
(Magma modular forms machinery); parametrizations for χ > 1 implementable.
**References.** Darmon–Granville 1995; Beukers 1998; Edwards 2004;
Poonen–Schaefer–Stoll 2007; Bennett–Chen–Dahmen–Yazdani, "Generalized Fermat
equations: a miscellany" (Int. J. Number Theory 11, 2015).


---

## 6. Polynomial–exponential equations

### `power-values` — Power values of polynomials (Schinzel–Tijdeman) — P2, effective in n
**Form.** f(x) = c·yⁿ, f fixed polynomial with ≥ 2 distinct roots, n ≥ 2 unknown.
**Status.** n is effectively bounded (Schinzel–Tijdeman 1976); for each fixed n it
is a superelliptic equation (effective). Includes `lebesgue-nagell` and, e.g.,
perfect powers among products of consecutive integers (Erdős–Selfridge: never a
perfect power).
**References.** Schinzel–Tijdeman 1976; Erdős–Selfridge 1975; Shorey–Tijdeman,
*Exponential Diophantine Equations* (1986) — the standard reference for this whole section.

### `exponential-diophantine` — Purely exponential equations — P2, effective (few terms)
**Form.** c₁·b₁^{n₁} + ⋯ + c_k·b_k^{n_k} = c (fixed bases, unknown exponents);
e.g. 2ᵃ + 3ᵇ = 5ᶜ, Jeśmanowicz conjecture instances. An unknown base raised to an
unknown exponent (Catalan, Goormaghtigh, Nagell–Ljunggren) is polynomial-exponential,
not purely exponential.
**Status.** Three-term cases: S-unit machinery gives finiteness and effectivity;
many-term cases hit the subspace theorem (ineffective). General exponential
Diophantine solvability is **undecidable** (Davis–Putnam–Robinson 1961) — the other
side of the decidability boundary.
**References.** Shorey–Tijdeman 1986; Evertse–Schlickewei–Schmidt 2002 (subspace);
Davis–Putnam–Robinson 1961.


---

## 7. Norm forms, index forms, unit equations

### `decomposable-form` — Decomposable form equations — P2, umbrella
**Form.** F(x₁,…,x_k) = m, F a product of linear forms over ℚ̄ (includes binary
forms, norm forms, index forms, discriminant forms).
**Status.** The unifying framework (Evertse–Győry): finiteness ↔ nondegeneracy
conditions via S-unit theory; effective in the norm-form/Thue cases of rank
conditions, ineffective in general (subspace theorem).
**References.** Evertse–Győry, *Discriminant Equations in Diophantine Number Theory*
(2017) and *Unit Equations* (2015).


---

## 8. Surfaces, higher dimension, and thin orbits

### `egyptian-fractions` — Unit fraction equations — P1 (Erdős–Straus), open
**Form.** 1/x₁ + ⋯ + 1/x_k = a/n.
**Status.** For fixed k, a, n: finitely many, enumerable by branch-and-bound.
As n varies: `erdos-straus` below; rich combinatorics (Erdős–Graham problems,
Bloom 2021 density result for unit fractions).
**References.** Graham's surveys; Bloom 2021.

### `diagonal-form` — Diagonal equations — P2, umbrella
**Form.** a₁x₁^k + ⋯ + a_sx_s^k = c.
**Status.** Umbrella for `waring`, `sum-of-three-cubes`, `equal-sums-like-powers`,
Fermat quartic surfaces (e.g. x⁴+y⁴+z⁴ = w⁴: Euler's conjecture, disproved by
Elkies 1988 — elliptic fibration method; minimal solution Frye);
local solvability decidable, global behavior varies wildly with (k, s).
**References.** Davenport–Lewis 1963; Elkies 1988.


---

## 9. The undecidability boundary


---

## 10. General references

- H. Cohen, *Number Theory, Vol. I: Tools and Diophantine Equations*, GTM 239 (2007) —
  the algorithmic compendium closest in spirit to this project.
- N. Smart, *The Algorithmic Resolution of Diophantine Equations*, LMS Student Texts 41 (1998).
- T. N. Shorey, R. Tijdeman, *Exponential Diophantine Equations* (1986).
- J.-H. Evertse, K. Győry, *Unit Equations in Diophantine Number Theory* (2015);
  *Discriminant Equations in Diophantine Number Theory* (2017).
- L. J. Mordell, *Diophantine Equations* (1969).
- L. E. Dickson, *History of the Theory of Numbers*, vol. II: Diophantine Analysis (1920).
- E. Bombieri, W. Gubler, *Heights in Diophantine Geometry* (2006).
- B. Poonen, *Rational Points on Varieties*, GSM 186 (2017).
- Y. Bilu, Y. Bugeaud, M. Mignotte, *The Problem of Catalan* (2014).
- B. M. M. de Weger, *Algorithms for Diophantine Equations* (1989).
