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


---

## 3. Cubic equations and genus one


---

## 4. Curves of higher genus and binary forms


---

## 5. Fermat-type equations


---

## 6. Polynomial–exponential equations


---

## 7. Norm forms, index forms, unit equations


---

## 8. Surfaces, higher dimension, and thin orbits


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
