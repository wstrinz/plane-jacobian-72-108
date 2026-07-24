# PRIOR_ART.md — novelty audit of the corner law vs. the published literature

2026-07-23. Deep read (full TeX where available), not keyword search. Executable
companion: `prior_art_postdiction_verify.py` (12 checks, in run_tests.sh).

## What was read

Full text: GGV1 (arXiv:1401.1784, bracket-ODE proposition L985–1114), GGV3
(1406.0886, §5 the (50,75) kill), GGV5 (1708.07936, chains / last-lower-corner),
GGHV22 (2204.14178, §Case (8,28) L1000–1311, §4 "(9,24) and (9,27)" L1399–1812,
§5 "(7,21)"), 1310.8249 (*A differential equation for polynomials related to the
JC*, all 999 lines), 2506.05697 (Ramirez–Valqui) main sections; 1409.6390,
1605.09430, 1708.09367, 1111.6100/1205.6827 at theorem level with targeted greps.
Abstract level: 2402.11135, Magnus-formula series I–IV (IV queried), Moh 1983
(scans only — see caveats), Horruitiner PUCP thesis 2018 (unfetchable; its arXiv
siblings 1708.07936/1708.09367 were read), Appelgate–Onishi 1985 (secondary).

## The decisive find: GGHV22 §4 solves ONE instance of the forcing ODE

At the (66,99)/(9,24) corner (C₃ = y⁸(y+1), a₀=9, q=8, r=0), GGHV22 L1571–1597
derives `6C₃f₁′ − 10C₃′f₁ = C₃²`, states uniqueness ("found using a CAS"), and
prints `f₁ = −(1/910)·y⁹(y+1)²·(243y⁴−81y³+54y²−42y+35)`, using the quartic's
separability/coprimality in the kill. This is exactly the parametric ODE at
(a,b,t,κ)=(2,3,3,1) — κ=t−2 holds — and the f-level corner law postdicts its
full signature with zero fitting freedom: ord=(e−1)q+1=9, mult=e=2,
deg=(e·a₀−q+1)+gap=15, unit cofactor degree gap+r=4, separable, coprime to
y(y+1). The published record contains its own out-of-sample confirmation of the
law (machine-checked in `prior_art_postdiction_verify.py`). GGHV22 L1612–1632
also publishes the D_k := C_kC₃^{5−2k} clearing transformation — the ancestor of
our D-transformation. None of this is generalized there: one CAS-solved case
inside one kill proof; no parametric family, no Φ, no N, no signature law.

## Verdicts per claim

1. **Closed-form Φ = f·C^N with explicit signatures: NOT FOUND — new.**
   No paper states Φ, any Φ closed form, or any y^A(y^k+1)^B formula for a
   tower's last element. Cite GGHV22 §4 for the one published closed-form
   last-element instance (f-level).
2. **Signature law (deg/ord/mult/cofactor; N = a[t(a+b−1)+1]−2b): NOT FOUND,
   even partially — new.** Nearest relatives, all distinct: GGV1's
   equal-multiplicity/separability proposition (cite as ancestor of the
   multiplicity bookkeeping); 1708.09367's intersection-number formulas
   (different quantity); GGV5's last-lower-corner lattice bounds (constraints,
   not signatures). No equivalent of N anywhere.
3. **κ=t−2 / parametric ODE / uniqueness / ramified-unramified dichotomy:
   PARTIALLY — instance layer published, law new.** Published: the
   bracket-to-ODE mechanism in general form (GGV1 eqq1); the solved instance
   with uniqueness (GGHV22 §4); κ=l−2 as computed instances (GGHV22 L1229 at
   l=4; 1310.8249 ψ₃ at l=3); ODE-reduction genre precedent (1310.8249 reduces
   B=16 to an unsolved Abel ODE). NOT published: the parametric (a,b,t,κ)
   family with general closed-form solution structure (forced cyclotomic
   residuals), general uniqueness, κ=t−2 as a structural theorem, and the
   ramified formula dg(e+N)−(dg−1) (no trace anywhere).
4. **Resonance gap (r=0 ⟹ unit cofactor; the (72,108) quartic, disc class 17):
   PARTIALLY — one instance in print (GGHV22 §4's quartic), no concept of gap,
   no r=0 statement, no cofactor formula, nothing on discriminant classes. The
   (72,108) case itself is explicitly left open in print (L268) with no ODE or
   quartic published.**

## Required citations for any writeup

- **GGHV22 (2204.14178) §4** — solved ODE instance, unique closed-form f₁, the
  quartic-cofactor instance, D_k transformation. The single most important
  citation; omitting it would be a real gap.
- **GGV1 (1401.1784)** Prop. L985ff — bracket ODE; equal-multiplicity/
  separability lemma.
- **1310.8249** — ODE-reduction genre precedent (Abel ODE, unsolved).
- **GGV3 (1406.0886) / 2506.05697** — C-expansion polynomial systems.
- Optionally the Magnus series (2201.06613 etc.) as parallel approximate-root
  framework.

## Caveats

- **Moh 1983** (Crelle 340, pp. 140–212): GDZ page scans only, not
  machine-readable here. Risk low (GGV papers cite him for case lists and
  reduction technique, and he had no chain/N formalism), but a per-case
  computation resembling a last element cannot be strictly excluded without
  reading the scans.
- **Horruitiner PUCP thesis 2018**: unfetchable repository handle; almost
  certainly matches 1708.07936/1708.09367 (read). Obtain if it surfaces.
- Two junk search hits (spam PDFs) discarded — not prior art.

## Bottom line

The framework is theirs (and cited); the laws are ours. Strongest
presentational move, supported by the record: the corner law *postdicts*
GGHV22's printed f₁ exactly — coefficients, signature, quartic — with zero
fitting freedom.
