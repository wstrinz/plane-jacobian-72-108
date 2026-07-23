# Lean certificate: the f37 ideal-membership identity, kernel-checked

**Date:** 2026-07-23
**Project:** `lean_certificates/` (self-contained Lean 4, no mathlib)
**Toolchain:** `leanprover/lean4:v4.32.1`

## What is proven

The program's headline fact (see `d2_plane_72_108/F37_SATURATION_REPORT.md`) is
that the degree-31 elimination polynomial `f31` lies in the pre-resultant ideal:

```
f31 = c1·G1 + c2·G2 + c3·G3 + c4·(G5body + Φ)     (exact, over ℚ)
```

in `ℚ[d̃2, d̃1, d̃0, d₋₁, d₋₂, d₋₃, d₋₄, Φ]`. This is what makes the entire
`{f37 = 0}` branch of the resultant a classical excess-factor artifact.

**What this build actually does (precise wording).** It does **not** "formally
prove the ideal-membership theorem" in the abstract-algebra sense. It
**Lean kernel-checks the emitted integer certificate under a small custom
sparse-polynomial implementation** (`Cert/Poly.lean`): the kernel verifies that
a specific integer polynomial identity between opaque `List (Nat × Int)` data
normalizes to equality. The trust-base notes below (including the honest caveat
that `add`/`mul` are tested, not semantically proven) delimit exactly what that
buys. Scope is **characteristic zero**: the certificate divides out
`Dmul = 46875 = 3·5⁶`, so the `f31 = 0` conclusion holds over any `ℚ`-algebra /
field of characteristic `≠ 3, 5` (see the field-scope note in
`d2_plane_72_108/F37_SATURATION_REPORT.md`).

The Lean theorem `Cert.f37_certificate` proves the equivalent **integer** identity

```
Dmul · f31  =  ĉ1·G1 + ĉ2·G2 + ĉ3·G3 + ĉ4·G4        (exact, over ℤ)
```

where `Gᵢ` are the denominator-cleared generators (`G4 = G5body + Φ`, cleared),
`ĉᵢ = Dmul·cᵢ` are the denominator-cleared cofactors, and `Dmul = 46875` is the
nonzero integer lcm of all cofactor denominators. Dividing by the nonzero
integer `Dmul` recovers the ℚ-identity above (elementary; the kernel-checked
content is the integer identity).

Statement (verbatim):

```lean
def certCombo : Poly :=
  add (add (add (mul chat1 G1) (mul chat2 G2)) (mul chat3 G3)) (mul chat4 G4)

theorem f37_certificate : certCombo = Df31 := by
  native_decide
```

Sizes: cofactors `ĉ1..ĉ4` have 5639 / 5228 / 5848 / 3418 terms; the generators
have 4 / 4 / 4 / 8 terms; `Df31 = Dmul·f31` has 102 terms; 8 variables; packing
radix `R = 32` (> max exponent 25).

## Trust story

Three layers, from most to least trusted:

1. **Lean kernel** certifies the polynomial arithmetic. The polynomials are
   opaque `List (Nat × Int)` data; every `mul` and `add` is performed by the
   small verified library `Cert/Poly.lean` (sparse, packed-key, sorted; all
   definitions structurally recursive). The kernel checks that the assembled
   right-hand side normalizes to exactly `Df31`. Nothing about the *result* is
   taken on faith from the generator, Singular, or sympy.

2. **`native_decide`** (used for the headline theorem) compiles the verified
   library to native code and runs it, then the kernel trusts that evaluation.
   This extends the trusted base to the Lean compiler + this machine's runtime.
   A pure-kernel `decide` cross-check of the same machinery is provided in
   `Cert/KernelCheck.lean` (see there for the honest scope of what `decide`
   reduces without the compiler).

3. **The generator** (`d2_plane_72_108/lean_export/export_certificate.py`)
   produces the *data* only. It imports the exact data-loading functions of the
   READ-ONLY reference verifier `f37_sat_verify.py` — `pre_resultant_generators`
   (parsed from the canonical `generators.json`; no pickle on this path),
   `load_cofactors`
   (`f37_sat_certificate.txt`), `load_f31` (`f31_deg31.txt`) — and never
   hand-copies a coefficient. It clears denominators, **re-verifies the integer
   identity in sympy independently**, then emits the sorted `List (Nat × Int)`
   literals. If the generator emitted wrong data, the Lean kernel identity would
   simply fail to hold — the kernel is the backstop.

**Axiom audit** (`#print axioms`):

- `Cert.f37_certificate` depends on axioms: `[…native_decide.ax…]` — i.e. ONLY
  the `native_decide` reflection axiom; no `sorry`, no `Classical.choice`.
- `Cert.mul_assoc_G123` (a `decide` cross-check on real generator data) depends
  on **no axioms** — pure kernel, zero trust extension.

The provenance chain is: Singular `lift()` → `f37_sat_certificate.txt` →
(re-checked by `f37_sat_verify.py` in sympy) → `export_certificate.py` (clears
denominators, re-checks in sympy) → `Cert/Data.lean` → **Lean kernel** verifies
`certCombo = Df31`.

**One honest caveat.** The theorem states that two `Poly` normal forms are equal.
That this *means* the polynomial identity `Dmul·f31 = Σ ĉᵢ·Gᵢ` relies on the
`Cert.Poly` `add`/`mul` being faithful polynomial arithmetic. Those operations
are **not machine-proven** correct against a semantic model here (they are simple,
and are exercised by the `KernelCheck.lean` commutativity/associativity lemmas on
real data, but that is testing, not a proof). A fully rigorous version would
prove `add`/`mul` correct w.r.t. an evaluation homomorphism (or reuse mathlib's
`MvPolynomial`); see "what remains" below. This is the one place trust rests on
un-verified (though tested and very simple) code rather than the kernel.

## How to rebuild

```
# 1. regenerate the Lean data from the exact certificate sources (sympy):
cd d2_plane_72_108/lean_export && python export_certificate.py

# 2. build + kernel-check the theorem:
cd ../../lean_certificates && lake build
```

`lake build` succeeding is the gate: it means the kernel accepted
`f37_certificate` (and the `decide` cross-check in `KernelCheck.lean`).

## Build outcome

- `lake build`: **SUCCESS** — theorem `Cert.f37_certificate` accepted by the kernel.
- Toolchain `leanprover/lean4:v4.32.1`; machine had ~24 GB RAM free during the run.
- Build times (cold): `Cert.Poly` 29s, `Cert.Data` 172s (elaborating the
  ~20k-term literal data), `Cert.Certificate` (the `native_decide`) 27s,
  `Cert.KernelCheck` (pure-kernel `decide` cross-checks) ~3s. Full `lake build`
  ≈ **3m20s**. Rebuilds with the olean cache are seconds.

## What remains to formalize for the full (72,108) program

This certificate formalizes **one node** of the claim graph: the f37 branch is
resultant excess because `f31 ∈ ⟨G1,G2,G3,G5body+Φ⟩`. The remaining program
(per `STATE.md`) is the `f31` window/cascade argument. A full formalization
would still need:

- **Provenance of the generators.** Here `G1..G4` are *given* integer data
  (matching `STATE.md` item 4, shipped as exact term lists in `generators.json`,
  originally regenerated from `t4_state.pkl`). Formalizing
  their derivation from the original (72,108) system (the `(D̃²)` linear
  substitutions producing `(D̃³)₋₁,₋₂,₋₃` and `(D̃³)₋₅ + Φ`) is not done here.
- **The specialization step.** The ℚ[Φ]-identity specializes to the genuine
  instance `Φ = f₁·C₄²⁸ ∈ ℚ[y]` by any ring hom `ℚ[Φ] → ℚ[y]`; formalizing that
  membership is preserved under ring homomorphisms (trivial in mathlib, not done
  in this self-contained build).
- **The `f31` branch itself.** The saturation report closes only the f37 branch;
  the `f31` window/cascade program remains the open task and is untouched here.
- **A mathlib bridge (optional).** Relating this ad-hoc `Poly` to
  `MvPolynomial (Fin 8) ℤ` would let the identity feed standard algebra lemmas
  (ideal membership, `Ideal.mem`), at the cost of a mathlib dependency.
```
