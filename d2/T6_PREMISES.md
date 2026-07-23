# T6_PREMISES — the two outline-only premises, closed to citations (2026-07-22)

Scope: the derivation audit (`T6_SELECTION_AUDIT.md`, `verify_derivation.py`)
reduced AUDIT.md Risk 3 to exactly two premises that were left as outlines
(`T6_SELECTION_AUDIT.md` §4). This file states each precisely, gives the full
argument, verifies every finite piece (`t6_premises_verify.py`, 14 checks, all
pass), and delivers a writeup-readiness verdict.

Sources (file + line ranges; "proved in the paper" / "proved here" / "assumed"
are distinguished throughout):

- **GGHV22** = `paper_src/2204.14178.tex` (Guccione–Guccione–Horruitiner–Valqui,
  "…from 100 to 108"), the paper whose open case is ours.
- **GGV1** = `paper_src/1401.1784_GGV1.tex` (Guccione–Guccione–Valqui, "On the
  shape of possible counterexamples to the Jacobian Conjecture", J. Algebra 471
  (2017) 13–74; arXiv:1401.1784 v3). Fetched this session (the one permitted
  addition under `paper_src/`); bib entry GGHV22 lines 2094–2105.

Both premises are **the same machinery GGHV22 runs, in full and published
detail, for the *twin* closed case (9,27)** (its §4, `seccion 4`,
GGHV22 lines 1399–1596), transported from `t=3` to our `t=4`. The paper leaves
*our* case (the (8,28) reduction, GGHV22 lines 1000–1010) open only "for lack
of computing power" at the final elimination (GGHV22 lines 252, 267–268, 276) —
**not** at this setup stage. So the premises are not our invention; they are
GGV1 Props 1.13 and 2.1 applied exactly as the paper applies them, plus finite
arithmetic that changes only through the numerology `t: 3→4`
(`v_{1,0}(P): 6→8`, `v_{1,0}(Q): 9→12`, `v_{1,0}(C): 3→4`, `v(F): −4→−5`).

The two GGV1 propositions, verbatim from the fetched source:

- **GGV1 Prop 1.13** = `pr v de un conmutador`, GGV1 lines 480–491
  (proved in GGV1, one line, from the homogeneous-decomposition Remark
  `re v de un conmutador`, GGV1 lines 453–478):
  > for `P,Q ∈ L^{(l)}\{0}`, `(ρ,σ)∈𝔙`,
  > `v_{ρ,σ}([P,Q]) ≤ v_{ρ,σ}(P)+v_{ρ,σ}(Q)−(ρ+σ)`, with **equality iff
  > `[ℓ_{ρ,σ}(P), ℓ_{ρ,σ}(Q)] ≠ 0`**, and in that case the leading bracket is
  > `ℓ_{ρ,σ}([P,Q])`.

- **GGV1 Prop 2.1** = `P y Q alineados`, GGV1 lines 514–563 (proved in GGV1;
  part (2a) is Joseph's Prop 2.1(2), the rest is elementary):
  > for `(ρ,σ)`-homogeneous `P,Q`, set `τ=v_{ρ,σ}(P)`, `μ=v_{ρ,σ}(Q)`. If
  > `[P,Q]=0` and `(μ,τ)≠(0,0)`, and `m,n` are coprime with `nτ=mμ`, then there
  > exist `R∈L^{(l)}` and `λ_P,λ_Q∈K^×` with **`P=λ_P R^m`, `Q=λ_Q R^n`**
  > (and `R∈L` when `P,Q∈L`).

---

## Premise 1 — the leading-form normalization `ℓ(P)=R²`, `ℓ(Q)=R³`, `C₄=y⁷(y+1)`

### 1.1 What the downstream derivation needs (precise statement)

For a subcase-(2) candidate — `P,Q∈K[x,y]`, `[P,Q]=x²`, with the Newton
polygons of GGHV22 Prop "Case (8,28)"(2) (GGHV22 line 1004),
`N(P)⊃{(8,14),(8,16)}`, `N(Q)⊃{(12,21),(12,24)}` as its `(1,0)`-leading
corners — there is a `(1,0)`-homogeneous `R` with

    ℓ_{1,0}(P) = R²,   ℓ_{1,0}(Q) = R³,   R = x⁴ C₄,   C₄ = y⁷(a₀+a₁y), a₀a₁≠0,

and, after a linear change of variables, `C₄ = y⁷(y+1)`. This is exactly what
`Proposition calculo de C` (its `t=4` analogue) consumes to launch the series
`C = x⁴C₄ + x³C₃ + …`, `C²=P`, `ℓ_{1,0}(C)=x⁴C₄` — the ansatz that
`verify_derivation.py` §C/§D takes as its starting point.

### 1.2 The argument

**(a) `[ℓP, ℓQ] = 0`.** Take `(ρ,σ)=(1,0)`, so `ρ+σ=1` and `v_{1,0}` is the
`x`-exponent. From the leading corners, `v_{1,0}(P)=8`, `v_{1,0}(Q)=12`
(*read off the polygons; proved here — the polygons themselves are GGHV22's
Prop "Case (8,28)", lines 1000–1010, proved in the paper*). Since `[P,Q]=x²`,
`v_{1,0}([P,Q])=2`. Then

    v_{1,0}([P,Q]) = 2  <  8 + 12 − 1 = 19 = v_{1,0}(P)+v_{1,0}(Q)−(ρ+σ),

a **strict** inequality, so by **GGV1 Prop 1.13** the equality case fails and
`[ℓ_{1,0}(P), ℓ_{1,0}(Q)] = 0`. *(Verified: `t6_premises_verify.py` P1a.)*

**(b) Common power `R`.** `ℓ_{1,0}(P)`, `ℓ_{1,0}(Q)` are `(1,0)`-homogeneous
with `[ℓP,ℓQ]=0`, `τ=8`, `μ=12`, `(μ,τ)≠(0,0)`. The coprime `(m,n)` with
`nτ=mμ` is `(m,n)=(2,3)` (since `3·8 = 2·12 = 24`). By **GGV1 Prop 2.1**,
there is a `(1,0)`-homogeneous `R` with `ℓ(P)=λ_P R²`, `ℓ(Q)=λ_Q R³`.
*(Verified: P1b.)*

**(c) `R = x⁴C₄`, `C₄=y⁷(a₀+a₁y)`, `a₀a₁≠0`.** `ℓ_{1,0}(P)` sits at
`x`-exponent 8, so `R²` sits at `x`-exponent 8 and `R=x⁴·(poly in y)`; write
`R=x⁴C₄`. Then `R²=x⁸C₄²` must equal `ℓ(P)=x⁸·(edge polynomial of N(P))`, whose
`y`-support is the segment `[14,16]` (corners `(8,14),(8,16)`); so `C₄²` has
`y`-degrees `14..16`, forcing `C₄` to have `y`-degrees `7..8`,
`C₄ = y⁷(a₀+a₁y)`. Both corners being present forces `a₀≠0` (else `ord_y C₄²`
would be 16, killing the `(8,14)` corner) and `a₁≠0` (else `deg_y C₄²` would be
14, killing the `(8,16)` corner). Consistently, `R³=x¹²C₄³` has `y`-support
`[21,24]`, matching `N(Q)`'s corners `(12,21),(12,24)`. *(Verified: P1c —
`R²` and `R³` reproduce all four corners, and the two degenerate collapses.)*
This step is *proved here* from the paper's polygons.

**(d) `R` is primitive.** For GGV1 Prop 2.1 to yield a *single* power (rather
than a Laurent polynomial in a smaller `R`), and for the strip range in Premise
2 to be `k<3`, one needs `R=x⁴C₄=x⁴y⁷(y+1)` to be no proper power `S^d`,
`d≥2`, of any element of `K[x,y]`. A `d`-th power requires `d` to divide the
multiplicity of every irreducible factor: `x` (mult 4), `y` (mult 7), `(y+1)`
(mult 1). Since `gcd(4,7,1)=1`, only `d=1`: `R` is primitive. *(Verified: P1d.
This is the `t=4` analogue of GGHV22's clause "`x³C₃` is not the positive power
of any element of `K[x,y]`", GGHV22 line 1513, `t=3`.)*

**(e) Normalization to `C₄=y⁷(y+1)`.** The substitution `y = (a₀/a₁)·ỹ` sends
`y⁷(a₀+a₁y)` to `(const)·ỹ⁷(ỹ+1)`; the overall constant and the leading
coefficient are absorbed into the definition of `R` (i.e. into `C₄`'s
normalization) and, together with a compensating scaling of `x`, into keeping
the bracket a scalar multiple of `x²`. *(Verified: P1e — the change of
variables produces exactly `ỹ⁷(ỹ+1)` up to a nonzero constant.)* This is the
`t=4` copy of GGHV22's clause "by a linear change of variables we may also
assume `C₃=y⁸(y+1)`", GGHV22 line 1414.

### 1.3 One honest caveat (attributed to the paper, not a new gap)

GGV1 Prop 2.1 supplies `ℓ(P)=λ_P R²`, `ℓ(Q)=λ_Q R³` with *two* scalars.
Absorbing `λ_P` into `R` leaves one residual scalar `λ' = λ_Q/λ_P^{3/2}`
(equivalently the constant `α` in the alignment identity `(ℓP)³ = α(ℓQ)²`);
the clean statement `ℓ(Q)=R³` is the case `λ'=1`. The downstream argument does
use `λ'=1`: Premise 2 opens with `ℓ(Q)=ℓ(C³)`, hence `v(Q−C³)<v(C³)=12`
(GGHV22 line 1508), which is exactly `λ'=1`. **GGHV22 folds this scalar
normalization into the same "there exists `R` with `ℓ(P)=R²`, `ℓ(Q)=R³`"
sentence and the "linear change of variables" clause (GGHV22 lines 1411–1414)
and uses it identically in its *published, closed* (9,27) proof** (and again
for the (9,24) case with base `x³y`, GGHV22 lines 1918–1927). It is a scalar
gauge-fix, not a mathematical gap; we rely on it exactly as the paper does for
its closed twin.

### 1.4 Verdict — Premise 1: **READY-WITH-CITATION**

Relies on GGV1 Prop 1.13 (GGV1 lines 480–491) and GGV1 Prop 2.1 (GGV1 lines
514–563) *used as stated*, plus finite arithmetic all verified in
`t6_premises_verify.py` (P1a–P1e). The transport `t=3→t=4` changes only the
valuations `6→8`, `9→12`; the structural inputs and the primitivity/normalization
clauses are the paper's own, applied to its own open case. The residual scalar
normalization (§1.3) is the paper's convention, identical to its closed-case
proof. No step is a vague "should follow."

---

## Premise 2 — the α-strip WLOG: `Q = C³ + λC⁻¹ + F`, `v_{1,0}(F) = −5`

### 2.1 What the downstream derivation needs (precise statement)

With `C` as in Premise 1 (`C²=P`, `ℓ_{1,0}(C)=x⁴C₄`), there exist
`α₂,α₁,α₀,α₋₁∈K` such that
`F := Q − C³ − α₂C² − α₁C − α₀ − α₋₁C⁻¹` has `v_{1,0}(F) = −5`; and after the
replacement `Q ↦ Q̃ = Q − α₂P − α₀`, `P ↦ P̃ = P + ⅔α₁` (which preserves
`[P,Q]=x²` and the polygons), one may **assume `α₂=α₁=α₀=0`**, i.e.

    Q = C³ + λ C⁻¹ + F,   λ∈K,   v_{1,0}(F) = −5.

This is precisely the input to `verify_derivation.py` §A (the forcing ODE for
`f₁=C₄³F₋₅`) and §B (λ enters only the dropped `j=4` slice).

### 2.2 The argument (= GGHV22 lines 1508–1546, transported `t=3→t=4`)

**(a) The strip terminates at `v(F)=−5`.** Because `ℓ(Q)=ℓ(C³)` (Premise 1,
§1.3), `v_{1,0}(Q−C³) < v_{1,0}(C³)=12`. Each strippable leading form is a pure
power `α_k R^k` (`R=x⁴C₄=ℓ(C)`), with `v_{1,0}=4k`. GGV1 Prop 1.13 + Prop 2.1
give it: whenever the current remainder `G` has
`v(G)+v(P)−1 > v(x²)=2`, i.e. `v(G) > −5`, the strict inequality forces
`[ℓP,ℓG]=0` (Prop 1.13) and hence `ℓ(G)=α_k R^k` (Prop 2.1, `R` primitive),
which we cancel by subtracting `α_k C^k` (since `ℓ(C^k)=R^k`), lowering `v`.
The admissible powers are `k` with `−4 ≤ 4k < 12`, i.e. `k∈{−1,0,1,2}`
(the range `−2<k<3` of GGHV22 line 1515), giving the four terms
`α₂C² + α₁C + α₀ + α₋₁C⁻¹`. The descent **halts** at `F` with `v(F)=−5`,
because `−5` is not a multiple of `v_{1,0}(R)=4`, so `ℓ(F)` cannot be any
`α_k R^k` and `[ℓP,ℓF]≠0`. *(Verified: P2c — the four powers `−4,0,4,8` lie in
`[−4,12)`, and `−5 ≢ 0 mod 4`.)* The construction of the `α_k` that makes
`[ℓP, ℓ(Q−C³−Σα_kC^k)]≠0` is GGHV22's citation to "the arguments of GGV3
Section 1" (GGHV22 line 1517); it is the elementary top-down cancellation just
described (*glue step filled here* — the substitution is now VERIFIED against
the source: `paper_src/1406.0886_GGV3.tex`, fetched 2026-07-22; GGV3's own
strip runs inside its Thm `principal` (1.8) proof via van den Essen Lemma
10.2.11 in the (1,1)-grading, equivalent to the Props 1.13/2.1 single-step
descent used here — see `GGV3_CITATION_CHECK.md`, verdict SUBSTITUTION
VERIFIED, no mathematical gap).

**(b) `v(F)=−5` exactly.** Each stripped term Poisson-commutes with `P=C²`:
`α₂C²=α₂P`, and `[C^k,C²]=0` for every `k` because powers of `C` are functions
of `C` (*verified here: P2a, `[C^k,C²]=0` for `k∈{3,1,0,−1}` with a generic
`C`*). Hence `[F,P] = [Q,P] = −x²`, so `v_{1,0}([F,P]) = 2`. Since
`[ℓF,ℓP]≠0`, GGV1 Prop 1.13's **equality** case gives
`2 = v(F) + v(P) − 1 = v(F) + 8 − 1`, i.e. `v(F) = −5`. *(Verified: P2b.)*
This is the `t=4` analogue of GGHV22 lines 1521–1525 (`t=3`: `v(F)=−4` from
`v(P)=6`).

**(c) The Remark WLOG (`α₂=α₁=α₀=0`).** Put `P̃ = P + ⅔α₁`, `Q̃ = Q − α₂P − α₀`.

- *Bracket and polygons preserved.* `[P̃,Q̃] = [P + ⅔α₁, Q − α₂P − α₀] =
  [P,Q] = x²` (constants bracket to 0; `[P,α₂P]=0`). *(Verified: P2e.)* Adding
  a constant to `P` and a scalar multiple of `P` plus a constant to `Q` leaves
  the top corners — hence the polygons and Premise 1 — untouched.
- *The `α₁C` term is absorbed.* `C̃ := √P̃ = √(C² + ⅔α₁) =
  C + ⅓α₁C⁻¹ + O(C⁻³)`, so `C̃³ = C³ + α₁C + O(C⁻¹)`. Thus
  `Q̃ = C³ + α₁C + λC⁻¹ + F = C̃³ + λC̃⁻¹ + F̃` for some `λ∈K`, the `α₁C` term
  now inside `C̃³` and the `C̃⁻¹` vs `C⁻¹` discrepancy absorbed into `F̃`
  (`v` unchanged). *(Verified: P2d — expanding `√(C²+⅔α₁)` gives coefficient
  `⅓α₁` on `C⁻¹`, and its cube gives coefficient exactly `1` on `C³` and
  `α₁` on `C`.)* The constant `⅔` is `t`-independent (it comes from the cube),
  so this Remark is *verbatim* the paper's (GGHV22 lines 1528–1546).

After the replacement we rename `(P̃,Q̃,C̃,F̃)→(P,Q,C,F)` and have
`Q = C³ + λC⁻¹ + F`, `v(F)=−5`, `α₂=α₁=α₀=0`, as required.

### 2.3 Verdict — Premise 2: **READY-WITH-CITATION**

The entire strip is GGHV22 lines 1508–1546 with the valuations shifted
`t=3→t=4` (`v(P):6→8`, `v(F):−4→−5`); its only structural inputs are GGV1
Props 1.13 and 2.1, already cited and used in Premise 1. Every finite piece —
the terminating `k`-range, `v(F)=−5`, the commutation `[C^k,C²]=0`, the `⅔α₁`
Remark algebra, and bracket preservation — is verified in
`t6_premises_verify.py` (P2a–P2e). The one prose citation ("GGV3 Section 1")
is an elementary descent filled here from Props 1.13/2.1. No vague step remains.

---

## Summary

| Premise | Statement | Verdict |
|---|---|---|
| 1 | `ℓ(P)=R²`, `ℓ(Q)=R³`, `R=x⁴C₄`, `C₄=y⁷(y+1)` | **READY-WITH-CITATION** (GGV1 Props 1.13, 2.1) |
| 2 | α-strip WLOG: `Q=C³+λC⁻¹+F`, `v(F)=−5` | **READY-WITH-CITATION** (GGHV22 §4 template) |

Both premises are the published GGV1/GGHV22 machinery applied to our open case,
with the only change being the numerology `t=3→t=4`. Neither is a mathematical
gap; both are "published proposition/argument used as stated" plus finite
arithmetic, and the finite arithmetic is verified exactly
(`t6_premises_verify.py`, 14 checks). The single honest attribution (§1.3, the
scalar `λ_Q=1` normalization) is the paper's own convention, used identically in
its peer-reviewed closed twin case (9,27).

Files this session: `T6_PREMISES.md` (this file), `t6_premises_verify.py`
(14 checks, all pass), `paper_src/1401.1784_GGV1.tex` (fetched GGV1 source).
