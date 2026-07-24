# The ζ-corrected tail theory: rigidity, the off-diagonal family, and the μ-ladder

Closes the boundary named by COMPOSITE_CHARTS.md §4 ("A0′=(2,0) breaks
the MODEL, not the chart"): what the forcing family becomes when the reduced
P-side element carries a pure-power defect x^ζ, and what that decides for
F12 and F13. Derivation: `zeta_tail.py`; checker: `zeta_tail_verify.py`
(**34/34, exit 0, --quiet**). Along the way this lane found a **new branch
structure at the UNCORRECTED family** that corrects one claim of
COMPOSITE_CHARTS.md and upgrades PHI_F7.md's judgment-6 observation to the
actual law (see §4).

## 1. Rigidity: the defect cannot live in the tail (one-slice theorem)

Model (COMPOSITE_CHARTS §4 commutator analysis, standard framework shape):
`P = x^ζ (C^a + tail)`, `Q = C^b + tail + F`, `v(F) < 0`, `[P,Q] = x^κ`,
`ell(C) = x^t c(y)`, κ = t−2 (fused-chart lemma). The head-head bracket's
top slice (the tails contribute nothing at this slice — that is the point;
terminology corrected per REVIEW_ZETA_MU.md item 1) is

```
[x^ζ C^a, C^b] -> ζ·b · x^(ζ+(a+b)t−1) · c^(a+b−1) c′   (nonzero for ζ>0),
```

it sits at κ + [ζ+1+(a+b−1)t] — strictly above the target — and nothing in
the model reaches it: any F-cross-term there would need
`v(F) = (a+b−j)t ≥ bt > 0`, contradicting `v(F) < 0`. So **a tail-carried
defect admits no solution at all**; the defect must be absorbed into the
algebraic element. (Verified symbolically at two truncations.)

## 2. The repair: D = x^η C and the off-diagonal family

GGV proportionality (st(P) ∥ st(Q)) distributes the defect evenly:
`D = x^η C`, η = ζ/a, restores commuting powers, and the forcing identity
`[D^a, x^s f c^(−b)] = x^κ` (s = κ+1−aT) yields the **ζ-corrected family**

```
a { T c f′ − [T(b−a) + κ + 1] c′ f } = c^(b−a+1),
T = t + η,     κ = t − 2   (chart-fixed: κ does NOT move with η).
```

This is exactly the standard two-parameter family evaluated **off the
κ = T−2 diagonal by η**. Fractional η lives in the refined chart x = wⁿ
(T_w = nT, κ_w = nκ+n−1, K_w = nK); the resonant degree d_res = K·a0/T is
chart-invariant. (All identities verified symbolically, including the
bracket chart rule [,]_w = n·w^(n−1)·[,]_x.)

**Universal η = −1 lemma:** T = t−1 ⟺ K = eT for every family; the ODE then
integrates in closed form, `(f/c^e)′ = 1/(aTc)`, so a polynomial f exists
iff **every residue of 1/c vanishes** — impossible in every root
configuration at dg ∈ {2,3} (any simple root has residue ≠ 0; tuning the
double-root residue to zero (β = 6α/5 at q=5) pushes the obstruction to the
remaining simple root; full-multiplicity roots give 15/α⁷ resp. −7/α⁸ ≠ 0).
**η = −1 is always dead.**

## 3. F12 and F13 sweeps

F12 (a=3, b=7, t=4, κ=2, a0=8, q=5, e=5, dg=3, gap=2, r=2):

| η | T_w | K_w | d_res | verdict |
|---|---|---|---|---|
| −1 | 3 | 15 | 40 | DEAD (K=eT log obstruction) |
| ±1/3, ±2/3, +1, +3 | — | — | non-integral | DEAD (min forced degree 34 > d_max 33) |
| 0 | 4 | 19 | 38 | standard family — the μ-ladder, §4 |
| **+2** | **6** | **27** | **36** | **CANONICAL COLLAPSE, §5** |
| −2, −3 | 2, 1 | 11, 7 | 44, 56 | arithmetically viable, **unmotivated — flagged OPEN** |

F13 (a=2, b=13, t=3, κ=1, a0=9, q=7, e=12, dg=2): **every** motivated defect
is dead (η=−1 by the lemma; ±1/2, ±1, +2 by degree count/non-integral
resonance); η=−2 viable-unmotivated; only η=0 survives (conditional; the
PHI_F7 dg=2 analysis applies). F13 j=1 is Orevkov's case.

**Verdict vs the directive's menu:** for every *motivated* nonzero defect the
answer is NO-POLYNOMIAL-SOLUTION — if the (8,24)/(9,21) reductions produce
any unit-size pure-power defect, F12/F13 carry **no standard-shape last
element at all** (they exit the tower class the corner law lives on). The
two model shapes that survive are η = 0 and (for F12 only) η = +2, both
law-consistent; selecting between them needs the actual polygon reduction
(judgment 5).

## 4. The μ-ladder (NEW; corrects COMPOSITE_CHARTS §5, realizes PHI_F7 judgment 6)

Indexing branches of the (un)corrected family by **μ = mult_{y+1}(g)**, the
F12 standard family (η=0, N=97, e+N=102) realizes **all three rungs**:

| μ | branch | signature |
|---|---|---|
| 1 | the μ=1 branch's two squarefree points (g₁ = (−1±√22)/8) | (814, 506, 102, 206) |
| 2 | g = (y+1)²(y−β), β a root of **195β⁴+120β³−40β²+32β−80** (2 real roots; system ≡ 0 exactly mod the quartic) | (814, 506, 203, 105) |
| 3 | g = (y+1)³, u = −(2048y⁴+2560y³+320y²−80y+35)/1155 | (814, 506, 304, 4) |

The μ=3 rung **is** the "ramified formulas (304, 4)" that COMPOSITE_CHARTS
§5 stated were *NOT realized* — that claim was an artifact of the
g^e-uniform ansatz, which cannot represent the ramified forced orders
(e−1)μ+1. **Correction:** at F12/η=0 the ramified signature IS realized,
alongside a genuinely new intermediate μ=2 rung.

The signatures obey the **μ-graded law**

```
deg = res + N·a0,   ord = ρ + N·q,
mult = μ(e+N) − (μ−1),   cof = gap + r(e+N) − (μ−1)(e+N−1),
```

which (via the identity **r = dg−1**, always) specializes *exactly* to the
old unramified law at μ=1 and to PHI_F7's ramified law (mult = dg(e+N)−(dg−1),
cof = gap+r) at μ=dg. PHI_F7's judgment-6 "mult_g-indexed unification,
observation only" is therefore now the actual branch structure, realized at
μ = 1, 2, 3 in one corner. **Parity refinement (scoped per
REVIEW_ZETA_MU.md item 3):** dg even kills μ=1 — PROVEN at dg=2 (PHI_F7's
eliminant enumeration) AND at dg=4 (MU_RUNGS_F10.md: all cubic-h partitions
eliminated by exact polynomials with zero real roots); even dg ≥ 6 remains
unenumerated (no gap>0,r>0 survey corner has dg ≥ 6). Odd dg keeps μ=1 available but does **not** exclude higher
rungs. Presentational notes from the same review: r = dg−1 is definitional
(r = a0−q−1, dg = a0−q), so the μ=dg specialization identity is exact
algebra, not an independent confirmation; and deg−ord = (e+N)dg+gap
identically, so the cofactor law is equivalent to the mult law — each
realized rung contributes ONE independent datum, not two.

## 5. F12 at η = +2: the canonical collapse

At T=6, K=27 the family collapses exactly like the seven standard corners:
f = A y²¹ g⁵ reduces the ODE to `9A(yg′ − 3g) = 1`, forcing g₂ = g₁ = 0,
g = y³+1 (gauge g(−1)=0), A = −1/27 — the forced-cyclotomic residual
H = y²−y+1 (same class as (75,125)/F2). With the off-diagonal N-formula
N = a[T(a+b)−(κ+1)]−2b = **157** and gap_eff = d_res − pure = **0** (the
defect shifts F12 into the clean gap=0 regime):

| μ | branch | signature |
|---|---|---|
| 1 | g = y³+1 (canonical) | (1292, 806, 162, 324) |
| 2 | — | **provably none** |
| 3 | g = (y+1)³, u = (8y²+4y−1)/27 | (1292, 806, 484, 2) |

Both rungs sit ON the μ-graded law. Note the μ-support differs between
η=0 ({1,2,3}) and η=+2 ({1,3}) — branch support is model-dependent.

## 6. What this means upstream

* **Corner law**: untouched at the twelve derived/audited points; the
  μ-graded law is a strict generalization that all of them obey (μ=1 for
  the unramified seven, μ=dg for the PHI_F7 four, μ=1 for (72,108)).
* **COMPOSITE_CHARTS.md §5**: the "NO ramified branch at dg=3" sentence and
  the "(814,506,304,4) NOT realized" parenthesis are corrected by §4 above;
  the conditional μ=1 signature itself stands.
* **Case compiler**: F12's conjectural flag should stay (model selection
  η ∈ {0, +2} open), but its reason can now cite this analysis: the surviving
  models and their full μ-rung signatures are enumerated.
* **PHI_F7 branch-completeness judgment**: dg=4 (F10) likely carries
  intermediate rungs μ=2,3 as well — named follow-up, one small solve each.

## `[judgment]` list

1. **[model form]** P = x^ζ(C-series), Q = C-series + F with v(F) < 0 is the
   COMPOSITE_CHARTS §4 commutator model (standard framework shape); the
   rigidity theorem is unconditional *within* it.
2. **[defect set]** The motivated candidates are |η| ≤ 1 (per-unit defect
   from A0′=(2,0) vs (1,0), either sign, fractional via ζ ∈ {±1,±2}); the
   scan extends to |η| ≤ 3. η = −2, −3 (F12) and η = −2 (F13) are
   arithmetically viable but correspond to no proposed reduction — flagged
   OPEN, nothing claimed. η = +2's canonical collapse is *discovered by
   scan*, post-hoc motivated by its unique restoration of the standard
   structure — model selection between η ∈ {0, +2} requires the actual
   (8,24) polygon reduction or a C-series build (no paper performs either).
3. **[chart carryover]** t = l and κ = t−2 for the escapes are the earlier
   fused-chart lemma + its judgment 3 (t = l unproven for the escapes).
4. **[branch completeness]** Partition analyses cover all root partitions of
   dg ≤ 3 with the (y+1)-place at a chosen rung; μ-rungs with the multiple
   root AWAY from −1 change no signature (mult is read at −1) and were not
   separately enumerated. dg=4 rungs (F10) not attempted.
5. **[N-formula]** The off-diagonal N = a[T(a+b)−(κ+1)]−2b extrapolates the
   published-corner bookkeeping into the D-model; same conditional status as
   every non-C-series-built point (PHI_CORNER4 judgment 3).

## Files

`zeta_tail.py` (derivation, all steps printed), `zeta_tail_verify.py`
(**34/34, exit 0, --quiet**; independent instances/truncations, fresh μ=2
system rebuild, exact mod-quartic residuals, explicit branch polynomials).
Nothing existing touched; run_tests.sh wiring left to the parent session.
