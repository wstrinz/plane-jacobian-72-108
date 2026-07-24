# R9_VALSPLIT — the divisibility-forced valuation split (deg_e = 10 wall)

**Status: COMPLETE (heavily truncated by an external outage) — verdict:
WALL SURVIVES the valuation split at both triage tiers.**
Files: `r9_valsplit.py` (case enumeration + checked exhaustiveness + sweep
driver), `r9_valsplit_results.json` (machine record), logs
`r9_valsplit_run.log`, `r9_valsplit_control.log`.

The named continuation after the dm4 elimination (`R9_SYMBOLIC.md` §3c):
Lane B proved the deg_e = 10 wall survives dm4 elimination and localized it
to the dm2/dm3 spare ansätze; this lane imposes the certified divisibility
lemma's full valuation content, splitting each state into finitely many
strictly smaller cases.

## 1. The split, and why it is sound

On the G-variety, `monic(e) | dm2·dm3` in Q[y] (certified identity,
`r9_symbolic_elim.py`).  With monic(e) = (y+1)⁹(y−r) for the R9 column and
(y+1)^a for the batch deg_e=10 T2 states, prime-factor valuation additivity
gives

    v_{y+1}(dm2) + v_{y+1}(dm3) ≥ 9,      v_{y−r}(dm2) + v_{y−r}(dm3) ≥ 1,

so every variety point (including dm2 ≡ 0 or dm3 ≡ 0, representable in any
case with A ≡ 0 / B ≡ 0) lies in some case (i, j), i = 0..9, j = 0..1:

    dm2 = (y+1)^i (y−r)^j A,        deg A ≤ 12 − i − j
    dm3 = (y+1)^{9−i} (y−r)^{1−j} B, deg B ≤ 14 − (9−i) − (1−j)

— **20 cases, 18 spare unknowns each** (down from 28).  Batch analogue:
i = 0..a, giving a+1 cases of 28−a unknowns.

**Exhaustiveness is a checked claim, not an assumption**
(`r9_valsplit.py check_exhaustive`): every achievable valuation profile
(α, β, ρ₂, ρ₃) with α+β ≥ n1, ρ₂+ρ₃ ≥ n2 and the degree caps is covered by
some case — verified exhaustively over the finite profile grid: R9 shape
(9,1): 5145 profiles; batch shapes a=7/8/9: 6958/6126/5295 profiles.  ALL
PASS.

A state is KILLED only if **all** its cases are exact-Q UNIT (each case
verdict is a sound necessary condition under its case hypothesis; the union
of the case varieties contains the state variety).  Per case: UNIT ⇒ case
dead; PROPER ⇒ inconclusive (still weaker than the full bridge — dm4
capped-degree polynomiality is not imposed); TIMEOUT ⇒ COST.

## 2. Implementation route (and a measured dead end)

The direct route — structured product ansätze substituted into the H-system
build — was **measured at >25 min/case** (z=1, i=5 pilot, killed): the
binomial-spread coefficients of (y+1)^i·A densify every convolution in the
per-term builder.  20 cases × 3 columns would have burned the budget on
sympy expansion alone.

The equivalent free route (`case_equations`): impose each case as **linear
equations on the generic spare coefficients** of the CACHED dm4-eliminated
build (`r9red_*.pkl`, validated by Lane B):

    v_{y+1}(p) ≥ m  ⟺  p(−1) = p′(−1) = … = p^{(m−1)}(−1) = 0   (char 0)
    v_{y−r}(p) ≥ 1  ⟺  p(r) = 0

Same case variety (A, B are exactly the generic coefficients modulo the
linear system); zero re-expansion; Singular eliminates the linear rows
instantly.  The base build's divisibility remainder rows stay (redundant
under the split, still sound).  Attack ladder, budgets, orphan-proof WSL
runner (timeout + ulimit -v 8G inside WSL): reused from Lane B verbatim.

## 3. Census (machine record `r9_valsplit_results.json`)

Exhaustiveness: **ALL PASS** (R9 shape (9,1): 5145 profiles; batch a=7/8/9:
6958/6126/5295).

| target | cases attempted | result |
|---|---|---|
| R9 z=1 (i0,j0) (i0,j1) (i1,j0) (i1,j1) | full 45s ×3-prime numroot triage | **COST ×4** — uniform TIMEOUT on every prime, ~790 s/case wall |
| R9 z=1 (i2,j0) | first prime TIMEOUT 45 s, then the driver was killed by an external outage | partial COST (log only) |
| R9 z=1 (i5,j0) | pilot: 2 primes TIMEOUT 45 s (killed at cutoff) | partial COST (log only) |
| R9 z=1 (i9,j1) **300 s CONTROL** — the maximally-constrained case (10 of dm2's 13 coefficients forced) at the full bridge's own control budget | p=10007 TIMEOUT 300 s; p=10009 TIMEOUT 302 s; third prime not completed (process gone at report time; `r9_valsplit_control.log`) | **COST at the 300 s tier** |
| R9 z=1 remaining 13 cases | NOT ATTEMPTED (outage + budget) | honest truncation |
| R9 z=2, z=3; batch 8 states | NOT ATTEMPTED | honest truncation |

Zero case kills, zero PROPER — every completed verdict is pure Gröbner
cost, no survival signal anywhere (same asymmetric semantics as Lane B).

## 4. Verdict

**WALL SURVIVES.** The valuation split — 10 linear constraints per case,
cutting the spare dimension 28 → 18 — does not move any attempted system
across the mod-p feasibility boundary: uniform 45 s TIMEOUTs across edge
(i=0,1) and middle (i=2,5) cases, and, decisively, the **structurally
easiest case (9,1) times out at the 300 s control budget** on two primes —
the same tier where the full bridge and Lane B's dm4-eliminated systems
failed.  This sharpens the localization: infeasibility persists even with
10 of dm2's 13 coefficients forced, so the swell is driven by the generic
dm3 tail coupled to the window-cap e — not by spare-dimension count.  Two
successive sound reductions (dm4 elimination; valuation split) have now
each strictly shrunk the system without denting the wall: the obstruction
is formulation-level for Gröbner engines, not size-level.  Named next
candidates, in order: (1) msolve/F4 on a split case (engine change on the
smallest formulation reached); (2) imposing the divisibility content of
e's unknown-coefficient tail factor (the flagged-omitted part, batch
states); (3) leaving GB entirely — exploit the linear structure of the
H-system's y-coefficient matrix in the forced-valuation coordinates.

## 5. [judgment] list

- **[inherited] k=6,7 window caps** for dm2, dm3 — same basis as every
  bridge kill (now PROVEN, `WINDOW_CAPS.md`).
- **[judgment] Batch split uses only the (y+1)^a factor** of e (as in Lane
  B): the unknown-coefficient tail of e also divides dm2·dm3 but is not
  imposed — sound to omit.
- **[judgment] Case coverage is a covering, not a partition** — cases
  overlap (A may itself vanish at −1); soundness needs only the covering
  direction, which is what exhaustiveness checks.
- All kills are **PENDING AUDIT**; each killed case stores integer-cleared
  generators + saturation factors (certificate-extractable,
  `kill_system` field).
