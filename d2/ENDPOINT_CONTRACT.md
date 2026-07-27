# ENDPOINT_CONTRACT.md — the reduced-polygon → bigraded-engine window contract

**Status:** specification (2026-07-24). Names the exact metadata every reduced
polygon must hand the bigraded window engine, per window coefficient, so that the
engine can **automatically** locate the coefficient whose forced vanishing kills
the case. Worked example: **(50,75)** fully instantiated.

Upstream producers: `polygon_reduction.py` (the reduction front end) and the
per-case corner models (`phi_corner4.py`, `phi_75_125.py`, `f2_tower.py`, …).
Downstream consumer: the bigraded `(u-weight, y-order)` window engine (the build
target of the campaign; the object F2_TOWER.md §3 calls the "period-84 /
two-coordinate window compiler").

---

## 1. Why the contract exists (the review's regression)

The (50,75) kill (GGV3 §5, reproduced in `f2_tower.py` / F2_TOWER.md §1) is a
**window-depth contradiction**: after the terminal window unknowns are
eliminated, the forced series `C_0` has support only down to `y^-6`, but corner
**primitivity** demands a nonzero `y^-10` coefficient. The kill *is* the mismatch
between one **required-nonzero** window coefficient and the forced support.

The engine can only rediscover such kills automatically if, for **every** window
coefficient, it is told the coefficient's **status**. The review's regression was
exactly: *"the annotator must rediscover the (50,75) required coefficient
automatically."* That is impossible without a contract that tags which
coefficient primitivity requires to be nonzero. This document is that contract.

---

## 2. The contract: one status per window coefficient

A **window coefficient** is indexed by `(series, y-order)`:
`c_{s, j}` = the coefficient of `y^j` in window series `s` (e.g. `s = 0` is the
corner series `C_0`; `s = -1, -2, …` are the capped/eliminated deep window
unknowns `C_{-1}, C_{-2}, …`; forcing data `F_{-k}, λ` enter as parameters).

For each window coefficient the reduced polygon must emit **exactly one** status:

| status | meaning | engine obligation |
|---|---|---|
| **`required-nonzero`** | corner **primitivity** (the corner condition, GGV "a6"/"b6") forces `c_{s,j} ≠ 0`. | If the equations force `c_{s,j}=0`, that is a **KILL** (contradiction). The engine must flag it. |
| **`forbidden`** | the reduction / window cap forces `c_{s,j} = 0` (the slot is outside the retained polygon support, i.e. below the forced window floor). | The engine must treat `c_{s,j}=0` as an identity; a solution needing it nonzero is inconsistent with the chart. |
| **`optional`** | a free / eliminable window unknown (slack); the `E`-equations determine it, and it may be zero or nonzero. | No obligation; carried as an unknown to be eliminated. |

**The kill predicate (what the engine computes).**

```
KILL  ⇔  ∃ coefficient c_{s,j} with status = required-nonzero
          whose y-order j lies strictly below the forced support floor
          of series s  (equivalently: the equations force c_{s,j} = 0).
```

I.e. a **`required-nonzero`** coefficient landing in a **`forbidden`** slot. The
engine finds the kill by scanning the required-nonzero list against the forced
floor — no case-specific human annotation.

**Required per reduced polygon (the emitted record).**

```
endpoint_contract = {
  "case":            <tag>,                       # e.g. "F2_j0_50_75"
  "chart":           <chart id>,                  # e.g. "gamma=3 reduced (GGV3 §5)"
  "corner":          {t, kappa, a0, q},           # the corner signature
  "window_caps":     { c_{s,j}: value_or_symbol },# capped leading window coeffs (a5)
  "coefficients": [
     { "series": s, "y_order": j, "status": <status>, "reason": <cite> }, ...
  ],
  "forced_floor":    { s: j_min },                # lowest y-order the equations support, per series
  "required_nonzero":[ (s, j), ... ],             # the primitivity list (redundant index)
}
```

The `forced_floor` is produced by the engine after elimination; everything else
is supplied by the reduction. The contract is satisfied when every coefficient in
`required_nonzero` has `y_order >= forced_floor[series]` — otherwise the case is
killed at exactly the offending `(s, j)`.

---

## 3. Worked example — (50,75), the γ=3 chart, fully instantiated

Source: F2_TOWER.md §1 (γ=3 window-depth kill), `f2_tower.py` / `f2_tower_verify.py`
§B (exact, PASS). Corner `(5,20)`: `t = ceil(20/5) = 4`, `kappa = 2`, `C = y` a
**monomial** so `q = ord_y C = 1`; reduced pair `(a,b) = (2,3)`.

> **REPAIRED 2026-07-26.** This line read ``Corner `(5,20) → (7/5,2)`:
> `t=5, kappa=3, a0=5, q=2` `` — the pre-repair signature, obtained by reading
> GGV5's final chain corner `(7\5,2)` as chart data via the dictionary
> `(t,q) = (l_final, b_final)`. That dictionary is valid **only** on the
> retraction shape `b_0 = t(a_0-1)`, and `(5,20)` fails it (`20 != 4*4`), so
> `polygon_reduction.final_corner_dictionary()` now raises here. See
> `PASSPORT_75_125_REPAIR.md`.
>
> **What does NOT move, and why that is the interesting part.** Every number in
> the rest of this section — the caps `C_{-1} = a y^3`, `C_{-2} = b y^4`, the
> forced floor `j_min(0) = -6`, the required-nonzero `c_{0,-10}` — is
> **unchanged**. Those live in GGV3 §5's own reduced γ-chart coordinates, which
> are supplied by GGV3's conditions (a1)–(a6), not by our polygon bookkeeping.
> That is precisely why the γ-chart kill survived a repair that moved `t`,
> `kappa`, `C`, `N`, `Phi` and `q_window`: the depth ledger and the u-weight
> window cone are **different objects**. Cf. F2_TOWER.md's "THE BRIDGE IS
> UNVERIFIED" banner, which is the same distinction seen from the other side.

**Window caps (a5).** The corner fixes the two leading deep window coefficients:

```
C_{-1} = a·y^3        # required-nonzero: a  (the corner leading coefficient)
C_{-2} = b·y^4        # required-nonzero: b
```

**Eliminated window unknowns.** `C_{-3}, C_{-4}, C_{-5}, C_{-6}, C_{-7}` are
solved out of `E_1..E_5 = 0`, `E_6 = -F_{-1}`, `E_7 = -F_{-2}` (yielding
`F_{-1} = -3 C_{-1} C_{-2}` and the forcing relation). Status: **`optional`**
(slack, eliminated).

**The forced corner series.** Substituting the caps gives (exact, checker §B):

```
C_0 = b²y²/a² + (2/3a²)( f_8 y² + f_6 + f_4 y^-2 + f_2 y^-4 + λ y^-6 )
```

so `C_0` has **support `y^{-6} .. y^{2}`** — the **forced floor of series 0 is
`j = -6`**.

**The primitivity requirement.** Corner primitivity (a6) demands a nonzero
`y^{-10}` coefficient of the corner series:

```
c_{0,-10}  status = required-nonzero   (reason: corner primitivity a6)
```

**The contract instance (abridged coefficient ledger for series 0):**

| coefficient | status | reason |
|---|---|---|
| `c_{-1,3}` (= a) | `required-nonzero` | window cap a5, corner leading coeff |
| `c_{-2,4}` (= b) | `required-nonzero` | window cap a5 |
| `c_{-3..-7, ·}` | `optional` | eliminated by `E_1..E_7` |
| `c_{0, 2}`, `c_{0,0}`, `c_{0,-2}`, `c_{0,-4}`, `c_{0,-6}` | `optional` | forced support of `C_0` |
| `c_{0, j}` for `j < -6` | `forbidden` | below the forced floor `j_min(0) = -6` |
| **`c_{0,-10}`** | **`required-nonzero`** | **corner primitivity (a6)** |

```
forced_floor      = { 0: -6, -1: 3, -2: 4 }
required_nonzero  = [ (-1,3), (-2,4), (0,-10) ]
```

**The kill, read straight off the contract.** `(0,-10)` is in
`required_nonzero`, but `-10 < forced_floor[0] = -6`, so `c_{0,-10}` is in a
`forbidden` slot. The kill predicate fires at `(s,j) = (0,-10)`:

> **(50,75) is killed because the coefficient primitivity marks
> `required-nonzero` — `c_{0,-10}` — is forced to zero by the window equations
> (the forced `C_0` bottoms out at `y^{-6}`).**

This is precisely the coefficient the review demanded the annotator rediscover
automatically; the contract supplies it as data, and the kill predicate locates
it with no case-specific reasoning.

**Cross-check (γ=2 chart).** The sibling γ=2 reconstruction (F2_TOWER.md §1) has
the same shape with a different required-nonzero slot: `e_{-10} ≠ 0` (b6),
forced to vanish once `C_0 = 0` collapses the support to `y^{-1}`. Same contract,
different `(s,j)` — confirming the schema is chart-general, not (50,75)-special.

---

## 4. Scope / judgment notes

1. **[required-nonzero source]** The `required-nonzero` tags come from the corner
   primitivity conditions (GGV "a6"/"b6"), which are part of the reduced-polygon
   data — audited for (8,28), reconstructed for (50,75) in the GGV3 §3 charts.
2. **[forbidden source]** The `forbidden` tags are the window-cap / retained-
   support decisions of `polygon_reduction.py`; they inherit that module's
   standing judgments (chart discharged at the polygon layer; residual-gauge
   **branch completeness REOPENED**, POLYGON_REDUCTION.md R3 2026-07-24 — so a
   ramified branch may reshuffle the `optional`/`forbidden` split of the residual
   slots, though not the `c_{0,-10}` primitivity requirement).
3. **[optional / elimination]** The `optional` slack unknowns are engine-side;
   the contract only needs to declare which coefficients are slack vs capped.
4. **[bigraded lift]** In the two-coordinate `(u-weight, y-order)` engine the
   `y_order` index becomes the second grading; the same three-status contract
   applies coordinate-wise. The period-jump obstruction of F2_TOWER.md (7 → 12,
   coprime/nonaligned) is a statement about how the `forbidden` floor moves
   between rungs — the refined-lattice engine reads it from `forced_floor`.

---

## Files

| file | role |
|---|---|
| `ENDPOINT_CONTRACT.md` | this spec — the per-window-coefficient status contract |
| `polygon_reduction.py` | upstream producer of the `forbidden` / cap metadata |
| `f2_tower.py` / `f2_tower_verify.py` | the (50,75) γ=2/γ=3 reconstruction the worked example is drawn from (exact, PASS) |
