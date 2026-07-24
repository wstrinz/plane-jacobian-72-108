# WINDOW_CAPS — the k=6,7,8 window caps, recited and proven (2026-07-23)

**Mission.** `FULL_SYSTEM_BRIDGE.md` §4 flags one [judgment] under every
bridge kill: the degree/order caps for the spare window variables
`dm2, dm3, dm4` (= `d_-2, d_-3, d_-4`, rows `k = 6,7,8`) were *used* by the
bridge but proven in `T3_WINDOW_AUDIT.md` only for the jetlift window
`k = 2..5`. This doc + `window_caps_verify.py` recite the full derivation
for the deeper rows. **Verdict: PROVEN — same premises, same induction,
no new conditionality. The [judgment] flag can be retired.**

## 1. The claims

Per window variable `d_{4-k} := C_{4-k}·C4^(2k-1)` (shifted, `d3 = 0`),
for `k = 6, 7, 8`:

| var | k | ord ≥ 12k | deg ≤ 15k (sub1) | deg ≤ 14k (sub2) | stripped deg ≤ (3k / 2k) |
|---|---:|---:|---:|---:|---|
| dm2 | 6 | 72 | 90  | 84  | 18 / 12 |
| dm3 | 7 | 84 | 105 | 98  | 21 / 14 |
| dm4 | 8 | 96 | 120 | 112 | 24 / 16 |

**Role.** These caps are the *entire content* the bridge adds beyond `f31`:
`full_system_bridge.augment()` introduces `dm2,dm3,dm4` as generic stripped
polynomials of exactly these degrees (`STRIP_DEGCAP`), and every bridge
`UNIT` kill says "no polynomials *within these caps* satisfy the G-system."
If the caps were wrong (too small), a kill could be a false negative. Hence
load-bearing under every one of the 35 bridge-sweep kills, the a8 closure,
and the R9 z=0 kill.

## 2. Premises (inherited — identical to the k=2..5 window)

- **[P1]** The Prop 4.3 Newton polygons of `P` (both subcases). Loaded from
  `paper_src/upstream_facts.json` (sha-pinned transcription; audited
  verbatim against the arXiv source in `T3_WINDOW_AUDIT.md` §1).
- **[P2]** `C4 = y^7(y+1)`: ord/deg forced by corners `(8,14),(8,16)`,
  coefficients normalised by the paper's linear-change freedom.
- **[P3]** The Laurent square root `C = x^4·C4 + C3·x^3 + …` of `P`
  (`C² = P`), the paper's template object.

Nothing else. In particular **no** new assumption enters at `k = 6,7,8`.

## 3. The derivation (each step mechanically verified)

**(a) Three valuation inductions close for ALL k ≥ 1** (`T3` §3, now
verified as symbolic identities in `k`, not just instances): with
`v_{a,b}(x^i y^j) = ai + bj` and the recursion
`C_{4-k} = -(1/2C4)·(P_{8-k} + Σ_{j=1}^{k-1} C_{4-j}C_{4-k+j})`,

| direction | v(C4) | hypothesis h(k) | slice bound | step |
|---|---:|---|---|---|
| (−1,1) sub1 deg | 8 | `v(C_{4-k}) ≤ 8−k` | `≤ 16−k` | `−8+(16−k) = 8−k` ✓ |
| (−2,1) sub2 deg | 8 | `≤ 8−2k` | `≤ 16−2k` | `−8+(16−2k) = 8−2k` ✓ |
| (2,−1) ord both | −7 | `≤ 2k−7` | `≤ 2k−14` | `7+(2k−14) = 2k−7` ✓ |

The product bound `h(j)+h(k−j)` is *j-free and equals the slice bound* in
every direction (exact closure — no slack), and the slice bounds come from
the corner maxima `max(j−i)=8`, `max(j−2i)=0`, `max(2i−j)=2`, **computed**
from [P1], not transcribed. The `k=6,7,8` steps consume the genuine slices
`P_2, P_1, P_0` (no vacuity).

**(b) D-transform.** `D_j := C_j·C4^(7−2j)` (`D_4 = 1`; `d_{4-k} = D_{4-k}`).
Symbolically in the x-exponent `j`:
`deg D_j ≤ (j+4) + 8(7−2j) = 60−15j` (sub1), `≤ 2j + 8(7−2j) = 56−14j`
(sub2), `ord D_j ≥ (2j−1) + 7(7−2j) = 48−12j`. Substituting `j = 4−k`
gives exactly `deg ≤ 15k / 14k`, `ord ≥ 12k` — for every `k`, in
particular the three rows of §1. (Polynomiality of `D_j` itself is the
already-audited `verify_derivation.py` §C exponent-cancellation, proven
down to `j = −13`.)

**(c) The d3-killing shift, made polynomial.** The unique shift
`x → x − s`, `s = c_3/(4C4)`, kills `c~_3`. New recitation (cleaner than
the valuation argument, and verified exactly): in D-coordinates the shift
is the **polynomial identity**

```
D~_j = Σ_{m=j}^{4} binom(m, m−j) · D_m · (−D_3/4)^(m−j)
```

(series-recomposition identity verified per row `j = 3..−4`; all `C4`
exponents cancel — so the *shifted* window variables are automatically
polynomial). It gives `D~_3 = D_3 − D_3 = 0`, and each term obeys
`deg[D_m·(−D_3/4)^(m−j)] ≤ cap(m) + (m−j)·cap(3) = cap(j)` because
`cap(3) = 15/14/12` is exactly the per-step slope — identically in `(m,j)`,
all three directions. **The caps survive the shift with no loss.**

**(d) Generic end-to-end corroboration.** For both subcases: random
exact-`Q` polynomials `P_0..P_7` supported on the *computed* corner hull
(`P_8 = C4²` forced), the division-free D-recursion run down to `D_{-4}`,
and every cap checked on unshifted and shifted variables for `k = 2..8`.
All hold — and the degree caps are **attained on all 7/7 rows** in both
regimes (the bounds are tight, matching T3's "empirically tight" and the
`Φ` weight-17 witness `deg 238 = 14·17`, `ord 204 = 12·17`).

## 4. Checker

`window_caps_verify.py` — 81 checks, ALL PASSED (`--quiet` supported;
exit 0 iff all pass). Sections W0–W5 mirror §3 plus consumer cross-checks:
`full_system_bridge.WEIGHT`/`STRIP_DEGCAP` and `jetlift.CONFIGS` agree with
the caps; spare-unknown totals 66 (sub1) / 45 (sub2); `f1` re-derived from
its forcing ODE for the tightness witness.

## 5. What this changes

- `FULL_SYSTEM_BRIDGE.md` §4's "[judgment — well-supported] Window caps for
  k=6,7,8" and `BRIDGE_SWEEP.md`'s "[judgment — inherited]" entry are now
  **recited and mechanically verified** — same trust tier as the k=2..5
  window (T3). Remaining conditionality is [P1]–[P3] only, shared by the
  entire program.
- Bonus lemma worth reusing: the D-coordinate shift identity §3(c) — it
  removes every appeal to "the shift preserves valuations" in favour of a
  term-by-term polynomial degree count.
