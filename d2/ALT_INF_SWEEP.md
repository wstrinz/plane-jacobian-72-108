# Alternate-regime max-plus degree sweep over the 27 open flipped branches

Date: 2026-07-22
Sweep: `alt_inf_sweep.py` -> `alt_inf_sweep.json`
Verifier: `alt_inf_sweep_verify.py` (PASSES)
Derivation transported: `ALT_REGIME_INF.md` (b) identities (I7)/(If)/(I0)

This runs the flipped-cascade **degree layer** (max-plus at infinity) of
`ALT_REGIME_INF.md` as an actual stratum-by-stratum enumeration -- the step
`ALT_REGIME_INF.md` (d) explicitly left open ("does NOT enumerate the 26
strata or claim new kills"). It sweeps every degree assignment in the sub1
sandwich for each open branch, runs the descending chain (D_t), and records
which degree states the degree layer refutes and which survive with what
obligations.

## 1. Branch list and source

Source: **`ALT_REGIME_L2.md` section 4 "Per-input-branch verdict"** -- the
`O` (open) rows remaining after the six h6/h5 cone kills of that document.
That table lists exactly **27 survivors: 13 T1 + 14 T2** (its section 1 header
states the count). These are the strata `a in {11,12,13,14}`; `a=15` and the
extra `a<=14` rows were already killed by `ALT_REGIME.md` / `ALT_REGIME_L2.md`
and are not in the denominator. The 27 (hardcoded in `alt_inf_sweep.py`
`OPEN_BRANCHES`, asserted 13+14):

- **T1 (13):** a11 `(0000)(1000)(1100)(1110)(1111)(3000)`;
  a12 `(0000)(1000)(1100)(1110)(3000)`; a14 `(0000)(1000)`.
- **T2 (14):** a11 `(0000)(1000)(1100)(1110)(1111)(3000)(3100)`;
  a12 `(0000)(1000)(1100)(1110)`; a13 `(0000)(1000)`; a14 `(0000)`.

(Cross-check: the same strata appear as `alternate_regime` rows `a in [11,15]`
in `split_place_ledger_sub1.json`; removing the `ALT_REGIME.md`/`_L2.md` kills
leaves this set.)

## 2. Semantics (mirrors the standard infinity layer, imported not edited)

`alt_inf_sweep.py` **imports** `cascade_engine.py` and reuses its exact code:
`deg_h_options` (h-degree drop rules + `Obligation` objects), the monomial
tables from `cascade_signature.load_levels()` (exponent order `(d2,d1,sigma,e)`,
the `e`-slot carrying `deg e`), and `DEG_U = 4`. No existing file is modified.

Degrees are the degree in the plane variable (valuation `-deg` at infinity);
`deg T = w = 3a-30`, `deg u = 4`, `deg_E = deg e - a`. Writing `H_f = deg h_f`,
`R_f = deg r_f`, the chain (D_t) gives (`ALT_REGIME_INF.md` (b)):

```
(I7) top anchor f=7 :  w + R_6      = H_7                                (unique)
(If) levels   f=6..1 :  w + R_{f-1} = max( 3(7-f) deg_E + H_f , 4 + R_f )
(I0) bottom  close f=0:  max( 21 deg_E + H_0 , 4 + R_0 )  must TIE
```

Rules, identical to `deg_h_options` / `descend_options_inf`:

- `H_f` = tropical max over surviving monomials. A **unique** achiever forces
  `H_f` (no drop, no vanishing). Several achievers permit a drop to any depth
  down to `0`, or identical vanishing (`NEG_INF`), each a recorded obligation
  (`degree_tie_drop` / `identical_vanishing`).
- In (If) the sum degree is the max of the two term degrees; it may drop below
  only when the two terms **tie** (`term1 == term2`), the drop carried as a
  `leading_cancellation` obligation. A unique max is forced; a forced negative
  `R_{f-1}` (nonzero polynomial of negative degree) is impossible -> dead path.
- (I0) can vanish **only** on a tie. A unique maximum at the close is a
  contradiction (nonzero leading term) -- the degree-layer kill.

Degree sandwich (sub1): `deg d2<=6, deg d1<=9, deg sigma<=12, deg e<=15`, and
`deg e >= a + sum(b)` (since `E = prod p_i^{b_i} F`, `deg F>=0`, so
`deg_E = deg e - a >= sum b`). `deg h_f <= 60-6f` is automatic and was
cross-checked (verifier PART 1). **No independent cap on `R_f`** is imposed
(`ALT_REGIME_INF.md` [judgment]): each `R_f` is propagated by the chain.

Exact integer arithmetic throughout; `NEG_INF` encodes an identically-zero
polynomial.

## 3. Results

All **27 branches are OPEN** on the degree layer alone. But the degree layer
is far from empty: it **refutes 33670 of 38360 (87.8%) of the enumerated
degree states**, leaving 4690 surviving states, each carrying explicit
leading-cancellation obligations that the finite-place (residue) layers must
discharge.

| verdict | branches | degree states | surviving | killed |
|:--|--:|--:|--:|--:|
| OPEN | 27 | 38360 | 4690 | 33670 (87.8%) |
| KILLED | 0 | | | |

Per-branch (full data in `alt_inf_sweep.json`):

| id | a | sum_b | br | verdict | states | surv | killed |
|:--|--:|--:|:--|:--|--:|--:|--:|
| a11_b0000_T1 | 11 | 0 | T1 | OPEN | 5600 | 591 | 5009 |
| a11_b1000_T1 | 11 | 1 | T1 | OPEN | 4480 | 487 | 3993 |
| a11_b1100_T1 | 11 | 2 | T1 | OPEN | 3360 | 395 | 2965 |
| a11_b1110_T1 | 11 | 3 | T1 | OPEN | 2240 | 319 | 1921 |
| a11_b1111_T1 | 11 | 4 | T1 | OPEN | 1120 | 259 | 861 |
| a11_b3000_T1 | 11 | 3 | T1 | OPEN | 2240 | 319 | 1921 |
| a12_b0000_T1 | 12 | 0 | T1 | OPEN | 4480 | 459 | 4021 |
| a12_b1000_T1 | 12 | 1 | T1 | OPEN | 3360 | 370 | 2990 |
| a12_b1100_T1 | 12 | 2 | T1 | OPEN | 2240 | 296 | 1944 |
| a12_b1110_T1 | 12 | 3 | T1 | OPEN | 1120 | 238 | 882 |
| a12_b3000_T1 | 12 | 3 | T1 | OPEN | 1120 | 238 | 882 |
| a14_b0000_T1 | 14 | 0 | T1 | OPEN | 2240 | 227 | 2013 |
| a14_b1000_T1 | 14 | 1 | T1 | OPEN | 1120 | 175 | 945 |
| a11_b0000_T2 | 11 | 0 | T2 | OPEN | 520 | 30 | 490 |
| a11_b1000_T2 | 11 | 1 | T2 | OPEN | 416 | 27 | 389 |
| a11_b1100_T2 | 11 | 2 | T2 | OPEN | 312 | 24 | 288 |
| a11_b1110_T2 | 11 | 3 | T2 | OPEN | 208 | 22 | 186 |
| a11_b1111_T2 | 11 | 4 | T2 | OPEN | 104 | 20 | 84 |
| a11_b3000_T2 | 11 | 3 | T2 | OPEN | 208 | 22 | 186 |
| a11_b3100_T2 | 11 | 4 | T2 | OPEN | 104 | 20 | 84 |
| a12_b0000_T2 | 12 | 0 | T2 | OPEN | 416 | 27 | 389 |
| a12_b1000_T2 | 12 | 1 | T2 | OPEN | 312 | 24 | 288 |
| a12_b1100_T2 | 12 | 2 | T2 | OPEN | 208 | 22 | 186 |
| a12_b1110_T2 | 12 | 3 | T2 | OPEN | 104 | 20 | 84 |
| a13_b0000_T2 | 13 | 0 | T2 | OPEN | 312 | 22 | 290 |
| a13_b1000_T2 | 13 | 1 | T2 | OPEN | 208 | 20 | 188 |
| a14_b0000_T2 | 14 | 0 | T2 | OPEN | 208 | 17 | 191 |

### Kill mechanisms (histogram in the JSON)

1. **Bottom-close unique maximum (23870 states).** For every reachable chain,
   `21 deg_E + H_0 != 4 + R_0`, so `E^21 h_0 + u r_0` has a nonzero leading
   term and cannot be `0`. This is the (I0) contradiction of
   `ALT_REGIME_INF.md` (b).
2. **T1 top-anchor negativity (~9500 states).** `2 deg d1 < w` forces
   `R_6 = H_7 - w = 2 deg d1 - w < 0`; `r_6 = h_7/T` is not a polynomial. The
   degree-layer image of `ALT_REGIME.md`'s `2 v_t(d1) >= w` anchor.
3. **Dead intermediate level (280 states).** Some (If) forces `R_{f-1} < 0`.

### Example surviving-state obligations

- Branch `a11_b0000_T1`, state `(deg d2,d1,sigma,e)=(5,2,10,11)`, `deg_E=0`:
  closes with the tie `21 deg_E + H_0 = 46 = 4 + R_0`. Obligations: an
  `h_0` `degree_tie_drop` of depth 4 (tied set `{sigma^2 d2^6, sigma^3 d2^4,
  sigma^4 d2^2, sigma^5}`) plus the closing `leading_cancellation` (depth 0)
  between `E^21 h_0` and `u r_0`.
- Branch `a12_b0000_T1`, tight-cap state `(6,9,12,15)`, `deg_E=3`: survives
  only via a **depth-17** `degree_tie_drop` on `h_0` (26-way tie) that pulls
  `21 deg_E + H_0` from `123` down to `106 = 4 + R_0`. This is precisely the
  instance `ALT_REGIME_INF.md`'s worked table calls "die": it dies **iff** the
  cancellation is disallowed. See [judgment] J1 below.

## 4. Verification (`alt_inf_sweep_verify.py`, PASSES)

Independent of the sweep (never imports `alt_inf_sweep`); recomputes `H_f` from
the monomial table and builds `r_f` as exact sympy rational functions.

- **PART 1** grounds the tropical `H_f` against exact sympy `deg h_f` for the
  generic `a=12` window (all 8 match, all `<= 60-6f`).
- **PART 2** kill A (top anchor): `a=12`, `deg d1=2`, `2*2=4 < w=6`;
  `r_6=h_7/T` has negative degree / non-trivial denominator.
- **PARTS 3-4** two bottom-close kills (`a=12` state `(-,3,-,12)`; `a=11`
  state `(-,2,-,11)`), re-derived by a hand chain with a **unique** maximum at
  every level (so no cancellation is available), ending `48 != 37` and
  `44 != 36`; each confirmed on an exact sympy window whose closing residual is
  nonzero of the dominant degree.
- **PART 5** open survivor `(5,2,10,11)`: an explicit witness chain is checked
  against every (I7)/(If)/(I0) identity, drops gated on genuine ties, close a
  tie.

## 5. Honest / ambiguous points -- [judgment]

- **[judgment] J1 -- "generic instance dies" is a narrower claim than this
  sweep's over-approximation.** `ALT_REGIME_INF.md`'s worked `a=12` table
  concludes "TIE fails -> die" for the tight-cap state. That holds **only in
  the no-cancellation sub-scenario**. This sweep adopts the policy stated in
  `ALT_REGIME_L2.md` section 1 -- "every tied leading form ... capable of
  cancelling to arbitrary depth ... a nonempty cone is recorded as an open
  residue problem" -- and the `deg_h_options` drop rule of the standard engine.
  Under that (sound, over-approximating) policy the tight-cap state **survives**
  with a depth-17 obligation. So the branches are OPEN, not killed, on the
  degree layer alone; the sweep's job is to catalog the surviving states and
  their obligations, which it does. A kill would need the finite-place residue
  layers to forbid the required cancellations.

- **[judgment] J2 -- T2 top anchor.** For T2 (`d1==0`), `h_7 = 8192 d1^2 == 0`,
  so the top anchor `T r_6 = h_7` gives `r_6 == 0` (`R_6 = NEG_INF`), exactly
  `ALT_REGIME.md` "Descending cascade" T2 paragraph ("For T2 (d1=0), r_6=0 and
  T r_5 = E^3 h_6"). I do **not** hardcode a separate T2 anchor: seeding
  `R_6 = NEG_INF` makes the general (If) at `f=6` collapse `4 + R_6` to
  `NEG_INF`, leaving the unique h-term `w + R_5 = 3 deg_E + H_6` -- structurally
  the same identity `ALT_REGIME.md` writes as `T r_5 = E^3 h_6`. The finer
  `T | r_5` divisibility strengthening in that paragraph is a **t-adic
  (finite-place) fact** and is deliberately not used here; the degree layer only
  reads `deg`.

- **[judgment] J3 -- zero flags enumerated, not chosen.** `sigma==0` and
  `d2==0` are enumerated as independent degree states (`deg=NEG_INF`), an
  over-approximation (never a silent restriction). The single combination
  excluded is **T3** (`d1==0` and `sigma==0`), which `split_place_ledger_sub1`
  states is "excluded globally"; those states are skipped and not counted.

- **[judgment] J4 -- no `R_f` cap; `r_f==0` only via a recorded obligation.**
  Per `ALT_REGIME_INF.md` (d) there is no proven cap on `deg r_f`, so none is
  imposed; `R_f` is purely chain-propagated. `r_f == 0` (`R_f = NEG_INF`) is
  reachable only through an `identical_vanishing`/`sum_vanish` obligation on the
  level above -- never assumed for free.

- **[judgment] J5 -- degree vs valuation layers are independent.** The
  top-anchor kill uses `2 deg d1 >= w` (a degree fact), which is weaker than the
  finite-place `2 v_t(d1) >= w` parity lemma of `ALT_REGIME.md`. The two layers
  are independent; this sweep is the degree (infinity) layer only and does not
  import the finite-place valuation bounds. Combining layers (intersecting
  surviving states with the residue/valuation constraints) is the next step and
  is left to the finite-place engine.
