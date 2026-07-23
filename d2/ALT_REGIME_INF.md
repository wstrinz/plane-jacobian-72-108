# Degree layer (max-plus at infinity) for the alternate regime

Date: 2026-07-22
Verifier: `alt_regime_inf_verify.py` (PASSES)
Scope: the 26 sub1 alternate-regime strata, `a in [11,15]`, `v = 30-3a < 0`.

This document transports the STANDARD-regime infinity layer
(`cascade_engine.py`: `deg_h_options` / `descend_options_inf` /
`inf_place_profiles`, `a <= 10`, `v = 30-3a >= 0`) to the FLIPPED regime. It is
a NEW derivation lane: it does not edit the standard engine, and every step is
source-linked to `ALT_REGIME.md`, `ALT_REGIME_L2.md`, `ALT_REGIME_AUDIT.md`,
`f31_graded.txt` and the sub1 caps in `sub1_cascade_verify.py`.

Notation matches `ALT_REGIME_L2.md` §1 verbatim:

```text
e = t^a E,   u = c q,   v = 30-3a < 0,   w = |v| = 3a-30 > 0,   T = t^w,
E = ehat (a unit at t and at every q-root),   deg_E = deg e - a.
```

Throughout, "degree" is the degree in the plane variable (equivalently the
valuation at infinity `v_inf = -deg`, `deg t = 1`, `deg T = w`, `deg u = 4`).
This is the INFINITY place. It is a different place from the t-adic and q-adic
valuation lemmas of `ALT_REGIME.md` §"Terminal plus first-level local lemmas"
(those are finite-place parity arguments); the two layers are independent.

---

## (a) Flipped per-level identities the alternate chain uses

`ALT_REGIME.md` "Flipped reduction" replaces the standard polynomial reduction
`F = t^(21a) G` (valid only as a Laurent equality when `v < 0`, table row 2) by
the genuine polynomial reduction, because `21a + 7v = 210` is the UNIQUE minimum
of the f-term t-order `21a + f v` (strictly decreasing in f for `v <= -3`,
`ALT_REGIME_AUDIT.md` C33-a):

```text
F = t^210 G',        G' = sum_{f=0..7} t^((7-f) w) u^f E^(21-3f) h_f(d~).   (G')
```

The chain is the DESCENDING cascade `(D_t)` of `ALT_REGIME_L2.md` §1
(equivalently `ALT_REGIME.md` "Descending cascade" with `r_l` in place of the
Laurent `g_l`). Writing the recursion in one closed form
(`ALT_REGIME_AUDIT.md` C33-b, `T r_{f-1} = E^(3(7-f)) h_f + u r_f`,
`r_7 = r_{-1} = 0`):

```text
level f=7 (TOP ANCHOR):  T r_6            = h_7                       (u r_7 = 0)
level f=6..1           :  T r_{f-1}        = E^(3(7-f)) h_f + u r_f
level f=0 (BOTTOM CLOSE):  0               = E^21 h_0 + u r_0
```

Correspondence to the standard chain `t^v g_{l+1} = ehat^3 g_l + u^l h_l`:

| standard-regime object (`inf_place_profiles`)          | alternate-regime replacement                                             |
|--------------------------------------------------------|--------------------------------------------------------------------------|
| ascending cofactors `g_1,...,g_terminal`               | DESCENDING cofactors `r_6, r_5, ..., r_0`                                 |
| level-0 anchor `t^v g_1 = h_0`                          | TOP anchor `T r_6 = h_7` and BOTTOM closing `E^21 h_0 + u r_0 = 0`        |
| LHS shift `v = 30-3a >= 0` (on `t^v g_{l+1}`)          | `+w = 3a-30 = |v|` (on `T r_{f-1}`, `deg T = w`)                          |
| `g`-side exponent `3*deg(ehat)` (constant in level)    | `3(7-f)*deg_E` (LEVEL-DEPENDENT: `E^(21-3f)` multiplies `h_f`)            |
| `h`-side shift `4l = DEG_U*level` (grows with level)   | constant `4 = DEG_U` (a SINGLE power of `u` on `r_f`; the `u`-powers      |
|                                                        | instead accumulate through the telescope, not per level)                 |
| terminal identity `ehat^3 g_7 + u^7 h_7 = G` (bottom-up)| survives verbatim via `g_1 = T h_0`, `g_{l+1}=T(E^3 g_l+u^l h_l)`,        |
|                                                        | `E^3 g_7 + u^7 h_7 = G'` (`ALT_REGIME.md` table; `ALT_REGIME_L2.md` (U))  |

The top anchor `T r_6 = h_7` follows from `h_7 = 8192 d1^2` carrying the least
explicit t-power (`ALT_REGIME_AUDIT.md` C33-b); `E^21 h_0 + u r_0 = 0` is the
new TERMINAL condition, "rather than a top equation reached by upward
t-divisibility" (`ALT_REGIME.md` "Descending cascade"). The recursion
telescopes exactly to `G' = 0` (`ALT_REGIME.md`/`ALT_REGIME_L2.md` [B];
`alt_regime_verify.py` check 3), so the degree bookkeeping below is bookkeeping
for `(D_t)` itself.

At a q-root `p` with `b = v_p(E) > 0` the same shape holds with `s = 3b-1`,
`A = T p^s` in place of `T` (`ALT_REGIME_L2.md` (D_p)); the degree layer is the
`t`-instance, where `T` is the only relevant pole factor.

---

## (b) Degree identities at infinity (max-plus)

Apply `deg(.)` to each line of `(D_t)`. With `deg T = w`, `deg u = 4`,
`deg(E^(3(7-f)) h_f) = 3(7-f) deg_E + deg h_f` (E is a nonzero polynomial
factor, so the product degree is exact), and writing `H_f := deg h_f`,
`R_f := deg r_f`:

```text
TOP ANCHOR   (f=7):   w + R_6 = H_7                                       (I7)
LEVELS f=6..1     :   w + R_{f-1} = max( 3(7-f) deg_E + H_f ,  4 + R_f )  (If)
BOTTOM CLOSE (f=0):   max( 21 deg_E + H_0 ,  4 + R_0 )   must be a TIE    (I0)
```

with the standard max-plus contract (`deg_h_options` /`descend_options_inf`
docstrings): the maximum on the right is ATTAINED (LHS equals it) unless the two
right-hand terms TIE, and a drop below the max is permitted only as a recorded
leading-coefficient cancellation (`Obligation` "leading_cancellation" /
"degree_tie_drop"). The degree of a sum never exceeds the max of the summand
degrees, with equality unless the leading forms cancel -- and cancellation can
only occur where the two degrees are equal. Concretely:

- `(I7)` has a UNIQUE term, so `R_6 = H_7 - w` is FORCED (no drop possible),
  the exact analogue of the standard level-0 anchor forcing `deg g_1`.
- `(If)`: if `3(7-f) deg_E + H_f != 4 + R_f`, the max is unique and
  `R_{f-1} = max - w` is FORCED. If the two are equal, `R_{f-1}` may be
  `max - w` or drop below it, the drop carried as a leading-cancellation
  obligation of depth `= (max - w) - R_{f-1}`.
- `(I0)`, the closing anchor, is where a CONTRADICTION is read off: the sum
  `E^21 h_0 + u r_0` must be identically `0`, so its leading term must cancel;
  this REQUIRES the tie `21 deg_E + H_0 = 4 + R_0`. If instead one term
  strictly dominates (unique maximum), the sum has a nonzero leading term of
  that degree and CANNOT be zero -- the would-be counterexample dies purely on
  the degree layer. (This is the flipped mirror of the standard exact-identity
  closing `ehat^3 g_l = -u^l h_l` in `descend_options_inf`, case `g_above_zero`.)

**Exact contribution of the flipped power.** `deg t^((7-f)|v|) = (7-f) w` is the
per-term prefactor in `(G')`; after factoring the global `t^210`, the residual
prefactor is exactly `t^((7-f)w)` (`ALT_REGIME.md` "Flipped reduction";
`alt_regime_verify.py` check 1 asserts `orders[f]-210 == (7-f) w`). In the
cascade `(D_t)` this prefactor is realised one power of `T = t^w` per descending
step -- hence the constant `+w` shift in every `(If)`, not a growing `(7-f)w`.
The `deg u = 4` contribution appears once per level (single `u` on `r_f`);
across the whole telescope the total `u`-power reaches `u^7` (matching `(G')`'s
`u^f`), but at INFINITY only the per-level `+4` enters each `(If)`.

**Worked table (verifier instance, `a=12`, `w=6`, `deg_E=3`, degrees at the
sub1 caps so `H_f = 60-6f` is tight).** Computed by `alt_regime_inf_verify.py`:

| f | line of `(D_t)`               | `3(7-f)deg_E + H_f` | `4 + R_f` | max | forced `R_{f-1}` |
|--:|-------------------------------|--------------------:|----------:|----:|-----------------:|
| 7 | `T r_6 = h_7`                  | `H_7 = 18` (unique) |    --     |  18 | `R_6 = 12`       |
| 6 | `T r_5 = E^3 h_6 + u r_6`      |                  33 |        16 |  33 | `R_5 = 27`       |
| 5 | `T r_4 = E^6 h_5 + u r_5`      |                  48 |        31 |  48 | `R_4 = 42`       |
| 4 | `T r_3 = E^9 h_4 + u r_4`      |                  63 |        46 |  63 | `R_3 = 57`       |
| 3 | `T r_2 = E^12 h_3 + u r_3`     |                  78 |        61 |  78 | `R_2 = 72`       |
| 2 | `T r_1 = E^15 h_2 + u r_2`     |                  93 |        76 |  93 | `R_1 = 87`       |
| 1 | `T r_0 = E^18 h_1 + u r_1`     |                 108 |        91 | 108 | `R_0 = 102`      |
| 0 | `0 = E^21 h_0 + u r_0`         |                 123 |       106 | 123 | TIE fails -> die |

In this generic (non-solution) window the `E^(3(7-f)) h_f` term strictly
dominates at every level (no ties, hence no drops), and the closing anchor has a
unique maximum `123 != 106`, so `E^21 h_0 + u r_0` has degree `123` and is
nonzero: the degree layer alone refutes the instance. A genuine counterexample
would have to realise the tie `21 deg_E + H_0 = 4 + R_0` at closing (and,
generally, ties at intermediate levels to lower the `R_f` enough), each tie a
leading-cancellation obligation handed to the finite-place layers.

---

## (c) Degree caps applicable in the alternate regime

From the sub1 window (`sub1_cascade_verify.py`: stripped caps
`(d2,d1,d0,e) = (6,9,12,15)`, `sigma <= 12`, `deg h_f <= 60-6f`, budget
`a + sum b <= 15`). Which survive into the flipped degree layer:

- **Source-variable window caps SURVIVE verbatim** — `deg d2 <= 6`,
  `deg d1 <= 9`, `deg d0 <= 12`, `deg sigma <= 12`, `deg e <= 15`, hence
  `deg_E = deg e - a <= 15 - a <= 4`. These are window facts independent of the
  regime (`ALT_REGIME.md`: q-place data and the graded identity are
  "window-independent" / "regime-independent"; `ALT_REGIME_AUDIT.md` lists them
  as the trusted ground truth). `alt_regime_inf_verify.py` builds its window
  inside these caps and asserts `deg_E = 15-a`.
- **`deg h_f <= 60-6f` SURVIVES** — it is a bound on the graded coefficient
  `h_f(d~)` from the same source caps; the checker asserts it for all `f`, tight
  (`= 60-6f`) at the cap window. This is what feeds `H_f` into `(If)`.
- **Terminal g-caps SURVIVE** — `ALT_REGIME.md` table: "`deg g7 <= 46`,
  `deg g6 <= 48` survive, sharpening in a stratum to `46 - 3 sum b_i`,
  `48 - 3 sum b_i`." These bound the bottom-up polynomial auxiliaries `g_l`
  (used for the q-terminal identities), and are the analogue of the standard
  `g_caps` consumed by `descend_options_inf`.
- **Old lower-level global `g_l` caps and the t-coupling DO NOT transfer** —
  `ALT_REGIME.md` table, last row: "negative `v` reverses the degree recursion
  and t-adic edge." The standard `inf_place_profiles` raises
  `ValueError("infinity layer requires the standard regime a <= 10")`, so its
  per-level `g_cap` schedule is not imported; the flipped layer uses the caps
  above plus the max-plus recursion `(If)` directly.

---

## (d) Honest gap list — [judgment]

Items I could NOT pin down verbatim from the repo sources and resolved by my own
reasoning (each marked and defended):

- **[judgment] No explicit degree cap on the descending cofactors `r_f`.** The
  repo docs give caps on `g_l` (terminal `deg g7 <= 46`, `deg g6 <= 48`) but not
  on the flipped `r_f`. I therefore do NOT impose an independent `r_f` cap; I
  only use `(If)` to DERIVE `R_f` from `H_f`, `deg_E`, `w`. This is conservative
  (it never asserts a kill from an unproven `r_f` bound). Deriving a rigorous
  `r_f` cap (e.g. via `T r_6 = h_7`, `R_6 = H_7 - w <= 60-42-w`) is left open.
- **[judgment] The infinity place is treated as strictly separate from the
  finite-place parity lemmas.** `ALT_REGIME.md`'s `v_P` lemmas (`P = t` or a
  q-root) are finite-place valuations; I read the DEGREE (infinity) layer off
  the same identities `(D_t)`. The repo does not state a standard-regime
  infinity layer for `a > 10` (the code guards against it), so the transport of
  the max-plus contract to the flipped chain is my derivation, justified by the
  fact that `(D_t)` is an exact polynomial identity and `deg` is a valuation.
- **[judgment] "Attained max unless a tie" is verified as a rational-function
  degree law, not as a per-stratum kill.** `alt_regime_inf_verify.py` forms
  `r_f` as exact rational functions from a random window and checks the max-plus
  identities and the closing contradiction; it does NOT enumerate the 26 strata
  or claim new kills. Turning the degree layer into a stratum-by-stratum kill
  engine (the analogue of `inf_place_profiles` Pareto enumeration, coupling
  `sum_p v_p(r_f) <= deg r_f`) is not attempted here.
- **[judgment] `deg u = 4` and the `21a+7v=210` minimum are used exactly as
  written in the sources**; the only inference is that a SINGLE `u` (not `u^l`)
  multiplies `r_f` in `(D_t)`, so the per-level `h`-shift is the constant `4`
  rather than the growing `4l` of the standard chain. This follows directly from
  the recursion `T r_{f-1} = E^(3(7-f)) h_f + u r_f` as printed in
  `ALT_REGIME_L2.md` §1 and confirmed by the exact telescope in
  `alt_regime_verify.py` check 3.
