# Phase F, work item F2 -- SUB2 divisor-reconstruction kill test

**Date:** 2026-07-23. **Status: PENDING AUDIT** (same-author layer over
`convolution_descent.py`, `cascade_engine.py`, `phase_f_defects.json`,
`phase_d_states_sub2.json`, `cascade_cones_qt_inf_rl.json`; no independent
audit). Every kill below is a *candidate* kill.
**New files (uncommitted):** `phase_f2_sub2.py` (runner), `phase_f2_sub2.json`
(census), this doc, `phase_f2_sub2_verify.py` (independent verifier,
**ALL CHECKS PASS**). READ-ONLY on every audited artifact; nothing committed.

## 0. What this closes

`PHASE_F2_SCALE.md` sec.2 ran the divisor-reconstruction idea on the sub2 front
only as a *demonstration + boundary*: it showed the mechanism on the single
`a10_b0000_T1` state and then **DEFERRED** the `b != 0000` forced states, noting
they "need a sub2 geometric-regime split reconstruction that isn't recorded" in
`phase_f_defects.json`. This item is that extension.

It differs from the ALT lane (`phase_f2_scale.py`) in two ways:

1. **The full master identity, not the level-0 tie.** Sub2 states live in the
   STANDARD regime, so the necessary condition is the whole master identity
   `f31 = sum_{f=0}^7 Phi^f e^(21-3f) h_f == 0` with `Phi = c (y+1)^30 q`,
   `c = -1/6630` -- exactly the `convolution_descent.py` machinery (imported
   read-only; the reconstruction *is* the ansatz). This is a strictly stronger
   object than the ALT lane's level-0 `h_0` tie at infinity: the very
   `a10_b0000_T1` state that sec.2 could only mark **NARROWED** (a single
   top-degree relation) is a **KILL at depth 2** under the full identity (two
   consecutive top coefficients are pure in `E` and mutually inconsistent).
2. **d2 is IMPOSED, not free.** Where a state's `d2` is itself forced defect-0,
   sub2 reconstructs it (`d2 = D (y+1)^{v_t} prod_j (y-r_j)^{v_j}`) rather than
   leaving it free; the choice is recorded per state (`d2_mode`). A defect-1
   `d2` (or `d1`, `sigma`) is reconstructed on its forced defect-0 part times a
   **free linear cofactor** `(y - u)` -- one extra, *unsaturated* unknown
   (conservative, the ALT free-d2 model).

## 1. Reconstruction (every step exact)

Places `S = { t : y=-1 ; r_1..r_4 : roots of q = 2048 y^4 - 512 y^3 + 320 y^2 -
240 y + 195 }`. For a fully-forced state at cell `(a, b, T1)` whose core
divisors (`d1, sigma, e`, and `d2` where present) are all defect 0
(`deg p == sum_places v_place(p)`):

```
e     = E (y+1)^a        prod_j (y - r_j)^{b_j}      (v_t(e)=a, v_{r_j}(e)=b_j)
d1    = X (y+1)^{v_t(d1)} prod_j (y - r_j)^{v_j(d1)}
sigma = S (y+1)^{v_t(sig)} prod_j (y - r_j)^{v_j(sig)}
d2    = D (y+1)^{v_t(d2)} prod_j (y - r_j)^{v_j(d2)}   (imposed; 0 if d2_zero)
```

The per-place valuations `v_place(p)` are read from the **exact** cascade-engine
Pareto place profiles (`cascade_engine.place_profiles` / `t_place_profiles`, the
same machinery `phase_f_defects.py` used), selecting the **unique** simultaneous
per-place split that reproduces the state's degrees. For `b = (1,0,0,0)` exactly
one root is marked; arithmetic is in `Q[r]/(q(r))` with `q(r)=0` adjoined.
`b = 0000` is over `Q`.

**The reconstructibility taxonomy (honest coverage).** "Defect 0" in
`phase_f_defects.json` is *per dimension* (`delta_p = deg p - max_j v_j(p)`,
each `p` maximised independently). A clean reconstruction additionally needs a
**single simultaneous split** of *all* polynomials at once, supported on `{t,
marked roots}`. Enumerating the joint splits sorts every forced state into:

| class | meaning | reconstructible? |
|---|---|---|
| **unique** (`nsol=1`, all divisors on `{t, marked}`) | one joint split | **yes -> run the kill test** |
| **ambiguous** (`nsol>1`) | divisor could sit on any of several Galois-conjugate unmarked roots | no (needs a 2nd marked root; heavy) |
| **no-simultaneous-split** (`nsol=0`) | each dimension hits its cap but no *joint* split does | no (divisor not determined) |
| **valuation-on-unmarked-root** | a forced valuation lands on an un-`e`-marked root | no (extra marked root) |

Only the **unique** class is reconstructed and tested; the rest are recorded as
`SKIPPED_*` -- honest non-coverage, not kills.

## 2. The tower

With every polynomial determined up to its leading scalar, the reconstructed
`(d2,d1,sigma,e)` are fed to `convolution_descent.ConvolutionDescent` and the
master coefficients are **walked from the top degree downward** (top = 250 for
`deg e = 10`, 242 for `deg e = 9`), each reduced mod `q(r)`, accumulated as an
ideal, and saturated by the leading scalars (`w * prod(scalars) - 1`; free
cofactors `u` are *not* saturated). The verdict is read off an exact grevlex
saturated Groebner basis: **unit ideal => KILL**. A per-state wall (spawned
worker, `F2SUB2_HARD` s) turns single-GB blowups into `PENDING_HARD_TIMEOUT`
instead of stalling the census.

## 3. Census (201 forced states, defect <= 1 core, verdicts PENDING AUDIT)

`phase_f2_sub2.json`. Target cells: the `a9_b1000_T1` family (both d2-flag
variants, sz variants) -- the deferred marked-root front -- plus the
`a10_b0000_T1` and `a9_b0000_T1` over-`Q` remainders for cross-check.

| verdict | count |
|---|---:|
| **KILLED** | **23** |
| SKIPPED_NO_SIMULTANEOUS_SPLIT | 71 |
| SKIPPED_AMBIGUOUS_SPLIT | 64 |
| PENDING_HARD_TIMEOUT | 26 |
| SKIPPED_VALUATION_ON_UNMARKED_ROOT | 17 |

Of the **23 kills**: 14 are defect-0, 9 use a free linear cofactor (defect-1).
By novelty:

| kills | count | field | note |
|---|---:|---|---|
| **NEW (no batch overlap)** | **16** | 15 over `Q[r]/(q)`, 1 over `Q` | the deferred `b1000` q-support states + 1 `b0000` deg-1-`d2` state the batch left UNRESOLVED |
| batch-overlap (cross-validation) | 7 | `Q` | reproduce existing `batch_convolution_sub2*` degree-tuple kills |

Kill depths of the NEW set: 10 at depth 2, 4 at depth 3, 1 at depth 4, 1 at
depth 7. The 15 marked-root kills are the headline: they are exactly the
`deg_e = 10` (`= a + b`) states, which are **absent from every batch kill list**
(the batch dropped their q-support), so the reconstruction is doing the work.

### The NEW marked-root kills (`a9_b1000_T1`, over `Q[r]/(q)`)

`(deg d1, deg sigma, deg d2, deg e)`; `-inf` = identically zero.

| cell | state | degs | kill depth | d2 |
|---|---|---|---:|---|
| `a9_b1000_T1_sz0_dz0` | 0,13,33,93 | (1,0,0,10),(1,2,3,10),(1,6,3,10),(3,0,3,10) | 2,3,4,3 | reconstructed |
| `a9_b1000_T1_sz0_dz0` | 1,5,6,18 | (1,0,1,10),(1,1,0,10),(1,1,1,10),(1,3,3,10) | 2,2,2,3 | recon + lin. cofactor |
| `a9_b1000_T1_sz0_dz1` | 0,2,6,18 | (1,0,-,10),(1,2,-,10),(1,6,-,10),(3,0,-,10) | 2 | d2 == 0 |
| `a9_b1000_T1_sz1_dz0` | 0,3 | (3,-,0,10),(3,-,3,10) | 2,3 | reconstructed |
| `a9_b1000_T1_sz1_dz1` | 0 | (3,-,-,10) | 2 | d2 == 0 |

Plus one over `Q`: `a9_b0000_T1_sz0_dz0` state4 `(2,0,4,9)`, `d2` defect-1 free
linear cofactor, KILL at depth 7 (a `b0000` tuple the batch left UNRESOLVED).

## 4. Overlap accounting

The `b0000` targets (`deg_e = 9` or `10` with no marked root) are over `Q` and
their **degree tuples coincide with the generic `batch_convolution_sub2*` runs**;
7 of my kills reproduce a batch kill exactly (a consistency cross-check, not new
coverage). The `b1000` targets carry `deg_e = 10` and q-support; **no `deg_e=10`
tuple appears in any batch kill list** (`batch_convolution_sub2.json`,
`_round2.json`, `batch_convolution_overnight.json`), confirming these 15
marked-root kills are the deferred states the batch could not reach. The single
NEW over-`Q` kill (`a9_b0000` state4) is a batch UNRESOLVED tuple that the
divisor reconstruction (defect-1 `d2` cofactor) closes.

## 5. Whole cell / branch closure

**No whole sub2 cell closes, and no branch closes.** Every targeted cell retains
`SKIPPED_*` and/or `PENDING` states (sec.3 table per cell in the JSON), and only
the defect-`<=1` fraction of any cell is even reconstructible -- the bulk of each
cell's states have core defect `>= 2` and are not touched. Example: the flagship
`a9_b1000_T1_sz0_dz0` has 270 states; 57 are defect-`<=1`, of which 8 KILL, 41
SKIP (ambiguous/no-split/unmarked-root), 8 PENDING. The kills are a genuine
new *frontier dent* on the marked-root front, not a cell kill.

## 6. Verification (`phase_f2_sub2_verify.py`, ALL CHECKS PASS)

Independent of `phase_f2_sub2.py` **and** of `convolution_descent.py`. It reads
the audited source `h_f` (`t5_90t1_verify.load_h`) and `Phi`, reconstructs the
`a9_b1000_T1` factors directly, and extracts the top two master coefficients by
an **independent reversed-series short convolution** (each factor reversed about
`y=infinity`, per-`f` products aligned by absolute degree), then reduces mod
`q(r)`. It certifies two kills and walks the explicit **hand-style chain**:

- **K1** `a9_b1000_T1`, `d2 == 0`, state `(1,0,10)`: the top coefficient (degree
  250) is `E^8 (A E^17 + B)` with `A = -6561`, `B` rational, **r-free** -- so
  `E != 0` forces `E^17` to a fixed rational. The next coefficient (degree 249),
  after that substitution, is **linear in `r` with nonzero `r`-coefficient**,
  forcing `r = 1/4 in Q`. But `q` is **irreducible over Q** (checked), so a
  marked root is never rational: **CONTRADICTION => KILL** (depth 2). The
  saturated ideal `(c_250, c_249, q(r), sat E X S) = (1)` is confirmed directly.
- **K2** same state with `d2 = D` a free nonzero constant: identical chain,
  `(c_250, c_249, q(r), sat E X S D) = (1)` -- d2-freedom at that level does not
  rescue it (parallel to `phase_f2_scale_verify` K2).

The independently computed `A = -6561`, `B`, and the r-coefficient agree with the
runner's coefficients (via `convolution_descent`) -- two engines, one answer.

## 7. Honest / ambiguous points -- [judgment]

- **[J1] The 23 kills are exact saturated-Groebner facts** over the stated field
  (`Q` or `Q[r]/(q)`), with every free cofactor left UNSATURATED. The master
  identity `f31 == 0` is a necessary condition on any genuine counterexample
  (it involves only `d0,d1,d2,e`, not the g-chain), so a unit ideal under the
  reconstructed divisors is a sound kill. Two are reproduced by the independent
  verifier with an explicit irreducibility contradiction.
- **[J2] Coverage is the unique-split stratum only.** 152 of 201 forced states
  are `SKIPPED` because their divisor is **not** cleanly determined: either no
  single simultaneous place split exists (`nsol=0`, per-dimension defect-0 is
  not joint-defect-0), or it is ambiguous across Galois-conjugate unmarked roots
  (`nsol>1`), or a forced valuation lands on an unmarked root. Reporting those as
  kills would be unsound; I do not. This is the real limit of the recorded
  finite-place data for `b != 0000`.
- **[J3] The full master identity is strictly stronger than the sec.2 tie.**
  `a10_b0000_T1` is NARROWED under the level-0 `h_0` tie (`PHASE_F2_SCALE.md`
  sec.2) but KILLED at depth 2 here. The sec.2 "single top relation" is only the
  `f=0` term; the deeper `f` terms of the master identity supply the second,
  inconsistent relation. So the sec.2 sub2 NARROWED verdict was an artifact of
  using the weaker object, not a real survival.
- **[J4] d2 is IMPOSED (recorded).** For defect-0 `d2` states I reconstruct `d2`
  rather than leaving it free (the sub2 distinction the task requested;
  `d2_mode = reconstructed_defect0`). K2 of the verifier shows that at deg-0 a
  *free* `d2 = D` dies anyway, so the imposed choice is not what drives those
  kills -- it only sharpens the higher-degree `d2` cases.
- **[J5] Defect-1 kills use a conservative free linear cofactor.** The 9
  defect-1 kills reconstruct the forced defect-0 part and multiply by a free
  `(y-u)` (one unsaturated unknown). A kill under that extra freedom is sound
  whatever the true extra root; a *large* free cofactor would restore freedom
  (the ALT/pilot phenomenon), which is why only defect `<= 1` is attempted.
- **[J6] 26 PENDING are cost, not mathematics.** High-degree marked-root states
  (`deg d1 = 6`, `deg sigma in {6,7}`, plus a cofactor) blow up the grevlex GB
  over `Q[r]/(q)`; they survive even at a 125 s wall with 4 coefficients. Guarded
  out, not resolved.
- **[J7] No cell or branch closes.** The kills dent the marked-root frontier
  (15 NEW `Q[r]/(q)` states + 1 NEW `Q` state) but every cell keeps unresolved
  states; this sharpens, and does not close, the sub2 `b != 0000` front.
