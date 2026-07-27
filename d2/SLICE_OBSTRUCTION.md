# SLICE_OBSTRUCTION.md — the positive-slice obstruction, compiled; the stacked P/Q audit; `a_t >= 9`

2026-07-25. Checker: `slice_obstruction_basis.py` (`--quiet --deep`, exit 0).
Read-only on every existing artifact: this lane wrote only
`SLICE_OBSTRUCTION.md`, `slice_obstruction_basis.py`,
`slice_obstruction_stage.json`. Pure sympy — no Singular, no msolve, no WSL,
no subprocess, no solver. Nothing was written to `state_kill_ledger.json`,
`proof_dag.json`, `phase_d_states*.json`, `frontier_rebuild.py`, or any
`positive_slice*` / `alt_*` / `frontier_*` file.

> **HEADLINE.** `POSITIVE_SLICE.md` emptied the last standard-sub2 cell with
> three conditions, all from `P = C^2`, all at `t = 0`. Its §9.3 flags the
> `Q = C^3 + lambda*C^-1 + F` side and the deeper jets as unused, and says the
> argument does not transfer to sub1. This lane builds the general obstruction
> calculus and runs it.
>
> **The Q-side slice formula is derived and exact:**
> `Q_M = y^(2M-3) * [u^(12-M)] H(u)^3 / t^(21-2M)`. The brief's index/sign
> claim is **CONFIRMED**, by derivation.
>
> **Stacking is the whole story.** `P` alone and `Q` alone each have cokernel
> **0** at every level, in both windows: a free fresh coefficient absorbs
> either support map completely. Stacked — because the *same* fresh coefficient
> must serve both — the cokernel is **`2n-3`**. Counting the two sides
> separately reports zero and misses everything.
>
> Because `L_n^P = 2` and `L_n^Q = 3`, the fresh coefficient cancels
> identically in `2*r_n - 3*p_n`, leaving the **cell-independent,
> window-independent** obstruction
> ```
> t^(2n-3)  |  [u^n]( 3*K^2 + 2*K^3 ),        K := H - 1
> ```
> (since `2H^3 - 3H^2 = -1 + 3K^2 + 2K^3`).
>
> **Cascading it forces `v_t(h_k) >= 2k-1` for `k = 1..5`.** The `d3`-killing
> shift is triangular across zero and does not move the spares, so
> `h_5 = dm1 = e` exactly, and `v_t(e) = a_t`. Hence
>
> **`a_t >= 9`, on both branches, in both windows.**
>
> | window | C08/C20 | cells | flag cases | states |
> |---|---|---|---|---|
> | **sub1** | ON (`rl`) | **34 -> 11** | **314 -> 134** | **7275 -> 3657** |
> | **sub1** | OFF (`norl`) | **34 -> 11** | **322 -> 142** | **8889 -> 4696** |
>
> **23 of the 34 standard-sub1 cells die.** Identical kill under both settings.
> The alternate regime is **untouched** (`a_t in {12,14} >= 9`) — a real and
> flagged negative result.
>
> **Two things a reader should attack first.** (i) `[QQ1]`, the `alpha`-strip
> WLOG `Q = C^3 + lambda*C^-1 + F` with `v_{1,0}(F) = -5`, which is
> `PROOF_INVENTORY.md` premise **C3 at confidence 2/4** — it is what makes the
> Q column exist at all, and it is not re-proved here. (ii) `[Q8]`/`S3.4`, the
> identification `h_5 = dm1`, which is what turns a valuation into a census
> delta.

---

## 0. What was asked, and what survived contact

| claim under test | verdict |
|---|---|
| unrestricted square-series slices are triangular with diagonal 2 (`m` for `H^m`) | **DERIVED** (S3.2). `p_n = 2h_n + q_n`, `r_n = 3h_n + q_n^Q`. |
| before support restrictions no slice is a consequence of earlier ones | **CONFIRMED** (S4.1/S4.2): P-only and Q-only cokernels are 0 at every level. |
| `#{new constraints} = dim coker(A_n L_n)` | **IMPLEMENTED LITERALLY** (S4), as an explicit matrix over `Q`; the canonical basis is the left kernel. |
| **the joint audit must STACK the two support maps** | **CONFIRMED, and it is the whole result** (S4.3). Separate counting gives `0 + 0 = 0`; stacked gives `2n-3`. The warning was not a technicality. |
| `t^(21-2i) \| [u^(12-i)]H^3` for `i <= 10` | **DERIVED and machine-checked** (S2.3.b), index and sign as stated. |
| the sub2 argument transfers to sub1 | **NO, and it did not need to.** The transferring object is not SPINE's forcing but the slice calculus itself, which is regime-parametric. `E_min` being vacuous in sub1 is irrelevant to it. |

Nothing in the brief's algebra was found wrong. The one correction of emphasis:
the obstruction that matters is **not** a bigger pile of P-slice conditions.
It is the single stacked family, and it needs **no cell input at all**.

---

## 1. The Q column, derived

`verify_derivation.py` §B establishes, and S2.1 recomputes from a generic unit
series, that `C = x^4 * (unit)`, so `C^-1` starts at `x^-4`; and `[QQ1]` gives
`v_{1,0}(F) = -5`, so `F` starts at `x^-5`. Hence

> for every `M >= -3`,  `Q_M = (C^3)_M`  exactly — no `lambda`, no `F`.

With `c_m = D_m*C4^(2m-7)`, `D_4 = 1`, `C4 = y^7*t` (the **same** in both
windows — the corners `(8,14),(8,16)` and `(12,21),(12,24)` are shared), and
the bridge's stripping `d_j = D_j / y^(12(4-j))`:

```
y :   7*(2M-21) + 12*(12-M)  =  2M-3
t :   2M-21
```

> ```
> Q_M  =  y^(2M-3) * [u^(12-M)] H(u)^3 / t^(21-2M)
> ```

`Q_M` is a polynomial and `gcd(y,t) = 1`, so `t^(21-2M) | [u^(12-M)]H^3` for
`M <= 10`. S2.3.b verifies the identity **exactly, on generic stripped `d`'s,
for every `M` in `[5,12]`** — the same standard `positive_slice.py` held itself
to on the P side.

In level coordinates `n = 12-M` (and `n = 8-M` on the P side) the two families
are, with `p_n := [u^n]H^2` and `r_n := [u^n]H^3`:

```
P :   t^(2n-2) | p_n     for n = 2..8 ;   p_n = 0 exactly for n >= 9
Q :   t^(2n-3) | r_n     for n = 2..15
```

`p_n = 0` for `n >= 9` because `P` has no negative `x`-powers. That is
*stronger* than any divisibility, and it is used as such.

---

## 2. Why stacking is not optional

`h_0 = 1`, so the fresh coefficient enters with

```
p_n = 2*h_n + q_n^P ,        r_n = 3*h_n + q_n^Q .
```

Both windows have `cap_n + 1 >= 2n-2` (`cap_n = lam*n`, `lam = 3` sub1 / `2`
sub2), so `h_n` surjects onto the forbidden P-jets **and** onto the forbidden
Q-jets. Each side alone is therefore **absorbable**:

| window | level `n` | `coker(P only)` | `coker(Q only)` | naive sum | **stacked** |
|---|---|---|---|---|---|
| sub1/sub2 | 2 | 0 | 0 | 0 | **1** |
| sub1/sub2 | 3 | 0 | 0 | 0 | **3** |
| sub1/sub2 | 4 | 0 | 0 | 0 | **5** |
| sub1/sub2 | 5..8 | 0 | 0 | 0 | **7,9,11,13** |

Stacked, `h_n` has already been spent satisfying the `2n-2` P-jets; the `2n-3`
Q-jets are then determined, and

```
2*3 - 3*2 = 0   =>   h_n cancels identically in  2*r_n - 3*p_n .
```

Since `2H^3 - 3H^2 = -1 + 3K^2 + 2K^3` with `K = H-1`, the canonical basis of
new obstructions at level `n` is the `2n-3` jets of

> ```
> t^(2n-3)  |  [u^n]( 3*K^2 + 2*K^3 )
> ```

involving `h_1..h_{n-1}` **only**. This is the entire new content, and it is
independent of the cell, the branch, and the window.

*Cheapest instance.* At `n = 2`: `[u^2](3K^2+2K^3) = 3*h_1^2`, so `t | h_1`,
i.e. **`d3(-1) = 0`**. Concretely, `P_6 = y^10*(2*d2+d3^2)/t^2` and
`Q_10 = 3*y^17*(d2+d3^2)/t` are both polynomials; subtracting gives `t | d3^2`.
`POSITIVE_SLICE.md` §3.3 explicitly leaves `eta := h(-1)` free — its controls
imposed P support only. The Q column pins it (S7.9).

---

## 3. The cascade

Every P condition is absorbable, so write

```
h_n = -q_n^P/2 + t^(2n-2)*g_n     (n <= 8) ,      h_n = -q_n^P/2   (n >= 9)
```

with `g_n` free. Every P condition then holds identically and the *only*
remaining content is the stacked family. Its lowest jet at each level is a
forced equation. `slice_obstruction_basis.py` S8 runs this and finds, at every
level, a jet that is a **unit times a single irreducible factor, linear in one
`g`-coefficient** — so each step is a forced consequence, never a choice of
component (S8.1):

| level `n` | required | lowest nonzero jet | forced |
|---|---|---|---|
| 2 | `t^1` | `t^0`: `3*g1_0^2` | `v_t(h_1) >= 1` |
| 3 | `t^3` | — (all vanish) | — |
| 4 | `t^5` | `t^4`: `(3/4)*(g1_1^2-2*g2_0)^2` | `v_t(h_2) >= 3` |
| 5 | `t^7` | — | — |
| 6 | `t^9` | `t^8`: `3*(...)^2` | `v_t(h_3) >= 5` |
| 7 | `t^11` | — | — |
| 8 | `t^13` | `t^12`: `(3/4)*(...)^2` | `v_t(h_4) >= 7` |
| 9 | `t^15` | — | — |
| **10** | `t^17` | `t^16`: `3*(... + g5_0)^2` | **`v_t(h_5) >= 9`** |

Uniformly: **level `2m` forces the `t^(2m-2)` coefficient of `h_m` to vanish**,
advancing `v_t(h_m)` from `2m-2` to `2m-1`. Odd levels contribute nothing
(S8.2). Each jet is a perfect square, which is why no case split arises.

`--deep` stops at level 10, the level that bounds `h_5 = e` and therefore
produces the census delta. `--deep12` additionally runs levels 11 and 12: level
11, like every odd level, contributes nothing; **level 12 fires in the same
shape and gives `v_t(h_6) >= 11`** (jet `t^20`, again a perfect square, again
linear in `g6_0`). Levels 14 and 16 (which would advance `h_7 = S`, `h_8 = T`)
were **not computed**; the bounds reported for those rows are the un-advanced
`v_t(h_7) >= 11`, `v_t(h_8) >= 12` and are not claimed sharp.

### 3.1 The bridge to the census

`window_caps_verify.py` W3's shift map uses **generalized** binomials
`binom(m, m-j)`, and `binom(m, m-j) = 0` whenever `m >= 0 > j`. So the
`d3`-killing shift is **triangular across zero**: no non-negative `D_m` feeds
any spare. In particular (S3.4)

> ```
> D*_{-1} = D~_{-1} ,     i.e.    h_5 = dm1 = e   exactly.
> ```

(`D*_{-2} = D~_{-2} - theta*D~_{-1}` already mixes, so this is a statement about
level 5 specifically, and it is the level that matters.)

Since `a_t = v_t(e)` by definition (`divisor_filter.py`: `e = gamma*t^a_t*prod
(y-r_i)^{b_i}` with the off-support factor a unit at `t`),

> ```
> a_t  >=  9 .
> ```

### 3.2 The cascade constrains; it does not contradict

Under `h_k = t^(2k-1)*(free)` the substitution `u = v/t^2` gives `K = Hhat(v)/t`
and

```
[u^n](3K^2+2K^3) = 3*t^(2n-2)*[v^n]Hhat^2 + 2*t^(2n-3)*[v^n]Hhat^3 ,
```

divisible by `t^(2n-3)` for **every** `n`. So the forced profile satisfies all
stacked conditions identically (S8.4b): the cascade pins valuations, it does not
empty the slice system. The bounds are also **sharp** — the joint control of §4
attains `v_t(h_k) = 2k-1` exactly for `k = 1..4` (S8.4).

---

## 4. Controls (non-negotiable — all pass)

**S6a — P side, genuine polygon-supported `P`.** Corners loaded from
`paper_src/upstream_facts.json`, `P_8 = C4^2` forced, every other slice random
inside the hull; four seeds. The D-recursion output meets the certified caps,
stripping is legal, the slice formula reproduces **all nine** slices
`P_0..P_8` exactly, and **every** P-side condition `t^(2n-2) | p_n` holds
identically for `n = 2..8` — not just the three `POSITIVE_SLICE.md` used.

**S6b — the JOINT control.** `C = x^4*C4 + a*y^5*x^3 + b*y^3*x^2 + c*y*x + e0`
is a Laurent **polynomial** in `x`, so `P = C^2` and `Q = C^3` are polynomials
in both variables with `lambda = F = 0`: a genuine point of the slice system for
every `(a,b,c,e0)`. Its stripped coordinates are `h_1 = a*t`, `h_2 = b*t^3`,
`h_3 = c*t^5`, `h_4 = e0*y*t^7` — polynomials, `h_0 = 1`. On it, **every**
obstruction functional — P-side, Q-side, and stacked — vanishes **identically in
`(a,b,c,e0)`, at every level `n = 2..8`**. The obstruction does not fire on
genuine data.

**S6c — MUTATION.** Perturbing the single support coefficient `h_1 = d3` by `+1`
makes the *predicted* level-2 stacked functional go nonzero
(`3*p_2-2*r_2 = -3*h_1^2`, residue `-3` mod `t`), while the same functional on
the unmutated instance is identically zero. The obstruction detects the mutation
and nothing else.

Both controls abort the run with exit 1 before the obstruction is reached.
Neither does.

**Limitation, stated plainly.** The joint control has `h_5..h_8 = 0` (any `C`
finite in `x` has zero spares, and any `C` with nonzero spares is an infinite
series whose `Q` conditions are exactly what is at issue). So the level-10 step
is controlled structurally — S8.4b (the profile satisfies it), S8.4c (`g5_0`
enters linearly with a unit coefficient, so it is a genuine equation on `e`) —
but **not** by an instance with `e != 0`. That is the weakest joint in this file
after `[QQ1]`.

---

## 5. Regression against `POSITIVE_SLICE.md`

Required: reproduce the three conditions and the constant-term contradiction on
`a10_b0000_T1`. Done, end to end, without importing `positive_slice.py` or
`spine.py`:

* the `n = 0` rows are re-derived from `generators.json` and factor exactly as
  SPINE's `g1,g2,g3,kbox` (S7.2);
* the forced `y = -1` values match `POSITIVE_SLICE.md` §5.2 exactly (S7.3);
* the inverse-shift formulas match §3.2 (S7.4), and levels `n = 2,3,4` are its
  slices `M = 6,5,4` (S7.5);
* **the count is re-derived, not assumed** (S7.6): with `V_n` free the depth-1
  P-only cokernel is `0`; with SPINE's `y=-1` value **pinned** it is `1` at each
  of `n = 2,3,4` — *exactly three* constant-term conditions;
* they reproduce `(A) = 7Y^2 - X(48X+8)` exactly (S7.7);
* their ideal, saturated at `gamma != 0`, is the **unit ideal** over `Q` (S7.8).

And the new machinery kills the same cell a **second, shorter way**: the stacked
level-2 functional gives `eta = 0` with no cell input, whereupon `(A)` collapses
to `-8X(6X+1) = 0` and dies against SPINE's forced `delta2` (S7.9–S7.11).

---

## 6. Frontier impact (READ-ONLY census)

`frontier_rebuild.STAGES` is read **by AST** — never imported, never executed,
never modified. The stage-2 universe reproduces `FRONTIER_REBUILD.md` on the
nose (34 / 314 / 7275 and 34 / 322 / 8889), which is itself a control.

Criterion: **`a_t < 9` => EMPTY**.

| window | C08/C20 | cells | flag cases | states |
|---|---|---:|---:|---:|
| sub1 | ON | 34 -> **11** | 314 -> **134** | 7275 -> **3657** |
| sub1 | OFF | 34 -> **11** | 322 -> **142** | 8889 -> **4696** |
| sub2 | ON | 5 -> 2 | 70 -> 34 | 2947 -> 1520 |
| sub2 | OFF | 5 -> 2 | 73 -> 37 | 3018 -> 1559 |

**sub1 killed (23):** `a2_b1111_T1`, `a3_b1000_T1`, `a3_b1110_T1`,
`a4_b0000_T1`, `a4_b1100_T1`, `a4_b1111_T1`, `a5_b1000_T1`, `a5_b1110_T1`,
`a6_b0000_T1`, `a6_b1100_T1`, `a6_b1111_T1`, `a6_b1111_T2`, `a7_b1000_T1`,
`a7_b1100_T1`, `a7_b1110_T1`, `a7_b1111_T1`, `a8_b0000_T1`, `a8_b0000_T2`,
`a8_b1000_T1`, `a8_b1100_T1`, `a8_b1110_T1`, `a8_b1111_T1`, `a8_b1111_T2`.

**sub1 surviving (11):** all of `a9_*` (5) and `a10_*` (6).

The sub1 row is a genuine frontier delta, because `stage3_spine` and
`stage4_positive_slice` kill nothing in sub1. **The sub2 row is not**: standard
sub2 is already EMPTY after stage 4. It is measured against the *stage-2*
universe on purpose, so that §6.1's corroboration is against cells this lane
could have contradicted.

The kill is **identical under C08/C20 ON and OFF**: `a_t >= 9` is a `t`-adic
valuation statement over `Q` with no square class, no splitting field and no
residue arithmetic, so the field-scope downgrade has no purchase on it. (The
net frontier count is not this file's to quote.)

### 6.1 Cross-corroboration — the strongest control available

Nothing here uses SPINE. Yet on sub2 the criterion independently re-kills
`a6_b1111_T1`, `a7_b1110_T1`, `a8_b1100_T1` — **three of SPINE's seven sub2
cells**, by a completely different mechanism (SPINE: a zero-slack degree count
on the G-rows; here: slice polynomiality alone). It kills nothing SPINE does not
also kill. *"No survivors" is exactly the shape a bug takes*, so agreeing with
kills nobody derived this way is real evidence. It also does not touch
`a9_b1000_*` (`a_t = 9`), which is the correct behaviour — the criterion is a
threshold, not a blanket.

Consistency with sub2's own arithmetic: there `deg e = 10 = a_t + sum(b_i)`
exactly, so `a_t >= 9` forces `(a,sum b) in {(9,1),(10,0)}` — precisely the two
sub2 cells that survive.

### 6.2 The alternate regime — a real negative result

All six surviving alternate T1 branches (`ALT_FRONTIER_V2.md`) have
`a_t in {12,14} >= 9`. **The criterion kills nothing there.** The `a_t` bound is
the wrong shape for `a_t >= 11`.

### 6.3 A conditional handover to the alternate lane — NOT claimed here

Level 12 **was** computed (`--deep12`) and gives `v_t(h_6) >= 11`. The inverse
shift gives `D*_{-2} = dm2 - (h/4)*dm1` with `v_t(h) >= 1` and `v_t(dm1) = a_t`,
so

```
v_t(R)  >=  min( 11 , 1 + a_t ) ,      i.e.   t^11 | R   whenever a_t >= 10.
```

On the six alternate branches `a_t in {12,14}`, so `t^11 | R`. `T1_BRANCH.md`'s
place trichotomy at `beta = -1` has two horns:

* **(H1)** `v_t(R) >= v_t(e) = a_t` — `12` or `14`, both `>= 11`: **no contradiction**;
* **(H2)** `30 = a_t + 2*v_t(R)`, i.e. `v_t(R) = (30-a_t)/2 = 9` (`a=12`) or `8`
  (`a=14`) — both **below 11**: **contradiction**.

So every alternate branch sitting on horn 2 is EMPTY by this bound. The ALT lane
reports that the odd-`a` branches die precisely because H2 needs `(30-a)/2` in
`Z`, which suggests H1 is already excluded there — **if that is so, all six
alternate branches die.**

**This lane does not claim it.** The horn premise and the H1 exclusion belong to
`ALT_FRONTIER_V2.md` / `T1_BRANCH.md`, are not re-derived here, and the
alternate entry of the stage record is deliberately left **empty**. This is a
handover, not a result.

---

## 7. PROVED / CHECKED / INFERRED

**PROVED here** — exact polynomial identity, exact linear algebra over `Q`, or
exact Gröbner/resultant, machine-checked, over the premises below:

* the Q-slice formula `Q_M = y^(2M-3)[u^(12-M)]H^3/t^(21-2M)`, and hence the
  Q conditions `t^(2n-3) | r_n` (§1);
* the cokernel table: P-only `0`, Q-only `0`, stacked `2n-3`, both windows, all
  levels (§2);
* the identity `2H^3-3H^2 = -1+3K^2+2K^3` and the exact cancellation of the
  fresh coefficient (§2);
* the shift's triangularity across zero, hence `h_5 = dm1` (§3.1);
* the cascade, every step a forced single-component deduction, hence
  `v_t(h_k) >= 2k-1` for `k = 1..5` (§3), and `v_t(h_6) >= 11` with `--deep12`;
* the profile's satisfiability and sharpness (§3.2);
* the controls and the mutation control (§4);
* the full `POSITIVE_SLICE.md` regression including the unit-ideal kill (§5);
* the census arithmetic, read-only (§6).

**CHECKED** — reproduced from an existing artifact without re-proving it:

| tag | statement | source |
|---|---|---|
| [Q1] | canonical `G1,G2,G3,G5body` | `generators.json` (loaded, never transcribed) |
| [Q2] | `Phi = c*t^30*q`, `c = -1/6630`, `q(-1) = 3315` | `verify_derivation.py` §A |
| [Q3] | `ord D_j >= 12k`, `deg D_j <= (12+lam)k`, `D_j = C_j*C4^(7-2j)`, `C4 = y^7*(y+1)` in **both** windows | `window_caps_verify.py` W2/W5 |
| [Q4] | the `d3`-killing shift and its D-coordinate form | `window_caps_verify.py` W3 |
| [Q7] | the Prop-4.3 sub1/sub2 corner sets | `paper_src/upstream_facts.json` |
| [Q8] | the G-system indeterminates are the **shifted** stripped `D~_j` | convention; `POSITIVE_SLICE.md` §3.3 |
| [QQ1] | `Q = C^3 + lambda*C^-1 + F`, `v_{1,0}(F) = -5` | `PROOF_INVENTORY.md` premise **C3 (2/4)**; `T6_PREMISES.md` §2 |
| [QC1] | `a_t = v_t(e)`, `e = dm1`, `deg e = a_t + sum b_i` | `divisor_filter.py`; `phase_d_states` schema |

**INFERRED** — nothing is asserted as a result. The one forward-looking remark
(the level-12 alternate-regime lead, §6.2) is labelled as not computed.

**The load-bearing imports, stated plainly.**

* **`[QQ1]` is the hinge — but it is not a new one.** The entire Q column, and
  therefore the stacking and everything here, rests on `Q_M = (C^3)_M` for
  `M >= -3`, i.e. on the `alpha`-strip WLOG (`T6_PREMISES.md` Premise 2,
  status *READY-WITH-CITATION*, GGHV22 §4 template; `PROOF_INVENTORY.md` C3 at
  **2/4**). Note the normalized form is the one that *already* underpins the
  whole pipeline: `STATE.md` §27 states it, and `verify_derivation.py` §D checks
  `regenerate_system.py`'s `D3(j)` rows against `(C^3)_{-j}` under exactly this
  premise. So `[QQ1]` is not an extra assumption this lane smuggles in — if it
  fell, the G-system itself would fall, not merely this file.

  What *is* new here is using the **non-negative** `x`-slices of `Q`; the
  pipeline only ever used `j = 4, 5` and below. That extension needs
  `v_{1,0}(F) = -5` to bound *all* of `F`, not just its leading form — and it
  does: `v_{1,0}` is the top `x`-degree (`verify_derivation.py` §A uses
  `ell(P) = x^8*C4^2` and `ell(2*C^3*F) = 2*x^7*C4^3*F_{-5}`), and `F` lives in
  `K[y,C4^-1]((x^-1))`, so `v_{1,0}(F) = -5` says `F` has **no** term of
  `x`-degree above `-5`. Likewise `C = x^4*(unit)` gives `C^-1` top degree
  `-4`. Both correction columns therefore start strictly below `M = -3`, which
  is exactly what S2.1 checks. The premise is used at its stated strength, not
  beyond it.
* **`[Q8]`/S3.4 is what makes it a census delta.** Without `h_5 = dm1` the
  cascade is a statement about unshifted coefficients and kills nothing.
  S3.4 proves the triangularity from the map `window_caps_verify.py` W3 itself,
  but the identification of the G-system variables with `D~` remains the
  convention `POSITIVE_SLICE.md` §3.3 flags.

---

## 8. What this does NOT do

1. **Nothing is entered into the ledger or the DAG.** No `state_kill_ledger.json`,
   no `proof_dag.json`, no `phase_d_states*.json`, no `frontier_rebuild.py` was
   written. §6 is a read-only census and has **not been audited**. Evidence
   grade is *exact-checked (same-author)*, **not** independently audited — the
   result needs a second, independently authored checker before it is load-bearing.
2. **The alternate regime is untouched** (§6.2), and the level-12 lead is a lead.
3. **The joint control has zero spares** (§4). The level-10 step, which is the
   one that produces `a_t >= 9`, is controlled structurally but not by an
   instance with `e != 0`.
4. **Only the stacked family is used.** The `y`-order conditions (`y^2 | p_8`,
   `y | r_11`, `y^3 | r_12`), the degree/upper-hull conditions, the exact
   equations `r_n = 0` for `n >= 16`, and the sharp relation
   `[u^17]H^3 = (1/6630)*t^30*q` (equivalent to `Phi = -y^204*[u^17]H^3`) are
   **all unused**. There is slack here that was never spent.
5. **The sub2 degree caps are not exploited.** In sub2 `deg g_n <= 2`, which is
   far tighter than the freedom the cascade was given. The cascade deliberately
   over-allows `g_n` (conservative), so its deductions hold a fortiori — but a
   sub2-specific run would be stronger.
6. **No Gröbner basis over the G-system, no modular arithmetic, no solver.**
   Nothing external was run, so there are no aborts, timeouts, or exit codes.

---

## 9. Reproduce

```
cd d2_plane_72_108
python -u slice_obstruction_basis.py --deep            # full derivation + report
python -u slice_obstruction_basis.py --quiet --deep    # exit 0 iff every check passes
python -u slice_obstruction_basis.py                   # levels <= 8 only (fast; no a_t bound)
python -u slice_obstruction_basis.py --deep12          # + levels 11,12 (slow; v_t(h_6) >= 11)
```

`--deep` (levels 2..10) is what produces `a_t >= 9` and the census delta;
`--deep12` adds only the §6.3 handover and changes no verdict here.

Read-only and pure sympy. The stage record is emitted as a **drop-in**,
`slice_obstruction_stage.json`, in `frontier_rebuild.STAGES`' schema —
`frontier_rebuild.py` is deliberately not edited, because another lane owns it:

```json
{"id": "stage5_slice_obstruction",
 "title": "Stacked P/Q positive-slice obstruction (a_t >= 9)",
 "checker": "python slice_obstruction_basis.py --quiet --deep",
 "dead": {"sub1": [23 cells], "sub2": ["a6_b1111_T1","a7_b1110_T1","a8_b1100_T1"]},
 "applies_after": "stage4_positive_slice"}
```

Until the owning lane appends it and re-runs `frontier_rebuild.py`, the frontier
artifacts do not reflect this result. That is a deliberate, flagged omission.
