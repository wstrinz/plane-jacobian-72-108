# Phase C residue-system worklist — f31 subcase-(2) survivors

> **STATUS (2026-07-22):** Headline counts below are superseded. The current sub2 frontier is **26 cells / 220 flag cases** (`cascade_cones_qt_inf_rl.json`, q+t+inf with residue kills and the T2 squeeze), not the "30-branch / 232-case" figure quoted here. The §4 level-5-squeeze claim was corrected by `T5_T2_COLUMN.md`: under partial q-support the conclusion weakens to `F²|G` (F the q-coprime remainder); 4 T2 cells killed, 8 open. Authoritative artifacts: `cascade_cones_qt_inf_rl.json`, `phase_d_states_sub2.json`, `T5_T2_COLUMN.md`, `STATE.md`. The pattern inventory and residue equations below remain valid.

> **Correction (2026-07-22, later same day).** Section 4's judgment that the level-5 squeeze clears all 12 T2 cells was falsified by execution: under partial q-support the coprimality conclusion weakens to `F^2 | G` and only 4 cells died (`T5_T2_COLUMN.md`); the remaining 8 were narrowed by `T5_T2_INFINITY.md`. The sub2 frontier is 26 cells, not 30. Pattern inventory and residue equations remain valid.

**Date:** 2026-07-22
**Scope:** the 30-branch / 232-case frontier of `cascade_cones_qt.json`
(q+t coupled, depth 4), the Phase C input described in
`CASCADE_ENGINE_REPORT.md` §"What the 30 survivors need".

**How to read this file.** Statements tagged **[data]** are mechanical
consequences of the witness records in `cascade_cones_qt.json` (re-derivable
by `phase_c_inventory.py`, which also re-checks every tied-monomial string
against the `h_l` tables). Statements tagged **[judgment]** are my
assessment of Phase C difficulty and lemma reuse, not forced by the data.
Every algebraic monomial is quoted from the `cascade_signature.py` /
`f31_graded.txt` `h_l` tables (reproduced in §0); obligation *semantics* are
quoted from `CASCADE_ENGINE_REPORT.md` §Semantics; local-descent *style* is
taken from `T5_90_T1.md` and `T5_STRATA_50_11.md`.

Provenance note **[data]**: `phase_c_inventory.py` confirms all 3621 tied
monomials across the 232 cases are exact members of the `h_l` tables — the
tied sets are always sub-multisets of a single `h_l` cut down by the case's
zero flags.

---

## 0. The `h_l` monomial tables used below

From `cascade_signature.py` (rewrite `sigma = 4*d0 - d2**2`, variables
`(d2,d1,sigma,e)`, `e = d̃₋₁`); levels 4/5/6/7 are the ones the depth-4
engine touches. Verbatim:

```
h_7 = 8192*d1^2
h_6 = 14336*d1^2*d2 + 8192*d1*e - 3072*sigma^2
h_5 = -12288*d1^2*d2^2 + 32256*d1^2*sigma + 18432*d1*d2*e - 9216*d2*sigma^2 + 2048*e^2
h_4 = -220752*d1^4 - 31232*d1^2*d2^3 - 23616*d1^2*d2*sigma - 3072*d1*d2^2*e
      + 34560*d1*e*sigma - 5184*d2^2*sigma^2 + 5632*d2*e^2 - 12096*sigma^3
```

The cascade identity (soundness ground truth, `CASCADE_ENGINE_REPORT.md` §Semantics,
`t5_multiplace_verify.py` checks 5–7) is, level by level,

```
t^v g_{l+1} = ehat^3 g_l + u^l h_l,   u = c q,   v = 30 - 3a,   c = -1/6630,
q = 2048y^4 - 512y^3 + 320y^2 - 240y + 195,   t = y+1,
T1 terminal: ehat^3 g_7 = -u^7 (8192 d1^2)
T2 terminal: ehat^3 g_6 = -u^6 (-3072 sigma^2)
deg g_l <= 10 + 3a.
```

Local reading at a q-root `p` (uniformizer `pi = y - p`): `t,q'` are units,
`v_p(u) = 1`, `v_p(e) = b`. At the place `t = y+1`: `v_t(u) = v_t(ehat) = 0`
(`q(-1)=3315`), the e-slot of each monomial costs `a = v_t(e)`.
For a quantity `X` with `v_p(X) = m` write `X* := (X / pi^m) mod pi` for its
leading residue.

Obligation kinds (`CASCADE_ENGINE_REPORT.md` §Semantics), i.e. the four ways
a level identity can be satisfied without a forced valuation:

- **`monomial_tie_rise`** — the tied minimum monomials of `h_l` share the
  minimum valuation `M`; to let `v_p(h_l)` rise to `M + depth`, the leading
  residues of the tied monomials must cancel to that depth.
- **`term_cancellation`** — the two *sides* `ehat^3 g_l` and `u^l h_l` are
  valuation-tied and must cancel to `depth` (uses the free leading residue of
  `g_l`; `tied=[]`).
- **`exact_identity`** — `g_{l+1} ≡ 0` forces the polynomial identity
  `ehat^3 g_l = -u^l h_l` (recorded `depth 0`; `tied=[]`).
- **`identical_vanishing`** — `h_l` is the zero polynomial at this place's
  surviving monomials (`depth 0`; needs ≥2 tied monomials).

---

## 1. Pattern inventory  [data]

Parsed by `phase_c_inventory.py` over the 232 survivor cases (30 branches).
Grouping key = `(place kind, cascade level, obligation kind, tied-monomial set)`.

- **41 distinct obligation patterns.**
- Raw obligation counts by kind: `term_cancellation` 1178, `monomial_tie_rise`
  635, `exact_identity` 435, `identical_vanishing` 260 (total 2508).
- By place: `q` 1925, `t` 583. By level: L4 921, L5 902, L6 685.

### 1.1 The 12 highest-frequency patterns

| # | place | L | kind | tied N | freq | cells | depth multiset |
|--:|:--|--:|:--|--:|--:|--:|:--|
| 1 | q | 6 | term_cancellation   | 0 | 344 | 17 | 1×344 |
| 2 | q | 5 | term_cancellation   | 0 | 299 | 29 | 1×299 |
| 3 | q | 4 | term_cancellation   | 0 | 241 | 29 | 1×241 |
| 4 | q | 5 | exact_identity      | 0 | 180 | 15 | 0×180 |
| 5 | q | 4 | exact_identity      | 0 | 168 | 19 | 0×168 |
| 6 | q | 6 | monomial_tie_rise   | 3 | 132 | 14 | 1×132 |
| 7 | t | 6 | term_cancellation   | 0 | 110 | 17 | 3–15 |
| 8 | t | 5 | term_cancellation   | 0 |  97 | 29 | 3–15 |
| 9 | t | 4 | term_cancellation   | 0 |  87 | 29 | 3–18 |
| 10| q | 5 | monomial_tie_rise   | 5 |  84 | 17 | 1×84 |
| 11| q | 4 | monomial_tie_rise   | 8 |  77 | 17 | 1×77 |
| 12| q | 4 | identical_vanishing | 8 |  68 | 17 | 0×68 |

Structural facts **[data]**:

- The `term_cancellation` patterns (#1–3, 7–9) carry no tied monomials: they
  constrain the *free leading residue of `g_l`*, not the window data. At `q`
  they are all depth 1; at `t` they run 3…18 (scaling with `a` and the
  monomial e-slot — see §5).
- `exact_identity` patterns (#4,5) are the `g_{l+1}≡0` book-keeping: every
  survivor case that zeroes a `g_l` carries one exact_identity per place at
  the level below.
- The *window-only* obligations — the ones that are genuine residue equations
  in `(d2,d1,sigma,e)` with no free `g` — are exactly the `monomial_tie_rise`
  and `identical_vanishing` families (895 of 2508). Their tied sets are always
  the surviving cut of a single `h_l`:
  - full `h_6` (3 mon), full `h_5` (5 mon), full `h_4` (8 mon) at b=0 T1
    places (nothing zeroed);
  - `d1≡0` cuts at T2 / local-`d1`-vanishing places:
    `h_6|_{d1≡0} = -3072 sigma^2` (a *single* monomial ⇒ never a tie),
    `h_5|_{d1≡0} = -9216 d2 sigma^2 + 2048 e^2`,
    `h_4|_{d1≡0} = -5184 d2^2 sigma^2 + 5632 d2 e^2 - 12096 sigma^3`.

### 1.2 Which cells each window-pattern touches  [data]

The three "full `h_l`" residue patterns are the backbone of the T1 column:

- **full `h_6` tie** `{14336 d1^2 d2, 8192 d1 e, -3072 sigma^2}` (pattern #6,
  132 occ) touches 14 T1 cells:
  `a10 b0000, a9 b0000, a9 b1000, a8 {0000,1000,1100}, a7 {1000,1100,1110,3000},
  a6 {1000,1100,1110}, a5 b1110` (all T1).
- **full `h_5` tie/vanish** (patterns #10,#12-analog, 84+61 occ) touches the
  same 17 T1 cells (adds `a6 {3000,3100}`, `a5 b3110`).
- **full `h_4` tie/vanish** (patterns #11,#12, 77+68 occ) touches those 17 T1
  cells.

The `d1≡0` cuts are the T2 column plus locally-`d1`-vanishing T1 cases:

- `h_5|_{d1≡0} = -9216 d2 sigma^2 + 2048 e^2` appears as an exact tie/vanish
  set in 10 cells (`a9 b1000 T2` and 9 T1 cells where `d1` vanishes at a
  root); the *mechanism* however applies to all 12 T2 cells (§4).
- `h_4|_{d1≡0}` cut `{-5184 d2^2 sigma^2, 5632 d2 e^2, -12096 sigma^3}`
  (30 occ) touches all 12 T2 cells.

---

## 2. Explicit residue equations for the top window patterns

Convention: an obligation of `depth D` on tied set `S ⊂ h_l` at place `p`
asserts that the polynomial `Σ_{m∈S} m`, evaluated at the leading residues
`(d2*,d1*,sigma*,e*)` of the window data at `p`, vanishes to order `D` in
`pi`. For a b=0 q-root every slot is a unit, so `depth 1` is literally
"the cut of `h_l` vanishes at the point `p`". This is the exact analogue of
the "alternating square" residue equations of `T5_90_T1.md` §2 (there the
successive `[x^k]N_l` coefficients are forced to be perfect squares that must
vanish); here the engine has pre-extracted the tied leading forms.

### P6 — full `h_6` tie, q-place, depth 1 (132 occ, 14 T1 cells)  [data]

At a b=0 root `p` (all slots units), the level-6 line needs `v_p(h_6) ≥ 1`:

```
14336 (d1*)^2 (d2*) + 8192 (d1*)(e*) - 3072 (sigma*)^2 = 0.          (P6)
```

One scalar equation per b=0 root; the same polynomial `h_6` is evaluated at
all four roots of the *fixed* `q`, so (P6) is really one condition on the
window polynomials `d2(y),d1(y),sigma(y),e(y)` reduced mod `q` — a resultant-
level coupling, not four independent points **[judgment]**.

### P10 — full `h_5` tie, q-place, depth 1 (84 occ, 17 T1 cells)  [data]

```
-12288 (d1*)^2(d2*)^2 + 32256 (d1*)^2(sigma*) + 18432 (d1*)(d2*)(e*)
   - 9216 (d2*)(sigma*)^2 + 2048 (e*)^2 = 0.                         (P10)
```

### P11 — full `h_4` tie, q-place, depth 1 (77 occ, 17 T1 cells)  [data]

```
-220752 (d1*)^4 - 31232 (d1*)^2(d2*)^3 - 23616 (d1*)^2(d2*)(sigma*)
   - 3072 (d1*)(d2*)^2(e*) + 34560 (d1*)(e*)(sigma*)
   - 5184 (d2*)^2(sigma*)^2 + 5632 (d2*)(e*)^2 - 12096 (sigma*)^3 = 0. (P11)
```

Pattern #12 is the same left-hand side as (P11) but with `identical_vanishing`
semantics: it fires when the case has already zeroed `g_5` (so the level-4
line is `ehat^3 g_4 = -u^4 h_4`) and forces `h_4 ≡ 0` as a polynomial at the
place, i.e. *all eight* leading residues cancel simultaneously, not merely the
sum. That is a strictly stronger demand than (P11) **[judgment]**.

### P-squeeze — `h_5|_{d1≡0}` tie, q-place, depth 1 (T2 backbone)  [data]

At a b=0 root of a T2 (or local-`d1`-vanishing) case:

```
-9216 (d2*)(sigma*)^2 + 2048 (e*)^2 = 0.                            (Psq)
```

This is exactly `h_5|_{d1≡0} = -9216 d2 sigma^2 + 2048 e^2` (`T5_STRATA_50_11.md`
[W2]), the input to the level-5 squeeze (§4). Globally, when `g_5 ≡ 0` the
level-5 line reads `t^v g_6 = u^5 h_5` with `u = cq`, forcing `q^5 | g_6` and
hence `q | (-9216 d2 sigma^2 + 2048 e^2)` — the polynomial form of (Psq)
across all four roots at once.

---

## 3. Depth histogram and cheapest targets

### 3.1 Per-cell difficulty  [data]

Ranking the 30 cells by the total obligation depth of their *cheapest*
survivor case (`phase_c_inventory.py`; `q_dep`/`t_dep` split; `d0` = number
of depth-0 exact_identity+identical_vanishing obligations in that case):

| rank | cell | total | q_dep | t_dep | d0 |
|--:|:--|--:|--:|--:|--:|
| 1 | a9 b=1000 T2  | 3  | 3 | 0  | 5  |
| 2 | a10 b=0000 T1 | 4  | 4 | 0  | 10 |
| 3 | a9 b=1000 T1  | 6  | 3 | 3  | 10 |
| 4 | a9 b=0000 T1  | 7  | 4 | 3  | 10 |
| 5 | a8 b=1100 T1  | 8  | 2 | 6  | 10 |
| 6 | a8 b=1000 T1  | 9  | 3 | 6  | 10 |
| 7 | a7 b=1110 T1  | 10 | 1 | 9  | 10 |
| 8 | a8 b=0000 T1  | 10 | 4 | 6  | 10 |
| … | (T1 a5–a7)    | 11–15 | 0–3 | 9–12 | 10 |
| 16+| T2 a5–a8 + T1 a5/a6 b=3xxx | 16–36 | 2–8 | 12–30 | 0/5 |

Depth histograms of individual obligations **[data]**:
`q`: depth 0 → 556, depth 1 → 1359, depth 2 → 4, depth 5 → 6.
`t`: depth 0 → 139, then 3→63, 6→101, 9→116, 11→3, 12→122, 15→33, 18→6.

So **all q-place residue depth is 0 or 1** (six depth-5 outliers aside); the
heavy numbers live entirely at the `t` place, and the T2 column's totals are
inflated by their deep `t` term-cancellations.

### 3.2 The depth metric is a proxy  [judgment]

Total depth *undercounts* `exact_identity`/`identical_vanishing` (both scored
0), which are the strongest obligations (a full polynomial identity, resp.
`h_l ≡ 0`). The `d0` column exposes this: the cheap T1 cells carry 10 depth-0
obligations each. Read the ranking as "cheap in the tie/cancellation sense";
weight the `d0` count when picking real effort.

### 3.3 Recommended first targets  [judgment]

1. **`a9 b=1000 T2`** (total 3, `t_dep`=0) — cheapest, and lands squarely on
   existing machinery (§4). Full sketch below.
2. **`a10 b=0000 T1`** (total 4, `t_dep`=0) — the degenerate top of the rigid
   tail; no root divides `e`, so it is pure "full-`h_l`" residue equations
   (P6/P10/P11) at the four q-roots with no b>0 complications. Its `t`
   obligations are all depth 0.
3. **`a9 b=1000 T1`** (total 6) — pairs with the `a=9` local-descent program
   already built in `T5_90_T1.md` / `t5_90t1_local_verify.py` (same stratum,
   T1 branch); the doc's nonconstant shape (4)/(5a) is the continuation.

### 3.4 Full residue system for the single easiest cell (`a9 b=1000 T2`)

All valuations below are **[data]** (the min-depth survivor case, `g_5≡0`,
`sigma≠0`, `d2≠0`); the reduction to a global divisibility is **[judgment]**
following the `T5_STRATA_50_11.md` §2 template.

Global setup: `a=9 ⇒ v=30-27=3`; T2 ⇒ `d1 ≡ 0`, `g_7 = 0`, terminal
`ehat^3 g_6 = 3072 u^6 sigma^2` (`u=cq`); `deg g_l ≤ 37`. `b=(1,0,0,0)`: one
q-root `p_1 | e` with `v_{p1}(e)=1`, roots `p_2,p_3,p_4 ∤ e`; place `t=y+1`
with `v_t(e)=9`. From `deg e ≤ 10`, `v_t(e)=9`, `v_{p1}(e)=1`:
`e = κ·t^9·(y-p_1)` with `p_1` a root of `q`.

Witness valuations (`v(g_4),g_5≡0,v(g_6)`) and obligations, per place:

| place | v(d2) | v(sigma) | v(g_4) | v(g_6) | obligations |
|:--|--:|--:|--:|--:|:--|
| `p_1` (b=1) | 0 | 2 | 3 | 7 | L4 exact_identity |
| `p_2,p_3,p_4` (b=0) | 0 | 0 | 4 | 6 | L5 tie_rise d1 **(Psq)**; L4 exact_identity |
| `t` (v_t e=9) | 3 | 0 | 0 | 0 | L4 exact_identity |

The residue system to refute:

- **(T)** terminal: `ehat^3 g_6 = 3072 c^6 q^6 sigma^2`, `deg g_6 ≤ 37`.
- **(E4)** `g_5 ≡ 0` ⇒ level-4 exact identity `ehat^3 g_4 = -u^4 h_4` at every
  place (`h_4|_{d1≡0} = -5184 d2^2 sigma^2 + 5632 d2 e^2 - 12096 sigma^3`).
- **(E5)** `g_5 ≡ 0` ⇒ level-5 line `t^3 g_6 = u^5 h_5`,
  `h_5|_{d1≡0} = -9216 d2 sigma^2 + 2048 e^2`. Since `v_p(g_6)=6` at each b=0
  root and `v_p(u^5)=5`, this forces `v_p(h_5)=1` there — the three depth-1
  (Psq) equations `-9216 d2(p_i) sigma(p_i)^2 + 2048 e(p_i)^2 = 0`.
- Consolidated: `t^3 g_6 = c^5 q^5 h_5` forces `q^5 | g_6`, hence
  **`q | (-9216 d2 sigma^2 + 2048 e^2)`** (the four-root form of (Psq)).

Refutation route **[judgment]**: substitute `sigma^2` from (T) into (E5) to
kill the `sigma^2` term (the T5_STRATA §2 "absorb into the ehat^3-part" move),
obtaining an `ehat^2 | ĝ`-type divisibility; with `deg e = 5`-class rigidity
this collapses `ehat` toward a constant and `sigma` toward a constant, at
which point Proposition E (infinity, §5) dominates. I.e. the *local* part of
this cell is completely covered by the level-5 squeeze; only the final degree
count is Phase D.

---

## 4. Reusable-lemma candidates

### L-squeeze (level-5 squeeze) — clears the entire T2 column  [judgment, precedent proven]

`T5_STRATA_50_11.md` §2 (Lemmas S1/S2) proves, for `d1≡0, sigma≢0`, the
level-5 line `t^v g_6 = ehat^3 g_5 + u^5 h_5` with
`h_5|_{d1≡0} = -9216 d2 sigma^2 + 2048 e^2` yields `ehat^2 | ĝ` after the
terminal `sigma^2` is absorbed. The doc's §4 states explicitly that "the
level-5 squeeze applies **verbatim** to every T2 branch (d1≡0, sigma≢0) of
every surviving stratum with a≤9". **Match [data]:** all 12 T2 survivor cells
(`a5 b1000`, `a6 {1000,1100,1110}`, `a7 {1000,1100,1110,3000}`,
`a8 {0000,1000,1100}`, `a9 b1000`) carry the (Psq) obligation (pattern
P-squeeze / its `h_4|_{d1≡0}` companion), and its hypotheses are met by
construction. One lemma, tuned only by the growing cap
`deg ĝ ≤ (10+3a) - 4(6-3a_q)`, reduces every T2 cell. This is the single
highest-leverage Phase C lemma.

### L-h6vanish — the b=0 level-6 residue across the T1 tail  [judgment]

Pattern P6 (full-`h_6` tie) fires identically at a b=0 root in 14 T1 cells.
A single lemma "no window tuple makes `14336 d1^2 d2 + 8192 d1 e - 3072 sigma^2`
vanish at a root of `q` while respecting the terminal `ehat^3 g_7 = -8192 u^7 d1^2`
budget" would discharge the level-6 obligation for all 14 at once. Likewise
P11/P12 (full `h_4`) is a common front across the same 17 T1 cells.

### L-Efib (Proposition E generalization) — infinity endgame, shared  [precedent proven for D=5; conjectured D≤9]

`T5_STRATA_50_11.md` §3 Proposition E kills any tuple reduced to
`d1≡0, sigma=s∈K^×, deg e = D` by an infinity-place Newton-polygon domination
(the `f=6` term wins uniquely, `H_6 = -3072 s^2 ≠ 0`). Proven at `D=5`
(Theorems 3–4); the `D≤9` extension is **CONJECTURED** in that file's §4 and
§5. Every T2 cell that L-squeeze reduces to `sigma=const` funnels into this
lemma. **This is a Phase D ingredient** (see §5), flagged here because it is
the shared closer of the whole T2 column.

---

## 5. What Phase C cannot see (Phase D / infinity)

Obligations whose refutation needs the place at infinity or global
leading-coefficient data, and must NOT be filed as local q-residue work:

1. **All `t`-place `term_cancellation` obligations (depths 3–18; 294 occ).**
   STATE.md §"T-PLACE COUPLING" records that "deep t-place cancellations are
   always available to zero-budget states — the t constraint is a residue-level
   constraint, as the by-hand a=9 proofs already showed". They carry the bulk
   of every T2 cell's total depth and cannot be closed by q-local ultrametrics.
   **[data on depths; Phase-D attribution is judgment backed by STATE.]**

2. **Every `exact_identity` (435 occ) and `identical_vanishing` (260 occ).**
   `exact_identity` asserts a genuine polynomial identity
   `ehat^3 g_l = -u^l h_l` up to `deg g_l ≤ 37`; refuting it is a degree/
   leading-coefficient argument, exactly the `T5_90_T1.md` §3 constant-E
   descent (degrees 238→226) and §2 infinity domination. These are depth-0 in
   the engine precisely because they are *not* local valuation obstructions.
   **[judgment]**

3. **The Proposition E closer (L-Efib above).** Infinity Newton polygon; and
   its `D≤9` form is still conjectural.

4. **The whole T1 rigid tail** `a ∈ {8,9,10}` (cells `a8/a9 b=0000/1000/…`,
   `a10 b0000`). `CASCADE_ENGINE_REPORT.md` §"What the 30 survivors need" and
   STATE.md §"CASCADE ENGINE PHASE B" both mark these as "the known rigid
   tail" to be "paired with the `a=9` local-descent machinery and infinity."
   Their cheap q-depth (§3) is real, but the closing argument is the
   `T5_90_T1.md` infinity descent, i.e. Phase D. **[judgment, per those docs.]**

**Net Phase-C-only surface [judgment]:** the depth-1 `q` `monomial_tie_rise`
residue equations (P6/P10/P11/Psq and their cuts, ≈ 480 occ across the T1
front) plus the level-5 squeeze reductions of the T2 column. Everything with
nonzero `t`-depth, every exact_identity, and every degree-domination closure
is Phase D.

---

## Appendix: reproduction

`phase_c_inventory.py` (standalone, read-only) regenerates §1's pattern table,
the tied-vs-`h_l` cross-check, and §3.1's cell ranking:

```
python phase_c_inventory.py          # human summary
python phase_c_inventory.py --json   # machine-readable
```

It imports `build_signature` from `cascade_signature.py` only to fetch the
`h_l` tables for the mismatch check; all other data is read from
`cascade_cones_qt.json`.
