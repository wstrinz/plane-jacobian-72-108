# Phase C residue-system worklist — f31 subcase-(1) survivors

**Date:** 2026-07-22
**Scope:** the 279-branch / 2253-case frontier of `cascade_cones_sub1_qt.json`
(subcase (1), q+t coupled, depth 4), the Phase C input described in
STATE.md §"SUB1 T-PLACE COUPLING". This document ports `PHASE_C_WORKLIST.md`
(the subcase-(2) precedent) to subcase (1) and measures pattern overlap
between the two windows.

**How to read this file.** Statements tagged **[data]** are mechanical
consequences of the witness records in `cascade_cones_sub1_qt.json`
(re-derivable by `phase_c_inventory_sub1.py`, which also re-checks every
tied-monomial string against the `h_l` tables). Statements tagged
**[judgment]** are my assessment of Phase C difficulty and lemma reuse, not
forced by the data. Every algebraic monomial is quoted from the
`cascade_signature.py` `h_l` tables (§0); obligation *semantics* are quoted
from `CASCADE_ENGINE_REPORT.md` §Semantics (reproduced in the sub2 doc §0).

Provenance note **[data]**: `phase_c_inventory_sub1.py` confirms all **26965**
tied monomials across the 2253 cases are exact members of the `h_l` tables
(0 mismatches) — the tied sets are always sub-multisets of a single `h_l` cut
down by the case's zero flags. The `h_l` tables are **window-independent**
(they are the same coefficients used in `PHASE_C_WORKLIST.md` §0), so every
residue equation below is literally the same polynomial as in subcase (2).

---

## 0. The `h_l` monomial tables used below  [data]

From `cascade_signature.py` (rewrite `sigma = 4*d0 - d2**2`, variables
`(d2,d1,sigma,e)`, `e = d̃₋₁`), verbatim, identical to the sub2 doc §0:

```
h_7 = 8192*d1^2
h_6 = 14336*d1^2*d2 + 8192*d1*e - 3072*sigma^2
h_5 = -12288*d1^2*d2^2 + 32256*d1^2*sigma + 18432*d1*d2*e - 9216*d2*sigma^2 + 2048*e^2
h_4 = -220752*d1^4 - 31232*d1^2*d2^3 - 23616*d1^2*d2*sigma - 3072*d1*d2^2*e
      + 34560*d1*e*sigma - 5184*d2^2*sigma^2 + 5632*d2*e^2 - 12096*sigma^3
```

Cascade identity, local reading at a q-root, and the four obligation kinds
(`monomial_tie_rise`, `term_cancellation`, `exact_identity`,
`identical_vanishing`) are exactly as in `PHASE_C_WORKLIST.md` §0 — not
repeated here. The only window-specific inputs are the **subcase-(1) caps**
(`sub1_cascade_verify.py`): `deg d2 ≤ 6`, `deg d1 ≤ 9`, `deg d0 ≤ 12`,
`deg sigma ≤ 12`, `deg e ≤ 15`; rigorous terminal caps `deg g_7 ≤ 46` (T1),
`deg g_6 ≤ 48` (T2); `deg h_f(d̃) ≤ 60 − 6f`. (These are looser than sub2's
`deg g_l ≤ 10+3a`, which is why sub1's frontier is wider — STATE.md §"SUB1
DEPTH-5 SWEEP".)

---

## 1. Pattern inventory  [data]

Parsed by `phase_c_inventory_sub1.py` over the 2253 survivor cases (279
branches). Grouping key = `(place kind, cascade level, obligation kind,
tied-monomial set)`.

- **67 distinct obligation patterns.** (sub2: 41.)
- Raw obligation counts by kind (total 21337): `term_cancellation` 9978,
  `monomial_tie_rise` 5209, `exact_identity` 4035, `identical_vanishing` 2115.
- By place: `q` 15580, `t` 5757. By level: L4 8382, L5 7569, L6 5386.

### 1.1 The 15 highest-frequency patterns

| # | place | L | kind | tied N | freq | cells | depth multiset |
|--:|:--|--:|:--|--:|--:|--:|:--|
| 1 | q | 6 | term_cancellation   | 0 | 2799 | 159 | 1×2799 |
| 2 | q | 5 | term_cancellation   | 0 | 2388 | 247 | 1×2388 |
| 3 | q | 4 | exact_identity      | 0 | 2108 | 269 | 0×2108 |
| 4 | q | 4 | term_cancellation   | 0 | 1720 | 247 | 1×1720 |
| 5 | t | 6 | term_cancellation   | 0 | 1263 | 168 | 3–12 |
| 6 | q | 5 | exact_identity      | 0 | 1120 |  89 | 0×1120 |
| 7 | t | 5 | term_cancellation   | 0 | 1046 | 258 | 3–14 |
| 8 | t | 4 | term_cancellation   | 0 |  762 | 258 | 3–18 |
| 9 | q | 6 | monomial_tie_rise   | 3 |  621 |  72 | 1×621 |
| 10| q | 5 | monomial_tie_rise   | 5 |  600 | 159 | 1×600 |
| 11| t | 4 | exact_identity      | 0 |  527 | 269 | 0×527 |
| 12| q | 4 | monomial_tie_rise   | 8 |  504 | 159 | 1×504 |
| 13| q | 4 | identical_vanishing | 8 |  402 | 126 | 0×402 |
| 14| q | 4 | monomial_tie_rise   | 3 |  383 | 159 | 1×383 |
| 15| q | 5 | monomial_tie_rise   | 2 |  379 | 131 | 1×263,2×90,5×8,10×18 |

Structural facts **[data]**, mirroring the sub2 doc:

- The `term_cancellation` patterns (#1,2,4,5,7,8) carry no tied monomials:
  they constrain the *free leading residue of `g_l`*, not the window data.
  **At `q` they are all depth 1** exactly as in sub2; at `t` they run 3…18.
- `exact_identity` patterns (#3,6,11) are the `g_{l+1}≡0` book-keeping.
- The window-only residue obligations — genuine equations in `(d2,d1,sigma,e)`
  — are the `monomial_tie_rise` and `identical_vanishing` families (7324 of
  21337). Their tied sets are always the surviving cut of a single `h_l`.
- Patterns #9/#10/#12 are the **full `h_6` / `h_5` / `h_4` ties**, i.e. the
  sub2 backbone patterns P6/P10/P11 recurring verbatim (§2). Pattern #13 is
  the full `h_4` `identical_vanishing`.

The q-depth histogram is again essentially two-valued: **every q obligation
has depth 0 or 1** except a thin tail (pattern #15's depths 2/5/10 and a
handful of level-4 depth-5/11 outliers). All the heavy depth lives at the `t`
place, whose depths reach 30 (see §3, §5). This is the same q/t split as sub2.

---

## 2. THE KEY METRIC — overlap with subcase (2)  [data]

Recomputed by `phase_c_inventory_sub1.py`, keying the sub2 survivors
(`cascade_cones_qt.json`) identically and comparing against the 67 sub1
patterns. **Definitive answer:**

| relation | count |
|:--|--:|
| sub2 patterns **absent** from sub1 (`sub2_only`) | **0** |
| sub1 patterns **identical** to a sub2 pattern (same key **and** depth-set) | **29** |
| sub1 patterns with the **same key, deeper/extra depths** | **12** |
| sub1 patterns **genuinely new** (key absent from sub2) | **26** |
| — of the 26 new, **reuse an existing sub2 tied-set** (only place/kind differs) | 20 |
| — of the 26 new, **genuinely new residue polynomial** | **6** |

Read this prominently:

> **Every one of the 41 subcase-(2) patterns recurs in subcase (1)**
> (`sub2 ⊆ sub1` at the pattern-key level). The sub1 frontier adds 26 new
> keys, but **20 of them reuse a residue polynomial already in the sub2
> library** and only **6 carry a residue polynomial not seen in sub2**.

**[judgment]** One residue-lemma library therefore serves both windows almost
entirely. Concretely, of the 67 sub1 pattern keys, 41 are shared with sub2 and
20 more re-use a shared tied-set — **61 of 67 (91%) touch no algebra beyond
the sub2 library**; only 6 need new lemmas, and (below) all 6 are sub-cuts of
the *same* `h_4`.

### 2.1 The 12 "same key, deeper depth" patterns  [data]

These share the *exact* residue polynomial with sub2; only the required order
of vanishing is larger, because sub1 spans the full `a = 0…10` (sub2's
survivors were the `a = 5…10` rigid tail). Examples:

- `t L6 term_cancellation` (tied=∅): sub2 depths `{3,6,9,12,15}`, sub1 depths
  `{2,3,4,6,7,9,10,12,13,14,15,17,18,20,21,22,24,26,27,30}`.
- `t L4 monomial_tie_rise` (`{-9216 d2 sigma^2, 2048 e^2}` cut): sub2
  `{3,6,9,12}`, sub1 `{3,6,…,27}`.
- `q L5 monomial_tie_rise` (tiedN=2): sub2 `{1,2,5}`, sub1 `{1,2,5,10}`.

The extra depths are all multiples/near-multiples of the term running down
from `v = 30 − 3a` (§3.3): the same lemma, evaluated at more `a`.

### 2.2 The 6 genuinely-new residue polynomials  [data]

All six are **cuts of `h_4`** produced by richer zero-flag combinations than
sub2's survivors reached (e.g. `d1≡0`+`sigma` partial, or `d2` partial):

| place | L | kind | tied set | freq |
|:--|--:|:--|:--|--:|
| q | 4 | identical_vanishing | `{-12096 σ³, 34560 d1 σ e}` | 16 |
| t | 4 | monomial_tie_rise   | `{-12096 σ³, 34560 d1 σ e}` | 14 |
| t | 4 | identical_vanishing | `{-3072 d2² d1 e, -31232 d2³ d1², 5632 d2 e²}` | 13 |
| q | 4 | identical_vanishing | `{-12096 σ³, 5632 d2 e²}` | 11 |
| t | 4 | monomial_tie_rise   | `{-3072 d2² d1 e, -31232 d2³ d1², 5632 d2 e²}` | 4 |
| t | 4 | identical_vanishing | `{-220752 d1⁴, 34560 d1 σ e}` | 2 |

**[judgment]** These are low-frequency (≤16 occ) two- and three-term subsets
of `h_4`; a single "no `h_4`-cut vanishes at a q-root of the fixed `q`" lemma
(the L-h4 family of §6) subsumes all of them together with P11/P12. **No new
lemma *family* is needed** — only new instantiations of the level-4 residue
lemma at additional cuts.

---

## 3. The 26-family anatomy

### 3.1 The family is a-parametric with an a-invariant cell set  [data]

`phase_c_inventory_sub1.py` confirms the survivor `(b-vector, branch)` set is
**identical for every `a = 0…8`** — exactly 26 pairs, listed below; and
`a = 9` drops to 24, `a = 10` to 21. STATE.md §"SUB1 DEPTH-4 SWEEP" records
the same count.

The 26 `(b, branch)` pairs (a ≤ 8):

```
(0000)T1 (0000)T2 (1000)T1 (1000)T2 (1100)T1 (1100)T2 (1110)T1 (1110)T2
(1111)T1 (1111)T2 (2000)T1 (2100)T1 (2110)T1 (3000)T1 (3000)T2 (3100)T1
(3100)T2 (3110)T1 (3110)T2 (3111)T1 (3300)T1 (3310)T1 (5000)T1 (5000)T2
(5100)T1 (5110)T1
```

**Split:** 17 T1 + 9 T2 at each `a ≤ 8`. The `b`-values are drawn from
`{0,1,2,3,5}` — the sub1 live-local set from `CASCADE_CONE_LEMMAS_SUB1.md`
(Family L: live `beta ∈ {0,1,2,3,5}` for T1, `{0,1,3,5}` for T2), never `4`.

**What dies as `a` grows [data]:** the family thins only through the
single-place Family-L kills, whose unconditional dead set is
`beta ∈ {4} ∪ [6, 15−a]` (`CASCADE_CONE_LEMMAS_SUB1.md`). As `a → 9,10` the
cap `15−a` descends, killing the high-`b` cells:

- `a = 9` loses `(3310)T1, (5110)T1` (2 cells) → 24;
- `a = 10` is missing 5 cells relative to the base family —
  `(3111)T1, (3300)T1, (3310)T1, (5100)T1, (5110)T1` — i.e. the two lost at
  `a=9` plus `(3111)T1, (3300)T1, (5100)T1` → 21.

So the *cell set* is a-invariant for `a ≤ 8` and then loses its highest-`b`
T1 members monotonically.

### 3.2 What depends on `a` inside a fixed cell  [data]

For a fixed `(b, branch)`, comparing across `a`:

- **The q-place obligation structure is completely a-invariant.** Every q
  witness carries the same patterns *and the same depths* `{0, 1}` at every
  `a`. (Verified for `(0000)T1`: q-depths `{0,1}` for all `a = 0…10`.) The
  q-residue equations (P6/P10/P11/Psq and the `h_l` cuts) do **not** depend
  on `a`.
- **The flag-case multiplicity varies with `a`** (13,12,19,19,21,25,22,22,23
  for `a = 0…8` in `(0000)T1`): the number of surviving zero-flag cases grows
  then plateaus as the valuation budget `v = 30 − 3a` shrinks.
- **The t-place depths are affine in `a`.** For the b=0, all-`g`-nonzero t
  witness the level-4/5/6 `term_cancellation` depth equals `v` exactly:

  | a | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
  |:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
  | v = 30−3a | 30 | 27 | 24 | 21 | 18 | 15 | 12 | 9 | 6 | 3 |
  | t term_cancel depth (L4=L5=L6) | 30 | 27 | 24 | 21 | 18 | 15 | 12 | 9 | 6 | 3 |

  i.e. **depth = v = 30 − 3a**, affine with slope −3. Deeper `monomial_tie_rise`
  t-depths are `v` minus the tied monomial's e-slot offset, hence also affine
  in `a`.

**[judgment]** The family is genuinely a-parametric, but the a-dependence is
confined to (i) flag-case count and (ii) the t-place depths, both governed by
`v = 30 − 3a`. The **residue polynomials and the entire q-place obligation
system are a-independent.** Consequently a Phase-C lemma proved at one `a`
(for the q-residue equations) proves it for all `a`; only the t-place /
infinity depth budget must be re-checked per `a`. This is exactly the
"single a-independent family, ideal for cone-lemma compression" of STATE.md.

---

## 4. Explicit residue equations for the top 3 residue patterns  [data]

The three highest-frequency *window* patterns (#9/#10/#12) are the full
`h_6`/`h_5`/`h_4` ties — identical polynomials to sub2's P6/P10/P11. At a
`b=0` q-root `p` (all slots units) the level-`l` line needs `v_p(h_l) ≥ 1`:

### P6 — full `h_6` tie, q-place, depth 1 (621 occ, 72 cells)

```
14336 (d1*)^2 (d2*) + 8192 (d1*)(e*) - 3072 (sigma*)^2 = 0.            (P6)
```

### P10 — full `h_5` tie, q-place, depth 1 (600 occ, 159 cells)

```
-12288 (d1*)^2(d2*)^2 + 32256 (d1*)^2(sigma*) + 18432 (d1*)(d2*)(e*)
   - 9216 (d2*)(sigma*)^2 + 2048 (e*)^2 = 0.                           (P10)
```

### P11 — full `h_4` tie, q-place, depth 1 (504 occ, 159 cells)

```
-220752 (d1*)^4 - 31232 (d1*)^2(d2*)^3 - 23616 (d1*)^2(d2*)(sigma*)
   - 3072 (d1*)(d2*)^2(e*) + 34560 (d1*)(e*)(sigma*)
   - 5184 (d2*)^2(sigma*)^2 + 5632 (d2*)(e*)^2 - 12096 (sigma*)^3 = 0. (P11)
```

Pattern #13 (full `h_4` `identical_vanishing`, 402 occ) is the same LHS as
(P11) but demands *all eight* leading residues cancel simultaneously (fires
when `g_5 ≡ 0`) — strictly stronger, as in sub2. The T2 backbone
`h_5|_{d1≡0}` cut (Psq) `-9216 (d2*)(sigma*)^2 + 2048 (e*)^2 = 0` and its
`h_4|_{d1≡0}` companion are the same as sub2 §2's (Psq); they carry all 99 T2
cells.

Because `q` is the *same fixed quartic* in both windows, each of P6/P10/P11
is again one condition on `(d2,d1,sigma,e) mod q` — a resultant-level coupling,
not four independent points **[judgment]**.

---

## 5. Depth histogram and cheapest targets

### 5.1 Per-cell difficulty  [data]

Ranking cells by total obligation depth of their *cheapest* survivor case
(`phase_c_inventory_sub1.py`). The very cheapest are degenerate:

| tier | cells | TOT | q_dep | t_dep | d0 | reading |
|:--|:--|--:|--:|--:|--:|:--|
| A | `a10 (1111)T1/T2` | 0 | 0 | 0 | 0 | **no obligations at all** — pure Phase D, no local handle |
| B | `a8/a9 (1111)`, `a9 (3111)T1`, `a6 (1111)T2` | 0 | 0 | 0 | 5–10 | 0 depth but 5–10 `exact_identity`/vanish constraints |
| C | `a8/a9/a10 (1110)/(3110)` | 1 | 1 | 0 | 5–10 | one depth-1 q-residue equation |
| D | `a10 (0000)T1`, `a9 (0000)T1` | 4 | 4 | 0 | 10 | four full-`h_l` q-residues, no `t`-depth |
| E | `a9/a10 (1000)T2` | 3 | 3 | 0 | 5 | T2 squeeze cell (§6) |

### 5.2 The metric is a proxy  [judgment]

As in sub2 §3.2, total depth **undercounts** the strongest obligations:
`exact_identity`/`identical_vanishing` both score 0. Tier-A cells
(`a10 (1111)`, `d0 = 0`) are *not* cheap to close — they carry **no local
obstruction whatsoever** and are pure infinity/Phase-D cells. Tier-B/C cells
carry many depth-0 constraints (the `d0` column). Weight `d0` and "lands on a
proven lemma", not raw total, when picking real effort.

### 5.3 Recommended three cheapest targets  [judgment]

Choosing cells that are cheap **and** land on existing machinery — the direct
analogues of the sub2 doc's three picks, which is itself evidence of reuse:

1. **`a9 (1000)T2` / `a10 (1000)T2`** (total 3, `t_dep = 0`, `g_5≡0`). Lands
   squarely on the **level-5 squeeze** (§6, `T5_STRATA_50_11.md`): the (Psq)
   residue `-9216 d2 sigma^2 + 2048 e^2` at the three `b=0` roots plus the
   consolidated `q | (-9216 d2 sigma^2 + 2048 e^2)`. Its `t`-obligations are
   depth 0. Identical in shape to sub2's cheapest cell.
2. **`a10 (0000)T1`** (total 4, `t_dep = 0`, `d0 = 10`). The degenerate top of
   the rigid tail: pure full-`h_l` residue equations (P6/P10/P11) at the four
   q-roots, no `b>0` complication, all `t`-obligations depth 0. Lands on the
   L-h6vanish / L-h4 lemma family (§6).
3. **`a9 (1000)T1`** (total 6, `q_dep = 3`, `t_dep = 3`). Pairs with the
   `a = 9` local-descent program of `T5_90_T1.md`; the T1 branch of the same
   stratum. Its q-part is three depth-1 residues; the depth-3 `t`-part is the
   `v = 30−27 = 3` term_cancellation of §3.2.

**[judgment]** All three land on lemmas that already exist for sub2 — the
level-5 squeeze, the full-`h_l` residue front, and the `a=9` machinery. The
T2-infinity narrowings of `T5_T2_INFINITY.md` apply to the sub1 T2 column too
(same `h_5|_{d1≡0}` collapse, same `F^5 W = 3072 c^6 Z^2` terminal), so the
`a9 (1000)T2` pick inherits that partial closure as well.

---

## 6. Phase C vs Phase D filing

Which obligations need only **local q-residue algebra** (Phase C) vs the place
at infinity / global leading coefficients (Phase D):

**Phase-C surface [judgment]:** the depth-1 `q` `monomial_tie_rise` residue
equations — P6/P10/P11/Psq and the six new `h_4` cuts (§2.2) — plus the
level-5 squeeze reductions of the T2 column. These are a-independent (§3.2),
so proving each once settles it for every `a`. Reusable lemmas:

- **L-squeeze** (level-5 squeeze, `T5_STRATA_50_11.md` §2, proven): clears the
  whole T2 column. **Match [data]:** all **99** T2 survivor cells (9 per `a`,
  `a = 0…10`) carry the (Psq) obligation and its `h_4|_{d1≡0}` companion; the
  `d1≡0, sigma≢0` hypotheses hold by construction. One lemma, retuned only by
  the sub1 terminal cap `deg g_6 ≤ 48`, reduces every T2 cell. Highest-leverage
  Phase-C lemma, exactly as in sub2.
- **L-h6vanish / L-h4** (`PHASE_C_WORKLIST.md` §4): a single "no window tuple
  makes `14336 d1^2 d2 + 8192 d1 e − 3072 sigma^2` (resp. any `h_4` cut) vanish
  at a root of the fixed `q` within the terminal budget" discharges P6 across
  the 72 T1 cells that carry it, and P11/P12 plus the six new `h_4` cuts of
  §2.2 across the level-4 front.

**Phase-D (must NOT be filed as local q-residue work) [judgment, per STATE.md
§"T-PLACE COUPLING" and `T5_T2_INFINITY.md`]:**

1. **All `t`-place `term_cancellation` depths (3…30; the bulk of the 5757
   t-obligations).** STATE.md: "deep t-place cancellations are always available
   to zero-budget states — the t constraint is a residue-level constraint".
   They carry each cell's total-depth inflation and are not closable by q-local
   ultrametrics. Their depths are affine in `a` (§3.2) and reach 30 at `a=0`.
2. **Every `exact_identity` (4035 occ) and `identical_vanishing` (2115 occ).**
   A genuine polynomial identity `ehat^3 g_l = −u^l h_l` up to `deg g_l ≤ 46/48`
   (T1/T2); refuting it is a degree/leading-coefficient argument
   (`T5_90_T1.md` §3 constant-`E` descent, §2 infinity domination). Depth-0 in
   the engine precisely because they are not local valuation obstructions.
3. **The Proposition-E closer** (`T5_STRATA_50_11.md` §3; `D≤9` still
   conjectural) and **the T2 infinity convolution** (`T5_T2_INFINITY.md`,
   eqs (6)–(12)): the shared closers of the whole T2/rigid-tail column.
4. **Tier-A cells** (`a10 (1111)`, `d0 = 0`): zero local obligations —
   entirely Phase D.

**Net [judgment]:** the sub1 Phase-C-only surface is the *same* residue
system as sub2 — P6/P10/P11/Psq, now with the six extra `h_4` cuts and the
99-cell (vs 12-cell) T2 column — and it is a-independent, so it is proved once
per pattern and quantified over `a`. Everything with nonzero `t`-depth, every
`exact_identity`/`identical_vanishing`, and every degree-domination closure is
Phase D.

---

## Appendix: reproduction

`phase_c_inventory_sub1.py` (standalone, read-only) regenerates §1's pattern
table, the tied-vs-`h_l` cross-check, §2's overlap metric, §3.1's family
a-invariance, and §5.1's cell ranking:

```
python phase_c_inventory_sub1.py          # human summary (incl. overlap block)
python phase_c_inventory_sub1.py --json    # machine-readable
```

It reads `cascade_cones_sub1_qt.json` (sub1) and `cascade_cones_qt.json`
(sub2, for the overlap metric), and imports `build_signature` from
`cascade_signature.py` only to fetch the window-independent `h_l` tables for
the mismatch check.
