# The Chain Natural-History Survey

A faithful re-implementation of GGV5's complete-chain enumeration algorithm, run
well past the published `v11(A_0) <= 35` tables, mapping the demography of the
problem out to `v11 <= 100`.

- **Enumerator:** `chain_survey.py` (pure integer/`Fraction` Python; no sympy).
- **Regression checker:** `chain_survey_verify.py` (`--quiet`; exit 0 on pass).
- **Full enumeration data:** `chain_survey_data.json` (3995 family rows at M=100).
- **Ground truth:** GGV5 = `paper_src/1708.07936_GGV5.tex`, the two family
  tables at lines 1674-1718, and the eight algorithms that generate them
  (`GetPossibleLastLowerCorners`, `GetStartingEdges`, `GetGeneratedCorners`,
  `GetCornerChildrenList`, `GetCompleteChains`, `GetIsAdmissible`,
  `GetmnFamilies`, `Main algorithm`).
- **Conventions** match `phi_corner4.py` / `composite_charts.py`: corner
  `(a\l,b)` = triple `(a,l,b)`; `t := l_final`, `kappa := t-2` (standard single
  chart), `a0 := A_0.x`, `q := b_final`, `dg := a0-q`, `r := a0-q-1`,
  `e := |m-n|+1`.

---

## 1. Reproduction verdict: EXACT (23/24 verbatim; the 24th is a paper typo)

Running the `Main algorithm` at `M = 35` reproduces the published tables to the
letter. The length distribution of the canonical admissible chains is
`{length 1: 17, length 2: 7}` — precisely the **17 length-1 families (F1-F17)**
and **7 length-2 families (F18-F24)** of the paper.

| Check | Result |
|---|---|
| All 24 published chains `(A_0, A_0', ..., final corner, k)` reproduced | **exact** — no missing, no extra |
| `(m,n)`-family parametrization matches the printed table verbatim | **23 of 24** |
| Diophantine identity `(m+n)bk - n(bl-a) = k` on every enumerated family | holds (M=35 and M=55) |
| Coprimality `gcd(m,n)=1` on every base pair and along every progression | holds |
| `kappa = t-2`, `dg = a0-q` on every family | holds |

**The single discrepancy is F6, and it is a typo in the paper, not in the
algorithm.** GGV5 prints `F6` as `(m,n) = (3j+4, 8j+10)`, i.e. base pair
`(4,10)` — but `gcd(4,10) = 2`, which violates the coprimality that Definition
"mn families" (line 1500) explicitly requires. `(4,10)` therefore cannot be a
member of `MN_k`. The enumerator returns the coprime correction
`(m,n) = (6j+7, 16j+18)`, base `(7,18)`, step `(6,16)`; it satisfies the same
Diophantine identity `(7+18)*4*2 - 18*11 = 2 = k`. F6 is the coprime-restricted
sub-progression that the algorithm's own `GetmnFamilies` produces (compare the
sister family F4, `(2j+3,12j+16)`, which the paper *does* print in its coprime
`k̄`-stepped form). `phi_corner4.py` transcribed the printed typo verbatim
`(4,3),(10,8)`; the algorithm corrects it.

Two structural notes on faithfully reproducing the table:

- **Branch-at-a-final-corner.** `GetChildrenAndFinalList` runs
  `GetCornerChildrenList` on *every* generated corner, including final-shaped
  ones. F22-F24 depend on this: the corner `(14/4,6)` is final-shaped
  (`l - a/b = 1.67 > 1`) yet carries `gcd(14,6)=2` children — it produces no
  `(m,n)`-family of its own (`I(14/4,6)` is empty because `gcd(6,10)=2`) but
  spawns the length-2 chains F22, F23, F24.
- **Canonical dedup.** The literal algorithm also emits redundant chains whose
  first edge is degenerate-vertical (`A_0'` with `b'=0`, direction `(1,0)`,
  refinement `rho=1`): its type-IIb "column" generation reproduces a corner
  `A_1` already reached directly as the type-IIa step `A_0' = A_1`. The
  published table keeps the canonical `A_0' = A_1` form (this is why F18-F21
  list `A_0'=(6,15)` rather than `(6,0)`). We dedup by the generated-corner
  sequence, preferring the `A_0'=A_1` representative; this yields exactly the 24
  rows. Both raw and canonical chain counts are recorded (`n_admissible_raw`,
  `n_canonical_chains`).

Aside: GGV5's prose (line 1671) says "2 admissible complete chains of length 2",
which is inconsistent with its own 7-row length-2 table. The algorithm yields
**7** length-2 family-chains, matching the table, not the prose.

---

## 2. Extension ceiling reached: v11 <= 100

The honest compute ceiling in the ~1.5 h budget was **M = 100** (95-100 s for
the top rung; the full sweep 35->100 runs in ~170 s). Growth of the complete-
chain tree is steep (~M^3.3): `complete chains` go 67 -> 33 567 over 35 -> 100.

| M | complete | adm (raw) | canonical | families | by length |
|---:|---:|---:|---:|---:|---|
| 35 | 67 | 27 | 23 | 24 | 1:17, 2:7 |
| 45 | 359 | 124 | 87 | 95 | 1:58, 2:37 |
| 55 | 914 | 310 | 214 | 251 | 1:123, 2:118, 3:10 |
| 65 | 2 581 | 838 | 528 | 619 | 1:235, 2:310, 3:74 |
| 75 | 6 077 | 1 678 | 977 | 1 092 | 1:342, 2:550, 3:200 |
| 85 | 12 029 | 2 732 | 1 668 | 1 890 | 1:601, 2:961, 3:328 |
| 100 | 33 567 | 5 706 | 3 403 | 3 995 | 1:1130, 2:1972, 3:869, 4:24 |

---

## 3. Demography

### 3a. Do chain types stabilize? — NO. New sporadic corners keep entering.

Every fine-grained census grows monotonically with the bound and shows **no
sign of flattening**:

| M | families | A_0 motifs | final-types (l,q,k) | corner-sigs (a0,t,q,k) | distinct t | max len |
|---:|---:|---:|---:|---:|---:|---:|
| 35 | 24 | 9 | 17 | 19 | 5 | 2 |
| 45 | 95 | 21 | 47 | 58 | 8 | 2 |
| 55 | 251 | 35 | 83 | 112 | 9 | 3 |
| 65 | 619 | 53 | 130 | 210 | 10 | 3 |
| 75 | 1 092 | 73 | 195 | 319 | 12 | 3 |
| 85 | 1 890 | 103 | 302 | 539 | 18 | 3 |
| 100 | 3 995 | 147 | 452 | 982 | 21 | 4 |

New corner types keep being born: the final-corner `l`-value (`= t`) first
reaches **13 at M=55, 16 at M=65, 25 at M=100**; the observed `t`-set at M=100 is
`{3..22, 25}`. The A_0-motif count grows roughly linearly in M. This is the
**anti-finite-grammar signal** at the level of literal corner signatures.

### 3b. Finite local grammars? — Yes at the COARSE level, no at the fine level.

Clustering the corners by a *coarse regime label* — `(sign of the resonance
gap, sign of r, dg parity, chain length)` — tells a different story from the raw
signatures:

| M | 35 | 45 | 55 | 65 | 75 | 85 | 100 |
|---|---:|---:|---:|---:|---:|---:|---:|
| coarse regimes | 8 | 10 | 17 | 20 | 20 | 20 | 26 |

The coarse-regime count **plateaus near ~20** across M = 65-85 and only steps up
(to 26) at M=100 — exactly when the first length-4 chains appear and add new
`(..., length=4)` combinations. Interpretation: **within a fixed maximum chain
length the number of local reduction grammars is bounded; the apparent
unbounded growth is driven by (i) unbounded numeric labels `(t, a0, q)` decorating
a bounded set of grammar shapes, and (ii) the opening of new chain-length
strata.** The finite-grammar hypothesis survives at the level of *reduction
shape*; it fails at the level of *labelled signature*.

This is corroborated by **heavy reuse of a few low-`t` final-corner templates**.
The most recurrent final reductions across distinct chains at M=100 are

```
(7/3,4) k1 : 100 chains      (11/3,8) k1 : 57
(8/3,5) k1 :  84             (13/3,7) k1 : 54
(10/3,7) k1 : 63             (7/4,3) k1 : 38
```

A handful of `t = 3` (and a few `t = 4`) templates are the durable "grammar
atoms", reused by up to a hundred different starting motifs; the ever-growing
tail of high-`t` templates each occur only once or twice.

### 3c. Maximum chain length grows very slowly.

`max length` climbs `2 (M<=45) -> 3 (M=55) -> 4 (M=100)`, i.e. logarithmically-
slowly, capped by `Lmax = bigOmega(gcd(b,(b-b')/rho)) + 1`. Length 2 is the
modal class (1972 of 3995 at M=100); length 4 is barely populated (24 rows, 14
distinct `(A_0, final)`), all descending from `A_0 = (18,72)` and `(24,72)`
through the shared prefixes `(18,72)->(18,54)->(12,30)->...` and
`(24,72)->(24,60)->(32/3,20)->...` — a single new depth stratum that opens right
at the M=100 boundary.

### 3d. Recurring motifs.

The named published motifs persist and recur: `(5,20)` carries 5 families,
`(8,24)` carries 4, `(9,24)` carries 4, `(6,15)` carries 2 at M=100. But they are
dwarfed by the high-`v11` motifs that dominate the tail:

```
(20,80): 323 families    (18,72): 189    (16,80): 142    (15,75): 99
(24,72): 318             (14,84): 156    (10,90): 114    ...
```

The "(5,20)-type" and "(8,24)-type" of the review are the small-`v11`
representatives of families `A_0 = (a0, 4*a0)` and `A_0 = (a0, 3*a0)` — the
`v_{1,-1}(A_0) < 0` cone fills in densely as the bound rises.

### 3e. Distributions (M = 100, 3995 families).

- **`k` (family index):** concentrated on `k=1` (2208), decaying `k=2:906,
  k=3:391, ... k=9:6`.
- **`t = l_final`:** peaks at `t=3 (850), t=5 (773), t=4 (688)`, long thin tail to
  `t=25`.
- **`dg = a0-q` (residual-divisor degree) and `r = dg-1` (spare count):** spread
  smoothly over `dg = 1..33` (`r = 0..32`); `r >= 1` for 3819/3995 families
  (`r = 0`, the "no-spare" corners, number 176). **`dg > 0` and `r >= 0` hold
  universally** — the negative-`r` regime is empty everywhere in `v11 <= 100`.
- **`e = |m-n|+1` at the base pair:** heavy-tailed; `e=2` (adjacent `m,n`) is the
  single largest class (858), tail out past `e = 390`.
- **Window-denominator law.** See the dated correction in **Section 3f** below —
  the proxy statistic `t*a - kappa` used here originally is superseded by the
  exact invariant `q_window = M/gcd(M,H)`, which changes the conclusion about
  `(72,108)`.

### 3f. CORRECTION (2026-07-24): the exact window invariant `q_window` supersedes the `t*a-kappa` proxy

> **This block replaces the proxy "window-denominator law" of Section 3e.** The
> original census scored each family by the reduced denominator of the proxy
> `t*a - kappa` (the F2 `5a-3` slope) and concluded that integrality of the
> window was "`(72,108)`'s private coincidence, not a family trait." That proxy
> is **not** the invariant that governs the bigraded window lattice, and the
> conclusion it supported is **wrong**. The exact object is proved and censused
> in `q_window_theorem.py` (checks appended to `chain_survey_verify.py`, block E).

**The invariant.** For a family with fixed corner data `(t, kappa, q)` and a
moving member `(a,b) = (m,n)(j)`,

```
q_window(a,b) = M / gcd(M, H),   M := t(a+b) - (kappa+1),   H := q(a+b) - 1.
```

**Identity (proved symbolically).** `t*H - q*M = q*(kappa+1) - t =: C`, a **fixed
corner integer** (the `(a+b)` terms cancel identically). Under `kappa=t-2`,
`C = q(t-1) - t`.

**Divisibility lemma.** `gcd(M,H) | C`. Verified to hold on **all 3995** census
family rows. Consequence — the **family-level statement**: along a fixed family
`(a+b)` grows linearly in `j`, so `M` grows linearly while `gcd(M,H)` stays
pinned to a divisor of the fixed `|C|`; hence **`q_window` grows ~linearly in
`M`, with cancellation bounded by the fixed corner integer `C`.** For F2
(`C=3`): `q_window = 7, 12, 17, 22, …` (step `5 = t(dm+dn)/gcd`).

**Known cases (exact).** `(72,108)`: `t=4,kappa=2,q=7,(a,b)=(2,3)` → `M=17,
H=34, gcd=17, q_window=1` (**integral**). `F2 a=2 (50,75)`: `M=21,H=9 → 7`;
`F2 a=3 (75,125)`: `M=36,H=15 → 12`; `F9 (56,84)`: `M=29,H=9 → 29`.

**The census verdict (what the proxy could not establish).** Evaluating
`q_window` at each distinct family's derived base member `(a,b)=(m0,n0)`:

| quantity | value |
|---|---:|
| distinct families (M=100) | 1848 |
| **integral (`q_window=1`) families** | **51** |
| distinct integral corner-shapes `(t,q,kappa,{a,b})` | **23** |
| `(72,108)` among them | yes |
| **`(72,108)` the UNIQUE integral case** | **NO** |

**`(72,108)` is not unique — it is one member of an arithmetic lattice of
integral windows.** The integral shapes fall into clean progressions, e.g. at
`t=3`: `q = 8,11,14,17,20,23` with base `(2,3),(3,4),(4,5),(5,6),(6,7),(7,8)`;
at `t=4`: `q = 7,11,15,18,19` — `(72,108)` sits at `(t=4,q=7,base{2,3})`. The
proxy census (`t*a-kappa`) genuinely **could not** resolve this — it reduced
against the wrong datum and read integrality as a coincidence. The exact
`q_window` settles it: integrality is a **structured, non-unique** phenomenon
governed by `M | H ⇔ M | C`.

---

## 4. Anomalies flagged

1. **F6 is a published typo** (Section 1): printed base pair `(4,10)` has
   `gcd = 2`, contradicting GGV5's own coprimality requirement. Correct coprime
   family: `(6j+7, 16j+18)`.

2. **The negative-resonance-gap regime is populated.** `composite_charts.py`
   (STEP 6) flagged "NEGATIVE gap: res < pure, regime unobserved anywhere". The
   extended survey **observes it in 588 of 3995 families** at M=100. It is the
   large-`|m-n|` tail: once `e = |m-n|+1` is large enough, the resonance degree
   `res = (t(b-a)+kappa+1)*a0/t` falls below the pure-ansatz degree
   `e*a0 - q + 1`. The regime is not exotic — it is simply invisible in the tiny
   `v11 <= 35` escape survey and unavoidable at scale. Any Phi theory that
   assumed `gap >= 0` needs a `gap < 0` branch.

3. **A new length stratum opens exactly at the ceiling.** The first length-4
   admissible chains appear at M=100 (none at M<=85), all from `A_0 in
   {(18,72),(24,72)}`. Whether length 5 first appears near M~130-160 is the
   natural next probe; the coarse-regime count is predicted to step up again when
   it does.

4. **`GetIsAdmissible` does not kill F18-F21.** Consistent with GGV5's own
   narrative, F18-F21 are admissible complete chains produced by the algorithm;
   they are excluded only by the *separate* last-lower-corner argument
   (GGV5 lines 1726-1786), reproduced independently in `composite_charts.py`.
   They appear in this survey's enumeration and are *not* filtered by the
   divisibility criteria.

5. **Prose/table mismatch in GGV5** ("2 admissible complete chains of length 2"
   vs. the 7-row length-2 table) — the algorithm sides with the table.

No family violated the Diophantine identity, coprimality, `kappa = t-2`, or
`dg = a0-q`. No corner outside the `r >= 0`, `dg >= 1` cone was found.

---

## 5. Bottom line

The problem's demography is **open, not closed**: past the published tables the
enumeration keeps minting new starting motifs, new final-corner reductions, and
new (unboundedly large) `t`-values, with the maximum chain length creeping up one
notch every ~40-50 units of `v11`. The finite-grammar picture holds only after
you quotient out the numeric labels: a small, apparently bounded set of coarse
reduction shapes (dominated by a few `t = 3,4,5` atoms) is reused across an
unbounded, linearly-growing population of labelled corners. The single hard
correction to the record is the F6 coprimality typo; the single conceptual
correction is that the negative-resonance-gap regime, previously "unobserved",
is in fact routine at scale.

---

## Files

| File | Role |
|---|---|
| `chain_survey.py` | the GGV5 enumerator + statistics + JSON export (`--sweep`, `--out`) |
| `q_window_theorem.py` | the exact `q_window = M/gcd(M,H)` invariant: identity + divisibility lemma (symbolic) + known-case table + census (Section 3f) |
| `chain_survey_verify.py` | regression: published table reproduced + invariants + `q_window` theorem checks block E (`--quiet`, exit 0) |
| `chain_survey_data.json` | full enumeration: per-M census + all 3995 family rows at M=100 |
| `CHAIN_SURVEY.md` | this report |
