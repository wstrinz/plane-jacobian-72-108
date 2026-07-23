# Depth-1 residue-congruence layer for the alternate regime

Date: 2026-07-22
Engine: `alt_residue_congruences.py` -> `alt_residue_congruences.json`
Verifier: `alt_residue_congruences_verify.py` (**PASSES**, 4 hand checks)
Inputs (all already audited, NONE edited here):
- surviving states: `alt_combined.json` (3102 states over 27 open branches; `ALT_COMBINED.md`)
- flipped chain: `alt_inf_sweep.py` `run_chain` (the `(D_t)` max-plus cascade; `ALT_REGIME_INF.md`)
- monomial tables + kill library: `cascade_engine.py` (`MONOMIALS`, `tropical_h_max_full`,
  `deg_h_options`, `FORBIDDEN_RISES`, `LC_U`), `RESIDUE_LEMMAS.md` (23 supports; C08/C20 the two kills)
- constant provenance: `cascade_inf_ties_verify.py` (`lc(u) = -1024/3315`, leading coeffs at
  infinity are base-field elements)

This is the step `ALT_COMBINED.md` [judgment] **J1** left open: closing a branch
"needs the residue congruences (the leading-cancellation coefficient conditions),
which are a strictly finer fact than any valuation lower bound and are out of
scope here." This document builds the **depth-1** slice of exactly those
congruences and classifies every state's system.

## 1. What the depth-1 layer is (every step cited)

Each surviving degree state fixes `(deg d2, deg d1, deg sigma, deg e)` and hence
`deg_E = deg e - a`. Running the flipped descending cascade
`(D_t)` (`ALT_REGIME_INF.md` (b), `T r_{f-1} = E^(3(7-f)) h_f + u r_f`,
`r_7 = 0`, top anchor `T r_6 = h_7`, bottom close `E^21 h_0 + u r_0 = 0`) with
the audited max-plus contract of `cascade_engine.deg_h_options`, each state emits
a minimal-witness list of **obligations**. An obligation is either a
`degree_tie_drop` at a level `f` (the leading form of `h_f` must cancel to some
depth `delta`) or a `leading_cancellation` (a cross-term `E^(3(7-f)) h_f` vs
`u r_f` cancellation). Writing `(D,X,S,E)` for the leading coefficients of
`(d2,d1,sigma,e)` at the place at infinity (same convention as
`RESIDUE_LEMMAS.md` §1; base-field elements by `cascade_inf_ties_verify.py` §C),
the **exact depth-one leading-coefficient equation** of each obligation is:

- `degree_tie_drop` at level `f`, tied support `T`, depth `delta`
  (`RESIDUE_LEMMAS.md` (IF)):

  ```
  sum_{coef * d2^k d1^x sigma^z e^b in T}  coef * D^k X^x S^z E^b  =  0.
  ```
  (This is the depth-1 equation; `delta >= 2` stacks `delta-1` deeper
  convolution equations on top of it, all containing this one.)

- bottom close `leading_cancellation` `E^21 h_0 + u r_0 = 0`
  (`ALT_REGIME_INF.md` (I0), `deg u = 4`, `lc(u) = -1024/3315`):

  ```
  lc(E)^21 * lc(h_0)  +  (-1024/3315) * lc(r_0)  =  0.
  ```

### The one structural fact that decides the whole layer

**For every one of the 3102 survivors the minimal obligation list is exactly two
records, BOTH at level 0** (the bottom close): one `degree_tie_drop` on the
`h_0` initial form, and the close `leading_cancellation`. There are **no
level-1..6 obligations**. The reason is the generic dominance already visible in
`ALT_REGIME_INF.md`'s worked table: at every intermediate level `f = 6..1` the
`g`-side `3(7-f) deg_E + H_f` strictly dominates the `h`-side `4 + R_f`
throughout the surviving window, so `h_f` is taken at its maximum (forced, no
drop, no obligation). Only the closing anchor `f = 0` forces a tie, producing the
two level-0 records.

## 2. Classification (soundness-first)

The **only** depth-1 kills proven anywhere in the library are C08 (level 5,
`6X^2D^2 - 9XDE - E^2`, square class 105) and C20 (level 4,
`61X^2D^2 + 6XDE - 11E^2`, square class 170); both are arithmetic — no all-nonzero
solution over `Q` or over the `q`-splitting field `Q(sqrt(17))`
(`RESIDUE_LEMMAS.md` §4, `residue_lemmas_verify.py`). At infinity the unknowns are
leading coefficients living in the base field, so those obstructions transport
(`cascade_inf_ties_verify.py` §C).

A state **dies at depth 1 iff every viable flipped chain is FORCED to drop `h_f`
on exactly a C08/C20 support** — equivalently, iff switching
`APPLY_RESIDUE_KILLS = True` (which makes `deg_h_options` refuse a drop on those
two supports) removes its last surviving chain. We decide this the sound way:
re-run `run_chain` with kills on and compare.

| class | rule | count |
|:--|:--|--:|
| **(a) KILL** | dies under kills-on (a required C08/C20 drop) | **0** |
| **(b) singleton obligatory tie** | `len(T)=1`, unique monomial cannot drop — a sweep bug | **0** (none; sweep is right) |
| **(c) CONSTRAINT** | keeps a chain under kills-on; each support a non-empty hypersurface | **3102** |

### The soundness subtlety (why 0 die even though C08/C20 appear)

The C08 (level 5) and C20 (level 4) forbidden supports **do** occur — as
tropical ties at their levels, in **523** and **180** state-occurrences
respectively (all in `sigma`-identically-zero T1 states, matching the
`RESIDUE_LEMMAS.md` §3 `sigma=0` column where C08/C20 live). **But the chain
takes `h_5`, `h_4` at their maximum there; it never needs those leading forms to
cancel.** A kill lemma forbids a *drop*; a tie sitting at the max with no drop
required is not an obstruction. Concretely `deg_h_options` at those levels returns
the max option `(maximum, ())` with or without kills, so kills-on `==` kills-off
for all 3102 states (verified: `kills_on_equals_kills_off = true`). Claiming a
kill here would violate the soundness direction the task demands — a
non-obligatory forbidden tie is **not** a kill.

The bottom-close cancellation is likewise never a kill: it is linear in the free
nonzero unknown `lc(r_0)` (`ALT_REGIME_INF.md` [judgment]: no cap on the
descending cofactors `r_f`), so `lc(r_0) = -lc(E)^21 lc(h_0)/lc(u)` is always a
nonzero solution.

## 3. Census

**States: 3102 killed=0 / constrained=3102. Whole-branch kills: 0.** The depth-1
residue-congruence layer closes **no** branch — an honest, sound negative result
that sharpens `ALT_COMBINED.md` J1.

Per branch (all 27 OPEN; `forbid_tie` = states carrying a non-obligatory C08/C20
tropical tie):

| id | states | constraint | kill | forbid_tie | id | states | constraint | kill | forbid_tie |
|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|
| a11_b0000_T1 | 408 | 408 | 0 | 99 | a12_b0000_T1 | 447 | 447 | 0 | 98 |
| a11_b1000_T1 | 301 | 301 | 0 | 71 | a12_b1000_T1 | 325 | 325 | 0 | 61 |
| a11_b1100_T1 | 202 | 202 | 0 | 43 | a12_b1100_T1 | 235 | 235 | 0 | 33 |
| a11_b1110_T1 | 133 | 133 | 0 | 24 | a12_b1110_T1 | 164 | 164 | 0 | 12 |
| a11_b1111_T1 | 80 | 80 | 0 | 9 | a12_b3000_T1 | 144 | 144 | 0 | 12 |
| a11_b3000_T1 | 106 | 106 | 0 | 16 | a14_b0000_T1 | 213 | 213 | 0 | 33 |
| a11_b0000_T2 | 27 | 27 | 0 | 0 | a14_b1000_T1 | 144 | 144 | 0 | 12 |
| a11_b1000_T2 | 22 | 22 | 0 | 0 | a12_b0000_T2 | 21 | 21 | 0 | 0 |
| a11_b1100_T2 | 17 | 17 | 0 | 0 | a12_b1000_T2 | 16 | 16 | 0 | 0 |
| a11_b1110_T2 | 13 | 13 | 0 | 0 | a12_b1100_T2 | 12 | 12 | 0 | 0 |
| a11_b1111_T2 | 9 | 9 | 0 | 0 | a12_b1110_T2 | 8 | 8 | 0 | 0 |
| a11_b3000_T2 | 12 | 12 | 0 | 0 | a13_b0000_T2 | 15 | 15 | 0 | 0 |
| a11_b3100_T2 | 8 | 8 | 0 | 0 | a13_b1000_T2 | 11 | 11 | 0 | 0 |
| | | | | | a14_b0000_T2 | 9 | 9 | 0 | 0 |
| **total** | | | | | | **3102** | **3102** | **0** | **703 occ / 523 states** |

The `h_0` initial forms collapse to **19 distinct hypersurfaces** in `(D,X,S,E)`.
Every one is multi-term (no singleton), matches **none** of C01-C23 (they live at
level 0, not 4/5/6), and carries an exhibited **all-nonzero rational point** — so
each is a genuine CONSTRAINT, not a kill. Selected supports (full list +
factorisations + witnesses in the JSON `support_catalog`):

| id | n | states | `h_0` initial form (factored) | rational point |
|--:|--:|--:|:--|:--|
| 15 | 2 | 987 | `-729 E (9E^3 + 8X^5)` | `(X,E)=(3,-6)` |
| 13 | 2 | 438 | `-E^2 (2560 D^5 + 6561 E^2)` | `(D,E)=(-10, 16000/81)` |
| 12 | 2 | 343 | `-2187 (3E^4 - 4S^5)` | `(S,E)=(12,24)` |
| 5 | 2 | 272 | `5832 X^5 (DX - E)` | `(D,X,E)=(1,-8,-8)` |
| 1 | 3 | 180 | `-2560 D^5 (DX - E)^2` | `(D,X,E)=(1,-8,-8)` |
| 9 | 3 | 193 | `-8 D X^2 (2D^3+27X^2)(160D^3-27X^2)` | `(D,X)=(-6,-4)` |
| 8 | 2 | 126 | `2187 S^2 (4S^3 + X^4)` | `(X,S)=(4,-4)` |
| 0 | 4 | 103 | `12 S^2 (4D^2+9S)^2 (5D^2+9S)` | `(D,S)=(3,-5)` |
| 18 | 26 | 13 | full `h_0` tie (26 monomials) | `(D,X,S,E)=(3/8,-8,-8,-6)` |

The `degree_tie_drop` depths run from 1 to 17 (`L0_tie_depth_histogram` in the
JSON; 2062 of 3102 states demand depth 17, i.e. the `h_0` leading form must
vanish to order 17). Depth 1 is only the first equation of each such tower — this
is precisely why the depth-1 layer, though sound, closes nothing.

## 4. Worked example systems

**CONSTRAINT, `a=11` T1, `(deg d2, deg d1, deg sigma, deg e) = (5,2,10,11)`**
(`deg_E = 0`). Two level-0 obligations. The `h_0` tie (tropical max 50, dropped to
46 to meet the close) has support `{d2^6 sigma^2, d2^4 sigma^3, d2^2 sigma^4,
sigma^5}`; depth-1 equation

```
960 D^6 S^2 + 6048 D^4 S^3 + 12636 D^2 S^4 + 8748 S^5 = 0
   = 12 S^2 (4D^2 + 9S)^2 (5D^2 + 9S) = 0,
```

satisfied by `(D,S) = (3,-5)` (`5D^2 + 9S = 45 - 45 = 0`). The close obligation is
`lc(E)^21 lc(h_0) + (-1024/3315) lc(r_0) = 0` (`21*0 + 46 = 4 + 42`, a tie), always
solvable for `lc(r_0)`. Depth `delta = 4`, so the full obligation additionally
demands the next three convolution coefficients vanish (beyond depth 1). Survives
kills-on -> CONSTRAINT.

**Soundness demonstrator, `a=11` T1, `(6,5,sigma=0,11)`** (`deg_E = 0`,
`sigma` identically zero). Here the level-5 tropical tie is **exactly the C08
support** `{d2^2 d1^2, d2 d1 e, e^2}` (`6X^2D^2 - 9XDE - E^2`, discriminant
`105`, not a square in `Q(sqrt(17))`) and the level-4 tie is **exactly C20**. If
the chain were forced to cancel `h_5` or `h_4` on those supports the state would
die. It is not: `deg_h_options(5,...)` and `deg_h_options(4,...)` both return the
maximum with no obligation, the chain uses `h_5, h_4` at full degree, and its only
obligations are the two level-0 records — identical with kills off and on. Hence
CONSTRAINT, not KILL. (Reproduced by `alt_residue_congruences_verify.py` check C.)

## 5. Honest / ambiguous points — [judgment]

- **[judgment] K1 — 0 depth-1 kills is a real, sound result, not an engine gap.**
  The kill test is the audited kills-on re-run over all 3102 states; it is a
  strict over-approximation of the finite-place layer (it forbids only the two
  arithmetically-proven drops). It changes nothing. Every claim of "no admissible
  solution" is confined to C08/C20, which are never a *required* obligation here.

- **[judgment] K2 — the forbidden supports appear but do not fire (the crux).**
  C08/C20 occur as tropical ties in 703 level-occurrences (523 states), all
  `sigma=0` T1. A kill lemma forbids a *drop*; these ties sit at the max with no
  drop required, so they impose nothing. Treating a non-obligatory forbidden tie
  as a kill would be unsound — this is the exact hazard the task flags, and I
  reject it. `n_with_forbidden_tie` is reported per branch so the audit is visible.

- **[judgment] K3 — all obligations are level-0, so C01-C23 (levels 4-6) never
  enter as obligations.** The intermediate `g`-side dominance (`ALT_REGIME_INF.md`
  worked table) forces `h_6..h_1` at their maxima across the whole surviving
  window; I verified no survivor carries a level-1..6 obligation. The residue
  congruences that matter here are the **level-0** `h_0` initial forms plus the
  close — a family (19 hypersurfaces) not previously catalogued, computed here
  directly from `MONOMIALS[0]`.

- **[judgment] K4 — CONSTRAINT is honest about depth.** Every survivor's `h_0`
  tie has depth up to 17; the depth-1 equation is necessary but far from
  sufficient. I record the depth so no reader mistakes "passes depth 1" for "is a
  counterexample." Closing a branch needs the deeper convolution coefficients of
  the level-0 tie **and** the `sigma = 4 d0 - d2^2` coherence
  (`RESIDUE_LEMMAS.md` §2) coupling `S` to `D` and `lc(d0)` — both strictly beyond
  this layer.

- **[judgment] K5 — the bottom close never kills.** With no cap on `r_f`
  (`ALT_REGIME_INF.md` [judgment]) `lc(r_0)` is a free nonzero unknown, so the
  linear close `lc(E)^21 lc(h_0) + lc(u) lc(r_0) = 0` always solves. It is retained
  as a CONSTRAINT label, never a kill.

- **[judgment] K6 — WHOLE-BRANCH KILLS: none (ENGINE+LEMMA-PROVEN).** No branch
  loses all its states at depth 1 (every branch is 100% CONSTRAINT). This is
  proven by the same sound kills-on test, cross-checked on a fresh 120-state
  sample by the verifier. **No whole-branch kill is claimed and none exists at
  depth 1.**

## 6. Verification (`alt_residue_congruences_verify.py`, PASSES)

Independent of the engine (re-derives from `cascade_engine`/`alt_inf_sweep`,
re-parses supports, never trusts the JSON forms):

- **A** — `lc(u) = c lc(q) = -1024/3315` (nonzero) from the source `q`; the close
  is linearly solvable, never a kill.
- **B** — CONSTRAINT state `(5,2,10,11)`: exactly two level-0 obligations; the
  recomputed `h_0` form equals `12 S^2 (4D^2+9S)^2 (5D^2+9S)`; point `(D,S)=(3,-5)`;
  survives kills-on.
- **C** — soundness state `(6,5,sigma=0,11)`: C08(L5)+C20(L4) forbidden ties
  present, C08 discriminant square class `105 != 17`; yet `deg_h_options` offers
  only the max there and the state survives identically kills off/on -> not a kill.
- **D** — census: `0` killed / `3102` constrained / `0` whole-branch kills; all 19
  hypersurface points re-verified; a fresh random 120-state kill-test matches the
  JSON classification.

All checks PASS.
