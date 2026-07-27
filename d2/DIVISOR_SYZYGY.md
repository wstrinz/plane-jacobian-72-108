# DIVISOR_SYZYGY.md — the universal K-syzygy, `e | Φ`, and the support collapse

2026-07-25. Checker: `divisor_syzygy.py` (7/7, `--quiet` in `run_tests.sh`).
Origin: external review (GPT), independently verified here before any use.

## 1. The identity

With `e = dm1`, `R = dm2`, `S = dm3`, `T = dm4`, the canonical G-system
generators satisfy the **exact** polynomial identity

```
2*(G5 + d2*G3 + d1*G2 + d0*G1)  ==  2*Phi - e*(d2*e^2 + 3*e*S + 3*R^2)
```

Every `dm4` term and every `d0/d1/d2` cross term cancels identically. Residual
is **exactly 0** (check C1). Writing `K` for the right-hand side,

```
<G1,G2,G3,G5>  ==  <G1,G2,G3,K>
```

This is **ideal equality**, not a necessary-condition weakening: `K` is an exact
combination of the generators, and conversely
`G5 = K/2 - d2*G3 - d1*G2 - d0*G1` (check C3). So the dense `G5` row can be
replaced by the sparse four-term `K` row with no saturation and no division
by `e`.

## 2. The consequence: `e | Φ`

On every genuine lift `K = 0`, i.e. `2*Phi = e*(d2*e^2 + 3*e*S + 3*R^2)`, hence

> **`e | Φ`**

with `Φ = -(1/6630)*(y+1)^30*q(y)` and `q` the fixed **squarefree** quartic.
Over the algebraic closure this forces `e = γ*(y+1)^a * (squarefree divisor of q)`:

* `e` has **no roots off** `{y=-1} ∪ {roots of q}` — the off-support factor is a
  unit;
* every simple `q`-root divides `e` to order **at most 1**, i.e. `b_i ∈ {0,1}`;
* `deg e = a + Σ b_i` — degree stops being a free Phase-D coordinate.

No degree stratum, no saturation, no `d1=0`, no genericity assumption. It holds
on the whole G-variety.

**Why this is the right shape.** Per `SESSION_HANDOFF.md` THE SPEC, Φ-depth kills
bottom-up and closes a cell only by reaching the top stratum — which happens in
**0 of 220** sub2 cells. This lemma is *not monotone in degree*: it deletes whole
**support cells**, top stratum included. That is the first mechanism in this
program with the shape the spec demands.

## 3. Sub2: the degree is forced

With the certified sub2 caps `deg d2 ≤ 4`, `deg R ≤ 12`, `deg S ≤ 14`,
`deg e ≤ 10`, the RHS has degree `E + max(4+2E, E+14, 24)`. Against
`deg Φ = 34`:

| `E` | 6 | 7 | 8 | 9 | **10** | 11 |
|---|---|---|---|---|---|---|
| max RHS deg | 30 | 31 | 32 | 33 | **34** | 37 |

`E ≤ 9` is impossible; `E = 10` fits exactly, and the cap gives `E ≤ 10`. So

> **`deg e = 10` exactly for every sub2 G-system solution**, and `a + Σ b_i = 10`.

### The open T2 columns collapse 8 → 3

| column | verdict | reason |
|---|---|---|
| `a9_b1000` | **ALIVE** | 9+1 = 10 |
| `a8_b1100` | **ALIVE** | 8+2 = 10 |
| `a7_b1110` | **ALIVE** | 7+3 = 10 |
| `a8_b0000` | dead | `a+Σb = 8 ≠ 10` |
| `a8_b1000` | dead | `a+Σb = 9 ≠ 10` |
| `a7_b1000` | dead | `a+Σb = 8 ≠ 10` |
| `a7_b1100` | dead | `a+Σb = 9 ≠ 10` |
| `a7_b3000` | dead | **`b_i ≥ 2` at a simple `q`-root** |

Note the two death modes are independent tests, not one condition.

**The live saturated pilot `a7_b3000_T2` is empty before any Gröbner run.** It
required a simple root of `q` to divide `e` three times. The in-flight run was
stopped on this basis. Consistent with its own controls, which returned UNIT.

## 3b. VERIFIED SINCE: the census, and a correction to the degree consequence

**`deg e = a + Σ b_i` is NOT automatic from the universe.** `phase_d_states_sub2.json`
records `deg_e = 10` for **3790** states whose `a_t + Σb < 10` — defect-1 `e` with
a free extra root (`PHASE_F2_SUB2.md`). The divisibility repairs this: `e | 2Φ`
confines every root of `e` to `{-1, r_1..r_4}`, forcing defect 0. So the
consequence survives, but it needs that extra step and must not be quoted as
immediate. (The external review gave 3844 for this count; 3790 is what the
authoritative artifact yields under every natural reading — `<10` and `!=10`
agree, and restricting to `b_i<=1` cells gives 3642. The structural point holds;
the number did not reproduce.)

> **CORRECTION, 2026-07-25 (frontier-rebuild lane; `divisor_filter.py` check C8,
> `FRONTIER_REBUILD.md` §6 row `EXTERNAL-DEFECT-COUNT`).** The external **3844
> DOES reproduce** and should no longer be recorded as unreproduced. The two
> numbers count different sets:
> * `deg_e != a_t + Σb` — **3844** states, the count of `e` carrying a free extra
>   root, which is the quantity this paragraph is about;
> * `deg_e == 10` **and** `a_t + Σb < 10` — **3790**, a strict subset.
>
> The 54-state gap is defect-1 `e` at `deg_e = 8` and `9`, in cells
> `a6_b1000_T1` (10 states, `deg_e` 8 over support 7) and `a8_b0000_T1` (44
> states, `deg_e` 9 over support 8). The parenthetical above is right that `<10`
> and `!=10` agree — but only for the `a_t + Σb` test at fixed `deg_e = 10`; it
> does not make 3790 the defect count. Recomputed by
> `python divisor_filter.py --quiet`.

**Census over all 220 sub2 cells** — counted independently from
`phase_d_states_sub2.json`, not from the lemma's own runner:

> **NOMENCLATURE CORRECTION, 2026-07-25.** The 220 objects this table calls
> "cells" are **flag cases**. The cell count `(a_t, b, branch)` in sub2 is **26**;
> the filter kills **18 of 26 cells / 140 of 220 flag cases / 4822 of 7888
> states**. Every state and flag-case figure below is confirmed by the compiler
> (`FRONTIER_REBUILD.md` §2, §6); only the word "cells" was wrong.

| death mode | cells | states |
|---|---|---|
| `b_j >= 2` at a simple `q`-root | 23 | 885 |
| `a + Σb != 10` | 117 | 3937 |
| **total, with NO Gröbner basis** | **140 of 220** | **4822 of 7888** |

Against **0 of 220** before. Surviving columns: `a10_b0000`, `a9_b1000`,
`a8_b1100`, `a7_b1110`, `a6_b1111` (T1), and the three T2 columns above.

**`t^a | R, S, T` is now VERIFIED on the T2 branch** (`d1 = 0`), not by the
sketch but by exhaustive enumeration of every order configuration
`(rho, s, delta2, delta0)` within the sub2 caps, accepting one only if both the
identity and `H3 = 0` can hold. True for `a = 7, 8, 9`. **T1 changes the case
analysis and is NOT covered.** For `a9_b1000_T2` the spare ansatz collapses
`{dm2:13, dm3:15, dm4:17} -> {4, 6, 8}`, i.e. **45 -> 18 = 45 - 3a**.

## 4. NOT verified here

* **`t^a | R,S,T` on the T1 branch.** T2 (`d1 = 0`) is settled in sec.3b by
  exhaustive enumeration; T1 changes the case analysis and is **open**.
* ~~The full frontier regeneration (cascade, Phase-D universes, rollup, proof DAG)
  under the new filter.~~ **DONE 2026-07-25** — `FRONTIER_REBUILD.md`,
  `divisor_filter.py`, `frontier_rebuild.py`. The filter is now a compiler stage
  (`phase_d_states.py --divisor-filter`), a ledger kill source
  (`state_kill_ledger.py --divisor-filter`) and a DAG input
  (`proof_dag.py --ledger state_kill_ledger_divfilter.json`). The sec.3b sub2
  numbers all reproduce; **sub1 is NOT analogous** — its caps make the degree
  forcing vacuous, so only `b_i ∈ {0,1}` and defect-0 act there.
* **The C08/C20 field-scope objection is untouched and independent of all of the
  above.** If residue kills are being read as geometric emptiness on the strength
  of "no point over Q or the q-splitting field", that is unsound — C is a
  Q-algebra, and the repo records real torus points involving sqrt(105),
  sqrt(170). That would REOPEN branches, whereas this lemma closes them; the net
  frontier count is not guessable until both land.

## 5. Independent cross-check that already existed

The generic-fiber lane, starting from the `dm4`-eliminated H-system at `d1 = 0`,
landed `2*(H5 + d2*H3) = dm1*K5` with
`K5 = 2Φ − 3·dm1·dm2² − d2·dm1³ − 3·dm1²·dm3`.

**That `K5` is literally this `K`** (check C4, residual 0). Two lanes, two
starting systems, one object — derived months apart by different routes. The
generic-fiber lane found the `d1=0` shadow; this is the universal form.

## 6. Standing guard

`G5 = Φ + G5body`. A stale `2Φ` transcription was a **real bug** in this repo,
corrected 2026-07-24. Check C2 asserts `coeff(G5, Φ) = 1` precisely because that
stale form would silently break C1 — the syzygy would fail to cancel and the
residual would be `Φ`, which looks like a derivation error rather than a
transcription error.
