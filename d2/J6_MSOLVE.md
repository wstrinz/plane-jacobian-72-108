# J6_MSOLVE — msolve pass on ALT_HUNT's 4 open [J6] states

**Date:** 2026-07-24. **Status: ALL 4 KILLED (PENDING AUDIT).**
**Headline: 49/49 fully-forced HUNT states are now dead and all 17 HUNT
cells are CLOSED(all0).** The [J6] wall was pure engine cost, as ALT_HUNT.md
§[J4] claimed: msolve's F4 + multi-modular solver decides in ~1 s the same
systems on which sympy Buchberger (600 s) and Singular `std` over Q (240 s,
re-checked here) both blow up.

New files (uncommitted): `j6_msolve.py` (runner), `j6_msolve_results.json`
(full record incl. mod-p corroboration), this doc. READ-ONLY on every
existing artifact; nothing committed.

## 1. What was attacked — exactly what Lane C recorded

The 4 OPEN states of `alt_hunt_results.json` (2 windows × {`a9_b1000_T1`
deg-6-d1 #state3, `a8_b1100_T1` deg-6-d1 #state0}, each with exactly one
admissible split after Galois dedup). Replay guarantees, enforced by
assertion before any msolve call:

- `reconstruct_general(case, combo)` (Lane C's own committed code path) was
  replayed on the recorded combo and its polys asserted **string-identical**
  to the recorded `polys` — 4/4 verbatim matches;
- where Lane C recorded accumulated master coefficients (`sub1:a9` depth 2),
  the replayed walk was asserted to reproduce those exact strings;
- the generator set per depth n is gens[..n] + class relations +
  `w·Π(scalars)−1` — identical to `alt_hunt_depth2.kill_test_record`; only
  the GB engine differs.

## 2. Census

| state | vars | verdict | kill depth | msolve walls | state wall |
|---|---|---|---|---|---|
| `sub2:a9_b1000_T1_sz1_dz1_gz-#state3` | E,X,c0_0,w | **KILLED** | 3 | 1.9 s (d2), 0.5 s (d3) | 18.9 s |
| `sub2:a8_b1100_T1_sz1_dz1_gz-#state0` | E,X,c0_0,c0_1,w | **KILLED** | 3 | 4.5 s (d2), 1.0 s (d3) | 23.2 s |
| `sub1:a9_b1000_T1_sz1_dz1_gz-#state3` | E,X,c0_0,w | **KILLED** | 3 | 2.7 s (d2), 0.6 s (d3) | 21.4 s |
| `sub1:a8_b1100_T1_sz1_dz1_gz-#state0` | E,X,c0_0,c0_1,w | **KILLED** | 3 | 4.1 s (d2), 0.6 s (d3) | 23.1 s |

msolve `[-1]` = no solution in the algebraic closure of Q — the saturated
system is empty, so **no Galois assignment of the admissible split
survives**: the same soundness statement as Lane C's unit-ideal kills.

## 3. Cell accounting

ALT_HUNT.md had 13/17 cells CLOSED(all0); the 4 open states kept the 4 T1
cells `sub{2,1}:a9_b1000_T1_sz1_dz1_gz-` and `sub{2,1}:a8_b1100_T1_sz1_dz1_gz-`
open. With all 4 states killed: **17/17 HUNT cells CLOSED(all0), 49/49
fully-forced states dead.** The s-unit BM-candidate residual layer
(`s_unit_results.json` candidate_rows) is fully closed. Per ALT_HUNT §[J5],
CLOSED(all0) does NOT close the enclosing branches (they carry
non-fully-forced states outside HUNT scope).

## 4. Structural finding: why [J6] blew up

The depth-2 truncations are **genuinely satisfiable** — msolve returns a
0-dimensional parametrization (isolated solutions) at depth 2 for all four
states. So sympy/Singular were not slow at detecting a unit ideal; they were
computing a full nontrivial Gröbner basis of a satisfiable system — the
expensive case. The depth-3 master coefficient is the killer. (Lane C's
depth cap was never the issue; its engine never got past depth 2's GB.)

## 5. Cross-engine corroboration (recorded in the results JSON)

- sympy grevlex Buchberger on the depth-3 killing system: **240 s timeout**
  (blowup is engine-family-wide, not verdict-dependent).
- Singular `std` over Q: **240 s timeout** (same).
- Singular `std` over F_p, p ∈ {10007, 10009, 100019} (good primes): **unit
  ideal, all four states, all three primes (12/12), ≤ 1 s each** — an
  independent second implementation corroborating msolve's char-0 verdict.
  (Campaign triage record: zero mod-p refutations of a char-0 claim, ever.)

## 6. Honest points — [judgment]

- **[J1] PENDING AUDIT**, like all ALT_HUNT kills: msolve is the deciding
  engine here; the named audit route is a Singular `lift` cofactor
  certificate over the same saturated system (`kill_certificate_tools.py`
  pipeline) — note the char-0 `std` timeout above means certificate
  extraction may need the msolve-assisted route or a longer budget.
- **[J2]** Soundness of "state dead" inherits ALT_HUNT [J2]'s split-variety
  argument verbatim (exhaustive splits, class-polynomial bijection, only
  leading scalars saturated).
- **[J3]** The mod-p checks corroborate but do not replace the char-0
  verdict; msolve `[-1]` (multi-modular + rational reconstruction) is the
  load-bearing statement.

## 7. Reproduction

```
python j6_msolve.py            # full pass (~90 s total)
# per-state .ms inputs staged as $HOME/j6_*.ms in WSL; outputs $HOME/j6_*.out
```

Orphan-proofing: msolve ran WSL-side under `timeout` + `ulimit -v 8G`;
final WSL process sweep: clean.
