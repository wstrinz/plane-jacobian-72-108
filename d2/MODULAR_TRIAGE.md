# Modular triage of the program's resistant ideals

**Status: RECONNAISSANCE. Predictions are NOT proofs** and are labelled as such
throughout. A mod-p Groebner verdict is evidence about emptiness/solvability
over the algebraic closure of F_p at a bounded coefficient/tie depth; it is not
a certificate over Q or R.

New, uncommitted files: `modular_triage.py` (builder + runner), this doc, and
`modular_triage.json` (the machine-readable verdict map). READ-ONLY on every
imported module and artifact.

## Method

For each system the exact generators are built in sympy (reusing the landed
`convolution_elim_qsupport`, `convolution_descent`, `phase_f2_scale`, and the
`phase_d_states_sub2` / `batch_convolution_sub2` metadata), coefficients are
reduced modulo p (denominators cleared to integers — never emitted as `/`), and
a Groebner basis is computed over F_p in Singular (WSL, `elim.lib`). Verdict per
system per prime:

* **UNIT** — reduced basis is `1`; the variety is empty over `\bar F_p`.
* **PROPER** — basis `!= 1`; Krull `dim(G)` recorded (cheap).

**Primes:** `10007, 10009, 100019`. The discriminant of
`q = 2048y^4-512y^3+320y^2-240y+195` factors as `2^36·3^2·5^2·13^2·17^3`, so the
**only** bad primes (dividing disc, `lc(q)=2048`, or the fixed denominator
`6630`) are `{2,3,5,13,17}`. All three chosen primes avoid these, keep q
separable, and each splits q with `>= 2` roots so that two-marked-root systems
can be specialized. **Marked roots are specialized to numeric roots of q mod p**
(the cleaner route the brief recommended) rather than adjoining `q(r)`.

**Per-prime Singular budget: 60 s** (TIMEOUT recorded, never silent).

### The one caveat that governs the read

A mod-p Groebner basis sees emptiness over an **algebraically closed** field. A
system that is empty only over **R** (a real / positive-definiteness obstruction)
still returns PROPER and is scored LIKELY-SOLVABLE here. **System 2 is exactly
this case** and must be read with that in mind: "LIKELY-SOLVABLE" means "has a
solution over some field", not "survives over R".

## Prediction rule

| prediction | definition |
|---|---|
| **LIKELY-EMPTY** | UNIT for every prime |
| **LIKELY-SOLVABLE** | PROPER for every prime, consistent dim |
| **MIXED** | verdicts disagree across primes (bad prime or genuine subtlety) |
| **INDETERMINATE** | timeout/skip — no verdict at the budget |

## Headline split (59 subsystems)

| prediction | count |
|---|---:|
| LIKELY-EMPTY | **24** |
| LIKELY-SOLVABLE | **31** |
| INDETERMINATE (timeout) | 4 |
| **MIXED** | **0** |

**No MIXED cases.** Every prime agreed wherever a verdict was returned — no
bad-prime artifacts, no cross-prime dimension disagreement.

## System 1 — R9 z=0..6 q-supported ideals

Ansatz `e = gamma(y+1)^9(y-r)`, `sigma = (y-r)^2 G(y)` (deg G = z), `d1=0`,
`deg d2 <= 4`; 8 master coefficients (deg 251 down), saturated by
`(gamma, g_z, G(r))`. (More generators over-determine a unit ideal, so it
*collapses faster* — 8 coefficients resolve z<=3 in seconds where 6 time out.)

| z | 10007 | 10009 | 100019 | prediction |
|---|---|---|---|---|
| 0 | UNIT | UNIT | UNIT | **LIKELY-EMPTY** |
| 1 | UNIT | UNIT | UNIT | **LIKELY-EMPTY** |
| 2 | UNIT | UNIT | UNIT | **LIKELY-EMPTY** |
| 3 | UNIT | UNIT | UNIT | **LIKELY-EMPTY** |
| 4 | TIMEOUT | TIMEOUT | TIMEOUT | INDETERMINATE |
| 5 | TIMEOUT | TIMEOUT | TIMEOUT | INDETERMINATE |
| 6 | TIMEOUT | TIMEOUT | TIMEOUT | INDETERMINATE |

z=4..6 time out purely on **saturation cost** (11–13 ring variables), not on any
sign of solvability: z=0..3 are cleanly UNIT and the construction is uniform in
z. **Read: the whole R9 column is LIKELY-EMPTY**, confirmed for z<=3, unconfirmed
(cost) for z=4..6.

## System 2 — a8 constant-E gauge stall states (24)

Gauge ansatz `e = gamma(y+1)^8` (gamma a nonzero parameter), generic
`d1,sigma,d2` per state; top ~3 accumulated master coefficients; saturated by
gamma. The accumulated residual is the catalogued `8192 b_i^2 + 9945 gamma^3 s_j^2`
(+ two lower coefficients).

**All 24: PROPER on all three primes → LIKELY-SOLVABLE**, with dimension set by
the free sigma tail: `deg_sigma 5 -> dim 2`, `6 -> 3`, `7,8 -> 4` (identical
across primes and across the d2/d1 sub-block).

**This is the caveat case.** The residual `8192 b^2 + 9945 gamma^3 s^2 = 0`
(gamma != 0) has honest solutions over F_p and over C, so mod-p can only ever say
PROPER. The real obstruction the batch doc conjectured — for `gamma > 0` a
positive-definite sum forcing a degree drop — is **invisible to a mod-p GB**.
Modular triage therefore neither confirms nor refutes these; it confirms they are
**not** killable by any characteristic-blind / closed-field argument.

## System 3 — alt NARROWED/UNOBSTRUCTED reconstruction tie-towers (18)

Reconstructed `d1,sigma,e` (defect-0 divisors), free d2, level-0 tie tower run to
**full tie depth** (mod-p is cheap, so pushed well past the depth-2 cap the
rational run stalled at), Rabinowitsch-saturated leading scalars, roots
specialized numerically.

| state (bid#sup) | field | 10007/10009/100019 | dim | prediction |
|---|---|---|---:|---|
| a11_b0000_T1 #1 (idx1) | Q | PROPER | 4 | LIKELY-SOLVABLE |
| a11_b0000_T1 #1 (idx2) | Q | PROPER | 5 | LIKELY-SOLVABLE |
| a11_b1000_T1 #1 | Q[r]/(q) | PROPER | 5 | LIKELY-SOLVABLE |
| a11_b1000_T1 #5 | Q[r]/(q) | PROPER | 1 | LIKELY-SOLVABLE |
| a11_b1100_T1 #1 | Q[r1,r2] | PROPER | 1 | LIKELY-SOLVABLE |
| a12_b0000_T1 #5 | Q | PROPER | 2 | LIKELY-SOLVABLE |
| a12_b0000_T1 #9 | Q | PROPER | 2 | LIKELY-SOLVABLE |
| a11_b1110_T1 #1 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a11_b3000_T1 #5 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a11_b3000_T1 #9 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a12_b1110_T1 #13 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a12_b3000_T1 #13 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a14_b1000_T1 #13 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a11_b1111_T2 #13 | Q | UNIT | — | **LIKELY-EMPTY** |
| a11_b3100_T2 #12 | Q[r1,r2] | UNIT | — | **LIKELY-EMPTY** |
| a12_b1110_T2 #12 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a13_b1000_T2 #11 | Q[r]/(q) | UNIT | — | **LIKELY-EMPTY** |
| a11_b1111_T1 #17 | Q | TIMEOUT | — | INDETERMINATE |

**10 LIKELY-EMPTY, 7 LIKELY-SOLVABLE, 1 INDETERMINATE.** This is the most
actionable finding: **pushing the tie tower to full depth mod-p flips 10 of the
18 "NARROWED" states to UNIT.** Their NARROWED label was an artifact of the
depth-2 cap the rational GB stalled at, not evidence of survival — they are
strong full-force kill candidates. The 7 that stay PROPER at full depth are the
genuine survivors: every `b0000` state (no q-support, free scalars) plus the low-
support `a11_b1000` and the two-root `a11_b1100` — these carry real residual
freedom and look genuinely UNOBSTRUCTED.

## System 4 — sub2 T2 pattern-B tie states (sample of 10, a7/a8 cells)

Master-coefficient systems built like R9 (T2: `d1=0`,
`e = gamma(y+1)^{a_t} prod(y-r_i)^{b_i}`, `sigma = prod(y-r_i)^{2b_i} G`), 6
master coefficients, saturated, single/no marked root specialized numerically.
Sample restricted to single-root (`b1000`,`b3000`) and no-root (`b0000`) cells to
keep the reconstruction faithful; multi-distinct-root pattern-B cells were not
sampled (deferred).

| state | root | 10007/10009/100019 | prediction |
|---|---|---|---|
| a8_b0000 dd2=0 dsig3 | none | UNIT | **LIKELY-EMPTY** |
| a8_b0000 dd2=1 dsig3 | none | UNIT | **LIKELY-EMPTY** |
| a8_b1000 dd2=0 dsig5 | 1 | UNIT | **LIKELY-EMPTY** |
| a8_b1000 dd2=1 dsig5 | 1 | UNIT | **LIKELY-EMPTY** |
| a8_b1000 dd2=2 dsig5 | 1 | UNIT | **LIKELY-EMPTY** |
| a7_b1000 dd2=0 dsig3 | 1 | UNIT | **LIKELY-EMPTY** |
| a7_b1000 dd2=1 dsig3 | 1 | UNIT | **LIKELY-EMPTY** |
| a7_b1000 dd2=2 dsig3 | 1 | UNIT | **LIKELY-EMPTY** |
| a7_b3000 dd2=0 dsig7 | 1 | UNIT | **LIKELY-EMPTY** |
| a7_b3000 dd2=1 dsig7 | 1 | UNIT | **LIKELY-EMPTY** |

**All 10: LIKELY-EMPTY.** The marked-root q-support master-coefficient systems
for the a7/a8 T2 cells look uniformly killable — the same signature as R9 z<=3.
(Ansatz is a faithful generalization of the R9 recipe, not a verbatim landed
construction; treat as reconnaissance.) Contrast with System 2: the q-support /
marked-root structure is what drives the emptiness, exactly the program thesis.

## Prioritized recommendations

### Hit with full rational force (LIKELY-EMPTY — kills look reachable)

1. **System 4 — the a7/a8 T2 pattern-B single-root cells (all 10 UNIT).** Highest
   value: a whole new region of T2 states signalling killable, previously
   uncharacterized. Extend the single-root exact q-support run (the R9 machinery)
   to a7/a8 `b1000`/`b3000`, then the multi-root `b1100`/`b1110` cells not
   sampled here.
2. **System 3 — the 10 NARROWED states that flipped to UNIT** (a11_b1110,
   a11_b3000×2, a12_b1110, a12_b3000, a14_b1000, and all four T2:
   a11_b1111/a11_b3100/a12_b1110/a13_b1000). Re-run the rational saturated tie
   tower **at full tie depth** (not the depth-2 cap) — mod-p says the unit ideal
   is there. Cheapest first: the over-Q `a11_b1111_T2`.
3. **System 1 — R9 z=0..3 (confirmed UNIT).** Finish/attach the exact rational
   CONTRADICTION certificate. For z=4..6 (timeout, not solvable-looking), raise
   the coefficient count and split saturation to beat the cost wall.

### Pivot to characterization (LIKELY-SOLVABLE — no closed-field kill exists)

4. **System 2 — all 24 a8 constant-E states.** Do **not** spend rational GB force
   here: mod-p proves no characteristic-blind kill exists. The only route is the
   **real / sign** argument (`8192 b^2 + 9945 gamma^3 s^2` positive-definite for
   gamma>0), i.e. the sign-split refinement the batch doc already flagged.
   Characterize, don't Groebner.
5. **System 3 — the 7 genuine survivors** (b0000 ×4, a11_b1000 low-support,
   a11_b1100 two-root, a12_b0000 ×2). PROPER at full depth with clean positive
   dimension → genuinely UNOBSTRUCTED at the tracked level; these need the deeper
   d2-vs-tie-depth argument or d2's own divisor, not more of this tie tower.

### Unresolved by this pass

* System 1 z=4..6 and System 3 `a11_b1111_T1 #17`: INDETERMINATE (60 s
  saturation timeout). A longer budget or a leaner saturation would likely
  resolve all four toward their column trend (EMPTY for z4..6, unknown for the
  a11_b1111_T1 flagship — its 7 Q-siblings were KILLED in phase_f2_scale, so
  EMPTY is the expectation).
