# CASE_COMPILER — pilot (Lane I, 2026-07-23)

**What this is.** The first working slice of the INDUCTIVE_PROGRAM.md
architecture item: a compiler from GGV5/GGHV chain data to a normalized CASE
DOSSIER. Input: a family id + `j` (GGV5 `v11<=35` length-1 table) or a special
externally-sourced tag (`GGHV_72_108`, `GGHV_108_144`). Output: canonical JSON
(`case_dossier_<tag>.json`) plus a human-readable rendering.

**What it is NOT.** No new mathematics. Every instantiated field restates a
landed, checked fact with its source; every schematic field describes structure
and names the computation that would instantiate it; every conditional
inference carries a raised flag. The compiler *packages* the seven-point corner
law (PHI_F14.md), the kappa = t−2 chart lemma (PHI_CORNER4.md), the Galois
transfer rule (GALOIS_LIBRARY.md §4), and the (72,108) machinery inventory —
it does not extend them.

Checker: `case_compiler_verify.py` — **71 checks** — validates the three pilot
dossiers against known landed values (audited (72,108); derived (75,125) = F2
j=1; derived F9 j=0 (56,84)), the law on all seven landed points, the Galois
routines on independent witnesses, survey-wide Diophantine identities and
conjectural-flag firing, and canonical-JSON determinism.

## Dossier schema (`case-dossier-v1`)

| section | contents | grounding |
|---|---|---|
| `case` | tag, family, j, (m,n), (a,b), degrees, v11, A0, A1, A0', k, chain length, source, Diophantine check | GGV5 table rows are Diophantine-checked (`(m+n)qk − n(ql−p) = k`); GGHV entries are marked externally sourced, identity N/A |
| `corner_signature` | a, b, t, kappa=t−2, a0, q, e=b−a+1, r=a0−q−1, gap=(q−1)−a0/t, dg=a0−q, N=a[t(a+b−1)+1]−2b | kappa=t−2 structural (PHI_CORNER4); gap/dg mini-lemmas (PHI_F14); non-integral gap is preserved as a fraction and flags the law undefined |
| `phi_prediction` | law signature (deg, ord_y, mult_{y+1}, cofactor deg); regime + grounded flag; conjectural flag with reasons; derived/audited reference when the point is landed | unified seven-point law (PHI_F14); the compiler ASSERTS agreement with every landed point it knows |
| `presentations` | **both** presentations with per-field `instantiated` markers (see boundary below) | INDUCTIVE_PROGRAM.md architecture amendment |
| `galois_transfer` | forcing-polynomial candidate + rationale; Galois label, disc, disc class; C08/C20 verdicts; conditionality status | GALOIS_LIBRARY.md §4 two-line rule |
| `transfer_inventory` | per-mechanism status: AS-IS / PARAMETRIC / METHOD-ONLY, with rationale and source | this file, kept with the compiler |
| `judgment` | the case's assembled conditionality list (unreduced polygon, N-formula slice status, forcing-poly identification, regime, chart scope) | mirrors the PHI_* docs' judgment style |

## The instantiated-vs-schematic boundary (honesty contract)

**Fully instantiated only for (72,108)** (the audited home case):
master-identity tower `(D,F,E,s) = (31,7,21,3)`; G-system generators
`G1,G2,G3,G5body+Phi` over `Q[d~2,d~1,d~0,e,r,s,dm4,Phi]`, weights
156/168/180/204, ~122-equation bridge; window caps 12k/15k/14k **PROVEN**
(WINDOW_CAPS.md); `Phi_stripped = c·t^30·q4` with `c = −1/6630` — and the
compiler records that this equals the corner-law object `Phi_full/(lc·y^204)`
exactly (FULL_SYSTEM_BRIDGE.md §1).

**Instantiated for every standard-chart case** (from corner data + the law):
the full corner signature; the predicted Phi signature; the stripped-Phi shape
`t^mult · u_gap · H^mult` with `H = (y^dg+1)/(y+1)` (dg odd) and
`deg u_gap = gap`; the forcing-polynomial candidate and its exact Galois data.

**Schematic everywhere else** — emitted as structure + the instantiating
computation, never as numeric guesses: the tower data `(D,F,E,s)` off (72,108)
(requires the case's cascade/C-series build; corner-144 recurrence is
"skeleton yes, numerics no"), G-weights and bridge equations (require the
case's D-transform), window-cap slopes (case degree arithmetic; the symbolic-
in-k induction METHOD transfers), and the cascade depth constant.

**Flags that fire automatically** (updated 2026-07-23 after the F7/composite/
ζ-tail landings): regime `gap>0 & r>0` is now GROUNDED — the F7 test came back
DIFFERS and the RAMIFIED law replaced the old conjecture (PHI_F7.md; refined
by the μ-ladder, ZETA_TAIL.md), so `law_signature` branches there instead of
flagging; `k=2` rows still flag (N-formula unverified); `A0' != (1,0)` still
flags with an updated reason (chart now SETTLED, κ=t−2 extends per
COMPOSITE_CHARTS.md; the residual conditionality is the ζ>0 commuting-tail
break, models enumerated in ZETA_TAIL.md); dg even no longer means UNKNOWN —
the residual is the ramified shape (PHI_F7.md; parity scoping per
REVIEW_ZETA_MU.md); non-integral gap (law undefined as stated) still flags.

## Galois transfer as compiled

The candidate forcing polynomial follows the INDUCTIVE_PROGRAM correspondence
(audited disc-17 quartic at (72,108); residual `H` — the 10th cyclotomic — at
(108,144)): gap-regime unit cofactor where landed (registry), residual `H`
for gap=0 with dg odd, UNKNOWN otherwise. Labels are computed exactly
(resolvent cubic + discriminant + factorization over `Q(sqrt(disc-class))`
for the C4/D4 split; constructive order-4 automorphism cross-check for the
cyclotomic in the verifier). Verdicts per the two-line rule: kill transfers
iff the splitting field misses `sqrt(105)` / `sqrt(170)`; D4/V4
witness-decided. Off (72,108) every verdict is stamped **CONDITIONAL** on the
case's residue layer reproducing the shape library (GALOIS_LIBRARY judgment
G1) — that residue analogue is underived off the home case.

Pilot results: (72,108) S4/17 → both kill (AUDITED); (75,125) C2/−3 → both
transfer-eligible (CONDITIONAL); F9 C4/5 → both transfer-eligible
(CONDITIONAL).

## Transfer inventory highlights

AS-IS: residue library + Galois descent (two-line check), corner law (within
the standard-chart length-1 class), modular triage. PARAMETRIC: cascade
engine, T2 squeeze, full-system bridge, window caps, divisor-lemma engine.
METHOD-ONLY: S-unit corner kills. Details + caveats in every dossier's
`transfer_inventory` section.

## Next steps this pilot makes concrete

1. **Instantiate a second case's deep structure** (the real T1/T2 test):
   build (75,125)'s C-series/D-transform from polygon data and fill the
   schematic fields — the dossier now states exactly which fields those are.
2. **F7 (Lane H)** flips the last conjectural regime row if it matches.
3. **Residue-layer analogue at one non-home case** would upgrade the Galois
   verdicts from CONDITIONAL and is the named blocker for layer-2 transfer.
