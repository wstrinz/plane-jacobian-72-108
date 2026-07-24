# R9_SYMBOLIC — symbolic spare-unknown elimination for the deg_e = 10 cost wall

**Status: COMPLETE — elimination certified; sweep verdict NEGATIVE (the
honest, structural kind).** Files: `r9_symbolic_elim.py` (one-time
elimination + `r9_eliminated_system.json`), `r9_symbolic_sweep.py` (per-state
driver + `r9_symbolic_sweep.json`, run log `r9_sweep_run.log`).

**Headline: eliminating dm4 (45 → 28 spare unknowns, + 9-10 divisibility
equations) does NOT break the deg_e = 10 cost wall — on either the
marked-root R9 column or the r-free batch T2 family. The wall is driven by
the dm2/dm3 spare ansatz coupled to a window-cap e, not by dm4. The named
next reduction (§5) is the divisibility-forced valuation split, which cuts
28 spare unknowns to 18 across a finite case family.**

## 1. What was eliminated, and the certificate

`BRIDGE_SWEEP.md` §6 named the cure for the R9 z≥1 / a10 deg_e=10 cost wall:
"a structurally reduced system (e.g. eliminating the spare unknowns
symbolically first), not a bigger budget." This lane builds exactly that.

Every G-system generator is **linear in dm4 = d_-4** (verified by
`r9_symbolic_elim.py`, pivot coefficients extracted, `pivot(G1) = 3·dm1`).
Since `dm1 = e` is a known nonzero polynomial in every cascade state, the
cross-multiplied combinations

```
H2 := dm1·G2 − dm2·G1          (weight 228, 5 terms)
H3 := dm1·G3 − dm3·G1          (weight 240, 6 terms)
H5 := dm1·G5 + (d0·dm1 + d1·dm2 + d2·dm3)·G1     (weight 264, 12 terms)
```

are **dm4-free elements of the ideal ⟨G1,G2,G3,G5⟩** — membership is by
construction, with the cofactors recorded in `r9_eliminated_system.json` and
the identities re-verified by exact expansion (ALL SYMBOLIC-ELIMINATION
CHECKS PASSED). Weighted homogeneity of each H is verified (228/240/264
under the T3 window weights).

**The divisibility lemma.** `G1 = 0` rearranges to
`3·dm2·dm3 = −dm1·(3/2·d1·dm1 + 3·d2·dm2 + 3·dm4)` (identity verified
exactly), all factors polynomials — hence on the G-variety

```
monic(e) | dm2·dm3   in Q[y].
```

Per state this contributes deg(monic(e)) further small quadratic equations
`rem(dm2·dm3, monic(e), y) = 0` — recovering most of G1's content without
dm4. For R9, monic(e) = (y+1)⁹(y−r); for the batch deg_e=10 T2 states,
(y+1)^a (the tail's unknown-coefficient factor is not exploited — flagged,
sound to omit).

## 2. The per-state reduced system and its verdict semantics

```
REDUCED(state) = [ y-coefficients of H2, H3, H5 on the stripped ansatz,
                   dm2, dm3 as bounded polys (sub2 caps 12/14; 28 unknowns) ]
               + [ rem(dm2·dm3, monic(e)) = 0 coefficient-wise ]
               + [ q(r) = 0 and the state saturations (marked-root states) ]
```

dm4's 17 spare coefficients never appear (45 → 28 spare unknowns).

**Asymmetric semantics (the honest part):** every REDUCED equation is a
proven necessary condition, so `UNIT` ⇒ **state killed** (candidate kill,
PENDING AUDIT). But REDUCED is *weaker* than the full bridge (dropping G1
forgets that dm4 must be a polynomial of capped degree ≤ 16), so `PROPER`
is **inconclusive** — never a survival signal. TIMEOUT remains pure cost.

The instantiation route (`r9_symbolic_sweep.build_reduced`) reuses
`fsb.augment`'s single expansion (its Gpolys) and forms the H's in
y-coefficient space by convolution with the small cofactor lists; the exact
cancellation of every dm4 coefficient M0..M16 is **asserted per
coefficient**, which independently re-checks the elimination identity on
every state. Emission/triage/exact-Q reuse `bridge_sweep`'s machinery
verbatim (same integer-cleared Singular, Rabinowitsch saturation, minpoly
and numroot formulations); the only substitution is an orphan-proof runner
(WSL-side `timeout` + `ulimit -v 8G` — a Windows relay timeout can no
longer leave a WSL Singular alive).

## 3. Results — census (machine record `r9_symbolic_sweep.json`)

| state | system (eqs/vars; full bridge was ~122/53+) | mod-p triage | exact-Q | verdict |
|---|---|---|---|---|
| R9 z=1 | 136 / 37 | numroot 3×TIMEOUT (45–60 s) | minpoly TIMEOUT 300 s (+ q(r) fallback TIMEOUT 300 s) | **COST** |
| R9 z=2 | 136 / 38 | numroot 3×TIMEOUT (45 s) | minpoly TIMEOUT 300 s | **COST** |
| R9 z=3 | 136 / 39 | numroot 3×TIMEOUT (45 s) | skipped by policy (no UNIT evidence) | **COST** |
| R9 z=4..6 | — | NOT ATTEMPTED (budget; monotone column trend after z=1,2,3) | — | honest truncation |
| batch a9 T2 deg_e=10, 8 cheapest (d2 ∈ {∅,0,1} × σ ∈ 2..4) | 134 / 33–36 | symbolic 3×TIMEOUT (45 s) each | skipped by policy | **COST ×8** |
| batch, remaining 82 of 90 | — | NOT ATTEMPTED (limit-8 cheapest-first sample; uniform COST) | — | honest truncation |

**The 300 s control (the decisive datum).** The full-bridge sweep's
deg_e=10 failures were at 300 s mod-p budgets; the sweep above used 45 s.
To compare like with like, the first batch state's reduced system was
re-triaged at **300 s on p=10007: TIMEOUT** (`control_300s` in the JSON).
The reduced system fails at the same budget where the full bridge failed —
the elimination shrank the system (17 fewer unknowns, comparable equation
count) but did **not** move it across the feasibility boundary. Zero kills,
zero PROPER: every verdict is pure Groebner cost, no survival signal
anywhere.

## 3b. What this negative result buys

1. **The wall is now localized.** dm4 was the cheapest spare unknown to
   eliminate (linear in every generator) — and removing it changes nothing.
   The swell therefore lives in the dm2/dm3 polynomial ansätze (28 scalars)
   coupled to the window-cap e (deg 10). Budget, engine, and formulation
   were already ruled out by BRIDGE_SWEEP §6; this lane rules out the easy
   structural reduction too.
2. **The elimination machinery is sound, cheap, and reusable** — the
   certified H-system + per-state convolution instantiation (validated
   byte-identical against the fsb.augment route on R9 z=1) is the platform
   any deeper reduction will instantiate through.

## 3c. The named next reduction (concrete, sound, not yet implemented)

The divisibility lemma is much stronger than the remainder equations used
here. For R9, `(y+1)^9 (y−r) | dm2·dm3` gives by valuation additivity

```
v_{y+1}(dm2) + v_{y+1}(dm3) >= 9     and     v_r(dm2) + v_r(dm3) >= 1,
```

so (dm2, dm3) splits into the finite case family

```
dm2 = (y+1)^i (y−r)^j · A,        deg A <= 12 − i − j
dm3 = (y+1)^{9−i} (y−r)^{1−j} · B,  deg B <= 14 − (9−i) − (1−j)
      i = 0..9,  j = 0..1
```

— **18 spare unknowns per case (down from 28), 20 cases**, each strictly
smaller than anything tried in this lane, and each a sound necessary
condition under its case hypothesis (the cases exhaust the variety). The
same split applies to batch states through `(y+1)^a`. This composes with
the H-system unchanged and is the direct continuation.

## 4. [judgment] list

- **[inherited] k=6,7 window caps** for dm2, dm3 — the same flagged
  `T3_WINDOW_AUDIT` extension every bridge kill rests on (k=8/dm4 is no
  longer used by this lane at all — one judgment row *retired* for these
  systems).
- **[judgment] Batch divisibility uses only the (y+1)^a factor** of e; the
  unknown-coefficient tail also divides dm2·dm3 but is not imposed
  (non-monic division would need lc-clearing; omitted, sound).
- All kills are **PENDING AUDIT**; each kill record stores the
  integer-cleared generators + saturation factors, sufficient to re-emit
  the exact system for a cofactor certificate
  (`kill_system` field of `r9_symbolic_sweep.json`).
