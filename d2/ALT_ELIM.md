# ALT_ELIM — the named cure applied to the alt-bridge wall (2026-07-24)

**Status: elimination CERTIFIED; the DECISIVE CONTROL still COSTs. The
dm4-eliminated alt bridge does NOT reproduce the known a12 kill — the wall is
DEEPER than spare count, exactly as `R9_SYMBOLIC.md` found for the standard
window. Zero kills; nothing to sweep, nothing to audit.** New files:
`alt_elim.py`, `alt_eliminated_system.json`, `alt_elim_results.json`, this doc.
Read-only on every landed module; nothing committed.

## 0. The job

`ALT_BRIDGE.md` §4a named the cure for the alt-bridge cost wall (both pilots
COST on four engines, the control not computationally reproduced): **eliminate
the 66 spare window unknowns symbolically FIRST**, the pattern
`R9_SYMBOLIC.md`/`r9_symbolic_elim.py` certified for the standard window
(dm4 eliminated via linear pivots, sub2 45→28). This lane applies it to the
alt (sub1) bridge and re-runs the two pilots on the reduced system.

## 1. The elimination + greedy pivot census (`alt_eliminated_system.json`)

The four pre-resultant generators `G1,G2,G3,G5(=G5body+Phi)` are **loaded once**
from `generators.json` and are **t-regime-independent**: the alt (sub1) bridge
and the standard (sub2) bridge share the *same* generators; only the spare-poly
degree caps differ (sub1 18/21/24 vs sub2 12/14/16). Hence the certified dm4
elimination transfers verbatim. Greedy pivot census (probe over the loaded
generators, `pivot_census` in the JSON):

| spare | linear in | pivot | guaranteed nonzero? |
|---|---|---|---|
| **dm4** | G1 (deg 1) | `3·dm1 = 3·e` | **YES** (e is a known nonzero state poly) |
| dm3 | G1 / G3 (deg 1) | `3·dm2` / `3·dm4` | no (spare-valued) |
| dm2 | G1 (deg 1) | `3·d2·dm1 + 3·dm3` | no (spare-valued) |

dm2, dm3 additionally appear only *quadratically* in G2, G5 (constant pivots but
degree 2). **⇒ dm4 is the UNIQUE spare admitting a guaranteed-nonzero-pivot
linear elimination.** Using `dm1 = e` as the multiplier,

```
H2 := dm1·G2 − dm2·G1            (weight 228)
H3 := dm1·G3 − dm3·G1            (weight 240)
H5 := dm1·G5 + (d0·dm1 + d1·dm2 + d2·dm3)·G1   (weight 264)
```

are **dm4-free elements of ⟨G1,G2,G3,G5⟩**. Cofactor certificate re-verified by
exact re-expansion (membership residuals `0`), weighted homogeneity checked, and
the G1 divisibility-lemma identity (`monic(e) | dm2·dm3` on the variety) verified
exactly — ALL SYMBOLIC-ELIMINATION CHECKS PASSED (shared verbatim with
`r9_eliminated_system.json`; the elimination is regime-independent).

**Elimination census (sub1).** Of the **66** spare scalars (dm2:19, dm3:22,
dm4:25), the **25 dm4 scalars are eliminated symbolically** at zero soundness
cost (**66 → 41**). No further guaranteed-nonzero-pivot linear elimination
exists (census above). Verdict semantics (asymmetric, from `R9_SYMBOLIC.md`):
`UNIT` ⇒ state killed (PENDING AUDIT); `PROPER` ⇒ inconclusive (the H-system is
strictly weaker than the full bridge — it drops G1's "dm4 is a capped
polynomial" coupling); `TIMEOUT` ⇒ COST.

## 2. Instantiation route (native Singular, not sympy)

The reduced H-system is **cubic in the 41 remaining spare symbols** (H2 carries
`dm2²·dm3`) and multiplies `e` up to the **4th power** (H3 has `dm1⁴/2`) over
the two marked roots — so the y-coefficients are cubic forms in 43 symbols with
r-symbolic coefficients. Expanding them **in sympy is intractable** (a single
state's build did not finish in 30 min CPU, RSS-bounded, pure symbolic swell).

The build is therefore pushed **entirely into Singular** (C-speed): the abstract
generators (loaded, `dm4:=0` — valid since dm4 cancels identically in every H)
and the small state polys are *defined* as Singular polynomials, the certified
combinations `H2 = dm1·G2 − dm2·G1` etc. are formed, and their y-coefficients
are extracted with `coeffs(H_i, y)` as ideal generators — the exact reduced
system, only the arithmetic moves. Emission clears the `y^N/M` parser trap
(memory: gm^8/N) by writing every rational coefficient *first* in its term.
The H-coefficients alone are the sound necessary system (187 equations); the G1
divisibility strengthening is omitted (its y-only division is itself costly and
does not affect kill-soundness). **msolve was not usable** here: its input is a
list of *pre-expanded* polynomials — exactly the intractable sympy step — so the
exact char-0 engine is Singular's own `std` over Q on the `coeffs()`-built ideal.

## 3. Verdicts — BOTH PILOTS COST; CONTROL NOT REPRODUCED (`alt_elim_results.json`)

Protocol per state: numroot mod-p triage (marked roots → numeric q-roots,
3 primes × 45 s, orphan-proof WSL `timeout` + `ulimit -v 8G`) → Singular-native
char-0 `std` over Q (300 s).

| state (deg d2=6) | reduced eqs | numroot mod-p (3p × 45 s) | exact char-0 Q std (300 s) | verdict |
|---|---|---|---|---|
| **CONTROL** `a12_b1110_T2` | 187 (dm4-elim, 41 spare) | TIMEOUT ×3 (45.3 s) | TIMEOUT (300.3 s) | **COST** |
| **PRIZE** `a11_b3100_T2` | 187 (dm4-elim, 41 spare) | TIMEOUT ×3 (45.3 s) | TIMEOUT (300.2 s) | **COST** |

**THE CONTROL VERDICT: still COST — the reduced bridge did NOT reproduce the
known a12 kill.** `a12_b1110_T2` deg d2=6 is a *known* kill (msolve closed it
exact char-0 on the direct/f31 route, `msolve_bridge_results.json`, making the
branch 8/8). The dm4-eliminated alt-bridge system for the same state times out
on every engine at the stated budgets. **The mechanism is NOT computationally
validated on the alt regime.** This is the identical shape of the
`R9_SYMBOLIC.md` negative result: eliminating the cheapest spare (dm4 — the only
one with a guaranteed-nonzero pivot) shrinks the system (25 fewer unknowns,
66→41) but does **not** move it across the feasibility boundary.

**PRIZE VERDICT: COST.** `a11_b3100_T2` deg d2=6 (the last open state of its 7/8
branch, two marked roots `r1≠r2`) also times out on the reduced system. The alt
frontier is **unchanged** by this lane; the deg d2=6 alt stragglers stay OPEN
with the tie-tower / d2-threshold lanes as status-quo owners.

**SWEEP CENSUS: not reached** (Step 3 runs only if the control kills; it did
not). 0 KILL, 0 PROPER/survival, 2 COST. The "kills PENDING AUDIT" discipline is
vacuous — nothing to hand an auditor.

## 4. What remains — the wall is deeper than spare count

The greedy census is exhaustive: **dm4 is the only sound symbolic linear
elimination** (unique guaranteed-nonzero pivot). Its removal leaves the swell
untouched, which localizes the wall precisely:

* The cost lives in the **dm2/dm3 cubic ansatz** (41 scalars, H is cubic in
  them) coupled to an `e` that **attains the sub1 window cap** (deg 15, the
  cost law `BRIDGE_SWEEP.md` §2/§6 flagged) and to the marked-root
  nonlinearity (one or two adjoined `q`-roots). Neither budget, engine,
  formulation, nor now spare-count moves it.
* The **named next reduction** is the divisibility-forced valuation split
  (`R9_SYMBOLIC.md` §3c): the lemma proven and certified here
  (`monic(e) | dm2·dm3`, cofactor-exact) splits `(dm2, dm3)` into a finite
  case family (per case ~18 spare unknowns instead of 41, each strictly
  smaller and a sound necessary condition under its case hypothesis). This
  lane sets it up (the divisibility identity is in `alt_eliminated_system.json`)
  but does not implement the case split — the direct continuation for any
  swell-tolerant follow-up.

## 5. [judgment] list

- **[J1] The elimination is regime-independent** and shares its certificate
  with `r9_eliminated_system.json`; `alt_eliminated_system.json` re-runs the
  exact-expansion checks and adds the sub1 scalar census + greedy pivot census.
- **[J2] The H-system alone is the sound necessary system** used here
  (`UNIT` ⇒ kill); the G1 divisibility strengthening is proven but omitted
  from the solve (soundness of the kill test is unaffected; a `PROPER` was
  already inconclusive).
- **[J3] Native-Singular instantiation is the exact reduced system** — the
  `coeffs(H_i,y)` route only relocates the arithmetic; the generators are
  loaded from `fsb.gsystem()` (never hand-copied) and `dm4:=0` is valid
  because dm4 cancels identically in every H (certified §1).
- **[J4] All kills PENDING AUDIT** — vacuous this lane (zero kills). Any future
  `UNIT` inherits the status; `kill_system` (the exact char-0 Singular program)
  is recorded in `alt_elim_results.json` when it occurs.
