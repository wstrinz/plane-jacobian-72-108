# ALT_BRIDGE — the full-system bridge analog for the alternate regime (2026-07-24)

**Status: construction EXISTS (DERIVED + mechanically spot-checked) — the
analog is the sub1-caps G-system. BOTH pilots returned COST (no kill; the
control was NOT computationally reproduced). Zero kills, so nothing to audit
yet. Verdicts + engine census in §4 / `alt_bridge_results.json`.** New files:
`alt_bridge.py`, `alt_bridge_results.json`, `alt_bridge_results_symbolic.json`
(first-run snapshot), `alt_bridge_run.log`, `alt_bridge_retry.log`, this doc.
Read-only on every landed module; nothing committed.

**Collect/finalize pass (2026-07-24) re-verified, cheaply and independently:**
`check_soundness()` → `{homogeneity G1..G5 = 156/168/180/204, phi_identical,
sub1_caps {dm2:18,dm3:21,dm4:24}, alt_states_in_sub1_window (deg e=15, sigma=12,
d2=6 — all at the 3k sub1 caps), complement_exact_mod_q}` all pass (0.8 s); and
`check_fast_builder()` → the fast y-convolution builder emits an equation set
**IDENTICAL** to `fsb.augment` on the synthetic control (5.5 s). The two
computational load-bearing steps of the construction hold. The pilot solves
were NOT re-run — they are documented multi-hundred-second timeouts across four
engine routes (§4); re-running only re-confirms COST.

## 0. The question

`BRIDGE_SWEEP.md` §3 skipped every alternate-regime target: "the bridge is a
standard-window object … outside its soundness domain." The alt frontier (24
open branches) was left to the tie-tower and msolve routes, which hit cost
walls exactly at the two `deg d2=6` stragglers. Does a sound bridge analog
exist for the alt window?

## 1. Answer: the analog exists, and it is the G-system with SUB1 caps

**Claim.** For an alternate-regime state (ALT_REGIME.md: subcase (1),
`a = v_t(e) ∈ 11..15`, `v = 30−3a < 0`), the full pre-resultant G-system
`(G1, G2, G3, G5body+Phi)` with the spare window unknowns `dm2, dm3, dm4`
bounded by the **sub1** stripped caps (`deg ≤ 18/21/24`, 66 spare scalars)
is a sound necessary system — verbatim `full_system_bridge.augment(state,
regime="sub1")` plus marked-root adjunction.

**Derivation** (each point tied to a landed, audited source):

1. **The generators are t-regime-independent.** `G1,G2,G3,G5body+Phi =
   (D~³)_{-1,-2,-3,-5}` after the `(D~²)` substitutions — x-level identities
   derived from `C² = P` and the d₃-killing shift (`T6_SELECTION_AUDIT.md`,
   `regenerate_system.py`). No step consumes the t-adic profile of `e` or
   the standard reduction `F = t^{21a}G`. The subcase enters only through
   `P`'s Newton polygon → the *caps*, not the identities.
2. **f31 ∈ ⟨G⟩ and f31 is regime-independent.** The cofactor certificate is
   pure algebra (`full_system_bridge_verify.py` V2), and `ALT_REGIME.md`'s
   own survival table, row 1: the graded identity "survives verbatim; it is
   window-independent."
3. **The stripped Phi is the identical object.** `ALT_REGIME.md` writes
   `Phi~ = t^30·u, u = cq` — precisely `phi_stripped() = c·t^30·q`
   (asserted identical in `alt_bridge.check_soundness`).
4. **The sub1 window caps hold for alt states.** `WINDOW_CAPS.md` proves
   `ord ≥ 12k`, `deg ≤ 15k` (sub1) for `k = 2..8` from premises
   [P1][P2][P3] alone — bidegree valuation inductions on the C-recursion,
   D-transform arithmetic, and the shift identity. **No t-adic input.**
   The alt layer itself already consumes the `k = 2..5` sub1 caps
   (`deg d1 ≤ 9`, `deg sigma ≤ 12`, `deg e ≤ 15`, `deg d2 ≤ 6` — the very
   caps in `ALT_REGIME.md`'s audit); the `k = 6,7,8` rows carry identical
   provenance and trust tier.
5. **Weighted-homogeneous stripping is exact** — an algebraic property of
   the generators (`G1..G5` weights 156/168/180/204, re-asserted at
   runtime), independent of any regime.

**Why BRIDGE_SWEEP §3's skip was right then and is superseded now.** Its
three objections, examined:

* *"alt states exceed the window cap (deg e ≤ 10)"* — that is the **sub2**
  cap; the sweep's own targets were sub2 states. Alt states are **sub1**
  states (`deg e ≤ 15`, attained). No alt state violates the sub1 window.
* *"the alt regime reduces at t^210, a different polynomial identity"* —
  that concerns the **cascade bookkeeping** (`F = t^210·G'`, the descending
  recursion). The bridge imports no cascade level; its Phi is the window
  variable, the same `c·t^30·q` in both regimes (point 3).
* *"the old caps / t-coupling do not transfer"* — `ALT_REGIME.md`'s row is
  about the cascade's **per-level `g_l` caps** and t-adic edge, which the
  bridge never uses. The **window** caps transfer because their derivation
  (point 4) never touched the t-side. (`WINDOW_CAPS.md` postdates
  `BRIDGE_SWEEP.md` — at sweep time the k=6,7,8 caps were themselves a
  [judgment] whose regime-scope was unknown; the skip was the honest call
  with the information then on the table.)

## 2. The augmented system for an alt state

State ansatz = the **audited alt witnesses** (`alt_combined.json`, defect-0;
reconstruction as in `d2_threshold.py`, audited by `d2_threshold_verify.py`):

```
a11_b3100_T2 : e = E(y+1)^11 (y−r1)^3 (y−r2),  sigma = S(y+1)^3 (y−r1)^7 (y−r2)^2
               field Q[r1,r2]/(q,q), r1 ≠ r2
a12_b1110_T2 : e = E(y+1)^12·comp,             sigma = S(y+1)^6·comp²
               comp = q/(2048(y−r)),           field Q[r]/(q)
d1 = 0 (T2),  d2 = generic deg-6 (7 free coeffs),  d0 = (d2²+sigma)/4
```

System = every y-coefficient of `G1,G2,G3,G5body+Phi` on this ansatz with
`dm2,dm3,dm4` generic stripped sub1 polynomials (66 scalars), coefficients
reduced mod `q(r_i)`, `q(r_i) = 0` adjoined, saturation by
`E·S·lc(d2)·(r1−r2)` (Rabinowitsch). Runner is orphan-proof (WSL-side
`timeout`). Verdict discipline: mod-p triage (3 primes, 45 s) → exact `Q`
(300 s); `UNIT` = kill, `PROPER` = loud survival signal, else `COST`.

## 3. Pilot design

* **Control: `a12_b1110_T2` deg d2=6** — already killed exact char-0 by
  msolve (`msolve_bridge_results.json`; that closure made the branch 8/8).
  A sound new mechanism must reproduce this known kill.
* **Prize: `a11_b3100_T2` deg d2=6** — the last open state of its 7/8
  branch (`FRONTIER_V2.md` "one state from closure"); tie-tower route: ~12
  GB RAM blowup; msolve route: TIMEOUT. A kill closes the branch.

## 4. Verdicts — BOTH PILOTS COST; NO KILL; CONTROL NOT REPRODUCED

Machine record: `alt_bridge_results.json`. The construction builds and both
pilot systems were driven through **four independent engine routes**; every
route timed out. **No `UNIT` (kill) and no `PROPER` (survival) was obtained on
either state.** Verdict for both: **COST**.

| state (deg d2=6) | system | build | Singular mod-p triage (3p, 45 s) | Singular exact `Q` (300 s) | numroot mod-p (3p, 90 s) | msolve char-0 (300 s) | verdict |
|---|---|---|---|---|---|---|---|
| **control** `a12_b1110_T2` | 182 eq / 76 vars (66 spare) | 151 s (fast) | TIMEOUT ×3 (49/45/45 s) | TIMEOUT (300.5 s) | TIMEOUT ×3 (90/90/91 s) | **TIMEOUT (306.3 s)** | **COST** |
| **prize** `a11_b3100_T2` | 183 eq / 77 vars (66 spare) | 338 s (fast) | TIMEOUT ×3 (~46 s) | (not reached; run cut) | TIMEOUT ×3 (91/90/91 s) | **TIMEOUT (300.4 s)** | **COST** |

The first run (`alt_bridge_run.log` → `alt_bridge_results_symbolic.json`, a12
only, a11 build+triage cut off) took the symbolic-`r` Singular route; the retry
(`alt_bridge_retry.log` → `alt_bridge_results.json`) took the numeric-root
mod-p + char-0 msolve route (the two engines that had beaten marked-root swell
elsewhere). Both routes, both states: TIMEOUT. The last recorded event before
Lane N died (API outage) was the `a11_b3100_T2` msolve TIMEOUT at 300.4 s,
written to `alt_bridge_results.json`; nothing was left genuinely mid-compute.

**The load-bearing negative finding.** The **control** `a12_b1110_T2` deg d2=6
is a *known* kill — msolve closed it exact char-0 on the direct/f31 route
(`msolve_bridge_results.json`; that closure made the branch 8/8). The bridge
system for the same state **times out on msolve** (and on every other engine).
**The bridge did NOT reproduce the control kill.** So the pilot does *not*
computationally validate the mechanism on the alt regime: the construction is
sound by derivation (§1, spot-checks pass), but its systems are strictly
*harder* to solve than the direct route that already worked. This is the
opposite of the standard-window pilot (`FULL_SYSTEM_BRIDGE.md` §5), where the
bridge collapsed an f31-survivor in ≈1–9 s. The +66 sub1 spare unknowns and the
marked-root (`q(r)=0`, up to two adjoined roots) nonlinearity make the deg d2=6
alt systems land squarely in the same swell family that already walled the
tie-tower (~12 GB) and direct-msolve routes — see `BRIDGE_SWEEP.md` §2/§5, the
"deg hitting the window cap" cost law (here `deg e = 15`, the sub1 cap, attained
with equality on both states).

**Recorded COST (msolve).** `a11_b3100_T2` deg d2=6, alt-bridge sub1 system
(183 eq / 77 vars / 66 spare, two marked roots `r1≠r2`): msolve exact char-0
**TIMEOUT at 300.4 s** (the emitted `.ms` was written to WSL `$HOME`; the run
also logged a benign WSL relay `chdir` warning — cosmetic, not the cause).

**Census: 0 KILL, 0 PROPER/survival, 2 COST.** No kill exists, so the
"all kills PENDING AUDIT" discipline (§5 [J4]) is currently *vacuous* — there is
nothing to hand a spec-only auditor. If a swell-tolerant follow-up ever produces
a `UNIT`, it inherits the PENDING-AUDIT status.

## 4a. Named next step

The cost wall is **engine-independent** (Singular std symbolic, Singular
numeric-root mod-p, msolve F4 char-0 all fail at the same budgets) — exactly the
signature `BRIDGE_SWEEP.md` §6 flagged for the R9/a10 deg_e=10 swell. The named
next step is therefore **not a bigger budget** but a *structurally reduced
system*: eliminate the 66 spare window unknowns `dm2,dm3,dm4` symbolically first
(they enter `G1,G2,G3` linearly in several blocks), reducing the bridge to a
much smaller system in the state coefficients + marked roots before handing it
to a solver. Secondary options: (i) specialize both marked roots to a common
numeric `q`-root pair mod a larger prime to strip the `r`-nonlinearity entirely
for a mod-p *reconnaissance* verdict (evidence only, not a `Q` certificate);
(ii) keep the deg d2=6 alt stragglers on the existing tie-tower / d2-threshold
lanes, which remain their status-quo owners. Until one of these lands, the two
deg d2=6 alt stragglers stay OPEN and the alt frontier is unchanged by this
lane.

## 5. [judgment] list

- **[J1] This doc overturns `BRIDGE_SWEEP.md` §3's applicability call** —
  on derivation, not computation. The overturn leans on `WINDOW_CAPS.md`
  (the caps' premises are t-regime-free), which postdates the sweep. The
  five soundness points are each cited to a landed artifact; the
  mechanically checkable ones are asserted in `alt_bridge.check_soundness`.
- **[J2] The G-system's subcase-independence** is inherited from the landed
  position that `augment(…, regime=…)` changes caps only (the generators are
  loaded once from `generators.json` for both regimes) — same trust basis as
  the sub2 bridge, `T6_SELECTION_AUDIT.md`.
- **[J3] State ansätze are the audited alt witnesses** (defect-0 sigma/e,
  free deg-6 d2, d1=0) — the same reconstruction `d2_threshold_verify.py`
  audits; no new valuation claims are imposed (notably none on `d2`,
  per `D2_THRESHOLD.md` [J1]).
- **[J4] All kills PENDING AUDIT** until a spec-only auditor consumes the
  recorded systems (certificate route: `kill_certificate_tools.py`).
