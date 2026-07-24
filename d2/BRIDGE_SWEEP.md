# BRIDGE_SWEEP — the full-system bridge across the remaining resistant core (2026-07-23)

**Status: candidate EXACT KILLS, PENDING AUDIT.** Machine record:
`bridge_sweep.json`. Verifier: `bridge_sweep_verify.py`, **ALL 4 CHECKS
PASSED**.

**Headline census — 35 exact-Q kills, 0 PROPER (no two-gate survivor
anywhere), 8 COST:**

| target group | states | EXACT KILL (UNIT/Q) | PROPER | COST |
|---|---:|---:|---:|---:|
| a8 constant-E class (Target 2, full) | 24 | **24** | 0 | 0 |
| a8 T2 deg_e=8 batch class (Target 4, full class) | 5 | **5** | 0 | 0 |
| a7 T2 deg_e=8 batch class (Target 4, full class) | 5 | **5** | 0 | 0 |
| R9 column z=0..6 (Target 1) | 7 | **1** (z=0) | 0 | 6 |
| a10 T1 constant-E deg_e=10 (Target 4, sampled) | 2 | 0 | 0 | 2 |
| alt-regime targets (Targets 3, 4-alt) | — | skipped: bridge not applicable (§3) | | |

**NOT ONE state anywhere returned PROPER under the bridge — mod-p or exact.
The first true two-gate survivor does not exist in anything this sweep could
afford to compute.** Every non-kill is a pure Groebner-cost verdict with an
explicit reason recorded, never a solvability signal.

**Mission.** `FULL_SYSTEM_BRIDGE.md` validated the endgame tool: augmenting a
resistant state with the full pre-resultant G-system (the window unknowns
`dm2,dm3,dm4` as bounded stripped polynomials) turns f31-resistant states into
~122-equation quadratic systems that go UNIT over exact `Q` in seconds. This
sweep runs that bridge across the entire remaining resistant core.

**New files (uncommitted):** `bridge_sweep.py` (the sweep + the marked-root
extension), `bridge_sweep.json` (machine census), `bridge_sweep_verify.py`
(independent re-derivation of two kills), this doc. READ-ONLY on every landed
module/artifact; nothing committed.

---

## 0. The marked-root extension (the one construction this sweep adds)

`full_system_bridge.augment` takes any `convolution_descent.Ansatz` (stripped
`d2,d1,sigma,e`) and returns the G-system equations with the spare window
unknowns. Its pilot handled only r-free states (gamma-saturation). The R9 /
T2-pattern-B states carry a **marked root** `r` of the fixed quartic
`q = 2048y^4-512y^3+320y^2-240y+195`. The extension (documented in
`bridge_sweep.py`, mirroring the pilot's stripping discipline exactly):

1. The state ansatz is the **landed** R9 construction
   `convolution_elim_qsupport.build_qsupport_ansatz(z)`:
   `e = gamma(y+1)^9(y-r)`, `sigma = (y-r)^2 G` (`deg G = z`), `d1 = 0`,
   `deg d2 <= 4`, `d0 = (d2^2+sigma)/4` — already stripped and already at the
   sub2 window caps (`deg e = 10 = 2*5`, `deg sigma <= 8 = 2*4`,
   `deg d2 = 4 = 2*2`), so `augment(..., regime="sub2")` applies verbatim; `r`
   rides along symbolically inside the coefficients.
2. Adjoin `q(r) = 0` (the exact `Q[r]/(q)` representation; equivalently, the
   exact-Q attempts also use Singular's native `minpoly` field `Q(r)` —
   mathematically identical since `q` is irreducible over `Q`).
3. Saturate by the state's genuine nonzero scalars `(gamma, lc(G), G(r))` via a
   single Rabinowitsch `w*prod-1` — verbatim the landed
   `build_qsupport_ansatz.saturation_factors` (G(r)!=0 makes `v_r(sigma) = 2`
   exact; lc(G)!=0 makes `deg G` exact; gamma!=0 the gauge).

Soundness is inherited unchanged from `FULL_SYSTEM_BRIDGE.md` §4 (the G-system
is the genuine window system; `f31` lies in the G-ideal, so the bridge is `>=`
the cascade; stripping is exact by weighted homogeneity; the `k=6,7,8` caps are
the flagged [judgment] extension of `T3_WINDOW_AUDIT`).

---

## 1. Target 2 — the a8 constant-E class: **24/24 EXACT KILLS over Q**

The entire a8 constant-E UNRESOLVED class of `batch_convolution_sub2.json`
(`a_t=8, T1, deg_e=8`; `deg_sigma in 5..8` x `d2 in {0, deg 0..4}`), swept with
the core bridge (`fsb.augment`, regime sub2, gamma-saturated):

| deg_sigma | states | mod-p (3 primes) | exact Q | wall (exact) |
|---|---:|---|---|---|
| 5 | 6 | UNIT all | **UNIT — KILLED** | 0.8–24 s |
| 6 | 6 | UNIT all | **UNIT — KILLED** | 0.6–6 s |
| 7 | 6 | UNIT all | **UNIT — KILLED** | 0.6–15 s |
| 8 | 6 | UNIT all | **UNIT — KILLED** | 0.9–15 s |

**Every state: mod-p UNIT on all three triage primes AND an exact `Q` UNIT
certificate** (Rabinowitsch gamma-saturation; integer-cleared emission). This
closes the two open a8 sub-questions at once:

* `MODULAR_TRIAGE.md` System 2 scored all 24 LIKELY-SOLVABLE and argued only a
  real/sign argument could kill them; `TRIAGE_HARVEST.md` reversed that mod-p
  for `deg_sigma` 5/6 but left **all 12 `deg_sigma` 7/8 states cost-bound with
  no exact-Q certificate**. The bridge kills all 24 **exactly over Q in
  seconds** — including the pilot's two states, now re-confirmed in-sweep.
* No real/sign mathematics is needed anywhere in the a8 class: the full window
  system is empty over the algebraic closure, per an exact rational
  certificate. The a8 chapter is closed (pending audit).

## 2. Target 1 — the R9 column z=0..6

**z=0: KILLED — exact `Q` UNIT in 3.4–27.8 s** (mod-p UNIT on all three triage
primes first, 17–26 s each). This is the **first exact-`Q` kill anywhere in the
R9 column**: the master-coefficient route was cost-bound over `Q` even at z=3
with four different attempts (TRIAGE_HARVEST.md Target 1). The bridge does
bypass the f31 swell at z=0.

**z=1..6: exact COST — and an honest structural finding.** For z >= 1 the
bridge system itself becomes the heavy object. At z=1, five formulations all
time out:

| formulation | engine | budget | result |
|---|---|---:|---|
| symbolic r, q(r) adjoined, char 0 | Singular std | 300 s | TIMEOUT |
| symbolic r, q(r) adjoined, char 0 | msolve | 600 s | TIMEOUT |
| number field Q(r) via `minpoly`, char 0 | Singular | 300 s | TIMEOUT |
| hybrid: + 7 f31 master coeffs, char 0 | Singular / msolve | 300 s each | TIMEOUT |
| numeric root of q mod p (reconnaissance) | Singular, 3 primes | 90 s/prime | TIMEOUT (z=1,2) |

The z=1..6 mod-p bridge verdicts are INDETERMINATE at this budget; z=0 numroot
mod-p is UNIT (18/9/13 s) consistent with the exact kill.

**Two honest conclusions, reported with care:**

1. **The R9 z>=1 timeouts are pure Groebner cost, not evidence of survival.**
   Soundness argument: `f31` lies in the G-ideal (`F37_SATURATION_REPORT` fact,
   re-verified by `full_system_bridge_verify.py` V2), so on the bridge variety
   every f31 master coefficient of the state vanishes. `MODULAR_TRIAGE.md`
   System 1 proved the 8-master-coefficient system UNIT mod p (all three
   primes) for z <= 3 under the identical `(gamma, lc G, G(r))` saturation.
   Hence the bridge system for z <= 3 is **logically empty over `F̄_p`** — its
   own GB simply doesn't terminate at these budgets. For z=4..6 both routes are
   cost-bound (as they were for the triage) and the column trend stands.
2. **The bridge's advantage is class-dependent — it does NOT bypass the R9
   marked-root swell for z >= 1.** The asymmetry inverts between classes: on
   constant-E / generic-tail states (a8, batch T2/T1) the bridge kills in
   seconds where deep-f31 was infeasible; on the marked-root R9 column the
   8-coefficient f31 system is the fast one mod p (seconds, z <= 3) while the
   122-equation bridge with 45 spare unknowns and the `G(r)`-saturation swells
   in every engine. The window constraints that collapse tail states do not
   cheaply collapse root-supported sigma. R9 z=1..6 exact-`Q` certificates
   remain the program's outstanding cost frontier (as before this sweep;
   nothing regressed — z=0 is strictly new).

## 3. Target 3 — the two deg-d2=6 alt stragglers: NOT APPLICABLE (honest skip)

`D2_THRESHOLD.md`'s two `deg d2=6` survivors (`a11_b3100_T2`, `a12_b1110_T2`)
live in the **alternate regime** (`a in 11..15`, `v = 30-3a < 0`,
`ALT_REGIME.md`). The bridge's G-system is a **standard-regime window object**:

* Its `Phi = c t^30 q` is the standard reduction's `Phi_full/y^204` at weight
  204; the alt regime reduces at `t^210` with per-f exponents `t^{(7-f)w}u^f` —
  a *different* polynomial identity (`F = t^210 G'`, not `F = t^{21a} G`).
* `ALT_REGIME.md`'s survival table states verbatim that the old lower-level
  global caps / t-coupling "**do not transfer**; negative `v` reverses the
  degree recursion and t-adic edge." The G-system generators
  `G1,G2,G3,G5body+Phi = (D~^3)_{-1,-2,-3,-5}` are exactly such standard-window
  levels, and the spare-unknown caps (`ord >= 12k`, `deg <= 14k/15k`) are
  T3-window facts for the standard sub1/sub2 windows.
* The alt states' `e = E(y+1)^{11..12}...` with `deg e = 14..15` **exceed the
  sub2 window cap** (`deg e <= 10` stripped) — the sub2/sub1 window bounds the
  bridge relies on simply do not hold for these divisors.

Applying the bridge there would impose standard-window necessary conditions on
states that are constructed to violate the standard window — unsound. **Skipped
honestly.** The same applies to the two cost-bound System-3 alt tie-tower
states (`a11_b1110_T1 sup1`, `a11_b3000_T1 sup9`) and the alt entries of the
batch census (Target 4): every alt-regime UNRESOLVED class is outside the
bridge's window. The alt stragglers remain with the (running, separate)
blowup-diagnosis lane and the d2-threshold route.

## 4. Verification — `bridge_sweep_verify.py`: **ALL 4 CHECKS PASSED**

Two kills re-derived with **fresh constructions** (generators loaded straight
from `t4_state.pkl`, spare unknowns re-declared with new names, Phi re-derived,
states rebuilt by hand — no `fsb.augment`, no `bridge_sweep` reuse):

| kill | LEG-1 | LEG-2 |
|---|---|---|
| **A: R9 z=0** (marked root) | Singular, fresh prime 32003: UNIT | numeric-root construction (`r -> root of q mod 65027`, no symbolic `r`, no `q(r)` generator), Singular fresh prime 65027: UNIT |
| **B: a8 dsig8 dd13 d2=0** | Singular, fresh prime 32003: UNIT | **msolve (different CAS), exact char 0: `[-1]` no solution over Qbar** |

Honest note on the different-CAS clause: for KILL A, msolve swells on the
marked-root G-system in *every* characteristic (char 0 > 400 s; even mod 65537
its F4 exceeded 400 s where Singular's std takes ~20 s), so KILL A's second leg
varies the construction and prime instead — flagged in the verifier source.
KILL B carries the true cross-CAS **exact rational** confirmation.

## 5. Target 4 — batch-census UNRESOLVED sample (standard-window part)

The `batch_convolution_sub2.json` UNRESOLVED census holds 130 states; 24 are
Target 2 (above); the rest were sampled cheapest-first (m = deg_e − a_t, then
total degree):

| class | states swept | verdict | exact-Q wall |
|---|---:|---|---|
| `a8 T2 deg_e=8` (m=0; deg_sigma 3, d2 in 0..4) — **entire class** | 5 | **5/5 KILLED** | 0.5–29 s |
| `a7 T2 deg_e=8` (m=1; deg_sigma 3, d2 in 0..4) — **entire class** | 5 | **5/5 KILLED** | 0.5–2.3 s |
| `a10 T1 deg_e=10` constant-E (deg_sigma 7, 8) | 2 | 2 COST (mod-p 300 s timeouts ×3 primes; exact 300 s timeout) | — |
| `a7/a8/a9 T2 deg_e=10` (m=1..3) | 0 | not reached (cheapest-first order; same deg_e=10 swell family as a10/R9, see §2) | — |

Two complete UNRESOLVED classes die outright (10 exact kills). The a10
deg_e=10 states are the constant-E cousins of the R9 swell: **bridge cost
tracks deg_e hitting the sub2 window cap (deg_e = 10)**, not the marked root
per se — a clean, previously-unknown cost law for the endgame tool. The
deg_e=10 T2 classes (a7 m=3, a8 m=2, a9 m=1; 90 states) were not reached and
remain for a swell-tolerant follow-up.

## 6. [judgment] list

- **[RETIRED 2026-07-23 — now PROVEN] The k=6,7,8 window caps** for
  `dm2,dm3,dm4` were the flagged extension of `T3_WINDOW_AUDIT`
  (`FULL_SYSTEM_BRIDGE.md` §4); every kill in this sweep rests on them.
  Now proven in `WINDOW_CAPS.md` + `window_caps_verify.py` (81 checks,
  in run_tests.sh) — same trust tier as the audited k=2..5 window.
- **[judgment] R9 z=1..3 mod-p emptiness by implication.** I claim the bridge
  variety is empty over `F̄_p` for z ≤ 3 by the chain (f31 ∈ G-ideal) +
  (MODULAR_TRIAGE System-1 8-coeff UNIT, same saturation) — not by a
  terminated bridge GB. The z=1..6 **exact-Q** verdicts remain COST, honestly
  open.
- **[judgment] Alt-regime inapplicability (§3).** The bridge is a
  standard-window object; the alt reduction (`F = t^210 G'`) and the alt
  states' deg_e > 10 put Targets 3 / 4-alt outside its soundness domain. This
  is an applicability argument from `ALT_REGIME.md`'s survival table, not a
  computation.
- **[judgment] Cheapest-first truncation.** The batch sample was stopped after
  the two a10 COST states rather than grinding the remaining deg_e=10 states
  at ~20 min/state; the two complete killed classes + two COST probes are
  representative, but 90 deg_e=10 T2 states are unswept.
- **[judgment — engine note] The deg_e=10 cost wall is engine-independent**
  (Singular std, msolve F4, minpoly field arithmetic, hybrid over-determination
  all fail at z=1) — future exact certificates for R9 z≥1 / a10 likely need a
  structurally reduced system (e.g. eliminating the spare unknowns
  symbolically first), not a bigger budget.

## 7. Reproduce

```
python bridge_sweep.py r9 0          # R9 z=0: triage + exact kill (~7 min)
python bridge_sweep.py a8            # all 24 a8 states (~1.5 h)
python bridge_sweep.py batch 20      # batch census sample, cheapest-first
python bridge_sweep_verify.py        # must end: ALL 4 ... CHECKS PASSED
```
