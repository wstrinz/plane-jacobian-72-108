# FACE_KILL_SWEEP — Milestone 1.5: the bigraded face detector, pointed at the frontier

> **STATUS (2026-07-24): 133 NEW STATE-LEVEL KILLS, BRANCH-JOIN DISCHARGED,
> PENDING INDEPENDENT AUDIT.** The criterion is derived, independently checked
> (`phi_depth_criterion_verify.py`, 36/36), linked to the exact symbolic
> computation on sampled states, and the branch-join gate is now closed (§5.1):
> the phase-D universes are already the post-cascade residual, so all 133 lie in
> still-open branches. They remain **same-author `claimed`** until audited.

Milestone 1 (`bigrade_annotator.py`) showed on three curated systems that an
extremal face of the (u-weight, y-order) lattice can expose a contradiction the
scalar Gröbner projection dissolves — R2 localised a home-case kill to a single
face equation with ZERO spare unknowns, in seconds, on 45 vars / 122 eqs. This
milestone points that detector at the frontier instead of at three hand-picked
systems.

Files (NEW): `face_kill_sweep.py` (per-state symbolic sweep),
`phi_depth_criterion.py` (closed-form O(1) specialisation),
`phi_depth_criterion_verify.py` (independent checker), `face_kill_sweep.json`,
`phi_depth_criterion.json`.

---

## 1. The mechanism

`G5 = G5body + Φ` is weighted-homogeneous of u-weight 17. Since u-weights add
under multiplication and the window floor of a weight-`w` symbol is `y^{12w}`,
**every** term of `G5` strips by exactly `y^{204}`:

| term | u-weight | floor |
|---|---:|---:|
| `d0·dm1·dm4` | 4+5+8 = 17 | 48+60+96 = 204 |
| `d0·dm2·dm3` | 4+6+7 = 17 | 48+72+84 = 204 |
| `d1·dm2·dm4` | 3+6+8 = 17 | 36+72+96 = 204 |
| `d1·dm3²` | 3+7+7 = 17 | 36+84+84 = 204 |
| `d2·dm3·dm4` | 2+7+8 = 17 | 24+84+96 = 204 |
| `dm1²·dm3` | 5+5+7 = 17 | 60+60+84 = 204 |
| `dm1·dm2²` | 5+6+6 = 17 | 60+72+72 = 204 |
| `Φ` | 17 | 204 |

> **Justification note (2026-07-24, from `monomial_lens.py`).** The seven
> `G5body` monomials **alone do not pin these weights**: their exponent lattice
> leaves a 1-dimensional null grading
> `{d0:-5, d1:-8, d2:-11, dm1:-2, dm2:1, dm3:4, dm4:7}`, so weight 17 is not
> forced by `G5` in isolation. Adding `G1,G2,G3` cuts the grading space to
> dimension 1, and its primitive generator is exactly the documented assignment
> — `(d0,d1,d2,dm1,dm2,dm3,dm4) = (4,3,2,5,6,7,8)`, `Phi = 17`, generator
> weights `13/14/15/17`. **The numbers below are confirmed; only their
> justification needed the wider generator set.**

so comparing them on one y-degree axis is legitimate — this is exactly the
`_symbol_homogeneity_and_cone` invariant, and it is re-derived from the symbol
weights in the verifier rather than assumed.

The stripped `Φ̃ = c·t^30·q`, `c = −1/6630`, has y-degree **exactly 34** with
`lc(Φ̃) = −1024/3315 ≠ 0` — and 34 is precisely the sub2 stripped cap `2·17`, so
Φ sits at its cap (the C6 tightness statement, independently recomputed).

**Therefore:** if for a given state every `G5body` term has stripped degree
`< 34` even with all three spares at their maximal admissible degree, then the
degree-34 coefficient of `G5` is `lc(Φ̃)`, the equation reads `−1024/3315 = 0`,
and **no admissible spares exist: the state is killed.**

This is the (50,75) window-depth kill mechanism (GGV3 §5) firing on the home
case: the state's data cannot reach the Φ window depth.

## 2. The closed-form criterion

The above is a pure DEGREE argument, so it needs no symbolic expansion:

> **KILL ⟺ max deg(G5body) < 34**, spares at caps

with `deg d0 = max(2·deg d2, deg σ)` and spare caps from
`full_system_bridge.STRIP_DEGCAP` (sub2 12/14/16, sub1 18/21/24 — inherited, not
recomputed, so this module cannot drift from the audited envelope).

That is O(1) per state. Full sweep of both universes:

| window | states | Φ-depth kills |
|---|---:|---:|
| sub2 | 7,888 | **195** |
| sub1 | 44,117 | **0** |

Per-`a_t` (sub2): `a_t=5` 0/332, `a_t=6` 7/1570, `a_t=7` 52/1843,
`a_t=8` 94/1892, `a_t=9` 42/1423, `a_t=10` 0/828.

**The criterion is regime-limited, and the limits are structural, not empirical:**

* `a_t = 10` is immune — `deg e ≥ 10` forces `dm1·dm2² ≥ 10+24 = 34`.
* **sub1 is entirely immune** — with caps 18/21/24, `dm1·dm2²` alone gives
  `deg e + 36 ≥ 36 > 34` for every state. No sub1 state is ever Φ-depth-killable.
* `a_t = 5` gets none because those states carry σ at its cap, lifting `d0`.

So this is coverage of the shallow end of sub2. It says nothing about the R9/alt
wall and cannot.

## 3. Verification

`phi_depth_criterion_verify.py` — **36/36 PASSED**, re-deriving from primitives
rather than importing the modules under test:

1. strip consistency (all 7 terms at u-weight 17 / floor 204);
2. `deg Φ̃ = 34 = 2·17`, `lc ≠ 0`, agrees with `_phi_stripped`;
3. caps equal `(14−12)k` / `(15−12)k` at `k = 6,7,8`;
4. the canonical `G5` carries Φ-coefficient **1** (see §4);
5. **load-bearing:** on 5 sampled KILL states, the exact symbolic degree-34
   coefficient of the canonical `G5` is *exactly* `−1024/3315`;
6. on 4 sampled non-kill states it carries free symbols (so the criterion is not
   merely failing to fire);
7. per-window roll-ups reconstruct `state_total` (7888 / 44117).

Independently, `face_kill_sweep.py` ran the full symbolic path on manifest
states and agrees, including reproducing the certified kill
`harvest:a8_dd2-inf_dd10_dsig5` and returning NO-FACE-CERTIFICATE on
`harvest:sub2T2_a7_b3000_dd2{0,1}_dsig7` — i.e. it discriminates.

## 4. A generator-normalisation bug found en route

The canonical generator is **`G5 = G5body + Φ`**:

* `full_system_bridge.py:107` — `"G5": st["G5body"] + PHI`
* `f37_sat_verify.py` — the load-bearing C11 membership certificate verifies
  `f31 == c1·G1 + c2·G2 + c3·G3 + c4·(G5body + Φ)`

But `bigrade_annotator.py:675` uses `G5 = 2·Φ + G5body`, transcribed from
**`FULL_SYSTEM_BRIDGE.md:62`**, which contradicts **line 114 of the same file**
(`(G5body+Phi)`). Line 62 is the erroneous one.

These differ by `Φ`, not by a nonzero scalar, so they are genuinely different
equations and conclusions do **not** transfer automatically. Impact:

* the window-depth kill **survives** (only `deg Φ = 34` and `lc(Φ) ≠ 0` are used);
* but the emitted certificate VALUE was wrong: canonical `−1024/3315`, not
  `−2048/3315`;
* **M1's R2 certificate value is wrong**, and R2's rank/consistency numbers were
  computed on a non-canonical `G5`. R3 is unaffected (it uses `_H_generators()`
  from `r9_eliminated_system.json`).

Both new modules use the canonical normalisation. `FULL_SYSTEM_BRIDGE.md:62` and
`bigrade_annotator.py:675` still need repair.

Found by an independent adversarial review that was specifically instructed to
attack the strip arithmetic; the strip arithmetic held, this did not.

## 5. What is NOT established (the honest gate)

### 5.1 Branch-join — DISCHARGED

Of the 195, **62 are already in `state_kill_ledger.json` and 133 are not**. The
concern was that the ledger is *state-level* while the cascade kills branches, so
some of the 133 might lie in branches already dead. **They do not.** Joining the
phase-D universes against the branch-audit artifacts:

| universe | cases | matched to a cascade branch | branch verdict |
|---|---:|---:|---|
| `phase_d_states_sub2.json` | 220 | 220 | **`survives` — all 220** |
| `phase_d_states_sub1.json` | 1145 | 1145 | **`survives` — all 1145** |

Zero unmatched, zero in a killed branch. This is by construction, not luck: the
sub1 universe self-describes as *"Complete residual degree states per surviving
flag case (sub1, q+t+inf, residue kills, T2 squeeze)"*
(`source_artifact: cascade_cones_sub1_qt_inf_rl.json`). The phase-D universes
ARE the open frontier at state level, already net of the q-cascade, the t-place
layer, the infinity layer, the residue kills and the T2 squeeze.

**So the 133 are genuine new state-level kills in still-open branches.**

### 5.2 Remaining gaps
* **[claimed, not audited]** The 133 are same-author at present. They are,
  however, produced by a mechanism with no shared code path with the exact-GB /
  msolve route — which is what makes them a candidate *independent audit* of the
  62 overlapping states (see §6).
* **[scope]** sub1 and the alt layer yield nothing here, structurally (§2). No
  claim is made about the R9/alt wall.
* **[modelling]** Spares are treated as arbitrary polynomials up to cap. This is
  the generous direction for a kill (more freedom ⇒ harder to kill), so it is the
  safe side, but it does assume no further structural constraint *lowers* the
  reachable degree in a way that would kill even more states.

## 5.3 Per-state symbolic sweep — partial run, stopped deliberately

`face_kill_sweep.py` (the full symbolic path, as distinct from the O(1)
criterion) was run over `kill_manifest.json` and **stopped at 26 of 49 entries**
after 4 hours. Census of what completed: **10 KILL, 7 NO-FACE-CERTIFICATE,
9 SKIPPED**.

This is a deliberate stop, not a crash, and the discarded remainder is
*informationally* cheap:

* every **sub2** entry completed — those are the ones the criterion can decide;
* every remaining entry is alt-layer or phase, i.e. **sub1 caps**, and sub1 is
  structurally immune to the Phi-depth kill (§2: `dm1*dm2^2 >= deg e + 36 > 34`
  for every sub1 state);
* the completed sub1-cap entries all returned NO-FACE-CERTIFICATE, consistent
  with that;
* the remaining entries average ~45 min each (66 spare unknowns vs 45), i.e.
  ~17 further hours.

The 7 NO-FACE-CERTIFICATE results are the useful negative datum: the detector
**discriminates** rather than firing on everything — notably
`harvest:sub2T2_a7_b3000_dd2{0,1}_dsig7`, where the larger sigma lifts `deg d0`
enough for `G5body` to reach degree 34.

## 6. The secondary product: an evidence-grade upgrade path

The 62 overlapping states are killed in the ledger by "triage/bridge exact GB
(UNIT ideal)" with **audit = PENDING**. The Φ-depth criterion re-derives them by
a *different mechanism* (window-depth degree argument, not a Gröbner unit-ideal
computation) with no shared code path — the repo's own definition of an
independent re-derivation. Wiring this as a ledger attribution source would
upgrade those states `PENDING → independently-audited`.

Given that the program's stated moat is audited state rather than throughput,
this may be the more valuable output of the two.

## 6.1 COMPLETENESS: the face detector's reach is exactly the Phi-depth criterion

The five bigrade functionals cannot find anything the Phi-depth criterion misses,
and the `--escalate` parameter-form branch is **unreachable in this
construction**. This is structural, not empirical:

* **Every one of the seven `G5body` terms contains at least one spare**
  (`dm2`, `dm3` or `dm4`). So the degree-34 coefficient of `G5` is either
  `lc(Phi)` ALONE (no `G5body` term reaches 34 -> a rational constant -> sound
  KILL) or it carries spare symbols (-> not spare-free -> nothing to classify).
  **Binary. There is no parameter-form case at the top face.**
* **`G1`, `G2`, `G3` can never present a spare-free face equation**, because each
  has a PURE-SPARE product sitting exactly at its cap:
  `dm2*dm3 = 12+14 = 26` (G1 cap), `dm2*dm4 = 12+16 = 28` and `dm3^2 = 28`
  (G2 cap), `dm3*dm4 = 14+16 = 30` (G3 cap). Both factors are at their own caps,
  so the product's leading coefficient (`R12*S14` etc.) is always present as a
  symbol.

Empirical confirmation over six sub2 states spanning the degree range
(`escalate_probe.log`):

| state | unknown-free face eqs | values |
|---|---:|---|
| small (Phi-depth KILL) | 2 | both CONSTANT `-1024/3315` |
| small-mid | 2 | both CONSTANT `-1024/3315` |
| mid | 2 | both CONSTANT `-1024/3315` |
| near-cap | **0** | (none) |
| full cap (all maximal) | **0** | (none) |
| high-degree | **0** | (none) |

**Consequence — this lane is CLOSED.** The face-detector milestone yields the
Phi-depth criterion and nothing else; there is no further kill hiding in the
other four functionals, and no reason to keep mining it. The `--escalate`
machinery below is correct and would matter for a construction with a
forced-inhomogeneous term in more than one generator, but for the (72,108)
G-system it is vacuous — which retroactively makes "off by default" the right
call for a stronger reason than caution.

## 6.2 CORRECTION (2026-07-24): the stratum-union objection was wrong

This file previously recorded that collapsing a flag-case's exact-degree strata
into one system is unpromising, on the grounds that "the relaxed system is less
constrained, hence a priori more likely consistent". **That reasoning applies to
the wrong object and is withdrawn.**

The union of the degree strata is NOT the relaxed affine variety `V(I)`. For a
residual polynomial written at its cap, `R = r_0 + r_1 y + ... + r_D y^D`, the
strata are `{r_D != 0}`, `{r_D = 0, r_{D-1} != 0}`, ... and their union is
exactly `{(r_0,...,r_D) != 0}` -- i.e. `R` not identically zero. That is a
**quasi-affine saturated locus** `V(I) \ V(B)`, tested by the saturation
`I : B^inf`, where `B` is the irrelevant / endpoint ideal. It preserves precisely
the open conditions that the per-state Rabinowitsch variables currently enforce;
it does not discard them.

So the correct formulation is a SATURATED flag-case scheme, and the objection
recorded here does not apply to it. The separate finding that the sub2 system has
a single 45-variable SCC still stands, but it bears on *triangularisation*, not on
whether the saturated union is decidable.

Two implementation caveats, worth stating before anyone builds it:
* it must be constructed BEFORE any state-specific division by a leading
  coefficient -- such a division is a chart condition whose zero branch must
  remain present;
* the right pilot is a small MIXED cell containing both killed and unresolved
  strata, not the largest cell.

## 7. Optional (c)-escalation — off by default

`face_kill_sweep.py --escalate` promotes a parameter-form CONSTRAINT to a KILL
when `val = 0` is inconsistent with the state's parameter ideal + saturation,
i.e. when `⟨val, q(r), w·∏sat − 1⟩ = ⟨1⟩`. It runs a Gröbner call in a worker
process under a variable-count guard (≤ 8) and a hard 60 s timeout; either guard
tripping downgrades to CONSTRAINT, never to KILL.

It is **off by default and never on the CI path**, because it reintroduces
exactly the cost the face method avoids and is only as sound as the modelled
constraint set. Escalated kills are `claimed`, not certified.
