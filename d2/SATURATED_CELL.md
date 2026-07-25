# SATURATED_CELL — one flag-case as ONE saturated scheme (pilot)

> **STATUS (2026-07-24): MECHANISM VALIDATED, PILOT CELL *NOT* CLOSED — COST.**
> The saturated flag-case scheme of `FACE_KILL_SWEEP.md` §6.2 is built, and its
> pipeline **reproduces every control**: both ledger kills of the pilot cell,
> both of them again as degree-EXACT strata, the **saturated union of the two**
> (the code path that the whole idea rests on), and a state recorded *not* empty
> comes back `PROPER`. On the target cell the exact-Q saturation closes the
> stratum union `deg d2 ≤ 2, deg sigma = 7` — **3 of the cell's 10 states in one
> computation, one of which is a new state-level kill** — and then hits a hard
> cost cliff one step further out (`deg d2 ≤ 3`: `TIMEOUT` at 900 s exact-Q,
> `UNIT` in 6.5 s mod 10007). The `deg sigma = 8` half of the cell does not
> finish **even mod p at 7 ring variables**, so the `f31` engine cannot close
> this flag-case at any width — while the *bigger* 50-variable G-system engine
> (`face_kill_sweep.build_state_system`) returns `UNIT` on both a state and a
> saturated union. **The next move is the G-system, not a bigger `f31` budget.**
>
> **The flag-case does NOT close here, and no timeout in this file is reported as
> a mathematical verdict.** What is established is that the saturated
> formulation is *correct and computable*, that it is *not free* (the union
> costs ~20× its own strata at `deg d2 ≤ 2`), and exactly where it stops. The
> mod-p rows are reconnaissance in the sense of `MODULAR_TRIAGE.md` — a
> prediction, never a kill.

Files (NEW, uncommitted): `saturated_cell.py`; `saturated_cell.json` (the
control battery, §4); and the run logs `sat_cell_{target,modp,gensig,n16,
gsys,rows}.log`, which are the primary evidence for §5/§6 and are quoted
verbatim below. READ-ONLY on every existing artifact and module.

**Deliberately stopped / still in flight at write-up** (stated so nobody reads a
missing row as a verdict): the `--target` exact-Q ladder, the `--ncoeff 16`
ladder and the `--sigma-generic` ladder were **stopped by hand** once each had
delivered its cost datum and every remaining row was provably more of the same
wall (precedent: `FACE_KILL_SWEEP.md` §5.3); the mod-10007 ladder, the `--gsys`
probe and a `--rows 1:2:2,2:2:2,1:1:2,2:1:2` run (the unions spanning *both*
`sigma` strata) were **still running** and will write `saturated_cell_modp.json`
/ `saturated_cell_gsys.json` / `saturated_cell_rows.json` when they finish.
Every number quoted in this file was observed; nothing is projected.

---

## 1. The object

A phase-D **flag case** is a cell `(a_t, b, branch, sigma_zero, d2_zero,
g_zero_levels)`; its states are the admissible degree tuples
`(deg d1, deg d2, deg sigma, deg e)`. Those states are **degree strata of one
ansatz**, not independent problems: `deg d2 = 3` means `a3 ≠ 0 AND a4 = 0`.

Write a state polynomial at the cell's maximum degree, `P = p_0 + … + p_D y^D`.
The strata `{deg P = d}`, `d ∈ [d_min, D]`, union to

```
{ (p_{d_min}, …, p_D) ≠ 0 }  =  A^{D+1} \ V(B),     B = <p_{d_min}, …, p_D>,
```

so the union of the cell's states is the **quasi-affine** locus `V(I) \ V(B)`,
whose emptiness is decided by the saturation `I : B^∞`:

```
V(I) \ V(B) = ∅   ⟺   I : B^∞ = <1>          (over an algebraically closed field)
```

With several such conditions (one per polynomial, plus the essential
`gamma ≠ 0`), the object to test is `I : (B_1 B_2 ⋯)^∞`, and
`I : (B_1B_2)^∞ = (I : B_1^∞) : B_2^∞`, so **iterated `sat` is exactly right** —
which is what `saturated_cell.saturated_program` emits (`sat(Is,B)[1]` in a
chain, `LIB "elim.lib"`).

Both implementation caveats of §6.2 are honoured:

* **(a) built before any leading-coefficient division.** Degree-exactness is
  imposed *only* as a saturation, never as a substitution or a chart. No branch
  where a leading coefficient vanishes is silently discarded.
* **(b) small MIXED pilot cell**, not the largest.

### 1.1 A free corollary: the Φ-depth criterion, read at cell level, closes nothing

The Φ-window-depth kill of `FACE_KILL_SWEEP.md` is monotone in the state
degrees, so a cell is *entirely* Φ-depth-killed **iff its top stratum is**. Run
over both universes (`phi_depth_criterion.classify_state`, O(1)/state):

| window | cells | fully Φ-depth-killed | partially |
|---|---:|---:|---:|
| sub2 | 220 | **0** | 30 |
| sub1 | 1145 | **0** | 0 |

So the cheap criterion cannot close a single flag-case, by a one-line argument.
That is precisely why the saturated scheme has to be run through a Gröbner
saturation, and it is the reason this pilot exists.

---

## 2. Target selection

Census of the 220 sub2 cells against `state_kill_ledger.json` (205 sub2 kills,
all 205 matched to a state key) and the Φ-depth criterion: **39 cells are mixed**
(some killed, some not); the smallest are 4, 10, 15, 17, 20 states.

The pilot is

> **`sub2 a7_b3000_T2`** — `a_t=7`, `b=(3,0,0,0)`, branch `T2`,
> `sigma_zero=False`, `d2_zero=False`, `g_zero_levels=()` — **10 states**.

```
states: 10   degree profile: {'deg_d2': [0, 1, 2, 3, 4], 'deg_d1': [None],
                              'deg_sigma': [7, 8], 'deg_e': [10]}
```

Chosen for four reasons, in order of weight:

1. **Its state set is an exact PRODUCT OF DEGREE INTERVALS** — `{0..4} × {7,8}` —
   so the stratum union is *one* saturation with **no over-cover**.
   `saturated_cell.cell_profile` asserts this and fails loud otherwise. (Across
   sub2, 88 of 220 cells have this shape, covering 2502 of 7888 states; sub1,
   254 of 1145 cells covering 640 of 44117. Non-product cells are still a finite
   union of boxes — just not a single saturation.)
2. **It is MIXED**: 2 of the 10 states are in `state_kill_ledger.json`
   (`triage/bridge exact GB (UNIT ideal)`, source `triage_harvest.json:system4`,
   audit `PENDING`), 8 are unresolved. Neither Φ-depth nor anything else touches
   the other 8.
3. **Those 2 kills were produced by a mechanism this pilot can rebuild
   verbatim** (`triage_harvest.build_sys4_8`), so the controls are real controls
   and not a different model dressed up as one.
4. It is small: `d1 ≡ 0` (T2), `e` forced by the cell flags (§3), leaving 5 + 3
   state coefficients.

**Why not the "known-closed cell" control the brief suggests.** There is **no
closed phase-D cell** — the headline of `SESSION_HANDOFF.md` is *0 of 220*. The
closed families it names (`a11_b1111_T1`, `a12_b1110_T2`, `a14_b0000_T2`) live in
the **alt defect-0 layer**, are `sub1`-cap / alternate regime, and their closures
come from three *different* lanes (`phase_f2_scale` defect-0 tie tower,
`D2_THRESHOLD.md` depth-8 saturated GB, exact-char-0 msolve) — none of which is
the master-coefficient saturated GB this pipeline runs. Reproducing them would
be reproducing three other lanes, not validating this one. The closure control
was therefore taken **at state level and at union level** instead (C1 / C1c
below), where it is exact and in-model. This is stated as a limitation, not as a
pass.

---

## 3. The construction

`saturated_cell.build_cell_system` (no leading-coefficient division anywhere):

```
e      = gamma · (y+1)^{a_t} · (y-r)^{Σb}          gamma ≠ 0,  q(r) = 0
d1     = 0                                          (T2 branch, cell flag)
d2     = a_0 + a_1 y + … + a_4 y^4                  endpoint ideal <a_0,…,a_4>
sigma  = (y-r)^{2Σb} · G,  G = g_0 + g_1 y + g_2 y^2 endpoint ideal <g_1,g_2>
                                                    plus G(r) ≠ 0
```

* `e` is **cell data, not state data**: `v_t(e) = a_t` and `v_{r_j}(e) = b_j` are
  the flag coordinates, and `deg e = 10 = a_t + Σb` for *every* state of this
  cell, so `e` is forced up to the scalar `gamma`. Imposing it does not narrow
  the flag-case. `gamma ≠ 0` is the "essential" saturation.
* the endpoint ideal `<a_0,…,a_4>` is exactly `deg d2 ∈ {0..4}`; `<g_1,g_2>` is
  exactly `deg G ∈ {1,2}`, i.e. `deg sigma ∈ {7,8}` — the cell, on the nose.
* `sigma = (y-r)^{2Σb} G` with `G(r) ≠ 0` is **inherited verbatim** from the lane
  that produced the two known kills. See §7 — it is the single most important
  caveat in this file.

Two engines, both necessary conditions, the second containing the first:

| engine | generators | ring |
|---|---|---|
| `f31` | top 8 nonzero y-coefficients of `f31 = Σ_f Φ^f e^{21-3f} h_f`, reduced mod `q(r)`, plus `q(r)` | 6–11 vars |
| `gsys` | every y-coefficient of `G1,G2,G3,G5 = G5body+Φ` on stripped spare ansaetze at the sub2 caps 12/14/16, via **`face_kill_sweep.build_state_system`** | 50–55 vars |

`<G-system> ∩ Q[d2,d1,d0,e,Φ] = <f31>` (`F37_SATURATION_REPORT` fact [5]), so
each y-coefficient of `f31(state)` lies in the ideal generated by the
y-coefficients of the `G_i` — the `gsys` ideal **contains** the `f31` ideal. An
`f31` kill therefore implies a `gsys` kill, and `gsys` is the stronger test.

Soundness of using only 8 master coefficients: the ideal is a **subset** of the
true one, so `UNIT` is still a kill; failure to be unit proves nothing.

---

## 4. Controls — 7 of 7 PASS

Quoted verbatim from the run (`python saturated_cell.py --controls`):

```
CONTROL C1 -- reproduce the two LEDGER KILLS of this cell
  (state_kill_ledger.json: 'triage/bridge exact GB (UNIT ideal)',
   source triage_harvest.json:system4 -> UNIT, 1.14s / 1.52s)
  C1_state_dd20_dsig7 (as-built strata)   UNIT  dim=-1   0.38s  vars=5  gens=9  [expected UNIT -> OK]
  C1_state_dd21_dsig7 (as-built strata)   UNIT  dim=-1   0.36s  vars=6  gens=9  [expected UNIT -> OK]

CONTROL C1b -- the SAME two strata, but DEGREE-EXACT (lc saturated)
  C1b_exact_dd20_dsig7                    UNIT  dim=-1   0.23s  vars=5  gens=9  [expected UNIT -> OK]
  C1b_exact_dd21_dsig7                    UNIT  dim=-1   0.33s  vars=6  gens=9  [expected UNIT -> OK]

CONTROL C1c -- the SATURATED UNION of those two known-empty strata
  (endpoint ideal <a0,a1>: exactly deg d2 in {0,1}; must be EMPTY)
  C1c_union_dd2le1_dsig7                  UNIT  dim=-1   0.34s  vars=6  gens=9  [expected UNIT -> OK]

CONTROL C2 -- a state recorded NOT empty must NOT come out UNIT
  (modular_triage.json system2 a8_dd20_dd10_dsig5: PROPER dim 2, 3 primes)
  a8_dd20_dd10_dsig5@p10007               PROPER dim=2   0.30s  vars=6  gens=4  [expected PROPER -> OK]
  a8_dd20_dd10_dsig5@p10009               PROPER dim=2   0.22s  vars=6  gens=4  [expected PROPER -> OK]
```

* **C1** reproduces the two ledger kills of the pilot cell exactly, by rebuilding
  their generator set — the mechanism check.
* **C1b** re-runs them as *degree-exact* strata (leading coefficient saturated,
  not divided out) — caveat (a) exercised.
* **C1c** is the load-bearing one: the **saturated union** `I : <a_0,a_1>^∞` of
  two strata known to be empty comes out empty. If the saturation code path were
  wrong, this is where it would show.
* **C2** takes a system the repo records as `PROPER` (dimension 2, three primes)
  and gets `PROPER`, dimension 2 — the pipeline is not manufacturing unit
  ideals.

---

## 5. The target ladder

`python saturated_cell.py --target` walks the stratum union outward. Exact Q
(`char 0`; a `UNIT` here is a kill) and mod 10007 (reconnaissance only):

| stratum union | strata covered | exact Q | mod 10007 |
|---|---|---|---|
| `deg d2 ≤ 1, deg sigma = 7` | 2 | `UNIT` 0.34 s | `UNIT` 0.31 s |
| `deg d2 ≤ 2, deg sigma = 7` | 3 | **`UNIT` 8.83 s** | `UNIT` 0.27 s |
| `deg d2 ≤ 3, deg sigma = 7` | 4 | **`TIMEOUT` 900 s** | `UNIT` 6.52 s |
| `deg d2 ≤ 4, deg sigma = 7` | 5 | not reached (§6, wall 2) | **`TIMEOUT` 600 s** |
| `deg d2 ≤ 1, deg sigma = 8` | 2 | not reached | **`TIMEOUT` 600 s** (7 vars!) |

```
TARGET -- saturated stratum unions in cell a7_b3000_T2  [satGr]
  T_dd2le1_dsig7_satGr    UNIT     dim=-1     0.34s  vars=6  gens=9
  T_dd2le2_dsig7_satGr    UNIT     dim=-1     8.83s  vars=7  gens=9
  T_dd2le3_dsig7_satGr    TIMEOUT  dim=None 900.00s  vars=8  gens=9

TARGET -- saturated stratum unions in cell a7_b3000_T2  [satGr_p10007]
  T_dd2le1_dsig7_satGr_p10007   UNIT     dim=-1     0.31s  vars=6  gens=9
  T_dd2le2_dsig7_satGr_p10007   UNIT     dim=-1     0.27s  vars=7  gens=9
  T_dd2le3_dsig7_satGr_p10007   UNIT     dim=-1     6.52s  vars=8  gens=9
  T_dd2le4_dsig7_satGr_p10007   TIMEOUT  dim=None 600.03s  vars=9  gens=9
  T_dd2le1_dsig8_satGr_p10007   TIMEOUT  dim=None 600.03s  vars=7  gens=9
```

The last two lines are the decisive ones, and they close the question for this
engine:

* at the cell's **full `d2` width** the `f31` engine does not finish **even
  mod p**;
* worse, the **`deg sigma = 8` half of the cell does not finish even at the
  NARROWEST `d2` width and with only 7 ring variables** — one variable *fewer*
  than a row that finished in 6.5 s. So the obstruction is **not variable
  count**: freeing `sigma`'s top coefficient (`deg G = 2` instead of `1`, i.e.
  `deg sigma = 8`) is what explodes the basis.

Since 5 of the cell's 10 states have `deg sigma = 8`, **the flag-case cannot be
closed by the `f31` engine at any width** — that is a measured cost statement,
not a claim about the geometry.

Two side experiments, both negative and both worth recording:

```
  T_dd2le1_dsig7_satGr_n16   UNIT     dim=-1     1.80s  vars=6   gens=17   (16 master coeffs)
  T_dd2le1_dsig7_genSig      TIMEOUT  dim=None 900.03s  vars=11  gens=9    (sigma fully generic)
```

* **Over-determination does not help here.** `MODULAR_TRIAGE.md`'s lesson ("more
  generators collapse the unit ideal faster") does not reproduce on this ladder:
  16 master coefficients made the narrowest row *slower* (1.80 s vs 0.34 s).
* **Dropping the pattern-B `sigma` ansatz is unaffordable.** With `sigma` fully
  generic (11 ring variables) even the **narrowest** union times out at 900 s
  exact-Q. So the ansatz-independence of §5's kill is not merely unproven — it is
  **measured to be out of reach with this engine**, which is a stronger and more
  useful statement than "not run".

**What this establishes.** Exact over Q, the quasi-affine locus

```
{ f31 ≡ 0 } ∩ { deg d2 ∈ {0,1,2} } ∩ { deg sigma = 7 } ∩ { gamma ≠ 0 } ∩ { G(r) ≠ 0 }
```

is **empty**, in one 8.83-second computation covering **3 of the cell's 10
states**. Two of the three are the ledger's existing kills; the third,

```
sub2|7|(3, 0, 0, 0)|T2|False|False|()|-inf|2|10|7
```

is **new** (`claimed`, pending audit — and note it is "new to the ledger" partly
because `modular_triage.build_system4` sampled this cell under a quota of 2, not
because the per-state route could not have reached it).

**What this does not establish.** The flag-case does not close. `deg d2 ≤ 3` is
where exact Q stops; the mod-p `UNIT` at that width is a *prediction* of
emptiness over `F̄_10007`, and by the discipline of `MODULAR_TRIAGE.md` that is
reconnaissance, not a kill.

---

## 6. Cost — where it actually breaks

Three separate walls, and none is the one that was expected:

**Wall 1 — the exact-Q Gröbner saturation.** `deg d2 ≤ 2 → 8.83 s`;
`deg d2 ≤ 3 → > 900 s`. One extra free coefficient in `d2` moves the exact
saturation by at least two orders of magnitude, while mod 10007 the same system
takes 6.52 s. This is the standard char-0 coefficient blow-up the repo has met
before, arriving at the first genuinely relaxed width.

**Wall 1b — and it is not variable count.** `deg d2 ≤ 1, deg sigma = 8` has
**7** ring variables and times out mod 10007 at 600 s, while `deg d2 ≤ 3,
deg sigma = 7` has **8** and finishes in 6.52 s. Freeing the top coefficient of
`sigma` is what explodes the basis, not adding an unknown.

**Wall 2 — building the `f31` generators at all.** The master-coefficient
convolution for `d2` of degree 4 (5 free coefficients feeding
`d0 = (d2²+sigma)/4`, then `h_f` of degree up to 7 in `d0`) runs for **tens of
minutes in sympy before Singular is even called** — the documented
cubic/higher-expansion trap, hit from a new direction.

**And the surprise: the *bigger* system is the cheaper one.** The `gsys` engine
(`face_kill_sweep.build_state_system`, 45 spare unknowns at the sub2 caps
12/14/16) builds the *whole* relaxed cell in **290 s** and a single stratum in
**61–70 s**, because the `G`-generators are **cubic** in the window variables
while `f31` is degree 21. Sizes and the first solve:

```
G-SYSTEM ENGINE (G1,G2,G3,G5body+Phi, spares at sub2 caps 12/14/16)
    build 69.94s, 122 equations, 45 spare unknowns
  GS_state_dd20_dsig7      UNIT  dim=-1  112.89s  vars=50  gens=123
    build 89.78s, 122 equations, 45 spare unknowns
  GS_union_dd2le1_dsig7    UNIT  dim=-1  200.89s  vars=51  gens=123
```

| gsys configuration | build | equations | ring vars | generator monomials |
|---|---:|---:|---:|---:|
| single stratum `dd2=0, dsig=7` | 61.3 s / 69.9 s | 122 | 50 | 16 019 |
| whole relaxed cell `dd2≤4, degG≤2` | 290.3 s | 122 | 55 | 40 667 |

`GS_state_dd20_dsig7` **reproduces the ledger kill through the full 50-variable
pre-resultant G-system** (mod 10007, 112.89 s) — the first time in this file's
lane that the G-system route has been driven to a verdict on a frontier state
with the spares present. **This is the direction worth pursuing**, and it is the
opposite of what the cost model in the queue assumed.

---

## 7. Honesty — what is NOT established

1. **The flag-case is not closed.** 0 of 220 sub2 cells were closed before this
   run and 0 are closed after it. The pilot closes a 3-stratum *sub*-union of one
   cell exactly, and predicts more mod p.
2. **The `sigma` ansatz is inherited, and it is load-bearing.** The shape
   `sigma = (y-r)^{2Σb} G` with `G(r) ≠ 0` comes from
   `modular_triage.build_system4` / `triage_harvest.build_sys4_8`, and
   `MODULAR_TRIAGE.md` §"System 4" flags it in its own words: *"Ansatz is a
   faithful generalization of the R9 recipe, not a verbatim landed construction;
   treat as reconnaissance."* **Every exact-Q kill in §5 is exactly as sound as
   the two ledger kills it extends — no more.** A `--sigma-generic` mode that
   drops the divisor structure entirely (fully generic `sigma`, endpoint ideal
   `<s_7,s_8>`) is implemented and was run; it **times out at 900 s on the
   narrowest union in the ladder** (`T_dd2le1_dsig7_genSig TIMEOUT 900.03s
   vars=11`). So **the ansatz-independence of the result is NOT established, and
   is measured to be beyond this engine's reach.**
3. **`G(r) ≠ 0` narrows the claim.** Saturating by it removes
   `v_r(sigma) > 2Σb` from the scheme, so the kill covers
   `cell ∩ {v_r(sigma) = 2Σb}`. `--no-sat-Gr` drops it; implemented, not run.
4. **Mod-p rows are not kills.** `UNIT` mod 10007 says empty over `F̄_10007`.
   Reported as prediction only, and never rolled into a kill count.
5. **Kills are `claimed`.** The one new state (`deg d2 = 2, deg sigma = 7`) is
   same-author, produced by a pipeline whose controls pass but which has not been
   independently re-derived. It is not entered in `state_kill_ledger.json` by
   this file.
6. **Timeouts are cost, not verdicts.** `T_dd2le3_dsig7_satGr TIMEOUT 900.00s`
   means the computation did not finish. Nothing is inferred from it.
7. **No closed-cell control exists to run.** §2 explains why; this is a real gap
   in the validation, and the C1/C1b/C1c/C2 battery is what stands in for it.
8. **A corollary that would follow, if the cell closed.** The sibling cell
   `sub2 a7_b3000_T2 g_zero_levels=[4]` (6 states) has degree tuples that are a
   *subset* of this cell's 10, and the `f31` master identity does not involve the
   `g`-chain at all (`phase_f2_sub2.py`: *"the master identity ignores the
   g-chain"*), so any `f31`-level kill of a degree tuple kills it in both cells.
   That corollary is recorded here because it is cheap, **not** because the
   premise has been discharged.

---

## 8. Verdict against the four outcomes the brief allowed

* not `UNIT saturation → the flag-case closes`;
* not `a few components`, and not `components reproduce nearly every state` — the
  saturation never returned a proper basis to decompose;
* **`COST (does not finish)`** — with the sizes, the timings, the three distinct
  walls, and a concrete, evidence-backed next step (§6: the `gsys` engine is
  cheaper to *build* than `f31` on the relaxed ansatz and has already returned
  `UNIT` on both a state and a saturated union at 50–51 variables).

The decisive negative is sharper than a bare timeout: the `f31` engine fails on
the `deg sigma = 8` half of the cell **at the narrowest width and at 7 ring
variables**, so no widening schedule, prime choice, or coefficient budget
rescues it. Closing this flag-case requires a different engine — and §6 says
which one, with numbers.

The idea itself survives this pilot in good shape: the saturated union is the
**right object** (C1c), it is **decidable at small width**, it delivered a state
the per-state ledger did not have, and the obstruction is arithmetic cost in a
place that has a known lever (over-determination, `msolve`, mod-p→exact lifting,
or the G-system route) rather than a conceptual defect.

---

## 9. Reproduce

```
cd d2_plane_72_108
python saturated_cell.py --controls                              # 7/7, ~11 min
python saturated_cell.py --target --timeout 900                  # exact-Q ladder
python saturated_cell.py --target --chars 10007,10009,100019     # recon ladder
python saturated_cell.py --target --rows 2:1:2 --chars 10007,0   # one union only
python saturated_cell.py --target --sigma-generic                # ansatz-independent variant
python saturated_cell.py --target --no-sat-Gr                    # drop v_r(sigma) exactness
python saturated_cell.py --gsys   --gsys-char 10007              # G-system engine
```

`--rows` selects ladder entries as `dmax:gmin:gmax` (`2:1:2` = `deg d2 ≤ 2` and
`deg sigma ∈ {7,8}`). `--chars` builds the generators once and solves at each
characteristic. `--out` lets parallel runs keep separate artifacts; a run merges
into an existing artifact by label rather than clobbering it.

Exit code is nonzero if any control mismatches its recorded expectation.

**Standing traps this run re-confirmed**, for whoever picks it up: Singular
parses `x^8/N` as `x^(8/N)`, so `poly_to_singular_exact` clears denominators and
now *raises* rather than truncating a non-integer coefficient; `sat(I,B)[1]`
needs `LIB "elim.lib"`; and a `subprocess` timeout kills `wsl.exe` but leaves the
`Singular` process alive inside WSL — check `wsl -d Ubuntu -- ps -eo pid,etime,comm`
for orphans after any timed-out run.
