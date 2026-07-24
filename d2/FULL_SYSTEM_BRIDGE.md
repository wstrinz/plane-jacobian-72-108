# FULL_SYSTEM_BRIDGE — from a cascade STATE to the FULL necessary system (2026-07-23)

**Mission.** The cascade / convolution engine tests one necessary condition,
`f31 = 0`. But `f31` is only the **elimination-ideal generator** of the
pre-resultant window ("G-") system in the variables `(d~2, d~1, d~0, e, Phi)`
(`F37_SATURATION_REPORT.md` fact [5]: `E := <G-system> ∩ Q[d~2,d~1,d~0,e,Phi]
= <f31>`). The G-system holds **strictly more information**: it constrains three
further window unknowns `(d_-2, d_-3, d_-4)` that `f31` has eliminated. A state
solvable under `f31`-alone is **not** a counterexample germ until it lifts
through the full G-system with those spare unknowns realised as honest
polynomials of the window-bounded degree. This doc + `full_system_bridge.py` +
`full_system_bridge_verify.py` build and verify that bridge, and pilot it.

**Headline result (pilot, verified).** On the `a8` constant-E stall class —
`MODULAR_TRIAGE.md` System 2, the program's likely-first f31-survivor
("LIKELY-SOLVABLE") — the **full G-system is EMPTY: an exact `UNIT` ideal over
`Q` in ≈ 1–9 s**. The full system **kills** the germ. See §5.

---

## 1. The variable dictionary (verbatim)

The cascade unknowns and the G-system indeterminates are the *same* window
d-variables `d_k := c_k·C4^(7−2k)` (`d4 = 1`, `d3` killed by the shift), read in
two notations. Verbatim from `full_system_bridge.py`:

```
cascade / window variable  |  G-system indeterminate  |  window meaning  |  k
  d2      (cascade d2)      |  d2                      |  d~2  = d_{4-2}  |  2
  d1      (cascade d1)      |  d1                      |  d~1  = d_{4-3}  |  3
  d0=(d2^2+sigma)/4         |  d0                      |  d~0  = d_{4-4}  |  4
  e       (cascade e)       |  dm1                     |  d_-1 = d_{4-5}  |  5
  --- spare window unknowns the bridge introduces ---
  r                         |  dm2                     |  d_-2 = d_{4-6}  |  6
  s                         |  dm3                     |  d_-3 = d_{4-7}  |  7
  (dm4)                     |  dm4                     |  d_-4 = d_{4-8}  |  8
  Phi = c t^30 q            |  Phi                     |  Phi/y^204       |  -
```

- **`e = dm1 = d_-1`** and **`sigma := 4·d~0 − d~2²`** are the cascade's
  reparametrisations (`t5_90t1_verify.py` lines 16–17). The engine works in
  `(d0,d1,d2,e)`; `sigma` merely re-coordinatises `d0`.
- **`r, s, dm4` are `d_-2, d_-3, d_-4`** — the spare unknowns the mission asked
  to identify. They are read *directly* from the G-system generators in
  `t4_state.pkl`: F37's report writes `e=d_-1, r=d_-2, s=d_-3`, and the fourth
  generator variable is `dm4=d_-4`. Confirmed against the ground-truth
  `regenerate_system.py` `S = 1 + d2·u² + d1·u³ + d0·u⁴ + Σ dm[k]·u^{4+k}`
  (so `dm[k] = d_{-k}`). **Not guessed — loaded and cross-checked
  (`full_system_bridge_verify.py` V1).**
- **`Phi = c·t^30·q`, `c = −1/6630`, `t = y+1`, `q = 2048y⁴−512y³+320y²−240y+195`**
  is the cascade's **stripped** Phi `= Phi_full/y^204` where
  `Phi_full = f1·C4^28` (deg 238, ord 204). Verified `Phi_full/y^204 = c·t^30·q`
  (`verify_derivation.py` §A; re-checked in `..._verify.py` V3).

The G-system generators, loaded from `t4_state.pkl` (cleared `/2` shown here as
in the source):

```
G1 = 3/2·d1·dm1² + 3·d2·dm1·dm2 + 3·dm1·dm4 + 3·dm2·dm3
G2 = −3/2·d0·dm1² + 3/2·d2·dm2² + 3·dm2·dm4 + 3/2·dm3²
G3 = −3·d0·dm1·dm2 − 3/2·d1·dm2² − 1/2·dm1³ + 3·dm3·dm4
G5 = Phi + G5body,     G5body = −3·d0·dm1·dm4 − 3·d0·dm2·dm3 − 3·d1·dm2·dm4
                                 − 3/2·d1·dm3² − 3·d2·dm3·dm4 − 3/2·dm1²·dm3
                                 − 3/2·dm1·dm2²
```

> **ERRATUM (2026-07-24).** This block previously read `G5 = 2·Phi + G5body`,
> contradicting the augmentation recipe below (`(G5body+Phi)`) and the canonical
> loader `full_system_bridge.py` (`st["G5body"] + PHI`). The authority is the C11
> membership certificate in `f37_sat_verify.py`, which verifies
> `f31 == c1·G1 + c2·G2 + c3·G3 + c4·(G5body + Phi)`. Corrected to `Phi + G5body`.
> The two forms differ by `Phi`, not by a nonzero scalar, so they are genuinely
> different equations. `bigrade_annotator.py` transcribed the erroneous form; see
> `FACE_KILL_SWEEP.md` §4 for the impact assessment (no landed kill changes; the
> M1 R2 certificate VALUE was wrong, `−1024/3315` not `−2048/3315`).

These are `(D~³)_{-1,-2,-3,-5}` after the `(D~²)` linear substitutions — exactly
`regenerate_system.py`'s system, validated by `T6_SELECTION_AUDIT.md`.

---

## 2. Degree / order caps (stripped coordinates)

**Full-window bounds (T3_WINDOW_AUDIT.md §3), per variable `d_{4-k}`:**

| var | window | k | ord ≥ | deg ≤ (sub1) | deg ≤ (sub2) |
|---|---|---:|---:|---:|---:|
| d~2  | d2  | 2 | 24 | 30  | 28  |
| d~1  | d1  | 3 | 36 | 45  | 42  |
| d~0  | d0  | 4 | 48 | 60  | 56  |
| d_-1 | dm1=e | 5 | 60 | 75  | 70  |
| **d_-2** | **dm2=r** | **6** | **72** | **90**  | **84**  |
| **d_-3** | **dm3=s** | **7** | **84** | **105** | **98**  |
| **d_-4** | **dm4**   | **8** | **96** | **120** | **112** |

`ord ≥ 12k`, `deg ≤ 15k` (sub1) / `14k` (sub2). T3 states these explicitly for
the window `k=2..5`; the `k=6,7,8` rows are the **[judgment]** extension of §4.

**Stripped caps** (the coordinates the cascade and the bridge actually use;
`V_stripped := V_full/y^{12k}`, subtracting `12k` from every bound → `ord ≥ 0`,
`deg ≤ 3k` (sub1) / `2k` (sub2)):

| spare var | k | stripped deg ≤ sub1 | #coeffs sub1 | stripped deg ≤ sub2 | #coeffs sub2 |
|---|---:|---:|---:|---:|---:|
| dm2 = r  | 6 | 18 | 19 | 12 | 13 |
| dm3 = s  | 7 | 21 | 22 | 14 | 15 |
| dm4      | 8 | 24 | 25 | 16 | 17 |
| **total spare unknowns** | | | **66** | | **45** |

**This is the ansatz cost of the bridge: +45 (sub2) or +66 (sub1) fresh scalar
unknowns per state**, on top of the state's own coefficients.

---

## 3. The augmentation recipe

Given a **state** = flags + degree assignment (a `convolution_descent.Ansatz`
giving stripped `d2, d1, sigma, e`; `d0=(d2²+sigma)/4`), and a regime
(`sub1`/`sub2`), the full necessary system is:

```
FULL(state) =
    [ optional: the state's f31 master-coefficient equations ]       (redundant*)
  + [ every y-coefficient of G1, G2, G3, (G5body+Phi) ]              (the content)
      with  Phi := c·t^30·q,   c = −1/6630,
      and   dm2, dm3, dm4  introduced as generic stripped polynomials
            of degree ≤ (2k or 3k),  their coefficients fresh unknowns
  + [ gamma ≠ 0 ]  (Rabinowitsch w·gamma−1, and any state nonzero-parameter )
```

`*` f31 lies in the G-ideal (§4), so its coefficients are **implied** by the
G-system coefficients; the bridge's default is `nf31=0` (**pure G-system**), and
the pilot confirms the pure G-system already carries the kill (§5).

**Why stripped substitution is exact.** Each `G_i` is **weighted-homogeneous**
under `w(d_{4-k}) = 12k`, `w(Phi) = 204`:

```
G1: weight 156   G2: weight 168   G3: weight 180   G5: weight 204
```

(verified in code and in `..._verify.py` V3). Hence, writing
`V_full = y^{w(V)}·V_stripped`, every monomial of `G_i` carries the same power
`y^{W_i}`, so

```
G_i(V_full) = y^{W_i} · G_i(V_stripped)   ⇒   G_i(full)=0  ⟺  G_i(stripped)=0.
```

The bridge therefore substitutes stripped ansätze and sets the y-coefficients to
zero — equivalent to the full-window equation, and matching the cascade's own
stripped `Phi` and stripped state ansatz.

**Implementation.** `full_system_bridge.augment(ansatz, regime)` returns the
equation list, unknowns, and sizes. `triage(...)` runs a 3-prime mod-p Singular
verdict (`UNIT`=empty over `F̄_p`); `exact_kill(...)` attempts an exact `UNIT`
over `Q`. The systems are small (≈122 eqs / ≈52–56 unknowns) and are emitted
**integer-cleared** (never a `/` coefficient), reusing `modular_triage.py`'s
emitter for the mod-p reductions.

---

## 4. Soundness — every added equation is a proven necessary condition

The bridge only strengthens the necessary condition the cascade already trusts.
Chain, each link cited:

1. **The G-system is the genuine window system.** `G1,G2,G3,G5body =
   (D~³)_{-1,-2,-3,-5}` after the `(D~²)` linear substitutions — bit-for-bit
   `regenerate_system.py`, whose equation-selection, λ-isolation, and clearing
   exponents are verified by `T6_SELECTION_AUDIT.md` / `verify_derivation.py`
   (48 checks). The spare unknowns `dm2,dm3,dm4` are exactly `d_-2,d_-3,d_-4`,
   window variables of a real counterexample.
2. **`Phi = c·t^30·q` is the genuine instance.** `= Phi_full/y^204`,
   `Phi_full = f1·C4^28` (`verify_derivation.py` §A; re-derived in `..._verify`
   V3). `c=−1/6630` is the forced ODE constant, not a free parameter here.
3. **The bridge is `≥` the cascade.** `f31 ∈ <G1,G2,G3,G5body+Phi>` with an
   explicit exact cofactor certificate (`F37_SATURATION_REPORT.md`;
   re-run in `..._verify` V2). So **every** zero of the G-system is a zero of
   `f31`: the bridge can only *add* kills, never remove a genuine solution.
4. **Stripping is exact** (§3, weighted-homogeneity) — no information lost or
   fabricated by working in the cascade's small-degree coordinates.
5. **The window caps are real constraints.** Over `F̄_p` / `C` with `dm2,dm3,dm4`
   *free*, `f31=0` ⇒ G-system solvable (elimination ideal `= <f31>`). Requiring
   `dm2,dm3,dm4` to be **polynomials of the bounded degree** is a strictly
   stronger, and genuinely necessary, condition (they *are* the window
   variables `d_-2,d_-3,d_-4`, forced polynomial with `ord≥12k`, `deg≤15k/14k`).
   This is precisely the extra content `f31` cannot see.

### [judgment] list (anything inferred / not line-by-line in a cited source)

- **[RETIRED 2026-07-23 — now PROVEN] Window caps for `k=6,7,8` (dm2,dm3,dm4).**
  Formerly [judgment — well-supported]: `T3_WINDOW_AUDIT.md` proved the caps
  explicitly only for the jetlift window `k=2..5`, and the extension to
  `k=6,7,8` was flagged as un-recited. Now recited and proven in
  `WINDOW_CAPS.md` + `window_caps_verify.py` (81 checks, wired into
  run_tests.sh): the three T3 §3 valuation inductions close as symbolic
  identities in `k`, the D-transform arithmetic gives `15k/14k/12k`
  identically in `k`, and the d₃-killing shift preserves the caps
  term-by-term. No conditionality beyond the program's inherited premises
  [P1][P2][P3]. The k=6,7,8 rows sit in the same trust tier as k=2..5.
- **[judgment] Choice `c` fixed vs symbolic.** The bridge fixes `c=−1/6630` (the
  genuine value). Keeping `c` a free parameter would test a one-parameter
  family; the fixed value is the actual (72,108) necessary condition, and is
  what the pilot uses.
- **[inherited premises] Everything below the cascade ansatz.** The GGV1
  `ℓ(P)=R², ℓ(Q)=R³` premises, the α-strip WLOG, and `C4=y⁷(y+1)` normalisation
  remain outline-only (`T6_SELECTION_AUDIT.md` §4) — unchanged by this bridge,
  which lives strictly above them.
- **[reconnaissance] mod-p verdicts.** A `UNIT`/`PROPER` mod-p read is evidence
  over `F̄_p`, not a proof over `Q` (`MODULAR_TRIAGE.md`). The pilot therefore
  escalates the mod-p `UNIT` to an **exact `Q`** `UNIT` (a real certificate).

---

## 5. PILOT — one a8 constant-E stall state

**State** (`batch_convolution_sub2.json`, `a_t=8, T1, deg_e=8, UNRESOLVED`;
index 2 = `a8_dd2-inf_dd1(2)_dsig(7)`; regime **sub2**):

```
e = gamma·(y+1)^8,   d2 = 0,   d1 = b0+b1 y+b2 y²,   sigma = s0+…+s7 y^7,
d0 = sigma/4.
f31-alone top residual (deg 242):  8192·b2² + 9945·gamma³·s7²   (× a nonzero const)
```

`MODULAR_TRIAGE.md` System 2 scores this class **LIKELY-SOLVABLE** from the
top-3 f31 coefficients — the residual `8192 b² + 9945 γ³ s² = 0` is honestly
solvable over a field (its only obstruction is real/positive-definite, invisible
to a closed-field GB). This is the program's likely-first f31-survivor.

**Verdicts (Singular, this run):**

| system | equations | unknowns | mod-p (3 primes) | exact `Q` |
|---|---:|---:|---|---|
| f31-alone, top-3 coeffs (`MODULAR_TRIAGE`) | 3 | 10 | **PROPER → LIKELY-SOLVABLE** | — |
| f31-alone, 16 coeffs (control) | 16 | 10 | UNIT (202/225/297 s) | (too heavy to run) |
| **full G-system (bridge, `nf31=0`)** | **122** | **56** | **UNIT (32/16/12 s)** | **UNIT (8.75 s)** |

Index 0 (`a8_dd2-inf_dd1(0)_dsig(5)`), pure G-system: exact `Q` **UNIT (1.4 s)**.

**Answer to THE QUESTION — the full system kills, decisively.**
The full G-system is an **exact `UNIT` ideal over `Q`** (empty, γ≠0) in seconds:
the a8 constant-E germ **does not lift**. Two honest nuances, reported with care:

1. The kill is **entirely the G-system's** — the pure G-system (`nf31=0`, *no*
   f31 coefficient at all, only `G1..G5` + the bounded `dm2,dm3,dm4` + γ≠0) is
   already `UNIT` over `Q`. The bridge's own content, not smuggled-in f31 data,
   closes the state.
2. f31-alone is **not literally unable** to kill this class — with 16 *deep*
   coefficients it also reaches `UNIT` mod-p. But that costs **~200–300 s per
   prime** on dense degree-238 objects and yields no affordable exact-`Q`
   certificate (`triage_harvest` confirmed only 2 of 24 over `Q`). The bridge
   reaches the **exact** kill on a ~120-equation, quadratic, small-coefficient
   system in **≈ 1–9 s**. The window-bounded `dm2,dm3,dm4` constraints are what
   make the full system collapse so much faster than deep f31.

So the pilot **validates the bridge as the endgame tool**: it converts an
apparent f31-survivor (MODULAR_TRIAGE "LIKELY-SOLVABLE") into a cheap **exact**
kill, and it does so through the *complete* necessary system — retiring the
conceptual risk that "f31 is only one necessary condition." No genuine surviving
germ was found in the a8 class.

---

## 6. Files (new, uncommitted; READ-ONLY on everything else)

- `FULL_SYSTEM_BRIDGE.md` — this spec.
- `full_system_bridge.py` — `gsystem()`, `augment()`, `triage()`, `exact_kill()`,
  `f31_alone_system()` (control), `run_pilot()`. Emits integer-cleared Singular.
- `full_system_bridge_verify.py` — 13 checks, **ALL PASSED**: V1 identical
  G-system (no hand-copy), V2 `f31 ∈ G-ideal` (soundness), V3
  weighted-homogeneity + `Phi` strip, V4 independent re-derivation of the pilot
  `UNIT` on a fresh prime.
- `full_system_bridge_pilot.json` — machine-readable pilot record.

### Reproduce

```
python full_system_bridge.py 2 sub2      # pilot: control + full + exact kill
python full_system_bridge_verify.py      # must end ALL 13 ... PASSED
```
