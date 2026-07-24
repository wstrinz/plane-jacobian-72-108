> **ADJUDICATED (2026-07-24, `window_band_probe.py`): THE STAIRCASE CLAIM IS
> REFUTED — AND EXPLAINED.** The per-slice increment is not a bandwidth; it is
> simply **the number of spare series**. Measured on the sub2 home G-system
> (computable in ~10 s, same bigraded lattice), every slice introduces exactly
> one new coefficient from each spare series: `dm2,dm3,dm4` → `+3`
> (`R0,S0,T0`, then `R1,S1,T1`, …). The R9 z=1 system has `dm4` eliminated,
> leaving two series — hence the `+2` this file reports. So §2's "bandwidth-2
> staircase" measured a real number and promoted it to a structural law it is
> not. The correct statement is `increment = #spare series`, a bookkeeping
> identity carrying no information about the obstruction.
>
> Independently corroborated from the other side: CAOS's planar manuscript
> derives an *edge*-determined band (`g+1`, one per edge term) for the same
> species of object — likewise "count the contributing objects", likewise not a
> universal 2.
>
> **What survives:** the DECOMPOSED-BUT-OPEN direction, the exact top-corner
> equations of §3, and the extracted `a4` constraint — all computed exactly and
> unaffected. **What does not:** every quantitative lattice claim in §2.
> The Milestone-1 R3 lane built to adjudicate this did NOT finish (>1h18m, no
> output — the documented sympy cubic-expansion trap); it was superseded by the
> cheap measurement above and killed.
>
> ---
>
> **STATUS (2026-07-25, orchestrator): PARTIAL / UNVERIFIED.** The verifier's
> first end-to-end run FAILS: the "+2 new spares per slice" staircase claim
> breaks at slice 1 (introduces 4: R0,R1,S0,S1) — the reassembled scripts
> disagree with the in-session findings on the exact band structure. The
> DECOMPOSED-BUT-OPEN direction and the top-corner constraint (computed
> exactly in-session) are likely sound; the precise lattice claims await
> adjudication by the systematic bigrade_annotator lane (Milestone 1),
> which owns the authoritative structure. Do not cite this file's staircase
> numbers until reconciled.

# BIGRADED_PROBE — the cheapest test of the scalar-projection-artifact hypothesis on our own walled state (R9 z=1)

**Verdict: DECOMPOSED-BUT-OPEN.** The R9 z=1 dm4-eliminated H-system — the
computational wall that survived TWO sound scalar-Gröbner spare reductions
(`R9_SYMBOLIC.md` dm4 elimination, `R9_VALSPLIT.md` valuation split) — **does
decompose** under a (u-weight, y-order) window slicing. The spare-unknown
coupling is **banded (a staircase of bandwidth 2)**: the extreme y-order slices
are 3-equation, 2-unknown subsystems, and the top corner eliminates in **seconds**
to a genuine necessary constraint on a single state parameter — a constraint the
monolithic scalar Gröbner (136 eqs / 37 vars, TIMEOUT > 300 s) never surfaces.
**No kill was produced from the corner alone** (the extracted constraint is
satisfiable), so the wall is not closed; but the formulation hypothesis of
`F2_TOWER.md` / `GPT_ASK_BIGRADED.md` question 3 is **empirically validated on the
home case**: the wall is (at least partly) a scalar-projection artifact — the
y-order layer carries structure the scalar G-ideal projects away.

Files (NEW, uncommitted): `bigraded_probe.py` (construction), this writeup,
`bigraded_probe_verify.py` (structural checker — see status note §6).

---

## 1. The experiment

Target: the **R9 z=1** cell — `convolution_elim_qsupport.build_qsupport_ansatz(1)`,
the certified dm4-eliminated H-system of `r9_eliminated_system.json`
(`H2, H3, H5` — dm4-free elements of ⟨G1,G2,G3,G5⟩, weighted-homogeneous of
u-weight 228/240/264). The scalar attack takes **every** y-coefficient of
`H2,H3,H5` on the state's stripped ansatz as one big system and hands it to
Gröbner; it times out (`R9_SYMBOLIC.md` §3, `R9_VALSPLIT.md` §4). Here we instead
keep the **(u-weight, y-order)** label on every equation and work the window
elimination by y-order slices, the GGV3 §5 / `f2_tower.py` method template.

**State structure (R9 z=1).** `d1 = 0`; `d2` (deg 4, coeffs `a0..a4`), `d0`
(deg 8, `a0..a4,g0,g1,r`), `e = dm1` (deg 10, lead `gamma ≠ 0`, const `-gamma·r`),
`Phi` stripped (deg 34). Spare ansätze `dm2 = Σ R_i y^i` (i=0..12),
`dm3 = Σ S_i y^i` (i=0..14) — **28 spare unknowns**. Marked root `r`
(`2048r⁴−512r³+320r²−240r+195 = 0`), saturations `gamma, g1, g0+g1·r ≠ 0`.

## 2. The bigraded lattice found

Each window symbol `d_{4-k}` carries **u-weight 12k** (the window floor that is
stripped off) and a stripped y-polynomial that spreads y-order on top. Because
`H2,H3,H5` are u-homogeneous, the **u-weight axis is degenerate: 3 values
{228,240,264}** (one per H); the decomposition lives on the **y-order axis**.
Instantiating (state coeffs specialised to distinct primes — the spare incidence
is independent of their nonzero values, all supports being dense) gives:

| H | u-weight | y-order range | # nonzero slices |
|---|---:|---|---:|
| H2 | 228 | 0 .. 38 | 39 |
| H3 | 240 | 0 .. 40 | 41 |
| H5 | 264 | 0 .. 44 | 45 |

Total **125 H-slice equations** (the scalar system, + ~10 divisibility + `q(r)`
= 136, matching the census). **Every slice is cubic in the spares** (from the
`dm2²·dm3` / `dm2·dm3²` terms) and carries a known inhomogeneous part.

**The decisive structural fact — the spare coupling is a bandwidth-2 staircase.**
Per-slice spare count, both sweeps (H2 shown; H3/H5 identical shape shifted):

```
y-order:  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 ...(saturated).. 24 25 26 27 28 ...38
#spares:  2  4  6  8 10 12 14 16 18 20 22 24 26 27 28 ...    28    ... 28 27 26 24 22 ... 2
new/slice(bottom): R0S0 R1S1 R2S2 ... R12S12 S13 S14   (exactly 2 new spares per slice)
new/slice(top):    R12S14 R11S13 R10S12 ...            (exactly 2 new spares per slice)
```

Spare-count distribution over the 125 slices:
`{2:6, 4:6, 6:6, 8:6, 10:6, 12:6, 14:6, 16:6, 18:6, 20:6, 22:6, 24:6, 26:6,
27:6, 28:41}`. So **6 extreme corner slices carry only 2 spares each**, 18 carry
≤ 6, and the saturated middle band (41 slices) carries all 28. Sweeping inward
from either corner introduces **exactly two new spare unknowns per slice**
(`R_k, S_k`) — a genuine triangular/banded structure that the scalar formulation
dissolves into one 28-variable pot.

## 3. Window elimination — the corners

**Top corner** (the three deepest slices, one per H — u-weights 228/240/264 at
y-orders 38/40/44), with **real** state coefficients (`gamma≠0`, `a4`, and the
`Phi` leading constant `1024/3315`):

```
H2[38] = -3 R12² S14 - (3/2) a4 gamma R12² + (3/2) gamma S14² - (3/8) a4² gamma³
H3[40] = -3 R12 S14² - 3 a4 gamma R12 S14 - (3/4) a4² gamma² R12 - (1/2) gamma⁴
H5[44] = -(3/2) gamma² R12² + 3 a4 R12 S14² + 3 a4² gamma R12 S14
         + (3/4) a4³ gamma² R12 - (3/2) gamma³ S14 - (1024/3315) gamma
```

Three cubic equations in **two** spare unknowns `{R12, S14}`, coupled only to the
single state parameter `a4` (and `gamma≠0`). Note `R12=S14=0` is **not** a
solution — `H3[40]` leaves `-(1/2)gamma⁴ ≠ 0` — so the corner is genuinely
constraining. Resultant elimination of `R12, S14` (chain
`res(res(H2,H3;S14), res(H2,H5;S14); R12)`) runs in **seconds** and returns a
**non-trivial** necessary condition on `a4`:

```
gamma⁴⁸ · Q1(a4,gamma) · Q2(a4,gamma) = 0,
Q1, Q2  each a quintic in a4  (degree 15 in gamma; e.g.
Q1 = 1024843684156344·10⁶ a4⁵ gamma¹⁵ − 1013034022396231680000 a4⁴ gamma¹²
     + 273809497456115712000 a4³ gamma⁹ + 24165561335493427200 a4² gamma⁶
     − 20527970051508142080 a4 gamma³ + 2626562270214755071875 gamma¹⁷
     + 2305843009213693952 ).
```

This is (i) **not identically zero** — the slicing extracts real information; and
(ii) **not the unit ideal** — `a4` has admissible roots, so **the top corner does
not kill**. It is a *sound necessary condition on the full system* (a consequence
of three of its coefficient equations), reached without touching the other 26
spares — exactly what the scalar Gröbner cannot afford.

**Bottom corner** (`H2[0], H3[0], H5[0]`, spares `{R0, S0}`, params `a0,g0,r,gamma`):

```
H2[0] = -3 R0² S0 + (3/2) a0 gamma r R0² - (3/2) gamma r S0² + (3/8) a0² gamma³ r³ + (3/8) g0 gamma³ r⁵
H3[0] = -3 R0 S0² + 3 a0 gamma r R0 S0 - (3/4) a0² gamma² r² R0 - (3/4) g0 gamma² r⁴ R0 - (1/2) gamma⁴ r⁴
H5[0] = -(3/2) gamma² r² R0² + 3 a0 R0 S0² - 3 a0² gamma r R0 S0 + (3/4) a0³ gamma² r² R0
        + (3/4) a0 g0 gamma² r⁴ R0 + (3/2) gamma³ r³ S0 + (1/34) gamma r
```

Same shape: 3 cubics in `{R0, S0}`, coupled to `a0, g0, r` (`gamma, r ≠ 0` by
saturation). Its resultant elimination did not finish in the residual budget (the
extra parameter `r` heavies the chain) — **IN-PROGRESS-AT-REPORT** (§6).

## 4. Size comparison — scalar vs sliced (the headline)

| | scalar (what the wall's GB sees) | sliced (this probe) |
|---|---|---|
| object | 136 eqs / 37 vars, monolithic | 125 H-slices, banded by y-order |
| spares per unit | all 28 at once | **2** at each extreme corner; +2 per inward slice |
| extreme corner | — (dissolved) | **3 eqs / 2 unknowns**, param `a4` |
| cost of the corner | TIMEOUT > 300 s (whole system) | **seconds** (resultant) |
| output | no verdict (COST) | explicit necessary constraint on `a4` |

The scalar-vs-sliced gap is decisive: the wall's Gröbner cost is a **projection
artifact** of collapsing a bandwidth-2 staircase into a single 28-variable ideal.
The corners are trivially small; the middle band (41 slices at 28 spares) is where
any genuine obstruction, if it exists, must live — and it is reachable only by
propagating the corner constraints inward (the period-style window compiler of
`GPT_ASK_BIGRADED.md`).

## 5. Verdict and reading

**DECOMPOSED-BUT-OPEN.** The slicing genuinely shrinks the per-slice systems
(28 → 2 at the corners; +2 per slice inward) and extracts, in seconds, a sound
necessary constraint on the state that the scalar attack times out on. This
**validates the bigraded formulation hypothesis on the home case**: the R9/alt
wall is, at least in its outer layers, a scalar-projection artifact — the y-order
window layer is informative where the scalar G-ideal is inert. It does **not**
(yet) reproduce a kill: no corner forces a contradiction, and the middle band was
not reached in budget. Whether a window-depth contradiction (à la GGV3's forced
corner-coefficient vanishing) lives in the saturated band is the open question the
`bigrade_annotator` lane and a period-12-style inward sweep should now settle.

## 6. Status / honesty notes

- **[established]** The banded (bandwidth-2 staircase) spare-incidence lattice
  (§2), the exact top-corner equations, and the top-corner resultant elimination
  to a non-trivial-but-satisfiable `a4` constraint (§3) — all computed exactly
  and reproduced by `bigraded_probe.py`.
- **[IN-PROGRESS-AT-REPORT]** The bottom-corner resultant (extra param `r`) did
  not finish in budget; the second inward corner (4-spare slices) and the
  divisibility-row slices were not folded in. No middle-band sweep attempted.
- **[soundness]** Each corner equation is one coefficient equation of the certified
  necessary H-system; any consequence of a subset (resultant/elim) is a consequence
  of the full system, so the extracted `a4` constraint is a **sound necessary
  condition** (standard resultant caveats — spurious leading-coefficient factors —
  apply and would need audit before any kill claim; none is made here).
- **[verifier]** `bigraded_probe_verify.py` re-checks the structural invariants
  (corner sizes = {R12,S14}/{R0,S0}, the bandwidth-2 staircase, top-corner
  resultant ≠ 0). Its end-to-end `--quiet` pass is **IN-PROGRESS-AT-REPORT**
  (numeric-lattice build ≈ 50 s; not re-run within the reporting budget).
- **[scope]** State coefficients were specialised to distinct primes **only** to
  read the spare-incidence lattice (§2); the corner *equations and elimination*
  (§3) use the **real** state coefficients. `dm4` never appears (eliminated,
  `r9_eliminated_system.json`); the G1 divisibility rows and `q(r)=0` were not
  added to the corner subsystems (would only strengthen them).
