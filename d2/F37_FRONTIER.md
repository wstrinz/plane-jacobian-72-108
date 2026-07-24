# F37_FRONTIER — scoping map of the f37 branch (2026-07-22)

> **STATUS (2026-07-22):** WHOLLY SUPERSEDED by `F37_SATURATION_REPORT.md`. The saturation this document calls "not tried" was executed and closed the entire f37 branch: `f31` lies in the pre-resultant ideal, so f37 (and `d₋₁²¹`) is a classical resultant artifact for both subcases. Every open question below is resolved; this file is retained as the scoping analysis that led there.

> **Superseded (2026-07-22, later same day).** The saturation program recommended here TERMINATED: `f31` lies in the pre-resultant ideal, so the entire f37 branch is a resultant artifact — see `F37_SATURATION_REPORT.md`. This document is retained as the scoping analysis that led there; its open questions are all resolved.

Scoping/mapping document for the `f37 ≡ 0` branch of the master identity
`f31 · f37 · d₋₁²¹ ≡ 0`. This file only reads and lightly verifies; it proves
nothing new. Provenance tags used throughout:

- **[P]** proven in the repo (with source file/section);
- **[N]** numerically supported only (not a proof);
- **[I]** my inference / gap I could not source;
- **[V]** freshly verified in this session with light sympy (command in text).

Notation follows `T5_F37_GRADED.md §0`: stripped windows `d̃k = dk/y^{12wk}`,
`e := d̃₋₁`, `t := y+1`, `q := 2048y⁴−512y³+320y²−240y+195` (irreducible over
`Q`, splits into 4 simple linear factors over `Q̄`), `Φ̃ = c·t³⁰·q`,
`c = −1/6630`, `v_t(Φ̃)=30`, `v_q(Φ̃)=1`, `deg Φ̃ = 34`, `σ := 4d̃0 − d̃2²`.

---

## 1. Component map of `{f37 ≡ 0}`

### 1.1 What `f37` is

`f37` is the 618-term second resultant factor
(`STATE.md` item 5, l.47; raw form `f37_deg37.txt`, one line, vars
`d2,d1,d0,m1,P` = `d₋₁,Φ`). **[V]** It is weighted-homogeneous of total
weight **134** under `w(d2,d1,d0,d₋₁,Φ) = (2,3,4,5,17)` — the direct analogue
of f31's weight 125 (`STATE.md` item 6, l.61). Verified this session:
`sympy` over `f37_deg37.txt` gives a single distinct monomial weight `134`
across all 618 terms. So the "grading" of §2 is *not* a homogeneity failure
(see §2).

### 1.2 Restriction to `Φ = 0` — the `dm1^18·h37` factorization **[P]**

`f37|_{Φ=0} = d₋₁¹⁸ · h37` (`T5_NOTES.md` step 2, l.58; `T5_F37_GRADED.md §1`
`[V2]`; `STATE.md` T5 entry l.257). Here `h37 = h_0`, the weight-44,
**145-term** pure-`d` factor (`f37_graded.txt` line for `h_0`;
`T5_F37_GRADED.md §1`). Contrast f31: `f31|_{Φ=0} = d₋₁²¹·h31`, `h31`
weight-20, 28 terms (`T5_NOTES.md` step 2). The `dm1` exponent 18 (vs f31's
21) is the `f=0` entry of the non-uniform power vector `p_f` (§2).

### 1.3 Is `h37` irreducible?

- Repo status: the notes assert irreducibility **only for `h31`**
  (`T5_NOTES.md` step 2, l.60: "with h31 irreducible … h37 weight 44
  (145 terms). Both checked by factorization" — the "irreducible" adjective is
  attached to `h31` only; `h37` irreducibility is *not* stated). **[I]** gap.
- **[V]** Resolved this session: `sympy.factor_list(h_0)` returns
  `27 · (irreducible 145-term, degree-19 factor)` — i.e. only the integer
  constant 27 splits off, and `h37` is **irreducible over `Q`**. (Command:
  parse `h_0` from `f37_graded.txt`, `factor_list` in `d2,d1,d0,dm1`.) This
  matches the f31 side and can be promoted to a repo fact if independently
  re-run.

### 1.4 The bigraded / Newton-polygon decomposition **[P]**

`T5_F37_GRADED.md §1` `[V2]`, mirrored in `f37_graded.txt` header:

```
f37 = Σ_{f=0}^{7} Φ^f · d₋₁^{p_f} · h_f(d2,d1,d0,d₋₁),
p_f = (18,15,12,9,6,4,2,0),   h_f weighted-homog of weight 134−17f−5p_f.
```

`h_f` term counts (145,124,106,88,78,51,25,1) and the terminal collapses
(`T5_F37_GRADED.md §5` `[V6]`):

```
h_7 = 221184·d̃1⁵                       (fifth power; f31: 8192·d1²)
h_6|_{d̃1=0} = −82944·d̃2·σ⁵             (NEW d̃2 factor; σ⁵ not σ²)
h_5|_{d̃1=0,σ=0} = 131072·d̃2²·e⁵         (f31: 2048·e²)
d̃2 | h_f|_{d̃1=0} for every f = 0…7
```

### 1.5 The free family `d2=d1=0` **[P]**

`f37|_{d̃1≡0, d̃2≡0} = 0` identically in `(d̃0,e,Φ)` — equivalently every
monomial of `f37` is divisible by `d1` or `d2` (`T5_F37_GRADED.md §6` `[V6]`;
`STATE.md` 2026-07-22 entry l.324). This is a genuine ~20-dimensional
solution family of the *bare* `f37` identity inside the sub2 windows. f31 has
**no** such family (f31's Lemma 1 forbids total degeneration;
`T5_F37_GRADED.md §6.2`). This is the single deepest structural difference
between the two branches.

### 1.6 The free family does NOT lift **[P]**

Restoring the three smallest pre-resultant equations (`H2,H3,H5` from
`t4_state.pkl`) on `d2=d1=0` yields the compact necessary system
`12rs(r²−es)=e⁵` and `3e(r²+es)=2Φ`; DVR-valuation + degree-window +
infinity arguments force `e=C·t¹⁰`, `R²−C·S=0`, contradicting `(P)`
(`F37_FREE_FAMILY_SYSTEM.md`, whole file; checker
`f37_free_family_verify.py`; `STATE.md` 2026-07-22 l.328). So the *bare*
resultant cannot empty its own free family, but the original system does.

### 1.7 The σ-locus is dead off the free family **[P]**

`{d̃1≡0, d̃0=d̃2²/4, d̃2≢0}`: `32·A′⁴·B′ = 27·e¹⁷` with `A′=2Φ̃+3d̃2e³`,
`B′=4Φ̃+3d̃2e³`; empty by Mason–Stothers (Theorem 2′,
`T5_F37_GRADED.md §7` `[V7]`; field-stable split-place version
`FIELD_SPLIT_AUDIT.md` "split-place sigma-locus theorem"). The `d̃2≡0`
sub-case of the σ-locus *is* the free family again (not a contradiction).

### 1.8 Numerics **[N]**

- All four factor×subcase window systems show genuine positive floors under
  central-difference polish; for f37: sub2 `~1.2e-6` (`STATE.md` T2 entry),
  sub1 `~4.7e-7` (`STATE.md` T3 entries — the tightest floor in the project).
- **Caveat now on record**: the T3 "f37 infeasible" verdict is *provably
  false* on the free-family locus `d̃1=d̃2=0` (those bases are sampled with
  probability 0 by `jetlift.py::make_base`), so it must be re-scoped to
  "infeasible off `{d̃1=d̃2=0}`" (`T5_F37_GRADED.md §6.4`).
- **Key discriminator [N]**: the 21 verified 60-digit numeric solutions of the
  *actual* system (3 random parameter sets) all satisfy `f31=0` and **never
  `f37=0`** (`STATE.md` item 5, l.56–60). f31 is the genuine generic
  component; f37 has never been hit by a real solution.

### 1.9 Exact-fact inventory (with source)

| fact | tag | source |
|---|---|---|
| f37 = 618-term factor, weight 134 homog. `(2,3,4,5,17)` | P/V | STATE l.47,61; f37_deg37.txt; this session |
| `f37|_{Φ=0}=d₋₁¹⁸·h37`, `h37`=145-term weight-44 | P | T5_NOTES l.58; T5_F37_GRADED §1 |
| `h37` irreducible over `Q` | V | this session (factor_list); NOT in notes |
| bigraded `p_f=(18,15,12,9,6,4,2,0)`, terminal collapses | P | T5_F37_GRADED §1,§5; f37_graded.txt |
| free family `f37|_{d1=d2=0}=0` | P | T5_F37_GRADED §6; STATE l.324 |
| free family does not lift to original system | P | F37_FREE_FAMILY_SYSTEM.md |
| σ-locus (`d̃2≢0`) empty (Mason–Stothers) | P | T5_F37_GRADED §7; FIELD_SPLIT_AUDIT |
| Lemma A′ non-uniform cascade exists | P | T5_F37_GRADED §3 |
| Lemma B′ (`q | h_0` when `q∤e`) | P | T5_F37_GRADED §4 |
| terminal degree-starvation kill set EMPTY (all 21 strata survive) | P | T5_F37_GRADED §8 |
| positive numeric floors, all 4 cells | N | STATE T2/T3 entries |
| real numeric solutions never satisfy f37 | N | STATE item 5 |

---

## 2. The grading obstruction, precisely

**f37 IS weighted-homogeneous** of weight 134 under the same 5-variable weights
`(2,3,4,5,17)` that make f31 weight-125 homogeneous ([V], §1.1). So the
obstruction is **not** a homogeneity failure. It is a *Newton-polygon*
(single-variable-substitution) failure, localized to the `d₋₁`-power vector.

**The precise failure** (`T5_F37_GRADED.md §2` `[V3][V4]`,
`T5_NOTES.md` l.186):

- f31 has `d₋₁`-powers `21,18,15,12,9,6,3,0` — a **uniform** drop of exactly 3.
  So `f31 = d₋₁²¹ · H(d̃, w)` with the *single* substitution `w := Φ/d₋₁³`, and
  the certificate becomes "`w` is a rational root of a degree-7 polynomial `H`
  with 8 small coefficients `h_f`" (`T5_NOTES.md` l.164–184). This linearizes
  the whole branch.
- f37 has `d₋₁`-powers `18,15,12,9,6,4,2,0` — the last three (`f=5,6,7`) drop by
  **2, not 3**. In the `(f, p_f)` plane the lower hull has vertices
  `(0,18),(4,6),(7,0)` with **two edges, slopes −3 and −2**; all 8 points sit
  on the boundary (`[V3]`). A single `w=Φ/d₋₁³` therefore does **not** clear
  `d₋₁` from the tail: `f37 = d₋₁¹⁸·H̃(d̃,d₋₁,w)` still carries `d₋₁·h5w⁵ +
  d₋₁²·h6w⁶ + d₋₁³·h7w⁷`. This is exactly what `T5_F37_GRADED.md` calls the
  "non-uniform grading" and what `STATE.md` (l.312) shortcuts as "f37
  (non-uniform grading)".

**A modified grading DOES exist** — two of them, both proven:

1. *Two-substitution normal form* (`§2` `[V4]`):
   `f37 = d₋₁¹⁸·A(d̃,w) + Φ⁵d₋₁⁴·B(d̃,z)`, `w=Φ/d₋₁³` (head, slope-3, `f≤4`),
   `z=Φ/d₋₁²` (tail, slope-2), `A=Σ_{f≤4}h_f W^f`, `B=h5+h6 Z+h7 Z²`. The
   deviation from an f31-style single-`w` object is exactly the tail
   `d₋₁·w⁵·B(z)`.
2. *Non-uniform cascade* Lemma A′ (§3, and see §3 below).

So the honest statement is: **f37 is weight-134 homogeneous like f31, but its
`d₋₁`-Newton-polygon has two slopes instead of one, so the clean scalar
`w`-root reformulation of f31 does not transport; a two-edge normal form or a
non-uniform-step cascade is required instead.**

---

## 3. Ledger requirements for the `(d2,d1)≠(0,0)` subcase

**Does f37 admit a sparse-companion cascade like f31? YES — already extracted
and proven** (`T5_F37_GRADED.md §3`, Lemma A′ `[V5][V8]`). With `e=t^a ê`,
`u:=Φ̃/t³⁰=c·q`, `Δ=(3,3,3,3,2,2,2)`, `δ_ℓ=30−a·Δ_ℓ`:

```
g_1     := h_0(d̃)/t^{δ1}                              (t^{δ1} | h_0 forced)
g_{ℓ+1} := (ê^{Δ_ℓ} g_ℓ + u^ℓ h_ℓ(d̃))/t^{δ_{ℓ+1}}     ℓ=1..6
ê² g_7 + u⁷ h_7(d̃) = 0                                (terminal)
```

The companions `h_ℓ` **are sparse** (`f37_graded.txt`): `h_7 = 221184·d1⁵`
(1 term), `h_6` (25 terms), `h_5` (51 terms), rising to `h_0` (145 terms) —
exactly the shape `cascade_signature.py` already consumes for f31 (it parses
`f31_graded.txt`, rewrites in `(d2,d1,σ,e)`, emits monomial/degree tables).

**Extraction sketch (cheap, mostly done):** point the existing
`cascade_signature.py` at `f37_graded.txt` instead of `f31_graded.txt` — the
file format is identical (header + `h_f (weight …, dm1-power …) = …` lines).
The signature object of `CASCADE_ENGINE_PLAN.md §1` populates as:

```
forcing_divisor : [(t,30,1),(p_i,1,1)×4]          # Φ̃ divisor after base change
master_exponents: p_f = (18,15,12,9,6,4,2,0)
cascade_step    : δ_ℓ = 30 − a·Δ_ℓ, Δ=(3,3,3,3,2,2,2)   # NON-uniform (f31: uniform 3)
window_caps     : deg d̃2≤4, d̃1≤6, d̃0≤8, e≤10  (sub2)
h_degree_caps   : deg h_f ≤ 2·wt_f = (88,84,80,76,72,58,44,30)
terminal_branches: T1 (ê²g7 = −221184 c⁷q⁷ d1⁵), T2 (ê²g6 = 82944 c⁶q⁶ d2 σ⁵), σ-locus
```

The split-place local state `(b,x,z,k,r_ℓ)` and the four-root
degree-budget coupling of `CASCADE_ENGINE_PLAN.md §2–3` transport verbatim;
only the transition arithmetic changes (`Δ_ℓ` non-uniform; terminal q-injection
is `ê²` not `ê³`).

**Two hard differences a split-place ledger for f37 must confront** (both
proven, `T5_F37_GRADED.md §8` `[V8]`):

1. **Terminal degree-starvation prunes NOTHING for f37.** The mechanism that
   killed 81 strata on the f31 ledger (`SPLIT_PLACE_LEDGER.md`) needs
   `4(7−2a_q) > 38+2a` (level 7) or `4(6−2a_q) > 48+2a` (level 6); max LHS is
   28 resp. 24, always below. Reason: f37's caps are ≈2× f31's and the
   terminal q-injection is `ê²` (not `ê³`). So the f37 ledger's *terminal
   pruning column would be all-zero* — every one of the 21 (or, after the
   split-place refinement, 327) strata survives to level-6/5 residue work.
2. **The free family is a permanently-live extra branch.** No f37 ledger built
   from the f37 identity + window caps alone can close the T2 sub-branch
   `d̃2≡0` (`T5_F37_GRADED.md §6.3`). An external constraint is mandatory
   (see §4).

**Smallest computation that would *decide* the cascade is usable:** it is
already decided (Lemma A′ proven, `t5_f37_verify.py` `[V5]`). The next
smallest decisive computation is to run the *head-only* window: below
`t`-order ≈150 the identity is the pure slope-3 system in `h_0…h_4`
(`T5_F37_GRADED.md §2, §11.3`) — four cascade blocks of 30 equations each at
`a=0`, objects of degree ≤88, an f31-shaped Gröbner target. Whether *those
five `h_f` admit a common `q`-adic vanishing* is the concrete gate.

---

## 4. The pre-resultant alternative — is all of f37 an excess component?

**Hypothesis (all-of-f37-fails-to-lift): plausible, and it is the best-founded
route.** Evidence:

- **[N]** The 21 verified numeric solutions of the original system all lie on
  `f31=0`, never on `f37=0` (`STATE.md` item 5). If f37 carried a genuine
  lifting component one would expect to have sampled it.
- **[P]** The one exactly-known solution family of the bare `f37` identity —
  the free family `d2=d1=0` — provably does **not** lift
  (`F37_FREE_FAMILY_SYSTEM.md`). The pattern "resultant component that is an
  artifact of elimination and dies when small pre-resultant equations are
  restored" is already demonstrated once.
- **[P]** No argument from the f37 identity + windows alone can close the
  branch (`T5_F37_GRADED.md §6.3`) — the branch *structurally requires*
  importing an original-system fact, which is precisely the pre-resultant
  method.

Together these make "f37 is entirely an excess/artifact component of the
resultant with no lift" the natural conjecture. It is **[I]** as a global
statement (only the free-family slice is proven).

**Computation that would test "f37 has NO lifting component at all":** the
resultant `f31·f37·d₋₁²¹` is a *necessary* condition extracted from the
G-system `⟨G1,G2,G3,G5+Φ⟩`. Every lift satisfies it, so a lift with `f31≠0`
must have `f37=0`. Hence f37 has no lift **iff** the G-system forces `f31=0`,
i.e.

```
V(G-system) ∩ {f31 ≠ 0} = ∅.
```

The decisive ideal-theoretic test (cheapest exact form):

1. saturate the pre-resultant ideal by `f31`:  `J := ⟨G-system⟩ : f31^∞`;
2. check `f37 ∈ √J`  (equivalently `1 ∈ J : f37^∞`, or `V(J)=∅`).

If `V(J)=∅`, then no solution has `f31≠0`, so the whole f37 branch is spurious
and the case reduces to f31 alone. This is a saturation + radical-membership
(or emptiness) test on the *small* G-system in `(d2,d1,d0,e,r,s,d₋₄)`, **not**
on the 618-term `f37`. `T5_NOTES.md` records that monolithic Gröbner on the
618-term window system times out, and `FIELD_SPLIT_AUDIT.md §6` warns off
materializing the global ideal — but the *saturation of the compact G-system*
is a much smaller object than any f37 computation and has not been tried. A
cheaper first pass: the *local/numeric* version — take one of the 21 verified
solutions, and Newton-search near it for a nearby G-system solution with
`f37=0, f31≠0`; a robust failure across strata is strong `[N]` corroboration
before committing to the exact saturation.

The free-family paper already isolates the restored small system on `d2=d1=0`
(`H2,H3,H5`); the general test is the same three equations kept with their
`d2,d1` terms, saturated instead of specialized.

---

## 5. Recommended attack order

Ranked, cheapest decisive computation first for each option.

### Option A (recommended primary) — pre-resultant "f37 is spurious"
Prove `V(⟨G-system⟩ : f31^∞) = ∅`, i.e. the original system forces `f31=0`,
collapsing the whole f37 branch at once.
- *Cheapest decisive computation:* saturate the compact G-system by `f31` in
  `(d2,d1,d0,e,r,s,d₋₄)` and test emptiness / `f37 ∈ √J` (§4). Small-variable,
  small-degree relative to any f37 object.
- *Cheap pre-screen `[N]`:* Newton-search from the 21 known f31-solutions for a
  nearby `f37=0, f31≠0` G-solution; expect failure.
- *Why first:* if it succeeds it closes f37 in one stroke; it is directly
  motivated by the two proven facts (free family doesn't lift; numeric
  solutions never hit f37); and it sidesteps the two structural blockers of the
  cone route (§3: empty terminal pruning, live free family).

### Option B (fallback / hybrid) — cone engine on `(d2,d1)≠0` + pre-resultant on the free family
Extend `CASCADE_ENGINE_PLAN.md`'s engine to the f37 signature (§3) for the
generic locus, and dispatch the free family separately to the pre-resultant
kill (already done, `F37_FREE_FAMILY_SYSTEM.md`).
- *Cheapest decisive computation:* run the head-only slope-3 block (levels 0–4,
  `a=0`, 120 equations in `h_0…h_4`, degree ≤88) and test for a common `q`-adic
  vanishing / residue contradiction (`T5_F37_GRADED.md §11.3`).
- *Why second:* the cascade exists and is sparse (proven), but §8 shows
  terminal starvation kills **zero** strata and §6.3 shows the identity alone
  can never finish — so this route *cannot* be self-contained; it must still
  import a pre-resultant fact for the free family. More engineering, more
  residue systems, larger caps than f31.

### Option C (not recommended alone) — cone engine as the sole tool
Provably insufficient: `T5_F37_GRADED.md §6.3` establishes that no argument
from the f37 identity + window caps can exclude the free family. Pursue only as
the `(d2,d1)≠0` half of Option B.

**Bottom line ordering:** A ≫ B > C. The strongest lever is that f37 has never
been hit by a real solution and its one known exact family is a proven
non-lift — so the highest-value cheap computation is the G-system saturation
that would show f37 is entirely an artifact of the resultant.

---

## Key findings (summary)

1. **f37 IS weighted-homogeneous** (weight 134, `(2,3,4,5,17)`; [V] this
   session, 618 terms one weight). The "non-uniform grading" is a
   *Newton-polygon* fact — the `d₋₁`-power vector `(18,15,12,9,6,4,2,0)` has
   two slopes (−3,−2), so f31's single `w=Φ/d₋₁³` root reformulation fails; a
   two-edge normal form / non-uniform cascade replaces it.
2. **A sparse-companion cascade exists** (Lemma A′, proven), so the cone engine
   *can* ingest f37 — but its terminal degree-starvation prunes **zero** strata
   (caps 2× f31, `ê²` vs `ê³`), and the free family is a permanently-live
   branch the identity alone cannot kill.
3. **"All of f37 fails to lift" is plausible** and is the best route: numeric
   solutions never satisfy f37, and the one proven exact family (`d2=d1=0`)
   provably does not lift. Decisive test = saturate the compact pre-resultant
   G-system by `f31` and check the variety `{f31≠0}` is empty.
4. Bonus verified fact not previously in the notes: **`h37` is irreducible over
   `Q`** (145 terms, degree 19).

File created: `C:\Users\wstri\dev\math-stuff\d2_plane_72_108\F37_FRONTIER.md`.
