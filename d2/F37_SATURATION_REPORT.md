# F37_SATURATION — the entire f37 branch is a resultant artifact (2026-07-22)

**Scope:** the f37 ≡ 0 branch of the master identity `f31 · f37 · d₋₁²¹ ≡ 0`
(`STATE.md` item 5).  The task was the "no-lift" saturation test of
`F37_FRONTIER.md` §4–5: show that no solution of the original (pre-resultant)
system lies on `{f37 = 0} \ {f31 = 0}`.

**Result (FULL, positive, exact).** The intended saturation is not merely
`(I_pre + ⟨f37⟩) : f31^∞ = (1)` — the stronger statement holds:

> **`f31` lies in the pre-resultant ideal itself:**
> `f31 ∈ ⟨G1, G2, G3, G5body + Φ⟩` over `Q[d̃2,d̃1,d̃0,d₋₁,d₋₂,d₋₃,d₋₄,Φ]`.

Consequently `f31` vanishes on the **entire** pre-resultant variety, over every
field and every specialization of `Φ` (including the genuine (72,108)
instance).  So every solution of the original system has `f31 = 0`; the whole
locus `{f37 = 0} \ {f31 = 0}` contains no solution.  The f37 factor — and the
`d₋₁²¹` factor — of the resultant are **classical resultant excess factors**,
not geometry.

This closes the f37 branch in one stroke (Option A of `F37_FRONTIER.md` §5),
superseding both the single proven `d̃2=d̃1=0` free-family slice
(`F37_FREE_FAMILY_SYSTEM.md`) and any stratum-by-stratum program on the f37 cone.

---

## 1. What was computed

The pre-resultant system is tiny (`STATE.md` item 4; regenerated from
`t4_state.pkl`), all polynomial over `Q`:

```
G1 = 3·d̃1·e² + 6·d̃2·e·r + 6·e·d₋₄ + 6·r·s
G2 = -3·d̃0·e² + 3·d̃2·r² + 6·r·d₋₄ + 3·s²
G3 = -6·d̃0·e·r - 3·d̃1·r² - e³ + 6·s·d₋₄
G5 = 2·Φ - 6·d̃0·e·d₋₄ - 6·d̃0·r·s - 6·d̃1·r·d₋₄ - 3·d̃1·s² - 6·d̃2·s·d₋₄ - 3·e²·s - 3·e·r²
```
(here `e=d₋₁, r=d₋₂, s=d₋₃`, denominator-cleared to integer coefficients).

Let `I = ⟨G1,G2,G3,G5⟩`.  Over `Q`, Gröbner (Singular 4.2.1, `dp`):

| # | statement | result |
|---|-----------|--------|
| 1 | `reduce(f31, std I)` | **0**  → `f31 ∈ I` |
| 2 | `f37 ∈ I` ? | **false** (f37 not forced — ideal is proper) |
| 3 | `d₋₁ ∈ I` ? | **false** |
| 4 | `reduce(f31·f37·d₋₁²¹, std I)` | **0** (master identity ∈ I, as expected) |
| 5 | `E := I ∩ Q[d̃2,d̃1,d̃0,d₋₁,Φ]` | **principal**, 1 generator, deg 31, 102 terms |
| 6 | `f31 ∈ E`, `E[1] ∈ ⟨f31⟩` | both **0** → **`E = ⟨f31⟩`** |

So the true elimination ideal of the pre-resultant system in the
`(d̃2,d̃1,d̃0,d₋₁,Φ)` variables is **exactly `⟨f31⟩`**.  The resultant returned
`f31·f37·d₋₁²¹`, a *multiple* of the generator; `f37` and `d₋₁²¹` are the
spurious factors that a resultant (as opposed to a Gröbner elimination) is known
to introduce.

**Independent cross-check (Part 2 of the script).** The same conclusion is
re-derived from the H-system `⟨H2,H3,H5⟩` (the actual resultant path of
`regenerate_system.py`, with `d₋₄` pre-eliminated), saturated by `d₋₁`:
`f31 ∈ ⟨H2,H3,H5⟩ : d₋₁^∞` and its elimination ideal is again `⟨f31⟩`
(1 gen / deg 31 / 102 terms).

**Explicit membership certificate.** Singular `lift()` produced cofactors
`c1..c4` (degrees 27–28, 3418–5848 terms, rational coefficients) with
`f31 = c1·G1 + c2·G2 + c3·G3 + c4·(G5body+Φ)`.  This identity is re-checked
**independently in sympy**, not trusting the Gröbner engine, by expanding the
combination and asserting it equals `f31` exactly (`f37_sat_verify.py`, PASS).

---

## 2. Why this settles the branch (logic)

On `V(I)`, the master identity gives `f31·f37·d₋₁²¹ = 0`.  A priori a solution
with `f31 ≠ 0` could sit on `{f37 = 0}` (or `{d₋₁ = 0}`).  Fact 1 removes the
possibility outright: `f31 ∈ I` forces `f31 = 0` at **every** point of `V(I)`.
Hence:

```
V(pre-resultant system)  ⊆  {f31 = 0},
```
so `V(pre-resultant) ∩ {f37 = 0} ⊆ {f31 = 0}` (the requested containment) holds
trivially, and more.  No `d₋₁ = 0` side condition is needed for the conclusion,
because `f31` vanishes identically on the variety; the separate proof that
`d₋₁ ≡ 0` is impossible (`STATE.md` item 5) is not even required here.

This matches, and explains, every prior signal:
- the 21 verified numeric solutions all had `f31 = 0`, never `f37 = 0`
  (`F37_FRONTIER.md` §1.8) — because the variety *is* `{f31 = 0}`;
- the one exact f37 family `d̃2=d̃1=0` provably did not lift
  (`F37_FREE_FAMILY_SYSTEM.md`) — a special case of "nothing on f37 lifts";
- no argument from the f37 identity + windows alone could close the branch
  (`T5_F37_GRADED.md` §6.3) — correct, because the obstruction is not intrinsic
  to f37 at all; it is that f37 is not in the ideal.

---

## 3. Strata of {f37=0} now proven not to lift

**All of them.**  The result is uniform over the whole f37 cone, not a stratum
list: every point of `{f37 = 0} \ {f31 = 0}` fails to solve the pre-resultant
system.  In particular this subsumes the previously-open strata of the f37
split-place ledger and the "permanently-live free family" branch of
`T5_F37_GRADED.md` §6.3 / `F37_FRONTIER.md` §3 — none of them lift, because none
of `{f37=0}\{f31=0}` lifts.  The (72,108) case therefore reduces to the `f31`
branch alone.

---

## 4. Computational cost / walls

None hit.  The decisive object is the pre-resultant ideal (4 generators, ≤ 12
terms each, 8 variables), far smaller than any f37 window object.  Every
Singular run (reductions, the `eliminate` for `E`, the H-system `sat`, the
`lift`) finishes in **seconds**; the sympy certificate check finishes in about a
minute.  The 45-minute cap was never approached.  The `LINALG_CERT_REPORT.md` /
`MSOLVE_REPORT.md` blow-ups were on the *post-resultant* window systems (f31/f37
as ~10⁵-term dense objects); the pre-resultant ideal never has that size,
which is exactly why the saturation route recommended in `F37_FRONTIER.md` §4 is
cheap.

Note: the *literal* saturation `(I+⟨f37⟩):f31^∞` was not needed and was not run
— computing the elimination ideal `E` and observing `E = ⟨f31⟩` (equivalently
`reduce(f31, std I) = 0`) is strictly cheaper and gives the stronger statement.

---

## 5. Logical scope and caveats

- The computation is over `Q` with `Φ` a **free indeterminate**; the genuine
  (72,108) instance is the specialization `Φ = f₁·C₄²⁸ ∈ Q[y]`, `d̃k ∈ Q[y]`, so
  it inherits `f31 = 0` from the ideal-membership identity (membership is
  preserved under any ring homomorphism `Q[Φ] → Q[y]`).
- The pre-resultant generators are exactly `STATE.md` item 4's
  `(D̃³)₋₁,₋₂,₋₃` and `(D̃³)₋₅ + Φ` after the `(D̃²)` linear substitutions —
  the same object `regenerate_system.py` builds and the T6 audit validated.
  The higher `(D̃²)/(D̃³)` slices only define `d₋₅…` and `λ` and do not touch
  `(d̃2,d̃1,d̃0,d₋₁,Φ)` (`T6_SELECTION_AUDIT.md`), so `E = ⟨f31⟩` is the complete
  relation among the f31/f37 variables.
- This is a statement about the *pre-resultant polynomial system*, i.e. the
  reduction the derivation actually produces.  It does not by itself prove the
  (72,108) case infeasible — it removes the f37 branch, leaving the `f31`
  window/cascade program (STATE.md, ongoing) as the sole remaining task.

---

## Files created (no existing files modified; nothing committed)

- `F37_SATURATION_REPORT.md` — this report.
- `f37_sat_confirm.sing` — self-contained reproducible Singular script; PART 1
  (G-system) facts [1]–[7], PART 2 (independent H-system) [8]–[10].
- `f37_sat_verify.py` — independent sympy verification of the explicit
  membership certificate `f31 = Σ cᵢ·Gᵢ` (does not trust the Gröbner engine),
  plus the master-identity excess-factor check.  PASS.
- `f37_sat_certificate.txt` — the four `lift()` cofactors `c1..c4` (data for the
  sympy checker).

### Reproduce

```
# sympy certificate (self-contained, no Singular):
python f37_sat_verify.py

# Gröbner facts (Singular; on this box pipe via WSL since /mnt/c is unavailable):
cat f37_sat_confirm.sing | wsl.exe -d Ubuntu -- bash -lc 'cd $HOME && Singular -q'
```
