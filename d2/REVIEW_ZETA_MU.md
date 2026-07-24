# REVIEW_ZETA_MU.md — adversarial review of ZETA_TAIL.md and the μ-ladder

Independent skeptical review pass (separately prompted, no code shared with the derivation), 2026-07-23. Targets: ZETA_TAIL.md / zeta_tail.py /
zeta_tail_verify.py; upstream COMPOSITE_CHARTS.md (incl. correction block),
PHI_F7.md, PHI_F14.md, PHI_CORNER4.md; the case compiler's F12 flag.
Companion computations: `review_zeta_mu.py` — **31/31 confirmations hold,
exit 0** (independent re-derivations; the F12 μ=2 quartic is re-obtained by a
*different* elimination route than the author's). Two of my own first-draft
test expectations were wrong and fixed (documented in §D below) — neither was
a defect in the reviewed work.

## Verdicts

### 1. One-slice rigidity theorem — **SOUND** (within its declared model), one cosmetic mislabel

The top slice of `[x^ζ(C^a + tail), C^b + tail + F]` at x-degree
`ζ+(a+b)t−1` receives exactly one contribution — the **head-head** bracket
`[x^ζC^a, C^b] = ζb·x^(ζ+(a+b)t−1)c^(a+b−1)c′` — because tail powers are
strictly lower (max cross index a+b−1, a full `t` below the slice) and an
F-cross-term would need `v(F) ≥ bt > 0` against `v(F) < 0`. I verified the
slice formula and its uniqueness at three (a,b) instances × two t values with
fully generic scalar tails (checks A), beyond the author's single instance.
The argument is complete *within* the model `P = x^ζ(C-series)` declared in
[judgment 1]; that model form is the real conditional boundary and is
honestly flagged. **Cosmetic:** both ZETA_TAIL.md §1 and zeta_tail.py call
this the "tail-tail bracket's top slice" — it is the head-head slice (the
tails contribute nothing there; that is the point). Wording should be fixed
to avoid confusing a future reader.

### 2. η=−1 universal death lemma — **SOUND**, scope matches use

`K = eT ⟺ T = t−1` is verified symbolically in both directions; the
integration identity `(f/c^e)′ = 1/(aTc)` holds (checked at e=5 and e=12,
the two invocation values). The residue argument does **not** need c
squarefree — it enumerates multiplicity configurations, and I reproduced
every residue independently via `sympy.residue`: any simple g-root has
residue `1/c′(α) ≠ 0`; the [2,1] tuning at q=5 kills Res_α only at
β = 6α/5, where Res_β = 78125/(7776·α⁷) ≠ 0; the [3] and [2] full-multiplicity
residues are 15/α⁷ and −7/α⁸. The lemma is claimed and used only at
dg ∈ {2,3} (F12, F13) — scope and invocation agree. (The y=0 residue never
needs checking: one nonvanishing residue suffices for the obstruction.)

### 3. μ-ladder specialization — **SOUND with two findings** (one presentational, one scope GAP)

- `r = dg−1` is **definitional** (r = a0−q−1, dg = a0−q), so "specializes
  exactly via the identity r = dg−1" is algebra, not evidence — verified
  symbolically both ways (μ=1 → old law; μ=dg → PHI_F7 law). No smuggling,
  but the phrasing dresses a triviality as a discovery.
- **Redundancy finding:** `deg − ord = (e+N)·dg + gap` identically, so the
  cof-law is *equivalent* to the mult-law given deg/ord (verified check C).
  Each realized rung therefore contributes exactly ONE new datum
  (mult = μ(e+N)−(μ−1)), not two. The law's presentation as four formulas
  overstates its per-rung content; worth one sentence in ZETA_TAIL.md.
- **GAP (scope overreach):** "Parity refinement: dg even kills μ=1
  (PHI_F7's theorem)". PHI_F7's own judgment 5 states branch-completeness is
  proven **only at dg=2**; at dg=4 (F10, F15) the μ=1 exclusion is
  conjecture, not theorem. COMPOSITE_CHARTS' new "dg-parity mechanism"
  paragraph inherits the same overreach ("at even dg the eliminant forces
  double or complex roots" — established at dg=2 only). Both sentences
  should be scoped "at dg=2; conjectured for higher even dg" until the F10
  rung enumeration (in flight) settles it. Not load-bearing for any kill.

### 4. F12 μ=2 rung — **CONFIRMED** by independent elimination

Rebuilt the branch system from the F12 ODE parameters; the residual quotient
has degree 4 (the naive top coefficient cancels because deg f = 38 = d_res —
T·38 = K·8 = 152), giving 5 equations affine-linear in the four u-coefficients.
Eliminating via **augmented-matrix maximal minors + gcd** (a different route
from the author's solve-then-substitute) reproduces exactly
`195β⁴+120β³−40β²+32β−80` up to a nonzero rational scalar and unit factors at
the excluded positions {0, −1}. The quartic has exactly 2 real roots, 0 and
−1 are not roots, the full 5-equation system vanishes **exactly** mod the
quartic, u(−1)·u(0)·lc(u) ≠ 0 mod the quartic (orders exact), and the
signature (814, 506, 203, 105) is exact bookkeeping. The μ=3 rung and its
u = −(2048y⁴+2560y³+320y²−80y+35)/1155 were also independently rebuilt and
match (checks D, E).

### 5. Cross-document consistency — **mostly consistent; three minor findings**

- **COMPOSITE_CHARTS correction block: accurate.** It faithfully reports what
  changed, quotes the correct quartic and signatures, and leaves COMPOSITE_CHARTS.md's
  μ=1 result standing (which is right). Its dg-parity paragraph shares the
  item-3 scope overreach (above).
- **Compiler F12 flag: consistent but incomplete.** It cites the η=0 rungs
  and the η=+2 collapse signature (1292,806,162,324), but omits (a) the
  η=+2 μ=3 rung (1292,806,484,2) and (b) the arithmetically-viable but
  unmotivated η ∈ {−2,−3} rows that ZETA_TAIL judgment 2 flags OPEN. The
  closing "model selection needs the actual polygon reduction" covers (b)
  loosely; adding "η=−2,−3 viable-unmotivated, OPEN" would make the flag
  faithful to the doc.
- **Latent sweep-logic gap (code robustness, not a current error):** in
  zeta_tail.py's sweep, a partition whose rung evades the forced-order bound
  by sub-resonance (`sub` integral and < om) is silently dropped from
  `min_forced` instead of blocking the DEAD verdict — if it ever fired, a
  DEAD row could be unsound. I verified programmatically it fires on **no**
  swept (family, η, partition) row (check F), so all current verdicts stand;
  the loop should still be patched to mark such rows "viable (sub-resonant
  evasion possible)" defensively.

## Overall assessment

**Needs the listed patches; no real hole.** The mathematics of the rigidity
theorem, the η=−1 lemma, the corrected family, and all concrete branch
computations survive adversarial re-derivation, including by an independent
elimination route. The required patches are: (i) scope the two "dg even kills
μ=1" sentences to dg=2, (ii) fix the "tail-tail" label, (iii) one sentence on
the cof-law redundancy, (iv) complete the compiler F12 flag, (v) the
defensive sweep-loop guard. Items (i) and (iv) matter for external readers;
the rest are hygiene.

## Files

`REVIEW_ZETA_MU.md` (this file), `review_zeta_mu.py` (31 confirmations,
exit 0, --quiet). No existing file touched; nothing committed.
