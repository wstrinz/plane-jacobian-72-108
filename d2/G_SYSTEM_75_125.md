# The (75,125) D-transform G-system: built, with a characterised a=3 boundary

## Verdict

**BUILT — the structure transfers; the physical-weight normalisation does not.**
Phase 2 of the TRANSFER TEST constructs the `(75,125)` analogue of the
`(72,108)` pre-resultant **G-system** (`FULL_SYSTEM_BRIDGE.md` §1,
`regenerate_system.py`) from phase 1's tower (`C_SERIES_75_125.md`). Applying the
parametric recipe — linear window `S^a`, forcing window `S^b`, with `(a,b,t) =
(3,5,5)` — yields a genuine, weight-homogeneous G-system:

```
10 quintic generators   G1..G9, G11    (forcing slices j = 1..11, skip j = 10)
 9 spare unknowns        dm2..dm10  =  d_-2 .. d_-10   (= d_-2 .. d_-(a-1)t)
u-grading weights (AP)   26,27,28,29,30,31,32,33,34,36   (Phi at 36; 35 = skipped G10)
```

against `(72,108)`'s 4 cubic generators, 3 spares, and weights `13,14,15,17`
(physically `156,168,180,204`). The construction is `g_system_75_125.py`; the
canonical machine-readable system is `g_system_75_125.json`; the exact checker is
`g_system_75_125_verify.py` (**all checks pass, exit 0**), whose primary control
rebuilds the `(72,108)` system from scratch and reproduces its published
generators `G1,G2,G3,(G5body+Phi)` **bit-for-bit**.

**The one genuine obstruction, characterised precisely.** The physical
(y-valuation) weight per window step,

```
W_step := ord_y(Phi) / M ,        M = b·t + jphi = 36  (the forcing slice),
```

is **`12` (integer)** at `(72,108)` (`204/17`) but **`201/36 = 67/12`
(NON-integral)** at `(75,125)`. [CORRECTED 2026-07-24: this is NOT an
`a ≥ 3` boundary — see the correction block below.] This was attributed to
`CORNER_144_COMPARISON.md` §5 already flagged as a **quasipolynomial** window cap
(`8w + ceil(w/5)` for `(108,144)`, denominator `5`; here the denominator is
`12`). The abstract G-system is weight-homogeneous under the intrinsic
**u-grading** (integer weights `26..36`), but the exact integer `y^W` **stripping**
that made the `(72,108)` bridge substitution valid (`FULL_SYSTEM_BRIDGE.md` §3)
does **not** carry over: the window variables have **quasi-affine**, not affine,
`y`-order. The weights are "arithmetic-progression-like" in the u-grading; the
physical normalisation is rational.

---

## 1. The recipe (parametric, from (72,108) + phase 1)

Write the normalised C-series in the tower parameter `u`:

```
S = Σ_m d_m u^(t-m),     d_t = 1,   d_{t-1} = 0  (x-shift),
    d_m = c_m · c^(a(t-m)-1)   (D-transform;  c = leading poly C).
```

`P = C^a` occupies the **linear window** `S^a`; the forcing term `F`, hence `Phi`,
occupies the **forcing window** `S^b` (phase 1, `C_SERIES_75_125.md` §4;
verified there that `(72,108)`/`(108,144)` use `S^2,S^3` / `S^3,S^4`). Two
families of slices:

- **Linear window** `La(k) := [u^(a t + k)] S^a`. The deepest new coefficient in
  `La(k)` is `dm_{(a-1)t+k}`, appearing linearly with coefficient `a` (from
  `d_t^{a-1}=1`). Solve `La(k)=0` sequentially to **eliminate** the deep window
  unknowns.
- **Forcing window** `Lb(j) := [u^(b t + j)] S^b`. After the linear
  substitutions these are the **generators** `G_j`. `Phi = f·C^N` is the
  `(D~^b)_{-jphi}` slice, `M = b t + jphi`, `jphi = -s = a t − κ − 1`.

For `(72,108)` (`a,b,t = 2,3,4`) this is exactly `regenerate_system.py`: linear
`S^2` slices `D2(k)`, forcing `S^3` slices `D3(j)`, cubic `G1,G2,G3,(G5body+Phi)`.
The builder reproduces those four polynomials **exactly** (checker §A).

---

## 2. What the (75,125) parameters make of the recipe

`(a,b,t,κ,q) = (3,5,5,3,2)`, `e = b−a+1 = 3`, `s = κ+1−a t = −11`,
`jphi = 11`, `M = b t + jphi = 36`, deepest window unknown `d_-31`.

| quantity | (72,108) `a=2` | (75,125) `a=3` | rule |
|---|---:|---:|---|
| linear window | `S^2` | `S^3` | `S^a` |
| forcing window | `S^3` (cubic) | `S^5` (quintic) | `S^b` |
| forcing slice `M` | `17` | `36` | `b t + jphi` |
| linear eliminations | `8` (k=1..9 skip 8) | `20` (k=1..21 skip 20) | `(b−a)t + jphi`, skip `(b−a)t + jphi−1` |
| **# generators** | `4` (j=1,2,3,5) | `10` (j=1..9,11) | `jphi − 1`, skip `jphi−1` |
| **# spare unknowns** | `3` (dm2,dm3,dm4) | `9` (dm2..dm10) | `(a−1)t − 1` |
| spare window | `d_-2..d_-4` | `d_-2..d_-10` | `d_-2 .. d_-(a-1)t` |

The **spare inventory jumps 3 → 9** purely because `a` goes `2 → 3`: the linear
window `S^a` only reaches `dm_{(a-1)t+k}` (its deepest linear coefficient), so the
un-eliminated spares are `dm_1..dm_{(a-1)t}`, i.e. `d_-1..d_-(a-1)t`; removing the
state's `e = dm1` leaves `(a-1)t − 1` spares. The **generator count rises 4 → 10**
because `Phi` sits at the deeper slice `jphi = 11` (vs `5`), and the forcing
slices below it (skipping `jphi−1`, whose window unknown is left undetermined,
matched to the skipped linear slice) all become generators.

The skip pairs `(forcing j = jphi−1) ↔ (linear k = (b−a)t + jphi−1)` share the
single undetermined window unknown `d_-((b-1)t+jphi-1)` (`d_-12` at `(72,108)`,
`d_-30` here), keeping the system closed in the spares — exactly as
`regenerate_system.py` drops `D2(8)`/`D3(4)` and their shared `dm12`.

Every generator is verified **u-grading homogeneous**: each `G_j` has all
monomials at weight `b t + j` (checker §C, re-parsing the stored strings
independently of the builder). Because `d_t = 1` is constant and the `a=3`
substitutions are cubic, the generators are **weight-homogeneous but NOT
total-degree homogeneous** — total-degree maxima grow `5,5,5,5,5,5,6,6,7,8` across
`G1..G9,G11` (at `(72,108)` the analogue `G1` already mixes degrees 2 and 3).

---

## 3. Weights, and the physical-weight obstruction

The generators are weight-homogeneous under the **intrinsic u-grading**
`w(d_m) = t − m`, `w(Phi) = M`. This grading is automatic (every monomial of a
`u`-slice shares its `u`-power) and is preserved by the linear substitutions.
Under it the forcing generators carry the arithmetic progression

```
w(G_j) = b t + j :   26, 27, 28, 29, 30, 31, 32, 33, 34, 36   (common diff 1),
```

with `w(Phi) = 36`. This is the "arithmetic-progression-like family" the mission
asked for; the value `35` is absent — it is the skipped generator `G10`.

The **physical** (y-order) weight is `W_step × (u-weight)`, with

```
W_step = ord_y(Phi) / M .
```

| case | `a` | `ord_y(Phi)` | `M` | `W_step` | integral? | physical G-weights |
|---|---:|---:|---:|---|---|---|
| `(72,108)` | 2 | `204` | `17` | `12` | **yes** | `156,168,180,204` |
| `(108,144)` | 3 | `205` | `25` | `41/5` | no | — (denominator 5) |
| **`(75,125)`** | **3** | **`201`** | **`36`** | **`67/12`** | **no** | — (denominator 12) |

At `(72,108)`, `W_step = 12` is an integer, so the u-grading scaled by 12 **is**
the `y`-valuation grading: the full-window generators satisfy
`G_i(V_full) = y^{W_i}·G_i(V_stripped)` exactly, which is what licenses the
bridge's stripped ansatz substitution (`FULL_SYSTEM_BRIDGE.md` §3). At
`(75,125)`, `W_step = 201/36 = 67/12` is **not** an integer: `M = 36` does **not**
divide `ord_y(Phi) = 201`. The window variables' `y`-order is therefore
**quasi-affine** (nearest-lattice, quasi-period `12`), the forcing generators are
homogeneous **only** in the u-grading, and the exact integer `y^W` stripping does
not transfer. This is precisely the `a ≥ 3` warning of
`CORNER_144_COMPARISON.md` §5 (quasipolynomial lower cap `8w + ceil(w/5)`,
denominator `5` = denominator of `(108,144)`'s `W_step = 41/5`); here the quasi
period is `12`, the denominator of `67/12`.

**This is the boundary, and it is a result.** The G-system as an algebraic ideal
(generators, spares, u-grading) transfers cleanly; the *window-cap layer* — the
piece that turns the ideal into a bridge with bounded polynomial spares — is
quasipolynomial, not the affine `12k/15k/14k` of `(72,108)`.

---

## 4. Φ-consistency

**Verdict: CONSISTENT (intrinsic grading).** In the forcing window `Phi` enters
the deepest generator `G11 = [u^36]S^5 + Phi` with the homogeneity-forced weight
`w(Phi) = M = 36`. This is the same slice phase 1 identified: `C_SERIES_75_125.md`
builds `Phi = f·C^98` as the `u^36` slice of `S^5`, with the slice-sum invariant
`clear = a·M − b = 103` giving `N = clear − b = 98`. The physical divisor
`Phi = −(1/9) y^201 (y^3+1)^101` has `ord_y = 201 = W_step·M = (67/12)·36` — the
Φ-consistency identity holds exactly as a rational statement. The **only**
departure from `(72,108)`'s master-identity check (`FULL_SYSTEM_BRIDGE.md` V3:
`Phi_full/y^204 = c·t^30·q`, an *integer* strip) is that here the strip power
`W_step·M = 201` is fine but the *per-step* normalisation `W_step` is
non-integral — the §3 obstruction, and nothing more.

So `Phi` plays the analogous role: it sits at the forcing slice with the weight
the G-system's homogeneity forces, and it reproduces phase 1's `N = 98`.

---

## 5. Controls

`g_system_75_125_verify.py` (exact sympy, `--quiet`, exit 0):

- **§A recipe control (independent ground truth).** Rebuild `(72,108)` from
  scratch: reproduces the published `G1,G2,G3,(G5body+Phi)` of
  `FULL_SYSTEM_BRIDGE.md` §1 **exactly**; spares `dm2,dm3,dm4`; u-weights
  `13,14,15,17`; `W_step = 12`; physical weights `[156,168,180,204]`. This
  validates the builder against externally-known values before it is trusted on
  `(75,125)`.
- **§B/§C JSON structure + homogeneity.** Corner signature, ring/variable order,
  9-spare inventory, 10-generator list, and (re-parsing the stored strings,
  builder-independently) u-grading homogeneity of every generator, weights the
  AP `26..36`.
- **§D obstruction.** `W_step = 201/36 = 67/12` non-integral; contrast with the
  integral `(72,108)` and the non-integral `(108,144)` `41/5`.
- **§E spot re-derivation (anti-fabrication).** Rebuild `(75,125)` `G1,G2` from
  scratch (truncated) and match the stored JSON generators.
- **§F Φ-consistency.** `M = 36`, `clear = 103`, `N = 98`, `ord_y(Phi) = 201 =
  W_step·M`.

---

## 6. What this instantiates for the case compiler

`g_system_75_125.json` (schema `g-system-v1`) fills exactly the fields
`case_compiler.py`'s `pre_resultant_G_system` block marks *schematic* for
`F2_j1_75_125`:

| dossier field | schematic → instantiated |
|---|---|
| `generators` | `null` → `["G1".."G9","G11"]` (10 quintic; full polynomials in the JSON) |
| `ring` | generic → `Q[d3,d2,d1,d0,dm1,dm2..dm10,Phi]` |
| `spare_dictionary` | generic → `dm2..dm10 = d_-2..d_-10` (9 spares) |
| `G_weights` | `null` → u-grading AP `[26..36]` (physical `67/12 × ·`, non-integral) |
| `window_caps` | "case arithmetic" → **OBSTRUCTED**: quasipolynomial, quasi-period 12 |
| `bridge` | schematic → 10 quintic generators, 9 spare unknowns + Phi over the state |

The honest status flip is: **structure INSTANTIATED, window-cap layer OBSTRUCTED
(a=3 boundary).**

---

## 7. `[judgment]` list — where this is conditional

1. **[inherited from phase 1]** The corner `(5,20)→(7/5,2)`, `(a,b,t,κ,q) =
   (3,5,5,3,2)`, `C = y^2(y^3+1)`, `Phi = −(1/9)y^201(y^3+1)^101` are phase 1's
   built objects (`C_SERIES_75_125.md`), themselves conditional on the standard
   unreduced-polygon reduction (`PHI_75_125.md` judgment 2). Unchanged here.
2. **[derived]** The G-system recipe (linear `S^a`, forcing `S^b`, slice indices,
   spare inventory, skip pairing) is derived from the `(72,108)` construction and
   validated by exact reproduction of its published generators — not assumed.
3. **[derived, this lane]** The `W_step` non-integrality (`67/12`) and the
   consequent quasipolynomial window-cap boundary are computed from the built
   `Phi` and slice `M`, and cross-checked against the `(108,144)` `a=3` datum
   (`41/5`, `CORNER_144_COMPARISON.md` §5). The abstract u-graded G-system is
   exact; the affine window caps of `(72,108)` provably do **not** exist here.

**Scope note.** This lane builds the G-system and characterises its weight
grading and the window-cap boundary. Turning the ideal into a running *bridge*
(the ~122-equation-analogue with bounded polynomial spares) would require the
quasipolynomial window caps to be pinned exactly — the pending item this boundary
result now scopes precisely.

---

## Files

- `G_SYSTEM_75_125.md` — this writeup.
- `g_system_75_125.py` — the parametric construction (linear/forcing windows,
  elimination, weights, obstruction); emits the JSON; run end-to-end (~90 s for
  the full `S^5` build). Reproduces the `(72,108)` control inline.
- `g_system_75_125.json` — canonical machine-readable G-system (documented
  variable order, generators, weights, spare inventory, obstruction) for
  `case_compiler.py`.
- `g_system_75_125_verify.py` — exact PASS/FAIL checker (`--quiet`, exit 0);
  control = `(72,108)` recipe reproducing the known G-weights and generators.


## CORRECTION (2026-07-24, GPT-Pro review 6; verified in f2_family_verify.py)

The non-integral W_step is **not** an `a ≥ 3` phase transition: already at
a=2 (F2 j=0, the (50,75) case) W_step = 75/21 = 25/7 is non-integral. The
meaningful invariant is the **window denominator**
q_window := denom(ord_y(Phi)/M), which for the F2 family equals exactly
**5a−3** (a=2: 7, a=3: 12, a=4: 17, ...; gcd(2a−1,5a−3)=1 always). The
integral W_step = 12 at (72,108) is a friendly coincidence of that corner,
not generic a=2 behavior. Consequence: the right generalization is an
engine for **rationally related gradings** (bigraded/two-coordinate window
lattice), not a separate `a ≥ 3` engine; for the fixed (75,125) case the
tactical route is a q_window=12-phase (period-12) window compiler.
