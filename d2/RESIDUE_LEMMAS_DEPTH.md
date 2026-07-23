# Depth-quantified residue lemmas — the a-parametric jet structure

**Date:** 2026-07-22
**Scope:** the surviving f31 subcase-(1) frontier
`cascade_cones_sub1_qt_inf_rl.json` (post-infinity, residue-kills fed back):
**171 branches / 1145 flag cases**, survivors at `a ∈ [2,10]`. This extends
`RESIDUE_LEMMAS.md` from the depth-1 initial forms (P6/P10/P11, the 23 tied
supports C01–C23) to the **a-quantified depth structure** — the jet tower that
the affine t-depth law `depth = 30−3a` demands.

**How to read this file.** Statements tagged **[data]** are mechanical
consequences of `cascade_cones_sub1_qt_inf_rl.json`, re-derivable by
`residue_lemmas_depth_verify.py`, which also re-parses every residue
coefficient from `f31_graded.txt` (via `cascade_signature.load_levels`; no
`h_l` coefficient is copied into the checker or this document). Statements
tagged **[judgment]** are assessments not forced by the data. Every jet
equation displayed below is reproduced symbolically by the checker from the
source tie polynomial — it is not hand-transcribed algebra.

---

## 0. The jet structure (RESIDUE_LEMMAS.md §1, made explicit)

Fix a place `p` with uniformizer `π` and a tied set `T ⊆ h_l` whose monomials
share one valuation-weight `m`. Expand each signature variable in local jets

```
d2 = π^k (D + D₁π + D₂π² + …),   d1 = π^x (X + X₁π + …),
σ  = π^z (S + S₁π + …),          e  = π^b (E + E₁π + …).
```

Because every monomial of `T` has weight `m`, the factor `π^m` divides out
uniformly and the tie polynomial `F = Σ_{M∈T} A_M M` expands as

```
F = π^m ( C₀ + C₁ π + C₂ π² + … ).
```

A **depth-δ** obligation (`monomial_tie_rise` of order δ, i.e. `v_p(F) ≥ m+δ`)
is exactly the system `C₀ = C₁ = … = C_{δ-1} = 0`. The checker proves, for each
pattern below, the exact identities

```
C₀ = IF(D,X,S,E)                         (the depth-1 initial form)
C₁ = ∇IF · (D₁,X₁,S₁,E₁)
C₂ = ∇IF · (D₂,X₂,S₂,E₂) + ½ (jet₁)ᵀ Hess(IF) (jet₁)
```

and in general `C_δ = ∇IF · (order-δ jets) + P_δ(lower jets)`. **The linear
coefficient of the order-δ jet is the fixed gradient `∇IF`, independent of δ.**

**The decisive structural fact [judgment, proved for the occurring supports in
§3–§4].** At the surviving t-place `t = y+1` (a *finite* place — `q(−1)=3315≠0`,
`t5_multiplace_verify.py`), the residues `(D,X,S,E)` and every jet
`(Dⱼ,Xⱼ,Sⱼ,Eⱼ)` are **free rational Taylor coefficients** of the unknown
polynomials `d2,d1,σ,e` about `y=−1` (σ’s jets are free through the free `d0`,
`σ=4d0−d2²`). They are *not* confined to the `q`-splitting field — that
confinement is the `q`-place phenomenon behind the only two kills C08/C20.
Consequently:

> **Smooth-point tower lemma.** If a residue point lies on `IF=0` with all
> required leading coefficients nonzero and `∇IF ≠ 0` (a *smooth* point of the
> hypersurface), then every `C_δ = 0` is a linear equation in one order-δ jet
> with nonzero coefficient, solvable for arbitrary lower data. Hence the full
> jet tower is solvable **to every depth δ**, i.e. the extended system is a
> **CONSTRAINT for every a** (each `a` demands `δ = 30−3a`), all at once.

A KILL could arise only if `∇IF ≡ 0` on the entire nonzero-leading locus of
`IF=0` (a hypersurface singular exactly there). §3–§4 show this never happens.

---

## 1. Census of the 1145 survivors  [data]

`residue_lemmas_depth_verify.py` §V2. Survivor obligations by place/kind:

| place | kind | count |
|:--|:--|--:|
| q | term_cancellation | 3413 |
| q | monomial_tie_rise | 2005 |
| q | exact_identity | 1656 |
| q | identical_vanishing | 968 |
| t | term_cancellation | 1441 |
| t | monomial_tie_rise | 734 |
| t | exact_identity | 414 |
| t | identical_vanishing | 242 |
| inf | leading_cancellation / degree_tie_drop / … | 7430 |

Every t-place tied monomial string is an exact term of the source `h_l`
(**0 mismatches**). Branches by `a`: `{2:2, 3:11, 4:14, 5:22, 6:25, 7:26,
8:26, 9:24, 10:21}`; cases by `a`: `{2:9, 3:59, 4:102, 5:135, 6:154, 7:165,
8:176, 9:168, 10:177}`.

**The affine t-depth law [data, §V3].** The dominant t-place
`term_cancellation` depth equals `30−3a` exactly at every `a`:

| a | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| mode t-depth | 24 | 21 | 18 | 15 | 12 | 9 | 6 | 3 |
| 30−3a | 24 | 21 | 18 | 15 | 12 | 9 | 6 | 3 |

(`a=10`: the level collapses to depth-0 `exact_identity`, `30−30=0`.)

### 1.1 What actually carries the a-growing depth  [data]

Two families carry the growing depth, and they are different in kind:

- **t-place `term_cancellation`, tied = ∅** (highest incidence: L6 598, L5
  481, L4 362 occ). These constrain the *free leading residue of `g_l`*, not
  the window residue polynomial. The jet tower of §0 does not attach an
  arithmetic obstruction to them — `g_l`’s coefficients are free — so they are
  Phase-D "always available" (STATE.md), never a local kill.
- **t-place `monomial_tie_rise`, tied = an `h_l` cut** (the window residue
  equations). These are the ones for which §0’s jet structure is defined. Ten
  distinct tied supports occur; the three highest-incidence, whose depth grows
  with `a`, are the subject of §3.

The a-growing window patterns, ranked by incidence:

| rank | place | L | support = `h_l` cut | RESIDUE_LEMMAS id | occ | cells | a range | depths |
|--:|:--|--:|:--|:--|--:|--:|:--|:--|
| 1 | t | 5 | `{d2²d1², d2σ², d1²σ}` | **C09** | 201 | 17 | 2–9 | 3,6,…,24 |
| 2 | t | 6 | `{σ², d2d1²}` | **C02** | 196 | 10 | 3–9 | 3,6,…,21 |
| 3 | t | 4 | `{σ³, d1⁴, d2d1²σ, d2³d1², d2²σ²}` | **C22** | 160 | 17 | 2–9 | 3,6,…,24 |
| 4 | t | 4 | `{σ³, d1⁴}` | C12 | 65 | 65 | 3–9 | 3,…,21 |
| 5 | t | 4 | `{σ³, d2²σ²}` | C14 | 44 | 44 | 5–9 | 3,…,15 |
| 6 | t | 6 | `{σ², d1e}` | C01 | 32 | 24 | 2–8 | 2,…,22 |
| 7 | t | 4 | `{d1⁴, d2³d1²}` | C16 | 21 | 21 | 4–9 | 3,…,18 |
| 8 | t | 6 | `{d2d1², d1e}` | C03 | 6 | 6 | 4–8 | 1,…,14 |
| 9 | t | 4 | `{σ³, d1σe}` | C11 | 6 | 6 | 6–8 | 1,5,9 |
| 10| t | 4 | `{σ³, d1⁴, d1σe}` | C18 | 3 | 3 | 5 | 13 |

The depth multiset of each pattern is `{30−3a : a in its range}` minus the
tied monomial’s e-slot offset — hence affine in `a` (§0). **Every occurring
support is a CONSTRAINT-class id (C01–C23 minus the kills); neither kill
support C08 nor C20 appears anywhere on the 1145 survivors** (§5).

---

## 2. Source identification of the top three  [data, §V4]

Rebuilt from `f31_graded.txt` (coefficients pulled by matching the support in
the parsed `h_l`, never copied):

```
C09 (L5):  IF = −12288 D²X² − 9216 D S² + 32256 S X²         = −1536·(8X²D²−21X²S+6DS²)
C02 (L6):  IF =  14336 D X² − 3072 S²                        =  1024·(14X²D−3S²)
C22 (L4):  IF = −31232 D³X² − 5184 D²S² − 23616 D S X²
                −12096 S³ − 220752 X⁴                        =   −16·(13797X⁴+1952X²D³
                                                                     +1476X²DS+324D²S²+756S³)
```

matching the RESIDUE_LEMMAS.md primitive relations for C09/C02/C22 up to the
displayed nonzero rational multiples.

---

## 3. Depth-2 and depth-3 jet equations, and the verdicts  [data, §V5]

For each pattern the checker substitutes the depth-2 jet series into the
source tie polynomial and extracts `C₀,C₁,C₂` symbolically, confirming the §0
identities `C₀=IF`, `C₁=∇IF·jet₁`, `C₂=∇IF·jet₂+½ jet₁ᵀH jet₁`.

### 3.1 Pattern 1 — C09 (L5), 201 occ, cells `{0000,1000,1100,1110,1111,2000,2100,2110,3000,3100,3110,3111,3300,3310,5000,5100,5110}T1`, a=2–9

```
C₁ = D₁(−24576 DX² − 9216 S²) + X₁(−24576 D²X + 64512 SX)
        + S₁(−18432 DS + 32256 X²)                                = 0    (depth 2)

C₂ = 32256 S₂X² − 9216 D₂S² − 24576 D₂X² + … + 32256 S X₁²
        − 12288 D₁²X² − 18432 D₁S S₁ − 9216 D S₁²                 = 0    (depth 3)
```

Smooth point (RESIDUE_LEMMAS witness, verified against source IF): `(D,X,S) =
(13/6, 13/6, 169/36)`, all nonzero, `IF=0`, `∂IF/∂D = −4077632/9 ≠ 0`.
Solving `C₁=0`, `C₂=0` for the pivot jet (forcing `X₁,S₁,X₂,S₂` to distinct
nonzero values) gives the **exact depth-3 solution** `D₁ = 954/377`,
`D₂ = 22464360/53582633`. **Verdict: CONSTRAINT at every depth, all a=2–9.**

### 3.2 Pattern 2 — C02 (L6), 196 occ, cells `{0000,1000,1100,1110,1111,3000,3100,3110,3111,5000}T1`, a=3–9

```
C₁ = 28672 DX·X₁ + 14336 X²·D₁ − 6144 S·S₁                       = 0    (depth 2)
C₂ = 28672 DX·X₂ + 14336 X²·D₂ − 6144 S·S₂
        + 14336 D X₁² + 28672 D₁X X₁ − 3072 S₁²                   = 0    (depth 3)
```

Smooth point `(D,X,S)=(3/14,1,1)`, `IF=0`, `∂IF/∂D = 14336 ≠ 0`. Exact
depth-3 solution `D₁ = −3/7`, `D₂ = 15/14`. **Verdict: CONSTRAINT at every
depth, all a=3–9.**

### 3.3 Pattern 3 — C22 (L4), 160 occ, cells `{0000,1000,1100,1110,1111,2000,2100,2110,3000,3100,3110,3111,3300,3310,5000,5100,5110}T1`, a=2–9

```
C₁ = D₁(−93696 D²X² − 10368 DS² − 23616 SX²)
        + X₁(−62464 D³X − 47232 DSX − 883008 X³)
        + S₁(−10368 D²S − 23616 DX² − 36288 S²)                   = 0    (depth 2)

C₂ = ∇IF·(D₂,X₂,S₂) + ½(D₁,X₁,S₁)ᵀ Hess(IF) (D₁,X₁,S₁)           = 0    (depth 3)
       [ = −93696 D²D₂X² − 62464 D³X X₂ − … − 1324512 X²X₁², full form in §V5 ]
```

Smooth point `(D,X,S) = (511/152, 511/152, −261121/17328)` (the `u=152/511`
witness), `IF=0`, `∂IF/∂D = −39274085745216/2476099 ≠ 0`. Exact depth-3 solution
`D₁ = −11134/1533`, `D₂ = −371919751915/14410745748`. **Verdict: CONSTRAINT at
every depth, all a=2–9.**

---

## 4. General no-jet-kill theorem  [data, §V6]

For **all ten** occurring t-place `monomial_tie_rise` supports (§1.1 table),
the checker rebuilds `IF` from source, confirms a consistent nonnegative
valuation-weight exists (the tie is realizable), and verifies a rational
nonzero-leading point on `IF=0` with `∇IF ≠ 0` (the RESIDUE_LEMMAS §4
certificates, checked against the source-derived `IF`). By the smooth-point
tower lemma (§0):

> **Every t-place `monomial_tie_rise` on the 1145 survivors is a CONSTRAINT at
> every depth `δ = 30−3a`, hence for every `a ∈ [2,10]`. No depth-2, depth-3,
> or any deeper jet equation turns a surviving t-place tie into a KILL.**

This is *one* proof per support covering all `a` at once: smoothness of a
single rational point discharges the entire affine depth tower. It rigorously
substantiates STATE.md’s "deep t-place cancellations are always available to
zero-budget states" for the backbone window ties.

---

## 5. Kill accounting  [data, §V7]

- **No new kills.** The depth extension produces no KILL: §4 shows every
  t-place tie is smooth, so no jet obstruction exists at the (free, finite)
  t-place. Arithmetic (square-class) obstructions live only at the `q`-place,
  where residues are confined to the `q`-splitting field `Q(√17)`.
- **The two program kills are q-place, depth-1, a-independent.** C08
  (`6X²D²−9XDE−E²`, disc square-class 105) and C20 (`11r²−6r−61`, square-class
  170) are the only kills (`RESIDUE_LEMMAS.md` §4). Their **supports occur 0
  times on the 1145 survivors** — their carrying branches were already removed
  upstream by the residue-kill feedback (`residue_kills: true`), so they sit in
  the 2007 `engine_killed_pending_audit` branches, not the frontier.
- **Instances any kill eliminates in this frontier: none** (they acted before
  it). Their upstream impact is quantified in `RESIDUE_LEMMAS.md` §6 (C08: 304
  obligations / 54 cells; C20: 17 / 8, sub1).

**Net.** On the surviving frontier the a-growing depth is entirely
CONSTRAINT-class: the affine law `depth = 30−3a` inflates the required jet
order, but smoothness makes the tower solvable at every order, uniformly in
`a`. The only obstructions in the whole program remain the two a-independent
depth-1 q-place kills.

---

## 6. [judgment] list — what is asserted but not fully forced

1. **t-place jets are free rational parameters.** The core of every CONSTRAINT
   verdict. It rests on: `t=y+1` being a finite place (verified nonzero
   `q(−1)`), and the Taylor coefficients of `d2,d1,σ,e` about `y=−1` being
   unconstrained *for the local obligation in isolation* (σ’s jets free via the
   free `d0`). This is the standard **local** scope of a residue lemma
   (`RESIDUE_LEMMAS.md` intro); it does not assert the jets are simultaneously
   free of the global degree caps or the q-place system — only that the t-place
   obligation is not a *local* kill. A global (Phase-D) closure may still
   impose more, exactly as `exact_identity`/`identical_vanishing` do.
2. **"depth = 30−3a governs the required δ."** Verified as the dominant law
   (§V3); a thin tail of offset depths (e.g. 25, 26, 20 at a=2,4) sits at
   `v` minus an e-slot offset and is also affine. The smooth-point lemma covers
   *all* depths regardless, so the verdict is insensitive to the exact δ.
3. **The general theorem is stated for the ten supports that occur**, not for
   hypothetical supports. Should the frontier shift, a new support would need
   its own smooth certificate (all C01–C23 except C08/C20 already have one).
4. **term_cancellation (tied=∅) t-obligations are treated as free `g_l`
   residues**, hence not jet-killable. This matches the engine semantics
   (`cascade_engine.py` `descend_options`) but is a Phase-D filing, not a
   theorem proved here.

---

## Appendix: reproduction

```
python residue_lemmas_depth_verify.py
```

Standalone, read-only. Re-parses `f31_graded.txt`, rebuilds the census over
`cascade_cones_sub1_qt_inf_rl.json` (0 tied-vs-source mismatches), checks the
affine law, identifies C09/C02/C22, proves the `C₀/C₁/C₂` jet identities,
exhibits each exact depth-3 solution, verifies the smooth certificate for all
ten t-place supports, and confirms the C08/C20 kill supports are absent from
the survivors. All checks PASS.
