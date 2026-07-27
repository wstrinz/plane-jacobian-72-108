# BELYI_PASSPORT — the top band of (72,108) is a Belyi/Hurwitz problem

2026-07-25. Repo `d2_plane_72_108`, task HEAD **`1e2d99b`** (the tree had advanced
to `877ffa4` when this landed). Checker: `belyi_passport.py`
(**90/90**, `--quiet` exit 0, ~6 min; **95/95** with `--singular`, **96/96**
with `--full`; `--fast` skips the number-field cross-check and needs no external
bundle). Read-only — no existing file in this repo or in Helali's artifact was
touched. New files only: this one and `belyi_passport.py`.

---

## 0. HEADLINE

> **The Hurwitz number of the passport `(2^10 1 | 3^7 | 17 1^4)` in degree 21 is
> exactly `5`. It is NOT 35.**

The task hypothesis — that the `35` of Helali's first block *is* the Hurwitz
number — is **false**. But the correspondence is real and is now closed exactly:

```
                  35   =   7    ×   5
                  |        |        |
        vdim of   |        |        +--  the Hurwitz number: 5 dessins
        his first |        +-----------  a residual mu_7 in his normalisation
        block     |                      a_1 = a_8 = 1  (acts FREELY)
                  +--------------------  our own from-scratch recomputation,
                                         vdim 35, DIM 0, and the SAME ideal
```

and the payoff is bigger than the original hypothesis would have been:

> **`L = Q[w]/(w^5 − w^4 + 3w^3 + 3w^2 + 26)`, the degree-5 field his whole
> endgame descends to, is the FIELD OF MODULI of a single Galois orbit of the
> 5 dessins — and `[L:Q] = 5` IS the Hurwitz number.**

Every number in his first block is now accounted for by representation theory,
with no Gröbner basis in the derivation. Concretely, his printed lex basis makes
the descent visible on its face: `H(a_7)` is *literally a polynomial in `a_7^7`*,
supported on exponents `{0,7,14,21,28,35}`, and each of `L[2]…L[6]` is
`a_7^{c_j}` times a polynomial in `a_7^7`. That is the `mu_7` grading, printed.

**Independent corroboration.** An independent party computed the same Hurwitz
number, 5, by the same route (Frobenius/Murnaghan–Nakayama), and independently
derived the same residual `mu_7`. My computation was performed before receiving
theirs and agrees to the number, to the freeness argument, and to the weights
`(1,2,3,4,5,6)` on `a_2..a_7`. **Two independent computations, same answer.**

### And a structural bonus

(72,108) is not a one-off. It is the member `k = 7` of a family, and the family's
Hurwitz numbers are the **Catalan numbers** (`J2`, verified for `k = 1..11`, and
`k = 13` off-line):

| `k = deg a` | 1 | 3 | 5 | **7** | 9 | 11 | 13 |
|---|---|---|---|---|---|---|---|
| `deg Φ` | 3 | 9 | 15 | **21** | 27 | 33 | 39 |
| `deg d` | 1 | 4 | 7 | **10** | 13 | 16 | 19 |
| `e_∞` | 2 | 7 | 12 | **17** | 22 | 27 | 32 |
| **Hurwitz** | 1 | 1 | 2 | **5** | 14 | 42 | 132 |
| **predicted `vdim`** `= k·h` | 1 | 3 | 10 | **35** | 126 | 462 | 1716 |
| **`vdim` measured** | 1 | **3** | **10** | **35** | — | — | — |

Passport at level `k` (odd): `( 2^{(3k−1)/2}, 1 | 3^k | (5k−1)/2, 1^{(k+1)/2} )`.
The Hurwitz number is `C_{(k−1)/2}`. The `vdim = k · Hurwitz` law is *measured
independently* at `k = 1, 3, 5, 7` — three of them are not (72,108), so this is a
real test of the correspondence and not a restatement of it.

---

## 1. PROVED / CHECKED / INFERRED

### PROVED (mathematical argument, machine-verified here, not machine-dependent)

- **P1. The passport is FORCED by (J4) — there is no genericity hypothesis.**
  Write `A = t·a`, `D = t²·d`; (J4) `2AD′ − 3A′D = t²` is equivalent to
  `a·d + 2t·a·d′ − 3t·a′·d = 1`. Reduce that identity mod `a`: every term but one
  carries a factor `a`, so at **any** root `t₀` of `a`,
  `−3 t₀ a′(t₀) d(t₀) = 1`. Hence `t₀ ≠ 0`, `a′(t₀) ≠ 0`, `d(t₀) ≠ 0`. Reduce mod
  `d`: at any root `t₁` of `d`, `2 t₁ a(t₁) d′(t₁) = 1`, so `t₁ ≠ 0`,
  `d′(t₁) ≠ 0`, `a(t₁) ≠ 0`. Therefore **`a` has 7 distinct nonzero roots, `d` has
  10 distinct nonzero roots, and `gcd(a,d) = 1` — proved, not assumed.** The
  `t⁰` coefficient is `a₁d₂ = 1`, so neither vanishes at 0 either. Checker `F1–F3`.
  *This is what makes the correspondence a bijection rather than a
  generic-fibre statement.*
- **P2. `β = D²/A³ = t·d²/a³` is a Belyi map of degree 21 with the stated
  passport.** `β′ = d·(ad + 2tad′ − 3ta′d)/a⁴ = d/a⁴` under (J4) (`E4`). With P1:
  ten double zeros from `d` plus the simple zero `t = 0` over `0`; seven triple
  poles from `a` over `∞`; and at `t = ∞`, `β(1/u) = d̃(u)²/ã(u)³` with
  `ã(0) = a₈`, `d̃(0) = d₁₂`, so `β(1/u) − β(∞)` vanishes to order exactly
  `deg a + deg d = 17`. Riemann–Hurwitz is then an *equality*
  (`10 + 14 + 16 = 40 = 2·21 − 2`), which forces the third fibre to be
  `(17, 1⁴)` and forces there to be **no fourth branch point**. `B1–B3`, `I9a–d`.
- **P3. Every triple in this passport is transitive.** Orbits are unions of
  `σ_1`-cycles, all of length 3, so every orbit size is divisible by 3. `σ_2` has
  a 17-cycle, so some orbit has size ≥ 17, hence 18 or 21. If 18, the complement
  is a single orbit of size 3 on which `σ_2 = id`, `σ_1` is a 3-cycle and `σ_0`
  an involution, and `σ_0σ_1 = id` would make `σ_0` a 3-cycle. Contradiction.
  `B4a`, `B4b` (the 3-set case enumerated exhaustively).
- **P4. Every cover here is rigid.** A deck automorphism is semiregular, so its
  order divides 21; it preserves cycle lengths of `σ_2`, hence fixes the unique
  17-cycle setwise and acts semiregularly on those 17 points, so its order also
  divides 17. `gcd(21,17) = 1`. Therefore `#dessins = N/21!` **exactly**, with no
  `1/|Aut|` weighting. `B5`.
- **P5. The residual group is `mu_7` and it acts freely.** `A ↦ μA(λt)`,
  `D ↦ νD(λt)` preserves (J4) iff `μνλ³ = 1`; preserving `a₁ = a₈ = 1` forces
  `μλ = 1` **and** `μλ⁸ = 1`, whose ratio gives `λ⁷ = 1`, then `μ = λ⁻¹`,
  `ν = λ⁻²`, and `a_k ↦ λ^{k−1} a_k` — weights `(1,2,3,4,5,6)` on `a₂..a₇`. None
  is `0 mod 7`, so a nontrivial fixed point forces `a₂ = … = a₇ = 0`; the origin
  is **not** a solution (residual `R₁₄` evaluates to `−8` there). Hence the action
  is free. `G4–G6`.
  *Machine-visible consequence:* every one of the six residual generators is
  `mu_7`-homogeneous, of character `n mod 7` for `n = 11..16` — characters
  `4,5,6,0,1,2`. `G4`.
- **P6. The general engine (see §4) and its Riemann–Hurwitz identity.** For any
  `(m, n, k, ℓ)`, `ℓ(m−1) + (nk−mℓ−1) + k(n−1) + (k+ℓ−1) = 2nk − 2` **identically**.
  `K1`. Also `p = ord A + ord D − 1` and `e_∞ = deg a + deg d` are forced, not
  fitted (`K2`, `K3`).

### CHECKED (exact machine computation, reproducible)

- **C1. Hurwitz number = 5.** Frobenius over all 792 irreps of `S_21`, with
  Murnaghan–Nakayama; only **19** irreps have all three characters nonzero (the
  17-cycle prunes hard). `Σ χ₀χ₁χ₂/dim = 31104/19019`;
  `N = 255 454 710 858 547 200 000`; `N/21! = 5` exactly. `C1–C3`.
- **C2. The character machinery is validated, not assumed.** `S_6` character
  table column orthogonality across all 11 classes; `Σ dim² = 8!`; and Frobenius
  vs **brute-force enumeration in `S_n`** on five genuinely branched passports
  including the `k = 1` member of our own family. `A1–A7`.
- **C3. Our own recomputation of the first block: `vdim = 35`, `DIM = 0`.**
  Rebuilt from scratch from the coefficient recursion
  `Σ_i α_i δ_{n−i} (1 + 2n − 5i) = [n = 0]`, `α₀ = α₇ = 1` — never read from his
  files. `modStd(I,1)` over **Q**. The `n = 17` equation is identically vacuous
  (17 equations, 19 unknowns, 2-dimensional symmetry group — consistent), and the
  11 pivots are `1,3,5,…,21`. `G1–G3`, `M1`.
- **C4. THE IDEALS ARE EQUAL.** My six from-scratch residuals evaluate to
  **exactly zero** at Helali's degree-35 lex point, computed in `Q[a₇]/(H)` on
  `python-flint`. So `V(his) ⊆ V(mine)`; `H` is irreducible of degree 35 so
  `V(his)` is 35 distinct points; `vdim(mine) = 35`; therefore `V(mine)` is those
  same 35 points, `I_mine` is radical, and `I_mine = I_his`. **`H6`.**
  *This closes the audit seam: no computation is being taken on trust in either
  direction.*
- **C5. `H(a_7) = G(a_7^7)`, `deg G = 5`, `G` irreducible, and `G` has a root in
  `L`.** Support of `H` is exactly `{0,7,14,21,28,35}` (`H2`). Helali's own
  `build_degree5.py` asserts `G(Φ) ≡ 0 mod F`; I re-verified that with my own
  `flint` arithmetic (`H4`). Hence `Q(a₇^7) = L`.
- **C6. THE INVARIANT TEST PASSES.** `a₂a₇`, `a₃a₆`, `a₄a₅` and `a₇^7` are all
  `mu_7`-invariant (weight sums `7, 7, 7, 42`, all `≡ 0 mod 7`), all have
  `a₇`-support contained in the multiples of 7, and each has characteristic
  polynomial over `Q(a₇)` equal to **(a degree-5 irreducible)^7**. Since each lies
  in `Q(a₇^7) = L` and has degree 5 over `Q`, and `[L:Q] = 5`, **each generates
  `L` exactly**. `H9`, `H10`. For `a₂a₇` the minimal polynomial is
  ```
  2768896 v^5 − 526303232 v^4 + 58496180029 v^3
            − 3550404654747 v^2 + 138280970186991 v − 1707884082093393
  ```
  with Galois group `S_5` — as have `w^5−w^4+3w^3+3w^2+26` and the other two
  invariants. (`S_5` means the 5 dessins carry the *full* symmetric Galois action;
  they are one orbit and nothing smaller.)
- **C7. Degeneracy is excluded.** `d₁₂ = δ₁₀` reduces to a **nonzero** element of
  `Q[a₇]/(H)`; `H` is irreducible, so that element is invertible, so
  `d₁₂ ≠ 0` at **all 35 points** and `deg d = 10` exactly (`H7`). Independently,
  in Singular the locus `V(I) ∩ {δ_max = 0}` is **empty** (`dim = −1`) at
  `k = 3, 5, 7` (`M2`). Together with P1 (`res(a,d) ≠ 0`, `disc a ≠ 0`,
  `disc d ≠ 0` are *proved*) every degeneracy that could fall outside the Belyi
  correspondence is excluded.
- **C8. End-to-end on an actual solution.** At a numeric root of `H` (90 digits):
  the six residuals vanish to `1e−85`; `a` has 7 distinct roots (min separation
  0.49), `d` has 10 (0.17), they share none (0.025), none is 0; `β(∞) = d₁₂²`.
  Then **exactly**, in `Q[a₇]/(H)`: the `u^1 … u^16` coefficients of `β(1/u)` are
  all zero and the `u^17` coefficient is not. Sixteen exact vanishings in a
  degree-35 number field. `I1–I9`.
- **C9. `vdim = k · Hurwitz` at `k = 1, 3, 5, 7`** → `1, 3, 10, 35` against
  Hurwitz `1, 1, 2, 5`. `M1`.

### INFERRED (stated as such; not proved here)

- **N1.** That the Hurwitz numbers of this family are *exactly* the Catalan
  numbers `C_{(k−1)/2}` — verified at `k = 1,3,5,7,9,11,13`, not proved. A proof
  would presumably come from the dessin being a plane trivalent tree-like map
  with one big face; I did not construct that bijection.
- **N2.** That the 5 dessins form a *single* Galois orbit. This follows from
  irreducibility of `H` over `Q` (35 = 7·5 with blocks of size 7 from `mu_7`,
  hence a transitive action on the 5 blocks), which I take from Helali's `H` —
  whose irreducibility I checked, but whose *provenance* as the eliminant is his
  computation plus my C4 (which proves the ideals agree, so this is solid; I
  label it INFERRED only because I did not re-run the elimination itself over Q).
- **N3.** The `(75,125)` band functional (§5). Arithmetically necessary, not
  sufficient — the polygon has not been derived.

---

## 2. The observation, re-derived from the actual polygons

Not taken from the task summary. From `N(P) = conv{(0,0),(1,0),(8,14),(8,16),(0,8)}`,
`N(Q) = conv{(0,0),(2,1),(12,21),(12,24),(0,12)}`, `[P,Q] = x²`:

- `t = xy²`, `z = y⁻¹` gives `x^i y^j = t^i z^{2i−j}` (`D1`), and
  `det ∂(t,z)/∂(x,y) = −1`, so **`[P,Q]_{t,z} = −x² = −t²z⁴`** (`D2`).
- Band census over the actual lattice points: `max(2i−j)` on `N(P)` is **2**,
  attained at `i = 1..8` (8 points ⇒ `A = t·a`, `deg a = 7`); on `N(Q)` it is
  **3**, attained at `i = 2..12` (11 points ⇒ `D = t²·d`, `deg d = 10`). `D3`, `D4`.
- The `z⁴` coefficient of `P_tQ_z − P_zQ_t` is `3A′D − 2AD′` (`D5`), so
  **(J4): `2AD′ − 3A′D = t²`** (`D6`).
- `d/dt(D^m/A^n) = (D^{m−1}/A^{n+1})(mAD′ − nA′D)` — verified at
  `(m,n) = (1,1),(2,3),(3,5),(4,7),(5,2)` (`E1`).
- Hence `β = D²/A³ = t d²/a³` and `β′ = d/a⁴`. `E2–E4`.

The Riemann–Hurwitz balance quoted in the brief is confirmed, but note §4: **it
was never a test.** It is an identity in `(m,n,k,ℓ)`.

---

## 3. Where 35 came from, exactly

The solution variety of (J4) in the 19 unknowns `a₁..a₈, d₂..d₁₂` carries a
2-dimensional group `(a,d) ↦ (α·a(λt), α⁻¹·d(λt))`. The `t^17` coefficient
equation is identically zero, leaving 17 equations — so the variety is
2-dimensional and the group acts with finite stabiliser. Helali's slice
`a₁ = a₈ = 1` meets each orbit in exactly **7** points (P5). Therefore

```
   vdim(first block) = 7 · #(orbits) = 7 · #(dessins) = 7 · 5 = 35 .
```

The `35` is a *normalisation artefact times* the Hurwitz number. And the block's
`mu_7` structure is not a hidden inference — it is printed in his own lex basis
(`H2`, `H5`) and used, unnamed, in his own `build_degree5.py`
(`G = fmpq_poly([H[7*k] for k in range(6)])`, plus an `assert` that every
coefficient array is supported on multiples of 7).

---

## 4. The general engine

**Inputs.** Reduced Newton polygons `N(P)`, `N(Q)`; a *primitive* integer
functional `f(i,j) = u·i − v·j` attaining its maximum on a face of **both**; and
the bracket exponent `κ` in `[P,Q] = x^κ`. Put `m = max_{N(P)} f`,
`n = max_{N(Q)} f`, and let `A`, `D` be the corresponding top-band coefficient
polynomials in `t`.

| step | statement |
|---|---|
| **top layer** | `m A D′ − n A′ D = γ t^p` |
| **THE GATE** | `γ ≠ 0` **iff** `u·κ = m + n − 1` |
| reduce | `g = gcd(m,n)`; the top object is `Φ = D^{m/g} / A^{n/g}` |
| forced | `p = ord_t A + ord_t D − 1` (automatic when `m·ord D ≠ n·ord A`) |
| forced | `m·deg D = n·deg A` (else `γ t^p` could not have low degree) |
| write | `A = t^{α_A} a`, `D = t^{α_D} d`, `k = deg a`, `ℓ = deg d` |
| **degree** | `deg Φ = n·k = (nk − mℓ) + m·ℓ` |
| **passport** | over `0`: `( m^ℓ , nk − mℓ )`; over `∞`: `( n^k )`; third: `( k+ℓ , 1^{nk−k−ℓ} )` |
| **RH** | **balances identically** — it is *not* a test |
| **rigidity** | `a`, `d` have simple, nonzero, coprime roots — *forced* by the equation |
| **count** | `vdim(normalised block) = (deg a) × (Hurwitz number)`, the `deg a` being a residual `μ_{deg a}` |

(72,108) is `u=2, v=1, κ=2, m=2, n=3, α_A=1, α_D=2, p=2, k=7, ℓ=10`:
`2·2 = 2+3−1` ✔, `e_∞ = 17` ✔, `deg Φ = 21` ✔, `vdim = 7·5 = 35` ✔.

### What this still needs to be a compiler

1. **Explicit reduced polygon vertices for an arbitrary corner.** This is the
   binding constraint. `polygon_reduction.py:436-440` emits them only for
   (72,108) and explicitly declines for the F2 family
   (`"the exact reduced vertex list is (m,n)-scaled {(0,0),(a0-q,0)? ...}"`,
   `red.pre_inversion = {}  # no published vertex list`). Without them,
   `deg A`, `deg D`, `ord A`, `ord D` are unknown and no passport can be written.
2. **A face-pairing rule.** The engine needs a functional maximised on a face of
   *both* polygons. At (72,108) both have an edge of direction `(1,2)`; there is
   no general theorem here yet.
3. **The gate, `u·κ = m + n − 1`, decided from corner data alone.** §5 shows the
   gate is a genuine, restrictive Diophantine condition — this is the real
   discriminator, and it is cheap once (1) and (2) are available.
4. **Freeness of the residual `μ_{deg a}`** — a one-line check per corner (is the
   origin a solution of the block?), but it must be checked, not assumed.

Given 1–4, the passport, the RH balance, the rigidity and the count are all
free — and one gets the field of the endgame *before running any CAS*.

---

## 5. (75,125) — an honest negative, plus a sharp prediction

**Inputs that are in the repo** (`L1`): `(a,b) = (3,5)` (`phi_75_125.py:27`),
`κ = 3` so `[P,Q] = x³` — *not* `x²* (`POLYGON_REDUCTION.md:108`,
`PHI_75_125.md:65`), `C = y²(y³+1)`, `R = x⁵C`, `l(P) = R³`, `l(Q) = R⁵`,
corner `(5,20) → (7/5,2)`, `t = l = 5`, `a₀ = 5`, `q = 2`.

**Result 1 — the obvious transplant is a tautology.** The face where the leading
forms literally are `R³` and `R⁵` is the `x`-degree face, `f = i`, `u = 1`. The
gate then reads `1·3 = 15 + 25 − 1 = 39` — **false**, so `γ = 0` and the top
equation degenerates to `3AD′ − 5A′D = 0`, i.e. `D³ = cA⁵`. In a UFD that says
`A = c₁E³`, `D = c₂E⁵`: it merely restates `l(P) = R³`, `l(Q) = R⁵`. **No Belyi
map, no finite classification, positive-dimensional solution set.** The same
computation at (72,108) on *its* `x`-degree face gives `2 ≠ 19` and
`D² = cA³` — equally vacuous. `L2`.

> So: writing down `Φ = D³/A⁵` at (75,125) because the reduced powers are `(3,5)`
> gives a true but empty statement. The brief's guess at the *shape* is right;
> the naive placement of it is not.

**Result 2 — the gate, and the minimal admissible face.** On a face derived from
`R` with `r = max_R f`, one has `m = a·r`, `n = b·r`, so the gate becomes

```
        u·κ  =  (a+b)·r − 1 .
```

At (72,108) (`κ=2`, `a+b=5`) this gives `2u = 5r−1`, solved by `r=1, u=2` —
i.e. `f = 2i − j`, **exactly the functional the coordinate change realises**
(`L3`). At (75,125) (`κ=3`, `a+b=8`) it reads `3u = 8r − 1`, and combining with
`r = 5u − 2v` (from `R`'s monomials `(5,2)`, `(5,5)`) forces `r ≡ 5 mod 6`:

| `r` | `u` | `v` | functional |
|---|---|---|---|
| **5** | **13** | **30** | **`f = 13i − 30j`** |
| 11 | 29 | 67 | `f = 29i − 67j` |
| 17 | 45 | 104 | `f = 45i − 104j` |

The minimal one gives `(m,n) = (15,25)`, `gcd = 5`, reduced `(3,5)` — so the top
object *is* `Φ = D³/A⁵`, but on the face `13i − 30j`, not on the `x`-degree face.
`L4`. **This is INFERRED (N3): arithmetically necessary, not sufficient.**
Whether `N(P)` and `N(Q)` actually have parallel edges of direction `(30,13)` is
a polygon question, and the polygon is not available.

**Result 3 — "does RH balance?" is not a test.** For `(m,n) = (3,5)` the balance
holds for **every** admissible `(k,ℓ)` — 484 pairs checked, all balance (`L5`),
as it must by the identity `K1`. The brief's success criterion "if it balances,
that is the beginning of a compiler" therefore cannot be met by balancing; the
discriminating condition is the gate `u·κ = m + n − 1`.

**Result 4 — what is blocked.** The concrete (75,125) passport is **not**
computed here, and I decline to guess it. `deg A`, `deg D`, `ord A`, `ord D`
require the reduced vertex list, which `polygon_reduction.py` declines to emit
(`L6`). That derivation — flip, shift to the foot `(−2,0)`, inversion
`(a,b) ↦ (5b−a, b)` applied to the F2 feet — is a self-contained next task.

This is consistent with, and independent of, `weight_lemma_75_125.py` (38/38):
that lane showed the **Φ-divisor** mechanism is not family-level. This lane shows
that the **top-band Belyi** mechanism is *also* not a transplant — for a
different and more elementary reason (the gate fails on the obvious face). The
two negatives are not the same negative.

---

## 6. Two limits, stated plainly

1. **This governs only the top band.** (J4) is one layer. `J3 … J0`, the negative
   bands, the `s = ±c` split, the `r`-elimination, the `h = 0 / h ≠ 0` branching
   and the 89 MB identity get **nothing** from this. The Belyi picture explains
   the *first block* and nothing below it.
2. **A nonzero Hurwitz number never kills a case.** It hands you a number field.
   That is exactly what happened here: the top band admits 5 perfectly good
   covers, and (72,108) died further down. What the Belyi layer buys is **exact,
   instant control over which field the endgame will live in, before any CAS
   runs** — historically the expensive thing to discover. At (72,108) that field
   is `L`, of degree 5, and one could have known `5` from a character sum in
   seconds.

---

## 7. Reproduction

```
python -u belyi_passport.py               # 90/90, ~6 min (96/96 with --full)
python -u belyi_passport.py --quiet       # exit 0 iff all pass
python -u belyi_passport.py --fast        # skips sections H, I (no bundle needed)
python -u belyi_passport.py --full        # all four mu_7 invariants
python -u belyi_passport.py --singular    # adds the vdim jobs (k=3,5,7) via WSL
```

Sections `H` and `I` read (never write) Helali's
`firstblock_Q_exact.out` from the extracted bundle; they self-skip if it is
absent. Nothing in his artifact or in this repo is modified.

**Dependencies:** `sympy`, `python-flint`, `mpmath`; `--singular` needs
Singular 4.2.1 via WSL.
