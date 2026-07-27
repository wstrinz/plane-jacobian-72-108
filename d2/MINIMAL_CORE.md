# MINIMAL_CORE.md — how short is the (72,108) proof, and what theorem did we prove?

Checker: `minimal_core.py` (60/60, `--quiet` exits 0 iff all pass; `--deep` extends
the character sums to k = 19). Companion scratch: `_scratch_hurwitz.py`,
`_scratch_trees.py`, `_scratch_vdim9.py`.

**This lane changes no verdict.** It adds no kill, touches no existing file, and
compiles nothing into the frontier. Everything below is either an *ablation* of an
already-audited fact or a *theorem* about the Belyi side.

---

## 0. Headline

**The minimal load-bearing spine.** Downstream of the upstream GGHV reduction, the
whole closure of (72,108) runs on **four polynomials, one polynomial identity, two
slice families, one valuation cascade, and nine integers.**

```
              ┌─ [U1] GGHV Prop 4.3: (72,108) → corner (8,28), two Newton configs
   UPSTREAM   │        ** CONDITIONAL: the uncited 7th-power sentence, L1132 **
              └─ [U2] C4 = y^7(y+1)
                              │
   ┌──────────────────────────┴──────────────────────────────────┐
   │  [P1]  G1, G2, G3, G5body           (= r_13, r_14, r_15, r_17)
   │        four polynomials, eight variables.  generators.json
   │  [P2]  Φ = −(1/6630)·y^204·t^30·q,  q = 2048y⁴−512y³+320y²−240y+195
   │        q irreducible squarefree, q(−1) = 3315,  v_t(Φ) = 30 EXACTLY, deg Φ = 34
   └─────────────────────────────┬───────────────────────────────┘
                                 │
                    [P3]  THE K-SYZYGY  (one expansion)
        2(G5 + d2·G3 + d1·G2 + d0·G1)  =  2Φ − e·(d2·e² + 3eS + 3R²)
                                 │
              ┌──────────────────┴──────────────────┐
        e = 0 │                                    │ e ≠ 0
              ▼                                    ▼
   [L1] 2Φ ∈ ideal ⟹ Φ = 0.              [P4] e | Φ ⟹ e = γ·t^a·(sqfree | q),
        CONTRADICTION. **ONE LINE.**            b_i ∈ {0,1};  sub2: deg e = 10
        C10 / subcase:dm1 EMPTY                       │
        (char ≠ 2)                                    │
                                                      ▼
   [P5] the d3-killing shift + its index-(−1) triangularity  D̃₋₁ = D₋₁
   [P6] the two slice families   (P<) t^(2n−2) | [uⁿ]H², = 0 for n ≥ 9
                                 (Q)  t^(2n−3) | [uⁿ]H³
                                 stacked: 2H³−3H² = −1+3K²+2K³  (h_n cancels)
                                                      │
                        ┌─────────────────────────────┴──────────────┐
                        ▼                                            ▼
   [P7] THE CASCADE  v_t(h_k) ≥ 2k−1 (k ≤ 5)          [P8] the dictionary
        ⟹ **a_t ≥ 9**                                      e = h_5, R = h_6+¼h_1h_5,
        + level 12 ⟹ v_t(h_6) ≥ 11                         S = h_7+½h_1h_6+1/16h_1²h_5,
                        │                                   d2 = h_2−⅜h_1²
                        └──────────────┬─────────────────────┘
                                       ▼
                        [P9] bracket collapse
                             B = h_2h_5² + 3h_5h_7 + 3h_1h_5h_6 + 3h_6²
                                       ▼
                        [L2] **a_t ≤ 9**   ⟵ FIVE INTEGERS (§2)
                                       ▼
                            **a_t = 9 EXACTLY**
                                       ▼
   [P10] a=9 reduction; cofactor identity  F·Z = (1/6)γ⁵t⁹Π⁴;  Π | B;
         gcd(Z,Π) = 1;  Z = ζ·t^z with 2 ≤ z ≤ 6
                                       ▼
   [L3] marked-support test: k=1,2,3 infeasible ∀z; k=4 only z=3   ← the ONE Gröbner step
   [L4] degree ledger: k=4 needs deg F = 22 > 17 ✗;  k=0 dies by the dichotomy
                                       ▼
              ** all five a9_b*_T1 cells EMPTY — FRONTIER EMPTY **
```

Everything not on that diagram is scaffolding. The repo has **227 `.py` files and
165 `.md` files**; the spine above touches **twelve** of them.

**Honest page estimate: 15–20 pages, with 9 machine steps, of which exactly ONE is
not hand-checkable in principle.** See §3.

**Theorem candidates: (a) PROVED. (b) PROVED, but as a *sharpening that refutes the
repo's own framing*. (c) NEGATIVE — the sandwich has one jaw.** See §4–§6.

---

## 1. Part 1 — ablation results

Method: drop an input, recompute the conclusion from a *re-implementation* (not the
target's own code), and record the range of values over which the conclusion
survives. Every ablation below is in `minimal_core.py` with a mutation control.

### 1.1 `deg e = 10` (sub2) — three of four caps are not at the margin

`divisor_syzygy` C5 derives `deg e = 10` from four certified sub2 caps. Which bind?

| cap | audited | conclusion `deg e = 10` holds for | slack down | slack up | verdict |
|---|---|---|---|---|---|
| `deg R ≤ 12` | 12 | `R ∈ [0,12]` | 12 | **0** | **THE binding cap** |
| `deg S ≤ 14` | 14 | `S ∈ [0,15]` | 14 | 1 | not at the margin |
| `deg d2 ≤ 4` | 4 | `d2 ∈ [0,6]` | 14 | 2 | not at the margin |
| `deg e ≤ 10` | 10 | `e ∈ [10,39]` | **0** | 29 | used only as the ceiling `E ≤ 10` |

At `E = 10` all three RHS branches `d2+2E`, `E+S`, `2R` equal **24 simultaneously** —
a three-way tie. That is why only the maximum matters and only one cap binds.
`R = 13` already re-admits `E = 8`. (`A1.2`–`A1.8`; mutation control `A1.7`.)

*Compression consequence:* a writeup needs `deg R ≤ 12` sharply and the other two
caps only to within 1–2. That removes a real burden, because the cap derivation is
one of the longer upstream pieces.

### 1.2 `a_t ≤ 9` — five integers

This is the largest compression found. Given the K-syzygy, `v_t(Φ) = 30` exactly,
the dictionary, and the bracket collapse, the syzygy reads `a_t + v_t(B) = 30`
**as an equality**. So `a_t` is refuted whenever *every* term of `B` has valuation
`> 30 − a_t`. That is four integer comparisons:

| term of `B` | valuation floor | at a=9 (bar 21) | at a=10 (bar 20) | at a=10, no level 12 |
|---|---|---|---|---|
| `h_2·h_5²` | `v₂ + 2a` | 21 | 23 | 23 |
| `3·h_5·h_7` | `a + v₇` | 21 | 22 | 21 |
| `3·h_1·h_5·h_6` | `v₁ + a + v₆` | 21 | 22 | 21 |
| `3·h_6²` | `2·v₆` | 22 | 22 | **20 ← survives** |

* **a = 9 survives**, and the binding term at the boundary is `h_2·h_5²` at exactly
  `21 = 30 − 9` — **not** `h_7`. (`A2.1`, non-vacuity.)
* **a = 10 dies with margin 2**, and every `a ≥ 10` up to 30 dies. (`A2.0`, `A2.2`.)
* **Cascade level 12 is load-bearing for exactly one term**, `3·h_6²`, and with
  **zero margin** — the admissible floor is precisely `v_t(h_6) ≥ 11`. This
  sharpens `SYZYGY_COLLISION` X11 from "load-bearing" to "load-bearing with zero
  margin". (`A2.3`, `A2.4`, `A2.10`.)

Per-input slack for the *full* conclusion ("a = 9 survives **and** every a ≥ 10 dies"):

| input | audited | admissible | verdict |
|---|---|---|---|
| `v_t(Φ)` | `= 30` | `{30, 31}` | **zero downward slack**, slack 1 up |
| `v_t(h_6) ≥ 11` (level 12) | 11 | `≥ 11` | **zero slack** |
| `v_t(h_3) + v_t(h_4)` | `5+7 = 12` | `≥ 11` | joint slack 1 |
| `v_t(h_2) ≥ 3` | 3 | `≥ 1` | slack 2 |
| `v_t(h_1) ≥ 1` | 1 | **anything** | **NOT NEEDED** |

`v_t(h_1) ≥ 1` is consumed by nothing in the collision — it is needed only to
*start* the cascade. `h_3` and `h_4` enter only through their sum, via `v_t(h_7)`.
(`A2.9`–`A2.14`.)

**Two "not needed" claims independently reconfirmed.** `v_t(h_7) ≥ 11` follows from
the `(P<)` convolution floor on the audited rows alone, with **no cascade level**.
Hence levels 14 and 16 are not needed, and `r_13 = 0` is not needed — matching
`AT_LE9_AUDIT` E9 and F8 by a different route. (`A2.5`, `A2.6`.)

> **A wrong guess, caught by mutation testing.** My first draft of `A2.9` asserted
> that `v_t(Φ) = 31` would revive `a = 10`. It does not — the kill has margin 2, so
> it tolerates `+1`. The check FAILED, and the truth (`{30,31}`, zero slack
> *downward* only) is the finding. Recorded because the discipline earned its keep.

### 1.3 What is scaffolding

Confirmed non-load-bearing, from the map plus these ablations:

* **`g4_row.py` (68/68) — kills nothing anywhere**, and §6 *proves* no counting test
  on the raw λ row can ever kill. Elimination-only. Not in the kill chain.
* **`slice_phi_yplace`'s `a_t ≤ 9`** — not established as written (shifted-chart
  `(P<)` is false on the repo's own control, `AT_LE9_AUDIT` C-2/F11). Not a second leg.
* **Cascade levels 14 and 16**; **`r_13 = 0`**; **`v_t(h_1) ≥ 1` for the collision**.
* **Level 12 for the five-cell closure** — provably inert there (`sub1_spine9` X1b,
  `spine9_audit` G4/G5), even though it is zero-margin load-bearing for `a_t ≤ 9`.
  The same input is critical in one place and inert in the next.
* **`SLICE_OBSTRUCTION` §8's whole unused inventory**: the y-order conditions, the
  degree/upper-hull conditions, the exact rows `r_n = 0`, the sharp relation
  `[u^17]H³ = (1/6630)t^30 q`, and the sub2 degree caps. Never spent.
* **The T2-only consequences** `R | e²`, `e·R | Φ`, `R = c(y+1)^ρ`, and SPINE's
  zero-slack degree count — all sub2/T2 accidents, and the header of
  `slice_obstruction_basis` already says "nothing below uses any of them".
* **`divisor_syzygy` C6/C7's 8 → 3 column collapse** — superseded by the later,
  stronger closure.

### 1.4 Zero-margin items a writeup must carry as such

Not compression, but the honest inverse. Five places have **no margin at all**:

1. `v_t(S) ≥ 11`, binding on `(1/16)h_1²h_5` at `2+9` — not on `h_7` (`at_le9_audit` D14).
2. `v_t(h_6) ≥ 11` (level 12) for `a_t ≤ 9` (this document, `A2.10`).
3. `v_t(Φ) = 30` on the low side (this document, `A2.9`).
4. `deg d2 ≤ 6` and `z ≤ 6`, both for `a9_b0000_T1` (`spine9_audit` F12/F13).
5. `deg R ≤ 12` for `deg e = 10` (this document, `A1.2`).

---

## 2. Is it short enough to write as a human-checkable proof?

**Yes.** Estimate:

| section | pages |
|---|---|
| Setup, notation, and the upstream reduction as a black box | 2 |
| Φ: the forcing ODE, uniqueness, `q` squarefree, `v_t(Φ) = 30` | 1.5 |
| The K-syzygy; the `e = 0` branch (one line); `e | Φ` | 1.5 |
| The d3-killing shift and the two slice families | 3 |
| The cascade `v_t(h_k) ≥ 2k−1`, and level 12 | 3 |
| The dictionary, the bracket collapse, `a_t ≤ 9` | 1.5 |
| The `a = 9` reduction and the cofactor identity | 2 |
| The five-cell closure: support test + degree ledger | 2.5 |
| **total** | **≈ 17** |

Call it **15–20 pages** plus a computational appendix.

**The nine steps that must stay machine-verified:**

| # | step | hand-checkable? |
|---|---|---|
| 1 | the K-syzygy expansion | yes, with patience (~30 terms) |
| 2 | the forcing-ODE solution and `q`'s irreducibility | yes (linear algebra, deg ≤ 14) |
| 3 | the d3-shift triangularity `D̃₋₁ = D₋₁` | yes, but index-specific and error-prone |
| 4 | slice-family cokernel ranks (`dim coker = 2n−3`) | tedious but yes |
| 5 | the cascade's per-level lowest jets to level 12 | tedious but yes |
| 6 | the bracket collapse (`−3/8 + 3/16 + 3/16 = 0`) | yes, trivially |
| 7 | the `a=9` exact divisions and `F·Z = (1/6)γ⁵t⁹Π⁴` | yes, with patience |
| 8 | **the marked-support feasibility test** (40 `(k,z)` pairs) | **NO — genuine Gröbner / number-field work** |
| 9 | the degree ledger | yes, trivially |

Only **#8** is irreducibly machine work, and even there `spine9_audit` E4 gives a
hand gcd route for `k = 1` and E7 shows the `k = 3` kill is structural rather than
arithmetic. Steps 3, 4, 5 are hand-checkable but are exactly where this repo's
historical errors lived, so they should stay gated regardless.

**Why the core is not smaller than ~15 pages.** The cascade (§P7) and the shift
(§P5) are ~6 pages together and are irreducibly inductive: the cascade is a
level-by-level argument whose *ordering* is load-bearing (`slice_obstruction_audit`
F7: dropping level 8 makes level 10 branch), and the shift's triangularity holds at
index −1 but demonstrably **fails** at index −2, so it cannot be stated as a clean
change of coordinates. No amount of compression removes those.

**Conditionality that survives into any writeup** — unchanged by this lane:
(i) GGHV Prop 4.3 L1132's uncited 7th-power edge claim (shared with Helali);
(ii) `[QQ1]`/α-strip, spent at its boundary `M = −5`;
(iii) the zero-margin list in §1.4;
(iv) the top `subcase → C0` DAG edge remains judgment-level.

---

## 3. Candidate (a) — the Catalan law: **PROVED**

`BELYI_PASSPORT.md` N1 states the law as INFERRED on four data points and says
exactly what was missing:

> "A proof would presumably come from the dessin being a plane trivalent tree-like
> map with one big face; I did not construct that bijection."

**That bijection is constructed here.**

### THEOREM
For odd `k`, put `m = (k−1)/2`. The genus-0 dessins with passport
`( 2^((3k−1)/2), 1 | 3^k | (5k−1)/2, 1^((k+1)/2) )` are in canonical bijection with
**rooted plane binary trees on `m` internal nodes**. Hence there are exactly
`C_m` of them, every one is connected, and every one has trivial automorphism group.

### Proof
Write `n = 3k = 6m+3`. Read the triple as a map: `σ₂` (type `3^k`) gives `k` trivalent
vertices, `σ₁` (type `2^((3k−1)/2) 1`) gives `(3k−1)/2` edges and **one leg** (its
unique fixed point), and the faces are the cycles of `σ₃`.

*Faces of degree 1 are self-loops.* `σ₁σ₂(d) = d` means `σ₂(d) = σ₁(d)`, i.e. `d`
and `σ₂(d)` are joined by an edge — a loop at a trivalent vertex whose two darts are
consecutive. The passport prescribes `(k+1)/2 = m+1` such monogons and one big face.

*Delete the monogons.* Turn the leg into an edge to a phantom univalent vertex, so
`V = 2m+2`, `E = 3m+2`, `F = m+2`, and Euler gives `2` — genus 0. Deleting the `m+1`
loop edges cannot disconnect, and leaves `V = 2m+2`, `E = 2m+1`, so `F = 1` and
`V − E = 1`: a **plane tree**. Its degrees are `m` vertices of degree 3 (the
non-loop-carrying `σ₂`-vertices) and `m+2` of degree 1, since
`(m+1)·1 + m·3 + 1·1 = 4m+2 = 2E`.

*Root it.* The leg is canonical (`σ₁` has exactly one fixed point), so the tree is
rooted at a leaf. Rooted plane trees with `m` degree-3 vertices and `m+2` leaves are
exactly rooted plane binary trees with `m` internal nodes: **`C_m` of them.**

*Invert.* Attach a loop at each non-root leaf and reopen the root edge into a leg.
The two possible cyclic orders `(a,b,c)` and `(a,c,b)` at such a leaf are conjugate
by the transposition of the two loop darts `a,b`, which is a dessin isomorphism, so
the inverse is well defined and the third partition's shape is forced:
`3k − (m+1) = (5k−1)/2`. ∎

### Two corollaries, and they REPAIR real gaps

`BELYI_PASSPORT.md` proves connectedness and rigidity **only at k = 7**.

**Corollary 1 (transitivity is automatic, every k).** An orbit not containing the
long cycle consists of `σ₃`-fixed points, so `σ₁σ₂ = id` there, so `σ₂ = σ₁⁻¹`. An
element of order dividing both 2 and 3 is trivial; but `σ₂` has no fixed points.
Hence the orbit is empty. *Two lines, no enumeration* — the repo needed an
exhaustive `S₃` search to kill the `(18,3)` pattern at `k = 7`. (`C1`, with a
discriminating mutation control `C2`.)

**Corollary 2 (`|Aut| = 1`, every k).** The centralizer of a transitive group acts
semiregularly. Any automorphism permutes `Fix(σ₁)`, which is a single point, so it
fixes that point, so it is trivial. (`C3`, and `D2` brute-forces `|Aut| = 1` on all
`1+1+2+5+14 = 23` dessins with `m ≤ 4`.)

> **Corollary 2 closes a gap that was load-bearing for the statement, not
> cosmetic.** The repo computes the Hurwitz number as `N / n!`, which is the dessin
> count *only* if `|Aut| = 1`. Its rigidity argument bounds `|Aut|` by
> `gcd(3k, (5k−1)/2)` — which is **3, not 1, at k = 5, 11, 17, 23, …** (`C4`). So at
> `k = 5` and `k = 11` — two of the seven data points the law was inferred from —
> the repo had no proof that the number it computed was the number of dessins.

### Verification status

| claim | grade | evidence |
|---|---|---|
| the bijection is well defined and injective | **PROVED** + CHECKED | `D1(m=0..7)`: every one of `1,1,2,5,14,42,132,429` trees yields a dessin with the **exact** passport, transitive, pairwise non-isomorphic under a canonicalizer verified sound (`D4`, 788 random relabellings) and discriminating (`D5`); mutation control `D3` |
| the count is `C_m` | **PROVED** | the theorem above |
| independent confirmation | **CHECKED** | Frobenius character sums, no trees: `k = 1…19` gives `1,1,2,5,14,42,132,429,1430,4862 = C_0…C_9`. Calibrated (`E0`) and negative-controlled (`E3`: a perturbed passport gives 15, not 5) |
| transitivity, `|Aut| = 1`, all `k` | **PROVED** | Corollaries 1–2 |

The law now holds **to k = 19 (degree 57)** by machine and **for all odd k** by
proof, against the four data points it was inferred from.

### Honest caveat on novelty
The proof is a standard plane-tree bijection: monogon deletion on a genus-0
trivalent constellation. In the dessins/constellations literature this is a routine
technique, and "Catalan numbers count rooted plane binary trees" is classical. So
this is a **clean and complete result, but almost certainly folklore-adjacent, and I
would not present it as a contribution to the Hurwitz-number literature.** Its real
value here is threefold: it upgrades the repo's own INFERRED headline to PROVED, it
closes the `|Aut|` gap described above, and it explains *why* the number is Catalan —
which is what makes the `35 = 7 × 5` decomposition intelligible rather than
coincidental.

### The VDIM ladder
`vdim = k × Hurwitz` is a **separate** law and remains INFERRED. Provenance of the
"measured at k = 1,3,5,7" claim, checked: `k = 1 → 1` is **trivial** (`build_block(1)`
returns an empty variable list, so `vdim = 1` in a zero-variable ring);
`k = 3 → 3` and `k = 5 → 10` are the only two genuinely independent Singular
measurements; `k = 7 → 35` needs either `--full` Singular or the external
`firstblock_Q_exact.out` bundle. So the honest count is **two nontrivial independent
data points**, not four. An attempt to extend the ladder to `k = 9`
(predicted `9 × 14 = 126`) via `modStd` on the 8-variable block is in
`_scratch_vdim9.py`. Its block is well formed (8 residuals `R[14..21]` in 8
unknowns, with `R[22] ≡ 0` — the exact analogue of the vacuous `R[17]` at `k = 7`),
but the Gröbner run **had not returned after ~50 minutes** and was still going when
this was written; see `_scratch_vdim9.log` for the outcome. Reported as an **open
item, not a result** — and note that a timeout is never a verdict. (With the Catalan
law now proved, the VDIM law is equivalent to the purely algebraic statement
`vdim = k·C_{(k−1)/2}`.)

---

## 4. Candidate (b) — why is the K-syzygy productive? **PROVED, and it refutes the repo's framing**

The brief asked whether the syzygy is "a feature of sporadic corners specifically."
**It is not, and sporadicity is the wrong axis entirely.**

`weight_lemma_75_125` explains the (72,108)/(75,125) difference by `q_window = 1`,
where `q_window` is the denominator of `ord_y(Φ)/M`: the y-order cap
`L(w) = ⌈αw/q⌉` is superadditive off the period, and at (75,125)
`L(6) + L(30) = 202 > 201 = L(36)`, so no relation `c·Φ = e·B` can exist. At
(72,108), `q_window = 1` so the carry vanishes identically.

**That is sufficient but not necessary.** Since `q | αM` (which is what it means for
`q` to be the denominator of `ord_y(Φ)/M`) and `gcd(α,q) = 1`, the two residues
`(−αw_e mod q)` and `(αw_e mod q)` sum to `q` or to `0`. Therefore:

> **CRITERION.** The y-order carry obstruction to a Φ-divisor relation vanishes
> **iff `q_window | w(e)`**. The carry is exactly 0 or 1, never more.

Verified exhaustively over **1 529 566** admissible splits with **zero
mismatches** (`F1`), and mutation-controlled: **30 326** admissible splits have
carry 0 with `q_window ≠ 1` (`F2`). At (75,125) the obstruction would vanish for
`w(e) ∈ {12, 24}` — exactly the multiples of `q_window` (`F3`). At (72,108),
`q_window = 1` divides every weight, so all 16 splits of `M = 17` are unobstructed
(`F4`).

**Consequences.**

* The correct statement is "(72,108) has the syzygy because `q_window | w(e)`",
  which `q_window = 1` implies but does not exhaust. Other corners can support the
  mechanism; the search space is larger than the repo believes.
* **Sporadicity is irrelevant.** `q_window | w(e)` is an arithmetic condition on one
  corner's `(ord_y Φ, M, w(e))`. Nothing connects it to membership in a GGV5
  `F_j` family. (72,108) *is* sporadic — it appears in no `F_j`, only in the
  length-1 chain `(8,28)|(11/4,7)|(3,2)|108` and the length-2 chain
  `(9,27)|(9,24)|(11/3,8)|(2,3)|108` — but that is a separate fact, and the two are
  logically unrelated.
* The repo's `q_window ≠ 1` family law still correctly rules out the whole F2
  family; it just rules it out for a *stronger* reason than stated. **REPAIRED
  2026-07-26:** the law is `q_window = 12a − 7 ∈ {17,29,41,…}`, not `5a − 3`, and
  `w(e) = t+1 = 5`, not 6 — both were pre-repair `(5,20)` chart values (`2adb92a`).
  The conclusion is now unconditional rather than accidental: every period in the
  family exceeds `w(e)`, so `w(e)` cannot be a nonzero multiple of any of them.

**Why C10 was missed, stated precisely** (this part of the brief is confirmed and is
worth keeping verbatim in a writeup): `e | Φ` silently presumes `e ≠ 0`, so the
syzygy covers **both** branches — one by division, one by contradiction. The K-syzygy
was already in the repo three times before anyone connected it to the `d₋₁ = 0` leaf.
The leaf's field scope is `char ≠ 2`, not char 0.

**Negative I could not overcome.** A sweep of `q_window` across all 34 GGV5
candidate corners is **not possible with the data in this repo**:
`paper_src/upstream_facts.json` carries Newton polygons for the **(8,28) corner
only**, so `ord_y(Φ)` and `M` are unavailable for the other 33. Testing the
criterion family-wide requires re-deriving each corner's polygon from GGV5 — work
never done here. I flag this rather than guess.

**Framing warning, from `LANDSCAPE_2026_07.md` §3.4 and endorsed here:** do not use
"syzygy" as a novelty hook. `1406.0886_GGV3.tex:1086-1109` applies the same
technique to a different relation, upstream of Abhyankar–Moh 1975. Scope any claim
to *this* relation and *this* conclusion, and lead with the `q_window | w(e)`
criterion, which is genuinely new and is a *classification* statement.

---

## 5. Candidate (c) — the two-sided sandwich: **NEGATIVE. It has one jaw.**

`LANDSCAPE_2026_07.md` §3.1 claims the novelty is *constructive pinning* (fixing
`a_t = 9` exactly) as against the literature's use of two-sided bounds only for
*contradiction*. `AT_LE9_AUDIT` C-1 already weakened this by showing the two proofs
of `a_t ≤ 9` share all four equations. **The ablation shows something stronger and
worse.**

The two *sides* of the sandwich are not independent either. Stripping each cascade
row and re-testing the upper bound:

| strip | `a_t ≤ 9` |
|---|---|
| `v_t(h_1) ≥ 1` | SURVIVES |
| `v_t(h_2) ≥ 3` | **FAILS** |
| `v_t(h_3) ≥ 5` | **FAILS** |
| `v_t(h_4) ≥ 7` | **FAILS** |
| level 12, `v_t(h_6) ≥ 11` | **FAILS** |

So the **upper** bound consumes three of the four cascade rows that produce the
**lower** bound, plus the cascade's level 12. There is no pincer. There is *one*
valuation cascade, read twice: once directly as a floor on `v_t(h_5) = v_t(e)`, and
once — through the K-syzygy, which converts `v_t(Φ)` into the **equality**
`a_t + v_t(B) = 30` — as a ceiling.

**Assessment.** The pin `a_t = 9` is real and is genuinely finer than anything in
the four manuscripts. But "two-sided sandwich" is a *misdescription of our own
method*, and it should not be sold as a reusable two-source technique. The reusable
unit is a **pair**:

> a valuation cascade that floors the coefficients, **plus** a divisor identity that
> turns a known exact valuation into an equation.

and the second half is available exactly when `q_window | w(e)` (§4). That is a
falsifiable statement about reusability, which "the sandwich is a new technique" is
not. Note also that `LANDSCAPE_2026_07.md` §4.3's own list of "still uncontested and
worth stating plainly" **does not include the sandwich** — the document already
declined to defend it.

Two prior-art caveats stand unresolved and both sit in this lane: **Makar-Limanov
2025** (Serdica 51, "On the shape of a counterexample to the two-dimensional
Jacobian conjecture") is Newton-polygon work from six months ago, metadata-verified
but **full text never obtained**; and **Borisov** (arXiv:1901.04073) does
combinatorial enumeration of necessary conditions for 2D Keller maps at degrees
(99,66), and `grep -ril borisov` over this repo returns nothing.

---

## 6. Summary — PROVED / CHECKED / INFERRED

| item | grade |
|---|---|
| `Hurwitz(passport_k) = C_{(k−1)/2}` for all odd k, via rooted plane binary trees | **PROVED** |
| transitivity automatic, all k | **PROVED** |
| `|Aut| = 1`, all k (repo had k=7 only; its argument fails at k = 5, 11, 17, 23) | **PROVED** |
| the carry obstruction vanishes iff `q_window \| w(e)` | **PROVED** |
| `q_window = 1` is sufficient, NOT necessary (30 326 counter-splits) | **PROVED** |
| the Catalan law holds to k = 19 by character sums | CHECKED |
| the bijection is exact for m ≤ 7 (429 dessins, canonicalizer sound + discriminating) | CHECKED |
| `deg e = 10` binds on one cap only, `deg R ≤ 12`, zero slack | CHECKED |
| `a_t ≤ 9` consumes 5 numbers; level 12 zero-margin; `v_t(h_1)` unused | CHECKED |
| `v_t(Φ) = 30` has zero downward slack, slack 1 upward | CHECKED |
| the sandwich has one jaw (upper bound consumes 3 of 4 cascade rows) | CHECKED |
| the writeup is 15–20 pages with 9 machine steps, 1 irreducible | **INFERRED** (my estimate) |
| `vdim = k × Hurwitz` | INFERRED — 2 nontrivial data points, not 4; k=9 inconclusive |
| `q_window \| w(e)` across the other 33 GGV5 corners | **UNTESTABLE with repo data** |

### Negatives, stated plainly

1. **The sandwich does not generalize as advertised** — it is not two independent
   bounds, it is one cascade plus one divisor identity (§5).
2. **The syzygy is not a sporadic-corner phenomenon.** `q_window | w(e)` has nothing
   to do with GGV5 family membership (§4).
3. **The Catalan proof is probably folklore.** Complete and gap-closing, but a
   routine bijection; do not oversell it (§3).
4. **The VDIM ladder rests on two data points, not four**, and my attempt to add a
   third (k = 9) had not returned after ~50 min — inconclusive, not negative (§3).
5. **The core is not smaller than ~15 pages**, because the cascade and the shift are
   irreducibly inductive and their ordering is load-bearing (§2).
6. **One of my own checks was wrong** about `v_t(Φ) = 31` and mutation testing caught
   it (§1.2).
