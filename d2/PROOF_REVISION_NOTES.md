# `PROOF_72_108.md` — revision 2 notes

**Date 2026-07-26. Base: the committed draft at `719e040` (837 lines, 15 735 words). Result: 992 lines, 19 746 words in the body.**

This file records exactly what changed and why, so that a reader who has the previous
draft can diff the *argument* rather than the prose. Two changes were mandated by an
external review; two more were opportunistic and are marked as such. Nothing in the
exclusion argument changed. Two repository claims were checked and found wrong; both
are recorded rather than quietly dropped.

---

## 1. The germ / admissible-germ split (DEFECT 1 — the blocking one)

### The defect

The previous draft defined a **germ** as a solution of `G1 = G2 = G3 = G5 = 0` in
`K[y]^7` with the fixed `Phi`, then stated "Theorem A′: there is no germ". But the
proof also consumes:

* the `P`- and `Q`-slice divisibility families (old §2.5),
* the *unshifted* coefficients `h_i` and the shift dictionary `h_5 = e` (old §§2.3, 7.1),
* one of the two configuration-specific degree cap profiles (old §2.5 table).

None of those follows from the four equations. So the stated theorem was **strictly
stronger than what was proved**. This is a statement-scope defect, not a hole: every
hypothetical counterexample does produce all that structure.

### How the split was made

**New §4 "Two objects: G-points and admissible germs"** now occupies the section slot
vacated by the demoted divisor-consequences section (see §2 below). It contains:

* **Definition 4.1 (`G`-point, a.k.a. algebraic germ)** — the bare four-equation object,
  with an explicit "nothing else is assumed" clause.
* **Definition 4.2 (admissible germ)** — the tuple
  `(d_2, d_1, d_0, e, R, S, T, H, c)` with `H(u) = sum h_k u^k`, `h_0 = 1`,
  `c in {(1),(2)}`, subject to four labelled hypotheses:
  * **(A1)** the four `G`-equations;
  * **(A2)** the slice divisibilities `t^{2n-2} | p_n` (n = 2..8), `p_n = 0` (n >= 9),
    `t^{2n-3} | r_n` (n = 2..15);
  * **(A3)** the shift dictionary (7.1.1), in particular `h_5 = e`;
  * **(A4)** a cap profile `deg d_j <= lambda(c) * w(d_j)`, `lambda = 3` / `2`.
* **Theorem 4.3 (transfer)** — a counterexample pair in either configuration determines
  an admissible germ. Its proof discharges (A1)–(A4) one at a time, naming which earlier
  subsection supplies each: (A1) from §§2.1–2.4 plus Lemma 2.4's *equivalence*; (A2) from
  Proposition 2.1 plus `P = C^2`; (A3) from the triangularity identity (2.3.1) plus
  premise `[I3]`; (A4) from Lemma 2.5 plus its shift-preservation clause.
* **Theorem A″: there is no admissible germ**, with the top-level logic boxed as
  `counterexample => admissible germ => contradiction`.
* **§4.3, a per-result ledger** stating for every lemma in the paper which object it is
  about and which of (A1)–(A4) it consumes.

### What the ledger revealed

Three facts that were invisible while everything was "a germ":

1. **Lemma 3.3 and Theorem 3.4 need only (A1).** The two results that structure the
   whole argument have the weakest hypotheses in the paper.
2. **(A4) is consumed by exactly two of the five cases** (`k = 0` and `k = 4`) plus §9.
   Three of the five kills are cap-free.
3. **§5 is about `H` alone** — no `G`-equation appears in it at all.

Also newly demoted-to-by-product, **Lemma 11.5 / Theorem 11.6 / Corollary 11.7 are the
only results in the paper that hold on a bare `G`-point**, and Corollary 11.7's
`a_t <= 10` is the only bound needing none of (A2), (A3), (A4).

### Consequential edits

* Old "Theorem A′" is **retired**. §4.2 and §13.4 state explicitly that "there is no
  germ" is *not* proved here and should not be claimed; §14 item 11 records
  "is `V(G1,G2,G3,G5)` empty?" as an **open** question, noting §10's quartic cover as the
  obvious tool since it too lives on a bare `G`-point.
* Every theorem statement retyped: `On every germ` -> `On every admissible germ`
  (Theorems 6.2, 7.2, 7.3, 8.6, 9.1, Corollary 8.3, Theorem 8.1's conclusion), or
  -> `On every G-point` (Lemma 3.3, Theorem 3.4, §11.3, §10's quartic cover).
* §1's skeleton now shows the two-step through the intermediate object; §1.1 and the
  header abstract updated.
* §2.4's old one-sentence definition of "germ" replaced by a forward pointer plus an
  explicit **"do not read `G`-point as germ"** warning.

---

## 2. Removing the old §4 from the spine (COMPRESSION 1 — verified before use)

### What moved

The previous §4 ("Consequences of the divisor structure") contained Lemma 4.1
(`e | S`, by Sylvester resultant + adjugate cofactors + integral closure), Theorem 4.2
(the place trichotomy), and Corollary 4.3 (`t^a | R,S,T` for `a <= 10`, plus a cap-free
`a_t <= 10`). Corollary 4.3 supplied §8's ansatz.

**All three are now §11.3**, renumbered Lemma 11.5, Theorem 11.6, Corollary 11.7. §11's
title changed from "Two by-products" to "Three by-products".

### What replaced it: Lemma 7.4 (new §7.5)

> `v_t(d_2) >= 2`, `v_t(d_1) >= 3`, `v_t(e) = 9`, `v_t(R) >= 10`, `v_t(S) >= 11`,
> `v_t(T) >= 12`; in particular `t^9 | R, S, T`.

Proof: four `min`s over the dictionary (7.1.1) evaluated on the cascade profile (6.2.1),
then one line of `G1 = 0`:

```
v_t(d_2) >= min( v(h_2), 2v(h_1) )                      = min(3, 2)          = 2
v_t(d_1) >= min( v(h_3), v(h_1)+v(h_2), 3v(h_1) )       = min(5, 4, 3)       = 3
v_t(R)   >= min( v(h_6), v(h_1)+v(h_5) )                = min(11, 1+9)       = 10
v_t(S)   >= min( v(h_7), v(h_1)+v(h_6), 2v(h_1)+v(h_5) )= min(12, 12, 11)    = 11
G1 = 0  =>  e*T = -(1/2) d_1 e^2 - d_2 e R - R S,  RHS terms at 21, 21, 21
         =>  v_t(eT) >= 21,  v_t(e) = 9  =>  v_t(T) >= 12
```

**Independently re-verified for this revision** (not taken from the review): over 200
random coefficient profiles obeying (6.2.1), the minima are exactly
`v_t(d_2) = 2`, `v_t(d_1) = 3`, `v_t(R) = 10`, `v_t(S) = 11`, `v_t(eT) = 21` — so the
bounds hold and are attained (hence sharp). Two findings beyond the review's brief:

* **The lemma does not need cascade level 12.** With only `v_t(h_{1..5}) >= (1,3,5,7,9)`
  and Lemma 6.1's consequences `v_t(h_6) >= 10`, `v_t(h_7) >= 11`, all four minima are
  *unchanged*. This independently corroborates `spine9_audit.py` **G4** / `sub1_spine9.py`
  **X1b** ("level 12 inert in §8"), and it means the zero-margin input of §7.3 buys
  nothing in §8.
* **The `v_t(T) >= 12` step is a three-way tie at 21**, so all four inputs are
  load-bearing — and the conclusion is identical under both readings of the `d`-ledger
  (`>= 5, >= 3` and `>= 3, >= 2`), which matters because the two readings disagree.

### What §8.4 now does

The old §8.4 hedged: "`2 <= z <= 6` under **every** reading of the shift convention",
with a parenthetical claiming that "the unshifted reading with level 12 gives the sharper
`z <= 5`". That hedge is gone. The window is now derived once, from Lemma 7.4's ledger:
`v_t(A) >= 1`, `v_t(v) >= 2`, `v_t(u) >= 2`, `v_t(w) >= 3`, hence `v_t(Z) >= 2` and
`v_t(F) >= 3`, hence `z` in `[2,6]`. The `z <= 5` claim is **not** repeated: with
`v_t(A) >= 2` the bound `v_t(F) >= min(2+2, 3) = 3` is unchanged, so it does not follow,
and `v_t(R) >= 11` is recorded as false at `spine9_audit.py` **G5**.

### What still depends on `e | S` — checked, and one review claim refuted

Named explicitly in §11.3's three closing notes and in §14 item 12:

| dependant | status |
|---|---|
| Theorem 11.6 (place trichotomy) | yes — it is proved *from* Lemma 11.5. Now a by-product itself. |
| §7.4(b), the alternate-regime corroboration (`alt_level12.py` L6.1, `alt_rebuild.py` `B4_trichotomy`) | yes, through horn 1. Corroboration only; Theorem 7.2 gives `{a >= 11}` empty uniformly without it. |
| §9 route (c) (`spine.py` **S19**/**S20**, the spare collapse `S = e*Sbar`) | yes. Routes (a) and (b) do not, and route (b) is machine-checked, so config (2) does not rest on it. |
| the older enumerative registry route (`g4_row.py`, `divisor_consequences.py`) | yes, but not on this spine (§1.5). |
| §10's quartic cover, *curve* form `v = S/e - u^2` | yes. The *polynomial* form `W^2 = R^4 + ...` with `W := eS - R^2` does not, and that is the form §10 states. |
| **Lemma 3.3 (the `e = 0` branch)** | **NO.** |

**The review's aside that `dm1_branch_verify` needed `e | S` is wrong, and I checked it
rather than repeating it.** `dm1_branch_verify.py` contains no `e|S`, no `S/e`, no
resultant and no integral-closure argument. Its **B1–B3** exist for the opposite reason:
to record that the post-elimination chain — the only place a division by `e` would enter
— is *unusable* on that branch, `sol4` having denominator exactly `2*dm1`, and that the
post-elimination pair `H2, H3` is blind there (both collapse to multiples of
`G1|_{dm1=0}/3 = dm2*dm3` and neither contains `Phi`). The branch is closed
**pre-elimination** by one integer identity (**C1**). This is now stated in §11.3 and
listed in §13.4 among claims the paper declines to repeat.

Also verified: **Lemma 7.4 does not use Corollary 11.7 either.** The compression is total
— after it, no Sylvester resultant, no integral-closure argument, and no place-by-place
valuation dichotomy appears anywhere on the spine. §10's "what the proof does not use"
list gained a bullet saying exactly that.

---

## 3. `[QQ1]` became an ordinary proposition — COMPLETED

`[QQ1]` is now **Proposition 2.1 (the alpha-strip normalisation)**, stated in §2.1 with
explicit hypotheses (the `(1,0)`-grading, `v(P) = 8`, `v(Q) = 12`, `l(P) = R^2`,
`l(Q) = R^3` with `R` primitive, `P = C^2`) and a five-step proof:

1. `v([P,Q]) = 2 < 19 = v(P)+v(Q)-1`, so **GGV1 Prop 1.13**'s equality clause forces
   `[lP, lQ] = 0`.
2. **GGV1 Prop 2.1** (aligned) with `R` primitive: anything commuting with `lP` is
   `alpha_k R^k`. Alignment `(m,n) = (2,3)` is forced by `3*8 = 2*12`.
3. `[C^k, C^2] = 0`, so stripping preserves the bracket; strippable levels are `v = 4k`,
   admissible range `k in {-1,0,1,2}`.
4. Termination: `[lF, lP] != 0` puts Prop 1.13 in its **equality** case, giving
   `2 = v(F) + 8 - 1`, i.e. `v_{1,0}(F) = -5` **exactly** — and `-5` is not a multiple of
   `v_{1,0}(R) = 4`, so the halt is *forced*, not assumed.
5. Clearing `alpha_2` (subtract `alpha_2 P`), `alpha_0` (subtract `alpha_0`), and
   `alpha_1` (replace `P` by `P + (2/3)alpha_1`, whose cube absorbs `alpha_1 C` into the
   `C^{-1}` column, i.e. into `lambda`). `k = -1` is retained as `lambda C^{-1}`.

Cited by content, not by number where the risk is real: the two GGV1 propositions are
identified by their statements (the bracket-degree inequality *with* its equality-iff
clause; the alignment/pure-power statement) as well as their numbers. Per the standing
discipline, theorem *numbers* in un-compiled TeX are not verifiable, so the content is
what carries the citation.

§12.2 no longer restates the proof; it now discusses standing only, and adds the
observation that the halt and the exponent **corroborate each other** rather than being
independent assertions. Trust tier is unchanged at **2/4** and said so in both places.
The provenance row for `t6_premises_verify.py` was expanded from a bare `P2a–P2e` to all
ten check IDs with what each one verifies.

---

## 4. The cap lemma — COMPLETED

The old §2.5 said "the certified degree caps, **read from the engine at runtime**".
That is now **Lemma 2.5** in a new **§2.6**, proved in three ingredients:

1. **Direction functionals off the polygons.** A linear form maxes at a vertex, so the
   corner sets determine `max(j-i) = 8` (config 1), `max(j-2i) = 0` (config 2),
   `max(2i-j) = 2` (both), giving per-slice bounds `deg P_i <= M - a*i`,
   `ord P_i >= 2i - M` for `i = 0..8`.
2. **One valuation induction, closed in `k`.** From `P = C^2`,
   `C_{4-k} = (P_{8-k} - sum C_{4-j}C_{4-(k-j)}) / (2 C_4)`; `deg_y` and `-ord_y` are
   valuations and the pivot `2C_4` never vanishes. The product bound is `j`-free and
   equals the slice bound with **no slack**. Result: `deg C_{jx} <= jx+4` / `<= 2 jx`,
   `ord C_{jx} >= 2jx - 1`.
3. **The `D`-transform and the strip.** `deg D_{jx} <= 60-15jx` / `56-14jx`,
   `ord D_{jx} >= 48-12jx`; at `jx = 4-k` this is `deg <= 15k` / `14k` and
   `ord >= 12k` (the window floor). Subtracting: `lambda = 15-12 = 3` / `14-12 = 2`.

Plus a shift-preservation clause: `cap(m) + (m-j)*cap(3) = cap(j)` identically in
`(m,j)` in all three directions, because `cap(3)` *is* the per-step slope
(15 / 14 / 12) — which is why the caps apply to the *shifted* variables the G-system
actually uses.

The evaluated table now shows `w = (2,3,4,5,6,7,8)` and fills in the two entries the old
table left as `—` (config (2) `d_1 <= 6`, `d_0 <= 8`). The two **zero-margin** entries
`deg d_2 <= 6` and `deg R <= 12` are boldfaced with their consumers named. Two sensitivity
facts were added:

* **Attainment forbids lowering a cap.** A generic series square root attains
  `deg D_k = lambda*k` and `ord D_k = 12k` exactly at every `k = 1..17` in both
  configurations. In particular lowering `deg e` to 9 would empty config (2) outright —
  an unearned proof — and attainment rules it out.
* **Raising `deg R` is the dangerous direction**: it enters the config-(2) test as
  `2 deg R`, so an off-by-one costs two degrees.

`caps_audit.py` was **re-run for this revision: 70/70 pass**, so its provenance mark went
from `●` (last logged suite run) to `✓` (re-run here), and the row now lists the real
check IDs — `A1–A5`, `B` (direction-keyed, not numbered), `C1–C9`, `D1–D4`, `E1–E4` per
regime, `F1–F8`, `G1–G3`, `H1–H4`, `I1–I3`, `J1–J6`, `K0–K5`.

---

## 5. Other edits, all bookkeeping

* §13.2 grew from "the nine machine steps" to **eleven**, with the cap lemma (#4) and the
  valuation ledger (#8) inserted. Both are *replacements*: #4 replaces a runtime cap
  read, #8 replaces a Sylvester resultant plus a place trichotomy. Both move work from
  "machine says so" toward hand-checkable. Still exactly **one** irreducible step, now
  numbered #10. §1's count updated to match.
* §2's lemmas renumbered: Lemma 2.1 -> **2.2** (the ODE), Lemma 2.2 -> **2.3** (the exact
  valuation), Lemma 2.3 -> **2.4** (nothing is lost), plus new Proposition 2.1 and
  Lemma 2.5. All six internal references chased.
* §13.3 rows retargeted: `divisor_consequences.py` `3.1, 4.1` -> `3.1, 11.3`;
  `t1_branch.py` `4, 10` -> `10, 11.3 — off the spine since Lemma 7.4`;
  `alt_level12.py` `4, 7.4` -> `7.4, 11.3`.
* §13.4 gained three entries: "there is no germ" as a statement of what is proved; that
  the `e = 0` branch consumes `e | S`; and that the caps are runtime reads.
* §14 gained items 11 (is the bare variety empty? **open**) and 12 (§9 route (c) consumes
  Lemma 11.5).
* §11.3's naming note now records `alt_level12.py` **L5.1**: the "place trichotomy" is a
  *dichotomy*; the genuine three-way trichotomy is the unrelated `T1 | T2 | T3` split.

### Constraints honoured

* `a_t <= 9` is still presented as **single-legged** (§1.3 unchanged). `slice_phi_yplace`
  is still not cited as corroboration. `AT_LE9_AUDIT`'s `C-1`/`C-2` are still flagged as
  narrative labels with machine correlates `E2` and `F11`/`F12`.
* Still conditional on GGHV22 Prop 4.3 exhaustiveness and on `[QQ1]` — now
  Proposition 2.1 — spent at its `M = -5` boundary.
* Helali's priority (2026-07-21, doi:10.5281/zenodo.21479814, verdict **SUBSUMES**)
  untouched and still §0.1.
* `proof_dag.json`'s recorded `C0` level is still stated as **`claimed`** with the
  recorded-vs-assessed gap intact.
* Horn-1 exclusion still cited as `alt_level12.py` **L6.1** corroborated by
  `alt_rebuild.py` **B4_trichotomy**; no `alt_regime.py` is cited anywhere.

---

## 6. Page estimate

| | previous draft | revision 2 |
|---|---|---|
| lines | 837 | 992 |
| body words | 15 735 | 19 746 |
| lines opening a display (`^$$`) | 58 | 71 |
| table rows | 85 | 101 |
| **estimated pages** (11pt article, ~500 words/page incl. displays and tables) | **~31** | **~39** |

Per-section, revision 2:

| section | pages |
|---|---|
| §0 status / priority / conditionality | 3.1 |
| §1 introduction | 1.4 |
| §2 conventions, `[QQ1]`, the cap lemma | 4.7 |
| §3 the K-syzygy | 1.6 |
| **§4 the two objects (new)** | **2.0** |
| §5 the stacked obstruction | 1.1 |
| §6 the cascade, `a_t >= 9` | 1.2 |
| §7 the collapse, `a_t <= 9`, **the ansatz (new §7.5)** | 3.5 |
| §8 the five cases | 4.5 |
| §9 configuration (2) | 1.0 |
| §10 what the proof does not use | 1.1 |
| §11 three by-products (incl. **demoted §11.3**) | 4.7 |
| §12 the imported boundary | 1.6 |
| §13 machine verification | 6.5 |
| §14 open items | 1.0 |
| references | 0.6 |
| **total** | **~39** |

The **mathematical core** (§§1–9 with §0's conditionality folded to a paragraph) is
about **21 pages**, against `MINIMAL_CORE.md`'s standing estimate of 15–20. The
difference is §4 (2 pages, new and non-negotiable — it is the theorem statement) and the
larger §2 (Proposition 2.1 and Lemma 2.5 together add ~1.5 pages of proof that was
previously a citation and a runtime read). Net: **the spine got one section shorter and
the setup got two proofs longer.**

---

## 7. What resisted compression, reported honestly

* **§§5–6 did not shrink at all.** The cascade's ordering is load-bearing — the level-12
  jet factors into two coprime non-constant factors *before* reduction by relations
  already derived, and only becomes a perfect square after — so the deduction cannot be
  restated as a change of coordinates. The ~6 pages stand; the draft says so in §6.1 and
  I did not attempt to argue otherwise.
* **§13 is 6.5 pages, a sixth of the paper, and grew.** The provenance table's job is to
  be checkable, and every entry added this round (the cap-lemma IDs, the ten
  `t6_premises_verify` IDs) made it longer while making it more useful. I see no honest
  way to compress it that does not amount to asking the reader to trust a summary. It
  should probably become an appendix at typesetting time, but it should not be cut.
* **§11.1 (the Catalan law) is 2+ pages for a by-product the paper itself calls
  folklore-adjacent.** It stays because it repairs the `|Aut| = 1` gap that the Hurwitz
  count depended on, but a referee may reasonably ask for it to be an appendix or a
  separate note.
* **§1.4's zero-margin table shrank by two rows' worth of hedging** (the two caps are now
  theorems) but the other three margins are irreducible facts about the arithmetic and
  cannot be compressed away.
* **The transfer Theorem 4.3's proof is a paragraph of pointers, not a self-contained
  argument.** That is deliberate — each clause points at the subsection that does the
  work — but it means §4 is only as strong as §§2.1–2.5, and a formalisation would have
  to expand it. It is the right shape for a paper and the wrong shape for Lean.

---

## 8. What a reviewer should attack next

In descending order of leverage, given this revision:

1. **`[QQ1]` / Proposition 2.1 at its `M = -5` boundary.** Still tier 2/4, still spent
   exactly where the `F` column first becomes nonempty. Now at least it is a proposition
   one can attack step by step; step 4 (the equality case of GGV1 Prop 1.13) is where a
   failure would land.
2. **(A2) as a hypothesis of Definition 4.2.** The slice families are the single largest
   thing an admissible germ carries beyond the four equations, and §5.3's joint control
   has zero spares — so the level-10 step is controlled structurally but not by an
   instance with `e != 0`.
3. **The single-legged `a_t <= 9`.** Unchanged by this revision, and still item 1 of §14.
4. **Cascade level 12's zero margin in §7.3** — now known to be *inert* in §8 and in
   Lemma 7.4, which localises the risk precisely to §7.3 and nowhere else.
