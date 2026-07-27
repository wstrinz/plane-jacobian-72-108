# The Cross-Corner Atlas

**GGV5's 34 possible counterexamples with `max(deg P, deg Q) <= 150`, run through
the five gates this campaign proved.**

- **Registry + checker:** `corner_atlas.py` (`--quiet`, exit 0 iff all pass) — **31/31**
- **Machine-readable atlas:** `corner_atlas.json` (34 rows, one per published case)
- **The transfer experiment:** `corner_atlas_toric_34.py` (`--quiet`) — **9/9**, ~80 s
- **Population source:** `paper_src/1708.07936_GGV5.tex`, §"Possible counterexamples
  with `max(deg P, deg Q) <= 150`" (tex:1792). Every datum carries its line.

---

## 0. The two answers, up front

### Q1 — "four or five recurring mechanisms, or 34 singletons?"

**A handful. The 34 rows collapse to SIX gate signatures, and four of them cover
32 of the 34.** The frontier is attackable in batches, not case by case.

| # | signature `G1｜G2｜G3｜G4｜G5` | n | what the class *is* |
|---|---|---|---|
| 1 | `FAIL｜FAIL｜FAIL｜FAIL｜UNK` | **19** | monomial corner, `t != 4` |
| 2 | `FAIL｜PASS｜FAIL｜FAIL｜UNK` | **9** | monomial corner, `t = 4` — **the (75,125) class** |
| 3 | `PASS｜FAIL｜UNK ｜FAIL｜PASS` | 2 | retracting, `t = 3`, `q_window = 1` |
| 4 | `PASS｜FAIL｜UNK ｜FAIL｜UNK` | 2 | retracting, `t = 3`, `q_window != 1` |
| 5 | `PASS｜PASS｜PASS｜PASS｜PASS` | **1** | **(72,108) — the case we closed** |
| 6 | `PASS｜PASS｜UNK ｜FAIL｜UNK` | 1 | (8,28)/(3,4)/144 |

And the classification is driven by **one integer test**:

> **28 of the 34 rows fail the retraction shape `b_0 = t·(a_0 − 1)` at `A_0`.**
> At such a corner `C` is a monomial, so `deg C = ord C = 1`; that single fact
> simultaneously (a) refuses GGV5's final-corner dictionary and (b) empirically
> empties the Belyi sweep. **One test, two mechanisms here, 28 rows.**
>
> ### ⚠ CORRECTED 2026-07-27 — this claim used to say THREE mechanisms
>
> The withdrawn third was "forces `lam = 0`, killing the slice cascade with no
> case-specific input at all". **That was wrong**, and it was wrong in a way that
> inflated this very headline. Gate `G3` computed `lam = (deg_y Φ − ord_y Φ)/M` —
> the **strip** object — and tested it against `D5`, whose `lam` is the **cap**
> object (the D-transform slope difference; `contact_lemma.py:539/545`).
> `lambda_two_objects.py` (9/9) proves `cap ≥ strip`, so a strip-based `FAIL` is
> **not** a cap-based `FAIL`: the substitution could only ever declare our own
> mechanism void where it is in fact available. See `g3_gate_defect.py` (16/16).
>
> Consequences, all now in the shipped artifact:
>
> * `G3` reads `{FAIL: 2, PASS: 3, UNKNOWN: 29}`. Twenty-nine rows are honestly
>   **UNKNOWN** — they need a computed reduced polygon, which exists in-repo only
>   at `(8,28)` and `(5,20)`.
> * Two `(5,20)` rows flip `FAIL → PASS`: **`F_2(2,3)/75`** and **`F_3(3,2)/75`**
>   (both `a = 2`, cap `2 ≥ 2`). The slice cascade **is** available there and we
>   had written it off.
> * `(75,125)` is **unaffected** — at `a = 3` the gate fails on the correct object
>   too. The flagship open case never moved.
> * The clustering LOOSENS: 6 → **8** signatures, top two 82% → **74%**. Part of
>   the old concentration was an artifact — a gate that wrongly answers `FAIL`
>   merges rows not actually known to behave alike.
>
> Independently corroborated: `yplace_transfer.py` (57/57) recomputed the cascade
> at a class row's `y`-place, levels 2→12, reproducing `PROOF` §6.1's exact shape.
> **The cascade does transfer.**
>
> A different third mechanism *does* survive, found the same day: monomiality also
> forces `dg = deg C − ord C = 0`, which makes the F2-family closed form
> `f = −1/(a·dg)·y^ρ(y^dg+1)^e` **undefined** (`−1/(a·0)`). So it remains "one
> test, three mechanisms" — but the third is the family closed form, not the slice
> cascade.

Only **18 distinct `A_0` corners** carry the 34 rows, and only **four** of them
retract: `(6,15)`, `(8,28)`, `(9,24)`, `(12,33)`. Of those four, **only `(8,28)`
has `t = 4`**; the other three all have `t = 3`.

### Q2 — does `(8,28)/(3,4)/144` admit the toric identity?

**NO — decisively, and by both halves of the mechanism.** See §5.
`(2,3)` specifically is doing the work; the `t = 4` corner is not enough.
This **hardens**, rather than softens, the "individual endgames" reading.

---

## 1. What these gates do and do NOT mean — read this first

**No row in this atlas is eliminated as a counterexample by this atlas.** These
are gates on **mechanisms**, not on cases. `G1 = FAIL` means *our chart
dictionary is not usable at that corner*; `G3 = FAIL` means *the slice-cascade
argument is void there*, not that the case is safe. GGV5's 34 are **not** 34 open cases.

> ### CORRECTED 2026-07-27 -- the open frontier is **27**, not 32
>
> This paragraph used to read "34 open cases minus the one GGV5 itself discards
> and the one this campaign closed", i.e. 32. **Five more were already ruled out
> before GGV5 was written.** GGV5 tex:1794: *"In [M] there are listed four cases
> (which correspond to six cases in our terminology) of possible counterexamples
> with `max(deg P, deg Q) <= 100`. **They are discarded by hand.**"* -- and the
> red pairs in its table are exactly those `<= 100` rows. tex:1818 then accounts
> for all six red rows: **five are Moh's**, and the sixth (starred) is `F_22`,
> which GGV5 discards itself. Moh's own sixth case never entered the 34 -- the
> algorithm filtered it on `(2,1)` not in `PLLC`.
>
> | | |
> |---|---|
> | already settled in the literature | **6** -- `F_1(3,4)/64`, `F_2(2,3)/75`, `F_3(3,2)/75`, `F_9(2,3)/84`, `F_17(2,3)/99` (Moh), `F_22(2,3)/96` (GGV5) |
> | closed by this campaign | 1 -- `(8,28)/(3,2)/108` |
> | **open** | **27** |
>
> Corroborated independently: `max_deg <= 100` partitions the 34 as exactly 6 / 28,
> matching the red count without reference to the colouring.
>
> **Evidence boundary.** The red partition is EXACT-CHECKED (`moh_discards.py`,
> 21/21). That Moh's five are *ruled out* is **CITATION-LEVEL** -- [M] has not
> been read here, and this rests on GGV5's characterisation of it, exactly as
> `prop43_audit.py` discharges the GGHV22 Prop 4.3 citation without re-deriving
> the mathematics.
>
> **Direction.** The error made the atlas claim *more* is open than is, so no
> case-level assertion was wrong -- but it mis-priced work. In particular the two
> rows whose `G3` flipped `FAIL -> PASS` in v0.4.1, `F_2(2,3)/75` and
> `F_3(3,2)/75`, are **both** Moh discards. They are not new open ground. What
> they are is the project's **first external control in the monomial regime** --
> two settled cases where the slice cascade is now known available, at a corner
> unlike `(8,28)`, where every mechanism here was calibrated.

The value of the atlas is the *converse* reading: it says which of the 34 could
possibly be attacked by which of our five tools, and it says that the answer is
highly structured.

---

## 2. The five gates

| | gate | test | status | evaluated |
|---|---|---|---|---|
| **G1** | chart exponent & dictionary | `t = ceil(b_0/a_0)`, retraction `b_0 = t(a_0−1)`, `l_final = t` | **INFERRED** (`t` rule not published) | 34/34 |
| **G2** | toric admissibility | `(t+1) \| (4t+9)` ⟺ `t = 4` | **PROVED** (`toric_general.py` B4) | 34/34 |
| **G3** | slice cascade | `gcd(m,n)=1`, `lam >= m`, `N_Q >= D_P+D_Q` | **PROVED** (`contact_lemma.py`) | 29/34 |
| **G4** | Belyi passport | `u·kappa = m_f+n_f−1` over the full branch sweep | **PROVED** gate, **CHECKED** engine | 34/34 |
| **G5** | divisor syzygy | `q_window \| w(e)` | **PROVED** criterion | 3/34 |

**PROVED / CHECKED / INFERRED breakdown.** G2, G3, G4's gate and G5's criterion
are each proved in-repo with passing test suites. G4's *reduction engine* is
CHECKED — `passport_75_125.Reduction` reproduces all five published GGHV22
reductions from `(a_0,b_0)` alone. **G1's `t = ceil(b_0/a_0)` is INFERRED and is
the single largest load-bearing assumption in the atlas**: 28 of the 34 kills
route through it. It is not citable as published (see `chart_exponent`'s
docstring and `CORNER_RESOLVENT.md` §5.1). See §4 for what this atlas adds to
its evidence base.

---

## 3. Evaluability — how much is actually decided

| gate | PASS | FAIL | UNKNOWN | the missing input |
|---|---|---|---|---|
| G1 | 6 | 28 | 0 | — |
| G2 | 11 | 23 | 0 | — |
| G3 | 1 | 28 | **5** | `ord_y(Φ)`, `deg_y(Φ)`, `M` at that corner |
| G4 | 1 | 33 | 0 | — |
| G5 | 3 | 0 | **31** | `w(e)` — the split-weight enumeration at that corner |

**G1, G2 and G4 are fully decided on all 34 rows.** G3 is decided on 29: the 28
monomial corners fail it with *no* case-specific input (the `lam ≡ 0` argument),
and `(8,28)/(3,2)` passes from its known Φ. The 5 UNKNOWNs are exactly the
retracting rows other than (72,108) — the ones where `C` is genuinely not a
monomial and `Φ` has simply never been derived. **That is the concrete, named
work item this atlas produces:** derive `Φ` at `(6,15)`, `(9,24)`, `(12,33)` and
at `(8,28)/(7\4,3)`, and G3 closes on the whole published population.

G5 is the weakest: `q_window` itself is now computed everywhere (see §6), but the
criterion `q_window | w(e)` is `w(e)`-dependent, so only `q_window = 1` yields an
unconditional verdict. Three rows have it.

---

## 4. G1, and what the atlas adds to the `t` rule

`chart_exponent(a_0,b_0) = ceil(b_0/a_0)` was validated on five published GGHV22
reductions and pinned at `(5,20)` by GGV3. **This atlas corroborates it on four
more corners, by a mechanism the original validation did not use:**

> At **every** corner of the 34 where the retraction shape holds — `(6,15)`,
> `(8,28)`, `(9,24)`, `(12,33)`, six rows — GGV5's independently printed
> `l_final` equals `ceil(b_0/a_0)` exactly. The converse fails (4 non-retracting
> rows also happen to agree), so this is a one-sided but non-trivial test that
> the rule passed six times without being fitted to.

Two further controls in the checker:

- **`G1-CTRL`**: `polygon_reduction.chart_exponent` (`ceil(b_0/a_0)`) and
  `passport_75_125`'s independent rule (r1) (`floor((b_0−1)/a_0)+1`) — two files,
  two derivations — agree on all 34.
- **`G1-MUT`**, the historic trap: the *existential* reading "∃ `l` with
  `b_0 = l(a_0−1)`" admits **13** rows; the guarded reading admits **6**. The
  seven extra rows are exactly the `l = 5` family — `F_1`, `F_2(2,3)`,
  `F_2(3,5)`, `F_3`, and both `(7,42)` rows. This is the bug that put `l = 5`
  into the repo, reproduced as a live discriminator.

**28 rows carry a `SUSPECT` dictionary flag in the JSON** (`dictionary_trust`),
including `(7,21)/F_9` — the published landmine — and every `(5,20)` row.

---

## 5. Deliverable 3: does the toric mechanism transfer?

### 5a. `(8,28)` is the unique corner with both preconditions

Verified and extended (`corner_atlas.py` checks D3a–D3d):

- **Among the SPORADIC rows**, `(8,28)` is the only corner passing both the
  retraction test and `t = 4`. **Your claim is CONFIRMED.**
- **Extended to the FAMILY rows: still only `(8,28)`.** No family corner passes
  both; the two family corners that retract, `(6,15)` and `(9,24)`, both have
  `t = 3`.
- **The conjunction is doing real work.** Retraction alone admits four corners;
  `t = 4` alone admits five corners / 11 rows. Neither alone isolates `(8,28)`,
  and `(12,33)` is a genuine near-miss (retracts, but `t = 3`).

Two of the 34 rows sit on that corner: `(8,28)/(11\4,7)/(3,2)/108` = (72,108),
and `(8,28)/(7\4,3)/(3,4)/144`.

### 5b. The verdict on `(8,28)/(3,4)/144`: **NO**

`corner_atlas_toric_34.py`, exact rational linear algebra, both halves tested,
each negative paired with a positive control through the identical code path:

| test | `(3,4)` at `t=4` | control `(2,3)` at `t=4` |
|---|---|---|
| weight-admissible exponents `k` | `{5,6,7,8}` — **5 is among them** | `{5}` |
| **T1** product-of-minors `≡ e^k` | **no**, at every admissible `k` (1/36/67/15 pairs tested) | **HIT** at `k=5` |
| **T2** cofactor-is-itself-a-minor (`win`) | **no**, `k = 3..7` (3/11/41/102/234 cofactor columns) | **IN** at `k=5`, minimal |
| T3 looser (`all`, state cofactors) | IN at `k=5` | IN at `k=4` |

**T2 is the condition you singled out** — the one that turns the relation into a
product and yields the divisor law and the `Π^4` contact order — and it fails.
The weight arithmetic does permit `k = 5` there (`toric_general` B4 is a fact
about `t` alone, as proved), so the prediction got a fair test and lost.

**Which answer is it: `(2,3)` is doing the work.** Same corner, same retraction,
same `t = 4`, same admissible exponent — and neither half of the mechanism
survives. `t = 4` and the retraction shape are **necessary but not sufficient**.

**Independent corroboration by a different mechanism.** The Belyi sweep says the
same thing on the same pair of rows: `(8,28)/(3,2)` produces 4 admissible faces
with valid ramification data (all on the `en-split k=1` branch, `s ∈ {2,3}`,
`split_e ∈ {None,4}`); `(8,28)/(3,4)` produces **zero**, over 56 legal branches.
Two unrelated mechanisms, one corner, and both single out `(m,n) = (2,3)`.

**Logged, not asserted.** The looser invariant T3 reads `4 / 5 / 6` at
`(2,3) / (3,4) / (3,5)`. Three points fit `max(m,n)+1`. Per `toric_general` D2's
own warning about one-point fits, this is **not** claimed as a law; it is logged
for the next corner to falsify.

---

## 6. Findings that were not asked for

**(i) `MINIMAL_CORE.md`'s "negative I could not overcome" is partly superseded.**
That document states a sweep of `q_window` across the 34 corners "is **not
possible with the data in this repo**", because `upstream_facts.json` carries
Newton polygons for `(8,28)` only. But `q_window_theorem.py` proves the **closed
form** `q_window = M/gcd(M,H)` with `M = t(a+b)−(kappa+1)`, `H = q(a+b)−1`, which
needs only `(t, kappa, ord C, m, n)` — all of which the guard supplies at every
corner. **`q_window` is now computed for all 34 rows.** Control `G5-CTRL`: fed
the *guarded* inputs it reproduces all three values the repo already knew —
`17` at (50,75), `29` at (75,125) (`= 12·3−7`, matching
`window_functions.q_window(3)`), and `1` at (72,108) — from a different
derivation path. What remains genuinely untestable is `w(e)`; `MINIMAL_CORE.md`'s
*other* statement, that `q_window | w(e)` is untestable off (8,28), stands.

**(ii) `q_window = 1` is not unique to (72,108), and not a sporadic phenomenon.**
Three rows have it: `(8,28)/(3,2)/108`, `(12,33)/(2,3)/135` — and **`F_17(2,3)/99`,
a FAMILY case.** This is direct data for `MINIMAL_CORE.md`'s claim that the
divisor syzygy is "not a sporadic-corner phenomenon" and is "unrelated to family
membership", which until now had one data point. At these three rows the syzygy's
carry obstruction vanishes on **every** split, unconditionally.

**(iii) The guard moves `q_window` a lot.** The JSON carries both
`q_window` (guarded) and `q_window_unguarded` (GGV5's final-corner dictionary
taken at face value). They differ on 28 of 34 rows — e.g. `F_9(2,3)/84` reads
`13` guarded vs `29` unguarded, and **eight** rows read `1` unguarded while their
guarded value is `13`, `17`, `21` or `25`. **Any landed result keyed to an
unguarded `q_window` off (72,108) should be re-checked**; the repo's landed `F_9`
signature is one of them.

**(iv) `(8,32)` and `(8,40)` are blow-ups of our corner.** Both length-2 sporadic
chains pass through the intermediate corner `(8,28)` — the (72,108) corner —
before reaching the same final corner `(11\4,7)`. Their `A_0` does not retract,
so they land in classes 2 and 1 respectively; but the repo applies the retraction
test at `A_0` only (`case_compiler` handles chain length 1). **Whether the guard
should be applied at `A_0` or at the last integer corner of the chain is an open
modelling question that changes the verdict on exactly these two rows.** Flagged,
not decided.

**(v) `contact_lemma.py:1115` — confirmed wrong, as you suspected.** The row
`(m, n, ell, lam) = (3, 5, 5, 3)` hardcodes the **pre-repair** chart parameter:
`ell` is the chart parameter and equals `t` (`MINIMAL_CORE.md` §2.0), `ell = 5`
is the refuted `l = 5`, and with the repaired `t = 4` and `lam = 0` the row
should read `(3, 5, 4, 0)`. **Reported, not fixed** — it did not block this
atlas, which derives `ell` from the guarded `t`. Note the repair is not cosmetic:
`lam = 0` means the `capok` assertion at that row is testing a cap that does not
exist.

---

## 7. Transcription integrity

The family `A_0` column is **not printed** in the counterexample table; it had to
be joined in from the chain-data table (tex:1678-1694, 1709-1715). Two published
identities are re-derived from the transcription and must hold:

- **X1** `max{deg P, deg Q} = v11(A_0)·max(m,n)` — **holds on all 34 rows.** This
  is what validates the join; a mis-joined `A_0` breaks it immediately.
  Mutation control: shifting `v11(A_0)` by 1 breaks it on every row.
- **X2** `(m+n)·q·k − n·(q·l − p) = k` (GGV5's own admissibility identity) —
  **holds on all 13 family rows**, which also confirms that the final corner of
  the two length-2 families `F_22`, `F_24` is `A_2`, not `A_1`.

Census controls: 13 FAMILY + 21 SPORADIC = 34; sporadic chain lengths split
9 / 11 / 1 exactly as the paper's prose states (tex:1825, 1839, 1862); every
family row resolves to a unique `j >= 0` in its family law.

---

## 8. Non-discriminating tests, declared

Per the project's standing rule that a trivially-true check is worse than no
check, the checker asserts these are non-discriminating and says so:

| sub-test | verdict on all 34 | why it is still kept |
|---|---|---|
| `gcd(m,n) = 1` | PASS ×34 | GGV5's own definition requires it. PROVED and load-bearing *in general* (`contact_lemma` F7 gives counterexamples at `(2,4),(3,6),(4,6),(2,6)`) — but it filters nothing here. |
| `N_Q >= D_P + D_Q` | PASS ×34 | min slack 6; `contact_lemma` F3 already notes it "only binds at small `ell`". |
| Belyi **proportional** class | FAIL ×34 | needs `kappa >= 2(m+n)−1`; min margin 5. Proved empty, so only **en-split** faces can carry a passport — which is what the sweep then tests. |

The remaining tests all discriminate; distributions are in §3 and in the JSON.

---

## 9. Where this leaves the frontier

1. **The next case is not "(75,125)". It is a class of nine.** Class 2 —
   monomial corner with `t = 4` — contains `F_2(3,5)/125` together with eight
   others: both `(9,36)/(17\9,4)` rows, `(9,36)/(11\3,8)`, `(8,32)/(3,2)/120`,
   both `(10,40)` rows, `F_2(2,3)/75` and `F_3(3,2)/75`. Everything proved about
   (75,125)'s monomial corner — `lam = 0`, `C` a monomial, dictionary refused —
   is a **class** fact, not a case fact. That is nine cases for the price of one.
2. **Class 1 (19 rows) is the same story with `t != 4`**, and additionally has no
   toric weight solution at all. 28 of 34 rows die on one integer test.
3. **The real remaining structure is six rows on four retracting corners.** Four
   of them (`F_7`, `F_8`, `F_17`, `(12,33)`) have `t = 3`, so no toric identity is
   possible there by G2 — but their `lam` is genuinely unknown and their cascade
   is genuinely open. **Deriving `Φ` at `(6,15)`, `(9,24)`, `(12,33)` is the
   single highest-leverage missing computation in the atlas**: it closes G3 on
   the entire published population.
4. **The `(8,28)` corner is exhausted.** Both its rows are now understood: one
   closed, one shown to carry neither the toric nor the Belyi mechanism.
5. **Caveat with teeth.** Items 1–2 rest on the INFERRED `t = ceil(b_0/a_0)` rule.
   If that rule is wrong the 28-row kill collapses. It now has ten independent
   corroborations (five published reductions, GGV3 at `(5,20)`, and the four
   retracting corners here) and zero contradictions, but it is not published.
   **Pinning it is worth more than any individual case.**
