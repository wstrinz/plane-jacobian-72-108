# GGV3_CITATION_CHECK — closing writeup gap G7 (the "GGV3 §1" strip-descent citation)

Date: 2026-07-22. Task: verify the citation `\cite{GGV3}*{Section 1}` that GGHV22
(arXiv:2204.14178) makes for the α-strip / lower-order-perturbation argument, which
`T6_PREMISES.md` fills from GGV1 (arXiv:1401.1784) Props 1.13 and 2.1 instead of from
GGV3 itself (GGV3 was not in `paper_src/` until this session).

**Verdict up front: SUBSTITUTION VERIFIED** (with one honest reconstruction nuance,
§4). GGV3 §1 was fetched, and it *does* contain the cited strip-descent argument.
GGHV22's "GGV3 §1" citation is narrow — it delegates only the *iteration* of a
single-step cancellation that GGHV22 has itself already set up two lines earlier with
GGV1 Props 1.13/2.1. T6_PREMISES's reconstruction reproduces exactly that iteration and
loses no content.

---

## 1. Identification of GGV3 — **CONFIRMED, high confidence**

| Field | Value |
|---|---|
| Repo alias | **GGV3** |
| arXiv id | **1406.0886** |
| Title | *A system of polynomial equations related to the Jacobian Conjecture* |
| Authors | Jorge A. Guccione, Juan José Guccione, Christian Valqui |
| Abstract | "We prove that the Jacobian conjecture is false if and only if there exists a solution to a certain system of polynomial equations. We analyse the solution set of this system. In particular we prove that it is zero dimensional." |

Evidence for the identification (three independent, mutually consistent):

1. `NEXT_CASES.md` §4 paper-chain table (L200): `[GGV3] arXiv:1406.0886 — "A system of
   polynomial equations related to the JC" — the polynomial-system machine; §5 discarded
   both 75 cases`.
2. GGHV22's own bibliography, `paper_src/2204.14178.tex` L2117–2123:
   `\bib{GGV3}{...} title={A system of polynomial equations related to the Jacobian
   Conjecture}, eprint={arXiv:1406.0886}`.
3. arXiv abstract page for 1406.0886 matches title/authors verbatim.

**Section-number cross-check (resolves the version-drift worry).** The fetched source is
GGV3 **v3 (posted 2024-04-07)**, which post-dates GGHV22 (2022). To confirm the section
numbering GGHV22 cited still matches, note GGHV22 cites GGV3 *twice*: `Section 1` (the
strip descent, L1517/L1923) **and** `section 5` (the 75-case kill, L298/L310/L311). In
the fetched v3 source, §1 is "The Jacobian Conjecture as a system of equations" (L249) and
§5 is "A modified system and an example" analysing exactly `(n,m)=(50,75)` — Moh's 75 case
(L235–247, L1708). Both citations land on the right section in the fetched version, so the
§1 numbering is stable across the versions GGHV22 saw. [judgment] Version drift is
therefore not a live risk for this citation.

---

## 2. What was fetched and saved

- Fetched `https://arxiv.org/e-print/1406.0886` (gzipped single-file LaTeX,
  original name `Polynomial_Equations_related_to_the_Jacobian_Conjecture_19_03_2024.tex`).
- Decompressed and saved as **`paper_src/1406.0886_GGV3.tex`** (2271 lines,
  `\end{document}` present at L2271 — complete), matching the existing convention
  (`1401.1784_GGV1.tex`, `1708.07936_GGV5.tex`, `2204.14178.tex`).

Section map of the fetched source:

| § | Line | Heading |
|---|---|---|
| Intro | 199 | Introduction |
| **1** | **249** | **The Jacobian Conjecture as a system of equations** ← cited |
| 2 | 762 | Properties of solutions of `S(n,m,(λ_i),F_{1-n})` |
| 3 | 1046 | The homogeneous system `S(n,m,F_{1-n})` |
| 4 | 1223 | Presentations of the solutions of `S(n,m,Y^{m+n-1})` |
| 5 | 1708 | A modified system and an example (the `(50,75)` kill) |

---

## 3. Locating the cited strip-descent inside GGV3 §1

GGHV22's citation, `2204.14178.tex` L1517 (identical copy at L1923 for the (9,24) base
`x³y` case):

> "…$\ell_{1,0}(Q-C^3)=\alpha_k (x^3 C_3)^k$, for some $k$ with $-2<k<3$. **Using the
> arguments of `\cite{GGV3}*{Section 1}`**, we find $\alpha_2,\alpha_1,\alpha_0,\alpha_{-1}\in K$
> such that $[\ell_{1,0}(P),\ \ell_{1,0}(Q-C^3-\alpha_2 C^2-\alpha_1 C-\alpha_0-\alpha_{-1}C^{-1})]\ne 0$."

The referenced argument is the **inductive "claim" inside the proof of GGV3 Theorem
`principal` (Theorem 1.8)**, `1406.0886_GGV3.tex` **L505–L535**. Verbatim skeleton:

- Start from `Q`, subtract the leading power `C^m` (`ℓ` of `Q` is a pure power of `C`).
- **Inductive step (L510–L535).** Given `λ_1,…,λ_i` with
  `deg(Q − C^m − λ_1 C^{m-1} − ⋯ − λ_i C^{m-i}) ≤ m−i−1`, set
  `Q̃ := Q − C^m − ⋯ − λ_i C^{m-i}`. If `deg(P)+deg(Q̃)−2 ≤ 0` stop (take remaining
  `λ = 0`); otherwise `deg(Q̃) > 2−n` and **Proposition `vdE 10.2.11` (L426)** yields
  `j`, `λ_j∈K^×` with `deg(Q̃ − λ_j C^{m-j}) < deg(Q̃)`. Subtract it; the exponent index
  advances `i+1 ≤ j ≤ m+n−3`. Iterate.
- **Termination (L535–L562).** The descent halts with a remainder
  `F = Q − C^m − Σ λ_i C^{m-i}` of `deg ≤ 2−n`, and then `F_+ = x^{1-n}(μ_0+μ_1 y)` with
  `μ_1 ≠ 0` by **Lemma `auxiliar` (L403)** — i.e. `[P,F]∈K^×` forces the remainder's
  leading form to be the *non-power* corner `x^{1-n}y`, so the descent cannot continue.

So GGV3 §1's strip descent = "peel leading powers `C^{m-i}` off `Q` one integer step at a
time until the leading bracket `[P,F]` is a nonzero scalar." That is *exactly* the
normal-form reduction `Q = Σ λ_i C^{m-i} + F` GGHV22 delegates.

**Its engine, however, is not GGV1 Props 1.13/2.1.** GGV3 §1 runs the descent in the
**(1,1)-grading** using van den Essen's **Lemma 10.2.11** (packaged as GGV3 Prop
`vdE 10.2.11`, L426) as the single-step degree-lowering lemma, plus GGV3 Lemma `auxiliar`
(L403) for termination. GGV3 §1 does **not** cite or restate GGV1 Props 1.13/2.1 (GGV3
uses GGV1 only in §5, for the reduction-of-degree technique — L236, L1716–1723). This is
the one substantive difference from T6_PREMISES's reconstruction; see §4–5.

---

## 4. Proposition-by-proposition comparison

Roles: **[E]** = single-step degree/valuation-lowering engine; **[T]** = termination
criterion; **[N]** = normalization of the final leading form.

| Role | GGV3 §1 (as cited) | GGV1 statement T6 uses instead | Where T6_PREMISES uses it | Same content? |
|---|---|---|---|---|
| Bracket-degree inequality (setup) | `wdeg([P,Q]) ≤ wdeg(P)+wdeg(Q)−|w|`, GGV3 L357 (stated, no equality clause) | **GGV1 Prop 1.13** (`pr v de un conmutador`, GGV1 L480–491): the *same* inequality **plus** the equality-iff-`[ℓP,ℓQ]≠0` clause | T6 P1a §1.2(a) L66–75; P2 §2.2(a) L157–169 | Yes — GGV1 Prop 1.13 is the sharper (with-equality) form; GGV3 gets the equality info implicitly through vdE 10.2.11 |
| **[E]** one strip step ("leading form is a pure power of `C`, so subtract `λ C^{m-i}`") | **GGV3 Prop `vdE 10.2.11`** (van den Essen Lemma 10.2.11), L426–435, (1,1)-grading | **GGV1 Prop 1.13 + Prop 2.1** (`P y Q alineados`, GGV1 L514–563), (1,0)-grading: `[ℓP,ℓG]=0` ⟹ `ℓ(G)=α_k R^k` for `R` primitive | T6 §2.2(a) L159–166 (single step); the iteration is the whole strip | **Equivalent, different toolset.** Both prove "the current leading form is a scalar × a power of `C`, cancel it." GGV3 via vdE degree drop; T6 via 1.13-forces-commuting-leading-forms + 2.1-forces-common-power |
| **[E]** iteration / bookkeeping | GGV3 Thm `principal` proof, inductive claim L505–535 | Repeated application of the single step above; range `−2<k<3` ⟹ `k∈{−1,0,1,2}` (four terms) | T6 §2.2(a) L164–169 | Yes — GGHV22's four `α_k` are the `t=4` copy of GGV3's `λ_i` list; the finite range is verified (`t6_premises_verify.py` P2c) |
| **[T]** termination | GGV3 Lemma `auxiliar` L403: `[P,F]∈K^×` ⟹ `F_+ = x^{1-n}(μ_0+μ_1y)`, `μ_1≠0` (leading form not a power ⟹ descent stops) | GGV1 Prop 1.13 **equality case**: `[ℓP,ℓF]≠0` ⟺ `v([F,P]) = v(F)+v(P)−1` | T6 §2.2(b) L175–182; halts at `v(F)=−5` because `−5 ≢ 0 (mod 4)` so `ℓ(F)` is not a power | Yes — both say "descent stops exactly when the leading bracket is a nonzero scalar"; T6's `v(F)=−5` ↔ GGV3's `deg(F)=2−n` (`t=4`: `2−? ` recomputed in the (1,0)-grading, `v(F)+v(P)−1 = 2` ⟹ `v(F)=−5`) |
| **[N]** final leading-form normalization | GGV3 L554–562: affine `y ↦ (y−μ_0)/μ_1` makes `F_+ = x^{1-n}y` | Change of variables `y=(a₀/a₁)ỹ` (Premise 1 §1.2(e)) + scalar gauge-fix (§1.3) | T6 §1.2(e) L104–111; §1.3 L113–126 | Yes — same normalization move (absorb the corner's two scalars), applied to `C₄=y⁷(y+1)` |

### (a) Does GGV3 §1 contain the same propositions (renumbered/generalized)?

**Partly.** GGV3 §1 contains the *same argument* (the iterated top-down strip) and the
*same output* (normal form `Q = Σλ_i C^{m-i} + F`, `F` at the minimal valuation with
`[P,F]∈K^×`). But it does **not** contain GGV1 Props 1.13/2.1 as such — its single-step
engine is van den Essen's Lemma 10.2.11 in the (1,1)-grading, and its termination lemma is
GGV3 Lemma `auxiliar`. So the propositions are **not** the literal 1.13/2.1; they are a
*different but equivalent* single-step-descent toolkit. [judgment]

### (b) Does the T6 substitution prove everything the citation was carrying?

**Yes.** The citation carries a strictly narrow obligation: *given* that each successive
leading form `ℓ_{1,0}(Q−C^3−Σα_k C^k)` is a scalar multiple of `R^k` (which **GGHV22
itself** establishes one line earlier, L1508–1516, via GGV1 Props 1.13/2.1 — not via
GGV3), produce the finite coefficient list `α_2,α_1,α_0,α_{-1}` that peels those forms off
until `[ℓP,ℓF]≠0`. That is pure iteration of the single step GGHV22 already has in hand.
T6_PREMISES (§2.2(a), L159–169) does exactly this iteration with the same Props 1.13/2.1,
plus the finite check that the admissible powers are `k∈{−1,0,1,2}` and that the descent
halts at `v(F)=−5` because `−5` is not a multiple of `v_{1,0}(R)=4`
(`t6_premises_verify.py` P2c). Nothing in GGV3 §1's argument is used by GGHV22 beyond this
iteration, so the substitution is complete.

---

## 5. Verdict: **SUBSTITUTION VERIFIED**

The `\cite{GGV3}*{Section 1}` reference resolves to the strip-descent inside GGV3 Theorem
`principal` (`1406.0886_GGV3.tex` L491–588, iteration at L505–535). That argument is
present, complete, and does exactly what T6_PREMISES Premise 2 reconstructs: strip leading
powers of `C` from `Q` one step at a time until the remainder's leading bracket is a
nonzero scalar. T6's reconstruction using **GGV1 Props 1.13/2.1** proves everything the
citation carries, because:

1. The obligation delegated to GGV3 §1 is only the *iteration* of a single cancellation
   step; GGHV22 sets up that single step itself with GGV1 Props 1.13/2.1 (L1508–1516),
   independent of GGV3.
2. The iteration is finite and explicitly bounded (`k∈{−1,0,1,2}`), and its termination
   (`v(F)=−5`, `−5 ∤ 4`) is verified arithmetic (`t6_premises_verify.py` P2c).
3. The `t=3→t=4` transport changes only valuations (`v(P): 6→8`, `v(F): −4→−5`); the
   structural engine is grading-general (GGV1 1.13/2.1 hold for any `(ρ,σ)∈𝔙`).

**Honest reconstruction nuance [judgment].** This is a *verified substitution*, not a
line-for-line match to GGV3's own proof. GGV3 §1 runs the descent in the **(1,1)-grading**
with **van den Essen Lemma 10.2.11** as the engine and GGV3 Lemma `auxiliar` for
termination; T6 runs it in the **(1,0)-grading** with **GGV1 Props 1.13/2.1**. The two are
mathematically equivalent (both peel one leading power per step and stop when the leading
bracket is a nonzero scalar), and the (1,0)/1.13-2.1 formulation is the one GGHV22 *itself*
uses for the single step — so the reconstruction is faithful to GGHV22's own machinery, not
a foreign import. No mathematical gap; the difference is purely which equivalent single-step
lemma names the descent.

---

## 6. What the writeup must say about this citation

1. **State the substitution explicitly** (do not silently rely on "GGV3 §1"). Recommended
   wording: *"GGHV22 cites `\cite{GGV3}*{Section 1}` for the construction of the strip
   coefficients `α_k`. GGV3 §1 (Theorem 1.8, the descent in its proof, arXiv:1406.0886
   L505–535) carries out this strip in the (1,1)-grading using van den Essen's Lemma
   10.2.11. We reconstruct the identical iteration in the (1,0)-grading from GGV1 Props
   1.13 and 2.1 — the same single-step tool GGHV22 uses one line earlier — so the
   substitution loses no content."*
2. **The source is now on hand**: `paper_src/1406.0886_GGV3.tex` (v3, 2024; complete).
   Section-1 numbering matches what GGHV22 cited (verified via the independent `GGV3 §5`
   ↔ 75-case cross-check).
3. **No gap remains.** The premise's verdict in `T6_PREMISES.md` ("READY-WITH-CITATION",
   §2.3) stands; the parenthetical "GGV3 is not in `paper_src/`" (T6_PREMISES L172–173) is
   now stale and should be updated to point at `paper_src/1406.0886_GGV3.tex` with the note
   that GGV3 §1's engine is van den Essen Lemma 10.2.11, reconstructed here from GGV1
   1.13/2.1.
4. **One caveat to keep [judgment]:** the fetched source is v3 (2024), later than the
   version GGHV22 cited. Section-number alignment is confirmed by the dual `§1`/`§5`
   cross-check, but the writeup should cite GGV3 by *theorem role* ("the JC-as-system strip
   descent, Thm `principal`") rather than a bare "Section 1", so the citation is robust to
   any future renumbering.
