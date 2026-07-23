# LINALG_CERT_REPORT — NulLA / Macaulay-matrix certificate scoping (2026-07-21)

Scoping question: is a linear-algebra Nullstellensatz certificate search feasible for the
a=0 (e_0 ≠ 0) stratum of the f31 subcase-(2) g-cascade system (T5_NP.md), and if a viable
instance exists, run it.

Setup: 32 unknowns — a_0..a_4 (wt 2), b_0..b_6 (wt 3), c_0..c_8 (wt 4), e_0..e_10 (wt 5);
window polys d̃2, d̃1, d̃0, e = d̃₋₁. Block 1 = the 30 t-coefficients t^0..t^29 of h_0(d̃),
each weighted-homogeneous of weight 20. Certificate ansatz: e_0^N = Σ λ_i g_i with λ_i
weighted-homogeneous of weight 5N − 20, over F_32003.

**VERDICT: NO-GO — and stronger than a sizing no-go.** The homogeneous block-1-only
certificate provably does not exist for ANY N (structural bigrading argument below,
confirmed by exact F_32003 rank computations for N = 4, 5, 6, 8). The mixed/affine
versions (blocks 1+2, or the 18-var pivot system) blow past 10^7 columns at the very
first useful degree. Details and tables follow.

---

## Task 1 — weighted monomial counts (32 vars, weights 2⁵ 3⁷ 4⁹ 5¹¹)

Generating function 1/((1−x²)⁵(1−x³)⁷(1−x⁴)⁹(1−x⁵)¹¹), exact DP.

| W | count | W | count | W | count |
|---|-------|---|-------|---|-------|
| 20 | 508,161 | 34 | 230,679,377 | 48 | 31,935,679,118 |
| 21 | 831,262 | 35 | 338,851,662 | 49 | 43,967,625,638 |
| 22 | 1,346,947 | 36 | 494,959,121 | 50 | 60,308,896,932 |
| 23 | 2,160,186 | 37 | 719,044,306 | 51 | 82,425,507,992 |
| 24 | 3,434,144 | 38 | 1,039,101,941 | 52 | 112,257,175,717 |
| 25 | 5,409,854 | 39 | 1,493,968,353 | 53 | 152,361,802,940 |
| 26 | 8,452,686 | 40 | 2,137,386,167 | 54 | 206,103,130,075 |
| 27 | 13,099,376 | 41 | 3,043,284,580 | 55 | 277,890,754,607 |
| 28 | 20,147,611 | 42 | 4,313,090,705 | 56 | 373,489,711,083 |
| 29 | 30,758,012 | 43 | 6,085,219,898 | 57 | 500,414,937,798 |
| 30 | 46,629,241 | 44 | 8,548,034,782 | 58 | 668,436,805,847 |
| 31 | 70,207,253 | 45 | 11,956,656,514 | 59 | 890,222,736,350 |
| 32 | 105,025,533 | 46 | 16,655,617,104 | 60 | 1,182,154,725,475 |
| 33 | 156,120,386 | 47 | 23,108,228,304 | | |

Low weights used below: cnt(5) = 46, cnt(10) = 1,789, cnt(15) = 36,256, cnt(18) = 183,569.
The weighted grading does grow far slower than the ordinary degree grading, but weight 20
is already a 508k-dimensional space, and each certificate step multiplies by ~6–8×.

## Block-1 equations, computed exactly (flint fmpz_mpoly, 4 s)

h_0(d̃(t)) has **508,161 terms total — it is DENSE in the weight-20 space** (support =
every weight-20 monomial, exactly cnt(20)). Terms per t-slice (palindromic in k ↔ 40−k):

| t^k | terms | t^k | terms | t^k | terms |
|-----|-------|-----|-------|-----|-------|
| 0 | 28 | 10 | 8,747 | 20 | 33,345 |
| 1 | 71 | 11 | 11,371 | 21 | 32,863 |
| 2 | 182 | 12 | 14,415 | 22 | 31,695 |
| 3 | 372 | 13 | 17,589 | 23 | 29,672 |
| 4 | 721 | 14 | 20,949 | 24 | 27,184 |
| 5 | 1,232 | 15 | 24,139 | 25 | 24,139 |
| 6 | 2,036 | 16 | 27,184 | 26 | 20,949 |
| 7 | 3,109 | 17 | 29,672 | 27 | 17,589 |
| 8 | 4,598 | 18 | 31,695 | 28 | 14,415 |
| 9 | 6,435 | 19 | 32,863 | 29 | 11,371 |

Block 1 (t^0..t^29) total nnz 480,630; all 41 slices 508,161. Average ~16k terms/equation.
h_1(d̃) similarly has 183,569 terms = cnt(18) (dense in weight-18 space).

## Task 2 — Macaulay matrix sizes for e_0^N = Σ λ_i g_i (30 block-1 eqs)

Rows = 30 × cnt(5N−20) multiplier monomials; cols = cnt(5N); full-matrix nnz =
cnt(5N−20) × 480,630 (monomial-shift preserves term counts).

| N | wt 5N | mult wt | #mult | rows | cols | nnz (full) |
|---|------|--------|-------|------|------|-----------|
| 4 | 20 | 0 | 1 | 30 | 508,161 | 4.8×10^5 |
| 5 | 25 | 5 | 46 | 1,380 | 5,409,854 | 2.2×10^7 |
| 6 | 30 | 10 | 1,789 | 53,670 | 46,629,241 | 8.6×10^8 |
| 7 | 35 | 15 | 36,256 | 1,087,680 | 338,851,662 | 1.7×10^10 |
| 8 | 40 | 20 | 508,161 | 15,244,830 | 2,137,386,167 | 2.4×10^11 |
| 9 | 45 | 25 | 5,409,854 | 1.6×10^8 | 1.2×10^10 | 2.6×10^12 |
| 10 | 50 | 30 | 46,629,241 | 1.4×10^9 | 6.0×10^10 | 2.2×10^13 |
| 11 | 55 | 35 | 3.4×10^8 | 1.0×10^10 | 2.8×10^11 | 1.6×10^14 |
| 12 | 60 | 40 | 2.1×10^9 | 6.4×10^10 | 1.2×10^12 | 1.0×10^15 |

Under the ~10^7-nnz budget: N = 4 comfortably, N = 5 marginally (2.2×10^7). N ≥ 6 out.
**But the full matrix is never needed — see the bigrading collapse below.**

## Task 3 — block 2 is NOT weighted-homogeneous (claim verified by substitution)

Block-2 numerator = 6630·e(t)³·g_1(t) + ũ(t)·h_1(d̃(t)) with
ũ(t) = −2048t⁴ + 8704t³ − 14144t² + 10608t − 3315 (= 6630·u, y = t−1; ũ(0) = −3315 ≠ 0).
Computed exactly: **all 30 slices t^0..t^29 have weight set {18, 35}** — the e³g_1 part
contributes weight 35 (= 15 + 20) and the ũh_1 part weight 18, in every slice. Sample
term counts: slice t^0: 8,769; t^5: 93,205; t^13: 376,447; t^29: 69,566. Max total degree
in the unknowns: 13 (matches the "degree ≤ 13" bookkeeping). The u-constants also break
the index grading (e³g_1 slices have index 30+k, ũh_1 slices index k−4..k). So the
homogeneous certificate can indeed only use block 1 (or the 41-slice h_0 ≡ 0 extension).

## The bigrading collapse: block-1 certificates provably fail for ALL N

**Observation (verified exactly on all 41 slices):** the 32 variables carry a second
grading, idx(a_i) = idx(b_i) = idx(c_i) = idx(e_i) = i (the subscript), and slice t^k of
h_0(d̃) is precisely the idx-degree-k component. (Numerical fingerprint: the 41 slices
have pairwise disjoint supports, 508,161 = Σ slice sizes.) Every block-1 generator g_k
is BIhomogeneous: (weight, idx) = (20, k). The 41-slice extension is the same ideal
family, bidegrees (20, 0..40).

**Theorem.** For every N and every subset of the 41 slices containing g_0,
e_0^N ∈ ideal(slices) ⟺ h_0(a_0,b_0,c_0,e_0) divides e_0^N in F_p[a_0,b_0,c_0,e_0]
— which is false for all N (h_0(a_0,b_0,c_0,e_0) has 28 terms, is not a monomial).

*Proof.* The ideal is bihomogeneous; e_0^N has bidegree (5N, 0). Projecting a putative
identity e_0^N = Σ λ_k g_k onto idx-degree 0 kills every term with k ≥ 1 (multiplier idx
degrees are ≥ 0), leaving e_0^N = (λ_0)₀ · g_0 inside the idx-0 subring F_p[a_0,b_0,c_0,e_0],
where g_0 = h_0(a_0,b_0,c_0,e_0). ∎

**Geometric witness (why no depth N can ever work):** any constant window
d̃2 = a_0, d̃1 = b_0, d̃0 = c_0, e = e_0 with h_0(a_0,b_0,c_0,e_0) = 0 kills ALL 41 slices
(slices k ≥ 1 vanish when the positive-index variables are 0), and since
h_0(0,0,0,e_0) = −6561e_0⁴, generic points of this hypersurface have e_0 ≠ 0. So the
block-1 variety (even the full h_0(d̃) ≡ 0 variety) contains a 3-parameter family of
points with e_0 ≠ 0 — over the algebraic closure and over F_p alike. Two consequences:

1. "Empty stratum" is FALSE for block 1 alone (and for the 41-slice extension): the
   negative certificate answer is *not* "certificate too deep" — the stratum is genuinely
   non-empty. Only block 2's u-mixing can kill these constant-window solutions.
2. dim(block-1 variety) ≥ 3, i.e. strictly above the expected 32 − 30 = 2 — block 1 is
   not a complete intersection (contains the constant-window hypersurface {h_0 = 0} ⊂ A⁴).

## Solve outcomes (exact, F_32003, flint nmod_mat)

By bihomogeneity the Macaulay matrix is block-diagonal by (weight, idx) bidegree and the
target e_0^N sees only the tiny bidegree-(5N, 0) block (idx-0 multipliers × g_0). Ran
those exactly, plus one full-support sanity check:

| instance | rows × support-cols | rank(A) | rank(A\|e_0^N) | verdict |
|----------|--------------------|---------|----------------|---------|
| N=4, bideg(20,0): {1}·g_0 | 1 × 28 | 1 | 2 | NOT in span |
| N=5, bideg(25,0): {e_0, a_0b_0}·g_0 | 2 × 42 | 2 | 3 | NOT in span |
| N=6, bideg(30,0): 7 idx-0 wt-10 mults · g_0 | 7 × 70 | 7 | 8 | NOT in span |
| N=8, bideg(40,0): 28 idx-0 wt-20 mults · g_0 | 28 × 144 | 28 | 29 | NOT in span |
| N=4 sanity, all 30 slices directly | 30 × 480,630 | 30 | 31 | NOT in span |

All consistent with the theorem (and the 30 slices are linearly independent, as the
disjoint supports force). No deeper N needs to be run — the answer is NO for all N.

## Task 3b/4 — affine (dehomogenized, e_0 = 1) sizing: blocks 1+2, and the 18-var pivot

Full-column-space counts for 1 = Σ λ_i g_i, deg λ_i ≤ D, equations degree ≤ 13:

| D | 31 vars, monos deg ≤ D+13 | 18 vars (pivot system), deg ≤ D+13 |
|---|---------------------------|-------------------------------------|
| 0 | 5.19×10^10 | 2.06×10^8 |
| 1 | 1.67×10^11 | 4.71×10^8 |
| 2 | 5.12×10^11 | 1.04×10^9 |
| 3 | 1.50×10^12 | 2.20×10^9 |
| 4 | 4.24×10^12 | 4.54×10^9 |

Both exceed 10^7 columns **before the first multiplier degree**: the 31-var system at
D = 0 is already 5×10^10, the 18-var pivot system at D = 0 already 2×10^8. (The pivot
system was sized only — t5_pivot_eqs.txt is being generated by another agent.)

Support-restricted nnz is no rescue: block-2 slices carry up to 376k terms each, so even
D = 1 (rows = 60 × 32) runs to ~10^8 nnz from block 2 alone. The only sub-budget affine
instance is D = 0 (60 rows, ~5×10^6 nnz): "is 1 an F_p-linear combination of the 60 raw
equations" — a certificate that shallow has essentially zero chance and its (near-certain)
negative outcome proves nothing, so it was not run.

## Conclusion (go/no-go)

* **Homogeneous block-1-only search: worse than no-go — provably fruitless.** e_0^N is
  never in the block-1 (or 41-slice) ideal, for any N, by the bigrading collapse; the
  stratum it would certify empty is in fact non-empty (constant-window family). This is
  the definitive version of the "block 1 alone has dim ≥ 2, hence no certificate"
  expectation — now a theorem with witness, and it upgrades "expected dim 2" to
  "dim ≥ 3, not a complete intersection".
* **Mixed/affine blocks 1+2: no-go on size.** 5×10^10 columns at D = 0; every increment
  multiplies by ~3. The weighted grading cannot be exploited there because block 2
  breaks BOTH gradings (verified: every slice mixes weights 18 and 35).
* **18-var pivot system: no-go on size** (2×10^8 columns at D = 0).
* **No certificate-search instance below ~10^7 nnz proves anything meaningful.** The
  meaningful statement (emptiness of blocks 1+2 on e_0 ≠ 0) inherently couples the two
  inhomogeneous weights; brute Macaulay linear algebra pays for that coupling in full
  monomial-space dimensions. The structure that makes the problem tractable (the t-adic
  cascade, the bigrading) is exactly what Gröbner-style elimination (Singular, as planned
  in T5_NP.md) exploits adaptively and flat linear algebra cannot.

Positive by-products of this scoping worth keeping:
1. The **index bigrading** of the h_f(d̃) slices (slice k = idx-degree-k component) —
   likely useful far beyond this report (e.g. it decomposes any block-1-internal linear
   algebra into 30 small independent pieces, and constrains Gröbner computations).
2. h_0(d̃) is dense in the weight-20 space (support = all 508,161 monomials); h_1(d̃)
   likewise in weight 18.
3. Block-1's variety strictly exceeds expected dimension (constant-window family), so
   any emptiness proof MUST use block 2 — quantifying exactly why T5_NP.md's "block 1 +
   block 2" pairing is the minimal meaningful target.

Scripts: scratchpad (sizing.py, expand.py, solve.py); flint fmpz_mpoly expansion ~200 s,
all rank tests < 1 s. No existing files modified; nothing committed.
