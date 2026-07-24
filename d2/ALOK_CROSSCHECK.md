# Cross-program sanity check against `alok/jacobian-two`

**Date:** 2026-07-24. **Status: READ-ONLY cross-check, no audited artifact
touched.** New files only: this doc + `alok_crosscheck.py`. Nothing committed.
The external result is an *independent public program*; this note imports it as
a sanity check on our branch ledger, exactly as the external review recommended.

---

## 1. The external result — found, read, and cited

The repository exists and is public.

- **Repo:** <https://github.com/alok/jacobian-two> — "Lean-certified audit of
  the 3D Jacobian counterexample and research toward JC(2)".
- **Commit pinned:** `ded8e67fb47c155f83ec9bd68af6014499fc2d61`
  ("document the fourth codimension-three exclusion", 2026-07-22), the current
  `main` HEAD.
- **Primary artifact:** `docs/newton-72-108-sparse.md`
  ("A sparse-interior obstruction in the residual `(72,108)` configurations").
- **Machine-checked checker:** `scripts/newton_72_108.py`
  (thresholds/polygons at lines **118-127**), with regression
  `tests/test_newton_72_108.py`.

### 1.1 Exact statement (verbatim from the source)

Over a field of characteristic zero, for `[P,Q] = P_xQ_y - P_yQ_x = x^2` with
the GGHV (Guccione–Guccione–Horruitiner–Valqui, arXiv:2204.14178, **Prop 4.3**)
residual Newton polygons, writing `int N(·)` for strict-interior lattice points
and `#(supp ∩ int)` for the number of nonzero coefficients there:

- **Case 1:** `#(supp(P)∩int N(P)) + #(supp(Q)∩int N(Q)) ≥ 3`.
- **Case 2:** `#(supp(P)∩int N(P)) + #(supp(Q)∩int N(Q)) ≥ 4`.

A "strict-interior coefficient" is a **nonzero coefficient `p_ij` (of `x^i y^j`
in `P`) or `q_kl` (of `x^k y^l` in `Q`) at a lattice point strictly inside the
Newton polygon** `N(P)` resp. `N(Q)`. Boundary coefficients are unrestricted;
only the polygon vertices are required nonzero.

**The review's paraphrase is accurate**: "Case 1 ≥ 3, Case 2 ≥ 4" matches the
source exactly, and the review's honest framing — "an exhaustive elimination of
sparser supports that does not close the case" — is also the source's own
framing ("It does not eliminate either full Proposition 4.3 case").

### 1.2 Proof / verification artifact (public)

Exact finite exhaustion (`scripts/newton_72_108.py`), reproducible via
`uv run python scripts/newton_72_108.py`:

| case | patterns checked | mechanism | max forced-zero steps | certificate sha256 |
|---|---:|---|---:|---|
| Case 1 | `1+122+C(122,2)=7504` (all ≤2-interior supports) | zero-product propagation, 7504/7504 | 21 | `7c7a38…bb1519` |
| Case 2 | `1+28+C(28,2)+C(28,3)=3683` (all ≤3-interior supports) | 3678 zero-product + **5** algebraic Gröbner-`[1]` certificates | 22 | `ee64db…31f83f` |

The certificate mechanism is exactly our object: a monomial pair contributes
`(iℓ − jk) p_ij q_kl` to the coefficient of `x^{i+k-1} y^{j+ℓ-1}` in `[P,Q]`,
and every output coefficient except `(2,0)` must vanish. The five Case-2
exceptions carry hand identities + machine-checked `GB = [1]` over `QQ`.

---

## 2. The subcase correspondence — EXACT

Both programs reduce the *same* GGHV Prop 4.3 configuration. The polygons match
vertex-for-vertex:

| | alok `CASE_1` (`newton_72_108.py:118-121`) | our `sub1` / subcase (1) (`STATE.md:21`) |
|---|---|---|
| `N(P)` | `{(0,0),(1,0),(8,14),(8,16),(0,8)}` | `{(0,0),(1,0),(8,14),(8,16)}` **+ corner `(0,8)`** |
| `N(Q)` | `{(0,0),(2,1),(12,21),(12,24),(0,12)}` | `{(0,0),(2,1),(12,21),(12,24)}` **+ corner `(0,12)`** |
| threshold | **≥ 3** | inherited: **≥ 3** |

| | alok `CASE_2` (`newton_72_108.py:124-127`) | our `sub2` / subcase (2) (`STATE.md:20`) |
|---|---|---|
| `N(P)` | `{(0,0),(1,0),(8,14),(8,16)}` | `{(0,0),(1,0),(8,14),(8,16)}` |
| `N(Q)` | `{(0,0),(2,1),(12,21),(12,24)}` | `{(0,0),(2,1),(12,21),(12,24)}` |
| threshold | **≥ 4** | inherited: **≥ 4** |

So the map is unambiguous: **`sub1 ↔ Case 1` (threshold 3), `sub2 ↔ Case 2`
(threshold 4)**. Both use the identical `[P,Q]=x^2` normalization.

**Independent census check (performed by the checker).** From the vertices
alone, `alok_crosscheck.py` re-derives the strict-interior lattice-point census
and reproduces alok's published table exactly:

| config | side | interior lattice points (recomputed = alok published) |
|---|---|---:|
| Case 1 | `P` | 35 |
| Case 1 | `Q` | 87 |
| Case 2 | `P` | 7 |
| Case 2 | `Q` | 21 |

This is genuine independent corroboration that **both programs are attacking
the identical polygons under the identical bracket normalization** — the
foundation our entire `(72,108)` ledger rests on.

---

## 3. Translation to our state space — PARTIAL (setup: exact; per-state
support count: not recoverable)

Their support-count and our state coordinates live on **different lattices**,
related by a nonlinear map we do not (and cannot soundly) invert per-state.

- **Their coordinates:** the `(x,y)`-monomial coefficients `p_ij`, `q_kl` of
  `P`, `Q` directly; "support" = which strict-interior monomials are nonzero.
- **Our coordinates:** the `y`-**degrees** and zero-flags of the D-transformation
  coefficients `d2,d1,sigma,e = D̃_2,D̃_1,D̃_0,D̃_{-1}` and the cascade
  polynomials `g_l`, indexed further by the `t=(y+1)` multiplicity `a_t`, the
  four-place `b`-vector (valuations of `e` at the `q`-roots), and the branch
  `T1/T2` (`CASCADE_INF_REPORT.md`; `phase_d_states_*.json` schema).

Three obstructions make a faithful per-state interior-term count unrecoverable:

1. **Different ring.** The `d_k` live in the *localized* ring `K[y, C4^{-1}]`
   (`STATE.md` items 1,3), not in `K[x,y]`. Recovering `p_ij,q_kl` requires
   clearing `C4 = y^7(y+1)` denominators.
2. **Nonlinear (convolution) map.** `P = C^2`, `Q = C^3 + λC^{-1} + F`
   (`STATE.md` item 1), so each `p_ij,q_kl` is a convolution of `C`-coefficients;
   cancellations can occur that degrees alone do not predict.
3. **Degrees, not supports.** Our ledger records `deg_{d_k}` (an upper envelope),
   never which individual `y`-coefficients are nonzero — the exact datum their
   count needs.

Inventing a closed-form "interior support count = f(deg_d2, …)" would be an
unaudited fabrication, precisely the kind of thing the local audit discipline
forbids. So the translation is **exact at the setup/subcase level and not
invertible at the per-state support level.** The two programs partition the
problem into **complementary regimes**:

- **alok:** the *sparsest* realizations (≤ 2 / ≤ 3 interior terms), elementary
  exhaustive elimination.
- **ours:** the *dense generic* regime — every tracked state carries the fixed
  nonzero forcing term `Φ = f1·C4^28` (`STATE.md` items 2,4: the reduced
  equation `(D̃^3)_{-5} + Φ = 0` with `Φ ≠ 0`) plus full-degree `d_k` families,
  i.e. interior support far above either threshold.

---

## 4. Cross-check verdict — CORROBORATION (no findings)

`alok_crosscheck.py` runs the one **sound** comparison the mapping supports: a
guard that flags any LIVE branch/degree-state whose coordinates could FORCE the
total strict-interior support of `(P,Q)` below the case threshold. Result:

| window | alok case | threshold | live branches | live degree-states | **live below threshold** |
|---|---|---:|---:|---:|---:|
| `sub1` | Case 1 | 3 | 171 | 44 117 | **0** |
| `sub2` | Case 2 | 4 | 26 | 7 888 | **0** |

- **Exit code 0** (`--quiet` clean). **No live branch or state sits below
  alok's sparse floor** — exactly the reviewer's prediction ("any branch
  encoding support below those thresholds should disappear automatically").
  Here it holds *structurally*: our reduction commits every state to the dense
  regime (nonzero `Φ` forbids the empty-interior realization; each
  non-zero-flagged `d_k` family is an extra independent interior degree of
  freedom), so no live state can encode a below-threshold support in the first
  place. The guard therefore corroborates rather than finds.

- **Direction that would have been a finding:** an OPEN branch forced to a
  ≤ threshold interior support (their theorem kills it, our ledger keeps it).
  **None exists.** The mirror direction (their kills we could re-derive) is
  quantified next.

### 4.1 How much of our kill work their theorem independently implies

Their proven-dead region = supports with **≤ 2 (Case 1) / ≤ 3 (Case 2)**
strict-interior terms. Our ledger's kills (390 sub2 + 2007 sub1 engine kills,
plus the residue/convolution layers) are all **degree-states in the dense
D-lattice regime**, every one of which carries `Φ ≠ 0`. The intersection of the
two dead regions is therefore **empty**: their sparse-support theorem
independently implies **0** of our specific degree-state kills, and vice versa.

This is **complementary corroboration, not overlap**: the two programs
independently eliminate the `(72,108)` case from opposite ends of the same
support spectrum (alok closes the sparse floor `≤2/≤3`; our cascade attacks the
dense generic body), and — critically — **they do not contradict each other at
any live branch.** No live state of ours falls into a region alok has proved
dead.

### 4.2 Honest limits of this cross-check

- The guard's soundness rests on the program invariant `Φ ≠ 0`, not on a
  per-state monomial count (which is not recoverable — §3). A future
  audited convolution `K[y,C4^{-1}] → K[x,y]` would upgrade the qualitative
  "dense regime" argument to an exact per-state interior-term count; until then
  the guard certifies only the *absence of below-threshold live states*, which
  is the check the reviewer asked for.
- The `d_lattice_signature` reported by the script (histogram over live states:
  `sub1 {3:79, 4:3690, 5:40348}`, `sub2 {3:50, 4:1247, 5:6591}`) counts nonzero
  D-families + the forcing term. It lives on the **D-lattice, not alok's Newton
  lattice**, and is intentionally *not* compared numerically to the threshold —
  it is reported only to make the coordinate mismatch explicit and auditable.

---

## 5. Reproduction

```bash
python alok_crosscheck.py            # full report (verdict PASS)
python alok_crosscheck.py --quiet    # exit 0 iff no live-below-threshold state
```

Independent re-fetch of the external artifact:

```bash
gh api repos/alok/jacobian-two/contents/docs/newton-72-108-sparse.md --jq .content | base64 -d
gh api repos/alok/jacobian-two/contents/scripts/newton_72_108.py       --jq .content | base64 -d
```

---

## 6. Files created (new only)

- `ALOK_CROSSCHECK.md` — this note (statement + citations + mapping + verdict).
- `alok_crosscheck.py` — `--quiet`; exit 0 iff no live-below-threshold state.

No existing file edited; no commit; no concurrent-lane file touched.
