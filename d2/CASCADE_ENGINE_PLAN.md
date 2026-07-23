# Lower-cascade valuation engine

**Date:** 2026-07-22
**Purpose:** replace one-stratum proofs by exact cone-level elimination across
the 420 live T1/T2 branch records in the geometric split-place ledger.

## Why this is now the primary program

The field repair enlarged the honest f31/subcase-(2) frontier from 21 uniform
rows to 327 geometric multiplicity vectors. Terminal tests and correctly
scoped proofs leave 235 strata and 420 live branches. Those records are not
420 unrelated systems: they share one cascade, one small list of `h_f`
monomials, four symmetric roots of `q`, and global degree caps.

The post-repair exact proofs already exhibit the reusable transitions:

```text
terminal UFD/parity
  -> minimum-valuation obstruction at level 6 or 5
  -> a square/residue condition when the minimum ties
  -> lower-level propagation
  -> infinity or degree contradiction.
```

The a_t=9 T1 work is the regression pilot. Its constant cell is closed by a
coefficient cascade; its nonconstant cell produces alternating square
conditions and is currently forced to `s^6|sigma` and `s^3|R`.

## 1. Cascade signature

The engine should consume a case description rather than a hand-expanded
resultant. For one master component record

```text
signature = {
  forcing_divisor: [(place, multiplicity, residue_degree)],
  master_exponents: p_f,
  cascade_step: v_t(Phi)-3*v_t(e),
  window_caps: deg(d2), deg(d1), deg(d0), deg(e),
  h_degree_caps: deg(h_f),
  h_monomials: h_f rewritten in (d2,d1,sigma,e),
  terminal_branches: T1, T2, T3,
}
```

For f31 the terminal ladder is

```text
h7 = 8192*d1^2,
h6 = -3072*sigma^2 + 14336*d1^2*d2 + 8192*d1*e,
h5 = -9216*d2*sigma^2 + 32256*d1^2*sigma
     -12288*d1^2*d2^2 + 18432*d1*d2*e + 2048*e^2.
```

`h4`, and later `h3`, are still sparse after replacing
`d0=(d2^2+sigma)/4`. The extractor must derive these monomials directly from
`f31_graded.txt`; no copied formula is trusted without a source assertion.

## 2. Local state and exact transition semantics

At a simple split place `p|q`, use a state

```text
(b, x, z, k, r_l, zero_flags, residue_constraints)
```

with

```text
b=v_p(E), x=v_p(d1), z=v_p(sigma), k=v_p(d2), r_l=v_p(g_l).
```

Zero polynomials use an explicit infinity value, not an artificial large
finite valuation. A local polynomial can vanish only if its minimum monomial
valuation occurs at least twice. Therefore each transition has two stages:

1. **tropical pruning:** reject a unique minimum;
2. **residue transition:** when minima tie, record the leading-coefficient
   equation and allow the valuation to rise only through a certified residue
   cancellation.

The cascade transition at level `l` is

```text
t^v*g_(l+1) = E^3*g_l + (c*q)^l*h_l.
```

At a root of `q`, `t` is a unit, so the candidate valuations are

```text
r_(l+1), 3*b+r_l, l+v_p(h_l).
```

The engine must not replace a tied minimum by an arbitrary larger valuation.
It should attach the small residue equation that permits the rise. The exact
squares found in `t5_90t1_local_verify.py` are the model behavior.

## 3. Global coupling

Four local states are coupled by polynomial degree:

```text
sum_i v_(p_i)(d1)    <= deg(d1),
sum_i v_(p_i)(sigma) <= deg(sigma),
sum_i v_(p_i)(d2)    <= deg(d2),
sum_i v_(p_i)(g_l)   <= deg(g_l).
```

Infinity is a fifth place. It supplies exact degree minima/maxima and detects
unique leading terms. The implementation should combine the four finite
places by dynamic programming on the accumulated degree budgets, quotienting
by permutation symmetry of the four roots. It should never enumerate labeled
four-tuples when their sorted valuation profile is sufficient.

## 4. Cone output, not row output

Each proof result should describe a region such as

```text
a in [a0,a1], sorted b satisfying affine inequalities,
branch=T1, local transition pattern=P,
contradiction=forced_degree > cap.
```

The human-facing artifact should contain ten or twenty such reusable lemmas,
plus a finite exceptional list with attached residue systems. The JSON ledger
can retain row-level witnesses, but every killed row must point to a cone
certificate and an independently checked inequality.

## 5. Implementation phases

### Phase A — signature extractor and regression suite

- Parse `h_0,...,h_7` from `f31_graded.txt`.
- Rewrite them in `(d2,d1,sigma,e)`.
- Emit monomial exponent tables and degree caps.
- Replay the existing a=7, a=9 T2, a=9 T1, and sigma-locus proofs.
- Treat zero-polynomial and cancellation edges explicitly.

### Phase B — levels 7 through 4 on the full ledger

- Start from every terminal-feasible branch in `split_place_ledger.json`.
- Propagate local valuation states through `h6`, `h5`, and `h4`.
- Couple the four roots using degree-budget dynamic programming.
- Report eliminated cones, surviving cones, and why each survivor needs a
  residue equation.

### Phase C — residue systems

- Normalize tied leading coefficients by units.
- Factor each residue equation before introducing new variables.
- Recognize repeated square/parity templates.
- Export only the exceptional systems that cannot be decided tropically.

### Phase D — levels 3 through 0 and infinity

- Continue only surviving cones.
- Add infinity-degree transitions at the earliest level where they are sharp.
- Use restored pre-resultant equations when a resultant component has excess
  solutions, following the f37 free-family method.

#### Phase D implementation design (2026-07-22, post-audit)

Infinity is a sixth place with `v_inf = -deg`, so the level identities get
the max-plus dual of the engine's min-plus semantics:

- `v_inf(t^v g_{l+1}) = -(v + deg g_{l+1})` exactly (product);
- for a sum: unique maximum degree forces the degree; a tie permits the
  degree to DROP only through leading-coefficient cancellation (dual of the
  tie-rise), recorded as an obligation with the tied leading terms;
- monomial degrees are exact given per-variable degree assignments
  `(deg d2, deg d1, deg sigma, deg e)`, enumerated within
  `sum_places v_p(x) <= deg x <= cap` — the valuation witnesses already
  bound each degree from below, closing the sandwich from both ends;
- `deg e` interacts with `a` and the b-vector: `deg e >= a + sum(b_i)`.

Implementation: reuse `descend_options` with negated valuations, or a
mirrored `descend_options_inf`; join infinity into the DFS with the shared
budget dims plus the degree-assignment choice. Regression targets, in
order: (i) replay the a=9 T2 uniform infinity kill (`T5_90_T2.md`), (ii)
replay the 43/50 constant-cell degree kills of `T5_90_T1.md`, (iii) the
degree-domination endgames of `T5_T2_COLUMN.md` (e.g. max(T0..T5)=234 <
deg(T6)=236), (iv) then sweep the 26 sub2 cells, the 279-branch sub1
family, and the 27 alternate branches. Kills remain subject to the
independent-audit gate (extend the spec-only checker with the same dual
semantics).

### Phase E — parametric signatures

After the current case closes, compare signatures for the first F_2 members
and the repeated `(8,28)` corner at degree 144. Keep exponents symbolic long
enough to determine whether the cone inequalities are affine in the family
parameter. This is the decision gate for a family/template program.

## 6. Correctness requirements

- Base-change before assigning prime-place valuations.
- Separate a zero polynomial from a high finite valuation.
- A tropical tie is necessary, not sufficient; every valuation rise needs a
  residue certificate.
- Source-link every `h_f` monomial table.
- Check cone certificates against brute row enumeration on the complete
  327-stratum ledger.
- Preserve an independent checker that does not call the generator's pruning
  routines.
- Keep infinity leading-coefficient cancellation separate from degree bounds.

## 7. Immediate target

Implement Phase A as a small `cascade_signature.py` module and tests. Then add
one level at a time to the ledger, using the a_t=9 local descent as the first
residue-transition regression. The first quantitative milestone is:

```text
terminal + levels 6/5/4 reduce 420 live branches to a certified cone list
with no row lacking either a contradiction or an explicit residue system.
```
