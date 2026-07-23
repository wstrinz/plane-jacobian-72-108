# Independent cascade-kill audit: subcase (1)

**Verdict:** all audited claims agree. The independent verifier confirms all
1,899 claimed depth-4 kills and all 279 claimed depth-4 survivors. The
standard-regime terminal check agrees on 2,614/2,614 branch records. No branch
was skipped or left undecided.

## Independence and inputs

`audit_cascade_kills_sub1.py` is a port of the from-scratch subcase-(2) audit.
It neither imports nor reads `cascade_engine.py` or any test of that engine.
The audit used `f31_graded.txt`, `split_place_ledger_sub1.json`,
`cascade_cones_sub1_depth4.json`, and `cascade_cones_sub1_depth5.json`. The
port was derived from `CASCADE_KILL_AUDIT.md`, `sub1_cascade_verify.py`,
`SPLIT_PLACE_LEDGER_SUB1.md`, and the named 2026-07-22 entries in `STATE.md`.

The verifier parses all eight `h_f` expressions itself, checks their weighted
homogeneity, substitutes `d0=(sigma+d2^2)/4`, expands over the rationals, and
extracts exponent vectors in `(d2,d1,sigma,e)`. It checks both the source and
rewritten monomials against `deg h_f <= 60-6f`.

## Cap derivation and recheck

For each standard-regime `a=0,...,10`, the code redoes the stated recurrences
line by line:

1. `v=30-3a` and `deg(ehat)<=15-a`.
2. `forward[1]=60-v`.
3. For increasing levels,
   `forward[l+1]=max(3(15-a)+forward[l],60-2l)-v`.
4. The backward terminal anchor is 46 at T1 level 7 and 48 at T2 level 6.
5. For decreasing levels,
   `backward[l]=max(v+backward[l+1],60-2l)`.
6. The search cap is `min(forward[l],backward[l])` independently at each level.

The resulting caps at the audited levels are:

| a | T1 g4 | T1 g5 | T1 g6 | T1 g7 | T2 g4 | T2 g5 | T2 g6 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 75 | 90 | 76 | 46 | 75 | 78 | 48 |
| 1 | 78 | 93 | 73 | 46 | 78 | 75 | 48 |
| 2 | 81 | 94 | 70 | 46 | 81 | 72 | 48 |
| 3 | 84 | 88 | 67 | 46 | 84 | 69 | 48 |
| 4 | 87 | 82 | 64 | 46 | 84 | 66 | 48 |
| 5 | 90 | 76 | 61 | 46 | 78 | 63 | 48 |
| 6 | 82 | 70 | 58 | 46 | 72 | 60 | 48 |
| 7 | 73 | 64 | 55 | 46 | 66 | 57 | 48 |
| 8 | 64 | 58 | 52 | 46 | 60 | 54 | 48 |
| 9 | 55 | 52 | 49 | 46 | 54 | 51 | 48 |
| 10 | 52 | 50 | 48 | 46 | 52 | 50 | 48 |

The explicit ledger terminal caps and every finite stored survivor-witness
usage in the depth-4 claims respect these derived caps. **Cap discrepancies:
none.** In particular, the audit never substitutes the unsupported uniform
cap `15+3a`.

## Search semantics and completeness

At each q-place the verifier uses `v_i(ehat)=b_i`, `v_i(u)=1`, and treats `t`
as a unit. It enumerates the global zero flags for `d2`, for `sigma` in T1,
and for every processed nonterminal `g_l`. T2 has `d1=0`; T1 has nonzero
`d1`; terminal `g` is nonzero.

The terminal identities are
`v_i(g7)=7+2v_i(d1)-3b_i` for T1 and
`v_i(g6)=6+2v_i(sigma)-3b_i` for T2. For each rewritten `h_l`, a unique
minimum fixes its order. A tied minimum permits every integer rise through
`60-6l` and identical vanishing. Structural zero forces vanishing. The
cascade equation uses the unchanged ultrametric rule for
`t^v g_(l+1)=ehat^3 g_l+u^l h_l`.

The relaxed tie rule is deliberately conservative: no residue equation is
required, and local `h_l` cancellations may occur independently. Thus the
feasible set is enlarged, which can only make a claimed kill harder to
confirm.

All nonzero orders have finite caps. Per-place transitions are exhaustive;
an exact interval test replaces literal enumeration of every tied `h_l` rise.
Partial paths are dominance-reduced only when they have the same current
`g_l`, so their lower continuation problem is identical. Completed local
profiles are Pareto-reduced by `(d2,d1,sigma,g_l,...)` resource use. The four
places interact only through upper bounds on sums, so this reduction and the
memoized four-place join preserve existence exactly.

## Terminal results

There are 1,307 ledger strata with `a<=10`, hence 2,614 T1/T2 terminal branch
records. Result: **2,614/2,614 agree**. Disagreements: **none**.

## Depth-4 results

Result: **2,178/2,178 agree**: 1,899 killed agreements and 279 survivor
agreements. Disagreements: **none**. Undecided/resource-limited branches:
**none**.

| a | T1 killed | T1 survives | T2 killed | T2 survives |
|---:|---:|---:|---:|---:|
| 0 | 145 | 17 | 209 | 9 |
| 1 | 141 | 17 | 194 | 9 |
| 2 | 130 | 17 | 170 | 9 |
| 3 | 115 | 17 | 141 | 9 |
| 4 | 94 | 17 | 111 | 9 |
| 5 | 74 | 17 | 85 | 9 |
| 6 | 53 | 17 | 62 | 9 |
| 7 | 36 | 17 | 44 | 9 |
| 8 | 21 | 17 | 29 | 9 |
| 9 | 12 | 15 | 18 | 9 |
| 10 | 6 | 12 | 9 | 9 |

The script's default mode prints the complete 2,178-row agreement table with
branch id, claimed verdict, independent verdict, and agreement marker.
`--quiet` suppresses only that table.

## Depth-5 consistency

Both claim files contain the same 2,178 branch ids. All **1,606/1,606**
claimed depth-5 kills are contained in the claimed depth-4 kill set. Subset
violations: **none**.

## Runtime and reproduction

The complete quiet run on native Windows Python 3/SymPy took **162.995
seconds** measured inside the verifier (164.7 seconds including command
wrapper overhead):

```powershell
python audit_cascade_kills_sub1.py --quiet
```

The process exited 0. Any terminal disagreement, depth-4 disagreement, cap
discrepancy, or depth-5 subset violation makes it exit nonzero.
