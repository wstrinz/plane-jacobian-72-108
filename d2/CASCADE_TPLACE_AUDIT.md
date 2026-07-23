# Independent t-place flag-case audit

**Date:** 2026-07-22  
**Checker:** `audit_tplace_cases.py`  
**Scope:** standard regime `a <= 10`, depth 4, subcase (2) and subcase (1)

## Result

The independent checker agrees with every flag-case claim in both q+t files.
It enumerated all legal global zero-flag cases on every open branch, including
cases absent from `survivor_cases` (which are claims of infeasibility).

| window | branches | all flag-cases | agreements | disagreements | q-only survivors | q+t claimed survivors | q+t audit survivors |
|---|---:|---:|---:|---:|---:|---:|---:|
| sub2 | 420 | 7,872 | 7,872 | 0 | 320 | 232 | 232 |
| sub1 | 2,178 | 41,592 | 41,592 | 0 | 2,519 | 2,253 | 2,253 |

Thus the coupled t place removes 88 of 320 q-feasible flag cases in sub2 and
266 of 2,519 in sub1, while leaving at least one feasible flag case on every
q-only surviving branch.

At branch level, the flag audit also agrees with the q+t files on 420/420
sub2 branches and 2,178/2,178 sub1 branches.

## Independent method

The checker was written from `CASCADE_ENGINE_REPORT.md`, the `T-PLACE
COUPLING` and `SUB1 T-PLACE COUPLING` entries in `STATE.md`, and the two prior
spec-only q-place checkers. It does not import or read `cascade_engine.py` or
its tests.

It independently parses `h_0,...,h_7` from `f31_graded.txt`, checks their
weighted homogeneity, substitutes

```text
d0 = (sigma + d2^2)/4,
```

and verifies the rewritten degree caps for both windows. The e-slot in every
rewritten h monomial is then charged `b_i` at a q-root and exactly `a` at t.

For each branch and every legal flag case, the checker exhausts integer local
orders for `d2`, `d1`, `sigma`, and all `g_l` from the terminal down through
level 4. The flag grids are:

- T1: both `d2_zero` values, both `sigma_zero` values, and all masks on
  `g_4,g_5,g_6`, giving 32 cases per branch;
- T2: both `d2_zero` values, `sigma_zero=False`, and all masks on `g_4,g_5`,
  giving 8 cases per branch.

The terminal g is never included in a zero mask. This is also why T2 does not
admit `sigma_zero=True`: its nonzero terminal is controlled by `sigma^2`.

At q, the checker uses

```text
r_(l+1) = ultrametric(3*b_i + r_l, l + v_i(h_l)).
```

At t, with `v=30-3a` and both `ehat` and `u` units, it uses

```text
v + s_(l+1) = ultrametric(s_l, w_l),
w_l = v_t(h_l).
```

The t terminals are imposed exactly as `s_7=2*v_t(d1)` for T1 and
`s_6=2*v_t(sigma)` for T2. The q terminals remain
`r_7=7+2*v_i(d1)-3*b_i` and `r_6=6+2*v_i(sigma)-3*b_i`.

All five places are joined in one DFS against the same global budget vector.
The variable caps are `(4,6,8)` for `(d2,d1,sigma)` in sub2 and `(6,9,12)`
in sub1. Sub2 uses the uniform g cap `10+3a`; the sub1 per-level caps are
re-derived from the same forward/backward recurrences as the prior independent
sub1 q-only checker. Pareto reduction only discards a local profile when
another profile for the same global flags uses no more of any shared budget.

The ultrametric policy is deliberately conservative: a unique minimum forces
the valuation; a tie may rise to any finite value up to the polynomial degree
cap or to infinity. If all monomials are killed by global zero flags, the h
value is exactly infinity. A tied minimum already above the degree cap can
only be infinity. This is the relaxed-tie convention used in the earlier
spec-only audits and can only make infeasibility harder to prove.

For each legal flag case, membership in the q+t file's `survivor_cases` list
is compared with the independent five-place decision in both directions:

- listed but audit-infeasible is a disagreement;
- absent but audit-feasible is a disagreement and prints a complete
  conservative five-place witness.

## Separate branch consistency check

The lightweight command

```text
python audit_tplace_cases.py --consistency-only
```

compares branch keys and verdicts in each q+t file with its already audited
q-only counterpart. Its separate output was:

```text
sub2_q_vs_qt_branch_consistency: PASS; agreements=420/420; mismatches=0
sub1_q_vs_qt_branch_consistency: PASS; agreements=2178/2178; mismatches=0
branch_consistency_overall: PASS; mismatches=0
```

No branch changes status under t coupling.

## Disagreements

None. There are no flag-case, branch-audit, file-structure, or q-versus-q+t
branch-consistency disagreements to enumerate.

## Reproduction and runtime

The end-to-end command was:

```text
python audit_tplace_cases.py --quiet
```

It exited 0 and reported:

```text
sub2_runtime_seconds: 32.620
sub1_runtime_seconds: 365.963
runtime_seconds: 398.998
```

The final measured runtime was **398.998 seconds** on the audit host. The
reported total includes parsing, cap checks, both exhaustive window audits,
all flag comparisons, and the in-run explicit branch consistency checks.
