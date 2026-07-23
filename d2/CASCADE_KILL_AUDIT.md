# Independent cascade-kill audit

**Verdict:** all claims agree. The independent verifier confirms all 390 claimed depth-4 kills and finds all 30 claimed survivors feasible under the deliberately relaxed tie semantics. The terminal-only sanity check agrees with all 654 ledger records.

## Inputs and independence

The verifier reads only `f31_graded.txt`, `split_place_ledger.json`, and `cascade_cones.json` at runtime. It was written from the semantics in `CASCADE_ENGINE_REPORT.md`, with the cascade algebra cross-checked against `t5_multiplace_verify.py` checks 5-7 and `T5_NP.md`. It neither imports nor reads `cascade_engine.py`, `test_cascade_engine.py`, or `cascade_signature.py`.

It parses all eight `h_l` expressions directly from the documented `h_N (weight W, dm1-power P) = expr` records. Before any valuation search, it checks every original monomial has weighted degree `20-2l` for weights `(2,3,4,5)` on `(d2,d1,d0,dm1)`. It then independently substitutes `d0=(sigma+d2^2)/4`, expands and combines terms over the rationals, and extracts exponent vectors in `(d2,d1,sigma,e)`. This self-check passed for `h_0,...,h_7`.

## Implemented semantics

At place `p_i`, `v_i(e)=b_i` and `v_i(u)=1`; `t` is a unit. Nonzero polynomial orders are nonnegative integers. Global zero flags are enumerated for `d2`, for `sigma` only in T1, and for every nonterminal `g_l`. T2 has `d1` identically zero; T1 has nonzero `d1`. Terminal `g` is always nonzero.

The terminal equations are implemented as `v_i(g_7)=7+2v_i(d1)-3b_i` for T1 and `v_i(g_6)=6+2v_i(sigma)-3b_i` for T2. At each lower level, every surviving rewritten monomial of `h_l` receives its exact linear valuation. A unique minimum fixes `v_i(h_l)`. A tied minimum permits every integer rise through `40-4l` and also identical vanishing. If zero flags remove every monomial, `h_l` is forced to vanish; if exactly one monomial survives, no rise or vanishing is allowed.

For `t^v g_{l+1}=e^3 g_l+u^l h_l`, the right-side term orders are `3b_i+v_i(g_l)` and `l+v_i(h_l)`. A unique minimum equals `v_i(g_{l+1})`; a tie can rise arbitrarily. If `g_{l+1}` is identically zero, the two term orders must agree exactly (including the case in which both terms vanish). If `g_l` is identically zero, the second term alone fixes `g_{l+1}`.

The four local chains are joined subject to `sum v_i(d2)<=4`, `sum v_i(d1)<=6`, `sum v_i(sigma)<=8`, and separately `sum v_i(g_l)<=10+3a` at each processed nonzero level. Identically zero polynomials consume no degree budget.

### Conservative interpretations

- Tie cancellation is accepted solely from valuation multiplicity; no residue equation is imposed. This is the relaxation requested for a kill audit.

- An `h_l` tied minimum is allowed to rise independently at each place, including to infinity. Literal global compatibility of an identity `h_l=0` is therefore not imposed across places. This enlarges the feasible set beyond actual polynomial assignments and can only make a kill harder to confirm.

- When all monomials are removed by zero flags, `h_l=0` is treated as forced. The "at least two surviving monomials" restriction is applied to cancellation of surviving terms, not to this structural zero case.

## Completeness of the search

All nonzero orders are finite-range integers: base-variable orders are bounded by their degree caps, every nonzero `g_l` order is bounded by `10+3a`, and every nonzero `h_l` order is bounded by `40-4l`. The verifier exhausts every global zero-flag combination and every local valuation transition. Completed local chains are Pareto-reduced only by their vector of consumed global degree budgets; this is exact because completed places interact solely through upper bounds on sums. A memoized depth-first join then exhausts the four places. Thus "killed" is reported only after every assignment in this relaxed finite space is exhausted.

## Terminal cross-check

The terminal checker independently exhausts the four auxiliary orders and evaluates the terminal formula and degree sums. Result: **654/654 agree**, comprising all 327 strata times T1/T2. Disagreements: **none**. This check runs before the descent and the program refuses to report cascade results if it fails.

## Depth-4 results

- Open records processed: **420**

- Agreement: **420/420**

- Confirmed claimed kills: **390/390**

- Confirmed claimed survivors: **30/30**

- Disagreements: **none**

### Full killed-agreement list (390)

- `a=0`: `a=0;b=2,1,1,1;T2`; `a=0;b=2,2,1,0;T2`; `a=0;b=2,2,1,1;T1`; `a=0;b=2,2,1,1;T2`; `a=0;b=2,2,2,0;T1`; `a=0;b=2,2,2,0;T2`; `a=0;b=2,2,2,1;T1`; `a=0;b=2,2,2,1;T2`; `a=0;b=3,1,1,1;T2`; `a=0;b=3,2,1,0;T2`; `a=0;b=3,2,1,1;T1`; `a=0;b=3,2,1,1;T2`; `a=0;b=3,2,2,0;T1`; `a=0;b=3,2,2,0;T2`; `a=0;b=3,2,2,1;T1`; `a=0;b=3,2,2,1;T2`; `a=0;b=3,2,2,2;T1`; `a=0;b=3,2,2,2;T2`; `a=0;b=3,3,1,1;T1`; `a=0;b=3,3,1,1;T2`; `a=0;b=3,3,2,0;T1`; `a=0;b=3,3,2,0;T2`; `a=0;b=3,3,2,1;T1`; `a=0;b=3,3,2,1;T2`; `a=0;b=3,3,2,2;T1`; `a=0;b=3,3,2,2;T2`; `a=0;b=3,3,3,0;T1`; `a=0;b=3,3,3,0;T2`; `a=0;b=3,3,3,1;T1`; `a=0;b=3,3,3,1;T2`; `a=0;b=4,1,1,1;T2`; `a=0;b=4,2,1,0;T2`; `a=0;b=4,2,1,1;T1`; `a=0;b=4,2,1,1;T2`; `a=0;b=4,2,2,0;T1`; `a=0;b=4,2,2,0;T2`; `a=0;b=4,2,2,1;T1`; `a=0;b=4,2,2,1;T2`; `a=0;b=4,2,2,2;T1`; `a=0;b=4,2,2,2;T2`; `a=0;b=4,3,1,0;T2`; `a=0;b=4,3,1,1;T1`; `a=0;b=4,3,1,1;T2`; `a=0;b=4,3,2,0;T1`; `a=0;b=4,3,2,0;T2`; `a=0;b=4,3,2,1;T1`; `a=0;b=4,3,2,1;T2`; `a=0;b=4,3,3,0;T1`; `a=0;b=4,3,3,0;T2`; `a=0;b=4,4,1,0;T2`; `a=0;b=4,4,1,1;T1`; `a=0;b=4,4,1,1;T2`; `a=0;b=4,4,2,0;T1`; `a=0;b=4,4,2,0;T2`; `a=0;b=5,1,1,1;T2`; `a=0;b=5,2,1,0;T2`; `a=0;b=5,2,1,1;T1`; `a=0;b=5,2,1,1;T2`; `a=0;b=5,2,2,0;T1`; `a=0;b=5,2,2,0;T2`; `a=0;b=5,2,2,1;T1`; `a=0;b=5,2,2,1;T2`; `a=0;b=5,3,1,1;T1`; `a=0;b=5,3,1,1;T2`; `a=0;b=5,3,2,0;T1`; `a=0;b=5,3,2,0;T2`; `a=0;b=5,4,1,0;T2`; `a=0;b=6,1,1,1;T2`; `a=0;b=6,2,1,0;T2`; `a=0;b=6,2,1,1;T1`; `a=0;b=6,2,1,1;T2`; `a=0;b=6,2,2,0;T1`; `a=0;b=6,2,2,0;T2`; `a=0;b=6,3,1,0;T2`; `a=0;b=7,1,1,1;T2`; `a=0;b=7,2,1,0;T2`

- `a=1`: `a=1;b=2,1,1,0;T2`; `a=1;b=2,1,1,1;T1`; `a=1;b=2,1,1,1;T2`; `a=1;b=2,2,0,0;T2`; `a=1;b=2,2,1,0;T1`; `a=1;b=2,2,1,0;T2`; `a=1;b=2,2,1,1;T1`; `a=1;b=2,2,1,1;T2`; `a=1;b=2,2,2,0;T1`; `a=1;b=2,2,2,0;T2`; `a=1;b=2,2,2,1;T1`; `a=1;b=2,2,2,1;T2`; `a=1;b=3,1,1,0;T2`; `a=1;b=3,1,1,1;T1`; `a=1;b=3,1,1,1;T2`; `a=1;b=3,2,0,0;T2`; `a=1;b=3,2,1,0;T1`; `a=1;b=3,2,1,0;T2`; `a=1;b=3,2,1,1;T1`; `a=1;b=3,2,1,1;T2`; `a=1;b=3,2,2,0;T1`; `a=1;b=3,2,2,0;T2`; `a=1;b=3,2,2,1;T1`; `a=1;b=3,2,2,1;T2`; `a=1;b=3,2,2,2;T1`; `a=1;b=3,2,2,2;T2`; `a=1;b=3,3,1,0;T1`; `a=1;b=3,3,1,0;T2`; `a=1;b=3,3,1,1;T1`; `a=1;b=3,3,1,1;T2`; `a=1;b=3,3,2,0;T1`; `a=1;b=3,3,2,0;T2`; `a=1;b=3,3,2,1;T1`; `a=1;b=3,3,2,1;T2`; `a=1;b=3,3,3,0;T1`; `a=1;b=3,3,3,0;T2`; `a=1;b=4,1,1,0;T2`; `a=1;b=4,1,1,1;T1`; `a=1;b=4,1,1,1;T2`; `a=1;b=4,2,0,0;T2`; `a=1;b=4,2,1,0;T1`; `a=1;b=4,2,1,0;T2`; `a=1;b=4,2,1,1;T1`; `a=1;b=4,2,1,1;T2`; `a=1;b=4,2,2,0;T1`; `a=1;b=4,2,2,0;T2`; `a=1;b=4,2,2,1;T1`; `a=1;b=4,2,2,1;T2`; `a=1;b=4,3,0,0;T2`; `a=1;b=4,3,1,0;T1`; `a=1;b=4,3,1,0;T2`; `a=1;b=4,3,1,1;T1`; `a=1;b=4,3,1,1;T2`; `a=1;b=4,3,2,0;T1`; `a=1;b=4,3,2,0;T2`; `a=1;b=4,4,0,0;T2`; `a=1;b=4,4,1,0;T1`; `a=1;b=4,4,1,0;T2`; `a=1;b=5,1,1,0;T2`; `a=1;b=5,1,1,1;T1`; `a=1;b=5,1,1,1;T2`; `a=1;b=5,2,0,0;T2`; `a=1;b=5,2,1,0;T1`; `a=1;b=5,2,1,0;T2`; `a=1;b=5,2,1,1;T1`; `a=1;b=5,2,1,1;T2`; `a=1;b=5,2,2,0;T1`; `a=1;b=5,2,2,0;T2`; `a=1;b=5,3,1,0;T1`; `a=1;b=5,3,1,0;T2`; `a=1;b=5,4,0,0;T2`; `a=1;b=6,1,1,0;T2`; `a=1;b=6,1,1,1;T1`; `a=1;b=6,1,1,1;T2`; `a=1;b=6,2,0,0;T2`; `a=1;b=6,2,1,0;T1`; `a=1;b=6,2,1,0;T2`; `a=1;b=6,3,0,0;T2`; `a=1;b=7,1,1,0;T2`; `a=1;b=7,2,0,0;T2`

- `a=2`: `a=2;b=1,1,1,0;T2`; `a=2;b=2,1,0,0;T2`; `a=2;b=2,1,1,0;T1`; `a=2;b=2,1,1,0;T2`; `a=2;b=2,1,1,1;T1`; `a=2;b=2,1,1,1;T2`; `a=2;b=2,2,0,0;T1`; `a=2;b=2,2,0,0;T2`; `a=2;b=2,2,1,0;T1`; `a=2;b=2,2,1,0;T2`; `a=2;b=2,2,1,1;T1`; `a=2;b=2,2,1,1;T2`; `a=2;b=2,2,2,0;T1`; `a=2;b=2,2,2,0;T2`; `a=2;b=2,2,2,1;T1`; `a=2;b=2,2,2,1;T2`; `a=2;b=3,1,0,0;T2`; `a=2;b=3,1,1,0;T1`; `a=2;b=3,1,1,0;T2`; `a=2;b=3,1,1,1;T1`; `a=2;b=3,1,1,1;T2`; `a=2;b=3,2,0,0;T1`; `a=2;b=3,2,0,0;T2`; `a=2;b=3,2,1,0;T1`; `a=2;b=3,2,1,0;T2`; `a=2;b=3,2,1,1;T1`; `a=2;b=3,2,1,1;T2`; `a=2;b=3,2,2,0;T1`; `a=2;b=3,2,2,0;T2`; `a=2;b=3,2,2,1;T1`; `a=2;b=3,2,2,1;T2`; `a=2;b=3,3,0,0;T1`; `a=2;b=3,3,0,0;T2`; `a=2;b=3,3,1,0;T1`; `a=2;b=3,3,1,0;T2`; `a=2;b=3,3,1,1;T1`; `a=2;b=3,3,1,1;T2`; `a=2;b=3,3,2,0;T1`; `a=2;b=3,3,2,0;T2`; `a=2;b=4,1,0,0;T2`; `a=2;b=4,1,1,0;T1`; `a=2;b=4,1,1,0;T2`; `a=2;b=4,1,1,1;T1`; `a=2;b=4,1,1,1;T2`; `a=2;b=4,2,0,0;T1`; `a=2;b=4,2,0,0;T2`; `a=2;b=4,2,1,0;T1`; `a=2;b=4,2,1,0;T2`; `a=2;b=4,2,1,1;T1`; `a=2;b=4,2,1,1;T2`; `a=2;b=4,2,2,0;T1`; `a=2;b=4,2,2,0;T2`; `a=2;b=4,3,0,0;T1`; `a=2;b=4,3,0,0;T2`; `a=2;b=4,3,1,0;T1`; `a=2;b=4,3,1,0;T2`; `a=2;b=4,4,0,0;T1`; `a=2;b=4,4,0,0;T2`; `a=2;b=5,1,0,0;T2`; `a=2;b=5,1,1,0;T1`; `a=2;b=5,1,1,0;T2`; `a=2;b=5,1,1,1;T1`; `a=2;b=5,1,1,1;T2`; `a=2;b=5,2,0,0;T1`; `a=2;b=5,2,0,0;T2`; `a=2;b=5,2,1,0;T1`; `a=2;b=5,2,1,0;T2`; `a=2;b=5,3,0,0;T1`; `a=2;b=5,3,0,0;T2`; `a=2;b=6,1,0,0;T2`; `a=2;b=6,1,1,0;T1`; `a=2;b=6,1,1,0;T2`; `a=2;b=6,2,0,0;T1`; `a=2;b=6,2,0,0;T2`; `a=2;b=7,1,0,0;T2`

- `a=3`: `a=3;b=1,1,0,0;T2`; `a=3;b=1,1,1,0;T1`; `a=3;b=1,1,1,0;T2`; `a=3;b=2,0,0,0;T2`; `a=3;b=2,1,0,0;T1`; `a=3;b=2,1,0,0;T2`; `a=3;b=2,1,1,0;T1`; `a=3;b=2,1,1,0;T2`; `a=3;b=2,1,1,1;T1`; `a=3;b=2,1,1,1;T2`; `a=3;b=2,2,0,0;T1`; `a=3;b=2,2,0,0;T2`; `a=3;b=2,2,1,0;T1`; `a=3;b=2,2,1,0;T2`; `a=3;b=2,2,1,1;T1`; `a=3;b=2,2,1,1;T2`; `a=3;b=2,2,2,0;T1`; `a=3;b=2,2,2,0;T2`; `a=3;b=2,2,2,1;T1`; `a=3;b=2,2,2,1;T2`; `a=3;b=3,0,0,0;T2`; `a=3;b=3,1,0,0;T1`; `a=3;b=3,1,0,0;T2`; `a=3;b=3,1,1,0;T1`; `a=3;b=3,1,1,0;T2`; `a=3;b=3,1,1,1;T1`; `a=3;b=3,1,1,1;T2`; `a=3;b=3,2,0,0;T1`; `a=3;b=3,2,0,0;T2`; `a=3;b=3,2,1,0;T1`; `a=3;b=3,2,1,0;T2`; `a=3;b=3,2,1,1;T1`; `a=3;b=3,2,1,1;T2`; `a=3;b=3,2,2,0;T1`; `a=3;b=3,2,2,0;T2`; `a=3;b=3,3,0,0;T1`; `a=3;b=3,3,0,0;T2`; `a=3;b=3,3,1,0;T1`; `a=3;b=3,3,1,0;T2`; `a=3;b=4,0,0,0;T2`; `a=3;b=4,1,0,0;T1`; `a=3;b=4,1,0,0;T2`; `a=3;b=4,1,1,0;T1`; `a=3;b=4,1,1,0;T2`; `a=3;b=4,1,1,1;T1`; `a=3;b=4,1,1,1;T2`; `a=3;b=4,2,0,0;T1`; `a=3;b=4,2,0,0;T2`; `a=3;b=4,2,1,0;T1`; `a=3;b=4,2,1,0;T2`; `a=3;b=4,3,0,0;T1`; `a=3;b=4,3,0,0;T2`; `a=3;b=5,0,0,0;T2`; `a=3;b=5,1,0,0;T1`; `a=3;b=5,1,0,0;T2`; `a=3;b=5,1,1,0;T1`; `a=3;b=5,1,1,0;T2`; `a=3;b=5,2,0,0;T1`; `a=3;b=5,2,0,0;T2`; `a=3;b=6,0,0,0;T2`; `a=3;b=6,1,0,0;T1`; `a=3;b=6,1,0,0;T2`; `a=3;b=7,0,0,0;T2`

- `a=4`: `a=4;b=1,0,0,0;T2`; `a=4;b=1,1,0,0;T1`; `a=4;b=1,1,0,0;T2`; `a=4;b=1,1,1,0;T1`; `a=4;b=1,1,1,0;T2`; `a=4;b=2,0,0,0;T1`; `a=4;b=2,0,0,0;T2`; `a=4;b=2,1,0,0;T1`; `a=4;b=2,1,0,0;T2`; `a=4;b=2,1,1,0;T1`; `a=4;b=2,1,1,0;T2`; `a=4;b=2,1,1,1;T1`; `a=4;b=2,1,1,1;T2`; `a=4;b=2,2,0,0;T1`; `a=4;b=2,2,0,0;T2`; `a=4;b=2,2,1,0;T1`; `a=4;b=2,2,1,0;T2`; `a=4;b=2,2,1,1;T1`; `a=4;b=2,2,1,1;T2`; `a=4;b=2,2,2,0;T1`; `a=4;b=2,2,2,0;T2`; `a=4;b=3,0,0,0;T1`; `a=4;b=3,0,0,0;T2`; `a=4;b=3,1,0,0;T1`; `a=4;b=3,1,0,0;T2`; `a=4;b=3,1,1,0;T1`; `a=4;b=3,1,1,0;T2`; `a=4;b=3,1,1,1;T1`; `a=4;b=3,1,1,1;T2`; `a=4;b=3,2,0,0;T1`; `a=4;b=3,2,0,0;T2`; `a=4;b=3,2,1,0;T1`; `a=4;b=3,2,1,0;T2`; `a=4;b=3,3,0,0;T1`; `a=4;b=3,3,0,0;T2`; `a=4;b=4,0,0,0;T1`; `a=4;b=4,0,0,0;T2`; `a=4;b=4,1,0,0;T1`; `a=4;b=4,1,0,0;T2`; `a=4;b=4,1,1,0;T1`; `a=4;b=4,1,1,0;T2`; `a=4;b=4,2,0,0;T1`; `a=4;b=4,2,0,0;T2`; `a=4;b=5,0,0,0;T1`; `a=4;b=5,0,0,0;T2`; `a=4;b=5,1,0,0;T1`; `a=4;b=5,1,0,0;T2`; `a=4;b=6,0,0,0;T1`; `a=4;b=6,0,0,0;T2`

- `a=5`: `a=5;b=1,0,0,0;T1`; `a=5;b=1,1,0,0;T1`; `a=5;b=1,1,0,0;T2`; `a=5;b=1,1,1,0;T2`; `a=5;b=2,0,0,0;T1`; `a=5;b=2,0,0,0;T2`; `a=5;b=2,1,0,0;T1`; `a=5;b=2,1,0,0;T2`; `a=5;b=2,1,1,0;T1`; `a=5;b=2,1,1,0;T2`; `a=5;b=2,1,1,1;T1`; `a=5;b=2,1,1,1;T2`; `a=5;b=2,2,0,0;T1`; `a=5;b=2,2,0,0;T2`; `a=5;b=2,2,1,0;T1`; `a=5;b=2,2,1,0;T2`; `a=5;b=3,0,0,0;T1`; `a=5;b=3,0,0,0;T2`; `a=5;b=3,1,0,0;T1`; `a=5;b=3,1,0,0;T2`; `a=5;b=3,1,1,0;T2`; `a=5;b=3,2,0,0;T1`; `a=5;b=3,2,0,0;T2`; `a=5;b=4,0,0,0;T1`; `a=5;b=4,0,0,0;T2`; `a=5;b=4,1,0,0;T1`; `a=5;b=4,1,0,0;T2`; `a=5;b=5,0,0,0;T1`; `a=5;b=5,0,0,0;T2`

- `a=6`: `a=6;b=2,0,0,0;T1`; `a=6;b=2,0,0,0;T2`; `a=6;b=2,1,0,0;T1`; `a=6;b=2,1,0,0;T2`; `a=6;b=2,1,1,0;T1`; `a=6;b=2,1,1,0;T2`; `a=6;b=2,2,0,0;T1`; `a=6;b=2,2,0,0;T2`; `a=6;b=3,0,0,0;T2`; `a=6;b=3,1,0,0;T2`; `a=6;b=4,0,0,0;T1`; `a=6;b=4,0,0,0;T2`

- `a=7`: `a=7;b=2,0,0,0;T1`; `a=7;b=2,0,0,0;T2`; `a=7;b=2,1,0,0;T1`; `a=7;b=2,1,0,0;T2`

- `a=8`: `a=8;b=2,0,0,0;T1`; `a=8;b=2,0,0,0;T2`

### Full survivor-agreement list (30)

- `a=5`: `a=5;b=1,0,0,0;T2`; `a=5;b=1,1,1,0;T1`; `a=5;b=3,1,1,0;T1`

- `a=6`: `a=6;b=1,0,0,0;T1`; `a=6;b=1,0,0,0;T2`; `a=6;b=1,1,0,0;T1`; `a=6;b=1,1,0,0;T2`; `a=6;b=1,1,1,0;T1`; `a=6;b=1,1,1,0;T2`; `a=6;b=1,1,1,1;T1`; `a=6;b=3,0,0,0;T1`; `a=6;b=3,1,0,0;T1`

- `a=7`: `a=7;b=1,0,0,0;T1`; `a=7;b=1,0,0,0;T2`; `a=7;b=1,1,0,0;T1`; `a=7;b=1,1,0,0;T2`; `a=7;b=1,1,1,0;T1`; `a=7;b=1,1,1,0;T2`; `a=7;b=3,0,0,0;T1`; `a=7;b=3,0,0,0;T2`

- `a=8`: `a=8;b=0,0,0,0;T1`; `a=8;b=0,0,0,0;T2`; `a=8;b=1,0,0,0;T1`; `a=8;b=1,0,0,0;T2`; `a=8;b=1,1,0,0;T1`; `a=8;b=1,1,0,0;T2`

- `a=9`: `a=9;b=0,0,0,0;T1`; `a=9;b=1,0,0,0;T1`; `a=9;b=1,0,0,0;T2`

- `a=10`: `a=10;b=0,0,0,0;T1`

### Full disagreement list

None. Consequently there are no disagreement traces or counter-witnesses to report.

## Runtime and reproduction

Clean full-table run on the stated native Windows Python 3/SymPy environment: **115.859 seconds** wall-clock as measured inside the verifier. An earlier quiet run gave 124.562 seconds.

Run from `d2_plane_72_108` with:

```powershell
python audit_cascade_kills.py
```

The default output contains the homogeneity result, terminal count, a verdict row for each of the 420 open branches, the agreement summary, any disagreements with conservative witnesses, and runtime. `--quiet` suppresses only the 420-row table.
