# Envelope bounds for the (72,108) D-transformation

This is the referee-facing proof of the coefficient windows used by the
(72,108) computation. Every finite calculation below is checked exactly by
envelope_bounds_verify.py; it uses integer/rational SymPy algebra only and
exits nonzero on failure.

## Theorem

Let

\[
C=x^4C_4+x^3C_3+\cdots,\quad C^2=P,\quad C_4=y^7(y+1),\qquad
D_r=C_rC_4^{7-2r}\in K[y],\quad w=4-r.
\]

For the polygons in GGHV22 Proposition “Case (8,28)” (lines 1000–1007), the
coefficients before and after x -> x-D_3/4 satisfy

\[
\begin{array}{c|cc}
 & \operatorname{ord}_yD_{4-w} & \deg_yD_{4-w}\\ \hline
\text{subcase (2)} & \ge12w & \le14w\\
\text{subcase (1)} & \ge12w & \le15w.
\end{array}
\]

The checker proves this for w=0,...,5. The w=5 step is needed because the
active resultant variables d2,d1,d0,d_-1 have weights 2,3,4,5.

After removing y^(12w), the degree caps are 2w in subcase (2) and 3w
in subcase (1): (4,6,8,10) and (6,9,12,15) at those four weights.

## Source and base data

Checks S1–S12 compare inputs directly with repository sources. GGHV22 lines
1000–1007 give [P,Q]=x^2 and both polygon lists; S2–S6 compare every corner
with T3_WINDOW_AUDIT.md lines 19–25. GGHV22 lines 1412–1414, 1428–1436,
and 1455–1467 give the published normalization and induction template
(S8–S10). S11–S12 check the t=4 data in STATE.md items 1–4 and
T3_WINDOW_AUDIT.md lines 52–55.

B1–B6 verify that R=x^4y^7(y+1) squared and cubed have leading corners
(8,14),(8,16) and (12,21),(12,24). B7–B10 compute

\[
\max_{N(P)_1}v_{-1,1}=8,\qquad
\max_{N(P)_2}v_{-2,1}=0,\qquad
\max_{N(P)_i}v_{2,-1}=2.
\]

The first is attained at (8,16),(0,8), the second at (0,0),(8,16), and
the third on the edge from (1,0) to (8,14). B11–B12 check that C4^-1
contributes -8 at infinity and +7 at zero.

## Recursion and valuation induction

Coefficient extraction from C^2=P gives

\[
C_{4-w}=\frac{P_{8-w}-\sum_{j=1}^{w-1}C_{4-j}C_{4-w+j}}{2C_4}.
\]

R1–R6 substitute this expression into the coefficient of x^(8-w) for
w=0,...,5. The three exact affine inductions I1.*, I2.*, I3.* give

\[
v_{-1,1}(C_{4-w})\le8-w,\quad
v_{-2,1}(C_{4-w})\le8-2w,\quad
v_{2,-1}(C_{4-w})\le2w-7.
\]

At each non-base step the polygon bound for P_(8-w) equals every product
bound, and the valuation of C4^-1 closes the induction. With r=4-w,

\[
\deg C_r\le r+4\ (\text{sub1}),\quad \deg C_r\le2r\ (\text{sub2}),
\quad \operatorname{ord}_yC_r\ge2r-1.
\]

## D-transform, magic weights, and shift

M8.* checks exponent cancellation in

\[
D_r=\tfrac12P_{r+4}C_4^{6-2r}-\tfrac12\sum D_iD_j\in K[y]
\]

downward from D4=1; each power on P_(r+4) is nonnegative. M1–M3 solve
the condition that the coefficient of r disappear from total valuation.
The unique solutions are

\[
v_{15,1}\ (\text{sub1}),\qquad v_{14,1}\ (\text{sub2}),\qquad
v_{-12,-1}\ (\text{order}).
\]

M4–M6 reduce the totals to 60, 56, and -48; M7.0–M7.5 rearrange them
to 15w, 14w, and 12w.

The coefficient of x^r in (x-a)^i is binomial(i,i-r)(-a)^(i-r).
T1 checks that a=D3/4 kills the x^3 term. Every contribution to the new
D_r has weight (i-r)+(4-i)=4-r; T2.* and T3.* verify preservation of all
bounds through r=-1,...,4.

## Downstream agreement

D1–D10 compare the theorem with consumers without importing them:

- regenerate_system.py assigns weights (2,3,4,5) (D3).
- jetlift.py has sizes (5,7,9,11) / (7,10,13,16) and slice counts
  (251,269) / (376,403) (D4–D5).
- sub1_cascade_verify.py consumes (6,9,12,15) and cascade_engine.py
  consumes subcase-(2) caps and both slopes (D6–D8).
- Source factors have exact weights 125 and 134 (D9), and weight-17 Phi
  attains deg=238=14*17, ord=204=12*17 (D10).

regenerate_system.py constructs weighted equations but contains no numeric
windows; those are consumed by the other files checked above.

## Mechanization findings

1. **F1_RECURSION_SIGN_TYPO.** GGHV22 lines 1462–1466 and
   T3_WINDOW_AUDIT.md line 55 print -(P+sum)/(2C4). This does not satisfy
   C^2=P; R7 detects its generic residual. The correct identity is above.
   The envelope proof survives because its valuation bounds are sign-insensitive.

2. **F2_WEIGHT_RANGE_STALENESS.** w=0,...,4 is not all the active program
   uses: d_-1 has weight 5. The checker covers w=0,...,5.

No envelope inequality failed mechanization.

## Named assumptions

1. **A1_PUBLISHED_POLYGON_REDUCTION.** GGHV22 Proposition “Case (8,28)”
   (lines 1000–1007) is used as stated. The checker verifies its transcription,
   not the proposition's proof.

2. **A2_COMMON_ROOT_AND_NORMALIZATION.** Cited GGV1 Propositions 1.13 and
   2.1, with GGHV22's scalar/linear normalization at lines 1411–1414, give
   ell(P)=R^2, ell(Q)=R^3, and R=x^4*y^7*(y+1). The checker verifies
   consequences, not those propositions or the WLOG gauge.

3. **A3_LAURENT_VALUATION_FRAMEWORK.** The field has characteristic zero;
   the embeddings in K((y^-1)) and K((y)) have standard multiplicative
   non-Archimedean valuations; and x -> x-D3/4 is valid formally.
