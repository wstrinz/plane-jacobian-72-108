# Fourth corner: t=7 confirms the corner-signature fit; kappa = t-2 is structural

## Verdict

Two new out-of-sample points, both derived by the exact corner-144 template and both
**matching the six-parameter fit of PHI_75_125.md with zero fitting freedom**:

```
(56,84)  F_9 j=0, corner (7,21) -> (11\7,2), (a,b)=(2,3), t=7:
    Phi = -(1/10) y^107 (y^5+1)^54
    signature (deg, ord_y, mult_(y+1), cofactor) = (377, 107, 54, 216)   MATCHES

(50,75)  F_2 j=0, corner (5,20) -> (7\5,2),  (a,b)=(2,3), t=5:
    Phi = -(1/6) y^75 (y^3+1)^38
    signature = (189, 75, 38, 76)                                        MATCHES
```

Every parameter (`a, b, t, kappa, a0, q`) is read off the corner data, so these are
genuine predictions, not fits. The fit now reproduces **five** cases exactly, at
`t in {4,5,7}`:

| case | corner | `(a,b)` | `t` | `(a0,q,e,N)` | signature |
|---|---|---|---|---|---|
| `(72,108)` | `(8,28)` | `(2,3)` | `4` | `(8,7,2,28)` | `(238,204,30,4)` r=0 * |
| `(108,144)` | `(8,28)` | `(3,4)` | `4` | `(8,3,2,67)` | `(550,205,69,276)` |
| `(75,125)` | `(5,20)` | `(3,5)` | `5` | `(5,2,3,98)` | `(504,201,101,202)` |
| **`(56,84)`** | **`(7,21)`** | **`(2,3)`** | **`7`** | **`(7,2,2,52)`** | **`(377,107,54,216)`** |
| **`(50,75)`** | **`(5,20)`** | **`(2,3)`** | **`5`** | **`(5,2,2,36)`** | **`(189,75,38,76)`** |

The pair `(56,84)`/`(50,75)`/`(72,108)` is particularly clean: **the same reduced pair
`(a,b)=(2,3)` on three different corners** — the signature moves exactly as the corner
data says it should. Conversely `(50,75)` vs `(75,125)` is the same corner with
different `(a,b)`. Corner-dependence and pair-dependence are now each isolated.

All checked by `phi_corner4_verify.py` — **40/40 pass** — with independent
routines (branch-complete nonlinear solve for the forced `g`, full linear solve for
uniqueness of `f`, `factor_list`-based signature extraction, the operator-bracket
identity re-derived at all four `(t,kappa)` tuples). Derivation: `phi_corner4.py`.

---

## kappa = t-2: FORCED on the standard-chart class (the correlation is structural)

The directive asked for a corner with `kappa != t-2`. Answer: **there is none in the
class every published reduction lives in, and that is a theorem, not bad luck.**

The final reduction chart is the Laurent substitution `(X,Y) -> (x^-1, x^l y)`. Its
Jacobian is

```
d(x^-1) ^ d(x^l y) = -x^(l-2) dx ^ dy        (symbolically verified for generic l)
```

so `[P,Q] = 1` in original coordinates becomes `[p,q] = -x^(l-2)` in the chart —
`kappa = l-2` with no reference to the polygon. And `t = l` because each linear factor
`(Y - r X^-l)` of `ell(C)` pulls back to `x^l (y - r)`. Hence **`kappa = t-2`
identically for every corner reduced by this chart**: all 15 length-1 families anchored
at `A_0' = (1,0)` in the GGV5 `v11 <= 35` tables (F1–F11, F14–F17), plus the `(8,28)`
corner of GGHV22 (whose `-x^2` Jacobian at `l=4` is the same computation, 2204 lines
1228–1230).

Consequence, stated precisely: within this class `(t, kappa)` is **one parameter, not
two**, and the fit should be written with `kappa` eliminated:

```
N = a[t(a+b) - (kappa+1)] - 2b   --(kappa=t-2)-->   N = a[t(a+b-1) + 1] - 2b .
```

The reduced formula reproduces `N = 28, 67, 98, 52, 36` at all five points (checked).
The five-parameter fit `(a, b, t, a0, q)` is now tested at three distinct `t` values
with two out-of-sample exact hits; it is overdetermined and validated, not
underdetermined as it was at three points.

### Where `kappa != t-2` could still live

Only two escapes exist in the `v11 <= 35` tables, and both sit outside every published
reduction:

1. **`A_0' = (2,0)` families (F12, F13).** The first polygon edge ends at height-2, so
   the final chart need not be the pure Laurent substitution. No paper reduces these.
2. **Length-2 chains (F18–F24).** The composed transformation stacks two corner charts;
   composing `x_prev = x_new^-1` through the earlier Jacobian factor `-x_prev^(l1-2)`
   suggests a total exponent like `l2 - l1` rather than `l2 - 2` (e.g. F22 would give
   `kappa = 0` vs `t-2 = 2`) — but this composition heuristic is **unverified**; the
   intermediate root shifts and the meaning of the second denominator need the papers'
   full machinery. This is the named follow-up if separating `t` from `kappa` ever
   becomes load-bearing; with `kappa = t-2` proven on the standard class, it currently
   is not.

---

## Candidate survey (smallest coprime j per family)

From `phi_corner4.py` (Diophantine identity `(m+n)qk - n(ql-p) = k` checked for all):

| fam | A0 | A1 | (m,n) | degrees | t | a0 | q | e | r | N | gap | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| F1 | (4,12) | (7\4,3) | (3,4) | (48,64) | 4 | 4 | 3 | 2 | 0 | 67 | 1 | r=0 class |
| F2 | (5,20) | (7\5,2) | (2,3) | (50,75) | 5 | 5 | 2 | 2 | 2 | 36 | 0 | **derived (5th pt)** |
| F3 | (5,20) | (8\5,3) | (3,2) | (75,50) | 5 | 5 | 3 | 2 | 1 | 36 | 1 | |
| F4 | (5,20) | (8\5,3) | (3,16) | (75,400) | 5 | 5 | 3 | 14 | 1 | 241 | 1 | k=2 |
| F5 | (5,20) | (9\5,4) | (9,5) | (225,125) | 5 | 5 | 4 | 5 | 0 | 312 | 2 | r=0 class |
| F6 | (5,20) | (9\5,4) | (7,18) | (175,450) | 5 | 5 | 4 | 12 | 0 | 811 | 2 | k=2 |
| F7 | (6,15) | (7\3,4) | (2,7) | (42,147) | 3 | 6 | 4 | 6 | 1 | 36 | 1 | t=3 |
| F8 | (6,15) | (8\3,5) | (3,7) | (63,147) | 3 | 6 | 5 | 5 | 0 | 70 | 2 | t=3 |
| **F9** | **(7,21)** | **(11\7,2)** | **(2,3)** | **(56,84)** | **7** | **7** | **2** | **2** | **4** | **52** | **0** | **derived (4th pt)** |
| F10 | (7,21) | (13\7,3) | (7,4) | (196,112) | 7 | 7 | 3 | 4 | 3 | 270 | 1 | t=7 |
| F11 | (7,21) | (13\7,3) | (2,5) | (56,140) | 7 | 7 | 3 | 4 | 3 | 76 | 1 | k=2 |
| F12 | (8,24) | (13\4,5) | (3,7) | (96,224) | ? | 8 | 5 | 5 | 2 | ? | — | A0'=(2,0) |
| F13 | (9,21) | (13\3,7) | (2,13) | (60,390) | ? | 9 | 7 | 12 | 1 | ? | — | A0'=(2,0) |
| F14 | (9,24) | (7\3,4) | (2,7) | (66,231) | 3 | 9 | 4 | 6 | 4 | 36 | 0 | t=3, gap 0 |
| F15 | (9,24) | (8\3,5) | (3,7) | (99,231) | 3 | 9 | 5 | 5 | 3 | 70 | 1 | t=3 |
| F16 | (9,24) | (10\3,7) | (3,5) | (99,165) | 3 | 9 | 7 | 3 | 1 | 56 | 3 | t=3 |
| F17 | (9,24) | (11\3,8) | (2,3) | (66,99) | 3 | 9 | 8 | 2 | 0 | 20 | 4 | r=0 class |

`gap` = resonant degree − pure-ansatz degree; `gap > 0` marks the `(72,108)`-style
extra-cofactor regime (F1, F5, F17, ... — a ready-made testbed for the resonance-gap
story, untouched here). F9 was chosen as primary: new `t=7`, `k=1`, adjacent pair,
`gap=0`, and `(b-1)/a = 1` integer so the corner-144 slice-index picture transfers
verbatim (an assumption `(75,125)` could NOT satisfy). **F14 (t=3, gap=0, k=1) is the
natural sixth point if more coverage is ever wanted.**

## The derivations (both one-line collapses)

F9: ODE `14 c f' - 26 c' f = c^2`, `c = y^2 g`, `deg g = 5`. Ansatz `f = A y^3 g^2`
collapses it; the coefficient system has a single valid branch forcing `g_1..g_4 = 0`,
top coefficient resonant, `g(-1)=0` + monic gives `g = y^5+1 = (y+1)(y^4-y^3+y^2-y+1)`
(the same 10th-cyclotomic residual as `(108,144)`), `A = -1/10`, `f = -(1/10) y^3
(y^5+1)^2` (deg 13 = resonant, unique polynomial solution of degree <= 13). Then
`N = 52` and `Phi = -(1/10) y^107 (y^5+1)^54`.

F2 j=0: `10 c f' - 18 c' f = c^2` gives `g = y^3+1`, `A = -1/6`,
`f = -(1/6) y^3 (y^3+1)^2` (deg 9, unique), `N = 36`,
`Phi = -(1/6) y^75 (y^3+1)^38`. Residual `H_2 = y^2-y+1`, same as `(75,125)` — the
residual is a **corner** invariant (`g = y^(a0-q)+1` throughout), independent of `(a,b)`.

## `[judgment]` list — where this is conditional

1. **`[judgment: chain data]`** Corner rows read from the GGV5 `v11 <= 35` tables
   (paper_src/1708.07936_GGV5.tex, final section) and re-checked by the Diophantine
   identity for all 17 length-1 families. Primary-source.
2. **`[judgment: unreduced polygon]`** As with `(75,125)`: the `(7,21)` corner's
   explicit reduction is performed in no paper; the derivation assumes the standard
   type-II.b root shift + Laurent chart, giving `t = l = 7`, `kappa = l-2 = 5`,
   `deg C = a0 = 7`, `q = 2`. Same conditional boundary as CORNER_144 / PHI_75_125.
3. **`[judgment: N formula]`** `N = a[t(a+b)-(kappa+1)] - 2b` extrapolated to
   `t=7`. Here `(b-1)/a = 1` is an integer, so the corner-144 forcing-slice index
   transfers verbatim — this point is *less* conditional than `(75,125)` (whose
   `(b-1)/a = 4/3` broke the slice picture). Note the two new matches are evidence for
   the formula, not independent of it: a from-scratch C-series build at one new corner
   remains the way to discharge this item entirely.
4. **`[judgment: kappa = t-2 scope]`** The forcing argument covers exactly the
   single-Laurent-chart class. F12/F13 and the length-2 chains are excluded honestly;
   the `kappa = l2 - l1` composition heuristic there is conjecture, not claim.
