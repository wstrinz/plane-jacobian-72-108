# T5 T2 survivor column: four cells infeasible, eight cells open

**Date:** 2026-07-22
**Branch:** `f31`, subcase (2), `d1=0`, `sigma=4d0-d2^2 != 0`.

This applies the level-5 squeeze to every survivor case in the twelve requested
split-place cells.  The result is deliberately narrower than the claim in
`PHASE_C_WORKLIST.md` section 4:

```text
PROVEN infeasible cells: a5 b1000;
                         a6 b1000, b1100, b1110.       (4 cells, 8 cases)
OPEN cells:              all four a7 cells, all three a8 cells,
                         a9 b1000.                     (8 cells, 24 cases)
```

All computational assertions below are checked by `t5_t2_column_verify.py`.
Its checks C1--C5 parse `f31_graded.txt` and `cascade_cones_qt.json`; no `h_l`
coefficient is entered independently.

## 0. Case notation and audited data

Over a splitting field write `q` (up to a unit) as `p1*p2*p3*p4`.  Put

```text
e = t^a R F,       R = product p_i^b_i,       gcd(F,tq)=1,
sigma = P Z,       P = product p_i^s_i,       gcd(Z,q)=1,
g6 = Q G,          Q = product p_i^m_i,       gcd(G,q)=1.
```

Here `s_i=v_i(sigma)` and `m_i=v_i(g6)` are the exact witness values.  Set
`B=sum b_i`, `S=sum s_i`, `M=sum m_i`, `f=deg F`, `z=deg Z`, and `g=deg G`.
The JSON gives the following common q-profile for every flag case in a cell.
The `F cap` column is the cap before the level-5 squeeze.

| `a` | `b` | `v=30-3a` | `s` | `m` | `B` | `F cap=10-a-B` | `G cap=10+3a-M` |
|---:|:---:|---:|:---:|:---:|---:|---:|---:|
| 9 | 1000 | 3 | 2000 | 7666 | 1 | 0 | 12 |
| 8 | 0000 | 6 | 0000 | 6666 | 0 | 2 | 10 |
| 8 | 1000 | 6 | 2000 | 7666 | 1 | 1 | 9 |
| 8 | 1100 | 6 | 2200 | 7766 | 2 | 0 | 8 |
| 7 | 1000 | 9 | 2000 | 7666 | 1 | 2 | 6 |
| 7 | 1100 | 9 | 2200 | 7766 | 2 | 1 | 5 |
| 7 | 1110 | 9 | 2220 | 7776 | 3 | 0 | 4 |
| 7 | 3000 | 9 | 7000 | `(11)666` | 3 | 0 | 2 |
| 6 | 1000 | 12 | 2000 | 7666 | 1 | 3 | 3 |
| 6 | 1100 | 12 | 2200 | 7766 | 2 | 2 | 2 |
| 6 | 1110 | 12 | 2220 | 7776 | 3 | 1 | 1 |
| 5 | 1000 | 15 | 2000 | 7666 | 1 | 4 | 0 |

For example, `7666` means `(7,6,6,6)`.  Every entry satisfies the terminal
law `3b_i+m_i=6+2s_i` [C3].  All 32 cases have `sigma != 0`.

Flag labels used below are:

* `N`: `d2 != 0`, no globally zero `g_l`;
* `G4`: `d2 != 0`, `g4=0`;
* `G5`: `d2 != 0`, `g5=0` (only `a9 b1000`);
* `D`: `d2=0`, no globally zero `g_l`.

The exact extra equations are

```text
G4:  t^v g5 = c^4 q^4 h4,
G5:  t^v g6 = c^5 q^5 h5,   F^3 R^3 g4 = -c^4 q^4 h4,
D:   h5 = 2048 e^2,          h4 = -12096 sigma^3.
```

They follow from the cascade and the source-parsed collapses [C1].

## 1. Split-support level-5 squeeze

The terminal equation is

```text
F^3 G = 3072 c^6 Z^2                                      (1)
```

after the prescribed `p_i` factors are cancelled; the split-root
normalizations only alter the displayed nonzero unit.  Since every `m_i>=1`,
`Q/q` is a polynomial.  Substituting the terminal relation into

```text
t^v g6 = (RF)^3 g5 + c^5 q^5(-9216 d2 sigma^2+2048 e^2)
```

gives the exact rearrangement [C2]

```text
t^v QG - 2048 c^5 q^5 t^(2a) R^2 F^2
    = R^3 F^3 (g5 + 19890 d2 (Q/q)G).                     (2)
```

This is the part of `T5_STRATA_50_11.md` that transfers.  Because
`gcd(F,tQ)=1`, reducing (2) modulo `F^2` gives

```text
F^2 | G.                                                   (3)
```

The modification is essential: it is the q-coprime remainder `F`, not the
whole split-supported factor `RF`, whose square divides `G`.  Thus the
uniform-quartic remainder bound quoted in the worklist does not apply.

Taking degrees in (1) and using (3) leaves exactly

```text
3f+g=2z,    2f<=g<=G cap,    f<=F cap,    z<=8-S.          (4)
```

All state lists below are the complete integer solutions of (4) [C4].

## 2. Infinity test and the four killed cells

Put `D=deg e=a+B+f` and `Sigma=deg sigma=S+z`, and substitute
`d0=(d2^2+sigma)/4`.  Let

```text
T_j = Phi^j e^(21-3j) H_j,   H_j=h_j(d2,0,(d2^2+sigma)/4,e).
```

The verifier parses every monomial of every `H_j` to bound `deg T_j` for
each possible `deg d2` in `{0,1,2,3,4}`, and separately for `d2=0` [C4].
The top term is exact because

```text
H6=-3072 sigma^2,   deg T6=204+3D+2Sigma.                (5)
```

For the four cells below, every state from (4) has `deg T6` strictly larger
than the source-derived caps for `T0,...,T5`:

| cell | `(f,z)` | `(D,Sigma)` | `max deg(T0..T5)` | `deg T6` |
|:--|:--|:--|--:|--:|
| `a5 b1000` | `(0,0)` | `(6,2)` | 218 | 226 |
| `a6 b1000` | `(0,0)` | `(7,2)` | 226 | 229 |
|  | `(0,1)` | `(7,3)` | 226 | 231 |
|  | `(1,3)` | `(8,5)` | 234 | 238 |
| `a6 b1100` | `(0,0)` | `(8,4)` | 234 | 236 |
|  | `(0,1)` | `(8,5)` | 234 | 238 |
| `a6 b1110` | `(0,0)` | `(9,6)` | 242 | 243 |

These are strict degree contradictions [C5].  They do not depend on whether
`g4` vanishes, so every authoritative flag case in these cells is killed:

| cell | case `N` | case `G4` | cell verdict |
|:--|:--|:--|:--|
| `a5 b1000` | **PROVEN infeasible** | **PROVEN infeasible** | **PROVEN infeasible** |
| `a6 b1000` | **PROVEN infeasible** | **PROVEN infeasible** | **PROVEN infeasible** |
| `a6 b1100` | **PROVEN infeasible** | **PROVEN infeasible** | **PROVEN infeasible** |
| `a6 b1110` | **PROVEN infeasible** | **PROVEN infeasible** | **PROVEN infeasible** |

This accounts for all 8 survivor cases in the four killed cells, not merely
one representative case.

## 3. Why the remaining infinity states do not close

Every residual below is one of two tied-leading-degree patterns:

* **A:** `(D,Sigma)=(8,3)`.  `T5` and `T6` both have degree 234.
* **B:** `D=10`, `2<=Sigma<=8`.  `T0` and `T5` share degree 250;
  when `deg d2=4`, further intermediate caps can also be 250, and when
  `Sigma=8`, `T6` is also degree 250.  If `d2=0`, the `T0`--`T5` tie remains.

Thus neither the Proposition-E constant-sigma argument nor the
`T5_90_T2.md` unique-top argument transfers to these states.  A tie is not a
contradiction: its leading-coefficient equation has not been refuted.

The t-place witnesses were also rechecked against
`v+s_(l+1)=ultrametric(s_l,w_l)` [C3].  They give no extra tropical kill:

| flag | t-place profile `(v_t(d2);s4,s5,s6)` | unresolved t condition |
|:--|:--|:--|
| `N` | `(0;0,0,0)` | level-5 and level-4 term cancellations, each depth `v` |
| `G4` | `(0;inf,0,0)` | level-5 term cancellation and level-4 monomial rise, depth `v` |
| `D` | `(inf;0,v,0)` | level-4 term cancellation, depth `2v` |
| `G5` (`a=9`) | `(3;0,inf,0)` | level-4 exact identity |

Accordingly the tables in the next sections state the surviving global
degree/factor problem together with the applicable flag equation from section
0.  They do not promote a residue obligation to a proof of infeasibility.

## 4. Priority cell: `a9 b1000`

Here `F` is already constant, `D=10`, and the residual states before using a
zero flag are

```text
(f,z,g;D,Sigma) = (0,z,2z;10,2+z),   z=0,1,2,3,4,5,6.    (R9)
```

All are pattern B.  The four JSON cases have the following individual
verdicts.

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | (1)--(4), (2), and one of the seven states `(R9)`; the degree-250 leading equation is unresolved. |
| `G5` | **open** | The exact line `t^3 g6=c^5q^5h5` eliminates `z<6`: its left degree is `28+2z`, while for `Sigma<8`, `deg h5=20` from `2048e^2`, so the right degree is 40.  Residual: `(f,z,g;D,Sigma)=(0,6,12;10,8)`, `deg h5=20` (no leading cancellation), plus `R^3g4=-c^4q^4h4`. |
| `G4` | **open** | One of `(R9)`, plus `t^3g5=c^4q^4h4`; the degree-250 tie remains. |
| `D` | **open** | One of `(R9)` with `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`; `T0` and `T5` still tie at degree 250. |

The `G5` narrowing, including `v_t(d2)=3`, is checked directly from that
case's witness and the parsed `h5` [C3,C4].  It does not kill the final
`Sigma=8` state.

## 5. The `a=8` group

### 5.1 Cell `a8 b0000`

After removing every uniquely dominated state, the complete residual is

```text
R80 = {(0,3,6;8,3),
       (2,5,4;10,5), (2,6,6;10,6),
       (2,7,8;10,7), (2,8,10;10,8)}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One state in `R80`; the first is pattern A and the other four are pattern B. |
| `G4` | **open** | One state in `R80`, together with `t^6g5=c^4q^4h4`. |
| `D` | **open** | One state in `R80`, with `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`. |

### 5.2 Cell `a8 b1000`

The squeeze and infinity pruning force `f=1`; the residual is

```text
R81 = {(1,3,3;10,5), (1,4,5;10,6),
       (1,5,7;10,7), (1,6,9;10,8)}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One of the four pattern-B states in `R81`. |
| `G4` | **open** | One state in `R81`, plus `t^6g5=c^4q^4h4`. |
| `D` | **open** | One state in `R81`, plus `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`. |

### 5.3 Cell `a8 b1100`

Here `F` is constant.  The residual is

```text
R82 = {(0,z,2z;10,4+z): z=0,1,2,3,4}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One of the five pattern-B states in `R82`. |
| `G4` | **open** | One state in `R82`, plus `t^6g5=c^4q^4h4`. |
| `D` | **open** | One state in `R82`, plus `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`. |

## 6. The `a=7` group

### 6.1 Cell `a7 b1000`

```text
R71 = {(0,1,2;8,3), (2,5,4;10,7), (2,6,6;10,8)}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One state in `R71`; the first is pattern A and the last two are pattern B. |
| `G4` | **open** | One state in `R71`, plus `t^9g5=c^4q^4h4`. |
| `D` | **open** | One state in `R71`, plus `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`. |

### 6.2 Cell `a7 b1100`

```text
R72 = {(1,3,3;10,7), (1,4,5;10,8)}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One of the two pattern-B states in `R72`. |
| `G4` | **open** | One state in `R72`, plus `t^9g5=c^4q^4h4`. |
| `D` | **open** | One state in `R72`, plus `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`. |

### 6.3 Cell `a7 b1110`

```text
R73 = {(0,0,0;10,6), (0,1,2;10,7), (0,2,4;10,8)}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One of the three pattern-B states in `R73`. |
| `G4` | **open** | One state in `R73`, plus `t^9g5=c^4q^4h4`. |
| `D` | **open** | One state in `R73`, plus `d2=0`, `h5=2048e^2`, `h4=-12096sigma^3`. |

### 6.4 Cell `a7 b3000`

```text
R74 = {(0,0,0;10,7), (0,1,2;10,8)}.
```

| case | verdict | exact residual statement |
|:--|:--|:--|
| `N` | **open** | One of the two pattern-B states in `R74`. |
| `G4` | **open** | One state in `R74`, plus `t^9g5=c^4q^4h4`. |

There is no `d2=0` survivor case in this cell's authoritative JSON.

## 7. Conclusions

The transferred and modified steps are:

1. **Transferred exactly:** terminal absorption of the `sigma^2` term and the
   constant `-9216c^5/(3072c^6)=19890`.
2. **Modified for partial q-support:** factor the prescribed split-root parts
   first.  Coprimality proves `F^2|G`, not `(RF)^2|G`; the cofactor cap uses
   the exact sum `M=sum_i(6+2s_i-3b_i)`.
3. **Transferred with a wider degree enumeration:** source-derived infinity
   domination kills all states in `a=5` and `a=6`.
4. **Does not transfer:** constant residual `F` does not imply constant
   `sigma` in the `a>=7` cells.  Their surviving degrees land on patterns A
   or B, where the infinity maximum is tied.
5. **Flag-sensitive refinement:** `g5=0` in `a9 b1000` forces the single
   degree state `(0,6,12;10,8)`, but does not contradict it.  `g4=0` and
   `d2=0` add the exact equations recorded above; no unsupported cancellation
   claim is made.

Final ledger:

```text
cells proven infeasible: 4 / 12       survivor cases killed:  8 / 32
cells open:               8 / 12       survivor cases open:   24 / 32
```

## 8. Verification map

Run from this directory:

```text
python t5_t2_column_verify.py
```

| check | content |
|:--|:--|
| C1 | Parses `f31_graded.txt`; checks the exact `h4`, `h5`, `h6`, `h7` collapses and `deg Phi=34`. |
| C2 | Checks the generalized split-support level-5 rearrangement and the coefficient 19890 symbolically over `Q`. |
| C3 | Reads all 12 cells and all 32 cases from `cascade_cones_qt.json`; checks flag inventory, every q-terminal law, and both t-place coupling lines. |
| C4 | Enumerates every solution of (4), parses all `H_j` monomials for infinity caps, checks every killed/open verdict, and checks the `G5` narrowing. |
| C5 | Checks and prints the seven strict `max(T0..T5)<deg(T6)` margins used for the four cell kills. |
