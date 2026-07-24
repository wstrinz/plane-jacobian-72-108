# The family grammar: a finite grammar of local obstruction types

Generalizes the F2-family symbolic closed form (`f2_family_verify.py`,
STATE.md "THE F2 FAMILY THEOREM") to **every** family of the GGV5 `v11 <= 35`
survey. For each of the 17 length-1 families (+ the 3 length-2 escapes) we
attempt the F2 trick — the ansatz `f = A y^rho g^e` driven by the collapse
identity — symbolically in the family parameter `j`, and classify the outcome.

Derivation: `family_grammar.py`. Exact checker: `family_grammar_verify.py`
(**210/210, --quiet, exit 0**); every closed form re-verified symbolically in
`j` AND by direct full-ODE substitution at `j = 0,1,2,3`, every landed derived
point reproduced exactly. Nothing existing edited; `run_tests.sh` untouched.

## RESCOPE (2026-07-24): CANONICAL FORCING BRANCHES, not exhaustive classification

> **What this theorem establishes is that every family has a CANONICAL MODELED
> BRANCH** — a distinguished forcing-ODE solution branch that falls into one of
> the three grammars (PURE / COFACTOR / RUNG-STRUCTURED), selected by the
> corner-data dichotomy on `(gap, r)` and the root-shift gauge. It does **NOT**
> establish an **exhaustive classification of every solution branch** of every
> family's forcing system. A genuine exhaustive statement would require primary
> decomposition of the forcing ideal and root-partition completeness of the
> residual divisor (which repeated-root / coexisting-multiplicity branches can
> violate — cf. the `mu=1,2,3` coexistence at `dg=3`, POLYGON_REDUCTION.md R3
> reopening, and the `dg`-even complex-branch caveat in §3's complex-scope
> discipline). Those are **not** ruled out here; the canonical branch is the
> representative we derive and check, not the only branch.
>
> **Read every "IRREGULAR 0" / "zero irregular" statement below as "zero families
> LACKING a canonical modeled branch"** — i.e. the grammar covers every family
> with a modeled branch, not "every family's every branch is one of these."

## 1. The one theorem behind the whole grammar

Corner data `(t, kappa=t-2, a0, q)` is **fixed** within a family; only the pair
`(a,b)(j)` moves, so `e = b-a+1` and `rho = (e-1)q+1` and `N` vary with `j` while
`dg = a0-q`, `r = a0-q-1`, `gap = (q-1) - a0/t` are **j-invariant**. Writing
`c = y^q g`, the forcing ODE is

```
a { t c f' - [t(b-a)+kappa+1] c' f } = c^(b-a+1).
```

**Collapse identity (all families, all dg):** `y (y^dg)' - dg (y^dg+1) = -dg`.

**CLOSED-FORM THEOREM (proved symbolically in j).** The pure F2 ansatz
`f = A y^rho (y^dg+1)^e` solves the ODE for *all j* **iff**

```
t - (kappa+1) q + dg = 0   <=>   a0 = t(q-1)   <=>   gap = 0,
```

and then `A = -1/(a·dg)`. The mechanism is two j-constant reductions plus the
collapse:

```
te - coef        = t - kappa - 1 = 1           (the kappa=t-2 identity in ODE form)
t*rho - coef*q   = t - (kappa+1)q              (constant in j; the "resonance" number)
```

so `a{...} = A a y^(qe) g^e {[t-(kappa+1)q + dg] g - dg}`, and the g-term dies
exactly at `gap = 0`, leaving `A a(-dg) = 1`. This is the F2 mechanism
(`y g' - 3 g = -3`, `A = -1/(3a)`) made general: F2 is the `dg=3, t=5, q=2`
instance. **`gap` is a pure corner-data invariant — no ODE needed to predict the
regime.**

The trichotomy that follows (with `r = dg-1` always):

| regime | condition | closed form in j |
|---|---|---|
| **PURE** | `gap = 0` | `f = -1/(a·dg) · y^rho (y^dg+1)^e` |
| **COFACTOR** | `r = 0` (`dg=1`), `gap>0` | `f = y^rho (y+1)^e · u`, `u` a UNIT of degree `gap` |
| **RUNG-STRUCTURED** | `r>0` and `gap>0` | unramified fails; `mu=dg` rung `f = y^rho (y+1)^(dg·e-(dg-1)) · u`, `deg u = gap+r` |
| **IRREGULAR** | `gap < 0` | none (res < pure); only among the composite escapes (F22) |

Both closed-form branches and the `mu=dg` ramified rung are **uniform in j** and
were re-derived symbolically; the RUNG family's `mu=dg` rung reproduces all four
PHI_F7 landed polynomials on the nose.

## 2. The grammar table

`gap = (q-1) - a0/t`; `dg = a0-q = r+1`; block size = `t`; `W_step` denominator
is the reduced denominator of `ord_y(Phi)/M`, `M = t(a+b)-(kappa+1)` (the F2 law
`5a-3` is the `q_window` analogue).

| fam | t | a0 | q | dg | r | gap | class | N_j (factored) | W_step denom | residual / Galois |
|---|---|---|---|---|---|---|---|---|---|---|
| **F2**  | 5 | 5 | 2 | 3 | 2 | 0 | PURE | `(3j+4)(5j+9)` | `5j+7` | `y^3+1`, H=Φ6 (C2, Q(√-3)) |
| **F9**  | 7 | 7 | 2 | 5 | 4 | 0 | PURE | `(3j+4)(7j+13)` | `21j+29` | `y^5+1`, H=Φ10 (C4) |
| **F14** | 3 | 9 | 4 | 5 | 4 | 0 | PURE | `(3j+4)(5j+9)` | `3j+5` | `y^5+1`, H=Φ10 (C4) |
| **F1**  | 4 | 4 | 3 | 1 | 0 | 1 | COFACTOR | `40j²+104j+67` | `4j+5` | `y+1` (rational) |
| **F5**  | 5 | 5 | 4 | 1 | 0 | 2 | COFACTOR | `220j²+525j+312` | `5j+6` | `y+1` (rational) |
| **F6**  | 5 | 5 | 4 | 1 | 0 | 2 | COFACTOR *(k=2)* | `165j²+402j+244` | `5j+6` | `y+1` (rational) |
| **F8**  | 3 | 6 | 5 | 1 | 0 | 2 | COFACTOR | `(6j+7)(7j+10)` | `3j+4` | `y+1` (rational) |
| **F17** | 3 | 9 | 8 | 1 | 0 | 4 | COFACTOR | `(13j+5)(15j+4)` | `3j+1` | `y+1` (rational) |
| **F3**  | 5 | 5 | 3 | 2 | 1 | 1 | RUNG (dg even) | `105j²+125j+36` | — | `(y+1)²` ramified; μ=1 real-empty |
| **F4**  | 5 | 5 | 3 | 2 | 1 | 1 | RUNG *(k=2)* | `140j²+368j+241` | — | `(y+1)²` ramified |
| **F7**  | 3 | 6 | 4 | 2 | 1 | 1 | RUNG (dg even) | `(3j+4)(5j+9)` | — | `(y+1)²` ramified; μ=1 real-empty |
| **F10** | 7 | 7 | 3 | 4 | 3 | 1 | RUNG (dg even) | `168j²+427j+270` | — | `(y+1)⁴` ramified; real μ∈{2,4} |
| **F11** | 7 | 7 | 3 | 4 | 3 | 1 | RUNG *(k=2)* | `28j²+93j+76` | — | `(y+1)⁴` ramified |
| **F15** | 3 | 9 | 5 | 4 | 3 | 1 | RUNG (dg even) | `(6j+7)(7j+10)` | — | `(y+1)⁴` ramified |
| **F16** | 3 | 9 | 7 | 2 | 1 | 3 | RUNG (dg even) | `(11j+8)(12j+7)` | — | `(y+1)²` ramified; μ=1 real-empty |
| **F12** | 4 | 8 | 5 | 3 | 2 | 2 | RUNG (dg ODD)* | `56j²+148j+97` | — | μ∈{1,2,3} realized (η=0) |
| **F13** | 3 | 9 | 7 | 2 | 1 | 3 | RUNG (dg even)* | `(3j+4)(8j+15)` | — | `(y+1)²` ramified |

`*` = conditional: F12/F13 are `A0'=(2,0)` (ζ-tail model, off-diagonal N —
ZETA_TAIL.md); F4/F6/F11 are `k=2` (N-formula unverified upstream). RUNG rows
have no clean `W_step` (the residual is fully ramified, not a graded unit).

**Census (length-1): CLOSED-FORM 8 (PURE 3 + COFACTOR 5), RUNG-STRUCTURED 9,
families LACKING a canonical modeled branch: 0.** (Per the 2026-07-24 rescope:
this counts canonical modeled branches, one per family — not an exhaustive
enumeration of all solution branches. The `gap<0` IRREGULAR cell is empty among
the length-1 families and appears only among the composite escapes, F22.)

## 3. Per-family closed forms

### PURE (gap = 0): F2, F9, F14
`f_j = -1/(a·dg) y^rho (y^dg+1)^e`, residual `H = Φ_{2dg}` (indexed by dg alone),
signature (μ=1) `(res+N a0, rho+N q, e+N, r(e+N))`. Block recurrence
`Phi_{j+1} = (a_j/a_{j+1}) y^(Δrho) (y^dg+1)^(Δe) c^(ΔN) Phi_j`; for F2 (Δe=1,
Δrho=q=2) this collapses to the STATE.md form `(a/(a+1)) C^(30a+3)`.

- **F2** `a=j+2, b=2j+3` (`e=a`): `f=-1/(3(j+2)) y^(2j+3)(y^3+1)^(j+2)`, `N=(3j+4)(5j+9)`.
  Landed: `j=0` (50,75) `(189,75,38,76)`; `j=1` (75,125) `(504,201,101,202)`. ✔
- **F9** `a=j+2, b=2j+3`: `f=-1/(5(j+2)) y^(2j+3)(y^5+1)^(j+2)`, `N=(3j+4)(7j+13)`.
  Landed: `j=0` (56,84) `(377,107,54,216)`. ✔
- **F14** `a=j+2, b=4j+7`: `f=-1/(5(j+2)) y^(3j+5)(y^5+1)^(3j+6)`, `N=(3j+4)(5j+9)`.
  Landed: `j=0` (66,231) `(375,165,42,168)`. ✔

### COFACTOR (r = 0, gap > 0): F1, F5, F6, F8, F17
`dg=1` is forced, `g=y+1`, and `f_j = y^rho (y+1)^e · u` with `u` a UNIT cofactor
of degree `gap`, coefficients rational in `j`, satisfying `u(0)≠0`, `u(-1)=-1/a`,
leading coeff `≠0`. Signature (μ=1) `(res+N a0, rho+N q, e+N, gap)` — the cofactor
degree **is** `gap`, the degree-`gap` analogue of the (72,108) quartic.

- **F1** `u=(4y-1)/(5(2j+3))`; `j=0` gives `(1/15)(4y-1)` → (48,64) `(275,205,69,1)`. ✔
- **F5** `u=(-25y²+5y-3)/(33(4j+5))`; **F8** `u=(-9y²+3y-2)/(14(2j+3))`;
  **F17** `deg u = 4` unit, all j-uniform (see `family_grammar.py`).

### RUNG-STRUCTURED (r > 0, gap > 0): F3,F4,F7,F10,F11,F12,F13,F15,F16
The unramified ansatz fails: for **dg even** `y^dg+1` has no root at `-1` (the
root-shift gauge needs `g(-1)=0`), forcing ramification; for **dg odd** (F12)
the residual is a squarefree branch variety, not `y^dg+1`. But the **`mu=dg`
fully-ramified rung** gives a uniform closed form in j:

```
f_j = y^rho (y+1)^(dg·e-(dg-1)) · u,     deg u = gap + r,   u coefficients rational in j.
```

This reproduces **all four PHI_F7 landed polynomials** exactly (checked by direct
ODE substitution):

- **F7** `j=0`: `(1/10) y^21 (y+1)^11 (9y²+3y-1)` → (42,147) `(250,165,83,2)`. ✔
- **F3** `j=0`: `(1/42) y^4 (y+1)^3 (25y²+15y-3)` → (75,50) `(189,112,75,2)`. ✔
- **F10** `j=0`: `(1/3740) y^10 (y+1)^13 (2401y⁴+5831y³+4165y²+595y-85)` → (196,112) `(1917,820,1093,4)`. ✔
- **F16** `j=0`: `(1/330) y^15 (y+1)^5 (243y⁴+81y³-27y²+15y-10)` → (99,165) `(528,407,117,4)`. ✔

The full branch structure is the **μ-graded law** (ZETA_TAIL.md / MU_RUNGS_F10.md),
which our checker confirms reproduces every published rung:

```
deg  = res + N·a0,                    ord  = rho + N·q,
mult = mu(e+N) - (mu-1),             cof  = gap + r(e+N) - (mu-1)(e+N-1).
```

Identically `cof = deg - ord - mult`, and (via `r=dg-1`) it specializes to the
unramified law at `mu=1` and PHI_F7's ramified law at `mu=dg`. Realized rungs:
F12 (η=0) `mu∈{1,2,3}` `(814,506,·,·)`; F10 real `mu∈{2,4}` `(1917,820,·,·)`.

**Complex-scope discipline (MU_RUNGS correction, adopted).** Branch schemes are
classified over Qbar. `dg` is **even** on every length-1 `k=1` `A0'=(1,0)` RUNG
row (F3,F7,F10,F15,F16); there the `mu=1` **real** locus is empty (PHI_F7's
complete factorization over C at `dg=2`; MU_RUNGS' real Sturm counts at `dg=4`),
but complex `mu=1` branches are **not excluded** at `dg≥4`. This is
branch-selection annotation — never a complex kill. The rational `mu=dg` rung is
the uniform representative we derive.

## 4. The escapes and the (72,108) exception

Length-2 composite escapes (charts derived in no paper — CONDITIONAL,
`composite_charts.py`):

| fam | (a,b)@j0 | t | a0 | q | dg | r | gap | verdict |
|---|---|---|---|---|---|---|---|---|
| **F22** | (2,3) | 4 | 8 | 2 | 6 | 5 | **−1** | **IRREGULAR** — res < pure, regime unobserved anywhere; the obstruction is a negative resonance gap (the ansatz's forced degree exceeds the tower window) |
| **F23** | (2,7) | 4 | 8 | 4 | 4 | 3 | 1 | RUNG (conditional) |
| **F24** | (3,4) | 8 | 8 | 3 | 5 | 4 | 1 | RUNG, dg ODD (conditional) |

**(72,108) is handled honestly as the resonance EXCEPTION, not a family member.**
It is the audited `(8,28)` GGHV corner (`a0=8,q=7,t=4,e=2,dg=1,r=0,gap=4`) — an
`r=0` COFACTOR-shaped point whose unit cofactor is the degree-4 quartic
`2048y⁴-512y³+320y²-240y+195`, i.e. `cof = gap = 4`. It obeys the same
`r=0`-amended law as F1 (`cof=gap`) but is its own corner; **no family formula in
this table generates it**, and none is claimed to. (Its neighbour `(108,144)`,
also `(8,28)`, is a gap=0 point `(550,205,69,276)` outside the 17 families.)

## 5. Surprises

1. **The whole grammar is a corner-data dichotomy on `(gap, r)`.** Whether the F2
   trick works in j is decided before any ODE is written: `gap=0` ⇒ pure,
   `r=0` ⇒ unit cofactor, both `>0` ⇒ ramified rungs, `gap<0` ⇒ irregular. The F2
   closed form was not special — it is the `gap=0` cell of a 2×2 table.
2. **`A = -1/(a·dg)` universally** for the pure cell — the `-1/(3a)` of F2 is just
   `dg=3`. And `u(-1) = -1/a` universally for the cofactor cell (a clean
   родовой invariant across F1/F5/F6/F8/F17).
3. **RUNG families are not formula-less.** The `mu=dg` ramified rung is a genuine
   uniform closed form in j (deg-`(gap+r)` unit cofactor with rational-in-j
   coefficients) that reproduces the PHI_F7 points — the "irregular" regime of
   the review is actually **rung-structured with a derivable representative**.
4. **`dg` even ⟺ RUNG on the standard rows**, `dg` odd ⟺ closed-form OR the F12
   exception — the parity mini-lemma of PHI_F7 falls out of `dg = a0-q` and the
   gauge `g(-1)=0`.

## 6. `[judgment]` list

1. **[chain data]** Corner rows transcribed from the GGV5 `v11≤35` tables exactly
   as in `phi_corner4.py`; `(a,b)(j)` are the survey's linear `(m,n)(j)` with
   `a=min`, order verified for `j≥0`. Primary-source.
2. **[unreduced polygon]** Non-`(8,28)` reductions are performed in no paper; the
   standard root-shift + fused-Laurent chart is assumed (`t=l`, `kappa=l-2`,
   `deg C=a0`, `q` from the table) — same boundary as CORNER_144 / PHI_75_125 /
   PHI_CORNER4 / PHI_F14 / PHI_F7.
3. **[branch selection]** The ODE admits multiple residual branches; the pure /
   cofactor / `mu=dg` branch is selected by continuity with the audited pattern
   (`(y+1) | C`) and the root-shift gauge. Real-locus emptiness of `mu=1` at even
   `dg` is annotation, not a complex kill (MU_RUNGS scope correction).
4. **[N-formula]** `N = a[t(a+b-1)+1]-2b` extrapolates the corner-144 bookkeeping;
   `(b-1)/a` integral (F2,F9,F14 less conditional). k=2 (F4,F6,F11) and
   `A0'=(2,0)` (F12,F13) are the more-conditional rows, flagged in the table.
5. **[escapes]** F22–F24 charts exist in no paper; their rows are conditional and
   the F22 IRREGULAR verdict is a within-model statement (negative gap).
