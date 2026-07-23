# NEXT_CASES: the candidate landscape above 125

Prepared 2026-07-21 from the full LaTeX sources of:

- **[GGHV22]** arXiv:2204.14178, Guccione–Guccione–Horruitiner–Valqui,
  *"Increasing the degree of a possible counterexample to the Jacobian Conjecture
  from 100 to 108"* — saved at `paper_src/2204.14178.tex` (2262 lines; recovered from a
  previous session's scratchpad, `\end{document}` present, complete).
- **[GGV5]** arXiv:1708.07936, Guccione–Guccione–Horruitiner–Valqui,
  *"Some algorithms related to the Jacobian Conjecture"* — fetched fresh from
  `https://arxiv.org/e-print/1708.07936`, saved at `paper_src/1708.07936_GGV5.tex`
  (2045 lines). This is the paper whose §§5–6 tables the 2204 paper's classification
  is read from.

Line numbers below refer to these two files.

---

## 1. Where 100, 108, 125 come from — exact quotes

**Abstract of [GGHV22]** (`2204.14178.tex` line 229):

> "We list all the pairs $(\deg(P),\deg(Q))$ with $\max\{\deg(P),\deg(Q)\}< 125$ for any
> hypothetical counterexample to the plane Jacobian Conjecture and discard them all,
> except the pair $(72,108)$ (and the symmetric pair $(108,72)$), thus we confirm the
> lower bound of 100 obtained by Moh and raise it up to 108."

**Introduction** (lines 251–252):

> "The only exception is the case $(\deg P,\deg Q)=(72,108)$, and so, if one manages to
> discard this case, it would increase the lower bound from $108$ up to $125$."

(lines 254–259):

> "In section 1 we list all the cases with $\max\{\deg(P),\deg(Q)\}<125$, following
> [GGV5]. There are 10 cases, from which we consider 5 to be discarded already by
> previous work of different authors. If we assume $\deg(P)<\deg(Q)$, then the remaining
> 5 cases satisfy $(\deg(P),\deg(Q))\in\{(56,84),(66,99),(72,108),(80,120)\}$, where
> there are two cases with $(\deg(P),\deg(Q))=(72,108)$."

(lines 267–268) — the open case:

> "For the other case with $(\deg(P),\deg(Q))=(72,108)$ we couldn't solve the
> corresponding system of polynomial equations, thus it is left open."

(lines 276–277) — the origin of "125 with enough computing power":

> "…up to 108. With enough computing power we would be able to raise it up from 108 to
> 125, since there is only one case left."

**The main theorem** (lines 286–288):

> "If $(P,Q)$ is a counterexample to the Jacobian Conjecture, then we have either
> $\max\{\deg(P),\deg(Q)\} \geq 125$, or $(\deg(P),\deg(Q))\in\{(72,108),(108,72)\}$."

**The 10-case table below 125** (lines 302–320), with the paper's "Discarded?" column.
Columns: corner $A_0$, $(m,n)$, $\max\deg$:

| $A_0$   | $(m,n)$ | max deg | status in [GGHV22] |
|---------|---------|---------|--------------------|
| (4,12)  | (3,4)   | 64      | [GGV4 §3.5], [Moh], [Heitmann] |
| (4,12)  | (5,7)   | 112     | [GGV4 §3.5] |
| (5,20)  | (2,3)   | 75      | [GGV3 §5]; "no detail in [M]" |
| (5,20)  | (3,2)   | 75      | [GGV3 §5]; "no detail in [M]" |
| (7,21)  | (2,3)   | 84      | "no detail in [M]" — discarded in [GGHV22] §2 (and again §5) |
| (8,24)  | (2,3)   | 96      | [GGV5 Prop 6.1] |
| **(8,28)** | **\*(3,2)** | **108** | **left open — our (72,108) case** |
| (8,32)  | (3,2)   | 120     | discarded in [GGHV22] §2 |
| (9,24)  | (2,3)   | 99      | "no detail in [M]" — discarded in [GGHV22] §4 |
| (9,27)  | (2,3)   | 108     | discarded in [GGHV22] §4 |

The two $(72,108)$ cases are $A_0=(9,27)$, $(m,n)=(2,3)$ (killed, §4) and
$A_0=(8,28)$, $(m,n)=(3,2)$ (open — the star in the table; this is the HANDOFF.md
"Prop 4.3 case (8,28)", reduced at `2204.14178.tex` lines 1000–1007 to
$[P,Q]=x^2$ with the two Newton-polygon subcases we call subcase 1 / subcase 2).

Provenance of the 10 cases (lines 291–292):

> "Following [GGV5], we take the smallest members of the families
> $F_1,F_2,F_3,F_9,F_{17},F_{22}$ in section 5 of [GGV5] and three additional cases in
> the tables of section 6 of [GGV5]."

---

## 2. The candidate-pair recipe (from GGV5)

All of this is algorithmic; [GGV5] even describes a C++/PostgreSQL/D3.js implementation
(`1708.07936_GGV5.tex` lines 1638–1666).

**Step 0 — degrees from corners.** If $(P,Q)$ is a counterexample there are coprime
$m,n>1$ and a corner $A_0=(a,b)$ with (GGV5 line 250):

> "the support of $P$ is contained in the rectangle with vertices
> $\{(0,0), m(a,0),m(a,b),m(0,b)\}$ … Note that $\deg(P)=m(a+b)$ and $\deg(Q)=n(a+b)$."

So $(\deg P,\deg Q) = (m\,v_{11}(A_0),\; n\,v_{11}(A_0))$ with $v_{11}(a,b)=a+b$.

**Step 1 — admissible complete chains.** Algorithm "GetCompleteChains" (GGV5 §"Main
algorithm") enumerates all *admissible complete chains*
$(\mathcal C_0,\dots,\mathcal C_j,\mathcal A_{j+1})$ of corners with $v_{11}(A_0)\le M$;
these encode the possible Newton-polygon corner geometry (Definitions "complete chain"
and "cond div" in GGV5). For $M=35$ there are 14 chains of length 1 and 2 of length 2
(line 1671). There are none with $v_{11}(A_0)<16$ (line 407 — reproving GGV1's
$\deg\ge 16\cdot\min(m,n)$-type bound).

**Step 2 — the Diophantine condition on $(m,n)$.** For a final corner
$\mathcal A=(a/l,b)$, define (GGV5 lines 1500–1513, Definition "mn families"):

$$I(\mathcal A) = \Big\{k\in\mathbb N:\ 1\le k<l-\tfrac ab,\ \gcd\Big(b,\tfrac{bl-a}{\gcd(k,bl-a)}\Big)=1\Big\},$$
$$\MN_k(\mathcal A)=\{(m,n):\ m,n>1,\ \gcd(m,n)=1,\ (m+n)bk-n(bl-a)=k\},$$

(the "ecuacion diofantica", lines 1449–1455, Proposition `extremosfinales1`; the
symmetric equation with $m\leftrightarrow n$ gives the swapped cases). Solutions form
arithmetic families $(m,n)=(m_0+j\Delta^{(1)},\,n_0+j\Delta^{(2)})$, $j\in\mathbb N_0$,
with $\Delta^{(1)}=\frac{bl-bk-a}{e_k}$, $\Delta^{(2)}=\frac{bk}{e_k}$,
$e_k=\gcd(k,bl-a)$ (lines 1519–1538, and Algorithm "GetmnFamilies", lines 1600–1635).

**Step 3 — tables.** Running this with $M=35$ gives 24 infinite families $F_1..F_{24}$
(GGV5 lines 1678–1717); families $F_{18}$–$F_{21}$ are then killed for *all* $(m,n)$ at
once by a structural argument (lines 1726–1786), and $F_{22}$'s smallest member by
Proposition "caso antisimetrico" (lines 1874–1927). Separately, GGV5 §"Possible
counterexamples with max ≤ 150" (lines 1792–1872) lists **all 34 cases with
$\max(\deg P,\deg Q)\le 150$**: 13 from the families plus 9 sporadic length-1 chains,
11 length-2 chains, and 1 length-3 chain (sporadic = corners with $v_{11}(A_0)>35$
that only enter below 150 for small $(m,n)$).

**Worked check of the recipe** (script: `paper_src/next_cases.py`; run output verified):

- F2's corner $(7/5,2)$: $bl-a=3$, $k=1$: $(m+n)\cdot2-3n=1 \Rightarrow n=2m-1$, giving
  $(2,3),(3,5),(4,7),\dots$ — exactly F2's $(j+2,2j+3)$. ✓
- The open corner $(8,28)$, final corner $\mathcal A_1=(11/4,7)$: $bl-a=17$,
  $k\in\{1,2\}$; $k=1$: $7m-10n=1 \Rightarrow (3,2),(13,9),(23,16),\dots$;
  $k=2$: $14m-3n=2 \Rightarrow (7,32),\dots$ So on this corner the open case $(3,2)$
  is followed by $(13,9)$, i.e. degrees $(468,324)$ — far away.
- $(8,28)$'s other final corner $(7/4,3)$: $3m-2n=1 \Rightarrow (3,4),(5,7),\dots$,
  giving the $(108,144)$ case at 144 and $(180,252)$ next. ✓

---

## 3. The computed next candidates at/above 125 (up to 150)

Recomputed independently with `paper_src/next_cases.py` from the family
parametrizations + sporadic chain tables; the script reproduces **exactly** the 10
cases below 125 of [GGHV22]'s table, and 24 cases in $[125,150]$ — consistent with
GGV5's count of 34 total ≤ 150. Sorted by max degree, in the $\deg P<\deg Q$
convention (each row is one *case*; equal degree pairs from different chains are
genuinely distinct cases):

| max | (deg P, deg Q) | source | $A_0$ | $(m,n)$ |
|-----|----------------|--------|-------|---------|
| **125** | **(75,125)** | family $F_2$, $j{=}1$ | (5,20) | (3,5) |
| 126 | (84,126) | chain len 1 | (7,35) | (2,3) |
| 126 | (84,126) | chain len 2 | (12,30) | (3,2) |
| 128 | (96,128) | family $F_{24}$, $j{=}0$ | (8,24) | (3,4) |
| 132 | (88,132) | chain len 1 | (11,33) | (2,3) |
| 135 | (90,135) | chain len 1 | (9,36) | (2,3) |
| 135 | (90,135) | chain len 1 | (9,36) | (3,2) |
| 135 | (90,135) | chain len 1 | (12,33) | (2,3) |
| 135 | (90,135) | chain len 2 | (9,36)→(9,24) | (2,3) |
| 140 | (56,140) | family $F_{11}$, $j{=}0$ | (7,21) | (2,5) |
| 140 | (84,140) | family $F_9$, $j{=}1$ | (7,21) | (3,5) |
| 144 | (96,144) | chain len 2 (×4 chains) | (12,36) | (2,3) |
| 144 | (96,144) | chain len 2 | (8,40) | (3,2) |
| 144 | (96,144) | chain len 3 | (12,36) | (3,2) |
| 144 | (108,144) | chain len 1 | **(8,28)** | (3,4) |
| 147 | (42,147) | family $F_7$, $j{=}0$ | (6,15) | (2,7) |
| 147 | (63,147) | family $F_8$, $j{=}0$ | (6,15) | (3,7) |
| 147 | (98,147) | chain len 1 (×2: (2,3) and (3,2)) | (7,42) | — |
| 150 | (100,150) | chain len 2 (×2 chains) | (10,40) | (3,2) |

Key observations:

1. **The next case after (72,108) is a single case: (75,125)** — family $F_2$,
   $A_0=(5,20)$, $(m,n)=(3,5)$. Its $j=0$ sibling is Moh's 75 case, discarded by the
   GGV3 polynomial-system method. So the first case above 125 is the "next rung" of an
   already-climbed ladder, a natural first test of any inductive/uniform argument.
2. **The open corner $(8,28)$ reappears at 144** with $(m,n)=(3,4)$ (degrees
   $(108,144)$). A corner-level structural kill of $(8,28)$ — as opposed to a kill of
   the specific $(3,2)$ system — would remove both the 108 and the 144 case.
3. **Case count grows quickly**: 10 cases below 125, 24 more in the next 25-degree
   band. Beyond 150 the tables of GGV5 stop; the families are infinite (their next
   members: $F_1(7,10)\to160$, $F_{22}(3,5)\to160$, $F_{16}(3,5)\to165$,
   $F_2(4,7)\to175$, $F_3(7,5)\to175$, …) and new sporadic corners keep entering, so an
   enumeration above 150 requires rerunning GGV5's Algorithm "Main algorithm"
   with a larger bound $M$ — mechanical, but not done in any published paper.

---

## 4. The paper chain (what each proved)

| Ref | Paper | Contribution |
|-----|-------|--------------|
| [K] Keller 1939 | Monatsh. Math. Phys. 47 | the conjecture |
| [A] Abhyankar 1977 | Tata Lectures 57 | expansion techniques background |
| [M] Moh 1983 | J. Reine Angew. Math. 340, 140–212 | no counterexample with $\max\deg\le 100$; "provides a detailed proof only for the smallest case" (2204 line 245–246); lists 4 cases (= 6 in GGV terminology, GGV5 line 1794) |
| [H] Heitmann 1990 | JPAA 64, 35–72 | independent computation of the 64 case; probably implicitly the $F_{22}$ 96 case "by symmetry reasons" (GGV5 line 1818) |
| [O] Orevkov 2001 | Tr. Mat. Inst. Steklova 235 | analyzed the $F_{13}$, $j=1$ possible counterexample (GGV5 lines 1788–1790); "counterexamples to the JC at infinity" |
| [LCW] Wang 2005 | Taiwanese J. Math. 9, 421–431 | case list cross-referenced for $F_{22}$ (GGV5 line 1818) |
| [GGV4] 2013 | Pro Mathematica 27 | differential-equation method; discarded 64 and $(80,112)$ |
| [GGV3] arXiv:1406.0886 | *A system of polynomial equations related to the JC* | the polynomial-system machine; §5 discarded both 75 cases |
| [GGV1] 2017 | J. Algebra 471, 13–74, *On the shape of possible counterexamples…* | support/corner theory: $(\rho,\sigma)$-leading forms, regular corners, $\deg\ge16$-type bound |
| [GGV2] arXiv:1605.09430 | *…Lower Side of the Newton Polygon* | lower-side constraints; Prop 3.29 / Rmk 3.31 corner exclusions (used for 120) |
| [GGV5] arXiv:1708.07936 | *Some algorithms related to the JC* | **the enumeration**: algorithms + tables of all 34 cases ≤ 150; killed $F_{18}$–$F_{21}$ (all $(m,n)$!) and the $F_{22}$ 96 case (Prop 6.1) |
| [GGV6] 2019 | Pro Mathematica 30, *Approximate roots and intersection numbers* | intersection-number inequality $I_M\ge I_m$ (Thm 7.3, used for 84) |
| [GGHV22] arXiv:2204.14178 | *from 100 to 108* | killed 120 (§2, corner argument), 84 (§2 via GGV6; §5 again by systems), 99 and one 108 (§4, systems + CAS elimination); reduced open (8,28) to $[P,Q]=x^2$ + two polygon subcases (§3); **left (72,108) open** |

---

## 5. How cases die, and prospects for an inductive/automatable pattern

**Observed kill mechanisms** (each case so far died by one of five):

1. *Corner exclusion / last-lower-corner arguments* (GGV2 Prop 3.29 & Rmk 3.31):
   killed 120 in half a page (`2204.14178.tex` lines 327–388) and $F_{18}$–$F_{21}$
   *for every* $(m,n)$ (GGV5 lines 1726–1786).
2. *Arithmetic of the Diophantine data*: GGV5's "caso antisimetrico" (lines
   1874–1927) kills $F_{22}(2,3)$ because $1=5\frac{k-m}{k}$ has no solution — pure
   number theory on $(k,m,n,b,l,a)$.
3. *Intersection-number inequality* (GGV6 Thm 7.3): killed 84 (lines 391–455).
4. *Newton-polygon reduction* (§3 of [GGHV22], following GGV1 automorphisms
   $x\!\leftrightarrow\! y$, $y\mapsto y+\lambda x^{-k}$, $x\mapsto x^{-1},y\mapsto yx^k$):
   does not kill by itself but shrinks supports drastically (e.g. (7,21) to corners
   $\{(0,0),(2,0),(3,1),(0,7)\}\times(2,3)$, lines 1313–1396).
5. *Polynomial systems + CAS elimination* (GGV3 method, §§4–5 of [GGHV22]): construct
   $C$ with $P=C^2$, $Q=C^3+\dots+F$ term by term (Prop "calculo de C", lines
   1416–1526), extract finitely many coefficient equations, eliminate with a CAS
   ("using a CAS (for example Mathematica) we eliminate the variables
   $d_{-10},\dots,d_{-2}$", lines 1689, 2056), then reach a univariate contradiction
   (degree/multiplicity mismatch, e.g. lines 1769–1785, 2085–2087). Killed 99, one 108,
   and (again) 84.

**Evidence for automatability.** The *enumeration* side is already fully automated
(GGV5 wrote C++ + SQL + a website). The *kill* side is semi-mechanical: mechanisms
(1)–(2) are cheap arithmetic checks on the chain data and could be run wholesale at any
bound; mechanism (5) is an algorithm in principle (polygon reduction → $C$-series →
finite system → elimination) whose only obstacle is computational: the paper's own
words for the surviving case are "we couldn't solve the corresponding system of
polynomial equations" and "With enough computing power we would be able to raise it up
from 108 to 125" (lines 268, 276). Our project's experience with (72,108) —
where the system's factors f31/f37 have 102/618 terms — shows the growth rate.

**Evidence for an inductive pattern.** Two genuine precedents for killing infinitely
many cases at once:

- GGV5's $F_{18}$–$F_{21}$ argument is uniform in $(m,n)$ — an entire infinite family
  dies by one corner computation.
- The Diophantine structure organizes all $(m,n)$ on a fixed corner into arithmetic
  progressions $(m_0+j\Delta^{(1)}, n_0+j\Delta^{(2)})$; any argument that only uses
  the corner data $(a/l,b,k)$ and congruences (mechanism 2) automatically covers all
  $j$.

**Against**: the cases that resist structural kills (Moh's cases, 96, 99, both 108s)
have so far each needed a bespoke system-of-equations endgame whose size grows with
$v_{11}(A_0)\cdot\max(m,n)$, and the per-band case count is growing (10 below 125 → 24
in the next band). The papers contain **no** stated conjecture of a uniform mechanism;
[GGHV22]'s own summary (lines 273–277) is explicitly the opposite of a single method:

> "all the techniques developed in each of the articles [GGV2], [GGV3] and [GGV5] are
> very useful in discarding some families of possible counterexamples. However, the
> Jacobian conjecture is a very hard problem, so we have to combine these techniques…"

**Bottom line for our program.** Closing (72,108) gives the clean statement "no
counterexample below degree 125". The wall above it is not one case but a thickening
band: (75,125) first (a family-$F_2$ ladder rung, plausibly amenable to the same GGV3
machinery that killed its $j=0$ sibling at 75), then 2 cases at 126, and 24 cases
before 150. A realistic post-(72,108) roadmap in increasing ambition: (i) rerun
mechanisms (1)–(2) wholesale above 125 to see which of the 24 die cheaply; (ii) attempt
a uniform-in-$j$ version of the GGV3 system argument on a full family ($F_2$ being the
test case); (iii) a corner-level kill of $(8,28)$ that covers both its 108 and its 144
incarnation would be the most direct payoff of whatever exact certificate closes our
case (T5), since our f31/f37 analysis lives exactly at that corner.
