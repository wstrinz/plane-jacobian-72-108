# CHORDAL_PROBE — treewidth / chordal-network diagnostic

**Date:** 2026-07-24. **Status: CLEAN NEGATIVE. THIS LANE IS CLOSED.**

**Headline: the constraint graph of every expanded coefficient system we care
about is the COMPLETE graph.** SUB2 is `K_45` (990 = 45·44/2 edges, treewidth
exactly 44). R9 is `K_28` (378 = 28·27/2 edges, treewidth exactly 27). A
chordal-elimination backend costs `n^{O(treewidth)}`; at treewidth `n−1` it is
*by definition* the dense algorithm with extra bookkeeping. There is no
"different elimination backend" hiding here.

New files: `chordal_probe.py` (the probe), `chordal_probe.json` (numbers),
`chordal_probe_generic.json`, `chordal_probe_symbol.json`,
`chordal_probe_run.log`, `chordal_probe_m2/*.m2` (Macaulay2 inputs), this doc.
READ-ONLY on every existing artifact; nothing else touched.

## 0. Why the question was worth asking, and what it actually asked

Earlier today the sub2 G-system's bipartite incidence / Dulmage–Mendelsohn +
SCC analysis produced **one** strongly connected component containing all 45
spare coefficients, and we concluded there is no topological elimination
ordering to exploit. SCC connectivity and **treewidth** are genuinely
different invariants — a path graph is connected and has treewidth 1 — so
"one SCC" does *not* by itself rule out a chordal decomposition. Macaulay2's
`Chordal` package (Cifuentes–Parrilo) exists precisely for connected systems
with small treewidth: it builds a chordal network and does elimination,
dimension, component and root counting along the clique tree at a cost
exponential in **treewidth**, not in the number of variables.

So the question had real content. It now has a real answer: **the treewidth
is maximal**, and it is maximal for a reason that is structural rather than
accidental (§5).

## 1. Method

`chordal_probe.py` builds, for each target, the **variable-interaction
(constraint / primal) graph**: vertices = variables, edge iff two variables
co-occur in the support of some polynomial. Every polynomial's support is
therefore a clique, which gives the cheap lower bound `tw ≥ maxsupport − 1`.

For upper bounds it runs **min-degree** and **min-fill** greedy elimination
orderings implemented in-file (eliminating `v` forms the clique `{v} ∪ N(v)`,
adds the fill edges, records `|N(v)|`); `max |N(v)|` is a treewidth upper
bound and `max(|N(v)|+1)` the largest maximal clique of the resulting chordal
completion. Cross-checks: `networkx.algorithms.approximation.treewidth_min_degree`
/ `treewidth_min_fill_in` (nx 3.4.2) and exact max-clique via `find_cliques`;
independently, Macaulay2 `Chordal`'s own `constraintGraph` / `chordalGraph` /
`treewidth` / `elimTree`. All four agree on every target.

Targets, and how each is built (all through the audited constructors, not
re-derived):

| tag | source | what it is |
|---|---|---|
| `J6` | `j6_msolve_results.json`, verbatim (gens + class relations + saturation) | the depth-3 saturated systems of the four ALT_HUNT `[J6]` states — a **known kill** (msolve `[-1]`, empty). Pipeline validation. |
| `SUB2` | `face_kill_sweep.build_state_system(..., "sub2")` with `d2` deg 4, `d1` deg 6, `sigma` deg 8, `e = gamma·(y+1)^10` | the sub2 home G-system at full caps: 45 spare unknowns, 22 state params, 122 equations |
| `R9` | `bigrade_annotator.build_R3()` (loads `r9_eliminated_system.json` via `_H_generators`) | the R9 z=1 dm4-eliminated H-system: 28 spare unknowns, 125 equations. The actual wall. |
| `SYMBOL` | `bigrade_annotator._G_generators()` / `_H_generators()` | the **pre-expansion** systems, as polynomials in the eight series symbols `d0,d1,d2,dm1,dm2,dm3,dm4,Phi` |

For SUB2 and R9 two variable models are reported: **spares** (vertices = spare
coefficients only, state parameters treated as generic coefficients — the model
in which "45 spare unknowns" is the variable count) and **full** (spares +
state parameters — the model a solver over ℚ actually faces).

**Fast path, and why it is sound for a negative.** The fully symbolic R9
expansion did not finish in 50 minutes of sympy. The probe therefore also
offers `target_R9_generic`, which runs `build_R3`'s construction with the state
parameters specialised at a random rational point (reusing
`bigrade_annotator._H_generators`, `_phi_stripped`, `_collect_y_equations`
unchanged). Specialisation can only ever **shrink** a support, never grow one,
so the specialised graph is a **subgraph** of the symbolic one — a complete
graph at a specialisation proves the symbolic graph is complete, with no
probabilistic caveat. The fast path was validated against the exact symbolic
SUB2 build, where the exact run is affordable: both give `V=45, E=990,
density 1.0, tw=44`, identical. (The exact symbolic R9 run was left going for
~1 h 50 min without completing `sp.expand` and was then killed; it is not
needed, and it is not quoted anywhere below.)

## 2. The numbers

`tw ≤` is the min-fill upper bound (min-degree agreed on every row); `tw ≥` is
the exact max-clique of the raw graph minus one. **Where the two meet, the
treewidth is exact, not bounded.**

| system | model | V | E | density | polys | max support | tw ≥ | tw ≤ | max clique in completion | fill edges |
|---|---|---|---|---|---|---|---|---|---|---|
| J6 `sub2:a9_b1000…#state3` | ring vars | 4 | 5 | 0.833 | 5 | 3 | **2** | **2** | 3 | 0 |
| J6 `sub2:a8_b1100…#state0` | ring vars | 5 | 8 | 0.800 | 6 | 4 | **3** | **3** | 4 | 0 |
| J6 `sub1:a9_b1000…#state3` | ring vars | 4 | 5 | 0.833 | 5 | 3 | **2** | **2** | 3 | 0 |
| J6 `sub1:a8_b1100…#state0` | ring vars | 5 | 8 | 0.800 | 6 | 4 | **3** | **3** | 4 | 0 |
| SYMBOL G-system | 8 series symbols | 8 | 28 | **1.000** | 4 | 8 | **7** | **7** | 8 | 0 |
| SYMBOL H-system | 8 series symbols | 8 | 21 | 0.750 | 3 | 7 | **6** | **6** | 7 | 0 |
| **SUB2 home G-system** (exact symbolic) | spares | **45** | **990** | **1.000** | 122 | 45 | **44** | **44** | **45** | 0 |
| SUB2 home G-system (exact symbolic) | full | 67 | 2211 | **1.000** | 122 | 67 | **66** | **66** | 67 | 0 |
| SUB2 (params specialised) — *fast-path cross-check* | spares | 45 | 990 | **1.000** | 122 | 45 | **44** | **44** | 45 | 0 |
| **R9 z=1 H-system** (params specialised) | spares | **28** | **378** | **1.000** | 125 | 28 | **27** | **27** | **28** | 0 |

`990 = 45·44/2` and `378 = 28·27/2` and `2211 = 67·66/2`: these are the
complete graphs `K_45`, `K_28`, `K_67` **exactly**. Density 1.000 is not a
rounded 0.9997. `fill edges = 0` because a complete graph is already chordal —
there is nothing to complete, and every elimination ordering is equally bad.

networkx cross-check, verbatim from `chordal_probe.json`:

```
SUB2 home G-system [spares] : {'treewidth_min_degree': 44, 'treewidth_min_fill_in': 44, 'exact_max_clique': 45}
SUB2 home G-system [full]   : {'treewidth_min_degree': 66, 'treewidth_min_fill_in': 66, 'exact_max_clique': 67}
R9 [spares]                 : {'treewidth_min_degree': 27, 'treewidth_min_fill_in': 27, 'exact_max_clique': 28}
```

## 3. Macaulay2 `Chordal` — available, and it agrees

`needsPackage "Chordal"` succeeds: **Chordal 0.2 under Macaulay2 1.19.1** in
WSL. So this is not an "unavailable, reporting graph numbers only" answer —
the package ran, and its own graph machinery reproduces the Python numbers on
all three systems.

J6 validation run (`chordal_probe_m2/j6_stage1b.m2`), verbatim:

```
key sub2:a9_b1000_T1_sz1_dz1_gz-#state3
nvars 4  ngens 5
constraintGraph graph ({E, X, c0_0, w}, {{X, E}, {c0_0, E}, {w, E}, {c0_0, X}, {w, X}})
chordalGraph chordalgraph ({c0_0, E, X, w}, {{c0_0, E}, {c0_0, X}, {E, X}, {E, w}, {X, w}})
treewidth(chordalGraph) 2
elimTree new ElimTree from {children => new HashTable from {null => {w}, c0_0 => {}, E => {c0_0}, w => {X}, X => {E}}, nodes => {c0_0, E, X, w}, parents => new HashTable from {c0_0 => E, E => X, w => null, X => w}, cliques => new HashTable from {c0_0 => {c0_0, E, X}, E => {E, X, w}, w => {w}, X => {X, w}}}
treewidth(elimTree) 2
suggestVariableOrder {c0_0, E, X, w}
```

R9 (`chordal_probe_m2/r9_graph.m2`, 125 generators, 28 spares, parameters
specialised) and SUB2 (`chordal_probe_m2/sub2_generic.m2`), verbatim:

```
nvars 28  ngens 125
constraintGraph edges 378  vertices 28
treewidth(chordalGraph) 27
treewidth(elimTree) 27
```
```
nvars 45  ngens 122
constraintGraph edges 990  vertices 45
treewidth(chordalGraph) 44
treewidth(elimTree) 44
```

Two further facts about the backend, both fail-loud:

**(a) `chordalTria` cannot run in this environment at all.** The package's
triangular-decomposition step shells out to Maple (or Epsilon), neither of
which is installed:

```
sh: 1: maple: Permission denied
/usr/share/Macaulay2/Chordal.m2:1087:9:(2):[10]: error: MapleInterface failed. Maple is needed unless ideal is binomial.
```

`Chordal.m2:1086-1089` confirms the fallback chain: `Maple` unless binomial,
`Epsilon` only in characteristic 0 if present. Our ideals are neither binomial
nor accompanied by a Maple licence. So `codimCount` / `rootCount` / the whole
`chordalTria` half of the package is **unavailable**, independently of the
treewidth verdict. `chordalNet` and `chordalElim` do work.

**(b) On J6 the network is built instantly but does not decompose.** Over
`ZZ/32003` (`chordal_probe_m2/j6_net_gf.m2`):

```
field ZZ/32003
chordalNet built, cpu .0028831
net ChordalNet{...5...}
structure None
chordalElim cpu .184216
net after elim ChordalNet{...5...}
```

`structure None` = neither binomial nor already triangular. The ℚ run
(`chordal_probe_m2/j6_net_qq2.m2`) makes the shape explicit and is the more
telling of the two:

```
chordalNet cpu .0030231
elimTree new ElimTree from {children => new HashTable from {null => {w}, E => {}, c0_0 => {X}, w => {c0_0}, X => {E}},
  nodes => {E, X, c0_0, w}, parents => new HashTable from {E => X, c0_0 => w, w => null, X => c0_0},
  cliques => new HashTable from {E => {E, X, c0_0, w}, c0_0 => {c0_0, w}, w => {w}, X => {X, c0_0, w}}}
--- net BEFORE elim ---
ChordalNet{ E => {{- 6561E^25 - 5832E^22X^5 + ... , - 164025E^25c0_0 - ... , - 1968300E^25c0_0^2 - ... , c0_0^4 + c0_0^3/4 + 5c0_0^2/32 + 15c0_0/128 + 195/2048, E*X*w - 1}}
            X => { }
            c0_0 => { }
            w => { }  }
```

**All five generators sit at the single node `E`; `X`, `c0_0`, `w` are empty.**
The network is a chain with the entire ideal in its first rank, so there is
nothing to propagate — small treewidth buys nothing when every generator
involves the deepest variable. (Note also that `chordalNet I` used the *ring's*
Lex order and got clique `{E,X,c0_0,w}`, i.e. width 3, rather than the width-2
tree `suggestVariableOrder` finds — even on the easy target the default is not
the optimum.)

`chordalElim` over ℚ then **did not return within 900 s** (`exit=124`, the
timeout kill), against 0.18 s over `ZZ/32003`. The cost is rational-coefficient
blowup on degree-25 polynomials — exactly the phenomenon msolve's multi-modular
F4 already disposes of in ~1 s (`J6_MSOLVE.md` §2, "0.5 s (d3)"). **The chordal
backend loses to the incumbent engine on the one target where its structural
precondition is actually satisfied.**

## 4. Is the density the fault of a few fat equations?

No. The probe rebuilds the graph from only the equations of support `≤ k`
(`layered` in the JSON). Even after discarding the fattest two thirds of the
system the graph is still near-complete:

**R9 (28 spares, 125 equations):**

| support cap | polys kept | edges | density | isolated vertices | tw ≤ |
|---|---|---|---|---|---|
| 2 | 6 | 2 | 0.005 | 24 | 1 |
| 7 | 18 | 30 | 0.079 | 16 | 5 |
| 14 | 42 | 182 | 0.481 | 1 | 13 |
| 21 | 60 | 314 | 0.831 | 0 | **19** |
| 28 (all) | 125 | 378 | 1.000 | 0 | **27** |

**SUB2 (45 spares, 122 equations):**

| support cap | polys kept | edges | density | isolated vertices | tw ≤ |
|---|---|---|---|---|---|
| 3 | 8 | 6 | 0.006 | 39 | 2 |
| 11 | 24 | 72 | 0.073 | 27 | 8 |
| 22 | 56 | 420 | 0.424 | 4 | 20 |
| 33 | 88 | 846 | 0.855 | 0 | **32** |
| 45 (all) | 122 | 990 | 1.000 | 0 | **44** |

Sparsification does not rescue the lane: throwing away 65 of R9's 125
equations still leaves treewidth ≥ 19 on 28 variables. And the discarded
equations are exactly the ones carrying the information — the whole system's
content is in the fat middle.

The support histograms say why. R9's is a perfect arithmetic ladder,
`{2:6, 4:6, 6:6, …, 26:6, 27:6, 28:41}` — six equations at each even support
from 2 up to 26, then **41 equations of full support 28**. SUB2's is the same
shape in steps of three, `{3:8, 6:8, 9:8, …, 33:8, …, 45:3}`.

## 5. Why it is complete — the structural reason, not an accident

The ladder is the signature of **polynomial convolution**. A generator like
`G1 = (3/2)·d1·dm1² + 3·d2·dm1·dm2 + 3·dm1·dm4 + 3·dm2·dm3` becomes, after the
spare ansätze are substituted and the result expanded in `y`, one equation per
`y`-order. The `y^j` coefficient of the product `dm2·dm3` is
`Σ_i R_i·S_{j−i}` — at the **extremes** of the `j` range only one or two terms
survive (support 2, 3), and in the **middle** every `R_i` and every `S_i` is
present at once. Add `dm2·dm4`, `dm3²`, `d0·dm2·dm3`, `d1·dm2·dm4` from `G2`,
`G3`, `G5body` and one middle-order equation touches all three spare series
simultaneously: support 45.

That is why the ladder rises linearly and then saturates, and it is why no
change of ansatz, degree cap or state will help. **Convolution of capped series
produces complete constraint graphs. The only sparse equations any such system
will ever have are the ones at the extremes of the `y`-order range.**

This also disposes of the third branch of the probe's interpretation — "small
cliques only *before* coefficient expansion, so the generic-`K(y)`
representation is the better abstraction". Measured, the symbol-level system is
*also* complete:

```
SYMBOL G-system (G1,G2,G3,G5)   V=8  E=28  density=1.0   tw = 7 (exact)
     G1   support  6: d1, d2, dm1, dm2, dm3, dm4
     G2   support  6: d0, d2, dm1, dm2, dm3, dm4
     G3   support  6: d0, d1, dm1, dm2, dm3, dm4
     G5   support  8: Phi, d0, d1, d2, dm1, dm2, dm3, dm4
SYMBOL H-system (H2,H3,H5)      V=8  E=21  density=0.75  tw = 6 (exact)
     H5   support  7: Phi, d0, d1, d2, dm1, dm2, dm3
```

`G5` alone is `K_8`. The symbol level is not *sparse*; it is merely *small*.
Its treewidth is `n−1` too — density just costs nothing when `n = 8`. So the
generic-`K(y)` abstraction is not recommended by *this* measurement: whatever
its other merits, it is not a treewidth story, and it should not be pitched as
one.

## 6. Verdict

Against the three interpretations set for this probe:

- **"small cliques on J6 AND R9 → build the backend."** ✗ Refuted. J6's
  cliques are small (3, 4) but only because J6 has 4–5 variables and density
  0.83; R9's maximum clique is **28 of 28 variables**.
- **"small cliques only before coefficient expansion → generic-K(y) is the
  better abstraction."** ✗ Refuted (§5). The symbol level is complete as well
  (`K_8`); it is small, not sparse.
- **"giant cliques everywhere → close this lane fast, and say so."** ✓ **This
  one.** `K_45`, `K_28`, `K_67`, `K_8`. Treewidth is exactly `n − 1` on every
  system with more than five variables.

**CLOSE THE LANE.** Chordal / treewidth-based elimination offers this campaign
nothing. The earlier one-SCC finding and this one are not the same result, but
they now point the same way, and this one is the stronger of the two: not only
is there no elimination ordering, there is no ordering that is even *better
than the worst*, because all `n!` of them produce the same maximal clique. A
`Chordal`-package backend on R9 would be `n^{O(27)}` on 28 variables — the
dense computation, with a clique tree bolted on. Do not build it.

**Anti-overclaim.** Three things this does NOT show. (i) It says nothing about
whether R9 is solvable — only that *this* structural handle is absent.
(ii) Treewidth is a property of the constraint graph, not of the ideal;
a different *formulation* (new variables, a different ansatz, a factored
representation that never expands the convolution) could have a different
graph. What is ruled out is applying a chordal backend to **these expanded
coefficient systems**, which is what was proposed. (iii) The `Chordal`
package's `chordalTria` half was never exercised, because Maple is absent —
but with treewidth `n−1` it could not have helped, so this is a footnote, not a
gap.

**Where the sparsity actually is, and why that is not a lead either.** The
layered tables and the support histograms locate every sparse equation these
systems have: they all live at the **extremes of the `y`-order range** (support
2 and 3 — six such equations in R9 leaving 24 variables isolated, eight in SUB2
leaving 39 isolated). That is exactly the extremal-face structure
`bigrade_annotator` / `face_kill_sweep` were built to attack, and this probe
independently confirms it is the *only* sparsity anywhere in the system.

But that is a closed door too, not an open one. `FACE_KILL_SWEEP.md` §6.1
(landed separately today) proves the face detector's reach is *exactly* the
Phi-depth criterion and nothing more, for a structural reason — every `G5body`
term contains a spare, and `G1`,`G2`,`G3` each carry a pure-spare product
sitting at their cap (`dm2·dm3 = 26`, `dm2·dm4 = dm3² = 28`, `dm3·dm4 = 30`),
so their extreme faces are never spare-free. So the correct joint reading is
harsher than "keep using the face detector": **the extremes are the only sparse
part of the system, and the extremes have already been fully mined.** This
probe adds the complementary half — that the non-extreme part is not merely
hard but *maximally* dense — and the two together say the sparsity structure of
these systems is exhausted.

## 7. Reproduce

```
python chordal_probe.py --only J6,SYMBOL,SUB2,GENERIC   # -> chordal_probe.json  (~6 min)
python chordal_probe.py --only R9                        # exact symbolic R9; did NOT finish in 1h50m
python chordal_probe.py --emit-m2                        # -> chordal_probe_m2/*.m2
```

`GENERIC` is the randomised-parameter path (§1) and covers R9 in ~2.5 min;
`--only R9` is the exact symbolic build and is not recommended.

Macaulay2 (WSL; `/mnt/c` is broken in this WSL, so pipe in via stdin and run
from `$HOME`):

```
wsl.exe -e bash -lc 'cat > $HOME/r9_graph.m2' < chordal_probe_m2/r9_graph.m2
wsl.exe -e bash -lc 'cd $HOME && M2 --script r9_graph.m2; cat r9_graph.txt'
```
