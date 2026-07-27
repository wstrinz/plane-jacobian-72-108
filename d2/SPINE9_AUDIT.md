# SPINE9_AUDIT.md — independent audit of the five-cell closure

> ## VERDICT: **CONFIRMED.**
>
> All five cells `a9_b{0000,1000,1100,1110,1111}_T1` are **EMPTY**, conditional
> on `a_t = 9`. Every load-bearing step of `SUB1_SPINE9.md` reproduced
> independently; **no step was found to be wrong, weakened, or vacuous**, and
> no cell survives.
>
> Checker: `spine9_audit.py` (`--quiet`, exit 0; **81/81**; ~12 s, pure sympy).
> Written from `generators.json`, `slice_obstruction_stage.json` and the runtime
> cap tables. `sub1_spine9.py` and `spine.py` were **not imported, not executed
> and not read** before the mathematics here was reconstructed; the audit's own
> ansatz, ideal formulation, second decision procedure and controls are its own.
> Read-only: `git status` shows **no tracked file modified**.
>
> **Two things this audit ADDS, and one of them matters a lot:**
>
> 1. **`SUB1_SPINE9.md` §10.2's "single integer needing a second party" now has
>    one, and it is not an integer that needs trusting.** `v_t(h_7) >= 11` — and
>    `v_t(h_6) >= 10` with it — is a two-line **consequence of the AUDITED rows**
>    `v_t(h_1..h_5) >= (1,3,5,7,9)` plus the P-side absorption form
>    `h_n = -q_n/2 + t^(2n-2) g_n`, `q_n = sum_{i=1..n-1} h_i h_{n-i}`, which
>    `SLICE_OBSTRUCTION_AUDIT.md` §4 itself lists as an imposed condition `(P<)`
>    and reproduces. `v_t(h_6) >= min(1+9, 3+7, 5+5, 10) = 10`,
>    `v_t(h_7) >= min(1+10, 3+9, 5+7, 12) = 11`, `v_t(h_8) >= 12`. Check `G10`.
>    So `a9_b0000_T1`'s cascade dependency rests on **audited** material, and
>    the exposure that `SUB1_SPINE9.md` §10.2 flags as the file's weakest point
>    is substantially smaller than it says.
> 2. **A non-vacuity control the audited lane does not have.** Its `X2`/`X2b`
>    exercise the support test only at `k = 1` and `k = 4`. `k = 2` is the one
>    value where "infeasible" is *not* structurally forced, so a silent bug
>    there would be invisible. Check `E6` builds a synthetic squarefree quartic
>    (`q(-1) = 3 != 0`) on which the **same** `k = 2` test returns **FEASIBLE**,
>    at the `z` the construction was designed for. Check `E7` shows, in the other
>    direction, that `k = 3` is infeasible on unrelated quartics too — i.e. that
>    kill is structural, not arithmetic, which is a weaker guarantee and is said
>    so below.
>
> **STANDING DEPENDENCY, stated as required.** Everything here is conditional on
> `a_t = 9`. Its lower half `a_t >= 9` is independently audited
> (`slice_obstruction_audit.py`); its **upper half `a_t <= 9` is same-author and
> under separate concurrent audit**. This verdict is conditional on that audit
> succeeding. If `a_t <= 9` falls, this file says nothing about the frontier.

---

## 0. What was audited, and against what

The claim: the last five enumerated `f31` sub1 frontier cells are empty. The
five are exactly the `a9_*` entries of `FRONTIER_REBUILD.md`'s 34-cell surviving
list, so "these five" is the right five.

Re-posed and re-decided here, from scratch:

| the lane's step | this audit's independent route | check |
|---|---|---|
| `P1` K-syzygy | recomputed `2*(G5body+Phi + d2*G3 + d1*G2 + d0*G1)` from `generators.json` | `A2`, `A3` |
| D1/D2 divisor filter | **PROVED here** from the syzygy + `Phi = c*t^30*q` factorisation, not imported from `divisor_filter.py` | `A13`, `A14` |
| `P3`/`P4` the `a=9` reduction | substitution + exact division in the **free ring** `Q[gamma,c,t,Pi,Qc,A,B,C,d0,d1,d2]` | `B1`–`B5` |
| `P7` the cofactor identity | independently re-derived; also re-derived with `t^9 -> t^a` free | `C3`–`C6` |
| `P7.c` on-variety confirmation | 12 genuine points of `{g1hat=g2hat=g3hat=0}` built by solving `g2hat` for `d0` and `g3hat` for `v` | `C7` |
| `P12` marked-support test | **two** decision procedures: a saturated ideal over `Q` posed my own way, **and** an exact rank-1 test over an explicit subfield of the splitting field | `E1`–`E7` |
| `P13`–`P15` degree ledger | re-derived from caps read at runtime; the `k=0` dichotomy re-posed as an explicit leading-term uniqueness search | `F1`–`F14` |
| `P9`/`P10` shift + `z` window | the shift re-derived from one generalized-binomial model that reproduces **all five** recorded forms at once | `G1`–`G11` |
| `P16` `t^9 \| R,S,T` | re-derived **here** from `e \| S` + the K-syzygy alone, cap-free; and again from the cascade profile | `G12`, `G13`, `G9` |
| `X1`–`X5` | all five re-posed and re-run | `B6`, `E3`, `E5`, `F12`–`F14`, `G4`, `G6` |

---

## 1. Where the brief said to push hardest

### 1.1 The cofactor identity — **CLEARED, decisively**

This was the highest-value check and it comes out clean. Re-derived from
`generators.json` alone (`C3`):

```
F*Z - (1/6)*gamma^5*t^9*Pi^4  ==  -gamma*A*g2hat + gamma^2*Pi*g3hat
```

residual exactly `0`, with `F = A*(u+2v)+w`, `Z = A^2 - gamma*Pi^2*v`,
`u = gamma*d2`, `w = (1/2)*gamma^2*d1*Pi`.

Three separate facts rule out the feared inheritance of SPINE's sub2 zero-slack
coincidence `(n+6)+(2n+4) = 3n+10`:

* **`C5`** — both sides have free symbols inside `{gamma, t, Pi, A, v, d0, d1, d2}`.
  It is an identity in a polynomial ring in those indeterminates. No degree, no
  cap, no slack count can enter a statement of that form. The `3n+10` count is
  a relation among *degrees*; there are no degrees here.
* **`C6`** — the identity survives replacing `t^9` by a free symbol `t^a`. It is
  the generic `a`-family identity. Nothing `a = 9`-specific, still less
  `a + n = 10`-specific, is being used.
* **`C4`** — `d0` is present in `g2hat` and in `g3hat` separately and cancels
  exactly in the combination. That cancellation is what makes the identity
  possible and it is a coefficient identity, not a degree count.

`SESSION_HANDOFF.md`'s warning is correct about SPINE and correctly **not
applicable** here.

### 1.2 `gcd(Z,Pi) = 1 => Z | t^9` — **CLEARED**

The chain, re-derived step by step:

* `kbox|_{Pi=0} = 3*A^2 - mu*t^3*Q` (`D1`). At a marked root
  `3*A(r)^2 = mu*(r+1)^3*Q(r)`, a product of three provably nonzero factors:
  `mu = 2c/gamma != 0` (`c = -1/6630`, `gamma = lc(e) != 0`);
  `(r+1) != 0` because `q(-1) = 3315` (`A6`, recomputed);
  `Q(r) != 0` because `q` is squarefree (`A5`, recomputed). So `A(r) != 0`.
* `g1|_{Pi=0} = A*B` (`D3`), so `B(r) = 0` at every marked root and `Pi | B`
  (`Pi` squarefree). `Z|_{Pi=0} = A^2 != 0` (`D4`), so `gcd(Z,Pi) = 1`.
* The right side of `(*)` is `(1/6)*gamma^5 * t^9 * Pi^4`. **`1/6` and
  `gamma^5` are units of the coefficient field** — `gamma` is the leading
  coefficient of `e`, a nonzero *scalar*, not a polynomial — and `t` is prime
  in the UFD `K[y]` with `gcd(t,Pi) = 1`. Nothing else hides in the right side.
  Hence `Z = zeta*t^z`, `0 <= z <= 9` (`D6`). **No factor is unaccounted for.**

At `k = 0` the marked-root step is vacuous and `gcd(Z,1) = 1` trivially, so the
chain is not silently assuming `k >= 1` (`D5`).

### 1.3 The saturated-ideal support test — **CLEARED by a second, disjoint decision procedure**

The lane's route is a Gröbner/saturation argument; a unit ideal over `Q` really
does mean no solution in any commutative `Q`-algebra, and `q` really is
irreducible over `Q` (`A4`, recomputed; `galois_group(q) = S4`, so a degree-1/2/3
factor genuinely lives in an extension and the ideal posing is the right one).

Rather than trust that, I decided the same question a **second, structurally
different way** and required the two to agree on all 40 `(k,z)` pairs (`E1`,
they do). The second route is exact linear algebra in an explicit subfield of
the splitting field:

* for fixed `Pi`, the condition `Pi^2 | (mu*t^3*Q - 3*zeta*t^z)` is **linear**
  in `(mu, zeta)`; write `M` for the `2k x 2` matrix whose columns are
  `rem(t^3*Q, Pi^2)` and `rem(-3*t^z, Pi^2)`;
* `col_zeta != 0` always (`gcd(Pi,t)=1`, `deg Pi^2 > 0`) and `col_mu != 0`
  always (`q` squarefree so `gcd(Pi,Q)=1`, hence `Pi^2` does not divide `Q`) —
  both **asserted at runtime**, so the criterion cannot be silently vacuous;
* therefore a solution with `mu != 0` **and** `zeta != 0` exists **iff
  `rank M = 1`**, i.e. iff every `2x2` minor vanishes. No saturation needed, no
  normalisation needed.

`Pi` is realised exactly: `k = 4` over `Q`; `k = 1, 3` over `Q[x]/(q)` (degree 4
— for `k=3`, `Pi = q/(2048(y-x))`); `k = 2` over `Q[x]/(sextic)` where the
sextic `32768p^6+24576p^5+16384p^4+5632p^3-10080p^2-2680p-495` is the resolvent
for `p = -(r_i+r_j)` and is **irreducible over `Q`** (verified), so its one
Galois orbit covers all six 2-subsets. In each case `q = Pi*Q` is asserted to
hold identically in the constructed field before the test runs.

Result, both routes agreeing:

| `k` | cell | feasible `z` in `[0,9]` |
|---:|---|---|
| 1 | `a9_b1000_T1` | **none** |
| 2 | `a9_b1100_T1` | **none** |
| 3 | `a9_b1110_T1` | **none** |
| 4 | `a9_b1111_T1` | `{3}` only |

A third, fully hand-checkable route for `k = 1` (`E4`) reduces the condition to
`q(r) = 0` and `(r+1)q''(r) = 2(z-3)q'(r)`, i.e. to
`gcd(q, (y+1)q'' - 2(z-3)q') = 1`, and that gcd is `1` for every `z` in `[0,9]`.

**Sweeping `z = 0..9` is complete, not a heuristic window**: `z <= 9` is forced
by `Z | t^9`, which costs no cap and no cascade.

**Non-vacuity, honestly graded.** `E3`: the test returns FEASIBLE at `(4,3)` on
the genuine quartic. `E5` (= `X2`): the `k=1` machinery returns FEASIBLE at
exactly `z = 3` on `y^4+y^3-9y^2+7`. `E6` (**new**): the `k=2` machinery returns
FEASIBLE on a synthetic squarefree quartic with `q(-1) != 0`. That covers
`k = 1, 2, 4`. For `k = 3` no such control is possible: with `deg Q = 1`,
`deg Pi^2 = 6` and `gcd(Pi,t) = 1`, **no** `z` can work for **any** squarefree
quartic with `q(-1) != 0` (`E7` demonstrates this on three unrelated quartics).
So the `k = 3` kill is structural rather than arithmetic — which is a *stronger*
kill, but it means the `k=3` "infeasible" cannot be distinguished from a bug by
a positive control. It is however cross-checked by two independent
implementations and by the hand degree count in `SUB1_SPINE9.md` §5.

### 1.4 The degree caps — **CLEARED, and every one traced to source**

Read at runtime from `cascade_engine.SUB1.aux_caps = (d1,sigma,d2) = (9,12,6)`,
`e_cap = 15`, `full_system_bridge.STRIP_DEGCAP["sub1"] = {dm2:18, dm3:21,
dm4:24}` (`A9`), cross-checked against `FRONTIER_REBUILD.md` §1's sub1 row
(`d2 = 6`, `R = 18`, `S = 21`, `deg e cap = 15`) (`A10`), and
`caps_audit.py --quiet` re-run by this auditor: **70/70, exit 0**.

Stripped consequences, re-derived (`F1`–`F4`):

| `k` | `deg A` | `deg u` | `deg v` | `deg w` | cap on `deg F` | `deg F = 9+4k-z` |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 9 | 6 | 12 | 9 | 21 | `9-z` |
| 1 | 9 | 6 | 11 | 10 | 20 | `13-z` |
| 2 | 9 | 6 | 10 | 11 | 19 | `17-z` |
| 3 | 9 | 6 | 9 | 12 | 18 | `21-z` |
| 4 | 9 | 6 | 8 | 13 | **17** | `25-z` |

* `k = 4`: `z = 3` (from §1.3) gives `deg F = 22 > 17`. **Contradiction**, with
  no cascade input (`F6`). Second route: `z <= 6` gives `deg F >= 19 > 17` (`F7`).
* `k = 0`: `deg F = 9-z <= 7 <= 21`, so the ledger says nothing; the kill is the
  separate dichotomy of §1.5.
* `F8` confirms the ledger is **not** a blanket argument — for each of
  `k = 0,1,2,3` it permits some admissible `z`.

`A15`: `deg e = 9+k <= 13 <= 15`, consistent with the `e` cap; `deg d0 = sigma`
is never used, as the lane says.

### 1.5 `a9_b0000_T1`'s zero margin — **CONFIRMED as zero, in three places**

`F9` re-derives the elimination (residual 0):

```
gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z ,      deg(gamma*u) <= 6 .
```

`F10` re-poses the dichotomy as an explicit leading-term search over
`deg A in {A=0, 0..39}` and `z in [0,6]`: a contradiction is registered **only**
when a *unique* contributor attains the maximal degree (so no cancellation could
rescue it). It closes: either `2*deg A >= 8` uniquely attains `>= 8 > 6`, or
`deg A <= 3` and the degree-`7` term `2048*mu != 0` uniquely survives; and
`2*deg A = 7` is impossible.

Three independent zero-margin dependencies, each verified to switch the kill
off when weakened by exactly one:

* `F12` (= `X3`): `deg d2 <= 7` instead of `6` — kill **off**.
* `F13` (= `X1`): `z <= 7` instead of `6` — kill **off** (at `z = 7` the term
  `3*zeta*t^7` can cancel `2048*mu*t^3*q`'s leading term).
* `G8`: `v_t(d1) >= 3` is what makes `v_t(F) >= 3`; it comes from `3*v_t(h_1)`,
  so `v_t(h_1) >= 1` is load-bearing too.

This is the one cell that consumes the slice cascade, and the lane says so.

### 1.6 `v_t(h_7) >= 11` — the weakest integer, and what happened to it

Confirmed load-bearing with **zero margin** (`G6`): `v_t(h_7) = 10` gives
`z_max = 7` and the `k=0` kill switches off; `= 9` gives `8`. `v_t(h_6) >= 10`
is load-bearing too (`G7`: `v_t(h_6) = 9` gives `z_max = 8`).

**But both reduce to audited material** (`G10`, `G11`). The cascade's own
absorption form, recorded at `SLICE_OBSTRUCTION.md` §3 line 164 and listed as
imposed condition `(P<)` in `SLICE_OBSTRUCTION_AUDIT.md` §4, is

```
p_n = 2*h_n + q_n ,  q_n = sum_{i=1}^{n-1} h_i*h_{n-i} ,  t^(2n-2) | p_n  (n = 2..8)
=>  h_n = -q_n/2 + t^(2n-2)*g_n  ,  g_n integral
=>  v_t(h_n) >= min( min_i (v(h_i) + v(h_{n-i})) , 2n-2 ) .
```

From the **audited** `v_t(h_1..h_5) >= (1,3,5,7,9)`
(`slice_obstruction_audit.py` F1–F5, `k = 1..5`):

```
v_t(h_6) >= min(1+9, 3+7, 5+5, 10) = 10
v_t(h_7) >= min(1+10, 3+9, 5+7, 12) = 11
v_t(h_8) >= min(1+11, 3+10, 5+9, 7+7, 14) = 12
```

— exactly the committed profile. (`ALT_LEVEL12.md` L3.7 already runs this
argument for `h_6`; what this audit adds is that the *same* argument delivers
`h_7` and `h_8`, so the un-advanced rows are not an independent unaudited
input.) Consequently `z <= 6` holds on the audited cascade rows alone.

Residual exposure at this point: the `(P<)` structure itself and the audited
`h_1..h_5` rows — i.e. exactly what `slice_obstruction_audit.py` (56/56) covers.
This is materially better than `SUB1_SPINE9.md` §10.2's own self-assessment.

### 1.7 The `z` window and the shift convention

`G1` re-derives the `d3`-killing shift from **one** generalized-binomial model
(`D~_j = sum_{i>=j} C(i,j)*theta^(i-j)*D_i` for `j >= 0`,
`D~_{-1-r} = sum_j C(r,j)*(-theta)^j*D_{-1-r+j}` for `r >= 0`, `theta = -h_1/4`)
and finds it reproduces **all five** recorded forms simultaneously, residual 0:

```
d2 = h2 - (3/8)h1^2      d1 = h3 - (1/2)h2h1 + (1/8)h1^3      e = h5
R  = h6 + (1/4)h1h5      S  = h7 + (1/2)h1h6 + (1/16)h1^2h5
```

That five-for-one fit is itself strong evidence the convention is being applied
correctly rather than fitted per-object.

`G3`: under **every** reading and level — A (shifted, `[I3]`) at level 10 and 12,
B (unshifted) at level 10 and 12 — the window is `2 <= z <= 6` (B/level-12 gives
the sharper `z <= 5`). `G4` (= `X1b`): level 12 is inert under reading A.
`G5`: **the lane's correction to the brief is right** —
`v_t(R) >= min(v_t(h_6), v_t(h_1)+a_t) = min(10,10) = 10` at `a_t = 9`, so
`v_t(R) >= 11` does **not** hold and the closure does not depend on level 12.

### 1.8 `t^9 | R,S,T` — re-derived here, cap-free

The four non-cascade cells need this, so it must not smuggle in a cap. `G12`,
`G13` derive it from `e | S` (imported; proved three independent ways in
`DIVISOR_CONSEQUENCES.md` §2) and the K-syzygy alone:

* `S = e*sbar` in `K = 0` gives `3*e*R^2 = 2*Phi - e^3*(d2 + 3*sbar)`, so
  `9 + 2*v(R) >= min(30, 27) = 27` and `v_t(R) >= 9`;
* `G1 = 0` gives `T = -R*(sbar + d2) - d1*e/2`, so
  `v(T) >= min(v(d1)+9, 9+v(R), v(R)+v(S)) - 9 >= 9`.

No degree cap, no `d1 = 0`, branch-independent. `G9` gives a second route
through the cascade profile (`v_t(R) >= 10`, `v_t(S) >= 11`, `v_t(T) >= 12`).

### 1.9 The controls — all seven fire, and none is vacuous

| control | this audit | fires? |
|---|---|---|
| `X1` (cascade weakening) | `F13`, `G6` | yes — `z_max` 6 → 7 → 8, kill off |
| `X1b` (level 12 inert) | `G4` | yes — `(10,11)` and `(11,11)` both give 6 |
| `X2` (`k=1` is arithmetic) | `E5`, two routes | yes — synthetic quartic FEASIBLE at exactly `z=3` |
| `X2b` (test not vacuous) | `E3` | yes — `(k,z) = (4,3)` FEASIBLE on the real quartic |
| `X3` (`deg d2` cap) | `F12` | yes — cap 7 kills the kill |
| `X4` (`deg R` cap) | `F14` | yes — off once `deg R` cap reaches 23 |
| `X5` (division falsifiable) | `B6` | yes — `t` does not divide `kbox` |
| **new: `k=2` non-vacuity** | `E6` | yes — synthetic quartic FEASIBLE |
| **new: `k=3` structural** | `E7` | yes — infeasible on unrelated quartics too |

---

## 2. Per-cell table

| cell | `k` | verdict | killed by | consumes the cascade? | margin |
|---|---:|---|---|---|---|
| `a9_b0000_T1` | 0 | **EMPTY** | degree dichotomy on `gamma*u = mu*t^3*q - 6A^2 + 3*zeta*t^z` | **YES** (`z <= 6`) | **ZERO**, in three independent places (`deg d2 <= 6`; `z <= 6`; `v_t(d1) >= 3`) |
| `a9_b1000_T1` | 1 | **EMPTY** | `Pi^2`-support infeasible for **every** `z in [0,9]` | no | n/a — three independent decision procedures agree |
| `a9_b1100_T1` | 2 | **EMPTY** | `Pi^2`-support infeasible for **every** `z in [0,9]` | no | n/a — two procedures agree; positive control `E6` present |
| `a9_b1110_T1` | 3 | **EMPTY** | `Pi^2`-support infeasible for **every** `z in [0,9]` | no | structurally forced by `q(-1) != 0` + degree count; no positive control possible |
| `a9_b1111_T1` | 4 | **EMPTY** | support pins `z = 3`, then `deg F = 22 > 17` | no | wide: cap would have to reach `deg R <= 23` (`F14`) |

**Conditional on `a_t = 9`, the enumerated `f31` frontier is EMPTY.**

---

## 3. PROVED / CHECKED / INFERRED

**PROVED here** — exact polynomial identity, exact linear algebra over an
explicitly constructed number field, or exact saturated Gröbner over `Q`,
machine-checked in characteristic 0:

* the K-syzygy from `generators.json`, hence `e | 2*Phi` (`A2`, `A3`);
* **D1 (`rad(e) | t*q`) and D2 (`b_i in {0,1}`)** — proved here from the syzygy
  and the factorisation of `Phi`, not imported (`A13`, `A14`);
* the arithmetic of `q`: irreducible, squarefree, `q(-1) = 3315`, `lc = 2048`,
  Galois group `S4` (`A4`–`A7`);
* the `a = 9` reduction `G1,G2,G3,K -> g1,g2,g3,kbox` with boxed power
  `t^(30-3a) = t^3`, as a **free-ring** identity (`B1`–`B5`), and its
  falsifiability (`B6`);
* the marked-root chain `A(r) != 0`, `Pi | B`, `gcd(Z,Pi) = 1` (`D1`–`D5`);
* the elimination of `C` and of `d0`, and the cofactor identity `(*)` — also in
  its generic-`a` form (`C1`–`C7`);
* the boxed row in `(u,v)` and the substituted form (5) (`D8`, `D9`);
* the marked-support classification, **all `k`, all `z in [0,9]`, by two
  independent decision procedures plus a third for `k=1`** (`E1`–`E7`);
* the degree ledger and both degree kills (`F1`–`F11`);
* the `d3`-killing shift and all five inverse-shift relations, from one model
  (`G1`, `G2`);
* the `z` window `2 <= z <= 6` under all four readings (`G3`–`G5`);
* **`v_t(h_6) >= 10`, `v_t(h_7) >= 11`, `v_t(h_8) >= 12` from the audited
  `h_1..h_5` rows plus `(P<)`** (`G10`, `G11`);
* `t^9 | R,S,T` at `a = 9`, from `e | S` + the syzygy, cap-free and
  branch-independent (`G12`, `G13`);
* every control (`B6`, `E3`, `E5`, `E6`, `E7`, `F12`–`F14`, `G4`, `G6`, `G7`).

**CHECKED** — reproduced or read from an existing artifact, not re-proved here:

| statement | source | how this audit used it |
|---|---|---|
| `a_t = 9` — **the standing premise** | lower half audited (`slice_obstruction_audit.py` 56/56); **upper half `a_t <= 9` same-author, under separate audit** | assumed throughout |
| `generators.json` faithfully encodes the `(72,108)` D2 system | de-pickled from `t4_state.pkl`, sha256 recorded in the file's own provenance | read as the foundation |
| `Phi = -(1/6630)*(y+1)^30*q`, `q = 2048y^4-512y^3+320y^2-240y+195` | appears identically in >10 independent verifier scripts (`alt_regime_verify.py`, `audit_gb_kills.py`, `audit_convolution_kills.py`, …) and `AUDIT.md` §4 | typed as a constant, arithmetic facts re-proved |
| sub1 caps `6/9/12/15/18/21/24` | `cascade_engine.SUB1`, `full_system_bridge.STRIP_DEGCAP`, read at runtime; `caps_audit.py --quiet` **re-run: 70/70** | read, not re-derived |
| the cascade rows `v_t(h_1..h_5) >= (1,3,5,7,9)` | `slice_obstruction_audit.py` F1–F5 | consumed; `h_6..h_8` then **derived** (`G10`) |
| the `(P<)` absorption `t^(2n-2) | p_n`, `n = 2..8` | `SLICE_OBSTRUCTION.md` §3, `SLICE_OBSTRUCTION_AUDIT.md` §4 | consumed for `G10` |
| `e \| S` | `DIVISOR_CONSEQUENCES.md` §2 (three independent proofs) | consumed for `G12`/`G13` |
| the five cells are the whole `a_t = 9` residue | `FRONTIER_REBUILD.md` §2b/§3, 34-cell list | read |
| `b`-vectors are per-place labels over `Qbar`, `e` need not be in `Q[y]` | `FIELD_SCOPE_AUDIT.md`, `FIELD_SPLIT_AUDIT.md` | this is why `k = 1,2,3` are non-vacuous branches and why the ideal/extension posing is the right one |

**INFERRED** — nothing. No step here is a plausibility argument.

---

## 4. What this audit does NOT settle

1. **`a_t <= 9`.** The whole result hangs off it and it is same-author. This
   audit does not touch it. Five cells now hang off that single premise.
2. **The upstream boundary.** `[QQ1]` (the `alpha`-strip WLOG,
   `PROOF_INVENTORY.md` C3 at 2/4) is upstream of the cascade and of the caps;
   untouched here.
3. **The caps' own proofs.** `caps_audit.py` was re-run (70/70) but its
   derivations were not re-proved from the geometry. `deg d2 <= 6` is a
   zero-margin input to `a9_b0000_T1`.
4. **`generators.json`'s provenance** from `t4_state.pkl` was not re-derived.
5. **`k = 3` has no positive control** and cannot have one: see §1.3. It is
   backed by two independent implementations plus a hand degree count.
6. **Nothing is written to any ledger, DAG or state file.** This audit adds
   `spine9_audit.py` and this file and nothing else; `git status` shows no
   tracked file modified. Wiring the closure into the frontier remains a
   separate decision, and remains conditional on `a_t <= 9`.

---

## 5. Reproduce

```
cd d2_plane_72_108
python -u spine9_audit.py            # full report, 81 checks, ~12 s
python -u spine9_audit.py --quiet    # exit 0 iff every check passes
```

Pure sympy. Imports `cascade_engine` and `full_system_bridge` only to read the
certified caps at runtime; reads `generators.json` and
`slice_obstruction_stage.json`. Does not run `slice_obstruction_basis.py`
(which would rewrite `slice_obstruction_stage.json` under `--deep12`), does not
run `run_tests.sh`, does not import or execute `sub1_spine9.py` or `spine.py`,
uses no Singular, no msolve, no WSL, no subprocess, no modular arithmetic.
