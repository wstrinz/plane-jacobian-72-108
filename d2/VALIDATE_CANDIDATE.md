# validate_candidate.py — usage & verdict glossary

An exact, standalone checker for a **claimed plane Jacobian counterexample**
(JC(2)): two polynomials `P, Q ∈ K[x,y]` alleged to have constant Jacobian yet
no polynomial inverse. It runs four staged checks, cheapest first, and reports a
verdict for each plus an overall verdict.

**Dependencies:** Python 3 + `sympy` only (no scipy; the convex hull is
implemented in-file). Run natively on Windows.

## Input formats

Pass one file (extension selects the parser):

**`.py`** — a plain Python module that assigns `P` and `Q` as sympy expressions.
The symbols `x`, `y` (and `sympy`, `symbols`) are pre-injected, so:

```python
# candidate.py
P = x + y**2
Q = y + (x + y**2)**2
```

**`.json`** — an object with string fields `"P"` and `"Q"`, each parsed by
`sympy.sympify` over the symbols `x, y`:

```json
{ "P": "x**2*y + x", "Q": "y**2" }
```

Both must be polynomials in `x, y` only; any other free symbol is rejected.

```
python validate_candidate.py candidate.py
python validate_candidate.py candidate.json
python validate_candidate.py --self-test
```

Exit code is `1` when the claim is **REFUTED**, else `0` (so a script can gate
on it). The full human-readable report always goes to stdout.

## What each stage checks

1. **BASIC** — computes `[P,Q] := P_x Q_y − P_y Q_x` exactly and checks it is a
   *nonzero constant*. Also reports `deg P, deg Q` and whether both exceed 1
   (a crude "non-invertible-looking" gate). If `[P,Q]` is non-constant or zero,
   the pair is refuted here.
2. **BOUND CHECK** — if `max(deg P, deg Q) < 125`, GGHV22 (arXiv:2204.14178,
   Prop 4.3) forces the shape `(deg P, deg Q) = (72,108)` up to symmetry for any
   genuine counterexample. A different shape is flagged (see caveat below).
3. **SHAPE CHECK** — computes the Newton-polygon corners (convex-hull vertices of
   the exponent set) of `P` and `Q` and compares them to the two Prop 4.3 subcase
   polygons. Those corner sets are **read at runtime from
   `T3_WINDOW_AUDIT.md` section 1** (which transcribes them verbatim from the
   paper, lines 1000–1007) — they are not hand-copied into the script. The
   reference describes the *reduced* `[P,Q]=x²` representative of total degree
   `(24,36)`; a raw `(72,108)` candidate is `3×` that and must be put through the
   paper's reduction (inversion morphism + normalization chain, the open **T6**
   debt) before the corners can match. So a non-match at full scale is **not** a
   refutation — it is reported as such.
4. **NEXT** — if the claim survives 1–3, the tool names (does **not** run) the
   deeper framework checks: the **f31** master-elimination necessary condition
   (`STATE.md` items 5–6) and the **frontier degree-states**
   (`phase_d_states_sub2.json` / `phase_d_states_sub1.json`, overview in
   `FRONTIER.md`).

## Verdict glossary

| Verdict | Meaning |
|---|---|
| **REFUTED** | A necessary condition (or a cited published theorem) is violated — non-constant/zero Jacobian, or a `<125` degree shape that isn't `(72,108)`. The pair is not a valid counterexample *as stated*. Exit 1. |
| **OUT OF SCOPE** | `[P,Q]` may be a nonzero constant, but a degree is `≤ 1`, so the map is (or looks like) an automorphism — not a counterexample claim at all. Exit 0. |
| **SURVIVES** | Passes every *implemented* necessary condition. This is **not** a certificate — see scope below. Exit 0. |

### The invertibility caveat (important)

The "non-invertible" gate is only the heuristic `deg P > 1 and deg Q > 1`. A
genuine plane automorphism can have *both* components of degree `> 1` (e.g.
`P = x + y²`, `Q = y + (x+y²)²`), pass that gate, and then trip the step-2 FLAG.
The tool therefore states the FLAG as a disjunction: *either* the pair is a real
non-invertible counterexample (in which case a `<125` off-shape degree is a hard
contradiction with GGHV22) *or* it is simply an automorphism and thus not a
counterexample. Deciding invertibility exactly is out of scope.

## Scope — honest statement

This tool validates **necessary** conditions only. It can **refute** a claim
(non-constant Jacobian, wrong degree shape, wrong Newton polygon). It **cannot
certify** that a surviving pair is a genuine counterexample: sufficiency would
require the full program — the f31 window-infeasibility argument, the frontier
degree-state exhaustion, and the still-open T6 reduction — none of which is
implemented here. "SURVIVES" means "not yet refuted by the cheap checks," nothing
more.

## Self-test

`python validate_candidate.py --self-test` runs:

- **(a)** `P=x, Q=y` — Jacobian `= 1` (constant) but degrees `≤ 1`: reported as
  **OUT OF SCOPE** (the identity/automorphism, not a counterexample).
- **(b)** `P=x³+y, Q=x+y²` — `[P,Q] = 6x²y − 1` is non-constant: **REFUTED** at
  step 1.
- a **shape-matcher unit demo** — synthetic polynomials whose Newton polygons
  equal Prop 4.3 subcase (2), calling the step-3 matcher directly to confirm the
  parser + hull reproduce the reference corners exactly. (They are not a Jacobian
  pair, so the full pipeline would stop them at step 1; hence the direct call.)

The self-test carries assertions, so it doubles as a regression check.
