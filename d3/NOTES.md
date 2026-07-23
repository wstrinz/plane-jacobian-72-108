# D3 — dimension-3 counterexample to the Jacobian Conjecture

## The statement
The Jacobian Conjecture JC(n): if `F: C^n -> C^n` is polynomial and
`det(JF)` is a nonzero constant, then `F` is invertible (in particular
injective). This map refutes it for `n = 3`:

```
F1 = (1 + xy)^3 z + y^2 (1 + xy)(4 + 3xy)
F2 = y + 3x(1 + xy)^2 z + 3x y^2 (4 + 3xy)
F3 = 2x - 3x^2 y - x^3 z
```

- `det(JF) = -2` identically (nonzero constant, so `F` is etale everywhere).
- `F` is **not injective**: the three distinct points
  `(0,0,-1/4)`, `(1,-3/2,13/2)`, `(-1,3/2,13/2)` all map to `(-1/4,0,0)`.
- The full fiber over `(-1/4,0,0)` is exactly those three points (verified by
  a Groebner/`solve` computation, no others).

A non-injective etale polynomial map is exactly what JC says cannot exist, so
this settles JC(3) in the negative.

## Provenance and status
Announced by L. Alpoge on 2026-07-19. The files here are our **independent
reconstruction and exact re-verification** of the published map — not a claim
of authorship. `verify.py` re-checks every stated property from scratch in
exact arithmetic and is wired into the repo's test runner.

## Relationship to the D2 work
D3 refutes JC in dimension 3. It does **not** resolve the plane case: the
degree-2 / plane Jacobian Conjecture JC(2) is a separate, still-open problem,
and the `d2_plane_72_108/` effort (raising the lowest possible counterexample
degree from 108 to 125) stands on its own regardless of D3.

## Reproduce
```
python3 verify.py      # exact checks; exits nonzero on any failure
```
