/-
  PURE-KERNEL cross-checks of the `Cert.Poly` library.

  The headline theorem (`Cert.f37_certificate`) is discharged by `native_decide`,
  which compiles the verified `mul`/`add` to native code -- fast, but it extends
  the trusted base to the Lean compiler.  The lemmas below are proved by `decide`,
  i.e. the **Lean kernel itself** reduces `mergeF`, `scaleShift`, `mul`, `add`
  (all structurally recursive) with no compiler in the loop.  They witness that
  the library's arithmetic is kernel-reducible; the only reason the full
  certificate uses `native_decide` is the sheer term count (~5–6k-term cofactors),
  not any reliance on `native_decide`'s semantics.

  Honest scope: these are small instances.  They do NOT re-prove the full
  certificate in the kernel; see `FEASIBILITY.md` for the measured cost of a
  pure-`decide` run of the full identity.
-/
import Cert.Poly
import Cert.Data

namespace Cert

/-- Hand-checkable literal identity in one variable `x` (key = exponent, radix
    irrelevant for a single field):  `(1 + x)·(1 + x) = 1 + 2x + x²`.
    Keys: `1 = x⁰`, `2` mislabel? no -- here key `n` denotes `xⁿ`. -/
example :
    mul [(0, 1), (1, 1)] [(0, 1), (1, 1)] = [(0, 1), (1, 2), (2, 1)] := by
  decide

/-- `(x - 1)·(x + 1) = x² - 1`, exercising sign cancellation of the constant
    term in `mergeF`. -/
example :
    mul [(0, -1), (1, 1)] [(0, 1), (1, 1)] = [(0, -1), (2, 1)] := by
  decide

/-- The library's `mul` commutes on the ACTUAL (small) generators `G1`, `G2` of
    the certificate -- a genuine multivariate instance decided purely by the
    kernel. -/
theorem mul_comm_G1_G2 : mul G1 G2 = mul G2 G1 := by decide

/-- Likewise `add` commutes on real generator data. -/
theorem add_comm_G3_G4 : add G3 G4 = add G4 G3 := by decide

/-- A real three-generator sub-identity, kernel-decided:
    `G1·G2·G3` is associative in the two obvious bracketings. -/
theorem mul_assoc_G123 : mul (mul G1 G2) G3 = mul G1 (mul G2 G3) := by decide

end Cert
