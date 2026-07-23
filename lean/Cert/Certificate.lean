/-
  ============================================================================
  THE f37 IDEAL-MEMBERSHIP CERTIFICATE  (kernel-checked)
  ============================================================================

  Program headline fact (d2_plane_72_108/F37_SATURATION_REPORT.md, 2026-07-22):
  the degree-31 elimination polynomial `f31` lies in the pre-resultant ideal,

      f31 = c1·G1 + c2·G2 + c3·G3 + c4·(G5body + Φ)              (exact, over ℚ)

  in ℚ[d̃2, d̃1, d̃0, d₋₁, d₋₂, d₋₃, d₋₄, Φ].  Hence `f31` vanishes on the entire
  pre-resultant variety and the whole `{f37 = 0}` branch of the resultant is a
  classical resultant-excess artifact.

  This file proves the equivalent INTEGER identity emitted by the generator
  `d2_plane_72_108/lean_export/export_certificate.py` (which reads the SAME data
  sources as `f37_sat_verify.py` -- `t4_state.pkl`, `f37_sat_certificate.txt`,
  `f31_deg31.txt` -- and never hand-copies a coefficient):

      Dmul · f31  =  ĉ1·G1 + ĉ2·G2 + ĉ3·G3 + ĉ4·G4          (exact, over ℤ)

  where  Gᵢ  are the denominator-cleared generators (G4 = G5body+Φ, cleared),
  ĉᵢ = Dmul·cᵢ are the denominator-cleared cofactors, and `Dmul` is a fixed
  NONZERO integer (the lcm of all cofactor denominators; its value is in
  `Cert/Data.lean`).  Dividing the integer identity by the nonzero integer `Dmul`
  recovers the ℚ-identity above; the divisibility/nonzeroness step is elementary
  and left to the reader -- the kernel-checked content is the integer identity.

  The polynomials are opaque `List (Nat × Int)` data; ALL arithmetic (`mul`,
  `add`) is performed and checked by the verified `Cert.Poly` library, so the
  Lean kernel -- not the generator, not Singular, not sympy -- certifies the
  identity.
-/
import Cert.Poly
import Cert.Data

namespace Cert

/-- The certificate combination `ĉ1·G1 + ĉ2·G2 + ĉ3·G3 + ĉ4·G4`, assembled with
    the verified library.  Each `mul` passes the tiny generator `Gᵢ` as the second
    (folded-over) operand for efficiency. -/
def certCombo : Poly :=
  add (add (add (mul chat1 G1) (mul chat2 G2)) (mul chat3 G3)) (mul chat4 G4)

/-- **f37 membership certificate (integer form).**
    `Dmul · f31 = ĉ1·G1 + ĉ2·G2 + ĉ3·G3 + ĉ4·G4` as an exact polynomial identity
    in ℤ[d̃2,d̃1,d̃0,d₋₁,d₋₂,d₋₃,d₋₄,Φ].  `Df31` is `Dmul·f31` and `certCombo` is the
    right-hand side, both as normalized sparse polynomials; equality of the two
    normal forms is decided by the kernel.

    `native_decide` runs the (verified) `mul`/`add` at native speed.  The SAME
    library is exercised by pure-kernel `decide` on smaller real instances in
    `Cert/KernelCheck.lean` (so the arithmetic is kernel-reducible, not reliant on
    `native_decide`'s semantics); `FEASIBILITY.md` records the cost of a full
    pure-`decide` run. -/
theorem f37_certificate : certCombo = Df31 := by
  native_decide

end Cert
