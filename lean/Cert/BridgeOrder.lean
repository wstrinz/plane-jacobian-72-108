import Cert.Poly

namespace Cert.BridgeOrder

/-- Negation in the tiny sparse integer-polynomial checker. -/
def neg (p : Poly) : Poly := scaleShift 0 (-1) p

/-- Formal derivative for the univariate encoding used by this certificate. -/
def deriv : Poly → Poly
  | [] => []
  | (0, _) :: p => deriv p
  | (k + 1, c) :: p => (k, Int.ofNat (k + 1) * c) :: deriv p

/-- Power computed through the certificate library's normalized product. -/
def pow : Poly → Nat → Poly
  | _, 0 => [(0, 1)]
  | p, n + 1 => mul (pow p n) p

/-- `c = y^3(y^5+1)`. -/
def c144 : Poly := [(3, 1), (8, 1)]

/-- `15 f = -y^4(y^5+1)^2`; clearing 15 does not change its y-order. -/
def f144Cleared : Poly := [(4, -1), (9, -2), (14, -1)]

/-- The denominator-cleared forcing ODE:
    `12 c (15f)' - 21 c' (15f) = 15 c^2`. -/
theorem forcing_ode_144 :
    add (scaleShift 0 12 (mul c144 (deriv f144Cleared)))
        (scaleShift 0 (-21) (mul (deriv c144) f144Cleared))
      = scaleShift 0 15 (mul c144 c144) := by
  decide

/-- Order at y=0 for a normalized nonzero sparse polynomial. -/
def ordY : Poly → Nat
  | [] => 0
  | (k, _) :: _ => k

/-- Slice weight `M=25` from the independently built tower. -/
def sliceM144 : Nat := 25

/-- The affine D-transform clearing sum `3M-b`. -/
def clearingExponent144 : Nat := 3 * sliceM144 - 4

/-- Removing the four powers already absorbed into `15f` leaves `N=67`. -/
def postClearingExponent144 : Nat := clearingExponent144 - 4

/-- The independently constructed cleared-slice candidate `f*c^N`. -/
def clearedPhi144 : Poly :=
  mul f144Cleared (pow c144 postClearingExponent144)

set_option maxRecDepth 100000 in
/-- Exact order of the cleared-slice candidate. This does not assert the
    remaining geometric identification with the global Phi. -/
theorem ordY_clearedPhi144 : ordY clearedPhi144 = 205 := by
  decide

/-- Tower arithmetic: clear=3M-b=71 and post-clearing N=71-b=67. -/
theorem tower_arithmetic_144 :
    clearingExponent144 = 71 ∧ postClearingExponent144 = 67 := by
  decide

/-- The bridge side agrees only after the independent construction. -/
theorem bridge_prediction_144 : 3 * 3 * 25 - 20 = 205 := by
  decide

end Cert.BridgeOrder

#print axioms Cert.BridgeOrder.forcing_ode_144
#print axioms Cert.BridgeOrder.ordY_clearedPhi144
#print axioms Cert.BridgeOrder.bridge_prediction_144