/-
  A tiny self-contained sparse multivariate polynomial library over `Int`.
  No `mathlib`, no external dependencies -- everything here reduces in the Lean
  kernel, so the certificate can be discharged by `decide` (and, for speed, by
  `native_decide`).

  REPRESENTATION.  A polynomial is a `List (Nat × Int)` of `(key, coeff)` pairs,
  kept strictly increasing in `key` with every `coeff` nonzero.  A `key` is a
  monomial's exponent vector packed in a fixed radix `R` (chosen by the generator
  to exceed every exponent that can occur):  key = Σ_i e_i · R^i.  Because no
  field ever reaches `R`, adding two monomials' exponent vectors is exactly
  ordinary `Nat` addition of their keys -- which is what makes multiplication a
  one-line `key ↦ key + k` shift.

  All definitions are *structurally* recursive (the merge uses an explicit fuel
  argument), so the kernel can evaluate them for `decide`.
-/

namespace Cert

/-- Sparse polynomial: `(packed-exponent key, coefficient)` pairs, kept strictly
    sorted by key with nonzero coefficients. -/
abbrev Poly := List (Nat × Int)

/-- Merge two key-sorted polynomials, adding coefficients on equal keys and
    dropping any coefficient that cancels to `0`.  `fuel`-driven to be
    structurally recursive (so it reduces in the kernel). -/
def mergeF : Nat → Poly → Poly → Poly
  | 0, _, _ => []
  | _, [], q => q
  | _, (k1, c1) :: p, [] => (k1, c1) :: p
  | fuel + 1, (k1, c1) :: p, (k2, c2) :: q =>
      if k1 < k2 then
        (k1, c1) :: mergeF fuel p ((k2, c2) :: q)
      else if k2 < k1 then
        (k2, c2) :: mergeF fuel ((k1, c1) :: p) q
      else
        let c := c1 + c2
        if c == 0 then mergeF fuel p q
        else (k1, c) :: mergeF fuel p q

/-- Sum of two key-sorted polynomials.  Fuel `= |p| + |q| + 1` always suffices
    (each recursive step drops at least one input element). -/
def add (p q : Poly) : Poly := mergeF (p.length + q.length + 1) p q

/-- Multiply a polynomial by a single monomial `(k, c)`.  Shifting every key by a
    constant `k` is monotone, so the result stays strictly sorted. -/
def scaleShift (k : Nat) (c : Int) (p : Poly) : Poly :=
  p.map (fun kc => (kc.1 + k, kc.2 * c))

/-- Polynomial product.  Folds over `q`, accumulating `Σ_t (p · t)` by `add`.
    Cost is `|p|·|q|` shifts plus linear merges, so **pass the smaller operand as
    `q`** (here always one of the tiny generators `Gᵢ`). -/
def mul (p q : Poly) : Poly :=
  q.foldl (fun acc kc => add acc (scaleShift kc.1 kc.2 p)) []

end Cert
