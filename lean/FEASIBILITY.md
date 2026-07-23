# Feasibility: kernel `decide` vs `native_decide` for the full identity

**Question.** Can the full f37 certificate `certCombo = Df31` be discharged by
pure-kernel `decide` (no compiler in the trusted base), or is `native_decide`
required?

## The computation

`certCombo` assembles `ĉ1·G1 + ĉ2·G2 + ĉ3·G3 + ĉ4·G4` with the verified library:

- cofactors `ĉ1..ĉ4`: 5639 / 5228 / 5848 / 3418 terms;
- generators `Gᵢ`: 4 / 4 / 4 / 8 terms;
- each `mul ĉᵢ Gᵢ` folds over the *small* `Gᵢ`, so it costs `|ĉᵢ|·|Gᵢ|`
  key-shift operations (~20k–47k) plus `|Gᵢ|` linear merges;
- 3 further `add`s merge the four ~40k-term products; final result 102 terms.

Total elementary `Nat`/`Int` operations: on the order of a few hundred thousand.

## Findings

- **`native_decide`:** **SUCCEEDS in ~27s** (compile + run of the verified
  `Cert.Poly` functions on the full data). This is the headline proof. The kernel
  trusts the resulting `Bool` via `Lean.ofReduceBool`. Trusted base = Lean kernel
  + Lean compiler + this machine's runtime.

- **pure `decide`:** **did NOT complete** on this machine. With default limits it
  first hits `maxRecDepth` (the ~40k-deep list spines exceed the 512 default).
  With `maxRecDepth 1000000` and `maxHeartbeats 0` it ran for ~5.5 minutes and
  passed 6.7 GB of resident memory in the elaborator's `whnf` reduction without
  finishing, at which point the probe was stopped to protect the machine. The
  same library IS kernel-reducible in principle — `Cert/KernelCheck.lean` proves
  real multivariate `mul`/`add`/associativity identities on the small generators
  purely by `decide` in ~3s. The obstruction for the *full* identity is only
  scale: the kernel's `whnf` evaluator carries large `List` spines and GMP `Nat`
  keys term-by-term, orders of magnitude slower (and heavier) than compiled code.

## What a pure-kernel proof of the full identity would take

The self-contained route is already "finite arithmetic," so no new mathematics is
needed — only kernel throughput. Options, in increasing trust-reduction value:

1. **Keep `native_decide`** (current). Fast, one extra trust assumption
   (the compiler), standard in the Lean community for large `decide`-shaped facts.
2. **`decide` with a bignum-packed evaluation.** Reformulate `Poly` so the kernel
   does fewer, cheaper reductions (e.g. encode a whole polynomial as one GMP
   `Nat` via a Kronecker substitution and reduce the identity to a single `Nat`
   equality the kernel checks natively). This can bring the full identity into
   pure-kernel range; it is the natural next step if removing the compiler from
   the trusted base is required.
3. **A mathlib `MvPolynomial` + `ring`/`Ideal.mem` proof.** Removes `decide`
   entirely in favor of `ring`-normalization, at the cost of a mathlib dependency
   and `ring` on ~40k-term intermediates (itself heavy).

The `KernelCheck.lean` lemmas demonstrate option-2 machinery works in-kernel at
small scale; scaling it is engineering, not new proof content.
