# Independent cross-author audit — f31/subcase-(1) alternate regime

Date: 2026-07-22
Auditor: independent (different author from `ALT_REGIME.md` / `ALT_REGIME_L2.md`)
Verifier: `alt_regime_audit_verify.py` (shares no code with the audited scripts
beyond the `f31_graded.txt` regex parse)
Scope: the weakest-provenance chain in the program (PROOF_INVENTORY gap G3,
C33–C34), claiming 25 of 52 alternate-regime branches (a ∈ [11,15]) killed.

## Bottom line

**Every audited claim is CONFIRMED.** My independent re-derivation of the
flipped reduction, the descending telescope, the q-terminal identities, the
first-level parity lemmas, and the deep h₅ cone reproduces the claimed
**25-branch kill list exactly (25/25 agreement: 13 T1 + 12 T2)**, and the
**27-branch residual list is honest** — no additional branch is killable with
the authors' own tools (nor with the stronger T3 σ-locus closure). No soundness
error found. Two harmless bookkeeping notes are recorded in §Notes.

The mathematics rests on the trusted ground truth only: the graded identity
`F = Σ_f Φ~^f e^(21-3f) h_f` (t5_multiplace_verify.py checks 2/5/6), the raw
`h_f` from `f31_graded.txt`, the sub1 caps `deg d1≤9, deg σ≤12, deg e≤15-a`
(sub1_cascade_verify.py), and the T3 σ-locus theorem (proven for deg e ≤ 15).

## Per-claim verdicts

### C33-a. Flipped reduction `F = t²¹⁰ G'` — CONFIRMED

t-order of the f-term is `30f + a(21-3f) = 21a + f·v`, `v = 30-3a`. For a ≥ 11,
`v ≤ -3`, so this is *strictly* decreasing in f; the minimum is at f=7 and
equals `21a + 7(30-3a) = 210`. The residual exponent of the f-term is exactly
`(7-f)w`, `w = 3a-30 > 0`, hence `G' = Σ_f T^(7-f) u^f E^(21-3f) h_f`, `T = t^w`,
is a genuine polynomial and `F = t²¹⁰ G'`. I verified this **exactly on my own
seeded random subcase-(1) windows at a=13 AND a=14** (the audited script only
checked a=12), for all eight f-terms and the full sum. `210` is the *exact*
t-order because f=7 is the unique minimizer of the explicit t-power (power 0).
[audit verifier parts 1–2]

### C33-b. Descending cascade / telescope — CONFIRMED

The recursion `T r_{f-1} = E^(3(7-f)) h_f + u r_f` (`r₇ = r₋₁ = 0`) is
*equivalent* to `G' = 0`. I proved this cleanly by solving the r's top-down and
obtaining the closed identity

```
G' = T⁷ (E²¹ h₀ + u r₀),
```

so `G' = 0 ⟺` the bottom equation `E²¹ h₀ + u r₀ = 0` (T ≠ 0). The bottom-up
auxiliaries reproduce the terminal law `E³ g₇ + u⁷ h₇ = G'`. The anchor `T | h₇`
(i.e. `2 v_t(d1) ≥ w`) is forced because h₇'s term carries the least explicit
t-power. [audit verifier part 3]

### C33-c. q-place invariance and terminal identities — CONFIRMED

At a q-root p, `t` (hence `T`) is a unit (q(-1) = 3315 ≠ 0), and
`f + b(21-3f) = 7 + (7-f)(3b-1)`, so dividing G' by p⁷ gives the *same* chain
with spacing `s = 3b-1`. The q-terminal identity is not "imported" — it follows
from `G' = 0` itself: `E³ g₇ = -u⁷ h₇` exactly, so taking p-valuations
(`v_p(E)=b, v_p(u)=1, v_p(h₇)=2 v_p(d1)`) gives

```
3b + v_p(g₇) = 7 + 2 v_p(d1),      and analogously   3b + v_p(g₆) = 6 + 2 v_p(σ).
```

This is genuinely regime-independent because only the q-local data and `G'=0`
enter; the (invalid) upward t-cascade is never used. At t the same identity
reads `v_t(g₇) = 2 v_t(d1)` (u, E are t-units). **The claim that the standard
terminal identity holds verbatim in the flipped regime is correct.**
[audit verifier part 5]

### C33-d. First-level parity lemmas — CONFIRMED

I re-implemented the parity argument from the raw h₆/h₅ monomials (source
re-derived in the (d1, d2, σ, e) basis: `h₆ = 14336 d1²d2 + 8192 d1 e - 3072 σ²`,
`h₅ = -12288 d1²d2² + 32256 d1²σ + 18432 d1 d2 e - 9216 d2 σ² + 2048 e²`):

- **T1:** odd s ⟹ `v(r₆) = 2x-s` is odd, unmatchable by any h₆ term
  (`2x+k`, `x+m` both `> 2x-s` for x<s; `2z` even), so `x ≥ s`; even s ⟹
  anchor `x ≥ s/2`. Gives min v(d1): t = {11:3,12:3,13:9,14:6,15:15},
  q = {1:1,2:5,3:4,4:11}.
- **T2 (d1=0):** `2m ≥ s ⟹ z ≥ s`; `2m < s` needs the e² term (order 2m) to
  match `2z-s`, impossible for odd s, giving `z=(s+2m)/2` for even s. Gives
  min v(σ): t = {11:3,12:6,13:9,14:12,15:15}, q = {1:2, 2:impossible, 3:7,
  4:impossible}.

All match the audited `T1_T/T1_Q/T2_T/T2_Q`. [audit verifier part 6]

### C34. Deep h₅ cone and the 6 new kills — CONFIRMED

I rebuilt the h₇→h₆→h₅ T1 tropical cone from scratch (ties allowed to cancel to
any depth — a safe over-approximation; empty cone = rigorous kill), using
**EXISTS-semantics** over the free valuations including k = v(d2). The tightened
per-place min v(d1) is `t = (a11:5, a12:3, a13:∅, a14:6)`, `q = (b1:1, b2:7,
b3:4, b4:∅)`, matching `ALT_REGIME_L2`'s T1F/T1Z tables. The 6 new kills follow
from `Σ min v(d1) > 9` or an empty local cone (e.g. a13 t-cone empty, b4 empty).
[audit verifier part 7]

### Kill list (all 52 branches) — CONFIRMED 25/25

A fully independent tropical + degree-budget computation over all 52 branches
(own cones, own convolution against `deg d1 ≤ 9`, `deg σ ≤ 12`, EXISTS-k)
reproduces the claimed kill list **exactly**:

| branch | independent kills | claimed | agree |
|---|---:|---:|:--:|
| T1 | 13 | 13 (7 L1 + 6 L2) | ✓ |
| T2 | 12 | 12 (all L1) | ✓ |
| total | **25** | **25** | ✓ |

Every T2 kill is either a q-parity impossibility (b ∈ {2,4}) or a forced
σ-order > 12 (⇒ σ = 0, closed by the T3 σ-locus theorem, which is proven for
deg e ≤ 15 and therefore covers this regime). Every T1 kill is a forced
d1-order > 9 or an empty local cone. [audit verifier part 8]

### Residual list (27 branches) — HONEST

27 open = **13 open T1 + 14 open T2**, matching `ALT_REGIME_L2` §3. Two
completeness probes with the authors' own tools:

1. The g₇/g₆ degree-order bound (their L1 mechanism) kills nothing beyond the
   cone kills — every g-bound kill is already in the 25.
2. **No** residual T1 branch survives only through the σ=0 route; each retains a
   genuine σ≠0 local solution in the over-approximated cone. So even invoking
   the stronger T3 σ-locus closure on σ=0 would not shrink the residual. The
   residual is not overclaimed. [audit verifier part 9]

## Notes (harmless; no action required)

- **N1 (k-quantifier).** `alt_regime_l2_verify.py`'s `project()` classifies a
  local (x,z) by `all(t1_local(...) for k)`. The mathematically-motivated
  over-approximation is `any` over k (the counterexample chooses d2, i.e. one
  good k). Using `all` is a *subset* of the correct cone and could in principle
  produce a false "empty ⇒ kill". I checked explicitly: **for every active
  place the ALL-k and ANY-k cones are identical**, so the audited kills are
  unaffected. Recommend switching to `any` for conceptual correctness, but it is
  not a bug in the current results.

- **N2 (fully-dead strata count).** `ALT_REGIME.md` reports "6 strata die in
  both branches"; that is the **pre-L2** figure. After the 6 new L2 T1 kills
  (four of which complete strata already T2-dead in L1), **10** strata are dead
  in both branches. `ALT_REGIME_L2.md` does not restate the "6", so there is no
  internal contradiction — only a stale cross-reference if the L1 number is read
  post-L2.

## Method

I verified the audited scripts also pass (`alt_regime_verify.py`,
`alt_regime_l2_verify.py` both green), but every verdict above rests on my own
derivations in `alt_regime_audit_verify.py`, which imports none of their logic.
Ground-truth inputs used: `f31_graded.txt` (raw, regex-parsed), the graded
identity and reduction from `t5_multiplace_verify.py`, the caps from
`sub1_cascade_verify.py`, and `split_place_ledger_sub1.json` (26 strata).
