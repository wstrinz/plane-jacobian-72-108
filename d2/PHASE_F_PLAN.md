# Phase F — the global residue obstruction algebra

**Date:** 2026-07-22 (late). **Origin:** external review (GPT) of the campaign
state, assessed and adopted here; this document is the working plan.
**Status:** PLAN — nothing below is claimed proven.

## The organizing observation

The program's two strongest structural results, read together:

1. **No-jet-kill** (`RESIDUE_LEMMAS_DEPTH.md`): the t-place jet towers lift to
   every order — survivors are LOCALLY flexible at every finite place.
2. **Caps attained** (`ALT_COMBINED.md` J1, and tight witnesses across both
   windows): forced valuation minima meet the degree caps exactly — survivors
   are GLOBALLY almost determined by their divisors.

The missing layer is global polynomiality itself: the local germs at
S = {-1, r_1..r_4} (t and the four q-roots) must be jets of the SAME low-degree
polynomial on P^1. The machinery couples ORDERS across places; it does not yet
couple the LEADING COEFFICIENTS and derivatives that global reconstruction
forces.

## The central quantity: divisor defect

For each polynomial p in {d2, d1, sigma, e, g_l} under a surviving state,

    delta_p = deg p - sum_{s in S} v_s(p)  (>= 0 by the audited coupling).

- delta = 0: p = lambda * prod (y-s)^{v_s} — determined up to ONE scalar; every
  local leading residue becomes a determined algebraic number
  lambda * prod_{i != j} (r_j - r_i)^{v_i} * (r_j+1)^{v_t} — structured
  elements of the S4 splitting field, exactly the shape behind the two proven
  arithmetic kills C08/C20 (square classes vs Q(sqrt(17))).
- delta = 1, 2: linear/quadratic cofactor only; generic coefficient vectors
  disappear.
- The jet-gluing formulation: H^0(P^1, O(d)) -> (+)_s O_s/m^{k_s}; requested
  jet length > d+1 forces confluent-Vandermonde linear compatibility; = d+1
  determines the polynomial (elementary Riemann-Roch / Hermite interpolation /
  CRT; determinants controlled by resultants of place polynomials).

## The pipeline (per canonical survivor signature)

1. Divisor extraction (from the audited witnesses/states).
2. Defect compression: replace each polynomial by forced divisor x cofactor of
   degree delta_p.
3. Global jet gluing: Hermite/CRT relations among all local residues.
4. Galois-aware coordinates: work over the field fixed by the stabilizer of
   the multiplicity pattern (marked root -> quartic field; unordered pair ->
   sextic resolvent; symmetric -> Q).
5. Saturation by all exact-order/nonzero-lc conditions.
6. Toric quotient: invariant ratios (e.g. b^2/gamma^3, s^2/gamma^5,
   c/gamma^17) + the fixed-c fiber kept AS AN EQUATION (this is the sound
   resolution of the gauge question both the codex and batch lanes hit).
7. Initial-ideal test: monomial in in_w(I : lc^infty) kills the cone (tropical
   VARIETY, not prevariety — the current engines test generators one at a
   time); else sparse resultant / component decomposition.
8. Certificate output: small identity, norm contradiction, or unit-ideal
   certificate — cone-level lemmas, not row kills.

Target theorem shape: every surviving valuation signature either fails global
jet interpolation or has a saturated residue initial ideal containing a
monomial.

## Immediate work items (priority order)

- **F1 (first compute, cheap):** defect histograms over every surviving state
  in both windows + alternate regime. If a large fraction of the frontier has
  all defects in {0,1,2}, this route is confirmed. -> `phase_f_defects.py`.
- **F2 (pilot, alt regime):** a branch where d1 and sigma both attain forced
  minima (J1 guarantees existence): substitute exact factor shapes into the
  (D_t)/closing-h_0 congruences.
- **F3 (pilot, R9):** redo the z=0 ideal over Q(alpha) in Singular (WSL; pipe
  via stdin — /mnt/c broken), computing only enough S-polynomials to test for
  a monomial in the saturated initial ideal (not a full basis).
- **F4 (block lemma):** toric-quotient parametric treatment of the a10
  constant-E block (the gamma^17 signal) — one lemma replacing ~100 descents.
- **F5 (S-unit layer):** at defect-0 states, normalized ratios are S-units in
  C(y); derive 3-4 term relations and apply Mason/Brownawell-Masser-type
  bounds; classify vanishing subsums as branches.
- **F6 (optional recon):** numerical irreducible decomposition of stubborn
  residue ideals to decide WHAT proof to seek (empty / 0-dim / positive-dim).

## Deprioritized (with reasons)

- Deeper t-jets (no-jet-kill theorem explains why they cannot help).
- Real sign-splits (the a8 audit: no real-embeddability available; complex
  witness exists).
- Longer generic sympy Groebner timeouts (the R9 ordering result shows the
  wins come from structure, not budget).
- Raw kill counts as the progress metric — the metric is now: canonical
  global residue ideals / whole branches CLOSED.
- Compactification/ramification geometry (F7-class long shot; anything strong
  enough there threatens all of JC(2) — separate campaign if ever).

## Relation to running work

The overnight batch/R9 jobs (see OVERNIGHT_RUN_MANIFEST.md) remain useful:
their kills shrink the frontier and their UNRESOLVED census is the input to
F1. But the next BUILD cycle is F1-F4, not another same-author forcing round.
