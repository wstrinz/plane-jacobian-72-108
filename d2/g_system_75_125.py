#!/usr/bin/env python3
"""g_system_75_125.py  (NEW; read-only over all existing artifacts)

THE TRANSFER TEST, phase 2 -- BUILD the (75,125) analogue of the D-transform
G-SYSTEM: the pre-resultant window presentation that made (72,108) tractable
(FULL_SYSTEM_BRIDGE.md sec.1, regenerate_system.py).

Recipe (parametric, derived from the (72,108) construction and phase-1's tower
C_SERIES_75_125.md, which fixed forcing window = S^b, linear window = S^a):

    S = sum_m d_m u^(t-m),   d_t = 1,  d_{t-1} = 0 (x-shift).
    Linear window  = S^a : the slices La(k) := [u^(a t + k)] S^a are LINEAR in the
        deepest new coefficient dm_{(a-1)t+k} (coefficient a from d_t^(a-1));
        solved sequentially to eliminate the deep window unknowns.
    Forcing window = S^b : the slices Lb(j) := [u^(b t + j)] S^b, after the linear
        substitutions, are the generators G_j (weight-homogeneous, degree b).
        Phi = f*C^N is the (D~^b)_{-jphi} slice at u^M, M = b t + jphi,
        jphi = -s = a t - kappa - 1 (phase-1 slice-sum invariant).

For (72,108): a,b,t = 2,3,4 -> linear S^2, forcing S^3, cubic generators
G1,G2,G3,(G5body+Phi) -- exactly regenerate_system.py / FULL_SYSTEM_BRIDGE.md.
For (75,125): a,b,t = 3,5,4 -> linear S^3, forcing S^5, QUINTIC generators.

*** REPAIRED 2026-07-26 (PASSPORT_75_125_REPAIR.md).  t was 5 and q was 2, read
off GGV5's final chain corner (7\\5,2) via the dictionary (t,q)=(l_final,b_final).
That dictionary holds only on the retraction shape b0 = l(a0-1), which (5,20)
fails; the chart exponent is ceil(20/5) = 4.  The (75,125) build is now
build_gsystem(3, 5, 4, 1, 80) and its inputs come through
polygon_reduction.corner_chart_data(), which RAISES rather than returning the
wrong pair.  Downstream: 8 generators (was 10), 7 spares dm2..dm8 (was 9,
dm2..dm10), M = 29 (was 36), W_step = 80/29 (was 201/36 = 67/12).  Both the
"structure transfers" conclusion and the "physical weight normalisation does NOT"
conclusion SURVIVE -- but every number in them changed, and the specific
period-12 / divisor-lattice content of the old story is refuted: 29 is prime. ***

Two structural facts this script establishes and the verifier re-checks:
  (1) The recipe transfers as an algebraic construction: the (75,125) G-system
      EXISTS, is weight-homogeneous under the intrinsic u-grading, and has a
      well-defined spare inventory.
  (2) The PHYSICAL (y-valuation) weight normalisation does NOT transfer cleanly:
      W_step := ord_y(Phi)/M is 12 (integer) at (72,108) but 80/29 (NON-integral)
      at (75,125).  The clean integer y^W stripping of the (72,108) bridge does
      not carry over; the grading is rational.  NOTE the repaired shape of this
      obstruction: because C = y is a monomial, ord_y(Phi) = deg_y(Phi) = 80, so
      the ord-slope and deg-slope COINCIDE (both 80/29) instead of being the two
      distinct slopes 67/12 and 14 of the superseded model.  q_window = 29 is
      PRIME, so the "divisor lattice of the period" story attached to 12 is gone.

Emits g_system_75_125.json (canonical, documented variable order) for the case
compiler to consume.  Independent checker: g_system_75_125_verify.py.
Exact sympy throughout.  Run end to end (~3-6 min for the full S^5 build).
"""
import json
import sys
import time

import sympy as sp

Phi = sp.Symbol("Phi")


# ---------------------------------------------------------------------------
# The parametric construction (works for any standard length-1 corner).
# ---------------------------------------------------------------------------
def build_gsystem(a, b, t, q, ordPhi, Nmax_override=None, jset=None, verbose=False):
    """Build the D-transform G-system for one corner.

    Returns the generators, the linear-elimination substitutions, the spare
    inventory, the intrinsic (u-grading) weights, and the physical-weight data.
    kappa = t-2 is structural (PHI_CORNER4.md); q = ord_y(C) (leading poly).
    ordPhi = ord_y(Phi_full) fixes the physical weight normalisation.
    """
    kappa = t - 2
    e = b - a + 1
    s = kappa + 1 - a * t
    jphi = -s                    # forcing-slice offset (= a t - kappa - 1)
    M = b * t + jphi             # forcing u-power slice (Phi lives here)
    deep = (b - 1) * t + jphi    # deepest window unknown reached by a generator

    # series coefficients indexed by u-power p:
    #   p=0 -> d_t = 1 ; p=1 -> d_{t-1} = 0 (shift) ;
    #   2<=p<=t -> d_{t-p} (shallow above-line) ; p>=t+1 -> dm_{p-t} = d_{-(p-t)}
    coeff = {0: sp.Integer(1), 1: sp.Integer(0)}
    dh, dm = {}, {}
    for p in range(2, t + 1):
        sym = sp.Symbol(f"d{t - p}")
        dh[t - p] = sym
        coeff[p] = sym
    for kk in range(1, deep + 1):
        sym = sp.Symbol(f"dm{kk}")
        dm[kk] = sym
        coeff[t + kk] = sym

    Nmax = M if Nmax_override is None else Nmax_override

    def conv(A, B):
        out = {}
        for pa, va in A.items():
            if pa > Nmax:
                continue
            for pb, vb in B.items():
                pp = pa + pb
                if pp > Nmax:
                    continue
                out[pp] = out.get(pp, sp.Integer(0)) + va * vb
        return out

    S = {p: c for p, c in coeff.items() if p <= Nmax}
    t0 = time.time()
    Spow = {1: S, 2: conv(S, S)}
    for n in range(3, b + 1):
        Spow[n] = conv(Spow[n - 1], S)
    Sa, Sb = Spow[a], Spow[b]
    if verbose:
        print(f"    convolution {time.time() - t0:.1f}s")

    def getc(P, p):
        return sp.expand(P.get(p, sp.Integer(0)))

    # ---- linear window: eliminate the deep window unknowns ----
    Klin = (b - a) * t + jphi              # deepest linear slice needed
    skiplin = (b - a) * t + (jphi - 1)     # slice tied to the skipped generator
    sub = {}
    for k in range(1, Klin + 1):
        if k == skiplin:
            continue
        slice_p = a * t + k
        if slice_p > Nmax:
            continue
        target = dm[(a - 1) * t + k]
        eq = sp.expand(getc(Sa, slice_p).xreplace(sub))
        c0 = eq.coeff(target, 1)           # linear coefficient (constant = a)
        assert c0 != 0 and not c0.free_symbols, (k, c0)
        sub[target] = sp.expand(-(eq - c0 * target) / c0)
    if verbose:
        print(f"    linear window {time.time() - t0:.1f}s ({len(sub)} eliminated)")

    # ---- forcing window: the generators ----
    Gs = {}
    jrange = range(1, jphi + 1) if jset is None else jset
    for j in jrange:
        if j == jphi - 1:                  # skipped generator (tied to skiplin)
            continue
        slice_p = b * t + j
        if slice_p > Nmax:
            continue
        raw = sp.expand(getc(Sb, slice_p).xreplace(sub))
        if j == jphi:
            raw = raw + Phi                # Phi enters the deepest forcing slice
        Gs[j] = raw
    if verbose:
        print(f"    forcing window {time.time() - t0:.1f}s ({len(Gs)} generators)")

    # ---- intrinsic (u-grading) weight of every symbol: w(d_m) = t - m ----
    def uweight(sym):
        nm = str(sym)
        if nm == "Phi":
            return M
        if nm.startswith("dm"):
            return t + int(nm[2:])
        return t - int(nm[1:])             # d{idx}: uweight = t - idx

    # spare inventory + state
    state = [dh[t - p] for p in range(2, t + 1)] + [dm[1]]   # d_{t-2}..d_0, e=dm1
    allsym = set()
    for g in Gs.values():
        allsym |= g.free_symbols
    allsym.discard(Phi)
    spares = sorted(allsym - set(state), key=lambda s: int(str(s)[2:]))

    # homogeneity: each generator's monomials share one u-weight = b t + j
    homog = {}
    for j, g in Gs.items():
        gg = g - (Phi if j == jphi else 0)
        P = sp.Poly(gg, *sorted(g.free_symbols - {Phi}, key=str))
        ws = {sum(uweight(v) * ex for v, ex in zip(P.gens, mon)) for mon, _ in P.terms()}
        homog[j] = ws

    # physical weight normalisation
    W_step = sp.Rational(ordPhi, M)

    return dict(a=a, b=b, t=t, q=q, kappa=kappa, e=e, s=s, jphi=jphi, M=M,
                deep=deep, Klin=Klin, skiplin=skiplin, Gs=Gs, sub=sub,
                dh=dh, dm=dm, state=state, spares=spares, homog=homog,
                uweight=uweight, W_step=W_step, ordPhi=ordPhi)


# ---------------------------------------------------------------------------
# The published (72,108) generators (FULL_SYSTEM_BRIDGE.md sec.1, /2-cleared),
# used as the recipe control: the parametric builder must reproduce them.
# ---------------------------------------------------------------------------
def published_72108():
    d0, d1, d2 = sp.symbols("d0 d1 d2")
    dm1, dm2, dm3, dm4 = sp.symbols("dm1 dm2 dm3 dm4")
    return {
        1: sp.Rational(3, 2) * d1 * dm1**2 + 3 * d2 * dm1 * dm2 + 3 * dm1 * dm4
           + 3 * dm2 * dm3,
        2: -sp.Rational(3, 2) * d0 * dm1**2 + sp.Rational(3, 2) * d2 * dm2**2
           + 3 * dm2 * dm4 + sp.Rational(3, 2) * dm3**2,
        3: -3 * d0 * dm1 * dm2 - sp.Rational(3, 2) * d1 * dm2**2
           - sp.Rational(1, 2) * dm1**3 + 3 * dm3 * dm4,
        5: (-3 * d0 * dm1 * dm4 - 3 * d0 * dm2 * dm3 - 3 * d1 * dm2 * dm4
            - sp.Rational(3, 2) * d1 * dm3**2 - 3 * d2 * dm3 * dm4
            - sp.Rational(3, 2) * dm1**2 * dm3 - sp.Rational(3, 2) * dm1 * dm2**2)
           + Phi,
    }


def canonical_varorder(r):
    """State (d_{t-2}..d_0, e=dm1), then spares dm2..dm_{(a-1)t}, then Phi."""
    return [str(v) for v in r["state"]] + [str(v) for v in r["spares"]] + ["Phi"]


def emit_json(r, path):
    a, b, t = r["a"], r["b"], r["t"]
    varorder = canonical_varorder(r)
    gens = {}
    for j in sorted(r["Gs"]):
        name = "G%d" % j
        body = r["Gs"][j] - (Phi if j == r["jphi"] else 0)
        P = sp.Poly(body, *sorted(body.free_symbols, key=str))
        degs = [sum(mon) for mon, _ in P.terms()]
        gens[name] = dict(
            slice_j=j,
            slice_upower=b * t + j,
            u_weight=b * t + j,
            has_phi=(j == r["jphi"]),
            total_degree_min=min(degs),
            total_degree_max=max(degs),
            num_terms=len(sp.Add.make_args(body)),
            poly=sp.sstr(sp.expand(r["Gs"][j])),
        )
    W = r["W_step"]
    dossier = dict(
        schema="g-system-v1",
        case=dict(tag="F2_j1_75_125", degrees=[75, 125],
                  corner="(5,20)  [GGV5 chain final corner (7\\5,2); NOT the chart]",
                  a=a, b=b, t=t, kappa=r["kappa"],
                  q=r["q"], e=r["e"], s=r["s"],
                  chart_exponent_rule="t = l = ceil(b0/a0) = ceil(20/5) = 4 "
                                      "[INFERRED rule, guarded in "
                                      "polygon_reduction.py sec.0b]",
                  repaired="2026-07-26: t 5->4, kappa 3->2, q 2->1, C "
                           "y^2(y^3+1)->y; see PASSPORT_75_125_REPAIR.md"),
        recipe=dict(
            linear_window="S^a = S^%d  (slices [u^(a t + k)], linear in dm_{(a-1)t+k})" % a,
            forcing_window="S^b = S^%d  (slices [u^(b t + j)], degree-%d generators)" % (b, b),
            forcing_slice_M=r["M"], jphi=r["jphi"],
            linear_eliminations=len(r["sub"]),
            linear_slice_range="k = 1..%d, skip k=%d" % (r["Klin"], r["skiplin"]),
            forcing_slice_range="j = 1..%d, skip j=%d" % (r["jphi"], r["jphi"] - 1),
        ),
        variable_order=varorder,
        ring="Q[%s]" % ",".join(varorder),
        state_variables=[str(v) for v in r["state"]],
        state_dictionary="d3,d2,d1,d0 = d_{t-2..t-5}; e = dm1 = d_{-1}",
        spare_variables=[str(v) for v in r["spares"]],
        spare_dictionary="dm2..dm%d = d_{-2}..d_{-%d} (= d_{-2}..d_{-(a-1)t})"
                         % ((a - 1) * t, (a - 1) * t),
        num_spares=len(r["spares"]),
        generators=gens,
        generator_names=["G%d" % j for j in sorted(r["Gs"])],
        num_generators=len(r["Gs"]),
        forcing_window_power=b,
        generator_degree_note=(
            "Generators are WEIGHT-homogeneous under the u-grading, NOT "
            "total-degree homogeneous (as at (72,108): G1 mixes deg 2 and 3). "
            "The forcing window is S^%d, but the a=3 linear substitutions "
            "(cubic, S^%d slices) inflate the total degree of the deeper "
            "generators: total-degree maxima run %s across G1..G%d,G%d."
            % (b, a, [max(sum(m) for m, _ in sp.Poly(
                r["Gs"][j] - (Phi if j == r["jphi"] else 0),
                *sorted((r["Gs"][j] - (Phi if j == r["jphi"] else 0)).free_symbols,
                        key=str)).terms())
               for j in sorted(r["Gs"])], r["jphi"] - 2, r["jphi"])),
        u_grading_weights=[b * t + j for j in sorted(r["Gs"])],
        phi_u_weight=r["M"],
        weight_grading_note=(
            "Generators are weight-homogeneous under the INTRINSIC u-grading "
            "w(d_m)=t-m, w(Phi)=M=%d.  Forcing-generator weights form the "
            "arithmetic progression %s (common difference 1; the value %d is "
            "absent -- the skipped generator G%d)."
            % (r["M"], [b * t + j for j in sorted(r["Gs"])],
               b * t + (r["jphi"] - 1), r["jphi"] - 1)),
        physical_weight=dict(
            W_step_num=int(W.p), W_step_den=int(W.q), W_step=str(W),
            ord_y_Phi=r["ordPhi"], forcing_slice_M=r["M"],
            integral=(W.q == 1),
            physical_generator_weights=(
                [int(W * (b * t + j)) for j in sorted(r["Gs"])] if W.q == 1
                else [str(W * (b * t + j)) for j in sorted(r["Gs"])]),
            note=("W_step = ord_y(Phi)/M is INTEGER: the physical y-order grading "
                  "coincides with the u-grading scaled by %d; exact y^W stripping "
                  "(FULL_SYSTEM_BRIDGE sec.3) transfers." % int(W)
                  if W.q == 1 else
                  "W_step = ord_y(Phi)/M = %s is NON-INTEGRAL (denominator %d): the "
                  "physical y-valuation grading is rational, not an integer multiple "
                  "of the u-grading.  The window variables have quasi-affine (not "
                  "affine) y-order, so the forcing generators are homogeneous only "
                  "in the u-grading; the exact integer y^W stripping of the (72,108) "
                  "bridge does NOT transfer.  This is the a>=3 boundary predicted by "
                  "CORNER_144_COMPARISON.md sec.5 (quasipolynomial window cap "
                  "ceil(w/5) for (108,144), denominator 5; denominator %d here)."
                  % (str(W), W.q, W.q))),
        deg_slope=dict(
            note=("Whether the y-DEGREE cap is affine depends on deg_y(Phi)/M "
                  "being an integer.  When C is a monomial, deg_y(Phi) = "
                  "ord_y(Phi) and the deg-slope equals W_step, so it is affine "
                  "only when W_step is -- i.e. never in this family.  The "
                  "(72,108) picture of TWO distinct slopes (ord 12, deg 14) has "
                  "no counterpart here."),
            deg_equals_ord=True),
        window_caps=dict(
            pattern=("ord y quasi-affine: physical step %d/%d = %s per window unit "
                     "-> quasipolynomial caps with quasi-period %d; NOT the affine "
                     "12k/15k/14k of (72,108).  And since C is a monomial the "
                     "y-DEGREE cap coincides with the y-order cap (same slope), so "
                     "there is no affine upper cap either."
                     % (r["ordPhi"], r["M"], str(W), W.q)),
            quasi_period=int(W.q),
            quasi_period_is_prime=bool(sp.isprime(W.q)),
            status="OBSTRUCTED: affine window-cap slopes do not exist "
                   "(non-integral W_step); the caps are quasipolynomial",
        ),
        phi_consistency=dict(
            verdict="CONSISTENT (intrinsic grading)",
            detail=("Phi enters the deepest forcing generator G%d at u-slice "
                    "M=%d with the homogeneity-forced weight %d, matching phase-1's "
                    "tower slice (C_SERIES_75_125.md: Phi = f*C^%d is the u^%d slice "
                    "of S^%d, clear = a*M-b = %d, N = %d).  ord_y(Phi)=%d=W_step*M "
                    "holds as a rational identity; the departure from (72,108) is "
                    "only the non-integral physical normalisation."
                    % (r["jphi"], r["M"], r["M"],
                       a * r["M"] - 2 * b, r["M"], b, a * r["M"] - b,
                       a * r["M"] - 2 * b, r["ordPhi"])),
        ),
        bridge=("~%d quintic generator equations (before y-coefficient expansion); "
                "%d spare window unknowns dm2..dm%d + Phi over the state "
                "(d3,d2,d1,d0,e)"
                % (len(r["Gs"]), len(r["spares"]), (a - 1) * t)),
        obstruction=("Structure transfers (slice equations, %d-generator forcing "
                     "system, %d-spare inventory, u-grading AP of weights). The "
                     "PHYSICAL weight normalisation does not: W_step = %d/%d = %s "
                     "is non-integral, so the affine window caps and exact "
                     "y-stripping of (72,108) are replaced by a quasipolynomial "
                     "(quasi-period %d) window layer.  Because C is a monomial the "
                     "y-degree cap coincides with the y-order cap, so unlike "
                     "(72,108) there is no affine upper slope at all."
                     % (len(r["Gs"]), len(r["spares"]), r["ordPhi"], r["M"],
                        str(W), W.q)),
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(json.dumps(dossier, sort_keys=True, indent=1) + "\n")
    return dossier


# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("G-SYSTEM TRANSFER, phase 2 :  build the (75,125) D-transform G-system")
    print("=" * 78)

    # ---- recipe control: reproduce the published (72,108) generators ----
    print("\n[control] (72,108)  a,b,t = 2,3,4  q=7  ord_y(Phi)=204")
    c = build_gsystem(2, 3, 4, 7, 204, verbose=True)
    pub = published_72108()
    for j in (1, 2, 3, 5):
        assert sp.expand(c["Gs"][j] - pub[j]) == 0, "control G%d mismatch" % j
    uw = [c["b"] * c["t"] + j for j in sorted(c["Gs"])]
    phys = [int(c["W_step"] * w) for w in uw]
    print("  reproduced G1,G2,G3,(G5body+Phi) EXACTLY (FULL_SYSTEM_BRIDGE sec.1)")
    print("  spares:", [str(s) for s in c["spares"]])
    print("  u-grading weights:", uw, "  W_step =", c["W_step"],
          "  physical weights:", phys)
    assert phys == [156, 168, 180, 204], phys
    print("  physical weights = [156,168,180,204]  MATCH known G-weights")

    # ---- the target: (75,125) ----
    # inputs DERIVED through the retraction guard, never transcribed
    import polygon_reduction as pr
    _cd = pr.corner_chart_data(5, 20, l_final=5, b_final=2, who="g_system_75_125")
    assert not _cd["retraction"] and _cd["monomial"], _cd
    _t, _q = _cd["t"], _cd["ord_C"]
    _a, _b = 3, 5
    _N = _a * (_t * (_a + _b) - (_cd["kappa"] + 1)) - 2 * _b
    _rho = (_b - _a + 1 - 1) * _q + 1
    _ordPhi = _rho + _N * _q                       # Phi = (1/a) y^(rho + N q)
    assert (_t, _q, _N, _ordPhi) == (4, 1, 77, 80), (_t, _q, _N, _ordPhi)
    print("\n[target] (75,125)  a,b,t = 3,5,%d  q=%d  ord_y(Phi)=%d   [REPAIRED: t "
          "was 5 and q was 2, from GGV5's final corner; the retraction guard now "
          "refuses that]" % (_t, _q, _ordPhi))
    r = build_gsystem(_a, _b, _t, _q, _ordPhi, verbose=True)
    uw = [r["b"] * r["t"] + j for j in sorted(r["Gs"])]
    print("  generators:", ["G%d" % j for j in sorted(r["Gs"])],
          "(quintic; slice j=1..%d skip %d)" % (r["jphi"], r["jphi"] - 1))
    print("  count law: #generators = jphi-1 = a*t-kappa-2 = %d;  #spares = "
          "(a-1)t-1 = %d  (both increment by t per rung, so t=4 gives +4/+4 where "
          "the superseded t=5 gave +5/+5)" % (r["jphi"] - 1, (r["a"] - 1) * r["t"] - 1))
    print("  spare inventory (%d):" % len(r["spares"]),
          [str(s) for s in r["spares"]], "= d_-2 .. d_-%d" % ((r["a"] - 1) * r["t"]))
    for j in sorted(r["Gs"]):
        assert r["homog"][j] == {r["b"] * r["t"] + j}, (j, r["homog"][j])
    print("  every generator is u-grading homogeneous; weights (AP):", uw)
    print("  Phi u-weight:", r["M"], " forcing slice M = b*t + jphi =", r["M"])
    W = r["W_step"]
    print("  physical normalisation  W_step = ord_y(Phi)/M = %d/%d ="
          % (r["ordPhi"], r["M"]), W, "(integral: %s)" % (W.q == 1))
    print("  because C = y is a monomial, deg_y(Phi) = ord_y(Phi) = %d, so the "
          "deg-slope EQUALS the ord-slope: there is no second, affine cap."
          % r["ordPhi"])
    print("  ->", "clean stripping transfers" if W.q == 1 else
          "NON-INTEGRAL: exact y-stripping does NOT transfer (a=3 boundary; "
          "quasi-period %d window caps)" % W.q)

    # slice-sum / N consistency with phase 1
    clear = r["a"] * r["M"] - r["b"]
    N = clear - r["b"]
    print("  phase-1 slice-sum: clear = a*M - b =", clear, " N = clear - b =", N,
          "(matches C_SERIES_75_125.md N=77:", N == 77, ")")

    out = "g_system_75_125.json"
    emit_json(r, out)
    print("\nwrote", out)

    print("\n" + "=" * 78)
    print("VERDICT: G-system BUILT (structure transfers) with a characterised")
    print("physical-weight OBSTRUCTION (W_step = %s non-integral)." % r["W_step"])
    print("  spares: %d  (dm2..dm%d)   generators: %d (quintic)"
          % (len(r["spares"]), (r["a"] - 1) * r["t"], len(r["Gs"])))
    print("  u-grading weights: %s   Phi-consistency: CONSISTENT"
          % [r["b"] * r["t"] + j for j in sorted(r["Gs"])])
    print("=" * 78)


if __name__ == "__main__":
    main()
