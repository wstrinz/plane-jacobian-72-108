#!/usr/bin/env python3
"""chain_survey.py  (NEW; read-only over all existing artifacts)

THE CHAIN NATURAL-HISTORY SURVEY -- a faithful re-implementation of GGV5's
complete-chain enumeration algorithm, run past the published v11<=35 tables to
map the demography of the problem.

WHAT THIS IS.  GGV5 (paper_src/1708.07936_GGV5.tex) gives an algorithm that,
for a bound M on v11(A_0), enumerates every ADMISSIBLE COMPLETE CHAIN and, for
each chain's final corner, the (m,n)-families that could realize it.  The
published output is the pair of tables at lines 1674-1718 (17 length-1 families
F1-F17, 7 length-2 families F18-F24).  This file reproduces those tables and
extends the enumeration.

The pipeline reproduced here (GGV5 algorithm names in brackets):
  * get_pllc            [GetPossibleLastLowerCorners, Alg. lines 345-383]
  * starting_edges      [GetStartingEdges,            Alg. lines 473-500]
  * generate_corners    [GetGeneratedCorners,         Alg. lines 708-741]
  * corner_children     [GetCornerChildrenList,       Alg. lines 786-819]
  * complete_chains     [GetCompleteChains,           Alg. lines 907-939]
  * get_is_admissible   [GetIsAdmissible,             Alg. lines 1296-1328]
  * mn_families         [GetmnFamilies,               Alg. lines 1599-1635]
  * enumerate_survey    [Main algorithm,              Alg. lines 1336-1362]

REGRESSION VERDICT (see chain_survey_verify.py for the exact checker):
  All 24 published CHAINS (A_0, A_0', ..., final corner, k) are reproduced
  exactly.  23 of the 24 (m,n)-family parametrizations match the printed table
  verbatim.  The sole exception is F6: the paper prints base pair (4,10), which
  has gcd 2 and therefore VIOLATES the coprimality that Definition "mn families"
  (line 1500) demands -- a typo.  The algorithm's coprime output (7,18) with
  step (6,16) is the correction; it satisfies the same Diophantine identity
  (m+n)bk - n(bl-a) = k (line 1480).

CANONICAL DEDUP.  The literal algorithm also emits redundant chains whose FIRST
edge is degenerate-vertical (A_0' with b'=0, direction (1,0), no lattice
refinement rho=1): its type-IIb "column" generation reproduces a corner A_1 that
is already reached directly as the type-IIa step A_0'=A_1.  The published table
keeps the canonical A_0'=A_1 form.  We dedup complete chains by their generated-
corner sequence (everything after the first corner), preferring the A_0'=A_1
representative -- this yields exactly the 24 published rows.  Both the raw and
canonical counts are recorded in the survey data.

Conventions match phi_corner4.py / composite_charts.py:
  corner (a\\l,b) stored as the integer triple (a,l,b), realization (a/l, b);
  t := l_final, kappa := t-2 (standard single-chart class), a0 := A_0[0],
  q := b_final, dg := a0-q, e := |m-n|+1, r := a0-q-1.

Exact integer/Fraction arithmetic only; no sympy needed for the enumeration.
Outputs chain_survey_data.json.  Checker: chain_survey_verify.py.
"""
import json
import sys
import time
from fractions import Fraction as Fr
from math import gcd, lcm, floor, ceil

# ---------------------------------------------------------------------------
# geometry helpers on rational points  (realizations (a/l, b))
# ---------------------------------------------------------------------------
def dir_of(P, Q):
    """dir(P-Q): primitive (rho,sigma), rho>0, with rho*dx + sigma*dy = 0
    (the (rho,sigma)-homogeneous direction of the edge P--Q)."""
    dx = P[0] - Q[0]
    dy = P[1] - Q[1]
    den = lcm(dx.denominator, dy.denominator)
    rx = int(dy * den)      # rho ~  dy
    ry = int(-dx * den)     # sigma ~ -dx
    g = gcd(abs(rx), abs(ry))
    rx //= g
    ry //= g
    if rx < 0:
        rx, ry = -rx, -ry
    if rx == 0:
        ry = 1 if ry > 0 else -1
    return (rx, ry)

def vrs(rho, sigma, P):
    return rho * P[0] + sigma * P[1]

def real(corner):
    a, l, b = corner
    return (Fr(a, l), Fr(b))

def slope_lt(d1, d2):
    """d1 < d2 by slope sigma/rho (rho>=0); rho=0 == slope -inf (smallest)."""
    r1, s1 = d1
    r2, s2 = d2
    if r1 == 0 and r2 == 0:
        return False
    if r1 == 0:
        return True
    if r2 == 0:
        return False
    return s1 * r2 < s2 * r1

def bigomega(n):
    """Number of prime factors counted with multiplicity."""
    n = abs(int(n))
    c = 0
    d = 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            c += 1
        d += 1
    if n > 1:
        c += 1
    return c

# ---------------------------------------------------------------------------
# GetPossibleLastLowerCorners  (GGV5 Alg., lines 345-383)
# ---------------------------------------------------------------------------
def get_pllc(xmax):
    PLLC = set()
    PFL = []                       # list of ((a,b),(rho,sigma)), insertion order
    for a in range(1, xmax + 1):
        b = 0
        # while b <= 1/2(2a - sqrt(4a-3) - 1); the paper's exact integer equivalent
        # is b < a and b <= (a-b-1)^2.
        while b < a and b <= (a - b - 1) ** 2:
            if b == 0:
                PFL.append(((a, b), (0, -1)))
                PLLC.add((a, b))
            elif a > 2 * b and b > 0:
                PFL.append(((a, b), (1, -2)))
                PLLC.add((a, b))
            else:
                rs_ab = (1, -1)
                for ((r, s), rs_rs) in PFL:
                    if not (r < a and s < b and (r - s) < (a - b)):
                        continue
                    N1 = gcd(a - r, b - s)
                    N2 = gcd(r, s)
                    rho = (b - s) // N1
                    sigma = (r - a) // N1
                    val = rho * a + sigma * b
                    g = gcd(abs(rho + sigma), abs(val))
                    if g == 0:
                        continue
                    vth = val // g
                    cond_dir = slope_lt(rs_rs, (rho, sigma)) and slope_lt((rho, sigma), rs_ab)
                    cond_v = val >= rho
                    cond_th = (vth <= N1) or (vth != 0 and N2 % vth == 0)
                    if cond_dir and cond_v and cond_th:
                        rs_ab = (rho, sigma)
                if slope_lt(rs_ab, (1, -1)):
                    PFL.append(((a, b), rs_ab))
                    PLLC.add((a, b))
            b += 1
    return PLLC

# ---------------------------------------------------------------------------
# valid-edge geometry
# ---------------------------------------------------------------------------
def is_simple(A, Ap):
    """The 'simple' predicate of Definition 'valid edges' (line 445), realised
    through Remark 'zazaza' (mu/d = (rho+sigma)/v_rho,sigma(A))."""
    a, l, b = A
    _, _, bp = Ap
    RA = real(A)
    RAp = real(Ap)
    rho, sigma = dir_of(RA, RAp)
    vA = vrs(rho, sigma, RA)
    if vA == 0:
        return False
    f2 = Fr(rho + sigma) / vA * b           # v01(enF)
    gap = rho // gcd(rho, l)
    return (f2 == gap + 1) and (gap > 1 or bp > 0)

def is_final(corner):
    a, l, b = corner
    if b == 0:
        return False
    return Fr(l) - Fr(a, b) > 1

# ---------------------------------------------------------------------------
# GetGeneratedCorners  (GGV5 Alg., lines 708-741)
# ---------------------------------------------------------------------------
def generate_corners(A, Ap):
    a, l, b = A
    _, _, bp = Ap
    RA = real(A)
    RAp = real(Ap)
    rho, sigma = dir_of(RA, RAp)
    if vrs(1, -1, RAp) < 0:
        return [Ap]                         # type-IIa: generated corner is A'
    l1 = lcm(rho, l)
    gap = rho // gcd(rho, l)
    gmax = int(min(Fr(b - bp, gap), Fr(b - 1)))

    def make(b1):
        a1 = Fr(a * l1, l) + Fr((b1 - b) * (-sigma) * l1, rho)
        return (int(a1), l1, b1)

    out = []
    if is_simple(A, Ap):
        A1 = make(gmax)
        a1, ll, bb = A1
        if bb != 0 and ((Fr(ll) - Fr(a1, bb) > 1) or gcd(a1, bb) > 1):
            out.append(A1)
    else:
        for b1 in range(bp + 1, gmax + 1):
            A1 = make(b1)
            a1, ll, bb = A1
            if vrs(1, -1, real(A1)) < 0 and ((Fr(ll) - Fr(a1, bb) > 1) or gcd(a1, bb) > 1):
                out.append(A1)
    return out

# ---------------------------------------------------------------------------
# GetCornerChildrenList  (GGV5 Alg., lines 786-819)
# ---------------------------------------------------------------------------
def corner_children(A, Ap, A1, PLLC):
    RA = real(A)
    RAp = real(Ap)
    rho, sigma = dir_of(RA, RAp)
    a1, l1, b1 = A1
    d1 = gcd(a1, b1)
    vA1 = vrs(rho, sigma, real(A1))
    lo = floor(1 + Fr(d1 * (rho + sigma), vA1))
    hi = d1
    if l1 > 1:
        hi = floor(l1 * (b1 * l1 - a1) + Fr(d1, b1))
    out = []
    for mu in range(lo, hi + 1):
        enF = (Fr(mu * a1, d1 * l1), Fr(mu * b1, d1))
        rho1, sigma1 = dir_of(enF, (Fr(1), Fr(1)))
        gapp = rho1 // gcd(rho1, l1)
        if gapp <= b1 and mu % d1 != 0:
            for j in range(1, b1 // gapp + 1):
                apx = Fr(a1, l1) + j * Fr(gapp * sigma1, rho1)
                apy = b1 - j * gapp
                num = apx * l1
                if num.denominator != 1:
                    continue
                A1p = (int(num), l1, apy)
                v = vrs(1, -1, (apx, Fr(apy)))
                ok = ((l1 > 1 and v != 0) or (l1 == 1 and v < 0) or
                      (l1 == 1 and v > 0 and (int(apx), apy) in PLLC))
                if ok:
                    out.append((A1, A1p))
    return out

def children_and_final(A, Ap, PLLC):
    """GetChildrenAndFinalList (lines 837-863): every generated corner is
    tested for finality AND handed to GetCornerChildrenList (the paper runs
    the latter unconditionally -- this is what produces the F22-F24 branch at
    the final-shaped corner (14/4,6))."""
    gens = generate_corners(A, Ap)
    children = []
    finals = []
    for A1 in gens:
        _, _, b1 = A1
        if is_final(A1):
            finals.append(A1)
        if b1 != 0:
            children.extend(corner_children(A, Ap, A1, PLLC))
    return children, finals

# ---------------------------------------------------------------------------
# GetStartingEdges  (GGV5 Alg., lines 473-500)
# ---------------------------------------------------------------------------
def starting_edges(a, b, PLLC):
    d = gcd(a, b)
    out = []
    for mu in range(1, d):
        enF = (Fr(mu * a, d), Fr(mu * b, d))
        rho, sigma = dir_of(enF, (Fr(1), Fr(1)))
        if rho <= 0:
            continue
        for i in range(1, b // rho + 1):
            apx = a + i * sigma
            apy = b - i * rho
            v = apx - apy
            if v < 0 or (v > 0 and (apx, apy) in PLLC):
                out.append(((a, 1, b), (apx, 1, apy)))
    seen = set()
    ded = []
    for e in out:
        if e not in seen:
            seen.add(e)
            ded.append(e)
    return ded

# ---------------------------------------------------------------------------
# GetCompleteChains  (GGV5 Alg., lines 907-939)
# ---------------------------------------------------------------------------
def complete_chains(C0, PLLC):
    A, Ap = C0
    _, _, b = A
    _, _, bp = Ap
    rho, _ = dir_of(real(A), real(Ap))
    Lmax = bigomega(gcd(b, (b - bp) // rho)) + 1
    open_chains = [[C0]]
    done = []
    j = 0
    while j < Lmax:
        popen = []
        for CH in open_chains:
            last = CH[-1]
            ch, fn = children_and_final(last[0], last[1], PLLC)
            for A1 in fn:
                done.append(CH + [('FINAL', A1)])
            for (A1, A1p) in ch:
                popen.append(CH + [(A1, A1p)])
        open_chains = popen
        j += 1
    return done

# ---------------------------------------------------------------------------
# GetIsAdmissible  (GGV5 Alg., lines 1296-1328; Definition cond div, line 1287)
# ---------------------------------------------------------------------------
def _q_of(edge):
    A, Ap = edge
    a, l, b = A
    rho, sigma = dir_of(real(A), real(Ap))
    val = rho * a + sigma * b * l              # = l * v_rho,sigma(A)
    g = gcd(abs((rho + sigma) * l), abs(val))
    return val // g                            # q_h = v/gcd(rho+sigma, v)

def _gap_of(edge):
    A, _ = edge
    _, l, _ = A
    rho, _ = dir_of(real(A), real(edge[1]))
    return rho // gcd(rho, l)

def get_is_admissible(CH):
    edges = [c for c in CH if c[0] != 'FINAL']
    j = len(edges) - 1
    if j <= 0:
        return True
    ok = True
    h = 0
    while h < j and ok:
        Ah, Aph = edges[h]
        ah, lh, bh = Ah
        aph, _, bph = Aph
        gap = _gap_of(edges[h])
        q = _q_of(edges[h])
        bhn = edges[h + 1][0][2]                 # b_{h+1}
        i = h + 1
        while i <= j and ok:
            Ai = edges[i][0]
            _, li, _ = Ai
            qp = _q_of(edges[i])
            D = gcd((bh - bph) // gap, bh)
            D = gcd(D, bhn)
            D = gcd(D, ah * li // lh)
            D = gcd(D, aph * li // lh)
            if bigomega(D) >= i - h and (qp != 0 and D % qp == 0) and (qp % q != 0):
                i += 1
            else:
                ok = False
        h += 1
    return ok

# ---------------------------------------------------------------------------
# BezoutCoefficients + GetmnFamilies  (GGV5 Alg., lines 1597-1635)
# ---------------------------------------------------------------------------
def bezout_min(x, y):
    """(M,N) positive with Mx - Ny = 1 and N minimal (>=1)."""
    assert gcd(x, y) == 1
    for N in range(1, x + 1):
        if (1 + N * y) % x == 0:
            return ((1 + N * y) // x, N)
    raise ValueError("no Bezout")

def mn_families(final):
    a, l, b = final
    bl_a = b * l - a
    kmax = ceil(Fr(l) - Fr(a, b)) - 1
    res = []
    for k in range(1, kmax + 1):
        e = gcd(k, bl_a)
        if gcd(b, bl_a // e) != 1:
            continue
        M, N = bezout_min(b, bl_a // e)
        n = N * k // e
        m = M - n
        D1 = (bl_a - b * k) // e
        D2 = b * k // e
        if m == 1 or n == 1:
            m += D1
            n += D2
        kbar = k // e
        if kbar == 1:
            res.append((k, 0, (m, n), (D1, D2)))
        else:
            for i in range(kbar):
                mi = m + i * D1
                ni = n + i * D2
                if gcd(mi, ni) == 1:
                    res.append((k, i, (mi, ni), (kbar * D1, kbar * D2)))
    return res

# ---------------------------------------------------------------------------
# Main algorithm + canonical dedup
# ---------------------------------------------------------------------------
def enumerate_complete(M):
    """All complete chains with v11(A_0) <= M (before admissibility)."""
    PLLC = get_pllc(M // 2)
    chains = []
    for a in range(2, M // 2 + 1):
        for b in range(a + 1, M - a + 1):
            for C0 in starting_edges(a, b, PLLC):
                chains.extend(complete_chains(C0, PLLC))
    return chains

def canonical_admissible(chains):
    """Filter to admissible, then dedup by generated-corner sequence (drop the
    redundant vertical-first-edge duplicates), preferring the A_0'=A_1 form."""
    adm = [CH for CH in chains if get_is_admissible(CH)]
    groups = {}
    for CH in adm:
        key = (CH[0][0],) + tuple(CH[1:])       # A_0 + everything after C_0's corner
        A0, A0p = CH[0]
        typeIIa = (len(CH) > 1 and CH[1][0] != 'FINAL' and A0p == CH[1][0])
        endsfinal = (len(CH) > 1 and CH[1][0] == 'FINAL')
        pref = 1 if (typeIIa or endsfinal) else 0
        cur = groups.get(key)
        if cur is None or pref > cur[0]:
            groups[key] = (pref, CH)
    return adm, [v[1] for v in groups.values()]

# ---------------------------------------------------------------------------
# per-family statistics  (conventions of phi_corner4.py)
# ---------------------------------------------------------------------------
def fmt_corner(c):
    a, l, b = c
    return f"({a}/{l},{b})" if l > 1 else f"({a},{b})"

def chain_signature(CH):
    """Return list of (corner, prime) pairs as readable strings + directions."""
    edges = [c for c in CH if c[0] != 'FINAL']
    final = CH[-1][1]
    corners = [fmt_corner(e[0]) for e in edges]
    primes = [fmt_corner(e[1]) for e in edges]
    dirs = [list(dir_of(real(e[0]), real(e[1]))) for e in edges]
    lvals = [e[0][1] for e in edges] + [final[1]]
    return corners, primes, dirs, lvals

def family_record(CH, k, i, base, step):
    A0, A0p = CH[0]
    final = CH[-1][1]
    a0 = A0[0]
    t = final[1]                 # l_final
    q = final[2]                 # b_final
    kappa = t - 2                # standard single-chart class (phi_corner4)
    dg = a0 - q                  # residual-divisor degree
    r = a0 - q - 1
    m0, n0 = base
    dm, dn = step
    a_small, b_big = sorted((m0, n0))
    e = b_big - a_small + 1
    # Diophantine identity  (m+n)bk - n(bl-a) = k   (GGV5 line 1480, final corner)
    af, lf, bf = final
    bl_a = bf * lf - af
    dio_lhs = (m0 + n0) * bf * k - n0 * bl_a
    # window-denominator analogue (F2's 5a-3 == t*a - (t-2) = t*a - kappa; see
    # f2_family_verify.py).  Reduced against the smaller family element.
    win_raw = t * a_small - kappa
    win_den = win_raw // gcd(win_raw, a_small) if win_raw else 0
    corners, primes, dirs, lvals = chain_signature(CH)
    length = len(corners)
    return {
        "A0": list(A0), "A0p": list(A0p),
        "final": [af, lf, bf], "final_str": fmt_corner(final),
        "length": length,
        "k": k, "i": i,
        "m0": m0, "n0": n0, "dm": dm, "dn": dn,
        "t": t, "kappa": kappa, "a0": a0, "q": q,
        "dg": dg, "r": r, "e_at_base": e,
        "win_denom_base": win_raw, "win_denom_reduced": win_den,
        "dio_residual_minus_k": dio_lhs - k,
        "motif": f"({A0[0]},{A0[2]})",
        "final_type": f"(l={lf},q={bf},k={k})",
        "corner_sig": [a0, t, q, k],            # j-independent corner signature
        "chain_corners": corners, "chain_primes": primes,
        "chain_dirs": dirs, "l_values": lvals,
    }

def survey(M):
    """Return (family_records, meta) for bound M."""
    chains = enumerate_complete(M)
    adm_raw, canon = canonical_admissible(chains)
    recs = []
    for CH in canon:
        final = CH[-1][1]
        for (k, i, base, step) in mn_families(final):
            recs.append(family_record(CH, k, i, base, step))
    meta = {
        "M": M,
        "n_complete_chains": len(chains),
        "n_admissible_raw": len(adm_raw),
        "n_canonical_chains": len(canon),
        "n_family_rows": len(recs),
    }
    return recs, meta

# ---------------------------------------------------------------------------
# demography aggregation over a sweep of bounds
# ---------------------------------------------------------------------------
def demography(sweep):
    per_M = []
    full = None
    t_start = time.time()
    for M in sweep:
        t0 = time.time()
        recs, meta = survey(M)
        meta["seconds"] = round(time.time() - t0, 1)
        # clustering statistics
        motifs = set(r["motif"] for r in recs)
        final_types = set(r["final_type"] for r in recs)
        corner_sigs = set(tuple(r["corner_sig"]) for r in recs)
        t_values = sorted(set(r["t"] for r in recs))
        max_len = max((r["length"] for r in recs), default=0)
        by_len = {}
        for r in recs:
            by_len[r["length"]] = by_len.get(r["length"], 0) + 1
        meta.update({
            "n_motifs_A0": len(motifs),
            "n_final_types_lqk": len(final_types),
            "n_corner_signatures": len(corner_sigs),
            "t_values": t_values,
            "max_chain_length": max_len,
            "family_rows_by_length": by_len,
        })
        per_M.append(meta)
        full = recs                     # keep the largest-M full record set
        print(f"[M={M:>3}] complete={meta['n_complete_chains']:>6} "
              f"adm={meta['n_admissible_raw']:>5} canon={meta['n_canonical_chains']:>5} "
              f"fams={meta['n_family_rows']:>4} motifs={len(motifs):>3} "
              f"final_types={len(final_types):>2} corner_sigs={len(corner_sigs):>3} "
              f"maxlen={max_len} t={t_values} ({meta['seconds']}s)", flush=True)
    print(f"total sweep time {time.time()-t_start:.1f}s", flush=True)
    return per_M, full

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", type=str, default="35,45,55,65,75,85,100")
    ap.add_argument("--out", type=str, default="chain_survey_data.json")
    args = ap.parse_args()
    sweep = [int(x) for x in args.sweep.split(",")]
    per_M, full = demography(sweep)
    data = {
        "generated_by": "chain_survey.py",
        "source": "GGV5 = paper_src/1708.07936_GGV5.tex (tables lines 1674-1718)",
        "sweep": sweep,
        "per_M": per_M,
        "families_at_max_M": full,
    }
    with open(args.out, "w") as f:
        json.dump(data, f, indent=1)
    print(f"wrote {args.out}: {len(full)} family rows at M={sweep[-1]}")

if __name__ == "__main__":
    main()
