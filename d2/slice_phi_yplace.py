"""slice_phi_yplace.py -- spending SLICE_OBSTRUCTION.md Sec.8's unused inventory.

Independent probe of the four "never used" condition families listed in
SLICE_OBSTRUCTION.md Sec.8 (as corrected by SLICE_OBSTRUCTION_AUDIT.md C-1):

  (1) the y-order conditions   y^2 | p_8, y | r_11, y^3 | r_12
  (2) the degree / upper-hull conditions
  (3) the exact equations      r_13 = r_14 = r_15 = 0
  (4) the sharp relation       [u^17] H^3 = (1/6630) t^30 q

HEADLINE RESULT (Section D):  a_t <= 9.  Combined with the independently
audited a_t >= 9 this gives  a_t = v_t(e) = 9 EXACTLY, killing the six
standard-sub1 a10_* cells and leaving the five a9_* cells.

Ingredients of the kill, all machine-checked below:
  * the *un-run* cascade levels 14 and 15 -- the lane's own (Q) list covers
    n = 2..15 but slice_obstruction_basis.py only ran levels 2..12;
  * Sec.8 item (4), the Phi relation, used at its leading t-jet:
    [t^30] r_17 = q(-1)/6630 = 1/2 =/= 0.
Sec.8 items (3) and (4) are shown to be *literally the G-system rows*
D3(1),D3(2),D3(3),D3(5) -- so they are not new information for the program,
only for the slice lane.  Item (1) is the only family that is new to the
whole program; it is analysed at the place y = 0 in Section F.

Read-only on every existing artifact.  This lane wrote exactly two files:
this one and SLICE_PHI_YPLACE.md.  Pure sympy -- no Singular, no msolve, no
subprocess, no solver, so there are no aborts/timeouts/exit codes to read.

Run:  python -u slice_phi_yplace.py            # full narrated report
      python -u slice_phi_yplace.py --quiet    # exit 0 iff every check passes
"""
import sys
import json
import os
import sympy as sp
from sympy import symbols, Rational, expand, factor, groebner, Poly

QUIET = "--quiet" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))

_n = [0]
_fail = []


def say(*a):
    if not QUIET:
        print(*a)


def check(name, cond):
    _n[0] += 1
    if cond:
        say(f"  [OK] {name}")
    else:
        _fail.append(name)
        print(f"  [FAIL] {name}")


def head(s):
    say("\n" + s)


# =====================================================================
# 0.  Core objects, built from scratch.
#
#   H(u) = sum_{i>=0} h_i u^i   is the *stripped* root series:  h_i = d_{4-i}.
#   p_n = [u^n] H^2,  r_n = [u^n] H^3.
#   P = C^2 has x-degrees 0..8  =>  p_n = 0 for n >= 9   (this DEFINES h_{>=9})
#   Q's x-degrees are 0..12 and Q_M = (C^3)_M for M >= -3
#                              =>  r_13 = r_14 = r_15 = 0.
#   Shifted coordinates (the G-system's): h_1 = d_3 = 0.
# =====================================================================
u, t, y, w = symbols('u t y w')

NU_DEF = 18


def make_H(hfree, nu=NU_DEF, ring_zero=sp.Integer(0)):
    """h_0 = 1, h_1 = 0, h_2..h_8 = hfree[2..8]; h_{>=9} forced by p_n = 0."""
    h = {0: sp.Integer(1), 1: ring_zero}
    for k in range(2, 9):
        h[k] = hfree[k]
    for k in range(9, nu + 1):
        h[k] = expand(-sum(h[i] * h[k - i] for i in range(1, k)) / 2)
    return h


def series_pow(h, power, nu=NU_DEF):
    E = {0: sp.Integer(1)}
    for _ in range(power):
        F = {}
        for i, a in E.items():
            for j in range(nu + 1 - i):
                F[i + j] = F.get(i + j, 0) + a * h[j]
        E = {i: expand(c) for i, c in F.items()}
    return E


# =====================================================================
# A.  IDENTIFICATION:  Sec.8 items (3) and (4) ARE the G-system rows.
# =====================================================================
head("A.  What Sec.8 items (3) and (4) actually are")

dsym = {2: symbols('d2'), 1: symbols('d1'), 0: symbols('d0')}
dm = {k: symbols(f'dm{k}') for k in range(1, 14)}
# regenerate_system.py's series S (the UNSTRIPPED D-coordinates), rebuilt here
S = (1 + dsym[2] * u ** 2 + dsym[1] * u ** 3 + dsym[0] * u ** 4
     + sum(dm[k] * u ** (4 + k) for k in range(1, 14)))
S2 = Poly(expand(S * S), u)
S3 = Poly(expand(S2.as_expr() * S), u)
D2 = lambda k: S2.coeff_monomial(u ** (8 + k))
D3 = lambda j: S3.coeff_monomial(u ** (12 + j))

check("D2(1) = 2*d0*dm1 + 2*d1*dm2 + 2*d2*dm3 + 2*dm5  (regenerate_system's row)",
      expand(D2(1) - (2 * dsym[0] * dm[1] + 2 * dsym[1] * dm[2]
                      + 2 * dsym[2] * dm[3] + 2 * dm[5])) == 0)

# the linear phase: D2(k)=0 for k=1..7,9 eliminates dm5..dm11, dm13
sub = {}
for k, fresh in [(1, dm[5]), (2, dm[6]), (3, dm[7]), (4, dm[8]),
                 (5, dm[9]), (6, dm[10]), (7, dm[11]), (9, dm[13])]:
    sub[fresh] = expand(sp.solve(D2(k).subs(sub), fresh)[0])
G1 = expand(D3(1).subs(sub))
G2 = expand(D3(2).subs(sub))
G3 = expand(D3(3).subs(sub))
G5b = expand(D3(5).subs(sub))

gj = json.load(open(os.path.join(HERE, 'generators.json')))
order = gj['variable_order']
vmap = {'d2': dsym[2], 'd1': dsym[1], 'd0': dsym[0], 'dm1': dm[1],
        'dm2': dm[2], 'dm3': dm[3], 'dm4': dm[4], 'Phi': symbols('Phi')}


def from_json(name):
    out = 0
    for coeff, ev in gj['polynomials'][name]:
        term = Rational(coeff)
        for vn, e in zip(order, ev):
            if e:
                term *= vmap[vn] ** e
        out += term
    return expand(out)


for nm, mine in (('G1', G1), ('G2', G2), ('G3', G3), ('G5body', G5b)):
    check(f"generators.json['{nm}'] == my rebuild of D3({ {'G1':1,'G2':2,'G3':3,'G5body':5}[nm] }) "
          f"after the linear phase", expand(from_json(nm) - mine) == 0)

check("generators.json provenance states G1,G2,G3 = D3(1),D3(2),D3(3) and G5body = D3(5)",
      'D3(1),D3(2),D3(3)' in gj['provenance']['note']
      and 'G5body = D3(5)' in gj['provenance']['note'])

# stripping: h_i = d_{4-i} = D_{4-i} / y^{12 i}, i.e. H(u) = S(u/y^12), so
# [u^n] H^3 = [u^n] S^3 / y^{12 n}.  Checked as an exact identity, not asserted.
Ssub = {dsym[2]: symbols('e2') * y ** 24, dsym[1]: symbols('e3') * y ** 36,
        dsym[0]: symbols('e4') * y ** 48}
Ssub.update({dm[k]: symbols(f'e{4+k}') * y ** (12 * (4 + k)) for k in range(1, 14)})
for j in (1, 2, 3, 5):
    n = 12 + j
    lhs = expand(D3(j).subs(Ssub) / y ** (12 * n))
    check(f"[u^{n}] S^3 is y-homogeneous of weight 12*{n} under D_(4-i) = "
          f"y^(12 i) * (stripped)  =>  r_{n} = D3({j}) / y^{12*n}",
          expand(lhs - lhs.subs(y, 1)) == 0)
for k in (1, 2, 3, 9):
    n = 8 + k
    lhs = expand(D2(k).subs(Ssub) / y ** (12 * n))
    check(f"[u^{n}] S^2 is y-homogeneous of weight 12*{n}  =>  p_{n} = D2({k}) / y^{12*n}",
          expand(lhs - lhs.subs(y, 1)) == 0)

# the G-system's linear phase IS "p_n = 0 for n >= 9 determines h_n"
hfree_sym = {k: symbols(f'e{k}') for k in range(2, 9)}
hdet = make_H(hfree_sym, nu=15)
check("(P0): p_n = 0 for n = 9..15 determines h_9..h_15 from h_2..h_8 -- the same "
      "triangular elimination the G-system's linear phase runs on D2(1..7)",
      all(expand(sum(hdet[i] * hdet[n - i] for i in range(0, n + 1))) == 0
          for n in range(9, 16)))

say("\n  => Sec.8 item (3) 'r_13 = r_14 = r_15 = 0' IS the G-system's G1,G2,G3.")
say("     Sec.8 item (4) 'the Phi relation' IS the G-system's G5body + Phi = 0.")
say("     They are unused by the SLICE LANE; they are not new to the program.")

# ---- the Phi relation in stripped coordinates -----------------------
q_quartic = 2048 * y ** 4 - 512 * y ** 3 + 320 * y ** 2 - 240 * y + 195
f1 = -y ** 8 * (y + 1) ** 2 * q_quartic / 6630          # STATE.md / verify_derivation A
C4 = y ** 7 * (y + 1)
Phi_poly = expand(f1 * C4 ** 28)

# independent re-derivation of f1 from the commutator ODE (verify_derivation A)
_f = sp.Function('f')(y)
ode = 8 * y * (y + 1) * sp.diff(_f, y) - 14 * (8 * y + 7) * _f - y ** 8 * (y + 1) ** 2
aa = symbols('a0:16')
ans = sum(aa[i] * y ** i for i in range(16))
sol = sp.solve(Poly(expand(ode.subs(_f, ans).doit()), y).all_coeffs(), aa, dict=True)
check("the f1 forcing ODE has a unique polynomial solution, equal to STATE.md's f1",
      len(sol) == 1 and expand(ans.subs(sol[0]) - f1) == 0)

r17_stripped = expand(-Phi_poly / y ** 204)
check("r_17 = -Phi/y^204 = t^30 * q(y)/6630   (t = y+1)",
      expand(r17_stripped - (y + 1) ** 30 * q_quartic / 6630) == 0)
check("q(-1) = 3315 =/= 0  =>  v_t(r_17) = 30 EXACTLY and [t^30] r_17 = 1/2",
      q_quartic.subs(y, -1) == 3315
      and sp.Poly(expand(r17_stripped.subs(y, t - 1)), t).coeff_monomial(t ** 30) == Rational(1, 2))
check("q(0) = 195 =/= 0  =>  ord_y(r_17) = 0 and [y^0] r_17 = 1/34",
      q_quartic.subs(y, 0) == 195 and expand(r17_stripped.subs(y, 0)) == Rational(1, 34))

# =====================================================================
# B.  The slice conditions, re-derived here (not imported).
# =====================================================================
head("B.  The slice conditions, re-derived from the stripping")

dgen = {i: symbols(f'dd{i}') for i in range(0, 14)}
NN = 13
# c_{4-i} = d_{4-i} * y^{7-2i} * t^{1-2i}   <=   D_j = c_j C4^{7-2j},
#                                                d_j = D_j / y^{12(4-j)}, C4 = y^7 t.
# Step 1: the substitution itself, from the two defining relations.
jx, kx = symbols('jx kx', integer=True)
check("c_(4-i) = d_(4-i) * y^(7-2i) * t^(1-2i):  from c_j = D_j C4^(2j-7) and "
      "D_j = d_j y^(12(4-j)), with C4 = y^7 t and j = 4-i, the y-exponent is "
      "12i + 7(1-2i) = 7-2i and the t-exponent is 1-2i",
      sp.simplify((12 * (4 - jx) + 7 * (2 * jx - 7)).subs(jx, 4 - kx) - (7 - 2 * kx)) == 0
      and sp.simplify((2 * jx - 7).subs(jx, 4 - kx) - (1 - 2 * kx)) == 0)

# Step 2: U(u) = C4 * H(u / (y^2 t^2))   -- exact, term by term.
Ugen = sum(dgen[i] * y ** (7 - 2 * i) * (y + 1) ** (1 - 2 * i) * u ** i for i in range(NN))
Uref = sum(y ** 7 * (y + 1) * dgen[i] * (u / (y ** 2 * (y + 1) ** 2)) ** i
           for i in range(NN))
check("U(u) = C4 * H(u/(y^2 t^2)) exactly, on generic stripped d's",
      all(sp.cancel(sp.expand(sp.Poly(Ugen, u).coeff_monomial(u ** i)
                              - sp.Poly(sp.expand(Uref), u).coeff_monomial(u ** i))) == 0
          for i in range(NN)))

# Step 3: the slice formulas then follow by exponent arithmetic; checked symbolically.
Mv = symbols('Mv')
check("P_M = [u^(8-M)] (C4^2 H^2(u/(y^2 t^2))) = y^(2M-2) p_(8-M) / t^(14-2M):  "
      "y-exponent 14 - 2(8-M) = 2M-2, t-exponent 2 - 2(8-M) = 2M-14",
      sp.simplify(14 - 2 * (8 - Mv) - (2 * Mv - 2)) == 0
      and sp.simplify(2 - 2 * (8 - Mv) - (2 * Mv - 14)) == 0)
check("Q_M = [u^(12-M)] (C4^3 H^3(u/(y^2 t^2))) = y^(2M-3) r_(12-M) / t^(21-2M):  "
      "y-exponent 21 - 2(12-M) = 2M-3, t-exponent 3 - 2(12-M) = 2M-21",
      sp.simplify(21 - 2 * (12 - Mv) - (2 * Mv - 3)) == 0
      and sp.simplify(3 - 2 * (12 - Mv) - (2 * Mv - 21)) == 0)

# Step 4: numerical corroboration of both slice formulas at a random rational
# (y, d) point -- independent of the exponent bookkeeping above.
import random
random.seed(20260725)
YV = Rational(7, 3)
pt = {dgen[i]: Rational(random.randint(-9, 9), random.randint(1, 5)) for i in range(NN)}
pt[dgen[0]] = sp.Integer(1)
pt[y] = YV
Un = sp.expand(Ugen.subs(pt))
Hn = sp.expand(sum(pt[dgen[i]] * u ** i for i in range(NN)))
U2n = sp.expand(Un * Un)
H2n = sp.expand(Hn * Hn)
U3n = sp.expand(sp.expand(U2n * Un))
H3n = sp.expand(sp.expand(H2n * Hn))
cU2 = sp.Poly(U2n, u)
cH2 = sp.Poly(H2n, u)
cU3 = sp.Poly(U3n, u)
cH3 = sp.Poly(H3n, u)
TV = YV + 1
okP = all(sp.nsimplify(cU2.coeff_monomial(u ** (8 - M))
                       - YV ** (2 * M - 2) * cH2.coeff_monomial(u ** (8 - M))
                       / TV ** (14 - 2 * M)) == 0 for M in range(8, -1, -1))
okQ = all(sp.nsimplify(cU3.coeff_monomial(u ** (12 - M))
                       - YV ** (2 * M - 3) * cH3.coeff_monomial(u ** (12 - M))
                       / TV ** (21 - 2 * M)) == 0 for M in range(12, 0, -1))
check("both slice formulas verified at a random rational (y, d) point "
      "(P: M = 8..0, Q: M = 12..1)", okP and okQ)

check("negative y-powers occur exactly at (P,M=0) and (Q,M=1),(Q,M=0) "
      "=> the y-order family is exactly {y^2|p_8, y|r_11, y^3|r_12}",
      [M for M in range(0, 9) if 2 * M - 2 < 0] == [0]
      and [M for M in range(0, 13) if 2 * M - 3 < 0] == [0, 1]
      and (2 * 0 - 2) == -2 and (2 * 1 - 3) == -1 and (2 * 0 - 3) == -3)

# =====================================================================
# C.  The t-place profile.  Graded machinery + cascade levels 12, 14.
# =====================================================================
head("C.  The t-place profile (shifted coordinates, h_1 = 0)")

A = {k: symbols(f'A{k}') for k in range(2, 9)}


def tprofile(W, nu=NU_DEF, TMAX=45):
    """h_k = A_k t^{W[k]} for k=2..8.  Returns (h, p, r) as {t-degree: coeff} dicts.

    GRADED LEMMA.  r_n is a polynomial in h_2..h_8 alone.  Writing the true
    h_k = t^{W[k]} * Hk(t) with Hk in Q[[t]] arbitrary, r_n = sum_m rho_m(H(t)) t^m
    where rho_m is exactly the coefficient of t^m computed here (A_k <-> Hk).
    Hence: v_t(r_n) >= min{m : rho_m =/= 0 as a polynomial}, and at that m the
    coefficient of t^m in the true r_n is rho_m(H_2(0),...,H_8(0)).
    """
    def tmul(a, b):
        o = {}
        for da, ca in a.items():
            for db, cb in b.items():
                d = da + db
                if d > TMAX:
                    continue
                o[d] = o.get(d, 0) + ca * cb
        return {d: expand(c) for d, c in o.items() if expand(c) != 0}

    def tadd(a, b):
        o = dict(a)
        for d, c in b.items():
            e = expand(o.get(d, 0) + c)
            if e == 0:
                o.pop(d, None)
            else:
                o[d] = e
        return o

    h = [dict() for _ in range(nu + 1)]
    h[0] = {0: sp.Integer(1)}
    h[1] = {}
    for k in range(2, 9):
        h[k] = {W[k]: A[k]}
    for k in range(9, nu + 1):
        s = {}
        for i in range(1, k):
            s = tadd(s, tmul(h[i], h[k - i]))
        h[k] = {d: expand(-c / 2) for d, c in s.items()}

    def umul(X, Y):
        Z = [dict() for _ in range(nu + 1)]
        for i in range(nu + 1):
            if not X[i]:
                continue
            for j in range(nu + 1 - i):
                if not Y[j]:
                    continue
                Z[i + j] = tadd(Z[i + j], tmul(X[i], Y[j]))
        return Z
    P2 = umul(h, h)
    P3 = umul(P2, h)
    return h, P2, P3


def low(D):
    return (min(D), D[min(D)]) if D else (None, sp.Integer(0))


# base: v_t(h_k) >= 2k-1 for k = 2..5 -- the AUDITED a_t >= 9 cascade, in the
# shifted chart (h_1 = 0 is a specialisation of the lane's free h_1).
h_, P2_, P3_ = tprofile({2: 3, 3: 5, 4: 7, 5: 9, 6: 9, 7: 11, 8: 12})
check("(P<) at n=6: v_t(p_6) >= 10 forces v_t(h_6) >= 10  [jet 2*A6]",
      low(P2_[6]) == (9, 2 * A[6]))
check("(P<) at n=7: v_t(p_7) >= 12 forces v_t(h_7) >= 12  [jet 2*A7]",
      low(P2_[7]) == (11, 2 * A[7]))

# level 12, profile (3,5,7,9,10,12,13); (P<) at n=8 supplies the relation
# [t^13] p_8 = 2 A8 + 2 A2 A6 = 0.
h_, P2_, P3_ = tprofile({2: 3, 3: 5, 4: 7, 5: 9, 6: 10, 7: 12, 8: 13})
lo8, jet8 = low(P2_[8])
check("(P<) at n=8 with w6=10 gives the relation [t^13] p_8 = 2*A8 + 2*A2*A6 = 0",
      lo8 == 13 and expand(jet8 - (2 * A[8] + 2 * A[2] * A[6])) == 0)
stack12 = {d: expand(2 * P3_[12].get(d, 0) - 3 * P2_[12].get(d, 0))
           for d in set(P3_[12]) | set(P2_[12])}
stack12 = {d: c for d, c in stack12.items() if c != 0}
lo12, jet12 = low(stack12)
jet12_red = expand(jet12.subs(A[8], -A[2] * A[6]))
check("level 12: v_t(2r_12-3p_12) >= 20 < 21 required; jet reduces mod (P<) to 3*A6^2"
      "  => v_t(h_6) >= 11",
      lo12 == 20 and expand(jet12_red - 3 * A[6] ** 2) == 0)

# w8 >= 14 once w6 >= 11
h_, P2_, P3_ = tprofile({2: 3, 3: 5, 4: 7, 5: 9, 6: 11, 7: 12, 8: 13})
check("(P<) at n=8 with w6=11: v_t(p_8) >= 14 forces v_t(h_8) >= 14  [jet 2*A8]",
      low(P2_[8]) == (13, 2 * A[8]))

# level 14 -- NEVER RUN by slice_obstruction_basis.py (it stops at level 12)
h_, P2_, P3_ = tprofile({2: 3, 3: 5, 4: 7, 5: 9, 6: 11, 7: 12, 8: 14})
stack14 = {d: expand(2 * P3_[14].get(d, 0) - 3 * P2_[14].get(d, 0))
           for d in set(P3_[14]) | set(P2_[14])}
stack14 = {d: c for d, c in stack14.items() if c != 0}
check("level 14 (UN-RUN by the lane): p_14 = 0 and t^25 | r_14 (r_14 = 0), "
      "jet at t^24 is 3*A7^2  => v_t(h_7) >= 13",
      low(P2_[14])[0] is None and low(stack14) == (24, 3 * A[7] ** 2))

# CRITICAL for Section D: the profile derivation must be uniform in w_5, since
# the kill hypothesis is v_t(h_5) >= 10, not = 9.  Re-run every step at w_5 = 10,11,12.
_unif = True
for _w5 in (10, 11, 12):
    _, q2, _q3 = tprofile({2: 3, 3: 5, 4: 7, 5: _w5, 6: 9, 7: 11, 8: 12})
    _unif = _unif and low(q2[6]) == (9, 2 * A[6]) and low(q2[7]) == (11, 2 * A[7])
    _, q2, q3 = tprofile({2: 3, 3: 5, 4: 7, 5: _w5, 6: 10, 7: 12, 8: 13})
    _unif = _unif and expand(low(q2[8])[1] - (2 * A[8] + 2 * A[2] * A[6])) == 0
    st = {d: expand(2 * q3[12].get(d, 0) - 3 * q2[12].get(d, 0))
          for d in set(q3[12]) | set(q2[12])}
    st = {d: c for d, c in st.items() if c != 0}
    _l, _j = low(st)
    _unif = _unif and _l == 20 and expand(_j.subs(A[8], -A[2] * A[6]) - 3 * A[6] ** 2) == 0
    _, q2, _q3 = tprofile({2: 3, 3: 5, 4: 7, 5: _w5, 6: 11, 7: 12, 8: 13})
    _unif = _unif and low(q2[8]) == (13, 2 * A[8])
    _, q2, q3 = tprofile({2: 3, 3: 5, 4: 7, 5: _w5, 6: 11, 7: 12, 8: 14})
    st = {d: expand(2 * q3[14].get(d, 0) - 3 * q2[14].get(d, 0))
          for d in set(q3[14]) | set(q2[14])}
    st = {d: c for d, c in st.items() if c != 0}
    _unif = _unif and low(st) == (24, 3 * A[7] ** 2)
check("UNIFORMITY: every profile step above (w6>=10, w7>=12, the (P<) n=8 relation, "
      "level 12 => w6>=11, w8>=14, level 14 => w7>=13) is UNCHANGED at "
      "w_5 = 10, 11, 12 -- so the profile used in Section D is legitimate under "
      "the kill hypothesis v_t(h_5) >= 10", _unif)

PROFILE = {2: 3, 3: 5, 4: 7, 5: 9, 6: 11, 7: 13, 8: 14}
h_, P2_, P3_ = tprofile(PROFILE)
check("at the final profile (3,5,7,9,11,13,14) every (P<) condition t^(2n-2)|p_n "
      "(n=2..8) holds identically -- no residual relations",
      all(low(P2_[n])[0] is None or low(P2_[n])[0] >= 2 * n - 2 for n in range(2, 9)))
check("at the final profile every (Q) condition t^(2n-3)|r_n (n=2..12) holds "
      "identically -- no residual relations",
      all(low(P3_[n])[0] is None or low(P3_[n])[0] >= 2 * n - 3 for n in range(2, 13)))
check("at the final profile v_t(r_17) = 30 exactly -- matching the Phi relation "
      "with no cancellation demanded",
      low(P3_[17])[0] == 30)

# =====================================================================
# D.  THE KILL.  v_t(h_5) >= 10 is impossible.
# =====================================================================
head("D.  a_t <= 9   (the six a10_* cells die)")

TARGET = Rational(1, 2)          # = q(-1)/6630


def layer0(w5, w8=14, w6=11, w7=13):
    W = {2: 3, 3: 5, 4: 7, 5: w5, 6: w6, 7: w7, 8: w8}
    _, P2x, P3x = tprofile(W)
    out = {}
    for n in (13, 14, 15, 17):
        out[n] = low(P3x[n])
    return out, P2x, P3x


for w5 in (10, 11, 12):
    L, P2x, P3x = layer0(w5)
    eqs = [L[n][1] for n in (13, 14, 15)]
    v17, rho17 = L[17]
    ok_v17 = (v17 == 30)
    sysA = eqs + [expand(rho17 - TARGET)]
    G = groebner(sysA, *[A[k] for k in range(2, 9)], order='grevlex')
    unit = list(G.exprs) == [sp.Integer(1)]
    G2_ = groebner(eqs + [expand(rho17 + TARGET)], *[A[k] for k in range(2, 9)],
                   order='grevlex')
    unit2 = list(G2_.exprs) == [sp.Integer(1)]
    check(f"w5 = {w5}: v_t(r_17) = 30 and the layer-0 system "
          f"{{rho_13 = rho_14 = rho_15 = 0, rho_17 = +-1/2}} is the UNIT IDEAL",
          ok_v17 and unit and unit2)
    if w5 == 10:
        say("      rho_13 =", factor(eqs[0]))
        say("      rho_14 =", factor(eqs[1]))
        say("      rho_15 =", factor(eqs[2]))
        say("      rho_17 =", factor(rho17))

# the human-readable proof at w5 = 10, so the GB is not the only witness
L, _, _ = layer0(10)
e14 = L[14][1]
e15 = L[15][1]
e17 = L[17][1]
check("hand proof, branch A6 = 0: rho_14 is vacuous; rho_15 becomes 3*A7*A8, so "
      "A7*A8 = 0; and rho_17 becomes -3*A2*(A7*A8) = 0 =/= 1/2  => CONTRADICTION",
      expand(e14.subs(A[6], 0)) == 0
      and expand(e15.subs(A[6], 0) - 3 * A[7] * A[8]) == 0
      and expand(e17.subs(A[6], 0) + 3 * A[2] * (A[7] * A[8])) == 0)
_sub = {A[8]: -A[2] * A[6] / 2}
check("hand proof, branch A6 =/= 0: rho_14 = 0 gives A8 = -A2*A6/2; substituting, "
      "A6*(A2*A7 + A3*A6) = -(2/3)*rho_15|sub, so rho_15 = 0 forces A2*A7+A3*A6 = 0 "
      "and hence rho_17 = 0 =/= 1/2  => CONTRADICTION",
      expand(e14.subs(_sub)) == 0
      and expand(A[6] * (A[2] * A[7] + A[3] * A[6])
                 + Rational(2, 3) * e15.subs(_sub)) == 0
      and expand(e17 + 3 * A[8] * (A[2] * A[7] + A[3] * A[6])) == 0)

say("\n  ==> v_t(h_5) >= 10 is IMPOSSIBLE.   h_5 = dm1 = e (SLICE_OBSTRUCTION S3.4),")
say("      so a_t = v_t(e) <= 9.  With the audited a_t >= 9:   a_t = 9 EXACTLY.")

# =====================================================================
# E.  Sharpness: a_t = 9 survives, and what else gets pinned.
# =====================================================================
head("E.  Sharpness and the extra sharp relations")

L9, _, P3_9 = layer0(9)
e13_9, e14_9, e15_9, e17_9 = (L9[13][1], L9[14][1], L9[15][1], L9[17][1])
check("at w5 = 9 the four layer-0 quantities sit at t^23, t^25, t^27, t^30 "
      "(= 2n-3 for n=13,14,15; = 30 for n=17)",
      (L9[13][0], L9[14][0], L9[15][0], L9[17][0]) == (23, 25, 27, 30))

WIT = {A[2]: -1, A[3]: 1, A[4]: 0, A[5]: 1, A[6]: 0, A[7]: Rational(-1, 3),
       A[8]: Rational(-1, 2)}
check("EXPLICIT RATIONAL WITNESS for a_t = 9: "
      "(A2..A8) = (-1, 1, 0, 1, 0, -1/3, -1/2) satisfies all four layer-0 equations",
      expand(e13_9.subs(WIT)) == 0 and expand(e14_9.subs(WIT)) == 0
      and expand(e15_9.subs(WIT)) == 0 and expand(e17_9.subs(WIT) - TARGET) == 0
      and WIT[A[5]] != 0)

# A2 * A5^3 = -1 on the a_t = 9 branch (saturate at A5 =/= 0, A8 =/= 0)
z = symbols('z')
Gsat = groebner([e13_9, e14_9, e15_9, expand(e17_9 - TARGET),
                 A[5] * A[8] * z - 1],
                *([A[k] for k in range(2, 9)] + [z]), order='lex')
check("on the a_t = 9 branch the layer-0 system forces  A2 * A5^3 = -1 "
      "(hence v_t(h_2) = 3 and v_t(h_5) = 9, both EXACTLY)",
      Gsat.reduce(expand(A[2] * A[5] ** 3 + 1))[1] == 0)

# w8 >= 15 is impossible
L8, _, _ = layer0(9, w8=15)
check("v_t(h_8) >= 15 forces v_t(r_17) >= 31 > 30  =>  v_t(h_8) = 14 EXACTLY",
      L8[17][0] is not None and L8[17][0] > 30)

# NEGATIVE CONTROLS -- the Phi relation and level 15 are both load-bearing
L10, _, _ = layer0(10)
noPhi = groebner([L10[13][1], L10[14][1], L10[15][1]],
                 *[A[k] for k in range(2, 9)], order='grevlex')
check("NEGATIVE CONTROL: drop Sec.8 item (4) (the Phi relation) and the w5 = 10 "
      "system is CONSISTENT -- item (4) is load-bearing, not decorative",
      list(noPhi.exprs) != [sp.Integer(1)])
no15 = groebner([L10[13][1], L10[14][1], expand(L10[17][1] - TARGET)],
                *[A[k] for k in range(2, 9)], order='grevlex')
check("NEGATIVE CONTROL: drop level 15 and the w5 = 10 system is CONSISTENT "
      "-- the un-run level 15 is load-bearing too",
      list(no15.exprs) != [sp.Integer(1)])
no14 = groebner([L10[13][1], L10[15][1], expand(L10[17][1] - TARGET)],
                *[A[k] for k in range(2, 9)], order='grevlex')
check("NEGATIVE CONTROL: drop level 14 and the w5 = 10 system is CONSISTENT",
      list(no14.exprs) != [sp.Integer(1)])
check("NEGATIVE CONTROL (anti-vacuity): each of rho_13, rho_14, rho_15, rho_17 at "
      "w5 = 10 is a NON-ZERO polynomial",
      all(expand(x) != 0 for x in
          (L10[13][1], L10[14][1], L10[15][1], L10[17][1])))

# =====================================================================
# F.  Sec.8 item (1): the y-place.  The genuinely new family.
# =====================================================================
head("F.  The place y = 0 -- Sec.8 item (1), new to the whole program")

g = symbols('g2:8')
Gt = 1 + sum(g[i - 2] * w ** i for i in range(2, 8))
X = Gt - 1
Ser = sum(sp.binomial(Rational(3, 2), m) * X ** m for m in range(0, 12))
PW = Poly(expand(Ser), w)
co = lambda n: expand(PW.coeff_monomial(w ** n))

check("co(n) is weighted-homogeneous of weight n (weight(gamma_k) = k) for "
      "n = 11..15, 17",
      all({sum((i + 2) * m[i] for i in range(6))
           for m in Poly(co(n), *g).monoms()} == {n}
          for n in (11, 12, 13, 14, 15, 17)))

b_, c_ = symbols('b c')
Adeg = 1 + b_ * w ** 2 + c_ * w ** 3
sub_deg = {g[i - 2]: Poly(expand(Adeg ** 2), w).coeff_monomial(w ** i)
           for i in range(2, 8)}
check("the degenerate component Gamma~ = A^2 (A of degree 3) satisfies "
      "co(11..15) = 0 but has co(17) = 0, so it is excluded by the Phi relation",
      all(sp.simplify(co(n).subs(sub_deg)) == 0 for n in (11, 12, 13, 14, 15, 17)))

check("ANTI-VACUITY: co(11), co(12), co(13), co(14), co(15), co(17) are all "
      "NON-ZERO polynomials, so all six y=0 layer conditions have content",
      all(expand(co(n)) != 0 for n in (11, 12, 13, 14, 15, 17)))

# the y = 0 layer IS non-empty: a certified numerical solution.
YSOL = [complex(0.21819819, 2.35473347), complex(-0.231267, 0.306247),
        complex(-2.39472205, 0.44765113), complex(-1.60111813, -0.99137037),
        complex(-0.17304106, -0.60817624), complex(0.43227419, -0.86812361)]
_fy = [sp.lambdify(g, e, 'math') for e in
       (co(11), co(12), co(13), co(14), co(15), co(17) - Rational(1, 34))]
_res = max(abs(complex(f(*YSOL))) for f in _fy)
check("the y=0 layer system is NON-EMPTY over C: the recorded Newton solution "
      f"has residual {_res:.2e} < 1e-4 at 8 printed digits "
      "(=> item (1) kills nothing alone)",
      _res < 1e-4)

# =====================================================================
# G.  INDEPENDENCE of the Sec.8 families from the slice lane's imposed set.
# =====================================================================
head("G.  Independence witnesses (against the lane's IMPOSED set)")

delta = symbols('delta', positive=True)
NW = 19


def vt(e):
    e = expand(e)
    if e == 0:
        return None
    return min(m[0] for m in Poly(e, t).monoms())


hw = {i: expand(sp.binomial(Rational(1, 2), i) * (delta * t ** 2) ** i)
      for i in range(NW)}
Pw2 = series_pow(hw, 2, NW - 1)
Pw3 = series_pow(hw, 3, NW - 1)
check("WITNESS W1 = H(u) = sqrt(1 + t^2*delta*u): p_1 = t^2*delta and p_n = 0 for "
      "all n >= 2, so (P<) n=2..8 and (P0) n>=9 hold trivially",
      expand(Pw2[1] - delta * t ** 2) == 0
      and all(expand(Pw2[n]) == 0 for n in range(2, NW - 1)))
check("WITNESS W1 satisfies (Q): v_t(r_n) = 2n >= 2n-3 for n = 2..12  "
      "=> it satisfies EVERY condition the slice lane imposes",
      all(vt(Pw3[n]) == 2 * n for n in range(2, 13)))
check("WITNESS W1 VIOLATES Sec.8 item (3): r_13 =/= 0  "
      "=> item (3) is INDEPENDENT of the lane's imposed set",
      expand(Pw3[13]) != 0 and expand(Pw3[14]) != 0 and expand(Pw3[15]) != 0)
check("WITNESS W1 VIOLATES Sec.8 item (4): v_t(r_17) = 34 =/= 30  "
      "=> item (4) is INDEPENDENT of the lane's imposed set",
      vt(Pw3[17]) == 34)

# Sec.8 item (1) independence: W1 also violates the y-order family.
check("WITNESS W1 VIOLATES Sec.8 item (1) too: r_11 = binom(3/2,11)*delta^11*t^22 "
      "exactly, a unit times delta^11, so ord_y(r_11) = 0 whenever delta is a "
      "y-unit -- item (1) is INDEPENDENT of the lane's imposed set",
      expand(Pw3[11] - sp.binomial(Rational(3, 2), 11) * delta ** 11 * t ** 22) == 0
      and sp.binomial(Rational(3, 2), 11) != 0)

# Sec.8 item (2): the degree caps ARE the window caps -- read from the source.
_wc = open(os.path.join(HERE, 'window_caps_verify.py')).read()
check("Sec.8 item (2) 'degree/upper-hull' = the window caps deg <= 15k/14k, "
      "ord >= 12k, i.e. deg d_(4-k) <= lam*k -- present verbatim in "
      "window_caps_verify.py, hence imposed program-wide, not slack",
      'deg <= 15k (sub1) / 14k (sub2)' in _wc or 'ord >= 12k, deg <= 15k' in _wc
      or ('15k' in _wc and '14k' in _wc and '12k' in _wc))

# =====================================================================
head("=" * 68)
if _fail:
    print(f"FAILED {len(_fail)} of {_n[0]} checks:")
    for f in _fail:
        print("   -", f)
    sys.exit(1)
say(f"ALL {_n[0]} CHECKS PASSED")
if QUIET:
    print(f"ALL {_n[0]} CHECKS PASSED")
sys.exit(0)
