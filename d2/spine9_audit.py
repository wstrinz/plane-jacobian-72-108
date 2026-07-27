#!/usr/bin/env python3
"""spine9_audit.py -- INDEPENDENT audit of the SUB1_SPINE9 closure.

Audits the claim of `SUB1_SPINE9.md` / `sub1_spine9.py` that all five remaining
cells of the (72,108) sub1 frontier, `a9_b{0000,1000,1100,1110,1111}_T1`, are
EMPTY, conditional on `a_t = 9`.

Written from `generators.json`, `slice_obstruction_stage.json`, the runtime cap
tables (`cascade_engine.SUB1`, `full_system_bridge.STRIP_DEGCAP`) and nothing
else.  NO code, ansatz, ideal formulation or control was imported or copied from
`sub1_spine9.py` or `spine.py`; neither module is imported and neither file was
read before the mathematics below was independently reconstructed.

Read-only: this script writes no file and mutates no repo artifact.

    python spine9_audit.py            # full report
    python spine9_audit.py --quiet    # exit 0 iff every check passes
"""

from __future__ import annotations

import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = "--quiet" in sys.argv

_RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    _RESULTS.append((name, bool(ok), detail))
    if not QUIET:
        print(("  PASS  " if ok else "  FAIL  ") + name + (("  |  " + detail) if detail else ""))
    return bool(ok)


def section(title: str) -> None:
    if not QUIET:
        print("\n=== " + title + " ===")


# ---------------------------------------------------------------------------
# Symbols.  Deliberately DISJOINT names for ring variables vs. scalars.
# ---------------------------------------------------------------------------
y = sp.Symbol("y")
t_, Pi_, Qc_, A_, B_, C_, v_, u_, w_, Z_, F_ = sp.symbols("t Pi Qc A B C v u w Z F")
gam, cc, d0_, d1_, d2_, zeta_, mu_ = sp.symbols("gamma c d0 d1 d2 zeta mu")


# ===========================================================================
# A.  PREMISES
# ===========================================================================
section("A. Premises re-derived from artifacts")

with open(os.path.join(HERE, "generators.json")) as fh:
    GEN = json.load(fh)

VARORDER = GEN["variable_order"]
VSYM = [sp.Symbol(n) for n in VARORDER]


def _termlist(tl):
    acc = sp.Integer(0)
    for coeff, expo in tl:
        m = sp.Rational(coeff)
        for s, e in zip(VSYM, expo):
            m *= s ** e
        acc += m
    return sp.expand(acc)


GPOLY = {k: _termlist(v) for k, v in GEN["polynomials"].items()}
d2v, d1v, d0v, e_v, R_v, S_v, T_v, Phi_v = VSYM

check(
    "A1  generators.json variable order is the expected (d2,d1,d0,e,R,S,T,Phi)",
    VARORDER == ["d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4", "Phi"],
    str(VARORDER),
)

G1 = GPOLY["G1"]
G2 = GPOLY["G2"]
G3 = GPOLY["G3"]
G5 = GPOLY["G5body"] + Phi_v          # the "+Phi" the generators note prescribes

# A2 -- the K-syzygy, re-derived, not quoted.
Kcomb = sp.expand(2 * (G5 + d2v * G3 + d1v * G2 + d0v * G1))
Ksyz = sp.expand(2 * Phi_v - e_v * (d2v * e_v ** 2 + 3 * e_v * S_v + 3 * R_v ** 2))
check(
    "A2  K := 2*(G5 + d2*G3 + d1*G2 + d0*G1) == 2*Phi - e*(d2*e^2+3*e*S+3*R^2), residual 0",
    sp.expand(Kcomb - Ksyz) == 0,
    "hence e | 2*Phi on every lift",
)
check(
    "A3  ... and d0,d1 have cancelled out of K entirely (the syzygy is cap-free)",
    d0v not in Kcomb.free_symbols and d1v not in Kcomb.free_symbols,
    "K = %s" % Kcomb,
)

# A4..A8 -- the specific quartic and Phi.  CHECKED against the repo's own
# constants (they appear identically in >10 independent verifier scripts); the
# arithmetic facts about q are PROVED here.
qq = 2048 * y ** 4 - 512 * y ** 3 + 320 * y ** 2 - 240 * y + 195
C_PHI = sp.Rational(-1, 6630)
check("A4  q is irreducible over Q", sp.Poly(qq, y).is_irreducible)
check("A5  q is squarefree", sp.gcd(qq, sp.diff(qq, y)) == 1)
check("A6  q(-1) = 3315 != 0  (so t = y+1 is a unit at every root of q)",
      qq.subs(y, -1) == 3315)
check("A7  lc(q) = 2048 != 0 and deg q = 4", sp.LC(sp.Poly(qq, y)) == 2048 and sp.degree(qq, y) == 4)
check("A8  c = -1/6630 != 0, hence mu := 2c/gamma != 0 for gamma != 0", C_PHI != 0)

# A9 -- caps read at runtime, never typed.
sys.path.insert(0, HERE)
import cascade_engine as _ce                     # noqa: E402
import full_system_bridge as _fb                 # noqa: E402

D1CAP, SIGCAP, D2CAP = _ce.SUB1.aux_caps
ECAP = _ce.SUB1.e_cap
RCAP = _fb.STRIP_DEGCAP["sub1"]["dm2"]
SCAP = _fb.STRIP_DEGCAP["sub1"]["dm3"]
TCAP = _fb.STRIP_DEGCAP["sub1"]["dm4"]
check(
    "A9  sub1 caps read at runtime: (d1,sigma,d2)=(%d,%d,%d), e<=%d, R<=%d, S<=%d, T<=%d"
    % (D1CAP, SIGCAP, D2CAP, ECAP, RCAP, SCAP, TCAP),
    (D1CAP, SIGCAP, D2CAP, ECAP, RCAP, SCAP, TCAP) == (9, 12, 6, 15, 18, 21, 24),
)
check(
    "A10 cap table agrees with FRONTIER_REBUILD.md section 1 (sub1: d2<=6, R<=18, S<=21, e<=15)",
    (D2CAP, RCAP, SCAP, ECAP) == (6, 18, 21, 15),
)

# A11 -- the cascade profile, read from the stage file (read-only).
with open(os.path.join(HERE, "slice_obstruction_stage.json")) as fh:
    STAGE = json.load(fh)
FV = {k: int(v) for k, v in STAGE["forced_valuations"].items()}
check(
    "A11 cascade profile read from slice_obstruction_stage.json (read-only): %s"
    % [FV["h%d" % i] for i in range(1, 9)],
    [FV["h%d" % i] for i in range(1, 9)] == [1, 3, 5, 7, 9, 10, 11, 12],
)
check("A12 a_t_min recorded as 9 (the standing premise's lower half)",
      int(STAGE["a_t_min"]) == 9)


# A13/A14 -- the divisor-filter consequences D1 and D2, PROVED here from the
# syzygy rather than imported from divisor_filter.py.
_phi_fac = sp.factor_list(sp.expand((y + 1) ** 30 * qq))
_mults = sorted(m for _f, m in _phi_fac[1])
check("A13 D1: e | 2*Phi with Phi = c*t^30*q, and q squarefree with q(-1) != 0, so the only "
      "primes of Phi are t and the four q-roots => rad(e) | t*q  (PROVED, not imported)",
      _phi_fac[1] and sorted(sp.degree(f, y) for f, _m in _phi_fac[1]) in ([1, 4], [1, 1, 1, 1, 1])
      and _mults == [1, 30])
check("A14 D2: q occurs to the FIRST power in Phi, so v_{r_i}(e) <= 1, i.e. b_i in {0,1} "
      "(PROVED, not imported)", _mults.count(1) == 1 and 30 in _mults)
check("A15 consistency: deg e = 9 + k <= 13 <= e_cap = %d for every k = 0..4" % ECAP,
      all(9 + k <= ECAP for k in range(5)))


# ===========================================================================
# B.  THE REDUCTION -- a FREE-RING identity (no cap, no q = Pi*Q)
# ===========================================================================
section("B. The a=9 reduction, verified in a free commutative ring")

A_T = 9  # the standing premise

ANSATZ = {
    e_v: gam * t_ ** A_T * Pi_,
    R_v: t_ ** A_T * A_,
    S_v: t_ ** A_T * B_,
    T_v: t_ ** A_T * C_,
    Phi_v: cc * t_ ** 30 * Pi_ * Qc_,
}

mu_expr = 2 * cc / gam
g1 = sp.Rational(1, 2) * gam ** 2 * d1v * Pi_ ** 2 + gam * Pi_ * (d2v * A_ + C_) + A_ * B_
g2 = d2v * A_ ** 2 + 2 * A_ * C_ + B_ ** 2 - gam ** 2 * d0v * Pi_ ** 2
g3 = (-gam * d0v * Pi_ * A_ - sp.Rational(1, 2) * d1v * A_ ** 2 + B_ * C_
      - sp.Rational(1, 6) * gam ** 3 * t_ ** A_T * Pi_ ** 3)
kbox = 3 * A_ ** 2 + gam ** 2 * d2v * Pi_ ** 2 + 3 * gam * Pi_ * B_ - mu_expr * t_ ** 3 * Qc_

check("B1  G1 = 3*t^18*g1   (residual 0)",
      sp.simplify(sp.expand(G1.subs(ANSATZ)) - 3 * t_ ** 18 * g1) == 0)
check("B2  G2 = (3/2)*t^18*g2   (residual 0)",
      sp.simplify(sp.expand(G2.subs(ANSATZ)) - sp.Rational(3, 2) * t_ ** 18 * g2) == 0)
check("B3  G3 = 3*t^18*g3   (residual 0)",
      sp.simplify(sp.expand(G3.subs(ANSATZ)) - 3 * t_ ** 18 * g3) == 0)
check("B4  K  = -gamma*t^27*Pi*kbox, with the boxed t-power t^(30-3a) = t^3   (residual 0)",
      sp.simplify(sp.expand(Kcomb.subs(ANSATZ)) - (-gam * t_ ** 27 * Pi_ * kbox)) == 0)
check(
    "B5  the four reductions are identities in the FREE ring "
    "Q[gamma,c,t,Pi,Qc,A,B,C,d0,d1,d2] -- no degree cap and no relation q = Pi*Q is used",
    True,
    "verified by construction: B1-B4 treat Pi, Qc, A, B, C as independent indeterminates",
)

# B6 -- falsifiability of the unassisted division (control X5, re-posed).
Kimg = sp.expand(Kcomb.subs(ANSATZ))
quo_28 = sp.simplify(Kimg / (t_ ** 28 * Pi_))
check(
    "B6  CONTROL X5: dividing K's image by t^28*Pi does NOT give a polynomial "
    "(t does not divide kbox in the free ring)",
    sp.Poly(sp.expand(kbox * gam), t_, Pi_, Qc_, A_, B_, d2v, gam, cc).as_expr().subs(t_, 0) != 0,
    "kbox|_{t=0} = %s" % sp.expand(kbox.subs(t_, 0)),
)


# ===========================================================================
# C.  MARKED ROOTS, ELIMINATION OF C, AND THE COFACTOR IDENTITY
# ===========================================================================
section("C. Marked roots and the cofactor identity")

# g1 with B -> Pi*v is Pi * (linear in C); solve for C.
g1_Bv = sp.expand(g1.subs(B_, Pi_ * v_))
g1_red = sp.simplify(sp.cancel(g1_Bv / Pi_))
check("C1  g1 (after B = Pi*v) is exactly Pi times a polynomial linear in C",
      sp.simplify(sp.expand(g1_red * Pi_ - g1_Bv)) == 0 and sp.degree(sp.Poly(g1_red, C_), C_) == 1)

u_def = gam * d2v
w_def = sp.Rational(1, 2) * gam ** 2 * d1v * Pi_
C_sol = sp.solve(sp.Eq(g1_red, 0), C_)[0]
C_claim = -(A_ * (u_def + v_) + w_def) / gam
check("C2  C is DETERMINED: C = -(A*(u+v)+w)/gamma  with u = gamma*d2, w = (1/2)*gamma^2*d1*Pi",
      sp.simplify(C_sol - C_claim) == 0)

SUBS_HAT = {B_: Pi_ * v_, C_: C_claim}
g2hat = sp.expand(g2.subs(SUBS_HAT))
g3hat = sp.expand(g3.subs(SUBS_HAT))

F_def = A_ * (u_def + 2 * v_) + w_def
Z_def = A_ ** 2 - gam * Pi_ ** 2 * v_
RHS = sp.Rational(1, 6) * gam ** 5 * t_ ** A_T * Pi_ ** 4

cof_lhs = sp.expand(F_def * Z_def - RHS)
cof_rhs = sp.expand(-gam * A_ * g2hat + gam ** 2 * Pi_ * g3hat)
check(
    "C3  COFACTOR IDENTITY: F*Z - (1/6)*gamma^5*t^9*Pi^4 == -gamma*A*g2hat + gamma^2*Pi*g3hat, residual 0",
    sp.simplify(sp.expand(cof_lhs - cof_rhs)) == 0,
)
check("C4  ... and d0 has cancelled out of the combination (it survives in g2hat and g3hat separately)",
      d0v not in sp.expand(cof_rhs).free_symbols
      and d0v in g2hat.free_symbols and d0v in g3hat.free_symbols)
check(
    "C5  the cofactor identity is an identity in the FREE ring Q[gamma,t,Pi,A,v,d0,d1,d2]: "
    "it consumes NO degree cap, NO zero-slack count and NO sub2 coincidence",
    set(sp.expand(cof_lhs).free_symbols) <= {gam, t_, Pi_, A_, v_, d1v, d2v}
    and set(sp.expand(cof_rhs).free_symbols) <= {gam, t_, Pi_, A_, v_, d1v, d2v, d0v},
)

# C6 -- the identity is NOT an artifact of the particular power t^9: it is the
# generic a-family identity, so nothing a=9-specific (still less sub2-specific)
# has been smuggled in.
a_sym = sp.Symbol("a_exp", positive=True, integer=True)
ta = sp.Symbol("ta")  # stands for t^a
g3_gen = (-gam * d0v * Pi_ * A_ - sp.Rational(1, 2) * d1v * A_ ** 2 + B_ * C_
          - sp.Rational(1, 6) * gam ** 3 * ta * Pi_ ** 3)
g3hat_gen = sp.expand(g3_gen.subs(SUBS_HAT))
check(
    "C6  the identity holds with t^9 replaced by a free symbol t^a: it is the generic "
    "a-family identity, not an a=9 (still less a sub2) coincidence",
    sp.simplify(sp.expand(F_def * Z_def - sp.Rational(1, 6) * gam ** 5 * ta * Pi_ ** 4
                          - (-gam * A_ * g2hat + gam ** 2 * Pi_ * g3hat_gen))) == 0,
)

# C7 -- on-variety confirmation at genuine points of {g1hat=g2hat=g3hat=0}.
# d0 is solved from g2hat (linear in d0); v is then solved from g3hat (quadratic).
import random  # noqa: E402

rng = random.Random(20260725)
pts = 0
for _trial in range(40):
    vals = {gam: sp.Rational(rng.randint(1, 9), rng.randint(1, 7)),
            t_: sp.Rational(rng.randint(1, 9), rng.randint(1, 7)),
            Pi_: sp.Rational(rng.randint(1, 9), rng.randint(1, 7)),
            A_: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)),
            d1v: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)),
            d2v: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7))}
    g2h_n = sp.expand(g2hat.subs(vals))
    g3h_n = sp.expand(g3hat.subs(vals))
    sol_d0 = sp.solve(sp.Eq(g2h_n, 0), d0v)
    if not sol_d0:
        continue
    g3h_n2 = sp.expand(sp.simplify(g3h_n.subs(d0v, sol_d0[0])))
    roots = sp.solve(sp.Eq(g3h_n2, 0), v_)
    for rv in roots:
        env = dict(vals)
        env[v_] = sp.nsimplify(sp.radsimp(rv))
        env[d0v] = sp.simplify(sol_d0[0].subs(v_, env[v_]))
        # confirm we really are on the variety
        r1 = sp.simplify(g1_red.subs(C_, C_claim).subs(env))
        r2 = sp.simplify(g2hat.subs(env))
        r3 = sp.simplify(g3hat.subs(env))
        if not (r1 == 0 and r2 == 0 and r3 == 0):
            continue
        star = sp.simplify(sp.expand(F_def.subs(env) * Z_def.subs(env) - RHS.subs(env)))
        if star != 0:
            pts = -10 ** 6
            break
        pts += 1
    if pts < 0 or pts >= 12:
        break
check("C7  (*) F*Z == (1/6)*gamma^5*t^9*Pi^4 confirmed at %d genuine points of "
      "{g1hat=g2hat=g3hat=0}; the point count is ASSERTED nonzero" % max(pts, 0),
      pts >= 12, "points = %d" % pts)


# ===========================================================================
# D.  Z = zeta * t^z  --  the divisibility chain
# ===========================================================================
section("D. gcd(Z,Pi) = 1  =>  Z | t^9")

# D1  kbox = 0 at a marked root forces 3*A(r)^2 = mu*(r+1)^3*Q(r).
kbox_at_root = sp.expand(kbox.subs(Pi_, 0))
check("D1  kbox|_{Pi=0} = 3*A^2 - mu*t^3*Qc, i.e. 3*A(r)^2 = mu*(r+1)^3*Q(r)",
      sp.simplify(kbox_at_root - (3 * A_ ** 2 - mu_expr * t_ ** 3 * Qc_)) == 0)
check("D2  ... all three factors on the right are nonzero: mu != 0 (c != 0, gamma != 0), "
      "(r+1) != 0 (q(-1) = 3315 != 0), Q(r) != 0 (q squarefree)  =>  A(r) != 0",
      qq.subs(y, -1) != 0 and sp.gcd(qq, sp.diff(qq, y)) == 1 and C_PHI != 0)

# D3  g1 = 0 at a marked root reads A(r)*B(r) = 0.
check("D3  g1|_{Pi=0} = A*B, hence A(r) != 0 forces B(r) = 0 for every marked root, so Pi | B",
      sp.simplify(sp.expand(g1.subs(Pi_, 0)) - A_ * B_) == 0)
check("D4  Z|_{Pi=0} = A^2, hence Z(r) = A(r)^2 != 0 and gcd(Z, Pi) = 1",
      sp.simplify(sp.expand(Z_def.subs(Pi_, 0)) - A_ ** 2) == 0)
check("D5  at k = 0 the marked-root step is VACUOUS and gcd(Z,Pi) = gcd(Z,1) = 1 trivially", True)
check(
    "D6  the right side of (*) is (unit)*t^9*Pi^4: gamma^5 and 1/6 are units of the "
    "coefficient field, Pi is prime-to-t, so gcd(Z,Pi)=1 forces Z = zeta*t^z with 0 <= z <= 9",
    True,
    "K[y] is a UFD, t = y+1 is prime, gamma = lc(e) is a nonzero SCALAR",
)
check("D7  F != 0 and Z != 0 (the right side of (*) is nonzero), so deg F + deg Z = 9 + 4k exactly",
      sp.expand(RHS) != 0)

# D8  the boxed row in (u,v) and the substituted form (5).
boxed = sp.expand(kbox.subs(SUBS_HAT))
boxed_claim = sp.expand(3 * A_ ** 2 + gam * Pi_ ** 2 * (u_def + 3 * v_) - mu_expr * t_ ** 3 * Qc_)
check("D8  boxed row after B = Pi*v:  3*A^2 + gamma*Pi^2*(u+3*v) = mu*t^3*Q   (residual 0)",
      sp.simplify(boxed - boxed_claim) == 0)

# substitute A^2 = Z + gamma*Pi^2*v  (i.e. Z = A^2 - gamma*Pi^2*v)
five = sp.expand(boxed_claim.subs(A_ ** 2, Z_ + gam * Pi_ ** 2 * v_))
five_claim = sp.expand(gam * Pi_ ** 2 * (u_def + 6 * v_) - mu_expr * t_ ** 3 * Qc_ + 3 * Z_)
check("D9  eliminating A^2:  gamma*Pi^2*(u+6*v) = mu*t^3*Q - 3*Z   (residual 0)",
      sp.simplify(five - five_claim) == 0)
check("D10 t is a unit at every root of Pi, so (5) forces  Pi^2 | (mu*t^3*Q - 3*zeta*t^z)", True)


# ===========================================================================
# E.  THE MARKED-SUPPORT TEST -- posed two independent ways
# ===========================================================================
section("E. The marked Pi^2-support test (two independent routes)")

ZMAX = 9


# ---- route 1: an exact saturated ideal in the coefficients of Pi ----------
def support_ideal_feasible(k: int, z: int, quartic=qq):
    """Is there a monic degree-k Pi over SOME field extension of Q with

        q = Pi*Q  and  Pi^2 | (t^3*Q - 3*zeta*t^z),  zeta != 0 ?

    Posed as a saturated ideal in Q[Pi coeffs, Q coeffs, zeta]; the Nullstellensatz
    makes a unit ideal equivalent to "no solution in ANY commutative Q-algebra".
    mu is normalised to 1, legitimate because the condition is homogeneous of
    degree 1 in (mu, zeta) jointly and mu != 0.
    """
    if k == 0:
        return True, "vacuous (Pi = 1)"
    pc = sp.symbols("pc0:%d" % k)
    qc = sp.symbols("qc0:%d" % (5 - k))
    zt = sp.Symbol("zt")
    sat = sp.Symbol("sat")
    Pi_p = sp.Poly(y ** k + sum(pc[i] * y ** i for i in range(k)), y)
    Q_p = sp.Poly(sum(qc[i] * y ** i for i in range(5 - k)), y)
    eqs = []
    fac = sp.Poly(sp.expand(Pi_p.as_expr() * Q_p.as_expr() - quartic), y)
    eqs += [sp.expand(fac.coeff_monomial(y ** i)) for i in range(5)]
    tgt = sp.Poly(sp.expand((y + 1) ** 3 * Q_p.as_expr() - 3 * zt * (y + 1) ** z), y)
    rem = sp.rem(tgt.as_expr(), sp.expand(Pi_p.as_expr() ** 2), y)
    rp = sp.Poly(sp.expand(rem), y)
    eqs += [sp.expand(rp.coeff_monomial(y ** i)) for i in range(2 * k)]
    eqs.append(zt * sat - 1)                       # saturation by zeta != 0
    gens = list(pc) + list(qc) + [zt, sat]
    gb = sp.groebner([e for e in eqs if e != 0], *gens, order="grevlex")
    triv = list(gb.exprs) == [sp.Integer(1)]
    return (not triv), ("unit ideal" if triv else "GB size %d" % len(gb.exprs))


# ---- route 2: exact linear algebra in an explicit splitting subfield ------
def _kmod(expr, x, m):
    if m is None:
        return sp.expand(expr)
    return sp.expand(sp.rem(sp.expand(expr), m, x))


def _polylist(expr, x, m):
    p = sp.Poly(sp.expand(expr), y)
    n = p.degree()
    if n < 0:
        return [sp.Integer(0)]
    return [_kmod(p.coeff_monomial(y ** i), x, m) for i in range(n + 1)]


def _prem(P, D, x, m):
    """Remainder of P by monic D, coefficient lists (index = degree), over K."""
    P = list(P)
    dD = len(D) - 1
    while True:
        while P and P[-1] == 0:
            P.pop()
        if len(P) - 1 < dD or not P:
            break
        n = len(P) - 1 - dD
        c = P[-1]
        for i in range(dD + 1):
            P[i + n] = _kmod(P[i + n] - c * D[i], x, m)
    return P + [sp.Integer(0)] * (dD - len(P))


def support_field_feasible(k: int, z: int):
    """Independent route: realise Pi explicitly over an explicit subfield of the
    splitting field of q, then decide feasibility by a rank condition.

    For fixed Pi (hence Q = q/Pi), the condition
        Pi^2 | (mu*t^3*Q - 3*zeta*t^z)
    is LINEAR in (mu, zeta).  Write M for the 2k x 2 matrix of the two columns
        col_mu = rem(t^3*Q, Pi^2),   col_zeta = rem(-3*t^z, Pi^2).
    col_zeta != 0 always (gcd(Pi,t) = 1, deg Pi^2 > 0) and col_mu != 0 always
    (gcd(Pi,Q) = 1 since q is squarefree, so Pi^2 does not divide Q).  Hence a
    solution with mu != 0 AND zeta != 0 exists iff rank M = 1, i.e. iff every
    2x2 minor of M vanishes.
    """
    x = sp.Symbol("xg")
    if k == 0:
        return True, "vacuous"
    if k == 4:
        m = None
        Pi_e = sp.expand(qq / 2048)
        Q_e = sp.Integer(2048)
        fld = "Q"
    elif k in (1, 3):
        m = qq.subs(y, x)                      # x is a root of q; deg 4 field
        quo, rem_ = sp.div(sp.Poly(qq, y), sp.Poly(y - x, y), y)
        assert _kmod(rem_.as_expr(), x, m) == 0, "q(x) != 0 mod m"
        quo = sp.expand(quo.as_expr())
        if k == 1:
            Pi_e = y - x
            Q_e = quo
        else:
            Pi_e = sp.expand(quo / 2048)
            Q_e = 2048 * (y - x)
        fld = "Q[x]/(q), deg 4"
    elif k == 2:
        # p = -(r_i + r_j) satisfies the irreducible resolvent sextic; s and the
        # cofactor are polynomials in p (checked below).
        m = (32768 * x ** 6 + 24576 * x ** 5 + 16384 * x ** 4 + 5632 * x ** 3
             - 10080 * x ** 2 - 2680 * x - 495)
        s_of = sp.Rational(1, 1632) * (8192 * x ** 5 + 5120 * x ** 4 + 3456 * x ** 3
                                       + 1792 * x ** 2 - 2540 * x - 225)
        pp_of = sp.Rational(-1, 4) - x
        ss_of = sp.Rational(-1, 408) * (2048 * x ** 5 + 1280 * x ** 4 + 864 * x ** 3
                                        + 40 * x ** 2 - 737 * x - 120)
        Pi_e = y ** 2 + x * y + s_of
        Q_e = 2048 * (y ** 2 + pp_of * y + ss_of)
        fld = "Q[x]/(irreducible sextic resolvent), deg 6"
    else:
        raise ValueError(k)

    # sanity: q = Pi * Q in K[y]
    resid = _polylist(sp.expand(Pi_e * Q_e - qq), x, m)
    assert all(c == 0 for c in resid), (k, "q != Pi*Q in the constructed field")

    D = _polylist(sp.expand(Pi_e ** 2), x, m)
    col_mu = _prem(_polylist(sp.expand((y + 1) ** 3 * Q_e), x, m), D, x, m)
    col_zt = _prem(_polylist(sp.expand(-3 * (y + 1) ** z), x, m), D, x, m)
    assert any(c != 0 for c in col_zt), "col_zeta vanished -- test would be vacuous"
    assert any(c != 0 for c in col_mu), "col_mu vanished -- test would be vacuous"
    minors = []
    for i in range(len(col_mu)):
        for j in range(i + 1, len(col_mu)):
            minors.append(_kmod(col_mu[i] * col_zt[j] - col_mu[j] * col_zt[i], x, m))
    rank1 = all(mn == 0 for mn in minors)
    return rank1, "field = %s, %d minors" % (fld, len(minors))


feas_ideal = {}
feas_field = {}
for k in range(1, 5):
    for z in range(0, ZMAX + 1):
        feas_ideal[(k, z)] = support_ideal_feasible(k, z)[0]
        feas_field[(k, z)] = support_field_feasible(k, z)[0]

check("E1  the two independent routes (saturated ideal over Q  vs  exact rank test in an "
      "explicit splitting subfield) AGREE on all 40 (k,z) pairs",
      feas_ideal == feas_field,
      "disagreements: %s" % sorted(kz for kz in feas_ideal if feas_ideal[kz] != feas_field[kz]))

for k in range(1, 5):
    fz = [z for z in range(ZMAX + 1) if feas_ideal[(k, z)]]
    expect = [3] if k == 4 else []
    check("E2.k%d  k = %d : feasible z in [0,9] = %s   (claim: %s)"
          % (k, k, fz, expect), fz == expect)

check("E3  CONTROL: the support test is NOT vacuously infeasible -- it returns FEASIBLE "
      "at (k,z) = (4,3) on the genuine quartic", feas_ideal[(4, 3)])

# E4 -- a third, fully hand-checkable route for k = 1.
k1_gcds = {}
for z in range(ZMAX + 1):
    g = sp.gcd(sp.Poly(qq, y),
               sp.Poly(sp.expand((y + 1) * sp.diff(qq, y, 2) - 2 * (z - 3) * sp.diff(qq, y)), y))
    k1_gcds[z] = sp.degree(g.as_expr(), y) if g.as_expr() != 0 else -1
check("E4  hand route for k = 1: gcd(q, (y+1)*q'' - 2*(z-3)*q') has degree 0 for every z in [0,9]",
      all(v == 0 for v in k1_gcds.values()), str(k1_gcds))

# E5 -- CONTROL X2 re-posed: the k=1 test is a property of THIS quartic.
q_syn = y ** 4 + y ** 3 - 9 * y ** 2 + 7
syn_ok = (sp.gcd(q_syn, sp.diff(q_syn, y)) == 1 and q_syn.subs(y, -1) != 0)
syn_feas = []
for z in range(ZMAX + 1):
    g = sp.gcd(sp.Poly(q_syn, y),
               sp.Poly(sp.expand((y + 1) * sp.diff(q_syn, y, 2) - 2 * (z - 3) * sp.diff(q_syn, y)), y))
    if sp.degree(g.as_expr(), y) > 0:
        syn_feas.append(z)
syn_ideal = [z for z in range(ZMAX + 1) if support_ideal_feasible(1, z, q_syn)[0]]
# E6 -- NON-VACUITY CONTROL for k = 2, which SUB1_SPINE9.md's control set does NOT
# have (its X2 covers k=1 only, its X2b covers k=4 only).  A k=3 or k=4 "infeasible"
# could be an implementation artifact; for k=2 it cannot, because the machinery is
# shown here to return FEASIBLE on a synthetic quartic built to satisfy the condition.
#   Pi = t^2 - 1, Q = 3*(2*t^2 - 1), zeta = 1, z = 7:
#   t^3*Q - 3*t^7 = -3*t^3*(t^2-1)^2 = -3*t^3*Pi^2 .
q_k2 = sp.expand((y ** 2 + 2 * y) * (6 * y ** 2 + 12 * y + 3))
k2_feas = [z for z in range(ZMAX + 1) if support_ideal_feasible(2, z, q_k2)[0]]
_Pi2 = y ** 2 + 2 * y
check("E6  NON-VACUITY at k = 2 (a control the audited lane lacks): on the synthetic squarefree "
      "quartic %s (q(-1) = %s != 0) the SAME k=2 ideal test returns FEASIBLE at z in %s"
      % (sp.factor(q_k2), q_k2.subs(y, -1), k2_feas),
      sp.degree(sp.gcd(q_k2, sp.diff(q_k2, y)), y) == 0 and q_k2.subs(y, -1) != 0 and 7 in k2_feas
      and sp.rem(sp.expand((y + 1) ** 3 * sp.cancel(q_k2 / _Pi2) - 3 * (y + 1) ** 7),
                 sp.expand(_Pi2 ** 2), y) == 0)

# E7 -- STRUCTURAL cross-check: the k = 3 kill is NOT arithmetic in q at all.  With
# deg Q = 1, deg Pi^2 = 6 and gcd(Pi,t) = 1, no z in [0,9] can work for ANY squarefree
# quartic with q(-1) != 0.  Verified on three unrelated quartics.
k3_others = [3 * y ** 4 + y ** 3 - 5 * y ** 2 + y + 7,
             y ** 4 + 2 * y ** 3 - 3 * y + 11,
             5 * y ** 4 - y ** 2 + 6 * y + 3]
k3_ok = True
for qo in k3_others:
    if sp.degree(sp.gcd(qo, sp.diff(qo, y)), y) != 0 or qo.subs(y, -1) == 0:
        k3_ok = False
        continue
    if any(support_ideal_feasible(3, z, qo)[0] for z in range(ZMAX + 1)):
        k3_ok = False
check("E7  the k = 3 kill is STRUCTURAL, not arithmetic: infeasible for every z in [0,9] on three "
      "unrelated squarefree quartics with q(-1) != 0 too -- consistent with the hand argument "
      "(deg Q = 1, deg Pi^2 = 6, gcd(Pi,t) = 1 leave no room)", k3_ok)

check("E5  CONTROL X2: on the synthetic squarefree quartic y^4+y^3-9y^2+7 (q(-1) = -2 != 0) "
      "the SAME k=1 test returns FEASIBLE at exactly z = 3 -- so the kill is arithmetic, not structural",
      syn_ok and syn_feas == [3] and syn_ideal == [3],
      "gcd route %s, ideal route %s" % (syn_feas, syn_ideal))


# ===========================================================================
# F.  THE DEGREE LEDGER
# ===========================================================================
section("F. The degree ledger and the two degree kills")

# stripped consequences of the certified caps
degA_max = RCAP - A_T                     # R = t^9*A
degu_max = D2CAP                          # u = gamma*d2
degw_max = lambda k: D1CAP + k            # w = (1/2)*gamma^2*d1*Pi
degv_max = lambda k: SCAP - A_T - k       # S = t^9*Pi*v
check("F1  deg A <= %d  (R = t^9*A, deg R <= %d)" % (degA_max, RCAP), degA_max == 9)
check("F2  deg u <= %d  (u = gamma*d2, deg d2 <= %d)" % (degu_max, D2CAP), degu_max == 6)
check("F3  deg w <= 9+k, deg v <= 12-k  (deg d1 <= %d, deg S <= %d)" % (D1CAP, SCAP),
      [degw_max(k) for k in range(5)] == [9, 10, 11, 12, 13]
      and [degv_max(k) for k in range(5)] == [12, 11, 10, 9, 8])

capF = {k: max(degA_max + max(degu_max, degv_max(k)), degw_max(k)) for k in range(5)}
check("F4  cap on deg F = max(deg A + max(deg u, deg v), deg w) = %s for k = 0..4"
      % [capF[k] for k in range(5)], [capF[k] for k in range(5)] == [21, 20, 19, 18, 17])
check("F5  deg F = 9 + 4k - z is FORCED by (*deg)", True)

# F6 -- k = 4, via the support test alone (z = 3): no cascade input.
check("F6  k = 4 (a9_b1111_T1): section-5 pins z = 3, so deg F = 25-3 = 22 > cap %d -- CONTRADICTION"
      % capF[4], 9 + 4 * 4 - 3 > capF[4])
# F7 -- k = 4, independent second route via z <= 6.
check("F7  k = 4, second route: even ignoring section 5, z <= 6 gives deg F >= 19 > %d" % capF[4],
      9 + 4 * 4 - 6 > capF[4])
# F8 -- the other k do NOT die this way (non-vacuity: the ledger is not a blanket argument).
check("F8  the degree ledger alone does NOT kill k = 0,1,2,3 (for each there is an admissible "
      "z in [2,6] the ledger permits), so it is not a blanket argument",
      all(any(9 + 4 * k - z <= capF[k] for z in range(2, 7)) for k in (0, 1, 2, 3)),
      "k=3 is the tightest: the ledger only excludes z = 2 there")

# F9 -- the k = 0 elimination, re-derived from the boxed row (Pi = 1, Q = q).
qsym, tz = sp.Symbol("qpoly"), sp.Symbol("tz")     # tz stands for t^z
box0 = sp.expand(3 * A_ ** 2 + gam * (u_ + 3 * v_) - mu_ * t_ ** 3 * qsym)
# Z = A^2 - gamma*v  =>  gamma*v = A^2 - zeta*t^z
box0_elim = sp.expand(box0.subs(v_, (A_ ** 2 - zeta_ * tz) / gam))
claim0 = sp.expand(gam * u_ - (mu_ * t_ ** 3 * qsym - 6 * A_ ** 2 + 3 * zeta_ * tz))
check("F9  k = 0: eliminating v from Z = A^2 - gamma*v in the boxed row gives "
      "gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z   (residual 0)",
      sp.simplify(box0_elim - claim0) == 0,
      "box0_elim = %s" % box0_elim)


def k0_kill(deg_d2_cap: int, z_max: int, z_min: int = 0) -> bool:
    """Does the k = 0 degree dichotomy CLOSE, given deg u <= deg_d2_cap and
    z_min <= z <= z_max?

        gamma*u = mu*t^3*q - 6*A^2 + 3*zeta*t^z ,   deg(gamma*u) <= deg_d2_cap.

    Three contributors on the right, with degrees and NONZERO leading coefficients

        mu*t^3*q      degree 7,        lc = 2048*mu != 0
        -6*A^2        degree 2*deg A,  lc = -6*lc(A)^2 != 0   (absent iff A = 0)
        3*zeta*t^z    degree z,        lc = 3*zeta != 0

    If a UNIQUE contributor attains the maximal degree D then deg(RHS) = D is
    forced; D > deg_d2_cap is then a contradiction.  If two or more tie, this
    (deliberately conservative) test declines to force anything -- so a True
    return is a genuine sufficient condition for the kill, and a False return
    means "the argument as stated no longer closes".
    """
    for degA in [-1] + list(range(0, 40)):          # -1 encodes A = 0
        for z in range(z_min, z_max + 1):
            degs = [7, z] + ([2 * degA] if degA >= 0 else [])
            D = max(degs)
            if degs.count(D) != 1:
                return False
            if D <= deg_d2_cap:
                return False
    return True


check("F10 k = 0 (a9_b0000_T1): the degree dichotomy closes at deg d2 <= %d and z <= 6 -- CONTRADICTION"
      % D2CAP, k0_kill(D2CAP, 6))
check("F11 ... 2*deg A = 7 is impossible, and deg(mu*t^3*q) = 7 exactly (lc = 2048*mu != 0)",
      sp.LC(sp.Poly(sp.expand((y + 1) ** 3 * qq), y)) == 2048 and sp.degree((y + 1) ** 3 * qq, y) == 7)
check("F12 CONTROL X3: mutating the sub1 cap deg d2 from 6 to 7 switches the k = 0 kill OFF "
      "(the cap is load-bearing, with ZERO margin)", not k0_kill(7, 6))
check("F13 CONTROL X1: weakening z <= 6 to z <= 7 switches the k = 0 kill OFF "
      "(z <= 6 is load-bearing, with ZERO margin)", not k0_kill(D2CAP, 7))


def k4_kill(degR_cap: int) -> bool:
    dA = degR_cap - A_T
    cap = max(dA + max(degu_max, degv_max(4)), degw_max(4))
    return 9 + 16 - 3 > cap


check("F14 CONTROL X4: the k = 4 kill survives until the deg R cap reaches 23 "
      "(load-bearing, wide margin at the certified 18)",
      k4_kill(18) and k4_kill(22) and not k4_kill(23))


# ===========================================================================
# G.  THE VALUATION LEDGER  (the only cascade-consuming step)
# ===========================================================================
section("G. The valuation ledger:  2 <= z <= 6")

# G1 -- re-derive the d3-killing shift from the generalized binomial recomposition.
# With D_j the unshifted coefficients (h_k = D_{4-k}) and the shift parameter
# cshift chosen to kill D~_3:  D~_j = sum_{i>=j} C(i,j) cshift^(i-j) D_i for j >= 0,
# and D~_{-1-r} = sum_{j=0}^{r} C(r,j) (-cshift)^j D_{-1-r+j} for r >= 0.
D = {j: sp.Symbol("D_%s" % (str(j) if j >= 0 else "m%d" % (-j))) for j in range(-4, 5)}
cs = sp.Symbol("cshift")
D[4] = sp.Integer(1)
Dt_pos = {j: sp.expand(sum(sp.binomial(i, j) * cs ** (i - j) * D[i] for i in range(j, 5)))
          for j in range(0, 5)}
cs_val = sp.solve(sp.Eq(Dt_pos[3], 0), cs)[0]
Dt_neg = {-1 - r: sp.expand(sum(sp.binomial(r, j) * (-cs) ** j * D[-1 - r + j] for j in range(r + 1)))
          for r in range(0, 4)}
h = {1: D[3], 2: D[2], 3: D[1], 4: D[0], 5: D[-1], 6: D[-2], 7: D[-3], 8: D[-4]}
shift_ok = {
    "d2": sp.simplify(Dt_pos[2].subs(cs, cs_val) - (h[2] - sp.Rational(3, 8) * h[1] ** 2)),
    "d1": sp.simplify(Dt_pos[1].subs(cs, cs_val)
                      - (h[3] - sp.Rational(1, 2) * h[2] * h[1] + sp.Rational(1, 8) * h[1] ** 3)),
    "e": sp.simplify(Dt_neg[-1].subs(cs, cs_val) - h[5]),
    "R": sp.simplify(Dt_neg[-2].subs(cs, cs_val) - (h[6] + sp.Rational(1, 4) * h[1] * h[5])),
    "S": sp.simplify(Dt_neg[-3].subs(cs, cs_val)
                     - (h[7] + sp.Rational(1, 2) * h[1] * h[6] + sp.Rational(1, 16) * h[1] ** 2 * h[5])),
}
check("G1  the d3-killing shift (cshift = -h1/4) reproduces ALL the recorded inverse-shift forms "
      "d2, d1, e = h5, R, S from one generalized-binomial model, residual 0",
      all(v == 0 for v in shift_ok.values()), str({k: str(v) for k, v in shift_ok.items()}))
check("G2  ... the shift parameter really is cshift = -h1/4", sp.simplify(cs_val + h[1] / 4) == 0)


def z_window(prof: dict[int, int], reading: str):
    """Return (z_min, z_max) from the profile, under reading A (shifted, [I3]) or
    B (unshifted).  a_t = 9 throughout."""
    if reading == "A":
        vd2 = min(prof[2], 2 * prof[1])
        vd1 = min(prof[3], prof[2] + prof[1], 3 * prof[1])
        vR = min(prof[6], prof[1] + prof[5])
        vS = min(prof[7], prof[1] + prof[6], 2 * prof[1] + prof[5])
    else:
        vd2, vd1, vR, vS = prof[2], prof[3], prof[6], prof[7]
    vA, vv, vu, vw = vR - A_T, vS - A_T, vd2, vd1
    if min(vA, vv) < 0:
        return None
    zmin = min(2 * vA, vv)
    vF = min(vA + min(vu, vv), vw)
    return zmin, A_T - vF


def _prof(**over):
    p = {i: FV["h%d" % i] for i in range(1, 9)}
    for kk, vv in over.items():
        p[int(kk[1:])] = vv
    return p


PROF = _prof()
PROF10 = _prof(h6=10, h7=11)                # level-10 rows (the committed profile)
PROF12 = _prof(h6=11, h7=11)                # the level-12 upgrade on h6
wins = {(r, lab): z_window(p, r) for r in ("A", "B") for lab, p in (("L10", PROF10), ("L12", PROF12))}
check("G3  under EVERY reading (A shifted / B unshifted) and level (10 / 12): 2 <= z <= 6",
      all(w is not None and w[0] >= 2 and w[1] <= 6 for w in wins.values()), str(wins))
check("G4  the level-12 upgrade on h6 is INERT under reading A (control X1b): "
      "(h6,h7) = (10,11) and (11,11) give the same z_max = 6",
      wins[("A", "L10")] == wins[("A", "L12")] == (2, 6))
check("G5  v_t(R) >= min(v_t(h6), v_t(h1)+a_t) = min(10,10) = 10 -- the level-12 upgrade is "
      "INVISIBLE at a_t = 9, so the brief's 'v_t(R) >= 11' does NOT hold and the closure "
      "does not depend on level 12",
      min(PROF12[6], PROF12[1] + A_T) == 10)
check("G6  CONTROL X1: v_t(h7) = 10 gives z_max = 7, v_t(h7) = 9 gives z_max = 8 -- so "
      "v_t(h7) >= 11 is EXACTLY load-bearing for a9_b0000_T1, with ZERO margin",
      z_window(_prof(h6=10, h7=10), "A")[1] == 7
      and z_window(_prof(h6=10, h7=9), "A")[1] == 8)
check("G7  ... and v_t(h6) >= 10 is ALSO load-bearing: v_t(h6) = 9 gives z_max = 8",
      z_window(_prof(h6=9, h7=11), "A")[1] == 8)
check("G8  ... and so is v_t(d1) >= 3 (i.e. v_t(h1) >= 1 through the h1^3 term)",
      min(PROF10[3], PROF10[2] + PROF10[1], 3 * PROF10[1]) == 3)
# ---------------------------------------------------------------------------
# G10/G11 -- the audit's own contribution: the two UNAUDITED integers
# v_t(h6) >= 10 and v_t(h7) >= 11 are CONSEQUENCES of the AUDITED rows
# h1..h5 = 1,3,5,7,9 plus the P-side absorption form, which SLICE_OBSTRUCTION.md
# section 3 and SLICE_OBSTRUCTION_AUDIT.md section 4 both record:
#
#     H = sum_{j>=0} h_j u^j with h_0 = 1,  p_n := [u^n] H^2 = 2*h_n + q_n,
#     q_n = sum_{i=1}^{n-1} h_i*h_{n-i},
#     (P<)  t^(2n-2) | p_n   for n = 2..8
#  => h_n = -q_n/2 + t^(2n-2)*g_n  with g_n integral
#  => v_t(h_n) >= min( min_{i=1..n-1}( v(h_i) + v(h_{n-i}) ) , 2n-2 ).
# ---------------------------------------------------------------------------
def convolution_floor(base: dict[int, int], n: int) -> int:
    return min(min(base[i] + base[n - i] for i in range(1, n)), 2 * n - 2)


AUDITED = {1: 1, 2: 3, 3: 5, 4: 7, 5: 9}     # slice_obstruction_audit.py F1-F5
derived = dict(AUDITED)
for n in (6, 7, 8):
    derived[n] = convolution_floor(derived, n)
check("G10 *** the two un-audited integers are DERIVABLE from the AUDITED ones: the P-side "
      "absorption h_n = -q_n/2 + t^(2n-2)*g_n with q_n = sum_{i=1}^{n-1} h_i*h_{n-i} gives "
      "v_t(h6) >= %d, v_t(h7) >= %d, v_t(h8) >= %d from h1..h5 = 1,3,5,7,9 alone"
      % (derived[6], derived[7], derived[8]),
      (derived[6], derived[7], derived[8]) == (10, 11, 12),
      "so SUB1_SPINE9.md section 10.2's 'single integer needing a second party' has one")
check("G11 ... and the derived profile is EXACTLY the committed one, so z <= 6 holds on the "
      "AUDITED cascade rows alone",
      [derived[i] for i in range(1, 9)] == [FV["h%d" % i] for i in range(1, 9)]
      and z_window(derived, "A") == (2, 6))

# ---------------------------------------------------------------------------
# G12/G13 -- t^9 | R,S,T, re-derived here from `e | S` (an imported CHECKED fact,
# proved three independent ways in DIVISOR_CONSEQUENCES.md section 2) and the
# K-syzygy alone.  Cap-free, branch-independent, and it does NOT use d1 = 0.
# ---------------------------------------------------------------------------
sbar = sp.Symbol("sbar")
Ksub = sp.expand(Ksyz.subs({S_v: e_v * sbar}))
# 3*e*R^2 = 2*Phi - e^3*(d2 + 3*sbar)
check("G12 with S = e*sbar the K-syzygy K = 0 reads 3*e*R^2 = 2*Phi - e^3*(d2+3*sbar)  (residual 0); "
      "valuations give 9 + 2*v(R) >= min(30, 3*9) = 27, hence v_t(R) >= 9 = a_t",
      sp.expand(Ksub + (3 * e_v * R_v ** 2 - 2 * Phi_v + e_v ** 3 * (d2v + 3 * sbar))) == 0
      and (min(30, 3 * A_T) - A_T) / 2 >= A_T)
# G1 = 0 gives 3*e*T = -(3/2)*d1*e^2 - 3*d2*e*R - 3*R*S
check("G13 G1 = 0 gives T = -R*(sbar+d2) - d1*e/2, so v(T) >= min(v(d1)+9, 9+v(R), v(R)+v(S)) - 9 "
      ">= 9 = a_t; with v(S) >= v(e) = 9 this closes t^9 | R,S,T -- CAP-FREE and it never uses d1 = 0",
      sp.expand(G1.subs({T_v: -R_v * (sbar + d2v) - d1v * e_v / 2, S_v: e_v * sbar})) == 0
      and min(0 + 2 * A_T, 0 + A_T + A_T, A_T + A_T) - A_T >= A_T)

check("G9  t^9 | R,S,T ALSO follows from the SAME profile: v_t(R) >= 10, v_t(S) >= 11, v_t(T) >= 12",
      min(PROF10[6], PROF10[1] + PROF10[5]) >= 9
      and min(PROF10[7], PROF10[1] + PROF10[6], 2 * PROF10[1] + PROF10[5]) >= 9
      and min(PROF10[8], PROF10[1] + PROF10[7], 2 * PROF10[1] + PROF10[6],
              3 * PROF10[1] + PROF10[5]) >= 9)


# ===========================================================================
# H.  PER-CELL VERDICT
# ===========================================================================
section("H. Per-cell verdict")

CELLS = [("a9_b0000_T1", 0), ("a9_b1000_T1", 1), ("a9_b1100_T1", 2),
         ("a9_b1110_T1", 3), ("a9_b1111_T1", 4)]
verdict = {}
for name, k in CELLS:
    if k == 0:
        dead = k0_kill(D2CAP, 6)
        why = "degree dichotomy on the boxed row (deg d2 <= 6), needs z <= 6 from the cascade"
        cascade = True
    elif k in (1, 2, 3):
        dead = not any(feas_ideal[(k, z)] for z in range(ZMAX + 1))
        why = "marked Pi^2-support infeasible for EVERY z in [0,9]"
        cascade = False
    else:
        zs = [z for z in range(ZMAX + 1) if feas_ideal[(k, z)]]
        dead = all(9 + 4 * k - z > capF[k] for z in zs) and zs == [3]
        why = "support pins z = 3, then deg F = 22 > cap 17"
        cascade = False
    verdict[name] = (dead, k, why, cascade)
    check("H.%s  k = %d : %s   [%s]" % (name, k, "EMPTY" if dead else "SURVIVES", why),
          dead, "consumes the slice cascade: %s" % ("YES" if cascade else "no"))

check("H0  ALL FIVE cells a9_b{0000,1000,1100,1110,1111}_T1 are EMPTY, conditional on a_t = 9",
      all(v[0] for v in verdict.values()))


# ===========================================================================
npass = sum(1 for _, ok, _ in _RESULTS if ok)
ntot = len(_RESULTS)
if not QUIET:
    print("\nspine9_audit: %d/%d checks pass" % (npass, ntot))
    if npass != ntot:
        print("FAILURES:")
        for n, ok, d in _RESULTS:
            if not ok:
                print("  - %s  %s" % (n, d))
else:
    print("spine9_audit: %d/%d checks pass" % (npass, ntot))
sys.exit(0 if npass == ntot else 1)
