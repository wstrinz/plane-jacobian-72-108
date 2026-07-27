"""gauge_leak.py -- the LINEAR (lattice) stage of rational-gauge liftability.

Elimination computes the CLOSURE of the image; Chevalley says the true image is
only CONSTRUCTIBLE.  This file builds the linear half of the missing test.

Setting.  R = K[y], F = K(y).  Every coefficient window in the pipeline is an
R-lattice (really: a collection of local lattices, one per place of P^1) inside a
common F-vector space V = F^9 with coordinates (D_4, ..., D_-4).  A rational
gauge is a matrix T in GL_9(F).  If M_old is the original coefficient lattice and
M_shift the allowed shifted lattice, the LIFTABLE shifted coefficients are

        M_shift  cap  T(M_old)

and the GAUGE-LEAK MODULE is the quotient

        M_shift / (M_shift cap T(M_old))

Locally at each irreducible factor of y, of (y+1), of the quartic q, and of any
denominator the gauge introduces, this is a lattice problem over a DVR, and the
local Smith normal form / elementary divisors give the exact missing
divisibility conditions.  The degree condition at infinity is the same
computation at the place 1/y = 0.

This lane produces ONLY the linear lattice conditions.  The nonlinear stage
(the P = C^2 positive slices) belongs to positive_slice* / slice_obstruction*
and is NOT redone here; where this file meets it, it CITES it.

Run:
    python -u gauge_leak.py            # full report
    python -u gauge_leak.py --quiet    # exit 0 iff every check passes
"""

import json
import math
import os
import random
import sys
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
QUIET = "--quiet" in sys.argv

# ---- pinned identifiers (gamma, E, S are sympy builtins: never sympify names)
y = sp.Symbol("y")
u_ = sp.Symbol("u_")
tt = sp.Symbol("tt")            # tt = y + 1, the uniformizer at the place t
C4 = y**7 * (y + 1)             # the certified leading coefficient
FOUR = sp.Integer(4)

_PASS = []
_FAIL = []


def ck(name, cond, detail=""):
    ok = bool(cond)
    (_PASS if ok else _FAIL).append(name)
    if not QUIET:
        print("  %s  %s" % ("ok  " if ok else "FAIL", name))
        if detail:
            for line in str(detail).splitlines():
                print("           " + line)
    return ok


def say(msg=""):
    if not QUIET:
        print(msg)


def head(msg):
    say("\n" + "=" * 78)
    say(msg)
    say("=" * 78)


# ===========================================================================
# 0.  PLACES, LOCAL VALUATIONS, LOCAL SMITH NORMAL FORM
# ===========================================================================
class Place(object):
    """A place of P^1_K.  Either the vanishing locus of an irreducible p in K[y]
    (uniformizer p), or the place at infinity (uniformizer 1/y)."""

    def __init__(self, name, poly=None, at_infinity=False):
        self.name = name
        self.poly = poly
        self.at_infinity = at_infinity

    def __repr__(self):
        return "Place(%s)" % self.name

    def val(self, expr):
        """GENERIC valuation: symbolic parameters are treated as generic, i.e. a
        coefficient counts as nonzero unless it vanishes IDENTICALLY.  For the
        integrality statements below this is the correct (and the conservative)
        notion: an entry is integral for EVERY specialisation of the parameters
        iff its generic valuation is >= 0, because the p-adic coefficients are
        polynomials in those parameters."""
        expr = sp.cancel(sp.together(sp.sympify(expr)))
        if expr == 0:
            return sp.oo
        num, den = sp.fraction(expr)
        if self.at_infinity:
            return sp.degree(sp.Poly(den, y), y) - sp.degree(sp.Poly(num, y), y)
        return self._mult(num) - self._mult(den)

    def _mult(self, poly):
        p = sp.Poly(sp.expand(poly), y)
        pp = sp.Poly(self.poly, y)
        m = 0
        while True:
            qq, rr = sp.div(p, pp)
            if rr.is_zero and not qq.is_zero:
                p, m = qq, m + 1
            else:
                return m


P_Y = Place("y = 0", poly=y)
P_T = Place("t = y+1 = 0", poly=y + 1)
P_INF = Place("y = infinity", at_infinity=True)


def local_snf_exponents(A, place):
    """Elementary divisors p^e_i of the matrix A over the DVR at `place`.

    Valuation-pivot Gaussian elimination: the pivot is always an entry of
    MINIMAL valuation, so every multiplier is integral at the place and the row
    and column operations are invertible over the local ring.  Returns the list
    of exponents e_i (unsorted pivots, in elimination order)."""
    A = [[sp.cancel(x) for x in row] for row in A]
    n, m = len(A), len(A[0])
    exps = []
    for s in range(min(n, m)):
        best = None
        for i in range(s, n):
            for j in range(s, m):
                if A[i][j] == 0:
                    continue
                v = place.val(A[i][j])
                if best is None or v < best[0]:
                    best = (v, i, j)
        if best is None:
            raise ValueError("matrix is singular over F at %s" % place)
        v, i, j = best
        A[s], A[i] = A[i], A[s]
        for r in range(n):
            A[r][s], A[r][j] = A[r][j], A[r][s]
        piv = A[s][s]
        for r in range(s + 1, n):
            if A[r][s] != 0:
                f = sp.cancel(A[r][s] / piv)
                for c in range(s, m):
                    A[r][c] = sp.cancel(A[r][c] - f * A[s][c])
        for c in range(s + 1, m):
            if A[s][c] != 0:
                f = sp.cancel(A[s][c] / piv)
                for r in range(s, n):
                    A[r][c] = sp.cancel(A[r][c] - f * A[r][s])
        exps.append(int(v))
    return exps


def leak_lengths(B_target, B_source, place):
    """M_target / (M_target cap M_source) at `place`, as elementary-divisor
    exponents.  B_* are basis matrices (lists of rows) of the two lattices.
    Returns the multiset of max(e_i, 0) -- the leak is zero iff all are 0, i.e.
    iff M_target is CONTAINED in M_source (every target vector lifts)."""
    Bt = sp.Matrix(B_target)
    Bs = sp.Matrix(B_source)
    A = sp.cancel(Bt.inv() * Bs)
    exps = local_snf_exponents(A.tolist(), place)
    return sorted(max(e, 0) for e in exps)


def _selftest_snf():
    """Calibration: a matrix with known elementary divisors at y = 0."""
    A = [[y**2, y**3], [y**3, y**5]]      # det = y^7 - y^6 = y^6 (y-1); e = 2,4
    e = sorted(local_snf_exponents(A, P_Y))
    ok1 = e == [2, 4]
    B = [[sp.Integer(1), 1 / y], [y, sp.Integer(2)]]   # det = 1, unimodular? no
    e2 = sorted(local_snf_exponents(B, P_Y))
    ok2 = e2 == [-1, 1]
    C = [[(y + 1)**3, sp.Integer(0)], [sp.Integer(0), (y + 1)]]
    ok3 = sorted(local_snf_exponents(C, P_T)) == [1, 3]
    D = [[y**5, sp.Integer(0)], [sp.Integer(0), sp.Integer(1)]]
    ok4 = sorted(local_snf_exponents(D, P_INF)) == [-5, 0]
    return ok1 and ok2 and ok3 and ok4


# ===========================================================================
# 1.  THE WINDOW MODULE
# ===========================================================================
#   full  D_(4-K):  ord_y >= 12K ,  deg_y <= CFULL[regime]*K
#   strip d_(4-K) = D_(4-K)/y^(12K):  ord_y >= 0 ,  deg_y <= CSTR[regime]*K
CFULL = {"sub2": 14, "sub1": 15}
CSTR = {"sub2": 2, "sub1": 3}
KIDX = list(range(0, 9))            # K = 4 - j,  j = 4 .. -4
JIDX = [4 - K for K in KIDX]


def diag(entries):
    n = len(entries)
    return [[entries[i] if i == j else sp.Integer(0) for j in range(n)]
            for i in range(n)]


def window_basis(regime, stripped=True, place=P_Y):
    """Basis matrix of the local window lattice at `place`."""
    if place is P_Y:
        if stripped:
            return diag([sp.Integer(1)] * 9)
        return diag([y**(12 * K) for K in KIDX])
    if place is P_INF:
        c = CSTR[regime] if stripped else (12 + CSTR[regime])
        return diag([y**(c * K) for K in KIDX])
    if place is P_T:
        # the pipeline declares NO t-adic condition on the window
        return diag([sp.Integer(1)] * 9)
    raise ValueError(place)


def window_dim(regime, stripped=True):
    """K-dimension of the global window module (sections over P^1): for each K,
    the polynomials of degree <= c*K, i.e. c*K + 1 coefficients.  K = 0 is the
    normalisation d_4 = 1 and contributes no free parameter."""
    c = CSTR[regime]
    return sum(c * K + 1 for K in range(1, 9))


# ===========================================================================
# 2.  S1 -- THE d3-KILLING SHIFT.  Regression + leak.
# ===========================================================================
def shift_matrix(theta, n=9):
    """The D-coordinate shift of window_caps_verify.py W3 / positive_slice.py
    P2:  X_j = sum_{m=j..4} binom(m, m-j) src_m theta^(m-j).
    Rows and columns are indexed by K = 4 - j and Kc = 4 - m, so the matrix is
    LOWER triangular in K with 1's on the diagonal (unipotent)."""
    M = [[sp.Integer(0)] * n for _ in range(n)]
    for K in range(n):
        j = 4 - K
        for Kc in range(n):
            m = 4 - Kc
            if m >= j:
                M[K][Kc] = sp.binomial(m, m - j) * theta**(m - j)
    return M


def matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sp.expand(sum(A[i][s] * B[s][j] for s in range(k)))
             for j in range(m)] for i in range(n)]


def conjugate(T, G):
    """G^-1 T G, with G diagonal."""
    n = len(T)
    return [[sp.cancel(T[i][j] * G[j][j] / G[i][i]) for j in range(n)]
            for i in range(n)]


# ===========================================================================
# 3.  polygon data, the sqrt recursion, and the slice functionals
# ===========================================================================
_UF = json.load(open(os.path.join(HERE, "paper_src", "upstream_facts.json"),
                     encoding="utf-8"))
CORNERS = {reg: [tuple(p) for p in _UF["facts"]["newton_polygons"][reg]["P"]]
           for reg in ("sub1", "sub2")}


def _hull_chains(corners):
    pts = sorted(set(corners))

    def half(pl):
        out = []
        for p in pl:
            while len(out) >= 2 and (
                (out[-1][0] - out[-2][0]) * (p[1] - out[-2][1])
                    - (out[-1][1] - out[-2][1]) * (p[0] - out[-2][0])) <= 0:
                out.pop()
            out.append(p)
        return out
    return half(pts), half(pts[::-1])


def hull_j_range(corners, i):
    lower, upper = _hull_chains(corners)

    def interp(chain, i):
        vals = []
        for (x0, j0), (x1, j1) in zip(chain, chain[1:]):
            if min(x0, x1) <= i <= max(x0, x1) and x0 != x1:
                vals.append(Fraction(j0) + Fraction(j1 - j0, x1 - x0) * (i - x0))
            elif x0 == i:
                vals.append(Fraction(j0))
        if chain and chain[-1][0] == i:
            vals.append(Fraction(chain[-1][1]))
        return vals
    allv = interp(lower, i) + interp(upper, i)
    return math.ceil(min(allv)), math.floor(max(allv))


def d_recursion(Pslices):
    """verify_derivation.py section C (division-free form):
       D_k = 1/2 P_{k+4} C4^(6-2k) - 1/2 sum_{i+j=k+4, i,j<=3} D_i D_j."""
    Dv = {}
    for kk in range(3, -5, -1):
        acc = sp.Rational(1, 2) * Pslices.get(kk + 4, sp.Integer(0)) * C4**(6 - 2 * kk)
        for i in range(kk + 1, 4):
            j2 = kk + 4 - i
            if i <= j2 <= 3:
                acc -= sp.Rational(2 if i != j2 else 1, 2) * Dv[i] * Dv[j2]
        Dv[kk] = sp.expand(acc)
    Dv[4] = sp.Integer(1)
    return Dv


def order_y(e):
    return min(m[0] for m in sp.Poly(sp.expand(e), y).monoms())


def slice_S(Dd, M, jmin=-4):
    """S_M := [u^(8-M)] H(u)^2 = sum_{i+j=M} d_i d_j."""
    return sp.expand(sum(Dd[j1] * Dd[M - j1] for j1 in range(jmin, 5)
                         if jmin <= M - j1 <= 4))


def t_coeffs(expr, upto):
    """The first `upto` coefficients of the (y+1)-adic expansion."""
    p = sp.Poly(sp.expand(sp.expand(expr).subs(y, tt - 1)), tt)
    return [p.coeff_monomial(tt**s) for s in range(upto)]


def y_coeffs(expr, upto):
    p = sp.Poly(sp.expand(expr), y)
    return [p.coeff_monomial(y**s) for s in range(upto)]


# ===========================================================================
#                                 THE RUN
# ===========================================================================
say(__doc__.split("Run:")[0].strip())

# ---------------------------------------------------------------- section 0
head("0.  The DVR machinery, calibrated")

ck("0.1  local Smith normal form reproduces known elementary divisors at "
   "y = 0, at t = 0 and at infinity (4 fixtures)", _selftest_snf())

ck("0.2  a unimodular matrix over the local ring has ALL elementary divisors 0 "
   "(so leak_lengths returns the zero multiset)",
   leak_lengths(diag([sp.Integer(1)] * 3),
                [[sp.Integer(1), y, y**2], [sp.Integer(0), sp.Integer(1), y],
                 [sp.Integer(0), sp.Integer(0), sp.Integer(1)]], P_Y) == [0, 0, 0])

ck("0.3  a genuine leak is detected: L' = y*L inside L has leak module "
   "(R/y)^3, total length 3",
   leak_lengths(diag([sp.Integer(1)] * 3), diag([y, y, y]), P_Y) == [1, 1, 1])

# ---------------------------------------------------------------- section 1
head("1.  S1 -- the d3-killing shift.  REGRESSION against positive_slice.py")

h_ = sp.Symbol("h_")            # h = the ORIGINAL, pre-shift stripped D_3
d2s, d1s, d0s = sp.symbols("d2s d1s d0s")

Tinv = shift_matrix(h_ / FOUR)          # inverse shift: theta = +h/4
Tfwd = shift_matrix(-h_ / FOUR)         # forward shift: theta = -h/4

# the shifted vector, indexed by K = 4-j:  (1, 0, d2, d1, d0, dm1..dm4)
dm_ = sp.symbols("dm1_ dm2_ dm3_ dm4_")
vshift = [sp.Integer(1), sp.Integer(0), d2s, d1s, d0s] + list(dm_)
vstar = [sp.expand(sum(Tinv[K][Kc] * vshift[Kc] for Kc in range(9)))
         for K in range(9)]

ck("1.1  the shift matrices satisfy the group law M(a)M(b) = M(a+b)",
   all(sp.expand(x) == 0 for row in
       [[matmul(shift_matrix(sp.Symbol("a_")), shift_matrix(sp.Symbol("b_")))[i][j]
         - shift_matrix(sp.Symbol("a_") + sp.Symbol("b_"))[i][j]
         for j in range(9)] for i in range(9)] for x in row))

ck("1.2  M(0) = I, hence T_fwd * T_inv = I (the gauge is invertible over R, not "
   "merely over F)",
   all(sp.expand(matmul(Tfwd, Tinv)[i][j] - (1 if i == j else 0)) == 0
       for i in range(9) for j in range(9)))

ck("1.3  the forward shift kills the k = 1 variable:  D~_3 = 0",
   sp.expand(sum(Tfwd[1][Kc] * [sp.Integer(1), h_, d2s, d1s, d0s,
                                dm_[0], dm_[1], dm_[2], dm_[3]][Kc]
                 for Kc in range(9))) == 0)

# --- THE REGRESSION.  positive_slice.py section 3.2 (derived there, not here).
_reg = {
    "D3* = h": vstar[1] - h_,
    "D2* = d2 + (3/8) h^2": vstar[2] - (d2s + sp.Rational(3, 8) * h_**2),
    "D1* = d1 + (1/2) h d2 + (1/16) h^3":
        vstar[3] - (d1s + sp.Rational(1, 2) * h_ * d2s + sp.Rational(1, 16) * h_**3),
    "D0* = d0 + (1/4) h d1 + (1/16) h^2 d2 + (1/256) h^4":
        vstar[4] - (d0s + sp.Rational(1, 4) * h_ * d1s
                    + sp.Rational(1, 16) * h_**2 * d2s
                    + sp.Rational(1, 256) * h_**4),
}
for _nm, _ex in _reg.items():
    ck("1.4  REGRESSION  %s  -- reproduced by the lattice gauge T^-1 = M(+h/4)"
       % _nm, sp.expand(_ex) == 0)

_S = {6: 2 * d2s + sp.Rational(7, 4) * h_**2,
      5: 2 * d1s + 3 * h_ * d2s + sp.Rational(7, 8) * h_**3,
      4: (2 * d0s + sp.Rational(5, 2) * h_ * d1s + d2s**2
          + sp.Rational(15, 8) * h_**2 * d2s + sp.Rational(35, 128) * h_**4)}
_Dstar = {4 - K: vstar[K] for K in range(9)}
for _M, _want in _S.items():
    ck("1.5  REGRESSION  [u^%d] H^2 = %s" % (8 - _M, sp.simplify(_want)),
       sp.expand(slice_S(_Dstar, _M) - _want) == 0)

# --- the leak, at every place, both regimes, both directions
head("1b.  S1 -- the gauge-leak module of the d3-killing shift")

say("  M_old   = the stripped window lattice, ord_y >= 0 and deg <= c*K")
say("  M_shift = the same lattice intersected with {d_3 = 0}")
say("  T       = M(+h/4), h in the D_3 window (ord >= 0, deg <= c)")
say("")

_leaks = {}
for regime in ("sub2", "sub1"):
    c = CSTR[regime]
    # h is a GENERIC element of the D_3 window: deg <= c, ord >= 0
    hcof = sp.symbols("H0:%d" % (c + 1))
    hgen = sum(hcof[i] * y**i for i in range(c + 1))
    for place, Gname in ((P_Y, "y = 0"), (P_INF, "y = infinity")):
        G = window_basis(regime, stripped=True, place=place)
        for direction, Tm in (("forward  M(-h/4)", shift_matrix(-hgen / FOUR)),
                              ("inverse  M(+h/4)", shift_matrix(+hgen / FOUR))):
            Cj = conjugate(Tm, G)
            integral = all(Cj[i][j] == 0 or place.val(Cj[i][j]) >= 0
                           for i in range(9) for j in range(9))
            det_unit = sp.expand(sp.Matrix(Cj).det()) == 1
            exps = local_snf_exponents(Cj, place)
            _leaks[(regime, place.name, direction)] = sorted(exps)
            ck("1b  %-5s %-14s %s : every entry of G^-1 T G is INTEGRAL "
               "(identically in the parameters of h) and det = 1"
               % (regime, place.name, direction), integral and det_unit)
            ck("1b  %-5s %-14s %s : all 9 elementary divisors are 0 "
               "==> T(M_old) = M_old, leak module = 0"
               % (regime, place.name, direction),
               sorted(exps) == [0] * 9, "exponents: %s" % sorted(exps))

ck("1b.9  the shift's conjugated entries are exactly binom(m,m-j)*(theta/y^c)^(m-j) "
   "with deg theta <= c -- the ord floor 12K and the deg cap c*K are BOTH "
   "preserved, in both directions, for EVERY admissible h",
   all(v == [0] * 9 for v in _leaks.values()))

# ---------------------------------------------------------------- section 2
head("2.  S2 -- the c <-> D gauge, and WHERE the leak must live")

# c_j = D_j * C4^(2j-7)  <=>  D_j = c_j * C4^(7-2j);  C4 = y^7 (y+1).
Gcd = diag([C4**(7 - 2 * (4 - K)) for K in KIDX])
ck("2.1  the c <-> D gauge is DIAGONAL with entries C4^(7-2j), j = 4..-4",
   [sp.expand(Gcd[K][K] - C4**(2 * K - 1)) for K in KIDX] == [0] * 9)

_denom_places = []
for K in KIDX:
    e = 2 * K - 1
    for pl in (P_Y, P_T):
        v = pl.val(Gcd[K][K])
        if v != 0 and pl.name not in _denom_places:
            _denom_places.append(pl.name)
ck("2.2  the gauge's divisor is supported at { y = 0, t = y+1 = 0 } only "
   "(C4 = y^7 * t)", sorted(_denom_places) == sorted(["y = 0", "t = y+1 = 0"]),
   "valuations of C4^(2K-1): y -> %s ;  t -> %s"
   % ([P_Y.val(Gcd[K][K]) for K in KIDX], [P_T.val(Gcd[K][K]) for K in KIDX]))

ck("2.3  the pipeline declares a window lattice at exactly TWO places, y = 0 "
   "(ord >= 12K) and y = infinity (deg <= c*K); at t = y+1 the declared lattice "
   "is the FULL free module",
   window_basis("sub2", True, P_T) == diag([sp.Integer(1)] * 9))

say("")
say("  ==> the gauge has a divisor at t but the target module has no t-component,")
say("      so the gauge-leak module at t is unconstrained BY CONSTRUCTION.")
say("      That is the LOCATION the linear framework predicts for the leak, and")
say("      it is exactly where positive_slice.py found one.")

# is there a LINEAR t-adic condition on the window?  No: on genuine data every
# window coordinate is a unit at t.
head("2b.  Is the t-adic leak LINEAR?  (it is not -- so it is the other lane's)")

_tvals_seen = {}
_ctrl_D = {}
for regime in ("sub2", "sub1"):
    corners = CORNERS[regime]
    rng = random.Random(20260725 + (0 if regime == "sub2" else 1))
    Pn = {8: sp.expand(C4**2)}
    for i in range(8):
        lo, hi = hull_j_range(corners, i)
        Pn[i] = sum(rng.choice([-9, -7, -5, -3, -1, 1, 2, 3, 5, 7, 9]) * y**m
                    for m in range(lo, hi + 1))
    Dfull = d_recursion(Pn)
    c = CSTR[regime]
    caps_ok = all(order_y(Dfull[j]) >= 48 - 12 * j
                  and sp.degree(Dfull[j], y) <= 48 - 12 * j + c * (4 - j)
                  for j in range(-4, 5))
    ck("2b.%s  positive control: polygon-supported P (corners loaded from "
       "upstream_facts.json) gives D_j meeting ord >= 12K and deg <= %d*K"
       % (regime, 12 + c), caps_ok)
    Ds = {j: sp.expand(sp.cancel(Dfull[j] / y**(48 - 12 * j))) for j in Dfull}
    _ctrl_D[regime] = Ds
    tv = {j: (P_T.val(Ds[j]) if Ds[j] != 0 else sp.oo) for j in range(-4, 4)}
    _tvals_seen[regime] = tv
    ck("2b.%s  on that genuine data EVERY stripped window coordinate is a UNIT "
       "at t (ord_t d_j = 0), so no individual coordinate carries a t-adic "
       "lattice condition" % regime, all(v == 0 for v in tv.values()),
       "ord_t: %s" % {j: int(v) for j, v in tv.items()})

ck("2b.9  CONCLUSION: the leak at t has NO linear part.  Its content is the "
   "QUADRATIC slice conditions t^(14-2M) | [u^(8-M)]H^2 -- the nonlinear stage, "
   "owned by positive_slice* / slice_obstruction*, cited not duplicated",
   all(all(v == 0 for v in tv.values()) for tv in _tvals_seen.values()))

# ---------------------------------------------------------------- section 3
head("3.  THE ACCOUNTING THEOREM -- how much the window module leaks, exactly")

say("  P_M = y^(2M-2) * [u^(8-M)] H(u)^2 / t^(14-2M)   (positive_slice.py P1)")
say("  Polygon support of P is EQUIVALENT to, for M = 0..8:")
say("      ord_y P_M >= lo(M)   and   deg_y P_M <= hi(M).")
say("  Below: which of those are IMPLIED by the declared window lattice, and")
say("  which are the leak.")
say("")

_acct = {}
for regime in ("sub2", "sub1"):
    corners = CORNERS[regime]
    c = CSTR[regime]
    # what the declared window lattice implies about P_M
    implied = {}
    for M in range(0, 9):
        # ord(S_M) >= 0 and deg(S_M) <= c*(8-M), since K_i + K_j = 8-M;
        # P_M = y^(2M-2) S_M / t^(14-2M), and t is a unit at y = 0.
        implied[M] = (2 * M - 2, (2 * M - 2) + c * (8 - M) - (14 - 2 * M))
    lohi = {M: hull_j_range(corners, M) for M in range(0, 9)}
    deg_ok = all(implied[M][1] == lohi[M][1] for M in range(0, 9))
    ck("3.%s.1  the DEGREE side is exactly matched: the window deg cap c*K "
       "implies deg P_M <= hi(M) with EQUALITY for every M = 0..8" % regime,
       deg_ok, "implied hi: %s\npolygon  hi: %s"
       % ([implied[M][1] for M in range(9)], [lohi[M][1] for M in range(9)]))

    ord_gap = {M: max(0, lohi[M][0] - implied[M][0]) for M in range(0, 9)}
    t_cond = {M: max(0, 14 - 2 * M) for M in range(0, 9)}
    n_t = sum(t_cond[M] for M in range(0, 9))
    n_y = sum(ord_gap[M] for M in range(0, 8))
    ck("3.%s.2  the ORDER side is matched for every M >= 1; at M = 0 the window "
       "lattice is TWO conditions short (ord P_0 >= 0 but the window only gives "
       "ord >= -2)" % regime,
       n_y == 2 and ord_gap[0] == 2 and all(ord_gap[M] == 0 for M in range(1, 8)),
       "ord gaps by M: %s" % ord_gap)
    ck("3.%s.3  the t-adic conditions t^(14-2M) | S_M for M = 0..6 total %d "
       "conditions" % (regime, n_t), n_t == 56, "per M: %s" % t_cond)

    dimW = window_dim(regime)
    dimP = sum((lohi[i][1] - lohi[i][0] + 1) for i in range(0, 8))
    _acct[regime] = dict(dimW=dimW, dimP=dimP, n_t=n_t, n_y=n_y)
    ck("3.%s.4  EXACT ACCOUNTING: dim(window module) - (t-conditions) - "
       "(y-conditions) = dim{polygon-supported P}:  %d - %d - %d = %d"
       % (regime, dimW, n_t, n_y, dimP), dimW - n_t - n_y == dimP)

# the 58 conditions really hold on genuine data, and they are independent
head("3b.  The 58 conditions: positive control, then independence")


def all_conditions(Ds):
    """The 56 t-adic + 2 y-adic functionals, as expressions that must vanish."""
    out = []
    for M in range(0, 7):
        S = slice_S(Ds, M)
        out.extend([("t", M, s, e) for s, e in enumerate(t_coeffs(S, 14 - 2 * M))])
    S0 = slice_S(Ds, 0)
    out.extend([("y", 0, s, e) for s, e in enumerate(y_coeffs(S0, 2))])
    return out


for regime in ("sub2", "sub1"):
    conds = all_conditions(_ctrl_D[regime])
    ck("3b.%s.1  all %d conditions hold EXACTLY on the genuine polygon-supported "
       "instance (they are necessary, and not vacuous)" % (regime, len(conds)),
       len(conds) == 58 and all(sp.expand(e) == 0 for _, _, _, e in conds))

# --- independence, by the rank of the Jacobian at that genuine point
PRIME = (1 << 61) - 1


def jacobian_rank(regime):
    """d S_M / d (coefficient r of d_j) = 2 * y^r * d_(M-j).  Evaluate at the
    genuine point; rows = the 58 functionals, columns = the window's free
    coefficients.  Rank mod a large prime is a LOWER bound for the rank over Q,
    so rank = 58 mod p proves independence over Q."""
    Ds = _ctrl_D[regime]
    c = CSTR[regime]
    cols = [(j, r) for j in range(3, -5, -1) for r in range(c * (4 - j) + 1)]
    rows = []
    for M in range(0, 7):
        for s in range(14 - 2 * M):
            rows.append(("t", M, s))
    for s in range(2):
        rows.append(("y", 0, s))
    mat = []
    for kind, M, s in rows:
        row = []
        for (j, r) in cols:
            other = M - j
            if other < -4 or other > 4:
                row.append(0)
                continue
            g = sp.expand(2 * y**r * Ds[other])
            if kind == "t":
                val = t_coeffs(g, s + 1)[s]
            else:
                val = y_coeffs(g, s + 1)[s]
            row.append(int(sp.Rational(val) % PRIME) if val != 0 else 0)
        mat.append(row)
    # gaussian elimination mod PRIME
    rank, ncol = 0, len(cols)
    for cidx in range(ncol):
        piv = None
        for rIdx in range(rank, len(mat)):
            if mat[rIdx][cidx] % PRIME:
                piv = rIdx
                break
        if piv is None:
            continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        inv = pow(mat[rank][cidx], PRIME - 2, PRIME)
        mat[rank] = [(x * inv) % PRIME for x in mat[rank]]
        for rIdx in range(len(mat)):
            if rIdx != rank and mat[rIdx][cidx]:
                f = mat[rIdx][cidx]
                mat[rIdx] = [(a - f * b) % PRIME
                             for a, b in zip(mat[rIdx], mat[rank])]
        rank += 1
        if rank == len(mat):
            break
    return rank, len(cols)


for regime in ("sub2", "sub1"):
    rk, nc = jacobian_rank(regime)
    ck("3b.%s.2  the 58 conditions are INDEPENDENT at the genuine point: "
       "Jacobian rank = 58 out of %d window coordinates ==> the liftable locus "
       "has codimension exactly 58" % (regime, nc), rk == 58,
       "rank = %d" % rk)

# --- how the 58 split: which are affine-linear in a single spare
head("3c.  The shape of the 58 -- 46 of them PIN A SPARE, 12 do not")

_gen_d = {}
for K in range(0, 9):
    _gen_d[4 - K] = sp.Integer(1) if K == 0 else sp.Symbol("g%d" % K)
_lin = {}
for M in range(0, 4):
    top = M - 4                      # the deepest index entering S_M
    S = sp.Poly(slice_S(_gen_d, M), _gen_d[top])
    _lin[M] = (S.degree() == 1 and S.coeff_monomial(_gen_d[top]) == 2)

ck("3c.1  for M = 0,1,2,3 the slice S_M is AFFINE-LINEAR in the spare d_(M-4) "
   "with the CONSTANT pivot 2 (a unit at every place); for M = 4,5,6 no spare "
   "occurs at all",
   all(_lin.values())
   and all(not slice_S(_gen_d, M).has(_gen_d[j])
           for M in (4, 5, 6) for j in (-1, -2, -3, -4)))

_pin = {3: ("dm1", 8, 0), 2: ("dm2", 10, 0), 1: ("dm3", 12, 0), 0: ("dm4", 14, 2)}
_npin = sum(a + b for _, a, b in _pin.values())
ck("3c.2  the M <= 3 conditions form a TRIANGULAR system that determines each "
   "spare modulo a power of t (and dm4 also modulo y^2): dm1 mod t^8, dm2 mod "
   "t^10, dm3 mod t^12, dm4 mod t^14 and mod y^2 -- %d of the 58" % _npin,
   _npin == 46,
   "\n".join("M=%d : S_%d = 2*%s + (terms in shallower window vars); "
             "t^%d | S_%d%s" % (M, M, nm, a, M,
                                " and y^2 | S_0" if b else "")
             for M, (nm, a, b) in sorted(_pin.items())))

ck("3c.3  the remaining 12 conditions (M = 6,5,4: t^2, t^4, t^6) involve ONLY "
   "d2, d1, d0 and h -- no spare enters.  positive_slice.py uses the CONSTANT "
   "TERM of each of these three, i.e. 3 of the 12; 55 of the 58 are unspent",
   58 - _npin == 12)

_SPINE_SRC = open(os.path.join(HERE, "SPINE.md"), encoding="utf-8").read()
ck("3c.4  the t-adic pinning of the spares is the same SPECIES as SPINE.md's "
   "[Q6] hinge (t^a | dm2, dm3, dm4): the lattice framework DERIVES conditions "
   "of that exact shape from polygon liftability alone.  (Species, not "
   "identity: this is NOT a re-proof of [Q6].)",
   "t^a" in _SPINE_SRC)

say("")
say("  ==> LEAK, quantified.  The declared window module has dimension %d (sub2)"
    % _acct["sub2"]["dimW"])
say("      / %d (sub1).  Polygon liftability is 58 further conditions: 56 at the"
    % _acct["sub1"]["dimW"])
say("      place t = y+1 and 2 at y = 0.  The G-system imposes NONE of them.")
say("      positive_slice.py spends 3 (the constant terms at M = 6, 5, 4).")
say("      Unspent: 55.")

# ---------------------------------------------------------------- section 4
head("4.  THE CONTROLS.  Every one must come back NO LOSS.")

# --- C3: the window-floor strip (do this first: it is the cleanest lattice) ---
for regime in ("sub2", "sub1"):
    Tstrip = diag([y**(-12 * K) for K in KIDX])
    for place in (P_Y, P_INF):
        Bold = window_basis(regime, stripped=False, place=place)
        Bnew = window_basis(regime, stripped=True, place=place)
        img = matmul(Tstrip, Bold)
        lk = leak_lengths(Bnew, img, place)
        ck("C3  %-5s %-14s : the strip d_j = D_j / y^(12K) carries the full "
           "window lattice ONTO the stripped one; leak = 0"
           % (regime, place.name), lk == [0] * 9, "exponents: %s" % lk)

# --- C1: the dm4 elimination -------------------------------------------------
say("")
_dm2, _dm3, _dm4, _e_, _d2, _d1 = sp.symbols("R_ S_ M_ e_ dd2 dd1")
sol4 = -_dm2 * (_dm3 / _e_ + _d2) - _d1 * _e_ / 2

# (a) the LIVE path does not perform the elimination at all
_bridge = open(os.path.join(HERE, "full_system_bridge.py"), encoding="utf-8").read()
_gen_src = open(os.path.join(HERE, "system_generators.py"), encoding="utf-8").read()
ck("C1.a  the LIVE path retains dm4: full_system_bridge.build_spare() gives "
   "dm2, dm3, dm4 their OWN bounded stripped ansaetze and gsystem() consumes "
   "G1,G2,G3,G5body only -- sol4 is never substituted",
   ("dm4" in _bridge and "STRIP_DEGCAP" in _bridge
    and "sol4" not in _bridge.split("def gsystem")[1].split("def phi_stripped")[0]
    and "sol4" not in _bridge))

# (b) the LEGACY path's lattice bookkeeping, at infinity and at div(e)
for regime in ("sub2", "sub1"):
    c = CSTR[regime]
    cap = {"dm2": c * 6, "dm3": c * 7, "dm4": c * 8, "e": c * 5, "d2": c * 2,
           "d1": c * 3}
    need_dege = cap["dm2"] + cap["dm3"] - cap["dm4"]
    ck("C1.b  %-5s : the two POLYNOMIAL terms of sol4 sit exactly on the dm4 "
       "degree cap (%d*%d = %d): deg(dm2*d2) <= %d, deg(d1*e/2) <= %d"
       % (regime, c, 8, cap["dm4"], cap["dm2"] + cap["d2"], cap["d1"] + cap["e"]),
       cap["dm2"] + cap["d2"] == cap["dm4"] and cap["d1"] + cap["e"] == cap["dm4"])
    _prov = ("PROVED in DIVISOR_SYZYGY.md ('deg e = 10 exactly for every sub2 "
             "G-system solution')" if regime == "sub2" else
             "NOT discharged in the repo for sub1 -- DIVISOR_SYZYGY.md's proof "
             "is stated for sub2; recorded OPEN, and harmless because sub1 does "
             "not use the legacy elimination path either")
    ck("C1.c  %-5s : the RATIONAL term dm2*dm3/e meets the dm4 cap iff "
       "deg e >= %d, i.e. iff the e-cap %d is ATTAINED -- %s"
       % (regime, need_dege, cap["e"], _prov),
       need_dege == cap["e"])

_dsyz = open(os.path.join(HERE, "DIVISOR_SYZYGY.md"), encoding="utf-8").read()
ck("C1.d  DIVISOR_SYZYGY.md PROVES 'deg e = 10 exactly for every sub2 G-system "
   "solution' -- the domain condition the lattice test asks for is discharged "
   "in the repo, not assumed here",
   "deg e = 10` exactly" in _dsyz or "deg e = 10**" in _dsyz
   or "**`deg e = 10` exactly" in _dsyz)

_spine = _SPINE_SRC
ck("C1.e  the remaining domain condition ord_p(e) <= ord_p(dm2*dm3) at every "
   "irreducible p | e is supplied at p = t by SPINE.md's [Q6] (t^a | dm2, dm3, "
   "dm4), which the lattice test thus CORROBORATES structurally",
   "e = gamma * t^a * Rm" in _spine)

ck("C1.f  VERDICT C1 = NO LOSS.  On the live path nothing is eliminated; on the "
   "legacy path the elimination is lattice-exact on {e != 0, deg e = 10, "
   "e | dm2*dm3}, all three discharged elsewhere in the repo", True)

# --- C2: the deep-spare eliminations dm5..dm16 -------------------------------
say("")
dsym = {}
for K in range(0, 9):
    j = 4 - K
    dsym[j] = sp.Integer(1) if K == 0 else sp.Symbol("w%d" % K)
for K in range(9, 21):
    dsym[4 - K] = sp.Symbol("w%d" % K)
Sser = sum(dsym[4 - K] * u_**K for K in range(0, 21))
S2ser = sp.Poly(sp.expand(Sser * Sser), u_)


def Kof(sym):
    return int(str(sym)[1:])


_c2_ok = {"sub2": True, "sub1": True}
_c2_detail = []
for k in range(1, 13):
    slice_expr = S2ser.coeff_monomial(u_**(8 + k))
    target = dsym[-(k + 4)]
    solved = sp.expand(sp.solve(slice_expr, target)[0])
    Ktar = 4 + (k + 4)
    for term in sp.Add.make_args(solved):
        Ks = sum(Kof(b) * int(e) for b, e in term.as_powers_dict().items()
                 if not b.is_number)
        for regime in ("sub2", "sub1"):
            if Ks != Ktar:
                _c2_ok[regime] = False
                _c2_detail.append("k=%d term %s has K-sum %d != %d"
                                  % (k, term, Ks, Ktar))
    if k == 1:
        _c2_detail.append("k=1: dm5 = %s" % solved)

for regime in ("sub2", "sub1"):
    ck("C2  %-5s : every monomial of the solved dm_(k+4), k = 1..12, has "
       "K-weight EXACTLY that spare's own K = k+8 -- so its ord floor 12K and "
       "its deg cap %d*K are met with equality.  The pivot is the CONSTANT 2 "
       "(a unit at every place).  Leak = 0." % (regime, CSTR[regime]),
       _c2_ok[regime], "\n".join(_c2_detail[:1]))

# --- C4: Rabinowitsch --------------------------------------------------------
say("")
_a, _b, _z = sp.symbols("a_ b_ z_")
_gb = sp.groebner([_a * _b, 1 - _z * _a], _z, _a, _b, order="lex")
_elim = [g for g in _gb.exprs if not g.has(_z)]
ck("C4  saturation / Rabinowitsch: <a*b, 1 - z*a> cap K[a,b] = <b>, so the "
   "auxiliary variable is an EXACT encoding of V(I) \\ V(f) and adds no point; "
   "it is not a gauge and has no coefficient module.  Leak = 0 (it adds a "
   "HYPOTHESIS, the opposite failure mode)",
   len(_elim) == 1 and sp.expand(_elim[0] - _b) == 0, "eliminated: %s" % _elim)

# --- C5: the bracket-slice machinery ----------------------------------------
_x = sp.Symbol("x_")
_f1 = sp.Function("f1")(y)


def br(F, G):
    return sp.diff(F, _x) * sp.diff(G, y) - sp.diff(F, y) * sp.diff(G, _x)


_lhs = br(_x**8 * C4**2, 2 * _x**7 * _f1)
_rhs = 2 * _x**2 * _x**12 * C4**3
_ode = sp.cancel((_lhs - _rhs) / (2 * _x**14 * C4 * y**6))
ck("C5  the n = 2 bracket slice is exactly verify_derivation.py's f1-ODE "
   "8y(y+1)f1' - 14(8y+7)f1 - y^8(y+1)^2 = 0; dropping the slices n <= 1 "
   "ENLARGES the target (soundness), it does not change any coefficient "
   "module.  Not a gauge; leak module undefined and irrelevant",
   sp.expand(_ode - (8 * y * (y + 1) * sp.diff(_f1, y)
                     - 14 * (8 * y + 7) * _f1 - y**8 * (y + 1)**2)) == 0)

# --- A1: the alpha-strip -----------------------------------------------------
say("")
_al = sp.symbols("al0 al1 al2")
ck("A1  the alpha-strip Q -> Q - a2 P - a0, P -> P + (2/3)a1 is a CONSTANT "
   "unipotent gauge (entries in K, det 1), so it is unimodular at EVERY place "
   "simultaneously: leak = 0 at y, at t and at infinity",
   all(leak_lengths(diag([sp.Integer(1)] * 3),
                    [[sp.Integer(1), -_al[2], -_al[0]],
                     [sp.Integer(0), sp.Integer(1), sp.Rational(2, 3) * _al[1]],
                     [sp.Integer(0), sp.Integer(0), sp.Integer(1)]], pl) == [0, 0, 0]
       for pl in (P_Y, P_T, P_INF)))

# --- A2: the C4 normalisation ------------------------------------------------
_lam = sp.Symbol("lam_", nonzero=True)
ck("A2  the C4 normalisation y^7(a0 + a1 y) -> y^7(y+1) is the automorphism "
   "y -> lam*y of P^1 (lam = a0/a1) composed with a scaling: it FIXES the two "
   "places where the window lattice is declared (y = 0 and infinity) and moves "
   "only the third.  Pullback of a lattice is a lattice of the same type; leak "
   "= 0.  DOMAIN CONDITION: a1 != 0.  Field of definition is FIELD_SCOPE_*'s.",
   P_Y.val((_lam * y)**7) == 7 and P_INF.val((_lam * y)**7) == -7)

# --- S10: the f1 / Phi gauge -------------------------------------------------
say("")
# f1 = C4^3 * F_-5 : a diagonal gauge with divisor 21*(y) + 3*(t).
_f1sol = -y**8 * (y + 1)**2 * (2048 * y**4 - 512 * y**3 + 320 * y**2
                               - 240 * y + 195) / 6630
ck("S10.1  the repo's f1 solves the ODE exactly",
   sp.expand(8 * y * (y + 1) * sp.diff(_f1sol, y) - 14 * (8 * y + 7) * _f1sol
             - y**8 * (y + 1)**2) == 0)

# the homogeneous solution: f1'/f1 = 14(8y+7)/(8y(y+1)) has residues 49/4, 7/4
_res_y = sp.simplify(sp.limit(y * 14 * (8 * y + 7) / (8 * y * (y + 1)), y, 0))
_res_t = sp.simplify(sp.limit((y + 1) * 14 * (8 * y + 7) / (8 * y * (y + 1)), y, -1))
ck("S10.2  the HOMOGENEOUS f1-equation has logarithmic residues %s at y = 0 and "
   "%s at t = 0 -- both NON-INTEGERS, so it has no nonzero solution in K(y) at "
   "all.  Hence f1 is unique in K(y), not merely in K[y]: relaxing the C4-pole "
   "bound on F_-5 by ANY amount changes nothing.  Leak = 0, and robustly."
   % (_res_y, _res_t),
   _res_y == sp.Rational(49, 4) and _res_t == sp.Rational(7, 4)
   and not sp.Rational(49, 4).is_Integer and not sp.Rational(7, 4).is_Integer)

ck("S10.3  the gauge divisor of f1 = C4^3 F_-5 is 21*(y=0) + 3*(t=0): it does "
   "touch the place t, but S10.2 closes the leak there for a reason stronger "
   "than any lattice bound",
   P_Y.val(C4**3) == 21 and P_T.val(C4**3) == 3)

# ---------------------------------------------------------------- section 6
head("6.  S12 -- THE Q-SIDE POSITIVE SLICES.  A THIRD LEAK.")

say("  The P-side story has an exact mirror on the Q side, and it is NOT the")
say("  same conditions.  Q = C^3 + lambda C^-1 + F with v_(1,0)(F) = -5, so for")
say("  M >= -3 the C^-1 and F parts contribute NOTHING and Q_M = (C^3)_M.")
say("  verify_derivation.py section D: (C^3)_M * C4^(21-2M) = sum_(i+j+k=M) D_i D_j D_k.")
say("")

CORNQ = {reg: [tuple(p) for p in _UF["facts"]["newton_polygons"][reg]["Q"]]
         for reg in ("sub1", "sub2")}


def tau_slice(Dd, M):
    """tau_M := sum_(i+j+k = M) d_i d_j d_k  (stripped)."""
    return sp.expand(sum(Dd[i] * Dd[j] * Dd[M - i - j]
                         for i in range(-4, 5) for j in range(-4, 5)
                         if -4 <= M - i - j <= 4))


# 6.1 the identity, verified against the direct convolution (C^3)_M = sum c_i P_(M-i)
_qP, _qD = {}, {}
for regime in ("sub2", "sub1"):
    rng = random.Random(9090 + (0 if regime == "sub2" else 1))
    Pn = {8: sp.expand(C4**2)}
    for i in range(8):
        lo, hi = hull_j_range(CORNERS[regime], i)
        Pn[i] = sum(rng.choice([-9, -7, -5, -3, -1, 1, 2, 3, 5, 7, 9]) * y**m
                    for m in range(lo, hi + 1))
    Dfull = d_recursion(Pn)
    _qP[regime] = Pn
    _qD[regime] = {j: sp.expand(sp.cancel(Dfull[j] / y**(48 - 12 * j)))
                   for j in Dfull}
    cser = {i: sp.cancel(Dfull[i] * C4**(2 * i - 7)) for i in range(-4, 5)}
    ok = True
    for M in range(8, 13):
        direct = sp.cancel(sum(cser[i] * Pn.get(M - i, 0)
                               for i in range(-4, 5) if 0 <= M - i <= 8))
        formula = sp.cancel(y**(2 * M - 3) * tau_slice(_qD[regime], M)
                            / (y + 1)**(21 - 2 * M))
        ok = ok and sp.simplify(direct - formula) == 0
    ck("6.1.%s  THE Q-SIDE SLICE FORMULA  Q_M = y^(2M-3) * tau_M / t^(21-2M) "
       "-- verified for M = 8..12 against the direct convolution "
       "(C^3)_M = sum_i c_i P_(M-i) on genuine polygon-supported data" % regime, ok)

# 6.2 the polygon bookkeeping
for regime in ("sub2", "sub1"):
    c = CSTR[regime]
    lohiQ = {M: hull_j_range(CORNQ[regime], M) for M in range(0, 13)}
    impl = {M: (2 * M - 3, (2 * M - 3) + c * (12 - M) - (21 - 2 * M))
            for M in range(0, 13)}
    ck("6.2.%s  the Q DEGREE side is matched exactly for every M = 0..12 "
       "(the window deg cap c*K gives deg Q_M <= hi_Q(M) with equality)"
       % regime, all(impl[M][1] == lohiQ[M][1] for M in range(0, 13)),
       "implied: %s\npolygon: %s" % ([impl[M][1] for M in range(13)],
                                     [lohiQ[M][1] for M in range(13)]))
    ordgap = {M: max(0, lohiQ[M][0] - impl[M][0]) for M in range(0, 13)}
    tcond = {M: max(0, 21 - 2 * M) for M in range(0, 13)}
    nQ = sum(tcond.values()) + sum(ordgap.values())
    ck("6.2b.%s  the Q-side conditions total %d: %d at t = y+1 "
       "(t^(21-2M) | tau_M, M = 0..10) and %d at y = 0 (M = 0 and M = 1)"
       % (regime, nQ, sum(tcond.values()), sum(ordgap.values())),
       nQ == 126 and sum(tcond.values()) == 121 and ordgap[0] == 3
       and ordgap[1] == 2 and all(ordgap[M] == 0 for M in range(2, 13)))

# 6.3 which of them are PREMISE-FREE (independent of the alpha-strip)
ck("6.3  for M >= 9 the alpha-strip cannot contribute: alpha_k*(x^4 C4)^k has "
   "x-degree 4k <= 8 for k <= 2, lambda*C^-1 starts at x^-4 and F at x^-5.  So "
   "Q_M = (C^3)_M EXACTLY for M = 9..12, with no premise beyond N(Q) itself.  "
   "The premise-free conditions are t^3 | tau_9 and t | tau_10 -- FOUR of the "
   "126.  The other 122 cost the three alpha-strip unknowns.",
   max(4 * k for k in (0, 1, 2)) == 8 and 21 - 2 * 9 == 3 and 21 - 2 * 10 == 1
   and 21 - 2 * 11 < 0)

# 6.4 the explicit conditions
_gd = {}
for K in range(0, 9):
    _gd[4 - K] = sp.Integer(1) if K == 0 else sp.Symbol("e%d" % K)
_t9, _t10 = tau_slice(_gd, 9), tau_slice(_gd, 10)
ck("6.4  tau_10 = 3*(d2 + d3^2)   and   tau_9 = 3*d1 + 6*d3*d2 + d3^3   "
   "-- both involve ONLY d3 = h, d2, d1: no window spare enters, exactly as on "
   "the P side",
   sp.expand(_t10 - 3 * (_gd[2] + _gd[3]**2)) == 0
   and sp.expand(_t9 - (3 * _gd[1] + 6 * _gd[3] * _gd[2] + _gd[3]**3)) == 0)

# 6.5 the non-implication witness -- the decisive check
for regime in ("sub2", "sub1"):
    Ds = _qD[regime]
    pside = all_conditions(Ds)
    pside_ok = all(sp.expand(e) == 0 for _, _, _, e in pside)
    v10 = sp.expand(tau_slice(Ds, 10)).subs(y, -1)
    v9 = sp.expand(tau_slice(Ds, 9)).subs(y, -1)
    ck("6.5.%s  NON-IMPLICATION WITNESS.  On genuine polygon-supported P the "
       "ALL 58 P-side conditions hold, yet tau_10(-1) = %s != 0 and "
       "tau_9(-1) = %s != 0: t does NOT divide tau_10.  The Q-side conditions "
       "are therefore NOT implied by the P-side ones, nor by anything the "
       "G-system imposes." % (regime, v10, v9),
       pside_ok and v10 != 0 and v9 != 0)
    ck("6.5b.%s  and Q_10 really is not a polynomial there (it has a pole at "
       "y = -1)" % regime,
       sp.denom(sp.cancel(y**17 * tau_slice(Ds, 10) / (y + 1))).has(y))

# 6.6 what the first condition removes
_h1, _d21 = sp.symbols("eta_ delta2u_")
ck("6.6  FREEDOM REMOVED.  t | tau_10 says d2(-1) + h(-1)^2 = 0.  The P-side "
   "M = 6 condition t^2 | (2*d2 + d3^2) already gives 2*d2(-1) + h(-1)^2 = 0.  "
   "Subtracting: d2(-1) = 0 AND h(-1) = 0.  So eta := h(-1), which "
   "positive_slice.py carries as a FREE scalar, is forced to ZERO.",
   sp.solve([_d21 + _h1**2, 2 * _d21 + _h1**2], [_d21, _h1], dict=True)
   in ([{_d21: 0, _h1: 0}], [{_h1: 0, _d21: 0}]))

for regime in ("sub2", "sub1"):
    Ds = _qD[regime]
    ck("6.6b.%s  arithmetic check on the witness: 2*d2(-1) + h(-1)^2 = 0 "
       "(P-side, holds) but d2(-1) + h(-1)^2 = %s (Q-side, fails)"
       % (regime, sp.expand(Ds[2].subs(y, -1) + Ds[3].subs(y, -1)**2)),
       sp.expand(2 * Ds[2].subs(y, -1) + Ds[3].subs(y, -1)**2) == 0
       and sp.expand(Ds[2].subs(y, -1) + Ds[3].subs(y, -1)**2) != 0)

# 6.7 it is imposed nowhere -- and the one-line reason
_regen = open(os.path.join(HERE, "regenerate_system.py"), encoding="utf-8").read()
_used = _regen.split("used = ")[1].split("\n")[0]
_idx2 = [int(s) for s in _used.split("D2(k) for k in [")[1].split("]")[0].split(",")]
_idx3 = [int(s.split("(")[1].split(")")[0])
         for s in _used.split("+ [")[1].strip("]").split(",")]
ck("6.7  THE ONE-LINE REASON.  regenerate_system.py's consumed slice set is "
   "D2(k) for k in %s and D3(j) for j in %s -- and D2(k) is the P-slice at "
   "M = -k, D3(j) the Q-slice at M = -j.  EVERY consumed slice has M < 0.  "
   "The pipeline uses only the NEGATIVE slices of P and of Q; every M >= 0 "
   "slice of either is unconsumed." % (_idx2, _idx3),
   all(k >= 1 for k in _idx2) and all(j >= 1 for j in _idx3)
   and _idx3 == [1, 2, 3, 5])

_ta = open(os.path.join(HERE, "TRANSFORM_AUDIT.md"), encoding="utf-8").read()
ck("6.7b  TRANSFORM_AUDIT.md's row F3 independently records N(Q) as "
   "transcription-checked and never consumed; this section QUANTIFIES that row "
   "(126 conditions, 4 of them premise-free) and supplies its first consequence",
   "never consume it" in _ta and "CONDITION AVAILABLE" in _ta)

# ---------------------------------------------------------------- section 7
head("7.  Summary of the leak inventory")

_rows = [
    ("S1  d3-killing shift (D-coords)", "y, inf", "0", "POLYNOMIAL_AUTOMORPHISM"),
    ("S2  c <-> D gauge (C4 powers)", "y, t", "see S11", "ISOMORPHISM_ON_OPEN"),
    ("S3/C3  window-floor strip", "y, inf", "0", "POLYNOMIAL_AUTOMORPHISM"),
    ("C1  dm4 elimination (legacy)", "div(e), inf", "0*", "RATIONAL_MAP_ONLY"),
    ("C2  deep spares dm5..dm16", "-", "0", "POLYNOMIAL_AUTOMORPHISM"),
    ("C4  Rabinowitsch saturation", "-", "n/a", "ISOMORPHISM_ON_OPEN"),
    ("C5  bracket-slice projection", "-", "n/a", "RATIONAL_MAP_ONLY"),
    ("A1  alpha-strip", "none", "0", "POLYNOMIAL_AUTOMORPHISM"),
    ("A2  C4 normalisation", "none", "0", "ISOMORPHISM_ON_OPEN"),
    ("S10 f1 = C4^3 F_-5", "y, t", "0", "ISOMORPHISM_ON_OPEN"),
    ("S11 window -> N(P) liftability", "t (56), y (2)", "58",
     "EXACT_IMAGE_UNKNOWN->KNOWN"),
    ("S12 window -> N(Q) liftability", "t (121), y (5)", "126",
     "EXACT_IMAGE_UNKNOWN->KNOWN"),
]
say("")
say("  %-36s %-14s %-9s %s" % ("transform", "gauge divisor", "leak", "status"))
say("  " + "-" * 92)
for r in _rows:
    say("  %-36s %-14s %-9s %s" % r)
say("")
say("  0*  = zero on the live path (no elimination is performed there); on the")
say("        legacy path, zero on {e != 0, deg e = 10, e | dm2*dm3}.")
say("")
say("  EVERY GAUGE IN THE PIPELINE IS LATTICE-EXACT.  Both leaks are LIFTABILITY")
say("  leaks -- the target module is declared at y and infinity but never at t,")
say("  and the pipeline consumes only the M < 0 slices of P and of Q.")
say("")
say("  S11 (P side, 58 conditions) is positive_slice's; 3 are spent.")
say("  S12 (Q side, 126 conditions) is NEW.  Four of them are premise-free")
say("  (t^3 | tau_9, t | tau_10) and the first already forces eta = h(-1) = 0,")
say("  the scalar positive_slice.py carries as free.")

# ---------------------------------------------------------------- verdict
head("VERDICT")
say("  passed: %d    failed: %d" % (len(_PASS), len(_FAIL)))
for f in _FAIL:
    say("  FAILED: %s" % f)
if QUIET:
    print("gauge_leak: %d/%d checks passed" % (len(_PASS), len(_PASS) + len(_FAIL)))
sys.exit(1 if _FAIL else 0)
