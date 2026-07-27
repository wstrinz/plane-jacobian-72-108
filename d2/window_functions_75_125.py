#!/usr/bin/env python3
"""window_functions_75_125.py  (REPAIRED 2026-07-26; see PASSPORT_75_125_REPAIR.md)

THE WINDOW FUNCTIONS for (75,125) -- the arithmetic layer any (75,125) window
compiler needs, derived exactly (pure sympy).

*** WHAT THIS FILE USED TO SAY, AND WHY IT IS WRONG. ***

It presented "the PERIOD-12 window functions": a quasipolynomial lower (y-order)
cap L(w) = ceil(67w/12) with quasi-period q_window = 5a-3 = 12, together with an
AFFINE upper (y-degree) cap U(w) = 14w -- two distinct slopes, exactly as at
(72,108) (ord 12, deg 14).  Both rested on the (75,125) corner data
t=5, kappa=3, C=y^2(y^3+1), which came from reading GGV5's final chain corner
(7\\5,2) as chart data.  That dictionary is valid only on the retraction shape
b0 = l(a0-1); (5,20) fails it (20 != 4*4), and the chart exponent is
l = ceil(20/5) = 4.  See polygon_reduction.py sec.0b for the guard.

*** WHAT IS TRUE. ***

With t=4, kappa=2 and C = y a MONOMIAL:

    Phi_a = (1/a) y^(12a^2-10a+2),   M_a = 12a-7,   ord_y(Phi_a) = deg_y(Phi_a)

so (a=3): ord = deg = 80, M = 29, W_step = 80/29, q_window = 29 -- and 29 is
PRIME.  Three consequences, each of which kills a piece of the old story:

  (R1) PERIOD 12 IS REFUTED.  q_window = 12a-7 (17 at a=2, 29 at a=3), not
       5a-3 (7, 12).  So there is no period-7 -> period-12 structure at all.
       Both periods are prime, so the "fractional-denominator classes
       {2,3,4,6,12}" / "divisor lattice of the period" reading has no
       counterpart: the denominator sets are just {1,17} and {1,29}.
  (R2) THE AFFINE DEGREE CAP DOES NOT EXIST.  deg_slope := deg_y(Phi)/M = 80/29
       is NOT an integer, so U(w) = deg_slope*w is not integral and there is no
       affine y-degree cap.  CAPS_AUDIT.md sec.5's "deg_slope = 504/36 = 14" is
       false, not merely tautological.
  (R3) THE TWO-SLOPE CONE COLLAPSES TO A RAY.  Because Phi is a monomial,
       ord_y(Phi) = deg_y(Phi), so the ord-slope and the deg-slope are the SAME
       number 80/29 and the stripped slope lambda = deg_slope - W_step = 0.  The
       (72,108) picture -- a cone between an ord-lower ray of slope 12 and a
       deg-upper ray of slope 14, with a 2-unit strip to strip away -- has NO
       (75,125) counterpart.  Under the extreme-ray premise the caps would pinch
       (ord >= (80/29)w and deg <= (80/29)w with ord <= deg), which is
       satisfiable only when 29 | w.  That is not a window system; it is a
       demonstration that the premise does not transfer.

WHAT SURVIVES, and it is worth keeping:

  (S1) The general window_law(ordPhi, M, degPhi) arithmetic -- W_step in lowest
       terms, quasi-period q = denom, beta_m = (-alpha m) mod q, and the
       identity floor((alpha w + beta_m)/q) = ceil(alpha w / q).  This is pure
       arithmetic, valid for any (alpha, q), and it is what the ord-side carry
       obstruction in weight_lemma_75_125.py sec.C consumes.
  (S2) The class-interaction / 1-cocycle structure: classes add mod q and the
       beta offsets satisfy beta_{m1} + beta_{m2} = beta_{(m1+m2) mod q} + q*carry
       with carry in {0,1}, the carry being exactly the ceil-superadditivity
       defect.  Also pure arithmetic.
  (S3) The (72,108) integral limit q_window = 1, ord >= 12w, deg <= 14w --
       untouched, and now visibly the exceptional case.

THE ONE REMAINING PREMISE, unchanged in status: that *Phi realises the extreme
(minimal ord/weight) ray of the window cone*, inherited verbatim from (72,108).
Under it the lower cap is ord >= (alpha/q) w, hence ord >= ceil(alpha w/q) = L(w),
and beta_m = (-alpha m) mod q is pinned.  Whether the actual (75,125) window cone
dips below the 80/29 ray at some non-Phi weight is still not decidable from the
u-grading alone.  Note the premise is now WEAKER-LOOKING than before, because at
(72,108) Phi sat strictly inside a cone (204 < 238) whereas here it is a monomial
with no interior at all -- so "extreme ray" has less content here, and (R3) is
the honest reading of that.
"""
from sympy import Rational, gcd, ceiling, floor, denom, isprime


# ---------------------------------------------------------------------------
# Family constants (the F2 family: a = j+2, b = 2a-1, t=4, kappa=2, C = y,
# Phi_a = (1/a) y^(12a^2-10a+2)).  Source: f2_family_verify.py (repaired),
# C_SERIES_75_125.md, and the retraction guard in polygon_reduction.py.
# ---------------------------------------------------------------------------
T_CHART, KAPPA, ORD_C, DEG_C = 4, 2, 1, 1


def family(a):
    """Return the exact (b,t,kappa,jphi,M,ordPhi,degPhi) for F2-family rung a."""
    b = 2 * a - 1
    t = T_CHART
    kappa = KAPPA
    jphi = a * t - kappa - 1             # = 4a-3
    M = b * t + jphi                     # = 12a-7
    N = a * (t * (a + b) - (kappa + 1)) - 2 * b        # = (3a-2)(4a-1)
    rho = (b - a + 1 - 1) * ORD_C + 1                  # = a
    ordPhi = rho + N * ORD_C             # ord_y(Phi_a) = 12a^2-10a+2
    degPhi = rho + N * DEG_C             # == ordPhi: C is a MONOMIAL
    return dict(a=a, b=b, t=t, kappa=kappa, jphi=jphi, M=M, N=N,
                ordPhi=ordPhi, degPhi=degPhi)


# ---------------------------------------------------------------------------
# The window-cap arithmetic, derived from (ordPhi, M, degPhi) for ANY case.
#   W_step = ord_y(Phi)/M  in lowest terms  =  alpha / q
#   L(w) = floor((alpha w + beta_m)/q) = ceil(alpha w / q)   (lower y-order cap)
#   deg_slope = deg_y(Phi)/M ; the upper cap U(w) = deg_slope*w is AFFINE only
#   when deg_slope is an integer.  It is NOT for the F2 family: see (R2)/(R3).
#   beta_m = (-alpha m) mod q,   m = w mod q
# ---------------------------------------------------------------------------
def window_law(ordPhi, M, degPhi):
    """Derive the window-cap law from a case's Phi signature at slice M.

    REPAIRED 2026-07-26: this function used to ASSERT that deg_y(Phi)/M is an
    integer ("deg cap not affine" was a hard error).  That assertion encoded the
    (72,108) shape as a universal law, and it is false at (5,20).  We now REPORT
    `deg_affine` instead of asserting it, and flag the degenerate case where the
    deg-slope equals the ord-slope (which happens exactly when Phi is a monomial,
    i.e. when C is).
    """
    W = Rational(ordPhi, M)              # W_step in lowest terms
    alpha, q = W.p, W.q                  # numerator / quasi-period
    degR = Rational(degPhi, M)
    deg_affine = (degR.q == 1)
    return dict(W_step=W, alpha=alpha, q=q,
                deg_slope=(int(degR) if deg_affine else degR),
                deg_affine=deg_affine,
                slopes_coincide=(degR == W),
                lam=degR - W,            # the "stripped slope"; 0 iff coincident
                q_is_prime=bool(isprime(q)),
                beta=[(-alpha * m) % q for m in range(q)])


def L(w, alpha, q, beta):
    """Lower y-order cap (quasipolynomial floor form)  ==  ceil(alpha w / q)."""
    return (alpha * w + beta[w % q]) // q


def L_ceil(w, alpha, q):
    """Same lower cap, ceiling form (the tight integer bound ord >= alpha w/q)."""
    return int(ceiling(Rational(alpha * w, q)))


def U(w, deg_slope):
    """Upper y-degree cap.  ONLY meaningful when deg_slope is an integer.

    Raises otherwise -- and it IS non-integral for every F2 rung, which is
    result (R2).  Callers that want a number in the non-affine case want
    U_ray() and should read (R3) first.
    """
    if Rational(deg_slope).q != 1:
        raise ValueError(
            "no affine y-degree cap: deg_slope = %s is not an integer.  For the "
            "F2 family this is always the case (deg_y(Phi) = ord_y(Phi) because "
            "C = y is a monomial), so the (72,108)-style affine upper cap has no "
            "counterpart -- see (R2)/(R3) in the module docstring." % deg_slope)
    return int(deg_slope) * w


def U_ray(w, alpha, q):
    """The degenerate 'upper cap' when the two slopes coincide: floor(alpha w/q).

    Reported for completeness, NOT as a usable cap: together with the lower cap
    L(w) = ceil(alpha w/q) it pinches (L(w) > U_ray(w) unless q | w), which is
    the content of (R3) -- the cone has collapsed to a ray and the extreme-ray
    premise cannot support two independent slopes here.
    """
    return (alpha * w) // q


def q_window(a):
    """The window-denominator law: q_window = denom(W_step) = 12a-7 = M_a."""
    f = family(a)
    return denom(Rational(f["ordPhi"], f["M"]))


# ---------------------------------------------------------------------------
# The (75,125) instance (a=3) -- the target.
# ---------------------------------------------------------------------------
def target_75_125():
    f = family(3)
    f.update(window_law(f["ordPhi"], f["M"], f["degPhi"]))
    return f


# ---------------------------------------------------------------------------
# u-weights of the G-system symbols (from g_system_75_125.py: w(d_m)=t-m,
# w(Phi)=M).  Used for the residue-class occupancy / interaction table.
# ---------------------------------------------------------------------------
def generator_uweights(a):
    """u-weights bt+j of the forcing generators G_j (j=1..jphi, skip jphi-1)."""
    f = family(a)
    b, t, jphi = f["b"], f["t"], f["jphi"]
    js = [j for j in range(1, jphi + 1) if j != jphi - 1]
    return [b * t + j for j in js]


def state_uweights(a):
    """u-weights of the state variables d_{t-2}..d_0, e=dm1 = d_{-1}."""
    t = family(a)["t"]
    return [t - m for m in range(t - 2, -2, -1)]   # m = t-2 .. -1


def spare_uweights(a):
    """u-weights of the spare window unknowns dm2..dm_{(a-1)t} = d_-2..d_-(a-1)t."""
    t = family(a)["t"]
    return [t - m for m in range(-2, -((a - 1) * t) - 1, -1)]


def cls(w, q):
    """Residue class of a u-weight."""
    return w % q


def compose(m1, m2, q):
    """Class-interaction: multiplication adds u-weights -> classes add mod q."""
    return (m1 + m2) % q


def carry(m1, m2, beta, q):
    """The 1-cocycle carry: beta_{m1}+beta_{m2} = beta_{(m1+m2)%q} + q*carry."""
    return (beta[m1] + beta[m2] - beta[(m1 + m2) % q]) // q


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    r = target_75_125()
    print("=" * 74)
    print("WINDOW FUNCTIONS for (75,125)  [F2 family, a=3]   -- REPAIRED 2026-07-26")
    print("=" * 74)
    print("  Phi_3 = (1/3) y^%d   (a MONOMIAL: ord_y = deg_y)" % r["ordPhi"])
    print("  M = 12a-7 =", r["M"], "  N = (3a-2)(4a-1) =", r["N"])
    print("  alpha =", r["alpha"], "   q = q_window = 12a-7 =", r["q"],
          "  (prime: %s)" % r["q_is_prime"])
    print("  W_step = ord_y(Phi)/M =", r["W_step"], "(non-integral -> ord cap quasipolynomial)")
    print("  deg_slope = deg_y(Phi)/M =", r["deg_slope"],
          " affine:", r["deg_affine"], "  slopes coincide:", r["slopes_coincide"],
          "  lambda =", r["lam"])
    print("  beta_m (m=0..%d)      = %s" % (r["q"] - 1, r["beta"]))
    print("  lower cap  L(w) = floor((%d w + beta_m)/%d) = ceil(%d w/%d)"
          % (r["alpha"], r["q"], r["alpha"], r["q"]))
    print("  upper cap  NONE: deg_slope is not an integer  ->  see (R2)/(R3)")
    print()
    a_, q_, be_ = r["alpha"], r["q"], r["beta"]
    print("  Phi at M=%d:  L(%d)=%d = ord_y(Phi) = deg_y(Phi)  [equality; M = q]"
          % (r["M"], r["M"], L(r["M"], a_, q_, be_)))
    print("  q_window == M exactly, so class(Phi) = 0 and NO 0 < w < M has carry 0")
    zc = [w for w in range(1, r["M"]) if L(w, a_, q_, be_) + L(r["M"] - w, a_, q_, be_)
          - L(r["M"], a_, q_, be_) == 0]
    print("  zero-carry splits of M:", zc, " (superseded model had [12,24])")
    print("  generator classes:", sorted(set(cls(w, q_) for w in generator_uweights(3))))
    print("  state classes    :", [cls(w, q_) for w in state_uweights(3)])
    print("  spare classes    :", [cls(w, q_) for w in spare_uweights(3)])
    print()
    print("  the family window-denominator law  q_window(a) = 12a-7 = M_a:")
    for aa in range(2, 8):
        print("    a=%d  M=%2d  ordPhi=%4d  W_step=%s  q_window=%2d (prime %s)"
              % (aa, family(aa)["M"], family(aa)["ordPhi"],
                 Rational(family(aa)["ordPhi"], family(aa)["M"]), q_window(aa),
                 isprime(q_window(aa))))
    print()
    print("  (72,108) integral limit, untouched:")
    l72 = window_law(204, 17, 238)
    print("    W_step=%s alpha=%d q=%d deg_slope=%s affine=%s lambda=%s"
          % (l72["W_step"], l72["alpha"], l72["q"], l72["deg_slope"],
             l72["deg_affine"], l72["lam"]))
    print("    -> ord >= 12w, deg <= 14w, lambda = 2: TWO distinct slopes and a")
    print("       2-unit strip.  That is what (72,108) has and (75,125) lacks.")
    print("=" * 74)
    print("HEADLINE: the 'period-12 window functions' are REFUTED.  q_window is")
    print("29 (prime), there is NO affine degree cap, and the two-slope window")
    print("cone COLLAPSES to a single ray because Phi is a monomial.")
    print("=" * 74)
