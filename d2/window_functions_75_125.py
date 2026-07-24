#!/usr/bin/env python3
"""window_functions_75_125.py  (NEW; read-only over all existing artifacts)

THE PERIOD-12 WINDOW FUNCTIONS for (75,125) -- the arithmetic layer any
(75,125) window compiler needs, derived exactly (pure sympy).

Context.  G_SYSTEM_75_125.md builds the (75,125) D-transform G-system and
characterises its physical-weight boundary: W_step = ord_y(Phi)/M = 201/36 =
67/12 is NON-integral, so the affine window caps of (72,108) (ord >= 12k,
deg <= 15k/14k) become QUASIPOLYNOMIAL with quasi-period q_window = 5a-3 = 12.
F2_TOWER.md locates the y-order lattice that carries this: period 7 at a=2,
period 12 at a=3, incommensurate (gcd=1).  This module DERIVES, exactly:

  (1) The floor/ceiling window-cap functions.  For a window object at intrinsic
      u-slice weight w:
          lower cap (y-order):   L(w) = floor((alpha*w + beta_m)/q)  = ceil(alpha*w/q)
          upper cap (y-degree):  U(w) = deg_slope * w                (affine)
      with  m = w mod q,  and for the F2 family at parameter a:
          alpha     = 10 a^2 - 8 a + 1     (a=3: 67,  a=2: 25)
          q         = 5 a - 3   (= q_window) (a=3: 12,  a=2: 7)
          deg_slope = 5 a - 1              (a=3: 14,  a=2: 9)
          beta_m    = (-alpha * m) mod q   (a=3: [0,5,10,3,8,1,6,11,4,9,2,7])
      Only the LOWER cap is quasipolynomial (W_step non-integral); the UPPER cap
      is affine (deg_slope integral for all a).  This is the exact content of
      "quasipolynomial window cap".

  (2) The class-interaction (composition) table.  Window objects at u-weight w
      carry residue class m = w mod q; multiplication adds u-weights, so classes
      compose additively mod q:  class(x*y) = (class x + class y) mod q.  The
      beta offsets are a group 1-cocycle for this: beta_{m1} + beta_{m2} =
      beta_{(m1+m2) mod q} + q * carry, carry in {0,1}, and the carry is exactly
      the +1 defect of ceil (superadditivity of L).  The G-generators occupy the
      classes {0,2,3,4,5,6,7,8,9,10} (all but {1,11}); Phi sits in class 0.

  (3) Consistency.  (a) the caps admit the Phi point (504,201,101,202) exactly at
      equality: L(M)=ord_y(Phi)=201, U(M)=deg_y(Phi)=504, because M=36=3q is a
      multiple of q so the floor is exact.  (b) under the family substitution the
      caps reduce to the a=2 period-7 analogue (alpha=25,q=7,deg_slope=9), the
      control against f2_tower.py's window table; and to the (72,108) integral
      limit (q=1, ord>=12w, deg<=14w) matching WINDOW_CAPS.md.  (c) q_window=5a-3
      with gcd(alpha,q)=1.

THE HONEST BOUNDARY.  Everything above is forced by the u-graded G-system, the
built Phi signature, and ONE structural premise inherited verbatim from
(72,108): *Phi realises the extreme (minimal ord/weight) ray of the window cone*
-- exactly the sense in which (72,108)'s Phi "sits at its caps" (ord 204 = 12*17,
deg 238 = 14*17).  Under that premise the lower cap is ord >= (alpha/q) w, hence
(integrality) ord >= ceil(alpha*w/q) = L(w), and the twelve beta_m are pinned as
beta_m = (-alpha*m) mod q -- fully determined.  The one thing NOT fixable from
the u-grading alone is whether the actual (75,125) window cone dips *below* the
67/12 ray at some non-Phi weight (which would raise some beta_m): that would need
the deeper Newton polygon of P -- the "unreduced polygon" data (C_SERIES_75_125.md
judgment 2) that only the actual bridge construction supplies.  So: SLOPES,
PERIOD, deg-cap, the class-0 (Phi) line, and the canonical extreme-ray beta_m are
DERIVED; the residual is the named, characterised boundary.  See §7 of the .md.
"""
from sympy import Rational, gcd, ceiling, floor, denom


# ---------------------------------------------------------------------------
# Family constants (the F2 family: a = j+2, b = 2a-1, t=5, kappa=3,
# C = y^2(y^3+1), Phi_a = -(1/(3a)) y^(30a^2-24a+3) (y^3+1)^(15a^2-12a+2)).
# Source: f2_family_verify.py (landed), C_SERIES_75_125.md.
# ---------------------------------------------------------------------------
def family(a):
    """Return the exact (b,t,kappa,jphi,M,ordPhi,degPhi) for F2-family rung a."""
    b = 2 * a - 1
    t = 5
    kappa = 3
    jphi = 5 * a - 4                     # = a t - kappa - 1
    M = b * t + jphi                     # = 15a - 9 = 3(5a-3)
    ordPhi = 30 * a ** 2 - 24 * a + 3    # ord_y(Phi_a)
    degPhi = 3 * (5 * a - 1) * (5 * a - 3)  # deg_y(Phi_a) = 75a^2-60a+9
    return dict(a=a, b=b, t=t, kappa=kappa, jphi=jphi, M=M,
                ordPhi=ordPhi, degPhi=degPhi)


# ---------------------------------------------------------------------------
# The window-cap arithmetic, derived from (ordPhi, M, degPhi) for ANY case.
#   W_step = ord_y(Phi)/M  in lowest terms  =  alpha / q
#   deg_slope = deg_y(Phi)/M                (must be integral: affine deg cap)
#   L(w) = floor((alpha w + beta_m)/q) = ceil(alpha w / q)   (lower y-order cap)
#   U(w) = deg_slope * w                                     (upper y-degree cap)
#   beta_m = (-alpha m) mod q,   m = w mod q
# ---------------------------------------------------------------------------
def window_law(ordPhi, M, degPhi):
    """Derive the window-cap law from a case's Phi signature at slice M."""
    W = Rational(ordPhi, M)              # W_step in lowest terms
    alpha, q = W.p, W.q                  # numerator / quasi-period
    degR = Rational(degPhi, M)
    assert degR.q == 1, "deg cap not affine (deg_slope non-integral): %s" % degR
    deg_slope = int(degR)
    beta = [(-alpha * m) % q for m in range(q)]
    return dict(W_step=W, alpha=alpha, q=q, deg_slope=deg_slope, beta=beta)


def L(w, alpha, q, beta):
    """Lower y-order cap (quasipolynomial floor form)  ==  ceil(alpha w / q)."""
    return (alpha * w + beta[w % q]) // q


def L_ceil(w, alpha, q):
    """Same lower cap, ceiling form (the tight integer bound ord >= alpha w/q)."""
    return int(ceiling(Rational(alpha * w, q)))


def U(w, deg_slope):
    """Upper y-degree cap (affine)."""
    return deg_slope * w


def q_window(a):
    """The window-denominator law: q_window = denom(W_step) = 5a-3."""
    f = family(a)
    return denom(Rational(f["ordPhi"], f["M"]))


# ---------------------------------------------------------------------------
# The (75,125) instance (a=3) -- the target.
# ---------------------------------------------------------------------------
def target_75_125():
    f = family(3)
    law = window_law(f["ordPhi"], f["M"], f["degPhi"])
    f.update(law)
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
    print("PERIOD-12 WINDOW FUNCTIONS for (75,125)  [F2 family, a=3]")
    print("=" * 74)
    print("  alpha = 10a^2-8a+1 =", r["alpha"],
          "   q = q_window = 5a-3 =", r["q"],
          "   deg_slope = 5a-1 =", r["deg_slope"])
    print("  W_step = ord_y(Phi)/M =", r["W_step"], "(non-integral -> ord cap quasipolynomial)")
    print("  beta_m (m=0..11)      =", r["beta"])
    print("  lower cap  L(w) = floor((%d w + beta_m)/%d) = ceil(%d w/%d)"
          % (r["alpha"], r["q"], r["alpha"], r["q"]))
    print("  upper cap  U(w) = %d w   (affine)" % r["deg_slope"])
    print()
    a, q, be, dsl = r["alpha"], r["q"], r["beta"], r["deg_slope"]
    print("  Phi at M=%d:  L(%d)=%d = ord_y(Phi),  U(%d)=%d = deg_y(Phi)  [equality]"
          % (r["M"], r["M"], L(r["M"], a, q, be), r["M"], U(r["M"], dsl)))
    print("  generator classes:", sorted(set(cls(w, q) for w in generator_uweights(3))),
          " (Phi in class 0; skipped G10 at w=35 -> class 11)")
    print("  spare classes     :", [cls(w, q) for w in spare_uweights(3)])
    print()
    print("  sample caps (w : L(w) .. U(w)):")
    for w in list(range(0, 13)) + [26, 36]:
        print("    w=%2d  ord>=%3d  deg<=%3d   [class %2d, beta=%2d]"
              % (w, L(w, a, q, be), U(w, dsl), cls(w, q), be[w % q]))
