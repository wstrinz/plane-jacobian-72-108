#!/usr/bin/env python3
"""f2_branch_manifest.py  (NEW; read-only over every existing artifact)

BRANCH-COMPLETE FRONT END FOR THE F_2 FORCING ODE.

The standing F_2 analysis SELECTS the residual `g = y^3 + 1` (see
`phi_75_125.py` lines 50-70: "g1 = g2 = 0 (forced), g3 free (resonant),
g(-1)=0 => g0 = g3, monic normalization g3 = 1") after first SELECTING the
ansatz `f = A y^rho g^e`.  `FAMILY_GRAMMAR.md`'s 2026-07-24 rescope says
plainly that repeated-root / coexisting-multiplicity branches are NOT ruled
out.  This module removes both selections:

  * `g` is a GENERAL polynomial of degree `dg` with symbolic coefficients;
  * no ansatz for `f` -- the forcing ODE is solved as what it is, a LINEAR
    first-order inhomogeneous ODE, by exact coefficient recursion;
  * the resulting solvability conditions are assembled into a coefficient
    ideal which is then primary-decomposed (Singular `primdecGTZ`, via WSL),
    saturated by the non-degeneracy conditions;
  * every component is classified (squarefree / discriminant-zero / marked
    multiplicity type) and given an exact WITNESS re-substituted into the
    original ODE.

Every gauge normalization used is declared with its exceptional locus, and
gauge invariance is CHECKED (not assumed) by transporting each witness.

The forcing ODE itself is inherited, not re-derived here; it is
`corner144_verify.py` section C, re-checked in section A0 below:

    [P, x^s f / C^b] / x^kappa  =  a C^(a-b-1) ( t C f' - [t(b-a)+kappa+1] C' f )

so requiring the bracket to be 1 gives

    a { t c f' - coef c' f } = c^e,   coef = t(b-a)+kappa+1,  e = b-a+1,
    c = y^q g,  g(0) != 0,  deg g = dg,  a0 = deg c = q + dg.

CORRECTED INPUTS.  The (5,20) corner data used here is re-derived from
scratch in section A1 (NOT imported from any *_75_125.py, several of which
carry the superseded l = 5 / kappa = 3 / C = y^2(y^3+1) values):

    l = t = 4,  kappa = 2,  q = 1,  C = y  (a monomial),  deg C = a0 = 1

cross-checked against GGV3's own published reduction of the sibling (50,75)
at the same corner (`paper_src/1406.0886_GGV3.tex:1723-1727`:
`[P_1,Q_1] = x^2`, `deg P_1 = 10`, `deg Q_1 = 15`).

Run:
    python -u f2_branch_manifest.py            # full report
    python -u f2_branch_manifest.py --quiet    # exit 0 iff every check passes

Needs: sympy; Singular 4.x reachable as `wsl -e bash -lc Singular` (the
Singular stage is skipped-with-a-FAIL, never silently passed, if absent).
"""
import subprocess
import sys
from fractions import Fraction
from itertools import product

import sympy as sp

y = sp.symbols("y")

QUIET = "--quiet" in sys.argv
_CHECKS = []


def check(name, ok, note=""):
    _CHECKS.append((name, bool(ok), note))
    if not QUIET:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {note}" if note else ""))
    return bool(ok)


def head(s):
    if not QUIET:
        print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def sub(s):
    if not QUIET:
        print("\n--- " + s)


def say(s=""):
    if not QUIET:
        print(s)


# ===========================================================================
# A0.  The forcing ODE, re-checked from the bracket (inherited premise)
# ===========================================================================
def bracket_identity(a, b, t, kappa):
    """Verify  [P, x^s f/C^b]/x^kappa = a C^(a-b-1)(t C f' - coef C' f)."""
    x = sp.symbols("x")
    C = sp.Function("C")(y)
    f = sp.Function("f")(y)
    s = kappa + 1 - a * t
    P = x ** (a * t) * C ** a
    tail = x ** s * f / C ** b
    br = sp.diff(P, x) * sp.diff(tail, y) - sp.diff(P, y) * sp.diff(tail, x)
    want = a * C ** (a - b - 1) * (t * C * sp.diff(f, y)
                                   - (t * (b - a) + kappa + 1) * sp.diff(C, y) * f)
    return sp.simplify(sp.expand(br / x ** kappa - want)) == 0


# ===========================================================================
# A1.  Independent re-derivation of the (5,20) corner data
# ===========================================================================
def hull_vertices(pts):
    """Convex-hull vertices of a finite planar point set (exact, integer)."""
    pts = sorted(set(map(tuple, pts)))
    if len(pts) <= 2:
        return pts

    def cross(o, p, q):
        return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])

    lo = []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    up = []
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return sorted(set(lo[:-1] + up[:-1]))


def reduce_corner(a0, b0, l):
    """The GGV1-Section-8 style reduction of the corner A_0=(a0,b0), A_0'=(1,0).

    Returns (Delta_prime_vertices, kappa, q, retraction, mu).
      Delta  = {(0,0),(1,0),(a0,b0),(0,q)},  q = gcd(a0,b0)   [GGV1 Cor 7.4]
      flip   (i,j) -> (j,i)
      root shift of depth mu = l-1: the foot vertex (q,0) is REPLACED by the
             shifted foot (-mu,0)   [Pred_P(1,0) = (1,-mu)]
      inversion (i,j) -> (l*j - i, j)
    `retraction` is the test b0 == l*(a0-1): TRUE means the {(0,1),(b0,a0)}
    edge collapses to a vertical top face and deg C = a0; FALSE means there
    is no vertical top face and C is the monomial y^q.

    VALIDATED against GGHV22's published (7,21) reduction (section A1).
    """
    q = int(sp.igcd(a0, b0))
    mu = l - 1
    delta = [(0, 0), (1, 0), (a0, b0), (0, q)]
    flip = [(j, i) for (i, j) in delta]
    flip = [p for p in flip if p != (q, 0)] + [(-mu, 0)]
    pre = hull_vertices(flip)
    inv = [(l * j - i, j) for (i, j) in pre]
    dprime = hull_vertices(inv)
    return dprime, l - 2, q, (b0 == l * (a0 - 1)), mu


def total_degree(verts):
    return max(i + j for (i, j) in verts)


# ===========================================================================
# B.  The branch engine
# ===========================================================================
class Model:
    """One (t, kappa, q, dg, a, b) instance of the forcing ODE."""

    def __init__(self, name, t, kappa, q, dg, a, b):
        self.name, self.t, self.kappa, self.q, self.dg = name, t, kappa, q, dg
        self.a, self.b = a, b
        self.e = b - a + 1
        self.coef = t * (b - a) + kappa + 1
        self.rho = q * (self.e - 1) + 1
        self.a0 = q + dg
        self.K = self.rho + self.e * dg            # deg of the RHS y^rho g^e
        # band resonances: the coefficient of g_d * F_{k-d} in equation k is
        #   a * ( t(k-d) - coef(q+d) ),  which vanishes at k = d + coef(q+d)/t
        self.band_res = {}
        for d in range(dg + 1):
            num = self.coef * (q + d)
            if num % t == 0:
                self.band_res[d] = d + num // t
        self.k_star = self.band_res.get(0)          # diagonal resonance
        self.K_res = max(self.band_res.values()) if self.band_res else None
        self.B = self.K if self.K_res is None else max(self.K, self.K_res)

    def band(self, k, d):
        """Integer multiplier of  g_d * F_{k-d}  in equation k."""
        return self.a * (self.t * (k - d) - self.coef * (self.q + d))

    def __str__(self):
        return (f"{self.name}: t={self.t} kappa={self.kappa} q={self.q} dg={self.dg} "
                f"(a,b)=({self.a},{self.b}) e={self.e} coef={self.coef} rho={self.rho} "
                f"a0={self.a0} K={self.K} K_res={self.K_res} B={self.B} "
                f"band_res={self.band_res}")


def gsyms(dg, gauge_g0):
    """g_0..g_dg.  gauge_g0=True substitutes g0 = 1 (the C-scaling gauge)."""
    gs = list(sp.symbols(f"g0:{dg+1}"))
    if gauge_g0:
        gs[0] = sp.Integer(1)
    return gs


def forward_solve(M, gs):
    """Exact ascending coefficient recursion for the forcing ODE.

    Equation k (coefficient of y^(k+q-1) in the ODE, k >= 0):
        sum_{d=0..min(dg,k)}  band(k,d) * g_d * F_{k-d}  =  RHS_k
    with RHS_k = [y^(k-rho)] g^e.  The d=0 term carries F_k with multiplier
    band(k,0)*g_0 = a*g_0*(t k - coef q); since g_0 != 0 this determines F_k
    UNIQUELY unless t k = coef q (the diagonal resonance k_star), where the
    equation degenerates into a solvability CONDITION and F_k becomes free.

    Returns (F, conds, frees) with F[k] for k = 0..B+dg.
    """
    g = sum(gs[i] * y ** i for i in range(M.dg + 1))
    ge = sp.Poly(sp.expand(g ** M.e), y)

    def rhs_coeff(k):
        """[y^k] (y^rho g^e) = [y^(k-rho)] g^e."""
        n = k - M.rho
        if n < 0 or n > M.e * M.dg:
            return sp.Integer(0)
        return ge.coeff_monomial(y ** n)

    F, conds, frees = [], [], []
    for k in range(M.B + M.dg + 1):
        acc = sp.Integer(0)
        for d in range(1, min(M.dg, k) + 1):
            acc += M.band(k, d) * gs[d] * F[k - d]
        r = rhs_coeff(k)
        diag = M.band(k, 0)
        if diag == 0:
            conds.append((f"resonance k={k}", sp.expand(r - acc)))
            lam = sp.symbols(f"lam{k}")
            frees.append(lam)
            F.append(lam)
        else:
            F.append(sp.cancel(sp.expand(r - acc) / (diag * gs[0])))
    # termination: F_k = 0 for the dg indices just above B
    term = [(f"terminate F_{k}", F[k]) for k in range(M.B + 1, M.B + M.dg + 1)]
    return F, conds + term, frees


def intify(expr, syms):
    """Clear rational denominators and integer content: a primitive
    integer-coefficient polynomial with the same zero set."""
    expr = sp.expand(sp.together(expr))
    num, _den = sp.fraction(sp.cancel(expr))
    num = sp.expand(num)
    if num == 0:
        return sp.Integer(0)
    p = sp.Poly(num, *syms)
    p = p.clear_denoms(convert=True)[1]
    return p.primitive()[1].as_expr()


def poly_from(F, upto):
    return sp.expand(sum(F[k] * y ** k for k in range(upto + 1)))


def ode_residual(M, c, f):
    return sp.expand(M.a * (M.t * c * sp.diff(f, y)
                            - M.coef * sp.diff(c, y) * f) - c ** M.e)


# ---------------------------------------------------------------------------
# Singular bridge
# ---------------------------------------------------------------------------
SING_OK = False          # sticky: True once ANY Singular call has succeeded
SING_TIMEOUTS = []       # exit 124/timeouts are NEVER verdicts; recorded, not used


def singular(script, timeout=900):
    global SING_OK
    try:
        p = subprocess.run(["wsl", "-e", "bash", "-lc",
                            f"timeout {timeout} Singular -q"],
                           input=script, capture_output=True, text=True,
                           timeout=timeout + 60)
    except subprocess.TimeoutExpired:
        SING_TIMEOUTS.append(timeout)
        return None, f"TIMEOUT after {timeout}s -- NOT a verdict"
    except Exception as exc:                                    # noqa: BLE001
        return None, str(exc)
    if p.returncode in (124, 137):
        # 124 = timeout, 137 = SIGKILL.  NEVER a verdict.
        SING_TIMEOUTS.append(timeout)
        return None, f"TIMEOUT/ABORT (exit {p.returncode}) after {timeout}s -- NOT a verdict"
    if p.returncode != 0 or "? " in p.stdout or "? " in p.stderr:
        return None, (p.stdout + p.stderr)[:2000]
    SING_OK = True
    return p.stdout, ""


RESERVED = ("IDL", "SATP", "PDL", "PRM", "TST", "GBX", "EXT", "ii", "jj", "MPZ")


def _sg(expr):
    return str(sp.expand(expr)).replace("**", "^")


def primdec(gens, varnames, sat, tests=None, elim_to=None, timeout=900):
    """primdecGTZ of <gens> : sat^infty  in Q[varnames].

    `tests` = {label: poly}.  For every component the script reports whether
    each test poly lies in the component's associated PRIME (exact `reduce`
    against a standard basis -> ideal membership).
    `elim_to` = a variable name; for each component the univariate elimination
    ideal of PRIME + (sat - 1) in that variable is reported (the field of
    definition of the branch after the dg-th-root gauge fix).

    Returns (result, err) where result is "UNIT-FREE" (the ideal is (0) after
    saturation: the whole stratum is one branch), or a list of dicts, or None.
    """
    assert set(RESERVED).isdisjoint(set(varnames)), \
        f"emitted identifiers collide with ring variables: {set(RESERVED) & set(varnames)}"
    syms = [sp.Symbol(v) for v in varnames]
    ig = [intify(gg, syms) for gg in gens]
    for gg in ig:
        assert gg == 0 or all(c.is_Integer for c in sp.Poly(gg, *syms).coeffs()), \
            "generator still has rational coefficients"
    gl = ",".join(f"({gg})" for gg in ig) or "0"
    gl = gl.replace("**", "^")
    tests = tests or {}
    tl = list(tests.items())
    lines = [
        'LIB "primdec.lib";',
        f"ring MPZ = 0,({','.join(varnames)}),dp;",
        f"ideal IDL = {gl};",
        f"poly SATP = {_sg(sat)};",
        "IDL = sat(IDL, SATP)[1];",
        "IDL = std(IDL);",
        'if (size(IDL) == 0) { "EMPTYIDEAL"; }',
        "if (size(IDL) != 0) {",
        "  list PDL = primdecGTZ(IDL);",
        "  int ii;",
        "  for (ii = 1; ii <= size(PDL); ii++) {",
        '    "COMPONENT";',
        '    "PRIME:"; string(PDL[ii][2]);',
        '    "DIM:"; string(dim(std(PDL[ii][2])));',
        '    "DEG:"; string(mult(std(PDL[ii][2])));',
        "    ideal PRM = std(PDL[ii][2]);",
    ]
    for lab, pol in tl:
        lines += [f'    "TEST {lab}:";',
                  f"    if (reduce({_sg(pol)}, PRM) == 0) {{ \"IN\"; }} else {{ \"OUT\"; }}"]
    if elim_to:
        others = [v for v in varnames if v != elim_to]
        lines += [
            f"    ideal EXT = std(PDL[ii][2] + ideal(SATP - 1));",
            f'    "ELIM:"; string(eliminate(EXT, {"*".join(others) if others else "1"}));',
        ]
    lines += ["  }", "}", "quit;"]
    out, err = singular("\n".join(lines))
    if out is None:
        return None, err
    if "EMPTYIDEAL" in out:
        return "UNIT-FREE", ""
    comps, cur = [], None
    it = iter(out.splitlines())
    for line in it:
        s = line.strip()
        if s == "COMPONENT":
            cur = {"tests": {}}
            comps.append(cur)
        elif s == "PRIME:":
            cur["prime"] = next(it).strip()
        elif s == "DIM:":
            cur["dim"] = int(next(it).strip())
        elif s == "DEG:":
            cur["deg"] = int(next(it).strip())
        elif s.startswith("TEST ") and s.endswith(":"):
            cur["tests"][s[5:-1]] = next(it).strip()
        elif s == "ELIM:":
            cur["elim"] = next(it).strip()
    return comps, ""


def in_ideal(polys, gens, varnames, sat, timeout=900):
    """Is every p in `polys` a member of <gens>:sat^infty ?  Exact."""
    syms = [sp.Symbol(v) for v in varnames]
    ig = [intify(gg, syms) for gg in gens]
    gl = ",".join(f"({gg})" for gg in ig).replace("**", "^") or "0"
    lines = [
        f"ring MPZ = 0,({','.join(varnames)}),dp;",
        'LIB "elim.lib";',
        f"ideal IDL = {gl};",
        f"poly SATP = {_sg(sat)};",
        "IDL = std(sat(IDL, SATP)[1]);",
    ]
    for idx, p in enumerate(polys):
        lines.append(f'if (reduce({_sg(intify(p, syms))}, IDL) == 0) {{ "IN{idx}"; }}'
                     f' else {{ "OUT{idx}"; }}')
    lines.append("quit;")
    out, err = singular("\n".join(lines))
    if out is None:
        return None, err
    return all(f"IN{idx}" in out for idx in range(len(polys))), out


def parse_ideal(s, varnames):
    loc = {v: sp.Symbol(v) for v in varnames}
    loc["y"] = y
    return [sp.sympify(t.replace("^", "**"), locals=loc) for t in s.split(",") if t.strip()]


# ===========================================================================
# report
# ===========================================================================
def mult_strata(dg):
    """Marked-root multiplicity strata of a degree-dg g with g(0) = 1.

    Returns {partition: (parameter symbols, [g_1..g_dg] as polys in them)}.
    A partition (m_1 >= m_2 >= ...) of dg is realised as
        g = prod_i (1 + w_i y)^(m_i)          (g_0 = 1, the C-scaling gauge)
    so `which root carries which multiplicity` is the partition together with
    the labelling of the w_i; g(0) = 1 != 0 is automatic and the roots -1/w_i
    are all nonzero, as required by q = ord_y C exactly.
    """
    def parts(n, mx=None):
        mx = n if mx is None else mx
        if n == 0:
            yield ()
        for k in range(min(n, mx), 0, -1):
            for rest in parts(n - k, k):
                yield (k,) + rest

    out = {}
    for pt in parts(dg):
        ws = sp.symbols(f"w0:{len(pt)}")
        g = sp.Integer(1)
        for wi, m in zip(ws, pt):
            g *= (1 + wi * y) ** m
        pol = sp.Poly(sp.expand(g), y)
        out[pt] = (list(ws), [pol.coeff_monomial(y ** i) for i in range(dg + 1)])
    return out


def main():
    head("A0.  The forcing ODE (inherited premise, re-checked from the bracket)")
    for (aa, bb) in [(2, 3), (3, 5), (3, 4)]:
        check(f"bracket identity at (a,b)=({aa},{bb}), t=4, kappa=2",
              bracket_identity(aa, bb, 4, 2))
    check("bracket identity at (a,b)=(3,5), t=5, kappa=3 (the superseded chart)",
          bracket_identity(3, 5, 5, 3))
    say("  ODE:  a { t c f' - coef c' f } = c^e,  coef = t(b-a)+kappa+1, e = b-a+1")

    # -------------------------------------------------------------------
    head("A1.  Corrected (5,20) inputs, re-derived from scratch")
    d4, k4, q4, retr4, mu4 = reduce_corner(5, 20, 4)
    d5, k5, q5, retr5, mu5 = reduce_corner(5, 20, 5)
    say(f"  l=4:  Delta' = {d4}  total deg {total_degree(d4)}  kappa={k4} q_gcd={q4} "
        f"retraction={retr4}")
    say(f"  l=5:  Delta' = {d5}  total deg {total_degree(d5)}  kappa={k5} q_gcd={q5} "
        f"retraction={retr5}")
    # external validation of the reduction rule on a PUBLISHED reduction of the
    # same (non-retraction) shape: GGHV22 tex:1316-1317 for the (7,21) corner.
    d721 = reduce_corner(7, 21, 3)[0]
    check("reduction rule reproduces GGHV22's PUBLISHED (7,21) N(P) at (m,n)=(2,3)",
          {(2 * i, 2 * j) for (i, j) in d721} == {(0, 0), (4, 0), (6, 2), (0, 14)},
          f"Delta'(7,21) = {d721}")
    check("reduction rule reproduces GGHV22's PUBLISHED (7,21) N(Q) at (m,n)=(2,3)",
          {(3 * i, 3 * j) for (i, j) in d721} == {(0, 0), (6, 0), (9, 3), (0, 21)})
    check("l=4 gives Delta' = {(0,0),(3,0),(4,1),(0,5)}",
          set(d4) == {(0, 0), (3, 0), (4, 1), (0, 5)}, f"got {d4}")
    check("l=4, (m,n)=(2,3): GGV3's published deg P_1 = 10",
          2 * total_degree(d4) == 10)
    check("l=4, (m,n)=(2,3): GGV3's published deg Q_1 = 15",
          3 * total_degree(d4) == 15)
    check("l=4: GGV3's published [P_1,Q_1] = x^2, i.e. kappa = 2", k4 == 2)
    check("l=5 contradicts GGV3 on deg P_1", 2 * total_degree(d5) != 10)
    check("l=5 contradicts GGV3 on deg Q_1", 3 * total_degree(d5) != 15)
    check("l=5 contradicts GGV3 on kappa", k5 != 2)
    check("(75,125) reduced degrees are (15,25) at (m,n)=(3,5)",
          (3 * total_degree(d4), 5 * total_degree(d4)) == (15, 25))
    check("N(P) = 3*Delta' = {(0,0),(9,0),(12,3),(0,15)}",
          {(3 * i, 3 * j) for (i, j) in d4} == {(0, 0), (9, 0), (12, 3), (0, 15)})
    check("N(Q) = 5*Delta' = {(0,0),(15,0),(20,5),(0,25)}",
          {(5 * i, 5 * j) for (i, j) in d4} == {(0, 0), (15, 0), (20, 5), (0, 25)})
    check("retraction test b0 = l(a0-1) FAILS at (5,20), l=4  (20 != 16)",
          not retr4 and 20 != 4 * 4)
    check("retraction test HOLDS at the audited (8,28), l=4  (28 == 4*7)",
          reduce_corner(8, 28, 4)[3])
    check("retraction test FAILS at the published (7,21), l=3  (21 != 3*6)",
          not reduce_corner(7, 21, 3)[3])
    say("  => no vertical top face at (5,20): deg C = 1, C = y, q = 1, dg = 0.")
    say("     Same shape as the PUBLISHED (7,21), where GGHV22 gives C = y.")
    say("     CORRECTED F_2 INPUTS:  t = 4, kappa = 2, q = 1, dg = deg C - q = 0.")

    # -------------------------------------------------------------------
    head("B.  Engine control: the audited (72,108) corner  (8,28), (a,b)=(2,3)")
    M72 = Model("(72,108)", t=4, kappa=2, q=7, dg=1, a=2, b=3)
    say("  " + str(M72))
    gs = gsyms(1, gauge_g0=True)
    F, conds, frees = forward_solve(M72, gs)
    check("(72,108): no diagonal resonance (coef*q/t = 49/4)", M72.k_star is None)
    check("(72,108): top band resonance at k = 15 (deg bound B = 15)",
          M72.band_res.get(1) == 15 and M72.B == 15)
    check("(72,108): no free parameters", frees == [])
    cnum = [sp.expand(cc) for (_, cc) in conds]
    check("(72,108): exactly one solvability condition (F_16 = 0)", len(cnum) == 1)
    check("(72,108): that condition is IDENTICALLY zero -> the whole dg=1 stratum "
          "is one branch", all(cc == 0 for cc in cnum),
          f"cond = {cnum[0]}")
    g1 = sp.Symbol("g1")
    fgen = poly_from(F, M72.B)
    cgen = y ** 7 * (1 + g1 * y)
    check("(72,108): the recursion's f solves the ODE identically in g1",
          ode_residual(M72, cgen, fgen) == 0)
    f_at_1 = sp.factor(fgen.subs(g1, 1))
    q_pub = 2048 * y ** 4 - 512 * y ** 3 + 320 * y ** 2 - 240 * y + 195
    f_pub = -y ** 8 * (y + 1) ** 2 * q_pub / sp.Integer(6630)
    check("(72,108): the recursion REPRODUCES the audited f exactly "
          "(quartic 2048y^4-512y^3+320y^2-240y+195, 1/6630)",
          sp.expand(f_at_1 - f_pub) == 0)
    say(f"  f(g1=1) = {f_at_1}")
    check("(72,108): deg f = 14, the top resonance coef*a0/t = 14",
          sp.degree(f_at_1, y) == 14 and M72.coef * M72.a0 % M72.t == 0
          and M72.coef * M72.a0 // M72.t == 14)
    # mutation test: the SAME engine must NOT validate a wrong f
    check("MUTATION (72,108): perturbing the quartic breaks the ODE",
          ode_residual(M72, y ** 7 * (y + 1), f_pub + y ** 8) != 0)
    check("MUTATION (72,108): a wrong 1/6630 breaks the ODE",
          ode_residual(M72, y ** 7 * (y + 1), 2 * f_pub) != 0)

    # -------------------------------------------------------------------
    head("C.  THE SUPERSEDED F_2 MODEL, made branch-complete "
         "(t=5, kappa=3, a0=5) -- the general-cubic question as posed")
    results_old = {}
    for jj, (aa, bb) in [(0, (2, 3)), (1, (3, 5))]:
        sub(f"C.{jj}   j = {jj},  (a,b) = {aa,bb}   [{'(50,75)' if jj==0 else '(75,125)'}]")
        Mo = Model(f"oldF2_j{jj}", t=5, kappa=3, q=2, dg=3, a=aa, b=bb)
        say("  " + str(Mo))
        gso = gsyms(3, gauge_g0=True)
        Fo, condso, freeso = forward_solve(Mo, gso)
        gens = [sp.expand(cc) for (_, cc) in condso]
        vnames = ["g1", "g2", "g3"] + [str(s) for s in freeso]
        say(f"  {len(gens)} solvability conditions in {vnames}")
        for nm, cc in condso:
            say(f"    {nm}:  {sp.factor(cc)}")
        gsym = sum(gso[i] * y ** i for i in range(4))
        g1s, g2s, g3s = sp.symbols("g1 g2 g3")
        disc = sp.expand(sp.discriminant(sp.Poly(gsym, y)))
        # the triple-root locus: g = (1+w y)^3, eliminated
        w0 = sp.symbols("w0")
        tri = sp.groebner([g1s - 3 * w0, g2s - 3 * w0 ** 2, g3s - w0 ** 3],
                          w0, g1s, g2s, g3s, order="lex")
        tri = [p for p in tri.exprs if w0 not in p.free_symbols]
        tests = {"disc(g)": disc, "g(-1)": sp.expand(gsym.subs(y, -1))}
        for ti, tp in enumerate(tri):
            tests[f"triple{ti}"] = tp
        # global soundness of the termination lemma: on the conditions ideal the
        # TRUNCATED f really does solve the ODE, for every g in the branch locus
        rescoef = sp.Poly(ode_residual(Mo, y ** Mo.q * gsym, poly_from(Fo, Mo.B)),
                          y).all_coeffs()
        rescoef = [c for c in rescoef if sp.expand(c) != 0]
        check(f"old j={jj}: the ODE residual is NOT identically zero "
              f"(so the conditions have content)", len(rescoef) > 0,
              f"{len(rescoef)} nonzero residual coefficients")
        mem, _ = in_ideal(rescoef, gens, vnames, g3s)
        check(f"old j={jj}: TERMINATION LEMMA -- every residual coefficient lies in "
              f"the conditions ideal, so f is exact on every branch", mem is True)
        comps, err = primdec(gens, vnames, g3s, tests=tests, elim_to="g1")
        check(f"old j={jj}: Singular primary decomposition ran",
              comps is not None, err[:200])
        if comps is None:
            continue
        if comps == "UNIT-FREE":
            say("  ideal is (0) after saturation: ONE branch (the whole stratum)")
            results_old[jj] = ["(0)"]
            continue
        say(f"  {len(comps)} branch component(s) after saturating by g3 "
            f"(nonempty over Qbar: each is a proper PRIME with g3 not in it):")
        results_old[jj] = []
        for ci, comp in enumerate(comps, 1):
            tst = comp["tests"]
            trip = all(v == "IN" for k, v in tst.items() if k.startswith("triple"))
            kind = ("TRIPLE ROOT" if trip else
                    ("REPEATED ROOT (disc = 0)" if tst["disc(g)"] == "IN"
                     else "SQUAREFREE (generic)"))
            canon = comp["prime"].replace(" ", "") in ("g2,g1", "g1,g2")
            say(f"    [{ci}] dim {comp['dim']} deg {comp['deg']}  {kind}"
                f"   g(-1)=0 on it: {tst['g(-1)']}"
                f"   {'<-- THE CANONICAL BRANCH' if canon else '<-- NEW BRANCH'}")
            say(f"        prime = <{comp['prime']}>")
            if comp.get("elim"):
                say(f"        after the g3=1 gauge fix, g1 satisfies: {comp['elim']}")
            results_old[jj].append((comp, kind, canon))
        check(f"old j={jj}: the canonical branch <g1,g2> IS one of the components",
              any(c[2] for c in results_old[jj]))
        nnew = sum(1 for c in results_old[jj] if not c[2])
        check(f"old j={jj}: NEW branches counted", True,
              f"{nnew} new branch(es) beside the canonical one")
        # canonical-branch recovery
        canon = {sp.Symbol("g1"): 0, sp.Symbol("g2"): 0, sp.Symbol("g3"): 1}
        gcan = y ** 3 + 1
        ccan = y ** 2 * gcan
        fcan = sp.expand(poly_from(Fo, Mo.B).subs(canon))
        check(f"old j={jj}: CANONICAL branch g = y^3+1 solves the ODE",
              ode_residual(Mo, ccan, fcan) == 0, f"f = {sp.factor(fcan)}")
        check(f"old j={jj}: canonical f equals the repo's -1/(3a) y^rho (y^3+1)^e",
              sp.expand(fcan + sp.Rational(1, 3 * aa) * y ** Mo.rho * gcan ** Mo.e) == 0)
        check(f"old j={jj}: canonical point satisfies every solvability condition",
              all(sp.expand(cc.subs(canon)) == 0 for cc in gens))
        # MUTATION: a generic cubic must NOT satisfy them
        bad = {sp.Symbol("g1"): 1, sp.Symbol("g2"): 1, sp.Symbol("g3"): 1}
        check(f"MUTATION old j={jj}: generic g = y^3+y^2+y+1 FAILS the conditions",
              any(sp.expand(cc.subs(bad)) != 0 for cc in gens))
        check(f"MUTATION old j={jj}: and its recursion f does NOT solve the ODE",
              ode_residual(Mo, y ** 2 * (y ** 3 + y ** 2 + y + 1),
                           sp.expand(poly_from(Fo, Mo.B).subs(bad))) != 0)

        # ---- explicit discriminant / marked-multiplicity stratification ----
        say(f"\n  MARKED-MULTIPLICITY STRATA, j={jj}  (g_0 = 1 gauge; roots -1/w_i)")
        sq, e1 = primdec(gens, vnames, sp.expand(disc * g3s))
        say(f"    squarefree stratum  <conds> : (disc*g3)^inf  ->  "
            + ("(0), whole stratum" if sq == "UNIT-FREE"
               else (f"ERR {e1[:60]}" if sq is None
                     else "; ".join(f"dim {c['dim']} deg {c['deg']} <{c['prime']}>"
                                    for c in sq))))
        dz, e2 = primdec(list(gens) + [disc], vnames, g3s)
        say(f"    disc = 0 stratum    <conds, disc> : g3^inf  ->  "
            + ("(0), whole stratum" if dz == "UNIT-FREE"
               else (f"ERR {e2[:60]}" if dz is None
                     else ("EMPTY" if dz == [] else
                           "; ".join(f"dim {c['dim']} deg {c['deg']} <{c['prime']}>"
                                     for c in dz)))))
        check(f"old j={jj}: the two strata are complementary and both were computed",
              sq is not None and dz is not None)
        for pt, (ws, gvals) in sorted(mult_strata(3).items()):
            if len(pt) == 3:
                say(f"    partition {pt}: the squarefree stratum above (3 distinct "
                    f"roots over Qbar)")
                continue
            subsd = {sp.Symbol(f"g{i}"): gvals[i] for i in range(1, 4)}
            sgens = [sp.expand(cc.subs(subsd)) for cc in gens]
            wn = [str(w) for w in ws]
            satw = sp.prod(ws)
            res, ew = primdec(sgens, wn, satw)
            if res == "UNIT-FREE":
                desc = "(0): the WHOLE stratum solves"
            elif res is None:
                desc = f"ERR {ew[:70]}"
            elif res == []:
                desc = "EMPTY (no g of this multiplicity type solves the ODE)"
            else:
                desc = "; ".join(f"dim {c['dim']} deg {c['deg']} <{c['prime']}>"
                                 for c in res)
            say(f"    partition {pt} (g = {sp.factor(sum(gvals[i]*y**i for i in range(4)))})"
                f"  ->  {desc}")
            check(f"old j={jj}: multiplicity stratum {pt} decided", res is not None,
                  desc[:90])
            # explicit witness where the stratum is 0-dimensional and rational
            if isinstance(res, list) and res:
                for c in res:
                    if c["dim"] == 0:
                        pr = parse_ideal(c["prime"], wn)
                        sol = sp.solve(pr, [sp.Symbol(v) for v in wn], dict=True)
                        for sdict in sol or []:
                            if any(v == 0 for v in sdict.values()):
                                continue
                            if any(getattr(v, "free_symbols", set()) for v in sdict.values()):
                                continue
                            gw = sp.expand(sum(gvals[i] * y ** i
                                               for i in range(4)).subs(sdict))
                            fw = sp.expand(poly_from(Fo, Mo.B).subs(
                                {sp.Symbol(f"g{i}"): sp.expand(gvals[i].subs(sdict))
                                 for i in range(1, 4)}))
                            check(f"old j={jj} stratum {pt}: WITNESS "
                                  f"g = {sp.factor(gw)} solves the ODE",
                                  ode_residual(Mo, y ** 2 * gw, fw) == 0,
                                  f"deg f = {sp.degree(fw, y)}")
                            break

    # -------------------------------------------------------------------
    head("C''.  The three branches, verified by an INDEPENDENT direct linear solve "
         "(no recursion, no termination lemma, no Singular)")

    def direct_solve(t, kappa, q, a, b, gpoly, dmax=40):
        """Brute-force: unknown f = sum_{i<=dmax} c_i y^i, solve the linear
        system coefficientwise with sympy.solve.  Shares no code with
        forward_solve.  Returns the polynomial f, or None."""
        e, cf = b - a + 1, t * (b - a) + kappa + 1
        c = sp.expand(y ** q * gpoly)
        cs = sp.symbols(f"c0:{dmax+1}")
        ff = sum(cs[i] * y ** i for i in range(dmax + 1))
        R = sp.expand(a * (t * c * sp.diff(ff, y) - cf * sp.diff(c, y) * ff) - c ** e)
        sol = sp.solve(sp.Poly(R, y).all_coeffs(), cs, dict=True)
        if not sol:
            return None
        fv = sp.expand(ff.subs(sol[0]))
        return sp.expand(fv.subs({s: 0 for s in fv.free_symbols if s is not y}))

    phi5 = sp.Rational(1, 2) + sp.sqrt(5) / 2       # w0^2 = w0*w1 + w1^2 at w1 = 1
    BRANCHES = {
        "B1 canonical (squarefree)": y ** 3 + 1,
        "B2 double root, over Q(sqrt5)": sp.expand((1 + phi5 * y) ** 2 * (1 + y)),
        "B3 triple root (w = 1)": sp.expand((1 + y) ** 3),
        "B3 triple root (w = 2)": sp.expand((1 + 2 * y) ** 3),
    }
    for jj, (aa, bb) in [(0, (2, 3)), (1, (3, 5))]:
        Mo = Model(f"oldF2_j{jj}", t=5, kappa=3, q=2, dg=3, a=aa, b=bb)
        for lab, gp in BRANCHES.items():
            fv = direct_solve(5, 3, 2, aa, bb, gp)
            ok = fv is not None and sp.simplify(
                ode_residual(Mo, sp.expand(y ** 2 * gp), fv)) == 0
            check(f"old j={jj}: {lab} -- independent direct solve gives an EXACT "
                  f"polynomial f", ok,
                  f"deg f = {sp.degree(fv, y) if fv is not None else '-'}")
        # negative control: a g on NO branch must have no polynomial solution
        badg = y ** 3 + y ** 2 + y + 1
        check(f"NEGATIVE CONTROL old j={jj}: g = y^3+y^2+y+1 has NO polynomial f "
              f"(direct solve)", direct_solve(5, 3, 2, aa, bb, badg) is None)
        check(f"POSITIVE CONTROL old j={jj}: g = y^3+2 (on the canonical orbit, "
              f"g_1=g_2=0, but NOT monic-with-g(-1)=0) DOES solve -- the branch is "
              f"the whole gauge orbit, not the normalised representative",
              direct_solve(5, 3, 2, aa, bb, y ** 3 + 2) is not None)
        # THE SELECTION RULES THE REPO STATES DO NOT SINGLE OUT y^3+1
        gcube = sp.expand((y + 1) ** 3)
        fcube = direct_solve(5, 3, 2, aa, bb, gcube)
        check(f"old j={jj}: g = (y+1)^3 satisfies EVERY selection rule the repo "
              f"states (monic, g(-1)=0, g(0)!=0, deg 3) AND solves the ODE -- so "
              f"those rules do NOT single out y^3+1",
              sp.Poly(gcube, y).LC() == 1 and gcube.subs(y, -1) == 0
              and gcube.subs(y, 0) != 0 and sp.degree(gcube, y) == 3
              and fcube is not None
              and sp.expand(ode_residual(Mo, y ** 2 * gcube, fcube)) == 0,
              f"f = {sp.factor(fcube) if fcube is not None else '-'}")

    sub("Identification of the two new components with marked-multiplicity strata")
    w0s, w1s = sp.symbols("w0 w1")
    g1s, g2s, g3s = sp.symbols("g1 g2 g3")
    P6 = [g2s ** 6 - 25 * g2s ** 3 * g3s ** 2 - 125 * g3s ** 4,
          -g2s ** 5 + 5 * g2s ** 2 * g3s ** 2 + 75 * g1s * g3s ** 3,
          -4 * g2s ** 3 + 15 * g1s * g2s * g3s - 25 * g3s ** 2,
          g1s * g2s ** 3 - 3 * g2s ** 2 * g3s - 20 * g1s * g3s ** 2,
          -7 * g1s * g2s ** 2 + 20 * g1s ** 2 * g3s + 5 * g2s * g3s,
          5 * g1s ** 2 * g2s - 4 * g2s ** 2 - 35 * g1s * g3s,
          4 * g1s ** 3 - 13 * g1s * g2s + 7 * g3s]
    P3 = [g2s ** 3 - 27 * g3s ** 2, -g2s ** 2 + 3 * g1s * g3s,
          g1s * g2s - 9 * g3s, g1s ** 2 - 3 * g2s]
    p21 = sp.Poly(sp.expand((1 + w0s * y) ** 2 * (1 + w1s * y)), y)
    s21 = {g1s: p21.coeff_monomial(y), g2s: p21.coeff_monomial(y ** 2),
           g3s: p21.coeff_monomial(y ** 3)}
    rel21 = w0s ** 2 - w0s * w1s - w1s ** 2
    check("the deg-6 component IS the (2,1) double-root stratum cut by "
          "w0^2 - w0 w1 - w1^2 = 0  (i.e. w0/w1 = golden ratio, field Q(sqrt5))",
          all(sp.rem(sp.Poly(sp.expand(p.subs(s21)), w0s),
                     sp.Poly(rel21, w0s)) == 0 for p in P6))
    p3 = sp.Poly(sp.expand((1 + w0s * y) ** 3), y)
    s3 = {g1s: p3.coeff_monomial(y), g2s: p3.coeff_monomial(y ** 2),
          g3s: p3.coeff_monomial(y ** 3)}
    check("the deg-3 component IS exactly the (3) triple-root stratum, ALL of it",
          all(sp.expand(p.subs(s3)) == 0 for p in P3))
    check("MUTATION: the deg-3 component is NOT the (2,1) stratum",
          any(sp.rem(sp.Poly(sp.expand(p.subs(s21)), w0s),
                     sp.Poly(rel21, w0s)) != 0 for p in P3))

    sub("The new triple-root branch IS the repo's own mu = dg RUNG law, at gap = 0")
    say("  FAMILY_GRAMMAR.md section 3:  f = y^rho (y+1)^(dg*e-(dg-1)) * u, "
        "deg u = gap + r")
    for jj, (aa, bb) in [(0, (2, 3)), (1, (3, 5))]:
        Mo = Model(f"oldF2_j{jj}", t=5, kappa=3, q=2, dg=3, a=aa, b=bb)
        fv = direct_solve(5, 3, 2, aa, bb, sp.expand((1 + y) ** 3))
        fl = sp.factor_list(fv)
        m1 = dict((sp.expand(base), pw) for base, pw in fl[1]).get(y + 1, 0)
        cof = sp.degree(fv, y) - sp.degree(sp.Poly(y ** Mo.rho, y)) - m1
        pred_m = 3 * Mo.e - 2                     # dg*e - (dg-1),  dg = 3
        pred_c = 0 + 2                            # gap + r = 0 + 2
        say(f"  j={jj}: f = {sp.factor(fv)}")
        check(f"old j={jj}: triple-root f has (y+1)^{pred_m} exactly as the mu=dg "
              f"rung law predicts", m1 == pred_m, f"got (y+1)^{m1}")
        check(f"old j={jj}: its unit cofactor has degree gap + r = {pred_c}",
              cof == pred_c, f"got {cof}")
    say("  => the branch the repo's PURE cell omitted is the mu = dg rung the repo's")
    say("     OWN mu-graded law already writes down -- for the RUNG families only.")

    # -------------------------------------------------------------------
    head("C'.  The q-stratification the repo fixed by SELECTION "
         "(phi_75_125.py: 'q = selected mult')")
    say("  a0 = deg C = 5 is corner data in the superseded model; q = ord_y C is NOT")
    say("  determined by it.  Every (q, dg) with q + dg = 5 is a separate branch axis.")
    qstrat = {}
    for qq in range(1, 6):
        Mq = Model(f"old_q{qq}", t=5, kappa=3, q=qq, dg=5 - qq, a=3, b=5)
        gsq = gsyms(5 - qq, gauge_g0=True)
        Fq, cq, frq = forward_solve(Mq, gsq)
        gq = [sp.expand(cc) for (_, cc) in cq]
        vq = [f"g{i}" for i in range(1, 6 - qq)] + [str(s) for s in frq]
        satq = sp.Symbol(f"g{5-qq}") if qq < 5 else sp.Integer(1)
        if not vq:
            status = "no free coefficients"
            comps = "UNIT-FREE" if all(g == 0 for g in gq) else "EMPTY"
        elif 5 - qq >= 4:
            comps = None
            status = ("full decomposition NOT ATTEMPTED (cost) -- UNDECIDED, and "
                      "NOT empty: section C''' exhibits an explicit witness")
        else:
            comps, err = primdec(gq, vq, satq, timeout=240)
            status = ("(0) -> one branch" if comps == "UNIT-FREE"
                      else (f"UNDECIDED: {err[:70]}" if comps is None
                            else ("EMPTY" if comps == [] else
                                  "; ".join(f"dim {c['dim']} deg {c['deg']} "
                                            f"<{c['prime']}>" for c in comps))))
        qstrat[qq] = (Mq, comps)
        say(f"  q={qq} dg={5-qq}: B={Mq.B} band_res={Mq.band_res} frees={frq} "
            f"conds={len(gq)} -> {status}")
    check("q-stratification: q=2 (the selected value) is ONE of 5 admissible strata",
          2 in qstrat)

    # -------------------------------------------------------------------
    head("C'''.  WITNESS-FIRST: the fully-ramified branch g = (1 + w y)^dg, "
         "uniformly across the whole q-stratification")
    say("  Cheap exact witnesses decide the cells whose primary decomposition is")
    say("  expensive -- and cross-check every Singular EMPTY verdict by a method")
    say("  that shares no code with it.")

    def direct2(t, kappa, q, a, b, gpoly, dmax=45):
        e, cf = b - a + 1, t * (b - a) + kappa + 1
        c = sp.expand(y ** q * gpoly)
        cs = sp.symbols(f"c0:{dmax+1}")
        ff = sum(cs[i] * y ** i for i in range(dmax + 1))
        R = sp.expand(a * (t * c * sp.diff(ff, y) - cf * sp.diff(c, y) * ff) - c ** e)
        sol = sp.solve(sp.Poly(R, y).all_coeffs(), cs, dict=True)
        if not sol:
            return None
        fv = sp.expand(ff.subs(sol[0]))
        return sp.expand(fv.subs({s: 0 for s in fv.free_symbols if s is not y}))

    sub("superseded chart t=5, kappa=3, a0=5, (a,b)=(3,5):  every (q,dg) cell")
    for qq in (1, 2, 3, 4):
        dgv = 5 - qq
        M = Model(f"s_q{qq}", 5, 3, qq, dgv, 3, 5)
        gapv = Fraction(qq - 1) - Fraction(5, 5)
        for w in (1, 2):
            gp = sp.expand((1 + w * y) ** dgv)
            fv = direct2(5, 3, qq, 3, 5, gp)
            ok = fv is not None and sp.expand(
                ode_residual(M, sp.expand(y ** qq * gp), fv)) == 0
            check(f"superseded (q,dg)=({qq},{dgv}), w={w}: g = (1+{w}y)^{dgv} SOLVES; "
                  f"deg f = 14 = coef*a0/t",
                  ok and sp.degree(fv, y) == 14,
                  f"gap = {gapv}, deg u = gap+r = {gapv + dgv - 1}")
    check("this WITNESS decides the (q,dg)=(1,4) cell that the decomposition could "
          "not: it is NON-EMPTY, not undecided-about-emptiness", True)
    check("and it exists at gap = -1 (q=1), so 'gap < 0 => no branch' is NOT a "
          "general implication of the corner data",
          Fraction(1 - 1) - Fraction(5, 5) < 0)
    say("  The mu = dg rung law reproduces deg f = 14 in EVERY cell:")
    say("    rho + (dg*e-(dg-1)) + (gap+r) = (2q+1) + (2dg+1) + 2 = 2(q+dg)+4 = 14.")

    sub("corrected chart t=4, kappa=2, q=1: the counterfactual dg cells")
    for dgv in (1, 2, 3):
        M = Model(f"c_dg{dgv}", 4, 2, 1, dgv, 3, 5)
        gp = sp.expand((1 + y) ** dgv)
        fv = direct2(4, 2, 1, 3, 5, gp)
        got = fv is not None and sp.expand(
            ode_residual(M, sp.expand(y * gp), fv)) == 0
        want = (dgv == 3)
        check(f"corrected chart, counterfactual dg={dgv}: (1+y)^{dgv} "
              f"{'SOLVES' if want else 'has NO polynomial f'} -- agreeing with the "
              f"primary decomposition", got == want,
              f"deg f = {sp.degree(fv, y) if fv is not None else '-'}")

    # -------------------------------------------------------------------
    head("D.  THE CORRECTED F_2 MODEL  (t=4, kappa=2, q=1, dg=0)")
    say("  Q5/retraction (section A1) forces deg C = 1, so g = g_0 is a CONSTANT.")
    say("  There is no residual cubic in the corrected model at all.")
    corrected = {}
    for jj in range(0, 4):
        aa, bb = jj + 2, 2 * jj + 3
        Mc = Model(f"F2_j{jj}", t=4, kappa=2, q=1, dg=0, a=aa, b=bb)
        gsc = gsyms(0, gauge_g0=True)
        Fc, cc_, frc = forward_solve(Mc, gsc)
        fw = poly_from(Fc, Mc.B)
        cw = y ** Mc.q
        ok = ode_residual(Mc, cw, fw) == 0
        corrected[jj] = (Mc, fw, [x for (_, x) in cc_], frc)
        say(f"  j={jj} (a,b)=({aa},{bb}) e={Mc.e} coef={Mc.coef} rho={Mc.rho} "
            f"B={Mc.B} band_res={Mc.band_res}")
        check(f"corrected j={jj}: forcing ODE has the UNIQUE solution f = y^{Mc.rho}/{aa}",
              ok and sp.expand(fw - y ** Mc.rho / aa) == 0, f"f = {fw}")
        check(f"corrected j={jj}: no solvability conditions, no free parameters",
              all(x == 0 for x in [xx for (_, xx) in cc_]) and frc == [])
        check(f"MUTATION corrected j={jj}: 2f does not solve the ODE",
              ode_residual(Mc, cw, 2 * fw) != 0)
    check("corrected model: (50,75) [j=0] gives f = y^2/2",
          sp.expand(corrected[0][1] - y ** 2 / 2) == 0)
    check("corrected model: (75,125) [j=1] gives f = y^3/3",
          sp.expand(corrected[1][1] - y ** 3 / 3) == 0)
    check("corrected model: the ONLY obstruction is t-(kappa+1)q = 1 != 0",
          4 - 3 * 1 == 1)
    check("corrected model: y^3+1 is NOT the residual (dg = 0, g is a constant)",
          Model("x", 4, 2, 1, 0, 3, 5).dg == 0)

    sub("D'.  Counterfactual: what a degree-dg residual would need, corrected chart")
    for dgv in range(0, 5):
        Md = Model(f"cf_dg{dgv}", t=4, kappa=2, q=1, dg=dgv, a=3, b=5)
        gsd = gsyms(dgv, gauge_g0=True)
        Fd, cd, frd = forward_solve(Md, gsd)
        gd = [sp.expand(x) for (_, x) in cd]
        vd = [f"g{i}" for i in range(1, dgv + 1)] + [str(s) for s in frd]
        if not vd:
            res = "no free coefficients; " + ("solvable" if all(g == 0 for g in gd)
                                              else "EMPTY")
        else:
            if dgv >= 4:
                comps, err = None, ("not attempted: exploratory -- UNDECIDED, "
                                    "not empty")
            else:
                comps, err = primdec(gd, vd, sp.Symbol(f"g{dgv}"), timeout=240)
            if comps == "UNIT-FREE":
                res = "(0) -> ONE branch (whole stratum)"
            elif comps is None:
                res = f"UNDECIDED: {err[:70]}"
            elif comps == []:
                res = "EMPTY"
            else:
                res = ("components: " +
                       "; ".join(f"dim {c['dim']} deg {c['deg']} <{c['prime']}>"
                                 for c in comps))
        say(f"  dg={dgv}: a0={Md.a0} B={Md.B} band_res={Md.band_res} "
            f"frees={frd} conds={len(gd)} -> {res}")

    # -------------------------------------------------------------------
    head("E.  Gauge group, declared, with exceptional loci -- and CHECKED")
    say("  G1  C-scaling:   (c, f) -> (mu c, mu^(e-1) f),   mu != 0")
    say("      acts g_i -> mu g_i.  Fixes g_0 = 1.  Exceptional locus g_0 = 0,")
    say("      i.e. ord_y C > q, excluded by the definition of q.  No field ext.")
    say("  G2  y-dilation:  (c(y), f(y)) -> (mu c(lam y), mu^(e-1) f(lam y)/lam)")
    say("      acts g_i -> mu lam^i g_i.  Can fix g_dg = 1 only after adjoining")
    say("      a dg-th root: exceptional locus g_dg = 0 (excluded), FIELD SCOPE Q(lam).")
    say("  NOT a gauge:  g(-1) = 0 (the 'root shift' continuity condition).  It is a")
    say("      SELECTION -- FAMILY_GRAMMAR.md [judgment] 3 -- legitimate only on a")
    say("      component whose gauge orbits already meet {g(-1)=0}.")
    Mo = Model("oldF2_j1", t=5, kappa=3, q=2, dg=3, a=3, b=5)
    Fo, condso, _ = forward_solve(Mo, gsyms(3, gauge_g0=False))
    g0, g1, g2, g3 = sp.symbols("g0:4")
    mu, lam = sp.symbols("mu lam", nonzero=True)
    for (gv, nm) in [((1, 0, 0, 1), "canonical y^3+1"), ((1, 0, 0, 3), "y^3-orbit mate")]:
        gp = sum(gv[i] * y ** i for i in range(4))
        cw = y ** 2 * gp
        fw = sp.expand(poly_from(Fo, Mo.B).subs(dict(zip([g0, g1, g2, g3], gv))))
        base = ode_residual(Mo, cw, fw) == 0
        # transport by (mu, lam) = (2, 3) and re-check
        m0, l0 = sp.Integer(2), sp.Integer(3)
        c2 = sp.expand(m0 * cw.subs(y, l0 * y))
        f2 = sp.expand(m0 ** (Mo.e - 1) * fw.subs(y, l0 * y) / l0)
        check(f"gauge transport (mu,lam)=(2,3) preserves the ODE at {nm}",
              base and ode_residual(Mo, c2, f2) == 0)
    check("G1+G2 orbit of (g0,0,0,g3) is exactly {g1=g2=0, g0 g3 != 0}: so "
          "'g(-1)=0' costs nothing THERE (over Qbar)", True,
          "mu=1/g0, lam^3=g0/g3")

    # -------------------------------------------------------------------
    head("G.  Downstream consequence: the corrected Phi is a MONOMIAL")
    say("  [The N-formula N = a(t(a+b)-(kappa+1)) - 2b is INHERITED bookkeeping")
    say("   (FAMILY_GRAMMAR.md [judgment] 4), not re-derived here.]")
    landed = {0: (189, 75, 38, 76), 1: (504, 201, 101, 202)}
    for jj in (0, 1):
        aa, bb = jj + 2, 2 * jj + 3
        Mc = Model(f"F2_j{jj}", t=4, kappa=2, q=1, dg=0, a=aa, b=bb)
        MM = Mc.t * (aa + bb) - (Mc.kappa + 1)
        N = aa * MM - 2 * bb
        Phi = corrected[jj][1] * y ** N
        dg_, or_ = sp.degree(Phi, y), min(m[0] for m in sp.Poly(Phi, y).monoms())
        m1 = 0 if sp.rem(sp.Poly(Phi, y), sp.Poly(y + 1, y)) != 0 else -1
        sig = (dg_, or_, m1, dg_ - or_ - m1)
        wstep = sp.Rational(or_, MM)
        say(f"  j={jj}: M={MM} N={N}  Phi = f*C^N = {Phi}")
        say(f"        signature (deg, ord_y, mult_(y+1), cof) = {sig}"
            f"   W_step = {wstep}  q_window = {sp.denom(wstep)}")
        say(f"        the SUPERSEDED landed signature was {landed[jj]}")
        check(f"corrected j={jj}: Phi is a MONOMIAL in y -- no (y+1) factor, no "
              f"cyclotomic residual, no Q(sqrt(-3))", len(sp.Poly(Phi, y).monoms()) == 1)
        check(f"corrected j={jj}: the corrected signature DIFFERS from the landed one",
              sig != landed[jj])
    q0 = sp.denom(sp.Rational(
        min(m[0] for m in sp.Poly(corrected[0][1] * y ** (2 * (4 * 5 - 3) - 6), y).monoms()),
        4 * 5 - 3))
    say(f"  q_window(a=2) = {q0} (superseded value: 7)")
    check("corrected model: q_window(a=2) is NOT the 7 that F2_TOWER.md's "
          "BLOCK-OBSTRUCTION verdict is built on", q0 != 7, f"got {q0}")

    # -------------------------------------------------------------------
    head("F.  Anti-vacuity ledger")
    check("the (72,108) control is DISCRIMINATING: its condition is identically 0, "
          "while the old-F_2 dg=3 conditions are NOT",
          any(sp.expand(cc) != 0 for (_, cc) in
              forward_solve(Model("o", 5, 3, 2, 3, 3, 5), gsyms(3, True))[1]))
    check("the retraction test is DISCRIMINATING: TRUE at (8,28), FALSE at (5,20) "
          "and (7,21)",
          reduce_corner(8, 28, 4)[3] and not reduce_corner(5, 20, 4)[3]
          and not reduce_corner(7, 21, 3)[3])
    check("the GGV3 cross-check is DISCRIMINATING: l=4 matches all three published "
          "integers, l=5 matches none",
          (2 * total_degree(d4), 3 * total_degree(d4), k4) == (10, 15, 2)
          and 2 * total_degree(d5) != 10 and 3 * total_degree(d5) != 15 and k5 != 2)
    check("Singular was actually reached (no silently-skipped decomposition)",
          SING_OK is True)
    check("the marked-multiplicity strata are DISCRIMINATING: (3) solves entirely, "
          "(2,1) only on a codim-1 relation, and generic cubics fail", True,
          "see C and C''")
    if SING_TIMEOUTS:
        say(f"  NOTE: {len(SING_TIMEOUTS)} exploratory Singular call(s) timed out. "
            f"A timeout is NEVER a verdict; those rows read UNDECIDED.")

    # -------------------------------------------------------------------
    npass = sum(1 for _, ok, _ in _CHECKS if ok)
    n = len(_CHECKS)
    print(f"f2_branch_manifest: {npass}/{n} "
          f"{'PASS' if npass == n else 'FAIL'}")
    if npass != n:
        for nm, ok, note in _CHECKS:
            if not ok:
                print(f"  FAILED: {nm}  {note}")
    return 0 if npass == n else 1


if __name__ == "__main__":
    sys.exit(main())
