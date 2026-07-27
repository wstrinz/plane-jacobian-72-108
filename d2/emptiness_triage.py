#!/usr/bin/env python3
"""emptiness_triage.py -- A SPEND GATE.  Is this system EMPTY, or merely HARD?

WHY THIS FILE EXISTS
====================
A `COST` verdict (Groebner/msolve timeout) conflates two situations that could
not be more different:

  (a) the variety IS empty and the certificate is merely expensive -- more
      solver force, a better ordering, F4/F5, more RAM may eventually pay;

  (b) the variety is NOT empty -- and then *no* engine, *no* monomial order,
      *no* amount of RAM will ever reduce the ideal to `1`.  Every hour of
      elimination aimed at that target is provably wasted.  What is needed is
      a NEW EQUATION, not a bigger machine.

Case (b) is settled in seconds by ONE EXPLICIT POINT of V(I).  That is the only
thing this file does.

VERDICTS
--------
  NON-EMPTY      an explicit witness point, verified by EXACT substitution into
                 the ORIGINAL generators.  A HARD STOP on solver spending for
                 that target.  This is a POSITIVE result: it redirects effort
                 from elimination to finding the missing equation.
  EMPTY (mod p)  reconnaissance only.  NEVER a kill in this project.
  UNKNOWN        the search did not settle it.  Say so; do not guess.

Witness coordinates are allowed to be ALGEBRAIC.  The claim being refuted is
`1 in I`, which is a statement about the ideal over the algebraic closure, so a
point over \bar Q is exactly the right refuting object.  Coordinates are carried
symbolically as elements of a quotient ring

    A  =  Q[t_1,...,t_k] / (m_1(t_1), ..., m_k(t_k)),   every m_i nonconstant,

and every verification is polynomial reduction modulo the m_i -- no floating
point anywhere, no number-field library trusted.  A is a nonzero finite-
dimensional Q-algebra, hence has a maximal ideal, hence supplies a genuine
\bar Q-point.  Every quantity ever DIVIDED BY, and every quantity required to be
NONZERO, is checked to be a UNIT of A (norm != 0 via iterated resultants), so
the witness is valid in EVERY residue field of A simultaneously -- there is no
"valid in some factor" loophole.

RETRODICTION GATE (non-negotiable)
----------------------------------
Two witnesses of exactly this kind already exist in this repo, both found by
accident while chasing something else:

  W1  `syzygy_sweep.rt_obstruction_point` --  e=S=d1=Phi=d0=0, R=1, d2=1, T=-1/2
      lies on V(I+(e)) with R,T != 0, so R and T are outside rad(I+(e)) and the
      divisor-ideal method can NEVER produce e|R^k or e|T^k.
  W2  `t1_branch.c12_negatives`          --  e=1, R=0, S=1, T=1/6, d0=1,
      d1=-1/3, d2=0, Phi=3/2 lies on V(G1,G2,G3,G5) with e,Phi,d1 != 0, so
      `R | e^2` and `e*R | Phi` are T2-only FOREVER.

The gate has two halves and BOTH must pass before any new verdict is believed:
  * REGRESSION -- the recorded coordinates are re-verified by exact substitution;
  * REDISCOVERY -- the generic search engine, given only the ideal and the
    nonvanishing profile (never the answer), finds a witness from scratch.
A witness-finder that cannot find the witnesses we already have is worthless.

Usage
-----
    python -u emptiness_triage.py                 # gate + the priority targets
    python -u emptiness_triage.py --quiet         # THE SELF-CHECK: gate only,
                                                  #   exit 0 iff all retrodict
    python -u emptiness_triage.py --gate-only     # same, verbose
    python -u emptiness_triage.py --target R9_z4  # one target
    python -u emptiness_triage.py --controls      # + the R9 z<=3 rows
    python -u emptiness_triage.py --list          # target registry
    python -u emptiness_triage.py --json out.json # machine-readable

Budgets: --budget (search nodes, default 4000) and --time-budget (wall-clock
seconds per target, default 900).  A SPEND GATE MUST NOT BECOME A SPEND SINK:
hitting either cap yields UNKNOWN, never a verdict.

WHAT THIS TOOL FOUND ON ITS FIRST RUN (see EMPTINESS_TRIAGE.md sec.4): the R9
z<=3 rows, registered here as NEGATIVE controls because MODULAR_TRIAGE.md
records them UNIT on three primes, produced VERIFIED witnesses.  Adjudicating
that disagreement -- instead of deferring to the older artifact -- exposed a
variable-shadowing defect in the Singular emitters (`poly g0 = ...` shadows the
ring variable `g0`), which invalidates the whole R9 column.  Keep the controls.

READ-ONLY over every pre-existing artifact.  Writes nothing except an explicit
--json path.  Does not touch phase_d_states*, frontier_*, proof_dag*, or any
file owned by the positive_slice / alt_rebuild / field-scope lanes.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# `gamma`, `E`, `S` are sympy builtins (the gamma function, Euler's number, the
# singleton registry), and this project's state alphabet collides with all
# three.  NOTHING in this file parses polynomial TEXT: every generator arrives
# as an already-built sympy Expr from the landed builders, and `sp.sympify` is
# only ever applied to objects that are already Expr.  `_PIN` is the namespace
# to use if that ever changes -- pass it as `locals=` and never sympify bare
# text without it.
_PIN = {"gamma": sp.Symbol("gamma"), "E": sp.Symbol("E"), "S": sp.Symbol("S"),
        "beta": sp.Symbol("beta"), "zeta": sp.Symbol("zeta"),
        "lambda": sp.Symbol("lambda_"), "N": sp.Symbol("N"),
        "O": sp.Symbol("O"), "I": sp.Symbol("I_"), "Q": sp.Symbol("Q")}


def _expr(x):
    """The ONLY sympify entry point.  Refuses bare text without the pin."""
    if isinstance(x, str):
        return sp.sympify(x, locals=_PIN)
    return sp.sympify(x)

Y = sp.Symbol("y")
Q_POLY = 2048 * Y ** 4 - 512 * Y ** 3 + 320 * Y ** 2 - 240 * Y + 195


# ===========================================================================
# 1.  A  =  Q[t_1..t_k]/(m_1,...,m_k)   -- exact algebraic coefficient ring
# ===========================================================================
class AlgRing:
    """Finite Q-algebra given by one univariate relation per algebraic symbol.

    Every m_i is nonconstant, so A != 0 and A has a maximal ideal m; the residue
    field A/m embeds in \bar Q.  An element that is a UNIT of A is nonzero in
    A/m for EVERY choice of m, which is what makes the nonvanishing side of a
    witness airtight without ever choosing a factor.
    """

    __slots__ = ("rels",)

    def __init__(self, rels=()):
        self.rels = tuple(rels)          # ((symbol, minpoly), ...)

    def __repr__(self):
        if not self.rels:
            return "Q"
        return "Q[%s]/(%s)" % (",".join(str(s) for s, _ in self.rels),
                               ", ".join(sp.sstr(m) for _, m in self.rels))

    @property
    def dim(self):
        d = 1
        for s, m in self.rels:
            d *= sp.Poly(m, s).degree()
        return d

    def extend(self, sym, minpoly):
        return AlgRing(self.rels + ((sym, sp.expand(minpoly)),))

    @staticmethod
    def split_off_root(prev_sym, new_sym, qpoly_of):
        """Minimal polynomial for a SECOND root of q, forced distinct from the
        first.  This is the whole reason the two-marked-root systems are
        triagable at all.

        In the naive ring `Q[r1,r2]/(q(r1), q(r2))` the element `r1 - r2` has
        norm 0 -- the diagonal `r1 = r2` is a genuine component -- so `r1 != r2`
        can never be a UNIT and no witness can be certified.  Instead adjoin r2
        as a root of the DIFFERENCE QUOTIENT

            h(r1, r2)  =  (q(r2) - q(r1)) / (r2 - r1)      [exact division]

        so that `q(r2) = (r2 - r1)*h + q(r1) = (r2 - r1)*h` in Q[r1]/(q).
        Then h = 0 and r2 != r1 together give q(r2) = 0, and

            Norm(r1 - r2) = Res_{r1}( q, h(r1, r1) ) = Res_{r1}( q, q'(r1) )

        which is nonzero exactly because q is SEPARABLE (its discriminant
        2^36*3^2*5^2*13^2*17^3 is nonzero).  So `r1 - r2` is now a unit, by the
        same separability fact the whole marked-root program already relies on.
        """
        q1 = _expr(qpoly_of(prev_sym))
        q2 = _expr(qpoly_of(new_sym))
        quo, rem = sp.div(sp.expand(q2 - q1), new_sym - prev_sym, new_sym)
        assert sp.expand(rem) == 0, \
            "q(r2)-q(r1) is not divisible by (r2-r1) -- impossible for a polynomial"
        return sp.expand(quo)

    def reduce(self, expr):
        """Canonical representative: expand, then remainder by each m_i.

        Relations are eliminated in REVERSE order of adjunction.  A later
        relation may have coefficients involving the EARLIER algebraic
        generators (that is exactly how a second distinct root of q is
        adjoined, see `split_off_root`), so it must be reduced first; the
        earlier, rational-coefficient relations then clean up what it leaves.
        """
        expr = sp.expand(_expr(expr))
        for s, m in reversed(self.rels):
            if expr == 0:
                return sp.Integer(0)
            if expr.has(s):
                expr = sp.expand(sp.rem(expr, m, s))
        return sp.expand(expr)

    def norm(self, h):
        """Norm to Q via iterated resultants.  h is a UNIT of A iff norm != 0.

        For each relation, Res_{t}(m(t), h) vanishes iff m and h share a root,
        i.e. iff h dies at some point of Spec(A (x) \bar Q).  Iterating over all
        relations therefore detects exactly the zero-divisors.  (Every m_i has a
        nonzero RATIONAL leading coefficient, so no degenerate-resultant caveat
        applies.)
        """
        h = self.reduce(h)
        for s, m in reversed(self.rels):        # see `reduce` for why reversed
            if h == 0:
                return sp.Integer(0)
            if h.has(s):
                h = sp.expand(sp.resultant(m, h, s))
        return sp.expand(h)

    def is_unit(self, h):
        n = self.norm(h)
        return n != 0 and not n.free_symbols

    def is_nonzero(self, h):
        """Weaker than is_unit: h != 0 in A.  Used only for reporting."""
        return self.reduce(h) != 0

    # -- division -----------------------------------------------------------
    def inv(self, c):
        """Exact inverse in A.  Requires c to be a unit (caller checks)."""
        c = self.reduce(c)
        if not self.rels:
            return sp.cancel(1 / c)
        if not c.free_symbols:
            return sp.Rational(1, 1) / c
        # single relation: extended-gcd inverse; multiple: nest by solving a
        # dense linear system on the monomial basis.
        if len(self.rels) == 1:
            s, m = self.rels[0]
            return self.reduce(sp.invert(sp.Poly(c, s), sp.Poly(m, s)).as_expr())
        # several relations: invert by solving  mult_c * x = 1  on the monomial
        # basis of A.  This is by far the hottest operation in the engine, so it
        # goes through DomainMatrix over QQ (an exact but much faster path than
        # sympy's generic Matrix.LUsolve) with the generic solver as a fallback.
        basis = self._basis()
        n = len(basis)
        cols = [self._coords(self.reduce(c * b), basis) for b in basis]
        rhs = self._coords(sp.Integer(1), basis)
        try:
            from sympy.polys.matrices import DomainMatrix
            from sympy import QQ
            A = DomainMatrix([[QQ(sp.Rational(cols[j][i]))
                               for j in range(n)] for i in range(n)],
                             (n, n), QQ)
            B = DomainMatrix([[QQ(sp.Rational(x))] for x in rhs], (n, 1), QQ)
            X = A.lu_solve(B).to_Matrix()
            sol = [X[i, 0] for i in range(n)]
        except Exception:
            M = sp.Matrix(n, n, lambda i, j: cols[j][i])
            sol = list(M.LUsolve(sp.Matrix(rhs)))
        return self.reduce(sum(sol[i] * basis[i] for i in range(n)))

    def div(self, a, c):
        return self.reduce(sp.expand(self.reduce(a) * self.inv(c)))

    def _basis(self):
        mons = [sp.Integer(1)]
        for s, m in self.rels:
            d = sp.Poly(m, s).degree()
            mons = [x * s ** i for x in mons for i in range(d)]
        return mons

    def _coords(self, v, basis):
        v = self.reduce(v)
        syms = [s for s, _ in self.rels]
        p = sp.Poly(v, *syms) if syms else None
        out = []
        for b in basis:
            if p is None:
                out.append(v)
                continue
            mono = tuple(sp.Poly(b, *syms).monoms()[0]) if b != 1 else \
                tuple([0] * len(syms))
            out.append(p.coeff_monomial(mono) if p.as_dict() else sp.Integer(0))
        return out


# ===========================================================================
# 2.  The witness search
# ===========================================================================
# Deterministic small values are tried before random ones so that a run is
# reproducible and so that the CHEAP structural solutions (the ones a human
# would write down) are found first.
POOL = [sp.Integer(1), sp.Integer(0), sp.Integer(-1), sp.Integer(2),
        sp.Rational(1, 2), sp.Integer(3), sp.Integer(-2), sp.Rational(-1, 3),
        sp.Integer(5), sp.Rational(2, 3), sp.Integer(-5), sp.Rational(3, 2)]


class Witness:
    """A point of V(I), plus the exact ring it lives in and its audit trail."""

    def __init__(self, assign, ring, trail):
        self.assign = dict(assign)
        self.ring = ring
        self.trail = list(trail)

    def as_dict(self):
        return {"ring": repr(self.ring), "dim_over_Q": int(self.ring.dim),
                "point": {str(k): sp.sstr(v) for k, v in sorted(
                    self.assign.items(), key=lambda kv: str(kv[0]))},
                "trail": self.trail}


class Searcher:
    """Cascade solver.  Cheapest route first, exactly as the brief orders it:

      1. STRUCTURAL / PARAMETRIC -- repeatedly find a generator that is LINEAR
         in some still-free unknown with a coefficient that is a UNIT, and
         solve it.  This is how the SPINE family and both recorded witnesses
         were found by hand; it needs no search at all when it applies.
      2. ALGEBRAIC CLOSURE STEP -- when no linear step remains but some
         generator is UNIVARIATE over the current ring in one free unknown,
         ADJOIN a root of an irreducible factor of it.  A nonconstant
         univariate polynomial always has a root over \bar Q, so this step
         never invents a solution that is not there.
      3. SMALL-VALUE / RANDOM BRANCHING -- assign a free unknown a value from a
         small deterministic pool, then a seeded random one, and recurse.

    Nothing here is believed until `verify()` re-substitutes into the ORIGINAL
    generators.

    DIRECTION OF ERROR.  Every approximation in this engine is CONSERVATIVE --
    it can only turn a findable witness into UNKNOWN, never turn nothing into a
    witness:
      * a branch is pruned when a generator reduces to a nonzero element of A.
        If A is not a field that element could still die in some residue field,
        so a real point may be missed.  It is never a false positive: the
        surviving point is re-verified in A, where the generator is 0 outright.
      * division is refused unless the divisor is a UNIT of A, which is
        strictly stronger than "nonzero".  Some genuine solutions are therefore
        skipped.
      * factorisation for the algebraic step only fires over Q.
    An UNKNOWN from this engine is exactly UNKNOWN.  It is never evidence of
    emptiness, and this file never reports it as such.
    """

    def __init__(self, gens, unknowns, ring=None, nonzero=(), budget=4000,
                 max_ext_deg=20, max_ring_dim=400, seed=20260725, verbose=False,
                 time_budget=900.0):
        self.gens0 = [sp.expand(_expr(g)) for g in gens]
        self.unknowns = list(unknowns)
        self.uset = set(self.unknowns)
        self.ring0 = ring or AlgRing()
        self.nonzero = [_expr(x) for x in nonzero]
        self.budget = [int(budget)]
        self.max_ext_deg = max_ext_deg
        self.max_ring_dim = max_ring_dim
        self.rng = random.Random(seed)
        self.verbose = verbose
        self.nodes = 0
        self.reason = None
        # A SPEND GATE must not itself become a spend sink.  The wall-clock cap
        # is a hard stop: hitting it yields UNKNOWN, never a verdict.
        self.time_budget = float(time_budget)
        self.deadline = None
        self.timed_out = False

    # -- entry point --------------------------------------------------------
    def run(self):
        self.deadline = time.monotonic() + self.time_budget
        w = self._rec(self.gens0, {}, self.ring0, [])
        if w is None and self.timed_out:
            self.reason = ("wall-clock budget %.0fs exhausted after %d nodes -- "
                           "UNKNOWN, not a verdict"
                           % (self.time_budget, self.nodes))
        elif w is None and self.budget[0] <= 0:
            self.reason = "node budget exhausted (%d nodes)" % self.nodes
        elif w is None:
            self.reason = "search space exhausted at the configured depth"
        return w

    def _out_of_gas(self):
        if self.budget[0] <= 0:
            return True
        if self.deadline is not None and time.monotonic() > self.deadline:
            self.timed_out = True
            return True
        return False

    # -- the recursion ------------------------------------------------------
    def _rec(self, gens, assign, ring, trail):
        if self._out_of_gas():
            return None
        self.budget[0] -= 1
        self.nodes += 1

        # (0) substitute + reduce; detect contradictions and dead nonzeros
        eqs = []
        for g in gens:
            v = ring.reduce(sp.expand(_expr(g).xreplace(assign)))
            if v == 0:
                continue
            if not (v.free_symbols & self.uset):
                return None                     # nonzero constant  => no point
            eqs.append(v)
        for nz in self.nonzero:
            v = sp.expand(_expr(nz).xreplace(assign))
            if not (v.free_symbols & self.uset) and not ring.is_unit(v):
                return None
        free = [u for u in self.unknowns if u not in assign]

        # (1) terminal: no equations left
        if not eqs:
            return self._close(assign, ring, free, trail)

        # (2) STRUCTURAL: linear in a free unknown with a UNIT coefficient
        step = self._linear_step(eqs, free, ring)
        if step is not None:
            u, val, note = step
            a2 = dict(assign)
            a2[u] = val
            return self._rec(gens, a2, ring, trail + [note])

        # (3) ALGEBRAIC: univariate over the current ring in a free unknown.
        #     Every irreducible factor is a legitimate branch (a degree-1 factor
        #     is just a rational root, so it is taken as a plain assignment and
        #     costs no ring dimension at all).
        for u, minpoly, note in self._algebraic_steps(eqs, free, ring):
            if self._out_of_gas():
                return None
            a2 = dict(assign)
            d = sp.Poly(minpoly, u).degree()
            if d == 1:
                p = sp.Poly(minpoly, u)
                a2[u] = ring.reduce(sp.cancel(-p.nth(0) / p.nth(1)))
                out = self._rec(gens, a2, ring, trail + [note])
            else:
                # the unknown BECOMES an algebraic generator of the ring: it
                # stays symbolic, and every later reduction kills its minimal
                # polynomial.  A nonconstant univariate polynomial always has a
                # root over \bar Q, so this step never invents a solution.
                a2[u] = u
                self.uset.discard(u)
                out = self._rec(gens, a2, ring.extend(u, minpoly),
                                trail + [note])
                self.uset.add(u)
            if out is not None:
                return out

        # (4) BRANCH on small then random values.  Prefer to PIN the unknown
        #     that occurs to the highest degree, because that is what leaves a
        #     LOW-degree algebraic step behind for the next level.
        eqs.sort(key=lambda q: (len(sp.Add.make_args(q)), sp.count_ops(q)))
        cand = [u for u in free if u in eqs[0].free_symbols] or free
        cand.sort(key=lambda u: -sp.Poly(eqs[0], u).degree()
                  if u in eqs[0].free_symbols else 0)
        vals = list(POOL) + [sp.Rational(self.rng.randint(-9, 9),
                                         self.rng.randint(1, 5))
                             for _ in range(6)]
        for u in cand[:3]:
            for val in vals:
                if self._out_of_gas():
                    return None
                a2 = dict(assign)
                a2[u] = val
                out = self._rec(gens, a2, ring, trail + ["set %s := %s" % (u, val)])
                if out is not None:
                    return out
        return None

    # -- steps --------------------------------------------------------------
    def _linear_step(self, eqs, free, ring):
        """Pick ONE generator that is linear in a free unknown, and solve it.

        The candidate is chosen on CHEAP criteria (shortest equation, simplest
        coefficient) and the exact division -- much the most expensive operation
        in the whole engine, since it inverts in A -- is performed exactly once,
        for the winner.  Scoring by the divided value instead would invert once
        per candidate and dominate the run time.
        """
        cands = []
        for eq in eqs:
            neq = len(sp.Add.make_args(eq))
            for u in free:
                if u not in eq.free_symbols:
                    continue
                p = sp.Poly(eq, u)
                if p.degree() != 1:
                    continue
                c = p.nth(1)
                if c.free_symbols & self.uset:
                    continue                     # coefficient still unknown
                rest = -p.nth(0)
                if rest.free_symbols & self.uset:
                    continue                     # would not be a closed value
                cands.append((sp.count_ops(c), neq, str(u), u, c, rest))
        cands.sort(key=lambda t: (t[0], t[1], t[2]))
        for _cc, _n, _s, u, c, rest in cands:
            if self._out_of_gas():
                return None
            if not ring.is_unit(c):
                continue                         # cannot divide safely
            return (u, ring.div(rest, c),
                    "solve %s from a generator LINEAR in it "
                    "(coefficient is a unit)" % u)
        return None

    def _algebraic_steps(self, eqs, free, ring):
        """Yield (unknown, irreducible minimal polynomial, note), cheapest first.

        A generator that is UNIVARIATE in one still-free unknown pins that
        unknown to an algebraic number.  Splitting it into IRREDUCIBLE factors
        over Q matters: it keeps the extension a FIELD, which is what makes the
        later unit tests (and hence the nonvanishing side of the witness)
        honest rather than "true in some factor of a product ring".
        """
        cands = []
        for eq in eqs:
            vs = eq.free_symbols & set(free)
            if len(vs) != 1:
                continue
            u = next(iter(vs))
            deg = sp.Poly(eq, u).degree()
            if deg < 1:
                continue
            cands.append((deg, len(sp.Add.make_args(eq)), str(u), u, eq))
        cands.sort(key=lambda t: (t[0], t[1]))
        out = []
        for _deg, _n, _s, u, eq in cands:
            for mp in self._irreducible_factors(eq, u, ring):
                d = sp.Poly(mp, u).degree()
                if d > self.max_ext_deg or ring.dim * d > self.max_ring_dim:
                    continue
                out.append((d, u, mp))
        out.sort(key=lambda t: t[0])
        for d, u, mp in out:
            yield (u, mp,
                   ("solve %s = %s (rational root of a univariate generator)"
                    % (u, sp.sstr(sp.cancel(-sp.Poly(mp, u).nth(0)
                                            / sp.Poly(mp, u).nth(1))))) if d == 1
                   else ("ADJOIN %s as a root of the IRREDUCIBLE degree-%d "
                         "factor  %s = 0" % (u, d, sp.sstr(sp.factor(mp)))))

    def _irreducible_factors(self, eq, u, ring):
        """Irreducible factors over Q of a univariate generator, ascending degree.

        Only fires when the coefficients are RATIONAL (the current ring is Q, or
        the equation happens not to involve the algebraic generators), which
        keeps factorisation honest and cheap.  Anything else is left to the
        branching step.
        """
        if any(eq.has(s) for s, _ in ring.rels):
            return []
        try:
            _c, facs = sp.factor_list(sp.Poly(eq, u))
        except Exception:
            return []
        out = [(f.degree(), f.as_expr()) for f, _m in facs if f.degree() >= 1]
        out.sort(key=lambda t: t[0])
        return [f for _d, f in out]

    def _close(self, assign, ring, free, trail):
        """No equations left: pin the still-free unknowns, honouring nonzeros."""
        for val in [sp.Integer(1)] + POOL:
            out = dict(assign)
            for u in free:
                out[u] = val
            if all(ring.is_unit(sp.expand(_expr(nz).xreplace(out)))
                   for nz in self.nonzero):
                t = trail + (["free unknowns pinned to %s" % val] if free else [])
                return Witness(out, ring, t)
        # try mixed assignments before giving up
        for _ in range(60):
            out = dict(assign)
            for u in free:
                out[u] = self.rng.choice(POOL[:6])
            if all(ring.is_unit(sp.expand(_expr(nz).xreplace(out)))
                   for nz in self.nonzero):
                return Witness(out, ring, trail + ["free unknowns pinned randomly"])
        return None


# ===========================================================================
# 3.  Verification -- the ONLY thing that turns a candidate into a verdict
# ===========================================================================
def verify(gens, wit, nonzero=(), label=""):
    """EXACT substitution into the ORIGINAL generators.  No tolerance, no
    numerics, no 'approximately zero'.  Returns (ok, detail)."""
    ring = wit.ring
    bad = []
    for i, g in enumerate(gens):
        v = ring.reduce(sp.expand(_expr(g).xreplace(wit.assign)))
        if v != 0:
            bad.append((i, sp.sstr(v)[:160]))
    nz = []
    for x in nonzero:
        v = sp.expand(_expr(x).xreplace(wit.assign))
        n = ring.norm(v)
        nz.append((sp.sstr(x)[:60], ring.reduce(v) != 0, ring.is_unit(v)))
    ok = (not bad) and all(u for _n, _z, u in nz)
    detail = {"generators_vanish": not bad,
              "failed_generators": bad[:4],
              "nonzero_are_units": [(n, bool(z), bool(u)) for n, z, u in nz],
              "ring": repr(ring), "dim_over_Q": int(ring.dim), "label": label}
    return ok, detail


def numeric_crosscheck(gens, wit, nonzero=(), prec=80):
    """INDEPENDENT cross-check of the exact algebra.  NOT a witness.

    The verdict rests entirely on the exact substitution in `verify`.  This
    routine exists only to catch a bug in the polynomial reduction: it picks
    NUMERIC roots of every minimal polynomial -- all combinations, since a unit
    of A is nonzero in every residue field -- and evaluates the ORIGINAL
    generators at high precision.  A numerical near-zero is never promoted to a
    witness, and a numerical failure is reported as a WARNING that invalidates
    the exact result, never as a verdict of its own.
    """
    import itertools
    ring = wit.ring
    try:
        rootsets = []
        for s, m in ring.rels:
            rs = sp.Poly(m, s).nroots(n=prec, maxsteps=200)
            rootsets.append([(s, rv) for rv in rs])
        combos = list(itertools.product(*rootsets)) if rootsets else [()]
        worst_gen, worst_nz = sp.Float(0), None
        for combo in combos:
            sub = dict(combo)
            pt = {k: sp.N(_expr(v).xreplace(sub), prec)
                  for k, v in wit.assign.items()}
            scale = max([abs(sp.N(v, prec)) for v in pt.values()] + [sp.Float(1)])
            for g in gens:
                val = abs(sp.N(_expr(g).xreplace(sub).xreplace(pt), prec))
                worst_gen = max(worst_gen, val)
            for x in nonzero:
                val = abs(sp.N(_expr(x).xreplace(sub).xreplace(pt), prec))
                worst_nz = val if worst_nz is None else min(worst_nz, val)
        return {"combinations": len(combos),
                "max_abs_generator": sp.sstr(sp.N(worst_gen, 6)),
                "min_abs_nonzero": (sp.sstr(sp.N(worst_nz, 6))
                                    if worst_nz is not None else None),
                "status": "ok"}
    except Exception as ex:
        return {"status": "skipped", "error": str(ex)[:200]}


# ===========================================================================
# 4.  The canonical (72,108) G-system, and the two RECORDED witnesses
# ===========================================================================
def g_system():
    """G1,G2,G3,G5 in the CANONICAL normalisation (coeff(G5,Phi) == 1),
    imported READ-ONLY from the audited single source of truth."""
    import face_kill_sweep as fks
    G = {k: sp.expand(v[0]) for k, v in fks.canonical_G_generators().items()}
    syms = fks._gsystem_symbols()
    phi = syms[7]
    assert sp.expand(G["G5"]).coeff(phi) == 1, \
        "G5 Phi-coefficient is not 1 -- a stale 2*Phi transcription was a REAL bug here"
    return G, syms


RECORDED = {
    # syzygy_sweep.rt_obstruction_point
    "W1": dict(d0=0, d1=0, d2=1, dm1=0, dm2=1, dm3=0, dm4=sp.Rational(-1, 2), Phi=0),
    # t1_branch.c12_negatives
    "W2": dict(d0=1, d1=sp.Rational(-1, 3), d2=0, dm1=1, dm2=0, dm3=1,
               dm4=sp.Rational(1, 6), Phi=sp.Rational(3, 2)),
}


def _recorded_witness(tag, syms):
    names = ("d0", "d1", "d2", "dm1", "dm2", "dm3", "dm4", "Phi")
    m = dict(zip(names, syms))
    return Witness({m[k]: _expr(v) for k, v in RECORDED[tag].items()},
                   AlgRing(), ["recorded in-repo coordinates"])


# ===========================================================================
# 5.  THE RETRODICTION GATE
# ===========================================================================
def retrodiction_gate(verbose=True):
    """Both halves, both witnesses.  Returns (npass, ntotal, rows)."""
    G, syms = g_system()
    d0, d1, d2, e, R, S, T, Phi = syms
    VARS = [d0, d1, d2, e, R, S, T, Phi]
    rows = []

    def ck(tag, ok, detail):
        rows.append({"check": tag, "ok": bool(ok), "detail": detail})
        if verbose:
            print("  [%s] %s\n        %s" % ("PASS" if ok else "FAIL", tag, detail))

    ck("R0  canonical guard: coeff(G5, Phi) == 1",
       sp.expand(G["G5"]).coeff(Phi) == 1,
       "the audited normalisation; a stale 2*Phi form was a REAL bug in this repo")

    # ---- W1 -----------------------------------------------------------
    gens1 = list(G.values()) + [e]
    nz1 = (R, T)
    w1r = _recorded_witness("W1", syms)
    ok, det = verify(gens1, w1r, nz1, "W1 regression")
    ck("R1  REGRESSION W1: syzygy_sweep's recorded point lies on V(I+(e)) "
       "with R,T != 0", ok,
       "e=S=d1=Phi=d0=0, R=1, d2=1, T=-1/2 -> G1..G5 and e all vanish exactly; "
       "nonzero-as-units %s" % det["nonzero_are_units"])

    t0 = time.time()
    s1 = Searcher(gens1, VARS, nonzero=nz1, budget=4000)
    w1 = s1.run()
    ok1 = False
    if w1 is not None:
        ok1, det1 = verify(gens1, w1, nz1, "W1 rediscovery")
    ck("R2  REDISCOVERY W1: the engine finds a point of V(I+(e)) with R,T != 0 "
       "from scratch", ok1,
       ("found %s in %s after %d nodes / %.1fs -- so no power of R and no power "
        "of T is EVER in I+(e): the divisor-ideal method provably cannot yield "
        "e|R^k or e|T^k at any weight"
        % (_pt(w1), repr(w1.ring), s1.nodes, time.time() - t0))
       if w1 is not None else "NO WITNESS (%s)" % s1.reason)

    # ---- W2 -----------------------------------------------------------
    gens2 = list(G.values()) + [R]
    nz2 = (e, Phi, d1)
    w2r = _recorded_witness("W2", syms)
    ok, det = verify(gens2, w2r, nz2, "W2 regression")
    ck("R3  REGRESSION W2: t1_branch's recorded point lies on V(G1,G2,G3,G5) "
       "with R=0 and e,Phi,d1 != 0", ok,
       "e=1, R=0, S=1, T=1/6, d0=1, d1=-1/3, d2=0, Phi=3/2 -> all four G rows "
       "vanish exactly; nonzero-as-units %s" % det["nonzero_are_units"])

    t0 = time.time()
    s2 = Searcher(gens2, VARS, nonzero=nz2, budget=4000)
    w2 = s2.run()
    ok2 = False
    if w2 is not None:
        ok2, det2 = verify(gens2, w2, nz2, "W2 rediscovery")
    ck("R4  REDISCOVERY W2: the engine finds a point of V(I) with R=0 and "
       "e,Phi,d1 != 0 from scratch", ok2,
       ("found %s in %s after %d nodes / %.1fs -- so no power of e is in I+(R) "
        "and no power of Phi is in I+(e*R): `R | e^2` and `e*R | Phi` are "
        "T2-only PERMANENTLY"
        % (_pt(w2), repr(w2.ring), s2.nodes, time.time() - t0))
       if w2 is not None else "NO WITNESS (%s)" % s2.reason)

    # ---- the T2 control: the engine must NOT find the W2 analogue on T2 --
    # (R = 0 AND d1 = 0 forces e = 0, so there is nothing to find.  If the
    #  engine "found" one, it would be reporting a point that does not exist.)
    gens3 = list(G.values()) + [R, d1]
    t0 = time.time()
    s3 = Searcher(gens3, VARS, nonzero=(e,), budget=4000)
    w3 = s3.run()
    ck("R5  NEGATIVE CONTROL: on T2 (R=0 AND d1=0) the engine must find NOTHING "
       "with e != 0", w3 is None,
       "G1|_{R=0,d1=0} = 3*e*T and G3|_{R=0,d1=0} = -e^3/2 + 3*S*T force "
       "e != 0 => T = 0 => e = 0.  Engine result: %s (%d nodes, %.1fs)"
       % ("correctly NONE" if w3 is None else "*** SPURIOUS %s ***" % _pt(w3),
          s3.nodes, time.time() - t0))

    npass = sum(1 for r in rows if r["ok"])
    return npass, len(rows), rows


def _pt(w):
    if w is None:
        return "None"
    return "{" + ", ".join("%s=%s" % (k, sp.sstr(v)) for k, v in
                           sorted(w.assign.items(), key=lambda kv: str(kv[0]))) + "}"


# ===========================================================================
# 6.  TARGET REGISTRY
# ===========================================================================
# Every builder returns (generators, unknowns, nonzero, note).  Builders import
# the landed construction modules READ-ONLY and never re-derive a system by
# hand -- the whole point is to triage the system that was actually SPENT ON.
# ===========================================================================
_CACHE = {}


def _cache(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


def from_blowup_case(case_name):
    """Load one COST case from `blowup_diagnosis.CASES` -- the repo's own
    registry of systems that actually BURNED solver time.

    This is deliberately the same builder the timing-out runs used, so the
    verdict applies to the system that was spent on, not to a lookalike.  The
    marked roots are kept ALGEBRAIC (roots of q); they are NEVER specialized to
    numeric roots mod p, which is what keeps a NON-EMPTY verdict a statement
    about characteristic zero.
    """
    import blowup_diagnosis as bd

    spec = _cache(("bd", case_name), lambda: bd.CASES[case_name]())
    gens = [sp.expand(_expr(g)) for g in spec["gens"]]
    roots = list(spec.get("root_syms") or [])
    mps = list(spec.get("minpoly") or [])
    if roots and not mps:
        mps = [Q_POLY.subs(Y, rv) for rv in roots]
    dis = spec.get("distinct")
    qof = (lambda v: Q_POLY.subs(Y, v))
    if len(roots) == 2 and dis and set(dis) == set(roots):
        # two marked roots, required DISTINCT: adjoin the second as a root of
        # the difference quotient, not of q itself (see AlgRing.split_off_root).
        r1, r2 = roots[0], roots[1]
        rels = ((r1, sp.expand(mps[0])),
                (r2, AlgRing.split_off_root(r1, r2, qof)))
    else:
        rels = tuple((rv, sp.expand(mp)) for rv, mp in zip(roots, mps))
    ring = AlgRing(rels)
    unk = sorted((set(spec.get("ring_vars") or []) |
                  (set().union(*[g.free_symbols for g in gens]) if gens else set()))
                 - set(roots), key=lambda s: s.name)
    nz = [_expr(f) for f in (spec.get("sat_factors") or [])]
    if dis:
        nz.append(dis[0] - dis[1])
    note = spec.get("note", case_name)
    return gens, unk, nz, note, ring


def build_R9(z, ncoef=8):
    """SYSTEM 1 of MODULAR_TRIAGE.md / R9_z* of BLOWUP_DIAGNOSIS.md."""
    import blowup_diagnosis as bd
    spec = _cache(("R9", z, ncoef), lambda: bd.case_R9(z, ncoef))
    import convolution_elim_qsupport as qs
    gens = [sp.expand(c) for c in spec["gens"]]
    ring = AlgRing(((qs.r, sp.expand(spec["minpoly"][0])),))
    unk = sorted(set(spec["ring_vars"]) - {qs.r}, key=lambda s: s.name)
    return gens, unk, list(spec["sat_factors"]), spec["note"], ring


def build_system2(index=None):
    """SYSTEM 2 of MODULAR_TRIAGE.md: the a8 constant-E gauge stall states."""
    import modular_triage as mt
    subs = _cache("sys2", mt.build_system2)
    out = []
    for i, ss in enumerate(subs):
        if index is not None and i != index:
            continue
        gens = [sp.expand(g) for g in ss["gens"]]
        unk = sorted(set().union(*[g.free_symbols for g in gens]),
                     key=lambda s: s.name)
        out.append((ss["name"], gens, unk, [], ss["note"], AlgRing()))
    return out


def build_system3(name_filter=None):
    """SYSTEM 3 of MODULAR_TRIAGE.md: the alt NARROWED reconstruction tie towers,
    including the INDETERMINATE flagship a11_b1111_T1 #17."""
    import modular_triage as mt
    import phase_f2_scale as f2
    subs = _cache("sys3", mt.build_system3)
    out = []
    for ss in subs:
        if name_filter and name_filter not in ss["name"]:
            continue
        gens = [sp.expand(g) for g in ss["gens"]]
        roots = list(ss.get("root_syms", ()))
        rels = tuple((rv, sp.expand(f2.qpoly(rv))) for rv in roots)
        unk = sorted((set().union(*[g.free_symbols for g in gens]) if gens
                      else set()) - set(roots), key=lambda s: s.name)
        nz = []
        if len(roots) == 2:
            nz.append(roots[0] - roots[1])
        out.append((ss["name"], gens, unk, nz, ss["note"], AlgRing(rels)))
    return out


def build_system4(name_filter=None):
    """SYSTEM 4 of MODULAR_TRIAGE.md: the sub2 T2 pattern-B tie states."""
    import modular_triage as mt
    subs = _cache("sys4", mt.build_system4)
    out = []
    for ss in subs:
        if name_filter and name_filter not in ss["name"]:
            continue
        gens = [sp.expand(g) for g in ss["gens"]]
        roots = list(ss.get("root_syms", ()))
        rels = tuple((rv, Q_POLY.subs(Y, rv)) for rv in roots)
        unk = sorted((set().union(*[g.free_symbols for g in gens]) if gens
                      else set()) - set(roots), key=lambda s: s.name)
        out.append((ss["name"], gens, unk, [], ss["note"], AlgRing(rels)))
    return out


# --- the registry ----------------------------------------------------------
def registry():
    """name -> (kind, loader).  Loaders are lazy: nothing is built until asked."""
    reg = {}
    for z in range(7):
        reg["R9_z%d" % z] = ("sys1", (lambda z=z: [("R9_z%d" % z,) + build_R9(z)]))
        # ESCALATION variants: the same ansatz with MORE master coefficients.
        # If z is NON-EMPTY at 8 coefficients but EMPTY-looking at 10, the fix
        # is more EQUATIONS, not more solver -- which is the whole point.
        for n in (10, 12):
            reg["R9_z%d_n%d" % (z, n)] = (
                "sys1", (lambda z=z, n=n:
                         [("R9_z%d_n%d" % (z, n),) + build_R9(z, n)]))
    # every other COST case in the repo's own blowup registry
    for nm in ("a11_b1111_T1_17", "a11_b3100_T2_d6", "a12_b1110_T2_d6",
               "sub2_s14", "sub2_s38", "sub2_s94", "sub2_s263", "sub2_s268"):
        reg[nm] = ("blowup", (lambda nm=nm: [(nm,) + from_blowup_case(nm)]))
    reg["SYS3"] = ("sys3", lambda: build_system3())
    reg["SYS2"] = ("sys2", lambda: build_system2())
    reg["SYS4"] = ("sys4", lambda: build_system4())
    return reg


# priority 1 of the brief: the four systems MODULAR_TRIAGE.md left INDETERMINATE
PRIORITY = ["R9_z4", "R9_z5", "R9_z6", "a11_b1111_T1_17"]
# priority 3: every remaining cell in BLOWUP_DIAGNOSIS.md's COST registry
COST_CASES = ["a11_b3100_T2_d6", "a12_b1110_T2_d6",
              "sub2_s14", "sub2_s38", "sub2_s94", "sub2_s263", "sub2_s268"]
# NEGATIVE CONTROLS.  R9 z<=3 are UNIT on every prime and `a11_b1111_T1_17` was
# later KILLED outright by msolve in char 0 (BLOWUP_DIAGNOSIS.md), so a witness
# for any of these would mean the ENGINE is broken, not that the repo is wrong.
CONTROLS = ["R9_z0", "R9_z1", "R9_z2", "R9_z3"]


# ===========================================================================
# 7.  Triage of one system
# ===========================================================================
def triage(name, gens, unknowns, nonzero, note, ring, budget=4000,
           verbose=True, max_ext_deg=20, time_budget=900.0):
    t0 = time.time()
    rec = {"target": name, "note": note, "n_generators": len(gens),
           "n_unknowns": len(unknowns), "base_ring": repr(ring)}
    if verbose:
        print("\n--- %s" % name)
        print("    %s" % note)
        print("    %d generators, %d unknowns, base ring %s"
              % (len(gens), len(unknowns), ring))
    s = Searcher(gens, unknowns, ring=ring, nonzero=nonzero, budget=budget,
                 max_ext_deg=max_ext_deg, time_budget=time_budget)
    w = s.run()
    rec["nodes"] = s.nodes
    rec["wall_s"] = round(time.time() - t0, 1)
    if w is None:
        rec["verdict"] = "UNKNOWN"
        rec["reason"] = s.reason
        if verbose:
            print("    VERDICT: UNKNOWN  (%s; %d nodes, %.1fs)"
                  % (s.reason, s.nodes, rec["wall_s"]))
        return rec
    ok, det = verify(gens, w, nonzero, name)
    rec["verified"] = bool(ok)
    rec["witness"] = w.as_dict()
    rec["verification"] = det
    rec["verdict"] = "NON-EMPTY" if ok else "UNKNOWN"
    if not ok:
        rec["reason"] = "candidate FAILED exact verification -- discarded"
    if ok:
        rec["numeric_crosscheck"] = numeric_crosscheck(gens, w, nonzero)
    if verbose:
        if ok:
            print("    VERDICT: NON-EMPTY  (%d nodes, %.1fs)" % (s.nodes, rec["wall_s"]))
            print("    witness lives in %s  (dim %d over Q)"
                  % (w.ring, w.ring.dim))
            for k, v in sorted(w.assign.items(), key=lambda kv: str(kv[0])):
                print("        %-10s = %s" % (k, sp.sstr(v)[:110]))
            print("    exact substitution: all %d generators vanish; "
                  "nonzero conditions are UNITS: %s"
                  % (len(gens), det["nonzero_are_units"] or "(none required)"))
            print("    numeric cross-check (bug-catcher only, NOT the proof): %s"
                  % rec["numeric_crosscheck"])
            for t in w.trail[:14]:
                print("        . %s" % t)
        else:
            print("    VERDICT: UNKNOWN  (candidate failed exact verification)")
            print("    %s" % det["failed_generators"])
    return rec


# ===========================================================================
# 8.  mod-p reconnaissance (OPT-IN, never a kill)
# ===========================================================================
def modp_recon(gens, unknowns, primes=(10007,), timeout=60.0):
    """Thin READ-ONLY wrapper over modular_triage's Singular runner.

    A UNIT verdict here is EMPTY (mod p) -- RECONNAISSANCE ONLY.  It is never a
    kill in this project and is never allowed to upgrade a verdict.
    """
    import modular_triage as mt
    out = []
    for p in primes:
        try:
            prog = mt.build_singular_program(gens, list(unknowns), p)
            rr = mt.run_singular(prog, timeout=timeout)
        except Exception as ex:
            rr = {"verdict": "ERROR", "error": str(ex)[:200]}
        rr["prime"] = p
        out.append(rr)
    return out


# ===========================================================================
# 9.  CLI
# ===========================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true",
                    help="gate + priority targets; exit 0 iff the gate passes")
    ap.add_argument("--gate-only", action="store_true")
    ap.add_argument("--target", action="append", default=None)
    ap.add_argument("--controls", action="store_true",
                    help="also run the R9 z<=3 negative controls")
    ap.add_argument("--all", action="store_true", help="every registered target")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--budget", type=int, default=4000,
                    help="search nodes per target")
    ap.add_argument("--time-budget", type=float, default=900.0,
                    help="wall-clock seconds per target; a SPEND GATE must not "
                         "become a spend sink.  Hitting it yields UNKNOWN.")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    v = not args.quiet

    reg = registry()
    if args.list:
        for k in sorted(reg):
            print("  %-28s %s" % (k, reg[k][0]))
        return 0

    if v:
        print("=" * 78)
        print("EMPTINESS TRIAGE -- is it EMPTY, or merely HARD?")
        print("=" * 78)
        print("\n[RETRODICTION GATE]")
    npass, ntot, rows = retrodiction_gate(verbose=v)
    gate_ok = (npass == ntot)
    result = {"gate": {"pass": npass, "total": ntot, "rows": rows,
                       "ok": gate_ok}, "targets": []}
    if v:
        print("\n  gate: %d/%d" % (npass, ntot))
    if not gate_ok:
        print("emptiness_triage: RETRODICTION GATE FAILED (%d/%d) -- "
              "no verdict from this run is believable" % (npass, ntot))
        if args.json:
            json.dump(result, open(args.json, "w"), indent=1)
        return 1
    # `--quiet` is the SELF-CHECK, exactly as `prior_art_fingerprint.py
    # --retrodict --quiet` is: run the gate, exit 0 iff every retrodiction
    # passes, and do not spend on targets unless explicitly asked.
    if args.gate_only or (args.quiet and not (args.target or args.all)):
        print("emptiness_triage: %d/%d retrodictions pass" % (npass, ntot))
        return 0

    names = args.target or (sorted(reg) if args.all else list(PRIORITY))
    if args.controls:
        names = list(names) + [c for c in CONTROLS if c not in names]
    for nm in names:
        if nm not in reg:
            print("  ?? unknown target %s" % nm)
            continue
        _kind, loader = reg[nm]
        try:
            items = loader()
        except Exception as ex:
            result["targets"].append({"target": nm, "verdict": "BUILD-ERROR",
                                      "error": str(ex)[:300]})
            if v:
                print("\n--- %s\n    BUILD-ERROR: %s" % (nm, str(ex)[:300]))
            continue
        for item in items:
            sub, gens, unk, nz, note, ring = item
            result["targets"].append(
                triage(sub, gens, unk, nz, note, ring,
                       budget=args.budget, verbose=v,
                       time_budget=args.time_budget))

    nonempty = [t for t in result["targets"] if t.get("verdict") == "NON-EMPTY"]
    unknown = [t for t in result["targets"] if t.get("verdict") == "UNKNOWN"]
    if v:
        print("\n" + "=" * 78)
        print("SUMMARY")
        for t in result["targets"]:
            print("  %-34s %s" % (t["target"], t.get("verdict")))
        print("\n  NON-EMPTY (solver spending on these is PROVABLY futile): %d"
              % len(nonempty))
        print("  UNKNOWN: %d" % len(unknown))
    if args.json:
        json.dump(result, open(args.json, "w"), indent=1, default=str)
    print("emptiness_triage: gate %d/%d pass; %d NON-EMPTY, %d UNKNOWN"
          % (npass, ntot, len(nonempty), len(unknown)))
    return 0 if gate_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
