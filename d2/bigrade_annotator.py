#!/usr/bin/env python3
"""bigrade_annotator.py  (NEW; read-only over all existing artifacts)

THE ENGINE BUILD, MILESTONE 1 -- the read-only bigrade annotator (diagnosis only;
no new solver).  For an existing coefficient system it:

  1. attaches to every scalar coefficient variable and every coefficient equation
     the bidegree (w, nu) = (u-weight, y-order), the integral defect
     (w, delta = Q*nu - P*w), and a place / root-algebra tag, then VERIFIES every
     monomial of every equation is bidegree-consistent (multiplication adds
     degrees);
  2. builds the bipartite incidence {equation blocks} <-> {variable blocks}
     grouped by bidegree and runs: connected components, a Dulmage-Mendelsohn /
     strongly-connected block-triangular decomposition, generic modular rank per
     block, extremal-face extraction under delta / w / mixed functionals, and a
     left-nullspace search on small terminal blocks;
  3. runs on THREE systems:
       R1 = the (50,75) window system  (regression: rediscover c_{0,-10});
       R2 = a known fast home-case control kill (full G-system, a8 pilot state);
       R3 = the walled R9 z=1 dm4-eliminated H-system (r9_eliminated_system.json).

Exact sympy throughout.  Sources (read-only): ENDPOINT_CONTRACT.md, F2_TOWER.md /
f2_tower.py, WINDOW_FUNCTIONS_75_125.md, FULL_SYSTEM_BRIDGE.md,
r9_eliminated_system.json, convolution_elim_qsupport.py, full_system_bridge_pilot.json.

Independent checker: bigrade_annotator_verify.py (--quiet, exit 0; the R1
rediscovery is the load-bearing check).
"""
from __future__ import annotations

import json
import os
import random
import sys
from collections import defaultdict

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
y = sp.symbols("y")


# =====================================================================
#  0.  Data model
# =====================================================================
class Var:
    __slots__ = ("sym", "name", "w", "nu", "place", "status", "block")

    def __init__(self, sym, w, nu, place, status):
        self.sym = sym
        self.name = str(sym)
        self.w = w
        self.nu = nu
        self.place = place
        self.status = status        # required-nonzero | forbidden | optional | param
        self.block = (w, nu)


class Eq:
    __slots__ = ("label", "w", "nu", "expr", "vs", "block")

    def __init__(self, label, w, nu, expr, vs):
        self.label = label          # (generator, y-order)
        self.w = w
        self.nu = nu
        self.expr = expr            # sympy expr in unknown symbols, = 0
        self.vs = vs                # frozenset of unknown symbols present
        self.block = (w, nu)


class System:
    def __init__(self, name, PQ, notes=""):
        self.name = name
        self.P, self.Q = PQ         # delta = Q*nu - P*w  (>=0 on/above extreme ray)
        self.vars = {}              # sym -> Var  (the scalar UNKNOWNS)
        self.params = {}            # sym -> Var  (nonzero state params / consts)
        self.eqs = []
        self.notes = notes
        self.contract = None        # optional endpoint contract (required_nonzero list)
        self.grade_axis = "both"    # "both" or "nu" (which axis the consistency check asserts)
        self.extra = {}             # per-system computed extras

    def delta(self, w, nu):
        return self.Q * nu - self.P * w

    def add_var(self, sym, w, nu, place, status):
        self.vars[sym] = Var(sym, w, nu, place, status)

    def add_param(self, sym, place, w=0, nu=0):
        self.params[sym] = Var(sym, w, nu, place, "param")

    def unknown_syms(self):
        return list(self.vars.keys())


# =====================================================================
#  1.  Bidegree annotation + monomial consistency
# =====================================================================
def _term_bigrade(term, grade_map):
    """(w, nu) of one monomial: sum of graded-symbol bigrades (with multiplicity)."""
    w = 0
    nu = 0
    saw = False
    for factor in sp.Mul.make_args(term):
        base, exp = factor.as_base_exp()
        if base in grade_map:
            v = grade_map[base]
            w += v.w * int(exp)
            nu += v.nu * int(exp)
            saw = True
    return w, nu, saw


# window-symbol u-weights  w(d_{4-k}) = k,  w(Phi) = 17   [(72,108) grading]
SYMBOL_WEIGHTS = None   # filled after _gsystem_symbols is defined


def check_bigrade_consistency(system):
    """Verify every monomial of every equation is bidegree-consistent.

    Two modes:
      * scalar mode (R1): every scalar coefficient variable AND every graded
        forcing/cap parameter carries a bidegree; each monomial of each scalar
        equation must sum to the equation's bidegree (multiplication adds
        degrees).  Asserts the nu (y-order) axis, the decisive grading.
      * symbol mode (R2/R3): the scalar equations are y-coefficients of a
        state-specialised generator, so the honest 'multiplication adds degrees'
        invariant is the u-weight homogeneity of the underlying window-symbol
        generator (checked pre-substitution), PLUS every emitted scalar equation
        must lie in the window cone delta = nu - 12*w >= 0.

    Returns (violations, checked_terms).  A LOUD finding if non-empty.
    """
    if system.extra.get("symbol_check"):
        return _symbol_homogeneity_and_cone(system)
    grade_map = {**system.vars, **system.params}
    violations = []
    checked = 0
    axis = system.grade_axis
    for eq in system.eqs:
        poly = sp.expand(eq.expr)
        for term in sp.Add.make_args(poly):
            w, nu, saw = _term_bigrade(term, grade_map)
            if not saw:
                w, nu = eq.w, eq.nu     # pure numeric inhomogeneity carries eq bidegree
            checked += 1
            ok = (nu == eq.nu) if axis == "nu" else (w == eq.w and nu == eq.nu)
            if not ok:
                violations.append((eq.label, str(term)[:60], (w, nu), (eq.w, eq.nu)))
    return violations, checked


def _symbol_homogeneity_and_cone(system):
    """R2/R3: u-weight homogeneity of each window-symbol generator + delta>=0 cone."""
    violations = []
    checked = 0
    for name, expr, w_G in system.extra["symbol_check"]:
        poly = sp.expand(expr)
        for term in sp.Add.make_args(poly):
            tw = 0
            for factor in sp.Mul.make_args(term):
                base, exp = factor.as_base_exp()
                if base in SYMBOL_WEIGHTS:
                    tw += SYMBOL_WEIGHTS[base] * int(exp)
            checked += 1
            if tw != w_G:
                violations.append(("SYMBOL", name, str(term)[:50], tw, w_G))
    # window-cone check on the emitted scalar equations
    for eq in system.eqs:
        checked += 1
        if system.delta(eq.w, eq.nu) < 0:
            violations.append(("CONE", eq.label, "delta<0", system.delta(eq.w, eq.nu)))
    return violations, checked


# =====================================================================
#  2.  Blocks, incidence, components, DM/SCC, rank, faces, nullspace
# =====================================================================
def variable_blocks(system):
    blocks = defaultdict(list)
    for v in system.vars.values():
        blocks[v.block].append(v.name)
    return blocks


def equation_blocks(system):
    blocks = defaultdict(list)
    for e in system.eqs:
        blocks[e.block].append(e.label)
    return blocks


class UnionFind:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def connected_components(system):
    """Union-find over the bipartite (equation, variable) incidence.

    Returns list of components, each = (set eq labels, set var names)."""
    uf = UnionFind()
    for e in system.eqs:
        etag = ("E", e.label)
        uf.find(etag)
        for s in e.vs:
            uf.union(etag, ("V", str(s)))
    comps = defaultdict(lambda: (set(), set()))
    for e in system.eqs:
        root = uf.find(("E", e.label))
        comps[root][0].add(e.label)
    for e in system.eqs:
        root = uf.find(("E", e.label))
        for s in e.vs:
            comps[root][1].add(str(s))
    # variables never appearing form singleton comps (ignored: not coupled)
    return [comps[k] for k in comps]


def _kuhn_matching(eqs, adj):
    """Maximum bipartite matching eq -> var (adj: eq_index -> set(var))."""
    match_var = {}      # var -> eq_index

    def try_aug(ei, seen):
        for v in adj[ei]:
            if v in seen:
                continue
            seen.add(v)
            if v not in match_var or try_aug(match_var[v], seen):
                match_var[v] = ei
                return True
        return False

    for ei in range(len(eqs)):
        try_aug(ei, set())
    match_eq = {ei: v for v, ei in match_var.items()}
    return match_eq, match_var


def _tarjan_scc(nodes, out_edges):
    index = {}
    low = {}
    onstack = {}
    stack = []
    result = []
    counter = [0]

    import sys as _sys
    _sys.setrecursionlimit(100000)

    def strong(v):
        index[v] = low[v] = counter[0]
        counter[0] += 1
        stack.append(v)
        onstack[v] = True
        for w in out_edges.get(v, ()):  # noqa
            if w not in index:
                strong(w)
                low[v] = min(low[v], low[w])
            elif onstack.get(w):
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                onstack[w] = False
                comp.append(w)
                if w == v:
                    break
            result.append(comp)

    for v in nodes:
        if v not in index:
            strong(v)
    return result       # reverse-topological order


def block_triangular(system):
    """Dulmage-Mendelsohn-style fine decomposition.

    Match equations to variables, orient var_i -> var_j when the equation matched
    to var_i also contains var_j, then Tarjan SCC => block-triangular blocks.
    Terminal blocks = sink SCCs (no edge leaving the SCC).  Returns dict of
    metrics + the SCC list (over matched variables)."""
    eqs = system.eqs
    var_names = [str(s) for s in system.vars]
    name_set = set(var_names)
    adj = []
    for e in eqs:
        adj.append({str(s) for s in e.vs if str(s) in name_set})
    match_eq, match_var = _kuhn_matching(eqs, adj)   # eq_index -> var, var -> eq_index

    # orient over matched variables
    matched_vars = list(match_var.keys())
    out = defaultdict(set)
    for v in matched_vars:
        ei = match_var[v]
        for w in adj[ei]:
            if w in match_var and w != v:
                out[v].add(w)
    sccs = _tarjan_scc(matched_vars, out)

    # map var -> scc id
    scc_of = {}
    for i, comp in enumerate(sccs):
        for v in comp:
            scc_of[v] = i
    # sink SCCs: no edge from comp to a *different* comp
    terminal = []
    for i, comp in enumerate(sccs):
        cs = set(comp)
        is_sink = True
        for v in comp:
            for w in out[v]:
                if scc_of.get(w) != i:
                    is_sink = False
                    break
            if not is_sink:
                break
        if is_sink:
            terminal.append(comp)

    sizes = sorted((len(c) for c in sccs), reverse=True)
    return {
        "n_matched": len(matched_vars),
        "n_unmatched_vars": len(name_set) - len(matched_vars),
        "n_unmatched_eqs": len(eqs) - len(match_eq),
        "scc_sizes": sizes,
        "largest_scc": sizes[0] if sizes else 0,
        "n_scc": len(sccs),
        "n_terminal": len(terminal),
        "terminal_sizes": sorted((len(c) for c in terminal), reverse=True),
        "sccs": sccs,
        "terminal": terminal,
    }


def generic_modular_rank(eqs, unknown_syms, param_syms, p=32003, seed=1):
    """Rank of the Jacobian d(eqs)/d(unknowns) at a random point mod p."""
    if not eqs or not unknown_syms:
        return 0, len(eqs), len(unknown_syms)
    rng = random.Random(seed)
    allsyms = list(unknown_syms) + list(param_syms)
    pt = {s: sp.Integer(rng.randrange(1, p)) for s in allsyms}
    rows = []
    for e in eqs:
        row = []
        expr = e.expr if hasattr(e, "expr") else e
        for s in unknown_syms:
            d = sp.diff(expr, s)
            row.append(0 if d == 0 else _rat_mod(d.xreplace(pt), p))
        rows.append(row)
    return _rank_mod(rows, p), len(eqs), len(unknown_syms)


def _rat_mod(expr, p):
    """Reduce a rational-number sympy expression to an int in F_p."""
    val = sp.nsimplify(sp.together(expr))
    num, den = sp.fraction(sp.Rational(val) if val.is_number else val)
    num, den = int(num), int(den)
    return (num % p) * pow(den % p, p - 2, p) % p


def _rank_mod(rows, p):
    M = [r[:] for r in rows]
    nrows = len(M)
    ncols = len(M[0]) if M else 0
    rank = 0
    col = 0
    for col in range(ncols):
        piv = None
        for i in range(rank, nrows):
            if M[i][col] % p != 0:
                piv = i
                break
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = pow(M[rank][col], p - 2, p)
        M[rank] = [(x * inv) % p for x in M[rank]]
        for i in range(nrows):
            if i != rank and M[i][col] % p != 0:
                f = M[i][col]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rank])]
        rank += 1
        if rank == nrows:
            break
    return rank


def extremal_face(system, functional, which="min"):
    """Extract the exposed face: the equations & variables minimising (or
    maximising) a linear functional f(w, nu) over the bigrade support."""
    vals_e = [(functional(e.w, e.nu), e) for e in system.eqs]
    if not vals_e:
        return None
    ext_e = (min if which == "min" else max)(v for v, _ in vals_e)
    face_e = [e for val, e in vals_e if val == ext_e]
    # the variables of the face are exactly the unknowns occurring in its equations
    face_v = set()
    for e in face_e:
        face_v |= {str(s) for s in e.vs}
    return {"value": ext_e, "eqs": face_e, "vars": sorted(face_v)}


def left_nullspace_certificate(face_eqs, unknown_syms, param_syms, p=32003, seed=3):
    """On a small block, test whether a left-null covector of the Jacobian hits a
    nonzero inhomogeneity -> a face-local kill certificate candidate.

    Returns dict {kill: bool, reason, n_eq, n_unk, rank}."""
    unk = [s for s in unknown_syms
           if any(e.expr.has(s) for e in face_eqs)]
    n_eq = len(face_eqs)
    if n_eq == 0:
        return {"kill": False, "reason": "empty face"}
    # inhomogeneous part: a face equation with NO unknowns at all whose value is a
    # forced-nonzero constant/parameter-expression => a window-depth kill certificate
    # (a required-nonzero forcing term stranded below the spare-window reach).
    for e in face_eqs:
        if not any(e.expr.has(s) for s in unk):
            val = sp.simplify(e.expr)
            if val != 0:
                nonzero = "nonzero constant" if not val.free_symbols else "generically-nonzero parameter form"
                return {"kill": True,
                        "reason": "inhomogeneous face equation, 0 unknowns (forced %s = 0)" % nonzero,
                        "eq": str(e.label), "value": str(val)[:48],
                        "n_eq": n_eq, "n_unk": 0, "rank": 0}
    # otherwise: is the face over-determined (more independent eqs than unknowns)?
    rank, ne, nu = generic_modular_rank(face_eqs, unk, param_syms, p=p, seed=seed)
    over = rank > len(unk)
    return {"kill": False, "reason": ("over-determined face (rank>unknowns)" if over
                                      else "consistent face (rank<=unknowns)"),
            "n_eq": n_eq, "n_unk": len(unk), "rank": rank, "over_determined": over}


# =====================================================================
#  3.  Reporting
# =====================================================================
def analyze(system, verbose=True):
    out = {}
    viol, checked = check_bigrade_consistency(system)
    out["bigrade_terms_checked"] = checked
    out["bigrade_violations"] = viol

    vblocks = variable_blocks(system)
    eblocks = equation_blocks(system)
    out["n_var_blocks"] = len(vblocks)
    out["n_eq_blocks"] = len(eblocks)
    out["n_vars"] = len(system.vars)
    out["n_eqs"] = len(system.eqs)

    comps = connected_components(system)
    csizes = sorted(((len(e), len(v)) for e, v in comps),
                    key=lambda t: (t[0] + t[1]), reverse=True)
    out["n_components"] = len(comps)
    out["component_sizes"] = csizes[:12]

    bt = block_triangular(system)
    out["block_triangular"] = {k: bt[k] for k in
                               ("n_matched", "n_unmatched_vars", "n_unmatched_eqs",
                                "largest_scc", "n_scc", "n_terminal", "terminal_sizes")}
    out["scc_size_hist"] = _hist(bt["scc_sizes"])

    rank, ne, nunk = generic_modular_rank(system.eqs, system.unknown_syms(),
                                          list(system.params))
    out["global_rank"] = rank
    out["global_rank_deficiency"] = min(ne, nunk) - rank

    # extremal faces
    faces = {}
    P, Q = system.P, system.Q
    functionals = {
        "delta(min)": (lambda w, nu: Q * nu - P * w, "min"),
        "delta(max)": (lambda w, nu: Q * nu - P * w, "max"),
        "w(min)": (lambda w, nu: w, "min"),
        "nu(max)": (lambda w, nu: nu, "max"),
        "mix 2w+3nu(min)": (lambda w, nu: 2 * w + 3 * nu, "min"),
    }
    for fname, (f, which) in functionals.items():
        fc = extremal_face(system, f, which)
        if fc is None:
            continue
        cert = left_nullspace_certificate(fc["eqs"], system.unknown_syms(),
                                          list(system.params))
        faces[fname] = {"value": fc["value"], "n_eqs": len(fc["eqs"]),
                        "n_vars": len(fc["vars"]), "certificate": cert}
    out["faces"] = faces

    if verbose:
        _print_report(system, out, bt)
    return out, bt


def _hist(sizes):
    h = defaultdict(int)
    for s in sizes:
        h[s] += 1
    return dict(sorted(h.items()))


def _print_report(system, out, bt):
    print("=" * 78)
    print("SYSTEM %s   (delta = %d*nu - %d*w)" % (system.name, system.Q, system.P))
    if system.notes:
        print("  " + system.notes)
    print("-" * 78)
    v = out["bigrade_violations"]
    print("  bigrade consistency: %d monomials checked, %d violations  %s"
          % (out["bigrade_terms_checked"], len(v), "OK" if not v else "*** LOUD FINDING ***"))
    for row in v[:6]:
        print("     VIOLATION", row)
    print("  scalars: %d unknown vars in %d bidegree-blocks; %d equations in %d blocks"
          % (out["n_vars"], out["n_var_blocks"], out["n_eqs"], out["n_eq_blocks"]))
    print("  connected components (bipartite incidence): %d" % out["n_components"])
    print("     largest (|eq|,|var|):", out["component_sizes"][:6])
    b = out["block_triangular"]
    print("  block-triangular / SCC: %d SCCs, largest=%d, terminal(sink)=%d %s"
          % (b["n_scc"], b["largest_scc"], b["n_terminal"], b["terminal_sizes"][:6]))
    print("     matched %d, unmatched vars %d, unmatched eqs %d"
          % (b["n_matched"], b["n_unmatched_vars"], b["n_unmatched_eqs"]))
    print("     SCC size histogram {size:count}:", out["scc_size_hist"])
    print("  generic modular rank: %d  (deficiency %d)"
          % (out["global_rank"], out["global_rank_deficiency"]))
    print("  extremal faces:")
    for fn, fc in out["faces"].items():
        c = fc["certificate"]
        tag = "KILL-CERT" if c.get("kill") else c.get("reason", "")
        print("     %-16s value=%-5s |eq|=%-3d |var|=%-3d  -> %s"
              % (fn, fc["value"], fc["n_eqs"], fc["n_vars"], tag))
    print()


# =====================================================================
#  4.  System builders
# =====================================================================
def _collect_y_equations(expr, weight_w, gen_name, unknown_syms, param_syms):
    """Expand a stripped generator in y; each y-coefficient is one equation with
    bidegree (w = weight_w, nu = 12*weight_w + j) [(72,108) grading]."""
    poly = sp.Poly(sp.expand(expr), y)
    eqs = []
    deg = poly.degree()
    for j in range(deg + 1):
        c = poly.nth(j)
        if c == 0:
            continue
        vs = frozenset(s for s in unknown_syms if c.has(s))
        eqs.append(Eq((gen_name, 12 * weight_w + j), weight_w, 12 * weight_w + j, c, vs))
    return eqs


# ---- R1 : the (50,75) window system (regression) --------------------
# Endpoint contract data (read-only, from ENDPOINT_CONTRACT.md sec.3); the
# required-nonzero list is contract INPUT, the forced_floor is COMPUTED here.
ENDPOINT_CONTRACT_50_75 = {
    "case": "F2_j0_50_75",
    "chart": "gamma=3 reduced (GGV3 sec.5)",
    "required_nonzero": [(-1, 3), (-2, 4), (0, -10)],   # (series, y_order)
    "caps": {(-1, 3): "a", (-2, 4): "b"},
}


def build_R1():
    """The (50,75) gamma=3 window system, at the scalar (series, y-order) level.

    The forced C_0 support (forced_floor) is COMPUTED by eliminating the linear
    window block (the coefficient-matching of the forcing relation), NOT read from
    the contract.  The required-nonzero list is the contract's data.
    """
    sysR = System("R1  (50,75) window", (25, 7),
                  notes="a=2 F2 chart; kill = y-order window-depth (F2_TOWER sec.1). "
                        "Decisive grading is nu (y-order); w = series index (block).")
    sysR.grade_axis = "nu"
    a, b, lam, f2, f4, f6, f8 = sp.symbols("a b lam f2 f4 f6 f8")
    # graded forcing/cap parameters (nu = the y-order each realises in the chart):
    #   a = C_{-1} lead @ y^3, b = C_{-2} lead @ y^4, lam @ y^0, f_{2k} @ y^{2k}.
    param_nu = {a: 3, b: 4, lam: 0, f2: 2, f4: 4, f6: 6, f8: 8}
    param_w = {a: -1, b: -2, lam: 0, f2: 0, f4: 0, f6: 0, f8: 0}
    for s in (a, b, lam, f2, f4, f6, f8):
        sysR.add_param(s, "forcing/cap", w=param_w[s], nu=param_nu[s])

    # C_0 = sum c0[j] y^j, unknown coefficients over a wide y-order window.
    JMIN, JMAX = -12, 2
    c0 = {j: sp.Symbol("c0_%d" % j if j >= 0 else "c0_m%d" % (-j)) for j in range(JMIN, JMAX + 1)}
    for j, sym in c0.items():
        status = "required-nonzero" if (0, j) in ENDPOINT_CONTRACT_50_75["required_nonzero"] else "optional"
        sysR.add_var(sym, 0, j, "corner series C_0", status)

    # The forcing relation (F2_TOWER sec.1, f2_tower.py a2_certificate, exact):
    #   3*(a y^3)^2 * C_0  =  3*(b y^4)^2 + 2*lam + 2*(f8 y^8+f6 y^6+f4 y^4+f2 y^2)
    Fm2 = f8 * y**8 + f6 * y**6 + f4 * y**4 + f2 * y**2
    C0series = sum(c0[j] * y**j for j in range(JMIN, JMAX + 1))
    lhs = sp.expand(3 * (a * y**3)**2 * C0series)
    rhs = sp.expand(3 * (b * y**4)**2 + 2 * lam + 2 * Fm2)
    relation = sp.expand(lhs - rhs)          # a Laurent poly in y; each power = 0

    # scalarize: multiply by y^12 to clear negatives, then take y-coefficients.
    cleared = sp.expand(relation * y**12)
    poly = sp.Poly(cleared, y)
    for pw in range(poly.degree() + 1):
        c = poly.nth(pw)
        if c == 0:
            continue
        yorder = pw - 12
        vs = frozenset(s for s in c0.values() if c.has(s))
        sysR.eqs.append(Eq(("forcing", yorder), 0, yorder, c, vs))

    sysR.contract = ENDPOINT_CONTRACT_50_75
    sysR.extra["c0"] = c0
    sysR.extra["relation"] = relation
    return sysR


def r1_forced_floor(sysR):
    """COMPUTE forced_floor[0] by solving the linear window block for c0[j].

    A coefficient c0[j] is forced-zero iff its solved value is identically 0
    (no free forcing parameter can make it nonzero)."""
    c0 = sysR.extra["c0"]
    a, b, lam, f2, f4, f6, f8 = sp.symbols("a b lam f2 f4 f6 f8")
    # linear system: solve for each c0[j].
    sol = sp.solve([e.expr for e in sysR.eqs], list(c0.values()), dict=True)
    assert sol, "window block did not solve"
    sol = sol[0]
    nonzero_orders = []
    for j, sym in c0.items():
        val = sp.simplify(sol.get(sym, sym))
        if val != 0:
            nonzero_orders.append(j)
    floor = min(nonzero_orders)
    return floor, sorted(nonzero_orders)


def r1_rediscover(sysR):
    """The load-bearing regression: scan the contract's required-nonzero list
    against the COMPUTED forced_floor and localise the kill automatically."""
    floor, support = r1_forced_floor(sysR)
    kills = []
    for (s, j) in sysR.contract["required_nonzero"]:
        if s == 0 and j < floor:            # required-nonzero coeff below forced floor
            kills.append((s, j))
    return {"forced_floor_series0": floor, "forced_support_series0": support,
            "required_nonzero": sysR.contract["required_nonzero"], "kills": kills}


# ---- (72,108) generators (read-only, verbatim) ----------------------
def _gsystem_symbols():
    return sp.symbols("d0 d1 d2 dm1 dm2 dm3 dm4 Phi")


_d0, _d1, _d2, _dm1, _dm2, _dm3, _dm4, _Phi = _gsystem_symbols()
SYMBOL_WEIGHTS = {_d0: 4, _d1: 3, _d2: 2, _dm1: 5, _dm2: 6, _dm3: 7, _dm4: 8, _Phi: 17}


def _G_generators():
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = _gsystem_symbols()
    R = sp.Rational
    G1 = R(3, 2) * d1 * dm1**2 + 3 * d2 * dm1 * dm2 + 3 * dm1 * dm4 + 3 * dm2 * dm3
    G2 = -R(3, 2) * d0 * dm1**2 + R(3, 2) * d2 * dm2**2 + 3 * dm2 * dm4 + R(3, 2) * dm3**2
    G3 = -3 * d0 * dm1 * dm2 - R(3, 2) * d1 * dm2**2 - R(1, 2) * dm1**3 + 3 * dm3 * dm4
    G5body = (-3 * d0 * dm1 * dm4 - 3 * d0 * dm2 * dm3 - 3 * d1 * dm2 * dm4
              - R(3, 2) * d1 * dm3**2 - 3 * d2 * dm3 * dm4 - R(3, 2) * dm1**2 * dm3
              - R(3, 2) * dm1 * dm2**2)
    # CANONICAL normalisation: G5 = G5body + Phi.  (Corrected 2026-07-24; this
    # read `2 * Phi + G5body`, transcribed from FULL_SYSTEM_BRIDGE.md:62, which
    # contradicted line 114 of that file and the canonical loader
    # full_system_bridge.py.  The authority is the C11 membership certificate in
    # f37_sat_verify.py: f31 == c1*G1 + c2*G2 + c3*G3 + c4*(G5body + Phi).
    # The forms differ by Phi, not by a scalar -- see FACE_KILL_SWEEP.md sec.4.)
    G5 = Phi + G5body
    return {"G1": (G1, 13), "G2": (G2, 14), "G3": (G3, 15), "G5": (G5, 17)}


def _H_generators():
    """Load the dm4-free H-system verbatim from r9_eliminated_system.json."""
    with open(os.path.join(HERE, "r9_eliminated_system.json")) as fh:
        d = json.load(fh)
    syms = dict(zip(("d0", "d1", "d2", "dm1", "dm2", "dm3", "dm4", "Phi"),
                    _gsystem_symbols()))
    local = {k: v for k, v in syms.items()}
    weight = {"H2": 19, "H3": 20, "H5": 22}   # 228/240/264 = 12*w
    Hs = {}
    for name, expr_str in d["H"].items():
        Hs[name] = (sp.sympify(expr_str, locals=local), weight[name])
    return Hs


def _phi_stripped():
    c = sp.Rational(-1, 6630)
    t = y + 1
    q = 2048 * y**4 - 512 * y**3 + 320 * y**2 - 240 * y + 195
    return sp.expand(c * t**30 * q)          # deg 34, ord 0


def _place_tag_72(w):
    return "class m=%d (Z/12); delta-floor nu>=%d" % (w % 12, 12 * w)


# ---- R2 : fast home-case control kill (full G-system, a8 pilot) ------
def build_R2():
    """The a8 constant-E pilot state (FULL_SYSTEM_BRIDGE sec.5, sub2), a documented
    exact-Q UNIT home-case kill.  Full G-system G1,G2,G3,G5=2Phi+G5body."""
    sysR = System("R2  a8 home-case kill (full G-system)", (12, 1),
                  notes="(72,108) home case; e=g*(y+1)^8, d2=0, d1 deg2, sigma deg7, "
                        "sub2 spares dm2/dm3/dm4 deg 12/14/16. Documented exact-Q UNIT.")
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = _gsystem_symbols()
    gamma = sp.Symbol("gamma")
    b = sp.symbols("b0 b1 b2")
    ssig = sp.symbols("s0 s1 s2 s3 s4 s5 s6 s7")
    for s in (gamma,) + b + ssig:
        sysR.add_param(s, "state param")

    d2v = sp.Integer(0)
    d1v = b[0] + b[1] * y + b[2] * y**2
    ev = gamma * (y + 1)**8
    sigma = sum(ssig[i] * y**i for i in range(8))
    d0v = sp.expand((d2v**2 + sigma) / 4)

    # sub2 spare ansaetze (stripped): deg caps 2k = 12,14,16 for k=6,7,8
    unknown = []
    spare_specs = [("R", 6, 12), ("S", 7, 14), ("T", 8, 16)]
    spare_series = {}
    for prefix, k, cap in spare_specs:
        coeffs = sp.symbols("%s0:%d" % (prefix, cap + 1))
        spare_series[k] = sum(coeffs[i] * y**i for i in range(cap + 1))
        for i, sym in enumerate(coeffs):
            sysR.add_var(sym, k, 12 * k + i, _place_tag_72(k), "optional")
            unknown.append(sym)

    subs = {d0: d0v, d1: d1v, d2: d2v, dm1: ev,
            dm2: spare_series[6], dm3: spare_series[7], dm4: spare_series[8],
            Phi: _phi_stripped()}
    gens = _G_generators()
    sysR.extra["symbol_check"] = [(name, gexpr, w) for name, (gexpr, w) in gens.items()]
    for name, (gexpr, w) in gens.items():
        inst = sp.expand(gexpr.xreplace(subs))
        sysR.eqs.extend(_collect_y_equations(inst, w, name, unknown, list(sysR.params)))
    return sysR


# ---- R3 : the walled R9 z=1 H-system --------------------------------
def build_R3():
    """R9 z=1 dm4-eliminated H-system (convolution_elim_qsupport.build_qsupport_ansatz(1)),
    reconstructed independently: d1=0, d2 deg4, e=g*(y+1)^9*(y-r), sigma=(y-r)^2*(g0+g1 y),
    d0=(d2^2+sigma)/4, spares dm2 deg12 / dm3 deg14 (28 spare), marked root q(r)=0."""
    sysR = System("R3  R9 z=1 H-system (walled)", (12, 1),
                  notes="dm4-free H2,H3,H5 (u-weight 19/20/22) on the R9 z=1 state; "
                        "marked root r in Q[r]/(q); the documented COST wall.")
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = _gsystem_symbols()
    gamma, r = sp.symbols("gamma r")
    a = sp.symbols("a0 a1 a2 a3 a4")           # d2 coeffs, deg 4
    g = sp.symbols("g0 g1")                    # G, deg z=1
    for s in (gamma, r) + a + g:
        sysR.add_param(s, "state param (root algebra Q[r]/q, deg 4)")

    d2v = sum(a[i] * y**i for i in range(5))
    d1v = sp.Integer(0)
    ev = gamma * (y + 1)**9 * (y - r)
    sigma = sp.expand((y - r)**2 * (g[0] + g[1] * y))
    d0v = sp.expand((d2v**2 + sigma) / 4)

    unknown = []
    spare_series = {}
    for prefix, k, cap in [("R", 6, 12), ("S", 7, 14)]:
        coeffs = sp.symbols("%s0:%d" % (prefix, cap + 1))
        spare_series[k] = sum(coeffs[i] * y**i for i in range(cap + 1))
        for i, sym in enumerate(coeffs):
            sysR.add_var(sym, k, 12 * k + i, _place_tag_72(k), "optional")
            unknown.append(sym)

    subs = {d0: d0v, d1: d1v, d2: d2v, dm1: ev,
            dm2: spare_series[6], dm3: spare_series[7], Phi: _phi_stripped()}
    hgens = _H_generators()
    sysR.extra["symbol_check"] = [(name, hexpr, w) for name, (hexpr, w) in hgens.items()]
    for name, (hexpr, w) in hgens.items():
        inst = sp.expand(hexpr.xreplace(subs))
        sysR.eqs.extend(_collect_y_equations(inst, w, name, unknown, list(sysR.params)))

    sysR.extra["spare_series"] = spare_series
    sysR.extra["monic_e"] = sp.expand((y + 1)**9 * (y - r))
    return sysR


def r3_band_structure(sysR):
    """Re-derive the TRUE (u,y) band structure of the R9 z=1 H-system (adjudication).

    - u-axis occupancy (distinct generator weights)
    - per y-order slice: the spare unknowns present, and how many are NEW as the
      slice index sweeps upward from 0 (the staircase increment)
    - the extreme corner: the top-y-order slice of each generator (eqs vs unknowns)
    """
    by_slice = defaultdict(set)               # nu -> set of spare names
    weights = set()
    per_gen_top = {}
    gen_ranges = defaultdict(lambda: [10**9, -10**9])
    for e in sysR.eqs:
        weights.add(e.w)
        spares = {s for s in map(str, e.vs)}
        by_slice[e.nu] |= spares
        gname = e.label[0]
        lo, hi = gen_ranges[gname]
        gen_ranges[gname] = [min(lo, e.nu), max(hi, e.nu)]

    # staircase: sweep nu ascending, count newly-introduced spares per occupied slice
    seen = set()
    staircase = []
    for nu in sorted(by_slice):
        new = by_slice[nu] - seen
        seen |= by_slice[nu]
        staircase.append((nu, len(by_slice[nu]), len(new), sorted(new)))

    # top corner per generator: highest-nu equation(s) and their unknowns
    corners = {}
    for gname, (lo, hi) in gen_ranges.items():
        top_eqs = [e for e in sysR.eqs if e.label[0] == gname and e.nu == hi]
        unk = set()
        for e in top_eqs:
            unk |= {s for s in map(str, e.vs)}
        corners[gname] = {"nu": hi, "n_eq": len(top_eqs), "unknowns": sorted(unk)}

    # merged extreme corner: the global top-nu face
    top_nu = max(by_slice)
    top_eqs = [e for e in sysR.eqs if e.nu == top_nu]
    top_unk = set()
    for e in top_eqs:
        top_unk |= {s for s in map(str, e.vs)}
    return {
        "u_axis_values": sorted(weights),
        "gen_yorder_ranges": {g: gen_ranges[g] for g in gen_ranges},
        "staircase_head": staircase[:6],
        "staircase_tail": staircase[-4:],
        "increments_seen": sorted(set(x[2] for x in staircase)),
        "corners_per_gen": corners,
        "global_top_corner": {"nu": top_nu, "n_eq": len(top_eqs),
                              "unknowns": sorted(top_unk)},
    }


# =====================================================================
#  5.  Main
# =====================================================================
def main():
    print("#" * 78)
    print("# BIGRADE ANNOTATOR  --  Milestone 1 (read-only diagnosis)")
    print("#" * 78)

    # ---- R1 : the decisive regression -------------------------------
    r1 = build_R1()
    out1, bt1 = analyze(r1)
    red = r1_rediscover(r1)
    print(">>> R1 REGRESSION -- automatic rediscovery from contract + block structure")
    print("    COMPUTED forced_floor[series 0] =", red["forced_floor_series0"],
          " (support", red["forced_support_series0"], ")")
    print("    contract required-nonzero       =", red["required_nonzero"])
    print("    KILLS (required-nonzero below computed floor):", red["kills"])
    assert red["kills"] == [(0, -10)], "R1 rediscovery FAILED"
    print("    => KILL localised at c_{0,-10}: required-nonzero coefficient forced")
    print("       into a forbidden slot (-10 < floor -6).  REGRESSION PASSED.\n")

    # ---- R2 : fast home-case control kill ---------------------------
    r2 = build_R2()
    out2, bt2 = analyze(r2)
    print(">>> R2 LOCALIZATION -- does the contradiction localise?")
    r2kill = None
    for fn, fc in out2["faces"].items():
        if fc["certificate"].get("kill"):
            r2kill = (fn, fc)
            break
    if r2kill:
        fn, fc = r2kill
        print("    localised: face %s -> %s (|eq|=%d, |var|=%d)"
              % (fn, fc["certificate"]["reason"], fc["n_eqs"], fc["n_vars"]))
    else:
        print("    no single-face inhomogeneous certificate; smallest terminal block:",
              out2["block_triangular"]["terminal_sizes"][:4])
    print()

    # ---- R3 : the walled R9 z=1 state -------------------------------
    r3 = build_R3()
    out3, bt3 = analyze(r3)
    band = r3_band_structure(r3)
    print(">>> R3 SIGNATURE + BAND-STRUCTURE ADJUDICATION")
    print("    u-axis values (distinct generator weights):", band["u_axis_values"])
    print("    generator y-order ranges:", band["gen_yorder_ranges"])
    print("    staircase (nu, #spares_on_slice, #NEW, new-names) head:")
    for row in band["staircase_head"]:
        print("       nu=%d  present=%d  new=%d  %s" % row)
    print("    per-slice NEW-spare increments observed:", band["increments_seen"])
    print("    top corner per generator:")
    for g, c in band["corners_per_gen"].items():
        print("       %s: nu=%d  %d eq  unknowns=%s" % (g, c["nu"], c["n_eq"], c["unknowns"]))
    print("    global top corner: nu=%d  %d eq  unknowns=%s"
          % (band["global_top_corner"]["nu"], band["global_top_corner"]["n_eq"],
             band["global_top_corner"]["unknowns"]))
    r3.extra["band"] = band

    print("\n" + "#" * 78)
    print("# comparison numbers (scalar-vs-bigraded):")
    for tag, out in (("R1", out1), ("R2", out2), ("R3", out3)):
        b = out["block_triangular"]
        print("#  %-4s vars=%-3d eqs=%-4d comps=%-3d largestSCC=%-3d terminals=%-3d rankdef=%d"
              % (tag, out["n_vars"], out["n_eqs"], out["n_components"],
                 b["largest_scc"], b["n_terminal"], out["global_rank_deficiency"]))
    print("#" * 78)
    return {"R1": (r1, out1, red), "R2": (r2, out2), "R3": (r3, out3, band)}


if __name__ == "__main__":
    main()
