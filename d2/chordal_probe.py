#!/usr/bin/env python3
"""chordal_probe.py -- NEW; READ-ONLY over every existing artifact.

STRUCTURAL DIAGNOSTIC: does any of our systems have exploitable CHORDAL
structure (small treewidth / small maximal cliques), which SCC connectivity
does not measure?

Motivation.  The sub2 G-system's bipartite incidence / Dulmage-Mendelsohn
analysis produced ONE strongly connected component containing all 45 spare
coefficients, so there is no topological elimination ordering.  But SCC
connectivity and TREEWIDTH are different invariants: a connected graph can
still admit a chordal completion whose maximal cliques are small, which is
exactly what Macaulay2's `Chordal` package (Cifuentes-Parrilo) exploits --
it builds a chordal network from a polynomial set and does elimination,
dimension, component and root counting along the clique tree, at a cost
exponential in TREEWIDTH rather than in the number of variables.

WHAT THIS SCRIPT MEASURES, for each target system:
  * the variable-interaction (constraint / primal) graph: vertices = ring
    variables, edge iff two variables co-occur in the support of some
    polynomial.  V, E, density, max/mean degree, #(complete-graph edges).
  * a treewidth UPPER BOUND from min-degree and min-fill elimination
    heuristics (implemented here; networkx only used for a cross-check when
    available).  Exact treewidth is NOT needed: an upper bound is the whole
    point -- a SMALL upper bound proves exploitable structure, and a LOWER
    bound (max clique of the raw graph, computed greedily/exactly for the
    small cases) certifies when no ordering can help.
  * the maximum clique size in the chordal completion produced by that
    ordering (= treewidth bound + 1).
  * a lower bound on treewidth: (a) the largest single-polynomial support
    (every polynomial's support is a clique, so tw >= maxsupport-1) and
    (b) the minimum-degree lower bound / degeneracy.

TARGETS
  J6      -- the depth-3 saturated system of the four ALT_HUNT [J6] states,
             read VERBATIM from `j6_msolve_results.json` (a KNOWN kill: msolve
             returns [-1] = empty at depth 3).  4-5 variables; pipeline
             validation.
  SUB2    -- the sub2 home G-system at full caps, built through
             `face_kill_sweep.build_state_system` (the audited constructor),
             with d2 deg 4, d1 deg 6, sigma deg 8, e = gamma*(y+1)^10.
             45 spare unknowns.
  R9      -- the R9 z=1 dm4-eliminated H-system, `bigrade_annotator.build_R3`,
             which loads `r9_eliminated_system.json` through `_H_generators`.
             28 spare unknowns.  The actual wall.

For SUB2 and R9 two variable models are reported, because the answer differs:
  spares   -- vertices = spare coefficient unknowns only (state parameters
              treated as generic coefficients).  This is the model in which
              "45 spare unknowns" is the variable count.
  full     -- vertices = spares + state parameters (the model an actual
              solver over Q would use).

Writes: chordal_probe.json (numbers), plus M2 input files under
        chordal_probe_m2/ .  Nothing else is touched.

Usage:
    python chordal_probe.py                 # all targets, graph numbers
    python chordal_probe.py --emit-m2       # also write the M2 scripts
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from collections import defaultdict

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
sys.path.insert(0, HERE)

OUT_JSON = os.path.join(HERE, "chordal_probe.json")
M2_DIR = os.path.join(HERE, "chordal_probe_m2")


# =====================================================================
#  1.  Graph machinery  (pure python; no dependency on networkx)
# =====================================================================
def build_graph(supports, vertices):
    """supports: iterable of iterables of vertex names.  Returns adjacency
    dict name -> set(names).  Every polynomial's support becomes a clique."""
    adj = {v: set() for v in vertices}
    for sup in supports:
        s = [v for v in sup if v in adj]
        for a, b in itertools.combinations(s, 2):
            adj[a].add(b)
            adj[b].add(a)
    return adj


def graph_stats(adj):
    V = len(adj)
    E = sum(len(s) for s in adj.values()) // 2
    maxE = V * (V - 1) // 2
    degs = sorted((len(s) for s in adj.values()), reverse=True)
    return {
        "vertices": V,
        "edges": E,
        "complete_graph_edges": maxE,
        "density": round(E / maxE, 6) if maxE else None,
        "max_degree": degs[0] if degs else 0,
        "min_degree": degs[-1] if degs else 0,
        "mean_degree": round(sum(degs) / V, 3) if V else 0,
        "isolated_vertices": sum(1 for d in degs if d == 0),
    }


def elimination_order(adj, rule="minfill"):
    """Greedy elimination-ordering heuristic.  Returns
    (order, width, fill_edges_added, clique_sizes).

    Eliminating v in the current (partially filled) graph forms the clique
    {v} u N(v); we record |N(v)| as the width contribution, add the fill
    edges among N(v), and delete v.  The maximum |N(v)| over the run is a
    TREEWIDTH UPPER BOUND and max(|N(v)|+1) is the largest maximal clique of
    the resulting chordal completion.

    rule = "mindeg"  : pick v minimising |N(v)|.
    rule = "minfill" : pick v minimising the number of fill edges it creates
                       (ties broken by degree).
    """
    g = {v: set(s) for v, s in adj.items()}
    order, width, fill_total, cliques = [], 0, 0, []
    remaining = set(g)
    while remaining:
        best, best_key = None, None
        for v in remaining:
            nb = g[v]
            d = len(nb)
            if rule == "mindeg":
                key = (d,)
            else:
                miss = 0
                nbl = sorted(nb)
                for i in range(len(nbl)):
                    ai = g[nbl[i]]
                    for j in range(i + 1, len(nbl)):
                        if nbl[j] not in ai:
                            miss += 1
                key = (miss, d)
            if best_key is None or key < best_key:
                best, best_key = v, key
        v = best
        nb = sorted(g[v])
        width = max(width, len(nb))
        cliques.append(len(nb) + 1)
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                if nb[j] not in g[nb[i]]:
                    g[nb[i]].add(nb[j])
                    g[nb[j]].add(nb[i])
                    fill_total += 1
        for u in nb:
            g[u].discard(v)
        del g[v]
        remaining.discard(v)
        order.append(v)
    return order, width, fill_total, cliques


def greedy_clique_lb(adj, tries=200, seed=0):
    """Cheap LOWER bound on the max clique of the RAW graph (hence on
    treewidth+1): repeated greedy growth from each vertex."""
    import random
    rng = random.Random(seed)
    verts = sorted(adj)
    best = 0
    best_set = []
    starts = verts if len(verts) <= tries else rng.sample(verts, tries)
    for v0 in starts:
        cand = set(adj[v0])
        clique = [v0]
        while cand:
            # pick the candidate with most connections inside cand
            u = max(cand, key=lambda x: len(adj[x] & cand))
            clique.append(u)
            cand &= adj[u]
        if len(clique) > best:
            best, best_set = len(clique), sorted(clique)
    return best, best_set


def degeneracy(adj):
    """Degeneracy = max over the min-degree peeling order.  It is a lower
    bound on treewidth."""
    g = {v: set(s) for v, s in adj.items()}
    k = 0
    while g:
        v = min(g, key=lambda x: len(g[x]))
        k = max(k, len(g[v]))
        for u in g[v]:
            g[u].discard(v)
        del g[v]
    return k


def analyse(name, supports, vertices, note=""):
    adj = build_graph(supports, vertices)
    st = graph_stats(adj)
    sup_sizes = sorted((len([v for v in s if v in adj]) for s in supports),
                       reverse=True)
    res = {
        "name": name,
        "note": note,
        "n_polynomials": len(supports),
        "max_support": sup_sizes[0] if sup_sizes else 0,
        "mean_support": round(sum(sup_sizes) / len(sup_sizes), 3) if sup_sizes else 0,
        "graph": st,
    }
    for rule in ("mindeg", "minfill"):
        t0 = time.monotonic()
        order, width, fill, cliques = elimination_order(adj, rule)
        res[rule] = {
            "treewidth_upper_bound": width,
            "max_clique_in_completion": max(cliques) if cliques else 0,
            "fill_edges_added": fill,
            "seconds": round(time.monotonic() - t0, 2),
        }
    lb, lbset = greedy_clique_lb(adj)
    res["clique_lower_bound_raw_graph"] = lb
    res["treewidth_lower_bound_from_clique"] = lb - 1
    res["treewidth_lower_bound_from_support"] = res["max_support"] - 1
    res["degeneracy"] = degeneracy(adj)
    res["networkx_crosscheck"] = _networkx_crosscheck(adj)
    res["support_histogram"] = _histogram(sup_sizes)
    res["layered"] = _layered(supports, vertices, sup_sizes)
    res["_adj"] = {k: sorted(v) for k, v in adj.items()}
    return res


def _histogram(sizes):
    h = defaultdict(int)
    for s in sizes:
        h[s] += 1
    return {str(k): h[k] for k in sorted(h)}


def _layered(supports, vertices, sup_sizes):
    """Is the density the fault of a handful of FAT equations?  Rebuild the
    graph from only the equations whose support is <= k, for a ladder of k,
    and report what treewidth bound survives.  If dropping the fat equations
    still leaves a near-complete graph, no sparsification can save the lane."""
    if not sup_sizes:
        return []
    out = []
    lo = sup_sizes[-1]
    hi = sup_sizes[0]
    ladder = sorted({lo, hi, max(lo, hi // 8), max(lo, hi // 4),
                     max(lo, hi // 2), max(lo, (3 * hi) // 4)})
    for k in ladder:
        keep = [s for s in supports if len(s) <= k]
        if not keep:
            continue
        adj = build_graph(keep, vertices)
        st = graph_stats(adj)
        _, w, _, _ = elimination_order(adj, "minfill")
        out.append({"support_cap": k, "n_polys_kept": len(keep),
                    "edges": st["edges"], "density": st["density"],
                    "isolated_vertices": st["isolated_vertices"],
                    "treewidth_upper_bound": w})
    return out


def _networkx_crosscheck(adj):
    """Independent treewidth upper bounds from networkx's own heuristics, and
    an EXACT max-clique of the raw graph when the graph is small enough."""
    try:
        import networkx as nx
        from networkx.algorithms.approximation import treewidth_min_degree, \
            treewidth_min_fill_in
    except Exception as exc:                       # pragma: no cover
        return {"available": False, "error": repr(exc)}
    G = nx.Graph()
    G.add_nodes_from(adj)
    for v, nb in adj.items():
        for u in nb:
            G.add_edge(v, u)
    md, _ = treewidth_min_degree(G)
    mf, _ = treewidth_min_fill_in(G)
    out = {"available": True, "nx_version": nx.__version__,
           "treewidth_min_degree": md, "treewidth_min_fill_in": mf}
    try:
        cl = nx.algorithms.clique.graph_clique_number(G)
        out["exact_max_clique"] = cl
    except Exception:
        best = 0
        for c in nx.find_cliques(G):
            best = max(best, len(c))
        out["exact_max_clique"] = best
    return out


# =====================================================================
#  2.  Target builders
# =====================================================================
def target_J6():
    """The depth-3 saturated systems of the four [J6] states, verbatim from
    j6_msolve_results.json (gens at all recorded degrees + class relations
    + saturation) -- i.e. exactly what msolve decided EMPTY."""
    path = os.path.join(HERE, "j6_msolve_results.json")
    d = json.load(open(path, encoding="utf-8"))
    out = []
    for r in d["results"]:
        assert r["verdict"] == "KILLED" and r["kill_depth"] == 3, r["key"]
        ring = r["ring_vars"]
        syms = {v: sp.Symbol(v) for v in ring}
        polys = []
        labels = []
        for g in r["gens"]:
            polys.append(sp.sympify(g["coefficient"], locals=syms))
            labels.append("mc_deg%d" % g["degree"])
        for i, s in enumerate(r["class_relations"]):
            polys.append(sp.sympify(s, locals=syms))
            labels.append("rel%d" % i)
        polys.append(sp.sympify(r["saturation"], locals=syms))
        labels.append("sat")
        supports = [sorted(str(s) for s in p.free_symbols & set(map(sp.Symbol, ring)))
                    for p in polys]
        out.append({
            "key": r["key"],
            "ring_vars": ring,
            "labels": labels,
            "polys": polys,
            "supports": supports,
        })
    return out


def _generic_state_polys():
    """The sub2 home G-system state asked for: d2 deg 4, d1 deg 6,
    sigma deg 8, e = gamma*(y+1)^10, all coefficients generic symbols."""
    y = sp.Symbol("y")
    a = sp.symbols("a0:5")
    b = sp.symbols("b0:7")
    s = sp.symbols("s0:9")
    gamma = sp.Symbol("gamma")
    return {
        "d2": sum(a[i] * y**i for i in range(5)),
        "d1": sum(b[i] * y**i for i in range(7)),
        "sigma": sum(s[i] * y**i for i in range(9)),
        "e": sp.expand(gamma * (y + 1)**10),
    }


def target_SUB2():
    import face_kill_sweep as fks
    polys = _generic_state_polys()
    sysR = fks.build_state_system(
        "CHORDAL-PROBE sub2 home G-system (full caps)", polys, "sub2",
        notes="d2 deg4, d1 deg6, sigma deg8, e=gamma*(y+1)^10; generic state.")
    return sysR


def target_R9():
    import bigrade_annotator as ba
    return ba.build_R3()


# ---------------------------------------------------------------------
#  Randomised-parameter specialisation (fast path).
#
#  The fully symbolic expansions of the R9 H-generators cost tens of minutes
#  in sympy.  For a SUPPORT question we do not need the symbolic coefficients:
#  specialising the state parameters at a random rational point makes `expand`
#  cheap and yields the GENERIC support -- a coefficient can only lose a spare
#  variable on a proper Zariski-closed subset of parameter space, so with
#  random parameters the support we see is, with probability 1, the support of
#  the symbolic system.  (Specialisation can only ever SHRINK a support, never
#  grow one, so the graph below is contained in the symbolic graph; if it is
#  already complete, so is the symbolic one -- no probabilistic caveat needed
#  for a NEGATIVE verdict.)
#
#  Cross-checked against the exact symbolic build on SUB2, where the exact run
#  is affordable.
# ---------------------------------------------------------------------
def _rand_rat(rng):
    return sp.Rational(rng.randrange(3, 5000), rng.randrange(1, 97))


def target_R9_generic(seed=20260724):
    """build_R3's construction with the state parameters specialised at a
    random rational point.  Uses bigrade_annotator's own _H_generators,
    _phi_stripped and _collect_y_equations -- only the parameter values
    differ from build_R3."""
    import random
    import bigrade_annotator as ba
    rng = random.Random(seed)
    yy = ba.y
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = ba._gsystem_symbols()

    gamma_v = _rand_rat(rng)
    r_v = _rand_rat(rng)
    a_v = [_rand_rat(rng) for _ in range(5)]      # d2 coeffs, deg 4
    g_v = [_rand_rat(rng) for _ in range(2)]      # G, deg z=1

    d2v = sum(a_v[i] * yy**i for i in range(5))
    d1v = sp.Integer(0)
    ev = sp.expand(gamma_v * (yy + 1)**9 * (yy - r_v))
    sigma = sp.expand((yy - r_v)**2 * (g_v[0] + g_v[1] * yy))
    d0v = sp.expand((d2v**2 + sigma) / 4)

    sysR = ba.System("R9 z=1 H-system (params specialised)", (12, 1),
                     notes="build_R3 construction, random rational state point")
    unknown, spare_series = [], {}
    for prefix, k, cap in [("R", 6, 12), ("S", 7, 14)]:
        coeffs = sp.symbols("%s0:%d" % (prefix, cap + 1))
        spare_series[k] = sum(coeffs[i] * yy**i for i in range(cap + 1))
        for i, sym in enumerate(coeffs):
            sysR.add_var(sym, k, 12 * k + i, ba._place_tag_72(k), "optional")
            unknown.append(sym)

    subs = {d0: d0v, d1: d1v, d2: d2v, dm1: ev,
            dm2: spare_series[6], dm3: spare_series[7],
            Phi: ba._phi_stripped()}
    for name, (hexpr, w) in ba._H_generators().items():
        inst = sp.expand(hexpr.xreplace(subs))
        sysR.eqs.extend(ba._collect_y_equations(inst, w, name, unknown, []))
    return sysR


def target_SUB2_generic(seed=20260724):
    """The same specialisation for the sub2 home G-system, used purely to
    CROSS-CHECK the fast path against the exact symbolic SUB2 build."""
    import random
    import face_kill_sweep as fks
    rng = random.Random(seed)
    yy = sp.Symbol("y")
    polys = {
        "d2": sum(_rand_rat(rng) * yy**i for i in range(5)),
        "d1": sum(_rand_rat(rng) * yy**i for i in range(7)),
        "sigma": sum(_rand_rat(rng) * yy**i for i in range(9)),
        "e": sp.expand(_rand_rat(rng) * (yy + 1)**10),
    }
    return fks.build_state_system(
        "CHORDAL-PROBE sub2 (params specialised)", polys, "sub2",
        notes="random rational state point; fast-path cross-check")


def target_SYMBOL():
    """The SYMBOL-LEVEL systems, BEFORE coefficient expansion: the G-system
    G1,G2,G3,G5 and the dm4-eliminated H-system H2,H3,H5 as polynomials in the
    eight series symbols d0,d1,d2,dm1,dm2,dm3,dm4,Phi.

    This is the third branch of the probe's interpretation: if the cliques are
    small only here, then the generic-K(y) / symbol-level representation is the
    right abstraction and the coefficient expansion is what destroys it."""
    import bigrade_annotator as ba
    syms = ba._gsystem_symbols()
    names = [str(s) for s in syms]
    out = []
    for tag, gens in (("G-system (G1,G2,G3,G5)", ba._G_generators()),
                      ("H-system (H2,H3,H5)", ba._H_generators())):
        sups = []
        for name, (expr, w) in sorted(gens.items()):
            sups.append(sorted({str(s) for s in expr.free_symbols} & set(names)))
        out.append((tag, sups, names,
                    {n: s for n, (s, _) in sorted(gens.items())}))
    return out


def system_supports(sysR, model):
    """model = 'spares' | 'full'."""
    spares = [str(v) for v in sysR.vars]
    params = [str(v) for v in sysR.params]
    verts = spares if model == "spares" else spares + params
    vset = set(verts)
    sups = []
    for e in sysR.eqs:
        fs = {str(s) for s in e.expr.free_symbols} & vset
        if fs:
            sups.append(sorted(fs))
    return verts, sups, spares, params


# =====================================================================
#  3.  M2 emission (Chordal package)
# =====================================================================
def m2_poly(expr, varmap):
    s = sp.sstr(sp.expand(expr))
    return s.replace("**", "^")


def emit_m2_j6(rec, path):
    ring = rec["ring_vars"]
    body = [
        'needsPackage "Chordal";',
        'R = QQ[%s, MonomialOrder=>Lex];' % ", ".join(ring),
        "I = ideal(",
    ]
    body.append(",\n".join("  " + m2_poly(p, ring) for p in rec["polys"]))
    body += [
        ");",
        'print("== ' + rec["key"] + '");',
        'print("nvars " | toString numgens R);',
        "G = constraintGraph I;",
        'print("constraintGraph " | toString G);',
        "CG = chordalGraph G;",
        'print("chordalGraph " | toString CG);',
        'print("treewidth " | toString treewidth CG);',
        "T = elimTree CG;",
        'print("elimTree treewidth " | toString treewidth T);',
        "N = chordalNet I;",
        'print("chordalNet " | toString N);',
        "elapsed chordalElim N;",
        'print("after chordalElim " | toString N);',
        "elapsed chordalTria N;",
        'print("after chordalTria " | toString N);',
        'print("isTriangular " | toString isTriangular N);',
        'print("codimCount " | toString codimCount N);',
        'print("dim(I) via net (codimCount above); ideal dim = " | toString dim I);',
    ]
    open(path, "w", encoding="utf-8").write("\n".join(body) + "\n")


def emit_m2_system(sysR, model, path, tag):
    verts, sups, spares, params = system_supports(sysR, model)
    if model == "spares":
        ringvars = spares
        coeffring = "frac(QQ[%s])" % ", ".join(params)
    else:
        ringvars = spares + params
        coeffring = "QQ"
    body = [
        'needsPackage "Chordal";',
        "kk = %s;" % coeffring,
        "R = kk[%s, MonomialOrder=>Lex];" % ", ".join(ringvars),
        "I = ideal(",
    ]
    polys = [e.expr for e in sysR.eqs]
    body.append(",\n".join("  " + m2_poly(p, ringvars) for p in polys))
    body += [
        ");",
        'print("== %s");' % tag,
        'print("nvars " | toString numgens R | "  ngens " | toString numgens I);',
        "G = constraintGraph I;",
        'print("constraintGraph " | toString G);',
        "CG = chordalGraph G;",
        'print("treewidth " | toString treewidth CG);',
        "T = elimTree CG;",
        'print("elimTree treewidth " | toString treewidth T);',
        'print("cliques " | toString cliques CG);',
    ]
    open(path, "w", encoding="utf-8").write("\n".join(body) + "\n")
    return len(polys), len(ringvars)


# =====================================================================
#  4.  Driver
# =====================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit-m2", action="store_true")
    ap.add_argument("--only", default=None, help="J6|SUB2|R9|GENERIC")
    ap.add_argument("--out", default=OUT_JSON)
    args = ap.parse_args()

    out = {"schema": 1,
           "item": "chordal / treewidth structural probe",
           "sympy": sp.__version__,
           "results": []}
    if args.emit_m2:
        os.makedirs(M2_DIR, exist_ok=True)

    sel = None if args.only is None else \
        {s.strip().upper() for s in args.only.split(",")}

    def want(t):
        return sel is None or t in sel

    # ---- J6 ---------------------------------------------------------
    if want("J6"):
        for rec in target_J6():
            r = analyse("J6 " + rec["key"], rec["supports"], rec["ring_vars"],
                        note="depth-3 saturated system, verbatim from "
                             "j6_msolve_results.json (msolve: EMPTY)")
            r["target"] = "J6"
            r["model"] = "ring vars (all)"
            print(fmt(r))
            out["results"].append(r)
            if args.emit_m2:
                fn = os.path.join(M2_DIR, "j6_%s.m2" %
                                  rec["key"].replace(":", "_").replace("#", "_"))
                emit_m2_j6(rec, fn)
                print("   emitted", fn)

    # ---- SUB2 -------------------------------------------------------
    if want("SUB2"):
        t0 = time.monotonic()
        sysR = target_SUB2()
        print("built SUB2 in %.1fs: %d eqs, %d spares, %d params"
              % (time.monotonic() - t0, len(sysR.eqs), len(sysR.vars),
                 len(sysR.params)))
        for model in ("spares", "full"):
            verts, sups, spares, params = system_supports(sysR, model)
            r = analyse("SUB2 home G-system [%s]" % model, sups, verts,
                        note="face_kill_sweep.build_state_system, sub2 caps; "
                             "d2 deg4 d1 deg6 sigma deg8 e=gamma*(y+1)^10")
            r["target"] = "SUB2"
            r["model"] = model
            r["n_spares"] = len(spares)
            r["n_params"] = len(params)
            print(fmt(r))
            out["results"].append(r)
            if args.emit_m2:
                fn = os.path.join(M2_DIR, "sub2_%s.m2" % model)
                np_, nv_ = emit_m2_system(sysR, model, fn, "SUB2 " + model)
                print("   emitted", fn, "(%d polys, %d vars)" % (np_, nv_))

    # ---- R9 ---------------------------------------------------------
    if want("R9"):
        t0 = time.monotonic()
        sysR = target_R9()
        print("built R9 in %.1fs: %d eqs, %d spares, %d params"
              % (time.monotonic() - t0, len(sysR.eqs), len(sysR.vars),
                 len(sysR.params)))
        for model in ("spares", "full"):
            verts, sups, spares, params = system_supports(sysR, model)
            r = analyse("R9 z=1 H-system [%s]" % model, sups, verts,
                        note="bigrade_annotator.build_R3 / r9_eliminated_system.json")
            r["target"] = "R9"
            r["model"] = model
            r["n_spares"] = len(spares)
            r["n_params"] = len(params)
            print(fmt(r))
            out["results"].append(r)
            if args.emit_m2:
                fn = os.path.join(M2_DIR, "r9_%s.m2" % model)
                np_, nv_ = emit_m2_system(sysR, model, fn, "R9 " + model)
                print("   emitted", fn, "(%d polys, %d vars)" % (np_, nv_))

    # ---- SYMBOL level (pre-expansion) --------------------------------
    if want("SYMBOL"):
        for tag, sups, names, gens in target_SYMBOL():
            r = analyse("SYMBOL " + tag, sups, names,
                        note="pre-expansion: vertices are the 8 series symbols")
            r["target"] = "SYMBOL"
            r["model"] = "series symbols"
            r["generator_supports"] = {k: sorted(str(s) for s in v.free_symbols)
                                       for k, v in gens.items()}
            print(fmt(r))
            for k, v in r["generator_supports"].items():
                print("     %-4s support %2d: %s" % (k, len(v), ", ".join(v)))
            out["results"].append(r)

    # ---- GENERIC (randomised-parameter fast path) --------------------
    if want("GENERIC"):
        for tag, builder in (("R9", target_R9_generic),
                             ("SUB2", target_SUB2_generic)):
            t0 = time.monotonic()
            sysR = builder()
            print("built %s-generic in %.1fs: %d eqs, %d spares"
                  % (tag, time.monotonic() - t0, len(sysR.eqs), len(sysR.vars)))
            verts, sups, spares, params = system_supports(sysR, "spares")
            r = analyse("%s [spares, params specialised]" % tag, sups, verts,
                        note="randomised-parameter fast path; support is a "
                             "SUBSET of the symbolic system's support")
            r["target"] = tag + "-generic"
            r["model"] = "spares (params specialised)"
            print(fmt(r))
            out["results"].append(r)

    for r in out["results"]:
        r.pop("_adj", None)
    json.dump(out, open(args.out, "w", encoding="utf-8"), indent=1)
    print("\nwrote", args.out)
    return 0


def fmt(r):
    g = r["graph"]
    return ("\n%-46s  V=%-4d E=%-6d density=%-8s maxdeg=%-4d\n"
            "   polys=%-5d maxsupport=%-4d  tw_lb(clique)=%-4d degeneracy=%-4d\n"
            "   mindeg : tw<=%-4d maxclique=%-4d fill=%d\n"
            "   minfill: tw<=%-4d maxclique=%-4d fill=%d"
            % (r["name"], g["vertices"], g["edges"], g["density"],
               g["max_degree"], r["n_polynomials"], r["max_support"],
               r["treewidth_lower_bound_from_clique"], r["degeneracy"],
               r["mindeg"]["treewidth_upper_bound"],
               r["mindeg"]["max_clique_in_completion"],
               r["mindeg"]["fill_edges_added"],
               r["minfill"]["treewidth_upper_bound"],
               r["minfill"]["max_clique_in_completion"],
               r["minfill"]["fill_edges_added"]))


if __name__ == "__main__":
    sys.exit(main())
