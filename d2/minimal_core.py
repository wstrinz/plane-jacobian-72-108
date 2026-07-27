#!/usr/bin/env python3
"""minimal_core.py -- the compression lane: what is load-bearing, and what theorem
did we prove?

This file adds NOTHING to the frontier and changes NO verdict.  It answers two
questions the (72,108) closure never answered:

  PART 1  (compression)  Which inputs are load-bearing?  Answered by ABLATION:
          drop a fact, recompute, see whether the conclusion survives.
          * A1  `deg e = 10` (sub2) rests on FOUR caps.  Exactly ONE of them,
                `deg R <= 12`, has zero slack; `deg d2 <= 4` has slack 2 and
                `deg S <= 14` has slack 1.  So three of the four caps could be
                loosened and the conclusion would stand.
          * A2  `a_t <= 9` is, downstream of the audited cascade base, PURE
                INTEGER ARITHMETIC: four term-valuations against one exact
                number.  Level 12 is load-bearing for exactly ONE of the four
                terms (`3*h_6^2`); levels 14/16 and `r_13 = 0` are not needed
                because `v_t(h_7) >= 11` already follows from the (P<)
                convolution floor with no cascade level at all.

  PART 2  (theorem)  The Catalan law of BELYI_PASSPORT.md is PROVED here, by the
          bijection its own N1 said was missing ("A proof would presumably come
          from the dessin being a plane trivalent tree-like map with one big
          face; I did not construct that bijection.").

          THEOREM.  For odd k, put m = (k-1)/2.  The genus-0 dessins with
          passport ( 2^((3k-1)/2), 1 | 3^k | (5k-1)/2, 1^((k+1)/2) ) are in
          canonical bijection with ROOTED PLANE BINARY TREES on m internal
          nodes.  Hence their number is the Catalan number C_m, every such
          dessin is connected, and every one has trivial automorphism group.

          Proof (three lines).  sigma_3 has (k+1)/2 fixed points, i.e. the map
          has (k+1)/2 - ... exactly m+1 MONOGON faces (self-loops at trivalent
          vertices) and one big face.  Deleting the m+1 loop edges leaves a
          connected planar map with V = 2m+2, E = 2m+1, hence (Euler) F = 1:
          a PLANE TREE, with m vertices of degree 3 and m+2 leaves.  The unique
          fixed point of sigma_1 (the leg) roots it at a leaf.  Rooted plane
          trees with m degree-3 vertices and m+2 leaves number C_m.  Reattaching
          a loop at each non-root leaf inverts the construction, and the two
          cyclic orders at such a leaf are conjugate by the transposition of the
          two loop darts, so the inverse is well defined.  []

          Two corollaries REPAIR real gaps in BELYI_PASSPORT.md, which proves
          connectedness and rigidity only at k = 7:
            * transitivity is automatic for every k (two lines, no enumeration);
            * |Aut| = 1 for every k -- whereas the repo's argument bounds |Aut|
              only by gcd(3k, (5k-1)/2), which is 3 (not 1) at k = 5, 11, 17,
              23, ...  Since the repo computes the Hurwitz number as N/n!, that
              quotient is the dessin count ONLY under |Aut| = 1, so the gap was
              load-bearing for the statement, not cosmetic.

  PART 2b (candidate b, sharpening)  WEIGHT_LEMMA_75_125's explanation of why
          the Phi-divisor syzygy exists at (72,108) and not at (75,125) is
          `q_window = 1`.  That is SUFFICIENT BUT NOT NECESSARY.  The exact
          criterion is proved here:

              the y-order carry obstruction vanishes  <=>  q_window | w(e).

          `q_window = 1` is the special case that divides every weight.  So the
          mechanism is not confined to corners with q_window = 1, and "(72,108)
          is sporadic" is NOT the explanation -- sporadicity and q_window are
          logically unrelated.

Read-only.  Writes nothing.  Usage:
    python minimal_core.py            # full report
    python minimal_core.py --quiet    # exit 0 iff every check passes
    python minimal_core.py --deep     # + Frobenius character sums to k = 19
"""
from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from functools import lru_cache
from math import factorial, gcd

sys.setrecursionlimit(100000)

# ==========================================================================
# check harness
# ==========================================================================
_OUT: list[str] = []
_NP = 0
_NT = 0
VERBOSE = True


def head(s: str) -> None:
    if VERBOSE:
        print("\n" + s)


def ck(cid: str, ok: bool, msg: str) -> bool:
    global _NP, _NT
    _NT += 1
    _NP += bool(ok)
    line = "  [%s] %-6s %s" % ("PASS" if ok else "FAIL", cid, msg)
    _OUT.append(line)
    if VERBOSE:
        print(line)
    return bool(ok)


def catalan(m: int) -> int:
    return factorial(2 * m) // (factorial(m) * factorial(m + 1))


# ==========================================================================
# PART 1, A1.  `deg e = 10` -- which of the four sub2 caps is load-bearing?
# ==========================================================================
DEG_PHI = 34
SUB2_CAPS = {"d2": 4, "R": 12, "S": 14, "e": 10}


def forced_deg_e(caps: dict) -> int | None:
    """Least E <= caps[e] with  E + max(d2 + 2E, E + S, 2R) >= deg Phi.

    This is `divisor_syzygy.forced_deg_e` re-implemented from the degree
    bookkeeping of  2*Phi = e*(d2*e^2 + 3*e*S + 3*R^2)  so that the ablation is
    not measuring the target's own code.
    """
    feas = [E for E in range(0, caps["e"] + 1)
            if E + max(caps["d2"] + 2 * E, E + caps["S"], 2 * caps["R"]) >= DEG_PHI]
    return min(feas) if feas else None


def section_A1() -> None:
    head("A1.  ablation: the four sub2 caps behind  deg e = 10")

    ok = forced_deg_e(SUB2_CAPS) == 10
    ck("A1.0", ok, "baseline caps %s force deg e = %s" % (SUB2_CAPS, forced_deg_e(SUB2_CAPS)))

    # cross-check against the target's own implementation (agreement, not import-of-truth)
    try:
        import divisor_syzygy as DS
        agree = DS.forced_deg_e(SUB2_CAPS)[0] == 10 and DS.DEG_PHI == DEG_PHI
        ck("A1.1", agree, "independent re-implementation agrees with divisor_syzygy "
                          "(E=%s, DEG_PHI=%s)" % (DS.forced_deg_e(SUB2_CAPS)[0], DS.DEG_PHI))
    except Exception as exc:                                          # noqa: BLE001
        ck("A1.1", False, "could not cross-check divisor_syzygy: %r" % exc)

    # per-cap slack
    slack = {}
    for cap in SUB2_CAPS:
        good = [v for v in range(0, 60)
                if forced_deg_e({**SUB2_CAPS, cap: v}) == 10]
        lo, hi = (min(good), max(good)) if good else (None, None)
        slack[cap] = (SUB2_CAPS[cap] - lo if lo is not None else None,
                      hi - SUB2_CAPS[cap] if hi is not None else None)

    ck("A1.2", slack["R"] == (12, 0),
       "deg R <= 12 has ZERO upward slack: R = 13 already re-admits E = 9.  "
       "This is the ONE binding cap.  (down,up) = %s" % (slack["R"],))
    ck("A1.3", slack["d2"][1] == 2,
       "deg d2 <= 4 has upward slack 2 (holds to d2 <= 6) -- NOT at the margin")
    ck("A1.4", slack["S"][1] == 1,
       "deg S <= 14 has upward slack 1 (holds to S <= 15) -- NOT at the margin")
    ck("A1.5", slack["e"] == (0, 49),
       "the cap deg e <= 10 is used only as the ceiling E <= 10; it cannot be "
       "lowered at all and may be raised freely.  (down,up) = %s" % (slack["e"],))

    # the triple coincidence at E = 10 -- and it is a coincidence, not a theorem
    at10 = (SUB2_CAPS["d2"] + 20, 10 + SUB2_CAPS["S"], 2 * SUB2_CAPS["R"])
    ck("A1.6", at10 == (24, 24, 24),
       "at E = 10 all three RHS branches equal 24 simultaneously: %s.  The "
       "forcing is a three-way tie, which is why only the max matters and only "
       "one cap binds" % (at10,))

    # MUTATION CONTROL: the check must be able to FAIL
    mut = forced_deg_e({**SUB2_CAPS, "R": 13})
    ck("A1.7", mut is not None and mut != 10,
       "mutation control: R -> 13 gives E = %s != 10, so A1.2 is a real test "
       "and not a tautology" % mut)
    mut2 = forced_deg_e({**SUB2_CAPS, "R": 11, "S": 13, "d2": 3})
    ck("A1.8", mut2 is None,
       "mutation control: shrinking all three caps makes the system INFEASIBLE "
       "(E = %s), i.e. the caps cannot be freely tightened either" % mut2)


# ==========================================================================
# PART 1, A2.  `a_t <= 9` is integer arithmetic
# ==========================================================================
V_PHI = 30                    # v_t(Phi), EXACT (divisor_consequences / at_le9_audit B7)
BASE = {1: 1, 2: 3, 3: 5, 4: 7}   # the AUDITED cascade rows v_t(h_k) >= 2k-1, k<=4


def convolution_floor(v: dict, n: int) -> int:
    """(P<) absorption  h_n = -q_n/2 + t^(2n-2) g_n  =>
       v(h_n) >= min( min_{i+j=n} (v(h_i)+v(h_j)), 2n-2 ).
    This is spine9_audit's G10/G11 floor, re-implemented."""
    conv = min((v[i] + v[n - i] for i in range(1, n) if i in v and (n - i) in v),
               default=10 ** 9)
    return min(conv, 2 * n - 2)


def profile(a: int, level12: bool, base: dict | None = None,
            L6: int = 11) -> dict:
    """v_t(h_1..h_7) as a function of a_t = a, with or without cascade level 12.
    `base` and `L6` are ABLATION knobs."""
    v = dict(BASE if base is None else base)
    v[5] = a
    v[6] = convolution_floor(v, 6)
    if level12:
        v[6] = max(v[6], L6)          # the single new input from cascade level 12
    v[7] = convolution_floor(v, 7)
    return v


def bracket_terms(v: dict, a: int) -> dict:
    """B = h_2*h_5^2 + 3*h_5*h_7 + 3*h_1*h_5*h_6 + 3*h_6^2   (syzygy_collision X6.1).
    Returns the t-valuation floor of each of the four terms."""
    return {"h2*h5^2": v[2] + 2 * a,
            "3*h5*h7": a + v[7],
            "3*h1*h5*h6": v[1] + a + v[6],
            "3*h6^2": 2 * v[6]}


def collision_kills(a: int, level12: bool, base: dict | None = None,
                    L6: int = 11, vphi: int | None = None) -> tuple[bool, dict, int]:
    """The K-syzygy gives 2*Phi = e*B with v_t(e) = a, so a + v_t(B) = v_t(Phi)
    EXACTLY.  If every term of B has valuation > v_t(Phi) - a then
    v_t(B) > v_t(Phi) - a: contradiction, so that a is impossible."""
    v = profile(a, level12, base, L6)
    T = bracket_terms(v, a)
    need = (V_PHI if vphi is None else vphi) - a
    return (min(T.values()) > need), T, need


def conclusion(base: dict | None = None, L6: int = 11,
               vphi: int | None = None) -> bool:
    """The FULL conclusion the collision must deliver: a_t = 9 survives (it must,
    the frontier's five cells sit there) AND every a_t >= 10 is refuted."""
    if not collision_kills(9, True, base, L6, vphi)[0]:
        return all(collision_kills(a, True, base, L6, vphi)[0]
                   for a in range(10, 31))
    return False


def section_A2() -> None:
    head("A2.  ablation: `a_t <= 9` downstream of the cascade base is 4 integers")

    # the kill, with level 12
    surv = [a for a in range(1, 31) if not collision_kills(a, True)[0]]
    ck("A2.0", surv == [1, 2, 3, 4, 5, 6, 7, 8, 9],
       "with level 12, every a_t >= 10 is refuted and a_t <= 9 survives.  "
       "surviving a: %s" % surv)

    kills9 = collision_kills(9, True)
    ck("A2.1", not kills9[0] and kills9[1]["h2*h5^2"] == 21 and kills9[2] == 21,
       "NON-VACUITY: at a = 9 the term h_2*h_5^2 lands on 21 = 30 - 9 exactly, "
       "so a = 9 is NOT refuted.  The binding term at the boundary is "
       "h_2*h_5^2, not h_7 (terms: %s, need > %d)" % (kills9[1], kills9[2]))

    k10 = collision_kills(10, True)
    ck("A2.2", k10[0] and min(k10[1].values()) == 22 and k10[2] == 20,
       "a = 10 dies with margin 2: terms %s all exceed %d" % (k10[1], k10[2]))

    # ABLATION: drop cascade level 12
    n10 = collision_kills(10, False)
    ck("A2.3", (not n10[0]) and n10[1]["3*h6^2"] == 20 == n10[2],
       "ABLATION -- drop cascade level 12 (v_t(h_6) >= 11 -> 10): a = 10 SURVIVES, "
       "and it survives on exactly ONE term, 3*h_6^2 at %d = %d.  Level 12 is "
       "load-bearing, and only for that term" % (n10[1]["3*h6^2"], n10[2]))
    others = {k: x for k, x in n10[1].items() if k != "3*h6^2"}
    ck("A2.4", all(x > n10[2] for x in others.values()),
       "...and the other three terms clear the bar without level 12: %s > %d"
       % (others, n10[2]))

    # v_t(h_7) >= 11 needs NO cascade level
    v_nolvl = profile(10, False)
    ck("A2.5", v_nolvl[7] == 11,
       "v_t(h_7) >= 11 follows from the (P<) convolution floor on the AUDITED "
       "rows (1,3,5,7) alone -- no cascade level, hence levels 14 and 16 are "
       "not needed and `r_13 = 0` is not needed.  v(h_7) = %d" % v_nolvl[7])
    ck("A2.6", profile(10, False)[6] == 10 and profile(10, True)[6] == 11,
       "the ONLY input level 12 supplies is v_t(h_6): 10 -> 11")

    # the alternate regime values
    ck("A2.7", all(collision_kills(a, True)[0] for a in (11, 12, 14)),
       "the alternate-regime values a_t = 12, 14 (and 11) are re-killed by the "
       "same four integers")

    # ---- the SLACK TABLE: per-input admissible range for the full conclusion
    ck("A2.8", conclusion(),
       "the full conclusion (a = 9 survives AND every a >= 10 dies) holds at the "
       "audited inputs -- this is the predicate every ablation below perturbs")

    vphi_ok = [v for v in range(20, 46) if conclusion(vphi=v)]
    ck("A2.9", vphi_ok == [30, 31],
       "ABLATION v_t(Phi): admissible range is %s.  ZERO downward slack (29 "
       "kills a = 9 too, destroying the frontier) and slack 1 upward.  So "
       "`v_t(Phi) = 30' must be known exactly on the low side; the earlier draft "
       "of this check guessed 31 would break it and was WRONG -- mutation "
       "testing caught that" % vphi_ok)

    l6_ok = [L for L in range(0, 20) if conclusion(L6=L)]
    ck("A2.10", min(l6_ok) == 11,
       "ABLATION cascade level 12: admissible v_t(h_6) floor is >= %d.  ZERO "
       "downward slack.  This sharpens SYZYGY_COLLISION X11 from 'load-bearing' "
       "to 'load-bearing with zero margin'" % min(l6_ok))

    rng = {}
    for j in (1, 2, 3, 4):
        rng[j] = [x for x in range(0, 16)
                  if conclusion(base={**BASE, j: x})]
    ck("A2.11", rng[1] == list(range(0, 16)),
       "ABLATION v_t(h_1) >= 1: EVERY value 0..15 still gives the conclusion.  "
       "The audited row v_t(h_1) >= 1 is NOT NEEDED by the collision at all -- "
       "it is scaffolding here (it is still needed to start the cascade)")
    ck("A2.12", min(rng[2]) == 1 and min(rng[3]) == 4 and min(rng[4]) == 6,
       "ABLATION the other audited rows: v_t(h_2) needs only >= %d (audited 3, "
       "slack 2); v_t(h_3) only >= %d (audited 5, slack 1); v_t(h_4) only >= %d "
       "(audited 7, slack 1)" % (min(rng[2]), min(rng[3]), min(rng[4])))

    # h_3 and h_4 enter only jointly, through v(h_7)
    joint = [(x, y) for x in range(0, 14) for y in range(0, 14)
             if conclusion(base={**BASE, 3: x, 4: y})]
    ck("A2.13", joint and min(x + y for x, y in joint) == 11,
       "...and they enter ONLY through the sum: the real requirement is "
       "v_t(h_3) + v_t(h_4) >= %d, against an audited 5 + 7 = 12.  Joint slack "
       "is 1, not 2" % min(x + y for x, y in joint))

    ck("A2.14", collision_kills(10, True)[1] and
       min(collision_kills(10, True)[1].values()) - collision_kills(10, True)[2] == 2,
       "MARGIN of the a = 10 kill is 2 in t-valuation.  NET RESULT: `a_t <= 9' "
       "consumes exactly FIVE numbers -- v_t(h_2) >= 1, v_t(h_3)+v_t(h_4) >= 11, "
       "v_t(h_6) >= 11, v_t(Phi) = 30 -- plus the bracket B and the syzygy")


# ==========================================================================
# PART 2.  the Catalan law
# ==========================================================================
def passport(k: int) -> tuple[tuple, tuple, tuple, int]:
    n = 3 * k
    P1 = tuple(sorted([2] * ((3 * k - 1) // 2) + [1], reverse=True))
    P2 = tuple([3] * k)
    P3 = tuple(sorted([(5 * k - 1) // 2] + [1] * ((k + 1) // 2), reverse=True))
    return P1, P2, P3, n


def binary_trees(m: int):
    """rooted plane binary trees with m internal nodes; leaf = None"""
    if m == 0:
        yield None
        return
    for i in range(m):
        for L in binary_trees(i):
            for Rt in binary_trees(m - 1 - i):
                yield (L, Rt)


def build_dessin(T, loop_mode: str = "ok"):
    """The bijection's inverse direction: tree -> dessin.
    loop_mode is a MUTATION knob; 'ok' is the real construction."""
    s1: dict[int, int] = {}
    s2: dict[int, int] = {}
    c = [0]

    def nd() -> int:
        c[0] += 1
        return c[0] - 1

    leg = nd()

    def visit(node, pd):
        if node is None:                        # leaf -> trivalent via a self loop
            a, b = nd(), nd()
            if loop_mode == "pair_parent":       # deliberately wrong
                s1[a] = pd
                s1[pd] = a
                s1[b] = b
            else:
                s1[a] = b
                s1[b] = a
            cyc = (pd, a, b)
        else:
            L, Rt = node
            du, dd_ = nd(), nd()
            s1[du], s1[dd_] = dd_, du
            eu, ed = nd(), nd()
            s1[eu], s1[ed] = ed, eu
            cyc = (pd, du, eu)
            visit(L, dd_)
            visit(Rt, ed)
        for i in range(3):
            s2[cyc[i]] = cyc[(i + 1) % 3]

    visit(T, leg)
    if loop_mode != "pair_parent":
        s1[leg] = leg                            # the unique fixed point
    return c[0], s1, s2, leg


def cycle_type(p, n: int) -> tuple:
    seen = [False] * n
    out = []
    for i in range(n):
        if seen[i]:
            continue
        L, j = 0, i
        while not seen[j]:
            seen[j] = True
            j = p[j]
            L += 1
        out.append(L)
    return tuple(sorted(out, reverse=True))


def is_transitive(s1, s2, n: int) -> bool:
    seen, st = {0}, [0]
    while st:
        x = st.pop()
        for p in (s1, s2):
            if p[x] not in seen:
                seen.add(p[x])
                st.append(p[x])
    return len(seen) == n


def canonical(n: int, s1, s2, seed: int):
    """BFS relabel from a canonically distinguished dart.  Deterministic, so
    isomorphic dessins have equal canonical forms (verified both ways below)."""
    lab = {seed: 0}
    order = [seed]
    i = 0
    while i < len(order):
        x = order[i]
        i += 1
        for p in (s2, s1):
            if p[x] not in lab:
                lab[p[x]] = len(order)
                order.append(p[x])
    if len(lab) != n:
        return None
    return (tuple(lab[s1[order[j]]] for j in range(n)),
            tuple(lab[s2[order[j]]] for j in range(n)))


def aut_order(n: int, s1, s2) -> int:
    """|Aut| by trying every image of dart 0 and propagating."""
    cnt = 0
    for tgt in range(n):
        pi = {0: tgt}
        order, i, ok = [0], 0, True
        while i < len(order) and ok:
            x = order[i]
            i += 1
            for p in (s1, s2):
                y, want = p[x], p[pi[x]]
                if y in pi:
                    if pi[y] != want:
                        ok = False
                        break
                else:
                    pi[y] = want
                    order.append(y)
        if ok and len(pi) == n and len(set(pi.values())) == n:
            cnt += 1
    return cnt


def section_B_family() -> None:
    head("B.  the passport family is a genus-0 family of degree 3k")
    bad_sum, bad_g = [], []
    for k in range(1, 40, 2):
        P1, P2, P3, n = passport(k)
        if not all(sum(P) == n for P in (P1, P2, P3)):
            bad_sum.append(k)
        if sum(n - len(P) for P in (P1, P2, P3)) != 2 * n - 2:
            bad_g.append(k)
    ck("B1", not bad_sum, "every partition sums to the degree n = 3k, for all odd "
                          "k <= 39 (deg = 3k, so k = 7 -> 21)")

    # B2 is NOT a test.  BELYI_PASSPORT.md K1 makes exactly this point about its
    # own RH check, and the brief for this lane names it as a known trap.  So
    # prove it is an identity, and supply a control showing the predicate CAN fail.
    import sympy as sp
    kk = sp.Symbol("kk")
    rh = (sp.Rational(3, 1) * kk - ((3 * kk + 1) / 2)      # n - #parts(P1)
          + 3 * kk - kk                                     # n - #parts(P2)
          + 3 * kk - ((kk + 3) / 2))                        # n - #parts(P3)
    ck("B2", sp.simplify(rh - (2 * (3 * kk) - 2)) == 0 and not bad_g,
       "[IDENTITY, NOT A TEST] Riemann-Hurwitz balances SYMBOLICALLY in k: "
       "the deficiency sum minus (2n-2) simplifies to %s for all k.  So 'RH "
       "balances' is vacuous as evidence -- it is recorded only to state that "
       "the family is genus 0 by construction" % sp.simplify(rh - (2 * (3 * kk) - 2)))

    perturbed_ok = sum(n - len(P) for P in
                       (passport(7)[0], passport(7)[1],
                        tuple(sorted([16, 2, 1, 1, 1], reverse=True)))) == 2 * 21 - 2
    ck("B3", not perturbed_ok,
       "CONTROL for B2: the perturbed third partition (16,2,1,1,1) does NOT "
       "balance RH, so the predicate in B2 is at least capable of failing -- "
       "it is the FAMILY that makes it automatic, not the code")

    P1, P2, P3, n = passport(7)
    ck("B4", (P1, P2, P3, n) == ((2,) * 10 + (1,), (3,) * 7, (17, 1, 1, 1, 1), 21),
       "k = 7 reproduces (72,108)'s passport (2^10,1 | 3^7 | 17,1^4), degree 21")
    ck("B5", all((3 * k - 1) % 2 == 1 and (5 * k - 1) % 2 == 1
                 for k in range(2, 40, 2)),
       "for EVERY even k both (3k-1)/2 and (5k-1)/2 are non-integral, so the "
       "family genuinely has no even members and the reindexing m = (k-1)/2 is "
       "onto the non-negative integers")


def section_C_lemmas() -> None:
    head("C.  the two structural lemmas (these REPAIR BELYI_PASSPORT.md, which "
         "proves both only at k = 7)")

    # LEMMA 1: transitivity is automatic.  An orbit missing the long cycle
    # consists of sigma_3-fixed points, so sigma_2 = sigma_1^{-1} there; an
    # element of order dividing both 2 and 3 is trivial; but sigma_2 has no
    # fixed points.  The only machine-checkable atom is:
    import itertools
    viol = []
    for nn in (3, 6):
        for p in itertools.permutations(range(nn)):
            t = cycle_type(list(p), nn)
            if t == tuple([3] * (nn // 3)):
                # is p an involution?  (that is what sigma_2 = sigma_1 forces)
                if all(p[p[i]] == i for i in range(nn)):
                    viol.append((nn, p))
    ck("C1", not viol,
       "LEMMA 1 atom: no permutation of cycle type 3^j is an involution "
       "(checked exhaustively on 3 and 6 points, %d violations).  Hence an "
       "orbit avoiding the long cycle is empty => EVERY triple with this "
       "passport is transitive, for every k" % len(viol))
    # mutation control for C1: type 2^j IS an involution, so the test discriminates
    inv2 = [p for p in itertools.permutations(range(4))
            if cycle_type(list(p), 4) == (2, 2) and all(p[p[i]] == i for i in range(4))]
    ck("C2", len(inv2) == 3,
       "mutation control: the same predicate finds %d involutions of type 2^2 on "
       "4 points, so C1 is discriminating and not vacuously true" % len(inv2))

    # LEMMA 2: |Aut| = 1.  |Fix(sigma_1)| = 1 in every member of the family.
    fixcounts = {k: passport(k)[0].count(1) for k in range(1, 40, 2)}
    ck("C3", set(fixcounts.values()) == {1},
       "LEMMA 2 hypothesis: sigma_1 has EXACTLY ONE fixed point for every odd k "
       "(the passport's trailing 1).  The centralizer of a transitive group is "
       "semiregular, so an automorphism permuting a single fixed point must fix "
       "it, hence is trivial => |Aut| = 1 for every k")

    # the gap this closes
    gapk = [k for k in range(1, 40, 2) if gcd(3 * k, (5 * k - 1) // 2) != 1]
    ck("C4", gapk and gapk[:4] == [5, 11, 17, 23],
       "GAP CLOSED: BELYI_PASSPORT.md's rigidity argument bounds |Aut| by "
       "gcd(3k, (5k-1)/2), which is 3 (not 1) at k = %s...  Since the repo "
       "computes the Hurwitz number as N/n!, |Aut| = 1 is required for that "
       "quotient to BE the dessin count, so the gap was load-bearing"
       % gapk[:4])


def section_D_bijection(mmax: int = 7) -> None:
    head("D.  THE BIJECTION  rooted plane binary trees (m nodes) -> dessins")
    allgood = True
    aut_seen: set[int] = set()
    for m in range(0, mmax + 1):
        k = 2 * m + 1
        P1, P2, P3, n = passport(k)
        seen, ntree, bad = set(), 0, 0
        for T in binary_trees(m):
            ntree += 1
            nn, s1, s2, leg = build_dessin(T)
            s1l = [s1[i] for i in range(nn)]
            s2l = [s2[i] for i in range(nn)]
            prod = [s1l[s2l[i]] for i in range(nn)]
            if nn != n:
                bad += 1
                continue
            if (cycle_type(s1l, nn), cycle_type(s2l, nn),
                    cycle_type(prod, nn)) != (P1, P2, P3):
                bad += 1
            if not is_transitive(s1, s2, nn):
                bad += 1
            cf = canonical(nn, s1, s2, leg)
            if cf is None:
                bad += 1
            else:
                seen.add(cf)
            if m <= 4:
                aut_seen.add(aut_order(nn, s1l, s2l))
        C = catalan(m)
        good = (ntree == C and len(seen) == C and bad == 0)
        allgood &= good
        ck("D1(m=%d)" % m, good,
           "k=%2d n=%2d : %3d trees -> %3d PAIRWISE NON-ISOMORPHIC dessins, all "
           "with the exact passport and transitive; Catalan_%d = %d"
           % (k, n, ntree, len(seen), m, C))
    ck("D2", aut_seen == {1},
       "|Aut| = 1 on every constructed dessin up to m = 4, by brute force over "
       "all n images of one dart (values seen: %s)" % sorted(aut_seen))

    # MUTATION CONTROL: a wrong loop attachment must break the passport
    m = 3
    P1, P2, P3, n = passport(2 * m + 1)
    survived = 0
    for T in binary_trees(m):
        try:
            nn, s1, s2, leg = build_dessin(T, loop_mode="pair_parent")
            s1l = [s1[i] for i in range(nn)]
            s2l = [s2[i] for i in range(nn)]
            prod = [s1l[s2l[i]] for i in range(nn)]
            if (cycle_type(s1l, nn), cycle_type(s2l, nn),
                    cycle_type(prod, nn)) == (P1, P2, P3):
                survived += 1
        except Exception:                                            # noqa: BLE001
            pass
    ck("D3", survived == 0,
       "mutation control: pairing the loop dart to the PARENT dart instead of "
       "its twin breaks the passport on all %d trees at m = 3 (%d survived) -- "
       "so D1 is a real test of the construction" % (catalan(m), survived))

    # SOUNDNESS of the canonical form: conjugation-invariant AND seed-sensitive
    import random
    rng = random.Random(20260726)
    nfail = 0
    ntest = 0
    for m in range(0, 7):
        for T in binary_trees(m):
            nn, s1, s2, leg = build_dessin(T)
            base = canonical(nn, s1, s2, leg)
            for _ in range(4):
                p = list(range(nn))
                rng.shuffle(p)
                q = {p[i]: p[s1[i]] for i in range(nn)}
                r = {p[i]: p[s2[i]] for i in range(nn)}
                ntest += 1
                if canonical(nn, q, r, p[leg]) != base:
                    nfail += 1
    ck("D4", nfail == 0 and ntest > 300,
       "canonical form is SOUND: invariant under %d random relabellings "
       "(%d failures), so D1's 'pairwise non-isomorphic' is not an artefact of "
       "labelling" % (ntest, nfail))
    nn, s1, s2, leg = build_dessin((None, (None, None)))
    alt = next(d for d in range(nn) if d != leg)
    ck("D5", canonical(nn, s1, s2, leg) != canonical(nn, s1, s2, alt),
       "canonical form is DISCRIMINATING: reseeding at a non-leg dart changes "
       "it, so D4 is not the trivial 'always equal' failure mode")


# ==========================================================================
# E.  independent count: Frobenius character sums
# ==========================================================================
def _partitions(n: int, maxp: int | None = None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in _partitions(n - p, p):
            yield (p,) + rest


def _to_beta(lam):
    l = len(lam)
    return tuple(lam[i] + (l - 1 - i) for i in range(l))


def _from_beta(beta):
    b = sorted(beta, reverse=True)
    l = len(b)
    return tuple(x for x in (b[i] - (l - 1 - i) for i in range(l)) if x > 0)


def _rim_hooks(lam, r: int):
    beta = list(_to_beta(lam))
    S = set(beta)
    for i, b in enumerate(beta):
        t = b - r
        if t < 0 or t in S:
            continue
        ht = sum(1 for x in beta if t < x < b)
        yield (-1) ** ht, _from_beta(beta[:i] + [t] + beta[i + 1:])


def _mk_chi(rho):
    @lru_cache(maxsize=None)
    def chi(lam, j):
        if j == len(rho):
            return 1 if sum(lam) == 0 else 0
        return sum(s * chi(mu, j + 1) for s, mu in _rim_hooks(lam, rho[j]))
    return lambda lam: chi(tuple(lam), 0)


def _dim(lam, n: int) -> int:
    conj = [sum(1 for x in lam if x > j) for j in range(lam[0])] if lam else []
    prod = 1
    for i, li in enumerate(lam):
        for j in range(li):
            prod *= (li - j) + (conj[j] - i) - 1
    return factorial(n) // prod


def _class_size(n: int, rho) -> int:
    from collections import Counter
    d = 1
    for L, mult in Counter(rho).items():
        d *= (L ** mult) * factorial(mult)
    return factorial(n) // d


def hurwitz(C1, C2, C3) -> Fraction:
    """#{triples in the three classes with product 1} / n!.  Equals the number of
    isomorphism classes because Lemmas 1 and 2 give transitivity and |Aut| = 1."""
    n = sum(C1)
    chi1, chi2, chi3 = (_mk_chi(tuple(sorted(C, reverse=True))) for C in (C1, C2, C3))
    tot = Fraction(0)
    for lam in _partitions(n):
        v3 = chi3(lam)
        if v3 == 0:
            continue
        v1 = chi1(lam)
        if v1 == 0:
            continue
        v2 = chi2(lam)
        if v2 == 0:
            continue
        tot += Fraction(v1 * v2 * v3, _dim(lam, n))
    N = Fraction(_class_size(n, C1) * _class_size(n, C2) * _class_size(n, C3),
                 factorial(n)) * tot
    assert N.denominator == 1
    return Fraction(int(N), factorial(n))


def section_E_frobenius(kmax: int) -> None:
    head("E.  independent count by Frobenius character sums (kmax = %d)" % kmax)
    # calibration on a case where brute force is possible
    C1, C2, C3, n = passport(1)
    ck("E0", hurwitz(C1, C2, C3) == 1,
       "calibration: the character machinery reproduces k = 1 (n = 3) -> 1")
    bad = []
    for k in range(1, kmax + 1, 2):
        m = (k - 1) // 2
        P1, P2, P3, n = passport(k)
        H = hurwitz(P1, P2, P3)
        if H != catalan(m):
            bad.append((k, H, catalan(m)))
        ck("E1(k=%d)" % k, H == catalan(m),
           "n=%2d  Hurwitz = %s  =  Catalan_%d = %d" % (n, H, m, catalan(m)))
    ck("E2", not bad,
       "the Catalan law holds by an INDEPENDENT method (character sums, no "
       "trees) for every odd k <= %d -- i.e. through C_%d = %d, well past the "
       "four data points BELYI_PASSPORT.md had"
       % (kmax, (kmax - 1) // 2, catalan((kmax - 1) // 2)))
    # a negative control: a neighbouring passport must NOT be Catalan
    P1, P2, P3, n = passport(7)
    P3b = tuple(sorted([16, 2, 1, 1, 1], reverse=True))     # perturb the third partition
    Hb = hurwitz(P1, P2, P3b)
    ck("E3", Hb != 5,
       "negative control: perturbing the third partition to %s gives Hurwitz = "
       "%s != 5, so E1 is a property of THIS passport family and not of the "
       "method" % (P3b, Hb))


# ==========================================================================
# F.  candidate (b): the exact carry criterion
# ==========================================================================
def Lcap(w: int, alpha: int, q: int) -> int:
    return -((-alpha * w) // q)


def carry(w1: int, w2: int, alpha: int, q: int) -> int:
    return Lcap(w1, alpha, q) + Lcap(w2, alpha, q) - Lcap(w1 + w2, alpha, q)


def section_F_carry() -> None:
    head("F.  candidate (b): q_window = 1 is SUFFICIENT, NOT NECESSARY")
    try:
        import weight_lemma_75_125 as W
        # REPAIRED 2026-07-26.  These calls used the SUPERSEDED (75,125)
        # signature (ord,M,deg) = (201,36,504) -> q_window = 12, and w(e) = 6.
        # The repaired signature is (80,29,80) -> q_window = 29, w(e) = t+1 = 5.
        # The old calls PASSED, because window_law is general arithmetic and
        # (67,12) is internally consistent -- it simply stopped describing this
        # case.  A gated check whose LABEL says "(75,125)" while its numbers do
        # not is the failure mode this whole audit exists to catch.
        ck("F0", (W.window_law(204, 17, 238)[1] == 1
                  and W.window_law(80, 29, 80)[1] == 29
                  and W.carry(5, 12, 12, 1) == 0
                  and W.carry(5, 24, 80, 29) == 1),
           "reproduce the repo's two data points from its own primitives: "
           "(72,108) q_window = 1, carry = 0; (75,125) q_window = 29, carry = 1 "
           "at w(e) = 5")
        _carry = W.carry
    except Exception as exc:                                          # noqa: BLE001
        ck("F0", False, "could not load weight_lemma_75_125: %r" % exc)
        _carry = carry

    # THE CRITERION.  Since q | alpha*M (that is what q_window being the
    # denominator of ordPhi/M means) and gcd(alpha,q) = 1, the two residues
    # (-alpha*w_e mod q) and (alpha*w_e mod q) sum to q or to 0.  Hence:
    bad, tested, zeros = [], 0, 0
    for q in range(1, 31):
        for alpha in range(1, 45):
            if gcd(alpha, q) != 1:
                continue
            for M in range(q, 160, q):
                for we in range(1, M):
                    c = _carry(we, M - we, alpha, q)
                    tested += 1
                    zeros += (c == 0)
                    if c != (0 if we % q == 0 else 1):
                        bad.append((q, alpha, M, we, c))
    ck("F1", not bad and tested > 10 ** 6,
       "CRITERION, exhaustive over %d admissible splits (%d with carry 0, %d "
       "mismatches): carry(w_e, M - w_e) = 0  <=>  q_window | w_e.  The carry "
       "is exactly 0 or 1, never more" % (tested, zeros, len(bad)))

    witnesses = sum(1 for q in range(2, 20) for alpha in range(1, 20)
                    if gcd(alpha, q) == 1
                    for M in range(q, 100, q) for we in range(1, M)
                    if _carry(we, M - we, alpha, q) == 0)
    ck("F2", witnesses > 1000,
       "MUTATION CONTROL on the repo's explanation: %d admissible splits have "
       "carry = 0 with q_window != 1.  So `q_window = 1` is sufficient but NOT "
       "necessary, and CANNOT be the characterisation of when the Phi-divisor "
       "syzygy can exist" % witnesses)

    ck("F3", carry(5, 24, 80, 29) == 1
       and all(carry(w, 29 - w, 80, 29) == 1 for w in range(1, 29)),
       "at (75,125) the obstruction is w(e) = 5 against q_window = 29, and it is "
       "STRICTLY STRONGER than the pre-repair reading: because q_window = 29 = M "
       "exactly, NO split with 1 <= w < M is a multiple of q_window, so the carry "
       "is 1 on EVERY admissible split and there is no weight at which the "
       "obstruction could vanish.  (The superseded chart gave q_window = 12 "
       "against M = 36, where it would have vanished at w in {12,24}.)  The "
       "corner is still ruled out by 'q_window does not divide w(e)', not by "
       "'q_window != 1'")

    ck("F4", all(carry(w, 17 - w, 12, 1) == 0 for w in range(1, 17)),
       "at (72,108), q_window = 1 divides every weight, so the obstruction "
       "vanishes for ALL 16 splits of M = 17 -- the syzygy is unobstructed for "
       "structural reasons that have nothing to do with sporadicity")


# ==========================================================================
def main() -> int:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--deep", action="store_true",
                    help="Frobenius sums to k = 19 (slow) instead of k = 13")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    section_A1()
    section_A2()
    section_B_family()
    section_C_lemmas()
    section_D_bijection(7)
    section_E_frobenius(19 if a.deep else 13)
    section_F_carry()

    if a.quiet:
        if _NP != _NT:
            print("minimal_core: %d/%d checks FAILED" % (_NT - _NP, _NT))
            for line in _OUT:
                if "[FAIL]" in line:
                    print(line)
            return 1
        print("minimal_core: %d/%d checks pass" % (_NP, _NT))
        return 0
    print("\n%d/%d checks pass" % (_NP, _NT))
    if _NP == _NT:
        print("""
  PROVED   Hurwitz(passport_k) = Catalan_{(k-1)/2}, by a bijection with rooted
           plane binary trees; plus transitivity and |Aut| = 1 for every k
           (BELYI_PASSPORT.md had both only at k = 7, and its rigidity
           argument provably fails at k = 5, 11, 17, 23, ...).
  PROVED   the y-order carry obstruction vanishes iff q_window | w(e).
           `q_window = 1' is the sufficient special case, not the criterion.
  CHECKED  `deg e = 10' binds on ONE cap, deg R <= 12 (zero slack).
  CHECKED  `a_t <= 9' is 4 term-valuations against 1 exact number; cascade
           level 12 is load-bearing for exactly one term, 3*h_6^2.""")
    return 0 if _NP == _NT else 1


if __name__ == "__main__":
    sys.exit(main())
