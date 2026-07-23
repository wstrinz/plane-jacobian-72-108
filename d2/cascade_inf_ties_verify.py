#!/usr/bin/env python3
"""Stage 2 of the Phase D infinity layer: exact tie equations at infinity.

Verifies that the max-plus tie obligations emitted by cascade_engine.py
carry the correct EXPLICIT leading-coefficient equations:

  A. the leading constants from the source: lc(q) = 2048, q(-1) = 3315,
     lc(u) = lc(Phi) = c*lc(q) = -1024/3315 (t = y+1 monic), matching the
     engine's LC_U label constant;
  B. the depth-one initial form of an infinity degree tie is the SAME
     polynomial in leading coefficients as the residue-lemma initial form
     (RESIDUE_LEMMAS.md section 1) — verified for the three backbone
     hypersurfaces P6 (level 6), P10 (level 5), P11 (level 4), whose full
     ties are all achieved at the degree state (deg d2, deg d1, deg sigma,
     deg e) = (2,3,4,5); the engine-side equation is reconstructed by
     parsing the tie labels, the reference side directly from the source
     tables, and both are compared to the published relations;
  C. the two arithmetic KILL lemmas C08 (level 5) and C20 (level 4) act as
     forbidden DROPS at infinity under --residue-kills: their depth-one
     equations have no all-nonzero solution over Q or the q-splitting
     field (proven in residue_lemmas_verify.py), and at infinity the
     unknowns are the LEADING coefficients, which live in the base field —
     so a tie on exactly that support can neither drop nor vanish;
  D. the engine-emitted leading_cancellation labels for the a_t=9 T2
     linear-E witness reconstruct the exact equation
     lc(ehat)^3*lc(g5) + (-1024/3315)^5*2048*lc(e)^2 = 0.

Run: python cascade_inf_ties_verify.py
"""

from __future__ import annotations

import re

import sympy as sp

import cascade_engine as ce
import t5_90t1_verify as base

NEG_INF = ce.NEG_INF
D, X, S, E = sp.symbols("D X S E")  # leading coeffs of d2, d1, sigma, e
LABEL = re.compile(r"^(-?\d+)\*d2\^(\d+)\*d1\^(\d+)\*sigma\^(\d+)\*e\^(\d+)$")


def parse_labels(labels: tuple[str, ...]) -> sp.Expr:
    total = sp.Integer(0)
    for label in labels:
        match = LABEL.match(label)
        assert match, label
        coefficient, kd, xd, zd, ed = (int(g) for g in match.groups())
        total += coefficient * D**kd * X**xd * S**zd * E**ed
    return total


def check_a_constants() -> None:
    y = base.y
    q = base.q
    c = sp.Rational(-1, 6630)
    assert sp.LC(q, y) == 2048
    assert q.subs(y, -1) == 3315
    lc_u = c * sp.LC(q, y)
    assert lc_u == sp.Rational(-1024, 3315)
    assert sp.Rational(ce.LC_U) == lc_u
    # Phi = c t^30 q with t monic: same top coefficient.
    phi = sp.expand(c * (y + 1) ** 30 * q)
    assert sp.LC(phi, y) == lc_u
    print("A. lc(q)=2048, q(-1)=3315, lc(u)=lc(Phi)=-1024/3315 = LC_U")


def check_b_backbone() -> None:
    # (2,3,4,5) puts every monomial of h6, h5, h4 on one common degree:
    # the full ties are exactly the P6/P10/P11 supports.
    degstate = (2.0, 3.0, 4.0, 5.0)
    flags = (False, False, False)
    published = {
        6: 14336 * X**2 * D + 8192 * X * E - 3072 * S**2,
        5: (
            -12288 * X**2 * D**2 + 32256 * X**2 * S + 18432 * X * D * E
            - 9216 * D * S**2 + 2048 * E**2
        ),
        4: (
            -220752 * X**4 - 31232 * X**2 * D**3 - 23616 * X**2 * D * S
            - 3072 * X * D**2 * E + 34560 * X * E * S - 5184 * D**2 * S**2
            + 5632 * D * E**2 - 12096 * S**3
        ),
    }
    for level, relation in published.items():
        maximum, labels, count, _ = ce.tropical_h_max_full(
            level, degstate, flags
        )
        assert len(labels) == count, (level, labels)  # full tie
        engine_side = parse_labels(labels)
        # Reference: recompute the tied initial form directly from source.
        rows = ce.MONOMIALS[level]
        values = [
            (ce.monomial_degree(exponents, degstate), exponents, coefficient)
            for exponents, coefficient in rows
        ]
        top = max(value for value, _, _ in values)
        assert top == maximum
        reference = sum(
            int(coefficient) * D**e[0] * X**e[1] * S**e[2] * E**e[3]
            for value, e, coefficient in values
            if value == top
        )
        assert sp.expand(engine_side - reference) == 0, level
        assert sp.expand(engine_side - relation) == 0, level
    print(
        "B. depth-1 infinity tie equations == residue-lemma initial forms "
        "(P6/P10/P11 at degstate (2,3,4,5))"
    )


def check_c_forbidden_drops() -> None:
    flags = (False, False, False)
    cases = [
        # (level, degstate achieving exactly the kill support as full max)
        (5, (2.0, 3.0, 0.0, 5.0)),   # C08: {a,c,e} tie at 10; z-terms low
        (4, (4.0, 5.0, 0.0, 9.0)),   # C20: {b,d,g} tie at 22; 4x=20 below
    ]
    for level, degstate in cases:
        maximum, labels, _, exponent_set = ce.tropical_h_max_full(
            level, degstate, flags
        )
        assert (level, exponent_set) in ce.FORBIDDEN_RISES, (
            level,
            exponent_set,
        )
        saved = ce.APPLY_RESIDUE_KILLS
        try:
            ce.APPLY_RESIDUE_KILLS = False
            options_off = ce.deg_h_options(level, degstate, flags)
            assert len(options_off) > 1  # drops offered without the lemma
            ce.APPLY_RESIDUE_KILLS = True
            options_on = ce.deg_h_options(level, degstate, flags)
            assert options_on == [(maximum, ())], options_on
            assert (
                ce.deg_h_options(level, degstate, flags, required=maximum - 1)
                == []
            )
            assert (
                ce.deg_h_options(level, degstate, flags, required=NEG_INF)
                == []
            )
        finally:
            ce.APPLY_RESIDUE_KILLS = saved
    print(
        "C. C08/C20 kill supports refuse infinity drops under "
        "--residue-kills (leading coefficients live in the base field)"
    )


def check_d_engine_labels() -> None:
    profiles = ce.inf_place_profiles(
        9, "T2", 37, 4, False, False, {6: False, 5: False, 4: False},
        (0.0, NEG_INF, 0.0, 10.0),
    )
    assert profiles
    cancellations = [
        obligation
        for profile in profiles
        for obligation in profile.obligations
        if obligation.kind == "leading_cancellation"
    ]
    assert cancellations
    obligation = cancellations[0]
    assert obligation.level == 5
    g_side, h_side = obligation.tied
    assert g_side == "lc(ehat)^3*lc(g5)"
    prefix = f"({ce.LC_U})^5*["
    assert h_side.startswith(prefix) and h_side.endswith("]"), h_side
    inner = parse_labels(tuple(h_side[len(prefix):-1].split(" + ")))
    assert sp.expand(inner - 2048 * E**2) == 0, inner
    print(
        "D. emitted a9-T2 tie labels reconstruct "
        "lc(ehat)^3*lc(g5) + (-1024/3315)^5*2048*lc(e)^2 = 0"
    )


def main() -> None:
    check_a_constants()
    check_b_backbone()
    check_c_forbidden_drops()
    check_d_engine_labels()
    print("infinity tie equations: PASS")


if __name__ == "__main__":
    main()
