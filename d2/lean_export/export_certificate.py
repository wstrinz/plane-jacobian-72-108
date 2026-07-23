#!/usr/bin/env python3
"""Generate kernel-checkable Lean 4 source for the f37 membership certificate.

MISSION.  Emit an INTEGER-coefficient polynomial identity, equivalent to the
program's headline fact

    f31 = c1*G1 + c2*G2 + c3*G3 + c4*(G5body + Phi)      (exact, over Q)

that a self-contained Lean library can verify by `decide` / `native_decide`.

PROVENANCE / TRUST.  This script does NOT hand-copy any coefficient.  It imports
the EXACT data-loading functions of ``f37_sat_verify.py`` (READ-ONLY) --
``pre_resultant_generators`` (the denominator-cleared integer generators G1,G2,G3,
G5body+Phi parsed from the canonical ``generators.json``), ``load_cofactors`` (the Singular
``lift()`` cofactors c1..c4 from ``f37_sat_certificate.txt``) and ``load_f31``
(``f31_deg31.txt``).  It then:

  1. clears the cofactor denominators: D = lcm of all denominators occurring in
     c1..c4; set  chat_i = D * c_i   (now integer),   L = D * f31 (integer).
     Then  L = sum_i chat_i * g_i  is an exact identity in Z[vars], and because
     D is a nonzero integer the original Q-identity  f31 = sum_i c_i*g_i
     follows by dividing by D.  D is recorded in the Lean statement.

  2. re-verifies  L - sum_i chat_i*g_i == 0  in sympy (independent of Singular).

  3. packs every monomial exponent-vector (8 vars, fixed order) into a single Nat
     key with radix R (asserted > every exponent that occurs in any product term,
     so key addition = exponent-vector addition without carry), and emits the
     polynomials as sorted Lean ``List (Nat x Int)`` literals plus the theorem.

Run:  python export_certificate.py   (writes into ../../lean_certificates/).
"""
from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
D2DIR = HERE.parent                      # d2_plane_72_108
ROOTREPO = D2DIR.parent
OUTDIR = ROOTREPO / "lean_certificates"
LIBDIR = OUTDIR / "Cert"

# Import the *exact* data sources of the reference verifier (never re-copy coeffs).
sys.path.insert(0, str(D2DIR))
import f37_sat_verify as ref  # noqa: E402

VARS = [ref.d2, ref.d1, ref.d0, ref.dm1, ref.dm2, ref.dm3, ref.dm4, ref.Phi]
VAR_NAMES = ["d2", "d1", "d0", "dm1", "dm2", "dm3", "dm4", "Phi"]
NV = len(VARS)


def poly_terms(e: sp.Expr) -> list[tuple[tuple[int, ...], int]]:
    """Expanded integer polynomial -> list of (exponent tuple, int coeff)."""
    p = sp.Poly(sp.expand(e), *VARS)
    out = []
    for monom, coeff in p.terms():
        c = sp.nsimplify(coeff)
        assert c == int(c), f"non-integer coeff {coeff}"
        out.append((tuple(int(x) for x in monom), int(c)))
    return out


def main() -> None:
    print("[1] loading exact sources via f37_sat_verify ...")
    gens = ref.pre_resultant_generators()      # cleared integer G1,G2,G3,G5body+Phi
    cof = ref.load_cofactors()                 # rational cofactors c1..c4
    f31 = ref.load_f31()                       # integer f31

    # --- clear cofactor denominators -------------------------------------------
    print("[2] clearing denominators ...")
    D = sp.Integer(1)
    for c in cof:
        pc = sp.Poly(sp.expand(c), *VARS)
        for coeff in pc.coeffs():
            D = sp.ilcm(D, sp.Rational(coeff).q)
    D = int(D)
    chat = [sp.expand(D * c) for c in cof]      # integer cofactors
    L = sp.expand(D * f31)                      # integer LHS

    # --- independent sympy re-verification of the integer identity -------------
    print(f"[3] verifying integer identity  D*f31 = sum chat_i*g_i   (D={D}) ...")
    combo = sp.expand(sum(ch * g for ch, g in zip(chat, gens)))
    residual = sp.expand(combo - L)
    assert residual == 0, f"INTEGER IDENTITY FAILED: residual {len(sp.Add.make_args(residual))} terms"
    print("    integer identity PASS (residual == 0)")

    # --- convert to term lists -------------------------------------------------
    g_terms = [poly_terms(g) for g in gens]
    c_terms = [poly_terms(c) for c in chat]
    L_terms = poly_terms(L)

    # radix: strictly greater than every per-variable exponent that can occur in
    # any product term chat_i*g_i (before cancellation).  A product monomial's
    # exponent in variable v is at most (max deg_v in chat_i) + (max deg_v in g_i);
    # take the max of that upper bound over i and v.  (Cheap: no big expansion.)
    def maxdeg_per_var(terms):
        m = [0] * NV
        for mono, _ in terms:
            for j in range(NV):
                if mono[j] > m[j]:
                    m[j] = mono[j]
        return m

    maxexp = 0
    for ct, gt in zip(c_terms, g_terms):
        mc, mg = maxdeg_per_var(ct), maxdeg_per_var(gt)
        maxexp = max(maxexp, max(a + b for a, b in zip(mc, mg)))
    for terms in (*g_terms, *c_terms, L_terms):
        maxexp = max(maxexp, max((max(m) for m, _ in terms), default=0))
    R = 1
    while R <= maxexp:
        R *= 2
    # R is a power of two strictly greater than maxexp.
    print(f"[4] max exponent = {maxexp}; radix R = {R}")
    assert maxexp < R

    def key(mono: tuple[int, ...]) -> int:
        k = 0
        for e in reversed(mono):          # mono[0] is lowest-order field
            k = k * R + e
        return k

    def to_pairs(terms):
        pairs = sorted(((key(m), c) for m, c in terms), key=lambda t: t[0])
        # keys are unique already (distinct monomials); assert strict sort
        for a, b in zip(pairs, pairs[1:]):
            assert a[0] < b[0], "duplicate key -- radix too small"
        return pairs

    g_pairs = [to_pairs(t) for t in g_terms]
    c_pairs = [to_pairs(t) for t in c_terms]
    L_pairs = to_pairs(L_terms)

    print("    term counts:")
    for i, gp in enumerate(g_pairs, 1):
        print(f"      G{i}: {len(gp)} terms")
    for i, cp in enumerate(c_pairs, 1):
        print(f"      chat{i}: {len(cp)} terms")
    print(f"      D*f31 (LHS): {len(L_pairs)} terms")

    # --- emit Lean -------------------------------------------------------------
    LIBDIR.mkdir(parents=True, exist_ok=True)

    def coeff_lean(c: int) -> str:
        # Render with the direct Int constructors so the elaborator does NOT run
        # OfNat/Neg typeclass synthesis per element (that is what timed out on the
        # multi-thousand-term cofactor literals).
        return f"Int.ofNat {c}" if c >= 0 else f"Int.negOfNat {-c}"

    def lean_poly_def(name: str, pairs, chunk=400) -> str:
        # emit as chunked ``(c0 ++ c1 ++ ...)`` to avoid one huge cons literal.
        chunks = [pairs[i:i + chunk] for i in range(0, len(pairs), chunk)]
        if not chunks:
            return f"def {name} : Poly := []\n"
        parts = []
        for ch in chunks:
            body = ", ".join(f"({k}, {coeff_lean(c)})" for k, c in ch)
            parts.append("[" + body + "]")
        joined = " ++ ".join(parts)
        return f"def {name} : Poly := {joined}\n"

    print(f"[5] writing Lean data to {LIBDIR / 'Data.lean'} ...")
    data_lines = [
        "-- AUTO-GENERATED by d2_plane_72_108/lean_export/export_certificate.py",
        "-- Do not edit by hand.  Source polynomials of the f37 membership certificate,",
        "-- variable order (packed low->high field): "
        + ", ".join(VAR_NAMES) + f";  radix R = {R}.",
        f"-- Denominator-clearing multiplier D = {D}.",
        "import Cert.Poly",
        "-- large literal data: lift the elaboration heartbeat/recursion limits.",
        "set_option maxHeartbeats 4000000",
        "set_option maxRecDepth 100000",
        "namespace Cert",
        "",
    ]
    for i, gp in enumerate(g_pairs, 1):
        data_lines.append(lean_poly_def(f"G{i}", gp))
    for i, cp in enumerate(c_pairs, 1):
        data_lines.append(lean_poly_def(f"chat{i}", cp))
    data_lines.append(lean_poly_def("Df31", L_pairs))
    data_lines.append("")
    data_lines.append(f"/-- Nonzero integer denominator-clearing multiplier. -/")
    data_lines.append(f"def Dmul : Int := {D}")
    data_lines.append("")
    data_lines.append("end Cert")
    (LIBDIR / "Data.lean").write_text("\n".join(data_lines), encoding="utf8")

    # small metadata for the .md / theorem header
    (OUTDIR / "cert_meta.txt").write_text(
        f"D={D}\nR={R}\nmaxexp={maxexp}\n"
        + "".join(f"G{i}={len(g_pairs[i-1])}\n" for i in range(1, 5))
        + "".join(f"chat{i}={len(c_pairs[i-1])}\n" for i in range(1, 5))
        + f"Df31={len(L_pairs)}\n",
        encoding="utf8",
    )
    print("[done] Data.lean + cert_meta.txt written.")


if __name__ == "__main__":
    main()
