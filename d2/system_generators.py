#!/usr/bin/env python3
"""Shared loader for the pre-resultant generators.

The verification path (f37_sat_verify.py, f37_free_family_verify.py,
lean_export/export_certificate.py, full_system_bridge.py) obtains the four
pre-resultant generators G1,G2,G3,G5body -- plus the derived H2,H3 and sol4 --
from the canonical, human-auditable ``generators.json``.  Nothing on the
mandatory path unpickles anything.

``generators.json`` was emitted ONCE from ``t4_state.pkl`` (see
extract step / regenerate_system.py) and its round-trip was checked exactly.
As an OPTIONAL provenance check, ``verify_pickle_provenance()`` re-loads the
pickle (if present) and confirms the JSON reproduces it; this is skipped
gracefully when the pickle is absent (e.g. a clean public clone)."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
GENERATORS_JSON = ROOT / "generators.json"

d2, d1, d0 = sp.symbols("d2 d1 d0")
dm1, dm2, dm3, dm4 = sp.symbols("dm1 dm2 dm3 dm4")
Phi = sp.symbols("Phi")
VARS = (d2, d1, d0, dm1, dm2, dm3, dm4, Phi)


def _expr_of(terms, var_syms) -> sp.Expr:
    acc = sp.Integer(0)
    for cs, exps in terms:
        c = sp.Rational(cs)
        mono = sp.Integer(1)
        for v, ex in zip(var_syms, exps):
            mono *= v ** ex
        acc += c * mono
    return sp.expand(acc)


def load_generators(path: Path | None = None) -> dict:
    """Parse generators.json -> dict of exact sympy expressions.

    Returns keys G1,G2,G3,G5body,H2,H3 (polynomials) and sol4 (a rational
    function).  No pickle is touched."""
    p = path or GENERATORS_JSON
    data = json.loads(p.read_text(encoding="utf-8"))
    order = data["variable_order"]
    var_syms = [sp.Symbol(n) for n in order]
    out: dict = {}
    for name, terms in data["polynomials"].items():
        out[name] = _expr_of(terms, var_syms)
    num = _expr_of(data["sol4"]["numerator"], var_syms)
    den = _expr_of(data["sol4"]["denominator"], var_syms)
    out["sol4"] = sp.together(num / den)
    out["_provenance"] = data.get("provenance", {})
    return out


def verify_pickle_provenance(gens: dict | None = None) -> str:
    """OPTIONAL provenance: confirm generators.json reproduces t4_state.pkl.

    Returns a human-readable status string.  Never required by the mandatory
    checks; skips gracefully when the pickle is absent (clean public clone)."""
    if gens is None:
        gens = load_generators()
    prov = gens.get("_provenance", {})
    pkl_name = prov.get("source_pickle", "t4_state.pkl")
    pkl_path = ROOT / pkl_name
    if not pkl_path.is_file():
        return (f"    [provenance SKIPPED] {pkl_name} absent; generators.json is the "
                f"self-contained source of record.")
    raw = pkl_path.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    expected = prov.get("source_pickle_sha256")
    tag = "matches recorded sha256" if sha == expected else f"WARNING sha256 {sha} != recorded {expected}"
    # NOTE: unpickling can execute arbitrary code; this optional path is only
    # taken when the trusted local pickle is present.  Deferred import so the
    # module never imports pickle on the mandatory JSON path.
    import pickle  # noqa: PLC0415
    st = pickle.loads(raw)  # trusted local blob, optional provenance only
    for k in ("G1", "G2", "G3", "G5body", "H2", "H3"):
        if sp.expand(gens[k] - sp.expand(st[k])) != 0:
            raise AssertionError(f"provenance MISMATCH on {k}")
    if sp.expand(gens["sol4"] - sp.together(st["sol4"])) != 0:
        raise AssertionError("provenance MISMATCH on sol4")
    return f"    [provenance OK] generators.json reproduces {pkl_name} exactly ({tag})."


if __name__ == "__main__":
    g = load_generators()
    print("loaded generators from", GENERATORS_JSON.name)
    print(verify_pickle_provenance(g))
