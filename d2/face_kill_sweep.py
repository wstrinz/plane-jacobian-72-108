#!/usr/bin/env python3
"""face_kill_sweep.py  --  THE ENGINE BUILD, MILESTONE 1.5

Sweep the bigraded FACE-FUNCTIONAL kill detector across frontier states.

Milestone 1 (`bigrade_annotator.py`) established, on three curated systems, that
an extremal face of the (u-weight, y-order) bigrade lattice can expose a
contradiction that the scalar Groebner projection dissolves: R2 localised a
home-case kill to a single face equation with ZERO spare unknowns, in seconds,
on a 45-var / 122-eq system.  That detector was pointed at three hand-picked
systems.  This module points it at the frontier.

METHOD.  For each entry of `kill_manifest.json` we call
`kill_certificate_msolve.resolve()` purely as a STATE RECONSTRUCTOR: we keep its
reconstructed `(d2, d1, sigma, e)` and its saturation/root data, and DISCARD its
scalar system.  We then re-instantiate the full multi-weight G-system
(G1,G2,G3,G5 at u-weights 13/14/15/17) on that state with stripped spare
ansaetze, exactly as `bigrade_annotator.build_R2` does, and run the face scan.

WHY NOT SWEEP THE MSOLVE SYSTEM DIRECTLY.  The msolve material is a single
graded identity (`f31_master` / `h0_tower`), so every equation carries the same
u-weight; the functional `delta = nu - 12*w` collapses to `nu`, the face
degenerates to "the top-degree coefficient", and the detector reduces to the
classical leading-coefficient argument the engines already exploit.  The face
detector has content only when several generator weights are present.

SOUNDNESS.  Each emitted equation is one y-coefficient of a generator of the
necessary system, so any consequence drawn from a SUBSET of them is a consequence
of the whole system.  The delicate step is `classify_face_value` -- see the long
note there; `bigrade_annotator.left_nullspace_certificate` is deliberately NOT
reused, because its parameter-form branch is unsound as a kill (it declares
"generically-nonzero parameter form" a KILL on symbolic non-vanishing alone,
which ignores that state parameters are constrained by the marked-root relation
and the saturation conditions).

Read-only over every existing artifact.  Writes `face_kill_sweep.json` only.

Usage:
    python face_kill_sweep.py                 # sweep all resolvable entries
    python face_kill_sweep.py --id ID [--id ID ...]
    python face_kill_sweep.py --quiet         # CI-shaped: summary + exit code
"""
from __future__ import annotations

import argparse
import json
import os
import re
import time
import traceback

import sympy as sp

import bigrade_annotator as ba
from bigrade_annotator import System, _gsystem_symbols, \
    _collect_y_equations, _phi_stripped, _place_tag_72, extremal_face, y

_PHI_SYMBOL = _gsystem_symbols()[7]          # the Phi window symbol

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
OUT_PATH = os.path.join(HERE, "face_kill_sweep.json")

# state polynomials arrive as `sp.sstr` text; every identifier in them must be
# pinned to a plain Symbol (see `state_polys_from_material`).
_IDENT_RE = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")

# Stripped spare degree caps.  Inherited from `full_system_bridge.STRIP_DEGCAP`
# -- the audited source of truth, checked by `window_caps_verify.py:306-311`
# against (14-12)k = 2k (sub2) and (15-12)k = 3k (sub1) -- rather than
# recomputed here, so this module cannot drift from the established envelope.
import full_system_bridge as _fsb

SPARE_KS = (6, 7, 8)
SPARE_PREFIX = {6: "R", 7: "S", 8: "T"}
SPARE_NAME = {6: "dm2", 7: "dm3", 8: "dm4"}
STRIP_DEGCAP = _fsb.STRIP_DEGCAP          # {"sub1": {...18,21,24}, "sub2": {...12,14,16}}

# Guards for the OPTIONAL (c)-escalation (`--escalate`).  It runs a Groebner
# call, so both a variable-count cap and a hard wall-clock timeout apply; either
# guard tripping downgrades the verdict to CONSTRAINT, never to KILL.
ESCALATION_MAX_VARS = 8
ESCALATION_TIMEOUT_S = 60

# The five bigrade functionals of Milestone 1, reused verbatim so that any kill
# found here is directly comparable to the R1/R2/R3 report.
FUNCTIONALS = {
    "delta(min)": (lambda P, Q: (lambda w, nu: Q * nu - P * w), "min"),
    "delta(max)": (lambda P, Q: (lambda w, nu: Q * nu - P * w), "max"),
    "w(min)": (lambda P, Q: (lambda w, nu: w), "min"),
    "nu(max)": (lambda P, Q: (lambda w, nu: nu), "max"),
    "mix 2w+3nu(min)": (lambda P, Q: (lambda w, nu: 2 * w + 3 * nu), "min"),
}


class Skip(Exception):
    """Entry cannot be swept soundly; the reason is recorded, never guessed past."""


def canonical_G_generators():
    """G1,G2,G3,G5 with the CANONICAL G5 normalisation `G5 = G5body + Phi`.

    `bigrade_annotator._G_generators` returns `G5 = 2*Phi + G5body`, which is
    NOT the canonical system: the authoritative sources all use `G5body + Phi` --
    `full_system_bridge.py:107` (`st["G5body"] + PHI`, the canonical
    generators.json loader) and `f37_sat_verify.py`, whose load-bearing C11
    membership certificate verifies `f31 == c1*G1 + c2*G2 + c3*G3 +
    c4*(G5body + Phi)`.  The `2*Phi` form was transcribed from
    `FULL_SYSTEM_BRIDGE.md:62`, which contradicts line 114 of the same file
    (`(G5body+Phi)`); line 62 is the erroneous one.

    `2*Phi + G5body` and `Phi + G5body` differ by `Phi`, not by a nonzero
    scalar, so they are genuinely different equations and conclusions do not
    transfer automatically.  (The window-depth kill happens to survive either
    way, since only `deg Phi = 34` and `lc(Phi) != 0` are used -- but the
    emitted certificate VALUE differs, -1024/3315 canonically vs -2048/3315.)

    NOTE (2026-07-24): `bigrade_annotator._G_generators` has since been CORRECTED
    to the canonical form, so this function now only ASSERTS the normalisation
    rather than repairing it.  It is kept as a fail-loud guard: if the `2*Phi`
    form is ever reintroduced upstream, every consumer here stops immediately
    instead of silently emitting certificates for the wrong equation.
    """
    gens = dict(ba._G_generators())
    g5, w5 = gens["G5"]
    coeff = sp.expand(g5).coeff(_PHI_SYMBOL)
    if coeff != 1:
        raise RuntimeError(
            "non-canonical G5: Phi coefficient is %s, expected 1. The canonical "
            "generator is G5body + Phi (full_system_bridge.py, f37_sat_verify.py); "
            "see FACE_KILL_SWEEP.md sec.4." % coeff)
    return gens


# =====================================================================
#  1.  Window resolution  (fail loud rather than guess a cap)
# =====================================================================
def resolve_window(entry):
    """Which envelope window (sub2 / sub1) governs this entry's spare caps?

    Guessing here would silently build the WRONG system -- a wrong cap changes
    the spare ansatz degree and therefore every face -- so anything not pinned
    by provenance is skipped with a reason instead of defaulted.
    """
    recipe = entry["recipe"]
    builder = recipe.get("builder")
    if builder == "harvest_sys4":
        name = recipe["name"]
        if name.startswith("sub2"):
            return "sub2"
        raise Skip("harvest_sys4 name does not pin a window: %s" % name)
    if builder == "harvest_a8":
        return "sub2"                      # batch_convolution_sub2.json
    if builder == "phase_f2_sub2":
        return "sub2"                      # phase_f2_sub2.py
    if builder == "blowup_case":
        if str(recipe.get("case", "")).startswith("sub2_s"):
            return "sub2"
        return "sub1"                      # alt layer -- see note below
    if builder in ("harvest_sys3", "d2_threshold"):
        # The alternate regime (a in 11..15, v = 30-3a < 0) is NOT a separate
        # window: it is a sub-stratum INSIDE subcase (1), so the sub1 column
        # governs.  `a = v_t(e)` is a multiplicity at the place t = y+1, while
        # the window caps are statements about ord at y=0 and total y-degree --
        # logically independent quantities.  Alt states sit AT the sub1 caps
        # with equality (`alt_bridge.py:128-133` asserts deg e = 15 = 3*5,
        # deg sigma = 12 = 3*4, deg d2 = 6 = 3*2), and `ALT_REGIME_INF.md:149`
        # records the caps as "window facts independent of the regime".
        # Using the sub2 column here would be an unsound (too small) cap --
        # exactly the error `BRIDGE_SWEEP.md` sec.3 made and `ALT_BRIDGE.md`
        # overturned.
        return "sub1"
    raise Skip("unknown builder %r" % builder)


# =====================================================================
#  2.  State reconstruction -> multi-weight G-system
# =====================================================================
def state_polys_from_material(material):
    """Recover (d2, d1, sigma, e) as sympy from the material's `polynomials`.

    `kill_certificate_msolve.expr_string` is `sp.sstr(sp.expand(.))`, so the
    round-trip is faithful as long as we sympify in a namespace where `y` and
    every state symbol are plain Symbols (no implicit function application).

    That namespace must be built EXPLICITLY: the state alphabet collides head-on
    with sympy's builtins -- `gamma` is the gamma function, `E` is Euler's
    number, `S` is the singleton registry, `beta`/`zeta`/`Q` likewise.  A bare
    `sympify` silently turns `gamma*(y+1)**8` into a function reference and the
    whole state is nonsense.  So every identifier in the text is pinned to a
    plain Symbol before parsing.
    """
    raw = dict(material.get("polynomials") or {})
    # `phase_f2_scale.reconstruct` names the sigma factor `sig`, while the
    # sub2 builders use `sigma`; normalise before checking completeness.
    if "sigma" not in raw and "sig" in raw:
        raw["sigma"] = raw.pop("sig")
    missing = [k for k in ("d2", "d1", "sigma", "e") if k not in raw]
    if missing:
        # d1 is legitimately absent/zero for T2 states (the branch forces d1==0);
        # a missing sigma or e would mean we did not reconstruct the state.
        hard = [k for k in missing if k != "d1"]
        if hard:
            raise Skip("material lacks state polynomials: %s" % ", ".join(hard))
    texts = [raw.get(k) for k in ("d2", "d1", "sigma", "e")]
    names = set()
    for text in texts:
        if text:
            names |= set(_IDENT_RE.findall(text))
    locals_ns = {name: sp.Symbol(name) for name in names}
    locals_ns["y"] = y

    out = {}
    for key, text in zip(("d2", "d1", "sigma", "e"), texts):
        if text is None:
            out[key] = sp.Integer(0)
            continue
        out[key] = sp.expand(sp.sympify(text, locals=locals_ns))
    return out


def build_state_system(label, polys, window, notes="", material=None):
    """Instantiate G1,G2,G3,G5 on a reconstructed state with stripped spares.

    Mirrors `bigrade_annotator.build_R2` (the validated construction) but takes
    the state from the ledger instead of hard-coding the a8 pilot.
    """
    caps = STRIP_DEGCAP[window]
    sysR = System(label, (12, 1), notes=notes)
    d0, d1, d2, dm1, dm2, dm3, dm4, Phi = _gsystem_symbols()

    d2v, d1v, sigmav, ev = polys["d2"], polys["d1"], polys["sigma"], polys["e"]
    d0v = sp.expand((d2v**2 + sigmav) / 4)

    # every free symbol of the state (bar y) is a PARAMETER, not an unknown:
    # we ask whether the spares exist over a generic point of this state family.
    state_syms = set()
    for expr in (d2v, d1v, sigmav, ev):
        state_syms |= {s for s in expr.free_symbols if s != y}
    for s in sorted(state_syms, key=str):
        sysR.add_param(s, "state param")

    unknown = []
    spare_series = {}
    for k in SPARE_KS:
        cap = caps[SPARE_NAME[k]]
        coeffs = sp.symbols("%s0:%d" % (SPARE_PREFIX[k], cap + 1))
        spare_series[k] = sum(coeffs[i] * y**i for i in range(cap + 1))
        for i, sym in enumerate(coeffs):
            sysR.add_var(sym, k, 12 * k + i, _place_tag_72(k), "optional")
            unknown.append(sym)

    subs = {d0: d0v, d1: d1v, d2: d2v, dm1: ev,
            dm2: spare_series[6], dm3: spare_series[7], dm4: spare_series[8],
            Phi: _phi_stripped()}
    gens = canonical_G_generators()
    sysR.extra["symbol_check"] = [(n, g, w) for n, (g, w) in gens.items()]
    for name, (gexpr, w) in gens.items():
        inst = sp.expand(gexpr.xreplace(subs))
        sysR.eqs.extend(_collect_y_equations(inst, w, name, unknown,
                                             list(sysR.params)))
    sysR.extra["window"] = window
    sysR.extra["spare_caps"] = {SPARE_NAME[k]: caps[SPARE_NAME[k]] for k in SPARE_KS}

    # The state's own parameter constraints, carried for the (c)-escalation:
    # marked-root relations q(r)=0 and the saturation factors that the ledger
    # declares nonzero.  Without these, a parameter form can never be promoted
    # to a kill (see `classify_face_value`).
    sysR.extra["root_variables"] = list((material or {}).get("root_variables") or [])
    sysR.extra["saturation_factors"] = list((material or {}).get("saturation_factors") or [])
    return sysR


# =====================================================================
#  3.  THE HARDENED FACE PREDICATE
# =====================================================================
def classify_face_value(val, system, escalate=False):
    """Decide what a face equation carrying ZERO spare unknowns actually proves.

    Context.  A face equation is one y-coefficient of a generator of the
    necessary system, restricted to the bigrade face exposed by a functional.
    If every spare unknown of the system lies at a strictly non-extremal
    bigrade, those unknowns drop out of the face and the equation degenerates
    to `val = 0`, where `val` is built from STATE PARAMETERS only.

    Three genuinely different situations hide in that one predicate:

      (a) `val` is a nonzero RATIONAL CONSTANT.  Then `val = 0` is false
          outright: the state admits no spares, KILL.  Unconditionally sound.
          (This is what R2 found, and why R2's verdict stands.)

      (b) `val` involves state parameters.  Then `val = 0` is a genuine
          NECESSARY CONDITION on the state -- new information, of exactly the
          kind the bigraded probe extracted for `a4` -- but it is NOT a kill:
          the parameters are free coordinates of the state family, and the
          hypersurface `val = 0` is generally a legal sub-family.
          `bigrade_annotator.left_nullspace_certificate` calls this case a KILL
          ("generically-nonzero parameter form"); that is the unsound branch
          this module exists to replace.

      (c) `val` involves state parameters AND is inconsistent with the state's
          own parameter ideal (the marked-root relation q(r) = 0) together with
          the saturation conditions (gamma != 0, leading coefficients != 0).
          Then it IS a kill -- but only an ideal-membership test can say so.

    Returns (verdict, detail) with verdict in {"KILL", "CONSTRAINT", "VACUOUS"}.
    """
    val = sp.simplify(val)
    if val == 0:
        return "VACUOUS", "face equation is identically satisfied"

    free = {s for s in val.free_symbols if s != y}
    if not free:
        # (a) unconditionally sound kill
        return "KILL", "forced nonzero rational constant %s = 0" % sp.nsimplify(val)

    # (b)/(c) -- parameter form.  DEFAULT POLICY: report as a necessary
    # constraint, never as a kill.  Strictly sound, possibly under-reporting.
    detail = "necessary condition on state params %s" % ", ".join(sorted(map(str, free)))
    if not escalate:
        return "CONSTRAINT", detail

    verdict, why = _escalate_constraint(val, system)
    if verdict == "KILL":
        return "KILL", "constraint incompatible with state ideal + saturation: %s" % why
    return "CONSTRAINT", "%s [escalation: %s]" % (detail, why)


def _escalate_constraint(val, system):
    """(c) Is `val = 0` INCONSISTENT with the state's own parameter constraints?

    The face gives the necessary condition `val = 0`.  The state's parameters are
    not free: they satisfy the marked-root relations `q(r) = 0` and must keep the
    declared saturation factors NONZERO.  So the state family survives this face
    iff the system

        { val = 0 } u { q(r) = 0 : r marked } u { w * prod(saturation) - 1 = 0 }

    has a solution (the last generator is the Rabinowitsch trick encoding
    `prod(saturation) != 0`).  If that ideal is the UNIT ideal, no admissible
    parameter point exists and the whole state family dies.

    This reintroduces a Groebner call -- the very cost the face method avoids --
    which is why it is flag-gated (`--escalate`) and never on the default/CI
    path.  It is nonetheless far smaller than the monolith: a handful of
    parameters rather than 45 spare unknowns.

    UNSOUNDNESS RISK, stated plainly: this is only as good as the modelled
    constraint set.  If the state carries a parameter relation we do NOT include
    here, the ideal can look like the unit ideal when the true (more constrained)
    parameter locus is simply empty for a different reason -- or, worse, we could
    miss a relation that makes a genuine solution invisible.  Adding constraints
    can only make the ideal larger, i.e. more likely to be declared unit, so an
    INCOMPLETE constraint set is the SAFE direction only for the CONSTRAINT
    verdict, not for KILL.  Treat every escalated kill as claimed-not-certified
    until independently audited.
    """
    import kill_certificate_msolve as kcm

    params = sorted({s for s in val.free_symbols if s != y}, key=str)
    root_names = set(system.extra.get("root_variables") or [])
    roots = [s for s in params if str(s) in root_names]

    gens = [sp.expand(val)]
    for rt in roots:
        gens.append(sp.expand(kcm.qpoly(rt)))       # canonical marked-root quartic

    sat_texts = system.extra.get("saturation_factors") or []
    variables = list(params)
    if sat_texts:
        # pin EVERY identifier to a plain Symbol -- `gamma`, `E`, `S` etc. are
        # sympy builtins and would otherwise parse as functions (same trap as
        # `state_polys_from_material`).
        ns = {}
        for text in sat_texts:
            ns.update({n: sp.Symbol(n) for n in _IDENT_RE.findall(text)})
        ns.update({str(s): s for s in params})
        ns["y"] = y
        sat = []
        for text in sat_texts:
            try:
                sat.append(sp.sympify(text, locals=ns))
            except Exception:
                return "CONSTRAINT", "unparseable saturation factor %r" % text[:40]
        product = sp.Integer(1)
        for factor in sat:
            product = product * factor
        wsym = sp.Symbol("w_sat")
        gens.append(sp.expand(wsym * product - 1))
        # saturation factors may involve state parameters absent from `val`;
        # they are genuine unknowns of the consistency question.
        extra_syms = set()
        for expr in sat:
            extra_syms |= {s for s in expr.free_symbols if s != y}
        variables = sorted(set(variables) | extra_syms, key=str) + [wsym]

    if len(variables) > ESCALATION_MAX_VARS:
        return "CONSTRAINT", ("skipped: %d variables exceeds escalation guard %d"
                              % (len(variables), ESCALATION_MAX_VARS))

    try:
        basis = _groebner_with_timeout(gens, variables, ESCALATION_TIMEOUT_S)
    except TimeoutError:
        return "CONSTRAINT", "escalation Groebner timed out at %ds" % ESCALATION_TIMEOUT_S
    except Exception as exc:
        return "CONSTRAINT", "escalation failed: %s" % type(exc).__name__

    if basis == ["1"]:
        return "KILL", "ideal <constraint, q(r), Rabinowitsch(sat)> = <1>"
    return "CONSTRAINT", "state ideal remains consistent (basis size %d)" % len(basis)


def _groebner_basis_strings(gen_strings, var_strings):
    """Worker entry point: Groebner basis, returned as strings (picklable)."""
    variables = [sp.Symbol(v) for v in var_strings]
    ns = {v: s for v, s in zip(var_strings, variables)}
    gens = [sp.sympify(g, locals=ns) for g in gen_strings]
    basis = sp.groebner(gens, *variables, order="grevlex")
    return [sp.sstr(e) for e in basis.exprs]


def _groebner_with_timeout(gens, variables, timeout_s):
    """Run the Groebner call in a worker process so a hang cannot wedge a sweep.

    sympy's `groebner` has no timeout and these ideals are adversarial by
    construction, so an in-process call risks stalling an entire run.
    """
    import concurrent.futures as cf

    gen_strings = [sp.sstr(g) for g in gens]
    var_strings = [str(v) for v in variables]
    with cf.ProcessPoolExecutor(max_workers=1) as pool:
        future = pool.submit(_groebner_basis_strings, gen_strings, var_strings)
        try:
            return future.result(timeout=timeout_s)
        except cf.TimeoutError:
            for proc in pool._processes.values():
                proc.kill()
            raise TimeoutError()


def face_scan(system, escalate=False):
    """Run every functional; classify each face with the hardened predicate."""
    results = {}
    for fname, (mk, which) in FUNCTIONALS.items():
        face = extremal_face(system, mk(system.P, system.Q), which)
        if face is None:
            continue
        unknowns = set(system.unknown_syms())
        free_faces = []
        for e in face["eqs"]:
            if not any(e.expr.has(s) for s in unknowns):
                verdict, detail = classify_face_value(e.expr, system, escalate=escalate)
                free_faces.append({"eq": str(e.label), "verdict": verdict,
                                   "detail": detail,
                                   "value": str(sp.simplify(e.expr))[:200]})
        results[fname] = {
            "value": face["value"],
            "n_eqs": len(face["eqs"]),
            "n_face_vars": len(face["vars"]),
            "unknown_free_equations": free_faces,
        }
    return results


def verdict_of(scan):
    """Roll the per-face classifications up to one state-level verdict."""
    seen = set()
    for fc in scan.values():
        for row in fc["unknown_free_equations"]:
            seen.add(row["verdict"])
    if "KILL" in seen:
        return "KILL"
    if "CONSTRAINT" in seen:
        return "CONSTRAINT"
    return "NO-FACE-CERTIFICATE"


# =====================================================================
#  4.  Sweep driver
# =====================================================================
def sweep(ids=None, verbose=True, escalate=False):
    import kill_certificate_msolve as kcm

    manifest = json.loads(open(kcm.MANIFEST_PATH, encoding="utf-8-sig").read())
    entries = manifest["entries"]
    if ids:
        wanted = set(ids)
        entries = [e for e in entries if e["id"] in wanted]
        missing = wanted - {e["id"] for e in entries}
        if missing:
            raise SystemExit("unknown manifest ids: " + ", ".join(sorted(missing)))

    rows = []
    for entry in entries:
        eid = entry["id"]
        row = {"id": eid, "category": entry.get("category"),
               "builder": entry["recipe"].get("builder")}
        t0 = time.time()
        try:
            window = resolve_window(entry)
            _members, material = kcm.resolve(entry)
            polys = state_polys_from_material(material)
            system = build_state_system(eid, polys, window,
                                        notes="reconstructed from kill_manifest",
                                        material=material)
            viol, checked = ba.check_bigrade_consistency(system)
            if viol:
                raise Skip("bigrade consistency FAILED (%d violations)" % len(viol))
            scan = face_scan(system, escalate=escalate)
            row.update({
                "status": "SWEPT", "window": window,
                "n_vars": len(system.vars), "n_eqs": len(system.eqs),
                "bigrade_terms_checked": checked,
                "spare_caps": {str(k): v for k, v in system.extra["spare_caps"].items()},
                "faces": scan, "verdict": verdict_of(scan),
            })
        except Skip as exc:
            row.update({"status": "SKIPPED", "reason": str(exc)})
        except Exception as exc:                      # loud, never silent
            row.update({"status": "ERROR", "reason": "%s: %s" % (type(exc).__name__, exc),
                        "traceback": traceback.format_exc()[-800:]})
        row["seconds"] = round(time.time() - t0, 2)
        rows.append(row)
        if verbose:
            tag = row.get("verdict") or row["status"]
            extra = "" if row["status"] == "SWEPT" else "  (%s)" % row.get("reason", "")[:70]
            print("  %-46s %-20s %6.1fs%s" % (eid[:46], tag, row["seconds"], extra))

    census = {}
    for r in rows:
        key = r.get("verdict") or r["status"]
        census[key] = census.get(key, 0) + 1
    return {"schema": "d2-face-kill-sweep-v1",
            "manifest": os.path.basename(str(kcm.MANIFEST_PATH)),
            "strip_degcap": STRIP_DEGCAP,
            "census": census, "rows": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", action="append", dest="ids")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--escalate", action="store_true",
                    help="promote parameter-form CONSTRAINTs to KILLs when they are "
                         "inconsistent with the state ideal + saturation (runs a "
                         "guarded Groebner call; OFF by default and never in CI -- "
                         "escalated kills are claimed, not certified)")
    args = ap.parse_args()

    if not args.quiet:
        print("=" * 78)
        print("FACE-FUNCTIONAL KILL SWEEP  (Milestone 1.5)")
        print("=" * 78)
    out = sweep(ids=args.ids, verbose=not args.quiet, escalate=args.escalate)
    out["escalation_enabled"] = bool(args.escalate)
    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)

    print("\nCENSUS:", ", ".join("%s=%d" % kv for kv in sorted(out["census"].items())))
    kills = [r["id"] for r in out["rows"] if r.get("verdict") == "KILL"]
    cons = [r["id"] for r in out["rows"] if r.get("verdict") == "CONSTRAINT"]
    if kills:
        print("KILLS:", ", ".join(kills))
    if cons:
        print("NEW NECESSARY CONSTRAINTS:", ", ".join(cons))
    print("wrote", os.path.basename(OUT_PATH))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
