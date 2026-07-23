#!/usr/bin/env python3
"""Lower-cascade valuation engine — Phase B of CASCADE_ENGINE_PLAN.md.

Propagates exact local valuation states at the four split places (roots of
q) down through cascade levels 6, 5, 4, coupled by global degree budgets,
starting from every open branch of split_place_ledger.json.

Algebraic ground truth (verified in t5_multiplace_verify.py, checks 5-7):

    t^v * g_{l+1} = ehat^3 * g_l + u^l * h_l,      u = c*q,  v = 30-3a,
    terminal T1:  ehat^3 * g_7 = -u^7 * h_7,       h_7 = 8192*d1^2,
    terminal T2:  ehat^3 * g_6 = -u^6 * h_6|_{d1=0},  h_6|_{d1=0} = -3072*sigma^2,
    deg g_l <= 10+3a,  deg d2 <= 4, deg d1 <= 6, deg sigma <= 8, deg e <= 10.

At a root p of q: t is a unit, v_p(u) = 1, v_p(e) = b exactly.  Monomial
valuations of h_l under a local state are therefore exact.  Ultrametric
semantics: a sum of terms with a unique minimum valuation has exactly that
valuation; a tie permits a rise only through residue cancellation, which is
recorded as an obligation, never silently granted.  h_l may vanish
identically only when at least two monomials survive the zero flags, and
that too is recorded as an obligation.

Soundness direction: the engine over-approximates the solution set.  A
branch reported killed has NO consistent valuation profile under exact
necessary conditions; a surviving branch carries explicit residue
obligations for Phase C.  New kills beyond the terminal layer are labelled
pending independent audit.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path

from cascade_signature import load_levels

ROOT = Path(__file__).resolve().parent
LEDGER_PATH = ROOT / "split_place_ledger.json"
JSON_OUT = ROOT / "cascade_cones.json"

# Exponent order matches cascade_signature.SIGNATURE_VARIABLES.
VAR_NAMES = ("d2", "d1", "sigma", "e")
GLOBAL_CAPS = {"d2": 4, "d1": 6, "sigma": 8}
INF = float("inf")
NEG_INF = float("-inf")

T1_TERMINAL, T2_TERMINAL = 7, 6
DEG_U = 4  # deg u = deg(c*q) = deg q
# Leading coefficient of u = c*q (equal to Phi's top coefficient, since
# Phi = c*t^30*q and t is monic): c*lc(q) = (-1/6630)*2048 = -1024/3315.
# Verified against the source q in cascade_inf_ties_verify.py.
LC_U = "-1024/3315"

# Proven forbidden rises (RESIDUE_LEMMAS.md, verified by
# residue_lemmas_verify.py): when the tied minimum-achiever set at a place
# is exactly one of these (level, exponent-set) shapes, the depth-1 residue
# equation has NO solution with all leading coefficients nonzero over Q or
# the q-splitting field, so neither a tie rise nor identical vanishing can
# occur there.  Exponent order matches VAR_NAMES = (d2, d1, sigma, e).
# C08: level 5, {d1^2 d2^2, d1 d2 e, e^2}: 6X^2D^2 - 9XDE - E^2 = 0 has
#      discriminant square class 105 (not a square in the S4 splitting
#      field, abelianization C2).
# C20: level 4, {d1^2 d2^3, d1 d2^2 e, d2 e^2}: 61X^2D^2 + 6XDE - 11E^2,
#      square class 170.
FORBIDDEN_RISES = frozenset(
    {
        (5, frozenset({(2, 2, 0, 0), (1, 1, 0, 1), (0, 0, 0, 2)})),
        (4, frozenset({(3, 2, 0, 0), (2, 1, 0, 1), (1, 0, 0, 2)})),
    }
)
# Default False so the audited artifacts remain reproducible; the
# lemma-assisted sweep enables it via --residue-kills.
APPLY_RESIDUE_KILLS = False


def load_monomials() -> dict[int, list[tuple[tuple[int, int, int, int], str]]]:
    """Exact (exponent, coefficient) tables per level, source-linked."""

    table = {}
    for index, level in load_levels().items():
        rows = []
        for record in level.monomial_records():
            exponents = tuple(record["exponents"][name] for name in VAR_NAMES)
            rows.append((exponents, record["coefficient"]))
        table[index] = rows
    return table


MONOMIALS = load_monomials()
H_DEGREE_CAP = {index: 40 - 4 * index for index in range(8)}


@dataclass(frozen=True)
class WindowConfig:
    """Window-dependent degree data; the h_f monomials are shared.

    sub2: proven caps (T3_WINDOW_AUDIT.md, t5_multiplace_verify.py check 7)
      aux (d1,sigma,d2) = (6,8,4), deg h_f <= 40-4f, uniform g-cap 10+3a
      (bottom-anchored induction closes exactly: 3(10-a)+(10+3a)=40).
    sub1: verified caps (sub1_cascade_verify.py) — NO uniform g-cap exists;
      per-level caps are min(forward, backward) recursions with terminal
      anchors deg g7 <= 46 (T1) / deg g6 <= 48 (T2).
    """

    name: str
    aux_caps: tuple[int, int, int]  # (d1, sigma, d2)
    e_cap: int  # total deg e cap; ledger budget a + sum(b) <= e_cap
    h_slope: int  # deg h_f <= h_slope * (20 - 2f) / 2 ... see h_cap

    def h_cap(self, level: int) -> int:
        # weighted-homogeneous of weight 20-2f with deg d~ <= m*weight
        return self.h_slope * (20 - 2 * level) // 2

    def g_caps(self, a: int, branch: str) -> dict[int, int]:
        if self.name == "sub2":
            return {level: 10 + 3 * a for level in range(1, 8)}
        v = 30 - 3 * a
        assert v >= 0, "standard regime only"
        ehat_cap = self.e_cap - a
        forward = {1: 60 - v}
        for level in range(1, 7):
            forward[level + 1] = max(
                3 * ehat_cap + forward[level], 60 - 2 * level
            ) - v
        backward: dict[int, int] = {}
        if branch == "T1":
            backward[7] = 46
            top = 6
        else:
            backward[7] = 46  # unused for T2 (g_7 == 0); kept for shape
            backward[6] = 48
            top = 5
        for level in range(top, 0, -1):
            backward[level] = max(v + backward[level + 1], 60 - 2 * level)
        return {
            level: min(forward[level], backward[level]) for level in range(1, 8)
        }


SUB2 = WindowConfig(name="sub2", aux_caps=(6, 8, 4), e_cap=10, h_slope=4)
SUB1 = WindowConfig(name="sub1", aux_caps=(9, 12, 6), e_cap=15, h_slope=6)
CONFIGS = {"sub2": SUB2, "sub1": SUB1}


@dataclass(frozen=True)
class Obligation:
    """A necessary residue condition attached to a surviving profile."""

    level: int
    kind: str  # "monomial_tie_rise" | "term_cancellation" | "exact_identity"
    #           | "identical_vanishing"
    depth: int
    tied: tuple[str, ...]


@dataclass(frozen=True)
class PlaceProfile:
    """One consistent local valuation chain at a single split place."""

    b: int  # e-slot valuation at this place: b_i at a q-root, a at t
    x: float  # v(d1); INF encodes d1 == 0 (branch T2)
    z: float  # v(sigma); INF encodes the sigma == 0 flag
    k: float  # v(d2); INF encodes the d2 == 0 flag
    r: tuple[tuple[int, float], ...]  # (level, v(g_level)); INF for g == 0
    obligations: tuple[Obligation, ...]
    place: str = "q"

    def budget_vector(self, levels: tuple[int, ...]) -> tuple[float, ...]:
        r_map = dict(self.r)
        return (
            0 if self.x == INF else self.x,
            0 if self.z == INF else self.z,
            0 if self.k == INF else self.k,
            *(0 if r_map[level] == INF else r_map[level] for level in levels),
        )


def specialized_monomials(
    level: int, flags: tuple[bool, bool, bool]
) -> list[tuple[tuple[int, int, int, int], str]]:
    """Drop monomials containing an identically-zero variable."""

    sigma_zero, d2_zero, d1_zero = flags
    kept = []
    for exponents, coefficient in MONOMIALS[level]:
        k_exp, x_exp, z_exp, _ = exponents
        if d1_zero and x_exp:
            continue
        if sigma_zero and z_exp:
            continue
        if d2_zero and k_exp:
            continue
        kept.append((exponents, coefficient))
    return kept


def monomial_valuation(
    exponents: tuple[int, int, int, int], k: float, x: float, z: float, b: int
) -> float:
    k_exp, x_exp, z_exp, e_exp = exponents
    total = e_exp * b
    for exponent, value in ((k_exp, k), (x_exp, x), (z_exp, z)):
        if exponent:
            if value == INF:
                return INF
            total += exponent * value
    return total


def tropical_h_full(
    level: int,
    state: tuple[float, float, float, int],
    flags: tuple[bool, bool, bool],
) -> tuple[float, tuple[str, ...], int, frozenset]:
    """Minimum valuation of h_level, achiever labels/exponents, term count."""

    k, x, z, b = state
    rows = specialized_monomials(level, flags)
    if not rows:
        return INF, (), 0, frozenset()
    best, achievers = INF, []
    for exponents, coefficient in rows:
        value = monomial_valuation(exponents, k, x, z, b)
        if value < best:
            best, achievers = value, [(exponents, coefficient)]
        elif value == best and value != INF:
            achievers.append((exponents, coefficient))
    if best == INF:
        return INF, (), len(rows), frozenset()
    labels = tuple(
        f"{coefficient}*d2^{e[0]}*d1^{e[1]}*sigma^{e[2]}*e^{e[3]}"
        for e, coefficient in achievers
    )
    exponent_set = frozenset(exponents for exponents, _ in achievers)
    return best, labels, len(rows), exponent_set


def tropical_h(
    level: int,
    state: tuple[float, float, float, int],
    flags: tuple[bool, bool, bool],
) -> tuple[float, tuple[str, ...], int]:
    """Exact minimum valuation of h_level, achievers, surviving-term count."""

    best, labels, count, _ = tropical_h_full(level, state, flags)
    return best, labels, count


def w_options(
    level: int,
    state: tuple[float, float, float, int],
    flags: tuple[bool, bool, bool],
    required: float | None = None,
    h_cap: int | None = None,
) -> list[tuple[float, tuple[Obligation, ...]]]:
    """Possible exact valuations w = v_p(h_level) with their obligations.

    Unique surviving monomial: w is forced to the monomial valuation and the
    polynomial cannot vanish identically.  Several monomials: w may exceed
    the tropical minimum only through residue cancellation among the tied
    minimum monomials, bounded per place by deg h_level; w = INF (h == 0
    identically) is possible only via total cancellation.  All granted rises
    carry obligations.  ``required`` restricts to a single value.
    """

    minimum, tied, term_count, exponent_set = tropical_h_full(
        level, state, flags
    )
    if minimum == INF:
        # Every surviving monomial already has infinite valuation (or none
        # survive): h_level vanishes at this place to any order.
        return [(INF, ())] if required in (None, INF) else []

    options: list[tuple[float, tuple[Obligation, ...]]] = []
    unique = len(tied) == 1
    cap = H_DEGREE_CAP[level] if h_cap is None else h_cap
    rise_forbidden = (
        APPLY_RESIDUE_KILLS and (level, exponent_set) in FORBIDDEN_RISES
    )

    def rise(w: float) -> tuple[Obligation, ...] | None:
        if w == minimum:
            return ()
        if unique or term_count < 2 or rise_forbidden:
            return None
        if w == INF:
            return (
                Obligation(level, "identical_vanishing", 0, tied),
            )
        if w > cap:
            return None
        return (
            Obligation(level, "monomial_tie_rise", int(w - minimum), tied),
        )

    if required is not None:
        if required < minimum:
            return []
        obligations = rise(required)
        return [] if obligations is None else [(required, obligations)]

    for w in range(int(minimum), cap + 1):
        obligations = rise(w)
        if obligations is not None:
            options.append((w, obligations))
    obligations = rise(INF)
    if obligations is not None:
        options.append((INF, obligations))
    return options


def descend_options(
    level: int,
    r_above: float,
    g_above_zero: bool,
    g_zero: bool,
    state: tuple[float, float, float, int],
    flags: tuple[bool, bool, bool],
    r_cap: int,
    g_shift: float | None = None,
    h_shift: float | None = None,
    h_cap: int | None = None,
) -> list[tuple[float, tuple[Obligation, ...]]]:
    """Solve t^v g_{l+1} = ehat^3 g_l + u^l h_l for the g_l valuation.

    Returns possible (r_l, obligations) pairs given that the LEFT side has
    valuation ``r_above`` at this place.  The two right-side candidates are
    ``g_shift + r_l`` and ``h_shift + w`` with w the h_l valuation.

    At a root of q: g_shift = 3b (from ehat^3), h_shift = level (from u^l),
    and r_above = v_p(g_{l+1}) since t is a unit — the defaults.
    At the place t = y+1: ehat and u are units (q(-1) = 3315 != 0, checked
    in t5_multiplace_verify.py), so g_shift = h_shift = 0, the caller
    passes r_above = v + s_{l+1} with v = 30-3a, and the state's e-slot
    valuation is a = v_t(e).

    Zero flags are global statements: g == 0 means the polynomial vanishes
    identically, giving exact identities instead of ultrametric minima.
    """

    b = state[3]
    if g_shift is None:
        g_shift = 3 * b
    if h_shift is None:
        h_shift = level
    results: list[tuple[float, tuple[Obligation, ...]]] = []

    if g_zero and g_above_zero:
        # 0 = 0 + u^l h_l: h_l must vanish identically.
        return [(INF, obligations) for _, obligations in w_options(
            level, state, flags, required=INF, h_cap=h_cap
        )]

    if g_zero:
        # t^v g_{l+1} = u^l h_l exactly.
        required = r_above - h_shift
        if required < 0:
            return []
        return [
            (INF, obligations)
            for _, obligations in w_options(
                level, state, flags, required=required, h_cap=h_cap
            )
        ]

    if g_above_zero:
        # Polynomial identity ehat^3 g_l = -u^l h_l: valuations match exactly.
        identity = Obligation(level, "exact_identity", 0, ())
        for r_l in range(0, r_cap + 1):
            required = g_shift + r_l - h_shift
            if required < 0:
                continue
            for _, obligations in w_options(
                level, state, flags, required=required, h_cap=h_cap
            ):
                results.append((r_l, obligations + (identity,)))
        return results

    if r_above == INF:
        raise ValueError("finite g_{l+1} flag with infinite valuation")

    options = w_options(level, state, flags, h_cap=h_cap)

    # Case (a): the g_l term dominates — g_shift + r_l = r_above < h_shift + w.
    r_l = r_above - g_shift
    if 0 <= r_l <= r_cap:
        for w, obligations in options:
            if w == INF or h_shift + w > r_above:
                results.append((r_l, obligations))
                break  # options are ordered; first admissible w is minimal

    # Case (b): the h_l term dominates — h_shift + w = r_above < g_shift + r_l.
    required = r_above - h_shift
    if required >= 0:
        for _, obligations in w_options(
            level, state, flags, required=required, h_cap=h_cap
        ):
            low = max(0, int(r_above - g_shift) + 1)
            for r_l in range(low, r_cap + 1):
                results.append((r_l, obligations))

    # Case (c): tie g_shift + r_l = h_shift + w <= r_above, rise by cancellation.
    for w, obligations in options:
        if w == INF:
            continue
        r_l = h_shift + w - g_shift
        if r_l < 0 or r_l > r_cap:
            continue
        tie_value = h_shift + w
        if tie_value > r_above:
            continue
        depth = int(r_above - tie_value)
        extra = (
            (Obligation(level, "term_cancellation", depth, ()),)
            if depth > 0
            else ()
        )
        candidate = (float(r_l), obligations + extra)
        if candidate not in results:
            results.append(candidate)

    return results


def pareto_minimal(
    items: list[tuple[tuple[float, ...], object]]
) -> list[tuple[tuple[float, ...], object]]:
    """Keep vector-minimal items (sound: cross-place coupling is <= caps)."""

    kept: list[tuple[tuple[float, ...], object]] = []
    for vector, payload in sorted(items, key=lambda item: (sum(item[0]), item[0])):
        if not any(
            all(previous[i] <= vector[i] for i in range(len(vector)))
            for previous, _ in kept
        ):
            kept.append((vector, payload))
    return kept


class ParetoAccumulator:
    """Bounded-memory Pareto frontier.

    Compacts the buffer whenever it grows past ``compact_at``.  Staged
    compaction is exact: Pareto(Pareto(A) + B) == Pareto(A + B), since an
    element dominated inside a chunk is dominated in the union and a true
    frontier element can never be dominated within its own chunk.  This
    bounds memory by the frontier size instead of the raw enumeration
    (the sub1 caps made the raw enumeration blow past physical RAM).
    """

    def __init__(self, compact_at: int = 20000) -> None:
        self.items: list[tuple[tuple[float, ...], object]] = []
        self.compact_at = compact_at

    def add(self, vector: tuple[float, ...], payload: object) -> None:
        self.items.append((vector, payload))
        if len(self.items) >= self.compact_at:
            self.items = pareto_minimal(self.items)

    def result(self) -> list[tuple[tuple[float, ...], object]]:
        return pareto_minimal(self.items)


# Cross-call memoization (the lru_cache was lost when g_zero became a
# dict; without it a ledger sweep recomputes every place set per stratum
# and flag case).  Keyed on canonical hashable forms; cleared between
# sweep chunks by clear_profile_caches().
_PROFILE_CACHE: dict = {}


def clear_profile_caches() -> None:
    _PROFILE_CACHE.clear()
    _INF_SIG_CACHE.clear()


def resolve_caps(
    branch: str, r_cap: int, config: WindowConfig | None, a: int | None
) -> tuple[dict[str, int], dict[int, int], dict[int, int | None]]:
    """(aux caps, per-level g caps, per-level h caps) for a place run.

    Default (config None or sub2): the audited uniform-cap behavior, with
    g cap = r_cap at every level and the module H_DEGREE_CAP.
    """

    if config is None or config.name == "sub2":
        return (
            dict(GLOBAL_CAPS),
            {level: r_cap for level in range(1, 8)},
            {level: None for level in range(1, 8)},
        )
    assert a is not None, "sub1 caps require the t-multiplicity a"
    d1_cap, sigma_cap, d2_cap = config.aux_caps
    return (
        {"d1": d1_cap, "sigma": sigma_cap, "d2": d2_cap},
        config.g_caps(a, branch),
        {level: config.h_cap(level) for level in range(1, 8)},
    )


def place_profiles(
    b: int,
    branch: str,
    r_cap: int,
    depth: int,
    sigma_zero: bool,
    d2_zero: bool,
    g_zero: dict[int, bool],
    config: WindowConfig | None = None,
    a: int | None = None,
) -> tuple[PlaceProfile, ...]:
    """All consistent local chains at one place, Pareto-reduced."""

    cache_key = (
        "q",
        APPLY_RESIDUE_KILLS,
        b,
        branch,
        r_cap,
        depth,
        sigma_zero,
        d2_zero,
        tuple(sorted(g_zero.items())),
        None if config is None else config.name,
        a,
    )
    if cache_key in _PROFILE_CACHE:
        return _PROFILE_CACHE[cache_key]

    d1_zero = branch == "T2"
    flags = (sigma_zero, d2_zero, d1_zero)
    terminal = T1_TERMINAL if branch == "T1" else T2_TERMINAL
    levels = tuple(range(terminal, depth - 1, -1))
    aux, g_caps, h_caps = resolve_caps(branch, r_cap, config, a)

    x_range: tuple[float, ...] = (
        (INF,) if d1_zero else tuple(range(aux["d1"] + 1))
    )
    z_range: tuple[float, ...] = (
        (INF,) if sigma_zero else tuple(range(aux["sigma"] + 1))
    )
    k_range: tuple[float, ...] = (
        (INF,) if d2_zero else tuple(range(aux["d2"] + 1))
    )

    accumulator = ParetoAccumulator()
    for x in x_range:
        for z in z_range:
            for k in k_range:
                state = (k, x, z, b)
                memo: dict[
                    tuple[int, float],
                    list[tuple[tuple[tuple[int, float], ...], tuple[Obligation, ...]]],
                ] = {}

                def tails(level_index: int, r_above: float):
                    """Completed descent tails below levels[level_index-1]."""

                    if level_index == len(levels):
                        return [((), ())]
                    key = (level_index, r_above)
                    if key in memo:
                        return memo[key]
                    level = levels[level_index]
                    collected = ParetoAccumulator(compact_at=20000)
                    for r_l, obligations in descend_options(
                        level,
                        r_above,
                        g_zero.get(level + 1, False),
                        g_zero.get(level, False),
                        state,
                        flags,
                        g_caps[level],
                        h_cap=h_caps[level],
                    ):
                        for tail, tail_obligations in tails(level_index + 1, r_l):
                            chain = ((level, r_l),) + tail
                            total = obligations + tail_obligations
                            collected.add(
                                tuple(
                                    0 if value == INF else value
                                    for _, value in chain
                                )
                                + (len(total),),
                                (chain, total),
                            )
                    memo[key] = [payload for _, payload in collected.result()]
                    return memo[key]

                # Terminal identity: ehat^3 g_T = -u^T h_T exactly.
                terminal_zero = g_zero.get(terminal, False)
                start_values = (
                    [INF] if terminal_zero else list(range(g_caps[terminal] + 1))
                )
                for r_t in start_values:
                    required = INF if r_t == INF else 3 * b + r_t - terminal
                    if required != INF and required < 0:
                        continue
                    starts = w_options(
                        terminal,
                        state,
                        flags,
                        required=required,
                        h_cap=h_caps[terminal],
                    )
                    if not starts:
                        continue
                    _, start_obligations = starts[0]
                    for tail, tail_obligations in tails(1, r_t):
                        profile = PlaceProfile(
                            b=b,
                            x=x,
                            z=z,
                            k=k,
                            r=((terminal, float(r_t)),) + tail,
                            obligations=start_obligations + tail_obligations,
                        )
                        accumulator.add(
                            profile.budget_vector(levels)
                            + (len(profile.obligations),),
                            profile,
                        )

    result = tuple(payload for _, payload in accumulator.result())
    _PROFILE_CACHE[cache_key] = result
    return result


def t_place_profiles(
    a: int,
    branch: str,
    r_cap: int,
    depth: int,
    sigma_zero: bool,
    d2_zero: bool,
    g_zero: dict[int, bool],
    config: WindowConfig | None = None,
) -> tuple[PlaceProfile, ...]:
    """All consistent chains at the place t = y+1, Pareto-reduced.

    At t: v_t(t^v) = v = 30-3a, v_t(ehat) = 0 (a = v_t(e) exactly), and
    v_t(u) = 0 since q(-1) = 3315 != 0.  The e-slot of every h_l monomial
    therefore costs a, the level identity reads
        v + s_{l+1} = ultrametric(s_l, w_l),
    and the terminals are s_7 = w_7 (T1) resp. s_6 = w_6 (T2).
    Requires the standard regime a <= 10 (v >= 0).
    """

    if a > 10:
        raise ValueError("t-place transition requires v = 30-3a >= 0")
    cache_key = (
        "t",
        APPLY_RESIDUE_KILLS,
        a,
        branch,
        r_cap,
        depth,
        sigma_zero,
        d2_zero,
        tuple(sorted(g_zero.items())),
        None if config is None else config.name,
    )
    if cache_key in _PROFILE_CACHE:
        return _PROFILE_CACHE[cache_key]
    v = 30 - 3 * a
    d1_zero = branch == "T2"
    flags = (sigma_zero, d2_zero, d1_zero)
    terminal = T1_TERMINAL if branch == "T1" else T2_TERMINAL
    levels = tuple(range(terminal, depth - 1, -1))
    aux, g_caps, h_caps = resolve_caps(branch, r_cap, config, a)

    x_range: tuple[float, ...] = (
        (INF,) if d1_zero else tuple(range(aux["d1"] + 1))
    )
    z_range: tuple[float, ...] = (
        (INF,) if sigma_zero else tuple(range(aux["sigma"] + 1))
    )
    k_range: tuple[float, ...] = (
        (INF,) if d2_zero else tuple(range(aux["d2"] + 1))
    )

    accumulator = ParetoAccumulator()
    for x in x_range:
        for z in z_range:
            for k in k_range:
                state = (k, x, z, a)
                memo: dict = {}

                def tails(level_index: int, s_above: float):
                    if level_index == len(levels):
                        return [((), ())]
                    key = (level_index, s_above)
                    if key in memo:
                        return memo[key]
                    level = levels[level_index]
                    collected = ParetoAccumulator(compact_at=20000)
                    r_above = INF if s_above == INF else v + s_above
                    for s_l, obligations in descend_options(
                        level,
                        r_above,
                        g_zero.get(level + 1, False),
                        g_zero.get(level, False),
                        state,
                        flags,
                        g_caps[level],
                        g_shift=0,
                        h_shift=0,
                        h_cap=h_caps[level],
                    ):
                        for tail, tail_obligations in tails(level_index + 1, s_l):
                            chain = ((level, s_l),) + tail
                            total = obligations + tail_obligations
                            collected.add(
                                tuple(
                                    0 if value == INF else value
                                    for _, value in chain
                                )
                                + (len(total),),
                                (chain, total),
                            )
                    memo[key] = [payload for _, payload in collected.result()]
                    return memo[key]

                # Terminal identity ehat^3 g_T = -u^T h_T at t: s_T = w_T.
                terminal_zero = g_zero.get(terminal, False)
                start_values = (
                    [INF] if terminal_zero else list(range(g_caps[terminal] + 1))
                )
                for s_t in start_values:
                    required = INF if s_t == INF else float(s_t)
                    starts = w_options(
                        terminal,
                        state,
                        flags,
                        required=required,
                        h_cap=h_caps[terminal],
                    )
                    if not starts:
                        continue
                    _, start_obligations = starts[0]
                    for tail, tail_obligations in tails(1, s_t):
                        profile = PlaceProfile(
                            b=a,
                            x=x,
                            z=z,
                            k=k,
                            r=((terminal, float(s_t)),) + tail,
                            obligations=start_obligations + tail_obligations,
                            place="t",
                        )
                        accumulator.add(
                            profile.budget_vector(levels)
                            + (len(profile.obligations),),
                            profile,
                        )

    result = tuple(payload for _, payload in accumulator.result())
    _PROFILE_CACHE[cache_key] = result
    return result


# ---------------------------------------------------------------------------
# Phase D: infinity as a sixth place (max-plus dual of the descent).
#
# With v_inf = -deg, the level identity t^v g_{l+1} = ehat^3 g_l + u^l h_l
# becomes, in degrees,
#     v + deg g_{l+1}  vs  max(3 deg ehat + deg g_l, 4l + deg h_l):
# a unique maximum forces the degree exactly; a tie permits the degree to
# DROP only through leading-coefficient cancellation (the dual of a residue
# tie-rise), recorded as an obligation with the tied leading monomials.
# Monomial degrees are exact given per-variable degree assignments
# (deg d2, deg d1, deg sigma, deg e); h_l is evaluated at the ORIGINAL
# variables, so the e-slot costs deg e while ehat^3 shifts by 3(deg e - a).
# Unlike the finite places the chain is run all the way down: levels
# terminal..1 plus the level-0 anchor t^v g_1 = h_0 (t5_multiplace_verify.py
# check 5), because the master-sum degree kills need the f=0 term.  NEG_INF
# encodes the degree of an identically-zero polynomial.


@dataclass(frozen=True)
class InfProfile:
    """One consistent degree chain at the place at infinity."""

    degs: tuple[float, float, float, float]  # (deg d2, deg d1, deg sigma, deg e)
    chain: tuple[tuple[int, float], ...]  # (level, deg g_level), terminal..1
    obligations: tuple[Obligation, ...]
    place: str = "inf"


def monomial_degree(
    exponents: tuple[int, int, int, int],
    degstate: tuple[float, float, float, float],
) -> float:
    k_exp, x_exp, z_exp, e_exp = exponents
    deg_k, deg_x, deg_z, deg_e = degstate
    total = e_exp * deg_e
    for exponent, value in ((k_exp, deg_k), (x_exp, deg_x), (z_exp, deg_z)):
        if exponent:
            if value == NEG_INF:
                return NEG_INF
            total += exponent * value
    return total


def tropical_h_max_full(
    level: int,
    degstate: tuple[float, float, float, float],
    flags: tuple[bool, bool, bool],
) -> tuple[float, tuple[str, ...], int, frozenset]:
    """Maximum degree of h_level, achiever labels/exponents, term count."""

    rows = specialized_monomials(level, flags)
    if not rows:
        return NEG_INF, (), 0, frozenset()
    best, achievers = NEG_INF, []
    for exponents, coefficient in rows:
        value = monomial_degree(exponents, degstate)
        if value > best:
            best, achievers = value, [(exponents, coefficient)]
        elif value == best and value != NEG_INF:
            achievers.append((exponents, coefficient))
    if best == NEG_INF:
        return NEG_INF, (), len(rows), frozenset()
    labels = tuple(
        f"{coefficient}*d2^{e[0]}*d1^{e[1]}*sigma^{e[2]}*e^{e[3]}"
        for e, coefficient in achievers
    )
    exponent_set = frozenset(exponents for exponents, _ in achievers)
    return best, labels, len(rows), exponent_set


def deg_h_options(
    level: int,
    degstate: tuple[float, float, float, float],
    flags: tuple[bool, bool, bool],
    required: float | None = None,
) -> list[tuple[float, tuple[Obligation, ...]]]:
    """Possible exact degrees of h_level with their obligations.

    Unique achieving monomial: the degree is forced to the monomial degree
    and the polynomial cannot vanish identically.  Several achievers: the
    degree may DROP below the tropical maximum only through leading-
    coefficient cancellation, floored at 0 (a nonzero polynomial has
    nonnegative degree); deg = NEG_INF (h == 0 identically) is possible
    only via total cancellation.  All granted drops carry obligations.
    Options are ordered maximum first (fewest obligations first).
    """

    maximum, tied, term_count, exponent_set = tropical_h_max_full(
        level, degstate, flags
    )
    if maximum == NEG_INF:
        return [(NEG_INF, ())] if required in (None, NEG_INF) else []

    unique = len(tied) == 1
    drop_forbidden = (
        APPLY_RESIDUE_KILLS and (level, exponent_set) in FORBIDDEN_RISES
    )

    def drop(value: float) -> tuple[Obligation, ...] | None:
        if value == maximum:
            return ()
        if unique or term_count < 2 or drop_forbidden:
            return None
        if value == NEG_INF:
            return (Obligation(level, "identical_vanishing", 0, tied),)
        if value > maximum or value < 0:
            return None
        return (
            Obligation(level, "degree_tie_drop", int(maximum - value), tied),
        )

    if required is not None:
        if required != NEG_INF and (required > maximum or required < 0):
            return []
        obligations = drop(required)
        return [] if obligations is None else [(required, obligations)]

    options: list[tuple[float, tuple[Obligation, ...]]] = [(maximum, ())]
    if not (unique or term_count < 2 or drop_forbidden):
        for value in range(int(maximum) - 1, -1, -1):
            options.append(
                (
                    value,
                    (
                        Obligation(
                            level,
                            "degree_tie_drop",
                            int(maximum) - value,
                            tied,
                        ),
                    ),
                )
            )
        options.append(
            (NEG_INF, (Obligation(level, "identical_vanishing", 0, tied),))
        )
    return options


def descend_options_inf(
    level: int,
    deg_above: float,
    g_zero: bool,
    degstate: tuple[float, float, float, float],
    flags: tuple[bool, bool, bool],
    g_cap: int,
    deg_ehat: int,
    v: int,
) -> list[tuple[float, tuple[Obligation, ...]]]:
    """Solve t^v g_{l+1} = ehat^3 g_l + u^l h_l for deg g_l (max-plus).

    ``deg_above`` is deg g_{l+1}; NEG_INF encodes g_{l+1} == 0.  Zero flags
    are global statements exactly as in the finite descent.  Returned pairs
    are (deg g_l, obligations); deg g_l = NEG_INF encodes g_l == 0.
    """

    g_above_zero = deg_above == NEG_INF
    h_shift = DEG_U * level
    g_shift = 3 * deg_ehat
    maximum, max_labels, _, _ = tropical_h_max_full(level, degstate, flags)

    def h_side(w: float) -> str:
        """Exact leading form of u^l h_l at degree ``h_shift + w``."""

        if w == maximum:
            return f"({LC_U})^{level}*[" + " + ".join(max_labels) + "]"
        return f"({LC_U})^{level}*lc(h{level}@deg={int(w)})"

    g_side = f"lc(ehat)^3*lc(g{level})"

    if g_zero and g_above_zero:
        # 0 = 0 + u^l h_l: h_l must vanish identically.
        return [
            (NEG_INF, obligations)
            for _, obligations in deg_h_options(
                level, degstate, flags, required=NEG_INF
            )
        ]

    if g_zero:
        # t^v g_{l+1} = u^l h_l exactly: degrees match.
        required = v + deg_above - h_shift
        return [
            (NEG_INF, obligations)
            for _, obligations in deg_h_options(
                level, degstate, flags, required=required
            )
        ]

    if g_above_zero:
        # Polynomial identity ehat^3 g_l = -u^l h_l: degrees match exactly
        # and the two leading coefficients must cancel.
        results: list[tuple[float, tuple[Obligation, ...]]] = []
        for deg_g in range(0, g_cap + 1):
            required = g_shift + deg_g - h_shift
            for w, obligations in deg_h_options(
                level, degstate, flags, required=required
            ):
                identity = Obligation(
                    level, "exact_identity", 0, (g_side, h_side(w))
                )
                results.append((float(deg_g), obligations + (identity,)))
        return results

    target = v + deg_above
    options = deg_h_options(level, degstate, flags)
    results = []

    # Case (a): the g_l term dominates — g_shift + deg_g = target > h-term.
    deg_g = target - g_shift
    if 0 <= deg_g <= g_cap:
        for w, obligations in options:
            if w == NEG_INF or h_shift + w < target:
                results.append((float(deg_g), obligations))
                break  # options are ordered; first admissible w is maximal

    # Case (b): the h_l term dominates — h_shift + w = target > g-term.
    required = target - h_shift
    for _, obligations in deg_h_options(
        level, degstate, flags, required=required
    ):
        hi = min(g_cap, int(target - g_shift) - 1)
        for deg_g in range(0, hi + 1):
            results.append((float(deg_g), obligations))

    # Case (c): tie g_shift + deg_g = h_shift + w >= target, drop by
    # leading-coefficient cancellation of the two sides.
    for w, obligations in options:
        if w == NEG_INF:
            continue
        deg_g = h_shift + w - g_shift
        if deg_g < 0 or deg_g > g_cap:
            continue
        tie_value = h_shift + w
        if tie_value < target:
            continue
        depth = int(tie_value - target)
        extra = (
            (
                Obligation(
                    level,
                    "leading_cancellation",
                    depth,
                    (g_side, h_side(w)),
                ),
            )
            if depth > 0
            else ()
        )
        candidate = (float(deg_g), obligations + extra)
        if candidate not in results:
            results.append(candidate)

    return results


_INF_SIG_CACHE: dict = {}


def inf_place_profiles(
    a: int,
    branch: str,
    r_cap: int,
    depth: int,
    sigma_zero: bool,
    d2_zero: bool,
    g_zero: dict[int, bool],
    degstate: tuple[float, float, float, float],
    config: WindowConfig | None = None,
) -> tuple[InfProfile, ...]:
    """All consistent degree chains at infinity, Pareto-reduced.

    The chain runs from the terminal level down to level 1 and closes with
    the level-0 anchor t^v g_1 = h_0.  Levels below ``depth`` carry no
    outer zero flag; their g_l == 0 branches are folded into the chain
    enumeration (deg NEG_INF) and recoverable from the witness chain.
    Pareto keeps chains with maximal shared-level degrees and minimal
    obligation counts: the join constrains sum_p v_p(g_l) <= deg g_l, so
    larger degrees are more permissive.
    """

    if a > 10:
        raise ValueError("infinity layer requires the standard regime a <= 10")
    v = 30 - 3 * a
    deg_e = degstate[3]
    deg_ehat = int(deg_e - a)
    if deg_ehat < 0:
        return ()
    d1_zero = branch == "T2"
    flags = (sigma_zero, d2_zero, d1_zero)
    terminal = T1_TERMINAL if branch == "T1" else T2_TERMINAL
    full_levels = tuple(range(terminal, 0, -1))
    shared_levels = tuple(range(terminal, depth - 1, -1))
    _, g_caps, _ = resolve_caps(branch, r_cap, config, a)

    h_signature = tuple(
        tropical_h_max_full(level, degstate, flags)
        for level in range(0, terminal + 1)
    )
    cache_key = (
        "inf",
        APPLY_RESIDUE_KILLS,
        branch,
        v,
        deg_ehat,
        r_cap,
        depth,
        tuple(sorted(g_zero.items())),
        None if config is None else config.name,
        h_signature,
    )
    if cache_key in _INF_SIG_CACHE:
        return tuple(
            InfProfile(
                degs=degstate,
                chain=chain,
                obligations=obligations,
            )
            for chain, obligations in _INF_SIG_CACHE[cache_key]
        )

    memo: dict = {}

    def tails(level_index: int, deg_above: float):
        """Completed chains below full_levels[level_index - 1]."""

        if level_index == len(full_levels):
            # Level-0 anchor: t^v g_1 = h_0 (deg_above is deg g_1 here).
            required = (
                NEG_INF if deg_above == NEG_INF else v + deg_above
            )
            return [
                ((), obligations)
                for _, obligations in deg_h_options(
                    0, degstate, flags, required=required
                )
            ]
        key = (level_index, deg_above)
        if key in memo:
            return memo[key]
        level = full_levels[level_index]
        zero_branches = (
            (g_zero.get(level, False),)
            if level >= depth
            else (False, True)
        )
        collected = ParetoAccumulator(compact_at=20000)
        for zero_flag in zero_branches:
            for deg_l, obligations in descend_options_inf(
                level,
                deg_above,
                zero_flag,
                degstate,
                flags,
                g_caps[level],
                deg_ehat,
                v,
            ):
                for tail, tail_obligations in tails(level_index + 1, deg_l):
                    chain = ((level, deg_l),) + tail
                    total = obligations + tail_obligations
                    collected.add(
                        tuple(
                            0 if value == NEG_INF else -value
                            for _, value in chain
                        )
                        + (len(total),),
                        (chain, total),
                    )
        memo[key] = [payload for _, payload in collected.result()]
        return memo[key]

    accumulator = ParetoAccumulator()
    # Terminal identity ehat^3 g_T = -u^T h_T: degrees match exactly.
    if not g_zero.get(terminal, False):
        for r_t in range(0, g_caps[terminal] + 1):
            required = 3 * deg_ehat + r_t - DEG_U * terminal
            starts = deg_h_options(
                terminal, degstate, flags, required=required
            )
            if not starts:
                continue
            _, start_obligations = starts[0]
            for tail, tail_obligations in tails(1, float(r_t)):
                chain = ((terminal, float(r_t)),) + tail
                total = start_obligations + tail_obligations
                accumulator.add(
                    tuple(
                        0 if value == NEG_INF else -value
                        for level, value in chain
                        if level in shared_levels
                    )
                    + (len(total),),
                    (chain, total),
                )

    result = tuple(payload for _, payload in accumulator.result())
    _INF_SIG_CACHE[cache_key] = result
    return tuple(
        InfProfile(degs=degstate, chain=chain, obligations=obligations)
        for chain, obligations in result
    )


def _dfs_budget_witness(per_place, vectors, suffix, caps):
    """First finite-place selection whose budget sums fit under ``caps``."""

    dims = len(caps)
    found: list[tuple[PlaceProfile, ...] | None] = [None]

    def dfs(index, sums, chosen):
        if found[0] is not None:
            return
        if index == len(per_place):
            found[0] = chosen
            return
        for profile, vector in zip(per_place[index], vectors[index]):
            new_sums = tuple(sums[i] + vector[i] for i in range(dims))
            bound = suffix[index + 1]
            if all(new_sums[i] + bound[i] <= caps[i] for i in range(dims)):
                dfs(index + 1, new_sums, chosen + (profile,))
                if found[0] is not None:
                    return

    dfs(0, tuple(0 for _ in range(dims)), ())
    return found[0]


def join_places_inf(
    b_vector: tuple[int, int, int, int],
    branch: str,
    r_cap: int,
    depth: int,
    sigma_zero: bool,
    d2_zero: bool,
    g_zero: dict[int, bool],
    a: int,
    t_place_a: int | None = None,
    config: WindowConfig | None = None,
    t2_squeeze: bool = False,
) -> tuple[tuple[PlaceProfile, ...], InfProfile] | None:
    """Join the finite places with the place at infinity.

    Degrees are first-class unknowns sandwiched between the finite-place
    valuation sums (enforced by running the budget DFS with caps tightened
    to the chosen degrees) and the window caps; deg e >= a + sum(b_i).

    ``t2_squeeze`` applies the proven level-5 squeeze F^2 | G
    (T5_T2_COLUMN.md section 1, claim C24): writing e = t^a R F with
    gcd(F, tq) = 1 and g6 = Q G with Q the prescribed split-root part, the
    squeeze gives deg g6 >= sum_p v_p(g6) + 2 deg F with
    deg F = deg e - a - sum(b_i) — a tightened cap on the g6 budget
    dimension.  Its hypothesis (every m_i = v_p(g6) >= 1, so Q/q is a
    polynomial) holds automatically from the terminal law
    m_i = 6 + 2 s_i - 3 b_i whenever no b_i equals 2 (b_i in {0,1}: >= 3;
    b_i in {3,5}: odd hence >= 1); the constraint is skipped otherwise.
    """

    terminal = T1_TERMINAL if branch == "T1" else T2_TERMINAL
    levels = tuple(range(terminal, depth - 1, -1))
    aux, g_caps, _ = resolve_caps(branch, r_cap, config, a)
    per_place = [
        place_profiles(
            b, branch, r_cap, depth, sigma_zero, d2_zero, g_zero, config, a
        )
        for b in b_vector
    ]
    if t_place_a is not None:
        per_place.append(
            t_place_profiles(
                t_place_a,
                branch,
                r_cap,
                depth,
                sigma_zero,
                d2_zero,
                g_zero,
                config,
            )
        )
    if any(not options for options in per_place):
        return None

    dims = 3 + len(levels)
    vectors = [
        [profile.budget_vector(levels) for profile in options]
        for options in per_place
    ]
    minima = [
        tuple(min(vector[i] for vector in place) for i in range(dims))
        for place in vectors
    ]
    suffix = [tuple(0 for _ in range(dims))]
    for place_min in reversed(minima):
        last = suffix[0]
        suffix.insert(0, tuple(place_min[i] + last[i] for i in range(dims)))
    min_sums = suffix[0]

    # Infinity only tightens the join: no finite-place witness under the
    # original caps means the branch is dead before degrees are chosen.
    base_caps = (
        aux["d1"],
        aux["sigma"],
        aux["d2"],
        *(g_caps[level] for level in levels),
    )
    if _dfs_budget_witness(per_place, vectors, suffix, base_caps) is None:
        return None

    e_cap = (SUB2 if config is None else config).e_cap
    e_low = a + sum(b_vector)
    x_domain: tuple[float, ...] = (
        (NEG_INF,)
        if branch == "T2"
        else tuple(range(aux["d1"], -1, -1))
    )
    z_domain: tuple[float, ...] = (
        (NEG_INF,) if sigma_zero else tuple(range(aux["sigma"], -1, -1))
    )
    k_domain: tuple[float, ...] = (
        (NEG_INF,) if d2_zero else tuple(range(aux["d2"], -1, -1))
    )

    for deg_e in range(e_cap, e_low - 1, -1):
        for x_deg in x_domain:
            if x_deg != NEG_INF and x_deg < min_sums[0]:
                continue
            for z_deg in z_domain:
                if z_deg != NEG_INF and z_deg < min_sums[1]:
                    continue
                for k_deg in k_domain:
                    if k_deg != NEG_INF and k_deg < min_sums[2]:
                        continue
                    squeeze_slack = (
                        2 * (deg_e - e_low)
                        if (
                            t2_squeeze
                            and branch == "T2"
                            and all(b != 2 for b in b_vector)
                        )
                        else 0
                    )
                    degstate = (k_deg, x_deg, z_deg, float(deg_e))
                    for inf_profile in inf_place_profiles(
                        a,
                        branch,
                        r_cap,
                        depth,
                        sigma_zero,
                        d2_zero,
                        g_zero,
                        degstate,
                        config,
                    ):
                        chain = dict(inf_profile.chain)
                        caps = (
                            aux["d1"] if x_deg == NEG_INF else x_deg,
                            aux["sigma"] if z_deg == NEG_INF else z_deg,
                            aux["d2"] if k_deg == NEG_INF else k_deg,
                            *(
                                g_caps[level]
                                if chain[level] == NEG_INF
                                else chain[level]
                                - (squeeze_slack if level == 6 else 0)
                                for level in levels
                            ),
                        )
                        if any(
                            min_sums[i] > caps[i] for i in range(dims)
                        ):
                            continue
                        witness = _dfs_budget_witness(
                            per_place, vectors, suffix, caps
                        )
                        if witness is not None:
                            return witness, inf_profile
    return None


def join_places(
    b_vector: tuple[int, int, int, int],
    branch: str,
    r_cap: int,
    depth: int,
    sigma_zero: bool,
    d2_zero: bool,
    g_zero: dict[int, bool],
    t_place_a: int | None = None,
    config: WindowConfig | None = None,
    a: int | None = None,
) -> tuple[PlaceProfile, ...] | None:
    """DFS over the coupled places under the global degree budgets.

    The four q-root places always participate; when ``t_place_a`` is given,
    the place t = y+1 joins as a fifth place sharing the same budgets.
    """

    terminal = T1_TERMINAL if branch == "T1" else T2_TERMINAL
    levels = tuple(range(terminal, depth - 1, -1))
    aux, g_caps, _ = resolve_caps(branch, r_cap, config, a)
    per_place = [
        place_profiles(
            b, branch, r_cap, depth, sigma_zero, d2_zero, g_zero, config, a
        )
        for b in b_vector
    ]
    if t_place_a is not None:
        per_place.append(
            t_place_profiles(
                t_place_a,
                branch,
                r_cap,
                depth,
                sigma_zero,
                d2_zero,
                g_zero,
                config,
            )
        )
    if any(not options for options in per_place):
        return None

    caps = (
        aux["d1"],
        aux["sigma"],
        aux["d2"],
        *(g_caps[level] for level in levels),
    )
    dims = len(caps)
    vectors = [
        [profile.budget_vector(levels) for profile in options]
        for options in per_place
    ]
    minima = [
        tuple(min(vector[i] for vector in place) for i in range(dims))
        for place in vectors
    ]
    suffix = [tuple(0 for _ in range(dims))]
    for place_min in reversed(minima):
        last = suffix[0]
        suffix.insert(0, tuple(place_min[i] + last[i] for i in range(dims)))

    found: list[tuple[PlaceProfile, ...] | None] = [None]

    def dfs(index: int, sums: tuple[float, ...], chosen: tuple[PlaceProfile, ...]):
        if found[0] is not None:
            return
        if index == len(per_place):
            found[0] = chosen
            return
        for profile, vector in zip(per_place[index], vectors[index]):
            new_sums = tuple(sums[i] + vector[i] for i in range(dims))
            bound = suffix[index + 1]
            if all(new_sums[i] + bound[i] <= caps[i] for i in range(dims)):
                dfs(index + 1, new_sums, chosen + (profile,))
                if found[0] is not None:
                    return

    dfs(0, tuple(0 for _ in range(dims)), ())
    return found[0]


def encode(value: float) -> object:
    return "inf" if value == INF else int(value)


def encode_deg(value: float) -> object:
    return "-inf" if value == NEG_INF else int(value)


def inf_witness_record(profile: InfProfile) -> dict[str, object]:
    return {
        "place": "inf",
        "deg_d2": encode_deg(profile.degs[0]),
        "deg_d1": encode_deg(profile.degs[1]),
        "deg_sigma": encode_deg(profile.degs[2]),
        "deg_e": encode_deg(profile.degs[3]),
        "deg_g": {
            str(level): encode_deg(value) for level, value in profile.chain
        },
        "obligations": [
            {
                "level": obligation.level,
                "kind": obligation.kind,
                "depth": obligation.depth,
                "tied": list(obligation.tied),
            }
            for obligation in profile.obligations
        ],
    }


def analyze_branch(
    a: int,
    b_vector: tuple[int, int, int, int],
    branch: str,
    depth: int,
    include_t: bool = False,
    include_inf: bool = False,
    config: WindowConfig | None = None,
    t2_squeeze: bool = False,
) -> dict[str, object]:
    """Full zero-flag case split for one (stratum, branch) record."""

    r_cap = 10 + 3 * a
    terminal = T1_TERMINAL if branch == "T1" else T2_TERMINAL
    below_terminal = [level for level in range(terminal - 1, depth - 1, -1)]

    survivors = []
    for sigma_zero in ((False, True) if branch == "T1" else (False,)):
        for d2_zero in (False, True):
            for zero_mask in range(2 ** len(below_terminal)):
                g_zero = {
                    level: bool(zero_mask >> i & 1)
                    for i, level in enumerate(below_terminal)
                }
                # The terminal g is a unit multiple of h_terminal / ehat^3 and
                # cannot vanish under the branch hypothesis (d1 != 0 resp.
                # sigma != 0).
                g_zero[terminal] = False
                inf_profile = None
                if include_inf:
                    joined = join_places_inf(
                        b_vector,
                        branch,
                        r_cap,
                        depth,
                        sigma_zero,
                        d2_zero,
                        g_zero,
                        a,
                        t_place_a=a if include_t else None,
                        config=config,
                        t2_squeeze=t2_squeeze,
                    )
                    witness = None if joined is None else joined[0]
                    inf_profile = None if joined is None else joined[1]
                else:
                    witness = join_places(
                        b_vector,
                        branch,
                        r_cap,
                        depth,
                        sigma_zero,
                        d2_zero,
                        g_zero,
                        t_place_a=a if include_t else None,
                        config=config,
                        a=a,
                    )
                if witness is not None:
                    survivors.append(
                        {
                            "sigma_zero": sigma_zero,
                            "d2_zero": d2_zero,
                            "g_zero_levels": sorted(
                                level for level, flag in g_zero.items() if flag
                            ),
                            "witness": [
                                {
                                    "place": profile.place,
                                    "b": profile.b,
                                    "v_d1": encode(profile.x),
                                    "v_sigma": encode(profile.z),
                                    "v_d2": encode(profile.k),
                                    "v_g": {
                                        str(level): encode(value)
                                        for level, value in profile.r
                                    },
                                    "obligations": [
                                        {
                                            "level": obligation.level,
                                            "kind": obligation.kind,
                                            "depth": obligation.depth,
                                            "tied": list(obligation.tied),
                                        }
                                        for obligation in profile.obligations
                                    ],
                                }
                                for profile in witness
                            ]
                            + (
                                [inf_witness_record(inf_profile)]
                                if inf_profile is not None
                                else []
                            ),
                            "obligation_count": sum(
                                len(profile.obligations) for profile in witness
                            )
                            + (
                                len(inf_profile.obligations)
                                if inf_profile is not None
                                else 0
                            ),
                        }
                    )
    return {
        "depth": depth,
        "status": (
            "engine_killed_pending_audit" if not survivors else "survives"
        ),
        "survivor_cases": survivors,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--depth",
        type=int,
        default=4,
        choices=(4, 5, 6, 7),
        help="lowest cascade level processed (7/6 = terminal only)",
    )
    parser.add_argument("--json-out", type=Path, default=JSON_OUT)
    parser.add_argument(
        "--with-t",
        action="store_true",
        help="couple the place t=y+1 as a fifth place (standard regime only)",
    )
    parser.add_argument(
        "--with-inf",
        action="store_true",
        help=(
            "couple infinity as a sixth place (max-plus degree layer, "
            "standard regime only)"
        ),
    )
    parser.add_argument(
        "--window",
        choices=("sub2", "sub1"),
        default="sub2",
        help="window configuration and ledger (sub1: standard regime a<=10)",
    )
    parser.add_argument(
        "--residue-kills",
        action="store_true",
        help="apply the proven forbidden-rise lemmas (RESIDUE_LEMMAS.md)",
    )
    parser.add_argument(
        "--t2-squeeze",
        action="store_true",
        help=(
            "apply the proven T2 level-5 squeeze F^2|G (T5_T2_COLUMN.md, "
            "C24) in the infinity join; auto-skipped when some b_i = 2"
        ),
    )
    parser.add_argument(
        "--max-rss-gb",
        type=float,
        default=6.0,
        help="abort cleanly (exit 3) if process RSS exceeds this many GiB",
    )
    args = parser.parse_args()

    try:
        import psutil

        process = psutil.Process()

        def rss_gb() -> float:
            return process.memory_info().rss / 2**30

    except ImportError:  # guard becomes a no-op without psutil

        def rss_gb() -> float:
            return 0.0

    global APPLY_RESIDUE_KILLS
    if args.residue_kills:
        APPLY_RESIDUE_KILLS = True
    config = CONFIGS[args.window]
    ledger_path = (
        LEDGER_PATH
        if args.window == "sub2"
        else ROOT / "split_place_ledger_sub1.json"
    )
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    results = []
    killed = surviving = 0
    skipped_alternate = 0
    resource_aborted = False
    # Process cheap strata first (high a => small b-ranges) and clear the
    # profile caches between a-chunks so memory stays bounded; checkpoint
    # the JSON after every chunk so a crash loses at most one chunk.
    rows = sorted(ledger["strata"], key=lambda row: -row["a_t"])
    current_a = None
    for row in rows:
        if args.window == "sub1" and row["a_t"] > 10:
            skipped_alternate += 1
            continue
        if resource_aborted:
            for branch in row["open_branches"]:
                results.append(
                    {
                        "a_t": row["a_t"],
                        "b": row["b"],
                        "branch": branch,
                        "depth": args.depth,
                        "status": "skipped_resource_limit",
                        "survivor_case_count": 0,
                        "survivor_cases": [],
                    }
                )
            continue
        if row["a_t"] != current_a:
            if current_a is not None:
                clear_profile_caches()
                write_payload(args, results, killed, surviving,
                              skipped_alternate, partial=True)
            current_a = row["a_t"]
        for branch in row["open_branches"]:
            outcome = analyze_branch(
                row["a_t"],
                tuple(row["b"]),
                branch,
                args.depth,
                include_t=args.with_t,
                include_inf=args.with_inf,
                config=config,
                t2_squeeze=args.t2_squeeze,
            )
            results.append(
                {
                    "a_t": row["a_t"],
                    "b": row["b"],
                    "branch": branch,
                    "depth": outcome["depth"],
                    "status": outcome["status"],
                    "survivor_case_count": len(outcome["survivor_cases"]),
                    "survivor_cases": outcome["survivor_cases"],
                }
            )
            if outcome["status"] == "engine_killed_pending_audit":
                killed += 1
            else:
                surviving += 1
            if rss_gb() > args.max_rss_gb:
                print(
                    f"RSS {rss_gb():.1f} GiB > {args.max_rss_gb} GiB after "
                    f"a={row['a_t']} b={row['b']} {branch}; aborting cleanly"
                )
                resource_aborted = True

    write_payload(args, results, killed, surviving, skipped_alternate,
                  partial=False)
    if resource_aborted:
        raise SystemExit(3)


def write_payload(args, results, killed, surviving, skipped_alternate,
                  partial: bool) -> None:
    skipped = sum(
        1 for r in results if r["status"] == "skipped_resource_limit"
    )
    payload = {
        "schema": 1,
        "description": (
            "Cascade-engine output: exact valuation descent "
            f"to level {args.depth} on the open split-place branches at "
            + ("the four q-places plus t" if args.with_t else "the four q-places")
            + (" plus infinity (max-plus degree layer)" if args.with_inf else "")
        ),
        "depth": args.depth,
        "places": ("q+t" if args.with_t else "q")
        + ("+inf" if args.with_inf else ""),
        "window": args.window,
        "residue_kills": APPLY_RESIDUE_KILLS,
        "t2_squeeze": bool(getattr(args, "t2_squeeze", False)),
        "partial_checkpoint": partial,
        "summary": {
            "open_branches_processed": killed + surviving,
            "engine_killed_pending_audit": killed,
            "surviving_branches": surviving,
            "alternate_regime_strata_skipped": skipped_alternate,
            "skipped_resource_limit": skipped,
        },
        "branches": results,
    }
    args.json_out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not partial:
        print(json.dumps(payload["summary"], indent=2))
        print(f"wrote {args.json_out.name}")


if __name__ == "__main__":
    main()
