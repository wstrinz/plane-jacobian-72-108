"""
t5_pivot_flint.py — pivoted t-adic lift for T5 branch (i) (a=0, f31 subcase-2)
reimplemented in python-flint with truncated-series arithmetic.

Replaces the failed t5_pivot.sing (repeated full re-substitution inside
qring S/(t^41) timed out at pivot slice 5). Here every window polynomial is a
list of nmod_mpoly coefficients (index = t-power) in the 18 FREE window
unknowns (a1..a4, b1..b6, e1..e8) over F_32003; products are truncated
convolutions, evaluation is Horner-structured so every big multiplication has
one SMALL factor.

Base point (from t5_pivot_gen.py, seed 20260721, kept for comparability):
    (a0,b0,c0,e0) = (20066, 6066, 31791, 22015)
with grad_d0 h0 and grad_dm1 h0 nonzero there (inverses 17370, 27284).

KEY STRUCTURAL FACT (discovered on the first full-truncation run): the lift is
weight-graded. Give a_i, b_i, e_i weight i (= the t-power they sit at); then
t-slice j of any window object is weighted-homogeneous of weight j, and the
computed pivots are DENSE in their weight spaces: term count of pivot j equals
    P(j) := #{monomials of weighted degree j in the 18 frees}
exactly (3, 9, 22, 51, 107, 217, 415, 771, 1379, 2407 for j=1..10 — matching
the Singular run's size(piv) for j<=5). Consequences, with
P(29) = 4,812,342 / P(30) = 6,535,586 / P(40) = 103,447,659 /
P(59) = 6,838,224,687:
  * block-1 equations Q_11..Q_29 top out at P(29) < 5M  -> computable;
  * g1 = H[30..40] needs entries of sizes P(30)..P(40)  -> OVER the 5M cap
    (the first run aborted at exactly 6,504,802 ~ P(30) terms in a weight-30
    power-table entry);
  * block-2 equations live at weight up to 59            -> ~10^9-10^10 terms.
So the pipeline computes block 1 completely (mod t^30 arithmetic — lossless
for slices 0..29) and STOPS, reporting the swell, per the task's safety rule.

Outputs: t5_pivot_eqs.txt (all 19 block-1 eqs, Singular syntax, chunked
lines), t5_pivot_eqs_small.txt (Q_11..Q_16 only, for fast GB staging),
t5_pivot_gb.sing (staged slimgb driver).
"""
import re
import sys
import time
from collections import defaultdict

from flint import nmod_mpoly_ctx, Ordering

P = 32003
TLEN = 30                 # block-1 truncation: slices 0..29 (t^30 ~ 0)
CAP = 5_000_000
BASE = (20066, 6066, 31791, 22015)          # a0 (d2), b0 (d1), c0 (d0), e0 (dm1)
EXPECT_INV_G0, EXPECT_INV_GM = 17370, 27284  # from t5_pivot_gen.py, sanity only

FREES = ([f'a{i}' for i in range(1, 5)] + [f'b{i}' for i in range(1, 7)]
         + [f'e{i}' for i in range(1, 9)])
WEIGHTS = [1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 8]

ctx = nmod_mpoly_ctx.get(tuple(FREES), P, Ordering.degrevlex)
GENS = ctx.gens()
ONE = ctx.constant(1)
Z = ctx.constant(0)


def cst(v):
    return ctx.constant(int(v) % P)


class TermSwell(Exception):
    pass


def chk(poly, where):
    n = len(poly)
    if n > CAP:
        raise TermSwell(f"{where}: {n} terms > cap {CAP}")
    return n


# dense weighted dimensions P(j)
PDIM = [0] * 80
PDIM[0] = 1
for w in WEIGHTS:
    for j in range(w, 80):
        PDIM[j] += PDIM[j - w]

# ---------------------------------------------------------------------------
# parse h_0..h_7 from f31_graded.txt into [(coeff mod p, (e_d2,e_d1,e_d0,e_dm1))]
# ---------------------------------------------------------------------------
SRC = open('f31_graded.txt').read()
H_SRC = {int(m.group(1)): m.group(2)
         for m in re.finditer(r'h_(\d) \(weight \d+, dm1-power \d+\) = (.+)', SRC)}
assert sorted(H_SRC) == list(range(8))

VAR_IDX = {'d2': 0, 'd1': 1, 'd0': 2, 'dm1': 3}


def parse_h(src):
    """Parse '8957952*d0**5 - 796...*d0**4*d2**2 + ...' exactly."""
    terms = []
    for piece in re.finditer(r'([+-]?)\s*(\d+)((?:\*[a-z]\w*(?:\*\*\d+)?)+)', src):
        sign, coef, mono = piece.group(1), int(piece.group(2)), piece.group(3)
        if sign == '-':
            coef = -coef
        exps = [0, 0, 0, 0]
        for vm in re.finditer(r'\*([a-z]\w*?)(?:\*\*(\d+))?(?=\*[a-z]|$)', mono):
            exps[VAR_IDX[vm.group(1)]] += int(vm.group(2) or 1)
        terms.append((coef, tuple(exps)))
    return terms


H_TERMS = {f: parse_h(H_SRC[f]) for f in range(8)}
assert len(H_TERMS[0]) == 28 and len(H_TERMS[1]) == 22 and len(H_TERMS[7]) == 1

# exact-integer sanity of the parser against a plain eval() of the source.
# Safe here: input is the trusted repo-local f31_graded.txt (pure arithmetic
# expressions), evaluated with builtins stripped; used only as a cross-check.
_ns = {'d2': 3, 'd1': 5, 'd0': 7, 'dm1': 11}
for f in range(8):
    want = eval(H_SRC[f], {'__builtins__': {}}, dict(_ns))
    got = sum(c * 3**e[0] * 5**e[1] * 7**e[2] * 11**e[3] for c, e in H_TERMS[f])
    assert want == got, f"parser mismatch at h_{f}"

# ---------------------------------------------------------------------------
# base verification + gradients (plain modular integer arithmetic)
# ---------------------------------------------------------------------------


def h_at(terms, pt):
    s = 0
    for c, e in terms:
        s += c * pow(pt[0], e[0], P) * pow(pt[1], e[1], P) \
               * pow(pt[2], e[2], P) * pow(pt[3], e[3], P)
    return s % P


def grad_at(terms, pt, k):
    s = 0
    for c, e in terms:
        if e[k] == 0:
            continue
        ee = list(e)
        ee[k] -= 1
        s += c * e[k] * pow(pt[0], ee[0], P) * pow(pt[1], ee[1], P) \
                      * pow(pt[2], ee[2], P) * pow(pt[3], ee[3], P)
    return s % P


assert h_at(H_TERMS[0], BASE) == 0, "base is not on h0=0"
G0 = grad_at(H_TERMS[0], BASE, 2)    # d(h0)/d(d0)
GM = grad_at(H_TERMS[0], BASE, 3)    # d(h0)/d(dm1)
assert G0 and GM
INV_G0 = pow(G0, P - 2, P)
INV_GM = pow(GM, P - 2, P)
assert INV_G0 == EXPECT_INV_G0 and INV_GM == EXPECT_INV_GM, \
    f"gradient mismatch vs t5_pivot_gen.py: {INV_G0},{INV_GM}"
NEG_INV_G0 = (P - INV_G0) % P
NEG_INV_GM = (P - INV_GM) % P
print(f"base {BASE} verified on h0=0; grad_d0={G0}, grad_dm1={GM}")

# ---------------------------------------------------------------------------
# truncated-series toolkit (series = python list of nmod_mpoly, index=t-power)
# ---------------------------------------------------------------------------


def ser_mul(A, B, tlen, tag="ser_mul"):
    C = [Z] * tlen
    na = min(len(A), tlen)
    for i in range(na):
        a = A[i]
        if not a:
            continue
        nb = min(len(B), tlen - i)
        for j in range(nb):
            b = B[j]
            if not b:
                continue
            C[i + j] += a * b
    for q in C:
        chk(q, tag)
    return C


class PowTab:
    """memoized powers of a window series at fixed truncation length."""

    def __init__(self, series, tlen, tag):
        self.tlen = tlen
        self.tag = tag
        one = [ONE] + [Z] * (tlen - 1)
        self.tab = [one, [series[i] if i < len(series) else Z
                          for i in range(tlen)]]

    def get(self, k):
        while len(self.tab) <= k:
            self.tab.append(ser_mul(self.tab[-1], self.tab[1], self.tlen,
                                    f"{self.tag}^{len(self.tab)}"))
        return self.tab[k]


def eval_h_horner(terms, D2s, D1s, D0s, Es, tlen, verbose=False):
    """Evaluate h(D2,D1,D0,E) as a truncated series, Horner in d0.

    Every large multiplication has one small factor: the E^d tables are the
    only big-by-medium products; the Horner steps multiply the (big)
    accumulator by the raw D0 window (<= 1+sum P(1..8) = 1856 terms).
    """
    t2 = PowTab(D2s, tlen, "D2")
    t1 = PowTab(D1s, tlen, "D1")
    te = PowTab(Es, tlen, "E")
    by_c = defaultdict(list)
    cmax = 0
    for coef, (x2, x1, x0, xm) in terms:
        by_c[x0].append((coef, x2, x1, xm))
        cmax = max(cmax, x0)
    R = None
    for c in range(cmax, -1, -1):
        t0 = time.perf_counter()
        A = [Z] * tlen
        for coef, x2, x1, xm in by_c.get(c, []):
            cc = cst(coef)
            cur = None
            for tab, e in ((t2, x2), (t1, x1)):
                if e:
                    s = tab.get(e)
                    cur = s if cur is None else ser_mul(cur, s, tlen, "small")
            if xm:
                s = te.get(xm)
                cur = s if cur is None else ser_mul(cur, s, tlen, "xE")
            if cur is None:
                A[0] += cc
            else:
                for i in range(tlen):
                    q = cur[i]
                    if q:
                        A[i] += q * cc
        if R is None:
            R = A
        else:
            R = ser_mul(R, D0s, tlen, "horner*D0")
            for i in range(tlen):
                if A[i]:
                    R[i] += A[i]
        if verbose:
            print(f"    horner c={c}: {time.perf_counter()-t0:.1f}s, "
                  f"acc sizes max {max(len(q) for q in R)}")
            sys.stdout.flush()
    for q in R:
        chk(q, "eval_h slice")
    return R


# ---------------------------------------------------------------------------
# windows
# ---------------------------------------------------------------------------
a0, b0, c0, e0 = BASE
D2 = [cst(a0)] + list(GENS[0:4])            # deg <= 4
D1 = [cst(b0)] + list(GENS[4:10])           # deg <= 6
D0 = [cst(c0)] + [Z] * 8                    # c1..c8 pivoted in
E = [cst(e0)] + list(GENS[10:18]) + [Z, Z]  # e9,e10 pivoted in

# ---------------------------------------------------------------------------
# phase 1: pivot lift, slices 1..10
# ---------------------------------------------------------------------------
SING_REF = [3, 9, 22, 51, 107]     # Singular size(piv) at slices 1..5
pivot_counts = []
print("\n--- pivot lift ---")
t_lift0 = time.perf_counter()
for j in range(1, 11):
    t0 = time.perf_counter()
    tl = j + 1
    S = eval_h_horner(H_TERMS[0], D2, D1, D0, E, tl)[j]
    if j <= 8:
        piv = S * NEG_INV_G0
        D0[j] = piv
        tgt = 'c%d' % j
    else:
        piv = S * NEG_INV_GM
        E[j] = piv
        tgt = 'e%d' % j
    # recheck: slice j must vanish after insertion
    S2 = eval_h_horner(H_TERMS[0], D2, D1, D0, E, tl)[j]
    assert not S2, f"slice {j} nonzero after pivot insertion — BUG"
    n = chk(piv, f"pivot {tgt}")
    pivot_counts.append(n)
    note = f", P({j})={PDIM[j]}{'=' if n == PDIM[j] else '!='}dense"
    ref = (f" (Singular ref {SING_REF[j-1]}"
           f"{' MATCH' if n == SING_REF[j-1] else ' MISMATCH!'})") \
        if j <= len(SING_REF) else ""
    print(f"slice {j:2d} -> pivot {tgt}: {n} terms{ref}{note}, "
          f"{time.perf_counter()-t0:.2f}s")
    sys.stdout.flush()
print(f"pivot lift total: {time.perf_counter()-t_lift0:.2f}s")

# ---------------------------------------------------------------------------
# phase 2: full H = h0(windows) mod t^30; block-1 consistency equations
# ---------------------------------------------------------------------------
print("\n--- block 1: h0(windows) mod t^30 (Horner) ---")
t0 = time.perf_counter()
H = eval_h_horner(H_TERMS[0], D2, D1, D0, E, TLEN, verbose=True)
print(f"h0(windows) mod t^30 done in {time.perf_counter()-t0:.1f}s")
for j in range(11):
    assert not H[j], f"H[{j}] nonzero — pivot lift inconsistent"
EQS = H[11:30]
print("block-1 eqs Q_11..Q_29 term counts (vs dense P(j)):")
for k, q in enumerate(EQS):
    j = 11 + k
    print(f"  Q_{j}: {len(q)} terms (P({j})={PDIM[j]}), "
          f"deg {q.total_degree()}")
sys.stdout.flush()

# ---------------------------------------------------------------------------
# phase 3: blocks 2..7 — infeasibility statement (term-swell stop)
# ---------------------------------------------------------------------------
print("\n--- blocks 2+: STOP (term swell, per safety rule) ---")
print("g1 = H[30..40] would need dense weighted components of sizes:")
print("  " + ", ".join(f"P({w})={PDIM[w]}" for w in range(30, 41)))
print(f"first entry P(30)={PDIM[30]} already exceeds the {CAP} cap")
print("(first full-truncation run aborted at 6,504,802 terms in a weight-30")
print(" power-table entry, matching P(30) up to near-density).")
print(f"block-2 equations would live at weight up to 59: P(59)={PDIM[59]}")

# ---------------------------------------------------------------------------
# phase 4: export
# ---------------------------------------------------------------------------
CHUNK = 20000  # terms per Singular source line


def write_eqs(fname, idx_list):
    t0 = time.perf_counter()
    with open(fname, 'w') as fh:
        fh.write(f"// t5_pivot_flint.py block-1 equations, base {BASE}, F_{P}\n")
        fh.write(f"// EQ[k] = t-slice (10+k) of h0(windows), k=1..{len(idx_list)}\n")
        for out_k, eq_i in enumerate(idx_list, 1):
            q = EQS[eq_i]
            s = str(q)
            terms = s.split(' + ')
            fh.write(f"EQ[{out_k}]=0;\n")
            for lo in range(0, len(terms), CHUNK):
                chunk = '+'.join(terms[lo:lo + CHUNK])
                fh.write(f"EQ[{out_k}]=EQ[{out_k}]+({chunk});\n")
    print(f"wrote {fname} ({time.perf_counter()-t0:.1f}s)")
    sys.stdout.flush()


write_eqs('t5_pivot_eqs_small.txt', list(range(0, 6)))     # Q_11..Q_16
write_eqs('t5_pivot_eqs.txt', list(range(len(EQS))))       # Q_11..Q_29

gb = f"""// generated by t5_pivot_flint.py — do not edit by hand
// block-1 pivoted-lift consistency equations, base {BASE}, 18 frees.
// Q_11..Q_29 are dense weighted-homogeneous of weight j (sizes P(11)..P(29),
// max 4,812,342 terms); blocks 2+ were infeasible (term swell > 5M).
// Staged: (A) degBound-capped std on Q_11..Q_16, (B) slimgb Q_11..Q_16,
// (C) slimgb on all 19.  1 in the ideal at ANY stage = per-base certificate
// (a subideal containing 1 certifies the full system).
option(prot);
ring R={P},({','.join(FREES)}),dp;
ideal EQ;
int j;
< "t5_pivot_eqs_small.txt";
"loaded small system, eqs:",ncols(EQ);

"=== stage A: std(Q_11..Q_16) with degBound=16 ===";
degBound=16;
ideal SGA=std(EQ);
degBound=0;
"stage A GB size:",size(SGA);
if (SGA[1]==1) {{ "*** 1 IN IDEAL (stage A) — PER-BASE CERTIFICATE ***"; quit; }}
"stage A (deg-capped, no unit found; not conclusive for consistency)";

"=== stage B: slimgb(Q_11..Q_16) ===";
ideal SGB=slimgb(EQ);
"stage B GB size:",size(SGB);
if (SGB[1]==1) {{ "*** 1 IN IDEAL (stage B) — PER-BASE CERTIFICATE ***"; quit; }}
"stage B dim:",dim(SGB);

"=== stage C: loading all 19 block-1 eqs (large file) ===";
ideal EQ2=EQ;   // keep
kill EQ;
ideal EQ;
< "t5_pivot_eqs.txt";
"loaded full block-1 system, eqs:",ncols(EQ);
ideal SGC=slimgb(EQ);
"stage C GB size:",size(SGC);
if (SGC[1]==1) {{ "*** 1 IN IDEAL (stage C) — PER-BASE CERTIFICATE ***"; }}
else {{ "stage C dim:",dim(SGC); "no unit — block-1 ideal PROPER for this base"; }}
quit;
"""
open('t5_pivot_gb.sing', 'w').write(gb)
print("wrote t5_pivot_gb.sing")
print("pivot term counts:", pivot_counts)
print("DONE")
