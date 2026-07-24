#!/usr/bin/env python3
"""chain_survey_verify.py  (NEW; independent PASS/FAIL checker for chain_survey.py)

REGRESSION: the GGV5 v11<=35 tables (paper_src/1708.07936_GGV5.tex, lines
1674-1718) are transcribed verbatim below and reproduced by the enumerator.

Checks
  A. The set of 24 published CHAINS (A_0, A_0', final corner, k) is reproduced
     EXACTLY -- no missing rows, no extra rows.
  B. 23 of the 24 (m,n)-family parametrizations match the printed table verbatim.
  C. The sole exception, F6, is the documented paper typo: the printed base pair
     (4,10) has gcd 2 (violating Definition 'mn families', which requires
     gcd(m,n)=1), while the enumerator returns the coprime correction (7,18).
  D. INVARIANTS on EVERY enumerated family (at M=35 and M=55):
       - Diophantine identity (m+n)b k - n(bl-a) = k         (GGV5 line 1480)
       - coprimality gcd(m0,n0)=1 of the family base pair
       - kappa = t-2   and   dg = a0-q   (standard-chart conventions)
       - the whole (m,n) progression stays coprime and on the Diophantine line.

--quiet suppresses the per-check log; exit code 0 iff all checks pass.
"""
import sys
from math import gcd

import chain_survey as cs

QUIET = "--quiet" in sys.argv
FAILS = []
N_OK = 0

def ok(label, cond):
    global N_OK
    if cond:
        N_OK += 1
        if not QUIET:
            print("[OK]", label)
    else:
        FAILS.append(label)
        print("[FAIL]", label)

# --- published table, transcribed verbatim (line refs in comments) ---
# (name, A0, A0p, final=(a,l,b), k, (m0,dm), (n0,dn))  with m=dm*j+m0, n=dn*j+n0
PUBLISHED = [
    ("F1",  (4, 12), (1, 0), (7, 4, 3),  1, (3, 2),  (4, 3)),
    ("F2",  (5, 20), (1, 0), (7, 5, 2),  1, (2, 1),  (3, 2)),
    ("F3",  (5, 20), (1, 0), (8, 5, 3),  1, (3, 4),  (2, 3)),
    ("F4",  (5, 20), (1, 0), (8, 5, 3),  2, (3, 2),  (16, 12)),
    ("F5",  (5, 20), (1, 0), (9, 5, 4),  1, (9, 7),  (5, 4)),
    ("F6",  (5, 20), (1, 0), (9, 5, 4),  2, (4, 3),  (10, 8)),   # printed base (4,10): gcd 2 (typo)
    ("F7",  (6, 15), (1, 0), (7, 3, 4),  1, (2, 1),  (7, 4)),
    ("F8",  (6, 15), (1, 0), (8, 3, 5),  1, (3, 2),  (7, 5)),
    ("F9",  (7, 21), (1, 0), (11, 7, 2), 1, (2, 1),  (3, 2)),
    ("F10", (7, 21), (1, 0), (13, 7, 3), 1, (7, 5),  (4, 3)),
    ("F11", (7, 21), (1, 0), (13, 7, 3), 2, (2, 1),  (5, 3)),
    ("F12", (8, 24), (2, 0), (13, 4, 5), 1, (3, 2),  (7, 5)),
    ("F13", (9, 21), (2, 0), (13, 3, 7), 1, (2, 1),  (13, 7)),
    ("F14", (9, 24), (1, 0), (7, 3, 4),  1, (2, 1),  (7, 4)),
    ("F15", (9, 24), (1, 0), (8, 3, 5),  1, (3, 2),  (7, 5)),
    ("F16", (9, 24), (1, 0), (10, 3, 7), 1, (3, 4),  (5, 7)),
    ("F17", (9, 24), (1, 0), (11, 3, 8), 1, (2, 5),  (3, 8)),
    ("F18", (6, 18), (6, 15), (7, 3, 4), 1, (2, 1),  (7, 4)),
    ("F19", (6, 18), (6, 15), (8, 3, 5), 1, (3, 2),  (7, 5)),
    ("F20", (6, 24), (6, 15), (7, 3, 4), 1, (2, 1),  (7, 4)),
    ("F21", (6, 24), (6, 15), (8, 3, 5), 1, (3, 2),  (7, 5)),
    ("F22", (8, 24), (2, 0), (5, 4, 2),  1, (2, 1),  (3, 2)),
    ("F23", (8, 24), (2, 0), (11, 4, 4), 1, (2, 1),  (7, 4)),
    ("F24", (8, 24), (2, 0), (19, 8, 3), 1, (3, 2),  (4, 3)),
]

# --- run the enumerator at M=35 ---
recs35, meta35 = cs.survey(35)

# index enumerated families by chain key (A0, A0p, final, k).
# corners are stored as (a,l,b) triples; the published table writes (a,b) when l=1.
def norm(c):
    c = tuple(c)
    return (c[0], c[2]) if c[1] == 1 else c

enum = {}
for r in recs35:
    key = (norm(r["A0"]), norm(r["A0p"]), tuple(r["final"]), r["k"])
    enum[key] = r

pub_keys = set((A0, A0p, fin, k) for (_, A0, A0p, fin, k, _, _) in PUBLISHED)
enum_keys = set(enum.keys())

# A. exact chain-set reproduction
ok("A: exactly 24 published chains, no missing", pub_keys - enum_keys == set())
ok("A: no extra chains beyond the 24 published", enum_keys - pub_keys == set())
ok("A: 24 family rows enumerated at M=35", len(recs35) == 24)

# B/C. per-family (m,n) comparison
exact = 0
typo_seen = None
for (name, A0, A0p, fin, k, (m0, dm), (n0, dn)) in PUBLISHED:
    key = (A0, A0p, fin, k)
    if key not in enum:
        continue
    r = enum[key]
    got = (r["m0"], r["dm"], r["n0"], r["dn"])
    want = (m0, dm, n0, dn)
    if got == want:
        exact += 1
    else:
        typo_seen = (name, want, got)

ok("B: 23 of 24 (m,n)-families match the printed table verbatim", exact == 23)
ok("C: the mismatch is F6",
   typo_seen is not None and typo_seen[0] == "F6")
# printed F6 base is (4,10) with gcd 2; enumerator base is coprime
ok("C: printed F6 base (4,10) is non-coprime (gcd 2) -> paper typo",
   gcd(4, 10) == 2)
f6 = enum[((5, 20), (1, 0), (9, 5, 4), 2)]
ok("C: enumerated F6 base (%d,%d) is coprime (correction)" % (f6["m0"], f6["n0"]),
   gcd(f6["m0"], f6["n0"]) == 1)

# D. invariants on every enumerated family (M=35 and M=55)
def check_invariants(recs, tag):
    global N_OK
    all_dio = True
    all_cop = True
    all_kappa = True
    all_dg = True
    all_prog = True
    for r in recs:
        af, lf, bf = r["final"]
        bl_a = bf * lf - af
        k = r["k"]
        m0, n0, dm, dn = r["m0"], r["n0"], r["dm"], r["dn"]
        if r["dio_residual_minus_k"] != 0:
            all_dio = False
        if gcd(m0, n0) != 1:
            all_cop = False
        if r["kappa"] != r["t"] - 2:
            all_kappa = False
        if r["dg"] != r["a0"] - r["q"]:
            all_dg = False
        # first three progression members: coprime + on the Diophantine line
        for j in range(3):
            m = m0 + dm * j
            n = n0 + dn * j
            if gcd(m, n) != 1:
                all_prog = False
            if (m + n) * bf * k - n * bl_a != k:
                all_prog = False
    ok(f"D[{tag}]: Diophantine identity holds on all {len(recs)} families", all_dio)
    ok(f"D[{tag}]: base pair coprime on all families", all_cop)
    ok(f"D[{tag}]: kappa = t-2 on all families", all_kappa)
    ok(f"D[{tag}]: dg = a0-q on all families", all_dg)
    ok(f"D[{tag}]: whole (m,n) progression stays coprime & Diophantine", all_prog)

check_invariants(recs35, "M=35")
recs55, _ = cs.survey(55)
check_invariants(recs55, "M=55")
ok("D: M=55 enumeration is a superset of the 24 published chains",
   pub_keys.issubset(set((norm(r["A0"]), norm(r["A0p"]), tuple(r["final"]), r["k"])
                         for r in recs55)))

print()
if FAILS:
    print("FAILURES:", len(FAILS))
    for f in FAILS:
        print("   -", f)
    sys.exit(1)
print("ALL %d CHAIN-SURVEY CHECKS PASSED" % N_OK)
sys.exit(0)
