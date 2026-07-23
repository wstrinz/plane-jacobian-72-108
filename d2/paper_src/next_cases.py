"""Recompute the GGV candidate (deg P, deg Q) pairs from the GGV5 recipe.

Recipe (GGV5 = arXiv:1708.07936):
 - deg P = m * v11(A0), deg Q = n * v11(A0), v11(a,b) = a+b   [GGV5 line ~250]
 - For a final corner (a/l, b): k in I(A) iff 1<=k<l-a/b and gcd(b,(b*l-a)/gcd(k,b*l-a))=1
 - MN_k(A) = {(m,n): m,n>1, gcd(m,n)=1, (m+n)*b*k - n*(b*l-a) = k}    ["ecuacion diofantica"]
   (plus the swapped equation with m<->n)
"""
from math import gcd

# --- Part 1: verify the Diophantine recipe regenerates the family (m,n) lists ---
def mn_solutions(a, l, b, maxmn=40):
    """All (k,(m,n)) with the diophantine equation, m,n>1 coprime, small."""
    out = []
    # k < l - a/b  <=>  k*b < l*b - a
    k = 1
    while k * b < l * b - a:
        e = gcd(k, b * l - a)
        if gcd(b, (b * l - a) // e) == 1:
            for m in range(2, maxmn):
                for n in range(2, maxmn):
                    if gcd(m, n) == 1 and (m + n) * b * k - n * (b * l - a) == k:
                        out.append((k, (m, n)))
        k += 1
    return out

print("check F2   corner (7/5,2):", mn_solutions(7, 5, 2, 12))
print("check (8,28) corner (11/4,7):", mn_solutions(11, 4, 7, 36))
print("check (8,28) corner (7/4,3):", mn_solutions(7, 4, 3, 12))
print("check F9   corner (11/7,2):", mn_solutions(11, 7, 2, 12))
print()

# --- Part 2: enumerate all candidate degree pairs with 125 <= max <= 150 ---
# Families from GGV5 tables (lines 1678-1717); F18-F21 discarded structurally
fams = {
 "F1": ((4,12), lambda j:(2*j+3, 3*j+4)),
 "F2": ((5,20), lambda j:(j+2,   2*j+3)),
 "F3": ((5,20), lambda j:(4*j+3, 3*j+2)),
 "F4": ((5,20), lambda j:(2*j+3, 12*j+16)),
 "F5": ((5,20), lambda j:(7*j+9, 4*j+5)),
 "F6": ((5,20), lambda j:(3*j+4, 8*j+10)),
 "F7": ((6,15), lambda j:(j+2,   4*j+7)),
 "F8": ((6,15), lambda j:(2*j+3, 5*j+7)),
 "F9": ((7,21), lambda j:(j+2,   2*j+3)),
 "F10":((7,21), lambda j:(5*j+7, 3*j+4)),
 "F11":((7,21), lambda j:(j+2,   3*j+5)),
 "F12":((8,24), lambda j:(2*j+3, 5*j+7)),
 "F13":((9,21), lambda j:(j+2,   7*j+13)),
 "F14":((9,24), lambda j:(j+2,   4*j+7)),
 "F15":((9,24), lambda j:(2*j+3, 5*j+7)),
 "F16":((9,24), lambda j:(4*j+3, 7*j+5)),
 "F17":((9,24), lambda j:(5*j+2, 8*j+3)),
 "F22":((8,24), lambda j:(j+2,   2*j+3)),
 "F23":((8,24), lambda j:(j+2,   4*j+7)),
 "F24":((8,24), lambda j:(2*j+3, 3*j+4)),
}
cases = []
for name, (A0, mn) in fams.items():
    v = A0[0] + A0[1]
    for j in range(0, 20):
        m, n = mn(j)
        if m > 1 and n > 1 and gcd(m, n) == 1:
            mx = max(m, n) * v
            if mx <= 150:
                cases.append((mx, (m*v, n*v), name, A0, (m, n), f"j={j}"))

# Sporadic chains from GGV5 section "<=150" tables (lines 1821-1872).
# (A0, (m,n), chain length); (8,32)->120 and (9,27)->108 are BELOW 125 (in 2204 table)
sporadic = [
 ((7,35),(2,3),1), ((7,42),(3,2),1), ((7,42),(2,3),1),
 ((8,28),(3,4),1), ((8,28),(3,2),1),                      # (3,2) is the open (108,72)
 ((9,36),(3,2),1), ((9,36),(2,3),1), ((11,33),(2,3),1), ((12,33),(2,3),1),
 ((8,32),(3,2),2), ((8,40),(3,2),2), ((9,27),(2,3),2), ((9,36),(2,3),2),
 ((10,40),(3,2),2), ((10,40),(3,2),2), ((12,30),(3,2),2),
 ((12,36),(2,3),2), ((12,36),(2,3),2), ((12,36),(2,3),2), ((12,36),(2,3),2),
 ((12,36),(3,2),3),
]
for A0, (m, n), ln in sporadic:
    v = A0[0] + A0[1]
    mx = max(m, n) * v
    cases.append((mx, (m*v, n*v), f"chain(len {ln})", A0, (m, n), ""))

below = [c for c in cases if c[0] < 125]
band  = [c for c in cases if 125 <= c[0] <= 150]
print(f"cases with max<125: {len(below)}  (paper says 10)")
for c in sorted(below): print("   ", c)
print(f"cases with 125<=max<=150: {len(band)}  (34-10=24 expected)")
for c in sorted(band): print("   ", c)
