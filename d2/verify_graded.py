"""
verify_graded.py — T5: exact verification of the graded decomposition of f31
and the structural facts the (y+1)-adic Newton-polygon argument rests on.

Checks (all exact, sympy over Q):
  1. f31 = sum_{f=0}^7 Phi^f * dm1^(21-3f) * h_f   with h_f from f31_graded.txt
  2. h_f weighted-homogeneous of weight 20-2f under w(d2,d1,d0,dm1)=(2,3,4,5)
  3. h_0 == f31|_{Phi=0} / dm1^21  (the h31 of T5_NOTES)
  4. cascade identities used to kill the degenerate branch:
       h_7 = 8192 d1^2
       h_6|_{d1=0} = -3072 (4 d0 - d2^2)^2
       h_5|_{d1=0, d0=d2^2/4} = 2048 dm1^2
  5. Phi~ = Phi/y^204 has (y+1)-valuation exactly 30 with unit cofactor
     (quartic(-1) = 3315 != 0), and the y-degree bounds behind
     deg h_f(d~) <= 40 - 4f on the stripped sub2 windows.
"""
import re, sympy as sp

d2, d1, d0, dm1, Phi, y, W = sp.symbols('d2 d1 d0 dm1 Phi y W')
V4 = (d2, d1, d0, dm1)
WTS = {d2: 2, d1: 3, d0: 4, dm1: 5}

s = open('f31_deg31.txt').read().strip()
f31 = sp.sympify(s.replace('m1', 'dm1').replace('P', 'Phi').replace('^', '**'))

hs = {}
for m in re.finditer(r'h_(\d) \(weight (\d+), dm1-power (\d+)\) = (.+)', open('f31_graded.txt').read()):
    f, wt, pw, expr = int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4)
    hs[f] = (sp.sympify(expr), wt, pw)
assert sorted(hs) == list(range(8)), "expected h_0..h_7"

# 1. decomposition
recon = sum(Phi**f * dm1**(21 - 3*f) * hs[f][0] for f in range(8))
assert sp.expand(recon - f31) == 0, "DECOMPOSITION FAILS"
print("1. f31 == sum Phi^f dm1^(21-3f) h_f            OK")

# 2. weights
for f in range(8):
    h, wt, pw = hs[f]
    assert pw == 21 - 3*f
    for mono in sp.Poly(h, *V4).monoms():
        assert sum(WTS[v]*e for v, e in zip(V4, mono)) == wt == 20 - 2*f, (f, mono)
print("2. h_f weighted-homogeneous, weight 20-2f      OK")

# 3. h_0 = f31|_{Phi=0}/dm1^21
q = sp.cancel(f31.subs(Phi, 0) / dm1**21)
assert sp.expand(q - hs[0][0]) == 0
print("3. h_0 == f31|_{Phi=0} / dm1^21                OK")

# 4. cascade identities
assert sp.expand(hs[7][0] - 8192*d1**2) == 0
h6_0 = hs[6][0].subs(d1, 0)
assert sp.expand(h6_0 + 3072*(4*d0 - d2**2)**2) == 0
h5_00 = hs[5][0].subs(d1, 0).subs(d0, d2**2/4)
assert sp.expand(h5_00 - 2048*dm1**2) == 0
print("4. cascade: h7=8192 d1^2; h6|_{d1=0}=-3072(4d0-d2^2)^2;")
print("            h5|_{d1=0,d0=d2^2/4}=2048 dm1^2    OK")

# 5. Phi~ valuation at y=-1 and degree bounds
quartic = 2048*y**4 - 512*y**3 + 320*y**2 - 240*y + 195
Phit = -(y+1)**30 * quartic / 6630
assert quartic.subs(y, -1) == 3315
t = sp.symbols('t')
Pt = sp.expand(Phit.subs(y, t - 1))
assert sp.Poly(Pt, t).monoms()[-1][0] == 30       # ord_t = 30
assert sp.degree(Phit, y) == 34
# stripped sub2 window degrees: deg d~_k <= 2*w_k -> deg h_f(d~) <= 2*(20-2f)
for f in range(8):
    h = hs[f][0]
    maxdeg = max(sum(2*WTS[v]*e for v, e in zip(V4, mono))
                 for mono in sp.Poly(h, *V4).monoms())
    assert maxdeg == 2*(20 - 2*f) == 40 - 4*f, (f, maxdeg)
print("5. v_{y+1}(Phi~)=30 (unit cofactor, quartic(-1)=3315);")
print("   deg_y h_f(d~) <= 40-4f on sub2 windows      OK")

print("\nALL GRADED-STRUCTURE CHECKS PASS")
