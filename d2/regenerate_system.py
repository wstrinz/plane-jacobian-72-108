"""Regenerate the t=4 system, linear phase, and Singular inputs from scratch.
Validates the generator against the published (7,21) t=3 system first."""
import sympy as sp, pickle
from sympy import symbols, expand, factor_list, resultant, Poly, together, numer

Phi = symbols('Phi'); u = symbols('u')
d2,d1,d0 = symbols('d2 d1 d0')
dm = {k: symbols(f'dm{k}') for k in range(1,14)}

# --- t=3 validation (must reproduce arXiv:2204.14178 sec.6 equations) ---
S3v = 1 + d1*u**2 + d0*u**3 + sum(dm[k]*u**(3+k) for k in range(1,11))
P2 = Poly(expand(S3v*S3v), u)
chk = P2.coeff_monomial(u**7)  # (D~^2)_{-1} for t=3
assert expand(chk - (2*d0*dm[1] + 2*d1*dm[2] + 2*dm[4])) == 0, "t=3 validation failed"
print("t=3 generator validation: OK")

# --- t=4 system ---
S = 1 + d2*u**2 + d1*u**3 + d0*u**4 + sum(dm[k]*u**(4+k) for k in range(1,14))
S2 = Poly(expand(S*S), u); S3 = Poly(expand(S2.as_expr()*S), u)
D2 = lambda k: S2.coeff_monomial(u**(8+k))
D3 = lambda j: S3.coeff_monomial(u**(12+j))
used = [D2(k) for k in [1,2,3,4,5,6,7,9]] + [D3(1),D3(2),D3(3),D3(5)]
assert not any(e.has(dm[12]) for e in used)
sub = {}
for k, fresh in [(1,dm[5]),(2,dm[6]),(3,dm[7]),(4,dm[8]),(5,dm[9]),(6,dm[10]),(7,dm[11]),(9,dm[13])]:
    sub[fresh] = expand(sp.solve(D2(k).subs(sub), fresh)[0])
G1 = expand(D3(1).subs(sub)); G2 = expand(D3(2).subs(sub))
G3 = expand(D3(3).subs(sub)); G5b = expand(D3(5).subs(sub))
sol4 = sp.solve(G1, dm[4])[0]
H = [expand(numer(together(e.subs(dm[4], sol4)))) for e in (G2, G3, G5b + Phi)]
A = expand(sp.prod(f for f,_ in factor_list(resultant(Poly(H[0],dm[3]),Poly(H[1],dm[3])))[1]))
B = expand(sp.prod(f for f,_ in factor_list(resultant(Poly(H[0],dm[3]),Poly(H[2],dm[3])))[1]))
pickle.dump(dict(G1=G1,G2=G2,G3=G3,G5body=G5b,sol4=sol4,H2=H[0],H3=H[1],A=A,B=B),
            open('t4_state.pkl','wb'))
Ah = [f for f,_ in factor_list(A)[1] if f.has(dm[2])][0]
Bh = [f for f,_ in factor_list(B)[1] if f.has(dm[2])][0]
def to_sing(e):
    return sp.sstr(expand(e)).replace('**','^').replace('dm1','m1').replace('dm2','x').replace('Phi','P')
open('Ain.txt','w').write(to_sing(Ah)); open('Bin.txt','w').write(to_sing(Bh))
print("Singular inputs written (Ain.txt, Bin.txt); state saved (t4_state.pkl)")
