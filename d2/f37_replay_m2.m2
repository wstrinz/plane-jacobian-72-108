-- Independent Macaulay2 replay of the f37 pre-resultant elimination theorem.
--
-- This script deliberately does not read generators.json, factor_*.txt,
-- f31_graded.txt, a pickle, a certificate, or generated CAS input.  It builds
-- the four generators by formal-series convolution and separately rebuilds
-- f31 through the classical resultant chain.

needsPackage "Elimination";

report = (label, ok) -> (
    if ok then print("PASS: " | label) else print("FAIL: " | label);
    ok
    );

-- -------------------------------------------------------------------------
-- 1. Build the pre-resultant generators from the mathematical definition.
-- -------------------------------------------------------------------------
A = QQ[d2,d1,d0,dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10,dm11,dm12,dm13,Phi];
dmVars = {dm1,dm2,dm3,dm4,dm5,dm6,dm7,dm8,dm9,dm10,dm11,dm12,dm13};

-- Coefficient of u^n in
-- S = 1 + d2*u^2 + d1*u^3 + d0*u^4 + sum_(k=1)^13 dm_k*u^(4+k).
sCoeff = n -> (
    if n == 0 then 1_A
    else if n == 2 then d2
    else if n == 3 then d1
    else if n == 4 then d0
    else if n >= 5 and n <= 17 then dmVars#(n-5)
    else 0_A
    );

coeffSquare = n -> sum(0..n, i -> sCoeff(i) * sCoeff(n-i));
coeffCube = n -> sum(0..n, i ->
    sum(0..(n-i), j -> sCoeff(i) * sCoeff(j) * sCoeff(n-i-j))
    );

-- Published t=3 checkpoint: [u^7]S3^2.  This checks the indexing convention
-- independently before the t=4 construction is used.
t3Coeff = n -> (
    if n == 0 then 1_A
    else if n == 2 then d1
    else if n == 3 then d0
    else if n >= 4 and n <= 13 then dmVars#(n-4)
    else 0_A
    );
t3u7 = sum(0..7, i -> t3Coeff(i) * t3Coeff(7-i));
t3Expected = 2*d0*dm1 + 2*d1*dm2 + 2*dm4;
t3OK = (t3u7 == t3Expected);

-- D2(k)=[u^(8+k)]S^2 is linear in the displayed fresh variable.  Solve the
-- eight equations successively, exactly as the mathematical linear phase.
linSubs = {};
linearOK = true;
steps = {
    {1,dm5},{2,dm6},{3,dm7},{4,dm8},
    {5,dm9},{6,dm10},{7,dm11},{9,dm13}
    };
scan(steps, stp -> (
    k := stp#0;
    z := stp#1;
    e := sub(coeffSquare(8+k), linSubs);
    c := diff(z,e);
    linearOK = linearOK and (c != 0_A) and (diff(z,c) == 0_A)
        and (diff(z,diff(z,e)) == 0_A);
    sol := -sub(e,{z=>0_A}) / c;
    linSubs = append(linSubs,z=>sol);
    ));

G1A = sub(coeffCube(13),linSubs);
G2A = sub(coeffCube(14),linSubs);
G3A = sub(coeffCube(15),linSubs);
G5bodyA = sub(coeffCube(17),linSubs);

-- A short human-readable formula check catches a wrong convolution index or
-- wrong linear slice without making these formulas the construction itself.
expectedG1 = 3/2*d1*dm1^2 + 3*d2*dm1*dm2 + 3*dm1*dm4 + 3*dm2*dm3;
expectedG2 = -3/2*d0*dm1^2 + 3/2*d2*dm2^2 + 3*dm2*dm4 + 3/2*dm3^2;
expectedG3 = -3*d0*dm1*dm2 - 3/2*d1*dm2^2 - 1/2*dm1^3 + 3*dm3*dm4;
expectedG5body = -3*d0*dm1*dm4 - 3*d0*dm2*dm3
    - 3*d1*dm2*dm4 - 3/2*d1*dm3^2 - 3*d2*dm3*dm4
    - 3/2*dm1^2*dm3 - 3/2*dm1*dm2^2;
constructionOK = t3OK and linearOK
    and (G1A == expectedG1) and (G2A == expectedG2)
    and (G3A == expectedG3) and (G5bodyA == expectedG5body);
report("formal-series construction of G1,G2,G3,G5body",constructionOK);

-- Retain only the eight variables used by the theorem.  The eliminated
-- variables come first in a native elimination order.
R = QQ[dm2,dm3,dm4,d2,d1,d0,dm1,Phi,MonomialOrder=>Eliminate 3];
toR = map(R,A,{
    d2,d1,d0,dm1,dm2,dm3,dm4,
    0_R,0_R,0_R,0_R,0_R,0_R,0_R,0_R,0_R,Phi
    });
G1 = toR G1A;
G2 = toR G2A;
G3 = toR G3A;
G5body = toR G5bodyA;
G5 = G5body + Phi;
I = ideal(G1,G2,G3,G5);

-- -------------------------------------------------------------------------
-- 2. Independently regenerate f31 by the historical resultant route.
-- -------------------------------------------------------------------------
-- These denominator-cleared combinations are obtained by solving G1 for dm4.
-- They eliminate dm4 without importing sol4 or a serialized H-system.
H2 = dm1*G2 - dm2*G1;
H3 = dm1*G3 - dm3*G1;
H5 = dm1*G5 + (d0*dm1+d1*dm2+d2*dm3)*G1;
assert(diff(dm4,H2) == 0_R);
assert(diff(dm4,H3) == 0_R);
assert(diff(dm4,H5) == 0_R);

factorBases = f -> (
    FF := factor f;
    apply(toList(0..(#FF-1)), i -> FF#i#0)
    );

print("INFO: computing the two dm3 resultants");
RA = resultant(H2,H3,dm3);
RB = resultant(H2,H5,dm3);
AhCandidates = select(factorBases RA, q -> diff(dm2,q) != 0_R);
BhCandidates = select(factorBases RB, q -> diff(dm2,q) != 0_R);
resultantInputsOK = (#AhCandidates == 1) and (#BhCandidates == 1);
report("unique dm2-bearing factors in the two intermediate resultants",
    resultantInputsOK);
if not resultantInputsOK then exit 1;
Ah = first AhCandidates;
Bh = first BhCandidates;

print("INFO: computing and factoring the final dm2 resultant");
master = resultant(Ah,Bh,dm2);
f31Candidates = select(factorBases master, q ->
    (first degree q == 31) and (diff(Phi,q) != 0_R));
f31Unique = (#f31Candidates == 1);
report("unique Phi-bearing total-degree-31 resultant factor",f31Unique);
if not f31Unique then exit 1;
f31 = first f31Candidates;

-- -------------------------------------------------------------------------
-- 3. Compute the true elimination ideal and replay the theorem.
-- -------------------------------------------------------------------------
print("INFO: computing I intersect QQ[d2,d1,d0,dm1,Phi]");
E = eliminate(I,{dm2,dm3,dm4});
EBasisIdeal = ideal(gens(gb E));
principalOK = (numgens EBasisIdeal == 1) and (EBasisIdeal_0 != 0_R);
eGenerator = if principalOK then EBasisIdeal_0 else 0_R;
matchOK = principalOK and (ideal(eGenerator) == ideal(f31));
membershipOK = ((f31 % (gb I)) == 0_R);

okA = report(
    "(a) elimination ideal in (d2,d1,d0,dm1,Phi) is principal",
    principalOK);
okB = report(
    "(b) its generator equals independently reconstructed f31 up to a QQ unit",
    matchOK);
okC = report(
    "(c) independently reconstructed f31 lies in the pre-resultant ideal",
    membershipOK);

print("INFO: elimination generators=" | toString(numgens EBasisIdeal)
    | ", f31 total degree=" | toString(first degree f31)
    | ", f31 terms=" | toString(numcols(monomials f31)));
allOK = constructionOK and resultantInputsOK and f31Unique and okA and okB and okC;
if allOK then (
    print("PASS: f37 theorem replay complete over QQ");
    exit 0
    ) else exit 1;

