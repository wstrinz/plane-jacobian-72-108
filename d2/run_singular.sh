#!/bin/bash
# Resultant + factorization of the master identity; exports every factor.
Singular -q -c '
ring r=0,(d2,d1,d0,m1,P,x),dp;
string sa=read("Ain.txt"); execute("poly A="+sa+";");
string sb=read("Bin.txt"); execute("poly B="+sb+";");
poly R=resultant(A,B,x);
"resultant terms:"; size(R);
list L=factorize(R);
int i;
for(i=2;i<=size(L[1]);i++){
  write(":w factor_"+string(i)+".txt", string(L[1][i]));
  "factor",i,": terms",size(L[1][i]),"deg",deg(L[1][i]),"mult",L[2][i];
}
quit;'
echo "Map exported factors to f31 (102 terms, deg 31) / f37 (618 terms, deg 37)"
echo "by size; the single-variable factor with multiplicity 21 is m1 (debris)."
