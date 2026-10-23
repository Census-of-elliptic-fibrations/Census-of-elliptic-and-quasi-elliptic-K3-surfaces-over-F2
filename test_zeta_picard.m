// test_zeta_picard.m
//
// Checks ZetaFunction and PicardNumber from functions_census.m on the first N curves of
// mini_census.m, and prints the distribution of the Picard numbers.  Each P2 must
//   - have degree 22 and constant term 1,
//   - satisfy the functional equation  u^22 P2(1/(q^2 u)) = +- P2(u) / q^22,
// and the Picard numbers must satisfy 2 <= rho(X) <= rho(X over F2bar) <= 22, with
// rho(X over F2bar) even (as for every K3 surface over a finite field).
//
// Run with  magma -b test_zeta_picard.m  from the repository root.  About 0.2 s a curve.

N := 100;

load "mini_census.m";
load "functions_census.m";

q := 2;
geometric := [];
arithmetic := [];
start := Cputime();
for i in [1..N] do
    E := curves[i];
    P2, Z := ZetaFunction(E);
    a := Coefficients(P2) cat [ 0 : j in [1..23 - #Coefficients(P2)] ];   // a[k+1] = coeff of u^k
    assert Degree(P2) eq 22 and a[1] eq 1;
    sign := a[23] div q^22;
    assert sign in {1, -1};
    assert forall{ k : k in [0..22] | a[23-k] eq sign * q^(22-2*k) * a[k+1] };

    g, r := PicardNumber(E);
    assert 2 le r and r le g and g le 22 and IsEven(g);
    Append(~geometric, g);
    Append(~arithmetic, r);
end for;
printf "%o curves checked in %o s\n", N, Cputime(start);

printf "rho over F2bar: %o\n", [ <g, #[ x : x in geometric | x eq g ]> : g in Sort(SetToSequence(Seqset(geometric))) ];
printf "rho over F2:    %o\n", [ <r, #[ x : x in arithmetic | x eq r ]> : r in Sort(SetToSequence(Seqset(arithmetic))) ];
printf "First curve: P2 = %o\n", ZetaP2(curves[1]);
