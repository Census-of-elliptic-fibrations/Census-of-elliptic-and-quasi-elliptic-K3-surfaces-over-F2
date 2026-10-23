// Reads lists [a1, a2, a3, a4, a6] of polynomials in s, t from a text file and,
// for each list, builds the Weierstrass equation
//     y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
// and the elliptic curve it defines over the function field of P^1.

datafile := "elliptic-fib-census.txt";

// s and t must be defined before the file is read, so that eval understands them
P<s,t> := PolynomialRing(GF(2), 2);      // the coefficients a_i live in Q[s,t]
R<x,y> := PolynomialRing(P, 2);                // polynomials in x, y over Q[s,t]

data := [];   // the coefficient lists [a1, a2, a3, a4, a6]
eqs  := [];   // the equations, each stored as a polynomial f meaning f = 0

for chunk in Split(Read(datafile), "]") do     // each "]" ends one list
    i := Position(chunk, "[");
    if i gt 0 then
        L := eval (Substring(chunk, i, #chunk - i + 1) cat "]");
        coeffs := [P | c : c in L];
        Append(~data, coeffs);
        a1, a2, a3, a4, a6 := Explode(coeffs);
        Append(~eqs, y^2 + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x + a6));
    end if;
end for;
printf "Read %o equations.\n", #eqs;

// The same curves as elliptic curves over the function field of P^1.
// That field is Q(t), where t stands for t/s (the chart s = 1), so each
// coefficient is dehomogenised by setting s = 1. The point s = 0 is not lost:
// it is the place at infinity of Q(t). From here on, t means this field's variable.
// Magma's EllipticCurve takes the coefficients in the order [a1, a2, a3, a4, a6].
// (Magma stops with an error if a curve is singular.)
K<t> := FunctionField(GF(2));
curves := [EllipticCurve([K | Evaluate(c, [K | 1, t]) : c in coeffs]) : coeffs in data];

// quit;   // uncomment if you run this from the terminal with output sent to a file
