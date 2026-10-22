// census.m
//
// An equation  y^2 + a1 xy + a3 y = x^3 + a4 x + a6  (a2 = 0) is stored as a vector in
// A^3 x A^7 x A^9 x A^13 over F2, i.e. an element of F2^32:
//
//     v = ( a1 | a3 | a4 | a6 ),   a_i = a_{i,0} s^(2i) + a_{i,1} s^(2i-1) t + ... + a_{i,2i} t^(2i),
//
// so that a_{i,j} (the coefficient of s^(2i-j) t^j) is entry j+1 of the block of a_i.
// Block positions in v:  a1 = 1..3,  a3 = 4..10,  a4 = 11..19,  a6 = 20..32.

F := GF(2);
P<s,t> := PolynomialRing(F, 2);
V := VectorSpace(F, 32);

DEG    := [2, 6, 8, 12];    // degrees of a1, a3, a4, a6
OFFSET := [0, 3, 10, 19];   // block k occupies OFFSET[k]+1 .. OFFSET[k]+DEG[k]+1

function FormFromSeq(c, d)
    // [a_{i,0}, ..., a_{i,d}]  ->  sum_j c[j+1] s^(d-j) t^j
    return &+[ P | c[j+1]*s^(d-j)*t^j : j in [0..d] ];
end function;

function SeqFromForm(f, d)
    return [ F | MonomialCoefficient(f, s^(d-j)*t^j) : j in [0..d] ];
end function;

function Forms(v)
    // v  ->  [a1, a3, a4, a6]
    return [ FormFromSeq([ v[OFFSET[k]+i] : i in [1..DEG[k]+1] ], DEG[k]) : k in [1..4] ];
end function;

function Vec(A)
    // [a1, a3, a4, a6]  ->  v
    return V ! &cat[ SeqFromForm(A[k], DEG[k]) : k in [1..4] ];
end function;

function ReduceA1A3Zero(v)
    // For a1 = a3 = 0: apply x -> x + w2^2, y -> y + w2 x + w6 with
    //   w2 = a_{4,0} s^2 + a_{4,4} st + a_{4,8} t^2   (removes s^8, s^4t^4, t^8 from a4: a4' = a4 + w2^4)
    //   w6 = sqrt of the even part of a6 + w2^2 a4 + w2^6   (removes s^12, s^10t^2, ..., t^12 from a6)
    // The coefficients of the result are polynomials in those of v.
    a1, a3, a4, a6 := Explode(Forms(v));
    error if a1 ne 0 or a3 ne 0, "ReduceA1A3Zero needs a1 = a3 = 0";
    c := SeqFromForm(a4, 8);
    w2 := c[1]*s^2 + c[5]*s*t + c[9]*t^2;
    b4 := a4 + w2^4;
    b6 := a6 + w2^2*a4 + w2^6;
    e := SeqFromForm(b6, 12);
    w6 := &+[ P | e[2*k+1]*s^(6-k)*t^k : k in [0..6] ];
    b6 := b6 + w6^2;
    return Vec([ P!0, P!0, b4, b6 ]);
end function;

K<T> := RationalFunctionField(F);   // F2(t), with T = t/s (the chart s = 1)

function Curve(v)
    // v  ->  y^2 + a1(1,T) xy + a3(1,T) y = x^3 + a4(1,T) x + a6(1,T)  over F2(T)
    a1, a3, a4, a6 := Explode([ Evaluate(f, [K!1, T]) : f in Forms(v) ]);
    return EllipticCurve([ K | a1, 0, a3, a4, a6 ]);
end function;

// ---------------------------------------------------------------------------------------
// Examples.  Set RUN_EXAMPLES := true (or paste the block into a Magma session) to run them.
// ---------------------------------------------------------------------------------------

RUN_EXAMPLES := true;

if RUN_EXAMPLES then

    // 1. Encoding: forms <-> vectors.
    A := [ s^2 + t^2, s^6 + s^3*t^3 + t^6, s^8 + s^5*t^3 + t^8, s^11*t + s^4*t^8 + t^12 ];
    v := Vec(A);
    assert Forms(v) eq A;
    printf "A = %o\n", A;
    printf "v = %o\n", v;

    // 2. An elliptic fibration (a1 ne 0): the generic fibre over F2(T).
    E := Curve(v);
    D := Discriminant(E);
    printf "E: %o\n", E;
    printf "j(E) = %o\n", jInvariant(E);

    // Singular fibres: LocalInformation(E) gives, for each bad place P of F2(T) (the place
    // at infinity included), the tuple <P, ord_P(Delta_min), ord_P(conductor), c_P, Kodaira type>.
    // A place of degree d accounts for d geometric fibres of the same type.
    eu := 0;
    for I in LocalInformation(E) do
        P, vD, vN, cP, KS := Explode(I);
        printf "  %o (deg %o): %o, ord(Delta) = %o, ord(N) = %o, c = %o\n",
               P, Degree(P), KS, vD, vN, cP;
        eu +:= Degree(P)*vD;
    end for;
    printf "Euler number = %o (24 for a K3)\n", eu;

    // 3. A quasi-elliptic fibration (a1 = a3 = 0): reduce a random vector to normal form.
    //    Curve() does not apply here, since y^2 = x^3 + a4 x + a6 is singular in char 2.
    w := V ! ([ F | 0 : i in [1..10] ] cat [ Random(F) : i in [11..32] ]);
    r := ReduceA1A3Zero(w);
    _, _, a4, a6 := Explode(Forms(r));
    assert forall{ m : m in [s^8, s^4*t^4, t^8] | MonomialCoefficient(a4, m) eq 0 };
    assert forall{ k : k in [0..6] | MonomialCoefficient(a6, s^(12-2*k)*t^(2*k)) eq 0 };
    assert ReduceA1A3Zero(r) eq r;      // the normal form is stable
    printf "w = %o\n  -> r = %o\n", w, r;
    printf "y^2 = x^3 + (%o) x + (%o)\n", a4, a6;

end if;
