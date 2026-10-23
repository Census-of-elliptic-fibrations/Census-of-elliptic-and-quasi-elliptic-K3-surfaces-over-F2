// functions_census.m

// Here are where we will program the functions that will output the data that we want to collect for the census of elliptic and quasi-elliptic K3 surfaces over F2.

// Useful existing functions already in the Magma environment, for E the generic fibre
// over K = F2(T) (e.g. curves[i] from mini_census.m).  Times are means per curve on the
// mini-census (timing_tests.m); they do not apply to the quasi-elliptic models, which
// are not elliptic curves in Magma's sense.
//
// LocalInformation(E):
//     Tate's algorithm at every place of bad reduction, the place at infinity included.
//     Returns a sequence of tuples <P, ord_P(Delta_min), ord_P(conductor), c_P, KS>:
//     the place P, the valuation of the minimal discriminant, the conductor exponent,
//     the Tamagawa number and the Kodaira symbol of the fibre.  A place of degree d
//     stands for d geometric fibres of that type, so the number of singular fibres is
//     sum Degree(P), and sum Degree(P)*ord_P(Delta_min) = 24 for a K3.  
//
// AnalyticInformation(E);
//     Information from the L-function of E.  Its first entry is the analytic rank, the
//     order of vanishing of L(E, s) at s = 1; this agreed with the multiplicity of the
//     root T = 1/2 of LFunction(E) on all of the mini-census.  It computes LFunction(E)
//     (about 0.2 s; the result is cached on E).  
//
// TorsionSubgroup(E);
//     The torsion subgroup of the Mordell-Weil group E(K), as an abstract abelian group
//     (Invariants(G) gives e.g. [2] for Z/2).  Over number fields it also returns a map
//     G -> E(K) giving the points (G, m := TorsionSubgroup(E)); not yet tried over F2(T).
//     About 8 ms.

// ---------------------------------------------------------------------------------------
// Zeta functions and Picard numbers of elliptic surfaces X -> P^1 over F_q.
//
// For the elliptic surface X with generic fibre E over K = F_q(T), b1 = b3 = 0 and
//
//     Z(X, u) = 1 / ( (1 - u) P2(u) (1 - q^2 u) ),
//
// where, by the Shioda-Tate decomposition of H^2,
//
//     P2(u) = (1 - q u)^2 * prod_v F_v(u) * L(E, u).
//
// The factor (1 - q u)^2 comes from the zero section and a fibre, L(E, u) is
// LFunction(E) (a polynomial in u = q^-s), and F_v(u) = det(1 - (q u)^deg(v) Frob_v)
// on the span of the non-identity components of the fibre at v.  Frobenius permutes
// those components; a cycle of length k contributes 1 - (q u)^(k deg v).  The
// permutation is fixed by the Kodaira type and the Tamagawa number c_v, which counts
// the components of multiplicity 1 defined over the residue field.
//
// The variable of the zeta function is u, to keep T free for F2(T) in mini_census.m.
// ---------------------------------------------------------------------------------------

Zu<u> := PolynomialRing(Integers());

function ComponentCycles(KS, c)
    // How Frobenius permutes the components of one singular fibre that do not meet the
    // zero section.
    //   Input:  KS, the Kodaira symbol as a string ("I4", "IV*", "I0*")
    //           c,  the Tamagawa number of the fibre (an integer, from LocalInformation)
    //   Output: a sequence of positive integers, one per Frobenius orbit on those
    //           components, giving the orbit's size; they add up to (#components - 1).
    //           E.g. split I4 -> [1, 1, 1], non-split I4 -> [1, 2], I1 or II -> [].
    ones := func< k | [ 1 : i in [1..k] ] >;
    if KS eq "II" then return []; end if;
    if KS eq "III" then return [1]; end if;
    if KS eq "IV" then return c eq 3 select [1, 1] else [2]; end if;
    if KS eq "IV*" then return c eq 3 select ones(6) else [1, 1, 2, 2]; end if;
    if KS eq "III*" then return ones(7); end if;
    if KS eq "II*" then return ones(8); end if;
    if KS[#KS] eq "*" then
        n := StringToInteger(Substring(KS, 2, #KS-2));  // I_n*: n + 4 non-identity components
        if n eq 0 then                             // the three outer ones are permuted
            if c eq 4 then return ones(4); end if;
            if c eq 2 then return [1, 1, 2]; end if;
            if c eq 1 then return [1, 3]; end if;
        else                                       // the two far ends may be swapped
            if c eq 4 then return ones(n+4); end if;
            if c eq 2 then return ones(n+2) cat [2]; end if;
        end if;
    else
        n := StringToInteger(Substring(KS, 2, #KS-1));  // I_n: an n-gon, split or reflected
        if n le 1 then return []; end if;
        if c eq n then return ones(n-1); end if;   // split (this also covers I2)
        if IsOdd(n) and c eq 1 then return [ 2 : i in [1..(n-1) div 2] ]; end if;
        if IsEven(n) and c eq 2 then return [1] cat [ 2 : i in [1..(n-2) div 2] ]; end if;
    end if;
    error "ComponentCycles: unexpected Tamagawa number", c, "for Kodaira type", KS;
end function;

function ZetaFunction(E)
    // The zeta function of the elliptic surface X with generic fibre E.
    //   Input:  E, an elliptic curve over F_q(T) (e.g. curves[i] from mini_census.m)
    //   Output: P2, the polynomial det(1 - u Frob | H^2(X)) in Z[u], of degree 22 for a K3
    //           Z,  the zeta function 1/((1 - u) P2(u) (1 - q^2 u)), a rational function in u
    //   Use P2, Z := ZetaFunction(E);  or just  P2 := ZetaFunction(E);
    q := #BaseRing(BaseRing(E));
    LI := LocalInformation(E);
    euler := &+[ Integers() | Degree(I[1])*I[2] : I in LI ];    // 12 chi(O_X); 24 for a K3

    P2 := (1 - q*u)^2 * (Zu ! Coefficients(LFunction(E)));
    for I in LI do
        d := Degree(I[1]);
        for k in ComponentCycles(Sprint(I[5]), I[4]) do
            P2 *:= 1 - (q*u)^(k*d);
        end for;
    end for;

    error if Degree(P2) ne euler - 2,
        "ZetaFunction: P2 has degree", Degree(P2), "but b2 =", euler - 2;
    Z := 1 / (FieldOfFractions(Zu) ! ((1 - u) * P2 * (1 - q^2*u)));
    return P2, Z;
end function;

function LFunctionSign(E)
    // The sign of the functional equation of the L-function of E.  Writing
    // L(E, u) = LFunction(E), of degree n, the functional equation is
    //     L(E, u) = eps * (q u)^n * L(E, 1/(q^2 u)),
    // and eps is the leading coefficient of L(E, u) divided by q^n.
    //   Input:  E, an elliptic curve over F_q(T) (e.g. curves[i] from mini_census.m)
    //   Output: eps, the integer 1 or -1
    q := #BaseRing(BaseRing(E));
    L := LFunction(E);
    eps := LeadingCoefficient(L) / q^Degree(L);
    error if eps notin {1, -1}, "LFunctionSign: leading coefficient is not +-q^n:", L;
    return Integers() ! eps;
end function;

function IsCyclotomicFactor(f)
    // Whether f is a cyclotomic polynomial (its roots are roots of unity).
    //   Input:  f, a monic irreducible polynomial over Q (a factor from Factorization)
    //   Output: true or false
    if not forall{ a : a in Coefficients(f) | IsIntegral(a) } then return false; end if;
    d := Degree(f);
    return exists{ n : n in [1..2*d^2+2] | EulerPhi(n) eq d and
                                           f eq Parent(f) ! CyclotomicPolynomial(n) };
end function;

function PicardNumber(E)
    // Picard numbers of the elliptic surface X with generic fibre E, read off from P2 by
    // the Tate conjecture (a theorem for elliptic K3 surfaces, by Artin and
    // Swinnerton-Dyer).
    //   Input:  E, an elliptic curve over F_q(T) (e.g. curves[i] from mini_census.m)
    //   Output: two integers,
    //           rho(X over the algebraic closure of F_q): the number of inverse roots of P2
    //               of the form q * (root of unity), with multiplicity
    //           rho(X over F_q): the multiplicity of the inverse root q itself
    //   Use geo, ar := PicardNumber(E);
    q := #BaseRing(BaseRing(E));
    P2 := ZetaFunction(E);
    QX<X> := PolynomialRing(Rationals());
    R := Evaluate(QX ! Coefficients(P2), X/q);          // inverse roots q*z  ->  z
    geometric := &+[ Integers() | Degree(F[1])*F[2] : F in Factorization(R) | IsCyclotomicFactor(F[1]) ];
    arithmetic := 0;
    while Evaluate(R, 1) eq 0 do
        R div:= X - 1;
        arithmetic +:= 1;
    end while;
    return geometric, arithmetic;
end function;


