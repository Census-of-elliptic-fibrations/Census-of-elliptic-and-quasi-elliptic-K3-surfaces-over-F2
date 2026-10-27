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

// ---------------------------------------------------------------------------------------
// Core functions.  They work from the L-function L = LFunction(E) and the local data
// LI = LocalInformation(E), so that a census run computes each of those once per surface
// (LFunction is the expensive one, about 0.2 s).  The functions further down that take E
// are wrappers around these, for interactive use.
// ---------------------------------------------------------------------------------------

function P2FromData(L, LI, q)
    // The H^2 factor of the zeta function, from the L-function and the singular fibres.
    //   Input:  L,  the L-function LFunction(E), a polynomial with integer coefficients
    //           LI, the local data LocalInformation(E)
    //           q,  the size of the constant field (2 for the census)
    //   Output: P2, the polynomial det(1 - u Frob | H^2(X)) in Z[u], of degree 22 for a K3
    euler := &+[ Integers() | Degree(I[1])*I[2] : I in LI ];    // 12 chi(O_X); 24 for a K3
    P2 := (1 - q*u)^2 * (Zu ! Coefficients(L));
    for I in LI do
        d := Degree(I[1]);
        for k in ComponentCycles(Sprint(I[5]), I[4]) do
            P2 *:= 1 - (q*u)^(k*d);
        end for;
    end for;
    error if Degree(P2) ne euler - 2,
        "P2FromData: P2 has degree", Degree(P2), "but b2 =", euler - 2;
    return P2;
end function;

function SignOfPolynomial(f, q)
    // The sign eps of the functional equation f(u) = eps (q u)^n f(1/(q^2 u)), n = deg f,
    // for f = L(E, u) or f = P2(u): the leading coefficient of f divided by q^n.
    //   Input:  f, a polynomial with integer coefficients whose inverse roots have
    //              absolute value q (L(E, u) or P2(u))
    //           q, the size of the constant field
    //   Output: eps, the integer 1 or -1
    eps := LeadingCoefficient(f) / q^Degree(f);
    error if eps notin {1, -1}, "SignOfPolynomial: leading coefficient is not +-q^n:", f;
    return Integers() ! eps;
end function;

function RootMultiplicityAtInverseQ(f, q)
    // How often the inverse root q occurs in f, i.e. the multiplicity of the root u = 1/q.
    //   Input:  f, a nonzero polynomial with integer or rational coefficients
    //           q, the size of the constant field
    //   Output: two values: the multiplicity r (an integer), and f / (1 - q u)^r as a
    //           polynomial over Q
    QX<X> := PolynomialRing(Rationals());
    g := QX ! Coefficients(f);
    r := 0;
    while Evaluate(g, 1/q) eq 0 do
        g div:= 1 - q*X;
        r +:= 1;
    end while;
    return r, g;
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

function PicardNumbersFromP2(P2, q)
    // Picard numbers read off from P2 by the Tate conjecture (a theorem for elliptic K3
    // surfaces, by Artin and Swinnerton-Dyer).
    //   Input:  P2, the H^2 factor of the zeta function (from P2FromData)
    //           q,  the size of the constant field
    //   Output: two integers,
    //           rho(X over the algebraic closure of F_q): the number of inverse roots of P2
    //               of the form q * (root of unity), with multiplicity
    //           rho(X over F_q): the multiplicity of the inverse root q itself
    QX<X> := PolynomialRing(Rationals());
    R := Evaluate(QX ! Coefficients(P2), X/q);          // inverse roots q*z  ->  z
    geometric := &+[ Integers() | Degree(F[1])*F[2] : F in Factorization(R) | IsCyclotomicFactor(F[1]) ];
    arithmetic := RootMultiplicityAtInverseQ(P2, q);
    return geometric, arithmetic;
end function;

function ZetaFunction(E)
    // The zeta function of the elliptic surface X with generic fibre E.
    //   Input:  E, an elliptic curve over F_q(T) (e.g. curves[i] from mini_census.m)
    //   Output: P2, the polynomial det(1 - u Frob | H^2(X)) in Z[u], of degree 22 for a K3
    //           Z,  the zeta function 1/((1 - u) P2(u) (1 - q^2 u)), a rational function in u
    //   Use P2, Z := ZetaFunction(E);  or just  P2 := ZetaFunction(E);
    q := #BaseRing(BaseRing(E));
    P2 := P2FromData(LFunction(E), LocalInformation(E), q);
    Z := 1 / (FieldOfFractions(Zu) ! ((1 - u) * P2 * (1 - q^2*u)));
    return P2, Z;
end function;

function ZetaP2(E)
    return ZetaFunction(E);
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
    P2 := ZetaP2(E);
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


// Håvard wrote the next two functions, so it's artisanal slop (not AI)
function ArtinTateWrapper(P2)
    // Returns rho (non-geometric), and Br(X) * #NS(X)_disc
    T := Parent(P2).1; // I added this line to avoid a warning about T being undefined (Alvaro)
    _, val := ArtinTateFormula(Parent(P2)!(T^22*Evaluate(P2, 1/T)), 2, 1);
    return val;
end function;

// returns height h 
// NB: return 0 if supersingular (but really, h = infty)
function HeightP2(P2)
    T := Parent(P2).1; // I added this line to avoid a warning about T being undefined (Alvaro)
    lambda := Slopes(NewtonPolygon(ChangeRing(P2, pAdicRing(2, 500))))[1];
    try
	return 1/(1-lambda);
    catch e
	return 0;
    end try;
end function;

// Elise wrote the next function
function PointCounts(equ , n)
// assuming we have a K3 over F_2
//equ should already be initalized as an elliptic curve over a univariate
//rational function field over GF(2)
  Qt := PolynomialRing(Rationals());
  R<T> := PowerSeriesRing(Rationals(), n + 5);
  P0 := 1 - T ; P1 := 1 ; P1 := Qt ! P1 ;
  P2, Z := ZetaFunction(equ) ;
  P3 := 1 ; P3 := Qt ! P3 ; P4 := 1-4*T ;
   // convert a polynomial into an exact element of R
  ser := func< f | &+[ R | Coefficient(f, i) * T^i : i in [0..Degree(f)] ] >;
  Z := ser(P1) * ser(P3) / ( ser(P0) * ser(P2) * ser(P4) );
  L := T * Derivative(Z) / Z;
  return [ Integers() ! Coefficient(L, r) : r in [1..n] ];
end function ;

// This is a version of the previous function that takes the zeta function as input, rather than the elliptic curve. That way, we can avoid recomputing the zeta function if we already have it.
function PointCountsFromZeta(P2, n)
  Qt := PolynomialRing(Rationals());
  R<T> := PowerSeriesRing(Rationals(), n + 5);
  P0 := 1 - T ; P1 := 1 ; P1 := Qt ! P1 ;
  P3 := 1 ; P3 := Qt ! P3 ; P4 := 1-4*T ;
   // convert a polynomial into an exact element of R
  ser := func< f | &+[ R | Coefficient(f, i) * T^i : i in [0..Degree(f)] ] >;
  Z := ser(P1) * ser(P3) / ( ser(P0) * ser(P2) * ser(P4) );
  L := T * Derivative(Z) / Z;
  return [ Integers() ! Coefficient(L, r) : r in [1..n] ];
end function ;

// ---------------------------------------------------------------------------------------
// Census rows.
//
// CensusRow turns one model of the census into one line of the census CSV file, with the
// columns of CensusHeader.  It calls the expensive Magma functions once each
// (LocalInformation, LFunction, TorsionSubgroup) and computes every other column from
// their results with the functions above.  Lists inside a column use spaces, and the
// singular fibres use ";" between places and ":" between fields, so no column contains
// a comma.  PlaceForm is not used by CensusRow; it recovers the place of a fibre, for
// interactive use.
// ---------------------------------------------------------------------------------------

function PlaceForm(P, facs)
    // A place of F_2(T), T = t/s, written as an irreducible binary form in s, t.
    //   Input:  P,    a place of bad reduction (the first entry of a LocalInformation tuple)
    //           facs, the irreducible factors in F_2[T] of the numerator of the
    //                 discriminant of the model (every finite bad place is a zero of one)
    //   Output: a string: the homogenisation s^d g(t/s) of the factor g vanishing at P,
    //           e.g. "t", "s + t", "s^2 + s*t + t^2"; or "s" for the place at infinity
    F := FunctionField(P);
    S<s,t> := PolynomialRing(GF(2), 2);
    for g in facs do
        if Valuation(F ! g, P) gt 0 then
            d := Degree(g);
            error if d ne Degree(P), "PlaceForm: degree of the factor and of the place differ";
            return Sprint(&+[ S | Coefficient(g, i) * s^(d-i) * t^i : i in [0..d] ]);
        end if;
    end for;
    error if Degree(P) ne 1, "PlaceForm: no factor vanishes at a place of degree", Degree(P);
    return "s";
end function;

function FibreRecords(LI)
    // The singular fibres of the surface, one record <deg P, f_P, c_P, KS> per bad place.
    // The places themselves are left out: they depend on the choice of coordinates on P^1,
    // while these records do not, so two isomorphic surfaces give the same string.  The
    // places can be recovered from the model by factoring its discriminant (PlaceForm).
    //   Input:  LI, the local data LocalInformation(E) of the generic fibre E
    //   Output: two strings,
    //           the records "deg:f_P:c_P:KS" joined by ";", sorted by degree, then Kodaira
    //               type, then f_P, then c_P, e.g. "1:2:1:II;2:1:2:I2"
    //           the Kodaira types over the algebraic closure, each type repeated deg P
    //               times, in the same order and separated by spaces, e.g. "II I2 I2"
    //   It also checks Ogg's formula ord_P(Delta) = f_P + m_P - 1 at each place, where
    //   m_P is the number of components of the fibre, and stops with an error if it fails.
    recs := [];
    for I in LI do
        P, vD, f, c, KS := Explode(I);
        ks := Sprint(KS);
        m := 1 + &+[ Integers() | k : k in ComponentCycles(ks, c) ];
        error if vD ne f + m - 1, "FibreRecords: Ogg's formula fails at a fibre of type", ks;
        Append(~recs, < Degree(P), ks, f, c >);
    end for;
    cmp := func< a, b | a lt b select -1 else (a eq b select 0 else 1) >;
    Sort(~recs, func< x, y | x[1] ne y[1] select x[1] - y[1] else
                             (x[2] ne y[2] select cmp(x[2], y[2]) else
                             (x[3] ne y[3] select x[3] - y[3] else x[4] - y[4])) >);
    records := Join([ Sprint(r[1]) cat ":" cat Sprint(r[3]) cat ":" cat Sprint(r[4]) cat ":" cat r[2]
                      : r in recs ], ";");
    types := Join(&cat[ [ r[2] : j in [1..r[1]] ] : r in recs ], " ");
    return records, types;
end function;

function TorsionLabel(G)
    // The torsion subgroup of the Mordell-Weil group as a short label.
    //   Input:  G, an abelian group (the first value of TorsionSubgroup(E))
    //   Output: a string, "trivial", "Z/2", "Z/2 x Z/4", ...
    inv := Invariants(G);
    return #inv eq 0 select "trivial" else Join([ "Z/" cat Sprint(n) : n in inv ], " x ");
end function;

function CensusHeader()
    // The first line of the census CSV file.
    //   Input:  none
    //   Output: a string, the column names separated by commas
    return "key,a1,mw_torsion,singular_fibres,L_function,L_sign," cat
           "analytic_rank,P2,zeta_sign,rho_F2,rho_geom,artin_height,br_disc," cat
           "pts_F2,pts_F4,pts_F8,pts_F16,pts_F32";
end function;

function CensusRow(key, A)
    // One line of the census CSV file, for one elliptic K3 surface.
    //   Input:  key, the integer orbit key of the model (from FibrationFiberOrbitKeysF2)
    //           A,   the model [a1, a2, a3, a4, a6] as binary forms in s, t over F_2,
    //                a2 = 0 (from FibrationFromKeyF2)
    //   Output: a string with the columns of CensusHeader, separated by commas:
    //           key                         the model (FibrationFromKeyF2 gives it back)
    //           a1                          its a1: 0, s^2, s*t or s^2 + s*t + t^2
    //           mw_torsion                  TorsionLabel of the torsion of E(F_2(T))
    //           singular_fibres             FibreRecords, "deg:f_P:c_P:KS;..."
    //           L_function, P2              coefficients, constant term first, by spaces
    //           L_sign, zeta_sign           signs of the functional equations of L and P2
    //           analytic_rank               multiplicity of the inverse root 2 of L
    //           rho_F2, rho_geom            Picard numbers over F_2 and its closure
    //           artin_height                height of the formal Brauer group, "inf" if
    //                                       supersingular (HeightP2)
    //           br_disc                     #Br(X) * |disc NS(X)| (ArtinTateWrapper)
    //           pts_F2, ..., pts_F32        #X(F_2^r) for r = 1..5 (PointCountsFromZeta)
    //   It stops with an error if a check fails (Ogg's formula, deg P2 = 22, the signs).
    K := RationalFunctionField(GF(2));
    E := EllipticCurve([ K | Evaluate(f, [K!1, K.1]) : f in A ]);

    LI := LocalInformation(E);
    L := LFunction(E);
    G := TorsionSubgroup(E);

    records := FibreRecords(LI);
    P2 := P2FromData(L, LI, 2);
    rank := RootMultiplicityAtInverseQ(L, 2);
    rho_geom, rho_F2 := PicardNumbersFromP2(P2, 2);
    h := HeightP2(P2);
    pts := PointCountsFromZeta(P2, 5);
    coeffs := func< f | Join([ Sprint(a) : a in Coefficients(f) ], " ") >;

    return Join([ Sprint(key), Sprint(A[1]),
                  TorsionLabel(G), records,
                  coeffs(L), Sprint(SignOfPolynomial(L, 2)), Sprint(rank),
                  coeffs(P2), Sprint(SignOfPolynomial(P2, 2)),
                  Sprint(rho_F2), Sprint(rho_geom),
                  h eq 0 select "inf" else Sprint(h), Sprint(ArtinTateWrapper(P2)) ]
                cat [ Sprint(n) : n in pts ], ",");
end function;
