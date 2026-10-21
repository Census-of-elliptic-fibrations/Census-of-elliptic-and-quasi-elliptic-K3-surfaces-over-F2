// Orbits of Weierstrass models y^2+a1*x*y+a3*y = x^3+a2*x^2+a4*x+a6 over
// F2[s,t], deg a_i = 2i, under GL(2,F2) on (s,t) and x -> x+q, y -> y+p*x+r.
// Minimality is not tested here; see k3_minimality_f2.m.
// census_driver.m shows usage.

function EF2BinaryBasis(R, d)
    s := R.1;
    t := R.2;
    return [ R | s^(d-i)*t^i : i in [0..d] ];
end function;

function EF2CheckModel(A)
    assert #A eq 5;
    R := Parent(A[1]);
    assert Type(R) eq RngMPol;
    assert Rank(R) eq 2;
    k := BaseRing(R);
    assert Type(k) eq FldFin;
    assert #k eq 2;
    degrees := [2,4,6,8,12];
    for i in [1..5] do
        assert Parent(A[i]) eq R;
        for m in Monomials(A[i]) do
            assert TotalDegree(m) eq degrees[i];
        end for;
    end for;
    return R;
end function;

function EF2Coefficients(f, d)
    return [ MonomialCoefficient(f, m)
             : m in EF2BinaryBasis(Parent(f), d) ];
end function;

function EF2BitKey(bits)
    key := 0;
    for b in bits do
        key := 2*key + (Integers()!b);
    end for;
    return key;
end function;

function EF2PairVector(c, d)
    k := BaseRing(Parent(c));
    return Vector(k, EF2Coefficients(c, 8) cat EF2Coefficients(d, 12));
end function;

function EF2PairFromVector(R, v)
    B8 := EF2BinaryBasis(R, 8);
    B12 := EF2BinaryBasis(R, 12);
    c := &+[ v[i]*B8[i] : i in [1..9] ];
    d := &+[ v[9+i]*B12[i] : i in [1..13] ];
    return c, d;
end function;

function NormalizeA2F2(A)
    // Substitute x -> x+a2, with y unchanged.
    R := EF2CheckModel(A);
    a1,a2,a3,a4,a6 := Explode(A);
    return [ R | a1, 0, a3+a1*a2, a4+a2^2, a6+a2*a4 ];
end function;

function FibrationDiscriminantF2(A)
    // A binary form of degree 24, nonzero iff the generic fibre is smooth.
    R := EF2CheckModel(A);
    a1,a2,a3,a4,a6 := Explode(A);
    return a3^4+a1^3*a3^3+a1^4*a4^2+a1^5*a3*a4
           +a1^6*a6+a1^4*a2*a3^2;
end function;

function TransformNormalizedF2(A, M, p, r)
    // First substitute (s,t)^T by M*(s,t)^T, then x -> x+q, y -> y+p*x+r
    // with p, r in the new base variables. Taking q = p^2+a1*p is exactly
    // the condition that the new a2 is 0.
    R := EF2CheckModel(A);
    assert A[2] eq 0;
    k := BaseRing(R);
    assert BaseRing(M) eq k;
    assert Nrows(M) eq 2 and Ncols(M) eq 2;
    assert Determinant(M) ne 0;
    p := R!p;
    r := R!r;
    for m in Monomials(p) do
        assert TotalDegree(m) eq 2;
    end for;
    for m in Monomials(r) do
        assert TotalDegree(m) eq 6;
    end for;
    s := R.1;
    t := R.2;
    st := [ M[1,1]*s+M[1,2]*t, M[2,1]*s+M[2,2]*t ];
    B := [ R | Evaluate(a, st) : a in A ];
    a1 := B[1];
    a3 := B[3];
    a4 := B[4];
    a6 := B[5];
    q := p^2+a1*p;
    b3 := a3+a1*q;
    b4 := a4+q^2+p*b3+a1*r;
    b6 := a6+a4*q+q^3+r^2+b3*r;
    return [ R | a1, 0, b3, b4, b6 ];
end function;

function CanonicalFibrationF2(A)
    // Returns the canonical model, with a2 = 0, and its key in [0,2^32). The
    // key is the coefficient bits of a1, a3, a4, a6, most significant first,
    // each in the order s^d, s^(d-1)*t, ..., t^d.
    E := NormalizeA2F2(A);
    R := Parent(E[1]);
    k := BaseRing(R);
    V22 := VectorSpace(k, 22);
    Vp := VectorSpace(k, 3);
    B2 := EF2BinaryBasis(R, 2);
    B6 := EF2BinaryBasis(R, 6);
    spaces := AssociativeArray(Integers());
    bestKey := 2^32;
    best := E;
    for g in GL(2, k) do
        M := Matrix(g);
        for v in Vp do
            p := &+[ v[i]*B2[i] : i in [1..3] ];
            B := TransformNormalizedF2(E, M, p, R!0);
            prefix := EF2Coefficients(B[1], 2)
                      cat EF2Coefficients(B[3], 6);
            prefixKey := EF2BitKey(prefix);
            found, S := IsDefined(spaces, prefixKey);
            if not found then
                // y -> y+r adds (a1*r, r^2+a3*r) to (a4,a6).
                // This map is GF(2)-linear, including its square term.
                S := sub< V22 |
                    [ V22!EF2PairVector(B[1]*r, r^2+B[3]*r)
                      : r in B6 ] >;
                spaces[prefixKey] := S;
            end if;
            w := ReduceVector(S, V22!EF2PairVector(B[4], B[5]));
            key := EF2BitKey(prefix cat Eltseq(w));
            if key lt bestKey then
                c, d := EF2PairFromVector(R, w);
                best := [ R | B[1], 0, B[3], c, d ];
                bestKey := key;
            end if;
        end for;
    end for;
    return best, bestKey;
end function;

function PartitionFibrationsF2(models)
    // classes[j] lists the indices of the input models equivalent to reps[j].
    reps := [* *];
    classes := [* *];
    byKey := AssociativeArray(Integers());
    for i in [1..#models] do
        E, key := CanonicalFibrationF2(models[i]);
        found, j := IsDefined(byKey, key);
        if found then
            indices := classes[j];
            Append(~indices, i);
            classes[j] := indices;
        else
            Append(~reps, E);
            Append(~classes, [i]);
            byKey[key] := #reps;
        end if;
    end for;
    return reps, classes;
end function;

// Staged enumeration: orbits of the prefix (a1,a3), then orbits of its
// stabilizer on the affine quotient of (a4,a6) by the sextic shifts. This
// follows Lemma 1.1 and Algorithm 1.2 of Auel-Kulkarni-Petok-Weinbaum,
// "A census of cubic fourfolds over F_2", arXiv:2306.09908.

function EF2BaseChange(f, M)
    R := Parent(f);
    s := R.1;
    t := R.2;
    return Evaluate(f,
        [ M[1,1]*s+M[1,2]*t, M[2,1]*s+M[2,2]*t ]);
end function;

function EF2FormFromMask(R, d, n)
    B := EF2BinaryBasis(R, d);
    return &+[ R!((n div 2^(d-i)) mod 2)*B[i+1] : i in [0..d] ];
end function;

function EF2VectorFromMask(V, n)
    d := Dimension(V);
    return V![ (n div 2^(d-i)) mod 2 : i in [1..d] ];
end function;

function FibrationFromKeyF2(R, key)
    // Accepts keys from CanonicalFibrationF2 and from the staged enumeration.
    R := EF2CheckModel([ R | 0,0,0,0,0 ]);
    assert Type(key) eq RngIntElt;
    assert key ge 0 and key lt 2^32;
    a1 := EF2FormFromMask(R, 2, key div 2^29);
    a3 := EF2FormFromMask(R, 6, (key div 2^22) mod 2^7);
    a4 := EF2FormFromMask(R, 8, (key div 2^13) mod 2^9);
    a6 := EF2FormFromMask(R, 12, key mod 2^13);
    return [ R | a1,0,a3,a4,a6 ];
end function;

function FibrationOrbitPrefixesF2(R)
    // One <a1,a3> per orbit of the prefix action. The four a1 are the
    // GL(2,F2) classes of quadratics: zero, square, split, irreducible.
    R := EF2CheckModel([ R | 0,0,0,0,0 ]);
    s := R.1;
    t := R.2;
    k := BaseRing(R);
    matrices := [ Matrix(g) : g in GL(2, k) ];
    quadratics := [ EF2FormFromMask(R, 2, n) : n in [0..7] ];
    prefixes := [* *];
    for a1 in [ R | 0, s^2, s*t, s^2+s*t+t^2 ] do
        H := [ M : M in matrices | EF2BaseChange(a1, M) eq a1 ];
        seen := [ false : i in [1..128] ];
        for n in [0..127] do
            if seen[n+1] then continue; end if;
            a3 := EF2FormFromMask(R, 6, n);
            Append(~prefixes, <a1,a3>);
            for M in H do
                b3 := EF2BaseChange(a3, M);
                for p in quadratics do
                    c3 := b3+a1*(p^2+a1*p);
                    j := EF2BitKey(EF2Coefficients(c3, 6));
                    seen[j+1] := true;
                end for;
            end for;
        end for;
    end for;
    return prefixes;
end function;

function FibrationFiberOrbitKeysF2(prefix : GenericSmoothOnly := true)
    // Keys of one representative per orbit over the given <a1,a3>, and
    // <quotient dimension, distinct affine maps, orbits, retained orbits>.
    // GenericSmoothOnly drops exactly the models with Delta = 0.
    assert #prefix eq 2;
    a1 := prefix[1];
    a3 := prefix[2];
    R := Parent(a1);
    R := EF2CheckModel([ R | a1,0,a3,0,0 ]);
    k := BaseRing(R);
    V22 := VectorSpace(k, 22);
    B6 := EF2BinaryBasis(R, 6);
    S := sub< V22 |
        [ V22!EF2PairVector(a1*r, r^2+a3*r) : r in B6 ] >;
    C := sub< V22 |
        [ ReduceVector(S, V22.i) : i in [1..22] ] >;
    CB := Basis(C);
    d := Dimension(C);
    assert d eq 22-Dimension(S);
    assert d le 16;
    V := VectorSpace(k, d);

    // Modulo the sextic shifts, the stabilizing (M,p) act on C as a finite
    // group of affine maps. Applying every map once to a point gives its orbit.
    maps := [* *];
    quadratics := [ EF2FormFromMask(R, 2, n) : n in [0..7] ];
    for g in GL(2, k) do
        M := Matrix(g);
        if EF2BaseChange(a1, M) ne a1 then continue; end if;
        b3 := EF2BaseChange(a3, M);
        for p in quadratics do
            q := p^2+a1*p;
            if b3+a1*q ne a3 then continue; end if;
            origin := ReduceVector(S,
                V22!EF2PairVector(q^2+p*a3, q^3));
            shift := V!Coordinates(C, origin);
            rows := [];
            for z in CB do
                c, e := EF2PairFromVector(R, z);
                c := EF2BaseChange(c, M);
                e := EF2BaseChange(e, M);
                w := ReduceVector(S, V22!EF2PairVector(c, e+q*c));
                Append(~rows, Coordinates(C, w));
            end for;
            linear := Matrix(k, d, d, &cat rows);
            duplicate := false;
            for T in maps do
                if T[1] eq linear and T[2] eq shift then
                    duplicate := true;
                    break;
                end if;
            end for;
            if not duplicate then Append(~maps, <linear,shift>); end if;
        end for;
    end for;
    assert #maps gt 0;

    seen := [ false : i in [1..2^d] ];
    prefixBits := EF2Coefficients(a1, 2) cat EF2Coefficients(a3, 6);
    keys := [ Integers() | ];
    total := 0;
    da6 := a1^6;
    da4 := a1^5*a3;
    da4sq := a1^4;
    dconstant := a1^3*a3^3+a3^4;
    for n in [0..2^d-1] do
        if seen[n+1] then continue; end if;
        v := EF2VectorFromMask(V, n);
        for T in maps do
            w := v*T[1]+T[2];
            j := EF2BitKey(Eltseq(w));
            seen[j+1] := true;
        end for;
        total +:= 1;
        z := V22!(&+[ v[i]*CB[i] : i in [1..d] ]);
        keep := true;
        if GenericSmoothOnly then
            if a1 eq 0 then
                keep := a3 ne 0;
            else
                c, e := EF2PairFromVector(R, z);
                keep := da6*e+da4*c+da4sq*c^2+dconstant ne 0;
            end if;
        end if;
        if keep then Append(~keys, EF2BitKey(prefixBits cat Eltseq(z))); end if;
    end for;
    return keys, <d,#maps,total,#keys>;
end function;

function EnumerateFibrationOrbitKeysF2(R :
    GenericSmoothOnly := true, PrefixIndices := [])
    // An empty PrefixIndices runs all prefixes. Each stats entry is the prefix
    // index followed by the info from FibrationFiberOrbitKeysF2. The keys
    // encode representatives but need not equal CanonicalFibrationF2 keys.
    prefixes := FibrationOrbitPrefixesF2(R);
    indices := PrefixIndices;
    if #indices eq 0 then indices := [1..#prefixes]; end if;
    assert #SequenceToSet(indices) eq #indices;
    for i in indices do assert i in [1..#prefixes]; end for;
    allkeys := [ Integers() | ];
    stats := [* *];
    for i in indices do
        keys, info := FibrationFiberOrbitKeysF2(prefixes[i] :
            GenericSmoothOnly := GenericSmoothOnly);
        allkeys cat:= keys;
        Append(~stats, <i,info[1],info[2],info[3],info[4]>);
    end for;
    return allkeys, stats;
end function;
