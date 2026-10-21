// Input: [a1,a2,a3,a4,a6] over F2[s,t], homogeneous of degrees 2i.
K3F2_U<u> := PolynomialRing(GF(2));
K3F2_primes := [u, u+1, u^2+u+1];
K3F2_powers := [[p^j : j in [1..12]] : p in K3F2_primes];
K3F2_q := [[K3F2_U | &+[K3F2_U | GF(2)!b[j+1]*u^j :
    j in [0..d-1]] : b in CartesianPower([0,1],d)] : d in [1,2]];
K3F2_v := [[K3F2_U | &+[K3F2_U | GF(2)!b[j+1]*u^j :
    j in [0..3*d-1]] : b in CartesianPower([0,1],3*d)] : d in [1,2]];

function K3F2_NonMinimal(B, k)
    a1,a2,a3,a4,a6 := Explode(B);
    pp := K3F2_powers[k];
    if a1 mod pp[1] ne 0 then return false; end if;
    d := Degree(pp[1]);
    // T(r,q,v)T(h,k,w) = T(r+h,q+k,v+q*h+w). The subgroup with
    // h in pi^2, k in pi, w in pi^3 preserves v(a_i) >= i, as does
    // its inverse. Right composition reduces r,q, then v+q*h modulo
    // pi^3. Thus these finite representatives meet every relevant coset.
    for q in K3F2_q[d] do
        r := (a2+q*a1+q^2) mod pp[2];
        if (a3+r*a1) mod pp[3] ne 0 then continue; end if;
        c4 := a4+q*a3+r*q*a1+r^2;
        c6 := a6+r*a4+r^2*a2+r^3;
        c3 := a3+r*a1;
        for v in K3F2_v[d] do
            if (c4+v*a1) mod pp[4] eq 0 and
               (c6+v*c3+v^2) mod pp[6] eq 0 then
                return true;
            end if;
        end for;
    end for;
    return false;
end function;

function IsEverywhereMinimalF2(A)
    R := Universe(A);
    s := R.1; t := R.2;
    U := K3F2_U; u := U.1;
    B := [U | Evaluate(a,[u,U!1]) : a in A];
    a1,a2,a3,a4,a6 := Explode(B);
    a1sq := a1^2; a3sq := a3^2;
    b8 := a1sq*a6+a1*a3*a4+a2*a3sq+a4^2;
    delta := a1sq^2*b8+a3sq^2+a1sq*a1*a3sq*a3;
    bad := [R | ];
    if delta eq 0 then return false, bad; end if;
    // delta = Delta(u,1); homogeneity recovers v_t(Delta) from its degree.
    if delta mod K3F2_powers[1][12] eq 0 and
       K3F2_NonMinimal(B,1) then Append(~bad,s); end if;
    if Degree(delta) le 12 then
        C := [U | Evaluate(a,[U!1,u]) : a in A];
        if K3F2_NonMinimal(C,1) then Append(~bad,t); end if;
    end if;
    if delta mod K3F2_powers[2][12] eq 0 and
       K3F2_NonMinimal(B,2) then Append(~bad,s+t); end if;
    if delta mod K3F2_powers[3][12] eq 0 and
       K3F2_NonMinimal(B,3) then Append(~bad,s^2+s*t+t^2); end if;
    return IsEmpty(bad), bad;
end function;
