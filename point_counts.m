load "functions_census.m" ;

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

function EqPointCounts(poly , n)
// assuming we have a K3 over F_2
//equ should already be initalized as an elliptic curve over a univariate
//rational function field over GF(2)
  Qt := PolynomialRing(Rationals());
  R<T> := PowerSeriesRing(Rationals(), n + 5);
  P0 := 1 - T ; P1 := 1 ; P1 := Qt ! P1 ;
  P2 := poly;
  P3 := 1 ; P3 := Qt ! P3 ; P4 := 1-4*T ;
   // convert a polynomial into an exact element of R
  ser := func< f | &+[ R | Coefficient(f, i) * T^i : i in [0..Degree(f)] ] >;
  Z := ser(P1) * ser(P3) / ( ser(P0) * ser(P2) * ser(P4) );
  L := T * Derivative(Z) / Z;
  return [ Integers() ! Coefficient(L, r) : r in [1..n] ];
end function ;
