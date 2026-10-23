load "functions_census.m" ;
load "MakeElliptic.m" ;

curves := MakeElliptic("mini_census.txt") ;
i := 900 ;

ZetaFunction(curves[i]) ;


exit ;
