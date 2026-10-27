// run_mini_census.m
//
// The census CSV for the 1000 models of mini_census.m: the same columns as run_census.m
// (CensusHeader in functions_census.m), written to mini_census_rows.csv, with failures in
// mini_census_errors.csv.  Each model gets the key it has in the full census, so its row
// can be compared with the full census later.
//
// Run from the repository root:
//     magma -b run_mini_census.m < /dev/null
// About 4 minutes.  For the 10000 models of small_census10k.m, change the line
// load "mini_census.m" below to load "small_census10k.m" (about 40 minutes).

load "elliptic_fibration_orbits_f2.m";
load "functions_census.m";
load "mini_census.m";

out_file := "mini_census_rows.csv";
err_file := "mini_census_errors.csv";

function KeyOfModel(A)
    // The orbit key of a census model: the coefficient bits of a1, a3, a4, a6, most
    // significant first, each in the order s^d, s^(d-1)*t, ..., t^d, as FibrationFromKeyF2
    // reads them.
    //   Input:  A, a model [a1, a2, a3, a4, a6] with a2 = 0, binary forms in s, t over F_2
    //   Output: the key, an integer in [0, 2^32)
    return EF2BitKey(EF2Coefficients(A[1], 2) cat EF2Coefficients(A[3], 6)
                     cat EF2Coefficients(A[4], 8) cat EF2Coefficients(A[5], 12));
end function;

// An error message on one line and without commas, so that it fits in one CSV cell.
OneCell := func< x | Join(Split(Join(Split(Sprint(x), "\n"), " "), ","), ";") >;

F := Open(out_file, "w");
Puts(F, CensusHeader());
Ferr := Open(err_file, "w");
Puts(Ferr, "index,key,error");

failed := 0;
start := Cputime();
for i in [1..#models] do
    A := models[i];
    key := KeyOfModel(A);
    try
        error if FibrationFromKeyF2(Parent(A[1]), key) ne A, "the key does not give back the model";
        Puts(F, CensusRow(key, A));
        Flush(F);
    catch err
        Puts(Ferr, Sprint(i) cat "," cat Sprint(key) cat "," cat OneCell(err`Object));
        Flush(Ferr);
        failed +:= 1;
    end try;
    if i mod 100 eq 0 then
        printf "%o / %o models, %o errors, %o s\n", i, #models, failed, Cputime(start);
    end if;
end for;

delete F;
delete Ferr;
printf "Done: %o rows in %o, %o errors in %o, %o s\n",
       #models - failed, out_file, failed, err_file, Cputime(start);
