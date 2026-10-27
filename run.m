// =====================================================================
// run.m  --  processes examples start..stop and writes a CSV file
//
// Launched by launch.sh as:
//     magma -b start:=1 stop:=375000 run.m
// For a quick test by hand:
//     magma -b start:=1 stop:=100 run.m
// =====================================================================

// ---------------------------------------------------------------------
// 1. Range of examples for this session (passed on the command line as
//    strings; defaults allow running the script without arguments)
// ---------------------------------------------------------------------
if not assigned start then start := "1";   end if;
if not assigned stop  then stop  := "100"; end if;
start := StringToInteger(start);
stop  := StringToInteger(stop);

outfile := "results/results_" cat Sprint(start) cat "_" cat Sprint(stop) cat ".csv";
batch   := 5000;     // rows held in memory before each write

// ---------------------------------------------------------------------
// 2. Shared structures -- create ONCE, outside the loop
// ---------------------------------------------------------------------
R<x,y,z> := PolynomialRing(Rationals(), 3);     // adjust to your ring

// ---------------------------------------------------------------------
// 3. Your functions (replace the placeholder bodies)
// ---------------------------------------------------------------------
load "functions_census.m" ;
load "point_counts.m";
GetExample := function(i)        // returns the i-th example
    return i;
end function;

Func1 := function(ex)            // builds the polynomial P2
    return x^2 + ex*y + z;
end function;

Func2 := function(P2)
    return TotalDegree(P2);
end function;

Func3 := function(P2)
    return #Terms(P2);
end function;

Func4 := function(P2)
    return IsIrreducible(P2);
end function;

// ---------------------------------------------------------------------
// 4. Helpers
// ---------------------------------------------------------------------
// Wrap a field in double quotes so commas inside it don't split columns
Q := func< s | "\"" cat s cat "\"" >;

// ---------------------------------------------------------------------
// 5. Open the output file (append mode; header only if file is new)
// ---------------------------------------------------------------------
isOld, _ := OpenTest(outfile, "r");
F := Open(outfile, "a");
if not isOld then
    Puts(F, "id,P2,result2,result3,result4");
    Flush(F);
end if;

// ---------------------------------------------------------------------
// 6. Main loop
// ---------------------------------------------------------------------
buf := [];
t0  := Cputime();

for i in [start..stop] do

    // Step 1: compute P2. If this fails, the other steps are skipped.
    ok := true;
    try
        ex := GetExample(i);
        P2 := Func1(ex);
        p2str := Q(Sprint(P2));
    catch e
        ok := false;
    end try;

    if ok then
        // Steps 2-4: each uses P2; one failure doesn't block the others
        try r2 := Sprint(Func2(P2)); catch e r2 := "ERROR"; end try;
        try r3 := Sprint(Func3(P2)); catch e r3 := "ERROR"; end try;
        try r4 := Sprint(Func4(P2)); catch e r4 := "ERROR"; end try;
        Append(~buf, Sprintf("%o,%o,%o,%o,%o", i, p2str, Q(r2), Q(r3), Q(r4)));
    else
        Append(~buf, Sprintf("%o,ERROR,NA,NA,NA", i));
    end if;

    // Write a batch to disk and report progress in the log
    if #buf ge batch then
        for s in buf do Puts(F, s); end for;
        Flush(F);
        buf := [];
        printf "done up to %o of %o  (%o s)\n", i, stop, Cputime(t0);
    end if;
end for;

// Write whatever is left
for s in buf do Puts(F, s); end for;
Flush(F);
delete F;                          // closes the file

printf "FINISHED %o..%o in %o s\n", start, stop, Cputime(t0);
quit;
