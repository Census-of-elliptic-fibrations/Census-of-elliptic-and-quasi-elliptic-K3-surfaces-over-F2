// timing_tests.m
//
// Times built-in Magma functions on the models in mini_census.txt, to estimate how long
// they would take on the whole census.  Each line of mini_census.txt is a model
// [a1, a2, a3, a4, a6] (binary forms in s, t over F2, a2 = 0); each is turned into its
// generic fibre over F2(T), T = t/s.
//
// Run with  magma -b timing_tests.m  from the repository root.  Output goes to the screen
// and to timing_report.txt.

CENSUS_SIZE := 1397598;
input_file  := "mini_census.txt";
report_file := "timing_report.txt";

K<T> := RationalFunctionField(GF(2));

// ---------------------------------------------------------------------------------------
// Read the sample.  Each line is parsed in a polynomial ring made inside eval, then
// evaluated at (1, T) to land in F2(T).
// ---------------------------------------------------------------------------------------

lines := [ l : l in Split(Read(input_file), "\n") | #l gt 0 ];
start := Cputime();
src := "Q<s,t> := PolynomialRing(GF(2),2); return [ "
       cat Join([ "[Q|" cat Substring(l, 2, #l-1) cat "]" : l in lines ], ", ") cat " ];";
models := eval src;
parse_time := Cputime(start);

// Magma caches results (local data, L-function, ...) on the curve, so every test gets
// freshly built curves; otherwise later tests would reuse the work of earlier ones.
Curves := func< | [ EllipticCurve([ K | Evaluate(f, [K!1, T]) : f in A ]) : A in models ] >;

start := Cputime();
curves := Curves();
build_time := Cputime(start);

// ---------------------------------------------------------------------------------------
// The tests.  Each is <name, function of E, number of curves to run it on>.  The slow ones
// run on a prefix of the sample only; the estimate is scaled from what was run.
// ---------------------------------------------------------------------------------------

tests := [*
    < "Discriminant",          func< E | Discriminant(E) >,          1000 >,
    < "jInvariant",            func< E | jInvariant(E) >,            1000 >,
    < "LocalInformation",      func< E | LocalInformation(E) >,      1000 >,
    < "KodairaSymbols",        func< E | KodairaSymbols(E) >,        1000 >,
    < "Conductor",             func< E | Conductor(E) >,             1000 >,
    < "MinimalModel",          func< E | MinimalModel(E) >,          1000 >,
    < "TorsionSubgroup",       func< E | TorsionSubgroup(E) >,       1000 >,
    < "LFunction",             func< E | LFunction(E) >,             100 >,
    < "AnalyticInformation",   func< E | AnalyticInformation(E) >,   100 >
*];

function Pad(x, w)
    // x as a string, right-padded with spaces to width w; reals get 4 significant digits
    str := Type(x) eq FldReElt select Sprint(RealField(4)!x) else Sprint(x);
    while #str lt w do
        str cat:= " ";
    end while;
    return str;
end function;

procedure Row(fields)
    widths := [22, 6, 12, 12, 12, 16];
    print &cat[ Pad(fields[i], widths[i]) : i in [1..#fields] ];
end procedure;

SetOutputFile(report_file : Overwrite := true);

printf "Sample: %o models from %o\n", #curves, input_file;
printf "Parsing: %o s, building curves: %o s\n\n", Pad(parse_time, 0), Pad(build_time, 0);
// Cputime ticks in steps of 10 ms, so the total is timed over the whole loop; the
// per-curve max is only meaningful when it is well above 10 ms.
Row(["function", "n", "total (s)", "mean (ms)", "max (ms)", "full census (h)"]);

for test in tests do
    name, f, n := Explode(test);
    curves := Curves();
    n := Min(n, #curves);
    worst := 0.0;
    failed := false;
    start := Cputime();
    for i in [1..n] do
        try
            t0 := Cputime();
            _ := f(curves[i]);
            worst := Max(worst, Cputime(t0));
        catch err
            printf "%o failed on model %o: %o\n", name, i, err`Object;
            failed := true;
        end try;
        if failed then
            break;
        end if;
    end for;
    total := Cputime(start);
    if not failed then
        Row([* name, n, total, 1000*total/n, 1000*worst, total/n*CENSUS_SIZE/3600 *]);
    end if;
end for;

UnsetOutputFile();
