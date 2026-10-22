//Modified by Elise

load "elliptic_fibration_orbits_f2.m";
load "k3_minimality_f2.m";

R<s,t> := PolynomialRing(GF(2),2);
start := Cputime();
prefixes := FibrationOrbitPrefixesF2(R);
a1_reps := [R|0,s^2,s*t,s^2+s*t+t^2];

// Minimal fibrations are written one per line to census.txt, gzipped at the end
census_file := "census.txt";
F := Open(census_file, "w");

for i in [1..#prefixes] do
    K,I := FibrationFiberOrbitKeysF2(prefixes[i]);
    for key in K do
        A := FibrationFromKeyF2(R,key);
        ok,bad := IsEverywhereMinimalF2(A);
        if ok then
            // Sprint(A) spreads the list over several lines, so join the forms by hand
            Puts(F, "[" cat Join([Sprint(f) : f in A], ", ") cat "]");
        end if;
    end for;
    assert #K eq I[4];
    a1 := prefixes[i][1];
    group := Index(a1_reps,a1);
    assert group ne 0;
end for;

Flush(F);
delete F;
System("gzip -f " cat census_file);

quit;
