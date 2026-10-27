// run_census.m
//
// Writes the census of elliptic K3 surfaces over F2 to a CSV file, one row per surface, with
// the columns of CensusHeader (functions_census.m).
//
// The work can be split into chunks that run at the same time: the orbit keys are dealt out
// in turn, so chunk c of n gets keys number c, c+n, c+2n, ... .  Run from the repository root:
//
//     magma -b chunk:=1 chunks:=16 run_census.m < /dev/null
//
// writes census_part_1_of_16.csv and census_errors_1_of_16.csv.  Without chunk and chunks it
// does the whole census as chunk 1 of 1.  Adding max_rows:=100 stops after 100 rows, for a test.
//
// The output file is appended to: keys whose row is already in it are skipped, so a run that
// was interrupted continues where it stopped.  A surface for which CensusRow fails goes to the
// errors file, with its key and the error message, and the run carries on.
//
// Time: about 0.2 s per surface, almost all of it in LFunction, so about 80 hours of CPU for the
// whole census, divided by the number of chunks running at the same time.  Each chunk also
// spends about 100 s enumerating the keys.

load "elliptic_fibration_orbits_f2.m";
load "k3_minimality_f2.m";
load "functions_census.m";

// ---------------------------------------------------------------------------------------
// Settings, from the command line (as strings) or from the session (as integers).
// ---------------------------------------------------------------------------------------

AsInteger := func< x | Type(x) eq MonStgElt select StringToInteger(x) else x >;

if assigned chunks then n_chunks := AsInteger(chunks); else n_chunks := 1; end if;
if assigned chunk then this_chunk := AsInteger(chunk); else this_chunk := 1; end if;
if assigned max_rows then limit := AsInteger(max_rows); else limit := 0; end if;     // 0: no limit
error if this_chunk lt 1 or this_chunk gt n_chunks, "run_census: need 1 <= chunk <= chunks";

suffix := Sprint(this_chunk) cat "_of_" cat Sprint(n_chunks) cat ".csv";
out_file := "census_part_" cat suffix;
err_file := "census_errors_" cat suffix;

// ---------------------------------------------------------------------------------------
// Resuming: the keys already in the output file, and whether it needs a header.
// ---------------------------------------------------------------------------------------

function ReadDone(file, columns)
    // The keys whose rows are already complete in a census file.
    //   Input:  file,    the name of a census CSV file (it need not exist)
    //           columns, the number of columns of a complete row
    //   Output: three values: the set of keys of complete rows, whether the file has any
    //           content, and whether its last line is unfinished (no final newline)
    try
        text := Read(file);
    catch e
        return {}, false, false;
    end try;
    done := { StringToInteger(Split(line, ",")[1]) : line in Split(text, "\n")
              | #Split(line, ",") eq columns and Split(line, ",")[1] ne "key" };
    return done, #text gt 0, #text gt 0 and text[#text] ne "\n";
end function;

columns := #Split(CensusHeader(), ",");
done, has_content, unfinished := ReadDone(out_file, columns);
_, err_has_content, err_unfinished := ReadDone(err_file, 2);

F := Open(out_file, "a");
if not has_content then Puts(F, CensusHeader()); end if;
if unfinished then Puts(F, ""); end if;           // end the cut-off line before appending
Ferr := Open(err_file, "a");
if not err_has_content then Puts(Ferr, "key,error"); end if;
if err_unfinished then Puts(Ferr, ""); end if;

// An error message on one line and without commas, so that it fits in one CSV cell.
OneCell := func< x | Join(Split(Join(Split(Sprint(x), "\n"), " "), ","), ";") >;

// ---------------------------------------------------------------------------------------
// The census: the keys of this chunk, the everywhere minimal models, one row each.
// ---------------------------------------------------------------------------------------

printf "Chunk %o of %o: writing %o (%o rows already there)\n", this_chunk, n_chunks, out_file, #done;

R<s,t> := PolynomialRing(GF(2), 2);
prefixes := FibrationOrbitPrefixesF2(R);

index := 0;        // position of the key among all keys, for dealing them out to chunks
written := 0;
failed := 0;
start := Cputime();
for prefix in prefixes do
    for key in FibrationFiberOrbitKeysF2(prefix) do
        index +:= 1;
        if (index - this_chunk) mod n_chunks ne 0 or key in done then
            continue;
        end if;
        A := FibrationFromKeyF2(R, key);
        if not IsEverywhereMinimalF2(A) then
            continue;
        end if;
        try
            Puts(F, CensusRow(key, A));
            Flush(F);
            written +:= 1;
            if written mod 1000 eq 0 then
                printf "%o rows, %o errors, %o s\n", written, failed, Cputime(start);
            end if;
        catch err
            Puts(Ferr, Sprint(key) cat "," cat OneCell(err`Object));
            Flush(Ferr);
            failed +:= 1;
        end try;
        if limit gt 0 and written ge limit then
            break;
        end if;
    end for;
    if limit gt 0 and written ge limit then
        break;
    end if;
end for;

delete F;
delete Ferr;
printf "Done: %o rows written, %o errors, %o s.  Rows in %o, errors in %o\n",
       written, failed, Cputime(start), out_file, err_file;
