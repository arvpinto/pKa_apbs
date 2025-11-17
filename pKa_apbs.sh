#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $(basename "$0") RESNAME RESNUMBER"
    echo "Ex.:   $(basename "$0") HIP 211"
    exit 1
fi

RESNAME="$1"
RESNUM="$2"

# initialize output files
: > ee_RES.dat
: > ee_RES_enz.dat
: > ee_no_RES.dat

for x in frame*.gro; do
    echo "Processing $x ..."

    # temporary working directories
    for d in RES_enz no_RES RES; do
        mkdir -p "$d"
    done

    # Convert input to PDB
    gmx editconf -f "$x" -o current.pdb >> apbs.log 2>&1

    # extract chain A
    pdb_chain -A current.pdb > current_CHAIN_A.pdb 

    # generate full enzyme PQR
    pdb2pqr --ff AMBER -k current_CHAIN_A.pdb RES_enz.pqr >> apbs.log 2>&1

    # generate "no_RES" version by zeroing the RESNAME RESNUM charges
    perl -pe '
        if (/'"$RESNAME"'\s+'"$RESNUM"'/) {
            s/^((?:\S+\s+){8})(\S+)/
              my $pad = length($2) - length("0.0000");
              $pad = 0 if $pad < 0;
              $1 . (" " x $pad) . "0.0000"
            /e;
        }
    ' RES_enz.pqr > no_RES.pqr

    # extract the single residue into RES.pqr
    grep "$RESNAME   $RESNUM" RES_enz.pqr > RES.pqr || {
        echo "Warning: Residue $RESNAME $RESNUM not found in $x"
    }
    echo -e "TER\nEND" >> RES.pqr

    # Move PQR files into their directories
    mv RES.pqr RES/
    mv RES_enz.pqr RES_enz/
    mv no_RES.pqr no_RES/

    # Prepare APBS input files
    sed "s/INPUT_NAME/RES.pqr/"      APBS_INPUT.in > RES/RES_apbs.in
    sed "s/INPUT_NAME/RES_enz.pqr/"  APBS_INPUT.in > RES_enz/RES_enz_apbs.in
    sed "s/INPUT_NAME/no_RES.pqr/"   APBS_INPUT.in > no_RES/no_RES_apbs.in

    # Run APBS (capture output)
    ( cd RES     && apbs RES_apbs.in     | tee RES_apbs.out ) >> apbs.log 2>&1
    ( cd RES_enz && apbs RES_enz_apbs.in | tee RES_enz_apbs.out ) >> apbs.log 2>&1
    ( cd no_RES  && apbs no_RES_apbs.in  | tee no_RES_apbs.out ) >> apbs.log 2>&1

    # Extract energies
    grep "Total electrostatic energy" RES_enz/RES_enz_apbs.out >> ee_RES_enz.dat
    grep "Total electrostatic energy" RES/RES_apbs.out        >> ee_RES.dat
    grep "Total electrostatic energy" no_RES/no_RES_apbs.out  >> ee_no_RES.dat

    # Cleanup safely
    rm -rf RES RES_enz no_RES
    rm -f current.pdb current_CHAIN_A.pdb

done

# --- Compute means, sample SDs and ΔE (explicitly using column 5) ---

get_stats_col5() {
    local file="$1"
    awk '
    # Only consider lines where field 5 looks like a number (supports scientific notation)
    ($5 ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
        x = $5 + 0
        sum += x
        sumsq += x*x
        n++
    }
    END {
        if (n > 0) {
            mean = sum / n
            if (n > 1) {
                sd = sqrt( (sumsq - n * mean * mean) / (n - 1) )
            } else {
                sd = 0
            }
            # print mean and sd in scientific notation
            printf "%.10e %.10e", mean, sd
        } else {
            # No valid numeric lines
            printf "nan nan"
        }
    }' "$file"
}

read avg_RES_enz sd_RES_enz <<< "$(get_stats_col5 ee_RES_enz.dat)"
read avg_RES     sd_RES     <<< "$(get_stats_col5 ee_RES.dat)"
read avg_no_RES  sd_no_RES  <<< "$(get_stats_col5 ee_no_RES.dat)"

# compute Δ = <RES_enz> - <no_RES> - <RES>
# handle nan values robustly in awk: result will be nan if any input is nan
delta=$(awk -v a="$avg_RES_enz" -v b="$avg_no_RES" -v c="$avg_RES" \
          'BEGIN {
               # check for non-numeric (nan) inputs
               if (a != a || b != b || c != c) { print "nan"; exit }
               printf "%.10e", (a - b - c)
           }')

echo ""
echo "-----------------------------------------------------------"
echo "RES_enz : mean = $avg_RES_enz   SD = $sd_RES_enz"
echo "RES     : mean = $avg_RES       SD = $sd_RES"
echo "no_RES  : mean = $avg_no_RES    SD = $sd_no_RES"
echo ""
echo "ΔE = <RES_enz> - <no_RES> - <RES> = $delta kJ/mol"
echo "-----------------------------------------------------------"
echo ""

echo "Done."

