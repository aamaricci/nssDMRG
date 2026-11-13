#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
dmrg_exe=${DMRG_EXE:-dmrg_kondo_1d}
ed_exe=${ED_EXE:-1d_Kondo_chain}
mpiexec=${MPIEXEC:-mpirun}
nproc=${NPROC:-2}
work_root=${WORK_ROOT:-${TMPDIR:-/tmp}/kondo-ed-dmrg-reference}
energy_tolerance=${ENERGY_TOLERANCE:-1e-10}
corr_tolerance=${CORR_TOLERANCE:-1e-9}
nsite=6

require_executable() {
    local executable=$1

    if [[ $executable == */* ]]; then
        [[ -x $executable ]] || {
            echo "compare_ed_dmrg: executable not found: $executable" >&2
            exit 2
        }
    else
        command -v "$executable" >/dev/null 2>&1 || {
            echo "compare_ed_dmrg: executable not found in PATH: $executable" >&2
            exit 2
        }
    fi
}

require_executable "$dmrg_exe"
require_executable "$ed_exe"
require_executable "$mpiexec"

if [[ -z $work_root || $work_root == / ]]; then
    echo "compare_ed_dmrg: unsafe WORK_ROOT: $work_root" >&2
    exit 2
fi

dmrg_dir=$work_root/dmrg
ed_dir=$work_root/ed
rm -rf "$work_root"
mkdir -p "$dmrg_dir" "$ed_dir"
cp "$script_dir/DMRG.conf" "$dmrg_dir/DMRG.conf"
cp "$script_dir/inputED.conf" "$ed_dir/inputED.conf"

echo "Running the DMRG reference in $dmrg_dir"
(
    cd "$dmrg_dir"
    "$mpiexec" -np "$nproc" "$dmrg_exe" FINPUT=DMRG.conf > dmrg.log
)

echo "Running the ED reference in $ed_dir"
(
    cd "$ed_dir"
    "$mpiexec" -np "$nproc" "$ed_exe" FINPUT=inputED.conf > ed.log
)

energy_file=$(find "$dmrg_dir" -maxdepth 1 \
    -name 'energyVSblock.length*' -type f -print -quit)
[[ -n $energy_file ]] || {
    echo "compare_ed_dmrg: DMRG energy file not found" >&2
    exit 3
}
[[ -s $ed_dir/egs.ed ]] || {
    echo "compare_ed_dmrg: ED energy file not found" >&2
    exit 3
}
[[ -s $dmrg_dir/SimpSzSeSzVSsite.dat ]] || {
    echo "compare_ed_dmrg: DMRG correlation file not found" >&2
    exit 3
}
[[ -s $ed_dir/Sz_corr.ed ]] || {
    echo "compare_ed_dmrg: ED correlation file not found" >&2
    exit 3
}

dmrg_energy=$(awk 'NF >= 2 {value=$2} END {printf "%.17g",value*6}' \
    "$energy_file")
ed_energy=$(awk 'NF >= 1 {value=$1} END {print value}' "$ed_dir/egs.ed")
energy_error=$(awk -v x="$dmrg_energy" -v y="$ed_energy" \
    'BEGIN {d=x-y; if(d<0)d=-d; printf "%.16e",d}')

awk -v nsite="$nsite" '
    NF >= nsite {for(i=1;i<=nsite;i++) value[i]=$i}
    END {for(i=1;i<=nsite;i++) print i,value[i]}
' "$ed_dir/Sz_corr.ed" > "$work_root/ed_czz.dat"

awk 'NF >= 2 {print $1,$2}' \
    "$dmrg_dir/SimpSzSeSzVSsite.dat" > "$work_root/dmrg_czz.dat"

corr_error=$(awk '
    NR==FNR {ed[$1]=$2; next}
    {
        d=$2-ed[$1]
        if(d<0)d=-d
        if(d>maximum)maximum=d
    }
    END {printf "%.16e",maximum}
' "$work_root/ed_czz.dat" "$work_root/dmrg_czz.dat")

printf '\n%-12s %22s %22s %14s\n' quantity DMRG ED abs_error
printf '%-12s %22.14f %22.14f %14.6e\n' \
    E0 "$dmrg_energy" "$ed_energy" "$energy_error"
printf '%-12s %22s %22s %14.6e\n' \
    max_Czz_error - - "$corr_error"

printf '\n%-6s %22s %22s %14s\n' site DMRG_Czz ED_Czz abs_error
awk '
    NR==FNR {ed[$1]=$2; next}
    {
        d=$2-ed[$1]
        if(d<0)d=-d
        printf "%-6d %22.14e %22.14e %14.6e\n",$1,$2,ed[$1],d
    }
' "$work_root/ed_czz.dat" "$work_root/dmrg_czz.dat"

awk -v error="$energy_error" -v tolerance="$energy_tolerance" \
    'BEGIN {exit !(error <= tolerance)}' || {
    echo "FAIL: energy error $energy_error exceeds $energy_tolerance" >&2
    exit 1
}
awk -v error="$corr_error" -v tolerance="$corr_tolerance" \
    'BEGIN {exit !(error <= tolerance)}' || {
    echo "FAIL: correlation error $corr_error exceeds $corr_tolerance" >&2
    exit 1
}

echo
echo "PASS: ED/DMRG reference is within the configured tolerances."
echo "Artifacts: $work_root"
