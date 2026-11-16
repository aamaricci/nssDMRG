#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
dmrg_exe=${DMRG_EXE:-dmrg_kondo_mobile_1d}
ed_exe=${ED_EXE:-1d_Kondo_chain}
mpiexec=${MPIEXEC:-mpirun}
nproc=${NPROC:-2}
work_root=${WORK_ROOT:-${TMPDIR:-/tmp}/kondo-mobile-ed-dmrg-reference}
energy_tolerance=${ENERGY_TOLERANCE:-1e-10}
density_tolerance=${DENSITY_TOLERANCE:-1e-10}
corr_tolerance=${CORR_TOLERANCE:-1e-9}
nsite=6
norb=2

require_executable() {
    local executable=$1

    if [[ $executable == */* ]]; then
        [[ -x $executable ]] || {
            echo "compare_ed_dmrg: executable not found:" \
                "$executable" >&2
            exit 2
        }
    else
        command -v "$executable" >/dev/null 2>&1 || {
            echo "compare_ed_dmrg: executable not in PATH:" \
                "$executable" >&2
            exit 2
        }
    fi
}

require_executable "$dmrg_exe"
require_executable "$ed_exe"
if ((nproc > 1)); then
    require_executable "$mpiexec"
fi

run_solver() {
    local executable=$1
    local input=$2

    if ((nproc == 1)); then
        "$executable" FINPUT="$input"
    else
        "$mpiexec" -np "$nproc" "$executable" FINPUT="$input"
    fi
}

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

echo "Running the mobile DMRG reference in $dmrg_dir"
(
    cd "$dmrg_dir"
    run_solver "$dmrg_exe" DMRG.conf > dmrg.log
)

echo "Running the mobile ED reference in $ed_dir"
(
    cd "$ed_dir"
    run_solver "$ed_exe" inputED.conf > ed.log
)

energy_file=$(find "$dmrg_dir" -maxdepth 1 \
    -name 'energyVSblock.length*' -type f -print -quit)
[[ -n $energy_file ]] || {
    echo "compare_ed_dmrg: DMRG energy file not found" >&2
    exit 3
}
for file in egs.ed dens.ed Sz_corr.ed; do
    [[ -s $ed_dir/$file ]] || {
        echo "compare_ed_dmrg: ED output not found: $file" >&2
        exit 3
    }
done
for file in NimpVSsite.dat SimpSzSeleSzVSsites.dat; do
    [[ -s $dmrg_dir/$file ]] || {
        echo "compare_ed_dmrg: DMRG output not found: $file" >&2
        exit 3
    }
done

# DMRG writes E/(L*Norb), whereas ED writes the total energy.
dmrg_energy=$(awk -v nsite="$nsite" -v norb="$norb" \
    'NF >= 2 {value=$2} END {printf "%.17g",value*nsite*norb}' \
    "$energy_file")
ed_energy=$(awk 'NF >= 1 {value=$1} END {print value}' \
    "$ed_dir/egs.ed")
energy_error=$(awk -v x="$dmrg_energy" -v y="$ed_energy" \
    'BEGIN {d=x-y; if(d<0)d=-d; printf "%.16e",d}')

# ED stores bath positions first and impurity positions second.
awk -v nsite="$nsite" '
    NF >= 2*nsite {
        for(i=1;i<=nsite;i++) value[i]=$(nsite+i)
    }
    END {for(i=1;i<=nsite;i++) print i,value[i]}
' "$ed_dir/dens.ed" > "$work_root/ed_nimp.dat"
awk 'NF >= 2 {print $1,$2}' \
    "$dmrg_dir/NimpVSsite.dat" > "$work_root/dmrg_nimp.dat"

density_error=$(awk '
    NR==FNR {ed[$1]=$2; next}
    {
        d=$2-ed[$1]
        if(d<0)d=-d
        if(d>maximum)maximum=d
    }
    END {printf "%.16e",maximum}
' "$work_root/ed_nimp.dat" "$work_root/dmrg_nimp.dat")

# Sz_corr.ed is the lower triangle of the complete bath+impurity matrix.
awk -v nsite="$nsite" '
    NR>nsite && NR<=2*nsite {
        i=NR-nsite
        for(j=1;j<=nsite;j++) print i,j,$j
    }
' "$ed_dir/Sz_corr.ed" > "$work_root/ed_czz.dat"
awk 'NF >= 3 {print $1,$2,$3}' \
    "$dmrg_dir/SimpSzSeleSzVSsites.dat" \
    > "$work_root/dmrg_czz.dat"

corr_error=$(awk '
    NR==FNR {ed[$1 FS $2]=$3; next}
    {
        d=$3-ed[$1 FS $2]
        if(d<0)d=-d
        if(d>maximum)maximum=d
    }
    END {printf "%.16e",maximum}
' "$work_root/ed_czz.dat" "$work_root/dmrg_czz.dat")

printf '\n%-18s %22s %22s %14s\n' \
    quantity DMRG ED abs_error
printf '%-18s %22.14f %22.14f %14.6e\n' \
    E0 "$dmrg_energy" "$ed_energy" "$energy_error"
printf '%-18s %22s %22s %14.6e\n' \
    max_Nimp_error - - "$density_error"
printf '%-18s %22s %22s %14.6e\n' \
    max_Czz_error - - "$corr_error"

check_tolerance() {
    local label=$1
    local error=$2
    local tolerance=$3

    awk -v error="$error" -v tolerance="$tolerance" \
        'BEGIN {exit !(error <= tolerance)}' || {
        echo "FAIL: $label error $error exceeds $tolerance" >&2
        exit 1
    }
}

check_tolerance energy "$energy_error" "$energy_tolerance"
check_tolerance density "$density_error" "$density_tolerance"
check_tolerance correlation "$corr_error" "$corr_tolerance"

echo
echo "PASS: mobile ED/DMRG reference is within all tolerances."
echo "Artifacts: $work_root"
