#!/usr/bin/env sh
set -eu
if [ "${TRACE-0}" = "1" ]; then set -x; fi

# Variables
FUNCS_LIB=$(realpath "$1")          # Path to the Function Library
. "$FUNCS_LIB"                      # Import Func Libraries
# XYZ=$(realpath "$2")            # Location of the .xyz S0 xtb opt file

# Header
NT=${NT:-40}                        # Number Threads
NP=${NP:-2}                         # Number of processors
NROOTS=${NROOTS:-5}                 # Orca NRoots
TDDFT_MAXDIM=${TDDFT_MAXDIM:-10}    # Orca MaxDim

CDIR=$(realpath "$(pwd)")           # Current Dir
WORKDIR="$SLURM_TMPDIR/$SLURM_JOB_ID" # Working Dir
XYZ="$CDIR/orca_freq.xyz"
checkFile "$XYZ"                # Check if the file exists

# Graceful Exit
cleanup() { cp -anr "$WORKDIR/." "$CDIR"; }
trap cleanup EXIT

# Body
base_dir="$WORKDIR"
mkdir -p "$base_dir"
tda_nac_dir="${base_dir}/S0_tdaNAC"
tda_soc_dir="${base_dir}/S0_tdaSOC"
calc_name='S0_TDA_SP'
(loadGnuParallelEnv &&
 export FUNCS_LIB XYZ NT NP calc_name NROOTS TDDFT_MAXDIM
 parallel -j 2 \
     '. "$FUNCS_LIB" \
     && mkdir -p {1} \
     && cp -n "$XYZ" {1} \
     && cd {1} || exit 1 \
     && buildOrcaInputBase "$TPL_ORCA_TDA_BODY" "$TPL_ORCA_TDA_HEADER" \
        "$TDDFT_MAXDIM" "$((NT/2))" "$XYZ" \
        | addOrcaInputRoots $NROOTS 1 \
        | {2} \
        > ${calc_name}.inp \
     && loadOrcaEnv $((NT/2)) $((NP/2)) \
     && runOrca ${calc_name}.inp 1> /dev/null \
     && checkNEcho ${calc_name}.out \
     && orca2mkl ${calc_name}' \
 ::: "$tda_nac_dir"  "$tda_soc_dir" \
 :::+ 'addOrcaInputS0Nac' 'addOrcaInputS0Soc'
)



# Output
echo "$base_dir"
