#!/usr/bin/env sh
set -eu
set -x

# Variables
FUNCS_LIB=$(realpath "$1"); shift   # Path to the Function Library
. "$FUNCS_LIB"                      # Import Func Libraries
XYZ_ARR="$*"                        # Path to XYZ files

# Header
NT=${NT:-40}                        # Number Threads
NP=${NP:-2}                         # Number of processors
NROOTS=${NROOTS:-5}                 # Orca NRoots
TDDFT_MAXDIM=${TDDFT_MAXDIM:-10}    # Orca MaxDim
DIRNAME='sp_comb'                   # Directory name

CDIR=$(realpath "$(pwd)")           # Current Dir
WORKDIR="$SLURM_TMPDIR/$SLURM_JOB_ID"  # Working Dir

# Graceful Exit
cleanup() { cp -anr "$WORKDIR/." "$CDIR"; }
trap cleanup EXIT

# Check .xyz Files
for opt_dir in $XYZ_ARR; do
    checkFile "$(realpath "$opt_dir")"
done

# Base Names of the OPT
# getNames () { for item in $1; do base=$(basename "$item"); echo "${base%.*}"; done }
getNames () { for item in $1; do base=$(basename "$item"); echo "${base%.*}"; done }
names=$(getNames "$XYZ_ARR")

# Body
base_dir="${WORKDIR}/$DIRNAME"
mkdir -p "$base_dir"
calc_name='TDA_SP_COMB'
(loadGnuParallelEnv &&
 export FUNCS_LIB TDDFT_MAXDIM NT NP NROOTS base_dir calc_name CDIR
parallel -j 5 \
     '. "$FUNCS_LIB" \
     && wdir="${CDIR}"/{1}_{3} \
     && mkdir -p "$wdir" \
     && cp -n {2} "$wdir"/orcaopt.xyz \
     && cd "$wdir" || exit 1 \
     && buildOrcaInputBase "$TPL_ORCA_TDA_BODY" "$TPL_ORCA_TDA_HEADER" \
        "$TDDFT_MAXDIM" $((NT/5)) orcaopt.xyz \
        | addOrcaInputComb "$NROOTS" {3} \
        > "${calc_name}".inp \
     && loadOrcaEnv $((NT/5)) $((NP/5)) \
     && runOrca "${calc_name}".inp 1> /dev/null \
     && checkNEcho "${calc_name}".out' \
 ::: $names \
 :::+ $XYZ_ARR \
 ::: 'Singlet' 'Triplet'
)

# Output
echo "$base_dir"
