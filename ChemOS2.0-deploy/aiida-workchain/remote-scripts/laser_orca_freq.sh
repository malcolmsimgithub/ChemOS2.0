#!/usr/bin/env sh
set -eu
if [ "${TRACE-0}" = "1" ]; then set -x; fi


# Variables
FUNCS_LIB=$(realpath "$1")          # Path to the Function Library
. "$FUNCS_LIB"                      # Import Func Libraries
# XYZ=$(realpath "$2")                # Location of the *.xyz file from xtb

# Header
NT=${NT:-40}                        # Number Threads
NP=${NP:-2}                         # Number of processors
NROOTS=${NROOTS:-5}                 # Orca NRoots
DIRNAME='freq'                      # Directory name

CDIR=$(realpath "$(pwd)")           # Current Dir
WORKDIR="$SLURM_TMPDIR/$SLURM_JOB_ID" # Working Dir
XYZ="$CDIR/S0_XTB.xyz"
checkFile "$XYZ"                    # Check if the file exists

# Graceful Exit
cleanup() { cp -anr "$WORKDIR/." "$CDIR"; }
trap cleanup EXIT

# Body
base_dir="${WORKDIR}/${DIRNAME}"
mkdir -p "$base_dir" && cp "$XYZ" "$base_dir"
cd "$base_dir" || exit 1

input_name=$(realpath "${base_dir}/FREQ_OPT")

buildOrcaInputFreq "$TPL_ORCA_FREQ_BODY" "$TPL_ORCA_TDA_HEADER" \
    "$NT" "$(basename "$XYZ")" \
    > "${input_name}.inp"
loadOrcaEnv "$NT" "$NP" && runOrca "${input_name}.inp" 1> /dev/null

# Output
echo "$base_dir"
