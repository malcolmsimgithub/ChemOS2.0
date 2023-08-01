#!/usr/bin/env sh
set -eu

#SBATCH --ntasks=40
#SBATCH --time=01:00:00

# if [ "${TRACE-0}" = "1" ]; then set -x; fi
set -x
# Variables
FUNCS_LIB="$1"         # Path to the Function Library
. "$FUNCS_LIB"                      # Import Func Libraries

# Parameters
STRUCTFILE="$2"
NT=${NT:-40}                        # Number Threads
NP=${NP:-2}                         # Number of processors
CHARGE=0
SPIN=1
SPIN_T=3

# Current DIR
CDIR=$(realpath "$(pwd)")

# Environment
MARVINDIR="/home/a/aspuru/$USER/bin/marvin/bin"
XTBDIR="/home/a/aspuru/$USER/bin/xtb"
WORKDIR="$SLURM_TMPDIR/$SLURM_JOB_ID"

# Graceful Exit
cleanup() { cp -anr "$WORKDIR/." "$CDIR"; }
trap cleanup EXIT


mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Get initial xyz
mkdir -p "./mol"

# Perform opt and crest
mkdir -p "./crest"

cd "./crest" || exit 1
xyz_opt=$(loadXtbEnv "$XTBDIR" "$NT" "$NP" \
    && runXtbOpt "$CDIR/structfile.xyz" "$CHARGE" "$SPIN")
xyz_crest=$(loadXtbEnv "$XTBDIR" "$NT" "$NP" \
    && runXtbCrest "$xyz_opt" "$CHARGE" "$SPIN")
cd "$WORKDIR"

# Compute S0 and T1
# mkdir -p "./opt/S0_XTB" "./opt/T1_XTB"
# cd "./opt" || exit 1
mkdir -p "./S0_XTB" "./T1_XTB"
xyz_best=$(getXtbCrestBest "$xyz_crest" "./crest_best.xyz")
(loadXtbEnv "$XTBDIR" "$NT" "$NP" &&
 loadGnuParallelEnv &&
 parallel -j 2 -N2 \
     "source $FUNCS_LIB \
     && loadXtbEnv $XTBDIR $((NT / 2)) $((NP / 2)) \
     && cd {1} \
     && runXtbOhess $xyz_best $CHARGE {2}" \
     ::: "S0_XTB" "$SPIN" "T1_XTB" "$SPIN_T"
)
checkNEcho "./S0_XTB/xtbopt.xyz"
checkNEcho "./T1_XTB/xtbopt.xyz"
