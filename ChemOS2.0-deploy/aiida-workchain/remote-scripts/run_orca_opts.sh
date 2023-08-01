#!/usr/bin/env sh
set -eu
set -x

# Header
getLast () { awk '{print $NF}' ; }
# FUNCS_LIB=$(dirname "$(realpath "$0")")'/lib/funcs.sh'
FUNCS_LIB="/home/a/aspuru/malcolms/laser/testing/lib/funcs.sh"
CDIR=$(realpath "$(pwd)")

. "$FUNCS_LIB"

# Load units
orca_opt="/home/a/aspuru/malcolms/laser/units/laser_orca_opt.sh"

# Expected Files
## From laser_xtb_crest.sh
S0_XTB_OPT_XYZ="${CDIR}/S0_XTB_OPT_XYZ.xyz"
T1_XTB_OPT_XYZ="${CDIR}/T1_XTB_OPT_XYZ.xyz"

## Orca Opt
# runOpt () { sbatch --dependency=afterok:"${JOB_ID_CREST}" \
#             "$orca_opt" "$FUNCS_LIB" "$1" "$2" "$3" \
#             | getLast ; }

runOpt () { "$orca_opt" "$FUNCS_LIB" "$1" "$2" "$3"; }
JOB_ID_OPT_S1=$(runOpt "$S0_XTB_OPT_XYZ" "S" "1")
JOB_ID_OPT_T1=$(runOpt "$T1_XTB_OPT_XYZ" "T" "1")
JOB_ID_OPT_T2=$(runOpt "$T1_XTB_OPT_XYZ" "T" "2")
JOB_ID_OPT_T3=$(runOpt "$T1_XTB_OPT_XYZ" "T" "3")
JOB_ID_OPT_T4=$(runOpt "$T1_XTB_OPT_XYZ" "T" "4")
