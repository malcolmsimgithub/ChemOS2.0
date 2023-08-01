#!/usr/bin/env sh


export OMP_NUM_THREADS=1
export OMP_STACKSIZE=200M
export KMP_STACKSIZE=200M
ulimit -s unlimited

#Multiwfn extract excited state information

natoms=$1

# python3 xtbhess2np.py $natoms

nele=$(awk '/Number of Electrons/{print $NF }' S0_tdaNAC_SP.out)
echo $nele >> orca_mo.out
grep "ORBITAL ENERGIES" S0_tdaNAC_SP.out -A$((nele + 3)) >> orca_mo.out
grep "CARTESIAN GRADIENT" S0_tdaNAC_SP.out -A$((natoms+2)) >> grads.out


multiwfn S0_tdaNAC_SP.molden -nt $np < ./ES_info.txt | tee ./ES_info.out

grep "Sr index" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "D index" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "RMSD of hole in" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "RMSD of electron in" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "H index" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep " t index" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "Delta_r =" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "lambda =" ./ES_info.out |nl >> ES_results.out; echo >> ES_results.out
grep "Transition electric dipole moment between ground state" ./ES_info.out -A7 >> ES_results.out
grep "Transition electric dipole moment between excited states" ./ES_info.out -A17 >> ES_results.out

echo "Multiwfn done"
echo $(date +%T)

echo "all done"
