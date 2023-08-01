#!/usr/bin/env python
import numpy as np
import sys

natoms = int(sys.argv[1])
hess_arr = np.empty(pow(natoms*3, 2))

with open('xtb_results/S0_hessian', 'r') as f_in:
    count = 0
    for i, line in enumerate(f_in):
        if i == 0:
            continue
        vibs = [float(x) for x in line.split()]
        for vib in vibs:
            hess_arr[count] = vib
            count += 1
#hess_arr = hess_arr.reshape((natoms*3, natoms*3))
#print(hess_arr)
np.save('hessian.npy', hess_arr.reshape((natoms*3, natoms*3)))



