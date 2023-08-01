#!/usr/bin/env python
import numpy as np

dips = np.zeros((6,6,3))
excs = np.zeros(6)

with open('final_step/ES_results.out','r') as f_in:
    for i, line in enumerate(f_in):
        if i in range(50, 55):
            dips[0][i-49] = [float(x) for x in line.split()[2:5]]
            dips[i-49][0] = [float(x) for x in line.split()[2:5]]
            excs[i-49] = float(line.split()[5])
        if i in range(59, 63):
            dips[1][i-57] = [float(x) for x in line.split()[2:5]]
            dips[i-57][1] = [float(x) for x in line.split()[2:5]]
        if i in range(64, 67):
            dips[2][i-61] = [float(x) for x in line.split()[2:5]]
            dips[i-61][2] = [float(x) for x in line.split()[2:5]]
        if i in range(68, 70):
            dips[3][i-64] = [float(x) for x in line.split()[2:5]]
            dips[i-64][3] = [float(x) for x in line.split()[2:5]]
        if i in range(71, 72):
            dips[4][i-66] = [float(x) for x in line.split()[2:5]]
            dips[i-66][4] = [float(x) for x in line.split()[2:5]]

np.save('dipmoments.npy',dips[:4,:4,:3])
np.save('excenergies.npy', excs[:4])

