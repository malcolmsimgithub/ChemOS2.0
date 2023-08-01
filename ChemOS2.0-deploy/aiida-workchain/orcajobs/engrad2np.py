import numpy as np
import sys

natoms = int(sys.argv[1])

grads = np.empty((4, natoms, 3))

with open('grads.out', 'r') as f_in:
    for i, line in enumerate(f_in):
        if 3 <= i < (natoms + 3):
            grads[1][i-3] = [float(x) for x in line.split()[-3:]]
        if (natoms + 7) <= i < (natoms * 2 + 7):
            grads[2][i - natoms -7] = [float(x) for x in line.split()[-3:]]
        if (natoms * 2 + 11) <= i < (natoms * 3 + 11):
            grads[3][i - natoms * 2 - 11] = [float(x) for x in line.split()[-3:]]
        if (natoms * 3 + 15) <= i < (natoms * 4 + 15):
            grads[0][i - natoms * 3 - 15] = [float(x) for x in line.split()[-3:]]

np.save('grads.npy', grads)
