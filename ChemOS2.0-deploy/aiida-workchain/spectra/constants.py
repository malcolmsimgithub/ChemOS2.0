#!/usr/bin/env python
from numpy import pi
# Physical constants used in computations
hbar = 0.6582119514            # fs eV
ev_nm = 1239.842
ev_wn = 8065.54
hartree_ev = 27.2113860217
bohr_ang = 0.52917721090380
kb = 8.617330350e-5            # eV
au_fs = hbar/hartree_ev
alpha_fc = 0.007297352569311  #  Fine structure constant
c_au = 1/alpha_fc             # Speed of light in au
ε_au = 1/(4*pi)
Na = 6.02214076 * 10**23
kcalmol_ev = 23.061

# Conversion factors are derived in units_calcs.py 
W_FACTOR =  0.6373391727317127
FLUO_RATE_FACTOR =  2.8571698162427876e-9# (e*bohr)^2/ε hbar c^3 in units of s^2
EXT_FACTOR =  154.64232186052422 # (e*bohr)^2 Na /ε hbar c in units of m^2
