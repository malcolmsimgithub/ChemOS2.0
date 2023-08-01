#!/usr/bin/env python
import numpy as np
import Units as units

U = units.SI()

print("Unit conversion factors")
print("-----------------------")
# Units of the Hessian
# --------------------
hessian = U("hartree / (amu * bohr * bohr)")

# w = np.sqrt(lambda) where lambda are the hessian eigenvalues
# hbar * w should have energy units
hbarw = U("hbar")* hessian**0.5
print("W_FACTOR = ", (hbarw/U("eV")).value)

# Units of extinction
# -------------------
au_dipole = U("bohr * charge_e") 
print("EXT_FACTOR = ", au_dipole**2 * U("Na/(vac_e*hbar*c)"))

# Fluorescence rate
# -----------------
print("FLUO_RATE_FACTOR = ", au_dipole**2/(U("vac_e * hbar * c**3")))
