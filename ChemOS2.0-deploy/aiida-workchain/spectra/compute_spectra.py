#!/usr/bin/env python
import numpy as np
# import cclib
import re
from constants import *
import periodictable

# Everything below is taken from the paper:
# --------------------------------------------------------------------
# On the theoretical prediction of fluorescence rates from first
# principles using the path integral approach
# by de Souza, Neese and Izsak
# J. Chem. Phys. 148, 034104 (2018); https://doi.org/10.1063/1.5010895
# --------------------------------------------------------------------
# Equation numbers refer to specific equations in said paper.

def compute_jt(temperature,
               nt, maxT,
               wi, Li, qi,
               wf, Lf, qf,
               freq_cutoff=10.0, emission=False):
    """Compute the response J(t) for a set of HOs.

    This is the main function that does all of the heavy lifting. It
    takes as arguments physical parameters (temperature, broadening
    etc.) and harmonic oscillator parameters (w,L,q) for two
    electronic surfaces and output the correlation function,
    
    J(t) = Tr [e^-i H_f t e^-i H_i t]
    
    where H_i and H_f are the hamiltonians of the two HO surfaces.
    This is done in a number of ways (not currently implemented)
    depending on the surfaces. Currently, this function only does the
    independent-mode displaced harmonic oscillator (IMDHO) approach,
    that is where Li = Lf and wi = wf.
    
    TODO:documentation

    Parameters
    ----------
    temperature:
        Temperature in Kelvin.
    nt:
        Number of time slices.
    maxT:
        Maximum time (in fs).
    wi, Li, qi:
        Initial surface parameters. wi and Li should be in the form
        return by normal_modes(), that is, as an array of
        frequencies^2 in au and a set of mass-weighted normal modes.
        qi is the equilibrium position of the initial surface in
        mass-weighted coordinates.
    wf, Lf, qf:
        Same as above but for the final electronic surface.
    
    Optional parameters
    -------------------
    freq_cutoff:
        A frequency cutoff in cm^-1 for wi and wf. Defaults to 10
        cm^-1. 

    Returns
    -------
    t: ndarray
        The time variable (in au).

    jt: ndarray
        The correlation function J(t) = rho^FC(t) in the paper.
    """

    t = np.linspace(0, maxT, nt) * au_fs # Time in au
    N = len(wf)

    # Frequencies of oscillators in eV
    wi_ev = np.sqrt(wi + 0.0j) * W_FACTOR
    wf_ev = np.sqrt(wf + 0.0j) * W_FACTOR

    whichi = wi_ev > freq_cutoff / ev_wn
    whichf = wf_ev > freq_cutoff / ev_wn

    wi_ev = wi_ev[whichi].real
    wf_ev = wf_ev[whichi].real

    # Frequencies of oscillators in inv time au
    wi_au = (wi_ev/hbar)/au_fs
    wf_au = (wf_ev/hbar)/au_fs

    # Compute beta = 1/kT (in au)
    beta = 1/(kb * temperature / hartree_ev) 

    # Compute Duchinsky rotation and vibrational displacement (eq. 7)
    J = Li[:,whichi].T.dot(Lf[:,whichf])
    K = np.sqrt(wi_au) * Li[:,whichi].T.dot(qf.flatten() - qi.flatten())
    # NOTE: extra factor of sqrt(wi)! Paper doesn't include but orca does
    # in its diagnostic message sum(K*K) = whatever. It's required for K
    # to be unitless.

    jt = np.zeros(len(t), dtype=np.complex)

    # Is J = I?
    is_J_eq_I = np.allclose(J , np.eye(len(J)))

    # is wf == wi?
    is_wf_eq_wi = np.allclose(wi, wf)

    # The article
    # JCTC 2013, 9, 9, 4097-4115
    # points out that the sign of time has to be reversed when
    # looking at emission compared to absoprtion.
    if emission:
        S = -1.0
    else:
        S =  1.0

    if is_J_eq_I and is_wf_eq_wi:
        # If J = I and wf = wi, we have the independent mode displaced
        # HO type of dynamics.
        jt_IMDHO(S*t,jt,beta,K, wi_au)
    elif is_J_eq_I:
        jt_DODF(S*t,jt,beta,K,wi_au,wf_au)
    else:
        jt_GEN(S*t,jt,beta,K,J,wi_au,wf_au)

    return t, jt


def broaden_jt(Gamma, Omega, t, jt):
    """Broaden the correlation function.

    see eq. 49 and 52 in JCTC 2013, 9, 9, 4097-4115.

    Parameters
    ----------
    Gamma: 
        Homogeneous (Lorentzian) broadening (in cm^-1). 
    Omega:
        Inhomogeneous (Gaussian) broadening (in cm^-1).
    t: ndarray
        The time variable (in au).
    jt: ndarray
        The correlation function J(t) = rho^FC(t) in the paper.

    Returns
    -------
    jt: ndarray
        The correlation function J(t) with broadening.
    """

    # Convert gamma and omega from cm^-1 to au
    gam = ((Gamma / ev_wn)/hbar) / au_fs
    ome = ((Omega / ev_wn)/hbar) / au_fs

    dt = t[1]-t[0]
    Nt = np.ones_like(jt)
    if ome > 0.0:
        Nt *=  np.exp(-t**2.0 * ome**2.0/2.0)

    if gam > 0.0:
        Nt *= np.exp(-gam * abs(t))

    return jt * Nt


def jt_GEN(t, jt, beta, K, J, wi_au, wf_au):
    """Inplace update of jt.

    See sec. IIB of paper."""
    raise NotImplementedError("Currently only IMDHO is implmented.")

def jt_DODF(t, jt, beta, K, wi_au, wf_au):
    """DODF inplace update of jt.

    See sec. IID of paper."""
    raise NotImplementedError("Currently only IMDHO is implmented.")

def jt_IMDHO(t, jt, beta, K, wi_au):
    """Independent-mode Displaced HO in-place update.

    See sec. IIE of paper."""
    for ti, time in enumerate(t):
        tau = -time - 1j * beta         # in au
        tau_bar = time                  # in au
        gamma1 = tau_bar * wi_au
        beta1 = tau * wi_au
        # eq.40 (except that K is = sqrt(wi_au) K)
        Sii = K**2  / (np.tan(np.pi/2 - beta1/2) + np.tan(np.pi/2 - gamma1/2))
        # eq. 39
        jt[ti] = np.exp(-1j * np.sum(Sii))


def jt2jw(t, jt):
    # Carrier frequency
    dt = t[1]-t[0]              # in au

    # eq. 13
    # J(w) = 2 Re int_0^inf dt e^{i w t} J(t)
    #      = Re int_-inf^inf dt e^{i w t} J(t)
    jw = np.fft.fft(jt) * dt/(2*np.pi)       # Units of au
    jw = np.fft.fftshift(jw)[::-1]                     # -iw -> iw
    # Normalization by 2 pi due to convention of DFT.

    # Make sure jw -> 0 as w -> inf
    jw -= jw[-1]

    # w in units of 1/au
    w = np.fft.fftshift(np.fft.fftfreq(len(t), dt) * 2 * np.pi)
    return w, jw

def compute_extinction(w, jw, n_solvent, transition_dipole, Ecarrier,
                       use_sb_formula=False):
    wcarrier = (Ecarrier/hbar) / au_fs

    # Final frequency axis in ev
    w_ev = w * au_fs  * hbar + Ecarrier
    F = np.sum(transition_dipole**2)
    if use_sb_formula:
        # Normalize absorption using strickler-berg expression.
        # see https://doi.org/10.1039/B412936A, eq 2 and 3

        # kr = 2.88e-9 n^2 <ν>^2 int ε(ν) dν = 0.668 <ν>^2 n^2 f

        # Therefore,
        # int ε(ν) dν = 0.668  f / 2.88e-9  = 231.9e6 f
        # let ε(ν) = A E(ν) where int E(ν) d ν = 1, then
        # A  = 231.9e6 f

        f = (2.0/3.0) * (Ecarrier/hartree_ev) * F
        nu = w_ev * ev_wn
        dnu = nu[1]-nu[0]

        A = 231.9e6 * f 
        extinct_m2_sb = A * jw.real/(np.sum(jw.real) * dnu)
        return w_ev, extinct_m2_sb

    else:
        # eq.1 in JCTC 2013, 9, 9, 4097-4115
        alpha = 2*((10.0/3.0) * np.pi /np.log(10)) * EXT_FACTOR * (w+wcarrier)

        # Extra factor of 2 is obtained from comparison with
        # strickler-berg expression.

        # Extinction in m^2
        extinct_m2 = alpha * jw.real * F

        # Convert to cm M^-1
        return w_ev, extinct_m2


def compute_fluorescence(w, jw, n_solvent, transition_dipole, Ecarrier,
                         use_sb_formula=False):
    wcarrier = (Ecarrier / hbar) / au_fs

    # Oscillator strength
    F = np.sum(transition_dipole**2)

    # Final frequency axis in ev
    w_ev = w * au_fs  * hbar - Ecarrier
    dw = w[1]-w[0]

    if use_sb_formula:
        # Compute rate using strickler-berg expression.
        # see https://doi.org/10.1039/B412936A, eq 3

        f = (2.0/3.0) * (-Ecarrier/hartree_ev) * F
        nu = w_ev * ev_wn
        dnu = nu[1]-nu[0]
        nu_avg = np.sum(nu * jw.real)/np.sum(jw.real)
        kr_sb = 0.668 * nu_avg**2 * n_solvent**2 * f


        # Correct fluo such that S fluo dw = kr_sb
        fluo = jw.real
        norm = np.sum(jw.real) * dw
        fluo = fluo/norm * kr_sb

        return w_ev, fluo, kr_sb


    else:
        # Use eq.1 in JCTC 2013, 9, 9, 4097-4115

        # alpha_OPE is for the emission intensity, the energy emitted by
        # one mole per unit time

        # alpha_OPE = 2 * Na/(3 * ε * c**3) * w**4

        # We convert that to the amount of photons emitted by dividing by
        # one factor of hbar w and remove the mole dependence

        # alpha = 2/(3 * ε * c**3) * w**3/hbar

        # We note that the equation is missing a factor of n^2 / pi from
        # eq. 10 in the main paper

        alpha = (2.0/(3.0*np.pi)) * FLUO_RATE_FACTOR \
            * ((w-wcarrier) * au_fs)**3 * n_solvent**2


        # # Fluo spectra in fs of time
        fluo = alpha * jw.real * F * 1e15
        kr = np.sum(fluo) * dw # in inverse s of time

        return w_ev, fluo, kr

def parse_xyz(
    xyz_str: str
):
    # _, _, *coords = xyz_str.strip().split("\n")
    # elements, coords = zip(*map(
    #     lambda xs: (xs[0], xs[1:4])
    #     , map(
    #         lambda s: s.split()
    #         , coords)
    # ))
    elements = []
    coords = []
    with open(xyz_str, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                natom = int(line)
            elif 2<= i <= natom+2:
                elements.append(line.split()[0])
                coords.append([float(x) for x in line.split()[1:4]])
    masses = np.array([getattr(periodictable, element).mass for element in elements])
    return (masses, np.asarray(coords, dtype=float))
