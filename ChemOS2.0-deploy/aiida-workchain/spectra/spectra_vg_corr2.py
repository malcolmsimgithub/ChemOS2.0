#!/usr/bin/env python

import numpy as np
import sys
from numpy import pi
import pandas as pd

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

def spectrum_vg():
    
    verbose = False
    masses, x0 = parse_xyz("xtb_results/S0_XTB_xyz.xyz") 
    H = np.load('final_step/hessian.npy')
    grad = np.load('final_step/grads.npy')
    dipo = np.load('final_step/dipmoments.npy')
    ee = np.load('final_step/excenergies.npy')

    # Smiles file
    try:
        f = open('molecule.smi')
        smiles_string = f.readline()[:-1]
        f.close()
    except FileNotFoundError:
        if verbose:
            print("No smiles string found.")
        smiles_string = ""
        # TODO: GENERATE FROM XYZ


    correlation_only = False
    if correlation_only:
        import uuid

    # Parameters for the spectra computation
    T = 298.15
    Gamma = 50.0
    Omega = 100.0
    cutoff = (50.0/ev_wn) /hartree_ev
    n_solvent = 1.344
    use_sb_formula = False

    # state_pairs = [[int(c) for c in s.split(',')] for s in args.state_pairs]
    # if len(args.state_pairs) == 1 and args.state_pairs[0] == 'all':
    state_pairs = [(0,i) for i in range(1,len(ee))] + [(1,0)]
    # else:
    #     state_pairs = [[int(c) for c in s.split(',')] for s in args.state_pairs]

    # TODO auto-magic computation of parameters
    # Parameters for correlation computations
    nt = 15000    # Number of time points
    maxT = 2000.0 # fs

    # Grid parameters
    on_grid = True
    if on_grid:
        g_low, g_high, g_n, g_unit = "0.5,8,2000,eV".split(',')
        if g_unit == 'nm':
            x_nm = np.linspace(float(g_low), float(g_high), int(g_n))
            x_ev = ev_nm/x_nm
        elif g_unit == 'eV':
            x_ev = np.linspace(float(g_low), float(g_high), int(g_n))
            x_nm = ev_nm/x_ev
        else:
            raise NotImplementedError("Grid string %s not understood"
                                    % "0.5,8,2000,eV".split(','))

        y_abs = np.zeros_like(x_ev)
        y_fluo = np.zeros_like(x_ev)


    # Start by transforming the Hessian to normal modes
    w,L = normal_modes(masses, H)
    if verbose:
        print("Normal modes and frequencies have been computed.")
        count = 0
        f = ""
        for wi in w:
            f += "%7.3f" % wi
            count +=1
            if count == 5:
                count =0
                f+= "\n"

        print("Frequencies^2 in au:\n" + f)


    # Move x0 to COM
    x0_com = move_to_com(masses,x0)

    # Starting point (mass-weighted)
    q0 = massweight(masses, x0_com.flatten())

    # Compute estimate of the minimum for all the surfaces
    qs = []
    for g in grad:
        qs += [vertical_gradient(masses, x0_com, w, L, g, freq_cutoff=cutoff)]

    # S0 energy (for conformer weighting)
    E_S0 = ee[0] - vib_potential(L, w, qs[0], q0)*27.2 

    # Start the computation
    for el_i, el_f in state_pairs:
        qi = qs[el_i].copy()
        qf = qs[el_f].copy()
        mu = dipo[el_i,el_f]

        # Eckart's condition
        # ------------------
        # Technically, I don't think thats necessary in the vertical
        # gradient calculations. It doesn't seem to do anything.
        #eT = eckart_T_matrix(qi,qf)
        #qi = qi.dot(eT)
        #Lt = nm_apply_coordinate_transform(L, eT)

        # qi is the minima of Ei but the transition energies are computed
        # at q0. We use the vibrational energy at q0 to correct the
        # surface energies at the minimums qi and qf.
        Ei = ee[el_i] - vib_potential(L, w, qi, q0) *27.2
        Ef = ee[el_f] - vib_potential(L, w, qf, q0) *27.2

        # The vertical excitation energy at qi = (V_f(qi) + E_f) - (V_i(qi) + Ei)
        carrier = (Ef + vib_potential(L, w, qf, qi)*27.2) - (Ei + vib_potential(L, w, qi, qi)*27.2)

        # Oscillator strength
        f = (2.0/3.0) * (carrier/hartree_ev) * np.sum(abs(mu)**2)

        output = dict(
            smiles = smiles_string,
            S0_energy = E_S0 * hartree_ev,
            state_i = el_i,
            state_f = el_f,
            mu_x = mu[0],
            mu_y = mu[1],
            mu_z = mu[2],
            osc_strength=f,
            electronic_gap=carrier,
            temperature=T,
        )

        if el_i < el_f:
            typ = "Absorption  "
            emi = False
        else:
            typ = "Fluorescence"
            emi = True

        if verbose:
            print("\n=========================================")
            print("%s  |%i>  --> |%i> " %(typ, el_i,el_f))
            print("μ²       = %10.5f au²" % np.sum(abs(mu)**2))
            print("Osc Str. = %10.5f" % f)
            print("ΔV       = %10.5f eV" % carrier)
            print("\nEckart T matrix")
            # print("  %10.5f %10.5f %10.5f" % tuple(eT[0,:]))
            # print("  %10.5f %10.5f %10.5f" % tuple(eT[1,:]))
            # print("  %10.5f %10.5f %10.5f" % tuple(eT[2,:]))
            sys.stdout.write("\n\nComputing correlation function...")

        # Compute correlation function
        t_au, corr_t = compute_jt(T, nt, maxT,
                                w, L, qi,
                                w, L, qf, emission=emi)

        if verbose:
            sys.stdout.write("  done!\n")

        if correlation_only:
            # If we were just computing the correlation functions, we save all
            # of them here. In the future, this will be default / only output.
            output.update(time=t_au, Jt=corr_t)
            np.savez("corr_%i_%i.npz" % ( el_i, el_f ), **output)


        else:
            # COMPUTE AND SAVE SPECTRAS
            if verbose:
                sys.stdout.write("Applying a Fourier transform...")

            corr_t_broadened = broaden_jt(Gamma, Omega, t_au, corr_t)
            freq_au, corr_w = jt2jw(t_au, corr_t_broadened) 
            output.update(
                w=freq_au,
                Jw=corr_w,
                homo_broadening=Gamma,
                inhomo_broadening=Omega)

            if verbose:
                sys.stdout.write("  done!\n")

            if typ == "Absorption  ":
                # Extinction in units of cm M^-1
                we1, abs1 = compute_extinction(freq_au, corr_w, n_solvent, mu, carrier,
                                               use_sb_formula=use_sb_formula)
                dw = we1[1]-we1[0]

                # Some other derived quantities
                maxe = np.argmax(abs1)
                nabs1 = abs1/np.sum(abs1)
                wcom = np.sum(nabs1 * we1)
                stw = np.sqrt(np.sum(nabs1 * we1**2) - wcom**2)

                if verbose:
                    print("Spectra computed!")
                    print("peak ε(ω)     %7.3f eV / %7.3f nm" % (we1[maxe], ev_nm/we1[maxe]))
                    print("ext at peak   %7.3f / M cm" % (abs1[maxe]))
                    print("spectral com  %7.3f eV / %7.3f nm" % (wcom, ev_nm/wcom))
                    print("         std  %7.3f eV" % (stw))
                    print("saving to file...")


                output.update(
                    n_solvent=n_solvent,
                    hw=we1,
                    ext_spectra=abs1,
                    ext_spectra_com=wcom,
                    ext_spectra_std_dev=stw,
                    ext_peak=abs1[maxe],
                    ext_peak_w=we1[maxe])

                np.savez("final_step/abs_%i_%i.npz" % ( el_i, el_f ), **output)

                if on_grid:
                    y_abs += np.interp(x_ev, we1, abs1)

            else:
                # Extinction in units of cm M^-1
                wf1, fluo1, rate = compute_fluorescence(freq_au, corr_w, n_solvent,
                                                        mu, carrier,
                                                        use_sb_formula=use_sb_formula)
                dw = wf1[1]-wf1[0]

                maxe = np.argmax(fluo1)
                nfluo1 = fluo1/(np.sum(fluo1))
                wcom = np.sum(wf1 * nfluo1)
                stw = np.sqrt(np.sum(nfluo1 * wf1**2) - wcom**2)

                if verbose:
                    print("Spectra computed!")
                    print("peak I(ω)     %7.3f eV / %7.3f nm" % (wf1[maxe], ev_nm/wf1[maxe]))
                    print("Fluo. k= %12.3f / ns  | τ= %12.3f ns " % (rate * 1e-9, 1/(rate*1e-9)))
                    print("spectral com  %7.3f eV / %7.3f nm" % (wcom, ev_nm/wcom))
                    print("         std  %7.3f eV" % (stw))
                    print("saving to file...")

                output.update(
                    n_solvent=n_solvent,
                    hw=wf1,
                    fluo_spectral_fluence = fluo1, # Units of photons / s freq
                    fluo_spectra_com=wcom,
                    fluo_spectra_std_dev=stw,
                    fluo_peak_w=wf1[maxe],
                    fluo_rate=rate)
                np.savez("final_step/fluo_%i_%i.npz" % ( el_i, el_f ), **output)

        
                if on_grid:
                    y_fluo += np.interp(x_ev, wf1, fluo1)

    if on_grid:
        out = pd.DataFrame(
            {"x_ev": x_ev,
            "x_nm": x_nm,
            "extinct":y_abs,
            "fluo":y_fluo}).set_index("x_ev")

        out.to_csv("final_step/specs.csv")


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


#!/usr/bin/env python
import numpy as np
# import cclib
import re
import periodictable


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

    jt = np.zeros(len(t), dtype=complex)

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


#!/usr/bin/env python
import numpy as np
from scipy import linalg as spl

def inv_massweight(masses,xi):
    return xi.flatten() / np.sqrt(np.repeat(masses,3))

def massweight(masses,xi):
    return xi.flatten() * np.sqrt(np.repeat(masses,3))

def get_eq_position(L, w, gq):
    # b = H^-1 g
    #   = L diag(w2^-1) L^T g
    b = np.zeros_like(gq)
    b = L.T.dot(gq)/w
    b = L.dot(b)
    return b

def move_to_com(masses,xi):
    # Make some copies
    x1 = xi.copy()

    # Translate to COM
    x1_com = masses.dot(x1)/np.sum(masses)
    x1[:,0] -= x1_com[0]
    x1[:,1] -= x1_com[1]
    x1[:,2] -= x1_com[2]
    return x1

def normal_modes(masses, H):
    """Compute normal modes from the Hessian matrix H.

    Parameters
    ----------
    masses : ndarray
        Atomic masses in amu.
    H : ndarray
        Hessian matrix in units of hartree/bohr^2. Not mass-weighted.
    
    Returns
    -------
    W2: ndarray
        Array of vibrational frequencies^2 in au.

    Lx: ndarray
        Normal modes (mass-weighted).
    """
    mcoord = np.repeat(masses, 3) # mass for each cart. coordinate

    # Compute the frequencies in a.u. and the normal mode matrix C
    W2, C = spl.eigh(H, np.diag(mcoord))

    # mass-weight C to obtain Lx
    Lx = np.einsum('i,ij->ij', np.sqrt(mcoord), C)

    return W2, Lx

def eckart_T_matrix(qi,qf):
    """Compute Eckart's T matrix for coordinates xi and xf.

    See paper at
    https://aip.scitation.org/doi/10.1063/1.1864872

    Parameters
    ----------
    qi, qf: ndarray
        Mass-weighted coordinates.

    Returns
    -------
    T: ndarray
        3x3 matrix corresponding to the Eckart transformation.
    """

    # eq. 3 in paper
    A = np.einsum('mi,mj->ij', qi, qf)

    A1 = A.dot(A.T)                  # above eq7
    A2 = A.T.dot(A)

    e1, xi = np.linalg.eigh(A1)     # eq8
    e2, eta = np.linalg.eigh(A2)    # eq9

    #  eq.7 
    T = np.einsum('ik,jk->ij', eta, xi)
    return T

def nm_apply_coordinate_transform(L, T):
    """Apply a coordinate transform T to normal modes L."""
    # Apply coordinate transform T to normal modes L
    new_L = np.zeros_like(L)
    natoms = L.shape[0]//3
    for j in range(L.shape[1]):
        v = L[:,j].reshape((natoms, 3))
        new_L[:,j] = v.dot(T).flatten()
    return new_L

def vib_potential(L, w, q0, q):
    """Compute the vibrational potential energy at point q.

    Returns the vibrational potential energy V(q) for the harmonic
    oscillators defined by L,w and q0.
    """
    return np.sum(np.abs(L.T.dot(q.flatten()-q0.flatten()))**2 * w)/2

def vertical_gradient(masses, xi, w, L, grad, freq_cutoff=1e-3):
    """Vertical gradient/IMDHO approximation to the eq. geom.

    This function computes the equilibrium geometry from gradient and
    normal mode information."""
    # Vertical gradient. Compute eq. geom. from...
    # x1 -> angstroms
    # xf = x1 - H^-1 grad
    # xf = x1 - b
    which = w>freq_cutoff

    gq = inv_massweight(masses,grad.flatten())
    b = inv_massweight(masses,get_eq_position(L[:,which], w[which], gq))
    q_out = massweight(masses, xi.flatten() - b)

    return q_out.reshape(xi.shape)
