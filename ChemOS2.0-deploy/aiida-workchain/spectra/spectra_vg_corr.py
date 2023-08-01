#!/usr/bin/env python

import numpy as np
import sys
import argparse
import pandas as pd

from constants import * 
from compute_spectra import *
from vibrations import *
# from rdkit.Chem import GetPeriodicTable

if __name__ == "__main__":
    # table = GetPeriodicTable()

    # Start parsing arguments
    parser = argparse.ArgumentParser(description="Compute the spectrum of a molecule using the Vertical Gradient method.")
    parser.add_argument("-v", help="Verbose.",action='store_true')

    # Input files
    parser.add_argument("-xyz", help=".xyz file for the molecular geometry (S0_opt.xyz).", default="S0_XTB_xyz.xyz", nargs="?")
    parser.add_argument("-hess", help="Hessian matrix file (hessian.npy).",
                        default='hessian.npy')
    parser.add_argument("-grad", help="Gradients file (grads.npy).",
                        default='grads.npy')
    parser.add_argument("-dipo", help="Transition dipoles (dipmoments.npy).",
                        default='dipmoments.npy')
    parser.add_argument("-ee", help="Excitation energies (excenergies.npy).",
                        default='excenergies.npy')

    # States to include
    parser.add_argument("state_pairs",
                        help="Pair of states to evaluate. "
                        + "Format is i,j k,l ... to compute transitions from i to j"
                        + " and k to l. Fluorescence properties are computed by default"
                        + " when j<i.", nargs='+')

    # Physical parameters for computation
    parser.add_argument("-T",
                        help="Temperature (Kelvin) --- defaults to 298.15K.",
                        default=298.15, type=float)
    parser.add_argument("-homo",
                        help="Homogeneous broadening (cm^-1) --- defaults to 50 cm^-1.",
                        default=50.0, type=float)
    parser.add_argument("-inhomo",
                        help="Inhomogeneous broadening (cm^-1) --- defaults to 100 cm^-1.",
                        default=100.0, type=float)
    parser.add_argument("-n",
                        help="Refractive index of solvent --- defaults to n=1.344 (Acetonitrile).",
                        default=1.344, type=float)
    parser.add_argument("--use-sb-formula",
                        help="Compute extinction and fluorescence spectra using a formula"
                        +" derived from the Strickler-Berg expression, as opposed to the"
                        +" \"first-principles\" approach. Results should be the same."
                        +" Defaults to false, as the Strickler-Berg formula is based on a"
                        +" slow-envelope approximation that does not hold for broad"
                        +" bands.",action="store_true")
    parser.add_argument("-cutoff",
                        help="Low frequency mode cutoff (cm^-1) --- defaults to 50 cm^-1.",
                        default=50.0, type=float)
    parser.add_argument("-corrN",
                        help="Number of time points for correlation function computation"
                        +" (defaults to 15000).",
                        default=15000, type=int)
    parser.add_argument("-corrT",
                        help="Maximum time for correlation function computation in fs"
                        +" (defaults to 2000 fs).",
                        default=2000.0, type=float)
    parser.add_argument("--correlations",
                        help="Only compute correlations, not spectra.",action="store_true")


    parser.add_argument("-grid",
                        help="If present, absorption and emission spectra will be integrated"
                        +" over an identical grid and exported as csv files.",
                        action='store_true')

    parser.add_argument("-grid-linspace",
                        help="Parameters defining the grid. These are a low cutoff, a high"
                        +" cutoff, the number of grid points and the unit (nm or eV), separated"
                        +" by spaces. Defaults to '0.5 8 2000 eV'.",
                        type=str, default="0.5,8,2000,eV".split(','), nargs=4)

    parser.add_argument("-grid-output",
                        help="Filename of csv gridded output. Defaults to specs.csv.",
                        type=str, default="specs.csv")

    args = parser.parse_args()
    verbose = args.v

    # XYZ file
    # cc_parser = cclib.io.ccopen(args.xyz)
    # data = cc_parser.parse()
    # x0 = data.atomcoords[-1,:,:] 
    # atomnos = data.atomnos
    # masses = np.array([table.GetMostCommonIsotopeMass(int(i)) for i in atomnos])

    masses, x0 = parse_xyz(args.xyz)
    #TODO: masses from periodic table



#     def parse_xyz(xyz_str: str):
#         _, _, *coords = xyz_str.strip().split("\n")
#         elements, coords = zip(*map(
#             lambda xs: (xs[0], xs[1:4])
#             , map(
#                 lambda s: s.split()
#                 , coords)
#         ))
#         return (elements, np.asarray(coords, dtype=float))

# with open("mol.xyz", 'r') as infile:
#     xyz = infile.read()

    # Other files
    H = np.load(args.hess)
    grad = np.load(args.grad)
    dipo = np.load(args.dipo)
    ee = np.load(args.ee)

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


    correlation_only = args.correlations
    if correlation_only:
        import uuid

    # Parameters for the spectra computation
    T = args.T
    Gamma = args.homo
    Omega = args.inhomo
    cutoff = (args.cutoff/ev_wn) /hartree_ev
    n_solvent = args.n
    use_sb_formula = args.use_sb_formula

    # state_pairs = [[int(c) for c in s.split(',')] for s in args.state_pairs]
    if len(args.state_pairs) == 1 and args.state_pairs[0] == 'all':
        state_pairs = [(0,i) for i in range(1,len(ee))] + [(1,0)]
    else:
        state_pairs = [[int(c) for c in s.split(',')] for s in args.state_pairs]

    # TODO auto-magic computation of parameters
    # Parameters for correlation computations
    nt = args.corrN    # Number of time points
    maxT = args.corrT # fs

    # Grid parameters
    on_grid = args.grid
    if on_grid:
        g_low, g_high, g_n, g_unit = args.grid_linspace
        if g_unit == 'nm':
            x_nm = np.linspace(float(g_low), float(g_high), int(g_n))
            x_ev = ev_nm/x_nm
        elif g_unit == 'eV':
            x_ev = np.linspace(float(g_low), float(g_high), int(g_n))
            x_nm = ev_nm/x_ev
        else:
            raise NotImplementedError("Grid string %s not understood"
                                    % args.grid_linspace)

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

                np.savez("abs_%i_%i.npz" % ( el_i, el_f ), **output)

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
                np.savez("fluo_%i_%i.npz" % ( el_i, el_f ), **output)

        
                if on_grid:
                    y_fluo += np.interp(x_ev, wf1, fluo1)

if on_grid:
    out = pd.DataFrame(
        {"x_ev": x_ev,
         "x_nm": x_nm,
         "extinct":y_abs,
         "fluo":y_fluo}).set_index("x_ev")

    out.to_csv(args.grid_output)
