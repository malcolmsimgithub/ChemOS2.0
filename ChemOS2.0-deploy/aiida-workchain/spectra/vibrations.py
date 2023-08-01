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
