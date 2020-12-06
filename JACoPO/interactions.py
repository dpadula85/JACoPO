#!/usr/bin/env python

# JACoPO.py: calculation of electronic couplings with various approaches.
# Copyright (C) 2016  Daniele Padula, Marco Campetella
# dpadula85@yahoo.it, marco.campetella82@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from JACoPO import tdc
from scipy.spatial.distance import cdist


def dipole_chgs(coords, qs):
    '''
    Function to compute the dipole due to a set of point charges.

    Parameters
    ----------
    coords: np.array (N,3).
        coordinates in au.
    qs: np.array (N).
        charges in au.

    Returns
    -------
    dipole: np.array (3).
        dipole in au.
    '''

    dipole = np.dot(qs, coords)

    return dipole


def dipole_den(cub):
    '''
    Function to compute the dipole due to an electronic density represented
    on a grid.

    Parameters
    ----------
    cub: Cube class instance.
        electronic density.

    Returns
    -------
    dipole: np.array (3).
        dipole in au.
    '''

    den = cub.data
    grid = cub.grid
    dV = cub.dV
    dipole = np.dot(den, grid) * dV

    return dipole


def coul_chgs(coords1, qs1, coords2, qs2, d=4.0):
    '''
    Function to compute coulombic interaction between two set of charges.

    Parameters
    ----------
    coords1: np.array (N,3).
        coordinates in au.
    qs1: np.array (N).
        charges in au.
    coords2: np.array (N,3).
        coordinates in au.
    qs2: np.array (N).
        charges in au.
    d: float (default: 4).
        distance threshold in au for gaussian smoothing of point charge.

    Returns
    -------
    coul: float.
        coulombic interaction in au.
    '''

    # Compute distance matrix
    D = cdist(coords1, coords2)

    # Compute interaction smoothing factors at short distances
    factors = np.ones_like(D)
    smooth = D[D < d]
    smooth = np.exp(-(smooth - d)**2 / 1.5)
    factors[D < d] = smooth

    # Compute coulombic interaction with smoothing prefactors
    R = factors / D
    coul = np.dot(np.dot(qs1, R), qs2)
    
    return coul


def coul_PDA(center1, dip1, center2, dip2):
    '''
    Function to compute coulombic interaction between two dipoles within the
    Point Dipole Approximation.

    Parameters
    ----------
    center1: np.array (3).
        coordinates in au.
    dip1: np.array (3).
        dipole in au.
    center2: np.array (3).
        coordinates in au.
    dip2: np.array (3).
        dipole in au.

    Returns
    -------
    coul: float.
        coulombic interaction in au.
    orifac: float.
        orientation factor between the two dipoles.
    '''

    r = center2 - center1
    rmod = np.linalg.norm(r)
    ur = r / rmod

    try:
        dip1mod = np.linalg.norm(dip1)
        udip1 = dip1 / dip1mod
        dip2mod = np.linalg.norm(dip2)
        udip2 = dip2 / dip2mod

    except RuntimeWarning:
        return (0, 0)

    orifac = np.dot(udip1, udip2) - 3 * (np.dot(udip1, ur) * np.dot(udip2, ur))
    coul = dip1mod * dip2mod * orifac / rmod**3

    return coul, orifac


def coul_TDC(cub1, cub2, thresh=1e-5):
    '''
    Function to compute coulombic interaction between electronic densities
    within the Transition Density Cube approximation. This function is a
    wrapper to a fortran call that actually carries out the calculation.

    Parameters
    ----------
    cub1: Cube class instance.
        electronic density.
    cub2: Cube class instance.
        electronic density

    Returns
    -------
    coul: float.
        coulombic interaction in au.
    '''

    den1 = cub1.data
    grid1 = cub1.grid
    dV1 = cub1.dV

    den2 = cub2.data
    grid2 = cub2.grid
    dV2 = cub2.dV

    coul = tdc.coultdc(den1, grid1, dV1, den2, grid2, dV2, thresh)

    return coul


if __name__ == '__main__':
    pass
