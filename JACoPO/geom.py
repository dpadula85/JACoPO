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
from scipy.spatial.distance import cdist


def centroid(coords, masses=None):
    '''
    Function to compute the centre (or the centre of mass) of a set of
    coordinates.

    Parameters
    ----------
    coord: np.array (N,3).
        coordinates.
    masses: np.array (N) (default: None).
        masses.

    Returns
    -------
    com: np.array (3).
        centre (or centre of mass) of the set of coordinates.
    '''

    com = np.average(coords, axis=0, weights=masses)

    return com


def rmse(x, y):
    '''
    Function to compute the Root Mean Square Error between two sets of points.

    Parameters
    ----------
    x: np.array (N,M).
        set of points.
    y: np.array (N,M).
        set of points.

    Returns
    -------
    err: float.
        Root Mean Squared Error.
    '''

    err = np.sqrt(np.mean( (x - y)**2))

    return err


# These come from https://github.com/charnley/rmsd
def kabsch(coords, ref):
    '''
    Function to transform a set of coordinates to match a reference.
    The transformation is done through the Kabsch algorithm. The order
    of points does not matter, as they get reordered to do calculations.

    Parameters
    ----------
    coords: np.array (N,D).
        set of points.
    ref: np.array (N,D).
        reference of points.

    Returns
    -------
    transf: np.array (N,D)
        set of points transformed to match the reference points.
    '''

    # Check for equal dimensions
    assert len(coords) == len(ref)

    P = coords
    Q = ref

    # # Reorder input coordinates according to their distance from their
    # # centroids
    # P, pidxs = _reorder_com(coords)
    # Q, qidxs = _reorder_com(ref)

    # Compute centroids
    com1 = centroid(P)
    com2 = centroid(Q)

    # Translate coordinates in the origin
    P -= com1
    Q -= com2

    # Get the optimal rotation matrix
    U = _kabsch(P, Q)

    # Rotate P unto Q
    P = np.dot(P, U)

    # Translate P unto ref
    P += com2

    # # Reverse P order to the originary one
    # transf = np.zeros_like(P)
    # transf[pidxs] = P
    transf = P

    return transf, U


def _kabsch(P, Q):
    '''
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters
    ----------
    P : np.array
        (N,D) matrix, where N is points and D is dimension.
    Q : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    U : np.array
        Rotation matrix (D,D)
    '''

    # Computation of the covariance matrix
    C = np.dot(P.T, Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def _reorder_com(P):
    '''
    Function to reorder xyz coordinates by distance of each atom from the
    centroid.

    Parameters
    ----------
    P : np.array

    Returns
    -------
    P : np.array
    idxs : np.array
    '''

    # Compute centroid
    com1 = centroid(P)

    # Calculate distance from centroid
    # Convert com1 to a 1,D matrix as required by cdist
    # Convert the result back to rank0 array
    D = cdist(P, com1.reshape(-1,com1.shape[0])).reshape(-1)

    # Get order from closer to further
    idxs = D.argsort()

    # Sort 
    P = P[idxs]

    return P, idxs


if __name__ == '__main__':
    pass
