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

try:
    import trden
    FModule = True

except ImportError:
    print(" WARNING!!!")
    print(" The Fortran Module could not be loaded.")
    print(" Coupling from Transition Densities will not be computed.")
    FModule = False

# Constants

au2ang = 0.5291771
au2wn = 2.194746e5

def dipole_chgs(struct, chgs):
    '''Calculates a dipole from a set of coordinates and atomic charges.'''

    return np.dot(struct.T, chgs)


def coup_chgs(struct1, chgs1, struct2, chgs2):
    '''Calculates Electronic Coupling from Transition Charges'''

    #
    # Convert Coordinates to au
    #
    coup = 0
    factor = 1
    for i in range(len(struct1)):
        for j in range(len(struct2)):
    
            a1 = struct1[i]
            chg1 = chgs1[i]
    
            a2 = struct2[j]
            chg2 = chgs2[j]
    
            d = np.linalg.norm(a1 - a2)

            #
            # Gaussian Blur factor for short distances
            #
            if d < 4:
                factor = np.exp(-(d - 4)**2/1.5)
    
            coup += factor * chg1 * chg2 / d
            factor = 1
    
    return coup * au2wn


def coup_PDA_OLD(struct1, atoms1, dip1, struct2, atoms2, dip2):
    '''Calculates Electronic Coupling according to the Point Dipole
    Approximation.'''

    com1 = np.dot(struct1.T, atoms1) / np.sum(atoms1)
    com2 = np.dot(struct2.T, atoms2) / np.sum(atoms2)

    r = com2 - com1
    rmod = np.linalg.norm(r)
    ur = r / rmod
    dip1mod = np.linalg.norm(dip1)
    udip1 = dip1 / dip1mod
    dip2mod = np.linalg.norm(dip2)
    udip2 = dip2 / dip2mod

    coup = dip1mod * dip2mod * (np.dot(udip1, udip2) - 3 * (np.dot(udip1, ur) * np.dot(udip2, ur))) / rmod**3
    orifac = (np.dot(udip1, udip2) - 3 * (np.dot(udip1, ur) * np.dot(udip2, ur)))

    return coup * au2wn, orifac


def coup_PDA(center1, dip1, center2, dip2):
    '''Calculates Electronic Coupling according to the Point Dipole
    Approximation.'''

    r = (center2 - center1)
    rmod = np.linalg.norm(r)
    ur = r / rmod

    try:
        dip1mod = np.linalg.norm(dip1)
        udip1 = dip1 / dip1mod
        dip2mod = np.linalg.norm(dip2)
        udip2 = dip2 / dip2mod

    # Needed for dipoles with norm = 0
    # I am able to catch this warning as an exception because in the imported
    # module QM_parser.util.util the following line is present:
    # warnings.filterwarnings("error")
    except RuntimeWarning:
        return (0, 0)

    orifac = np.dot(udip1, udip2) - 3 * (np.dot(udip1, ur) * np.dot(udip2, ur))
    coup = dip1mod * dip2mod * orifac / rmod**3

    return coup * au2wn, orifac


def kabsch(struct1, struct2):
    '''Returns the RMSD calculated with Kabsch's algorithm.'''

    # to np.array
    # struct1 = np.array([ [atom[1], atom[2], atom[3]] for atom in struct1 ])    
    # struct2 = np.array([ [atom[1], atom[2], atom[3]] for atom in struct2 ])    

    # check for consistency in number of atoms
    assert len(struct1) == len(struct2)
    L = len(struct1)
    assert L > 0

    # Center the two fragments to their center of coordinates
    com1 = np.sum(struct1, axis=0) / float(L)
    com2 = np.sum(struct2, axis=0) / float(L)
    struct1 -= com1
    struct2 -= com2

    # Initial residual, see Kabsch.
    E0 = np.sum(np.sum(struct1 * struct1, axis=0), axis=0) + \
         np.sum(np.sum(struct2 * struct2, axis=0), axis=0)

    # This beautiful step provides the answer. V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = np.linalg.svd(np.dot(np.transpose(struct2), struct1))

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation. V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    # The rotation matrix U is simply V*Wt
    U = np.dot(V, Wt)
 
    # rotate and translate the molecule
    # struct2 = np.dot((struct2), U)
    # struct2 = struct2 + com1

    return RMSD, U, com2, com1


if __name__ == '__main__':
    pass
