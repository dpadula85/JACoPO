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

import os
import sys
import time
import numpy as np
import argparse as arg
from elements import ELEMENTS

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

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='Calculates Electronic Coupling from Transition Charges and Densities.', formatter_class=arg.ArgumentDefaultsHelpFormatter)

    # Optional arguments
    parser.add_argument('--chg1', default=None, type=str, help='''File with coordinates and charges for monomer 1.''')
    parser.add_argument('--chg2', default=None, type=str, help='''File with coordinates and charges for monomer 2.''')

    parser.add_argument('--cub1', default='mon1.cub', type=str, help='''Transition Density Cube for monomer 1.''')
    parser.add_argument('--selcub1', default=None, nargs='+', type=str,
            help='''Atom Selection for Transition Density Cube for monomer 1. This can either be a list or a file.''')

    parser.add_argument('--geo1', default=None, type=str,
            help='''Geometry on which the Transition Density Cube of monomer 1 will be projected. (Units: Angstrom)''')
    parser.add_argument('--selgeo1', default=None, nargs='+', type=str,
            help='''Atom Selection for geometry 1. This can either be a list or a file, which should contain the list.''')

    parser.add_argument('--cub2', default='mon2.cub', type=str, help='''Transition Density Cube for monomer 2.''')
    parser.add_argument('--selcub2', default=None, nargs='+', type=str,
            help='''Atom Selection for Transition Density Cube for monomer 2. This can either be a list or a file, which should contain the list.''')

    parser.add_argument('--geo2', default=None, type=str,
            help='''Geometry on which the Transition Density Cube of monomer 2 will be projected. (Units: Angstrom)''')
    parser.add_argument('--selgeo2', default=None, nargs='+', type=str,
            help='''Atom Selection for geometry 2. This can either be a list or a file.''')

    parser.add_argument('--thresh', default=1e-5, type=float, help='''Threshold for Transition Density Cubes.''')

    parser.add_argument('--coup', default=None, type=str, choices=['chgs', 'den'],
            help='''Method of Calculation of the Electronic Coupling. If no method is specified, both methods will
            be used.''')

    parser.add_argument('-o', '--output', default='Coup.out', type=str, help='''Output File.''')

    args = parser.parse_args()

    return args


class CUBE:
    def __init__(self, fname):

        f = open(fname, 'r')
        for i in range(2): f.readline() # echo comment
        tkns = f.readline().split() # number of atoms included in the file followed by the position of the origin of the volumetric data
        self.natoms = int(tkns[0])
        self.origin = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])

        # The next three lines give the number of voxels along each axis (x, y, z) followed by the axis vector.
        tkns = f.readline().split() #
        self.NX = int(tkns[0])
        self.X = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])
        tkns = f.readline().split() #
        self.NY = int(tkns[0])
        self.Y = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])
        tkns = f.readline().split() #
        self.NZ = int(tkns[0])
        self.Z = np.array([float(tkns[1]),float(tkns[2]),float(tkns[3])])

        # The last section in the header is one line for each atom consisting of 5 numbers, the first is the atom number, second (?), the last three are the x,y,z coordinates of the atom center.
        self.atoms = []
        for i in range(self.natoms):
            tkns = map(float, f.readline().split())
            self.atoms.append([tkns[0], tkns[2], tkns[3], tkns[4]])

        # Volumetric data
        self.data = np.zeros((self.NX,self.NY,self.NZ))
        i=0
        for s in f:
            for v in s.split():
                self.data[i/(self.NY*self.NZ), (i/self.NZ)%self.NY, i%self.NZ] = float(v)
                i+=1
        if i != self.NX*self.NY*self.NZ: raise NameError, "FSCK!"


    def dump(self, f):

        # output Gaussian cube into file descriptor "f".
        # Usage pattern: f=open('filename.cube'); cube.dump(f); f.close()
        print >>f, "CUBE file\ngenerated by piton _at_ erg.biophys.msu.ru"
        print >>f, "%4d %.6f %.6f %.6f" % (self.natoms, self.origin[0], self.origin[1], self.origin[2])
        print >>f, "%4d %.6f %.6f %.6f"% (self.NX, self.X[0], self.X[1], self.X[2])
        print >>f, "%4d %.6f %.6f %.6f"% (self.NY, self.Y[0], self.Y[1], self.Y[2])
        print >>f, "%4d %.6f %.6f %.6f"% (self.NZ, self.Z[0], self.Z[1], self.Z[2])
        for atom in self.atoms:
            print >>f, "%s %d %s %s %s" % (atom[0], 0, atom[1], atom[2], atom[3])
        for ix in xrange(self.NX):
            for iy in xrange(self.NY):
                for iz in xrange(self.NZ):
                    print >>f, "%.5e " % self.data[ix,iy,iz],
                    if (iz % 6 == 5): print >>f, ''
                print >>f,  ""


    def mask_sphere(self, R, Cx,Cy,Cz):

        # produce spheric volume mask with radius R and center @ [Cx,Cy,Cz]
        # can be used for integration over spherical part of the volume
        m=0*self.data
        for ix in xrange( int(ceil((Cx-R)/self.X[0])), int(floor((Cx+R)/self.X[0])) ):
            ryz=np.sqrt(R**2-(ix*self.X[0]-Cx)**2)
            for iy in xrange( int(ceil((Cy-ryz)/self.Y[1])), int(floor((Cy+ryz)/self.Y[1])) ):
                rz=np.sqrt(ryz**2 - (iy*self.Y[1]-Cy)**2)
                for iz in xrange( int(ceil((Cz-rz)/self.Z[2])), int(floor((Cz+rz)/self.Z[2])) ):
                    m[ix,iy,iz]=1
        return m


def parse_TrDen(cubfile):

    TrDen1 = CUBE(cubfile)
    
    TrD1 = np.asfortranarray(TrDen1.data)
    
    # structure
    struct1 = np.array(TrDen1.atoms)
    
    # calculate the volume element
    dVx1 = TrDen1.X[0]
    dVy1 = TrDen1.Y[1]
    dVz1 = TrDen1.Z[2]
    
    # Grid points
    NX1 = TrDen1.NX
    NY1 = TrDen1.NY
    NZ1 = TrDen1.NZ
    
    # Origin of the cube
    O1 = TrDen1.origin

    return TrD1, dVx1, dVy1, dVz1, NX1, NY1, NZ1, O1, struct1


def dip_TrDen(TrDen, dVx, dVy, dVz, O):
    '''Calculates a dipole from a Transition Density cube.'''

    NX, NY, NZ = TrDen.shape
    dV = dVx * dVy * dVz
    mu = np.zeros(3)
    TrDen = TrDen * dV

    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                p0 = O[0] + i * dVx
                p1 = O[1] + j * dVy
                p2 = O[2] + k * dVz
                mu[0] += p0 * TrDen[i][j][k]
                mu[1] += p1 * TrDen[i][j][k]
                mu[2] += p2 * TrDen[i][j][k]

    return mu


def dipole_chgs(struct, chgs):
    '''Calculates a dipole from a set of coordinates and atomic charges.'''

    return np.dot(struct.T, chgs)


def coup_chgs(struct1, chgs1, struct2, chgs2):
    '''Calculates Electronic Coupling from Transition Charges as described in
    J. Phys. Chem. B, 2006, 110, 17268.'''

    coup = 0
    for i in range(len(struct1)):
        for j in range(len(struct2)):
    
            a1 = struct1[i]
            chg1 = chgs1[i]
    
            a2 = struct2[j]
            chg2 = chgs2[j]
    
            d = np.linalg.norm(a1 - a2)
    
            coup += chg1 * chg2 / d
    
    return coup * au2wn


def coup_PDA(struct1, atoms1, dip1, struct2, atoms2, dip2):
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

    return coup * au2wn


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


def process_selection(string):

    string =  ','.join(string).replace(',,',',')

    try:
        f = open(string, 'r')
        string = f.readlines()
        f.close()
        string =  ','.join(string).replace(',,',',')
        string = string.replace(',', ' ')
        string = map(lambda x: x - 1, extend_compact_list(string))

    except IOError:
        string = string.replace(',', ' ')
        string = map(lambda x: x - 1, extend_compact_list(string))

    return string


def format_selection(intlist):

    s = ''
    for i in intlist:
        s += '%3d ' % (i + 1) 

    return s


def read_geo(geofile):

    atgeo = np.loadtxt(geofile, usecols=[0], dtype="|S5")
    structgeo = np.loadtxt(geofile, usecols=[1,2,3]) / au2ang

    return atgeo, structgeo


def extend_compact_list(idxs):

    extended = []

    # Uncomment this line if idxs is a string and not a list
    idxs = idxs.split()

    for idx in idxs:

        to_extend = idx.split('-')

        if len(to_extend) > 1:

            sel =  map(int, to_extend)
            extended += range(sel[0],sel[1]+1,1)

        else:
        
            extended.append(int(idx))
    
    return extended


def banner(text=None, ch='=', length=78):
    """Return a banner line centering the given text.
    
        "text" is the text to show in the banner. None can be given to have
            no text.
        "ch" (optional, default '=') is the banner line character (can
            also be a short string to repeat).
        "length" (optional, default 78) is the length of banner to make.

    Examples:
        >>> banner("Peggy Sue")
        '================================= Peggy Sue =================================='
        >>> banner("Peggy Sue", ch='-', length=50)
        '------------------- Peggy Sue --------------------'
        >>> banner("Pretty pretty pretty pretty Peggy Sue", length=40)
        'Pretty pretty pretty pretty Peggy Sue'
    """
    if text is None:
        return ch * length

    elif len(text) + 2 + len(ch)*2 > length:
        # Not enough space for even one line char (plus space) around text.
        return text

    else:
        remain = length - (len(text) + 2)
        prefix_len = remain / 2
        suffix_len = remain - prefix_len
    
        if len(ch) == 1:
            prefix = ch * prefix_len
            suffix = ch * suffix_len

        else:
            prefix = ch * (prefix_len/len(ch)) + ch[:prefix_len%len(ch)]
            suffix = ch * (suffix_len/len(ch)) + ch[:suffix_len%len(ch)]

        return prefix + ' ' + text + ' ' + suffix


def checkfile(filename):

    if not os.path.isfile(filename):
        print(banner(text='ERROR', ch='#', length=80))
        print("File %s not found!" % filename)
        sys.exit()


if __name__ == '__main__':

    args = options()
    calctype = args.coup
    outfile = args.output

    start = time.time()
    #
    # Process Transition Charges
    #
    ChgsDone = False
    if not calctype or calctype == 'chgs':

        chg1file = args.chg1
        chg2file = args.chg2

        if chg1file and chg2file:

            checkfile(chg1file)
            checkfile(chg2file)

            data1 = np.genfromtxt(chg1file)
            data2 = np.genfromtxt(chg2file)

            atoms1 = np.genfromtxt(chg1file,usecols=[0], dtype="|S5")
            atoms1 = map(lambda x: ELEMENTS[x].mass, atoms1)
            struct1 = data1[:,1:4] / au2ang
            chgs1 = data1[:,-1]
            
            atoms2 = np.genfromtxt(chg2file,usecols=[0], dtype="|S5")
            atoms2 = map(lambda x: ELEMENTS[x].mass, atoms2)
            struct2 = data2[:,1:4] / au2ang
            chgs2 = data2[:,-1]

            coupchgs = coup_chgs(struct1, chgs1, struct2, chgs2)
            dip1chgs = dipole_chgs(struct1, chgs1)
            dip1chgsmod = np.linalg.norm(dip1chgs)
            dip2chgs = dipole_chgs(struct2, chgs2)
            dip2chgsmod = np.linalg.norm(dip2chgs)

            coup_PDA_chgs = coup_PDA(struct1, atoms1, dip1chgs, struct2, atoms2, dip2chgs)
            ChgsDone = True

    #
    # Process Transition Cubes
    #
    TrDenDone = False
    if FModule and (not calctype or calctype == 'den'):

        cub1file = args.cub1
        cub2file = args.cub2
        thresh = args.thresh
        selcub1 = args.selcub1
        selcub2 = args.selcub2
        geo1 = args.geo1
        geo2 = args.geo2
        selgeo1 = args.selgeo1
        selgeo2 = args.selgeo2

        if cub1file and cub2file:

            checkfile(cub1file)
            checkfile(cub2file)

            # Parse .cub files
            TrDenD, dVxD, dVyD, dVzD, NXD, NYD, NZD, OD, structD = parse_TrDen(cub1file)
            atomsD = structD[:,0]
            structD = structD[:,1:]

            TrDenA, dVxA, dVyA, dVzA, NXA, NYA, NZA, OA, structA = parse_TrDen(cub2file)
            atomsA = structA[:,0]
            structA = structA[:,1:]

            # Generate parameters for grid generation
            dVD = dVxD * dVyD * dVzD
            ND = NXD * NYD * NZD
            dVA = dVxA * dVyA * dVzA
            NA = NXA * NYA * NZA

            # Reshape TrDen 3D array to 1D array and rescale it to have its integral zero
            TrDenD = TrDenD.reshape(ND)
            TrDenA = TrDenA.reshape(NA)
            TrDenD -= sum(TrDenD) / len(TrDenD)
            TrDenA -= sum(TrDenA) / len(TrDenA)

            # Generate 4D array of grid points and reshape it to a 2D array
            gridD = trden.gengrid(OD, dVxD, dVyD, dVzD, NXD, NYD, NZD)
            gridD = gridD.reshape(ND,3)
            gridD = np.asfortranarray(gridD)
            gridA = trden.gengrid(OA, dVxA, dVyA, dVzA, NXA, NYA, NZA)
            gridA = gridA.reshape(NA,3)
            gridA = np.asfortranarray(gridA)

            # Rotate grids if reference geometries are provided
            if geo1:

                structD_rmsd = np.copy(structD)

                if selcub1:

                    selcub1 = process_selection(selcub1)
                    structD_rmsd = structD_rmsd[selcub1]

                # Get structure from final geometry file
                checkfile(geo1)
                atgeoD, structgeoD = read_geo(geo1)
                structgeoD_rmsd = np.copy(structgeoD)

                if selgeo1:

                    selgeo1 =  process_selection(selgeo1)
                    structgeoD_rmsd = structgeoD_rmsd[selgeo1]

                # Transform grid and structure according to the RMSD
                RMSDD, MD, T1D, T2D = kabsch(structgeoD_rmsd, structD_rmsd)

                structD = structD - T1D
                structD = np.dot(structD, MD)
                structD = structD + T2D

                gridD = gridD - T1D
                gridD = np.dot(gridD, MD)
                gridD = gridD + T2D

            if geo2:

                structA_rmsd = np.copy(structA)

                if selcub2:

                    selcub2 = process_selection(selcub2)
                    structA_rmsd = structA_rmsd[selcub2]

                # Get structure from final geometry file
                checkfile(geo2)
                atgeoA, structgeoA = read_geo(geo2)
                structgeoA_rmsd = np.copy(structgeoA)

                if selgeo2:

                    selgeo2 =  process_selection(selgeo2)
                    structgeoA_rmsd = structgeoA_rmsd[selgeo2]


                # Transform grid and structure according to the RMSD
                RMSDA, MA, T1A, T2A= kabsch(structgeoA_rmsd, structA_rmsd)

                structA = structA - T1A
                structA = np.dot(structA, MA)
                structA = structA + T2A

                gridA = gridA - T1A
                gridA = np.dot(gridA, MA)
                gridA = gridA + T2A


            # Dipoles
            dip1den = trden.diptrde(TrDenD, gridD, dVD)
            dip1denmod = np.linalg.norm(dip1den)
            dip2den = trden.diptrde(TrDenA, gridA, dVA)
            dip2denmod = np.linalg.norm(dip2den)

            # Coupling
            coupden = trden.couptrde(TrDenA, gridA, dVA, TrDenD, gridD, dVD, thresh)
            coup_PDA_den = coup_PDA(structD, atomsD, dip1den, structA, atomsA, dip2den)
            TrDenDone = True

    elapsed = (time.time() - start)
    elapsed = time.strftime("%H:%M:%S", time.gmtime(elapsed))
    #
    # Write a logfile with results
    #
    with open(outfile, 'w') as f:
        f.write('##############\n')
        f.write('##  JACoPO  ##\n')
        f.write('##############\n')
        f.write('\n')
        f.write('JACoPO: Just Another COupling Program, Obviously\n')
        f.write('JACoPO.py  Copyright (C) 2016  Daniele Padula, Marco Campetella\n')

        if ChgsDone:
            f.write('\n\n')
            f.write('###################################\n')
            f.write('##  Coupling Transition Charges  ##\n')
            f.write('###################################\n')
            f.write('\n')
            f.write('Coupling calculated from transition charges according to\n')
            f.write('J. Phys. Chem. B, 2006, 110, 17268\n')
            f.write('\n')
            f.write('Donor Structure and Charges:\n')
            f.write(open(chg1file).read())
            f.write('\n')
            f.write('Donor Electric Transition Dipole Moment from Transition Charges in a.u.:\n')
            f.write('%8s %8s %8s %15s\n' % ('x', 'y', 'z', 'norm'))
            f.write('%8.4f %8.4f %8.4f %15.4f\n' % (dip1chgs[0], dip1chgs[1], dip1chgs[2], dip1chgsmod))
            f.write('\n')
            f.write('\n')
            f.write('Acceptor structure and Charges:\n')
            f.write(open(chg2file).read())
            f.write('\n')
            f.write('Acceptor Electric Transition Dipole Moment from Transition Charges in a.u.:\n')
            f.write('%8s %8s %8s %15s\n' % ('x', 'y', 'z', 'norm'))
            f.write('%8.4f %8.4f %8.4f %15.4f\n' % (dip2chgs[0], dip2chgs[1], dip2chgs[2], dip2chgsmod))
            f.write('\n')
            f.write('\n')
            f.write('Electronic Coupling according to PDA from Dipoles from Transition Charges in cm-1:\n')
            f.write('%-10.2f\n' % coup_PDA_chgs)
            f.write('\n')
            f.write('Electronic Coupling in cm-1:\n')
            f.write('%-10.2f\n' % coupchgs)

        if TrDenDone:
            f.write('\n\n')
            f.write('#####################################\n')
            f.write('##  Coupling Transition Densities  ##\n')
            f.write('#####################################\n')
            f.write('\n')
            f.write('Coupling calculated from transition densities according to\n')
            f.write('J. Phys. Chem. B, 1998, 102, 5378\n')
            f.write('\n')
            f.write('\n')

            if geo1:
                f.write('Structure and Transition Density Cube in %s moved to match geometry in %s\n' % (cub1file, geo1))

                if selcub1 and selgeo1:
                    f.write('Atom correspondence between the cub and the geometry file, respectively:\n')
                    f.write('%s' % format_selection(selcub1))
                    f.write('\n')
                    f.write('%s' % format_selection(selgeo1))
                    f.write('\n')

                f.write('RMSD (Ang): %8.4f\n' % (RMSDD * au2ang))
                f.write('\n')


            f.write('Donor Electric Transition Dipole Moment from Transition Density in a.u.:\n')
            f.write('%8s %8s %8s %15s\n' % ('x', 'y', 'z', 'norm'))
            f.write('%8.4f %8.4f %8.4f %15.4f\n' % (dip1den[0], dip1den[1], dip1den[2], dip1denmod))
            f.write('\n')
            f.write('\n')

            if geo2:
                f.write('Structure and Transition Density Cube in %s moved to match geometry in %s\n' % (cub2file, geo2))

                if selcub2 and selgeo2:
                    f.write('Atom correspondence between the cub and the geometry file, respectively:\n')
                    f.write('%s' % format_selection(selcub2))
                    f.write('\n')
                    f.write('%s' % format_selection(selgeo2))
                    f.write('\n')

                f.write('RMSD (Ang): %8.4f\n' % (RMSDA * au2ang))
                f.write('\n')

            f.write('Acceptor Electric Transition Dipole Moment from Transition Density in a.u.:\n')
            f.write('%8s %8s %8s %15s\n' % ('x', 'y', 'z', 'norm'))
            f.write('%8.4f %8.4f %8.4f %15.4f\n' % (dip2den[0], dip2den[1], dip2den[2], dip2denmod))
            f.write('\n')
            f.write('\n')
            f.write('Electronic Coupling according to PDA from Dipoles from Transition Densities in cm-1:\n')
            f.write('%-10.2f\n' % coup_PDA_den)
            f.write('\n')
            f.write('Electronic Coupling in cm-1:\n')
            f.write('%-10.2f\n' % coupden)

        if not FModule:
            f.write('\n\n')
            f.write(" WARNING!!!\n")
            f.write(" The Fortran Module could not be loaded.\n")
            f.write(" Coupling from Transition Densities was not computed.\n")

        f.write('\n')
        f.write('\n')
        f.write('#####################\n')
        f.write('##  Results Table  ##\n')
        f.write('#####################\n')
        f.write('\n')
        f.write(' Calculation time: %s\n' % elapsed)
        f.write('\n')
        f.write('# Method          Coupling (cm-1)\n')
        f.write('#--------------------------------\n')

        if ChgsDone:
            f.write(' Tr Chgs         %8.2f \n' % coupchgs)
            f.write(' PDA Dip Chgs    %8.2f \n' % coup_PDA_chgs)

        if TrDenDone:
            f.write(' Tr Den          %8.2f \n' % coupden)
            f.write(' PDA Dip Den     %8.2f \n' % coup_PDA_den)
