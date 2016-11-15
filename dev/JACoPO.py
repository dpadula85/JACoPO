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

import Coup
from Opts import *
import ParseInput as PI
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
wn2eV = 8065.73

def calc_com(coords, atoms):

    masses = np.array([ ELEMENTS[x].mass for x in atoms ])
    com = np.dot(coords.T, masses) / np.sum(masses)

    return com


def format_selection(intlist):

    s = ''
    for i in intlist:
        s += '%3d ' % (i + 1) 

    return s


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


def print_dict(opts_dict, title=None, outstream=None):

    if outstream:
        sys.stdout = outstream

    if not title:
        title = "Options"

    print(banner(ch="=", length=60))
    print(title)
    print
    fmt = "%-20s %-20s"
    # print(fmt % ("# Option", "Value"))
    for k, v in sorted(opts_dict.iteritems()):

        if type(v) is str:
            pass

        if type(v) is int or type(v) is float:
            v = str(v)

        if type(v) is list:
            v = ', '.join(map(str, v))

        print(fmt % (k, v))

    print(banner(ch="=", length=60))
    print

    return


if __name__ == '__main__':

    Opts = options()

    if Opts['OutFile']:
        sys.stdout = open(Opts['OutFile'], 'w')

    start = time.time()
    name = ' ' * 24 + 'JACoPO' + ' ' * 24
    print(banner(ch='#', length=60))
    print(banner(text=name, ch='#', length=60))
    print(banner(ch='#', length=60))
    print
    print('JACoPO: Just Another COupling Program, Obviously')
    print('JACoPO Copyright (C) 2016 Daniele Padula, Marco Campetella')
    
    if Opts['Verb'] > 2:
        print
        print_dict(Opts)

    #
    # Parse Input Files
    #
    inigeo1 = None
    inigeo2 = None
    fingeo1 = None
    fingeo2 = None

    # Get geometries in dimer
    if Opts['FinGeo1File']:
        at1, fingeo1 = PI.read_geo(Opts['FinGeo1File'])

    if Opts['FinGeo2File']:
        at2, fingeo2 = PI.read_geo(Opts['FinGeo2File'])


    if Opts['Coup'] == 'chgs':

        # Monomer 1
        if Opts['IniGeo1File']:
            dum1, inigeo1 = PI.read_geo(Opts['IniGeo1File'])

        chgs1 = PI.read_chg(Opts['Chgs1File'])

        # Monomer 2
        if Opts['IniGeo2File']:
            dum2, inigeo2 = PI.read_geo(Opts['IniGeo2File'])

        chgs2 = PI.read_chg(Opts['Chgs2File'])

    elif Opts['Coup'] == 'tdc':

        # Monomer 1
        cub1 = PI.Cube(Opts['Cub1File'])
        at1 = np.array([ ELEMENTS[x].symbol for x in cub1.atoms[:,0] ])
        inigeo1 = cub1.atoms[:,1:]
        TrDen1 = cub1.data
        grid1 = cub1.grid
        dV1 = cub1.dV

        # Monomer 2
        cub2 = PI.Cube(Opts['Cub2File'])
        at2 = np.array([ ELEMENTS[x].symbol for x in cub2.atoms[:,0] ])
        inigeo2 = cub2.atoms[:,1:]
        TrDen2 = cub2.data
        grid2 = cub2.grid
        dV2 = cub2.dV


    #
    # Assign final geometry
    #
    if fingeo1 is None:
        fingeo1 = inigeo1
        inigeo1 = None

    if fingeo2 is None:
        fingeo2 = inigeo2
        inigeo2 = None

    #
    # Determine transformation matrices
    #
    if Opts['Verb'] > 2:
        print
        print(banner(ch="=", length=60))
        print("Monomer 1")
        print
        print
        print("%d" % len(at1))
        print
        for i in range(len(at1)):
            coor = fingeo1[i] * au2ang
            atom = [ at1[i], coor[0], coor[1], coor[2] ]
            print("%-5s %14.8f %14.8f %14.8f" % tuple(atom))
 
 
        print(banner(ch="=", length=60))
        print

        print
        print(banner(ch="=", length=60))
        print("Monomer 2")
        print
        print
        print("%d" % len(at2))
        print
        for i in range(len(at2)):
            coor = fingeo2[i] * au2ang
            atom = [ at2[i], coor[0], coor[1], coor[2] ]
            print("%-5s %14.8f %14.8f %14.8f" % tuple(atom))
 
 
        print(banner(ch="=", length=60))
        print

    # Monomer 1
    if inigeo1 is not None and fingeo1 is not None:

        inigeo1_rmsd = np.copy(inigeo1)
        fingeo1_rmsd = np.copy(fingeo1)

        if Opts['Sel1Geo']:
            sel1geo =  PI.read_sel(Opts['Sel1Geo'])
            inigeo1_rmsd = inigeo1[sel1geo]

        if Opts['Sel1Cub']:
            sel1cub =  PI.read_sel(Opts['Sel1Cub'])
            inigeo1_rmsd = inigeo1[sel1cub]

        RMSD1, M1, T11, T21 = Coup.kabsch(fingeo1_rmsd, inigeo1_rmsd)

        if Opts['Verb'] > 1:

            if Opts['IniGeo1File']:
                print
                print(banner(ch="=", length=60))
                print('Transformation of geometry from %s to %s' % (Opts['IniGeo1File'], Opts['FinGeo1File']))
                print
                print('RMSD (Ang): %14.8f' % (RMSD1 * au2ang))
                print(banner(ch="=", length=60))

            elif Opts['Cub1File']:
                print
                print(banner(ch="=", length=60))
                print('Transformation of geometry from %s to %s' % (Opts['Cub1File'], Opts['FinGeo1File']))
                print
                print('RMSD (Ang): %14.8f' % (RMSD1 * au2ang))
                print(banner(ch="=", length=60))

        #
        # Transform properties
        #
        if grid1 is not None:
            grid1 = grid1 - T11
            grid1 = np.dot(grid1, M1)
            grid1 = grid1 + T21

        if Opts['Dip1File']:
            dip1ext = np.loadtxt(Opts['Dip1File'])
            dip1ext = np.dot(dip1ext, M1)
            dip1extmod = np.linalg.norm(dip1ext)

        if Opts['SaveCub'] and Opts['Cub1File']:
        
            transfcub = PI.Cube(Opts['Cub1File'])
        
            new_x = np.atleast_2d(transfcub.X)
            new_x = np.dot(new_x, M1)
            transfcub.X = new_x.reshape(3).tolist()
        
            new_y = np.atleast_2d(transfcub.Y)
            new_y = np.dot(new_y, M1)
            transfcub.Y = new_y.reshape(3).tolist()
        
            new_z = np.atleast_2d(transfcub.Z)
            new_z = np.dot(new_z, M1)
            transfcub.Z = new_z.reshape(3).tolist()
        
            new_origin = np.atleast_2d(transfcub.origin)
            new_origin = new_origin - T11
            new_origin = np.dot(new_origin, M1)
            new_origin = new_origin + T21
            transfcub.origin = new_origin.reshape(3).tolist()
        
            for i, atom in enumerate(transfcub.atoms):
                atom[1] = fingeo1[i,0]
                atom[2] = fingeo1[i,1]
                atom[3] = fingeo1[i,2]
        
            with open("transf1.cub", "w") as f:
                transfcub.dump(f)


    # Monomer 2
    if inigeo2 is not None and fingeo2 is not None:

        inigeo2_rmsd = np.copy(inigeo2)
        fingeo2_rmsd = np.copy(fingeo2)

        if Opts['Sel2Geo']:
            sel2geo =  PI.read_sel(Opts['Sel2Geo'])
            inigeo2_rmsd = inigeo2[sel2geo]

        if Opts['Sel2Cub']:
            sel2cub =  PI.read_sel(Opts['Sel2Cub'])
            inigeo2_rmsd = inigeo2[sel2cub]

        RMSD2, M2, T12, T22 = Coup.kabsch(fingeo2_rmsd, inigeo2_rmsd)

        if Opts['Verb'] > 1:

            if Opts['IniGeo2File']:
                print
                print(banner(ch="=", length=60))
                print('Transformation of geometry from %s to %s' % (Opts['IniGeo2File'], Opts['FinGeo2File']))
                print
                print('RMSD (Ang): %14.8f' % (RMSD2 * au2ang))
                print(banner(ch="=", length=60))

            elif Opts['Cub2File']:
                print
                print(banner(ch="=", length=60))
                print('Transformation of geometry from %s to %s' % (Opts['Cub2File'], Opts['FinGeo2File']))
                print
                print('RMSD (Ang): %14.8f' % (RMSD2 * au2ang))
                print(banner(ch="=", length=60))

        #
        # Transform properties
        #
        if grid2 is not None:
            grid2 = grid2 - T12
            grid2 = np.dot(grid2, M2)
            grid2 = grid2 + T22

        if Opts['Dip2File']:
            dip2ext = np.loadtxt(Opts['Dip2File'])
            dip2ext = np.dot(dip2ext, M2)
            dip2extmod = np.linalg.norm(dip2ext)

        if Opts['SaveCub'] and Opts['Cub2File']:
        
            transfcub = PI.Cube(Opts['Cub2File'])
        
            new_x = np.atleast_2d(transfcub.X)
            new_x = np.dot(new_x, M2)
            transfcub.X = new_x.reshape(3).tolist()
        
            new_y = np.atleast_2d(transfcub.Y)
            new_y = np.dot(new_y, M2)
            transfcub.Y = new_y.reshape(3).tolist()
        
            new_z = np.atleast_2d(transfcub.Z)
            new_z = np.dot(new_z, M2)
            transfcub.Z = new_z.reshape(3).tolist()
        
            new_origin = np.atleast_2d(transfcub.origin)
            new_origin = new_origin - T12
            new_origin = np.dot(new_origin, M2)
            new_origin = new_origin + T22
            transfcub.origin = new_origin.reshape(3).tolist()
        
            for i, atom in enumerate(transfcub.atoms):
                atom[1] = fingeo2[i,0]
                atom[2] = fingeo2[i,1]
                atom[3] = fingeo2[i,2]
        
            with open("transf2.cub", "w") as f:
                transfcub.dump(f)

    #
    # Calculate Couplings
    #
    if not Opts['SkipCoup']:

        #
        # Charges
        #
        if Opts['Coup'] == 'chgs':

            # Charges
            coup = Coup.coup_chgs(fingeo1, chgs1, fingeo2, chgs2)
            dip1 = Coup.dipole_chgs(fingeo1, chgs1)
            dip1mod = np.linalg.norm(dip1)
            dip2 = Coup.dipole_chgs(fingeo2, chgs2)
            dip2mod = np.linalg.norm(dip2)

            # PDA
            center1 = calc_com(fingeo1, at1)
            center2 = calc_com(fingeo2, at2)
            coup_PDA, orifac = Coup.coup_PDA(center1, dip1, center2, dip2)

            if Opts['Dip1File'] and Opts['Dip2File']:
                coup_PDA_ext, orifac_ext = Coup.coup_PDA(center1, dip1ext, center2, dip2ext)

        #
        # TDC
        #
        if Opts['Coup'] == 'tdc' and FModule:

            # Dipoles
            dip1 = trden.diptrde(TrDen1, grid1, dV1)
            dip1mod = np.linalg.norm(dip1)
            dip2 = trden.diptrde(TrDen2, grid2, dV2)
            dip2mod = np.linalg.norm(dip2)

            # TDC
            coup = trden.couptrde(TrDen1, grid1, dV1, TrDen2, grid2, dV2, Opts['Thresh'])

            # PDA
            center1 = calc_com(fingeo1, at1)
            center2 = calc_com(fingeo2, at2)
            coup_PDA, orifac = Coup.coup_PDA(center1, dip1, center2, dip2)

            if Opts['Dip1File'] and Opts['Dip2File']:
                coup_PDA_ext, orifac_ext = Coup.coup_PDA(center1, dip1ext, center2, dip2ext)


        #
        # Print Results
        #
        print
        print(banner(ch="=", length=60))
        print('Results of the calculation with the %s method' % Opts['Coup'].upper())
        print

        if Opts['Verb'] > 0:
            print
            print(banner(ch="-", length=60))
            print(banner(text='Transition Dipole Moments (au)', ch=' ', length=60))
            print
            print(banner(text='  Monomer 1  ', ch="-", length=60))
            print('Method        x            y            z          norm')
            print(banner(ch="-", length=60))
            print('%s   %12.6f %12.6f %12.6f %12.6f' % (Opts['Coup'].upper(), dip1[0], dip1[1], dip1[2], dip1mod))

            if Opts['Dip1File']:
                print('%s   %12.6f %12.6f %12.6f %12.6f' % ('Ext.', dip1ext[0], dip1ext[1], dip1ext[2], dip1extmod))

            print
            print(banner(text='  Monomer 2  ', ch="-", length=60))
            print('Method        x            y            z          norm')
            print(banner(ch="-", length=60))
            print('%s   %12.6f %12.6f %12.6f %12.6f' % (Opts['Coup'].upper(), dip2[0], dip2[1], dip2[2], dip2mod))

            if Opts['Dip2File']:
                print('%s   %12.6f %12.6f %12.6f %12.6f' % ('Ext.', dip2ext[0], dip2ext[1], dip2ext[2], dip2extmod))

            print(banner(ch="-", length=60))
            print
            print('Orientation Factor %s: %8.4f' % (Opts['Coup'].upper(), orifac))

            if Opts['Dip1File'] and Opts['Dip2File']:
                print('Orientation Factor %s: %8.4f' % ('Ext.', orifac_ext))

        print
        print(banner(ch="-", length=60))
        print(banner(text='Couplings', ch=' ', length=60))
        print
        print('Method                         cm-1           eV')
        print(banner(ch="-", length=60))
        print('%s                %16.8f %16.8e' % (Opts['Coup'].upper(), coup, coup / wn2eV))
        print('PDA Dip %s        %16.8f %16.8e' % (Opts['Coup'].upper(), coup_PDA, coup_PDA / wn2eV))

        if Opts['Dip1File'] and Opts['Dip2File']:
            print('PDA Dip %s        %16.8f %16.8e' % ('Ext.', coup_PDA_ext, coup_PDA_ext / wn2eV))

        print(banner(ch="=", length=60))

        elapsed = (time.time() - start)
        print
        print("Calculation Time: %s" % time.strftime("%H:%M:%S", time.gmtime(elapsed)))
