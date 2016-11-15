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

import argparse as arg

# Constants

au2ang = 0.5291771
au2wn = 2.194746e5

def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(description='''Calculates Electronic Coupling
                                            from Transition Charges and Densities.''',
                                            formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input files
    #
    inp = parser.add_argument_group("Input Data")

    # Monomer 1
    inp.add_argument('--chg1', default=None, type=str, dest='Chgs1File',
                     help='''Charges for monomer 1.''')

    inp.add_argument('--geo1', default=None, type=str, dest='FinGeo1File',
                     help='''Geometry of monomer 1 in the dimer. (Units: Angstrom)''')

    inp.add_argument('--selgeo1', default=None, nargs='+', type=str, dest='Sel1Geo',
                     help='''Atom Selection for monomer 1. This can either be a list
                     or a file.''')

    inp.add_argument('--refgeo1', default=None, type=str, dest='IniGeo1File',
                     help='''Geometry of monomer 1. (Units: Angstrom)''')

    inp.add_argument('--cub1', default=None, type=str, dest='Cub1File',
                     help='''Transition Density Cube for monomer 1.''')

    inp.add_argument('--selcub1', default=None, nargs='+', type=str, dest='Sel1Cub',
                     help='''Atom Selection for Transition Density Cube for monomer 1.
                     This can either be a list or a file.''')

    inp.add_argument('--fac1', default=1.0, type=float, dest='Fac1',
                     help='''Scaling factor for monomer 1.''')

    inp.add_argument('--dip1', default=None, type=str, dest='Dip1File',
                     help='''Dipole file for monomer 1.''')


    # Monomer 1
    inp.add_argument('--chg2', default=None, type=str, dest='Chgs2File',
                     help='''Charges for monomer 2.''')

    inp.add_argument('--geo2', default=None, type=str, dest='FinGeo2File',
                     help='''Geometry of monomer 2 in the dimer. (Units: Angstrom)''')

    inp.add_argument('--selgeo2', default=None, nargs='+', type=str, dest='Sel2Geo',
                     help='''Atom Selection for monomer 2. This can either be a list
                     or a file.''')

    inp.add_argument('--refgeo2', default=None, type=str, dest='IniGeo2File',
                     help='''Geometry of monomer 2. (Units: Angstrom)''')

    inp.add_argument('--cub2', default=None, type=str, dest='Cub2File',
                     help='''Transition Density Cube for monomer 2.''')

    inp.add_argument('--selcub2', default=None, nargs='+', type=str, dest='Sel2Cub',
                     help='''Atom Selection for Transition Density Cube for monomer 2.
                     This can either be a list or a file.''')

    inp.add_argument('--fac2', default=1.0, type=float, dest='Fac2',
                     help='''Scaling factor for monomer 2.''')

    inp.add_argument('--dip2', default=None, type=str, dest='Dip2File',
                     help='''Dipole file for monomer 2.''')


    #
    # Calculations Options
    #
    calc = parser.add_argument_group("Calculation Options")

    calc.add_argument('--coup', default=None, type=str, choices=['chgs', 'tdc'],
                      help='''Method of Calculation of the Electronic Coupling.
                      The choice is exclusive.''', dest='Coup', required=True)

    calc.add_argument('--thresh', default=1e-5, type=float, dest='Thresh',
                      help='''Threshold for Transition Density Cubes.''')

    calc.add_argument('--nocoup', default=False, action="store_true",
                      help='''Skip coupling calculation.''', dest='SkipCoup')


    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--output', default=None, type=str, dest='OutFile',
                     help='''Output File.''')
    out.add_argument('--savecub', default=False, action="store_true", dest='SaveCub',
                     help='''Save Transition Density Cubes after transformation
                     in space.''')

    out.add_argument('-v', '--verbosity',
                     default=0, action="count", dest="Verb",
                     help='''Verbosity level''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts

if __name__ == '__main__':
    pass
