#!/usr/bin/env python

import numpy.f2py.f2py2e as f2py2e
import sys

sys.argv +=  "--opt='-O3' -lgomp --f90flags='-fopenmp' -c -m trden trden.f90 ".split()
f2py2e.main()
