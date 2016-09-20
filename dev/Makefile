.PHONY: all clean 

all: trden.so

trden.so: trden.f90
	f2py  --opt='-O3 -funroll-loops' -lgomp --f90flags='-fopenmp' -c -m trden trden.f90 
	
alternative:
	f2py  --opt='-03 -funroll-loops' -lgomp --f90flags='-fopenmp' --f90exec=/usr/bin/gfortran -c -m trden trden.f90

clean:
	rm -rf trden.so
