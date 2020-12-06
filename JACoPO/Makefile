.PHONY: all clean 

all: tdc.so

tdc.so: tdc.f90
	f2py  --opt='-O3 -funroll-loops' -lgomp --f90flags='-fopenmp' -c -m tdc tdc.f90 
	
clean:
	rm -rf tdc.*.so
