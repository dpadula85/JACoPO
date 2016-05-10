.PHONY: all clean 

all: trden.so

trden.so: trden.f90
	f2py  --opt='-O3' -lgomp --f90flags='-fopenmp' -c -m trden trden.f90 
	
alternative:
	./compile.py

clean:
	rm -rf trden.so
