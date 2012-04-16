# this is for vanilla gcc/gfortran/mpicc/mpif90 to be adjusted on target system

CFLAGS = -fpic
FCFLAGS = -fpic -fno-second-underscore -fcray-pointer

CC=gcc
FC=gfortran
MPICC = mpicc
MPIFC = mpif90

# end of vanilla gcc/gfortran/mpicc/mpif90 defaults

all: lib

clean:
	rm -f *.o *.a *.mod

lib: fortran_code c_code
	ar rcv libtest.a *.o

modules:
	$(MPIFC) $(FCFLAGS) -c RPN_COMM_mod.f90

fortran_code: modules
	$(MPIFC) $(FCFLAGS) -c RPN_COMM*.f

f77name:
	cd ../tools ;\
        gcc -c ftn-mangling-detect.c ;\
        $(FC) $(FCFLAGS) ftn-mangling.f ftn-mangling-detect.o
	../tools/a.out >f77name.h
	cd ../tools ;\
	rm -f a.out *.o

c_code: f77name
	$(CC) $(CFLAGS) -c RPN_COMM*.c
