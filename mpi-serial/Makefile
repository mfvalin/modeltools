FC=gfortran
CC=cc

lib:	libmpi-serial.a

fort.o:	fort.F90
	$(FC) -c fort.F90

libmpi-serial.a: collective.o comm.o group.o handles.o list.o mpi.o pack.o recv.o req.o send.o time.o fort.o
	ar rcv libmpi-serial.a *.o

clean:
	rm -f *.o *.a
