FC=gfortran
CC=cc

lib:	libmpi-serial.a

fort.o:	fort.F90
	$(FC) -c fort.F90

libmpi-serial.a: collective.o comm.o group.o handles.o list.o mpi.o pack.o recv.o req.o send.o time.o fort.o
	$(CC) -c collective.c comm.c group.c handles.c list.c mpi.c pack.c recv.c req.c send.c time.c
	$(FC) -c fort.F90
	ar rcv libmpi-serial.a *.o

clean:
	rm -f *.o *.a
