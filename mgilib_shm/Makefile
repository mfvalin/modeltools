
all: mgilib2.o atm.Abs oce.Abs test.Abs mpi.Abs mgi_shm.Abs libmgi.a


mgi_shm.Abs: mgi_shm.c
	s.cc -I. -O 2 -Dmgi_shm=main mgi_shm.c -o mgi_shm.Abs

mpi.Abs: mgilib2.o f_mgi_test.F90
	s.f90 -mpi -DMODEL f_mgi_test.F90 mgilib2.o -o mpi.Abs

atm.Abs: mgilib2.o f_mgi_test.F90
	s.f90 -mpi -DMODEL -DATM f_mgi_test.F90 mgilib2.o -o atm.Abs

oce.Abs: mgilib2.o f_mgi_test.F90
	s.f90 -mpi -DMODEL -DOCE f_mgi_test.F90 mgilib2.o -o oce.Abs

test.Abs: mgilib2.o f_mgi_test.F90
	s.f90 f_mgi_test.F90 mgilib2.o -o test.Abs

mgilib2.o: mgilib2.c mgi.h
	s.cc -I. -c -DWITHOUT_GOSSIP -O 2 mgilib2.c -o mgilib2.o

libmgi.a: mgilib2.o
	ar rcv libmgi.a mgilib2.o

clean:
	rm -f mgilib2.o atm.Abs oce.Abs mgi_shm.Abs mpi.Abs libmgi.a

