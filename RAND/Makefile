
SOURCES = randomgeneric.c random_gaussian.c   random_mt19937.c  random_r250.c  \
          random_shr3.c random_xorshiftrotate.c  random_xsr128.c  random_xsr128r.c

all: ftest ctest

interfaces: randomfunctions.h randomfunctions.inc

randomfunctions.h: ${SOURCES}
	echo '#include <randomgeneric.h>' >randomfunctions.h
	cat  ${SOURCES} | grep -w InTc | sed 's:[\t ]*//[\t ]*!InTc!.*:;:' >> randomfunctions.h

randomfunctions.inc: ${SOURCES}
	cat  ${SOURCES} | grep -w InTf | sed 's:[\t ]*!InTf!.*::' > randomfunctions.inc

ftest:	librandom.a demo_rand.F90
	rm -f *.o
	s.f90 -openmp demo_rand.F90 -L. -lrandom -o ftest
	rm -f *.o

ctest: librandom.a random_test.c
	rm -f *.o
	s.cc -mpi random_test.c -L. -lrandom -o ctest
	rm -f *.o

librandom.a: randomfunctions.h ${SOURCES}
	rm -f *.o
	s.cc -O3 -c ${SOURCES}
	ar rcv librandom.a *.o
	rm -f *.o

clean:
	rm -f ctest ftest *.o *.a

