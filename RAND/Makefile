
SOURCES = randomgeneric.c random_gaussian.c   random_mt19937.c  random_r250.c  \
          random_shr3.c random_xorshiftrotate.c  random_xsr128.c  random_xsr128r.c \
          randomgeneric.h

all: ftest ctest

self_tests: random_r250.Abs random_mt19937.Abs random_shr3.Abs random_xsr128.Abs random_xsr128r.Abs

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

random_r250.Abs: librandom.a random_generic_test.c
	s.cc -mpi -O2 -ftree-vectorize -mfma -mavx2 -DTEST_R250 random_generic_test.c -L. -lrandom -o $@
	./$@
	rm -f $@

random_mt19937.Abs: librandom.a random_generic_test.c
	s.cc -mpi -O2 -ftree-vectorize -mfma -mavx2 -DTEST_MT19937 random_generic_test.c -L. -lrandom -o $@
	./$@
	rm $@

random_shr3.Abs: librandom.a random_generic_test.c
	s.cc -mpi -O2 -ftree-vectorize -mfma -mavx2 -DTEST_SHR3 random_generic_test.c -L. -lrandom -o $@
	./$@
	rm $@

random_xsr128.Abs: librandom.a random_generic_test.c
	s.cc -mpi -O2 -ftree-vectorize -mfma -mavx2 -DTEST_XSR128 random_generic_test.c -L. -lrandom -o $@
	./$@
	rm $@

random_xsr128r.Abs: librandom.a random_generic_test.c
	s.cc -mpi -O2 -ftree-vectorize -mfma -mavx2 -DTEST_XSR128R random_generic_test.c -L. -lrandom -o $@
	./$@
	rm $@

clean:
	rm -f ctest ftest *.o *.a

