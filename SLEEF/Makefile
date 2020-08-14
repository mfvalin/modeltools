
CC = clang
FC = flang
CFLAGS = -march=core-avx2
NP = 1005

tests : sleef_wrapper.o ftn_sleef_test.F90 sleef.inc
	$(FC) ftn_sleef_test.F90 sleef_wrapper.o -I. -L/opt/sleef/lib -lsleef -lm -O3 ${LD_SET_RPATH} -DNP=${NP} -o test.Abs
	gfortran ftn_sleef_test.F90 sleef_wrapper.o -I. -L/opt/sleef/lib -lsleef -lm -O3 ${LD_SET_RPATH} -DNP=${NP} -o test_gnu.Abs

sleef.inc : functions.h
	grep '^//[ !]' functions.h | sed 's:^//::' >sleef.inc

functions.h : FUNCTION.sh  FUNCTION_1_2.sh  FUNCTION_2_1.sh FUNCTIONS.sh
	./FUNCTIONS.sh >functions.h

sleef_wrapper.o : sleef_wrapper.c functions.h
	$(CC) -c sleef_wrapper.c -I/opt/sleef/include -I. -O3 $(CFLAGS)

clean :
	rm -f test.Abs functions.h sleef_wrapper.o sleef.inc
