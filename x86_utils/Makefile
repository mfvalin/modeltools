
all: fp_x87_xmm_control.o fp_x87_xmm_control.h fp_x87_xmm_control.inc

fp_x87_xmm_control.o : fp_x87_xmm_control.c
	gcc -c fp_x87_xmm_control.c

fp_x87_xmm_control.inc: fp_x87_xmm_control.c
	gcc -o a.out fp_x87_xmm_control.c -DPRINT_DEFS -lm
	./a.out f >fp_x87_xmm_control.inc
	rm -f a.out fp_x87_xmm_control.o
	gcc -c fp_x87_xmm_control.c
	rm -f a.out
	gfortran -fno-range-check fp_x87_xmm_control_test.F90 fp_x87_xmm_control.o -lm
	./a.out
	rm -f fp_x87_xmm_control_test.o a.out

fp_x87_xmm_control.h: fp_x87_xmm_control.c
	gcc -o a.out fp_x87_xmm_control.c -DPRINT_DEFS -lm
	./a.out f >fp_x87_xmm_control.h
	rm -f a.out

