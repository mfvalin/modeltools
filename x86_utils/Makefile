
ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

ifeq (/intel,$(findstring /intel,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := -no-wrap-margin
EXTRA_C := -march=core-avx2
endif

ifeq (gfortran,$(findstring gfortran,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C := -mavx2 -mfma
endif

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

intrp_bicub_yx_f.F90:
	ln -sf intrp_bicub_yx_f.F90 intrp_bicub_yx.c

intrp_bicub_yx.o: intrp_bicub_yx.c
	s.cc -O2 -DTIMING -c ${EXTRA_C} intrp_bicub_yx.c -DWITH_MONO

intrp_bicub_yx: intrp_bicub_yx.o intrp_bicub_yx_f.F90
	s.f90 ${EXTRA_FTN} -DF_TEST -o intrp_bicub_yx intrp_bicub_yx_f.F90 intrp_bicub_yx.o

clean:
	rm -f intrp_bicub_yx intrp_bicub_yx.o a.out fp_x87_xmm_control.o 
