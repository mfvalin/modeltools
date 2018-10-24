
ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

ifneq (,$(findstring intel,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := -no-wrap-margin
EXTRA_C   := -march=core-avx2
COMPILER = s.cc
FCOMPILER = s.f90
endif

ifneq (,$(findstring gfortran,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C   := -mavx2 -mfma
COMPILER = s.cc
FCOMPILER = s.f90
endif

ifneq (,$(findstring gnu,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C   := -mavx2 -mfma
COMPILER = gcc
FCOMPILER = gfortran
endif

ifneq (,$(findstring llvm,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C   := -mavx2 -mfma
COMPILER = clang
FCOMPILER = flang
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

intrp_bicub_yx.inc: intrp_bicub_yx.c
	grep -w InTf intrp_bicub_yx.c >intrp_bicub_yx.inc

intrp_bicub_yx_f.F90:
	ln -sf  intrp_bicub_yx.c intrp_bicub_yx_f.F90

intrp_bicub_yx.o: intrp_bicub_yx.c
	$(COMPILER) -O2 -DTIMING -c ${EXTRA_C} intrp_bicub_yx.c

intrp_bicub_yx: intrp_bicub_yx.c intrp_bicub_yx_f.F90 intrp_bicub_yx.inc
	$(COMPILER) -O2 -DTIMING -c ${EXTRA_C} intrp_bicub_yx.c
	$(FCOMPILER) ${EXTRA_FTN} -DF_TEST -o intrp_bicub_yx intrp_bicub_yx_f.F90 intrp_bicub_yx.o
	rm -f intrp_bicub_yx.o
	$(COMPILER) -O3 -DTIMING -c intrp_bicub_yx.c
	$(FCOMPILER) ${EXTRA_FTN} -DF_TEST -o intrp_bicub_yx0 intrp_bicub_yx_f.F90 intrp_bicub_yx.o
	rm -f intrp_bicub_yx.o

test_search_list: search_list.c test_search_list.c
	$(COMPILER) -O3 -DSELF_TEST -c ${EXTRA_C} search_list.c
	$(COMPILER) -O3 ${EXTRA_C} test_search_list.c search_list.o -o test_search_list
	rm -f search_list.o

clean:
	rm -f intrp_bicub_yx0 intrp_bicub_yx intrp_bicub_yx.o intrp_bicub_yx.s a.out fp_x87_xmm_control.o intrp_bicub_yx.inc
	rm -f test_search_list search_list.o 
