
SHELL = /bin/bash

ifeq "$(BASE_ARCH)" "$(EC_ARCH)"
$(error FATAL: EC_ARCH is equal to BASE_ARCH, no compiler architecture is defined, ABORTING)
endif

ifneq (,$(findstring intel,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := -no-wrap-margin
EXTRA_C   := -march=core-avx2
CCOMPILER = s.cc
FCOMPILER = s.f90
endif

ifneq (,$(findstring gfortran,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C   := -mavx2 -mfma
CCOMPILER = s.cc
FCOMPILER = s.f90
endif

ifneq (,$(findstring gnu,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C   := -mavx2 -mfma
CCOMPILER = gcc
FCOMPILER = gfortran
endif

ifneq (,$(findstring llvm,$(EC_ARCH)))
# $(error FATAL: intel compiler  detected)
EXTRA_FTN := 
EXTRA_C   := -mavx2 -mfma
CCOMPILER = clang
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

intrp_bicub_yx.o: intrp_bicub_yx.c fast_interp.h
	$(CCOMPILER) -O3 -c ${EXTRA_C} intrp_bicub_yx.c

test_bicub: intrp_bicub_yx.c intrp_bicub_yx_f.F90 fast_interp.h
	rm -f test_bicub
	$(CCOMPILER) -O3 -DTIMING -c ${EXTRA_C} intrp_bicub_yx.c
	$(FCOMPILER) -DF_TEST intrp_bicub_yx_f.F90  intrp_bicub_yx.o -o test_bicub

intrp_bicub_yx: intrp_bicub_yx.c intrp_bicub_yx_f.F90 intrp_bicub_yx.inc
	$(CCOMPILER) -O2 -DTIMING -c ${EXTRA_C} intrp_bicub_yx.c
	$(FCOMPILER) ${EXTRA_FTN} -DF_TEST -o intrp_bicub_yx intrp_bicub_yx_f.F90 intrp_bicub_yx.o
	rm -f intrp_bicub_yx.o
	$(CCOMPILER) -O3 -DTIMING -DNO_SIMD -c intrp_bicub_yx.c
	$(FCOMPILER) ${EXTRA_FTN} -DF_TEST -o intrp_bicub_yx_nosimd intrp_bicub_yx_f.F90 intrp_bicub_yx.o
	rm -f intrp_bicub_yx.o
	if [[ -n $(RUN_TESTS) ]] ; then ./intrp_bicub_yx_nosimd | head -8 ; echo === ; ./intrp_bicub_yx | head -8 ; echo === ; fi
	if [[ -n $(RUN_TESTS) ]] ; then ./intrp_bicub_yx_nosimd | grep E- ; echo === ; ./intrp_bicub_yx | grep E- ; echo === ; fi

test_search_list: search_list.c test_search_list.c
	$(CCOMPILER) -O3 -DSELF_TEST -c ${EXTRA_C} search_list.c
	$(CCOMPILER) -O3 ${EXTRA_C} test_search_list.c search_list.o -o test_search_list
	rm -f search_list.o

tricublin_f90.inc: tricublin_beta.c
	grep InTf tricublin_beta.c | sed 's/[ ]*!InTf![ ]*//' > tricublin_f90.inc

test_tricublin.Abs: tricublin_beta_test.F90 tricublin_f90.inc tricublin_beta.c
	s.cc -c -O2 -march=core-avx2 tricublin_beta.c -DTIMING
	s.f90 $(MPI) -O2 -march=core-avx2 tricublin_beta_test.F90 tricublin_beta.o -o test_tricublin.Abs

boundary_clip.inc: boundary_clip.F90
	build_f90_interface.sh boundary_clip.F90 >boundary_clip.inc

memory_arena.inc: memory_arena.c
	grep InTf memory_arena.c | sed -e 's:^// ::' -e 's/[ ]*!InTf.*//' >memory_arena.inc

circular_buffer.inc: circular_buffer.c
	grep InTf circular_buffer.c | sed -e 's://::' -e 's/[ !]*InTf.*//' >circular_buffer.inc

test_arena.Abs: test_memory_arena.F90 memory_arena.c memory_arena.h memory_arena.inc circular_buffer.inc circular_buffer.c circular_buffer.h
	rm -f test_arena.Abs
	s.cc -c -O 2 memory_arena.c circular_buffer.c 
	s.f90 -mpi test_memory_arena.F90 memory_arena.o circular_buffer.o -o test_arena.Abs

clean:
	rm -f intrp_bicub_yx0 intrp_bicub_yx intrp_bicub_yx.o intrp_bicub_yx.s a.out fp_x87_xmm_control.o intrp_bicub_yx.inc
	rm -f test_search_list search_list.o test_bicub test_arena.Abs memory_arena.o xxhash_testc.Abs xxhash_testf.Abs xxhash.o

xxhash.inc: xxhash.h
	grep InTf xxhash.h >xxhash.inc

xxhash.o: xxhash.c xxhash.h
	$(CCOMPILER) -I. ${EXTRA_C} xxhash.c -O3 -c

xxhash_test: test_xxhash.c testf_xxhash.F90 xxhash.c xxhash.inc xxhash.h
	$(CCOMPILER) -I. ${EXTRA_C} xxhash.c -O3 -c -DTIMING
	$(CCOMPILER) test_xxhash.c xxhash.o -o xxhash_testc.Abs
	$(FCOMPILER) ${EXTRA_FTN} testf_xxhash.F90 xxhash.o -o xxhash_testf.Abs
