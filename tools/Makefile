
#include Makefile.common
CC = s.cc
F90C = s.f90

TESTS = test001 test002

all: $(TESTS) u.re_tag_scale
	echo making all

runtest: test_missing

test001:
	cat f_glob_fn.c | sed 's://:/;/:g' |cpp -P -DFTEST | sed 's:/;/://:g' >f_glob_fn2.f90
	$(CC) -c f_glob_fn.c
	$(F90C) f_glob_fn2.f90 f_glob_fn.o -o test001
	rm -f f_glob_fn2.f90

test002:
	cat f_glob_fn.c | sed 's://:/;/:g' |cpp -P -DFTEST -DFULL | sed 's:/;/://:g' >f_glob_fn2.f90
	$(CC) -c f_glob_fn.c
	$(F90C) f_glob_fn2.f90 f_glob_fn.o -o test002
	rm -f f_glob_fn2.f90

testfst: fstouv_byname.o test_fstouv_by_name.o
	touch file1.tfs file2.tfs file3.tfs
	mkdir -p mydir mydir2
	touch mydir/file1.fst mydir/file2.fst mydir/file3.fst mydir/file4.fst
	touch mydir/file1.Fst mydir/file2.Fst mydir/file3.Fst mydir/file4.Fst
	touch mydir2/file1.aaa mydir2/file2.aaa mydir2/file3.aaa
	touch mydir2/file1.bbb mydir2/file2.bbb mydir2/file3.bbb
	rm -f mydir2/.dir  mydir2/.Dir
	echo file1.bbb > mydir2/.Dir
	echo file2.bbb >> mydir2/.Dir
	echo file2.bbb >> mydir2/.Dir
	echo file1.aaa > mydir2/.dir
	echo file2.aaa >> mydir2/.dir
	echo file3.aaa >> mydir2/.dir
	$(CC) -c f_glob_fn.c
	$(F90C) fstouv_byname.o test_fstouv_by_name.o f_glob_fn.o -o testfst
	rm fstouv_byname.o test_fstouv_by_name.o f_glob_fn.o

runtestfst: testfst
	./testfst
	rm -rf mydir mydir2
	rm -f file1.tfs file2.tfs file3.tfs testfst

test_missing:
	$(CC) -c fstd98.c fst_missing.c DlInterface.c
	$(F90C) -o a.out bidon_001.f90 test_missing.f90 fstd98.o fst_missing.o DlInterface.o  -lrmn_015
	./a.out
	echo "missing.fst created, must be removed manually"
	rm -f a.out test_missing.o fstd98.o DlInterface.o bidon_001.o fst_missing.o

u.re_tag_scale: u.re_tag_scale.F90
	$(F90C) -o $@ $< -lrmn_015

clean:
	rm -f $(TESTS) u.re_tag_scale ./a.out *.o
