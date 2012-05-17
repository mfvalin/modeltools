include Makefile.common

TESTS = test001 test002

ALL: $(TESTS)
	echo making all

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

testfst: fstouv_byname.o test_fstouv_byname.o
	touch file1.tfs file2.tfs file3.tfs
	mkdir -p mydir mydir2
	touch mydir/file1.fst mydir/file2.fst mydir/file3.fst mydir/file4.fst
	touch mydir2/file1.aaa mydir2/file2.aaa mydir2/file3.aaa
	rm -f mydir2/.dir
	echo file1.aaa > mydir2/.dir
	echo file2.aaa >> mydir2/.dir
	echo file3.aaa >> mydir2/.dir
	$(CC) -c f_glob_fn.c
	$(F90C) fstouv_byname.o test_fstouv_byname.o f_glob_fn.o -o testfst
	rm fstouv_byname.o test_fstouv_byname.o f_glob_fn.o

runtestfst: testfst
	./testfst
	rm -rf mydir mydir2
	rm -f file1.tfs file2.tfs file3.tfs testfst
