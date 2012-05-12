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
