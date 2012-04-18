include RPN_COMM_version.inc
$(info rpn_comm VERSION IS $(RPN_COMM_version_s))
ifneq (,$(PLAIN))
  include Makefile.default
else
  ifdef EC_ARCH
    ifneq (,$(wildcard Makefile.$(BASE_ARCH)))
      include $(wildcard Makefile.$(BASE_ARCH))
    endif
    ifneq (,$(wildcard $(subst /,.,Makefile.$(EC_ARCH))))
      include $(wildcard $(subst /,.,Makefile.$(EC_ARCH)))
    endif
    include Makefile.ECsetup
  endif
endif
LIB = rpncomm

all: lib

clean:
	rm -f *.o

veryclean:
	rm -f *.o *.a *.mod f77name.h

lib: fortran_code c_code
	ar rcv lib$(LIB)_$(RPN_COMM_version).a *.o

modules:
	$(MPIF90C) $(FCFLAGS) -c RPN_COMM_mod.f90

fortran_code: modules
	$(MPIFC) $(FCFLAGS) -c RPN_COMM*.f

f77name:
	cd ../tools ;\
        $(CC) $(CFLAGS) -c ftn-mangling-detect.c ;\
        $(FC) $(FCFLAGS) ftn-mangling.f ftn-mangling-detect.o
	../tools/a.out >f77name.h
	cd ../tools ;\
	rm -f a.out *.o

c_code: f77name
	$(CC) $(CFLAGS) -c RPN_COMM*.c
