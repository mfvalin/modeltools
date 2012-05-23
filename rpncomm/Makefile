include RPN_COMM_version.inc
$(info rpn_comm VERSION IS $(RPN_COMM_version_s))

LIB = rpn_comm
CLEAN = rpn_comm_fortran_stubs.f rpn_comm_c_stubs.c $(STUB_LIBRARY) $(LIBRARY)

include Makefile.common

TESTS = TEST_001 TEST_002
FOBJECTS = $(subst .f,.o,$(wildcard RPN_COMM*.f))
COBJECTS = $(subst .c,.o,$(wildcard RPN_COMM*.c))
FMODULES = RPN_COMM_mod.o
LIBRARY = $(LIBDIR)/lib$(LIB)_$(RPN_COMM_version).a
STUB_LIBRARY = $(LIBDIR)/lib$(LIB)stubs_$(RPN_COMM_version).a


ALL: lib stublib $(TESTS)

lib: $(LIBRARY)

stublib: $(STUB_LIBRARY)

rpn_comm_fortran_stubs.f: rpn_comm_stubs.sh
	$(SHELL) rpn_comm_stubs.sh fortran

rpn_comm_c_stubs.c: rpn_comm_stubs.sh
	$(SHELL) rpn_comm_stubs.sh c

$(STUB_LIBRARY): rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	mkdir -p $(LIBDIR)
	ar rcv $(STUB_LIBRARY) rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o

$(subst .f,.o,$(wildcard RPN_COMM*.f)):	RPN_COMM_mod.o

$(LIBRARY): $(FMODULES) f77name.h $(FOBJECTS) $(COBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY) $(FOBJECTS) $(COBJECTS) $(FMODULES)

TEST_001: $(LIBRARY) TEST_001.f
	$(MPIFC) $(FCFLAGS) TEST_001.f -L$(LIBDIR) -l$(LIB)_$(RPN_COMM_version) -o TEST_001

TEST_002: $(LIBRARY) TEST_002.f90
	$(MPIF90C) $(FCFLAGS) TEST_002.f90 -L$(LIBDIR) -l$(LIB)_$(RPN_COMM_version) -o TEST_002
