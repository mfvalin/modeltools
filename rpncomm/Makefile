include RPN_COMM_version.inc
$(info rpn_comm VERSION IS $(RPN_COMM_version_s))

LIB = rpncomm

include Makefile.common

TESTS = TEST_001 TEST_002

ALL:	lib $(TESTS)

lib:	$(LIBDIR)/lib$(LIB)_$(RPN_COMM_version).a

FOBJECTS = $(subst .f,.o,$(wildcard RPN_COMM*.f))
COBJECTS = $(subst .c,.o,$(wildcard RPN_COMM*.c))
FMODULES = RPN_COMM_mod.o
LIBRARY = $(LIBDIR)/lib$(LIB)_$(RPN_COMM_version).a

$(subst .f,.o,$(wildcard RPN_COMM*.f)):	RPN_COMM_mod.o

$(LIBRARY): $(FMODULES) f77name.h $(FOBJECTS) $(COBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY) *.o

TEST_001:	$(LIBRARY) TEST_001.f
	$(MPIF90C) $(FCFLAGS) TEST_001.f -L$(LIBDIR) -l$(LIB)_$(RPN_COMM_version) -o TEST_001

TEST_002:	$(LIBRARY) TEST_002.f90
	$(MPIF90C) $(FCFLAGS) TEST_002.f90 -L$(LIBDIR) -l$(LIB)_$(RPN_COMM_version) -o TEST_002
