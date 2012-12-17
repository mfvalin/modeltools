SHELL=/bin/bash

ifneq (,$(realpath ./Makefile))

.SUFFIXES:

ifeq (,$(BASE_ARCH))
 OBJDIR := custom/unknown
else
 OBJDIR := $(EC_ARCH)
endif
$(info Phase one, OBJDIR=$(OBJDIR))

MAKETARGET=$(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile \
	SRCDIR=$(CURDIR) $(MAKECMDGOALS)

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[[ -d $@ ]] || mkdir -p $@
	+@ echo "MAKETARGET=$(MAKETARGET)"
	+@$(MAKETARGET)
	
Makefile : ;
%.mk :: ;
% :: $(OBJDIR) ; :

.PHONY: fullclean
fullclean:
	rm -rf $(EC_ARCH)

else

$(info Phase two)
VPATH=$(SRCDIR)

#$(info TOP IS $(TOP)=$(realpath $(TOP)), VPATH=$(VPATH))
#ifeq (,$(realpath $(VPATH)/Makefile.common))
#  $(error Makefile.common not found)
#endif

include $(VPATH)/RPN_COMM_version.inc
$(info this is RPN_COMM version $(RPN_COMM_version_s))

include $(VPATH)/Makefile.common
include $(VPATH)/dependencies.mk

LIB = rpn_comm
CLEAN = rpn_comm_fortran_stubs.f rpn_comm_c_stubs.c $(STUB_LIBRARY) $(LIBRARY)
TESTS = TEST_000.Abs TEST_001.Abs TEST_002.Abs TEST_003.Abs
FMODULES = RPN_COMM_mod.o
LIBRARY = $(LIBDIR)/lib$(LIB)_$(RPN_COMM_version).a
STUB_LIBRARY = $(LIBDIR)/lib$(LIB)stubs_$(RPN_COMM_version).a

ALL:  lib tests

tests:	$(TESTS)

stublib: $(STUB_LIBRARY)

lib: $(LIBRARY)

rpn_comm_fortran_stubs.f: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh fortran

rpn_comm_c_stubs.c: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh c

$(STUB_LIBRARY): rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	mkdir -p $(LIBDIR)
	ar rcv $(STUB_LIBRARY) rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o

$(LIBRARY): $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY) $(OBJECTS)

TEST_000.Abs: $(LIBRARY)
TEST_001.Abs: $(LIBRARY)
TEST_002.Abs: $(LIBRARY)
TEST_003.Abs: $(LIBRARY)

endif