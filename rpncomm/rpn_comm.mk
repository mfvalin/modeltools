include $(VPATH)/RPN_COMM_version.inc
$(info this is RPN_COMM version $(RPN_COMM_version_s))

# general building rules
include $(VPATH)/Makefile.common
# sources specific (mechanically generated) dependencies
include $(VPATH)/dependencies.mk

LIB      = rpn_comm
CLEAN    = rpn_comm_fortran_stubs.f rpn_comm_c_stubs.c \
           $(STUB_LIBRARY) $(LIBRARY) $(VPATH)/rpn-comm_$(RPN_COMM_version_s)_multi.ssm
CLEANDIRS= $(VPATH)/rpn-comm_$(RPN_COMM_version_s)_multi $(LIBDIR)
TESTS    = TEST_000.Abs TEST_001.Abs TEST_002.Abs TEST_003.Abs
FMODULES = RPN_COMM_mod.o
LIBNAME  = $(LIB)_$(RPN_COMM_version)
LIBRARY  = $(LIBDIR)/lib$(LIBNAME).a
STUB_LIBRARY = $(LIBDIR)/lib$(LIB)stubs_$(RPN_COMM_version).a
SOURCES  = $(INCDECKS) $(CDECKS) $(FDECKS) $(HDECKS) $(F90DECKS)

ALL:  lib tests

tests:	$(TESTS)

stublib: $(STUB_LIBRARY)

lib: $(LIBRARY)

$(VPATH)/dependencies.mk:
	(cd $(VPATH) ; find . -maxdepth 1 -type f | ../tools/mk.dependencies.pl >dependencies.mk )

ssm-package:
	rm -rf $(VPATH)/rpn-comm_${RPN_COMM_version_s}_multi
	(cd $(VPATH) ; tar zxf ssmtemplate_1.0_all.ssm ; mv ssmtemplate_1.0_all rpn-comm_$(RPN_COMM_version_s)_multi )
	(cd $(VPATH) ; cp rpn_comm_stubs.sh $(SOURCES) rpn-comm_$(RPN_COMM_version_s)_multi/src/.)
	(cd $(VPATH) ; cp Makefile Makefile.common Makefile.default Makefile.ECsetup *.mk \
	    rpn-comm_$(RPN_COMM_version_s)_multi/src/.)
	(cd $(VPATH) ; tar zcf rpn-comm_$(RPN_COMM_version_s)_multi.ssm rpn-comm_$(RPN_COMM_version_s)_multi)

rpn_comm_fortran_stubs.f: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh fortran

rpn_comm_c_stubs.c: $(VPATH)/rpn_comm_stubs.sh
	$(SHELL) $(VPATH)/rpn_comm_stubs.sh c

$(STUB_LIBRARY): rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	mkdir -p $(LIBDIR)
	ar rcv $(STUB_LIBRARY) rpn_comm_fortran_stubs.o rpn_comm_c_stubs.o
	(cd $(LIBDIR) ; ln -sf lib$(LIB)stubs_$(RPN_COMM_version).a lib$(LIB)stubs.a)

$(LIBRARY): $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar rcv $(LIBRARY) $(OBJECTS)
	(cd $(LIBDIR) ; ln -sf lib$(LIB)_$(RPN_COMM_version).a  lib$(LIB).a)

TEST_000.Abs: $(LIBRARY)
TEST_001.Abs: $(LIBRARY)
TEST_002.Abs: $(LIBRARY)
TEST_003.Abs: $(LIBRARY)

