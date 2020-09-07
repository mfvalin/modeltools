
# environment variable ZFP_HOME is expected to point to the install directory of the zfp compression package

CC = cc

LIBDIR = ../lib
INCDIR = ../include

F90_INCLUDES = $(INCDIR)/model_tools.inc
F90_PART_INCLUDES = model_tools.inc.000 \
                    model_tools.inc.001

LIBTOOLS = $(LIBDIR)/libtools.a

OBJECTS = zfp_compress_3d.o analyze_errors.o

all: inc test_zfp.Abs

inc : $(F90_INCLUDES)

$(F90_INCLUDES) : $(F90_PART_INCLUDES)
	cd $(INCDIR) ; \
	ar x includes.a ; \
	cat model_tools.inc.* >model_tools.inc ; \
	rm -f model_tools.inc.*

model_tools.inc.000: zfp_compress_3d.c
	grep InTf $< | sed -e 's:^//: :' -e 's/[ ]*!InTf.*//'  >$@
	ar rcv $(INCDIR)/includes.a $@
	rm $@

model_tools.inc.001: analyze_errors.c
	grep InTf $< | sed -e 's:^//: :' -e 's/[ ]*!InTf.*//'  >$@
	ar rcv $(INCDIR)/includes.a $@
	rm $@

$(LIBTOOLS): $(OBJECTS)
	ar rcv $(LIBTOOLS) $(OBJECTS)

analyze_errors.o: analyze_errors.c
	${CC} -O2 -c $<

zfp_compress_3d.o: zfp_compress_3d.c
	${CC} -O2 -I${ZFP_HOME}/include -c $<

test_zfp.Abs: $(LIBTOOLS) $(F90_INCLUDES) test_zfp.F90
	s.f90 test_zfp.F90 -I$(INCDIR) -L$(LIBDIR) -L${ZFP_HOME}/lib -ltools -lrmn -lzfp -lm -o $@

clean:
	rm -f *.o *.mod *.Abs a.out *.inc.*

veryclean: clean
	ar d $(LIBTOOLS) $(OBJECTS)
	ar d $(INCDIR)/includes.a $(F90_PART_INCLUDES)
	cd $(INCDIR) ; \
	rm -f $(F90_PART_INCLUDES) ; \
	touch $(F90_PART_INCLUDES) ; \
	ar x includes.a ; \
	cat model_tools.inc.* >model_tools.inc ; \
	rm -f $(F90_PART_INCLUDES)
