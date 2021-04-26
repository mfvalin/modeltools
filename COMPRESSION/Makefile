
# environment variable ZFP_HOME is expected to point to the install directory of the zfp compression package

CC = cc
FC = gfortran
FCFLAGS = -ffree-line-length-0 -fcray-pointer -fconvert=big-endian -frecord-marker=4 -fno-second-underscore

LIBDIR = ../lib
INCDIR = ../include

F90_INCLUDES = $(INCDIR)/model_tools.inc
F90_PART_INCLUDES = model_tools.inc.000 \
                    model_tools.inc.001 \
                    model_tools.inc.002 \
                    model_tools.inc.003

LIBTOOLS = $(LIBDIR)/libtools.a

OBJECTS = zfp_compress_3d.o analyze_errors.o cdf97.o cdf53.o dwthaar.o

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

model_tools.inc.002: cdf97.c
	grep InTf $< | sed -e 's:^//: :' -e 's/[ ]*!InTf.*//'  >$@
	ar rcv $(INCDIR)/includes.a $@
	rm $@

model_tools.inc.003: cdf53.c
	grep InTf $< | sed -e 's:^//: :' -e 's/[ ]*!InTf.*//'  >$@
	ar rcv $(INCDIR)/includes.a $@
	rm $@

$(LIBTOOLS): $(OBJECTS)
	ar rcv $(LIBTOOLS) $(OBJECTS)

dwthaar.o: dwthaar.c
	${CC} -I. -O2 -c $<

cdf53.o: cdf53.c
	${CC} -I. -O2 -c $<

cdf97.o: cdf97.c
	${CC} -I. -O2 -c $<

analyze_errors.o: analyze_errors.c
	${CC} -I. -O2 -c $<

zfp_compress_3d.o: zfp_compress_3d.c
	${CC} -I. -O2 -I${ZFP_HOME}/include -c $<

test_zfp.Abs: $(LIBTOOLS) $(F90_INCLUDES) test_zfp.F90
	${FC} ${FCFLAGS} -I. test_zfp.F90 -I$(INCDIR) -L$(LIBDIR) -L${ZFP_HOME}/lib -ltools -lrmn -lzfp -lfse -lm -o $@

clean:
	rm -f *.o *.mod *.Abs a.out *.inc.* diag_file

veryclean: clean
	ar d $(LIBTOOLS) $(OBJECTS)
	ar d $(INCDIR)/includes.a $(F90_PART_INCLUDES)
	cd $(INCDIR) ; \
	rm -f $(F90_PART_INCLUDES) ; \
	touch $(F90_PART_INCLUDES) ; \
	ar x includes.a ; \
	cat model_tools.inc.* >model_tools.inc ; \
	rm -f $(F90_PART_INCLUDES)
