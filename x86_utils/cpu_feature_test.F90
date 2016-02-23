program test_feature
  use ISO_C_BINDING
  implicit none
#define IN_FORTRAN_CODE
#include <cpu_type.h>
  if(processor_has_feature(FLAG_AVX2) .ne. 0) print *,'AVX2'
  if(processor_has_feature(FLAG_AVX) .ne. 0) print *,'AVX'
  if(processor_has_feature(FLAG_FMA) .ne. 0) print *,'FMA'
  if(processor_has_feature(FLAG_BMI) .ne. 0) print *,'BMI'
  if(processor_has_feature(FLAG_SSE4) .ne. 0) print *,'SSE4'
  if(processor_has_feature(FLAG_SSE3) .ne. 0) print *,'SSE3'
  if(processor_has_feature(FLAG_SSE2) .ne. 0) print *,'SSE2'
  stop
end
