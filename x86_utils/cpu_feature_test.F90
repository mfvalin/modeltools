program test_feature
  use ISO_C_BINDING
  implicit none
#define IN_FORTRAN_CODE
#include <cpu_type.h>
  real(C_DOUBLE) :: s1, s2
  integer(C_INT64_T) :: t1, t2
  if(cpu_has_feature(FLAG_AVX2) .ne. 0) print *,'AVX2'
  if(cpu_has_feature(FLAG_AVX) .ne. 0) print *,'AVX'
  if(cpu_has_feature(FLAG_FMA) .ne. 0) print *,'FMA'
  if(cpu_has_feature(FLAG_BMI) .ne. 0) print *,'BMI'
  if(cpu_has_feature(FLAG_SSE4) .ne. 0) print *,'SSE4'
  if(cpu_has_feature(FLAG_SSE3) .ne. 0) print *,'SSE3'
  if(cpu_has_feature(FLAG_SSE2) .ne. 0) print *,'SSE2'
  print *,'CPU apicid =',get_cpu_id()
  s1 = rdtsc_seconds()
  s2 = rdtsc_seconds()
  print *,'rdtsc_seconds overhead=',s2
  s1 = rdtscp_seconds()
  s2 = rdtscp_seconds()
  print *,'rdtscp_seconds overhead=',s2
  t1 = rdtsc()
  t2 = rdtsc()
  print *,'rdtsc overhead =',t2-t1,wall_clock_seconds(t2-t1)
  t1 = rdtscp()
  t2 = rdtscp()
  print *,'rdtscp overhead =',t2-t1,wall_clock_seconds(t2-t1)
  stop
end
