program test_feature
  use ISO_C_BINDING
  implicit none
#define IN_FORTRAN_CODE
#include <cpu_type.h>
  real(C_DOUBLE) :: s1, s2
  real :: r
  integer(C_INT64_T) :: t1, t2
  integer(C_INT) :: status
  if(cpu_has_feature(FLAG_AVX2) .ne. 0) print *,'AVX2'
  if(cpu_has_feature(FLAG_AVX) .ne. 0) print *,'AVX'
  if(cpu_has_feature(FLAG_FMA) .ne. 0) print *,'FMA'
  if(cpu_has_feature(FLAG_BMI) .ne. 0) print *,'BMI'
  if(cpu_has_feature(FLAG_SSE4) .ne. 0) print *,'SSE4'
  if(cpu_has_feature(FLAG_SSE3) .ne. 0) print *,'SSE3'
  if(cpu_has_feature(FLAG_SSE2) .ne. 0) print *,'SSE2'

  status = get_fp_status_ctl()
  print 100,'FP status = ',status
  if(iand(FP_STATUS_PE,status) .ne. 0) print *,"Precision ON"
  if(iand(FP_STATUS_UE,status) .ne. 0) print *,"Underflow ON"
  if(iand(FP_STATUS_OE,status) .ne. 0) print *,"Overflow ON"
  if(iand(FP_STATUS_ZE,status) .ne. 0) print *,"Zero divide ON"
100 format(A,z8.8)
  r = 1.0E-30; r = r*r; /* underflow */
  r = 1.0/r;            /* zero divide */
  r = 1.0E30; r = r*r;  /* overflow */
  status = get_fp_status_ctl()
  print *,"forcing underflow, zero divide, overflow"
  print 100,'FP status = ',status
  if(iand(FP_STATUS_PE,status) .ne. 0) print *,"Precision ON"
  if(iand(FP_STATUS_UE,status) .ne. 0) print *,"Underflow ON"
  if(iand(FP_STATUS_OE,status) .ne. 0) print *,"Overflow ON"
  if(iand(FP_STATUS_ZE,status) .ne. 0) print *,"Zero divide ON"

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
