program test_sleef  ! timing test for vector subroutines calling libsleef wrapper
  use ISO_C_BINDING
  implicit none
#define USE_LIBM_cbrt
#define USE_LIBM_sincos
#define USE_LIBM_pow
#include <libm.hf>
  include 'sleef.inc'

  integer, parameter :: NP = 100000
  real(C_FLOAT), dimension(NP) :: f4
  real(C_FLOAT), dimension(NP) :: r4, r4b
  real(C_DOUBLE), dimension(NP) :: f8
  real(C_DOUBLE), dimension(NP) :: r8, r8b
  integer(C_LONG_LONG) :: t0, t1, t2, t3, t4, t(4)
  integer :: i, j

  do i = 1, NP
    f4(i) = 10.0 * i / NP
    f8(i) = 10.0 * i / NP
  enddo
! V_rdtsc() reads the time stamp counter (X86 only at this moment)
  print *,' function    ftn-r4  sleef-r4    ftn-r8  sleef-r8 (ns/result)'
  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = sin(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_sin(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = sin(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_sin(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'sin',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = sin(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_sin35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = sin(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_sin35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'sin35',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = cos(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_cos(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = cos(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_cos(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'cos',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = cos(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_cos35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = cos(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_cos35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'cos35',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
#if defined(USE_LIBM_sincos)
      call sincos(f4(i), r4(i), r4b(i))
#else
      r4(i)  = sin(f4(i))
      r4b(i) = cos(f4(i))
#endif
    enddo
    t1 = V_rdtsc()
    call v_sincos(f4, r4, r4b, NP)
    t2 = V_rdtsc()
    do i = 1, NP
#if defined(USE_LIBM_sincos)
      call sincos(f8(i), r8(i), r8b(i))
#else
      r8(i)  = sin(f8(i))
      r8b(i) = cos(f8(i))
#endif
    enddo
    t3 = V_rdtsc()
    call v_sincos(f8, r8, r8b, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'sincos',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
#if defined(USE_LIBM_sincos)
      call sincos(f4(i), r4(i), r4b(i))
#else
      r4(i)  = sin(f4(i))
      r4b(i) = cos(f4(i))
#endif
    enddo
    t1 = V_rdtsc()
    call v_sincos35(f4, r4, r4b, NP)
    t2 = V_rdtsc()
    do i = 1, NP
#if defined(USE_LIBM_sincos)
      call sincos(f8(i), r8(i), r8b(i))
#else
      r8(i)  = sin(f8(i))
      r8b(i) = cos(f8(i))
#endif
    enddo
    t3 = V_rdtsc()
    call v_sincos35(f8, r8, r8b, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'sincos35',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = tan(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_tan(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = tan(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_tan(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'tan',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = tan(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_tan35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = tan(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_tan35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'tan35',t/3.7/NP

  do i = 1, NP
    f4(i) = 1.0 - (2.0*i)/np
    f8(i) = 1.0 - (2.0*i)/np
  enddo

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = asin(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_asin(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = asin(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_asin(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'asin',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = asin(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_asin35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = asin(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_asin35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'asin35',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = acos(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_acos(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = acos(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_acos(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'acos',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = acos(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_acos35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = acos(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_acos35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'acos35',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = atan(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_atan(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = atan(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_atan(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'atan',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = atan(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_atan35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = atan(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_atan35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'atan35',t/3.7/NP

  do i = 1, NP
    f4(i) = (2.0*i)/np
    f8(i) = (2.0*i)/np
  enddo

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = exp(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_exp(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = exp(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_exp(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'exp',t/3.7/NP

#define pow(x,y) x**y
  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = pow(f4(i),f4(i))
    enddo
    t1 = V_rdtsc()
    call v_pow(f4, f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = pow(f8(i),f8(i))
    enddo
    t3 = V_rdtsc()
    call v_pow(f8, f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'pow',t/3.7/NP

  do i = 1, NP
    f4(i) = (2.0*i)
    f8(i) = (2.0*i)
  enddo

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = log(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_log(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = log(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_log(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'log',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = log(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_log35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = log(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_log35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'log35',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = sqrt(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_sqrt(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = sqrt(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_sqrt(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'sqrt',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = cbrt(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_cbrt(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = cbrt(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_cbrt(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'cbrt',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = cbrt(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_cbrt35(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = cbrt(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_cbrt35(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'cbrt35',t/3.7/NP


  do i = 1, NP
    f4(i) = (10.0*i)/np
    f8(i) = (10.0*i)/np
  enddo
  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = sinh(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_sinh(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = sinh(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_sinh(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'sinh',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = cosh(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_cosh(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = cosh(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_cosh(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'cosh',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = asinh(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_asinh(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = asinh(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_asinh(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'asinh',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = acosh(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_acosh(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = acosh(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_acosh(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'acosh',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = tanh(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_tanh(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = tanh(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_tanh(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'tanh',t/3.7/NP

  t = 999999999999_8
  do j = 1, 10
    t0 = V_rdtsc()
    do i = 1, NP
      r4(i) = atanh(f4(i))
    enddo
    t1 = V_rdtsc()
    call v_atanh(f4, r4, NP)
    t2 = V_rdtsc()
    do i = 1, NP
      r8(i) = atanh(f8(i))
    enddo
    t3 = V_rdtsc()
    call v_atanh(f8, r8, NP)
    t4 = V_rdtsc()
    t(1) = min(t(1),t1-t0)
    t(2) = min(t(2),t2-t1)
    t(3) = min(t(3),t3-t2)
    t(4) = min(t(4),t4-t3)
  enddo
  print 1,'atanh',t/3.7/NP


1 format(A10,4F10.2)
end program