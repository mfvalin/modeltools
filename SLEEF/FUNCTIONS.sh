#!/bin/bash
# Copyright Recherche en Prevision Numerique and contributors 2020.
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
# wrappers for libsleef https://sleef.org
#
# M.Valin 2020/07/31
#
set -x
#prefixes for the various interfaces
export Vlibm=${Vlibm:-Vm_}            # functions calling libm
export Vlibsleef=${Vlibsleef:-Vsl_}   # functions calling libsleef
export Vfortran=${Vfortran:-V_}       # Fortran interface to functions calling libsleef
#
# the above environment variables are used by scripts FUNCTION.sh, FUNCTION_1_2.sh, FUNCTION_2_1.sh
#
# some possibly missing prototypes for libm funcions
cat <<EOT
//!   Copyright Recherche en Prevision Numerique and contributors 2020.
//!   Distributed under the Boost Software License, Version 1.0.
//!   https://www.boost.org/LICENSE_1_0.txt
//!   wrappers for libsleef https://sleef.org
//!   M.Valin 2020/07/31
#if ! defined(_GNU_SOURCE)
  void sincos(double x, double *sin, double *cos);
  void sincosf(float x, float *sin, float *cos);
  double exp10(double x);
  float exp10f(float x);
#endif
EOT
#
# timing function and its Fortran interface (returns clock cycles)
# for architectures other than X86-64(__x86_64__) and ARM-64(__aarch64__), a dymmy function is generated
cat <<EOT
// interface 
//   function ${Vfortran}rdtsc() result(t) BIND(C,name='Vsl_rdtsc')
//   import :: C_LONG_LONG
//   integer(C_LONG_LONG) :: t
//   end function ${Vfortran}rdtsc
// end interface
static uint64_t time0 = 0;
uint64_t ${Vlibsleef}rdtsc(void) {   // serialized version on X86
#if defined(__x86_64__) || defined(USE_RDTSCP)
  uint32_t lo, hi, rcx;
  __asm__ volatile ("rdtscp" :  "=a" (lo), "=d" (hi), "=c" (rcx) );
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#endif
#if defined(__aarch64__)
  asm volatile ("isb; mrs %0, cntvct_el0" : "=r" (time0));
  return time0;
#endif
#if !defined(__x86_64__) && !defined(__aarch64__)
  return time0++;  // dummy for the time being
#endif
}

EOT
# one input, one output
./FUNCTION.sh sin10 sin35 cos10 cos35 tan10 tan35 \
              log10 log35 log210 log235 log1p10 log1010 exp10 exp210 exp1010 expm110 \
              sqrt05 sqrt35 cbrt10 cbrt35 \
              asin10 asin35 acos10 acos35 atan10 atan35 \
              sinh10 sinh35 cosh10 cosh35 tanh10 tanh35 \
              asinh10 acosh10 atanh10 \
              erf10 erfc15 tgamma10 lgamma10
# two inputs, one output
./FUNCTION_1_2.sh atan210 atan235 pow10 hypot05 hypot35
# one input, twu outputs
./FUNCTION_2_1.sh sincos10 sincos35
