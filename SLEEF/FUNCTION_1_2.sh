#!/bin/bash
# Copyright Recherche en Prevision Numerique and contributors 2020.
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
# wrappers for libsleef https://sleef.org
#
# M.Valin 2020/07/31
#
# functions with two inputs and one output
#
# environment variables ${Vlibm}, ${Vfortran}, ${Vlibsleef} are expected to be supplied by script FUNSTIONS.sh
# prefixes for the various interfaces
export Vlibm=${Vlibm:-Vm_}            # functions calling libm
export Vlibsleef=${Vlibsleef:-Vsl_}   # functions calling libsleef
export Vfortran=${Vfortran:-V_}       # Fortran interface to functions calling libsleef
#
LaSt=""
for FN in $* ; do
FuNcTiOn=${FN%??}
PrEcIsIoN=${FN#${FuNcTiOn}}
PoStFiX=""
((PrEcIsIoN >  10)) && PoStFiX="${PrEcIsIoN}"

# ${Vlibm} function calling libm (only done if precision <= 10)
((PrEcIsIoN <= 10)) && cat <<EOT
void ${Vlibm}${FuNcTiOn}_f(float *f1, float *f2, float *r, int n){   // 32 bit floating point
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}f(f1[i], f2[i]);
}
//===========================================================================================
void ${Vlibm}${FuNcTiOn}_d(double *f1, double *f2, double *r, int n){   // 64 bit floating point
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}(f1[i], f2[i]);
}
//===========================================================================================
EOT

# ${Vlibsleef} version calling libsleef, with specified precision
# ${Vfortran}  Fortran interface to ${Vlibsleef} version above
# if __SLEEF_H__ is not defined, the ${Vlibsleef} version calls the ${Vlibm} function
cat <<EOT
// interface ${Vfortran}${FuNcTiOn}${PoStFiX}               ! generic interface
//   subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}(f1, f2, r, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_f${PoStFiX}')
//     import :: C_INT, C_FLOAT
//     integer(C_INT), intent(IN), value :: n
//     real(C_FLOAT), dimension(n), intent(IN)  :: f1, f2   ! 32 bit floating point inputs
//     real(C_FLOAT), dimension(n), intent(OUT) :: r        ! 32 bit floating point result
//   end subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}
//   subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}(f1, f2, r, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_d${PoStFiX}')
//     import :: C_INT, C_DOUBLE
//     integer(C_INT), intent(IN), value :: n
//     real(C_DOUBLE), dimension(n), intent(IN)  :: f1, f2   ! 64 bit floating point inputs
//     real(C_DOUBLE), dimension(n), intent(OUT) :: r        ! 64 bit floating point result
//   end subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}
// end interface
//===========================================================================================
void ${Vlibsleef}${FuNcTiOn}_f${PoStFiX}(float *f1, float *f2, float *r, int n){
  int i = 0;
#if defined(__SLEEF_H__) && defined(__x86_64__)
#if defined(__x86_64__)

#if defined(__AVX2__)
  for( ; i<n-7 ; i+=8){  // blocks of 8 values if AVX2 available
    _mm256_storeu_ps(r+i, Sleef_finz_${FuNcTiOn}f8_u${PrEcIsIoN}avx2(_mm256_loadu_ps(f1+i), _mm256_loadu_ps(f2+i))) ;
  }
  for( ; i<n-3 ; i+=4){  // blocks of 4 values if AVX2 available
    _mm_storeu_ps(r+i, Sleef_finz_${FuNcTiOn}f4_u${PrEcIsIoN}avx2128(_mm_loadu_ps(f1+i), _mm_loadu_ps(f2+i))) ;
  }

#else    // __AVX2__

  for( ; i<n-3 ; i+=4){  // blocks of 4 values if SSE2/SSE4.2 available
#if defined(__SSE4_2__)
    _mm_storeu_ps(r+i, Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse4(_mm_loadu_ps(f1+i), _mm_loadu_ps(f2+i))) ;
#else     // __SSE4_2__
    _mm_storeu_ps(r+i, Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse2(_mm_loadu_ps(f1+i), _mm_loadu_ps(f2+i))) ;
#endif    // __SSE4_2__
  }

#endif    // __AVX2__

  for( ; i<n ; i++){     // one value at a time for leftovers or all if no SIMD available
#if defined(__AVX2__)
    r[i] = Sleef_finz_${FuNcTiOn}f1_u${PrEcIsIoN}purecfma(f1[i], f2[i]) ;
#else     // __AVX2__
    r[i] = Sleef_cinz_${FuNcTiOn}f1_u${PrEcIsIoN}purec(f1[i], f2[i]) ;
#endif    // __AVX2__
  }

#else     // __x86_64__
  for( ; i<n ; i++){     // one value at a time for leftovers or all if no SIMD available
    r[i] = Sleef_${FuNcTiOn}f_u${PrEcIsIoN}(f1[i], f2[i]) ;                // sleef other architectures
  }
#endif    // __x86_64__

#else     // __SLEEF_H__
  ${Vlibm}${FuNcTiOn}_f(f1+i, f2+i, r+i, n);   // not using SLEEF, call vector libm function directly
#endif    // __SLEEF_H__
}
//===========================================================================================
void ${Vlibsleef}${FuNcTiOn}_d${PoStFiX}(double *f1, double *f2, double *r, int n){
  int i = 0;
#if defined(__SLEEF_H__) && defined(__x86_64__)
#if defined(__x86_64__)

#if defined(__AVX2__)
  for( ; i<n-3 ; i+=4){  // blocks of 4 values if AVX2 available
    _mm256_storeu_pd(r+i, Sleef_finz_${FuNcTiOn}d4_u${PrEcIsIoN}avx2(_mm256_loadu_pd(f1+i),  _mm256_loadu_pd(f2+i))) ;
  }
  for( ; i<n-1 ; i+=2){  // blocks of 2 values if SSE2/SSE4.2 available
    _mm_storeu_pd(r+i, Sleef_finz_${FuNcTiOn}d2_u${PrEcIsIoN}avx2128(_mm_loadu_pd(f1+i), _mm_loadu_pd(f2+i))) ;
  }
  for( ; i<n ; i++){     // one value at a time for leftovers or all if no SIMD available
    r[i] = Sleef_finz_${FuNcTiOn}d1_u${PrEcIsIoN}purecfma(f1[i], f2[i]) ;
  }

#else     // __AVX2__

  for( ; i<n-1 ; i+=2){  // blocks of 2 values if SSE2/SSE4.2 available
#if defined(__SSE4_2__)
    _mm_storeu_pd(r+i, Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse4(_mm_loadu_pd(f1+i), _mm_loadu_pd(f2+i))) ;
#else     // __SSE4_2__
    _mm_storeu_pd(r+i, Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse2(_mm_loadu_pd(f1+i), _mm_loadu_pd(f2+i))) ;
#endif    // __SSE4_2__
  }
  for( ; i<n ; i++){     // one value at a time for leftovers or all if no SIMD available
    r[i] = Sleef_cinz_${FuNcTiOn}d1_u${PrEcIsIoN}purec(f1[i], f2[i]) ;
  }
#endif    // __AVX2__

#else    // __x86_64__
  for( ; i<n ; i++){     // one value at a time for leftovers or all if no SIMD available
    r[i] = Sleef_${FuNcTiOn}_u${PrEcIsIoN}(f1[i], f2[i]) ;                // sleef other architectures
  }
#endif    // __x86_64__

#else     // __SLEEF_H__
  ${Vlibm}${FuNcTiOn}_d(f1, f2, r, n);   // not using SLEEF, call vector libm function directly
#endif    // __SLEEF_H__
}
//===========================================================================================
EOT
done
