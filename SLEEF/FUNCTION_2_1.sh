#!/bin/bash
# Copyright Recherche en Prevision Numerique and contributors 2020.
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
# wrappers for libsleef https://sleef.org
#
# M.Valin 2020/07/31
#
# functions with one input and two outputs
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
void ${Vlibm}${FuNcTiOn}_f(float *f, float *r1, float *r2, int n){   // 32 bit floating point
  int i;
  for(i=0 ; i<n ; i++) ${FuNcTiOn}f(f[i], r1+i, r2+i);
}
//===========================================================================================
void ${Vlibm}${FuNcTiOn}_d(double *f, double *r1, double *r2, int n){   // 64 bit floating point
  int i;
  for(i=0 ; i<n ; i++) ${FuNcTiOn}(f[i], r1+i, r2+i);
}
//===========================================================================================
EOT

# ${Vlibsleef} version calling libsleef, with specified precision
# ${Vfortran}  Fortran interface to ${Vlibsleef} version above
# if __SLEEF_H__ is not defined, the ${Vlibsleef} version calls the ${Vlibm} function
cat <<EOT
// interface ${Vfortran}${FuNcTiOn}${PoStFiX}               ! generic interface
//   subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}(f, r1, r2, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_f${PoStFiX}')
//     import :: C_INT, C_FLOAT
//     integer(C_INT), intent(IN), value :: n
//     real(C_FLOAT), dimension(n), intent(IN)  :: f        ! 32 bit floating point input
//     real(C_FLOAT), dimension(n), intent(OUT) :: r1, r2   ! 32 bit floating point results
//   end subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}
//   subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}(f, r1, r2, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_d${PoStFiX}')
//     import :: C_INT, C_DOUBLE
//     integer(C_INT), intent(IN), value :: n
//     real(C_DOUBLE), dimension(n), intent(IN)  :: f        ! 64 bit floating point input
//     real(C_DOUBLE), dimension(n), intent(OUT) :: r1, r2   ! 64 bit floating point results
//   end subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}
// end interface
//===========================================================================================
void ${Vlibsleef}${FuNcTiOn}_f${PoStFiX}(float *f, float *r1, float *r2, int n){
  int i = 0;
  Sleef_float_2  rst;
#if defined(__SLEEF_H__) && defined(__x86_64__)
#if defined(__x86_64__)

  Sleef___m128_2 dst128;

#if defined(__AVX2__)
  Sleef___m256_2 dst256;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values if AVX2 available
    dst256 = Sleef_finz_${FuNcTiOn}f8_u${PrEcIsIoN}avx2(_mm256_loadu_ps(f+i)) ;
    _mm256_storeu_ps(r1+i, dst256.x) ;
    _mm256_storeu_ps(r2+i, dst256.y) ;
  }
  for(i=0 ; i<n-3 ; i+=4){  // blocks of 4 values if AVX2 available
    dst128 = Sleef_finz_${FuNcTiOn}f4_u${PrEcIsIoN}avx2128(_mm_loadu_ps(f+i)) ;
    _mm_storeu_ps(r1+i, dst128.x) ;
    _mm_storeu_ps(r2+i, dst128.y) ;
  }

  for( ; i<n ; i++){     // one value at a time for leftovers
    rst   = Sleef_finz_${FuNcTiOn}f1_u${PrEcIsIoN}purecfma(f[i]) ;
    r1[i] = rst.x ; r2[i] = rst.y ;
    i++;
  }
#else     // __AVX2__

  for(i=0 ; i<n-3 ; i+=4){  // blocks of 4 values if SSE2/SSE4.2 available
#if defined(__SSE4_2__)
    dst128 = Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse4(_mm_loadu_ps(f+i)) ;
#else     // __SSE4_2__
    dst128 = Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse2(_mm_loadu_ps(f+i)) ;
#endif    // __SSE4_2__
    _mm_storeu_ps(r1+i, dst128.x) ;
    _mm_storeu_ps(r2+i, dst128.y) ;
  }

  for( ; i<n ; i++){     // one value at a time for leftovers
    rst   = Sleef_cinz_${FuNcTiOn}f1_u${PrEcIsIoN}purec(f[i]) ;
    r1[i] = rst.x ; r2[i] = rst.y ;
    i++;
  }
#endif    // __AVX2__

#else     // __x86_64__
  for( ; i<n ; i++){     // one value at a time for leftovers
    rst = Sleef_${FuNcTiOn}f_u${PrEcIsIoN}(f[i]) ;                  // sleef other architectures
    r1[i] = rst.x ; r2[i] = rst.y ;
  }
#endif    // __x86_64__

#else     // __SLEEF_H__
  ${Vlibm}${FuNcTiOn}_f(f, r1, r2, n);   // not using SLEEF, call libm function directly
#endif    // __SLEEF_H__
}
//===========================================================================================
void ${Vlibsleef}${FuNcTiOn}_d${PoStFiX}(double *f, double *r1, double *r2, int n){
  int i = 0;
  Sleef_double_2  rst;
#if defined(__SLEEF_H__) && defined(__x86_64__)
#if defined(__x86_64__)

  Sleef___m128d_2 dst128;

#if defined(__AVX2__)
  Sleef___m256d_2 dst256;

  for( ; i<n-3 ; i+=4){  // blocks of 4 values if AVX2 available
    dst256 = Sleef_finz_${FuNcTiOn}d4_u${PrEcIsIoN}avx2(_mm256_loadu_pd(f+i)) ;
    _mm256_storeu_pd(r1+i, dst256.x) ;
    _mm256_storeu_pd(r2+i, dst256.y) ;
  }
  for( ; i<n-1 ; i+=2){  // blocks of 2 values if AVX2 available
    dst128 = Sleef_finz_${FuNcTiOn}d2_u${PrEcIsIoN}avx2128(_mm_loadu_pd(f+i)) ;
    _mm_storeu_pd(r1+i, dst128.x) ;
    _mm_storeu_pd(r2+i, dst128.y) ;
  }

  for( ; i<n ; i++){     // one value at a time for leftovers
    rst   = Sleef_finz_${FuNcTiOn}d1_u${PrEcIsIoN}purecfma(f[i]) ;
    r1[i] = rst.x ; r2[i] = rst.y ;
  }

#else     // __AVX2__

  for( ; i<n-1 ; i+=2){  // blocks of 2 values if SSE2/SSE4.2 available
#if defined(__SSE4_2__)
    dst128 = Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse4(_mm_loadu_pd(f+i)) ;
#else     // __SSE4_2__
    dst128 = Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse2(_mm_loadu_pd(f+i)) ;
#endif    // __SSE4_2__
    _mm_storeu_pd(r1+i, dst128.x) ;
    _mm_storeu_pd(r2+i, dst128.y) ;
  }

  for( ; i<n ; i++){     // one value at a time for leftovers
    rst   = Sleef_cinz_${FuNcTiOn}d1_u${PrEcIsIoN}purec(f[i]) ;
    r1[i] = rst.x ; r2[i] = rst.y ;
  }

#endif    // __AVX2__

#else    // __x86_64__
  for( ; i<n ; i++){     // one value at a time for leftovers
    rst = Sleef_${FuNcTiOn}_u${PrEcIsIoN}(f[i]) ;                  // sleef other architectures
    r1[i] = rst.x ; r2[i] = rst.y ;
  }
#endif    // __x86_64__

#else     // __SLEEF_H__
  ${Vlibm}${FuNcTiOn}_d(f, r1, r2, n);   // not using SLEEF, call libm function directly
#endif    // __SLEEF_H__
}
//===========================================================================================
EOT
done
