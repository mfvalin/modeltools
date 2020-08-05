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
LaSt=""
for FN in $* ; do
FuNcTiOn=${FN%??}
PrEcIsIoN=${FN#${FuNcTiOn}}
PoStFiX=""
((PrEcIsIoN >  10)) && PoStFiX="${PrEcIsIoN}"

# ${Vlibm} function calling libm (only done if precision <= 10)
((PrEcIsIoN <= 10)) && cat <<EOT

void ${Vlibm}${FuNcTiOn}_f(float *f, float *r1, float *r2, int n){
  int i;
  for(i=0 ; i<n ; i++) ${FuNcTiOn}f(f[i], r1+i, r2+i);
}

void ${Vlibm}${FuNcTiOn}_d(double *f, double *r1, double *r2, int n){
  int i;
  for(i=0 ; i<n ; i++) ${FuNcTiOn}(f[i], r1+i, r2+i);
}
EOT

# ${Vlibsleef} version calling libsleef, with specified precision
# ${Vfortran}  Fortran interface to ${Vlibsleef} version above
# if __SLEEF_H__ is not defined, the ${Vlibsleef} version calls the ${Vlibm} function
cat <<EOT
// interface ${Vfortran}${FuNcTiOn}${PoStFiX}
//   subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}(f, r1, r2, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_f${PoStFiX}')
//     import :: C_INT, C_FLOAT
//     integer(C_INT), intent(IN), value :: n
//     real(C_FLOAT), dimension(n), intent(IN)  :: f
//     real(C_FLOAT), dimension(n), intent(OUT) :: r1, r2
//   end subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}
//   subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}(f, r1, r2, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_d${PoStFiX}')
//     import :: C_INT, C_DOUBLE
//     integer(C_INT), intent(IN), value :: n
//     real(C_DOUBLE), dimension(n), intent(IN)  :: f
//     real(C_DOUBLE), dimension(n), intent(OUT) :: r1, r2
//   end subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}
// end interface

void ${Vlibsleef}${FuNcTiOn}_f${PoStFiX}(float *f, float *r1, float *r2, int n){
#if defined(__SLEEF_H__)
  int i;
  Sleef_float_2  rst;
#if defined(__x86_64__)
#if defined(__AVX2__)
  __m256 src;
  Sleef___m256_2 dst;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values
    src = _mm256_loadu_ps(f+i) ;
    dst = Sleef_finz_${FuNcTiOn}f8_u${PrEcIsIoN}avx2(src) ;
    _mm256_storeu_ps(r1+i, dst.x) ;
    _mm256_storeu_ps(r2+i, dst.y) ;
  }
  while(i<n){
    rst   = Sleef_finz_${FuNcTiOn}f1_u${PrEcIsIoN}purecfma(f[i]) ;
    r1[i] = rst.x ;
    r2[i] = rst.y ;
    i++;
  }
#else
  __m128 src;
  Sleef___m128_2 dst;

  for(i=0 ; i<n-3 ; i+=4){  // blocks of 4 values
    src = _mm_loadu_ps(f+i) ;
#if defined(__SSE4_2__)
    dst = Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse4(src) ;
#else
    dst = Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse2(src) ;
#endif
    _mm_storeu_ps(r1+i, dst.x) ;
    _mm_storeu_ps(r2+i, dst.y) ;
  }
  while(i<n){
    rst   = Sleef_cinz_${FuNcTiOn}f1_u${PrEcIsIoN}purec(f[i]) ;
    r1[i] = rst.x ;
    r2[i] = rst.y ;
    i++;
  }
#endif
#endif
#else
  ${Vlibm}${FuNcTiOn}_f(f, r1, r2, n);
#endif
}

void ${Vlibsleef}${FuNcTiOn}_d${PoStFiX}(double *f, double *r1, double *r2, int n){
#if defined(__SLEEF_H__)
  int i;
  Sleef_double_2  rst;
#if defined(__x86_64__)
#if defined(__AVX2__)
  __m256d src;
  Sleef___m256d_2 dst;

  for(i=0 ; i<n-3 ; i+=4){  // blocks of 4 values
    src = _mm256_loadu_pd(f+i) ;
    dst = Sleef_finz_${FuNcTiOn}d4_u${PrEcIsIoN}avx2(src) ;
    _mm256_storeu_pd(r1+i, dst.x) ;
    _mm256_storeu_pd(r2+i, dst.y) ;
  }
  while(i<n){
    rst   = Sleef_finz_${FuNcTiOn}d1_u${PrEcIsIoN}purecfma(f[i]) ;
    r1[i] = rst.x ;
    r2[i] = rst.y ;
    i++;
  }
#else
  __m128d src;
  Sleef___m128d_2 dst;

  for(i=0 ; i<n-1 ; i+=2){  // blocks of 2 values
    src = _mm_loadu_pd(f+i) ;
#if defined(__SSE4_2__)
    dst = Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse4(src) ;
#else
    dst = Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse2(src) ;
#endif
    _mm_storeu_pd(r1+i, dst.x) ;
    _mm_storeu_pd(r2+i, dst.y) ;
  }
  while(i<n){
    rst   = Sleef_cinz_${FuNcTiOn}d1_u${PrEcIsIoN}purec(f[i]) ;
    r1[i] = rst.x ;
    r2[i] = rst.y ;
    i++;
  }
#endif
#endif
#else
  ${Vlibm}${FuNcTiOn}_d(f, r1, r2, n);
#endif
}
EOT
done
