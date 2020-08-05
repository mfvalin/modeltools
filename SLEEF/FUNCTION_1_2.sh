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
LaSt=""
for FN in $* ; do
FuNcTiOn=${FN%??}
PrEcIsIoN=${FN#${FuNcTiOn}}
PoStFiX=""
((PrEcIsIoN >  10)) && PoStFiX="${PrEcIsIoN}"

# ${Vlibm} function calling libm (only done if precision <= 10)
((PrEcIsIoN <= 10)) && cat <<EOT

void ${Vlibm}${FuNcTiOn}_f(float *f1, float *f2, float *r, int n){
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}f(f1[i], f2[i]);
}

void ${Vlibm}${FuNcTiOn}_d(double *f1, double *f2, double *r, int n){
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}(f1[i], f2[i]);
}
EOT

# ${Vlibsleef} version calling libsleef, with specified precision
# ${Vfortran}  Fortran interface to ${Vlibsleef} version above
# if __SLEEF_H__ is not defined, the ${Vlibsleef} version calls the ${Vlibm} function
cat <<EOT
// interface ${Vfortran}${FuNcTiOn}${PoStFiX}
//   subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}(f1, f2, r, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_f${PoStFiX}')
//     import :: C_INT, C_FLOAT
//     integer(C_INT), intent(IN), value :: n
//     real(C_FLOAT), dimension(n), intent(IN)  :: f1, f2
//     real(C_FLOAT), dimension(n), intent(OUT) :: r
//   end subroutine ${Vfortran}${FuNcTiOn}_f${PoStFiX}
//   subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}(f1, f2, r, n) BIND(C,name='${Vlibsleef}${FuNcTiOn}_d${PoStFiX}')
//     import :: C_INT, C_DOUBLE
//     integer(C_INT), intent(IN), value :: n
//     real(C_DOUBLE), dimension(n), intent(IN)  :: f1, f2
//     real(C_DOUBLE), dimension(n), intent(OUT) :: r
//   end subroutine ${Vfortran}${FuNcTiOn}_d${PoStFiX}
// end interface

void ${Vlibsleef}${FuNcTiOn}_f${PoStFiX}(float *f1, float *f2, float *r, int n){
#if defined(__SLEEF_H__)
  int i;
#if defined(__x86_64__)
#if defined(__AVX2__)
  __m256 src1, src2, dst;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values
    src1 = _mm256_loadu_ps(f1+i) ;
    src2 = _mm256_loadu_ps(f2+i) ;
    dst  = Sleef_finz_${FuNcTiOn}f8_u${PrEcIsIoN}avx2(src1, src2) ;
    _mm256_storeu_ps(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_finz_${FuNcTiOn}f1_u${PrEcIsIoN}purecfma(f1[i], f2[i]) ;
    i++;
  }
#else
  __m128 src1, src2, dst;

  for(i=0 ; i<n-3 ; i+=4){  // blocks of 4 values
    src1 = _mm_loadu_ps(f1+i) ;
    src2 = _mm_loadu_ps(f2+i) ;
    dst  = Sleef_cinz_${FuNcTiOn}f4_u${PrEcIsIoN}sse2(src1, src2) ;
    _mm_storeu_ps(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_cinz_${FuNcTiOn}f1_u${PrEcIsIoN}purec(f1[i], f2[i]) ;
    i++;
  }
#endif
#endif
#else
  ${Vlibm}${FuNcTiOn}_f(f1, f2, r, n);
#endif
}

void ${Vlibsleef}${FuNcTiOn}_d${PoStFiX}(double *f1, double *f2, double *r, int n){
#if defined(__SLEEF_H__)
  int i;
#if defined(__x86_64__)
#if defined(__AVX2__)
  __m256d src1, src2, dst;

  for(i=0 ; i<n-3 ; i+=4){  // blocks of 4 values
    src1 = _mm256_loadu_pd(f1+i) ;
    src2 = _mm256_loadu_pd(f1+i) ;
    dst  = Sleef_finz_${FuNcTiOn}d4_u${PrEcIsIoN}avx2(src1, src2) ;
    _mm256_storeu_pd(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_finz_${FuNcTiOn}d1_u${PrEcIsIoN}purecfma(f1[i], f2[i]) ;
    i++;
  }
#else
  __m128d src1, src2, dst;

  for(i=0 ; i<n-1 ; i+=2){  // blocks of 2 values
    src1 = _mm_loadu_pd(f1+i) ;
    src2 = _mm_loadu_pd(f1+i) ;
    dst  = Sleef_cinz_${FuNcTiOn}d2_u${PrEcIsIoN}sse2(src1, src2) ;
    _mm_storeu_pd(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_cinz_${FuNcTiOn}d1_u${PrEcIsIoN}purec(f1[i], f2[i]) ;
    i++;
  }
#endif
#endif
#else
  ${Vlibm}${FuNcTiOn}_d(f1, f2, r, n);
#endif
}
EOT
done
