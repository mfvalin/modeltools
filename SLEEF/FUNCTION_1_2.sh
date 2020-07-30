#!/bin/bash
# functions with two inputs and one output
LaSt=""
for FN in $* ; do
FuNcTiOn=${FN%??}
PrEcIsIoN=${FN#${FuNcTiOn}}
PoStFiX=""
((PrEcIsIoN >  10)) && PoStFiX="${PrEcIsIoN}"
((PrEcIsIoN <= 10)) && cat <<EOT

void ${FuNcTiOn}_f(float *f1, float *f2, float *r, int n){
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}f(f1[i], f2[i]);
}

void ${FuNcTiOn}_d(double *f1, double *f2, double *r, int n){
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}(f1[i], f2[i]);
}
EOT

cat <<EOT

void v_${FuNcTiOn}_f${PoStFiX}(float *f1, float *f2, float *r, int n){
  int i;
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
}

void v_${FuNcTiOn}_d${PoStFiX}(double *f1, double *f2, double *r, int n){
  int i;
  __m256d src1, src2, dst;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values
    src1 = _mm256_loadu_pd(f1+i) ;
    src2 = _mm256_loadu_pd(f1+i) ;
    dst  = Sleef_finz_${FuNcTiOn}d4_u${PrEcIsIoN}avx2(src1, src2) ;
    _mm256_storeu_pd(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_finz_${FuNcTiOn}d1_u${PrEcIsIoN}purecfma(f1[i], f2[i]) ;
    i++;
  }
}
EOT
done
