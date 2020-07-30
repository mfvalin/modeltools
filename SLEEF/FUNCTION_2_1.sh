#!/bin/bash
LaSt=""
for FN in $* ; do
FuNcTiOn=${FN%??}
PrEcIsIoN=${FN#${FuNcTiOn}}
PoStFiX=""
((PrEcIsIoN >  10)) && PoStFiX="${PrEcIsIoN}"
((PrEcIsIoN <= 10)) && cat <<EOT

void ${FuNcTiOn}_f(float *f, float *r1, float *r2, int n){
  int i;
  for(i=0 ; i<n ; i++) ${FuNcTiOn}f(f[i], r1+i, r2+i);
}

void ${FuNcTiOn}_d(double *f, double *r1, double *r2, int n){
  int i;
  for(i=0 ; i<n ; i++) ${FuNcTiOn}(f[i], r1+i, r2+i);
}
EOT

cat <<EOT

void v_${FuNcTiOn}_f${PoStFiX}(float *f, float *r1, float *r2, int n){
  int i;
  __m256 src;
  Sleef_float_2  rst;
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
}

void v_${FuNcTiOn}_d${PoStFiX}(double *f, double *r1, double *r2, int n){
  int i;
  __m256d src;
  Sleef_double_2  rst;
  Sleef___m256d_2 dst;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values
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
}
EOT
done
