#!/bin/bash
LaSt=""
for FN in $* ; do
FuNcTiOn=${FN%??}
PrEcIsIoN=${FN#${FuNcTiOn}}
PoStFiX=""
((PrEcIsIoN >  10)) && PoStFiX="${PrEcIsIoN}"
((PrEcIsIoN <= 10)) && cat <<EOT

void ${FuNcTiOn}_f(float *f, float *r, int n){
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}f(f[i]);
}

void ${FuNcTiOn}_d(double *f, double *r, int n){
  int i;
  for(i=0 ; i<n ; i++) r[i] = ${FuNcTiOn}(f[i]);
}
EOT

cat <<EOT

void v_${FuNcTiOn}_f${PoStFiX}(float *f, float *r, int n){
  int i;
  __m256 src, dst;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values
    src = _mm256_loadu_ps(f+i) ;
    dst = Sleef_finz_${FuNcTiOn}f8_u${PrEcIsIoN}avx2(src) ;
    _mm256_storeu_ps(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_finz_${FuNcTiOn}f1_u${PrEcIsIoN}purecfma(f[i]) ;
    i++;
  }
}

void v_${FuNcTiOn}_d${PoStFiX}(double *f, double *r, int n){
  int i;
  __m256d src, dst;

  for(i=0 ; i<n-7 ; i+=8){  // blocks of 8 values
    src = _mm256_loadu_pd(f+i) ;
    dst = Sleef_finz_${FuNcTiOn}d4_u${PrEcIsIoN}avx2(src) ;
    _mm256_storeu_pd(r+i, dst) ;
  }
  while(i<n){
    r[i] = Sleef_finz_${FuNcTiOn}d1_u${PrEcIsIoN}purecfma(f[i]) ;
    i++;
  }
}
EOT
done
