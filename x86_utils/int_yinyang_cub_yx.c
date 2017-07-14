static float cp133 = 0.166666666667;
static float cm133 = -0.16666666667;
static float cp5 = .5;
static float cm5 = -.5;
static float one = 1.0;
static float two = 2.0;

#include <stdint.h>
static uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}

#include <immintrin.h>

void int_yinyang_cub_yx(float *f, float *r, int ni, int ninj, int nk, int np, double x, double y){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3 ;
  __m128  fr0, fr1, fr2, fr3, ftr ;
  __m128d ft0, ft1;
  __m128  fmi, fma;
#else
  double fd0[4], fd1[4], fd2[4], fd3[0] ;
#endif
  double  wx[4], wy[4] ;
  int ni2 = ni + ni;
  int ni3 = ni2 + ni;
  int i, k;

  wx[0] = cm133*x*(x-one)*(x-two);
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;
  fwy0 = _mm256_set1_pd(wy[0]) ;
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
#endif
  for(k=0 ; k<nk ; k++){
#if defined(__AVX2__) && defined(__x86_64__)
    fr0 = _mm_loadu_ps(f) ;
    fmi = fr0 ;
    fma = fr0 ;
    fd0 = _mm256_cvtps_pd(fr0) ;
    fr1 = _mm_loadu_ps(f+ni) ;
    fmi = _mm_min_ps(fmi,fr1) ;
    fma = _mm_max_ps(fma,fr1) ;
    fd1 = _mm256_cvtps_pd(fr1) ;
    fr2 = _mm_loadu_ps(f+ni2) ;
    fmi = _mm_min_ps(fmi,fr2) ;
    fma = _mm_max_ps(fma,fr2) ;
    fd2 = _mm256_cvtps_pd(fr2) ;
    fr3 = _mm_loadu_ps(f+ni3) ;
    fmi = _mm_min_ps(fmi,fr3) ;
    fma = _mm_max_ps(fma,fr3) ;
    fd3 = _mm256_cvtps_pd(fr3) ;

    fd0 = _mm256_mul_pd(fd0,fwy0) ;
    fd0 = _mm256_fmadd_pd(fd1,fwy1,fd0) ;
    fd0 = _mm256_fmadd_pd(fd2,fwy2,fd0) ;
    fd0 = _mm256_fmadd_pd(fd3,fwy3,fd0) ;
    fd0 = _mm256_mul_pd(fd0,fwx) ;

    fr0 = _mm_permute_ps(fmi,0xEE) ;
    fmi = _mm_min_ps(fmi,fr0) ;
    fr0 = _mm_permute_ps(fmi,0x55) ;
    fmi = _mm_min_ps(fmi,fr0) ;
    fr1 = _mm_permute_ps(fma,0xEE) ;
    fma = _mm_max_ps(fma,fr1) ;
    fr1 = _mm_permute_ps(fma,0x55) ;
    fma = _mm_max_ps(fma,fr1) ;

    ft0 = _mm256_extractf128_pd(fd0,0) ;
    ft1 = _mm256_extractf128_pd(fd0,1) ;
    ft0 = _mm_add_pd(ft0,ft1) ;
    ft1 = _mm_permute_pd(ft0,0x55) ;
    ft0 = _mm_add_pd(ft0,ft1) ;
    ftr = _mm_cvtsd_ss(ftr,ft0) ;
#if defined(MONO)
    ftr = _mm_min_ss(ftr,fmi) ;
    ftr = _mm_max_ss(ftr,fma) ;
#endif
    _mm_store_ss(r,ftr) ;
#else
    for(i=0 ; i<4 ; i++){
//       fd0[i] = f[i];
//       fd1[i] = f[i+ni];
//       fd2[i] = f[i+ni2];
//       fd3[i] = f[i+ni3];
//       fd0[i] = fd0[i]*wy[0] + fd1[i]*wy[1] + fd2[i]*wy[2] + fd3[i]*wy[3];
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
//       fd0[i] = fd0[i]*wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
#endif
    f+= ninj;
    r += np;
  }
}
