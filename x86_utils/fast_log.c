/* 
 * Copyright (C) 2018  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#if !defined SELF_TEST

#include <stdint.h>
#include <immintrin.h>

typedef union{
    float    f;
    uint32_t i;
  }fi;

static float two_m_23 = 1.1920928955078125e-07f;  // 2 ** -23
static float two_p_23 = 8388608.0f;               // 2 ** +23
static float c127 = -127.0f;                      // compensation for exponent bias
// static float c127 = -126.96f;

static inline __m256 simd_inline_fake_log2(float *src, __m256 vm23, __m256 v127){
  return _mm256_fmadd_ps( _mm256_cvtepi32_ps( (__m256i) _mm256_loadu_ps(src) ), vm23, v127);
}

void simd_v_fake_log2(float *rz, float *vz, int n){  // inverse function is simd_v_fake_exp2
  int i;
  __m256  vf, vm23, v127;

  vm23 = _mm256_broadcast_ss(&two_m_23);
  v127 = _mm256_broadcast_ss(&c127);
  for(i=0; i<n-7;i++){ // blocks of 8 elements
    _mm256_storeu_ps(rz + i, simd_inline_fake_log2(vz + i , vm23, v127));
  }
  if(i < n){            // last 8 elements (may endup redoing work already done)
    _mm256_storeu_ps(rz + n - 8, simd_inline_fake_log2(vz + n - 8 , vm23, v127));
  }
}

static inline __m256i simd_inline_fake_exp2(float *src, __m256 vp23, __m256 v127){
  return _mm256_cvtps_epi32( _mm256_mul_ps( _mm256_sub_ps(  _mm256_loadu_ps(src), v127), vp23) ); // (src _ 127.0) * 2**23
}

void simd_v_fake_exp2(float *rz, float *vz, int n){    // undo what simd_v_fake_exp2 did
  int i;
  __m256  vf, vp23, v127;

  vp23 = _mm256_broadcast_ss(&two_p_23);
  v127 = _mm256_broadcast_ss(&c127);
  for(i=0; i<n-7;i++){ // blocks of 8 elements
    _mm256_storeu_ps(rz + i, (__m256) simd_inline_fake_exp2(vz + i , vp23, v127));
  }
  if(i < n){            // last 8 elements (may endup redoing work already done)
    _mm256_storeu_ps(rz + n - 8, (__m256) simd_inline_fake_exp2(vz + n - 8 , vp23, v127));
  }
}

float fake_log(float z) {  /* |x.f| must be  >= 1.0 (inverse function of fake_exp) */
  int exp;
  int sign;
  fi x;

  x.f = z;
  sign = x.i & (1 << 31);                   /* extract sign */
  exp = 0xFF & (x.i >> 23);                 /* IEEE exponent */
  exp = exp - 128;                          /* get log2(x) -1  from IEEE exponent */
  x.i = (x.i & 0x3FFFFFFF) | 0x3f800000;    /* 1.0 <= X < 2.0 (force exponent to 127) */
  x.f = x.f + exp;                          /* pseudo log  */
  x.i = x.i | sign;                         /* restore sign */
  return (x.f);
}

float fake_exp(float f){   /* |result| will be >= 1.0 (inverse function of fake_log) */
  fi x;
  int exp;
  int sign;

  x.f = f;
  sign = x.i & (1 << 31);                         /* save sign */
  x.i = x.i & 0x7FFFFFFF;                         /* strip sign */
  exp = x.f;                                      /* get log2-1  */
  x.f = x.f - exp + 1;                            /*  1.0 <= X < 2.0 */
  exp += 127;                                     /* IEEE exponent with proper offset */
  x.i = (x.i & 0x7FFFFF ) | (exp << 23) | sign;   /* rebuild float and restore sign */
  return (x.f);
}

float fake_log2(float z){
  fi xfi;
  float x;

  xfi.f = z;
  x = xfi.i;
  return (x * two_m_23) - 127.0f;
}

float fake_exp2(float z){
  fi xfi;

  xfi.f = z;
  xfi.i = (xfi.f + 127.0f) * two_p_23;
  return xfi.f;
}

float fake_log2a(float z){
  fi xfi;
  float x;

  xfi.f = z;
  x = xfi.i;
//   return (x * two_m_23.f) - 126.9569643f;
  return (x * two_m_23) - 126.96f;  // seems to give the lowest max error
}

static inline float  fasterlog2 (float x)
{
  union { float f; unsigned int i; } vx = { x };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;
//   return y - 126.94269504f;
  return y - 126.96f;
}

float fake_exp2a(float z){
  fi xfi;

  xfi.i = (z + 126.96f) * two_p_23;
  return xfi.f;
}

static inline float fastlog2 (float x)
{
  union { float f; uint32_t i; } vx = { x };
  union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
           - 1.498030302f * mx.f 
           - 1.72587999f / (0.3520887068f + mx.f);
}

float fasterlog2_(float *x) { float y=*x; return (fasterlog2(y)) ;}
float fastlog2_(float *x) { float y=*x; return (fastlog2(y)) ;}

#define cast_uint32_t (uint32_t)
static inline float
fasterpow2 (float p)
{
  float clipp = (p < -126) ? -126.0f : p;
  union { uint32_t i; float f; } v = { cast_uint32_t ( (1 << 23) * (clipp + 126.94269504f) ) };
  return v.f;
}

static inline float
fastpow2 (float p)
{
  float offset = (p < 0) ? 1.0f : 0.0f;
  float clipp = (p < -126) ? -126.0f : p;
  int w = clipp;
  float z = clipp - w + offset;
  union { uint32_t i; float f; } v = { cast_uint32_t ( (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) ) };

  return v.f;
}

float entropyfast_(int *val, int *N, int *Mask)
{
  int n=*N;
  int mask=*Mask;
  int i;
  int tab[65536];
  float k=1.0/n;
  float sum=0.0;

  for(i=0;i<65536;i++) tab[i] = 0;
  for(i=0;i<n;i++) tab[val[i]&mask] += 1;
//  for(i=0;i<65536;i++) if(tab[i]>0) { float temp=k*tab[i] ; sum -= temp * fastlog2(temp); } ;
  for(i=0;i<65536;i++)  { float temp=k*tab[i] ; sum -= temp * fastlog2(temp); } ;
  return sum;
}

float entropyfaster_(int *val, int *N, int *Mask)
{
  int n=*N;
  int mask=*Mask;
  int i;
  int tab[65536];
  float k=1.0/n;
  float sum=0.0;

  for(i=0;i<65536;i++) tab[i] = 0;
  for(i=0;i<n;i++) tab[val[i]&mask] += 1;
//  for(i=0;i<65536;i++) if(tab[i]>0) { float temp=k*tab[i] ; sum -= temp * fasterlog2(temp); } ;
  for(i=0;i<65536;i++) { float temp=k*tab[i] ; sum -= temp * fasterlog2(temp); } ;
  return sum;
}

#else

#include <stdio.h>
#include <stdint.h>
#include <math.h>

uint64_t Nanocycles(void) {
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
//   __asm__ volatile ("mfence");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#define ABS(a) ( (a>0) ? a : -(a) )
#define MIN(a,b) ( (a<b) ? a : b )
#define MAX(a,b) ( (a>b) ? a : b )
#define NINTRV 1005
#define NREP 100

float fake_log2(float x);
float fake_exp2(float x);
float fake_log2a(float x);
void simd_v_fake_log2(float *r, float *f, int n);
void simd_v_fake_exp2(float *r, float *f, int n);

int main(){
  float x = 0.000045;
  float mult = 1.0 + .7/NINTRV;
  float log_10_2;
//   float log_2_10;
  float y, z0, z1, z2, z3;
  float err1, err2, err3, errmin1, errmax1, errmin2, errmax2, errmin3, errmax3;
  int i, j;
  float f[NINTRV], r[NINTRV], e[NINTRV], g[NINTRV], h[NINTRV];
  uint64_t t0, t1;
  int t[NREP], best;

  y = 2.0;
//   y = 0.5;
  for(i=0 ; i<NINTRV; i++){
    r[i] = 0;
    f[i] = y;
    e[i] = -1.0;
    y *= mult;
//     y /= mult;
  }
  x *= .01;
  log_10_2 = log10f(2.0f);
//   log_2_10 = 1/log_10_2;
  errmax1 = 0; errmin1 = 1000000;
  errmax2 = 0; errmin2 = 1000000;
  errmax3 = 0; errmin3 = 1000000;
  err1 = 0; err2 = 0; err3 = 0;

  printf("log10f T(ns) =");
  best = 999999;
  for(j=0 ; j<NREP ; j++){
    t0 = Nanocycles();
    for(i=0 ; i<NINTRV ; i++) r[i] = log10f(f[i]);
    t1 =  Nanocycles();
    t[j] = t1 - t0;
    best = MIN(best,t[j]);
  }
  for(j=0 ; j<NREP ; j+=NREP/10){
    printf(" %5d",t[j]);
  }  
  printf(" best = %5d\n\n",best);

  printf("Scalar T(ns) =");
  best = 999999;
  for(j=0 ; j<NREP ; j++){
    t0 = Nanocycles();
    for(i=0 ; i<NINTRV ; i++) r[i] = fake_log2(f[i]);
    t1 =  Nanocycles();
    t[j] = t1 - t0;
    best = MIN(best,t[j]);
  }
  for(j=0 ; j<NREP ; j+=NREP/10){
    printf(" %5d",t[j]);
  }  
  printf(" best = %5d\n\n",best);

  printf("VecLog T(ns) =");
  best = 999999;
  for(j=0 ; j<NREP ; j++){
    t0 = Nanocycles();
    simd_v_fake_log2(r, f, NINTRV);
    t1 =  Nanocycles();
    t[j] = t1 - t0;
    best = MIN(best,t[j]);
  }
  for(j=0 ; j<NREP ; j+=NREP/10){
    printf(" %5d",t[j]);
  }  
  printf(" best = %5d\n\n",best);

  printf("VecExp T(ns) =");
  best = 999999;
  for(j=0 ; j<NREP ; j++){
    t0 = Nanocycles();
    simd_v_fake_exp2(e, r, NINTRV);
    t1 =  Nanocycles();
    t[j] = t1 - t0;
    best = MIN(best,t[j]);
  }
  for(j=0 ; j<NREP ; j+=NREP/10){
    printf(" %5d",t[j]);
  }  
  printf(" best = %5d\n\n",best);

  simd_v_fake_log2(r, f, NINTRV);
  simd_v_fake_exp2(e, r, NINTRV);  // e should be close to f

  simd_v_fake_log2(g, e, NINTRV);  // reversibility test
  simd_v_fake_exp2(h, g, NINTRV);  // h should be equal to e
  simd_v_fake_log2(g, h, NINTRV);
  simd_v_fake_exp2(h, g, NINTRV);
  simd_v_fake_log2(g, h, NINTRV);
  simd_v_fake_exp2(h, g, NINTRV);
  simd_v_fake_log2(g, h, NINTRV);
  simd_v_fake_exp2(h, g, NINTRV);
  simd_v_fake_log2(g, h, NINTRV);
  simd_v_fake_exp2(h, g, NINTRV);  // final h should be equal to e
  printf("        f       h=e(log(f))      f-h      1-h/e(log(h))  errScal        err2        errVec      log10(f)\n");
  for(i=0 ; i<NINTRV; i++){
    z1 = fake_log2(f[i]) * log_10_2;    // scalar version, straight
    z2 = fake_log2a(f[i]) * log_10_2;   // scalar version, compensated
//     z3 = fasterlog2(y) * log_10_2;
//     z3 = fake_log(y) * log_10_2;
    z3 = r[i] * log_10_2;               // vector version, compensated
    z0 = log10f(f[i]);
    err1 += ABS(1.0-z1/z0);
    err2 += ABS(1.0-z2/z0);
    err3 += ABS(1.0-z3/z0);
    errmin1 = MIN( ABS(1.0-z1/z0),errmin1);
    errmax1 = MAX( ABS(1.0-z1/z0),errmax1);
    errmin2 = MIN( ABS(1.0-z2/z0),errmin2);
    errmax2 = MAX( ABS(1.0-z2/z0),errmax2);
    errmin3 = MIN( ABS(1.0-z3/z0),errmin3);
    errmax3 = MAX( ABS(1.0-z3/z0),errmax3);
    if(i<5 || i >NINTRV-5) 
      printf(" %12.6f %12.6f %12.6f %12.6g %12.6f %12.6f %12.6f %12.6f \n",f[i],e[i],f[i]-e[i],1.0-h[i]/e[i],ABS(1.0-z1/z0),ABS(1.0-z2/z0),ABS(1.0-z3/z0),z0);
  }
  printf("MIN MAX AVG Scal = %12.5f %12.5f %12.5f\n",errmin1,errmax1,err1/NINTRV);
  printf("MIN MAX AVG Vect = %12.5f %12.5f %12.5f\n",errmin3,errmax3,err3/NINTRV);
  printf("MIN MAX AVG Scl2 = %12.5f %12.5f %12.5f\n",errmin2,errmax2,err2/NINTRV);
//   for (i = 1 ; i < 50000 ; i++ ){
//     y = fasterlog2(x);
//     j = (y*4096) + .5;
//     y = j;
//     y /= 4096;
//     z = fasterpow2(y);
//     errmax = (errmax < ( z/x)) ? ( z/x) : errmax;
//     errmin = (errmin > ( z/x)) ? ( z/x) : errmin;
//     if( i%1000 == 1) printf("%12.12E   %12.12E  %12F  %12.12E\n", x, y, 1.0/(1.0 - z/x), 1.0 - z/x);
//     x = x * mult;
//   }
//   printf("%12.12F  %12.12F %12.12F  %12.12F \n",errmax-1,1.0-errmin,1.0/(errmax-1),1.0/(1.0-errmin));
  return 0;
}
#endif