#include <stdio.h>
#include <stdint.h>
#if defined(NO_SSE)
#undef __SSE__
#undef __AVX2__
#endif
#if defined(__SSE__)
#include <x86intrin.h>
#endif

#define VLEN 4
#define VSHIFT 3
#define VMASK 7

// IEEE 754 components (float)
#define IEEE32_EXP(a) ( ((uint32_t) (a) >> 23) & 0xFF )
#define IEEE32_SGN(a) ( ((uint32_t) (a) >> 31) )
// rebuild IEEE 754 float from components
#define IEEE32(sign,exp,mantissa) ( (sign << 31) | ( ((exp+127) & 0xFF) << 23) | ((mantissa) & 0x7FFFFF) )
// minmax
#define MAX(a,b) ( (a > b) ? (a) : (b) )
#define MIN(a,b) ( (a < b) ? (a) : (b) )
// build immediate operand for shuffle control, 2 bits per selector (0-3)
#define SELECT(a,b,c,d) ( a + (b<<2) + (c<<4) + (d<<6) )

#if defined(__AVX2__)
// endian swap of bytes in 32 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
static uint8_t swapindex_8_32[] = {  3,  2,  1,  0,  7,  6,  5,  4, 11, 10,  9,  8, 15, 14, 13, 12};
#endif

// endian swap of 8 bit tokens in 32 bit tokens (s -> d)
void swap_8_in_32(void *s, void *d, int n){
#if defined(__AVX2__)
  __m256i ix, vs0, vs1;
  uint32_t n2;
#endif
  uint32_t i;
  uint32_t t;
  uint32_t *s0, *s1;
  uint32_t *d0, *d1;

  s0 = (uint32_t *) s;
  s1 = s0;
  d0 = (uint32_t *) d;
  d1 = d0;
  i = 0;
#if defined(__AVX2__)
  n2 = (n >> 4);        // number of 16 token chunks
  s1 = s0 + (n2 << 3);  // "halfway" point
  d1 = d0 + (n2 << 3);  // "halfway" point
  ix = _mm256_loadu_si256((__m256i const *) swapindex_8_32);
  for( ; i < n - 15 ; i += 16){                        // endian swap of 16 32 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 8 32 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs0 = _mm256_shuffle_epi8 (vs0, ix);               // shuffle (vector bswap equivalent)
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 8 32 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s0 += 8;
    s1 += 8;
    d0 += 8;
    d1 += 8;
  }
#endif
  for( ; i < n ; i++){    // loop over remainder
    t = *s1++;
    t = (t >> 24) | (t << 24) | ((t >> 8) & 0xFF00) | ((t << 8) & 0xFF0000);    // bswap
    *d1++ = t;
  }
}

// endian swap of 32 bit tokens in 64 bit tokens (s -> d)
void swap_32_in_64(void *s, void *d, int n){
#if defined(__AVX2__)
  __m256i vs0, vs1;
  uint32_t n2;
#endif
  uint32_t i ;
  uint64_t t;
  uint64_t *s0, *s1;
  uint64_t *d0, *d1;

  s0 = (uint64_t *) s;
  s1 = s0;
  d0 = (uint64_t *) d;
  d1 = d0;
  i = 0;
#if defined(__AVX2__)
  n2 = (n >> 3);        // number of 8 token chunks
  s1 = s0 + (n2 << 2);  // "halfway" point
  d1 = d0 + (n2 << 2);  // "halfway" point
  for( ; i < n - 7 ; i +=8){                           // endian swap of 8 64 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 4 64 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 4 64 bit tokens
    vs0 = _mm256_shuffle_epi32(vs0 , 0x1B);
    vs1 = _mm256_shuffle_epi32(vs1 , 0x1B);
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 4 64 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 4 64 bit tokens
    s0 += 4;
    s1 += 4;
    d0 += 4;
    d1 += 4;
  }
#endif
  for( ; i < n ; i++){    // loop over remainder
    t = *s1++;
    t = (t >> 32) | (t << 32);
    *d1++ = t;
  }
}

// get min value, max value, very rough average of 32 bit float array
// the average is only ~OK if n is a multiple of 2 * VLEN
// otherwise some points are processed (added) twice (we sum 2* np2 points)
// then quantize z_in into quant (16 bit unsigned tokens)
void FloatFastQuantizeLinear(float *z_in, uint16_t *quant, int n, int nbits, float *Max, float *Min, float *Avg, float *Rescl)
{
#if defined(__AVX2__) && defined(__FMA__)
  __m256  y0, yymin, yysca, point5;
  __m256i iy0, yymsk;
  __m128i ix0, ix1;
  int mask;
#else
  int j;
#endif
  int i, np, np2, offset, maxexp, imax, imin;
  int izero;
  float scale, fzero, fmin, quantum, sign, temp;
  float *z_in0 = z_in;
#if defined(__SSE__)
  __m128  x0, x1, x2, xxmin, xxmax, xxsum;
#else
  float tmax, tmin, tsum, vsum[VLEN], vmax[VLEN], vmin[VLEN];
#endif
  union {
    unsigned int i;
    float f;
  } xmax, xmin, xrng, xscl, q;
  float half = 0.5;
// initialize min, max, sum
#if defined(__SSE__)
  xxmin  = _mm_loadu_ps(&z_in[ 0]) ;  // min = max = 8 first elements of array
  xxmax  = xxmin;
  xxsum  = xxmin;
  xxsum  = _mm_xor_ps(xxsum,xxsum) ;  // xxsum = 0
#else
  for(j=0 ; j<VLEN ; j++){
    vmax[j] = z_in[j] ;         // min = max = VLEN first elements of array
    vmin[j] = vmax[j] ;
    vsum[j] = 0.0 ;             // xxsum = 0
  }
#endif
  np     = ((n+VMASK) >> VSHIFT ) << VSHIFT ; // next higher multiple of 2 * VLEN
  np2    = np >> 1;                           // np/2
  offset =  n - np2 ;                         // offset to "upper" part with length a multiple of VLEN
  z_in0 = z_in;
  for (i=0 ; i<np2 ; i+=VLEN){
#if defined(__SSE__)
    // SIMD code (only needs SSE instruction set)
    x0    = _mm_loadu_ps(&z_in[i]) ;          // stream from beginning
    x1    = _mm_loadu_ps(&z_in[i+offset]) ;   // stream from "midpoint"
    xxsum = _mm_add_ps(xxsum,x0);
    xxmax = _mm_max_ps(xxmax,x0);
    xxmin = _mm_min_ps(xxmin,x0);
    xxsum = _mm_add_ps(xxsum,x1);
    xxmax = _mm_max_ps(xxmax,x1);
    xxmin = _mm_min_ps(xxmin,x1);
#else
    // straight C code
    for(j=0 ; j<VLEN ; j++){                //    MAX(a,b) ( (a > b) ? (a) : (b) )
      vmax[j] = MAX(vmax[j],z_in0[j]) ;
      vmin[j] = MIN(vmin[j],z_in0[j]) ;
      vsum[j] = vsum[j] + z_in0[j] ;
      vmax[j] = MAX(vmax[j],z_in0[j+offset]) ;
      vmin[j] = MIN(vmin[j],z_in0[j+offset]) ;
      vsum[j] = vsum[j] + z_in0[j+offset] ;
    }
    z_in0 += VLEN;
#endif
  }
#if defined(__SSE__)
// final reduction operations on last 4 elements using vector instructions 
  x0 = _mm_shuffle_ps(xxmin,xxmin,SELECT(2,3,2,3)) ;  //  put elements 2 and 3 on top of 0 and 1
  x1 = _mm_shuffle_ps(xxmax,xxmax,SELECT(2,3,2,3)) ;
  x2 = _mm_shuffle_ps(xxsum,xxsum,SELECT(2,3,2,3)) ;
  xxmin = _mm_min_ps(xxmin,x0);   // f(0,2) , f(1,3), ?? , ?? ( f = min/max/add)
  xxmax = _mm_max_ps(xxmax,x1);
  xxsum = _mm_add_ps(xxsum,x2);
  x0 = _mm_shuffle_ps(xxmin,xxmin,SELECT(1,1,1,1)) ;  // put f(1,3) on top of f(0,2) (element 0 of vector)
  x1 = _mm_shuffle_ps(xxmax,xxmax,SELECT(1,1,1,1)) ;
  x2 = _mm_shuffle_ps(xxsum,xxsum,SELECT(1,1,1,1)) ;
  xxmin = _mm_min_ps(xxmin,x0);   // f( f(0,2) , f(1,3) ) , ?? , ?? , ?? ( f = min/max/add)
  xxmax = _mm_max_ps(xxmax,x1);
  xxsum = _mm_add_ps(xxsum,x2);
  _mm_store_ss(Min,xxmin);
  _mm_store_ss(Max,xxmax);
  _mm_store_ss(Avg,xxsum);
  *Avg /= np;
#else
// final wrapup, scalar loop over VLEN terms, only in non SSE case
  tmax = vmax[0] ;
  tmin = vmin[0] ;
  tsum = vsum[0] ;
  for(j=1 ; j<VLEN ; j++){
    tmax = MAX(tmax , vmax[j]) ;
    tmin = MIN(tmin , vmin[j]) ;
    tsum = tsum + vsum[j] ;    
  }
//   for(j=0;j<2;j++) {vsum[j] += vsum[j+2]; }   // force same summation order as SSE version (VLEN = 4)
//   tsum = vsum[0] + vsum[1];
  *Max  = tmax;
  *Min  = tmin;
  *Avg  = tsum / np ;
#endif
  if(*Max < 0) {       // both min and max are negative, we switch to absolute values, and swap min/max
      temp = *Min;
      *Min = - (*Max);
      *Max = - temp;
      sign = -1.0;     // everything is negative
  }else{
      sign = +1.0;
  }
// code to quantize floating poing into 16 bit tokens (nbits significant bits)
  xmax.f = *Max;
  xmin.f = *Min;
  xrng.f = xmax.f - xmin.f;
printf("range = %f, min = %f, max = %f\n",xrng.f,xmin.f,xmax.f);
//   maxexp = MAX( IEEE32_EXP(xmax.i) , IEEE32_EXP(xmin.i) );
//   maxexp = MAX( IEEE32_EXP(xrng.i) , maxexp);
  maxexp = IEEE32_EXP(xrng.i);
  xscl.i = IEEE32(0,(nbits-1) - (maxexp-127),0);
  scale = xscl.f;        // scale factor, used to convert floats to range  0 <= X < 2**nbits -1 
  xscl.i = IEEE32(0,(maxexp-127) - (nbits-1),0); // inverse scaling factor
  quantum = xscl.f ;   // quantization interval = 1.0/scale
printf("old scale = %f, old quantum = %f\n",scale,quantum);
  // make sure that quantum is not considered zero when added to min or max
  while ((quantum + xmax.f == xmax.f) && (quantum + xmin.f == xmin.f)) { quantum = quantum * 2.0 ; scale = scale / 2.0 ; }

  if (*Min > 0) {                        // all numbers positive or all numbers negative
    izero = (0 - *Min) * scale - .5 ;
  }else{                                 // max > 0, min < 0
    izero = (0 - *Min) * scale + .5 ;
  }
  q.f = quantum;
printf("new scale = %f, new quantum = %f %8x, new min = %f, izero = %d\n",scale, quantum, q.i, *Min,izero);
  if((*Min < 0)  || (*Min < xrng.f * 128.0)){
    fmin = -(izero * quantum);                                  // quantized minimum
    if( (*Min - fmin) > (.5 * quantum) ) fmin += quantum ;      // bring to closest quantum
    imax = (*Max - fmin) * scale + .5;                          // quantized max
    imin = (*Min - fmin) * scale + .5;                          // quantized min
    printf("final min = %f, deltamin = %g, coded max,min = %x %x\n",fmin,fmin-(*Min),imax,imin);
    *Min = fmin;
    izero = (0 - fmin) * scale;
    fzero = *Min + izero * quantum;
    printf("rescaled fzero = %f, izero = %d, new min = %f\n",fzero,izero,*Min);
  }
  scale   = scale * sign;      // if sign is negative, everything is negated, restored = min + Rescl * quantized
  *Rescl  = quantum * sign;
  quantum = quantum * sign;
  *Min    = *Min * sign;
  *Max    = *Max * sign;
#if defined(__AVX2__) && defined(__FMA__)
  mask = ~ (-1 << nbits) ;
  // uses AVX2, SSE4, FMA instructions
  yysca  = _mm256_set1_ps(scale) ;
  point5 = _mm256_set1_ps(half) ;
  yymin  = _mm256_set1_ps(*Min) ;
  ix0    = (__m128i) _mm_broadcast_ss( (void *) &mask ) ;
  yymsk  = (__m256i) yymin;   // useless action to silence compile time warning
  yymsk  = _mm256_inserti128_si256(yymsk,ix0,0) ;
  yymsk  = _mm256_inserti128_si256(yymsk,ix0,1) ;
#endif
  z_in0 = z_in;
  for (i=0 ; i<np2 ; i+=VLEN){
#if defined(__AVX2__) && defined(__FMA__)
    x0    = _mm_loadu_ps(&z_in0[i]) ;          // stream from beginning
    x1    = _mm_loadu_ps(&z_in0[i+offset]) ;   // stream from "midpoint"
    y0    = _mm256_insertf128_ps(y0,x0,0) ;
    y0    = _mm256_insertf128_ps(y0,x1,1) ;
    y0    = _mm256_sub_ps(y0,yymin) ;          // x - min
    y0    = _mm256_fmadd_ps(y0,yysca,point5) ; // * scale + .5
    iy0   = _mm256_cvttps_epi32(y0) ;          // to integer (lower 16 bits used only)
    iy0   = _mm256_min_epi32(iy0,yymsk) ;
    ix0   = _mm256_extracti128_si256(iy0,0) ;  // lower 128 bits
    ix1   = _mm256_extracti128_si256(iy0,1) ;  // upper 128 bits
    ix0   = _mm_packus_epi32(ix0,ix1) ;        // pack 32 -> 16 ix0 and ix1 , then stuff both back into ix0 (ix1 in upper part)
    _mm_storel_epi64((void *)&quant[i],ix0) ;              // store lower 64 bits
    _mm_storeh_pi((void *)&quant[i+offset],(__m128)ix0) ;  // store upper 64 bits
#else
    for(j=0 ; j<VLEN ; j++){
      quant[j] = (z_in0[j] - *Min) * scale + half;
      quant[j+offset] = (z_in0[j+offset] - *Min) * scale + half;
    }
    z_in0 += VLEN;
    quant += VLEN;
#endif
  }  // for (i=0 ; i<np2 ; i+=VLEN)  (both for AVX2 and normal code)
}

#define NBITS 16
// inverse of linear quantizer, short to float with bias restoration
void Short2Float(uint16_t *restrict q, float *restrict s, int n, float bias, float scale){
#if ! defined(__AVX2__)
  uint32_t i;
#endif

#if defined(__AVX2__) && defined(__FMA__)
  __m256  y0, y1, y2, y3, yscal, ybias;
  __m128i ix0, ix1, ix2, ix3;
  __m256i iy0, iy1, iy2, iy3;
  uint16_t *q1;
  float *s1;
  int n1;
  yscal  = _mm256_set1_ps(scale) ;
  ybias  = _mm256_set1_ps(bias) ;
  n1 = (n + 31) & 0xFFFFFFE0;      // multiple 0f 32 
  n1 >>= 1;                        // multiple of 16
  q1 = q + n - n1;                 // "beginning"
  s1 = s + n - n1;                 // "mid point" (may overlap tail end of "beginning")
#endif

#if defined(__AVX2__) && defined(__FMA__)
  while(n1 > 0){   // 16 elements per iteration, 16 at "beginning", 16 at "mid point"
    ix0   = _mm_loadu_si128((const __m128i *) &q[0] ) ; // "beginning"
    iy0   = _mm256_cvtepu16_epi32(ix0) ;                // unpack 16 to 32 bits
    ix1   = _mm_loadu_si128((const __m128i *) &q[8] ) ; // "beginning" + 8
    iy1   = _mm256_cvtepu16_epi32(ix1) ;                // unpack 16 to 32 bits
    ix2   = _mm_loadu_si128((const __m128i *) &q1[0]) ; // "mid point"
    iy2   = _mm256_cvtepu16_epi32(ix2) ;                // unpack 16 to 32 bits
    ix3   = _mm_loadu_si128((const __m128i *) &q1[8]) ; // "mid point" + 8
    iy3   = _mm256_cvtepu16_epi32(ix3) ;                // unpack 16 to 32 bits

    y0    = _mm256_cvtepi32_ps(iy0) ;                   // convert to float
    y1    = _mm256_cvtepi32_ps(iy1) ;
    y2    = _mm256_cvtepi32_ps(iy2) ;
    y3    = _mm256_cvtepi32_ps(iy3) ;

    y0    = _mm256_fmadd_ps(y0,yscal,ybias);            // (value * scale) + bias
    y1    = _mm256_fmadd_ps(y1,yscal,ybias);
    y2    = _mm256_fmadd_ps(y2,yscal,ybias);
    y3    = _mm256_fmadd_ps(y3,yscal,ybias);

    _mm256_storeu_ps(s   ,y0);                          // store at "beginning"
    _mm256_storeu_ps(s+ 8,y1);                          // store at "beginning" + 8
    _mm256_storeu_ps(s1  ,y2);                          // store at "mid point"
    _mm256_storeu_ps(s1+8,y3);                          // store at "mid point" + 8
    
    q  += 16;
    s  += 16;
    q1 += 16;
    s1 += 16;
    n1 -= 16;
  }
#else
  while(n >= 32){     // 32 elements per iteration
    for(i=0 ; i<8 ; i++){
      s[i]    = (q[i]    * scale) + bias;
      s[i+8]  = (q[i+8]  * scale) + bias;
      s[i+16] = (q[i+16] * scale) + bias;
      s[i+24] = (q[i+24] * scale) + bias;
    }
    q += 32;
    s += 32;
    n = n - 32;
  }
  while(n >= 8){     // 8 elements per iteration
    for(i=0 ; i<8 ; i++){
      s[i]    = (q[i]    * scale) + bias;
    }
    q += 8;
    s += 8;
    n = n - 8;
  }
  while(n-->0) {     // leftovers, one by one
    *s++ = (*q++ * scale) + bias;
  }
#endif
}

// linear quantizer, float to short with bias removal
void Float2Short(uint16_t *restrict q0, float *restrict s0, unsigned int n, float *bias, float *rscl){
  int n1, i, maxexp, izero;
#if defined(__AVX2__) && defined(__FMA__)
  float *restrict s1;
  uint16_t *restrict q1;
  __m128  x0, x1, x2, x3, xxmin, xxmax;
  __m256  y0, y1, yymin, yymax, yysca, point5;
  __m128i ix0, ix1, ix2, ix3;
  __m256i iy0, iy1;
#endif
  union {
    unsigned int i;
    float f;
  } xmin, xmax, xrng, xscl;
  float fzero, scale;

  n1 = (n + 15) & 0xFFFFFFF0 ; // multiple of 16
  n1 = n1 >> 1 ;               // multiple of 8
#if defined(__AVX2__) && defined(__FMA__)
  s1 = s0 + n - n1 ;           // "midpoint" with possible small overlap
  q1 = q0 + n - n1 ;           // "midpoint" with possible small overlap
  xxmin = _mm_loadu_ps(s0) ;
  xxmax = xxmin;
  yymin = _mm256_loadu_ps(&s0[0]);     // useless code to suppress warning
  yymin = _mm256_insertf128_ps(yymin,xxmin,0) ;
  yymin = _mm256_insertf128_ps(yymin,xxmin,1) ;
  yymax = yymin ;
  for(i=0 ; i<n1 ; i+=8){                   // loop 1, get max and min (16 items per pass)
    x0    = _mm_loadu_ps(&s0[i]) ;          // stream from beginning
    y0    = _mm256_insertf128_ps(y0,x0,0) ;
    x1    = _mm_loadu_ps(&s1[i]) ;          // stream from "midpoint"
    y0    = _mm256_insertf128_ps(y0,x1,1) ;
    yymin = _mm256_min_ps(yymin,y0);
    yymax = _mm256_max_ps(yymax,y0);
    x2    = _mm_loadu_ps(&s0[i+4]) ;        // stream from beginning + 4
    y1    = _mm256_insertf128_ps(y1,x2,0) ;
    x3    = _mm_loadu_ps(&s1[i+4]) ;        // stream from "midpoint" + 4
    y1    = _mm256_insertf128_ps(y1,x3,1) ;
    yymin = _mm256_min_ps(yymin,y1);
    yymax = _mm256_max_ps(yymax,y1);
  }
  xxmin = (__m128) _mm256_extracti128_si256((__m256i) yymin,0);
  x0    = (__m128) _mm256_extracti128_si256((__m256i) yymin,1);
  xxmin = _mm_min_ps(xxmin,x0);
  xxmax = (__m128) _mm256_extracti128_si256((__m256i) yymax,0);
  x1    = (__m128)_mm256_extracti128_si256((__m256i) yymax,1);
  xxmax = _mm_max_ps(xxmax,x1);
  // scalar wrapup for min/max from 4 way vector
  x0    = _mm_shuffle_ps(xxmin,xxmin,SELECT(2,3,2,3)) ;  //  put elements 2 and 3 on top of 0 and 1
  x1    = _mm_shuffle_ps(xxmax,xxmax,SELECT(2,3,2,3)) ;
  xxmin = _mm_min_ps(xxmin,x0);   // min(0,2) , min(1,3), 2 , 3
  xxmax = _mm_max_ps(xxmax,x1);   // max(0,2) , max(1,3), 2 , 3
  x0    = _mm_shuffle_ps(xxmin,xxmin,SELECT(1,1,1,1)) ;  // put 1 on top of 0
  x1    = _mm_shuffle_ps(xxmax,xxmax,SELECT(1,1,1,1)) ;
  xxmin = _mm_min_ss(xxmin,x0);
  xxmax = _mm_max_ss(xxmax,x1);
  _mm_store_ss(&xmin.f,xxmin);
  _mm_store_ss(&xmax.f,xxmax);
#else
  xmin.f = s0[0];
  xmax.f = s0[0];
  for(i=1 ; i<n ; i++){
    xmin.f = MIN(xmin.f,s0[i]);
    xmax.f = MAX(xmax.f,s0[i]);
  }
#endif
  xrng.f = xmax.f - xmin.f;                  // compute scaling factors from min and max
  maxexp = IEEE32_EXP(xrng.i);
  maxexp = maxexp - 127;                     // remove IEEE exponent bias
  xscl.i = IEEE32(0,(NBITS-1) - (maxexp),0); // scaling factor (float to int)
  scale  = xscl.f;
  xscl.i = IEEE32(0,(maxexp) - (NBITS-1),0); // inverse scaling factor (int to float)
  *rscl  = xscl.f;
  if(IEEE32_SGN(xmax.i) + IEEE32_SGN(xmin.i) == 1 ) {   // min < 0 and max >0
    izero = (0 - xmin.f) * scale + .5;
    fzero = xmin.f + izero * xscl.f;
    xmin.f -= fzero ;
  }
  *bias  = xmin.f;

#if defined(__AVX2__) && defined(__FMA__)
  yysca  = _mm256_set1_ps(scale) ;
  point5 = _mm256_set1_ps(0.5) ;
  yymin  = _mm256_set1_ps(xmin.f) ;
  for(i=0 ; i<n1 ; i+=8){      // quantization loop (remove bias, then quantize) (16 items per pass)
    x0    = _mm_loadu_ps(&s0[i]) ;             // stream from beginning
    y0    = _mm256_insertf128_ps(y0,x0,0) ;
    x1    = _mm_loadu_ps(&s1[i]) ;             // stream from "midpoint"
    y0    = _mm256_insertf128_ps(y0,x1,1) ;

    x2    = _mm_loadu_ps(&s0[i+4]) ;           // stream from beginning + 4
    y1    = _mm256_insertf128_ps(y1,x2,0) ;
    x3    = _mm_loadu_ps(&s1[i+4]) ;           // stream from "midpoint" + 4
    y1    = _mm256_insertf128_ps(y1,x3,1) ;

    y0    = _mm256_sub_ps(y0,yymin) ;          // x - min
    y1    = _mm256_sub_ps(y1,yymin) ;          // x - min
    y0    = _mm256_fmadd_ps(y0,yysca,point5) ; // * scale + .5
    y1    = _mm256_fmadd_ps(y1,yysca,point5) ; // * scale + .5
    iy0   = _mm256_cvttps_epi32(y0) ;          // to integer (lower 16 bits used only)
    iy1   = _mm256_cvttps_epi32(y1) ;          // to integer (lower 16 bits used only)

    ix0   = _mm256_extracti128_si256(iy0,0) ;  // lower 128 bits
    ix1   = _mm256_extracti128_si256(iy0,1) ;  // upper 128 bits
    ix0   = _mm_packus_epi32(ix0,ix1);         // 32 bit to 16 bit with unsigned saturation
    _mm_storel_pi((void *)&q0[i],(__m128) ix0) ;         // store lower 64 bits at beginning
    _mm_storeh_pi((void *)&q1[i],(__m128) ix0) ;   // store upper 64 bits at "midpoint"

    ix2   = _mm256_extracti128_si256(iy1,0) ;  // lower 128 bits
    ix3   = _mm256_extracti128_si256(iy1,1) ;  // upper 128 bits
    ix2   = _mm_packus_epi32(ix2,ix3);         // 32 bit to 16 bit with unsigned saturation
    _mm_storel_pi((void *)&q0[i+4],(__m128) ix2) ;       // store lower 64 bits at beginning + 4
    _mm_storeh_pi((void *)&q1[i+4],(__m128) ix2) ; // store upper 64 bits at "midpoint" + 4
  }
#else
  while(n >= 8){
    for(i=0 ; i<8 ; i++){
      q0[i] = (s0[i] - xmin.f) * scale;
    }
    q0 += 8;
    s0 += 8;
    n  -= 8;
  }
  while(n-->0) *q0++ = (*s0++ - xmin.f) * scale;
#endif
}

#if defined(SELF_TEST)
#define ASIZE 1026
#define RND 0.5
#define NBITS 16
int main()
{
  int i;
  float a[ASIZE+10], A[ASIZE+10];
  short unsigned int ia[ASIZE+10];
  int nbits = 16;
  float mi, ma, av, range, scal, mi0, fac, rescl, minval;
  int ima, imi, exp;
  union {
    unsigned int i;
    float f;
  } m;
//   int p2m16 = 0x37800000 ;  // 2.0 ** -16
//   int p2p16 = 0x47800000 ;  // 2.0 ** +16
//   int one   = 0x3f800000 ;  // 1.0
//   int onep  = 0x3f80000f ;  // 1.0 + epsilon
  int round = 0x000080;
  int mask  = 0x00FFFF00 ;

  mi = 1.0;
  for (i=0;i<23;i++) mi = mi * 2.0;
  m.f = mi;
  printf("2 ** 23   = %8.8x , exp = %d\n\n",m.i,m.i>>23);

  m.f = 1.000001;
//   printf("%15f = %8.8x \n",m.f,((m.i) & 0x7FFFFF) | 0x800000);
  imi = m.f * mi;
//   printf(" imi            = %8.8x \n",imi);

  mi  = 1048574 ;
//   mi = mi * 8 + 7;
  m.f = mi;
  printf("%15g = %8.8x \n",mi,((m.i) & 0x7FFFFF) | 0x800000);
  exp = m.i >> 23;
  m.i = (127 + (150 - exp)) << 23;
  fac = m.f;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  printf(" fac = %g\n",fac);
  printf(" imi            = %8.8x, %12d \n\n",imi,imi);

  mi = 65535;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);

  mi = 4094.5;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);

  mi = 254.9;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);

  mi = 255;
  imi = fac * mi + round;
  imi = (imi >> 8) << 8;
  if (imi > mask) imi = mask;
  m.f = mi;
  printf("%15f,  %8.8x\n",mi,(m.i & 0x00FFFFFF ) | 0x800000);
  printf(" imi            = %8.8x, %12d \n\n",imi, imi);
//   exit (0);
  
  a[0] = 1.0;
  for (i=1 ; i<ASIZE ; i++){
    a[i] = a[i-1]*1.01;
//     b[i] = a[i];
  }
  minval = -8.3 ;
  a[4] = minval ; a[7] = 65537.99 + a[4] ;
  printf("a[1,4,7,ASIZE-1] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  FloatFastQuantizeLinear(&a[1],&ia[1],ASIZE-2,nbits,&ma,&mi,&av,&rescl);
  for(i=1 ; i<ASIZE-1 ; i++) { A[i] = mi + ia[i] * rescl ; }
  printf("mi=%f, av=%f, ma=%f, quantum=%f\n",mi,av,ma,rescl);
  av = 0;
  for (i=1 ; i<ASIZE-1 ; i++) { av = av + a[i] ; } ;
  av = av / (ASIZE-2) ;
  printf("expected mi=%f, expected av=%f, expected ma=%f\n",minval,av,a[7]);
  printf("a[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  printf("A[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",A[1],A[4],A[7],A[ASIZE-2]);
  printf("-[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",A[1]-a[1],A[4]-a[4],A[7]-a[7],A[ASIZE-2]-a[ASIZE-2]);
  printf("ia[1,4,7,ASIZE-2] %8hu %8hu %8hu %8hu\n\n\n\n",ia[1],ia[4],ia[7],ia[ASIZE-2]);
  
  a[0] = 10.5;
  for (i=1 ; i<ASIZE ; i++){
    a[i] = a[i-1]+1.01;
  }
  minval = 8.3 ;
  a[4] = minval ; a[7] = 6553000.99 + a[4] ;
  printf("a[1,4,7,ASIZE-1] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  FloatFastQuantizeLinear(&a[1],&ia[1],ASIZE-2,nbits,&ma,&mi,&av,&rescl);
  for(i=1 ; i<ASIZE-1 ; i++) { A[i] = mi + ia[i] * rescl ; }
  printf("mi=%f, av=%f, ma=%f, quantum=%f\n",mi,av,ma,rescl);
  av = 0;
  for (i=1 ; i<ASIZE-1 ; i++) { av = av + a[i] ; } ;
  av = av / (ASIZE-2) ;
  printf("expected mi=%f, expected av=%f, expected ma=%f\n",minval,av,a[7]);
  printf("a[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  printf("A[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",A[1],A[4],A[7],A[ASIZE-2]);
  printf("-[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",A[1]-a[1],A[4]-a[4],A[7]-a[7],A[ASIZE-2]-a[ASIZE-2]);
  printf("ia[1,4,7,ASIZE-2] %8hu %8hu %8hu %8hu\n\n\n\n",ia[1],ia[4],ia[7],ia[ASIZE-2]);

  a[0] = 10.5;
  for (i=1 ; i<ASIZE ; i++){
    a[i] = a[i-1]+1.0123;
  }
  a[4] = 5.12345;
  for (i=1 ; i<ASIZE-1 ; i++) {a[i] = a[i] - 2000000.0;}
  minval = a[4];
  printf("a[1,4,7,ASIZE-1] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  FloatFastQuantizeLinear(&a[1],&ia[1],ASIZE-2,nbits,&ma,&mi,&av,&rescl);
  for(i=1 ; i<ASIZE-1 ; i++) { A[i] = mi + ia[i] * rescl ; }
  printf("mi=%f, av=%f, ma=%f, quantum=%f\n",mi,av,ma,rescl);
  av = 0;
  for (i=1 ; i<ASIZE-1 ; i++) { av = av + a[i] ; } ;
  av = av / (ASIZE-2) ;
  printf("expected mi=%f, expected av=%f, expected ma=%f\n",minval,av,a[7]);
  printf("a[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",a[1],a[4],a[7],a[ASIZE-2]);
  printf("A[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",A[1],A[4],A[7],A[ASIZE-2]);
  printf("-[1,4,7,ASIZE-2] %8f %8f %8f %8f\n",A[1]-a[1],A[4]-a[4],A[7]-a[7],A[ASIZE-2]-a[ASIZE-2]);
  printf("ia[1,4,7,ASIZE-2] %8hu %8hu %8hu %8hu\n\n\n\n",ia[1],ia[4],ia[7],ia[ASIZE-2]);

  m.f = ma;
  ima = m.i >> 23;
  m.f = mi;
  imi = m.i >> 23;
  imi &= 0xFF;
  if(ima > imi) m.i = (m.i >> (ima-imi)) << (ima-imi) ;
  mi0 = m.f;

  range = (ma - mi0) ;
  m.f = range ;                     // range
  m.i = m.i >> 23 ;                 // exponent of range
  m.i = m.i + 1 ;                   // range = next power of 2 > range
  m.i = (m.i - NBITS) -127 ;        // range / 2**NBITS, remove exponent bias
  m.i = 127 - m.i;                  // 1 / range , add exponent bias
  m.i = m.i << 23;                  // shift exponent into position
  scal =  m.f;                      // 1.0 / discretization interval
  printf("m.f = %f, m.i=%8.8x, range = %f, discr = %f, mi = %f, mi0 = %f\n",m.f,m.i,range,1.0/scal,mi,mi0);

//   m.i = p2m16;
//   scal = range * m.f;
  imi = (mi - mi0)*scal + RND;
  ima = (ma - mi0)*scal + RND;
  printf("range = %f, %f = %d, %f = %d\n",range,mi,imi,ma,ima);

  return 0;
}
#endif
