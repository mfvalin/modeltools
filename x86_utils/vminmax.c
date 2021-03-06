//  useful routines for C and FORTRAN programming
//  Copyright (C) 2020  Environnement Canada
// 
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
// 
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
// 

#include <stdint.h>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#define IS_SIGNED   0x80000000U
#define PLUS_MINUS  0x60000000U
#define HAS_PLUS    0x40000000U
#define HAS_MINUS   0x20000000U
#define HAS_ZEROS   0x10000000U

typedef union{
  uint32_t u;
  float    f;
}u_float;

typedef struct{
  uint32_t flags;   // version + packing kind + flags
  uint32_t nelem;   // size in elements of array q
  uint32_t shift;   // shift count (upper 16 bits) + nbits (lower 16 bits)
  u_float  maxza;   // largest absolute value (1 if all values are 0)
  u_float  minza;   // smallest nonzero absolute value
  float    avgz;    // rough average
  float    avga;    // rough average of absolute values
  float    maxz;    // maximum value (signed)
  float    minz;    // minimum value (signed)
  } packed_meta;
  
#if defined(__AVX2__) && defined(WITH_SIMD)
static int32_t mask16[] = { -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 };
#endif

// get from a float array of n elements :
//   min, max, average values
//   largest and smallest (nonzero) absolute values
//   smallest (possibly zero) absolute value
//   "refloated" "integer" averaged absolute value
//   HAS_PLUS, HAS_MINUS, IS_SIGNED flags
packed_meta QuantizeHeaderMinMax(float *fz, int n){
  int i = 0 ;
//   int32_t s0, s1 ;
  uint32_t minzb, minza, maxza ;
  float minz, maxz, avgz ;
  u_float zu ;
  packed_meta pm ;
  uint64_t avgu ;
  int32_t *iz = (int32_t *) fz;
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256  vsumz, vminz, vmaxz, vz, vt ;
//   __m256i v0, vmask, vs0, vs1, vmina, vminb, vmaxa, va, vsuma, vti ;
  __m256i v0, vmask, vs1, vmina, vminb, vmaxa, va, vsuma, vti ;
  __m128  tf0, tf1, tf2 ;
  __m128i ti0, ti1, ti2 ;
#endif

  avgz  = 0.0 ;
  avgu  = 0 ;
//   s0    = 0 ;                                // sign bit will be 0 only if all signs are 0
//   s1    = -1 ;                               // sign bit will be 1 only if all signs are 1
  minza = 0x7FFFFFFF  ;                      // all bits on, largest unsigned value
  minzb = 0x7FFFFFFF  ;                      // all bits on, largest unsigned value
  maxza = 1 ;                                // will stay this way only if all z values are 0
  minz  = fz[0] ;
  maxz  = fz[0] ;
  i     = 0 ;

#if defined(__AVX2__) && defined(WITH_SIMD)
  if(n >= 8) {                               // the AVX2 vectorized part will not handle n < 8
  i = (n & 7) ?  (n & 7) : 8 ;
  vt    = _mm256_loadu_ps((float *) (mask16 + (8-i)) ) ;
  vz  = _mm256_loadu_ps(fz) ;
  vsumz = _mm256_and_ps(vz, vt) ;                       // mask unwanted elements from floating sum
  v0  = _mm256_xor_si256((__m256i) vz, (__m256i) vz) ;  // set to zero
//   vs0 = v0 ;                                            // set to all zeroes  (s0)
  vs1 = _mm256_cmpeq_epi32(v0, v0) ;                    // set to all ones    (s1)
  vmask = _mm256_srli_epi32(vs1, 1) ;                   // 0x7FFFFFFF
  va    = _mm256_and_si256((__m256i) vz, vmask) ;       // absolute value in integer register
  vmaxz = vz ;
  vminz = vz ;
  vmaxa = va ;
  vminb = va ;
  vti   = _mm256_and_si256(va, (__m256i) vt) ;           // mask unwanted elements from sum
  vsuma = _mm256_blend_epi32(v0, vti, 0x55) ;            // odd elements of va zeroed
  vsuma = _mm256_add_epi64(vsuma, _mm256_blend_epi32(v0, _mm256_bsrli_epi128(vti, 4), 0x55) ) ;  // even elements of va zeroed
  vmina = vmask ;
  va    = _mm256_blendv_epi8 (va, vmina, _mm256_cmpeq_epi32(va, v0) ) ;  // use vmina where va == 0
  vmina = _mm256_min_epu32(vmina, va) ;           // smallest absolute value
  for( ; i < n ; i+=8){
    vz    = _mm256_loadu_ps(fz + i) ;               // value
    vmaxz = _mm256_max_ps(vmaxz, vz) ;
    vminz = _mm256_min_ps(vminz, vz) ;
//     vs0   = _mm256_or_si256(vs0, (__m256i) vz)  ;   // will remain 0 only if all signs are 0
//     vs1   = _mm256_and_si256(vs1, (__m256i) vz) ;   // will remain 1 only if all signs are 1
    vsumz = _mm256_add_ps(vz, vsumz) ;              // sum of floating values
    va    = _mm256_and_si256((__m256i) vz, vmask) ; // absolute value in integer register (treated as an integer)
    vmaxa = _mm256_max_epu32(vmaxa, va) ;           // largest aboslute value
    vminb = _mm256_min_epu32(vminb, va) ;           // smallest absolute value
    vsuma = _mm256_add_epi64(vsuma, _mm256_blend_epi32(v0, va, 0x55) ) ;                          // odd elements of va zeroed
    vsuma = _mm256_add_epi64(vsuma, _mm256_blend_epi32(v0, _mm256_bsrli_epi128(va, 4), 0x55) ) ;  // even elements of va zeroed
    va    = _mm256_blendv_epi8 (va, vmina, _mm256_cmpeq_epi32(va, v0) ) ;  // use vmina where va == 0
    vmina = _mm256_min_epu32(vmina, va) ;           // smallest absolute value
  }
//   ti0 = _mm256_extracti128_si256(vs0, 0) ;                      // fold vs0
//   ti0 = _mm_or_si128(ti0, _mm256_extracti128_si256(vs0, 1) ) ;  // 4 elements 0|4 1|5 2|6 3|7
//   ti0 = _mm_or_si128(ti0, _mm_srli_si128(ti0, 8) ) ;            // 2 elements 0|4|2|6 1|5|3|7  x  x
//   ti0 = _mm_or_si128(ti0, _mm_srli_si128(ti0, 4) ) ;            // 1 element
//   _mm_store_ss((float *)&s0, (__m128) ti0) ;

//   ti1 = _mm256_extracti128_si256(vs1, 0) ;                      // fold vs1
//   ti1 = _mm_and_si128(ti1, _mm256_extracti128_si256(vs1, 1) ) ; // 4 elements 0&4 1&5 2&6 3&7
//   ti1 = _mm_and_si128(ti1, _mm_srli_si128(ti1, 8) ) ;           // 2 elements 0&4&2&6 1&5&3&7  x  x
//   ti1 = _mm_and_si128(ti1, _mm_srli_si128(ti1, 4) ) ;           // 1 element
//   _mm_store_ss((float *)&s1, (__m128) ti1) ;

// note :
// _mm_store_ss / _mm_store_sd used instead of _mm_storeu_si32 / _mm_storeu_si64
// because of missing Intel intrinsics for some compilers
// another alternative would be to define _mm_storeu_si32 / _mm_storeu_si64 as macros
// #define _mm_storeu_si32(addr, what)  _mm_store_ss((float  *)addr, (__m128)  what)
// #define _mm_storeu_si64(addr, what)  _mm_store_sd((double *)addr, (__m128d) what)
//
  ti2 = _mm_add_epi64(_mm256_extracti128_si256(vsuma, 0), _mm256_extracti128_si256(vsuma, 1)) ;  // 0+2 1+3
  ti2 = _mm_add_epi64(ti2,  _mm_srli_si128(ti2, 8) ) ;          // 1 element
  _mm_store_sd((double *)&avgu, (__m128d) ti2) ;                // store folded sum of absolute values

  tf0 = _mm_add_ps(_mm256_extractf128_ps(vsumz, 0), _mm256_extractf128_ps(vsumz, 1)) ;  // 0+4 1+5 2+6 3+7
  tf0 = _mm_add_ps(tf0, _mm_shuffle_ps(tf0, tf0, 0x0E)) ;       // 4 elements 0+4+2+6 1+5+3+7  x  x
  tf0 = _mm_add_ps(tf0, _mm_shuffle_ps(tf0, tf0, 0x01)) ;       // 2 elements 0+4+2+6+1+5+3+7  x  x  x
  _mm_store_ss(&avgz, tf0) ;                                    // store folded sum

  tf1 = _mm_max_ps(_mm256_extractf128_ps(vmaxz, 0), _mm256_extractf128_ps(vmaxz, 1)) ;  // 0+4 1+5 2+6 3+7
  tf1 = _mm_max_ps(tf1, _mm_shuffle_ps(tf1, tf1, 0x0E)) ;       // 4 elements 0+4+2+6 1+5+3+7  x  x
  tf1 = _mm_max_ps(tf1, _mm_shuffle_ps(tf1, tf1, 0x01)) ;       // 2 elements 0+4+2+6+1+5+3+7  x  x  x
  _mm_store_ss(&maxz, tf1) ;                                    // store folded maximum value

  tf2 = _mm_min_ps(_mm256_extractf128_ps(vminz, 0), _mm256_extractf128_ps(vminz, 1)) ;  // 0+4 1+5 2+6 3+7
  tf2 = _mm_min_ps(tf2, _mm_shuffle_ps(tf2, tf2, 0x0E)) ;       // 4 elements 0+4+2+6 1+5+3+7  x  x
  tf2 = _mm_min_ps(tf2, _mm_shuffle_ps(tf2, tf2, 0x01)) ;       // 2 elements 0+4+2+6+1+5+3+7  x  x  x
  _mm_store_ss(&minz, tf2) ;                                    // store folded minimum value

  ti0 = _mm_max_epu32(_mm256_extracti128_si256(vmaxa, 0), _mm256_extracti128_si256(vmaxa, 1)) ; // 0,4 1,5 2,6 3,7
  ti0 = _mm_max_epu32(ti0, _mm_srli_si128(ti0, 8) ) ;            // 2 elements 0,4,2,6 1,5,3,7
  ti0 = _mm_max_epu32(ti0, _mm_srli_si128(ti0, 4) ) ;            // 1 element
  _mm_store_ss((float *)&maxza, (__m128) ti0) ;                  // stote folded largest absolute value

  ti1 = _mm_min_epu32(_mm256_extracti128_si256(vmina, 0), _mm256_extracti128_si256(vmina, 1)) ; // 0,4 1,5 2,6 3,7
  ti1 = _mm_min_epu32(ti1, _mm_srli_si128(ti1, 8) ) ;            // 2 elements 0,4,2,6 1,5,3,7
  ti1 = _mm_min_epu32(ti1, _mm_srli_si128(ti1, 4) ) ;            // 1 element
  _mm_store_ss((float *)&minza, (__m128) ti1) ;                  // stote folded smallest nonzero absolute value

  ti2 = _mm_min_epu32(_mm256_extracti128_si256(vminb, 0), _mm256_extracti128_si256(vminb, 1)) ; // 0,4 1,5 2,6 3,7
  ti2 = _mm_min_epu32(ti2, _mm_srli_si128(ti2, 8) ) ;            // 2 elements 0,4,2,6 1,5,3,7
  ti2 = _mm_min_epu32(ti2, _mm_srli_si128(ti2, 4) ) ;            // 1 element
  _mm_store_ss((float *)&minzb, (__m128) ti2) ;                  // stote folded smallest absolute value
  }
#endif
  for( ; i < n ; i++){
    zu.u = iz[i] ;
//     s0 |= zu.u ;                             // sign bit will be 0 only if all signs are 0
//     s1 &= zu.u ;                             // sign bit will be 1 only if all signs are 1
    avgz = avgz + zu.f ;
    minz = (zu.f < minz) ? zu.f : minz ;     // minimum float value
    maxz = (zu.f > maxz) ? zu.f : maxz ;     // maximum float value
    zu.u = (zu.u & 0x7FFFFFFF) ;             // abs(value)
    maxza= ( zu.u > maxza ) ? zu.u : maxza ; // largest absolute value
    minzb= ( zu.u < minzb ) ? zu.u : minzb ; // smallest (possibly zero) absolute value
    avgu = avgu + zu.u ;
    zu.u = ( zu.u == 0 ) ? minza  : zu.u  ;  // replace 0 with minza for minza computation
    minza= ( zu.u < minza ) ? zu.u : minza ; // smallest non zero absolute value
  }
  if(minza > maxza) minza = maxza - 1 ;      // if all values are 0, maxza = 1, minza = 0
  pm.flags = 0 ;
  if( minzb == 0 ) pm.flags |= HAS_ZEROS ;
  pm.minza.u = minza ;
  pm.maxza.u = maxza ;
  pm.avgz    = avgz / n  ;
  pm.minz    = minz  ;
  pm.maxz    = maxz  ;
  zu.u       = avgu / n ;
  pm.avga    = zu.f ;

  if(minz <  0) pm.flags |= HAS_MINUS ;      // there are negative numbers
  if(maxz >  0) pm.flags |= HAS_PLUS  ;      // there are positive numbers or zeroes
  if( (maxz > 0) && (minz < 0) ){            // there are positive and negative numbers
    pm.flags |= IS_SIGNED ;
  }
  return pm ;
}

static uint64_t rdtscp(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#if defined(SELF_TEST)

#include <stdio.h>
#include <stdlib.h>

#if ! defined(NP)
#define NP 400
#endif
#if ! defined(NJ)
#define NJ 25000
#endif

int main(int argc, char **argv){
  float z[NP] ;
  int i, j ;
  packed_meta pm ;
  uint64_t t0, t1;
  float tim, fact, clock;

  fact  = 1.0;
  clock = 1.0;
  if(argc > 2) {
    fact = atof(argv[1]) / atof(argv[2]) ;
    clock =  atof(argv[2]) ;
  }

//   for(i=0 ; i<NP ; i++) z[i] = i - 10.5 ;
  z[0] = -1.5 ; z[1] = 0.0  ; z[2] = 1.0 ; for(i=3 ; i<NP ; i++) z[i] = z[i-1] * (1.0f + 11.0f/NP) ;

  t0 = rdtscp() ;
  for(j=0 ; j<NJ; j++) pm = QuantizeHeaderMinMax(z, NP) ;
  t1 = rdtscp() ;
  printf(" Z[0,1,2,NP] = %12.6g %12.6g %12.6g %12.6g\n",z[0],z[1],z[2],z[NP-1]);
  if(pm.flags & HAS_ZEROS) printf(" HAS_ZEROS\n");
  printf("min = %12.6f, max = %12.7g, avg = %12.7g, avga = %12.7g, mid = %12.7g, mina = %12.6f, maxa = %12.7g, time = %lu\n", 
         pm.minz, pm.maxz, pm.avgz, pm.avga, (pm.maxz+pm.minz)*.5f, pm.minza.f, pm.maxza.f, t1-t0) ;
  tim = t1-t0;
  tim = tim / NP / NJ ;
  printf("cycles | ns per element = %8.2f | %8.2f \n",tim*fact,tim*fact/clock) ;
//   printf("cycles per element = %8.2f, nanoseconds per element = %8.2f \n",tim*fact,tim*fact/clock) ;
}
#endif
