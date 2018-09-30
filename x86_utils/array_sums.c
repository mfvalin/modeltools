/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2018  Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <immintrin.h>
#include <stdint.h>

static int32_t mask0[] = { -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0} ;

// non compensated sum of float array, result in double
double FloatArraySum(float *f, int32_t n){
  int32_t n0 = (n & 7);  // modulo 7
  int32_t n1, n2;
  int32_t i;
  double sum1 = 0.0;
  float *f1, *f2;
#if defined(__AVX__) && defined(__x86_64__)
  __m256  v0;      // first slice (0 -> 7 elements used)
  __m256  m0;      // mask used to zero unwanted elements in first slice
  __m256d v1, v2;  // data vectors
  __m256d s1, s2;  // running sum vectors
  __m128d s0;      // final wrapup
#else
  double sum2 = 0.0;
#endif

  n1 = n - n0;    // remaining length is now a multiple of 8
  n2 = n1 >> 1;
  f1 = f + n0;
  f2 = f1 + n2;   // mid point of what is left

#if defined(__AVX__) && defined(__x86_64__)
  v0 = _mm256_loadu_ps( f );                             // get 8 elements
  m0 = _mm256_loadu_ps ( (float *) &(mask0[8-n0]) );
  v0 = _mm256_and_ps (v0, m0);                           // zero unneeded ones
  s1 = _mm256_cvtps_pd( _mm256_extractf128_ps(v0,0) );   // lower part of v0 promoted to double
  s2 = _mm256_cvtps_pd( _mm256_extractf128_ps(v0,1) );   // upper part of v0 promoted to double
  for (i = 0 ; i < n2 ; i += 4){                         // 8 elements processed per iteration, 2 streams of 4 elements
    v1 =  _mm256_cvtps_pd( _mm_loadu_ps( &(f1[i]) ) );   // from beginning
    v2 =  _mm256_cvtps_pd( _mm_loadu_ps( &(f2[i]) ) );   // from midpoint
    s1 = _mm256_add_pd(s1, v1);
    s2 = _mm256_add_pd(s2, v2);
  }
  s1 = _mm256_add_pd(s1, s2);                            // sum of partial sums ( 4 elements)
  s0 = _mm_add_pd( _mm256_extractf128_pd(s1,0) , _mm256_extractf128_pd(s1,1) );  // sum reduced to 2 elements
  _mm_store_sd(&sum1, _mm_hadd_pd(s0, s0));               // sum of the 2 elements from vector s0
#else
  for(i = 0 ; i < n0 ; i++) {  // first slice (0 -> 7 elements used)
    sum1 += f[i];
  }
  for(i = 0 ; i < n2 ; i++) {  // two add streams, one from beginning, one from midpoint
    sum1 += f1[i];
    sum2 += f2[i];
  }
  sum1 += sum2;
#endif
  return sum1;
}

// compensated sum of float array, result in double
double FloatArraySumCompensated(float *f, int32_t n){
  int32_t n0 = (n & 7);  // modulo 7
  int32_t n1, n2;
  int32_t i;
  double sum1 = 0.0;
  float *f1, *f2;
#if defined(__AVX__) && defined(__x86_64__)
  __m256  v0;      // first slice (0 -> 7 elements used)
  __m256  m0;      // mask used to zero unwanted elements in first slice
  __m256d v1, v2;  // data vectors
  __m256d s1, s2;  // running sum vectors
  __m256d c1, c2;  // running error vectors
  __m256d y1, y2;  // temporary vectors
  __m128d s0;      // final wrapup
#else
  double sum2 = 0.0;
#endif

  n1 = n - n0;    // remaining length is now a multiple of 8
  n2 = n1 >> 1;
  f1 = f + n0;
  f2 = f1 + n2;   // mid point of what is left

#if defined(__AVX__) && defined(__x86_64__)
  v0 = _mm256_loadu_ps( f );                             // get 8 elements
  m0 = _mm256_loadu_ps ( (float *) &(mask0[8-n0]) );
  v0 = _mm256_and_ps (v0, m0);                           // zero unneeded ones
  s1 = _mm256_cvtps_pd( _mm256_extractf128_ps(v0,0) );   // lower part of v0 promoted to double
  s2 = _mm256_cvtps_pd( _mm256_extractf128_ps(v0,1) );   // upper part of v0 promoted to double
  c1 =  _mm256_xor_pd(c1, c1);                           // set compensation terms to zero
  c2 =  _mm256_xor_pd(c2, c2);
  for (i = 0 ; i < n2 ; i += 4){                         // 8 elements processed per iteration, 2 streams of 4 elements
    v1 =  _mm256_cvtps_pd( _mm_loadu_ps( &(f1[i]) ) );   // from beginning  t = f
    v2 =  _mm256_cvtps_pd( _mm_loadu_ps( &(f2[i]) ) );   // from midpoint
    y1 = _mm256_sub_pd(v1,c1);                           // y = t - c
    y2 = _mm256_sub_pd(v2,c2);
    v1 = s1;                                             // t = s
    v2 = s2;
    s1 = _mm256_add_pd(s1,y1);                           // s = s + y
    s2 = _mm256_add_pd(s2,y2);
    c1 = _mm256_sub_pd(s1,v1);                           // c = (s - t)
    c2 = _mm256_sub_pd(s2,v2);
    c1 = _mm256_sub_pd(c1,y1);                           // c = (s - t) - y
    c2 = _mm256_sub_pd(c2,y2);
  }
  s1 = _mm256_add_pd(s1, s2);                            // sum of partial sums ( 4 elements)
  s0 = _mm_add_pd( _mm256_extractf128_pd(s1,0) , _mm256_extractf128_pd(s1,1) );  // sum reduced to 2 elements
  _mm_store_sd(&sum1, _mm_hadd_pd(s0, s0));               // sum of the 2 elements from vector s0
#else
  for(i = 0 ; i < n0 ; i++) {  // first slice (0 -> 7 elements used)
    sum1 += f[i];
  }
  for(i = 0 ; i < n2 ; i++) {  // two add streams, one from beginning, one from midpoint
    sum1 += f1[i];
    sum2 += f2[i];
  }
  sum1 += sum2;
#endif
  return sum1;
}

// non compensated sum of double array, result in double
double DoubleArraySum(double *f, uint32_t n){
  int32_t n0 = (n & 3);  // modulo 3
  int32_t n1, n2;
  int32_t i;
  double sum1 = 0.0;
  double *f1, *f2;
#if defined(__AVX__) && defined(__x86_64__)
  __m256d m0;      // mask used to zero unwanted elements in first slices
  __m256d v1, v2;  // data vectors
  __m256d s1, s2;  // running sum vectors
  __m128d s0;      // final wrapup
#else
  double sum2 = 0.0;
#endif
#if defined(__AVX__) && defined(__x86_64__)
  n0 = (n & 4);                      // modulo 7 larger than 4 ?
  s1 = _mm256_loadu_pd( f );
  m0 = _mm256_loadu_pd( (double *) &(mask0[8-n0*2]) );   // will be all zeroes if modulo 7 smaller than 4
  s1 = _mm256_and_pd(s1, m0);

  s2 = _mm256_loadu_pd ( f + n0 );   // if modulo 7 smaller than 4, fetch same as before
  n0 = (n & 3);                                          // modulo 3
  m0 = _mm256_loadu_pd( (double *) &(mask0[8-n0*2]) );
  s2 = _mm256_and_pd(s2, m0);

  n0 = (n & 7);
  n1 = n - n0;    // remaining length is now a multiple of 8
  n2 = n1 >> 1;
  f1 = f + n0;
  f2 = f1 + n2;   // mid point of what is left
  for (i = 0 ; i < n2 ; i += 4){         // 8 elements processed per iteration, 2 streams of 4 elements
    v1 =  _mm256_loadu_pd( &(f1[i]) );   // from beginning
    v2 =  _mm256_loadu_pd( &(f2[i]) );   // from midpoint
    s1 = _mm256_add_pd(s1, v1);
    s2 = _mm256_add_pd(s2, v2);
  }
  s1 = _mm256_add_pd(s1, s2);                            // sum of partial sums ( 4 elements)
  s0 = _mm_add_pd( _mm256_extractf128_pd(s1,0) , _mm256_extractf128_pd(s1,1) );  // sum reduced to 2 elements
  _mm_store_sd(&sum1, _mm_hadd_pd(s0, s0));               // sum of the 2 elements from vector s0
#else
  n0 = (n & 3);
  n1 = n - n0;    // remaining length is now a multiple of 4
  n2 = n1 >> 1;
  f1 = f + n0;
  f2 = f1 + n2;   // mid point of what is left
  for(i = 0 ; i < n0 ; i++) {  // first slice (0 -> 3 elements used)
    sum1 += f[i];
  }
  for(i = 0 ; i < n2 ; i++) {  // two add streams, one from beginning, one from midpoint
    sum1 += f1[i];
    sum2 += f2[i];
  }
  sum1 += sum2;
#endif
  return sum1;
}

