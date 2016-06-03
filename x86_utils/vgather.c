/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2015  Environnement Canada
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

static int I0123[] = {0,1,2,3};

#if defined(USE4x4x4)
void VGather4x4x4(int *src, int *dst, int ni, int ninj){
  int *base1;
  int ni2 = ni+ni;
#if defined(__AVX2__) && defined(__x86_64__)
  __m128i v128a, v128b;
  __m256i vindx;
  __m256i v1, v2, v3, v4, v5, v6, v7, v8;
#endif

#if defined(__AVX2__) && defined(__x86_64__)
  v128a = _mm_loadu_si128((__m128i const*) I0123);
  v128b = _mm_set1_epi32(ni);
  v128b = _mm_add_epi32(v128a,v128b);
  vindx = _mm256_inserti128_si256(vindx,v128a,0);
  vindx = _mm256_inserti128_si256(vindx,v128b,1);

  base1 = src;
  v1 = _mm256_i32gather_epi32((int const*) base1    , vindx, 4);
  v2 = _mm256_i32gather_epi32((int const*) base1+ni2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst +  0), v1);
  _mm256_storeu_si256((__m256i *) (dst +  8), v2);

  base1 += ninj;
  v3 = _mm256_i32gather_epi32((int const*) base1    , vindx, 4);
  v4 = _mm256_i32gather_epi32((int const*) base1+ni2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst + 16), v3);
  _mm256_storeu_si256((__m256i *) (dst + 24), v4);

  base1 += ninj;
  v5 = _mm256_i32gather_epi32((int const*) base1    , vindx, 4);
  v6 = _mm256_i32gather_epi32((int const*) base1+ni2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst + 32), v5);
  _mm256_storeu_si256((__m256i *) (dst + 40), v6);

  base1 += ninj;
  v7 = _mm256_i32gather_epi32((int const*) base1    , vindx, 4);
  v8 = _mm256_i32gather_epi32((int const*) base1+ni2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst + 48), v7);
  _mm256_storeu_si256((__m256i *) (dst + 56), v8);

#else
  dst[0] = src[0];
#endif
}
#endif

#if defined(USE4x4x4rd)
void VGather4x4x4rd(float *src, double *dst, int ni, int ninj, double *abcd){
  float *base1;
#if defined(__AVX2__) && defined(__x86_64__)
  __m128i vindx;
  __m128  t1, t2, t3, t4;
  __m256d v0, v1, v2, v3, v4, va, vb, vc, vd;
#endif

#if defined(__AVX2__) && defined(__x86_64__)
  vindx = _mm_loadu_si128((__m128i const*) I0123);

  base1 = src;
  t1 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  va = _mm256_set1_pd(abcd[0]);
  v1 = _mm256_cvtps_pd(t1);
  v0 = _mm256_mul_pd(v1,va);        // v0 = v1*a
  base1 = base1+ninj;
  t2 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  vb = _mm256_set1_pd(abcd[1]);
  v2 = _mm256_cvtps_pd(t2);
  v0 = _mm256_fmadd_pd(vb,v2,v0);   // v0 += v2*b
  base1 = base1+ninj;
  t3 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  vc = _mm256_set1_pd(abcd[2]);
  v3 = _mm256_cvtps_pd(t3);
  v0 = _mm256_fmadd_pd(vc,v3,v0);   // v0 += v3*c
  base1 = base1+ninj;
  t4 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  vd = _mm256_set1_pd(abcd[3]);
  v4 = _mm256_cvtps_pd(t4);
  v0 = _mm256_fmadd_pd(vd,v4,v0);   // v0 += v4*d
  _mm256_storeu_pd((double *) (dst +  0), v0);

  src += ni; base1 = src;
  t1 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v1 = _mm256_cvtps_pd(t1);
  v0 = _mm256_mul_pd(v1,va);        // v0 = v1*a
  base1 = base1+ninj;
  t2 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v2 = _mm256_cvtps_pd(t2);
  v0 = _mm256_fmadd_pd(vb,v2,v0);   // v0 += v2*b
  base1 = base1+ninj;
  t3 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v3 = _mm256_cvtps_pd(t3);
  v0 = _mm256_fmadd_pd(vc,v3,v0);   // v0 += v3*c
  base1 = base1+ninj;
  t4 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v4 = _mm256_cvtps_pd(t4);
  v0 = _mm256_fmadd_pd(vd,v4,v0);   // v0 += v4*d
  _mm256_storeu_pd((double *) (dst +  4), v0);

  src += ni; base1 = src;
  t1 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v1 = _mm256_cvtps_pd(t1);
  v0 = _mm256_mul_pd(v1,va);        // v0 = v1*a
  base1 = base1+ninj;
  t2 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v2 = _mm256_cvtps_pd(t2);
  v0 = _mm256_fmadd_pd(vb,v2,v0);   // v0 += v2*b
  base1 = base1+ninj;
  t3 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v3 = _mm256_cvtps_pd(t3);
  v0 = _mm256_fmadd_pd(vc,v3,v0);   // v0 += v3*c
  base1 = base1+ninj;
  t4 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v4 = _mm256_cvtps_pd(t4);
  v0 = _mm256_fmadd_pd(vd,v4,v0);   // v0 += v4*d
  _mm256_storeu_pd((double *) (dst +  8), v0);

  src += ni; base1 = src;
  t1 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v1 = _mm256_cvtps_pd(t1);
  v0 = _mm256_mul_pd(v1,va);        // v0 = v1*a
  base1 = base1+ninj;
  t2 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v2 = _mm256_cvtps_pd(t2);
  v0 = _mm256_fmadd_pd(vb,v2,v0);   // v0 += v2*b
  base1 = base1+ninj;
  t3 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v3 = _mm256_cvtps_pd(t3);
  v0 = _mm256_fmadd_pd(vc,v3,v0);   // v0 += v3*c
  base1 = base1+ninj;
  t4 = _mm_i32gather_ps((float const*) base1    , vindx, 4);
  v4 = _mm256_cvtps_pd(t4);
  v0 = _mm256_fmadd_pd(vd,v4,v0);   // v0 += v4*d
  _mm256_storeu_pd((double *) (dst + 12), v0);

#else
  dst[0] = src[0];
#endif
}
#endif

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#define DATA_SIZE 9000
main(){
  int my_data[DATA_SIZE];
  int i, j;
  int ni, ninj, limit;
  int collect[64];

  for (i=0 ; i<DATA_SIZE ; i++) my_data[i] = i;
  ni = 100 ;
  ninj = 1000;
  limit = DATA_SIZE-4*ninj+1;
  printf("limit=%d\n",limit);
  for (j=0 ; j<100000 ; j++) {
    for (i=0 ; i<limit ; i++) VGather4x4x4(my_data+i, collect, ni, ninj);
  }
//   for (i=0 ; i<64 ; i++) collect[i] = -1;
//   VGather4x4x4(my_data, collect, 100, 10000);
//   for (i=0 ; i<64 ; i++) printf("I = %d, G = %d \n",i,collect[i]);
}
#endif
