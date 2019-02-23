/* Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2019  Environnement Canada
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
#include <stdint.h>
#if defined(NO_SSE)
#undef __SSE__
#undef __AVX2__
#else
#include <immintrin.h>
#endif

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
  n2 = (n >> 3);        // number of full 8 token chunks
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
  for( ; i < n ; i++){    // loop over leftovers (everything if not AVX2
    t = *s1++;
    t = (t >> 32) | (t << 32);
    *d1++ = t;
  }
}
