#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <immintrin.h>

// NOTE: Fortran order assumed for matrices
/*
  basic block, 8x8 transpose of 32 bit elements within 256 bit AVX/AVX2 registers
  3 pass shuffle, first 2 operations are within 128 bit lane, last one is across 128 bit lanes

     low part   high part                   low part   high part                   low part   high part
   <---------> <--------->  32 bit op     <---------> <--------->  64 bit op     <---------> <--------->  128 bit op    <---------> <--------->
0  a0 a1 a2 a3 a4 a5 a6 a7  unpacklo(0,1) a0 b0 a1 b1 a4 b4 a5 b5  unpacklo(0,2) a0 b0 c0 d0 a1 b1 c1 d1  perm(0,4)LL   a0 b0 c0 d0 e0 f0 g0 h0
1  b0 b1 b2 b3 b4 b5 b6 b7  unpackhi(0,1) a2 b2 a3 b3 a6 b6 a7 b7  unpackhi(0,2) a4 b4 c4 d4 a5 b5 c5 d5  perm(0,4)HH   a1 b1 c1 d1 e1 f1 g1 h1
2  c0 c1 c2 c3 c4 c5 c6 c7  unpacklo(2,3) c0 d0 c1 d1 c4 d4 c5 d5  unpacklo(1,3) a2 b2 c2 d2 a3 b3 c3 d3  perm(2,6)LL   a2 b2 c2 d2 e2 f2 g2 h2
3  d0 d1 d2 d3 d4 d5 d6 d7  unpackhi(2,3) c2 d2 c3 d3 c6 d6 c7 d7  unpackhi(1,3) a6 b6 c6 d6 a7 b7 c7 d7  perm(2,6)HH   a3 b3 c3 d3 e3 f3 g3 h3
4  e0 e1 e2 e3 e4 e5 e6 e7  unpacklo(4,5) e0 f0 e1 f1 e4 f4 e5 f5  unpacklo(4,6) e0 f0 g0 h0 e1 f1 g1 h1  perm(1,5)LL   a4 b4 c4 d4 e4 f4 g4 h4
5  f0 f1 f2 f3 f4 f5 f6 f7  unpackhi(4,5) e2 f2 e3 f3 e6 f6 e7 f7  unpackhi(4,6) e4 f4 g4 h4 e5 f5 g5 h5  perm(1,5)HH   a5 b5 c5 d5 e5 f5 g5 h5
6  g0 g1 g2 g3 g4 g5 g6 g7  unpacklo(6,7) g0 h0 g1 h1 g4 h4 g5 h5  unpacklo(5,7) e2 f2 g2 h2 e3 f3 g3 h3  perm(3,7)LL   a6 b6 c6 d6 e6 f6 g6 h6
7  h0 h1 h2 h3 h4 h5 h6 h7  unpackhi(6,7) g2 h2 g3 h3 g6 h6 g7 h7  unpackhi(5,7) e6 f6 g6 h6 e7 f7 g7 h7  perm(3,7)HH   a7 b7 c7 d7 e7 f7 g7 h7

permLL uses _mm256_permute2f128_si256 (a,b, 0x20), bot128 = bot128 from a, top128 = bot128 from b
permHH uses _mm256_permute2f128_si256 (a,b, 0x31), bot128 = top128 from a, top128 = top128 from b
*/
#if defined(__AVX2__) && defined(__x86_64__)
void Transpose32_4x4(uint32_t *a00, int la1, int la2, int la3, uint32_t *b00, int lb1, int lb2, int lb3){
  __m128i r0, r1, r2, r3;
  __m128i d0, d1, d2, d3;
  __m128i s0, s1, s2, s3;

  r0 = _mm_loadu_si128((__m128i const *) &a00[0  ]);  // a0 a1 a2 a3
  r1 = _mm_loadu_si128((__m128i const *) &a00[la1]);  // b0 b1 b2 b3
  r2 = _mm_loadu_si128((__m128i const *) &a00[la2]);  // c0 c1 c2 c3
  r3 = _mm_loadu_si128((__m128i const *) &a00[la3]);  // d0 d1 d2 d3

  d0 = _mm_unpacklo_epi32(r0,r1);  // a0 b0 a1 b1
  d1 = _mm_unpackhi_epi32(r0,r1);  // a2 b2 a3 b3
  d2 = _mm_unpacklo_epi32(r2,r3);  // c0 d0 c1 d1
  d3 = _mm_unpackhi_epi32(r2,r3);  // c2 d2 c3 d3

  s0 = _mm_unpacklo_epi64(d0,d2);  // a0 b0 c0 d0
  s1 = _mm_unpackhi_epi64(d0,d2);  // a1 b1 c1 d1
  s2 = _mm_unpacklo_epi64(d1,d3);  // a2 b2 c2 d2
  s3 = _mm_unpackhi_epi64(d1,d3);  // a3 b3 c3 d3

  _mm_storeu_si128((__m128i *) &b00[  0],s0);  // a0 b0 c0 d0
  _mm_storeu_si128((__m128i *) &b00[lb1],s1);  // a1 b1 c1 d1
  _mm_storeu_si128((__m128i *) &b00[lb2],s2);  // a2 b2 c2 d2
  _mm_storeu_si128((__m128i *) &b00[lb3],s3);  // a3 b3 c3 d3
}

void Transpose32_4x8(uint32_t *a00, int la1, int la2, int la3, uint32_t *b00, int lb1, int lb2, int lb3, int lb4){
  uint32_t *b01 = b00+lb4;
  __m256i y0, y1, y2, y3;
  __m256i t0, t1, t2, t3;

  y0 = _mm256_loadu_si256 ((__m256i const *) &a00[  0]);  // a0 a1 a2 a3 a4 a5 a6 a7
  y1 = _mm256_loadu_si256 ((__m256i const *) &a00[la1]);  // b0 b1 b2 b3 b4 b5 b6 b7
  y2 = _mm256_loadu_si256 ((__m256i const *) &a00[la2]);  // c0 c1 c2 c3 c4 c5 c6 c7
  y3 = _mm256_loadu_si256 ((__m256i const *) &a00[la3]);  // d0 d1 d2 d3 d4 d5 d6 d7

  t0 = _mm256_unpacklo_epi32(y0,y1); // a0 b0 a2 b2 a4 b4 a6 b6
  t1 = _mm256_unpackhi_epi32(y0,y1); // a1 b1 a3 b3 a5 b5 a7 b7
  t2 = _mm256_unpacklo_epi32(y2,y3); // c0 d0 c2 d2 c4 d4 c6 d6
  t3 = _mm256_unpackhi_epi32(y2,y3); // c1 d1 c3 d3 c5 d5 c7 d7

  y0 = _mm256_unpacklo_epi64(t0,t2); // a0 b0 c0 d0 a4 b4 c4 d4
  y1 = _mm256_unpackhi_epi64(t0,t2); // a2 b2 c2 d2 a6 b6 c6 d6
  y2 = _mm256_unpacklo_epi64(t1,t3); // a1 b1 c1 d1 a5 b5 c5 d5
  y3 = _mm256_unpackhi_epi64(t1,t3); // a3 b3 c3 d3 a7 b7 c7 d7

  _mm_storeu_si128((__m128i *) &b00[  0], _mm256_extracti128_si256 (y0, 0) ); // a0 b0 c0 d0
  _mm_storeu_si128((__m128i *) &b01[  0], _mm256_extracti128_si256 (y0, 1) ); // a4 b4 c4 d4
  _mm_storeu_si128((__m128i *) &b00[lb1], _mm256_extracti128_si256 (y1, 0) ); // a2 b2 c2 d2
  _mm_storeu_si128((__m128i *) &b01[lb1], _mm256_extracti128_si256 (y1, 1) ); // a6 b6 c6 d6
  _mm_storeu_si128((__m128i *) &b00[lb2], _mm256_extracti128_si256 (y2, 0) ); // a1 b1 c1 d1
  _mm_storeu_si128((__m128i *) &b01[lb2], _mm256_extracti128_si256 (y2, 1) ); // a5 b5 c5 d5
  _mm_storeu_si128((__m128i *) &b00[lb3], _mm256_extracti128_si256 (y3, 0) ); // a3 b3 c3 d3
  _mm_storeu_si128((__m128i *) &b01[lb3], _mm256_extracti128_si256 (y3, 1) ); // a7 b7 c7 d7
  
}

void Transpose32_8x4(uint32_t *a00, int la1, int la2, int la3, int la4, uint32_t *b00, int lb1, int lb2, int lb3){
  uint32_t *a01 = a00+la4;
  __m256i y0, y1, y2, y3;
  __m256i t0, t1, t2, t3, tt;

  tt = _mm256_xor_si256(tt,tt);  // dummy statement for the following insertf128 instructions

  t0 = _mm256_insertf128_si256(tt, _mm_loadu_si128((__m128i const *) &a00[0  ]) , 0); // a0 a1 a2 a3
  t0 = _mm256_insertf128_si256(t0, _mm_loadu_si128((__m128i const *) &a01[0  ]) , 1); // a0 a1 a2 a3 e0 e1 e2 e3
  t1 = _mm256_insertf128_si256(tt, _mm_loadu_si128((__m128i const *) &a00[la1]) , 0); // b0 b1 b2 b3
  t1 = _mm256_insertf128_si256(t1, _mm_loadu_si128((__m128i const *) &a01[la1]) , 1); // b0 b1 b2 b3 f0 f1 f2 f3
  t2 = _mm256_insertf128_si256(tt, _mm_loadu_si128((__m128i const *) &a00[la2]) , 0); // c0 c1 c2 c3
  t2 = _mm256_insertf128_si256(t2, _mm_loadu_si128((__m128i const *) &a01[la2]) , 1); // c0 c1 c2 c3 g0 g1 g2 g3
  t3 = _mm256_insertf128_si256(tt, _mm_loadu_si128((__m128i const *) &a00[la3]) , 0); // d0 d1 d2 d3
  t3 = _mm256_insertf128_si256(t3, _mm_loadu_si128((__m128i const *) &a01[la3]) , 1); // d0 d1 d2 d3 h0 h1 h2 h3

  y0 = _mm256_unpacklo_epi32(t0,t1); // a0 b0 a1 b1 e0 f0 e1 f1
  y1 = _mm256_unpackhi_epi32(t0,t1); // a2 b2 a3 b3 e2 f2 e3 f3
  y2 = _mm256_unpacklo_epi32(t2,t3); // c0 d0 c1 d1 g0 h0 g1 h1
  y3 = _mm256_unpackhi_epi32(t2,t3); // c2 d2 c3 d3 g2 h2 g3 h3

  t0 = _mm256_unpacklo_epi64(y0,y2); // a0 b0 c0 d0 e0 f0 g0 h0
  t1 = _mm256_unpackhi_epi64(y0,y2); // a1 b1 c1 d1 e1 f1 g1 h1
  t2 = _mm256_unpacklo_epi64(y1,y3); // a2 b2 c2 d2 e2 f2 g2 h2
  t3 = _mm256_unpackhi_epi64(y1,y3); // a3 b3 c3 d3 e3 f3 g3 h3

  _mm256_storeu_si256((__m256i *) &b00[  0],t0); // a0 b0 c0 d0 e0 f0 g0 h0
  _mm256_storeu_si256((__m256i *) &b00[lb1],t1); // a1 b1 c1 d1 e1 f1 g1 h1
  _mm256_storeu_si256((__m256i *) &b00[lb2],t2); // a2 b2 c2 d2 e2 f2 g2 h2
  _mm256_storeu_si256((__m256i *) &b00[lb3],t3); // a3 b3 c3 d3 e3 f3 g3 h3
}

void Transpose32_8x8(uint32_t *a00, int la1, int la2, int la3, int la4, uint32_t *b00, int lb1, int lb2, int lb3, int lb4){
  uint32_t *a01 = a00+la4;
  uint32_t *b01 = b00+lb4;
  __m256i t0, t1, t2, t3, t4, t5, t6, t7;
  __m256i y0, y1, y2, y3, y4, y5, y6, y7;
  __m256i x0, x1, x2, x3, x4, x5, x6, x7;

  y0 = _mm256_loadu_si256 ((__m256i const *) &a00[  0]); // a0 a1 a2 a3 a4 a5 a6 a7
  y1 = _mm256_loadu_si256 ((__m256i const *) &a00[la1]); // b0 b1 b2 b3 b4 b5 b6 b7
  y2 = _mm256_loadu_si256 ((__m256i const *) &a00[la2]); // c0 c1 c2 c3 c4 c5 c6 c7
  y3 = _mm256_loadu_si256 ((__m256i const *) &a00[la3]); // d0 d1 d2 d3 d4 d5 d6 d7
  y4 = _mm256_loadu_si256 ((__m256i const *) &a01[  0]); // e0 e1 e2 e3 e4 e5 e6 e7
  y5 = _mm256_loadu_si256 ((__m256i const *) &a01[la1]); // f0 f1 f2 f3 f4 f5 f6 f7
  y6 = _mm256_loadu_si256 ((__m256i const *) &a01[la2]); // g0 g1 g2 g3 g4 g5 g6 g7
  y7 = _mm256_loadu_si256 ((__m256i const *) &a01[la3]); // h0 h1 h2 h3 h4 h5 h6 h7

  t0 = _mm256_unpacklo_epi32(y0,y1); // a0 b0 a1 b1 a4 b4 a5 b5
  t1 = _mm256_unpackhi_epi32(y0,y1); // a2 b2 a3 b3 a6 b6 a7 b7
  t2 = _mm256_unpacklo_epi32(y2,y3); // c0 d0 c1 d1 c4 d4 c5 d5
  t3 = _mm256_unpackhi_epi32(y2,y3); // c2 d2 c3 d3 c6 d6 c7 d7
  t4 = _mm256_unpacklo_epi32(y4,y5); // e0 f0 e1 f1 e4 f4 e5 f5
  t5 = _mm256_unpackhi_epi32(y4,y5); // e2 f2 e3 f3 e6 f6 e7 f7
  t6 = _mm256_unpacklo_epi32(y6,y7); // g0 h0 g1 h1 g4 h4 g5 h5
  t7 = _mm256_unpackhi_epi32(y6,y7); // g2 h2 g3 h3 g6 h6 g7 h7

  y0 = _mm256_unpacklo_epi64(t0,t2); // a0 b0 c0 d0 a1 b1 c1 d1
  y1 = _mm256_unpackhi_epi64(t0,t2); // a4 b4 c4 d4 a5 b5 c5 d5
  y2 = _mm256_unpacklo_epi64(t1,t3); // a2 b2 c2 d2 a3 b3 c3 d3
  y3 = _mm256_unpackhi_epi64(t1,t3); // a6 b6 c6 d6 a7 b7 c7 d7
  y4 = _mm256_unpacklo_epi64(t4,t6); // e0 f0 g0 h0 e1 f1 g1 h1
  y5 = _mm256_unpackhi_epi64(t4,t6); // e4 f4 g4 h4 e5 f5 g5 h5
  y6 = _mm256_unpacklo_epi64(t5,t7); // e2 f2 g2 h2 e3 f3 g3 h3
  y7 = _mm256_unpackhi_epi64(t5,t7); // e6 f6 g6 h6 e7 f7 g7 h7

  t0 = _mm256_permute2f128_si256 (y0,y4, 0x20); // a0 b0 c0 d0 e0 f0 g0 h0
  t4 = _mm256_permute2f128_si256 (y0,y4, 0x31); // a1 b1 c1 d1 e1 f1 g1 h1
  t1 = _mm256_permute2f128_si256 (y1,y5, 0x20); // a2 b2 c2 d2 e2 f2 g2 h2
  t5 = _mm256_permute2f128_si256 (y1,y5, 0x31); // a3 b3 c3 d3 e3 f3 g3 h3
  t2 = _mm256_permute2f128_si256 (y2,y6, 0x20); // a4 b4 c4 d4 e4 f4 g4 h4
  t6 = _mm256_permute2f128_si256 (y2,y6, 0x31); // a5 b5 c5 d5 e5 f5 g5 h5
  t3 = _mm256_permute2f128_si256 (y3,y7, 0x20); // a6 b6 c6 d6 e6 f6 g6 h6
  t7 = _mm256_permute2f128_si256 (y3,y7, 0x31); // a7 b7 c7 d7 e7 f7 g7 h7

  _mm256_storeu_si256((__m256i *) &b00[  0],t0); // a0 b0 c0 d0 e0 f0 g0 h0
  _mm256_storeu_si256((__m256i *) &b01[  0],t4); // a1 b1 c1 d1 e1 f1 g1 h1
  _mm256_storeu_si256((__m256i *) &b00[lb1],t1); // a2 b2 c2 d2 e2 f2 g2 h2
  _mm256_storeu_si256((__m256i *) &b01[lb1],t5); // a3 b3 c3 d3 e3 f3 g3 h3
  _mm256_storeu_si256((__m256i *) &b00[lb2],t2); // a4 b4 c4 d4 e4 f4 g4 h4
  _mm256_storeu_si256((__m256i *) &b01[lb2],t6); // a5 b5 c5 d5 e5 f5 g5 h5
  _mm256_storeu_si256((__m256i *) &b00[lb3],t3); // a6 b6 c6 d6 e6 f6 g6 h6
  _mm256_storeu_si256((__m256i *) &b01[lb3],t7); // a7 b7 c7 d7 e7 f7 g7 h7
  b00 = b00 + lb4 + lb4;   // bump by 8 rows
}
#endif
#if ! defined(TEST_ONLY)
int TransposeBy4bytes(void *a, int la1, void *b, int lb1, int ni, int nj){
  int *a00, *a01, *b00, *b01;
  int *a000, *b000;
  int i, j, j0, ja0, ja00, jb0, jb00;
  int la2 = la1 + la1;
  int la3 = la2 + la1;
  int la4 = la3 + la1;
  int lb2 = lb1 + lb1;
  int lb3 = lb2 + lb1;
  int lb4 = lb3 + lb1;
  int ni8 ;
  int nj8 ;
  __m256i t0, t1, t2, t3, t4, t5, t6, t7;
  __m256i y0, y1, y2, y3, y4, y5, y6, y7;
  __m256i x0, x1, x2, x3, x4, x5, x6, x7;
  __m128i r0, r1, r2, r3;
  __m128i d0, d1, d2, d3;
  __m128i s0, s1, s2, s3;

  if (ni <= 0) ni = la1 ;
  if (nj <= 0) nj = lb1;
  ni8 = ni & 0x7FFFFFF8 ;   // lower multiple of 8
  nj8 = nj & 0x7FFFFFF8 ;   // lower multiple of 8
  if ( (la1 < ni) || (lb1 < nj) ) {
    fprintf(stderr,"ERROR: TransposeBy4bytes, la1=%d < ni=%d or lb1=%d < nj=%d\n",la1,ni,lb1,nj);
    return(1);
  }
// b[j,i] = a[i,j]
  a000 = (int *)a;
  b000 = (int *)b;
// basic block for transpose is 8x8
#if defined(FETCH128)
  ni8 = ni & 0x7FFFFFFC ;   // lower multiple of 4
  nj8 = nj & 0x7FFFFFFC ;   // lower multiple of 4
  for (j=0 ; j<nj8 ; j+= 4){
    a00 = a000;        // base address for row
//     a01 = a00 + la4;   // base of next block of 4x8
    b00 = b000;   // base address for column
    for ( i=0 ; i<ni8 ; i+=4){     // transpose a[i:i+7,J:j+7] into b[j:j+7,i:i+7] using 4 x 4 blocks

      r0 = _mm_loadu_si128((__m128i const *) &a00[0  ]);    // a0 a1 a2 a3
      r1 = _mm_loadu_si128((__m128i const *) &a00[la1]);    // b0 b1 b2 b3
      r2 = _mm_loadu_si128((__m128i const *) &a00[la2]);    // c0 c1 c2 c3
      r3 = _mm_loadu_si128((__m128i const *) &a00[la3]);    //d0 d1 d2 d3

      s0 = _mm_unpacklo_epi32(r0,r1);    // a0 b0 a1 b1
      s1 = _mm_unpackhi_epi32(r0,r1);    // a2 b2 a3 b3
      s2 = _mm_unpacklo_epi32(r2,r3);    // c0 d0 c1 d1
      s3 = _mm_unpackhi_epi32(r2,r3);    // c2 d2 c3 d3

      d0 = _mm_unpacklo_epi64(s0,s2);    // a0 b0 c0 d0
      d1 = _mm_unpackhi_epi64(s0,s2);    // a1 b1 c1 d1
      d2 = _mm_unpacklo_epi64(s1,s3);    // a2 b2 c2 d2
      d3 = _mm_unpackhi_epi64(s1,s3);    // a3 b3 c3 d3

      _mm_storeu_si128((__m128i *)&b00[  0] , d0);
      _mm_storeu_si128((__m128i *)&b00[lb1] , d1);
      _mm_storeu_si128((__m128i *)&b00[lb2] , d2);
      _mm_storeu_si128((__m128i *)&b00[lb3] , d3);

      a00 = a00 + 4;           // bump by 8 cols
//       a00 = a00 + 8;           // bump by 8 cols
//       a01 = a00 + la4;         // base of next block of 4x8

//       b01 = b00 + lb4;         // base of next block of 4x8
//       b00 = b00 + lb4 + lb4;   // bump by 8 rows
      b00 = b00 + lb4;          // bump by 4 rows
    } // for i
    for ( ; i<ni ; i++){  // leftovers along i , j -> j+7
      b00[0] = a00[  0] ;  // b[j  ,i] = a[i,j  ]
      b00[1] = a00[la1] ;  // b[j+1,i] = a[i,j+1]
      b00[2] = a00[la2] ;  // b[j+2,i] = a[i,j+2]
      b00[3] = a00[la3] ;  // b[j+3,i] = a[i,j+3]
//       b00[4] = a00[la4    ] ;  // b[j+4,i] = a[i,j+4]
//       b00[5] = a00[la4+la1] ;  // b[j+5,i] = a[i,j+5]
//       b00[6] = a00[la4+la2] ;  // b[j+6,i] = a[i,j+6]
//       b00[7] = a00[la4+la3] ;  // b[j+7,i] = a[i,j+7]
      b00 += lb1;
      a00++;
    }
//     a000 = a000 + la4 + la4;   // bump by 8 rows
//     b000 = b000 + 8;           // bump by 8 cols
    a000 = a000 + la4 ;        // bump by 4 rows
    b000 = b000 + 4;           // bump by 4 cols
  }  // for j
#else
  for (j=0 ; j<nj8 ; j+= 8){
    a00 = a000;   // base address for row
    a01 = a00 + la4;   // base of next block of 4x8
    b00 = b000;   // base address for column
    for ( i=0 ; i<ni8 ; i+=8){     // transpose a[i:i+7,J:j+7] into b[j:j+7,i:i+7]
      y0 = _mm256_loadu_si256 ((__m256i const *) &a00[  0]);
      y1 = _mm256_loadu_si256 ((__m256i const *) &a00[la1]);
      y2 = _mm256_loadu_si256 ((__m256i const *) &a00[la2]);
      y3 = _mm256_loadu_si256 ((__m256i const *) &a00[la3]);
      y4 = _mm256_loadu_si256 ((__m256i const *) &a01[  0]);
      y5 = _mm256_loadu_si256 ((__m256i const *) &a01[la1]);
      y6 = _mm256_loadu_si256 ((__m256i const *) &a01[la2]);
      y7 = _mm256_loadu_si256 ((__m256i const *) &a01[la3]);
      a00 = a00 + 8;           // bump by 8 cols
      a01 = a00 + la4;   // base of next block of 4x8

      t0 = _mm256_unpacklo_epi32(y0,y1); // _mm_prefetch (a00    , 3) ;
      t1 = _mm256_unpackhi_epi32(y0,y1); // _mm_prefetch (a00+la1, 3) ;
      t2 = _mm256_unpacklo_epi32(y2,y3); // _mm_prefetch (a00+la2, 3) ;
      t3 = _mm256_unpackhi_epi32(y2,y3); // _mm_prefetch (a00+la3, 3) ;
      t4 = _mm256_unpacklo_epi32(y4,y5); // _mm_prefetch (a01    , 3) ;
      t5 = _mm256_unpackhi_epi32(y4,y5); // _mm_prefetch (a01+la1, 3) ;
      t6 = _mm256_unpacklo_epi32(y6,y7); // _mm_prefetch (a01+la2, 3) ;
      t7 = _mm256_unpackhi_epi32(y6,y7); // _mm_prefetch (a01+la3, 3) ;

      y0 = _mm256_unpacklo_epi64(t0,t2);
      y1 = _mm256_unpackhi_epi64(t0,t2);
      y2 = _mm256_unpacklo_epi64(t1,t3);
      y3 = _mm256_unpackhi_epi64(t1,t3);
      y4 = _mm256_unpacklo_epi64(t4,t6);
      y5 = _mm256_unpackhi_epi64(t4,t6);
      y6 = _mm256_unpacklo_epi64(t5,t7);
      y7 = _mm256_unpackhi_epi64(t5,t7);

      t0 = _mm256_permute2f128_si256 (y0,y4, 0x20);
      t4 = _mm256_permute2f128_si256 (y0,y4, 0x31);
      t1 = _mm256_permute2f128_si256 (y1,y5, 0x20);
      t5 = _mm256_permute2f128_si256 (y1,y5, 0x31);
      t2 = _mm256_permute2f128_si256 (y2,y6, 0x20);
      t6 = _mm256_permute2f128_si256 (y2,y6, 0x31);
      t3 = _mm256_permute2f128_si256 (y3,y7, 0x20);
      t7 = _mm256_permute2f128_si256 (y3,y7, 0x31);

      b01 = b00 + lb4;   // base of next block of 4x8
      _mm256_storeu_si256((__m256i *) &b00[  0],t0);
      _mm256_storeu_si256((__m256i *) &b01[  0],t4);
      _mm256_storeu_si256((__m256i *) &b00[lb1],t1);
      _mm256_storeu_si256((__m256i *) &b01[lb1],t5);
      _mm256_storeu_si256((__m256i *) &b00[lb2],t2);
      _mm256_storeu_si256((__m256i *) &b01[lb2],t6);
      _mm256_storeu_si256((__m256i *) &b00[lb3],t3);
      _mm256_storeu_si256((__m256i *) &b01[lb3],t7);
      b00 = b00 + lb4 + lb4;   // bump by 8 rows
    } // for i
    for ( ; i<ni ; i++){  // leftovers along i , j -> j+7
      b00[0] = a00[  0] ;  // b[j  ,i] = a[i,j  ]
      b00[1] = a00[la1] ;  // b[j+1,i] = a[i,j+1]
      b00[2] = a00[la2] ;  // b[j+2,i] = a[i,j+2]
      b00[3] = a00[la3] ;  // b[j+3,i] = a[i,j+3]
      b00[4] = a00[la4    ] ;  // b[j+4,i] = a[i,j+4]
      b00[5] = a00[la4+la1] ;  // b[j+5,i] = a[i,j+5]
      b00[6] = a00[la4+la2] ;  // b[j+6,i] = a[i,j+6]
      b00[7] = a00[la4+la3] ;  // b[j+7,i] = a[i,j+7]
      b00 += lb1;
      a00++;
    }
    a000 = a000 + la4 + la4;   // bump by 8 rows
    b000 = b000 + 8;           // bump by 8 cols
  }  // for j
#endif
  a000 = (int *)a;
  b000 = (int *)b;
  jb00 = j;
  ja00 = j*la1;
  for (i = 0 ; i<ni ; i++){    // leftovers along j, i = 0 -> nj-1
    ja0 = ja00;
    jb0 = jb00;
    for(j0=j ; j0<nj ; j0++){
      b000[jb0] = a000[ja0];   // b[j,i] = a[i,j]
      ja0 += la1;
      jb0++;
    }
    jb00 += lb1;
    ja00++;
  }
  return(0);
}

int TransposeBy8bytes(void *a, int la1, void *b, int lb1, int ni, int nj){
  // after transpose, b[j,i] = a[i,j]
  int la2 = la1 + la1;
  int la3 = la2 + la1;
  int la4 = la3 + la1;
  int la8 = la4 + la4;
  int lb2 = lb1 + lb1;
  int lb3 = lb2 + lb1;
  int lb4 = lb3 + lb1;
  int i, j, jj, ni8, nj8;
  long long *a000, *a00, *a01, *a10, *a11, *b000, *b00, *b01, *b10, *b11;
  __m256i y0, y1, y2, y3;
  __m256i t0, t1, t2, t3;
  __m256i x0, x1, x2, x3;
  __m128i p0, p1;

  if (ni <= 0) ni = la1 ;
  if (nj <= 0) nj = lb1;
  if ( (la1 < ni) || (lb1 < nj) ) {
    fprintf(stderr,"ERROR: TransposeBy8bytes, la1=%d < ni=%d or lb1=%d < nj=%d\n",la1,ni,lb1,nj);
    return(1);
  }
//   printf("la1,ni,lb1,nj= %d %d %d %d\n",la1,ni,lb1,nj);
  a000 = (long long *)a;
  b000 = (long long *)b;
  ni8 = ni & 0x7FFFFFF8 ;   // lower multiple of 8
  nj8 = nj & 0x7FFFFFF8 ;   // lower multiple of 8

#if defined(FETCH128)
// following 2 macros missing in gcc primitives
// make 256 bit item (ymm) from 2 x 128 bit items (xmm)
#define _mm256i_from_2m128i( lo, hi) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)
// fetch 2 128 bit items and load them into a 256 bit register (ymm)
#define _mm256_loadu_lh__m128i( loaddr, hiaddr) \
    _mm256i_from_2m128i(_mm_loadu_si128(loaddr), _mm_loadu_si128(hiaddr))

// faster if ni * nj <= 2048 ?
  for ( i=0 ; i<ni8 ; i+=4){
    a00 = &a000[i] ;
    b00 = &b000[i*lb1] ;
    for (j=0 ; j<nj8 ; j+= 4){
      y0 = _mm256_loadu_lh__m128i( (__m128i const*) &a00[    0] , (__m128i const*) &a00[la2  ] ) ;   // a0 a1 c0 c1
      y1 = _mm256_loadu_lh__m128i( (__m128i const*) &a00[la1  ] , (__m128i const*) &a00[la3  ] ) ;   // b0 b1 d0 d1
      y2 = _mm256_loadu_lh__m128i( (__m128i const*) &a00[    2] , (__m128i const*) &a00[la2+2] ) ;   // a2 a3 c2 c3
      y3 = _mm256_loadu_lh__m128i( (__m128i const*) &a00[la1+2] , (__m128i const*) &a00[la3+2] ) ;   // b2 b3 d2 d3

      x0 = _mm256_unpacklo_epi64(y0,y1);              // a0 b0 c0 d0
      x1 = _mm256_unpackhi_epi64(y0,y1);              // a1 b1 c1 d1
      x2 = _mm256_unpacklo_epi64(y2,y3);              // a2 b2 c2 d2
      x3 = _mm256_unpackhi_epi64(y2,y3);              // a3 b3 c3 d3

      _mm256_storeu_si256((__m256i *)&b00[  0],x0);   // store 4 x 4 into b00 (from a00)
      _mm256_storeu_si256((__m256i *)&b00[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b00[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b00[lb3],x3);

      a00 += la4 ;
      b00 += 4;
    }
  }
#else
// faster if ni * nj > 2048
  for ( i=0 ; i<ni8 ; i+=8){
    a00 = &a000[i] ;
    b00 = &b000[i*lb1] ;
    for (j=0 ; j<nj8 ; j+= 8){
      a01 = a00 + 4 ; a10 = a00 + la4 ; a11 = a10 + 4; 
      b01 = b00 + 4 ; b10 = b00 + lb4 ; b11 = b10 + 4;
      // basic 8 x 8 transpose (a00, a01, a10, a11 into b00, b10, b01, b11) ( 4 4x4 sub matrices)
      // fetch 4 x 4 from a00
      y0 = _mm256_loadu_si256((__m256i const *)&a00[  0]);   // a0 a1 a2 a3
      y1 = _mm256_loadu_si256((__m256i const *)&a00[la1]);   // b0 b1 b2 b3
      y2 = _mm256_loadu_si256((__m256i const *)&a00[la2]);   // c0 c1 c2 c3
      y3 = _mm256_loadu_si256((__m256i const *)&a00[la3]);   // d0 d1 d2 d3

      // shuffle pass 1 on a00 -> b00
      t0 = _mm256_unpackhi_epi64(y0,y1);                     // a1 b1 a3 b3
      t1 = _mm256_unpackhi_epi64(y2,y3);                     // c1 d1 c3 d3
      t2 = _mm256_unpacklo_epi64(y0,y1);                     // a0 b0 a2 b2
      t3 = _mm256_unpacklo_epi64(y2,y3);                     // c0 d0 c2 d2

      y0 = _mm256_loadu_si256((__m256i const *)&a01[  0]);   // prefetch 4 x 4 from a01
      y1 = _mm256_loadu_si256((__m256i const *)&a01[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a01[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a01[la3]);

      // shuffle pass 2 on a00 -> b00
      x0 = _mm256_permute2f128_si256(t2,t3,0x20);            // a0 b0 c0 d0
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);            // a1 b1 c1 d1
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);            // a2 b2 c2 d2
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);            // a3 b3 c3 d3

      _mm256_storeu_si256((__m256i *)&b00[  0],x0);   // store 4 x 4 into b00 (from a00)
      _mm256_storeu_si256((__m256i *)&b00[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b00[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b00[lb3],x3);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a01 -> b10
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      y0 = _mm256_loadu_si256((__m256i const *)&a10[  0]);   // prefetch 4 x 4 from a10
      y1 = _mm256_loadu_si256((__m256i const *)&a10[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a10[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a10[la3]);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a01 -> b10
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b10[  0],x0);   // store 4 x 4 into b10 (from a01)
      _mm256_storeu_si256((__m256i *)&b10[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b10[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b10[lb3],x3);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a10 -> b01
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      y0 = _mm256_loadu_si256((__m256i const *)&a11[  0]);   // prefetch 4 x 4 from a11
      y1 = _mm256_loadu_si256((__m256i const *)&a11[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a11[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a11[la3]);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a10 -> b01
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b01[  0],x0);   // store 4 x 4 into b01 (from a10)
      _mm256_storeu_si256((__m256i *)&b01[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b01[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b01[lb3],x3);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a11 -> b11
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a11 -> b11
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b11[  0],x0);   // store 4 x 4 into b11 (from a11)
      _mm256_storeu_si256((__m256i *)&b11[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b11[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b11[lb3],x3);
      a00 += la8 ;
      b00 += 8;
    }
  }
#endif
  for (i=ni8 ; i<ni ; i++){
    for ( jj=0 ; jj<nj ; jj++){   // b[jj,i] = a[i,jj]
      b000[jj + lb1*i] = a000[i + jj*la1];
    }
  }
  for (j=nj8 ; j<nj ; j++){
    for (i=0 ; i<ni ; i++){   // b[j,i] = a[i,j]
      b000[j + lb1*i] = a000[i + j*la1];
    }
  }
  return(0);
}
#endif

#if defined(UNIT_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define NI 15
#define NJ 13

#define PRINT for(j=0 ; j<9 ; j++){ for(i=0 ; i<9 ; i++){fprintf(stderr,"%2d ", s[j][i]);} ; for(i=0 ; i<9 ; i++){fprintf(stderr," %2d", d[j][i]);} fprintf(stderr,"\n");} fprintf(stderr,"\n\n")

void reset(void *p){
  uint32_t *t = memset(p, 0, sizeof(uint32_t)*NI*NJ);
}

main(){
  uint32_t s[NJ][NI], d[NI][NJ];  // s will get transposed into D
  uint32_t *t;
  uint32_t i, j, k;

  reset(s); reset(d);
  k=1;
  for(j=0 ; j<8 ; j++){
    for(i=0 ; i<8 ; i++){
      s[j][i] = k++;
    }
  }
//   PRINT; 
  fprintf(stderr,"================ Transpose32_8x8 ================\n");
  Transpose32_8x8(&s[0][0], NI, NI*2, NI*3, NI*4, &d[0][0], NJ, NJ*2, NJ*3, NJ*4);
  PRINT;
  reset(d);
  fprintf(stderr,"================ Transpose32_4x4 ================\n");
  Transpose32_4x4(&s[0][0], NI, NI*2, NI*3, &d[0][0], NJ, NJ*2, NJ*3);
  PRINT;
  reset(d);
  fprintf(stderr,"================ Transpose32_4x8 ================\n");
  Transpose32_4x8(&s[0][0], NI, NI*2, NI*3, &d[0][0], NJ, NJ*2, NJ*3, NJ*4);
  PRINT;
  reset(d);
  fprintf(stderr,"================ Transpose32_8x4 ================\n");
  Transpose32_8x4(&s[0][0], NI, NI*2, NI*3, NI*4, &d[0][0], NJ, NJ*2, NJ*3);
  PRINT;
}
#endif

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
static double cvt = 1.0 / 3.7E+9 ;  // assumes 3.7GHz clock

double MPI_Wtime(){
  return cvt * rdtsc();
}

main(int argc, char **argv){
  double *mtx1_4x4;
  double *mtx2_4x4;
  int *mat1_4x4;
  int *mat2_4x4;
  int i, j, k;
  double MPI_Wtime();
  double t0, t1;
  double bytes;
  int NI, NJ, ni, nj;
  int nni=atoi(argv[1]);
  int nnj=atoi(argv[2]);
  int nk=atoi(argv[5]);
  int errors;

  NI = atoi(argv[3]);
  NJ = atoi(argv[4]);
  ni = nni;
  nj = nnj;
  if(ni <= 0) ni = NI;
  if(nj <= 0) nj = NJ;
  mtx1_4x4 = (double *) malloc(sizeof(double)*NI*NJ);
  mtx2_4x4 = (double *) malloc(sizeof(double)*NI*NJ);
  mat1_4x4 = (int *) malloc(sizeof(int)*NI*NJ);
  mat2_4x4 = (int *) malloc(sizeof(int)*NI*NJ);
  printf("ni = %d, nj = %d, nk= %d\n",ni,nj,nk);
  for (i=0 ; i<NI*NJ ; i++){
    mtx1_4x4[i] = i;
    mtx2_4x4[i] = 999;
    mat1_4x4[i] = i;
    mat2_4x4[i] = 999;
  }
  printf("\nInput Matrix\n");
  if(NI*NJ < 10) {
    for(j=NJ-1 ; j>= 0 ; j--){
      printf("row %2d ",j);
      for (i=0 ; i<NI ; i++){
	printf("%4.0f",mtx1_4x4[i+NI*j]);
      }
      printf("\n");
    }
    printf("Col  < ");
    for (i=0 ; i<NI ; i++) printf("%4d",i) ;
    printf(">\n");
  }

  if(TransposeBy8bytes(mtx1_4x4, NI, mtx2_4x4, NJ, nni, nnj)) exit(1);

  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk/10 ; k++) {
    for(j=0 ; j<nj ; j++){
      for(i=0 ; i<ni ; i++){
	mat2_4x4[j + i*NJ] = mat1_4x4[i + NI*j] ;
      }
    }
    bytes = bytes + ni * nj * 16;
  }
  t1 = MPI_Wtime() ;
  printf("\nT8ref = %f, %f GB/s, %f ns/pt, %f us/transpose\n",
	 (t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*10*1000*1000*1000,(t1-t0)/nk*10*1000*1000);

  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk ; k++) {
    if(TransposeBy8bytes(mtx1_4x4, NI, mtx2_4x4, NJ, nni, nnj)) exit(1);
    bytes = bytes + ni * nj * 16;
  }
  t1 = MPI_Wtime() ;
  printf("\nT8 = %g, %f GB/s, %f ns/pt, %f us/transpose\n",
	 (t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*1000*1000*1000,(t1-t0)/nk*1000*1000);

  errors = 0;
  for(j=0 ; j<nj ; j++){
    for(i=0 ; i<ni ; i++){
      if( mtx1_4x4[i + NI*j] != mtx2_4x4[j + i*NJ] ) errors++;
    }
  }
  printf("\nOutput Matrix (double), errors = %d\n\n",errors);
  if(NI*NJ < 10) {
    for(j=NI-1 ; j>= 0 ; j--){
      printf("row %2d ",j);
      for (i=0 ; i<NJ ; i++){
	printf("%4.0f",mtx2_4x4[i+NJ*j]);
      }
      printf("\n");
    }
    printf("Col  < ");
    for (i=0 ; i<NJ ; i++) printf("%4d",i) ;
    printf(">\n");
  }

  if(TransposeBy4bytes(mat1_4x4, NI, mat2_4x4, NJ, nni, nnj)) exit(1);

  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk/10 ; k++) {
    for(j=0 ; j<nj ; j++){
      for(i=0 ; i<ni ; i++){
	mat2_4x4[j + i*NJ] = mat1_4x4[i + NI*j] ;
      }
    }
    bytes = bytes + ni * nj * 8;
  }
  t1 = MPI_Wtime() ;
  printf("\nT4ref = %f, %f GB/s, %f ns/pt, %f us/transpose\n",
	 (t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*10*1000*1000*1000,(t1-t0)/nk*10*1000*1000);

  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk ; k++) {
    if(TransposeBy4bytes(mat1_4x4, NI, mat2_4x4, NJ, nni, nnj)) exit(1);
    bytes = bytes + ni * nj * 8;
  }
  t1 = MPI_Wtime() ;

  errors = 0;
  for(j=0 ; j<nj ; j++){
    for(i=0 ; i<ni ; i++){
      if( mat1_4x4[i + NI*j] != mat2_4x4[j + i*NJ] ) errors++;
    }
  }
  printf("\nOutput Matrix (int), errors = %d\n",errors);
  printf("T4 = %g, %f GB/s, %f ns/pt, %f us/transpose\n",(t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*1000*1000*1000,(t1-t0)/nk*1000*1000);
  if(NI*NJ < 10) {
    for(j=NI-1 ; j>= 0 ; j--){
      printf("row %2d ",j);
      for (i=0 ; i<NJ ; i++){
	printf("%4d",mat2_4x4[i+NJ*j]);
      }
      printf("\n");
    }
    printf("Col  < ");
    for (i=0 ; i<NJ ; i++) printf("%4d",i) ;
    printf(">\n");
  }
}
#endif
