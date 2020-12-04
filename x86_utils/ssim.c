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

// process 4 elements of x and y on 4 lines (16 values)
// ni = distance between 2 lines
// output 5 values
// sum(x), sum(x*x), sum(y), sum(y*y), sum(x*y)
// building block for SSIM score
void sums16(float *r, float *x, float *y, int ni){
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m128 vx0, vy0, sxy, vt0, vt1, vts;
  __m256 vxy, vs1, vs2, vst;
  int *l = (int *) r;
#else
  int i, j;
  float x2[4], xy[4], sx[4], sy[4], y2[4];
#endif

#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  vx0 = _mm_loadu_ps(x);   // row 0
  vy0 = _mm_loadu_ps(y);
  vxy = _mm256_insertf128_ps(vxy,vx0,0);   // _mm256_set_m128(vy0,vx0);
  vxy = _mm256_insertf128_ps(vxy,vy0,1);   // x[0,1,2,3] y[0,1,2,3]
  sxy = _mm_mul_ps(vy0,vx0);               // x*y[0,1,2,3]
  vs1 = vxy;                               // sums of x and y
  vs2 = _mm256_mul_ps(vxy,vxy);            // sums of x*x and y*y

  x += ni ; y += ni;
  vx0 = _mm_loadu_ps(x);   // row 1
  vy0 = _mm_loadu_ps(y);
  vxy = _mm256_insertf128_ps(vxy,vx0,0);   // _mm256_set_m128(vy0,vx0);
  vxy = _mm256_insertf128_ps(vxy,vy0,1);
  sxy = _mm_fmadd_ps(vx0,vy0,sxy);         // running sums of x*y
  vs1 = _mm256_add_ps(vxy,vs1);            // running sums of x and y
  vs2 = _mm256_fmadd_ps(vxy,vxy,vs2);      // running sums of x*x and y*y

  x += ni ; y += ni;
  vx0 = _mm_loadu_ps(x);   // row 2
  vy0 = _mm_loadu_ps(y);
  vxy = _mm256_insertf128_ps(vxy,vx0,0);   // _mm256_set_m128(vy0,vx0);
  vxy = _mm256_insertf128_ps(vxy,vy0,1);
  sxy = _mm_fmadd_ps(vx0,vy0,sxy);         // running sums of x*y
  vs1 = _mm256_add_ps(vxy,vs1);            // running sums of x and y
  vs2 = _mm256_fmadd_ps(vxy,vxy,vs2);      // running sums of x*x and y*y

  x += ni ; y += ni;
  vx0 = _mm_loadu_ps(x);   // row 3
  vy0 = _mm_loadu_ps(y);
  vxy = _mm256_insertf128_ps(vxy,vx0,0);   // _mm256_set_m128(vy0,vx0);
  vxy = _mm256_insertf128_ps(vxy,vy0,1);
  sxy = _mm_fmadd_ps(vx0,vy0,sxy);        // xy = sum(x*y)
  vs1 = _mm256_add_ps(vxy,vs1);           // x1 = sum(x)  , y1 = sum(y)
  vs2 = _mm256_fmadd_ps(vxy,vxy,vs2);     // x2 = sum(x*x), y2 = sum(y*y)

  // now fold sums
  vst = _mm256_hadd_ps(vs1,vs2);          // x1[0+1] x1[2+3]  x2[0+1] x2[2+3] y1[0+1] y1[2+3]  y2[0+1] y2[2+3] 
  vt0 = _mm256_extractf128_ps(vst,0);     // x1[0+1] x1[2+3]  x2[0+1] x2[2+3]
  vt1 = _mm256_extractf128_ps(vst,1);     // y1[0+1] y1[2+3]  y2[0+1] y2[2+3] 
  vts = _mm_hadd_ps(vt0,vt1);             // x1[0+1+2+3] x2[0+1+2+3] y1[0+1+2+3] y2[0+1+2+3]
  _mm_storeu_ps(r,vts);                   // sum(x), sum(x*x), sum(y), sum(y*y)
  vt0 = _mm_hadd_ps(sxy,sxy);             // xy[0+1] xy[2+3] xy[0+1] xy[2+3]
  vt1 = _mm_hadd_ps(vt0,vt0);             // xy[0+1+2+3] xy[0+1+2+3] xy[0+1+2+3] xy[0+1+2+3]
  l[4] = _mm_extract_epi32( (__m128i) vt1, 0);   // store first value, sum(x*y)
#else
  for(i=0 ; i<4 ; i++){        // row 0
    x2[i] = x[i] * x[i];
    y2[i] = y[i] * y[i];
    sx[i] = x[i];
    sy[i] = y[i];
    xy[i] = x[i] * y[i];
  }
  for(j=0 ; j<3 ; j++){        // next 3 rows
    x += ni ; y += ni;
    for(i=0 ; i<4 ; i++){
      x2[i] += x[i] * x[i];
      y2[i] += y[i] * y[i];
      sx[i] += x[i];
      sy[i] += y[i];
      xy[i] += x[i] * y[i];
    }
  }
  r[0] = sx[0] + sx[1] + sx[2] + sx[3];   // fold sums
  r[1] = x2[0] + x2[1] + x2[2] + x2[3];
  r[2] = sy[0] + sy[1] + sy[2] + sy[3];
  r[3] = y2[0] + y2[1] + y2[2] + y2[3];
  r[4] = xy[0] + xy[1] + xy[2] + xy[3];
#endif
}

// process 8 elements of x and y on 8 lines (64 values)
// ni = distance between 2 lines
// output 5 values
// sum(x), sum(x*x), sum(y), sum(y*y), sum(x*y)
// building block for SSIM score
void sums64(float *r, float *x, float *y, int ni){
  int j;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m128 vx0, vy0, sxy, vt0, vt1, vts;
  __m256 vxy, vs1, vs2, vst;
  int *l = (int *) r;
#else
  int i;
  float x2[8], xy[8], sx[8], sy[8], y2[8];
#endif

#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  vx0 = _mm_loadu_ps(x);                   // row 0
  vy0 = _mm_loadu_ps(y);
  vxy = _mm256_insertf128_ps(vxy,vx0,0);   // concatenate x and y
  vxy = _mm256_insertf128_ps(vxy,vy0,1);   // x[0,1,2,3] y[0,1,2,3]
  sxy = _mm_mul_ps(vy0,vx0);               // x*y[0,1,2,3]
  vs1 = vxy;                               // sums of x and y
  vs2 = _mm256_mul_ps(vxy,vxy);            // sums of x*x and y*y

  vx0 = _mm_loadu_ps(x+4);
  vy0 = _mm_loadu_ps(y+4);
  vxy = _mm256_insertf128_ps(vxy,vx0,0);   // concatenate x and y
  vxy = _mm256_insertf128_ps(vxy,vy0,1);   // x[0,1,2,3] y[0,1,2,3]
  sxy = _mm_fmadd_ps(vx0,vy0,sxy);         // running sums of x*y
  vs1 = _mm256_add_ps(vxy,vs1);            // running sums of x and y
  vs2 = _mm256_fmadd_ps(vxy,vxy,vs2);      // running sums of x*x and y*y

  for(j = 0 ; j < 7 ; j++){                  // next 7 rows
    x += ni ; y += ni;
    vx0 = _mm_loadu_ps(x);                   // next row
    vy0 = _mm_loadu_ps(y);
    vxy = _mm256_insertf128_ps(vxy,vx0,0);   // concatenate x and y
    vxy = _mm256_insertf128_ps(vxy,vy0,1);   // x[0,1,2,3] y[0,1,2,3]
    sxy = _mm_fmadd_ps(vx0,vy0,sxy);         // running sums of x*y
    vs1 = _mm256_add_ps(vxy,vs1);            // running sums of x and y
    vs2 = _mm256_fmadd_ps(vxy,vxy,vs2);      // running sums of x*x and y*y

    vx0 = _mm_loadu_ps(x+4);
    vy0 = _mm_loadu_ps(y+4);
    vxy = _mm256_insertf128_ps(vxy,vx0,0);   // concatenate x and y
    vxy = _mm256_insertf128_ps(vxy,vy0,1);   // x[0,1,2,3] y[0,1,2,3]
    sxy = _mm_fmadd_ps(vx0,vy0,sxy);         // running sums of x*y
    vs1 = _mm256_add_ps(vxy,vs1);            // running sums of x and y
    vs2 = _mm256_fmadd_ps(vxy,vxy,vs2);      // running sums of x*x and y*y
  }

  // now fold sums
  vst = _mm256_hadd_ps(vs1,vs2);          // x1[0+1] x1[2+3]  x2[0+1] x2[2+3] y1[0+1] y1[2+3]  y2[0+1] y2[2+3] 
  vt0 = _mm256_extractf128_ps(vst,0);     // x1[0+1] x1[2+3]  x2[0+1] x2[2+3]
  vt1 = _mm256_extractf128_ps(vst,1);     // y1[0+1] y1[2+3]  y2[0+1] y2[2+3] 
  vts = _mm_hadd_ps(vt0,vt1);             // x1[0+1+2+3] x2[0+1+2+3] y1[0+1+2+3] y2[0+1+2+3]
  _mm_storeu_ps(r,vts);                   // sum(x), sum(x*x), sum(y), sum(y*y)
  vt0 = _mm_hadd_ps(sxy,sxy);             // xy[0+1] xy[2+3] xy[0+1] xy[2+3]
  vt1 = _mm_hadd_ps(vt0,vt0);             // xy[0+1+2+3] xy[0+1+2+3] xy[0+1+2+3] xy[0+1+2+3]
  l[4] = _mm_extract_epi32( (__m128i) vt1, 0);   // store first value, sum(x*y)
#else
  for(i=0 ; i<8 ; i++){     // row 0
    x2[i] = x[i] * x[i];
    y2[i] = y[i] * y[i];
    sx[i] = x[i];
    sy[i] = y[i];
    xy[i] = x[i] * y[i];
  }
  for(j=0 ; j<7 ; j++){     // next 7 rows
    x += ni ; y += ni;
    for(i=0 ; i<8 ; i++){
      x2[i] += x[i] * x[i];
      y2[i] += y[i] * y[i];
      sx[i] += x[i];
      sy[i] += y[i];
      xy[i] += x[i] * y[i];
    }
  }
  r[0] = sx[0];             // fold sums
  r[1] = x2[0];
  r[2] = sy[0];
  r[3] = y2[0];
  r[4] = xy[0];
  for(i = 1 ; i<8 ; i++){
    r[0] += sx[i];
    r[1] += x2[i];
    r[2] += sy[i];
    r[3] += y2[i];
    r[4] += xy[i];
  }
#endif
}
