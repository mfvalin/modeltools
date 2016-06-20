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

static int ni=0;
static int ninj=0;
static float cp167 =  0.1666666667;
static float cm167 = -0.1666666667;
static float cp5 = .5;
static float cm5 = -.5;
static float one = 1.0;
static float two = 2.0;

void setninj(int ini, int ininj) { ni = ini; ninj = ininj;}

float tricub_x86_f(float *src, float *abcd, float x, float y){
  float *s;
  float x0, x1, x2, x3, y0, y1, y2, y3;
  float dst[4];
#if defined(__AVX2__) && defined(__x86_64__)
  __m256 v1, v2, v3, v4;
  __m256 va, vb, vc, vd;
  __m128 va4, vb4, vc4, vd4;
  __m128 v128a, v128b;
  __m128 vy0, vy1, vy2, vy3;
#else
  int i, ni2, ni3, ninj2, ninj3;
  float va4[4], vb4[4], vc4[4], vd4[4];
  ninj2 = ninj + ninj;
  ninj3 = ninj2 + ninj;
  ni2 = ni + ni;
  ni3 = ni2 + ni;
#endif
#if defined(__AVX2__) && defined(__x86_64__)

// ==== interpolation along Z, vector length is 16 (2 vectors of length 8 per plane) ====

  va = _mm256_broadcast_ss(abcd);   // promote constants to vectors
  vb = _mm256_broadcast_ss(abcd+1);
  vc = _mm256_broadcast_ss(abcd+2);
  vd = _mm256_broadcast_ss(abcd+3);

  s = src;                          // rows 0 and 1, 4 planes (Z0, Z1, Z2, Z3)
  v128a = _mm_loadu_ps(s);          // Z0 row 0
  v1 = _mm256_insertf128_ps(v1,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z0 row 1
  v1 = _mm256_insertf128_ps(v1,v128b,1);
  v1 = _mm256_mul_ps(v1,va);        // v1 = v1*va

  s += ninj;
  v128a = _mm_loadu_ps(s);          // Z1 row 0
  v2 = _mm256_insertf128_ps(v2,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z1 row 1
  v2 = _mm256_insertf128_ps(v2,v128b,1);
  v1 = _mm256_fmadd_ps(v2,vb,v1);   // v1 += v2*vb

  s += ninj;
  v128a = _mm_loadu_ps(s);          // Z2 row 0
  v3 = _mm256_insertf128_ps(v3,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z2 row 1
  v3 = _mm256_insertf128_ps(v3,v128b,1);
  v1 = _mm256_fmadd_ps(v3,vc,v1);   // v1 += v3*vc

  s += ninj;
  v128a = _mm_loadu_ps(s);          // Z3 row 0
  v4 = _mm256_insertf128_ps(v4,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z3 row 1
  v4 = _mm256_insertf128_ps(v4,v128b,1);
  v1 = _mm256_fmadd_ps(v4,vd,v1);   // v1 += v4*vd
                                    // split vector of length 8 into 2 vectors of length 4
  vy0 = _mm256_extractf128_ps(v1,0);// Y0 : row 0 (v1 low)
  vy1 = _mm256_extractf128_ps(v1,1);// Y1 : row 1 (v1 high)

  s = src + 2*ni;                   // rows 2 and 3, 4 planes (Z0, Z1, Z2, Z3)
  v128a = _mm_loadu_ps(s);          // Z0 row 2
  v1 = _mm256_insertf128_ps(v1,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z0 row 3
  v1 = _mm256_insertf128_ps(v1,v128b,1);
  v1 = _mm256_mul_ps(v1,va);        // v1 = v1*va

  s += ninj;
  v128a = _mm_loadu_ps(s);          // Z1 row 2
  v2 = _mm256_insertf128_ps(v2,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z1 row 3
  v2 = _mm256_insertf128_ps(v2,v128b,1);
  v1 = _mm256_fmadd_ps(v2,vb,v1);   // v1 += v2*vb

  s += ninj;
  v128a = _mm_loadu_ps(s);          // Z2 row 2
  v3 = _mm256_insertf128_ps(v3,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z2 row 3
  v3 = _mm256_insertf128_ps(v3,v128b,1);
  v1 = _mm256_fmadd_ps(v3,vc,v1);   // v1 += v3*vc

  s += ninj;
  v128a = _mm_loadu_ps(s);          // Z3 row 2
  v4 = _mm256_insertf128_ps(v4,v128a,0);
  v128b = _mm_loadu_ps(s+ni);       // Z3 row 3
  v4 = _mm256_insertf128_ps(v4,v128b,1);
  v1 = _mm256_fmadd_ps(v4,vd,v1);   // v1 += v4*vd
                                    // split vector of length 8 into 2 vectors of length 4
  vy2 = _mm256_extractf128_ps(v1,0);// Y2 : row 2  (v1 low)
  vy3 = _mm256_extractf128_ps(v1,1);// Y3 : row 3  (v1 high)

// ==== interpolation along Y, vector length is 4 (4 rows) ====

  y0 = cm167*y*(y-one)*(y-two);
  y1 = cp5*(y+one)*(y-one)*(y-two);
  y2 = cm5*y*(y+one)*(y-two);
  y3 = cp167*y*(y+one)*(y-one);

  va4 = _mm_broadcast_ss(&y0);      // promote constants to vectors
  vb4 = _mm_broadcast_ss(&y1);
  vc4 = _mm_broadcast_ss(&y2);
  vd4 = _mm_broadcast_ss(&y3);

  vy0 = _mm_mul_ps(vy0,va4);        //    vy0 * va4
  vy0 = _mm_fmadd_ps(vy1,vb4,vy0);  // += vy1 * vb4
  vy0 = _mm_fmadd_ps(vy2,vc4,vy0);  // += vy2 * vc4
  vy0 = _mm_fmadd_ps(vy3,vd4,vy0);  // += vy3 * vd4
  
  _mm_storeu_ps(dst,vy0);           // store 4 values along X
#else
  y0 = cm167*y*(y-one)*(y-two);
  y1 = cp5*(y+one)*(y-one)*(y-two);
  y2 = cm5*y*(y+one)*(y-two);
  y3 = cp167*y*(y+one)*(y-one);
  for (i=0 ; i<4 ; i++){
    va4[i] = src[i    ]*abcd[0] + src[i    +ninj]*abcd[1] +  src[i    +ninj2]*abcd[2] + src[i    +ninj3]*abcd[3];
    vb4[i] = src[i+ni ]*abcd[0] + src[i+ni +ninj]*abcd[1] +  src[i+ni +ninj2]*abcd[2] + src[i+ni +ninj3]*abcd[3];
    vc4[i] = src[i+ni2]*abcd[0] + src[i+ni2+ninj]*abcd[1] +  src[i+ni2+ninj2]*abcd[2] + src[i+ni2+ninj3]*abcd[3];
    vd4[i] = src[i+ni3]*abcd[0] + src[i+ni3+ninj]*abcd[1] +  src[i+ni3+ninj2]*abcd[2] + src[i+ni3+ninj3]*abcd[3];
    dst[i] = va4[i]*y0 + vb4[i]*y1 + vc4[i]*y2 + vd4[i]*y3;
  }
#endif

// ==== interpolation along x, scalar ====

  x0 = cm167*x*(x-one)*(x-two);
  x1 = cp5*(x+one)*(x-one)*(x-two);
  x2 = cm5*x*(x+one)*(x-two);
  x3 = cp167*x*(x+one)*(x-one);

  return(dst[0]*x0 + dst[1]*x1 + dst[2]*x2 + dst[3]*x3);
}
