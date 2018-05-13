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
static float cp5 = 0.5;
static float cm5 = -0.5;
static float one = 1.0;
static float two = 2.0;

void setninj(int ini, int ininj) { ni = ini; ninj = ininj;}

void tricub_x86_f4(float *d, float *f, float *abcd, float x, float y){
  float *s;
  float x0, x1, x2, x3, y0, y1, y2, y3;
#if defined(__AVX2__) && defined(__x86_64__)
  __m256 cz0, cz1, cz2, cz3, cya;
  __m256 za, ya, ta;
  __m256 zb, yb, tb;
  __m128 vx0, vx1, vx2, vx3;
  __m128 cxt;
#else
#endif

  y0 = cm167*y*(y-one)*(y-two);        // coefficients for interpolation along y
  y1 = cp5*(y+one)*(y-one)*(y-two);
  y2 = cm5*y*(y+one)*(y-two);
  y3 = cp167*y*(y+one)*(y-one);

  x0 = cm167*x*(x-one)*(x-two);        // coefficients for interpolation along x
  x1 = cp5*(x+one)*(x-one)*(x-two);
  x2 = cm5*x*(x+one)*(x-two);
  x3 = cp167*x*(x+one)*(x-one);

#if defined(__AVX2__) && defined(__x86_64__)

// ==== interpolation along Z, vector length is 16 (2 vectors of length 8 per plane) ====

  cz0 = _mm256_broadcast_ss(abcd);   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm256_broadcast_ss(abcd+1);
  cz2 = _mm256_broadcast_ss(abcd+2);
  cz3 = _mm256_broadcast_ss(abcd+3);

  ya =_mm256_xor_ps(ya,ya)  ; yb = _mm256_xor_ps(yb,yb);    // zero row accumulator

  s = f;                          // row 0, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y0);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y1);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + 2*ni;                          // row 2, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y2);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + 3*ni;                          // row 3, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y3);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+8);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  vx0 = _mm256_extractf128_ps(ya,0) ; cxt =  _mm_broadcast_ss(&x0); vx0 = _mm_mul_ps(vx0,cxt);
  vx1 = _mm256_extractf128_ps(ya,1) ; cxt =  _mm_broadcast_ss(&x1); vx0 = _mm_fmadd_ps(vx1,cxt,vx0);
  vx2 = _mm256_extractf128_ps(yb,0) ; cxt =  _mm_broadcast_ss(&x2); vx0 = _mm_fmadd_ps(vx2,cxt,vx0);
  vx3 = _mm256_extractf128_ps(yb,1) ; cxt =  _mm_broadcast_ss(&x3); vx0 = _mm_fmadd_ps(vx3,cxt,vx0);
  _mm_storeu_ps(d,vx0);
   
#else
#endif

}
#if defined(SELF_TEST)
#define NI 30
#define NJ 20
#define NK 85
main(){
  float array[NI*NJ*NK];
  float dest[NI*NJ*NK];
  float a[4];
  float avg=0.0;
  int i,j;
  int nijk = NI*NJ*NK - 5*NI*NJ;
  for (i=0 ; i<NI*NJ*NK ; i++) array[i] = 1.0;
  a[0] = 1.0 ; a[1] = 1.1 ; a[2] = 1.2 ; a[3] = 1.3 ;
  setninj(NI,NI*NJ);
  for (j=0;j<1000;j++) {
    for (i=0 ; i<nijk ; i++) dest[i] = tricub_x86_f(&array[i], &a[0], .3, .7);
  }
  for (i=0 ; i<nijk ; i++) avg = avg + dest[i];
  printf("%f %f\n",avg,dest[1]);
}
#endif
