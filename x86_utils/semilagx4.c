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
  int i, ni2, ni3, ninj2, ninj3, ninjl;
#if defined(__AVX2__) && defined(__x86_64__)
  __m256 cz0, cz1, cz2, cz3, cya;
  __m256 za, ya, ta;
  __m256 zb, yb, tb;
  __m128 vx0, vx1, vx2, vx3;
#else
  float va4[4], vb4[4], vc4[4], vd4[4];
  float dst[16];
#endif
  ni2 = ni + ni;
  ni3 = ni2 + ni;

  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be 0
  if(abcd[2] == 0.0 && abcd[3] == 0.0) ninjl = 0;

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

  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz0,za);   // plane 0
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz0,zb);

  s += ninj;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz1,za);   // plane 1
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz1,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz2,za);   // plane 2
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz2,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz3,za);   // plane 3
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y1);      // promote constant to vector

  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz0,za);   // plane 0
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz0,zb);

  s += ninj;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz1,za);   // plane 1
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz1,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz2,za);   // plane 2
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz2,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz3,za);   // plane 3
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + 2*ni;                          // row 2, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y2);      // promote constant to vector

  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz0,za);   // plane 0
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz0,zb);

  s += ninj;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz1,za);   // plane 1
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz1,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz2,za);   // plane 2
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz2,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz3,za);   // plane 3
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + 3*ni;                          // row 3, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set to zero
  cya = _mm256_broadcast_ss(&y3);      // promote constant to vector

  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz0,za);   // plane 0
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz0,zb);

  s += ninj;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz1,za);   // plane 1
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz1,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz2,za);   // plane 2
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz2,zb);

  s += ninjl;
  za  = _mm256_fmadd_ps(_mm256_load_ps(s  ),cz3,za);   // plane 3
  zb  = _mm256_fmadd_ps(_mm256_load_ps(s+8),cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  vx0 =   _mm_mul_ps(_mm_broadcast_ss(&x0),_mm256_extractf128_ps(ya,0));      //    x0 * ya(lo)  (x0 * column 0)
  vx0 = _mm_fmadd_ps(_mm_broadcast_ss(&x1),_mm256_extractf128_ps(ya,1),vx0);  // += x1 * ya(hi)  (x1 * column 1)
  vx0 = _mm_fmadd_ps(_mm_broadcast_ss(&x2),_mm256_extractf128_ps(yb,0),vx0);  // += x2 * yb(lo)  (x2 * column 2)
  vx0 = _mm_fmadd_ps(_mm_broadcast_ss(&x3),_mm256_extractf128_ps(yb,1),vx0);  // += x3 * yb(hi)  (x3 * column 3)

  _mm_storeu_ps(d,vx0);   // store 4 results  (column 0)
   
#else
  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  for (i=0 ; i<16 ; i++){
    va4[i] = s[i    ]*abcd[0] + s[i    +ninj]*abcd[1] +  s[i    +ninj2]*abcd[2] + s[i    +ninj3]*abcd[3];
    vb4[i] = s[i+ni ]*abcd[0] + s[i+ni +ninj]*abcd[1] +  s[i+ni +ninj2]*abcd[2] + s[i+ni +ninj3]*abcd[3];
    vc4[i] = s[i+ni2]*abcd[0] + s[i+ni2+ninj]*abcd[1] +  s[i+ni2+ninj2]*abcd[2] + s[i+ni2+ninj3]*abcd[3];
    vd4[i] = s[i+ni3]*abcd[0] + s[i+ni3+ninj]*abcd[1] +  s[i+ni3+ninj2]*abcd[2] + s[i+ni3+ninj3]*abcd[3];
    dst[i] = va4[i]*y0 + vb4[i]*y1 + vc4[i]*y2 + vd4[i]*y3;
  }
  d[0] = dst[0]*x0 + dst[4]*x1 + dst[ 8]*x2 + dst[12]*x3;
  d[1] = dst[1]*x0 + dst[5]*x1 + dst[ 9]*x2 + dst[13]*x3;
  d[2] = dst[2]*x0 + dst[6]*x1 + dst[10]*x2 + dst[14]*x3;
  d[3] = dst[3]*x0 + dst[7]*x1 + dst[11]*x2 + dst[15]*x3;
#endif

}

void tricub_x86_4f(void *dd, void *F1, void *F2, void *F3, void *F4, float *abcd, float x, float y, int NI, int NJ){
  float *s1, *s2, *s3, *s4;
  float x0, x1, x2, x3, y0, y1, y2, y3;
  float dst[13];
  int ni = NI;
  int ninj = ni*NJ;
  int ninjl;
  float *d = (float *) dd;
  float *f1 = (float *) F1;
  float *f2 = (float *) F2;
  float *f3 = (float *) F3;
  float *f4 = (float *) F4;
  int ni2 = ni + ni;
  int ni3 = ni2 + ni;

#if defined(__AVX__) && defined(__x86_64__)
  __m128 cz0, cz1, cz2, cz3, cz4, cya;
  __m128 za, ya, ta;
  __m128 zb, yb, tb;
  __m128 zc, yc, tc;
  __m128 zd, yd, td;
  __m128 vx0, vx1, vx2, vx3;
  __m128 cx0, cx1, cx2, cx3, cxt;
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

  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be 0
  if(abcd[2] == 0.0 && abcd[3] == 0.0) ninjl = 0;

#if defined(__AVX__) && defined(__x86_64__)
  // ==== interpolation along Z, vector length is 12 (3 vectors of length 4 per plane) ====
  cz0 = _mm_broadcast_ss(abcd);   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm_broadcast_ss(abcd+1);
  cz2 = _mm_broadcast_ss(abcd+2);
  cz3 = _mm_broadcast_ss(abcd+3);

  s1 = f1; s2 = f2; s3 = f3; s4 = f4;                // row 0, 4 planes (Z0, Z1, Z2, Z3)

  cya = _mm_broadcast_ss(&y0);      // promote constant to vector

  za  = _mm_mul_ps(_mm_load_ps(s1),cz0); zb  = _mm_mul_ps(_mm_load_ps(s2),cz0);             // plane 0
  zc  = _mm_mul_ps(_mm_load_ps(s3),cz0); zd  = _mm_mul_ps(_mm_load_ps(s4),cz0);

  s1 += ninj; s2 += ninj; s3 += ninj; s4 += ninj;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz1,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz1,zb);   // plane 1
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz1,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz1,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz2,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz2,zb);   // plane 2
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz2,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz2,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz3,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz3,zb);   // plane 3
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz3,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz3,zd);

  ya  = _mm_mul_ps(za,cya); yb  = _mm_mul_ps(zb,cya);               // accumulation along y
  yc  = _mm_mul_ps(zc,cya); yd  = _mm_mul_ps(zd,cya);

  s1 = f1 + ni; s2 = f2 + ni; s3 = f3 + ni; s4 = f4 + ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)

  cya = _mm_broadcast_ss(&y1);      // promote constant to vector

  za  = _mm_mul_ps(_mm_load_ps(s1),cz0); zb  = _mm_mul_ps(_mm_load_ps(s2),cz0);             // plane 0
  zc  = _mm_mul_ps(_mm_load_ps(s3),cz0); zd  = _mm_mul_ps(_mm_load_ps(s4),cz0);

  s1 += ninj; s2 += ninj; s3 += ninj; s4 += ninj;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz1,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz1,zb);   // plane 1
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz1,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz1,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz2,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz2,zb);   // plane 2
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz2,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz2,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz3,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz3,zb);   // plane 3
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz3,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz3,zd);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb);     // accumulation along y
  yc  = _mm_fmadd_ps(zc,cya,yc); yd  = _mm_fmadd_ps(zd,cya,yd);

  s1 = f1 + ni2; s2 = f2 + ni2; s3 = f3 + ni2; s4 = f4 + ni2;                          // row 2, 4 planes (Z0, Z1, Z2, Z3)

  cya = _mm_broadcast_ss(&y2);      // promote constant to vector

  za  = _mm_mul_ps(_mm_load_ps(s1),cz0); zb  = _mm_mul_ps(_mm_load_ps(s2),cz0);             // plane 0
  zc  = _mm_mul_ps(_mm_load_ps(s3),cz0); zd  = _mm_mul_ps(_mm_load_ps(s4),cz0);

  s1 += ninj; s2 += ninj; s3 += ninj; s4 += ninj;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz1,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz1,zb);   // plane 1
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz1,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz1,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz2,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz2,zb);   // plane 2
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz2,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz2,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz3,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz3,zb);   // plane 3
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz3,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz3,zd);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb);     // accumulation along y
  yc  = _mm_fmadd_ps(zc,cya,yc); yd  = _mm_fmadd_ps(zd,cya,yd);

  s1 = f1 + ni3; s2 = f2 + ni3; s3 = f3 + ni3; s4 = f4 + ni3;                 // row 3, 4 planes (Z0, Z1, Z2, Z3)

  cya = _mm_broadcast_ss(&y3);      // promote constant to vector

  za  = _mm_mul_ps(_mm_load_ps(s1),cz0); zb  = _mm_mul_ps(_mm_load_ps(s2),cz0);             // plane 0
  zc  = _mm_mul_ps(_mm_load_ps(s3),cz0); zd  = _mm_mul_ps(_mm_load_ps(s4),cz0);

  s1 += ninj; s2 += ninj; s3 += ninj; s4 += ninj;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz1,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz1,zb);   // plane 1
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz1,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz1,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz2,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz2,zb);   // plane 2
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz2,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz2,zd);

  s1 += ninjl; s2 += ninjl; s3 += ninjl; s4 += ninjl;
  za  = _mm_fmadd_ps(_mm_load_ps(s1),cz3,za); zb  = _mm_fmadd_ps(_mm_load_ps(s2),cz3,zb);   // plane 3
  zc  = _mm_fmadd_ps(_mm_load_ps(s3),cz3,zc); zd  = _mm_fmadd_ps(_mm_load_ps(s4),cz3,zd);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb);     // accumulation along y
  yc  = _mm_fmadd_ps(zc,cya,yc); yd  = _mm_fmadd_ps(zd,cya,yd);

  vx0 =   _mm_mul_ps(ya,_mm_broadcast_ss(&x0));       // accumulation along x
  vx0 = _mm_fmadd_ps(yb,_mm_broadcast_ss(&x1),vx0);
  vx0 = _mm_fmadd_ps(yc,_mm_broadcast_ss(&x2),vx0);
  vx0 = _mm_fmadd_ps(yd,_mm_broadcast_ss(&x3),vx0);

  _mm_storeu_ps(d,vx0);   // store 4 results
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
