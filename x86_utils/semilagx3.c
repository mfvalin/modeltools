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

#if defined(__AVX2__) && defined(__x86_64__)
#include <immintrin.h>
#endif

#include <stdint.h>

static float cp167 =  0.1666666667;
static float cm167 = -0.1666666667;
static float cp5 = 0.5;
static float cm5 = -0.5;
static float one = 1.0;
static float two = 2.0;

// void setninj(int ini, int ininj) {
//   static int ni=0;
//   static int ninj=0;
//   ni = ini; ninj = ininj;}

// Fortran indexed arrays of dimension (3,ni,nj,nk) are expected :
//    element i+1,j,k is 3 reals       from element i,j,k
//    element i,j+1,k is 3*NI reals    from element i,j,k (1 line)
//    element i,j,k+1 is 3*NI*NI reals from element i,j,k (1 plane)
// d     : dimension 3 array to receive result of interpolation
// f     : pointer to lower left corner of 4 x 4 x 4 subarray used for interpolation (CUBIC-CUBIC-CUBIC case)
//         pointer to lower left corner of 4 x 4 x 2 subarray used for interpolation (CUBIC-CUBIC-LINEAR case)
// abcd  : coefficients to use for the interpolation polynomial along z (dimension(4))
// x     : fractional index along x in middle cell (0 <= x <= 1.0)
// y     : fractional index along y in middle cell (0 <= y <= 1.0)
// NI    : number of real triplets on a line
// NJ    : number of lines in a plane
// NOTE:  in the linear case, a[2] and a[3] MUST BE 0.0, a[0] = 1 - z, a[1] = z , where z is the
//        fractional index along z in middle cell (0 <= z <= 1.0)
#if defined(NEVER_EVER_TRUE)
interface
  subroutine tricub_x86_f3(dd, ff, abcd, x, y, NI, NJ) bind(C,name='tricub_x86_f3')
    import :: C_FLOAT, C_INT
    real(C_FLOAT),  intent(OUT), dimension(3)   :: dd
    real(C_FLOAT),  intent(IN),  dimension(3,*) :: ff
    real(C_FLOAT),  intent(IN),  dimension(4)   :: abcd
    real(C_FLOAT),  intent(IN),  value :: x, y
    integer(C_INT), intent(IN),  value :: NI, NJ
  end subroutine tricub_x86_f3
  subroutine tricub_x86_3f(dd, ff1, f2, ff3, abcd, x, y, NI, NJ) bind(C,name='tricub_x86_3f')
    import :: C_FLOAT, C_INT
    real(C_FLOAT),  intent(OUT), dimension(3)   :: dd
    real(C_FLOAT),  intent(IN),  dimension(*)   :: ff1, ff2, ff3
    real(C_FLOAT),  intent(IN),  dimension(4)   :: abcd
    real(C_FLOAT),  intent(IN),  value :: x, y
    integer(C_INT), intent(IN),  value :: NI, NJ
  end subroutine tricub_x86_3f
end interface
#endif
void tricub_x86_f3(void *dd, void *ff, float *abcd, float x, float y, int NI, int NJ){
  float *s;
  float x0, x1, x2, x3, y0, y1, y2, y3;
  float dst[13];
  int ni = 3*NI;
  int ninj = ni*NJ;
  int ninjl;
  int *d = (int *) dd;
  float *f = (float *) ff;

#if defined(__AVX2__) && defined(__x86_64__)
  __m256 cz0, cz1, cz2, cz3, cya;
  __m256 za, ya, ta;
  __m256 zb, yb, tb;
  __m128 vx0, vx1, vx2, vx3;
  __m128 cx0, cx1, cx2, cx3, cxt;
#else
  float va4[4], vb4[4], vc4[4], vd4[4];
  int i, ni2, ni3, ninj2, ninj3;
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

#if defined(__AVX2__) && defined(__x86_64__)

// ==== interpolation along Z, vector length is 16 (2 vectors of length 8 per plane) ====

  cz0 = _mm256_broadcast_ss(abcd);   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm256_broadcast_ss(abcd+1);
  cz2 = _mm256_broadcast_ss(abcd+2);
  cz3 = _mm256_broadcast_ss(abcd+3);

  ya =_mm256_xor_ps(ya,ya)  ; yb = _mm256_xor_ps(yb,yb);    // set y accumulator to zero

  s = f;                          // row 0, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set z accumulator to zero
  cya = _mm256_broadcast_ss(&y0);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set z accumulator to zero
  cya = _mm256_broadcast_ss(&y1);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + 2*ni;                          // row 2, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set z accumulator to zero
  cya = _mm256_broadcast_ss(&y2);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

  s = f + 3*ni;                          // row 3, 4 planes (Z0, Z1, Z2, Z3)
  za  = _mm256_xor_ps(za,za); zb = _mm256_xor_ps(zb,zb);    // set z accumulator to zero
  cya = _mm256_broadcast_ss(&y3);      // promote constant to vector

  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 0
  za  = _mm256_fmadd_ps(ta,cz0,za); zb  = _mm256_fmadd_ps(tb,cz0,zb);

  s += ninj;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 1
  za  = _mm256_fmadd_ps(ta,cz1,za); zb  = _mm256_fmadd_ps(tb,cz1,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 2
  za  = _mm256_fmadd_ps(ta,cz2,za); zb  = _mm256_fmadd_ps(tb,cz2,zb);

  s += ninjl;
  ta  = _mm256_load_ps(s)   ; tb = _mm256_load_ps(s+4);  // plane 3
  za  = _mm256_fmadd_ps(ta,cz3,za); zb  = _mm256_fmadd_ps(tb,cz3,zb);

  ya  = _mm256_fmadd_ps(za,cya,ya); yb  = _mm256_fmadd_ps(zb,cya,yb);    // accumulation along y

// dst : x(0,0) x(1,0) x(2,0)   x(0,1) x(1,1) x(2,1)   x(0,2) x(1,2) x(2,2) x(0,3)   x(1,3) x(2,3) 0
  _mm256_storeu_ps(dst,ya); _mm256_storeu_ps(dst+4,yb); dst[12] = 0.0;  // store result of accumulation along y

// cheating load : only the first 3 values are of any interest in vx0 ... vx3
//          vx0 : x(0,0) x(1,0) x(2,0) x(0,1)   (low part of ya)
//          vx1 : x(0,1) x(1,1) x(2,1) x(0,2)
//          vx2 : x(0,2) x(1,2) x(2,2) x(0,3)
//          vx3 : x(0,3) x(1,3) x(2,3) 0   ( instead of load, could have used a shuffle on yb)
  vx0 = _mm256_extractf128_ps(ya,0); cxt =  _mm_broadcast_ss(&x0); vx0 = _mm_mul_ps(vx0,cxt);
  vx1 = _mm_loadu_ps(dst+3)        ; cxt =  _mm_broadcast_ss(&x1); vx0 = _mm_fmadd_ps(vx1,cxt,vx0);
  vx2 = _mm_loadu_ps(dst+6)        ; cxt =  _mm_broadcast_ss(&x2); vx0 = _mm_fmadd_ps(vx2,cxt,vx0);
  vx3 = _mm_loadu_ps(dst+9)        ; cxt =  _mm_broadcast_ss(&x3); vx0 = _mm_fmadd_ps(vx3,cxt,vx0);
//   _mm_storeu_ps(dst,vx0);    // store 4 values, move only the first 3 into output array d
//   d[0] = dst[0];             // could also use _mm_maskstore_ps (float * mem_addr, __m128i mask, __m128 a)
//   d[1] = dst[1];
//   d[2] = dst[2];
  d[0] = _mm_extract_ps (vx0, 0);  // extract element 0 and store
  d[1] = _mm_extract_ps (vx0, 1);  // extract element 1 and store
  d[2] = _mm_extract_ps (vx0, 2);  // extract element 2 and store
   
#else
  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  for (i=0 ; i<12 ; i++){
    va4[i] = s[i    ]*abcd[0] + s[i    +ninj]*abcd[1] +  s[i    +ninj2]*abcd[2] + s[i    +ninj3]*abcd[3];
    vb4[i] = s[i+ni ]*abcd[0] + s[i+ni +ninj]*abcd[1] +  s[i+ni +ninj2]*abcd[2] + s[i+ni +ninj3]*abcd[3];
    vc4[i] = s[i+ni2]*abcd[0] + s[i+ni2+ninj]*abcd[1] +  s[i+ni2+ninj2]*abcd[2] + s[i+ni2+ninj3]*abcd[3];
    vd4[i] = s[i+ni3]*abcd[0] + s[i+ni3+ninj]*abcd[1] +  s[i+ni3+ninj2]*abcd[2] + s[i+ni3+ninj3]*abcd[3];
    dst[i] = va4[i]*y0 + vb4[i]*y1 + vc4[i]*y2 + vd4[i]*y3;
  }
  d[0] = dst[0]*x0 + dst[3]*x1 + dst[6]*x2 + dst[ 9]*x3;
  d[1] = dst[1]*x0 + dst[4]*x1 + dst[7]*x2 + dst[10]*x3;
  d[2] = dst[2]*x0 + dst[5]*x1 + dst[8]*x2 + dst[11]*x3;
#endif

}
void tricub_x86_3f(void *dd, void *F1, void *F2, void *F3, float *abcd, float x, float y, int NI, int NJ){
  float *s1, *s2, *s3;
  float x0, x1, x2, x3, y0, y1, y2, y3;
  float dst[13];
  int ni = NI;
  int ninj = ni*NJ;
  int ninjl;
  int *d = (int *) dd;
  float *f1 = (float *) F1;
  float *f2 = (float *) F2;
  float *f3 = (float *) F3;

#if defined(__AVX2__) && defined(__x86_64__)
  __m128 cz0, cz1, cz2, cz3, cya;
  __m128 za, ya, ta;
  __m128 zb, yb, tb;
  __m128 zc, yc, tc;
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

#if defined(__AVX2__) && defined(__x86_64__)
// ==== interpolation along Z, vector length is 16 (2 vectors of length 8 per plane) ====
  cz0 = _mm_broadcast_ss(abcd);   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm_broadcast_ss(abcd+1);
  cz2 = _mm_broadcast_ss(abcd+2);
  cz3 = _mm_broadcast_ss(abcd+3);

  ya =_mm_xor_ps(ya,ya) ; yb = _mm_xor_ps(yb,yb) ; yc = _mm_xor_ps(yc,yc) ;    // set y accumulator to zero

  s1 = f1; s2 = f2; s3 = f3;                // row 0, 4 planes (Z0, Z1, Z2, Z3)

  za  = _mm_xor_ps(za,za); zb = _mm_xor_ps(zb,zb); zc = _mm_xor_ps(zc,zc);    // set z accumulator to zero
  cya = _mm_broadcast_ss(&y0);      // promote constant to vector

  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 0
  za  = _mm_fmadd_ps(ta,cz0,za); zb  = _mm_fmadd_ps(tb,cz0,zb); zc  = _mm_fmadd_ps(tc,cz0,zc);

  s1 += ninj; s2 += ninj; s3 += ninj;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 1
  za  = _mm_fmadd_ps(ta,cz1,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 2
  za  = _mm_fmadd_ps(ta,cz2,za); zb  = _mm_fmadd_ps(tb,cz2,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 3
  za  = _mm_fmadd_ps(ta,cz3,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb); yc  = _mm_fmadd_ps(zc,cya,yc);    // accumulation along y

  s1 = f1 + ni; s2 = f2 + ni; s3 = f3 + ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)

  za  = _mm_xor_ps(za,za); zb = _mm_xor_ps(zb,zb); zc = _mm_xor_ps(zc,zc);    // set z accumulator to zero
  cya = _mm_broadcast_ss(&y1);      // promote constant to vector

  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 0
  za  = _mm_fmadd_ps(ta,cz0,za); zb  = _mm_fmadd_ps(tb,cz0,zb); zc  = _mm_fmadd_ps(tc,cz0,zc);

  s1 += ninj; s2 += ninj; s3 += ninj;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 1
  za  = _mm_fmadd_ps(ta,cz1,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 2
  za  = _mm_fmadd_ps(ta,cz2,za); zb  = _mm_fmadd_ps(tb,cz2,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 3
  za  = _mm_fmadd_ps(ta,cz3,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb); yc  = _mm_fmadd_ps(zc,cya,yc);    // accumulation along y

  s1 = f1 + 2*ni; s2 = f2 + 2*ni; s3 = f3 + 2*ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)

  za  = _mm_xor_ps(za,za); zb = _mm_xor_ps(zb,zb); zc = _mm_xor_ps(zc,zc);    // set z accumulator to zero
  cya = _mm_broadcast_ss(&y2);      // promote constant to vector

  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 0
  za  = _mm_fmadd_ps(ta,cz0,za); zb  = _mm_fmadd_ps(tb,cz0,zb); zc  = _mm_fmadd_ps(tc,cz0,zc);

  s1 += ninj; s2 += ninj; s3 += ninj;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 1
  za  = _mm_fmadd_ps(ta,cz1,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 2
  za  = _mm_fmadd_ps(ta,cz2,za); zb  = _mm_fmadd_ps(tb,cz2,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 3
  za  = _mm_fmadd_ps(ta,cz3,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb); yc  = _mm_fmadd_ps(zc,cya,yc);    // accumulation along y

  s1 = f1 + 3*ni; s2 = f2 + 3*ni; s3 = f3 + 3*ni;                          // row 1, 4 planes (Z0, Z1, Z2, Z3)

  za  = _mm_xor_ps(za,za); zb = _mm_xor_ps(zb,zb); zc = _mm_xor_ps(zc,zc);    // set z accumulator to zero
  cya = _mm_broadcast_ss(&y3);      // promote constant to vector

  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 0
  za  = _mm_fmadd_ps(ta,cz0,za); zb  = _mm_fmadd_ps(tb,cz0,zb); zc  = _mm_fmadd_ps(tc,cz0,zc);

  s1 += ninj; s2 += ninj; s3 += ninj;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 1
  za  = _mm_fmadd_ps(ta,cz1,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 2
  za  = _mm_fmadd_ps(ta,cz2,za); zb  = _mm_fmadd_ps(tb,cz2,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  s1 += ninjl; s2 += ninjl; s3 += ninjl;
  ta  = _mm_load_ps(s1)   ;      tb = _mm_load_ps(s2);          tc = _mm_load_ps(s3);  // plane 3
  za  = _mm_fmadd_ps(ta,cz3,za); zb  = _mm_fmadd_ps(tb,cz1,zb); zc  = _mm_fmadd_ps(tc,cz2,zc);

  ya  = _mm_fmadd_ps(za,cya,ya); yb  = _mm_fmadd_ps(zb,cya,yb); yc  = _mm_fmadd_ps(zc,cya,yc);    // accumulation along y

// dst : x(0,0) x(1,0) x(2,0)   x(0,1) x(1,1) x(2,1)   x(0,2) x(1,2) x(2,2)   x(0,3) x(1,3) x(2,3)   0
  _mm_storeu_ps(dst,ya); _mm_storeu_ps(dst+4,yb);_mm_storeu_ps(dst+8,yc); dst[12] = 0.0;  // store result of accumulation along y

// cheating load : only the first 3 values are of any interest in vx0 ... vx3
//          vx0 : x(0,0) x(1,0) x(2,0) x(0,1)   (low part of ya)
//          vx1 : x(0,1) x(1,1) x(2,1) x(0,2)
//          vx2 : x(0,2) x(1,2) x(2,2) x(0,3)
//          vx3 : x(0,3) x(1,3) x(2,3) 0   ( instead of load, could have used a shuffle on yc)
  vx0 = ya                         ; cxt =  _mm_broadcast_ss(&x0); vx0 = _mm_mul_ps(vx0,cxt);
  vx1 = _mm_loadu_ps(dst+3)        ; cxt =  _mm_broadcast_ss(&x1); vx0 = _mm_fmadd_ps(vx1,cxt,vx0);
  vx2 = _mm_loadu_ps(dst+6)        ; cxt =  _mm_broadcast_ss(&x2); vx0 = _mm_fmadd_ps(vx2,cxt,vx0);
  vx3 = _mm_loadu_ps(dst+9)        ; cxt =  _mm_broadcast_ss(&x3); vx0 = _mm_fmadd_ps(vx3,cxt,vx0);
//   _mm_storeu_ps(dst,vx0);    // store 4 values, move only the first 3 into output array d
//   d[0] = dst[0];             // could also use _mm_maskstore_ps (float * mem_addr, __m128i mask, __m128 a)
//   d[1] = dst[1];
//   d[2] = dst[2];
  d[0] = _mm_extract_ps (vx0, 0);  // extract element 0 and store
  d[1] = _mm_extract_ps (vx0, 1);  // extract element 1 and store
  d[2] = _mm_extract_ps (vx0, 2);  // extract element 2 and store
#else
#endif
}

#if defined(SELF_TEST)
#define NI 300
#define NJ 200
#define NK 85
#include <stdio.h>
#include <sys/time.h>

main(){
  float array[NI*NJ*NK];
  float dest[3*NI*NJ*NK];
  float a[4];
  float avg=0.0;
  float dx, dy, dz, expected;
  struct timeval t1, t2;
  long long tm1, tm2;
  int i,j,k,ijk;
  float *p = array;
  int nijk = NI*NJ*NK - 5*NI*NJ;
  for (k=0 ; k<NK ; k++ ){
    for (j=0 ; j<NJ ; j++ ){
      for (i=0 ; i<NI ; i++ ){
	p[0] = i + j + k;
	p[1] = p[0] + 10;
	p[2] = p[0] + 100;
	p += 3;
      }
    }
  }
//   for (i=0 ; i<3*NI*NJ*NK ; i++) {
//     array[i] = 1.0;
//   }
//   a[0] = 1.0 ; a[1] = 1.1 ; a[2] = 1.2 ; a[3] = 1.3 ;
//   setninj(3*NI,3*NI*NJ);

  dz = .333;
  dx = .555;
  dy = .777;
#if defined(LINEAR)
  a[0] = 1.0 - dz;
  a[1] = dz;
  a[2] = 0.0;
  a[3] = 0.0;
#else
  a[0] = cm167*dz*(dz-one)*(dz-two);        // coefficients for interpolation along z
  a[1] = cp5*(dz+one)*(dz-one)*(dz-two);
  a[2] = cm5*dz*(dz+one)*(dz-two);
  a[3] = cp167*dz*(dz+one)*(dz-one);
#endif
#if defined(LINEAR)
  expected = 1.0 + dx + 1.0 + dy + 0.0 + dz ;
#else
  expected = 1.0 + dx + 1.0 + dy + 1.0 + dz ;
#endif
  tricub_x86_f3(&dest[0],&array[0], &a[0], dx, dy, NI, NJ);
  printf("got %10.6f, %10.6f, %10.6f, expected %10.6f\n",dest[0],dest[1],dest[2],expected);

#if defined(LINEAR)
  expected = 1.0 + dx + 2.0 + dy + 0.0 + dz ;
#else
  expected = 1.0 + dx + 2.0 + dy + 1.0 + dz ;
#endif
  tricub_x86_f3(&dest[0],&array[3*NI], &a[0], dx, dy, NI, NJ);
  printf("got %10.6f, %10.6f, %10.6f, expected %10.6f\n",dest[0],dest[1],dest[2],expected);

#if defined(LINEAR)
  expected = 2.0 + dx + 2.0 + dy + 1.0 + dz ;
#else
  expected = 2.0 + dx + 2.0 + dy + 2.0 + dz ;
#endif
  tricub_x86_f3(&dest[0],&array[3*NI*NJ + 3 + 3*NI], &a[0], dx, dy, NI, NJ);
  printf("got %10.6f, %10.6f, %10.6f, expected %10.6f\n",dest[0],dest[1],dest[2],expected);

  for (i=0 ; i<nijk ; i++) tricub_x86_f3(&dest[i*3], &array[i], &a[0], .3, .7, NI, NJ);
  gettimeofday(&t1,NULL);
  for (j=0;j<100;j++) {
    for (i=0 ; i<nijk ; i++) tricub_x86_f3(&dest[i*3], &array[i], &a[0], .3, .7, NI, NJ);
  }
  gettimeofday(&t2,NULL);
  tm1 = t1.tv_usec + t1.tv_sec * 1000000;
  tm2 = t2.tv_usec + t2.tv_sec * 1000000;
  tm2 = (tm2 - tm1) / j;
  printf("time per iteration = %ld, npts = %d\n",tm2,3*NI*NJ*(NK-5));
  for (i=0 ; i<nijk ; i++) avg = avg + dest[i];
  printf("%f %f\n",avg,dest[1]);
}
#endif
