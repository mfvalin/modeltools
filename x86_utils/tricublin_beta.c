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

#if defined(TIMING)
uint64_t rdtscp_(void) {   // version "in order" avec "serialization"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
#endif

#if defined(__AVX2__) && defined(__x86_64__)
static double cp133 =  1.0/3.0;
#endif
static double cp167 =  1.0/6.0;
static double cm167 = -1.0/6.0;
static double cm5 = -0.5;
static double cp5 = 0.5;
static double one = 1.0;
static double two = 2.0;

// Fortran dimensions: f(NI,NJ,NK)
// NI   : length of a line
// NINJ : length of a plane
// f : address of lower left corner of the 4 x 4 x 4 or 4 x 4 x 2 box
// interpolation is done along z, then along y, then along x
// values are interpolated using the cubic/linear polynomial coefficients
// pz(1) and pz(4) both == 0.0 mean linear interpolation along z
// pz(2) = 1.0 - dz, pz(3) = dz are expected in the linear case (0.0 <= dz <= 1.0)
// the above might get adjusted in the future
// in the z linear case, planes 0 and 1 are the same as are planes 2 and 3
void Tricublin_zyxf_beta(float *d, float *f1, double *px, double *py, double *pz, int NI, int NINJ){
  int ni = NI;
  int ninj = NINJ;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  int ni2, ni3;
#if defined(__AVX2__) && defined(__x86_64__)
  float *s1;
  int32_t *l = (int32_t *)d;
  __m256d cz0, cz1, cz2, cz3, cy, cx;
  __m256d za;
  __m256d ya;
  __m128d va;
  __m128  ta;
#else
  double va4[4], vb4[4], vc4[4], vd4[4], dst[4];
  int i, ninj2, ninj3;
#endif
  
  ni2 = ni + ni;
  ni3 = ni2 + ni;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(pz[0] == 0.0 && pz[3] == 0.0) ninjl = 0;

#if defined(__AVX2__) && defined(__x86_64__)
  // ==== interpolation along Z, vector length is 12 (3 vectors of length 4 per plane) ====
  cz0 = _mm256_broadcast_sd(pz  );   // coefficients for interpolation along z, promote to vectors
  cz1 = _mm256_broadcast_sd(pz+1);
  cz2 = _mm256_broadcast_sd(pz+2);
  cz3 = _mm256_broadcast_sd(pz+3);

  s1 = f1;      // row 0, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_mul_pd(za,cy);

  s1 = f1 + ni;        // row 1, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+1);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_fmadd_pd(za,cy,ya);

  s1 = f1 + ni2;     // row 2, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+2);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_fmadd_pd(za,cy,ya);

  s1 = f1 + ni3;     // row 3, 4 planes (Z0, Z1, Z2, Z3)
  cy  = _mm256_broadcast_sd(py+3);

  za  = _mm256_mul_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)), cz0);           // plane 0
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz1,za);       // plane 1
  s1 += ninj;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz2,za);       // plane 2
  s1 += ninjl;
  za  = _mm256_fmadd_pd( _mm256_cvtps_pd(_mm_loadu_ps(s1)),cz3,za);       // plane 3

  ya  = _mm256_fmadd_pd(za,cy,ya);

  // interpolation along x
  cx = _mm256_loadu_pd(px) ;     // make sure to use "unaligned" load to avoid segfault
  ya = _mm256_mul_pd(ya, cx );   // multiply by px
  // use folding to sum 4 terms
  va = _mm_add_pd( _mm256_extractf128_pd(ya,0) , _mm256_extractf128_pd(ya,1) ); va = _mm_hadd_pd(va,va);
  ta = _mm_cvtpd_ps(va); l[0] = _mm_extract_epi32((__m128i) ta, 0);  // convert results to float and store
#else
  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  // field 1
  for (i=0 ; i<4 ; i++){
    va4[i] = f1[i    ]*pz[0] + f1[i    +ninjl]*pz[1] +  f1[i    +ninj2]*pz[2] + f1[i    +ninj3]*pz[3];
    vb4[i] = f1[i+ni ]*pz[0] + f1[i+ni +ninjl]*pz[1] +  f1[i+ni +ninj2]*pz[2] + f1[i+ni +ninj3]*pz[3];
    vc4[i] = f1[i+ni2]*pz[0] + f1[i+ni2+ninjl]*pz[1] +  f1[i+ni2+ninj2]*pz[2] + f1[i+ni2+ninj3]*pz[3];
    vd4[i] = f1[i+ni3]*pz[0] + f1[i+ni3+ninjl]*pz[1] +  f1[i+ni3+ninj2]*pz[2] + f1[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[1]*px[1] + dst[2]*px[2] + dst[3]*px[3];
#endif
}

void Tricubic_coeffs_d(double *px, double *py, double *pz, double x, double y, double z){

  pz[0] = cm167*z*(z-one)*(z-two);        // coefficients for interpolation along z
  pz[1] = cp5*(z+one)*(z-one)*(z-two);
  pz[2] = cm5*z*(z+one)*(z-two);
  pz[3] = cp167*z*(z+one)*(z-one);
  pz[4] = 0.0;
  pz[5] = 1.0 - z;
  pz[6] = z;
  pz[7] = 0.0;

  py[0] = cm167*y*(y-one)*(y-two);        // coefficients for interpolation along y
  py[1] = cp5*(y+one)*(y-one)*(y-two);
  py[2] = cm5*y*(y+one)*(y-two);
  py[3] = cp167*y*(y+one)*(y-one);
  py[4] = 0.0;
  py[5] = 1.0 - y;
  py[6] = y;
  py[7] = 0.0;

  px[0] = cm167*x*(x-one)*(x-two);        // coefficients for interpolation along x
  px[1] = cp5*(x+one)*(x-one)*(x-two);
  px[2] = cm5*x*(x+one)*(x-two);
  px[3] = cp167*x*(x+one)*(x-one);
  px[4] = 0.0;
  px[5] = 1.0 - x;
  px[6] = x;
  px[7] = 0.0;
}

#if defined(SELF_TEST)
#define NI 300
#define NJ 200
#define NK 85
#include <stdio.h>
#include <sys/time.h>

int main(int argc, char **argv){
  float array[3*NI*NJ*NK];
  float a1[NI*NJ*NK];
  float a2[NI*NJ*NK];
  float a3[NI*NJ*NK];
  float dest[3*NI*NJ*NK];
  float mi[3];
  float ma[3];
  float lin[3];
  double pz[16], px[16], py[16];
  float avg=0.0;
  double dx, dy, dz;
  double expected1, expected2;
  struct timeval t1, t2;
  long long tm1, tm2;
  int i,j,k,ijk;
  float *p = array;
  float *p1 = a1;
  float *p2 = a2;
  float *p3 = a3;
  int nijk = NI*NJ*NK - 5*NI*NJ;
  int II, JJ, KK;

  for (k=0 ; k<NK ; k++ ){
    for (j=0 ; j<NJ ; j++ ){
      for (i=0 ; i<NI ; i++ ){
        p[0] = i + j + k;
        p[1] = p[0] + 10;
        p[2] = p[0] + 100;
        *p1++ = p[0];
        *p2++ = p[1];
        *p3++ = p[2];
        p += 3;
      }
    }
  }
  dz = .333;
  dx = .555;
  dy = .777;
  II = 10;
  JJ = 10;
  KK = 10;
  Tricubic_coeffs_d(px, py, pz, dx, dy, dz);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  Tricublin_zyxf_beta(&dest[0], &a1[II + JJ*NI + KK*NI*NJ], px, py, pz, NI, NI*NJ);
  Tricublin_zyxf_beta(&dest[1], &a2[II + JJ*NI + KK*NI*NJ], px, py, pz, NI, NI*NJ);
  Tricublin_zyxf_beta(&dest[2], &a3[II + JJ*NI + KK*NI*NJ], px, py, pz, NI, NI*NJ);
  expected2 = II + 1 + dx + JJ + 1 + dy + KK + 1 + dz ;
  printf("Cubic case  :\ngot %10.6f, %10.6f, %10.6f, expected %10.6f %10.6f %10.6f\n",dest[0],dest[1],dest[2],expected2,expected2+10,expected2+100);

  dest[0] = -1.0 ; dest[1] = -1.0 ; dest[2] = -1.0 ;
  Tricublin_zyxf_beta(&dest[0], &a1[II + JJ*NI + KK*NI*NJ], px, py, pz+4, NI, NI*NJ);
  Tricublin_zyxf_beta(&dest[1], &a2[II + JJ*NI + KK*NI*NJ], px, py, pz+4, NI, NI*NJ);
  Tricublin_zyxf_beta(&dest[2], &a3[II + JJ*NI + KK*NI*NJ], px, py, pz+4, NI, NI*NJ);
  expected2 = II + 1 + dx + JJ + 1 + dy + KK + 0 + dz ;
  printf("Linear case :\ngot %10.6f, %10.6f, %10.6f, expected %10.6f %10.6f %10.6f\n",dest[0],dest[1],dest[2],expected2,expected2+10,expected2+100);

}
#endif
