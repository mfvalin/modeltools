#if ! defined(F_TEST)
/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2017  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
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
static float cp133 =  0.166666666666666667E0;
static float cm133 = -0.166666666666666667E0;
static float cp5 =  .5;
static float cm5 = -.5;
static float one = 1.0;
static float two = 2.0;

#if defined(TIMING)
#include <stdio.h>
#include <stdint.h>
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}
#endif
#include <immintrin.h>
/*
           interpolate a column from a 3D source array f, put results in array r
  
   f       3D source array, fortran dimension(ni,nj,nk),  ni*nj = ninj
   ni      distance between f(i,j,k) and f(i,j+1,k)
   ninj    distance between f(i,j,k) and f(i,j,k+1)
   r       2D array, fortran dimension(np,nk)
   nk      number of levels
   xx      i coordinate in i j fractional index space of desired column
   yy      j coordinate in i j fractional index space of desired column
  
   f is assumed to point to f(1,1,1)
  
   xx = 2.5, yy = 2.5 would be the center ot the square formed by
   f(2,2,k) f(3,2,k) f(2,3,k) f(3.3.k)  (where 1 <= k <= nk)
  
   to call from FORTRAN, the following interface is needed
  interface
    subroutine int_yinyang_cub_yx(f, r, ni, ninj, nk, np, x, y) bind(C,name='int_yinyang_cub_yx')
      import :: C_INT, C_FLOAT, C_DOUBLE
      real(C_FLOAT), dimension(*), intent(IN) :: f
      real(C_FLOAT), dimension(*), intent(OUT) :: r
      real(C_DOUBLE), intent(IN), value :: x, y
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np
    end subroutine int_yinyang_cub_yx
  end interface
 */
void int_yinyang_cub_yx(float *f, float *r, int ni, int ninj, int nk, int np, double xx, double yy){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  fr0, fr1, fr2, fr3, frt ;
  __m128d ft0, ft1;
#else
  double fd0[4], fd1[4], fd2[4], fd3[0] ;
#endif
  double  wx[4], wy[4] ;
  double  x, y;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
// printf("DEBUG: f = %p, r = %p \n",f,r);
// printf("DEBUG: f[0] = %f, r[0] = %f, ni = %d, ninj = %d, nk = %d, np = %d, xx = %f, yy = %f\n",f[0],r[0],ni,ninj,nk,np,xx,yy);
  x = xx - 1.0 ; y = yy - 1.0; // xx and yy are in "ORIGIN 1"
  ix = xx ; ix = ix - 1;   // xx and yy are in "ORIGIN 1"
  iy = yy ; iy = iy - 1;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;
//   printf("DEBUG: ix=%d, iy=%d, ni=%d\n",ix, iy, ni);
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  // fetch 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    // interpolation along J level k, prefetch 4 rows for level k+1
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
  // interpolation along j , level nk
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // interpolation along i: multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    f+= ninj;
    r += np;
  }
#endif
}

#if defined(MONO)
void int_yinyang_cub_yx_mono(float *f, float *r, int ni, int ninj, int nk, int np, double xx, double yy){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  fr0, fr1, fr2, fr3, frt, fmi, fma ;
  __m128d ft0, ft1;
#else
  double fd0[4], fd1[4], fd2[4], fd3[0] ;
#endif
  double  wx[4], wy[4] ;
  double  x, y;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
// printf("DEBUG: f = %p, r = %p \n",f,r);
// printf("DEBUG: f[0] = %f, r[0] = %f, ni = %d, ninj = %d, nk = %d, np = %d, xx = %f, yy = %f\n",f[0],r[0],ni,ninj,nk,np,xx,yy);
  x = xx - 1.0 ; y = yy - 1.0; // xx and yy are in "ORIGIN 1"
  ix = xx ; ix = ix - 2;   // xx and yy are in "ORIGIN 1"
  iy = yy ; iy = iy - 2;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  // fetch 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,0]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,0]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,0]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,0]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;                               // point to next level
    // interpolation along J level k, prefetch 4 rows for level k+1
    fmi = _mm_min_ps(fr0,fr1) ;             // min of first 2 rows
    fma = _mm_max_ps(fr0,fr1) ;             // max of first 2 rows
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k+1]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fmi = _mm_min_ps(fmi,fr2) ;             // min of first 3 rows
    fma = _mm_max_ps(fma,fr2) ;             // max of first 3 rows
    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k+1]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fmi = _mm_min_ps(fmi,fr3) ;             // min of the 4 rows
    fma = _mm_max_ps(fma,fr3) ;             // max of the 4 rows
    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k+1]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    // get minimum of 4 vector elements
    frt = _mm_permute_ps(fmi,0xEE) ;        // fmi[2]              fmi[3]             fmi[2]  fmi[3] 
    fmi = _mm_min_ps(fmi,frt) ;             // min(fmi[0],fmi[2])  min(fmi[1],fmi[3]) fmi[2]  fmi[3]
    frt = _mm_permute_ps(fmi,0x55) ;        // fmi[1]              fmi[1]             fmi[1]  fmi[1]
    fmi = _mm_min_ss(fmi,frt) ;             // min(fmi[0],fmi[1])

    // get maximum of 4 vector elements
    frt = _mm_permute_ps(fma,0xEE) ;        // fma[2]              fma[3]             fma[2]  fma[3] 
    fma = _mm_max_ps(fma,frt) ;             // max(fma[0],fma[2])  max(fma[1],fma[3]) fma[2]  fma[3]
    frt = _mm_permute_ps(fma,0x55) ;        // fma[1]              fma[1]             fma[1]  fma[1]
    fma = _mm_max_ss(fma,frt) ;             // max(fma[0],fma[1])

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k+1]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x05) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    frt = _mm_max_ss(frt,fmi) ;             // mono min constraint
    frt = _mm_min_ss(frt,fma) ;             // mono max constraint
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
  // interpolation along j , level nk
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fmi = _mm_min_ps(fr0,fr1) ;                // min of first 2 rows
  fma = _mm_max_ps(fr0,fr1) ;                // max of first 2 rows
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fmi = _mm_min_ps(fmi,fr2) ;                // min of first 3 rows
  fma = _mm_max_ps(fma,fr2) ;                // max of first 3 rows
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fmi = _mm_min_ps(fmi,fr3) ;                // min of the 4 rows
  fma = _mm_max_ps(fma,fr3) ;                // max of the 4 rows
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // interpolation along i: multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x05) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]

  // get minimum of 4 vector elements
  frt = _mm_permute_ps(fmi,0xEE) ;        // fmi[2]              fmi[3]             fmi[2]  fmi[3] 
  fmi = _mm_min_ps(fmi,frt) ;             // min(fmi[0],fmi[2])  min(fmi[1],fmi[3]) fmi[2]  fmi[3]
  frt = _mm_permute_ps(fmi,0x55) ;        // fmi[1]              fmi[1]             fmi[1]  fmi[1]
  fmi = _mm_min_ss(fmi,frt) ;             // min(fmi[0],fmi[1])

  // get maximum of 4 vector elements
  frt = _mm_permute_ps(fma,0xEE) ;        // fma[2]              fma[3]             fma[2]  fma[3] 
  fma = _mm_max_ps(fma,frt) ;             // max(fma[0],fma[2])  max(fma[1],fma[3]) fma[2]  fma[3]
  frt = _mm_permute_ps(fma,0x55) ;        // fma[1]              fma[1]             fma[1]  fma[1]
  fma = _mm_max_ss(fma,frt) ;             // max(fma[0],fma[1])

  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
  frt = _mm_max_ss(frt,fmi) ;             // mono min constraint
  frt = _mm_min_ss(frt,fma) ;             // mono max constraint
  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    f+= ninj;
    r += np;
  }
#endif
}
#endif

#if defined(MONO8)
void int_yinyang_cub_yx_mono(float *f, float *r, int ni, int ninj, int nk, int np, double xx, double yy){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt, fmi, fma ;
  __m128  fr0, fr1, fr2, fr3, frt ;   // frt is used as a scalar, fr0->fr3 are 128 bit aliases for fd0->fd3
  __m128d ft0, ft1, smi, sma ;        // smi, sma are used as scalars (reduction of fmi, fma)
  double dd0[4], dd1[4], dd2[4], dd3[0] ;
#else
  double fd0[4], fd1[4], fd2[4], fd3[0] ;
#endif
  double  wx[4], wy[4] ;
  double  x, y, minval, maxval;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
// printf("DEBUG: f = %p, r = %p \n",f,r);
// printf("DEBUG: f[0] = %f, r[0] = %f, ni = %d, ninj = %d, nk = %d, np = %d, xx = %f, yy = %f\n",f[0],r[0],ni,ninj,nk,np,xx,yy);
  x = xx - 1.0 ; y = yy - 1.0; // xx and yy are in "ORIGIN 1"
  ix = xx ; ix = ix - 2;   // xx and yy are in "ORIGIN 1"
  iy = yy ; iy = iy - 2;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;
//   printf("DEBUG: ix=%d, iy=%d, ni=%d\n",ix, iy, ni);
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx  = _mm256_loadu_pd(wx) ;            // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;          // scalar * vector not available,
  fwy1 = _mm256_set1_pd(wy[1]) ;          // promote scalars to vectors
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
  frt  = _mm_xor_ps(frt,frt) ;            // set frt to zero
  // prefetch and promote to double 4 rows, level 1
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
//   _mm256_storeu_pd (dd0,fd0);
//   _mm256_storeu_pd (dd1,fd1);
//   _mm256_storeu_pd (dd2,fd2);
//   _mm256_storeu_pd (dd3,fd3);
//   printf("DEBUG: %f %f %f %f\n",dd3[0],dd3[1],dd3[2],dd3[3]);
//   printf("DEBUG: %f %f %f %f\n",dd2[0],dd2[1],dd2[2],dd2[3]);
//   printf("DEBUG: %f %f %f %f\n",dd1[0],dd1[1],dd1[2],dd1[3]);
//   printf("DEBUG: %f %f %f %f\n",dd0[0],dd0[1],dd0[2],dd0[3]);
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    // interpolation along J, level k,
    // prefetch and promote to double 4 rows, level k+1
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fmi = _mm256_min_pd(fd1,fd2) ;          // min of rows 1 and 2 (ignore 0 and 3)
    fma = _mm256_max_pd(fd1,fd2) ;          // max of rows 1 and 2 (ignore 0 and 3)

    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    // get minimum of fmi vector elements 1 and 2 (ignore 0 and 3)
    ft0 = _mm256_extractf128_pd(fmi,0) ;    // fmi[0]              fmi[1]
    ft1 = _mm256_extractf128_pd(fmi,1) ;    // fmi[2]              fmi[3]
    ft0 = _mm_permute_pd(ft0,0x1)      ;    // fmi[1]              fmi[1]
    smi = _mm_min_sd(ft0,ft1) ;             // min(fmi[1],fmi[2])

    // get maximum of fma vector elements 1 and 2 (ignore 0 and 3)
    ft0 = _mm256_extractf128_pd(fma,0) ;    // fma[0]              fma[1]
    ft1 = _mm256_extractf128_pd(fma,1) ;    // fma[2]              fma[3]
    ft0 = _mm_permute_pd(ft0,0x1) ;         // fma[1]              fma[1]
    sma = _mm_max_sd(ft0,ft1) ;             // max(fma[1],fma[2])

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    ft0 = _mm_max_sd(ft0,smi) ;             // max(result, min value)
    ft0 = _mm_min_sd(ft0,sma) ;             // min(result, max value)
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert result to float
//   _mm_store_sd(&minval,smi) ;
//   _mm_store_sd(&maxval,sma) ;
//   if(k == 0) printf("DEBUG: LIMITS(1) = %f %f\n",minval,maxval);
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
  // interpolation along j , level nk
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fmi = _mm256_min_pd(fd1,fd2) ;             // min of rows 1 and 2 (ignore 0 and 3)
  fma = _mm256_max_pd(fd1,fd2) ;             // max of rows 1 and 2 (ignore 0 and 3)
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;

  // get minimum of fmi vector elements 1 and 2 (ignore 0 and 3)
  ft0 = _mm256_extractf128_pd(fmi,0) ;    // fmi[0]              fmi[1]
  ft1 = _mm256_extractf128_pd(fmi,1) ;    // fmi[2]              fmi[3]
  ft0 = _mm_permute_pd(ft0,0x1)      ;    // fmi[1]              fmi[1]
  smi = _mm_min_sd(ft0,ft1) ;             // min(fmi[1],fmi[2])

  // get maximum of fma vector elements 1 and 2 (ignore 0 and 3)
  ft0 = _mm256_extractf128_pd(fma,0) ;    // fma[0]              fma[1]
  ft1 = _mm256_extractf128_pd(fma,1) ;    // fma[2]              fma[3]
  ft0 = _mm_permute_pd(ft0,0x1) ;         // fma[1]              fma[1]
  sma = _mm_max_sd(ft0,ft1) ;             // max(fma[1],fma[2])

  // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]

  ft0 = _mm_max_sd(ft0,smi) ;             // max(result, min value)
  ft0 = _mm_min_sd(ft0,sma) ;             // min(result, max value)
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert result to float
//   _mm_store_sd(&minval,smi) ;
//   _mm_store_sd(&maxval,sma) ;
//   printf("DEBUG: LIMITS(NK) = %f %f\n",minval,maxval);
  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=0 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    minval = f[ni+1] ;
    maxval = f[ni+1] ;
    minval = (minval > f[ni+2] ) ? f[ni+2] : minval ;
    maxval = (maxval < f[ni+2] ) ? f[ni+2] : maxval ;
    minval = (minval > f[ni2+1] ) ? f[ni2+1] : minval ;
    maxval = (maxval < f[ni2+1] ) ? f[ni2+1] : maxval ;
    minval = (minval > f[ni2+2] ) ? f[ni2+2] : minval ;
    maxval = (maxval < f[ni2+2] ) ? f[ni2+2] : maxval ;
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    r[0] = (r[0] < minval ) ? minval : r[0] ;
    r[0] = (r[0] > maxval ) ? maxval : r[0] ;
    f+= ninj;
    r += np;
  }
#endif
}
#endif

#if defined(MONO_)
void int_yinyang_cub_yx_mono_(float *f, float *r, int ni, int ninj, int nk, int np, double xx, double yy){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3, fdt ;
  __m128  fr0, fr1, fr2, fr3, frt, fmi, fma ;
  __m128d ft0, ft1;
#else
  double fd0[4], fd1[4], fd2[4], fd3[0] ;
  float the_min[4], the_max[4];
  float fmi, fma, rt;
#endif
  double  wx[4], wy[4] ;
  double  x, y;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;
  int ix, iy;
// printf("DEBUG: f = %p, r = %p \n",f,r);
// printf("DEBUG: f[0] = %f, r[0] = %f, ni = %d, ninj = %d, nk = %d, np = %d, xx = %f, yy = %f\n",f[0],r[0],ni,ninj,nk,np,xx,yy);
  x = xx - 1.0 ; y = yy - 1.0; // xx and yy are in "ORIGIN 1"
  ix = x ; ix = ix - 1;   // xx and yy are in "ORIGIN 1"
  iy = y ; iy = iy - 1;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm133*x*(x-one)*(x-two);       // polynomial coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);       // polynomial coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;

  for(k=0 ; k<nk ; k++){
    // fetch 4 rows, level k
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
    f+= ninj;
    // interpolation along J level k, prefetch 4 rows for level k+1
    fmi = _mm_min_ps(fr0,fr1) ;             // min of first 2 rows
    fma = _mm_max_ps(fr0,fr1) ;             // max of first 2 rows
    fdt = _mm256_mul_pd(fd0,fwy0) ;         // sum of row[j] * coefficient[j]

    fmi = _mm_min_ps(fmi,fr2) ;             // min of first 3 rows
    fma = _mm_max_ps(fma,fr2) ;             // max of first 3 rows
    fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;

    fmi = _mm_min_ps(fmi,fr3) ;             // min of the 4 rows
    fma = _mm_max_ps(fma,fr3) ;             // max of the 4 rows
    fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;

    // get minimum of 4 vector elements
    frt = _mm_permute_ps(fmi,0xEE) ;        // fmi[2]              fmi[3]             fmi[2]  fmi[3] 
    fmi = _mm_min_ps(fmi,frt) ;             // min(fmi[0],fmi[2])  min(fmi[1],fmi[3]) fmi[2]  fmi[3]
    frt = _mm_permute_ps(fmi,0x55) ;        // fmi[1]              fmi[1]             fmi[1]  fmi[1]
    fmi = _mm_min_ss(fmi,frt) ;             // min(fmi[0],fmi[1])

    // get maximum of 4 vector elements
    frt = _mm_permute_ps(fma,0xEE) ;        // fma[2]              fma[3]             fma[2]  fma[3] 
    fma = _mm_max_ps(fma,frt) ;             // max(fma[0],fma[2])  max(fma[1],fma[3]) fma[2]  fma[3]
    frt = _mm_permute_ps(fma,0x55) ;        // fma[1]              fma[1]             fma[1]  fma[1]
    fma = _mm_max_ss(fma,frt) ;             // max(fma[0],fma[1])

    fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x05) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fdt[0]+fdt[2]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[0]+fdt[2]+fdt[1]+fdt[3] to float
    frt = _mm_max_ss(frt,fmi) ;             // mono min constraint
    frt = _mm_min_ss(frt,fma) ;             // mono max constraint
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
#else
  for(k=1 ; k<nk ; k++){
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
      the_min[i] = (f[i] < f[i+ni]) ? f[i] ? f[i+ni] ;
      the_max[i] = (f[i] > f[i+ni]) ? f[i] ? f[i+ni] ;
      the_min[i] = (f[i+ni2] < the_min[i]) ? f[i+ni2] : the_min[i]) ;
      the_max[i] = (f[i+ni2] > the_max[i]) ? f[i+ni2] : the_max[i]) ;
      the_min[i] = (f[i+ni3] < the_min[i]) ? f[i+ni3] : the_min[i]) ;
      the_max[i] = (f[i+ni3] > the_max[i]) ? f[i+ni3] : the_max[i]) ;
    }
    fmi = (the_min[0] < the_min[1]) ? the_min[0] : the_min[1] ;
    fma = (the_min[0] > the_min[1]) ? the_min[0] : the_min[1] ;
    fmi = (the_min[2] < fmi) ? the_min[2] : fmi ;
    fma = (the_min[2] > fma) ? the_min[2] : fmi ;
    fmi = (the_min[3] < fmi) ? the_min[3] : fmi ;
    fma = (the_min[3] > fma) ? the_min[3] : fmi ;
    rt = fd0[0] + fd0[1] + fd0[2] + fd0[3];
    rt = (rt < fma) ? rt : fma ;
    rt = (rt > fmi) ? rt : fmi ;
    r[0] = rt;
    f+= ninj;
    r += np;
  }
#endif
}
#endif

#if defined(C_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define NI 65
#define NJ 27
#define NK 80
#define NP 20

int main(int argc,char **argv){
  float f[NK][NJ][NI] ;
  float r[NK][NP];
  double x[NP], y[NP];
  int i, j, k;
  uint64_t t1, t2;

  for(i=0 ; i<NI ; i++){
    for(j=0 ; j<NJ ; j++){
      for(k=0 ; k<NK ; k++){
	f[k][j][i] = i + j + 2;
      }
    }
  }
  for(i=0 ; i<NP ; i++){
    x[i] = 2 + i ;
    y[i] = 2 + i ;
  }
  t1 = rdtsc();
  for(i=0 ; i<NP ; i++){
    int_yinyang_cub_yx(&f[0][0][0], &r[0][i], NI, NI*NJ, NK, NP, x[i], y[i]) ;
  }
  t2 = rdtsc();
  k = t2 - t1;
  printf(" r = %f %f\n",r[0][0],r[0][1]);
  printf("time = %d clocks for %d values, %d flops\n",k,NP*NK,NP*NK*35);
}
#endif
#else
program test_interp
  use ISO_C_BINDING
  implicit none
  integer, parameter :: NI=65
  integer, parameter :: NJ=27
  integer, parameter :: NK=81
  integer, parameter :: NP=4
  integer, parameter :: HX=2
  integer, parameter :: HY=2
  integer, parameter :: NR=25
  real(C_FLOAT), dimension(1-HX:NI+HX , 1-HY:NJ+HY , NK) :: f
  real(C_FLOAT), dimension(NP,NK) :: r
  real(C_DOUBLE), dimension(NP) :: x, y, xmin, xmax
  integer :: i, j, k
  integer :: i0, j0
  integer*8, external :: rdtsc
  integer*8 :: t1, t2, tmg1(NR), tmg2(nr)
  integer :: nidim, ninjdim
  interface
    subroutine int_yinyang_cub_yx(f, r, ni, ninj, nk, np, x, y) bind(C,name='int_yinyang_cub_yx') !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np                                       !InTf!
    end subroutine int_yinyang_cub_yx                                                             !InTf!
    subroutine int_yinyang_cub_yx_mono(f, r, ni, ninj, nk, np, x, y) bind(C,name='int_yinyang_cub_yx_mono') !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np                                       !InTf!
    end subroutine int_yinyang_cub_yx_mono                                                        !InTf!
    subroutine int_yinyang_cub_yx_mono_(f, r, ni, ninj, nk, np, x, y) bind(C,name='int_yinyang_cub_yx_mono_') !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np                                       !InTf!
    end subroutine int_yinyang_cub_yx_mono_                                                       !InTf!
  end interface
#define FXY(A,B,C) (1.3*(A)**3 + 1.4*(A)**2 + (A)*1.5 + 2.3*(B)**3 + 2.4*(B)**2 + (B)*2.5 + (C))

  r = 9999.99
  do k = 1 , NK
    do j = 1-HY , NJ+HY
      do i = 1-HX , NI+HX
!        f(i,j,k) = i + j  + k
        f(i,j,k) = FXY(i*1.0 , j*1.0 , k*1.0)
      enddo
    enddo
!    print *,f(1,1,k),f(2,2,k)
  enddo
  do i = 1 , NP
    x(i) = i - .1
    y(i) = i - .1
  enddo
  nidim = NI + 2*HX
  ninjdim = nidim * (NJ + HY*2)
!  print *,'nidim=',nidim,' , ninjdim=',ninjdim
!  print *,'x=',x,' y=',y
!  print *,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),2),f(nint(x(2)),nint(y(2)),2)
!  print 101,loc(f(1,1,1)), loc(r(1,1))
101 format(2Z17)
!  print *,f(1,1,1), r(1,1)
  do j = 1, NR
    t1 = rdtsc()
    do i = 1 , NP
      call int_yinyang_cub_yx( f(1,1,1), r(i,1), nidim, ninjdim, NK, NP, x(i), y(i) )
    enddo
    t2 = rdtsc()
    tmg1(j) = t2 - t1
!    print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'
    t1 = rdtsc()
    do i = 1 , NP
      call int_yinyang_cub_yx_mono( f(1,1,1), r(i,1), nidim, ninjdim, NK, NP, x(i), y(i) )
    enddo
    t2 = rdtsc()
    tmg2(j) = t2 - t1
!    print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'
  enddo
  print 100,'direct =',tmg1
  print 100,'mono   =',tmg2
  print 100,'flops  =',NP*NK*35
100 format(A,40I6)
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,'X coordinates'
  print 102,x(:)
  print *,'Y coordinates'
  print 102,y(:)
  print *,'F matrix'
  do j = 4, 1-hy, -1
    print 102,f(1-hx:5,j,1)
  enddo
  print *,'MONO: limit (min, max)'
  do i = 1 , NP
    i0 = x(i)
    j0 = y(i)
    xmin(i) = i0 + j0
    xmax(i) = i0 + j0 + 2
  enddo
  print 102,xmin(:) + 1, xmin(:) + NK
  print 102,xmax(:) + 1, xmax(:) + NK
  print *,'MONO: expected'
  print 102,FXY(x(:),y(:),1), FXY(x(:),y(:),NK)
!  print 102,x(:)+y(:)+1, x(:)+y(:)+NK
  print *,' got'
  print 102,r(:,1),r(:,NK)
  print *,' delta'
  print 102,(r(:,1)-FXY(x(:),y(:),1)),(r(:,NK)-FXY(x(:),y(:),NK))
  print 103,(r(:, 1)-FXY(x(:),y(:), 1))  / FXY(x(:),y(:), 1) , & 
            (r(:,NK)-FXY(x(:),y(:),NK))  / FXY(x(:),y(:),NK)

  do i = 1 , NP
    call int_yinyang_cub_yx( f(1,1,1), r(i,1), nidim, ninjdim, NK, NP, x(i), y(i) )
  enddo
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,'DIRECT: expected'
  print 102,FXY(x(:),y(:),1), FXY(x(:),y(:),NK)
  print *,' got'
  print 102,r(:,1),r(:,NK)
  print *,' delta'
  print 102,(r(:,1)-FXY(x(:),y(:),1)),(r(:,NK)-FXY(x(:),y(:),NK))
  print 103,(r(:, 1)-FXY(x(:),y(:), 1))  / FXY(x(:),y(:), 1) , & 
            (r(:,NK)-FXY(x(:),y(:),NK))  / FXY(x(:),y(:),NK)

102 format(16F15.6)
103 format(16E15.2)
end
#endif
