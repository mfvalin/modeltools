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
#if defined(NEVER_TRUE)
 to test (Intel compilers) :
 rm -f intrp_biquint_yx intrp_biquint_yx.o intrp_biquint_yx_f.F90
 s.cc -O2 -DTIMING -c -march=core-avx2 intrp_biquint_yx.c
 ln -sf intrp_biquint_yx.c intrp_biquint_yx_f.F90
 s.f90 -no-wrap-margin -DF_TEST -o intrp_biquint_yx intrp_biquint_yx_f.F90 intrp_biquint_yx.o
 ./intrp_biquint_yx
#endif

static double one   = 1.0;
static double two   = 2.0;
static double three = 3.0;
static double cp120 = (double)  0.00833333333333333333E0;
static double cm120 = (double) -0.00833333333333333333E0;
static double cp24  = (double)  0.041666666666666664E0;
static double cm24  = (double) -0.041666666666666664E0;
static double cp12  = (double)  0.0833333333333333333E0;
static double cm12  = (double) -0.0833333333333333333E0;

#if defined(TIMING)
#include <stdio.h>
#include <stdint.h>
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
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
  
   f       3D Fortran indexed source array, real*4, dimension(mini:maxi,minj:maxj,nk)
           ni = (maxi-mini-1)
           ninj = (maxi-mini-1)*(maxj-minj-1)
   ni      distance between f(i,j,k) and f(i,j+1,k)
   ninj    distance between f(i,j,k) and f(i,j,k+1)
   r       2D array, fortran dimension(np,nk)
   nk      number of 2D planes
   xx      i coordinate in i j fractional index space of desired column ( xx(l) is assumed >= mini+1 )
   yy      j coordinate in i j fractional index space of desired column ( yy(l) is assumed >= minj+1 )
  
   f is assumed to point to f(1,1,1)
  
   xx = 2.5, yy = 2.5 would be the center ot the square formed by
   f(2,2,k) f(3,2,k) f(2,3,k) f(3.3.k)  (where 1 <= k <= nk)
  
   to call from FORTRAN, the following interface is used

    subroutine intrp_biquint_yx(f, r, ni, ninj, nk, np, x, y)      bind(C,name='intrp_biquint_yx')
    subroutine intrp_biquint_yx_mono(f, r, ni, ninj, nk, np, x, y) bind(C,name='intrp_biquint_yx_mono')
      real(C_FLOAT), dimension(*), intent(IN) :: f
      real(C_FLOAT), dimension(*), intent(OUT) :: r
      real(C_DOUBLE), intent(IN), value :: x, y
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np
 */
void intrp_biquint_yx(float *f, float *r, int ni, int ninj, int nk, int np, double xx, double yy){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fd4, fd5, fe0, fe1, fe2, fe3, fe4, fe5;
  __m256d fdw, few, fwy, fdt, fet ;
  __m128  fr0, fr1, fr2, fr3, fr4, fr5, frt ;
  __m128d ft0, ft1;
#else
  double fd0[6], fd1[6], fd2[6], fd3[6] ;
#endif
  float *f2;
  double  wx[6], wy[6] ;
  double  x, y;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int ni4 = ni3 + ni;   // + 4 rows
  int ni5 = ni4 + ni;   // + 5 rows
  int i, k;
  int ix, iy;
// printf("DEBUG: f = %p, r = %p \n",f,r);
// printf("DEBUG: f[0] = %f, r[0] = %f, ni = %d, ninj = %d, nk = %d, np = %d, xx = %f, yy = %f\n",f[0],r[0],ni,ninj,nk,np,xx,yy);
  x = xx - 1.0 ; y = yy - 1.0; // xx and yy are in "ORIGIN 1"
  ix = xx ; if(ix > xx) ix = ix -1 ; ix = ix - 2;   // xx and yy are in "ORIGIN 1"
  iy = yy ; if(iy > yy) iy = iy -1 ; iy = iy - 2;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  ix = ix -1 ; iy = iy -1 ; // cubic to quintic addressing
  f = f + ix + iy * ni;
//   printf("DEBUG: ix=%d, iy=%d, ni=%d\n",ix, iy, ni);
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm120*(x+one)*(x)*    (x-one)*(x-two)*(x-three);       // polynomial coefficients along i
  wx[1] = cp24* (x+two)*(x)*    (x-one)*(x-two)*(x-three);
  wx[2] = cm12* (x+two)*(x+one)*(x-one)*(x-two)*(x-three);
  wx[3] = cp12* (x+two)*(x+one)*(x)*    (x-two)*(x-three);
  wx[4] = cm24* (x+two)*(x+one)*(x)*    (x-one)*(x-three);
  wx[5] = cp120*(x+two)*(x+one)*(x)*    (x-one)*(x-  two);

  wy[0] = cm120*(y+one)*(y)*    (y-one)*(y-two)*(y-three);       // polynomial coefficients along j
  wy[1] = cp24* (y+two)*(y)*    (y-one)*(y-two)*(y-three);
  wy[2] = cm12* (y+two)*(y+one)*(y-one)*(y-two)*(y-three);
  wy[3] = cp12* (y+two)*(y+one)*(y)*    (y-two)*(y-three);
  wy[4] = cm24* (y+two)*(y+one)*(y)*    (y-one)*(y-three);
  wy[5] = cp120*(y+two)*(y+one)*(y)*    (y-one)*(y-  two);

#if defined(__AVX2__) && defined(__x86_64__)
  // fetch 6 rows, level 1
  f2 = f + 2;
  fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
  fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr0 = _mm_loadu_ps(f2) ;                // row 1 : f[i+2:i+5 , j   ,k]
  fe0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
  fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
  fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr1 = _mm_loadu_ps(f2+ni) ;             // row 2 : f[i+2:i+5 , j+1 ,k]
  fe1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
  fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
  fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr2 = _mm_loadu_ps(f2+ni2) ;            // row 3 : f[i+2:i+5 , j+2 ,k]
  fe2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
  fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
  fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  fr3 = _mm_loadu_ps(f2+ni3) ;            // row 4 : f[i+2:i+5 , j+3 ,k]
  fe3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
  fr4 = _mm_loadu_ps(f+ni4) ;             // row 4 : f[i:i+3 , j+4 ,k]
  fd4 = _mm256_cvtps_pd(fr4) ;            // promote row 4 to double
  fr4 = _mm_loadu_ps(f2+ni4) ;            // row 4 : f[i+2:i+5 , j+4 ,k]
  fe4 = _mm256_cvtps_pd(fr4) ;            // promote row 4 to double
  fr5 = _mm_loadu_ps(f+ni5) ;             // row 4 : f[i:i+3 , j+5 ,k]
  fd5 = _mm256_cvtps_pd(fr5) ;            // promote row 4 to double
  fr5 = _mm_loadu_ps(f2+ni5) ;            // row 4 : f[i+2:i+5 , j+5 ,k]
  fe5 = _mm256_cvtps_pd(fr5) ;            // promote row 4 to double
  for(k=0 ; k<nk-1 ; k++){
    f+= ninj;
    f2 = f + 2;
    // interpolation along J level k, prefetch k rows for level k+1
    fwy = _mm256_set1_pd(wy[0]) ;
    fdt = _mm256_mul_pd(fd0,fwy) ;          // sum of row[j] * coefficient[j]
    fet = _mm256_mul_pd(fe0,fwy) ;          // sum of row[j] * coefficient[j]
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double
    fr0 = _mm_loadu_ps(f2) ;                // row 1 : f[i+2:i+5 , j   ,k]
    fe0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fwy = _mm256_set1_pd(wy[1]) ;
    fdt = _mm256_fmadd_pd(fd1,fwy,fdt) ;
    fet = _mm256_fmadd_pd(fe1,fwy,fet) ;
    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double
    fr1 = _mm_loadu_ps(f2+ni) ;             // row 2 : f[i+2:i+5 , j+1 ,k]
    fe1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fwy = _mm256_set1_pd(wy[2]) ;
    fdt = _mm256_fmadd_pd(fd2,fwy,fdt) ;
    fet = _mm256_fmadd_pd(fe2,fwy,fet) ;
    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double
    fr2 = _mm_loadu_ps(f2+ni2) ;            // row 3 : f[i+2:i+5 , j+2 ,k]
    fe2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    fwy = _mm256_set1_pd(wy[3]) ;
    fdt = _mm256_fmadd_pd(fd3,fwy,fdt) ;
    fet = _mm256_fmadd_pd(fe3,fwy,fet) ;
    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double
    fr3 = _mm_loadu_ps(f2+ni3) ;            // row 4 : f[i+2:i+5 , j+3 ,k]
    fe3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    fwy = _mm256_set1_pd(wy[4]) ;
    fdt = _mm256_fmadd_pd(fd4,fwy,fdt) ;
    fet = _mm256_fmadd_pd(fe4,fwy,fet) ;
    fr4 = _mm_loadu_ps(f+ni4) ;             // row 5 : f[i:i+3 , j+4 ,k]
    fd4 = _mm256_cvtps_pd(fr4) ;            // promote row 5 to double
    fr4 = _mm_loadu_ps(f2+ni4) ;            // row 4 : f[i+2:i+5 , j+4 ,k]
    fe4 = _mm256_cvtps_pd(fr4) ;            // promote row 4 to double

    fwy = _mm256_set1_pd(wy[5]) ;
    fdt = _mm256_fmadd_pd(fd5,fwy,fdt) ;
    fet = _mm256_fmadd_pd(fe5,fwy,fet) ;
    fr5 = _mm_loadu_ps(f+ni5) ;             // row 6 : f[i:i+3 , j+5 ,k]
    fd5 = _mm256_cvtps_pd(fr5) ;            // promote row 5 to double
    fr5 = _mm_loadu_ps(f2+ni5) ;            // row 4 : f[i+2:i+5 , j+5 ,k]
    fe5 = _mm256_cvtps_pd(fr5) ;            // promote row 4 to double

    // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
    fdw = _mm256_loadu_pd(&wx[0]) ;
    few = _mm256_loadu_pd(&wx[2]) ;
    fdt = _mm256_mul_pd(fdt,fdw) ;          // fff[0:3]
    fet = _mm256_mul_pd(fet,few) ;          // fff[2:5]
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fff[2]                      fff[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fff[0]                      fff[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fff[0]+fff[2]               fff[1]+fff[3]
    ft1 = _mm256_extractf128_pd(fet,1) ;    // fff[4]                      fff[5]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fff[0]+fff[2]+fff[4]        fff[1]+fff[3]+fff[5]
    ft1 = _mm_permute_pd(ft0,0x3) ;         // fff[1]+fff[3]+fff[5]        fff[0]+fff[2]+fff[4]
    ft0 = _mm_add_sd(ft0,ft1) ;             // fff[0]+fff[2]+fff[1]+fff[3]+fff[4]+fff[5]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert to float
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
  // interpolation along j , level nk
  fwy = _mm256_set1_pd(wy[0]) ;
  fdt = _mm256_mul_pd(fd0,fwy) ;            // sum of row[j] * coefficient[j]
  fet = _mm256_mul_pd(fe0,fwy) ;            // sum of row[j] * coefficient[j]
  fwy = _mm256_set1_pd(wy[1]) ;
  fdt = _mm256_fmadd_pd(fd1,fwy,fdt) ;
  fet = _mm256_fmadd_pd(fe1,fwy,fet) ;
  fwy = _mm256_set1_pd(wy[2]) ;
  fdt = _mm256_fmadd_pd(fd2,fwy,fdt) ;
  fet = _mm256_fmadd_pd(fe2,fwy,fet) ;
  fwy = _mm256_set1_pd(wy[3]) ;
  fdt = _mm256_fmadd_pd(fd3,fwy,fdt) ;
  fet = _mm256_fmadd_pd(fe3,fwy,fet) ;
  fwy = _mm256_set1_pd(wy[4]) ;
  fdt = _mm256_fmadd_pd(fd4,fwy,fdt) ;
  fet = _mm256_fmadd_pd(fe4,fwy,fet) ;
  fwy = _mm256_set1_pd(wy[5]) ;
  fdt = _mm256_fmadd_pd(fd5,fwy,fdt) ;
  fet = _mm256_fmadd_pd(fe5,fwy,fet) ;

  // interpolation along i: multiply by coefficients along x , then sum elements (using vector folding)
  fdw = _mm256_loadu_pd(&wx[0]) ;
  few = _mm256_loadu_pd(&wx[2]) ;
  fdt = _mm256_mul_pd(fdt,fdw) ;          // fff[0:3]
  fet = _mm256_mul_pd(fet,few) ;          // fff[2:5]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fff[2]                      fff[3]
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fff[0]                      fff[1]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fff[0]+fff[2]               fff[1]+fff[3]
  ft1 = _mm256_extractf128_pd(fet,1) ;    // fff[4]                      fff[5]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fff[0]+fff[2]+fff[4]        fff[1]+fff[3]+fff[5]
  ft1 = _mm_permute_pd(ft0,0x3) ;         // fff[1]+fff[3]+fff[5]        fff[0]+fff[2]+fff[4]
  ft0 = _mm_add_sd(ft0,ft1) ;             // fff[0]+fff[2]+fff[1]+fff[3]+fff[4]+fff[5]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert to float
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
    intrp_biquint_yx(&f[0][0][0], &r[0][i], NI, NI*NJ, NK, NP, x[i], y[i]) ;
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
  real(C_DOUBLE), dimension(NP) :: x, y, xmin, xmax, xmin2, xmax2
  integer :: i, j, k
  integer :: i0, j0
  integer*8, external :: rdtsc
  integer*8 :: t1, t2, tmg1(NR), tmg2(nr)
  integer :: nidim, ninjdim
  interface
    subroutine intrp_biquint_yx(f, r, ni, ninj, nk, np, x, y) bind(C,name='intrp_biquint_yx')         !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np                                       !InTf!
    end subroutine intrp_biquint_yx                                                                 !InTf!
    subroutine intrp_biquint_yx_mono(f, r, ni, ninj, nk, np, x, y) bind(C,name='intrp_biquint_yx_mono') !InTf!
      import :: C_INT, C_FLOAT, C_DOUBLE                                                          !InTf!
      real(C_FLOAT), dimension(*), intent(IN) :: f                                                !InTf!
      real(C_FLOAT), dimension(*), intent(OUT) :: r                                               !InTf!
      real(C_DOUBLE), intent(IN), value :: x, y                                                   !InTf!
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np                                       !InTf!
    end subroutine intrp_biquint_yx_mono                                                            !InTf!
  end interface
#define FXY(A,B,C) (1.1*(A)**5 + 1.2*(A)**4 + 1.3*(A)**3 + 1.4*(A)**2 + (A)*1.5 + 2.1*(B)**5 + 2.2*(B)**4 + 2.3*(B)**3 + 2.4*(B)**2 + (B)*2.5 + (C))

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
    x(i) = i - .999
    y(i) = i - .999
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
      call intrp_biquint_yx( f(1,1,1), r(i,1), nidim, ninjdim, NK, NP, x(i), y(i) )
    enddo
    t2 = rdtsc()
    tmg1(j) = t2 - t1
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
  print *,'F matrix (1)'
  do j = 5, 1-hy, -1
    print 102,f(1-hx:5,j,1)
  enddo
  print *,'F matrix (nk)'
  do j = 5, 1-hy, -1
    print 102,f(1-hx:5,j,NK)
  enddo
  print *,'MONO: limit (min, max)'
  do i = 1 , NP
    i0 = x(i)
    j0 = y(i)
    xmin(i) = minval(f(i0:i0+1,j0:j0+1,1))
    xmax(i) = maxval(f(i0:i0+1,j0:j0+1,1))
    xmin2(i) = minval(f(i0:i0+1,j0:j0+1,NK))
    xmax2(i) = maxval(f(i0:i0+1,j0:j0+1,NK))
  enddo
  print 102,xmin(:) , xmin2(:)
  print 102,xmax(:) , xmax2(:)
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
    call intrp_biquint_yx( f(1,1,1), r(i,1), nidim, ninjdim, NK, NP, x(i), y(i) )
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
