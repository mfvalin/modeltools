#if ! defined(FORTRAN_CODE)
static float cp133 =  0.166666666666666667E0;
static float cm133 = -0.166666666666666667E0;
static float cp5 =  .5;
static float cm5 = -.5;
static float one = 1.0;
static float two = 2.0;

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

#include <immintrin.h>

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
  ix = x ; ix = ix - 1;   // xx and yy are in "ORIGIN 1"
  iy = y ; iy = iy - 1;
  x  = x - 1 - ix;
  y  = y - 1 - iy;
  f = f + ix + iy * ni;
// printf("DEBUG: f[0] = %f, ix = %d, iy = %d, x = %f, y = %f\n",f[0],ix,iy,x,y);
  wx[0] = cm133*x*(x-one)*(x-two);
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
#endif
#if defined(__AVX2__) && defined(__x86_64__)
  // fetch 4 rows, level 0
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
    // interpolation along J level k, prefetch 4 rows level k+1
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

    // interpolation along i, multiply by coefficients along x , then sum elements
    fdt = _mm256_mul_pd(fdt,fwx) ;
    ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
    ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
    ft1 = _mm_permute_pd(ft0,0x55) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[1]+fdt[3]+fdt[1]+fdt[3] fdt[1]+fdt[3]+fdt[1]+fdt[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[1]+fdt[3]+fdt[1]+fdt[3] to float
    _mm_store_ss(r,frt) ;                   // store float
    r += np;
  }
  // interpolation along J , last value along k
  fdt = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
  fdt = _mm256_fmadd_pd(fd1,fwy1,fdt) ;
  fdt = _mm256_fmadd_pd(fd2,fwy2,fdt) ;
  fdt = _mm256_fmadd_pd(fd3,fwy3,fdt) ;
  // interpolation along i, multiply by coefficients along x , then sum elements
  fdt = _mm256_mul_pd(fdt,fwx) ;
  ft0 = _mm256_extractf128_pd(fdt,0) ;    // fdt[0]                      fdt[1]
  ft1 = _mm256_extractf128_pd(fdt,1) ;    // fdt[2]                      fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[0]+fdt[2]               fdt[1]+fdt[3]
  ft1 = _mm_permute_pd(ft0,0x55) ;        // fdt[1]+fdt[3]               fdt[1]+fdt[3]
  ft0 = _mm_add_pd(ft0,ft1) ;             // fdt[1]+fdt[3]+fdt[1]+fdt[3] fdt[1]+fdt[3]+fdt[1]+fdt[3]
  frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fdt[1]+fdt[3]+fdt[1]+fdt[3] to float
  _mm_store_ss(r,frt) ;                   // store float
#else
  for(k=1 ; k<nk ; k++){
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
void int_yinyang_cub_yx_mono(float *f, float *r, int ni, int ninj, int nk, int np, double x, double y){
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d fd0, fd1, fd2, fd3, fwx, fwy0, fwy1, fwy2, fwy3 ;
  __m128  fr0, fr1, fr2, fr3, frt ;
  __m128d ft0, ft1;
  __m128  fmi, fma;
#else
  double fd0[4], fd1[4], fd2[4], fd3[0] ;
#endif
  double  wx[4], wy[4] ;
  int ni2 = ni + ni;    // + 2 rows
  int ni3 = ni2 + ni;   // + 3 rows
  int i, k;

  wx[0] = cm133*x*(x-one)*(x-two);
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm133*y*(y-one)*(y-two);
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

#if defined(__AVX2__) && defined(__x86_64__)
  fwx = _mm256_loadu_pd(wx) ;              // vector of coefficients along i
  fwy0 = _mm256_set1_pd(wy[0]) ;           // scalar * vector not available, promote scalar to vector
  fwy1 = _mm256_set1_pd(wy[1]) ;
  fwy2 = _mm256_set1_pd(wy[2]) ;
  fwy3 = _mm256_set1_pd(wy[3]) ;
#endif
  for(k=0 ; k<nk ; k++){
#if defined(__AVX2__) && defined(__x86_64__)
    fr0 = _mm_loadu_ps(f) ;                 // row 1 : f[i:i+3 , j   ,k]
    fmi = fr0 ;                             // min = first row
    fma = fr0 ;                             // max = first row
    fd0 = _mm256_cvtps_pd(fr0) ;            // promote row 1 to double

    fr1 = _mm_loadu_ps(f+ni) ;              // row 2 : f[i:i+3 , j+1 ,k]
    fmi = _mm_min_ps(fmi,fr1) ;
    fma = _mm_max_ps(fma,fr1) ;
    fd1 = _mm256_cvtps_pd(fr1) ;            // promote row 2 to double

    fr2 = _mm_loadu_ps(f+ni2) ;             // row 3 : f[i:i+3 , j+2 ,k]
    fmi = _mm_min_ps(fmi,fr2) ;
    fma = _mm_max_ps(fma,fr2) ;
    fd2 = _mm256_cvtps_pd(fr2) ;            // promote row 3 to double

    fr3 = _mm_loadu_ps(f+ni3) ;             // row 4 : f[i:i+3 , j+3 ,k]
    fmi = _mm_min_ps(fmi,fr3) ;             // min of 4 rows
    fma = _mm_max_ps(fma,fr3) ;             // max of 4 rows
    fd3 = _mm256_cvtps_pd(fr3) ;            // promote row 4 to double

    // interpolation along J
    fd0 = _mm256_mul_pd(fd0,fwy0) ;            // sum of row[j] * coefficient[j]
    fd0 = _mm256_fmadd_pd(fd1,fwy1,fd0) ;
    fd0 = _mm256_fmadd_pd(fd2,fwy2,fd0) ;
    fd0 = _mm256_fmadd_pd(fd3,fwy3,fd0) ;

    // get minimum of 4 vector elements
    frt = _mm_permute_ps(fmi,0xEE) ;        // fmi[2]              fmi[3]             fmi[2]  fmi[3] 
    fmi = _mm_min_ps(fmi,frt) ;             // min(fmi[0],fmi[2])  min(fmi[1],fmi[3]) fmi[2]  fmi[3]
    frt = _mm_permute_ps(fmi,0x55) ;        // fmi[1]              fmi[1]             fmi[1]  fmi[1]
    fmi = _mm_min_ps(fmi,frt) ;             // min(fmi[0],fmi[1])

    // get maximum of 4 vector elements
    frt = _mm_permute_ps(fma,0xEE) ;        // fma[2]              fma[3]             fma[2]  fma[3] 
    fma = _mm_max_ps(fma,frt) ;             // max(fma[0],fma[2])  max(fma[1],fma[3]) fma[2]  fma[3]
    frt = _mm_permute_ps(fma,0x55) ;        // fma[1]              fma[1]             fma[1]  fma[1]
    fma = _mm_max_ps(fma,frt) ;             // max(fma[0],fma[1])

    // interpolation along i, multiply by coefficients along x , then sum elements
    fd0 = _mm256_mul_pd(fd0,fwx) ;
    ft0 = _mm256_extractf128_pd(fd0,0) ;    // fd0[0]                      fd0[1]
    ft1 = _mm256_extractf128_pd(fd0,1) ;    // fd0[2]                      fd0[3]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fd0[0]+fd0[2]               fd0[1]+fd0[3]
    ft1 = _mm_permute_pd(ft0,0x55) ;        // fd0[1]+fd0[3]               fd0[1]+fd0[3]
    ft0 = _mm_add_pd(ft0,ft1) ;             // fd0[1]+fd0[3]+fd0[1]+fd0[3] fd0[1]+fd0[3]+fd0[1]+fd0[3]
    frt = _mm_cvtsd_ss(frt,ft0) ;           // convert fd0[1]+fd0[3]+fd0[1]+fd0[3] to float
    frt = _mm_min_ss(frt,fmi) ;
    frt = _mm_max_ss(frt,fma) ;
    _mm_store_ss(r,frt) ;                   // store float
#else
    for(i=0 ; i<4 ; i++){                   // easily vectorizable form
      fd0[i] = ( f[i]*wy[0] + f[i+ni]*wy[1] + f[i+ni2]*wy[2] + f[i+ni3]*wy[3] ) * wx[i];
    }
    r[0] = fd0[0] + fd0[1] + fd0[2] + fd0[3];
#endif
    f+= ninj;
    r += np;
  }
}
#endif

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define NI 65
#define NJ 27
#define NK 80
#define NP 2

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
  integer, parameter :: NK=80
  integer, parameter :: NP=4
  integer, parameter :: HX=2
  integer, parameter :: HY=2
  real(C_FLOAT), dimension(1-HX:NI+HX , 1-HY:NJ+HY , NK) :: f
  real(C_FLOAT), dimension(NP,NK) :: r
  real(C_DOUBLE), dimension(NP) :: x, y
  integer :: i, j, k
  integer*8, external :: rdtsc
  real *8 :: t1, t2
  integer :: nidim, ninjdim
  interface
    subroutine int_yinyang_cub_yx(f, r, ni, ninj, nk, np, x, y) bind(C,name='int_yinyang_cub_yx')
      import :: C_INT, C_FLOAT, C_DOUBLE
      real(C_FLOAT), dimension(*), intent(IN) :: f
      real(C_FLOAT), dimension(*), intent(OUT) :: r
      real(C_DOUBLE), intent(IN), value :: x, y
      integer(C_INT), intent(IN), value :: ni, ninj, nk, np
    end subroutine int_yinyang_cub_yx
  end interface

  r = 9999.99
  do k = 1 , NK
    do j = 1-HY , NJ+HY
      do i = 1-HX , NI+HX
        f(i,j,k) = i + j  + k
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
  t1 = rdtsc()
  do i = 1 , NP
    call int_yinyang_cub_yx( f(1,1,1), r(i,1), nidim, ninjdim, NK, NP, x(i), y(i) )
  enddo
  t2 = rdtsc()
  print *,'time=',t2-t1,' cycles for',NP*NK*35,' values'
!  print 102,f(nint(x(1)),nint(y(1)),1),f(nint(x(2)),nint(y(2)),1),f(nint(x(1)),nint(y(1)),NK),f(nint(x(2)),nint(y(2)),NK)
!  print 102,x(1)+y(1)+1,x(2)+y(2)+1,x(1)+y(1)+NK,x(2)+y(2)+NK
  print *,' expected'
  print 102,x(:)+y(:)+1, x(:)+y(:)+NK
  print *,' got'
  print 102,r(:,1),r(:,NK)
  print *,' delta'
  print 102,(r(:,1)-(x(:)+y(:)+1)),(r(:,NK)-(x(:)+y(:)+NK))
  print 103,(r(:,1)-(x(:)+y(:)+1))/(x(:)+y(:)+1) , (r(:,NK)-(x(:)+y(:)+NK))/(x(:)+y(:)+NK)
102 format(16F15.6)
103 format(16E15.2)
end
#endif
