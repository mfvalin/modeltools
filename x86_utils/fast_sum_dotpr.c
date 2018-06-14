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
/* Fortran interface definitions
   interface
    ! void fast_dot_product_f3(const float *fa, const float *fb, const float *fc, double *r, int n)
    subroutine fast_dot_product_f3(fa, fb, fc, r, n) bind(C,name='fast_dot_product_f3')
      import :: C_FLOAT, C_DOUBLE, C_INT
      integer(C_INT), intent(IN), value :: n
      real(C_FLOAT), dimension(n), intent(IN) :: fa, fb, fc
      real(C_DOUBLE), dimension(3), intent(OUT) :: r
    end subroutine fast_dot_product_f3
    ! double fast_dot_product_f(const float *fa, const float *fa, int n)
    function fast_dot_product_f(fa, fb, n) result(r) bind(C,name='fast_dot_product_f')
      import :: C_FLOAT, C_DOUBLE, C_INT
      integer(C_INT), intent(IN), value :: n
      real(C_FLOAT), dimension(n), intent(IN) :: fa, fb
      real(C_DOUBLE) :: r
    end function fast_dot_product_f 
    ! double fast_normdot_f(const float *fa, int n)
    function fast_normdot_f(fa, n) result(r) bind(C,name='fast_normdot_f')
      import :: C_FLOAT, C_DOUBLE, C_INT
      integer(C_INT), intent(IN), value :: n
      real(C_FLOAT), dimension(n), intent(IN) :: fa
      real(C_DOUBLE) :: r
    end function fast_normdot_f 
   end interface
 */
#include <immintrin.h>
#include <stdint.h>
typedef int xsum_length;
typedef double xsum_flt;

// triple dot product fa*fa, fa*fb, fb*fc  (real arrays)
// result in r[0], r[1], r[2]              (double)
void fast_dot_product_f3(const float *fa, const float *fb, const float *fc, double *r, int n){
  __m256d da0, da1, db0, db1, dc0, dc1;
  __m256d tda, tdb, tdc;
  __m256  tt0;
  __m128  tf0;
  float t1[8], t2[8], t3[8];
  double ta[4], tb[4], tc[4];
  int i, j;
  tt0 = _mm256_xor_ps(tt0,tt0);
  _mm256_storeu_ps((float *) &t1[ 0], tt0);
  _mm256_storeu_ps((float *) &t2[ 0], tt0);
  _mm256_storeu_ps((float *) &t3[ 0], tt0);
  for(i=0 ; i<(n&0x7); i++) { t1[i] = fa[i] ; t2[i] = fb[i];  t3[i] = fc[i]; }

  tf0 = _mm_loadu_ps((float *) &t1[0]);
  tda = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t2[0]);
  tdb = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t3[0]);
  tdc = _mm256_cvtps_pd(tf0);
  da0 = _mm256_mul_pd(tda, tda);
  db0 = _mm256_mul_pd(tda, tdb);
  dc0 = _mm256_mul_pd(tdb, tdc);

  tf0 = _mm_loadu_ps((float *) &t1[4]);
  tda = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t2[4]);
  tdb = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t3[4]);
  tdc = _mm256_cvtps_pd(tf0);
  da1 = _mm256_mul_pd(tda, tda);
  db1 = _mm256_mul_pd(tda, tdb);
  dc1 = _mm256_mul_pd(tdb, tdc);

  for(j=i ; j<n ; j+=8){      // 8 values processed per iteration
    tf0 = _mm_loadu_ps((float *) &fa[j + 0]);
    tda = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &fb[j + 0]);
    tdb = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &fc[j + 0]);
    tdc = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    da0 = _mm256_fmadd_pd (tda, tda, da0);
    dc0 = _mm256_fmadd_pd (tdb, tdc, dc0);
    db0 = _mm256_fmadd_pd (tda, tdb, db0);
#else
    tdc = _mm256_mul_pd (tdb, tdc);
    dc0 = _mm256_add_pd (tdc, dc0);
    tdb = _mm256_mul_pd (tda, tdb);
    db0 = _mm256_add_pd (tdb, db0);
    tda = _mm256_mul_pd (tda, tda);
    da0 = _mm256_add_pd (tda, da0);
#endif

    tf0 = _mm_loadu_ps((float *) &fa[j + 4]);
    tda = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &fb[j + 4]);
    tdb = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &fc[j + 4]);
    tdc = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    da1 = _mm256_fmadd_pd (tda, tda, da1);
    db1 = _mm256_fmadd_pd (tda, tdb, db1);
    dc1 = _mm256_fmadd_pd (tdb, tdc, dc1);
#else
    tdc = _mm256_mul_pd (tdb, tdc);
    dc1 = _mm256_add_pd (tdc, dc1);
    tdb = _mm256_mul_pd (tda, tdb);
    db1 = _mm256_add_pd (tdb, db1);
    tda = _mm256_mul_pd (tda, tda);
    da1 = _mm256_add_pd (tda, da1);
#endif
  }
  da0 =  _mm256_add_pd (da0, da1);   // add subtotal pairs into da0, db0, dc0
  db0 =  _mm256_add_pd (db0, db1);
  dc0 =  _mm256_add_pd (dc0, dc1);
  _mm256_storeu_pd ((double *) &ta[0], da0) ;
  _mm256_storeu_pd ((double *) &tb[0], db0) ;
  _mm256_storeu_pd ((double *) &tc[0], dc0) ;
  r[0] = ta[0] + ta[1] + ta[2] + ta[3];
  r[1] = tb[0] + tb[1] + tb[2] + tb[3];
  r[2] = tc[0] + tc[1] + tc[2] + tc[3];
}

// dot product f1*f2 (real arrays)
// function result   (double)
double fast_dot_product_f(const float *f1, const float *f2, int n){
  __m256d d0, d1, d2, d3;
  __m256d td0, td1;
  __m256  tt0;
  __m128  tf0;
  float t1[16], t2[16];
  double td[4];
  int i, j;
  tt0 = _mm256_xor_ps(tt0,tt0);
  _mm256_storeu_ps((float *) &t1[ 0], tt0);
  _mm256_storeu_ps((float *) &t1[ 8], tt0);
  _mm256_storeu_ps((float *) &t2[ 0], tt0);
  _mm256_storeu_ps((float *) &t2[ 8], tt0);
  for(i=0 ; i<(n&0xF); i++) { t1[i] = f1[i] ; t2[i] = f2[i]; }

  tf0 = _mm_loadu_ps((float *) &t1[ 0]);
  d0 = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t2[ 0]);
  td0 = _mm256_cvtps_pd(tf0);
  d0 = _mm256_mul_pd(d0, td0);

  tf0 = _mm_loadu_ps((float *) &t1[ 4]);
  d1 = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t2[ 4]);
  td0 = _mm256_cvtps_pd(tf0);
  d1 = _mm256_mul_pd(d1, td0);

  tf0 = _mm_loadu_ps((float *) &t1[ 8]);
  d2 = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t2[ 8]);
  td0 = _mm256_cvtps_pd(tf0);
  d2 = _mm256_mul_pd(d2, td0);

  tf0 = _mm_loadu_ps((float *) &t1[12]);
  d3 = _mm256_cvtps_pd(tf0);
  tf0 = _mm_loadu_ps((float *) &t2[12]);
  td0 = _mm256_cvtps_pd(tf0);
  d3 = _mm256_mul_pd(d3, td0);

  for(j=i ; j<n ; j+=16){      // 16 values processed per iteration
    tf0 = _mm_loadu_ps((float *) &f1[j + 0]);
    td0 = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &f2[j + 0]);
    td1 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d0 = _mm256_fmadd_pd (td0, td1, d0);
#else
    td0 = _mm256_mul_pd(td1, td0);
    d0 = _mm256_add_pd (td0, d0);
#endif

    tf0 = _mm_loadu_ps((float *) &f1[j + 4]);
    td0 = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &f2[j + 4]);
    td1 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d1 = _mm256_fmadd_pd (td0, td1, d1);
#else
    td0 = _mm256_mul_pd(td1, td0);
    d1 = _mm256_add_pd (td0, d1);
#endif

    tf0 = _mm_loadu_ps((float *) &f1[j + 8]);
    td0 = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &f2[j + 8]);
    td1 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d2 = _mm256_fmadd_pd (td0, td1, d2);
#else
    td0 = _mm256_mul_pd(td1, td0);
    d2 = _mm256_add_pd (td0, d2);
#endif

    tf0 = _mm_loadu_ps((float *) &f1[j +12]);
    td0 = _mm256_cvtps_pd(tf0);
    tf0 = _mm_loadu_ps((float *) &f2[j +12]);
    td1 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d3 = _mm256_fmadd_pd (td0, td1, d3);
#else
    td0 = _mm256_mul_pd(td1, td0);
    d3 = _mm256_add_pd (td0, d3);
#endif
  }
  d0 =  _mm256_add_pd (d0, d1);   // add 4 subtotals into d0
  d2 =  _mm256_add_pd (d2, d3);
  d0 =  _mm256_add_pd (d0, d2);
  _mm256_storeu_pd ((double *) &td[0], d0) ;
  return td[0] + td[1] + td[2] + td[3];
}

// dot product f1*f1 (real array)
// function result   (double)
double fast_normdot_f(const float *f1, int n){
  __m256d d0, d1, d2, d3;
  __m256d td0, td1;
  __m256  tt0;
  __m128  tf0;
  float t1[16];
  double td[4];
  int i, j;
  tt0 = _mm256_xor_ps(tt0,tt0);
  _mm256_storeu_ps((float *) &t1[ 0], tt0);
  _mm256_storeu_ps((float *) &t1[ 8], tt0);
  for(i=0 ; i<(n&0xF); i++) { t1[i] = f1[i] ; }

  tf0 = _mm_loadu_ps((float *) &t1[ 0]);
  d0 = _mm256_cvtps_pd(tf0);
  d0 = _mm256_mul_pd(d0, d0);

  tf0 = _mm_loadu_ps((float *) &t1[ 4]);
  d1 = _mm256_cvtps_pd(tf0);
  d1 = _mm256_mul_pd(d1, d1);

  tf0 = _mm_loadu_ps((float *) &t1[ 8]);
  d2 = _mm256_cvtps_pd(tf0);
  d2 = _mm256_mul_pd(d2, d2);

  tf0 = _mm_loadu_ps((float *) &t1[12]);
  d3 = _mm256_cvtps_pd(tf0);
  d3 = _mm256_mul_pd(d3, d3);

  for(j=i ; j<n ; j+=16){      // 16 values processed per iteration
    tf0 = _mm_loadu_ps((float *) &f1[j + 0]);
    td0 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d0 = _mm256_fmadd_pd (td0, td0, d0);
#else
    td0 = _mm256_mul_pd(td0, td0);
    d0 = _mm256_add_pd (td0, d0);
#endif

    tf0 = _mm_loadu_ps((float *) &f1[j + 4]);
    td0 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d1 = _mm256_fmadd_pd (td0, td0, d1);
#else
    td0 = _mm256_mul_pd(td0, td0);
    d1 = _mm256_add_pd (td0, d1);
#endif

    tf0 = _mm_loadu_ps((float *) &f1[j + 8]);
    td0 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d2 = _mm256_fmadd_pd (td0, td0, d2);
#else
    td0 = _mm256_mul_pd(td0, td0);
    d2 = _mm256_add_pd (td0, d2);
#endif

    tf0 = _mm_loadu_ps((float *) &f1[j +12]);
    td0 = _mm256_cvtps_pd(tf0);
#if defined(WITH_FMA)
    d3 = _mm256_fmadd_pd (td0, td0, d3);
#else
    td0 = _mm256_mul_pd(td0, td0);
    d3 = _mm256_add_pd (td0, d3);
#endif
  }
  d0 =  _mm256_add_pd (d0, d1);   // add 4 subtotals into d0
  d2 =  _mm256_add_pd (d2, d3);
  d0 =  _mm256_add_pd (d0, d2);
  _mm256_storeu_pd ((double *) &td[0], d0) ;
  return td[0] + td[1] + td[2] + td[3];
}

// sum of contents of double array vec, length n, using tweaked version of Kahan's algorithm
// function result (double)
xsum_flt xsum_sum_kahan (const xsum_flt *vec, xsum_length n){
  xsum_flt s, t, c, y;
  xsum_length i, j;
#if defined(USE_ORIGINAL_CODE)
  int always = 1;
#else
  int always = 0;
#endif
  if(always || n < 16){
    s = 0.0;
    c = 0.0;
    { for (j = 1; j < n; j += 2) 
      { y = vec[j-1] - c;
	t = s;
	s += y;
	c = (s - t) - y;
	y = vec[j] - c;
	t = s;
	s += y;
	c = (s - t) - y;
      }
      for (j = j-1; j < n; j++) 
      { y = vec[j] - c;
	t = s;
	s += y;
	c = (s - t) - y;
      }
    }
    return (xsum_flt) s;
  }
#if ! defined(USE_ORIGINAL_CODE)
  xsum_flt vs[8], vt[4], vc[4], vy[4];
  __m256d vs0, vt0, vc0, vy0;
  __m256d vs1, vt1, vc1, vy1;

  vc0 = _mm256_xor_pd(vc0,vc0);                // error = 0
  vc1 = _mm256_xor_pd(vc1,vc1);
  _mm256_storeu_pd ((double *) &vs[0], vc0) ;  // initialize vs to 0
  _mm256_storeu_pd ((double *) &vs[4], vc0) ;
  for(i=0 ; i<(n&0x7); i++) vs[i] = vec[i];    // copy n mod 7 elements

//     for(i=0 ; i<4; i++){
//       vy[i] = vec[j+i] - vc[i];
//       vt[i] = vs[i];
//       vs[i] += vy[i];
//       vc[i] = (vs[i] - vt[i]) - vy[i];
//     }

  vs0 = _mm256_loadu_pd ((double *) &vs[0]) ;   // first 4 "sums"
  vs1 = _mm256_loadu_pd ((double *) &vs[4]) ;   // next 4 "sums"
  for(j=i ; j<n ; j+=8){
    vt0 = _mm256_loadu_pd((double *) &vec[j+0]);  // vec[j+i]
    vt1 = _mm256_loadu_pd((double *) &vec[j+4]);
    vy0 = _mm256_sub_pd(vt0,vc0);                // vy[i] = vec[j+i] - vc[i]
    vy1 = _mm256_sub_pd(vt1,vc1);
    vt0 = vs0;                                   // vt[i] = vs[i]
    vt1 = vs1;
    vs0 = _mm256_add_pd(vs0,vy0);                // vs[i] += vy[i]
    vs1 = _mm256_add_pd(vs1,vy1);
    vc0 = _mm256_sub_pd(vs0,vt0);                // vc[i] = (vs[i] - vt[i])
    vc1 = _mm256_sub_pd(vs1,vt1);
    vc0 = _mm256_sub_pd(vc0,vy0);                // vc[i] = (vs[i] - vt[i]) - vy[i]
    vc1 = _mm256_sub_pd(vc1,vy1);
  }
//   s = vs[0]; c = 0;
//   for (j = 1; j < 4; j++) {
//     y = vs[j] - c;
//     t = s;
//     s += y;
//     c = (s - t) - y;
//   }
//   return s + c;
  vs0 = _mm256_add_pd(vs0,vc0);  // add error terms [0,1,2,3]
  vs1 = _mm256_add_pd(vs1,vc1);  // add error terms [4,5,6,7]
  vt0 = vs0;                     // vt[i] = vs[i]
  vy0 = vs1;
  vs0 = _mm256_add_pd(vs0,vs1);  // vs[i] += vy[i]   ( vy = vs1 )
  vc0 = _mm256_sub_pd(vs0,vt0);  // vc[i] = (vs[i] - vt[i])
  vc0 = _mm256_sub_pd(vc0,vy0);  // vc[i] = (vs[i] - vt[i]) - vy[i]
  vs0 = _mm256_add_pd(vs0,vc0);  // add error terms

  vs1 = _mm256_shuffle_pd(vs0,vs0,0x5) ;  // vs0[1,0,3,2]

  // vs0 = [0,1,2,3]
  // vs1 = [1,0,3,2]
  vt0 = vs0;                     // vt[i] = vs[i]
  vy0 = vs1;                     // vs1 - 0 because  error = 0
  vs0 = _mm256_add_pd(vs0,vs1);  // vs[i] += vy[i]   ( vy = vs1 )
  vc0 = _mm256_sub_pd(vs0,vt0);  // vc[i] = (vs[i] - vt[i])
  vc0 = _mm256_sub_pd(vc0,vy0);  // vc[i] = (vs[i] - vt[i]) - vy[i]
  vs0 = _mm256_add_pd(vs0,vc0);  // add error terms  [0+1,1+0,2+3,3+2]

  _mm256_storeu_pd ((double *) &vs[0], vs0) ;   // store for scalar wrapup
  return vs[0] + vs[2] ;                        // no fancy footwork for last addition
#endif
}

#if defined(SELF_TEST)

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#define NP 1000
#define NPF 1000000
#define REP 10000
#include <stdio.h>
int main(){
  int i, j;
  double z[NP], r[3];
  float f1[NPF], f2[NPF], f3[NPF];
  double sum, sum2, sum3, sum4;
  uint64_t t0, t1;
  sum2 = 0;
  sum3 = 0;
  sum4 = 0;
  for(i=0 ; i<NP; i++) z[i] = (i+1) * 111.111;
  for(i=0 ; i<NPF; i++) { f1[i] = 1.1  ; f2[i] = 1.5 ; f3[i] = 1.0 ;}
  for(i=0 ; i<NP; i++) sum2 += z[i];
  for(i=0 ; i<NPF; i++) sum3 += f1[i]*f2[i];
  for(i=0 ; i<NPF; i++) sum4 += f1[i]*f1[i];
  sum = xsum_sum_kahan (z , NP);
  printf("sum: expected = %g, got = %g\n\n",sum2,sum);

  t0 = rdtsc();
  for (j=0; j<REP; j++) sum = fast_dot_product_f(f1, f2, NPF);
  t1 = rdtsc();
  printf("dot time : %g cycles/value\n",1.0*(t1-t0)/NPF/REP);
  printf("dot: expected = %g, got = %g\n\n",sum3,sum);

  t0 = rdtsc();
  for (j=0; j<REP; j++) sum = fast_dot_product_f(f1, f1, NPF);
  t1 = rdtsc();
  printf("self dot time : %g cycles/value\n",1.0*(t1-t0)/NPF/REP);
  printf("dot: expected = %g, got = %g\n\n",sum4,sum);

  t0 = rdtsc();
  for (j=0; j<REP; j++) sum = fast_normdot_f(f1, NPF);
  t1 = rdtsc();
  printf("normdot time : %g cycles/value\n",1.0*(t1-t0)/NPF/REP);
  printf("dot: expected = %g, got = %g\n\n",sum4,sum);

  t0 = rdtsc();
  for (j=0; j<REP; j++) fast_dot_product_f3(f1, f2, f3, r, NPF);
  t1 = rdtsc();
  sum2 = 0; sum3 = 0; sum4 = 0;
  for(i=0 ; i<NPF; i++) sum2 += f1[i]*f1[i];
  for(i=0 ; i<NPF; i++) sum3 += f1[i]*f2[i];
  for(i=0 ; i<NPF; i++) sum4 += f2[i]*f3[i];
  printf("dot3 time : %g cycles/value\n",1.0*(t1-t0)/NPF/REP/3);
  printf("dot3:      got = %g, %g, %g\n",r[0],r[1],r[2]);
  printf("dot3: expected = %g, %g, %g\n\n",sum2,sum3,sum4);
  return 0;
}
#endif
