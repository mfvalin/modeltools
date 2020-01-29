/* useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2015  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdint.h>
#include <immintrin.h>

static int64_t __attribute__((aligned(64))) mask64[] = { -1, -1, -1, -1, 0, 0, 0, 0 };
static int32_t *mask32 = (int32_t *) mask64;

#if defined(OLD_CODE)

double dot_product_s3(int ni, int nj, int li, float *a, float *b, float *c){
  double result = 0;
  register double aa, bb, cc;
  int i;
  for( ; nj > 0 ; nj--){  // loop over rows
    for(i=0 ; i<ni ; i++){
      aa = a[i] ;
      bb = b[i] ;
      cc = c[i] ;
      result +=  aa * bb * cc;
    }
    a += li ; b += li ; c += li ; // jump between rows
  }
  return result;
}

double dot_product_v3(int ni, int nj, int li, float *a, float *b, float *c){
  double result = 0;
  int i;

#if defined(__AVX2__) & defined(__FMA__)
  int res = ni & 3;
  __m256d d0, d1, d2, d3, td0, td1, td2, td3;
  __m128d f0, f1;
  __m128  tf0;

  td0 = _mm256_xor_pd(td0,td0);
  td1 = _mm256_xor_pd(td1,td1);

  for( ; nj > 0 ; nj--){  // loop over rows
    if(res > 0){
      tf0 = _mm_load_ps((float *) (mask32 + 4 + 4 - res) );   // first ni modulo 4 elements if modulo nonzero
      d0 = _mm256_cvtps_pd(_mm_and_ps(tf0, _mm_load_ps(a)));  // mask unused elements
      d1 = _mm256_cvtps_pd(_mm_and_ps(tf0, _mm_load_ps(b)));
      d2 = _mm256_cvtps_pd(_mm_and_ps(tf0, _mm_load_ps(c)));
      d0 = _mm256_mul_pd(d0, d1);
      td0 = _mm256_fmadd_pd(d0, d2, td0);
    }

    a+=res ; b+= res ; c+= res ;
    for(i = res ; i+7<ni; i+=8){                    // 8 elements per loop
      d0 = _mm256_cvtps_pd(_mm_load_ps(a+0));
      d1 = _mm256_cvtps_pd(_mm_load_ps(b+0));
      d2 = _mm256_cvtps_pd(_mm_load_ps(c+0));
      d0 = _mm256_mul_pd(d0, d1);
      td0 = _mm256_fmadd_pd(d0, d2, td0);

      d0 = _mm256_cvtps_pd(_mm_load_ps(a+4));
      d1 = _mm256_cvtps_pd(_mm_load_ps(b+4));
      d2 = _mm256_cvtps_pd(_mm_load_ps(c+4));
      d0 = _mm256_mul_pd(d0, d1);
      td1 = _mm256_fmadd_pd(d0, d2, td1);

      a+=8 ; b+=8 ; c+=8 ;
    }
    td0 = _mm256_add_pd(td0, td1);

    if(i+3 < ni){                  // last 4 elements if need be
      d0 = _mm256_cvtps_pd(_mm_load_ps(a));
      d1 = _mm256_cvtps_pd(_mm_load_ps(b));
      d2 = _mm256_cvtps_pd(_mm_load_ps(c));
      d0 = _mm256_mul_pd(d0, d1);
      td0 = _mm256_fmadd_pd(d0, d2, td0);
      a+=4 ; b+=4 ; c+=4 ;
    }
    a += (li - ni) ; b += (li - ni) ; c += (li - ni) ; // jump between rows
  }

  f0 = _mm256_extractf128_pd(td0, 0);        // f0 = td0[0] td0[1]
  f1 = _mm256_extractf128_pd(td0, 1);        // f1 = td0[2] td0[3]
  f0 = _mm_add_pd(f0, f1);                   // f0 = (td0[0]+td0[2]) (td0[1]+td0[3])
  f1 = _mm_permute_pd(f0, 3);                // f1 = f0[1] f0[1]
  f0 = _mm_add_sd(f0, f1);                   // f0[0] = f0[0] + f1[0] = f0[0] + f0[1]
  _mm_store_sd(&result,f0);  // convert to single precision
#else
  for( ; nj > 0 ; nj--){  // loop over rows
    for(i=0 ; i<ni ; i++){
      result += a[i] * b[i] * c[i];
    }
    a += li ; b += li ; c += li ; // jump between rows
  }
#endif
  return result;
}

#endif

// aa, bb are 3D arrays (i,j,k)
// cc     is a 2D array (i,j)
double dot_product_332d(int ni, int nj, int nk, int li, int lij, float *aa, float *bb, float *cc){
  double result = 0;
  int i, k;
  float *a, *b, *c, *a0, *b0, *c0;

#if defined(__AVX2__) & defined(__FMA__)
  int res = ni & 3;
  __m256d d0, d1, d2, d3, td0, td1, td2, td3;
  __m128d f0, f1;
  __m128  tf0;

  td0 = _mm256_xor_pd(td0,td0);
  td1 = _mm256_xor_pd(td1,td1);
  tf0 = _mm_load_ps((float *) (mask32 + 4 + 4 - res) );   // first ni modulo 4 elements if modulo nonzero

  for( ; nj > 0 ; nj--){          // loop over rows
    a0 = aa ; b0 = bb ; c0 = cc;  // row base
    for(k = 0 ; k < nk ; k++ ){   // loop over planes
      a = a0 ; b = b0 ; c = c0 ;  // base addresses for the row

//       if(res > 0){
      d0 = _mm256_cvtps_pd(_mm_and_ps(tf0, _mm_load_ps(a)));  // mask unused elements (all if moduo is zero)
      d1 = _mm256_cvtps_pd(_mm_and_ps(tf0, _mm_load_ps(b)));
      d2 = _mm256_cvtps_pd(_mm_and_ps(tf0, _mm_load_ps(c)));
      d0 = _mm256_mul_pd(d0, d1);
      td0 = _mm256_fmadd_pd(d0, d2, td0);
//       }

      a += res ; b += res ; c += res ;                 // skip over modulo(ni,4) elements
      for(i = res ; i+15 < ni; i += 16){               // 16 elements per loop (4 way SIMD unroll)
	d0 = _mm256_cvtps_pd(_mm_load_ps(a+0));
	d1 = _mm256_cvtps_pd(_mm_load_ps(b+0));
	d2 = _mm256_cvtps_pd(_mm_load_ps(c+0));
	d0 = _mm256_mul_pd(d0, d1);                    // a[i]*b[i]
	td0 = _mm256_fmadd_pd(d0, d2, td0);            // += (a[i]*b[i]) * c[i]

	d0 = _mm256_cvtps_pd(_mm_load_ps(a+4));
	d1 = _mm256_cvtps_pd(_mm_load_ps(b+4));
	d2 = _mm256_cvtps_pd(_mm_load_ps(c+4));
	d0 = _mm256_mul_pd(d0, d1);
	td1 = _mm256_fmadd_pd(d0, d2, td1);

	d0 = _mm256_cvtps_pd(_mm_load_ps(a+8));
	d1 = _mm256_cvtps_pd(_mm_load_ps(b+8));
	d2 = _mm256_cvtps_pd(_mm_load_ps(c+8));
	d0 = _mm256_mul_pd(d0, d1);
	td0 = _mm256_fmadd_pd(d0, d2, td0);

	d0 = _mm256_cvtps_pd(_mm_load_ps(a+12));
	d1 = _mm256_cvtps_pd(_mm_load_ps(b+12));
	d2 = _mm256_cvtps_pd(_mm_load_ps(c+12));
	d0 = _mm256_mul_pd(d0, d1);
	td1 = _mm256_fmadd_pd(d0, d2, td1);

	a+=16 ; b+=16 ; c+=16 ;                       // i loop (next group of 16)
      }
      td0 = _mm256_add_pd(td0, td1);

      for(i = res ; i+3<ni; i+=4){                    // remainder at 4 elements per loop
	d0 = _mm256_cvtps_pd(_mm_load_ps(a));
	d1 = _mm256_cvtps_pd(_mm_load_ps(b));
	d2 = _mm256_cvtps_pd(_mm_load_ps(c));
	d0 = _mm256_mul_pd(d0, d1);
	td0 = _mm256_fmadd_pd(d0, d2, td0);
	a+=4 ; b+=4 ; c+=4 ;                          // i loop (last groups of 4 if needed)
      }
      a0 += lij ; b0 += lij ;        // bump row base to next plane (k loop)
    }
    aa += li ; bb += li ; cc += li ; // bump base row to next row (j loop)
  }
  f0 = _mm256_extractf128_pd(td0, 0);        // f0 = td0[0] td0[1]
  f1 = _mm256_extractf128_pd(td0, 1);        // f1 = td0[2] td0[3]
  f0 = _mm_add_pd(f0, f1);                   // f0 = (td0[0]+td0[2]) (td0[1]+td0[3])
  f1 = _mm_permute_pd(f0, 3);                // f1 = f0[1] f0[1]
  f0 = _mm_add_sd(f0, f1);                   // f0[0] = f0[0] + f1[0] = f0[0] + f0[1]
  _mm_store_sd(&result,f0);  // convert to single precision
#else
  for( ; nj > 0 ; nj--){  // loop over rows
    for(i=0 ; i<ni ; i++){
      result += a[i] * b[i] * c[i];
    }
    a += li ; b += li ; c += li ; // jump between rows
  }
#endif
  return result;
}
