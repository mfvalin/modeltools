/*
 * Copyright (C) 2021  Environnement et Changement climatique Canada
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
 * Author:
 *     M. Valin,   Recherche en Prevision Numerique, 2021
 */

#include <stdint.h>

static float c_mini = 0.0f ;
static float c_maxi = 0.0f ;
static float c_minj = 0.0f ;
static float c_maxj = 0.0f ;
static float g_mini = 0.0f ;
static float g_maxi = 0.0f ;
static float g_minj = 0.0f ;
static float g_maxj = 0.0f ;

void SetClippingLimits(float cmini, float cmaxi, float cminj, float cmaxj , float gmini, float gmaxi, float gminj, float gmaxj){
  c_mini = cmini ;
  c_maxi = cmaxi ;
  c_minj = cminj ;
  c_maxj = cmaxj ;
  g_mini = gmini ;
  g_maxi = gmaxi ;
  g_minj = gminj ;
  g_maxj = gmaxj ;
}

#define VL 16

void ClipTrajectories(float *alpha, float *beta, int *l1, int *l2, int ni, int *indx1, int *indx2){
  int i, k, ix1, ix2 ;
  int clip1[VL], clip2[VL] ;
  int vlmod = ni & (VL - 1) ;
  float *a = alpha ;
  float *b = beta ;

  ix1 = *indx1 ;
  ix2 = *indx2 ;
  for(i=0 ; i+VL-1<ni ; i+=VL){
    for(k=0 ; k<VL ; k++){
      clip2[k]   = (a[k] < g_mini) ? -1 : 0 ;  // a < g_mini
      clip2[k]  |= (a[k] > g_maxi) ? -1 : 0 ;  // a > g_maxi
      clip2[k]  |= (b[k] < g_minj) ? -1 : 0 ;  // b < g_minj
      clip2[k]  |= (b[k] > g_maxj) ? -1 : 0 ;  // b > g_maxj

      clip1[k]   = (a[k] < c_mini) ? -1 : 0 ;  // a < c_mini
      clip1[k]  |= (a[k] > c_maxi) ? -1 : 0 ;  // a > c_maxi
      clip1[k]  |= (b[k] < c_minj) ? -1 : 0 ;  // b < c_minj
      clip1[k]  |= (b[k] > c_maxj) ? -1 : 0 ;  // b > c_maxj
      // if outside both cell and domain bounds, mark as outside domain
      clip1[k] = clip1[k] & (~clip2[k]) ;
    }
    for(k = 0 ; k < VL ; k++){
      l1[ix1] = i + k ;  // store index, just in case
      ix1 -= clip1[k] ;  // bump pointer if there was clipping
      l2[ix2] = i + k ;  // store index, just in case
      ix2 -= clip2[k] ;  // bump pointer if there was clipping
    }
    a += VL ;
    b += VL ;
  }
  for(k=0 ; k<vlmod ; k++){      // leftovers (at most VL -1 points)
    clip2[k]   = (a[k] < g_mini) ? -1 : 0 ;  // a < g_mini
    clip2[k]  |= (a[k] > g_maxi) ? -1 : 0 ;  // a > g_maxi
    clip2[k]  |= (b[k] < g_minj) ? -1 : 0 ;  // b < g_minj
    clip2[k]  |= (b[k] > g_maxj) ? -1 : 0 ;  // b > g_maxj

    clip1[k]   = (a[k] < c_mini) ? -1 : 0 ;  // a < c_mini
    clip1[k]  |= (a[k] > c_maxi) ? -1 : 0 ;  // a > c_maxi
    clip1[k]  |= (b[k] < c_minj) ? -1 : 0 ;  // b < c_minj
    clip1[k]  |= (b[k] > c_maxj) ? -1 : 0 ;  // b > c_maxj
    // if outside both cell and domain bounds, mark as outside domain
    clip1[k] = clip1[k] & (~clip2[k]) ;
    l1[ix1] = i + k ;  // store index, just in case
    ix1 -= clip1[k] ;  // bump pointer if there was clipping
    l2[ix2] = i + k ;  // store index, just in case
    ix2 -= clip2[k] ;  // bump pointer if there was clipping
  }
  *indx1 = ix1 ;
  *indx2 = ix2 ;
}


#if defined(__AVX2__) && defined(__x86_64__) && defined(USE_SIMD)
#include <immintrin.h>
#define CMP_LT_OS 17 
#define CMP_GT_OS 14

static int mask[16] = { -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 } ;

int ClipTrajectoriesAVX2(float *alpha, float *beta, int *l1, int *l2, int ni, int *indx1, int *indx2){
  int i, k, ix1, ix2, count, active ;
  __m256 msk ;
  __m256 a, b, c1, c2 ;
  __m256 cmini, cmaxi, cminj, cmaxj ;
  __m256 gmini, gmaxi, gminj, gmaxj ;
  int xc1[8], xc2[8] ;
  int ni7 = ni & 7 ;

  ix1 = 0 ;
  ix2 = 0 ;
  count = 0 ;
  cmini = _mm256_broadcast_ss(&c_mini) ;
  cmaxi = _mm256_broadcast_ss(&c_maxi) ;
  cminj = _mm256_broadcast_ss(&c_minj) ;
  cmaxj = _mm256_broadcast_ss(&c_maxj) ;
  gmini = _mm256_broadcast_ss(&g_mini) ;
  gmaxi = _mm256_broadcast_ss(&g_maxi) ;
  gminj = _mm256_broadcast_ss(&g_minj) ;
  gmaxj = _mm256_broadcast_ss(&g_maxj) ;

  // take care of mod(ni,8) first points if necessary
  if(ni7 > 0){
    msk   = _mm256_loadu_ps((float *) &mask[8 - ni7]) ;      // first ni7 points only
    a     = _mm256_loadu_ps(alpha) ;         // fetch alpah
    a     = _mm256_and_ps(a,msk) ;           // ignore masked points (force to 0)
    b     = _mm256_loadu_ps(beta) ;          // fetch beta
    b     = _mm256_and_ps(b, msk) ;          // ignore masked points (force to 0)
    // check against tile bounds
    c1    =                  _mm256_cmp_ps(a, cmini, CMP_LT_OS)  ;  // a < c_mini
    c1    = _mm256_or_ps(c1, _mm256_cmp_ps(b, cminj, CMP_LT_OS)) ;  // b < c_minj
    c1    = _mm256_or_ps(c1, _mm256_cmp_ps(a, cmaxi, CMP_GT_OS)) ;  // a > c_maxi
    c1    = _mm256_or_ps(c1, _mm256_cmp_ps(b, cmaxj, CMP_GT_OS)) ;  // b > c_maxj
    c1    = _mm256_and_ps(c1, msk) ;   // ignore masked points (force result to 0, no clipping)
    // check against domain bounds
    c2    =                  _mm256_cmp_ps(a, gmini, CMP_LT_OS)  ;  // a < g_mini
    c2    = _mm256_or_ps(c2, _mm256_cmp_ps(b, gminj, CMP_LT_OS)) ;  // b < g_minj
    c2    = _mm256_or_ps(c2, _mm256_cmp_ps(a, gmaxi, CMP_GT_OS)) ;  // a > g_maxi
    c2    = _mm256_or_ps(c2, _mm256_cmp_ps(b, gmaxj, CMP_GT_OS)) ;  // b > g_maxj
    c2    = _mm256_and_ps(c2, msk) ;   // ignore masked points (force result to 0, no clipping)
    active = _mm256_movemask_ps (_mm256_or_ps(c2,c1) ) ;
    // if outside both cell and domain bounds, mark as outside domain
    c1    = _mm256_andnot_ps (c2, c1) ;   // c1 = ~c2 & c1  set to zero if c2 is one

    if( active ){  // at least one condition is true
      _mm256_storeu_ps ((float *) xc1, c1);
      _mm256_storeu_ps ((float *) xc2, c2);
      for(k=0 ; k<8 ; k++){
        l1[ix1] = k ;
        l2[ix2] = k ;
        ix1 = ix1 - xc1[k] ;
        ix2 = ix2 - xc2[k] ;
      }
    }
    alpha += ni7 ;  // bump source pointers by effective loop count
    beta  += ni7 ;
  }

  for(i=ni7 ; i+7< ni ; i+=8){

    a     = _mm256_loadu_ps(alpha) ;    // fetch alpha
    b     = _mm256_loadu_ps(beta) ;     // fetch beta
    // check against tile bounds
    c1    =                  _mm256_cmp_ps(a, cmini, CMP_LT_OS)  ;  // a < c_mini
    c1    = _mm256_or_ps(c1, _mm256_cmp_ps(b, cminj, CMP_LT_OS)) ;  // b  < c_minj
    c1    = _mm256_or_ps(c1, _mm256_cmp_ps(a, cmaxi, CMP_GT_OS)) ;  // a > c_maxi
    c1    = _mm256_or_ps(c1, _mm256_cmp_ps(b, cmaxj, CMP_GT_OS)) ;  // b  > c_maxj
    // check against domain bounds
    c2    =                  _mm256_cmp_ps(a, gmini, CMP_LT_OS)  ;  // a < g_mini
    c2    = _mm256_or_ps(c2, _mm256_cmp_ps(b, gminj, CMP_LT_OS)) ;  // b  < g_minj
    c2    = _mm256_or_ps(c2, _mm256_cmp_ps(a, gmaxi, CMP_GT_OS)) ;  // a > g_maxi
    c2    = _mm256_or_ps(c2, _mm256_cmp_ps(b, gmaxj, CMP_GT_OS)) ;  // b  > g_maxj
    active = _mm256_movemask_ps (_mm256_or_ps(c2,c1) ) ;
    // if outside both cell and domain bounds, mark as outside domain
    c1    = _mm256_andnot_ps (c2, c1) ;   // c1 = ~c2 & c1  set to zero if c2 is one

    if( active ){  // at least one condition is true
      _mm256_storeu_ps ((float *) xc1, c1);
      _mm256_storeu_ps ((float *) xc2, c2);
      count++ ;
      for(k=0 ; k<8 ; k++){      // process lists
        l1[ix1] = i+k ;          // cell clip list
        l2[ix2] = i+k ;          // global clip list
        ix1 = ix1 - xc1[k] ;
        ix2 = ix2 - xc2[k] ;
      }
    }
    alpha += 8 ;  // bump source pointers by loop count
    beta  += 8 ;
  }
  *indx1 = ix1 ;
  *indx2 = ix2 ;
  return count ;
}
#endif

#if defined(SELF_TEST)

#include <stdio.h>
#include <time.h>

static double LOCAL_time(){
  struct timespec tv ;
  double value ;
  clock_gettime(CLOCK_MONOTONIC, &tv);
  value = tv.tv_nsec * 1E-9 ;
  value += tv.tv_sec ;
  return value ;
}

#define NPTS 1280007

int main(int argc, char **argv){
  float alpha[NPTS], beta[NPTS] ;
  int i, ii ;
  int indx1, indx2 ;
  int l1[NPTS], l2[NPTS] ;
  double t0, t1 ;

  ii = 0 ;
  for(i=0 ; i<NPTS ; i++){
    alpha[ii] = (1.0f * i) / (NPTS - 1) ;
    beta[ii] = (1.0f * i) / (NPTS - 1) ;
    ii = ii + 41 ;                     // 41 is alpha prime number
    if(ii >= NPTS) ii = ii - NPTS ;
  }
//   printf(" %f %f \n",alpha[0],alpha[NPTS-1]);
//   for(i=0 ; i<NPTS ; i++) printf(" %4.2f",alpha[i]) ; printf("\n") ;

  SetClippingLimits(0.06f, 0.94f, 0.06f, 0.94f, 0.03f, .97f, 0.03f, .97f) ;
  t0 = LOCAL_time() ;
  indx1 = 0 ;
  indx2 = 0 ;
#if defined(USE_SIMD)
  ii = ClipTrajectoriesAVX2(alpha, beta, l1, l2, NPTS, &indx1, &indx2) ;
#else
  ClipTrajectories(alpha, beta, l1, l2, NPTS, &indx1, &indx2) ;
#endif
  t1 = LOCAL_time() ;
  printf("indx1 = %d, indx2 = %d, count = %d, t = %8.3f\n",indx1, indx2, ii, (t1 - t0)/NPTS * 1.0E9);
//   for(i=0 ; i<10 && i < indx1 ; i++) printf(" %8d",l1[i]) ; printf(" indx1\n") ;
//   for(i=0 ; i<10 && i < indx2 ; i++) printf(" %8d",l2[i]) ; printf(" indx2\n") ;
}

#endif
