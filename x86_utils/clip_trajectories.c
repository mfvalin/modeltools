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

// interface
//   subroutine SetClippingLimits(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj) bind(C,name='SetClippingLimits')
//     import :: C_FLOAT
//     implicit none
//     real(C_FLOAT), intent(IN), value :: cmini, cmaxi, cminj, cmaxj  ! cell (tile) limits
//     real(C_FLOAT), intent(IN), value :: gmini, gmaxi, gminj, gmaxj  ! grid limits
//   end subroutine SetClippingLimits
// end interface
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
// interface
//   subroutine ClipTrajectories(alpha, beta, l1, l2, ni, indx1, indx2) bind(C,name='ClipTrajectories')
//     import :: C_INT, C_FLOAT
//     implicit none
//     real(C_FLOAT), dimension(ni), intent(IN) :: alpha   ! tartgets along x (i)
//     real(C_FLOAT), dimension(ni), intent(IN) :: beta    ! tartgets along y (j)
//     integer(C_INT), dimension(*), intent(OUT) :: l1     ! clipping list for tile (cell)
//     integer(C_INT), dimension(*), intent(OUT) :: l2     ! clipping list for grid
//     integer(C_INT), intent(IN), value :: ni             ! dimension of alpha and beta
//     integer(C_INT), intent(OUT) :: indx1                ! number of points clipped in l1 list
//     integer(C_INT), intent(OUT) :: indx2                ! number of points clipped in l2 list
//   end subroutine ClipTrajectories
// end interface
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

#if ! defined(NPTS)
#define NPTS 1280007
#endif

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
  ClipTrajectories(alpha, beta, l1, l2, NPTS, &indx1, &indx2) ;
  t1 = LOCAL_time() ;
  printf("indx1 = %d, indx2 = %d, count = %d, t = %8.3f\n",indx1, indx2, ii, (t1 - t0)/NPTS * 1.0E9);
  if(indx1 < 10){
    for(i=0 ; i<10 && i < indx1 ; i++) printf(" %8d",l1[i]) ; printf(" indx1\n") ;
  }
  if(indx2 < 10){
    for(i=0 ; i<10 && i < indx2 ; i++) printf(" %8d",l2[i]) ; printf(" indx2\n") ;
  }
}

#endif
