/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2019  Division de Recherche en Prevision Numerique
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
#include <stdint.h>

static double cp133 =  1.0/3.0;
static double cp167 =  1.0/6.0;
static double cm167 = -1.0/6.0;
static double cp5 =  .5;
static double cm5 = -.5;
static double one = 1.0;
static double two = 2.0;

#if defined(__x86_64__)
#include <immintrin.h>
// use separate multiply and add instructions if fused multiply-add not available
#if defined(__AVX__) && ! defined(__FMA__)
#define _mm256_fmadd_ps(a,b,c) _mm256_add_ps(_mm256_mul_ps(a,b),c)
#define _mm256_fmadd_pd(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b),c)
#define _mm_fmadd_ps(a,b,c) _mm_add_ps(_mm_mul_ps(a,b),c)
#define _mm_fmadd_pd(a,b,c) _mm_add_pd(_mm_mul_pd(a,b),c)
#endif
#endif

// compute 2D X Y cubic interpolation coefficients (Origin 0 internally)
// assuming CONSTANT spacing along x, CONSTANT spacing along y
//
// xx    : position along x in index space of the interpolation target, input (origin 1)
// yy    : position along x in index space of the interpolation target, input (origin 1)
// wx    : output array, 4 elements, coefficients for cubic interpolation along x
// wy    : output array, 4 elements, coefficients for cubic interpolation along y
// qminx, qmaxx : limits along x used as box corner constraints
// qminy, qmaxy : limits along y used as box corner constraints
// offsetx  : offset applied to xx (global to local optional mapping) (usually 0)
// offsety  : offset applied to yy (global to local optional mapping) (usually 0)
//
// function result : displacement from point(1,1) in array to lower left corner of 4 x 4 subarray
//                   used for the bicubic interpolation
//
// inline version, used by multiple routines, the optimizer will perform the appropriate inlining
//
static inline int bicub_coeffs_inline0(double xx, double yy, int ni, double *wx, double *wy,
                                  int qminx, int qmaxx, int qminy, int qmaxy, int offsetx, int offsety){
  double x, y;
  int ix, iy;

  ix = xx ; 
  if(ix > xx) ix = ix -1 ;              // ix > xx if xx is negative
  ix = ix - 2;
  ix = (ix < qminx) ? qminx : ix;       // apply boundary limits on ix and iy
  ix = (ix > qmaxx) ? qmaxx : ix;
  x  = xx - 2 - ix;

  iy = yy ; 
  if(iy > yy) iy = iy -1 ;              // iy > yy if yy is negative
  iy = iy - 2;
  iy = (iy < qminy) ? qminy : iy;
  iy = (iy > qmaxy) ? qmaxy : iy;
  y  = yy - 2 - iy;

  wx[0] = cm167*x*(x-one)*(x-two);      // cubic interpolation coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm167*y*(y-one)*(y-two);      // cubic interpolation coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

  return ix + iy * ni;   // displacement to lower left corner of 4 x 4 subarray from point(1,1)
}

// compute 2D X Y cubic interpolation coefficients (Origin 1 intenally)
// assuming CONSTANT spacing along x, CONSTANT spacing along y
//
// xx    : position along x in index space of the interpolation target, input (origin 1)
// yy    : position along x in index space of the interpolation target, input (origin 1)
// wx    : output array, 4 elements, coefficients for cubic interpolation along x
// wy    : output array, 4 elements, coefficients for cubic interpolation along y
// qminx, qmaxx : limits along x used as box corner constraints
// qminy, qmaxy : limits along y used as box corner constraints
// offsetx  : offset applied to xx (global to local optional mapping) (usually 0)
// offsety  : offset applied to yy (global to local optional mapping) (usually 0)
//
// function result : displacement from point(1,1) in array to lower left corner of 4 x 4 subarray
//                   used for the bicubic interpolation
//
// inline version, used by multiple routines, the optimizer will perform the appropriate inlining
//
static inline int bicub_coeffs_inline1(double xx, double yy, int ni, double *wx, double *wy,
                                  int qminx, int qmaxx, int qminy, int qmaxy, int offsetx, int offsety){
  double x, y;
  int ix, iy;

  ix = xx ; 
  ix = (ix < qminx) ? qminx : ix;       // apply boundary limits on ix
  ix = (ix > qmaxx) ? qmaxx : ix;
  x  = xx - ix;
  ix - ix + offsetx;

  iy = yy ; 
  iy = (iy < qminy) ? qminy : iy;       // apply boundary limits on iy
  iy = (iy > qmaxy) ? qmaxy : iy;
  y  = yy - iy;
  iy = iy + offsety

  wx[0] = cm167*x*(x-one)*(x-two);      // cubic interpolation coefficients along i
  wx[1] = cp5*(x+one)*(x-one)*(x-two);
  wx[2] = cm5*x*(x+one)*(x-two);
  wx[3] = cp133*x*(x+one)*(x-one);

  wy[0] = cm167*y*(y-one)*(y-two);      // cubic interpolation coefficients along j
  wy[1] = cp5*(y+one)*(y-one)*(y-two);
  wy[2] = cm5*y*(y+one)*(y-two);
  wy[3] = cp133*y*(y+one)*(y-one);

  return (ix-1) + (iy-1) * ni;   // displacement to lower left corner of 4 x 4 subarray from point(1,1)
}

// compute 3D X Y Z cubic/linear interpolation coefficients (Origin 1)
// assuming CONSTANT spacing along x, CONSTANT spacing along y, VARIABLE spacing along z
// function value is true is interpolation along z will be linear, false otherwise
static inline int tricub_coeffs3_inline(double *wxyz, int *offset, float px8, float py8, float pz8, ztab *lv){
  int ix, iy, iz, ijk, zlinear;
  double pxy[2], *base, *pos;
  double zza, zzb, zzc, zzd, zzab, zzcd, dz, px, py, pz;
  int i, j;

    px = px8 ;                      // fractional index positions along x, y, z (float to double)
    py = py8 ;
    pz = pz8 ;

    ix = px ;
    px = px - ix;                   // px is now deltax (fractional part of px)
    ix = ix + lv->offi;             // global to local grid remapping
    ijk = ix - 2;                   // x displacement (elements), ix assumed to always be >1 and < ni-1

    wxyz[12] = 0.0;                 // linear interpolation coefficients along x
    wxyz[13] = 1.0 - px;            // 4 coefficients needed because of SIMD interpolation along x
    wxyz[14] = px;                  // first and last coefficients MUST be 0.0
    wxyz[15] = 0.0;

    iy = py ;
    py = py - iy;                   // py is now deltay (fractional part of py)
    iy = iy + lv->offj;             // global to local grid remapping
    ijk = ijk + (iy - 2) * lv->ni;  // add y displacement (rows), ix assumed to always be >1 and < nj-1

    wxyz[16] = 1.0 - py;            // linear interpolation coefficients along y
    wxyz[17] = py;

    iz = pz ; 
    if(iz<1) iz = 1; 
    if(iz>lv->nk-1) iz = lv->nk-1;  // iz < 1 or iz > nk-1 will result in linear extrapolation
    dz = pz - iz;                   // dz is now "fractional" part of pz  (may be <0 or >1 if extrapolating)
    ijk = ijk + (iz -1) * lv->nij;  // add z displacement (2D planes)
    wxyz[18] = 1.0 - dz;            // linear interpolation coefficients along z
    wxyz[19] = dz;

    iz--;                           // iz needs to be in "origin 0" (C index from Fortran index)
    zlinear = (iz - 1) | (lv->nk - 3 - iz); 
    zlinear >>= 31;                 // nonzero only if iz < 1 or iz > nk -3 (top and bottom intervals)
    if(! zlinear) ijk = ijk - lv->nij;  // not the linear case, go down one 2D plane to get lower left corner of 4x4x4 cube
    *offset = ijk;

    // now we can compute the coefficients along z using iz and dz
    if(zlinear){
      wxyz[ 8] = 0.0;                    // coefficients for linear interpolation along z
      wxyz[ 9] = 1.0 - dz;
      wxyz[10] = dz;
      wxyz[11] = 0.0;
    }else{
      base  = &(lv->ocz[4*iz]);  // precomputed inverses of denominators
      pos   = &(lv->z[iz]);
      pz  = dz * pos[1] + (1.0 - dz) * pos[0];   // pz is now an absolute position
      zza = pz - pos[-1] ; zzb = pz - pos[0] ; zzc = pz - pos[1] ; zzd = pz - pos[2] ; 
      zzab = zza * zzb ; zzcd = zzc * zzd;
      // printf("target = %8.5f, dz = %8.5f, levels = %8.5f %8.5f\n",*PZ,dz,pos[0],pos[1]);
      // printf("target = %8.5f, base = %8.5f %8.5f %8.5f %8.5f\n",pz,base[0],base[1],base[2],base[3]);
      // printf("dz     = %8.5f, zz   = %8.5f %8.5f %8.5f %8.5f\n",dz,zza,zzb,zzc,zzd);
      // printf("iz     = %d, lv->odz[iz] = %8.5f\n",iz,lv->odz[iz]);
      wxyz[ 8] = zzb * zzcd * base[0];   //   wxyz[16] = TRIPRD(pz,pos[1],pos[2],pos[3]) * base[0];
      wxyz[ 9] = zza * zzcd * base[1];   //   wxyz[17] = TRIPRD(pz,pos[0],pos[2],pos[3]) * base[1];
      wxyz[10] = zzd * zzab * base[2];   //   wxyz[18] = TRIPRD(pz,pos[0],pos[1],pos[3]) * base[2];
      wxyz[11] = zzc * zzab * base[3];   //   wxyz[19] = TRIPRD(pz,pos[0],pos[1],pos[2]) * base[3];
    }
//     alternative formulation using independent FMAs and products
//     wxyz[ 0] = (px*px  - px)   * (-(cp167*px) + cp133);
//     wxyz[ 1] = (cp5*px - one)  * (px*px       - one);
//     wxyz[ 2] = (px*px  + px)   * (-(cp5*px)   + one);
//     wxyz[ 3] = (px*px  + px)   * (cp167*px    - cp167);
// 
//     wxyz[ 8] = (py*py  - py)   * (-(cp167*py) + cp133);
//     wxyz[ 9] = (cp5*py - one)  * (py*py       - one);
//     wxyz[10] = (py*py  + py)   * (-(cp5*py)   + one);
//     wxyz[11] = (py*py  + py)   * (cp167*py    - cp167);

    wxyz[ 0] = cm167*px*(px-one)*(px-two);        // coefficients for cubic interpolation along x
    wxyz[ 1] = cp5*(px+one)*(px-one)*(px-two);
    wxyz[ 2] = cm5*px*(px+one)*(px-two);
    wxyz[ 3] = cp167*px*(px+one)*(px-one);

    wxyz[ 4] = cm167*py*(py-one)*(py-two);        // coefficients for cubic interpolation along y
    wxyz[ 5] = cp5*(py+one)*(py-one)*(py-two);
    wxyz[ 6] = cm5*py*(py+one)*(py-two);
    wxyz[ 7] = cp167*py*(py+one)*(py-one);

  return zlinear;  // linear / cubic flag
}
