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

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <immintrin.h>
#include <float.h>

typedef union{
  double  d;
  int64_t l;
}d_l_p;

typedef struct{
  double t0[ 8];   // search accelerator for t1 (64 entries max in t1)
  double t1[64];   // search accelerator for t2 (at most 256 levels in t2)
  double t2m1;     // normally equal to t2[0] (guard value)
  double t2[257];  // table for levels (max 256 usable), entry nk is a guard
  double top;      // t2[nk-2], next to last value used for research clipping
  double *odz;     // pointer to inverse of deltaz array [nk doubles]
  double *ocz;     // pointer to inverse of coefficient denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of levels (max 255)
  uint32_t nij;    // ni*nj, distance between levels
}lvtab;

// constants used for bicubic interpolation
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
static int64_t ONE = 1;
static double cp133 =  1.0/3.0;
#endif
static double cp167 =  1.0/6.0;
static double cm167 = -1.0/6.0;
static double cm5 = -0.5;
static double cp5 = 0.5;
static double one = 1.0;
static double two = 2.0;

// triple product used for Lagrange cubic polynomials coefficients in the non constant case
#define TRIPRD(x,a,b,c) ((x-a)*(x-b)*(x-c))

// inverse denominators r from positions a b c d
static inline void denominators(double *r, double a, double b, double c, double d){
  r[0] = 1.0 / TRIPRD(a,b,c,d);
  r[1] = 1.0 / TRIPRD(b,a,c,d);
  r[2] = 1.0 / TRIPRD(c,a,b,d);
  r[3] = 1.0 / TRIPRD(d,a,b,c);
}

// coefficients c from inverse denominators d, levels z, and position t
static inline void zcoeffs(double *c, double t, double *z, double *d){
  c[0] = TRIPRD(t,z[1],z[2],z[3]) * d[0];
  c[1] = TRIPRD(t,z[0],z[2],z[3]) * d[1];
  c[2] = TRIPRD(t,z[0],z[1],z[3]) * d[2];
  c[3] = TRIPRD(t,z[0],z[1],z[2]) * d[3];
}

// allocate lookup table set and
// return pointer to filled table set
// targets are expected to be positive, and monotonically increasing or decreasing
lvtab * Vsearch_setup(double *targets, int n, int ni, int nj){
  int i;
  lvtab *lv;
  double pad;

  if(n > 256) return NULL;
  if(0 != posix_memalign( (void **) &lv, 32, sizeof(lvtab) ) ) return NULL ;

  pad = DBL_MAX;                                        // increasing values, pad with large positive number
  if(targets[1] < targets[0]) pad = -pad;               // decreasing values, pad with large negative number
  for(i=0 ; i<n   ; i++) { lv->t2[i] = targets[i] ; }   // z coordinate table
  for(i=n ; i<256 ; i++) { lv->t2[i] = pad ; }          // pad at end of table
  lv->t2[n] = 2*lv->t2[n-1] - lv->t2[n-2];
  lv->t2m1 = 2*targets[0] - targets[1];

  for(i=0  ; i<64 ; i++) { lv->t1[i] = lv->t2[(i<<2)] ; } ;   // lookup for t2 (every fourth entry)
  for(i=0  ; i<8  ; i++) { lv->t0[i] = lv->t1[(i<<3)] ; } ;   // lookup for t1 (every eighth entry)
  // inverse of intervals ( 0 to n-2 only make sense )
  lv->odz = malloc(n * sizeof(double));
  if(NULL == lv->odz){    // malloc failed
    free(lv);             // deallocate lv
    return NULL;
  }
  for(i=0 ; i<n-1   ; i++) { lv->odz[i] = (targets[i+1] - targets[i]) ; }
  lv->odz[n - 1] = 0.0;   // NEVER USED
  // denominators for Lagrange cubic polynomials coefficients ( entries 1 to n-2 make sense )
  // entries 0 and n-1 are fudged
  lv->ocz = malloc(4 * n * sizeof(double));
  if(NULL == lv->ocz){    // malloc failed
    free(lv->odz);        // deallocate odz
    free(lv);             // deallocate lv
    return NULL;
  }
//   denominators( &(lv->ocz[0]) , targets[0], targets[1], targets[2], targets[3]);                 // level 0 coeffs are normally not used
  for(i=0 ; i<n-1   ; i++) {
    denominators( &(lv->ocz[4*i]) , lv->t2[i-1], lv->t2[i  ], lv->t2[i+1], lv->t2[i+2]);
  }
// i = n-2;
// printf("DE : %8.5f %8.5f %8.5f %8.5f \n",lv->t2[i-1], lv->t2[i  ], lv->t2[i+1], lv->t2[i+2]);
// i = 4*(n-3);
// printf("OV : %8.5f %8.5f %8.5f %8.5f \n",lv->ocz[i-1], lv->ocz[i  ], lv->ocz[i+1], lv->ocz[i+2]);
// i = 4*(n-2);
// printf("OV : %8.5f %8.5f %8.5f %8.5f \n",lv->ocz[i-1], lv->ocz[i  ], lv->ocz[i+1], lv->ocz[i+2]);
//   denominators( &(lv->ocz[4*(n-1)]) , targets[n-4], targets[n-3], targets[n-2], targets[n-1]);   // level n-1 coeffs are normally not used

  lv->top = lv->t2[n-2]; // next to last value along z
  lv->ni = ni;           // nb of points along x
  lv->nj = nj;           // nb of points along y
  lv->nk = n;            // nb of points along z
  lv->nij = ni*nj;       // ni * nj
  return lv;             // return pointer to filled table
}

// free lookup table set created by Vsearch_setup
int Vsearch_free(lvtab *lv){
  int errors = 0;
  if(lv == NULL) return -1;
  if(lv->odz == NULL) { 
    errors++;
  }else{ 
    free(lv->odz);
  }
  if(lv->ocz == NULL) { 
    errors++;
  }else{ 
    free(lv->ocz);
  }
  free(lv);
  return errors;
}

// NOTES: 
//      min(x,y) :  y + ( (x - y) & ( (x - y) >> 63 ) )
//      max(x,y) :  x - ( (x - y) & ( (x - y) >> 63 ) )
//      _mm256_blendv_pd (__m256d a, __m256d b, __m256d mask)  a if masksign == 0 , b if == 1
//      max(x,y) = blendv(x , y, x-y)
//      min(x,y) = blendv(y , x, x-y)

// version for monotonically increasing positive table values
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
// returns position ix of target in level list lv->t2  (ix has origin 0)
// such that lv->t2[ix] <= target <= lv->t2[ix+1]
// ix will never be negative nor larger than nk - 2
//
// inline function used by Vsearch_list_inc and Vsearch_list_inc_n
static inline int Vsearch_list_inc_inline(double target, lvtab *lv){
  int ix;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m256i v0, v1, vc, t;
  int m0, m1;

  t     = (__m256i) _mm256_broadcast_sd(&target);
  vc    = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  if(target < lv->t0[0]) t =(__m256i)  _mm256_broadcast_sd(&lv->t0[0]);    // target < first element in table
  if(target > lv->top  ) t =(__m256i)  _mm256_broadcast_sd(&lv->top  );    // target > next to last element in table

  t    = _mm256_add_epi64((__m256i) t, vc);   // we want target < tbl[i] so we add 1 to target

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[0]));         // 8 values to scan from table t0
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);  // t - tbl >= 0 if (target-1) < tbl
  v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);            // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  ix   = (_mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) << 3;         // index into t1 table for 8 values scan

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix  ]));      // 8 values to scan from table t1
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);
  v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  ix   = (ix + _mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) <<2 ;    // index into t2 table for final 4 value scan

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t2[ix  ]));      // only 4 values to scan from table t2
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);
  ix = ix + _mm_popcnt_u32(m0) - 1 ;
#else
  int i, j;
  d_l_p dlt, dlr, dlm;

  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l < dlr.l) target = dlr.d;  // target > first element in table
  if(dlt.l > dlm.l) target = dlm.d;  // target < next to last element in table
//   if(target < lv->t0[0]) target = lv->t0[0];  // target < first element in table
//   if(target > lv->top  ) target = lv->top  ;    // target > next to last element in table

  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t0[i]) j = i - 1; }
  ix = j << 3;
  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
  j = 3;
  for(i=3 ; i>0 ; i--) { if(target < lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
#endif
  return ix;
}
int Vsearch_list_inc(double target, lvtab *lv){
  return Vsearch_list_inc_inline(target, lv);
}
int Vsearch_list_inc_n(int *ix, double *target, lvtab *lv, int n){
  int i, ii;
  for(i=0 ; i<n ; i++){
    ii = Vsearch_list_inc_inline(target[i], lv);
    ix[i] = ii;
  }
  return ii;
}

// used to determine calling overhead
// double Vsearch_list_inc_do(double target, lvtab *lv){
//   return target;
// }

// same as Vsearch_list_inc but returns the position as a fractional index
// may be < 0 or > nk-1 if target is beyond table lv->t2 extreme values
double Vsearch_list_inc_d(double target, lvtab *lv){
  double temp;
  int ix = Vsearch_list_inc_inline(target, lv);
  temp = ix + ( (target - lv->t2[ix]) * lv->odz[ix] );  // (target - lowerlevel) / deltaz
  return temp;
}
double Vsearch_list_inc_dn(double *list, double *target, lvtab *lv, int n){
  int i, ix;
  double temp;
  for(i=0 ; i<n ; i++){
    temp = target[i];
    ix = Vsearch_list_inc_inline(temp, lv);
    temp = ix + ( (temp - lv->t2[ix]) * lv->odz[ix] );
    list[i] = temp;
  }
  return list[n-1];
}
#if defined(OLD_CODE)
double Vsearch_list_inc_d_orig(double target, lvtab *lv){
  int ix;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m256i v0, v1, vc, t;
  int m0, m1;

  t     = (__m256i) _mm256_broadcast_sd(&target);
  vc    = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  if(target < lv->t0[0]) t =(__m256i)  _mm256_broadcast_sd(&lv->t0[0]);    // target < first element in table
  if(target > lv->top  ) t =(__m256i)  _mm256_broadcast_sd(&lv->top  );    // target > next to last element in table

  t    = _mm256_add_epi64((__m256i) t, vc);   // we want target < tbl[i] so we add 1 to target

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[0]));                  // 8 values to scan from table t0
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);  // t - tbl >= 0 if (target-1) < tbl
  v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);              // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  ix   = (_mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) << 3;         // index into t1 table for 8 values scan

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix  ]));               // 8 values to scan from table t1
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);
  v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  ix   = (ix + _mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) <<2 ;             // index into t2 table for final 4 value scan

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t2[ix  ]));               // only 4 values to scan from table t2
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);
  ix = ix + _mm_popcnt_u32(m0) - 1 ;
#else
  int i, j;

  if(target < lv->t0[0]) target = lv->t0[0];  // target < first element in table
  if(target > lv->top  ) target = lv->top  ;    // target > next to last element in table

  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t0[i]) j = i - 1; }
  ix = j << 3;
  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
  j = 3;
  for(i=3 ; i>0 ; i--) { if(target < lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
#endif
  return ix + ( (target - lv->t2[ix]) * lv->odz[ix] );  // (target - lowerlevel) / deltaz
}

// used to measure call sequence overhead
int Vcoef_xyz_inc_o(uint32_t *ixyz, double *cxyz, double *px, double *py, double *pz, lvtab *lv, int n){
  return 0;
}
#endif

// calculate coefficients for x, y, z interpolation from lv table (see Vsearch_setup) and positions
// PX    : array of positions along x (dimension n) (fractional index space, "origin 1")
// PY    : array of positions along y (dimension n) (fractional index space, "origin 1")
// PZ    : array of positions along z (dimension n) (positions in same units as lv->t2)
// cxyz  : coefficients for tricubic or bicubic-linear interpolation (dimension 24*n)
//         cxyz[ 0- 7,I] : along x for point I
//         cxyz[ 8-15,I] : along y for point I
//         cxyz[16-23,I] : along z for point I
// ixyz  : array of indices along z
// lv    : x y z axis description  (see Vsearch_setup)
// function return : index of last point along z
int Vcoef_xyz_inc(uint32_t *ixyz, double *cxyz, double *PX, double *PY, double *PZ, lvtab *lv, int n){
  int ix, irep;
  double pxy[2], *base, *pos;
  double zza, zzb, zzc, zzd, zzab, zzcd, dz, px, py, pz;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m256i v0, v1, vc, t, ttop, tbot, ttemp;
  __m128d vxy, vc1, vc5, vc3, vp6, xxmx, xxpx, xxm1, x5p1, x6m3, x5m1, x6m6, dtmp;
  __m128d vr0, vr1, vr2, vr3;
  __m128i itmp;
  int m0, m1;
#else
  int i, j;
#endif
  for(irep=0 ; irep <n ; irep++){
    px = PX[irep] ; py = PY[irep] ; pz = PZ[irep] ;    // target positions along x, y, z

    ix = px ; px = px - ix;  // px is now dx (fractional part of px)
    ix = py ; py = py - ix;  // py is now dy (fractional part of py)

#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
    pxy[0] = px;
    pxy[1] = py;
    vxy  = _mm_loadu_pd(pxy);              // coefficients for x and y will be interleaved

    vc1  = _mm_set1_pd(one);               //  one    (1.0)
    vc5  = _mm_set1_pd(cp5);               //  cp5    (0.5)
    vc3  = _mm_set1_pd(cp133);             //  cp133  (1/3)
    vp6  = _mm_set1_pd(cp167);             //  cp167  (1/6)

    t     = (__m256i) _mm256_broadcast_sd(&pz);             // target
    vc    = (__m256i) _mm256_broadcast_sd((double *) &ONE); // integer constant 1
    tbot  = (__m256i) _mm256_broadcast_sd(&(lv->t0[0]));    // first element in table
    ttop  = (__m256i) _mm256_broadcast_sd(&(lv->top  ));    // next to last element in table
    if(pz < lv->t0[0]) t = tbot ;    // pz < first element in table, set to first element in table
    if(pz > lv->top  ) t = ttop;     // pz > next to last element in table et to next to last element
    
    t    = _mm256_add_epi64((__m256i) t, vc);   // we want pz < tbl[i] so we add 1 to pz (integer mode)

    v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[0]));    // 8 values to scan from table t0
    v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[4]));
    v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);            // t - tbl >= 0 if (pz-1) < tbl
    v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
    m0   = _mm256_movemask_pd((__m256d) v0);                      // transform subtract result signs into mask
    m1   = _mm256_movemask_pd((__m256d) v1);
    ix   = (_mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) << 3;    // index into t1 table for 8 values scan

    // alternate formula to compute coefficients taking advantage of independent FMAs
    //   pa = cm167*xy*(xy-one)*(xy-two)     = xy*(xy-one) * (-cp167)*(xy-two)  =  (xy*xy  - xy)   * (-(cp167*xy) + cp133)
    //   pb = cp5*(xy+one)*(xy-one)*(xy-two) = cp5*(xy-two) * (xy+one)*(xy-one) =  (cp5*xy - one)  * (xy*xy       - one)
    //   pc = cm5*xy*(xy+one)*(xy-two)       = xy*(xy+one) * (-cp5)*(xy-two)    =  (xy*xy  + xy)   * (-(cp5*xy)   + one)
    //   pd = cp167*xy*(xy+one)*(xy-one)     = xy*(xy+one) * cp167*(xy-one)     =  (xy*xy  + xy)   * (cp167*xy    - cp167)
    // STEP1 : 7 independent terms using 3 different FMAs  a*b+c, a*b-c, (-a*b)+c
    xxmx = _mm_fmsub_pd(vxy,vxy,vxy);      //  xy*xy       - xy
    x6m3 = _mm_fnmadd_pd(vxy,vp6,vc3);     //  -(cp167*xy) + cp133
    x5p1 = _mm_fmsub_pd(vc5,vxy,vc1);      //  cp5*xy      - one
    xxm1 = _mm_fmsub_pd(vxy,vxy,vc1);      //  xy*xy       - one
    xxpx = _mm_fmadd_pd(vxy,vxy,vxy);      //  xy*xy       + xy
    x5m1 = _mm_fnmadd_pd(vxy,vc5,vc1);     //  -(cp5*xy)   + one
    x6m6 = _mm_fmsub_pd(vxy,vp6,vp6);      //  cp167*xy    - cp167

    v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix  ]));               // 8 values to scan from table t1
    v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix+4]));
    v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);
    v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
    m0   = _mm256_movemask_pd((__m256d) v0);
    m1   = _mm256_movemask_pd((__m256d) v1);
    ix   = (ix + _mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) <<2 ;             // index into t2 table for final 4 value scan

    // STEP2 : multiply STEP 1 terms (4 independent operations)
    vr0  = _mm_mul_pd(xxmx,x6m3);          // coefficients for x and y are interleaved
    vr1  = _mm_mul_pd(x5p1,xxm1);
    vr2  = _mm_mul_pd(xxpx,x5m1);
    vr3  = _mm_mul_pd(xxpx,x6m6);

    v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t2[ix  ]));               // only 4 values to scan from table t2
    v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);
    m0   = _mm256_movemask_pd((__m256d) v0);
    ix = ix + _mm_popcnt_u32(m0) - 1 ;

    // final unshuffle to separate even terms (px) and odd terms (py) before storing them (independent operations)
    _mm_storeu_pd(cxyz   ,_mm_unpacklo_pd(vr0,vr1));  // cxyz[ 0: 1] = cx[0], cx[1]
    _mm_storeu_pd(cxyz+ 2,_mm_unpacklo_pd(vr2,vr3));  // cxyz[ 2: 3] = cx[2], cx[3]
    _mm_storeu_pd(cxyz+ 8,_mm_unpackhi_pd(vr0,vr1));  // cxyz[ 8: 9] = cy[0], cy[1]
    _mm_storeu_pd(cxyz+10,_mm_unpackhi_pd(vr2,vr3));  // cxyz[10:11] = cy[2], cy[3]

#else

    if(pz < lv->t0[0]) pz = lv->t0[0];  // pz < first element in table
    if(pz > lv->top  ) pz = lv->top  ;  // pz > next to last element in table

    j = 7;
    for(i=7 ; i>0 ; i--) { if(pz < lv->t0[i]) j = i - 1; }
    ix = j << 3;
    j = 7;
    for(i=7 ; i>0 ; i--) { if(pz < lv->t1[ix+i]) j = i - 1; }
    ix = (ix + j) <<2 ;
    j = 3;
    for(i=3 ; i>0 ; i--) { if(pz < lv->t2[ix+i]) j = i - 1; }
    ix = ix + j;           // we have the index along z, we can compute coefficients along z

    cxyz[ 0] = cm167*px*(px-one)*(px-two);        // coefficients for cubic interpolation along x
    cxyz[ 1] = cp5*(px+one)*(px-one)*(px-two);
    cxyz[ 2] = cm5*px*(px+one)*(px-two);
    cxyz[ 3] = cp167*px*(px+one)*(px-one);

    cxyz[ 8] = cm167*py*(py-one)*(py-two);        // coefficients for cubic interpolation along y
    cxyz[ 9] = cp5*(py+one)*(py-one)*(py-two);
    cxyz[10] = cm5*py*(py+one)*(py-two);
    cxyz[11] = cp167*py*(py+one)*(py-one);
#endif
    *ixyz++ = ix;
    // now we compute the coefficients along z using ix
    base  = &(lv->ocz[4*ix]);  // precomputed inverses of denominators
    pos   = &(lv->t2[ix]);
    pz = *PZ;
    zza = pz - pos[-1] ; zzb = pz - pos[0] ; zzc = pz - pos[1] ; zzd = pz - pos[2] ; 
    dz = lv->odz[ix] * zzb;   // (pz - pos[2]) / (pos[2] - pos[1])
    zzab = zza * zzb ; zzcd = zzc * zzd;
    cxyz[16] = zzb * zzcd * base[0];   //   cxyz[16] = TRIPRD(pz,pos[1],pos[2],pos[3]) * base[0];
    cxyz[17] = zza * zzcd * base[1];   //   cxyz[17] = TRIPRD(pz,pos[0],pos[2],pos[3]) * base[1];
    cxyz[18] = zzd * zzab * base[2];   //   cxyz[18] = TRIPRD(pz,pos[0],pos[1],pos[3]) * base[2];
    cxyz[19] = zzc * zzab * base[3];   //   cxyz[19] = TRIPRD(pz,pos[0],pos[1],pos[2]) * base[3];

    cxyz[ 4] = 0.0;                    // coefficients for linear interpolation along x
    cxyz[ 5] = 1.0 - px;
    cxyz[ 6] = px;
    cxyz[ 7] = 0.0;
    cxyz[12] = 0.0;                    // coefficients for linear interpolation along y
    cxyz[13] = 1.0 - py;
    cxyz[14] = py;
    cxyz[15] = 0.0;
    cxyz[20] = 0.0;                    // coefficients for linear interpolation along z
    cxyz[21] = 1.0 - dz;
    cxyz[22] = dz;
    cxyz[23] = 0.0;
    cxyz += 24;
  }
  return ix;
}

// calculate coefficients for x, y, z interpolation from lv table (see Vsearch_setup) and positions
// PX    : array of positions along x (dimension n) (fractional index space, "origin 1")
// PY    : array of positions along y (dimension n) (fractional index space, "origin 1")
// PZ    : array of positions along z (dimension n) (fractional index space, "origin 1")
// cxyz  : coefficients for tricubic or bicubic-linear interpolation (dimension 24*n)
//         cxyz[I][ 0- 7] : along x for point I   Fortran cxyz( 0: 7,I)
//         cxyz[I][ 8-15] : along y for point I   Fortran cxyz( 8:15,I)
//         cxyz[I][16-23] : along z for point I   Fortran cxyz(16:23,I)
//         in each group of 8 coefficients, the first 4 are for the cubic case, the last 4 for the linear case
// ixyz  : array of indices along z
// lv    : x y z axis description  (see Vsearch_setup) (OPAQUE OBJECT)
// n     : number of points to process (PX[n], PY[n], PZ[n], ixys[n]) ( cxyz[n][24] )
//
// Fortran dimensions : PX(n), PY(n), PZ(n), ixys(n), cxyz(24,n)
// function return : index for last point along z
typedef struct{
  float px;    // position along x in index space
  float py;    // position along y in index space
  float pz;    // position along z in index space
  float z;     // absolute position along z 
} pxpypzz;
// int Vcoef_xyz_incr(uint32_t *ixyz, double *cxyz, double *PX, double *PY, double *PZ, lvtab *lv, int n){
int Vcoef_ixyz8_pxyz8(uint32_t *ixyz, double *cxyz, pxpypzz *PXYZ, lvtab *lv, int n){
  int ix, irep, ijk, linear;
  double pxy[2], *base, *pos;
  double zza, zzb, zzc, zzd, zzab, zzcd, dz, px, py, pz;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m256i v0, v1, vc, t, ttop, tbot, ttemp;
  __m128d vxy, vc1, vc5, vc3, vp6, xxmx, xxpx, xxm1, x5p1, x6m3, x5m1, x6m6, dtmp;
  __m128d vr0, vr1, vr2, vr3;
  __m128i itmp;
  int m0, m1;
#else
  int i, j;
#endif
  for(irep=0 ; irep <n ; irep++){                      // loop over points
    px = PXYZ[irep].px ;            // fractional index positions along x, y, z (float to double)
    py = PXYZ[irep].py ;
    pz = PXYZ[irep].pz ;

    ix = PXYZ[irep].px ;
    px = px - ix;                   // px is now deltax (fractional part of px)
    ijk = ix - 2;                   // x displacement (elements), ix assumed to always be >1 and < ni-1

    ix = PXYZ[irep].py ;
    py = py - ix;                   // py is now deltay (fractional part of py)
    ijk = ijk + (ix -2) * lv->ni;   // add y displacement (rows), ix assumed to always be >1 and < nj-1

    ix = PXYZ[irep].pz ; 
    if(ix<1) ix = 1; 
    if(ix>lv->nk-1) ix = lv->nk-1;  // ix < 1 or ix > nk-1 will result in linear extrapolation
    dz = pz - ix;                   // dz is now "fractional" part of pz  (may be <0 or >1 if extrapolating)
    ijk = ijk + (ix -1) * lv->nij;  // add z displacement (2D planes)

    ix--;                           // ix needs to be in "origin 0" (C index from Fortran index)
    linear = (ix - 1) | (lv->nk - 3 - ix); 
    linear >>= 31;                  // nonzero only if ix < 1 or ix > nk -3 (top and bottom intervals)
    if(! linear) ijk = ijk - lv->nij;  // not the linear case, go down one 2D plane to get lower left corner of 4x4x4 cube
    ixyz[irep] = ijk;               // store collapsed displacement
    // now we compute the coefficients along z using ix and dx
    if(linear){
      cxyz[16] = 0.0;                    // coefficients for linear interpolation along z
      cxyz[17] = 1.0 - dz;
      cxyz[18] = dz;
      cxyz[19] = 0.0;
    }else{
      base  = &(lv->ocz[4*ix]);  // precomputed inverses of denominators
      pos   = &(lv->t2[ix]);
//     pz  = dz * pos[1] + (1.0 - dz) * pos[0];
      pz = PXYZ[irep].z  ;       // absolute position along z is used to compute coefficients
      zza = pz - pos[-1] ; zzb = pz - pos[0] ; zzc = pz - pos[1] ; zzd = pz - pos[2] ; 
      zzab = zza * zzb ; zzcd = zzc * zzd;
      // printf("target = %8.5f, dz = %8.5f, levels = %8.5f %8.5f\n",*PZ,dz,pos[0],pos[1]);
      // printf("target = %8.5f, base = %8.5f %8.5f %8.5f %8.5f\n",pz,base[0],base[1],base[2],base[3]);
      // printf("dz     = %8.5f, zz   = %8.5f %8.5f %8.5f %8.5f\n",dz,zza,zzb,zzc,zzd);
      // printf("ix     = %d, lv->odz[ix] = %8.5f\n",ix,lv->odz[ix]);
      cxyz[16] = zzb * zzcd * base[0];   //   cxyz[16] = TRIPRD(pz,pos[1],pos[2],pos[3]) * base[0];
      cxyz[17] = zza * zzcd * base[1];   //   cxyz[17] = TRIPRD(pz,pos[0],pos[2],pos[3]) * base[1];
      cxyz[18] = zzd * zzab * base[2];   //   cxyz[18] = TRIPRD(pz,pos[0],pos[1],pos[3]) * base[2];
      cxyz[19] = zzc * zzab * base[3];   //   cxyz[19] = TRIPRD(pz,pos[0],pos[1],pos[2]) * base[3];
    }

    cxyz[ 4] = 0.0;                    // coefficients for linear interpolation along x
    cxyz[ 5] = 1.0 - px;
    cxyz[ 6] = px;
    cxyz[ 7] = 0.0;
    cxyz[12] = 0.0;                    // coefficients for linear interpolation along y
    cxyz[13] = 1.0 - py;
    cxyz[14] = py;
    cxyz[15] = 0.0;
    cxyz[20] = 0.0;                    // coefficients for linear interpolation along z
    cxyz[21] = 1.0 - dz;
    cxyz[22] = dz;
    cxyz[23] = 0.0;
// printf("pz,dz,ix,nk,linear = %8.5f %8.5f %d %d %d\n",pz,dz,ix,lv->nk,linear);
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
    pxy[0] = px;                           // vector of length 2 to do x and y in one shot
    pxy[1] = py;
    vxy  = _mm_loadu_pd(pxy);              // coefficients for x and y will be interleaved

    vc1  = _mm_set1_pd(one);               //  one    (1.0)
    vc5  = _mm_set1_pd(cp5);               //  cp5    (0.5)
    vc3  = _mm_set1_pd(cp133);             //  cp133  (1/3)
    vp6  = _mm_set1_pd(cp167);             //  cp167  (1/6)

    // alternate formula to compute coefficients taking advantage of independent FMAs
    //   pa = cm167*xy*(xy-one)*(xy-two)     = xy*(xy-one) * (-cp167)*(xy-two)  =  (xy*xy  - xy)   * (-(cp167*xy) + cp133)
    //   pb = cp5*(xy+one)*(xy-one)*(xy-two) = cp5*(xy-two) * (xy+one)*(xy-one) =  (cp5*xy - one)  * (xy*xy       - one)
    //   pc = cm5*xy*(xy+one)*(xy-two)       = xy*(xy+one) * (-cp5)*(xy-two)    =  (xy*xy  + xy)   * (-(cp5*xy)   + one)
    //   pd = cp167*xy*(xy+one)*(xy-one)     = xy*(xy+one) * cp167*(xy-one)     =  (xy*xy  + xy)   * (cp167*xy    - cp167)
    // STEP1 : 7 independent terms using 3 different FMAs  a*b+c, a*b-c, (-a*b)+c
    xxmx = _mm_fmsub_pd(vxy,vxy,vxy);      //  xy*xy       - xy
    x6m3 = _mm_fnmadd_pd(vxy,vp6,vc3);     //  -(cp167*xy) + cp133
    x5p1 = _mm_fmsub_pd(vc5,vxy,vc1);      //  cp5*xy      - one
    xxm1 = _mm_fmsub_pd(vxy,vxy,vc1);      //  xy*xy       - one
    xxpx = _mm_fmadd_pd(vxy,vxy,vxy);      //  xy*xy       + xy
    x5m1 = _mm_fnmadd_pd(vxy,vc5,vc1);     //  -(cp5*xy)   + one
    x6m6 = _mm_fmsub_pd(vxy,vp6,vp6);      //  cp167*xy    - cp167

    // STEP2 : multiply STEP 1 terms (4 independent operations)
    vr0  = _mm_mul_pd(xxmx,x6m3);          // coefficients for x and y are interleaved
    vr1  = _mm_mul_pd(x5p1,xxm1);
    vr2  = _mm_mul_pd(xxpx,x5m1);
    vr3  = _mm_mul_pd(xxpx,x6m6);

    // final unshuffle to separate even terms (px) and odd terms (py) before storing them (independent operations)
    _mm_storeu_pd(cxyz   ,_mm_unpacklo_pd(vr0,vr1));  // cxyz[ 0: 1] = cx[0], cx[1]
    _mm_storeu_pd(cxyz+ 2,_mm_unpacklo_pd(vr2,vr3));  // cxyz[ 2: 3] = cx[2], cx[3]
    _mm_storeu_pd(cxyz+ 8,_mm_unpackhi_pd(vr0,vr1));  // cxyz[ 8: 9] = cy[0], cy[1]
    _mm_storeu_pd(cxyz+10,_mm_unpackhi_pd(vr2,vr3));  // cxyz[10:11] = cy[2], cy[3]

#else
//     cxyz[ 0] = (px*px  - px)   * (-(cp167*px) + cp133);
//     cxyz[ 1] = (cp5*px - one)  * (px*px       - one);
//     cxyz[ 2] = (px*px  + px)   * (-(cp5*px)   + one);
//     cxyz[ 3] = (px*px  + px)   * (cp167*px    - cp167);
// 
//     cxyz[ 8] = (py*py  - py)   * (-(cp167*py) + cp133);
//     cxyz[ 9] = (cp5*py - one)  * (py*py       - one);
//     cxyz[10] = (py*py  + py)   * (-(cp5*py)   + one);
//     cxyz[11] = (py*py  + py)   * (cp167*py    - cp167);

    cxyz[ 0] = cm167*px*(px-one)*(px-two);        // coefficients for cubic interpolation along x
    cxyz[ 1] = cp5*(px+one)*(px-one)*(px-two);
    cxyz[ 2] = cm5*px*(px+one)*(px-two);
    cxyz[ 3] = cp167*px*(px+one)*(px-one);

    cxyz[ 8] = cm167*py*(py-one)*(py-two);        // coefficients for cubic interpolation along y
    cxyz[ 9] = cp5*(py+one)*(py-one)*(py-two);
    cxyz[10] = cm5*py*(py+one)*(py-two);
    cxyz[11] = cp167*py*(py+one)*(py-one);
#endif
    cxyz += 24;
  }
  return ix;  // vertical index for last point
}

// version for monotonically decreasing positive table values
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
static inline int Vsearch_list_dec_inline(double target, lvtab *lv){
  int ix, j;
#if defined(__AVX2__) && defined(__x86_64__) && defined(SIMD)
  __m256i v0, v1, vc, t;
  int m0, m1, pop0, pop1;

  t    = (__m256i) _mm256_broadcast_sd(&target);
  vc   = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  if(target > lv->t0[0])  t = (__m256i) _mm256_broadcast_sd(&lv->t0[0]);   // target > first element in table
  if(target < lv->top  )  t = (__m256i) _mm256_broadcast_sd(&lv->top  );   // target < next to last element in table

  t    = _mm256_sub_epi64((__m256i) t, vc);   // we want target > tbl[i] so we subtract 1 from target

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[0]));                  // 8 values to scan from table t0
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) v0);  // t - tbl >= 0 if (target+1) > tbl
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) v1);
  m0   = _mm256_movemask_pd((__m256d) v0);              // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);                            // count bits on in mask
  pop1 = _mm_popcnt_u32(m1);
  j    = pop0+pop1-1;
  ix   = j << 3;                                        // index into t1 table for 8 values scan

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix  ]));              // 8 entries to scan from table t1
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) v0);
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) v1);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = _mm_popcnt_u32(m1);
  j    = pop0+pop1-1;
  ix   = (ix + j) <<2 ;                                 // index into t2 table for final 4 value scan

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t2[ix  ]));              // only 4 entries to scan from table t2
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) v0);
  m0   = _mm256_movemask_pd((__m256d) v0);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = 0;
  j = pop0+pop1-1;
  ix = (ix + j) ;
#else
  int i;
  d_l_p dlt, dlr, dlm;

  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l > dlr.l) target = dlr.d;  // target > first element in table
  if(dlt.l < dlm.l) target = dlm.d;  // target < next to last element in table
//   if(target > lv->t0[0]) target = lv->t0[0];  // target > first element in table
//   if(target < lv->top  ) target = lv->top;  // target < next to last element in table

  j = 7;
  for(i=7 ; i>0 ; i--) { if(target > lv->t0[i]) j = i - 1; }
  ix = j << 3;
  j = 7;
  for(i=7 ; i>0 ; i--) { if(target > lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
  j = 3;
  for(i=3 ; i>0 ; i--) { if(target > lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
#endif
  return ix;
}
int Vsearch_list_dec(double target, lvtab *lv){
  return Vsearch_list_dec_inline(target, lv);
}
void Vsearch_list_dec_n(int *ix, double *target, lvtab *lv, int n){
  int i;
  for(i=0 ; i<n ; i++){
    ix[i] = Vsearch_list_dec_inline(target[i], lv);
  }
}

#if defined(SELF_TEST)
// straight C version for comparison purposes
int search_list_dec(double target, lvtab *lv){
  int ix;
  int i, j;
  d_l_p dlt, dlr, dlm;

  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l > dlr.l) target = dlr.d;  // target > first element in table
  if(dlt.l < dlm.l) target = dlm.d;  // target < next to last element in table
  j = 7;

  for(i=7 ; i>0 ; i--) { if(target > lv->t0[i]) j = i - 1; }
  ix = j << 3;

  j = 7;
  for(i=7 ; i>0 ; i--) { if(target > lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;

  j = 3;
  for(i=3 ; i>0 ; i--) { if(target > lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;

  return ix;
}

int search_list_inc(double target, lvtab *lv){
  int ix;
  int i, j;

  d_l_p dlt, dlr, dlm;
  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l < dlr.l) target = dlr.d;  // target > first element in table
  if(dlt.l > dlm.l) target = dlm.d;  // target < next to last element in table

  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t0[i]) j = i - 1; }
  ix = j << 3;

  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;

  j = 3;
  for(i=3 ; i>0 ; i--) { if(target < lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
  return ix;
}
#endif
