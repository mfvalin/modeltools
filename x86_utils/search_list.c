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
  double t0[ 8];
  double t1[64];   // will support at most 256 levels
  double t2[257];
  double top;      // t2[nk-2]
  double x0;       // x origin
  double odx;      // inverse of deltax
  double y0;       // y origin
  double ody;      // inverse of deltay
  double *odz;     // pointer to inverse of deltaz array [nk doubles]
  double *ocz;     // pointer to inverse of coefficient denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of planes (max 255)
  uint32_t nij;    // ni*nj, distance between planes
}lvtab;

#if defined(__AVX2__) && defined(__x86_64__)
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

// return pointer to filled table
// targets are expected to be positive, and monotonically decreasing
lvtab * Vsearch_setup(double *targets, int n){
  int i;
  lvtab *lv;
  double pad = DBL_MAX;
printf("Vsearch_setup: n = %d, %10.3f %10.3f %10.3f %10.3f \n",n,targets[0],targets[1],targets[n-2],targets[n-1]);
  if(n > 256) return NULL;
  if(0 != posix_memalign( (void **) &lv, 32, sizeof(lvtab) ) ) return NULL ;

  if(targets[1] < targets[0]) pad = -pad;    // decreasing values
  for(i=0 ; i<n   ; i++) { lv->t2[i] = targets[i] ; }   // z coordinate table
  for(i=n ; i<256 ; i++) { lv->t2[i] = pad ; }          // pad at end of table

  for(i=0  ; i<64 ; i++) { lv->t1[i] = lv->t2[(i<<2)] ; } ;   // lookup for t2 (every fourth entry)
  for(i=0  ; i<8  ; i++) { lv->t0[i] = lv->t1[(i<<3)] ; } ;   // lookup for t1 (every eighth entry)
  // inverse of intervals ( 0 to n-2 only make sense )
  lv->odz = malloc(n * sizeof(double));
  if(NULL == lv->odz){
    free(lv);
    return NULL;
  }
  for(i=0 ; i<n-1   ; i++) { lv->odz[i] = (targets[i+1] - targets[i]) ; }
  lv->odz[n - 1] = 0.0;   // NEVER USED
  // denominators for Lagrange cubic polynomials coefficients ( 1 to n-2 only make sense )
  lv->ocz = malloc(4 * n * sizeof(double));
  if(NULL == lv->ocz){
    free(lv->odz);
    free(lv);
    return NULL;
  }
  denominators( &(lv->ocz[0]) , targets[0], targets[1], targets[2], targets[3]);   // level 0 coeffs are normally not used
  for(i=1 ; i<n-1   ; i++) {
    denominators( &(lv->ocz[4*i]) , targets[i-1], targets[i  ], targets[i+1], targets[i+2]);
  }
  denominators( &(lv->ocz[n-1]) , targets[n-4], targets[n-3], targets[n-2], targets[n-1]);   // level n-1 coeffs are normally not used

  lv->top = lv->t2[n-2]; // next to last value along z
  lv->x0 = 0.0;          // base coordinate along x
  lv->odx = 0.0;         // deltax
  lv->y0 = 0.0;          // base coordinate along y
  lv->ody = 0.0;         // deltay
  lv->ni = 0;            // nb of points along x
  lv->nj = 0;            // nb of points along y
  lv->nk = n;            // nb of points along z
  lv->nij = 0;           // ni * nj
  return lv;
}

// NOTE: 
//      min(x,y) :  y + ( (x - y) & ( (x - y) >> 63 ) )
//      max(x,y) :  x - ( (x - y) & ( (x - y) >> 63 ) )
//      _mm256_blendv_pd (__m256d a, __m256d b, __m256d mask)  a if masksign == 0 , b if == 1
//      max(x,y) = blendv(x , y, x-y)
//      min(x,y) = blendv(y , x, x-y)
// version for monotonically increasing positive table values
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
int Vsearch_list_inc(double target, lvtab *lv){
  int ix;
#if defined(__AVX2__) && defined(__x86_64__)
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
  return ix;
}

double Vsearch_list_inc_d(double target, lvtab *lv){
  int ix;
#if defined(__AVX2__) && defined(__x86_64__)
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
  return ix + ( (target - lv->t2[ix]) * lv->odz[ix] );
}

// px    : position along x
// py    : position along y
// pz    : position along z
// cxyz  : coefficients for tricubic or bicubic-linear interpolation
//         cxyz[ 0- 7] : along x
//         cxyz[ 8-15] : along y
//         cxyz[16-23] : along z
// lv    : x y z axis description
int Vcoef_xyz_inc(double *cxyz, double px, double py, double pz, lvtab *lv){
  int ix;
  double pxy[2], *base, *denom, *pos;
#if defined(__AVX2__) && defined(__x86_64__)
  __m256i v0, v1, vc, t, ttop, tbot, ttemp;
  __m128d vxy, vc1, vc5, vc3, vp6, xxmx, xxpx, xxm1, x5p1, x6m3, x5m1, x6m6, dtmp;
  __m128d vr0, vr1, vr2, vr3;
  __m128i itmp;
  int m0, m1;

  pxy[0] = (px - lv->x0) * lv->odx;      // position to fractional index
  pxy[1] = (py - lv->y0) * lv->ody;      // position to fractional index
  vxy  = _mm_loadu_pd(pxy);              // coefficients for x and y will be interleaved
  dtmp = _mm_floor_pd(vxy);              // convert positions x y into interval dx dy
  vxy  = _mm_sub_pd(vxy,dtmp);           // subtract integer part

  vc1  = _mm_set1_pd(one);               //  one    (1.0)
  vc5  = _mm_set1_pd(cp5);               //  cp5    (0.5)
  vc3  = _mm_set1_pd(cp133);             //  cp133  (1/3)
  vp6  = _mm_set1_pd(cp167);             //  cp167  (1/6)

  t     = (__m256i) _mm256_broadcast_sd(&pz);
  vc    = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  tbot  = (__m256i) _mm256_broadcast_sd(&(lv->t0[0]));    // first element in table
  ttop  = (__m256i) _mm256_broadcast_sd(&(lv->top  ));    // next to last element in table
  if(pz < lv->t0[0]) t = tbot ;    // pz < first element in table, set to first element in table
  if(pz > lv->top  ) t = ttop;     // pz > next to last element in table et to next to last element
//   ttemp = (__m256i) _mm256_sub_epi64(t,tbot);
//   t     = (__m256i) _mm256_blendv_pd((__m256d) t, (__m256d) tbot, (__m256d) ttemp);
//   ttemp = (__m256i) _mm256_sub_epi64(ttop,t);
//   t     = (__m256i) _mm256_blendv_pd((__m256d) t, (__m256d) ttop, (__m256d) ttemp);
//   if(dlt.l < dlr.l) t =(__m256i)  _mm256_broadcast_sd(&dlr.d);    // pz < first element in table
//   if(dlt.l > dlm.l) t =(__m256i)  _mm256_broadcast_sd(&dlm.d);    // pz > next to last element in table
  
  t    = _mm256_add_epi64((__m256i) t, vc);   // we want pz < tbl[i] so we add 1 to pz

  v0   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[0]));                  // 8 values to scan from table t0
  v1   = _mm256_lddqu_si256((__m256i const *) &(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) v0,(__m256i) t);  // t - tbl >= 0 if (pz-1) < tbl
  v1   = _mm256_sub_epi64((__m256i) v1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);              // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  ix   = (_mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1) << 3;         // index into t1 table for 8 values scan

  // alternate formula to compute coefficients to take advantage of FMAs
  //   p = cm167*xy*(xy-one)*(xy-two)     = xy*(xy-one) * (-cp167)*(xy-two)  =  (xy*xy  - xy)   * (-(cp167*xy) + cp133)
  //   p = cp5*(xy+one)*(xy-one)*(xy-two) = cp5*(xy-two) * (xy+one)*(xy-one) =  (cp5*xy - one)  * (xy*xy       - one)
  //   p = cm5*xy*(xy+one)*(xy-two)       = xy*(xy+one) * (-cp5)*(xy-two)    =  (xy*xy  + xy)   * (-(cp5*xy)   + one)
  //   p = cp167*xy*(xy+one)*(xy-one)     = xy*(xy+one) * cp167*(xy-one)     =  (xy*xy  + xy)   * (cp167*xy    - cp167)
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
  int i, j;

  if(pz < lv->t0[0]) pz = lv->t0[0];  // pz < first element in table
  if(pz > lv->top  ) pz = lv->top  ;    // pz > next to last element in table

  j = 7;
  for(i=7 ; i>0 ; i--) { if(pz < lv->t0[i]) j = i - 1; }
  ix = j << 3;
  j = 7;
  for(i=7 ; i>0 ; i--) { if(pz < lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
  j = 3;
  for(i=3 ; i>0 ; i--) { if(pz < lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
  pxy[ 0] = px;
  pxy[ 1] = py;

  cxyz[ 0] = cm167*pxy[0]*(pxy[0]-one)*(pxy[0]-two);        // coefficients for cubic interpolation along x
  cxyz[ 1] = cp5*(pxy[0]+one)*(pxy[0]-one)*(pxy[0]-two);
  cxyz[ 2] = cm5*pxy[0]*(pxy[0]+one)*(pxy[0]-two);
  cxyz[ 3] = cp167*pxy[0]*(pxy[0]+one)*(pxy[0]-one);

  cxyz[ 8] = cm167*pxy[1]*(pxy[1]-one)*(pxy[1]-two);        // coefficients for cubic interpolation along y
  cxyz[ 9] = cp5*(pxy[1]+one)*(pxy[1]-one)*(pxy[1]-two);
  cxyz[10] = cm5*pxy[1]*(pxy[1]+one)*(pxy[1]-two);
  cxyz[11] = cp167*pxy[1]*(pxy[1]+one)*(pxy[1]-one);
#endif
  // now we can compute the coefficients along z
  base  = &(lv->ocz[4*ix]);
  denom = &(lv->odz[ix]);
  pos   = &(lv->t2[ix]);
  cxyz[16] = TRIPRD(pz,pos[1],pos[2],pos[3]) * denom[0];
  cxyz[17] = TRIPRD(pz,pos[0],pos[2],pos[3]) * denom[1];
  cxyz[18] = TRIPRD(pz,pos[0],pos[1],pos[3]) * denom[2];
  cxyz[19] = TRIPRD(pz,pos[0],pos[1],pos[2]) * denom[3];
  return ix;
}

// version for monotonically decreasing positive table values
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
int Vsearch_list_dec(double target, lvtab *lv){
  int ix, j;
#if defined(__AVX2__) && defined(__x86_64__)
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

  if(target > lv->t0[0]) target = lv->t0[0];  // target > first element in table
  if(target < lv->top  ) target = lv->top;  // target < next to last element in table

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
//   printf("target = %f\n",target);
  for(i=7 ; i>0 ; i--) { if(target > lv->t0[i]) j = i - 1; }
  ix = j << 3;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t0[j-1],lv->t0[j],lv->t1[ix],lv->t1[ix+7]);
  j = 7;
  for(i=7 ; i>0 ; i--) { if(target > lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j],lv->t2[ix],lv->t2[ix+7]);
  j = 3;
  for(i=3 ; i>0 ; i--) { if(target > lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
//   printf("j = %d, ix = %d, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j]);
  return ix;
}

int search_list_inc(double target, lvtab *lv){
  int ix;
  int i, j;
//   d_l_p dlt, dlr, dlm;

//   dlt.d = target;
//   dlr.d = lv->t0[0];
//   dlm.d = lv->top;
  if(target < lv->t0[0]) target = lv->t0[0];  // target < first element in table
  if(target > lv->top  ) target = lv->top;    // target > next to last element in table

  j = 7;
//   printf("target = %f\n",target);
  for(i=7 ; i>0 ; i--) { if(target < lv->t0[i]) j = i - 1; }
  ix = j << 3;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t0[j-1],lv->t0[j],lv->t1[ix],lv->t1[ix+7]);
  j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j],lv->t2[ix],lv->t2[ix+7]);
  j = 3;
  for(i=3 ; i>0 ; i--) { if(target < lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
//   printf("j = %d, ix = %d, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j]);
  return ix;
}
#endif
