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
  double t2[257];  // table for levels (max 256 usable)
  double top;      // t2[nk-2]
  double x0;       // x origin
  double odx;      // inverse of deltax
  double y0;       // y origin
  double ody;      // inverse of deltay
  double *odz;     // pointer to inverse of deltaz array [nk doubles]
  double *ocz;     // pointer to inverse of coefficient denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of levels (max 255)
  uint32_t nij;    // ni*nj, distance between levels
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

// allocate lookup table and
// return pointer to filled table
// targets are expected to be positive, and monotonically increasing or decreasing
lvtab * Vsearch_setup(double *targets, int n){
  int i;
  lvtab *lv;
  double pad;

//   printf("Vsearch_setup: n = %d, %10.3f %10.3f %10.3f %10.3f \n",n,targets[0],targets[1],targets[n-2],targets[n-1]);
  if(n > 256) return NULL;
  if(0 != posix_memalign( (void **) &lv, 32, sizeof(lvtab) ) ) return NULL ;

  pad = DBL_MAX;                                        // increasing values, pad with large positive number
  if(targets[1] < targets[0]) pad = -pad;               // decreasing values, pad with large negative number
  for(i=0 ; i<n   ; i++) { lv->t2[i] = targets[i] ; }   // z coordinate table
  for(i=n ; i<256 ; i++) { lv->t2[i] = pad ; }          // pad at end of table

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
  lv->ocz = malloc(4 * n * sizeof(double));
  if(NULL == lv->ocz){    // malloc failed
    free(lv->odz);        // deallocate odz
    free(lv);             // deallocate lv
    return NULL;
  }
  denominators( &(lv->ocz[0]) , targets[0], targets[1], targets[2], targets[3]);                 // level 0 coeffs are normally not used
  for(i=1 ; i<n-1   ; i++) {
    denominators( &(lv->ocz[4*i]) , targets[i-1], targets[i  ], targets[i+1], targets[i+2]);
  }
  denominators( &(lv->ocz[4*(n-1)]) , targets[n-4], targets[n-3], targets[n-2], targets[n-1]);   // level n-1 coeffs are normally not used

  lv->top = lv->t2[n-2]; // next to last value along z
  lv->x0 = 0.0;          // origin of x coordinates
  lv->odx = 0.0;         // deltax
  lv->y0 = 0.0;          // origin of y coordinates
  lv->ody = 0.0;         // deltay
  lv->ni = 0;            // nb of points along x
  lv->nj = 0;            // nb of points along y
  lv->nk = n;            // nb of points along z
  lv->nij = 0;           // ni * nj
  return lv;             // return pointer to filled table
}

// free lookup table created by Vsearch_setup
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

// used to determine calling overhead
double Vsearch_list_inc_do(double target, lvtab *lv){
  return target;
}

// same as Vsearch_list_inc but returns the position as a fractional index
// may be < 0 or > nk-1 if target is beyond table lv-t2 extreme values
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
  return ix + ( (target - lv->t2[ix]) * lv->odz[ix] );  // (target - lowerlevel) / deltaz
}

// used to measure call sequence overhead
int Vcoef_xyz_inc_o(uint32_t *ixyz, double *cxyz, double *px, double *py, double *pz, lvtab *lv, int n){
  return 0;
}

// calculate coefficients for x, y, z interpolation from lv table (see Vsearch_setup) and positions
// PX    : array of positions along x (dimension n)
// PY    : array of positions along y (dimension n)
// PZ    : array of positions along z (dimension n)
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
#if defined(__AVX2__) && defined(__x86_64__)
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

#if defined(__AVX2__) && defined(__x86_64__)
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
  zza = pz - pos[0] ; zzb = pz - pos[1] ; zzc = pz - pos[2] ; zzd = pz - pos[3] ; 
  dz = lv->odz[ix] * zzc;   // (pz - pos[2]) / (pos[2] - pos[1])
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
