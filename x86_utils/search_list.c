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

typedef union{
  double  d;
  int64_t l;
  double *p;
}d_l_p;

typedef struct{
  double x0;       // x origin
  double ovdx;     // inverse of deltax
  double y0;       // y origin
  double ovdy;     // inverse of deltay
  d_l_p  odz;      // pointer to inverse of deltaz array [nk doubles]
  d_l_p  ocz;      // pointer to inverse of denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of planes (max 255)
  uint32_t nij;    // ni*nj, distance between planes
  double t0[8];
  double t1[64];   // will support at most 255 levels
  double t2[256];
  double top;      // t2[nk-2]
}lvtab;

// typedef struct{
//   uint64_t dummy[8];
//   uint64_t t0[8];
//   uint64_t t1[64];
//   uint64_t t2[256];
// }lvtabi;

static int64_t ONE = 1;

// return pointer to filled table
// targets are expected to be positive, and monotonically decreasing
lvtab * Vsearch_setup(double *targets, int n){
  int i;
  lvtab *lv;
printf("Vsearch_setup: n = %d, %10.3f %10.3f %10.3f %10.3f \n",n,targets[0],targets[1],targets[n-2],targets[n-1]);
  if(n > 256) return NULL;
  if(0 != posix_memalign( (void **) &lv, 32, sizeof(lvtab) ) ) return NULL ;

//   for(i=0 ; i<n   ; i++) { lv->t2[i] = targets[i] ; }
  memcpy(lv->t2 , targets, sizeof(double)*n);
  for( i=n ; i<256 ; i++) { lv->t2[i] = lv->t2[n-1] ; }
//   if(n < 256) memset(&lv->t2[n],0,sizeof(double)*(256-n));

  for(i=0 ; i<64 ; i++) { lv->t1[i] = lv->t2[(i<<2)] ; } ;
  for(i=0 ; i<8  ; i++) { lv->t0[i] = lv->t1[(i<<3)] ; } ;
  lv->odz.p = NULL;
  lv->ocz.p = NULL;
  lv->top = lv->t2[n-2]; // (lv->t2[n-1] + lv->t2[n-2]) * .5;
  return lv;
}

// version for monotonically increasing positive table values
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
int Vsearch_list_inc(double target, lvtab *lv){
  d_l_p dlt, dlr, dlm;
  int ix, i, j;
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d t;
  __m256d tbl0, tbl1;
  __m256i v0, v1, vc;
  int m0, m1, pop0, pop1;

  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  t     = _mm256_broadcast_sd(&target);
  if(dlt.l < dlr.l) dlt.l = dlr.l;    // target < first element in table
  if(dlt.l > dlm.l) dlt.l = dlm.l;    // target > next to last element in table
  t = _mm256_broadcast_sd(&dlt.d);
//   printf("target = %f\n",dlt.d);

  vc   = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  t    = (__m256d) _mm256_add_epi64((__m256i) t, vc);   // we want target < tbl[i] so we add 1 to target

  tbl0 = _mm256_load_pd(&(lv->t0[0]));                  // 8 values to scan from table t0
  tbl1 = _mm256_load_pd(&(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) tbl0,(__m256i) t);  // t - tbl >= 0 if (target-1) < tbl
  v1   = _mm256_sub_epi64((__m256i) tbl1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);              // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  j    = _mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1;   // count bits on in mask
  ix   = j << 3;                                        // index into t1 table for 8 values scan

  tbl0 = _mm256_load_pd(&(lv->t1[ix  ]));               // 8 values to scan from table t1
  tbl1 = _mm256_load_pd(&(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) tbl0,(__m256i) t);
  v1   = _mm256_sub_epi64((__m256i) tbl1,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  j    = _mm_popcnt_u32(m0) + _mm_popcnt_u32(m1) - 1;   // count bits on in mask
  ix   = (ix + j) <<2 ;                                 // index into t2 table for final 4 value scan

  tbl0 = _mm256_load_pd(&(lv->t2[ix  ]));               // only 4 values to scan from table t2
  v0   = _mm256_sub_epi64((__m256i) tbl0,(__m256i) t);
  m0   = _mm256_movemask_pd((__m256d) v0);
  ix = ix + _mm_popcnt_u32(m0) - 1 ;
#else
  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l < dlr.l) target = dlr.d;  // target < first element in table
  if(dlt.l > dlm.l) target = dlm.d;    // target > next to last element in table

//   j = 7;
//   printf("target = %f\n",target);
  for(i=7 ; i>0 ; i--) { if(target < lv->t0[i]) j = i - 1; }
  ix = j << 3;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t0[j-1],lv->t0[j],lv->t1[ix],lv->t1[ix+7]);
//   j = 7;
  for(i=7 ; i>0 ; i--) { if(target < lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j],lv->t2[ix],lv->t2[ix+7]);
//   j = 3;
  for(i=3 ; i>0 ; i--) { if(target < lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
//   printf("j = %d, ix = %d, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j]);
#endif
  return ix;
}

// version for monotonically decreasing positive table values
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
int Vsearch_list_dec(double target, lvtab *lv){
  d_l_p dlt, dlr, dlm;
  int ix, i, j;
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d t;
  __m256d tbl0, tbl1;
  __m256i v0, v1, vc;
  int m0, m1, pop0, pop1;

  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  t    = _mm256_broadcast_sd(&target);
  if(dlt.l > dlr.l) dlt.l = dlr.l;    // target > first element in table
  if(dlt.l < dlm.l) dlt.l = dlm.l;    // target < next to last element in table
  t = _mm256_broadcast_sd(&dlt.d);
// printf("target = %f\n",dlt.d);

  vc   = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  t    = (__m256d) _mm256_sub_epi64((__m256i) t, vc);   // we want target > tbl[i] so we subtract 1 from target

  tbl0 = _mm256_load_pd(&(lv->t0[0]));                  // 8 values to scan from table t0
  tbl1 = _mm256_load_pd(&(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);  // t - tbl >= 0 if (target+1) > tbl
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);              // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);                            // count bits on in mask
  pop1 = _mm_popcnt_u32(m1);
  j    = pop0+pop1-1;
  ix   = j << 3;                                        // index into t1 table for 8 values scan
// printf("ix = %d %f %f %f %f\n",ix,lv->t1[ix],lv->t1[ix+1],lv->t1[ix+6],lv->t1[ix+7]);
  tbl0 = _mm256_load_pd(&(lv->t1[ix  ]));              // 8 entries to scan from table t1
  tbl1 = _mm256_load_pd(&(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = _mm_popcnt_u32(m1);
  j    = pop0+pop1-1;
// printf("j = %d, %f %f %f \n",j,lv->t1[ix+j],lv->t1[ix+j+1],lv->t1[ix+j+2]);
  ix   = (ix + j) <<2 ;                                 // index into t2 table for final 4 value scan
// printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t2[ix],lv->t2[ix+1],lv->t2[ix+2],lv->t2[ix+3]);
  tbl0 = _mm256_load_pd(&(lv->t2[ix  ]));              // only 4 entries to scan from table t2
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);
  m0   = _mm256_movemask_pd((__m256d) v0);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = 0;
  j = pop0+pop1-1;
  ix = (ix + j) ;
// printf("j = %d, ix = %d\n",j,ix);
#else
  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l > dlr.l) target = dlr.d;  // target > first element in table
  if(dlt.l < dlm.l) target = dlm.d;  // target < next to last element in table
//   j = 7;
//   printf("target = %f\n",target);
  for(i=7 ; i>0 ; i--) { if(target > lv->t0[i]) j = i - 1; }
  ix = j << 3;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t0[j-1],lv->t0[j],lv->t1[ix],lv->t1[ix+7]);
//   j = 7;
  for(i=7 ; i>0 ; i--) { if(target > lv->t1[ix+i]) j = i - 1; }
  ix = (ix + j) <<2 ;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j],lv->t2[ix],lv->t2[ix+7]);
//   j = 3;
  for(i=3 ; i>0 ; i--) { if(target > lv->t2[ix+i]) j = i - 1; }
  ix = ix + j;
//   printf("j = %d, ix = %d, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j]);
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
  d_l_p dlt, dlr, dlm;

  dlt.d = target;
  dlr.d = lv->t0[0];
  dlm.d = lv->top;
  if(dlt.l < dlr.l) target = dlr.d;  // target < first element in table
  if(dlt.l > dlm.l) target = dlm.d;    // target > next to last element in table

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

uint64_t rdtscp_(void) {   // version "in order" avec "serialization"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#include <sys/time.h>
#define NT 251
int main(){
  int i, status, ix;
  double T;
  struct timeval tv0, tv1;
  int64_t t0, t1;
  int32_t count, cnt;
  lvtab *lv, *lv2;
  double levels[256];
  double levels2[256];

  for(i=0 ; i<NT ; i++) { levels[i] = NT - i ; levels2[i] = i + 1 ;}
  lv  = Vsearch_setup(levels, NT);
  lv2 = Vsearch_setup(levels2,NT);

//   printf("table0 =\n"); for(i=0 ; i<8 ; i++) printf("%10.3f",lv->t0[i]); printf("\n");
//   for(T=255.0 ; T>253.9 ; T=T-.5){
//     ix = search_list_dec(T,lv);
//     printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//     ix = Vsearch_list_dec(T,lv);
//     printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n\n\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//   }
//   exit(0);
//   T = 1.5;
//   ix = search_list_dec(T,lv);
//   printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//   T = 1.0;
//   ix = search_list_dec(T,lv);
//   printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);

  T = NT + .5 ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);

  T = .5 ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);

  count = (NT - 1) * 1000000;
  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(Vsearch_list_dec(T,lv) == -1) exit(1);
//     count++;
    ix = Vsearch_list_dec(T,lv); ; // cnt += ix;
    if(T > lv->t2[ix] || T < lv->t2[ix+1]){
      printf("ERROR(Vdecr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
      printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
      exit(1);
    }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDd  : %ld microseconds for %d iterations, %f us/iter (%d)\n",t1-t0,count,T,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(Vsearch_list_inc(T,lv2) == -1) exit(1);
    ix = Vsearch_list_inc(T,lv2); // cnt += ix;
    if(T < lv2->t2[ix] || T > lv2->t2[ix+1]){
      printf("ERROR(Vincr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
      printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv2->t2[ix], T-lv2->t2[ix+1]);
      exit(1);
    }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDi  : %ld microseconds for %d iterations, %f us/iter (%d)\n",t1-t0,count,T,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(search_list_dec(T,lv) == -1) exit(1);
    ix = search_list_dec(T,lv); // cnt += ix;
    if(T > lv->t2[ix] || T < lv->t2[ix+1]){
      printf("ERROR(Vdecr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
      printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
      exit(1);
    }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SCALARd : %ld microseconds for %d iterations, %f us/iter\n",t1-t0,count,T);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(search_list_inc(T,lv2) == -1) exit(1);
    ix = search_list_inc(T,lv2); // cnt += ix;
    if(T < lv2->t2[ix] || T > lv2->t2[ix+1]){
      printf("ERROR(incr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
      printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv2->t2[ix], T-lv2->t2[ix+1]);
      exit(1);
    }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SCALARi : %ld microseconds for %d iterations, %f us/iter\n",t1-t0,count,T);
}
#endif
