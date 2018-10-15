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

typedef struct{
  double t0[8];
  double t1[64];
  double t2[256];
}lvtab;

typedef struct{
  uint64_t t0[8];
  uint64_t t1[64];
  uint64_t t2[256];
}lvtabi;

typedef union{
  double  d;
  int64_t l;
}d_or_l;

static int64_t ONE = 1;

// return pointer to filled table
// targets are expected to be positive, and monotonically decreasing
lvtab * Vsearch_setup(double *targets, int n){
  int i;
  lvtab *lv;

  if(n > 256) return NULL;
  if(0 != posix_memalign( (void **) &lv, 64, sizeof(lvtab) ) ) return NULL ;

//   for(i=0 ; i<n   ; i++) { lv->t2[i] = targets[i] ; }
  memcpy(lv->t2 , targets, sizeof(double)*n);
//   for(    ; i<256 ; i++) { lv->t2[i] = 0.0 ; 
  if(n < 256) memset(&lv->t2[n],0,sizeof(double)*(256-n));

  for(i=0 ; i<64 ; i++) { lv->t1[i] = lv->t2[(i<<2)] ; } ;
  for(i=0 ; i<8  ; i++) { lv->t0[i] = lv->t2[(i<<5)] ; } ;
  return lv;
}

// SIMD version
// target and tables are of type double but processed as if they were 64 bit positive integers
// 8 8 4 scan pattern
int Vsearch_list(double target, lvtab *lv){
  d_or_l dlt, dlr;
  int ix, i, j;
#if defined(__AVX2__) && defined(__x86_64__)
  __m256d t;
  __m256d tbl0, tbl1;
  __m256i v0, v1, vc;
  int m0, m1, pop0, pop1;

//   printf("target = %f\n",target);
  dlt.d = target;
  dlr.d = lv->t0[0];
  if(dlt.l > dlr.l) return 0;  // target > first element in table

  t    = _mm256_broadcast_sd(&target);
  vc   = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  t    = (__m256d) _mm256_sub_epi64((__m256i) t, vc);   // we want target > tbl[i] so we subtract 1 from target

  tbl0 = _mm256_load_pd(&(lv->t0[0]));                 // 8 values to scan from table t0
  tbl1 = _mm256_load_pd(&(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);  // t - tbl >= 0 if (target+1) > tbl
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);              // transform subtract result signs into mask
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);                            // count bits on in mask
  pop1 = _mm_popcnt_u32(m1);
  j    = pop0+pop1-1;
  ix   = j << 3;                                        // index into t1 table for 8 values scan

  tbl0 = _mm256_load_pd(&(lv->t1[ix  ]));              // 8 entries to scan from table t1
  tbl1 = _mm256_load_pd(&(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = _mm_popcnt_u32(m1);
  j    = pop0+pop1-1;
  ix   = (ix + j) <<2 ;                                 // index into t2 table for final 4 value scan

  tbl0 = _mm256_load_pd(&(lv->t2[ix  ]));              // only 4 entries to scan from table t2
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);
  m0   = _mm256_movemask_pd((__m256d) v0);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = 0;
  j = pop0+pop1-1;
  ix = (ix + j) ;
#else
  dlt.d = target;
  dlr.d = lv->t0[0];
  if(dlt.l > dlr.l) return 0;  // target > first element in table

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
#endif
  return ix;
}

#if defined(SELF_TEST)
// straight C version for comparison purposes
int search_list(double target, lvtab *lv){
  int ix;
  int i, j;
  d_or_l dlt, dlr;

  dlt.d = target;
  dlr.d = lv->t0[0];
  if(dlt.l > dlr.l) return 0;  // target > first element in table

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
int main(){
  int i, status, ix;
  double T;
  struct timeval tv0, tv1;
  int64_t t0, t1, count;
  lvtab *lv;
  double levels[256];

  for(i=0 ; i<256 ; i++) { levels[i] = 255.0 - i ; }
  lv = Vsearch_setup(levels,250);

//   printf("table0 =\n"); for(i=0 ; i<8 ; i++) printf("%10.3f",lv->t0[i]); printf("\n");
  for(T=255.0 ; T>253.9 ; T=T-.5){
    ix = search_list(T,lv);
//     printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
    ix = Vsearch_list(T,lv);
//     printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n\n\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  }
//   exit(0);
//   T = 1.5;
//   ix = search_list(T,lv);
//   printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//   T = 1.0;
//   ix = search_list(T,lv);
//   printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  T = 256.0 ;
  ix = Vsearch_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  T = 0.000005 ;
  ix = Vsearch_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  count = 0;
  gettimeofday(&tv0,NULL);
  for(T=0.000005 ; T<255.0 ; T=T+.000001){
//     ix = Vsearch_list(T,lv); 
    if(Vsearch_list(T,lv) == -1) exit(1);
    count++;
//     if(T > lv->t2[ix] || T < lv->t2[ix+1]){
//       printf("ERROR(Vsearch): T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//       printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
//       exit(1);
//     }
//     ix = search_list(T,lv);
//     if(T > lv->t2[ix] || T < lv->t2[ix+1]){
//       printf("ERROR(search): T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//       printf("               T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
//       exit(1);
//     }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMD   : %ld microseconds for %ld iterations, %f us/iter\n",t1-t0,count,T);
  gettimeofday(&tv0,NULL);
  for(T=0.000005 ; T<255.0 ; T=T+.000001){
    if(search_list(T,lv) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SCALAR : %ld microseconds for %ld iterations, %f us/iter\n",t1-t0,count,T);
}
#endif
