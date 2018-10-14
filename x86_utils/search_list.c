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
#include <immintrin.h>

typedef struct{
  double t0[8];
  double t1[64];
  double t2[256];
}lvtab;

lvtab *lv;

typedef struct{
  uint64_t t0[8];
  uint64_t t1[64];
  uint64_t t2[256];
}lvtabi;

static int64_t ONE = 1;

// SIMD version
int Vsearch_list(double target, lvtab *lv){
  __m256d t;
  __m256d tbl0, tbl1;
  __m256i v0, v1, vc;
  int m0, m1, pop0, pop1, ix, j;
  lvtabi *lvi = (lvtabi *) lv;

//   printf("target = %f\n",target);
  t    = _mm256_broadcast_sd(&target);
  vc   = (__m256i) _mm256_broadcast_sd((double *) &ONE);
  t    = (__m256d) _mm256_sub_epi64((__m256i) t, vc);  // we want target > tbl so we subtract 1 from t

  tbl0 = _mm256_loadu_pd(&(lv->t0[0]));
  tbl1 = _mm256_loadu_pd(&(lv->t0[4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);  // t - tbl >= 0 if (t+1) > tbl
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = _mm_popcnt_u32(m1);
  j = pop0+pop1-1;
  ix = j << 3;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t0[j-1],lv->t0[j],lv->t1[ix],lv->t1[ix+7]);

  tbl0 = _mm256_loadu_pd(&(lv->t1[ix  ]));
  tbl1 = _mm256_loadu_pd(&(lv->t1[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);
  v1   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);
  m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);
  pop1 = _mm_popcnt_u32(m1);
  j = pop0+pop1-1;
  ix = (ix + j) <<2 ;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j],lv->t2[ix],lv->t2[ix+7]);

  tbl0 = _mm256_loadu_pd(&(lv->t2[ix  ]));
//   tbl1 = _mm256_loadu_pd(&(lv->t2[ix+4]));
  v0   = _mm256_sub_epi64((__m256i) t,(__m256i) tbl0);
//   v1   = _mm256_sub_epi32((__m256i) t,(__m256i) tbl1);
  m0   = _mm256_movemask_pd((__m256d) v0);
//   m1   = _mm256_movemask_pd((__m256d) v1);
  pop0 = _mm_popcnt_u32(m0);
//   pop1 = _mm_popcnt_u32(m1);
  pop1 = 0;
  j = pop0+pop1-1;
  ix = (ix + j) ;
//   printf("j = %d, ix = %d, %f %f, %f %f\n",j,ix,lv->t1[j-1],lv->t1[j],lv->t2[ix],lv->t2[ix+7]);

  return ix;
}

// straight C version
int search_list(double target, lvtab *lv){
  int ix = 0;
  int i, j;

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

#if defined(SELF_TEST)
main(){
  int i, status, ix;
  double T;

  status = posix_memalign((void **) &lv, 64, sizeof(lvtab));
  for(i=0 ; i<256 ; i++) { lv->t2[i] = 255.0 - i ; }
  for(i=0 ; i<64 ; i++) { lv->t1[i] = lv->t2[(i<<2)] ; } ;
  for(i=0 ; i<8  ; i++) { lv->t0[i] = lv->t2[(i<<5)] ; } ;

//   printf("table0 =\n"); for(i=0 ; i<8 ; i++) printf("%10.3f",lv->t0[i]); printf("\n");
  for(T=255.0 ; T>253.9 ; T=T-.5){
    ix = search_list(T,lv);
    printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
    ix = Vsearch_list(T,lv);
    printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n\n\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  }
//   exit(0);
  T = 1.5;
  ix = search_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  T = 1.0;
  ix = search_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  T = 0.5;
  ix = search_list(T,lv);
  printf("T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  for(T=0.05 ; T<255.0 ; T=T+.000001){
    ix = Vsearch_list(T,lv); if(ix == -1) exit(1);
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
}
#endif
