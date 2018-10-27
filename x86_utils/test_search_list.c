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
  double t1[64];   // will support at most 255 levels
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

lvtab * Vsearch_setup(double *targets, int n);
// int Vsearch_list_inc_2(double target, lvtab *lv);
int Vsearch_list_inc(double target, lvtab *lv);
double Vsearch_list_inc_d(double target, lvtab *lv);
// int Vsearch_list_dec_2(double target, lvtab *lv);
int Vsearch_list_dec(double target, lvtab *lv);
// straight C version for comparison purposes
int search_list_dec(double target, lvtab *lv);
int search_list_inc(double target, lvtab *lv);
int Vcoef_xyz_inc(double *cxyz, double px, double py, double pz, lvtab *lv);

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
  int i, ix;
  double T;
  struct timeval tv0, tv1;
  int64_t t0, t1;
  int32_t count, cnt;
  lvtab *lv, *lv2;
  double levels[256];
  double levels2[256];
  double cxyz[24];
//   double tg[4];

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

  T = NT  ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//   ix = Vsearch_list_dec_2(T,lv);
//   printf("Tdec2= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
//   ix = Vsearch_list_inc_2(T,lv2);
  ix = Vcoef_xyz_inc(cxyz, .5, .5, T, lv2);
  printf("Tinc2= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);

  T = NT+.5  ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//   ix = Vsearch_list_dec_2(T,lv);
//   printf("Tdec2= %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
//   ix = Vsearch_list_inc_2(T,lv2);
  ix = Vcoef_xyz_inc(cxyz, .5, .5, T, lv2);
  printf("Tinc2= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);

  T = 0.5 ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//   ix = Vsearch_list_dec_2(T,lv);
//   printf("Tdec2= %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
//   ix = Vsearch_list_inc_2(T,lv2);
  ix = Vcoef_xyz_inc(cxyz, .5, .5, T, lv2);
  printf("Tinc2= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  printf("\n");
// exit (0);
  count = (NT - 1) * 1000000;
  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if(Vsearch_list_dec(T,lv) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDd  : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
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
  printf("SIMDd+ : %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

//   cnt = 0;
//   gettimeofday(&tv0,NULL);
//   for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(Vsearch_list_dec_2(T,lv) == -1) exit(1);
//   }
//   gettimeofday(&tv1,NULL);
//   t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
//   t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
//   T = t1-t0; T /= count;
//   printf("SIMDd2 : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);
// 
//   cnt = 0;
//   gettimeofday(&tv0,NULL);
//   for(T=1.000000 ; T<NT ; T=T+.000001){
//     ix = Vsearch_list_dec_2(T,lv); ; // cnt += ix;
//     if(T > lv->t2[ix] || T < lv->t2[ix+1]){
//       printf("ERROR(Vdecr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//       printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
//       exit(1);
//     }
//   }
//   gettimeofday(&tv1,NULL);
//   t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
//   t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
//   T = t1-t0; T /= count;
//   printf("SIMDd2+: %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if((ix = Vsearch_list_inc(T,lv2)) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDi  : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
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
  printf("SIMDi+ : %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if((ix = Vsearch_list_inc_d(T,lv2) ) == -1) exit(1);
//     if((ix = Vcoef_xyz_inc(cxyz, .5, .5, T, lv2) ) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDi2 : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    ix = Vsearch_list_inc_d(T,lv2);
//     ix = Vcoef_xyz_inc(cxyz, .5, .5, T, lv2);
    if(T < lv2->t2[ix] || T > lv2->t2[ix+1]){
      printf("ERROR(Vincr2) : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
      printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv2->t2[ix], T-lv2->t2[ix+1]);
      exit(1);
    }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDi2+: %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if(search_list_dec(T,lv) == -1) exit(1);
//     ix = search_list_dec(T,lv); // cnt += ix;
//     if(T > lv->t2[ix] || T < lv->t2[ix+1]){
//       printf("ERROR(Vdecr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
//       printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
//       exit(1);
//     }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SCALARd : %ld microseconds for %d iterations, %f ns/iter\n",t1-t0,count,T*1000.0);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if(search_list_inc(T,lv2) == -1) exit(1);
//     ix = search_list_inc(T,lv2); // cnt += ix;
//     if(T < lv2->t2[ix] || T > lv2->t2[ix+1]){
//       printf("ERROR(incr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
//       printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv2->t2[ix], T-lv2->t2[ix+1]);
//       exit(1);
//     }
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SCALARi : %ld microseconds for %d iterations, %f ns/iter\n",t1-t0,count,T*1000.0);
}
