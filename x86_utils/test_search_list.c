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
  double t2m1;     // normally equal to t2[0]
  double t2[257];
  double top;      // t2[nk-2]
  double *odz;     // pointer to inverse of deltaz array [nk doubles]
  double *ocz;     // pointer to inverse of coefficient denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of planes (max 255)
  uint32_t nij;    // ni*nj, distance between planes
}lvtab;

lvtab * Vsearch_setup(double *targets, int n, int ni, int nj);
// int Vsearch_list_inc_2(double target, lvtab *lv);
int Vsearch_list_inc(double target, lvtab *lv);
double Vsearch_list_inc_d(double target, lvtab *lv);
// int Vsearch_list_dec_2(double target, lvtab *lv);
int Vsearch_list_dec(double target, lvtab *lv);
// straight C version for comparison purposes
int search_list_dec(double target, lvtab *lv);
int search_list_inc(double target, lvtab *lv);
int Vcoef_xyz_inc(uint32_t *ixyz, double *cxyz, double *PX, double *PY, double *PZ, lvtab *lv, int n);
int Vcoef_xyz_inc_o(uint32_t *ixyz, double *cxyz, double *PX, double *PY, double *PZ, lvtab *lv, int n);
int Vcoef_xyz_incr(uint32_t *ixyz, double *cxyz, double *PX, double *PY, double *PZ, lvtab *lv, int n);
int Vsearch_list_inc_n(int *ix, double *target, lvtab *lv, int n);
double Vsearch_list_inc_dn(double *list, double *target, lvtab *lv, int n);

uint64_t rdtscp_(void) {   // version "in order" avec "serialization"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

void Tricubic_coeffs_d(double *px, double *py, double *pz, double x, double y, double z){
static double cp167 =  1.0/6.0;
static double cm167 = -1.0/6.0;
static double cm5 = -0.5;
static double cp5 = 0.5;
static double one = 1.0;
static double two = 2.0;

  pz[0] = cm167*z*(z-one)*(z-two);        // coefficients for interpolation along z
  pz[1] = cp5*(z+one)*(z-one)*(z-two);
  pz[2] = cm5*z*(z+one)*(z-two);
  pz[3] = cp167*z*(z+one)*(z-one);
  pz[4] = 0.0;
  pz[5] = 1.0 - z;
  pz[6] = z;
  pz[7] = 0.0;

  py[0] = cm167*y*(y-one)*(y-two);        // coefficients for interpolation along y
  py[1] = cp5*(y+one)*(y-one)*(y-two);
  py[2] = cm5*y*(y+one)*(y-two);
  py[3] = cp167*y*(y+one)*(y-one);
  py[4] = 0.0;
  py[5] = 1.0 - y;
  py[6] = y;
  py[7] = 0.0;

  px[0] = cm167*x*(x-one)*(x-two);        // coefficients for interpolation along x
  px[1] = cp5*(x+one)*(x-one)*(x-two);
  px[2] = cm5*x*(x+one)*(x-two);
  px[3] = cp167*x*(x+one)*(x-one);
  px[4] = 0.0;
  px[5] = 1.0 - x;
  px[6] = x;
  px[7] = 0.0;
}

#include <sys/time.h>
#define NT 251
int main(){
  int i, ix, iy, iz;
  double T;
  struct timeval tv0, tv1;
  int64_t t0, t1;
  int32_t count, cnt;
  lvtab *lv, *lv2;
  double levels[256];
  double levels2[256];
  double cxyz[240];
  uint32_t ixyz[10];
  double point5 = .5;
  double pxs[10], pys[10], pzs[10];
  double xt, yt, zt;
  double rxyz[24];      // reference coeffs
  int ixs[10];

  for(i=0 ; i<10 ; i++) { pxs[i] = 2.5 ; pys[i] = 2.4 ; pzs[i] = 1.3 ;}
  for(i=0 ; i<NT ; i++) { levels[i] = NT - i ; levels2[i] = i + 1 ;}
  lv  = Vsearch_setup(levels, NT, 10, 20);
  lv2 = Vsearch_setup(levels2,NT, 10, 20);   // ascending levels

#if defined(CHECK)
  for(i = 0 ; i < 3 ; i++){
    ix = pxs[0] ; xt = pxs[0] - ix;
    iy = pys[0] ; yt = pys[0] - iy;
    iz = pzs[0] ; zt = pzs[0] - iz;
    printf("ix, iy, iz, ixyz = %d %d %d %d\n",ix,iy,iz,(ix-1) + (iy-1)*10 + (iz-1) * 200);
    printf("Ta : %8.5f %8.5f %8.5f, dxyz = %8.5f %8.5f %8.5f\n\n",pxs[0],pys[0],pzs[0],xt,yt,zt);
    Tricubic_coeffs_d(rxyz, rxyz+8, rxyz+16, xt, yt, zt);
    printf("rx : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[ 0],rxyz[ 1],rxyz[ 2],rxyz[ 3],rxyz[ 0]+rxyz[ 1]+rxyz[ 2]+rxyz[ 3]);
    printf("ry : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[ 8],rxyz[ 9],rxyz[10],rxyz[11],rxyz[ 8]+rxyz[ 9]+rxyz[10]+rxyz[11]);
    printf("rz : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[16],rxyz[17],rxyz[18],rxyz[19],rxyz[16]+rxyz[17]+rxyz[18]+rxyz[19]);
    printf("\n");
    ix = Vcoef_xyz_inc(ixyz, rxyz, pxs, pys, pzs, lv2, 1);
    printf("ix, ixyz[0] = %d %d\n",ix,ixyz[0]);
    printf("cx : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[ 0],rxyz[ 1],rxyz[ 2],rxyz[ 3],rxyz[ 0]+rxyz[ 1]+rxyz[ 2]+rxyz[ 3]);
    printf("cy : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[ 8],rxyz[ 9],rxyz[10],rxyz[11],rxyz[ 8]+rxyz[ 9]+rxyz[10]+rxyz[11]);
    printf("cz : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[16],rxyz[17],rxyz[18],rxyz[19],rxyz[16]+rxyz[17]+rxyz[18]+rxyz[19]);
    printf("czl: %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[20],rxyz[21],rxyz[22],rxyz[23],rxyz[20]+rxyz[21]+rxyz[22]+rxyz[23]);
    printf("\n");
    ix = Vcoef_xyz_incr(ixyz, rxyz, pxs, pys, pzs, lv2, 1);
    printf("ix, ixyz[0] = %d %d\n",ix,ixyz[0]);
    printf("cx : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[ 0],rxyz[ 1],rxyz[ 2],rxyz[ 3],rxyz[ 0]+rxyz[ 1]+rxyz[ 2]+rxyz[ 3]);
    printf("cy : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[ 8],rxyz[ 9],rxyz[10],rxyz[11],rxyz[ 8]+rxyz[ 9]+rxyz[10]+rxyz[11]);
    printf("cz : %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[16],rxyz[17],rxyz[18],rxyz[19],rxyz[16]+rxyz[17]+rxyz[18]+rxyz[19]);
    printf("czl: %8.5f %8.5f %8.5f %8.5f %8.5f \n",rxyz[20],rxyz[21],rxyz[22],rxyz[23],rxyz[20]+rxyz[21]+rxyz[22]+rxyz[23]);
    printf("========================================================================\n");
    pxs[0] += 1.33;
    pys[0] += 5.44;
    pzs[0] += 124.5;
  }
#endif
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

// exit(0);
  T = NT-0.5  ;
#if defined(CHECK)
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = Vcoef_xyz_inc(ixyz, cxyz, pxs, pys, &T, lv2, 1);
  printf("Tcoef= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);

  T = NT+.5  ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = Vcoef_xyz_inc(ixyz, cxyz, pxs, pys, &T, lv2, 1);
  printf("Tcoef= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);

  T = 0.5 ;
  ix = Vsearch_list_dec(T,lv);
  printf("Tdec = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = search_list_dec(T,lv);
  printf("Td   = %f, index = %d , tbl[index] = %g, tbl[index+1] = %g\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
  ix = Vsearch_list_inc(T,lv2);
  printf("Tinc = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = Vcoef_xyz_inc(ixyz, cxyz, pxs, pys, &T, lv2, 1);
  printf("Tcoef= %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
  ix = search_list_inc(T,lv2);
  printf("Ti   = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n\n",T, ix, lv2->t2[ix], lv2->t2[ix+1]);
// exit (0);
#endif
  count = (NT - 1) * 1000000;
//   cnt = 0;
//   gettimeofday(&tv0,NULL);
//   for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(Vsearch_list_dec(T,lv) == -1) exit(1);
//   }
//   gettimeofday(&tv1,NULL);
//   t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
//   t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
//   T = t1-t0; T /= count;
//   printf("SIMDd  : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);
// 
//   cnt = 0;
//   gettimeofday(&tv0,NULL);
//   for(T=1.000000 ; T<NT ; T=T+.000001){
//     ix = Vsearch_list_dec(T,lv); ; // cnt += ix;
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
//   printf("SIMDd+ : %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

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
  printf("SCALARi+: %ld microseconds for %d iterations, %f ns/iter\n",t1-t0,count,T*1000.0);

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
  for(T=1.000000 ; T<NT ; T=T+.00001){
    for(i=0 ; i<10 ; i++) {pzs[i] = T + i*.000001 ; }
    if((ix = Vsearch_list_inc_n(ixs,pzs,lv2,10)) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDin : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if(Vsearch_list_inc_d(T,lv2) < 0) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDid : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.00001){
    for(i=0 ; i<10 ; i++) {pzs[i] = T + i*.000001 ; }
    if(Vsearch_list_inc_dn(pzs, pzs, lv2, 10) < 0) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDidn: %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

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
    if((ix = Vcoef_xyz_inc(ixyz, cxyz,  pxs, pys, &T, lv2, 1) ) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDi2 : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    ix = Vcoef_xyz_inc(ixyz, cxyz, pxs, pys, &T, lv2, 1);
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
  printf("SIMDi2+: %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    if((ix = Vcoef_xyz_incr(ixyz, cxyz,  pxs, pys, &T, lv2, 1) ) == -1) exit(1);
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDp2 : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.000001){
    ix = Vcoef_xyz_incr(ixyz, cxyz, pxs, pys, &T, lv2, 1);
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
  printf("SIMDp2+: %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.00001){
    for(i=0 ; i<10 ; i++) {pzs[i] = T + i*.000001 ; }
    if((ix = Vcoef_xyz_inc(ixyz, cxyz, pxs, pys, pzs, lv2, 10)) == -1) exit(1);;
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDx2 : %ld microseconds for %d iterations, %f ns/iter (%d)\n",t1-t0,count,T*1000.0,cnt);

  cnt = 0;
  gettimeofday(&tv0,NULL);
  for(T=1.000000 ; T<NT ; T=T+.00001){
    for(i=0 ; i<10 ; i++) {pzs[i] = T + i*.000001 ; }
    if((ix = Vcoef_xyz_incr(ixyz, cxyz, pxs, pys, pzs, lv2, 10)) == -1) exit(1);;
  }
  gettimeofday(&tv1,NULL);
  t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
  t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
  T = t1-t0; T /= count;
  printf("SIMDp2x: %ld microseconds for %d iterations, %f ns/iter (%d)\n\n",t1-t0,count,T*1000.0,cnt);

//   cnt = 0;
//   gettimeofday(&tv0,NULL);
//   for(T=1.000000 ; T<NT ; T=T+.000001){
//     if(search_list_dec(T,lv) == -1) exit(1);
// //     ix = search_list_dec(T,lv); // cnt += ix;
// //     if(T > lv->t2[ix] || T < lv->t2[ix+1]){
// //       printf("ERROR(Vdecr)  : T = %f, index = %d , tbl[index] = %f, tbl[index+1] = %f\n",T, ix, lv->t2[ix], lv->t2[ix+1]);
// //       printf("                T - tbl[index] = %f, T - tbl[index+1] = %f\n",T-lv->t2[ix], T-lv->t2[ix+1]);
// //       exit(1);
// //     }
//   }
//   gettimeofday(&tv1,NULL);
//   t0 = tv0.tv_sec; t0 *= 1000000 ; t0 += tv0.tv_usec ;
//   t1 = tv1.tv_sec; t1 *= 1000000 ; t1 += tv1.tv_usec ;
//   T = t1-t0; T /= count;
//   printf("SCALARd : %ld microseconds for %d iterations, %f ns/iter\n",t1-t0,count,T*1000.0);
}
