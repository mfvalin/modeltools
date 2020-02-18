//     functions for C and FORTRAN programming
//     Copyright (C) 2020  Recherche en Prevision Numerique
// 
//     This software is free software; you can redistribute it and/or
//     modify it under the terms of the GNU Lesser General Public
//     License as published by the Free Software Foundation,
//     version 2.1 of the License.
// 
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//     Lesser General Public License for more details.
//
//****P* librkl/RPN kernel library smoothing routines
// DESCRIPTION
//  set of Fortran/C callable routines to perform in place (1 2 1)/4 1 or 2 dimensional smoothing
//  smoothed value = (previous_point + 2 * current_point + next_point) / 4
//  in the 2D case, FORTRAN storage order is assumed
//  one smoothing pass will remove any signal with a wavelength of 2 grid points,
//  attenuate by 50% any signal with a wavelength of 4 grid points
//  attenuate by 25% any signal with a wavelength of 8 grid points
//****
#include <stdint.h>
static inline uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

// 1 D horizontal smoothing pass used by lissage_121_2D
// src : data source
// dst : result
// n   : number of points
// when the back point is not available, (i-1=0 or i+1=ni index) a mirror condition is assumed
// i.e. src[-1] is assumed equal to src[1], src[ni] is assumed equal to src[ni-2]
static void inline lissage_121_1D_inline(float *dst, float *src, int n){  // 1 2 1 smoothing pass along i
  int i;
  dst[0] = .5f * (src[0] + src[1]);
  for(i=1 ; i<n-1 ; i++) dst[i] = .25f * src[i-1] + .5f * src[i] + .25f * src[i+1];
  dst[n-1] = .5f * (src[n-2] + src[n-1]);
}

//****f* librkl/2D linear smoothing
// Synopsis
// smoothing is performed IN-PLACE, (1 2 1)*.25 smoothing
// one pass from bottom of array towards the top to minimize memory traffic
// the horizontal smoothing pass is performed when storing result
// two temporary rows are needed
// storage of result into source array is delayed until replaced row is no longer needed for computations
// when the back/forward point is not available, (-1, ni, or nj index) a mirror condition is assumed
// row 1 is used where row -1 would be needed, row nj-2 is used where row nj would be needed
//
// Fortran interface
//   subroutine lissage_121_2d(Z, ni, nj, li) bind(C,name='lissage_121_2D')  !InTf!
//     import :: C_INT, C_FLOAT                               !InTf!
//     real(C_FLOAT), dimension(nj,nj), intent(INOUT) :: z    !InTf!
//     integer(C_INT), intent(IN), value :: ni                !InTf!
//     integer(C_INT), intent(IN), value :: nj                !InTf!
//     integer(C_INT), intent(IN), value :: li                !InTf!
//   end subroutine lissage_121_2D                            !InTf!
// ARGUMENTS
// z  : real(float) array of dimension (ni,nj) (Fortran order assumed)
// li : number of useful points in row
// ni : row dimension
// nj : number of rows
void lissage_121_2D(float *z, int ni, int nj, int li){
//****
  float t1[li], t2[li];  // dynamic local arrays 
  float *bot, *row, *top, *new, *old, *t;
  int i, j;
  old = t1 ; new = t2;
  row = z ; top = z + ni;                                  // rows 0 and 1
  for(i=0 ; i<li ; i++) old[i] = .5f * (row[i] + top[i]);  // bottom, row 0, smooth along j

  for(j=1; j<nj-1 ; j++){              // middle rows, smoothing along j, then i
    bot = row; row = top ; top += ni;  // current set of rows (bot, row, top : j-1 , j , j+1)
    for(i=0 ; i<li ; i++) new[i] = .25f * bot[i] + .5f * row[i] + .25f * top[i];  // row j, smooth along j
    lissage_121_1D_inline(bot, old, li);      // row j-1, smooth along i and store (row j-1 no longer needed as a source)
    t = old ; old = new ; new = t;     // swap old and new
  }
  for(i=0 ; i<li ; i++) new[i] = .5f * row[i] + .5f * top[i];   // top, row nj-1, smooth along j
  lissage_121_1D_inline(row, old, li);        // row nj-2 , smooth along i and store
  lissage_121_1D_inline(top, new, li);        // row nj-1 , smooth along i and store
}

//****f* librkl/1D linear smoothing
// Synopsis
// smoothing is performed IN-PLACE, (1 2 1)*.25 smoothing
// when the back/forward point is not available, (-1, niindex) a mirror condition is assumed
// point 1 is used where point -1 would be needed, point ni-2 is used where point nj would be needed
//
// Fortran interface
//   subroutine lissage_121_1d(Z, ni) bind(C,name='lissage_121_1D')  !InTf!
//     import :: C_INT, C_FLOAT                               !InTf!
//     real(C_FLOAT), dimension(nj), intent(INOUT) :: z       !InTf!
//     integer(C_INT), intent(IN), value :: ni                !InTf!
//   end subroutine lissage_121_1D                            !InTf!
// ARGUMENTS
// z  : real(float) array of dimension (ni)
// ni : row dimension
void lissage_121_1D(float *z, int ni){
//****
  int i;
  for(i=0 ; i<ni-1; i++) z[i] = 0.5f * (z[i] + z[i+1]);
  z[ni-1] = z[ni-2];
  for(i=ni-2 ; i>0 ; i--) z[i] = 0.5f * (z[i] + z[i-1]);
}

#if defined(SELF_TEST)
#include <stdio.h>
int main(){
  float z1[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f };  // 2 deltax
  float z2[] = {1.0f, 0.0f,-1.0f, 0.0f, 1.0f, 0.0f,-1.0f, 0.0f, 1.0f, 0.0f,-1.0f, 0.0f };  // 4 deltax
  float z3[] = {1.0f, 1.0f, 0.0f, 0.0f,-1.0f,-1.0f, 0.0f, 0.0f, 
                1.0f, 1.0f, 0.0f, 0.0f,-1.0f,-1.0f, 0.0f, 0.0f };   // 8 deltax
  float z4[] = {0.0f,-1.0f,-1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f,
                0.0f,-1.0f,-1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f };   // 8 deltax (z3 with offset)
  int i;
  for(i=0 ; i<12 ; i++) printf(" %4.2f",z1[i]) ; printf("\n");
  lissage_121_1D(z1, 12);
  for(i=0 ; i<12 ; i++) printf(" %4.2f",z1[i]) ; printf("\n\n");
  for(i=0 ; i<12 ; i++) printf(" %4.2f",z2[i]) ; printf("\n");
  lissage_121_1D(z2, 12);
  for(i=0 ; i<12 ; i++) printf(" %4.2f",z2[i]) ; printf("\n\n");
  for(i=0 ; i<16 ; i++) printf(" %4.2f",z3[i]) ; printf("\n");
  lissage_121_1D(z3, 16);
  for(i=0 ; i<16 ; i++) printf(" %4.2f",z3[i]) ; printf("\n\n");
  for(i=0 ; i<16 ; i++) printf(" %4.2f",z4[i]) ; printf("\n");
  lissage_121_1D(z4, 16);
  for(i=0 ; i<16 ; i++) printf(" %4.2f",z4[i]) ; printf("\n\n");
}
#endif

