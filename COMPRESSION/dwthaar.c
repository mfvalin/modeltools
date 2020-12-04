/* 
 * Copyright (C) 2020  Recherche en Prevision Numerique
 *                     Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 */
//****P* librkl/wavelet-transforms
// Synopsis
//
// Haar wavelets
// https://en.wikipedia.org/wiki/Haar_wavelet
//
// this code is using a lifting implementation
// https://en.wikipedia.org/wiki/Lifting_scheme
//
// 1 dimensional transform, "in place", "even/odd split", or "in place with even/odd split"
//   original data
//   +--------------------------------------------------------+
//   |                  N data                                |
//   +--------------------------------------------------------+
//
//   transformed data (in place, no split, even number of data)
//   +--------------------------------------------------------+
//   |   N data, even/odd, even/odd, ..... , even/odd         +
//   +--------------------------------------------------------+
//
//   transformed data (in place, no split, odd number of data)
//   +--------------------------------------------------------+
//   |   N data, even/odd, even/odd, ..... , even/odd, even   +
//   +--------------------------------------------------------+
//
//   transformed data, in place with even/odd split
//   +--------------------------------------------------------+
//   | (N+1)/2 even data            |    (N/2) odd data       |
//   +--------------------------------------------------------+
//
//   original data                     transformed data (2 output arrays)
//   +------------------------------+  +-------------------+  +------------------+
//   |             N data           |  | (N+1)/2 even data |  |   N/2 odd data   |
//   +------------------------------+  +-------------------+  +------------------+
//
//   even data are the "approximation" terms ("low frequency" terms)
//   odd data are the "detail" terms         ("high frequency" terms)
//
// 2 dimensional in place with 2 D split
//   original data                               transformed data (in same array)
//   +------------------------------------+      +-------------------+----------------+
//   |                  ^                 |      +                   |                |
//   |                  |                 |      +   even i/odd j    |  odd i/odd j   |
//   |                  |                 |      +                   |                |
//   |                  |                 |      +                   |                |
//   |                  |                 |      +-------------------+----------------+
//   |               NJ data              |      +                   |                |
//   |                  |                 |      +                   |                |
//   |                  |                 |      +   even i/even j   |  odd i/even j  |
//   |<----- NI data ---|---------------->|      +                   |                |
//   |                  v                 |      +                   |                |
//   +------------------------------------+      +-------------------+----------------+
//   the process can be applied again to the even/even transformed part to achieve a multi level transform
//
//****

#include <stdio.h>
#include <dwthaar.h>

#define A      1.0f
#define B      0.5f

//****f* librkl/F_DWTHAAR_1D_split_N_even
// Synopsis
//
// Forward DWT transform (analysis)
// n           : number of data points (MUST BE EVEN)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_split_N_even(x, e, o, n) bind(C,name='F_DWTHAAR_1D_split_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_DWTHAAR_1D_split_N_even               !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_split_N_even(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int nodd = n >> 1;

  for(i = 0 ; i < nodd ; i++) {         // n is even, nodd = neven, process i, i+1 pair
    o[i] = x[i+i+1] - A * x[i+i] ;      // predict odd terms
    e[i] = x[i+i] + B * o[i] ;          // update even terms
  }
}

//****f* librkl/F_DWTHAAR_1D_inplace_N_even
// Synopsis
//
// Forward DWT transform (analysis) (in place, x is overwritten)
// n           : number of data points (MUST BE EVEN)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_inplace_N_even(x, e, o, n) bind(C,name='F_DWTHAAR_1D_inplace_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_DWTHAAR_1D_inplace_N_even             !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_inplace_N_even(float *x, int n){    // InTc
//****
  int i;
  
  for (i = 0; i < n - 1; i += 2) {   // n is even, nodd = neven, process i, i+1 pair
    x[i+1] -= A * x[i  ] ;           // predict odd terms
    x[i  ] += B * x[i+1] ;           // update even terms
  }
}

//****f* librkl/F_DWTHAAR_1D_split_N_odd
// Synopsis
//
// Forward DWT transform (analysis)
// n           : number of data points (MUST BE ODD)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_split_N_odd(x, e, o, n) bind(C,name='F_DWTHAAR_1D_split_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_DWTHAAR_1D_split_N_odd                !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_split_N_odd(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;   // neven should be nodd + 1
  int nodd  = n >> 1;

  for(i = 0 ; i < nodd ;    i++) {       // n is odd, nodd = neven-1, process i, i+1 pair
    o[i] = x[i+i+1] - A * x[i+i] ;       // predict odd terms
    e[i] = x[i+i] + B * o[i] ;         // update even terms
  }
  e[neven-1] = x[n-1] + B * o[nodd-1];   // last even term
}

//****f* librkl/F_DWTHAAR_1D_inplace_N_odd
// Synopsis
//
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (MUST BE ODD)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_inplace_N_odd(x, e, o, n) bind(C,name='F_DWTHAAR_1D_inplace_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_DWTHAAR_1D_inplace_N_odd              !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_inplace_N_odd(float *x, int n){    // InTc
//****
  int i;

  for (i = 0; i < n - 1; i += 2) {     // n is odd, nodd = neven-1, process i, i+1 pair
    x[i+1] -= A * x[i  ] ;             // predict odd terms
    x[i  ] += B * x[i+1] ;             // update even terms
  }
  x[n - 1] += B * x[n - 2];            // last term is even
}
  
//****f* librkl/F_DWTHAAR_1D_split
// Synopsis
//
// Forward DWT transform (analysis)
// n           : number of data points (even or odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_split(x, e, o, n) bind(C,name='F_DWTHAAR_1D_split')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_DWTHAAR_1D_split                        !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_split(float *x, float *e, float *o, int n){    // InTc
//****
  if(n < 3) {
    if(n > 0) e[0] = x[0];
    if(n > 1) o[0] = x[1];
    return;
  }
  if(n & 1){
    F_DWTHAAR_1D_split_N_odd(x, e, o, n);
  }else{
    F_DWTHAAR_1D_split_N_even(x, e, o, n);
  }
}
  
//****f* librkl/F_DWTHAAR_1D_inplace
// Synopsis
//
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_inplace(x, e, o, n) bind(C,name='F_DWTHAAR_1D_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_DWTHAAR_1D_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_inplace(float *x, int n){    // InTc
//****
  if(n < 3) return;
  if(n & 1){
    F_DWTHAAR_1D_inplace_N_odd(x, n);
  }else{
    F_DWTHAAR_1D_inplace_N_even(x, n);
  }
}
  
//****f* librkl/F_DWTHAAR_1D_split_inplace
// Synopsis
//
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_DWTHAAR_1D_split_inplace(x, e, o, n) bind(C,name='F_DWTHAAR_1D_split_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_DWTHAAR_1D_split_inplace              !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_1D_split_inplace(float *x, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int lasteven = neven - 1;
  int nodd  = n >> 1;
  int lastodd = nodd - 1;
  float odd[n];

  if(n < 3) return;  
  for (i = 0; i < nodd; i ++) { 
    odd[i] = x[i+i+1] - A * x[i+i] ;         // predict odd terms and copy them to temporary odd array
    x[i]   = x[i+i  ] + B * odd[i] ;         // update even terms and copy them to beginning of array x
  }
  x[lasteven] = x[lasteven+lasteven] + B * odd[lastodd] ;  // neven may be one greater than nodd

  for (i = 0; i < nodd; i ++) x[i+neven] = odd[i]  ;    // copy odd terms to upper part of x
}

//****f* librkl/I_DWTHAAR_1D_split_N_even
// Synopsis
//
// Inverse DWT transform (synthesis)
// n           : number of data points (MUST BE EVEN)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_split_N_even(x, e, o, n) bind(C,name='I_DWTHAAR_1D_split_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_DWTHAAR_1D_split_N_even               !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_split_N_even(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = neven;

  for (i = 0; i < nodd ; i ++) {           // n is even, nodd = neven, process i, i+1 pair
    x[i+i  ] = e[i] - B * o[i] ;           // unupdate even terms
    x[i+i+1] = o[i] + A * x[i+i] ;         // unpredict odd terms
  }
}

//****f* librkl/I_DWTHAAR_1D_inplace_N_even
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (MUST BE EVEN)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_inplace_N_even(x, e, o, n) bind(C,name='I_DWTHAAR_1D_inplace_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_DWTHAAR_1D_inplace_N_even             !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_inplace_N_even(float *x, int n){    // InTc
//****
  int i;
  
  for (i = 0; i < n-1; i += 2) {       // n is even, nodd = neven, process i, i+1 pair
    x[i  ] -= B * x[i+1] ;             // unupdate even terms
    x[i+1] += A * x[i] ;               // unpredict odd terms
  }
}

//****f* librkl/I_DWTHAAR_1D_split_N_odd
// Synopsis
//
// Inverse DWT transform (synthesis)
// n           : number of data points (MUST BE ODD)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_split_N_odd(x, e, o, n) bind(C,name='I_DWTHAAR_1D_split_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_DWTHAAR_1D_split_N_odd                !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_split_N_odd(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;
  for (i = 0; i < nodd; i ++) {            // n is odd, nodd = neven-1, process i, i+1 pair
    x[i+i  ] = e[i] - B * o[i] ;             // unupdate even terms
    x[i+i+1] = o[i] + A * x[i+i] ;           // unpredict odd terms
  }
  x[n-1] = e[neven-1] - B * o[nodd-1];       // last term is even
}

//****f* librkl/I_DWTHAAR_1D_inplace_N_odd
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (MUST BE ODD)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_inplace_N_odd(x, e, o, n) bind(C,name='I_DWTHAAR_1D_inplace_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_DWTHAAR_1D_inplace_N_odd              !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_inplace_N_odd(float *x, int n){    // InTc
//****
  int i;
  float t = x[n-2] ;                     // save last odd value

  for (i = 0; i < n - 2; i += 2) {       // n is odd, nodd = neven-1, process i, i+1 pair
    x[i  ] -= B * x[i+1] ;               // unupdate even terms
    x[i+1] += A * x[i] ;                 // unpredict odd terms
  }
  x[n - 1] -= B * t ;                    // last term is even, unupdate it
}

//****f* librkl/I_DWTHAAR_1D_split
// Synopsis
//
// Inverse DWT transform (synthesis)
// n           : number of data points (even or odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_split(x, e, o, n) bind(C,name='I_DWTHAAR_1D_split')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_DWTHAAR_1D_split                        !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_split(float *x, float *e, float *o, int n){    // InTc
//****
  if(n < 3) {   // 3 points minimum
    if(n > 0) x[0] = e[0] ;
    if(n > 1) x[1] = o[0] ;
    return;
  }
  if(n & 1){
    I_DWTHAAR_1D_split_N_odd(x, e, o, n);
  }else{
    I_DWTHAAR_1D_split_N_even(x, e, o, n);
  }
}

//****f* librkl/I_DWTHAAR_1D_split_inplace
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_split_inplace(x, e, o, n) bind(C,name='I_DWTHAAR_1D_split_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_DWTHAAR_1D_split_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_split_inplace(float *x, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;
  float even[n];
  float *odd = x + neven;

  if(n < 3) return;   // 3 points minimum
  for(i = 0 ; i < neven ; i++) even[i] = x[i] ;             // copy even terms to temporary array
  for(i = 0 ; i < nodd ; i++) x[i+i+1] = odd[i] ;           // unshuffle odd terms
  for(i = 0 ; i < neven ; i++) x[i+i]  = even[i] ;          // unshuffle even terms
  if(n & 1){
    I_DWTHAAR_1D_inplace_N_odd(x, n);
  }else{
    I_DWTHAAR_1D_inplace_N_even(x, n);
  }
}

//****f* librkl/I_DWTHAAR_1D_inplace
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_DWTHAAR_1D_inplace(x, e, o, n) bind(C,name='I_DWTHAAR_1D_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_DWTHAAR_1D_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_1D_inplace(float *x, int n){    // InTc
//****
  if(n < 3) return;   // 3 points minimum
  if(n & 1){
    I_DWTHAAR_1D_inplace_N_odd(x, n);
  }else{
    I_DWTHAAR_1D_inplace_N_even(x, n);
  }
}

#define VCONTRIB(DEST,SCALE,SRC1,SRC2,N) { int i; for(i=0 ; i<N ; i++) {DEST[i] += SCALE *(SRC1[i] + SRC2[i]) ; } }
#define VSCALE(WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHAT[i] *= FAC ; } } 
#define VSCALE2(WHERE,WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHERE[i] = WHAT[i]; } } 

//****f* librkl/F_DWTHAAR_2D_split_inplace
// Synopsis
//  
// Forward DWT transform (analysis) (in place, x is overwritten)
// ni     : number of point in rows
// nj     : number of points in columns
// lni    : storage length of rows
// x[nj][lni]  [C]        array to transform
// x(lni,lnj)  [Fortran]  array to transform
//
// FORTRAN interface
// interface        !InTf
// subroutine F_DWTHAAR_2D_split_inplace(x, ni, lni, nj) BIND(C,name='F_DWTHAAR_2D_split_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine F_DWTHAAR_2D_split_inplace                  !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_2D_split_inplace(float *x, int ni, int lni, int nj){    // InTc
//****
  int jeven = (nj+1) >> 1;
  int jodd = nj >> 1;
  int ieven = (ni+1) >> 1;
  int i, j;
  float *rowe, *rowo, *rowt, *rowx;
  int lni2 = lni+lni;
  float t[((nj>>1)*ni)] ;  // temporary storage for odd rows

  rowx = x ;          // running pointer for even row storage
  rowe = x ;          // even rows from x
  rowo = x + lni ;    // odd rows from x
  rowx = x ;          // running pointer for even row storage
  rowt = &(t[0]) ;    // odd rows in t
  for(j = 0 ; j < jodd ; j++){
    F_DWTHAAR_1D_split(rowo, rowt, rowt + ieven, ni) ;    // first odd row
    if(j == 0){                                           // first even row is done in place
      F_DWTHAAR_1D_split_inplace(rowx, ni) ; 
    }else{                                                // other even rows, packed at bottom of x
      F_DWTHAAR_1D_split(rowe, rowx, rowx + ieven, ni) ; 
    }
    for(i = 0 ; i < ni ; i++) {
      rowt[i] -= A * rowx[i] ;     // predict odd row
      rowx[i] += B * rowt[i] ;     // update even row
    }
    rowe += lni2 ;               // next even row to be fetched
    rowx += lni ;                // next even row to be stored
    rowo += lni2 ;               // next odd row to be fetched
    rowt += ni ;                 // next odd row to be stored
  }
  if(jodd != jeven){                          // nj odd
    rowt -= ni ;       // last odd row
    F_DWTHAAR_1D_split(rowe, rowx, rowx + ieven, ni) ;     // last row is even
    for(i = 0 ; i < ni ; i++) rowx[i] += B * rowt[i] ;     // update last even row
  }
  rowo = x + lni * jeven ;                         // first odd row in x, above even rows
  rowt = &(t[0]) ;                                 // first odd row in t
  for(j = 0 ; j < jodd ; j++){
    for(i = 0 ; i < ni ; i++) rowo[i] = rowt[i] ;  // copy odd rows back to x
    rowt += ni ;
    rowo += lni ;
  }
}

//****f* librkl/I_DWTHAAR_2D_split_inplace
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// ni     : number of point in rows
// nj     : number of points in columns
// lni    : storage length of rows
// x[nj][lni]  [C]        array to transform
// x(lni,lnj)  [Fortran]  array to transform
//
// FORTRAN interface
// interface        !InTf
// subroutine I_DWTHAAR_2D_split_inplace(x, ni, lni, nj) BIND(C,name='I_DWTHAAR_2D_split_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine I_DWTHAAR_2D_split_inplace                        !InTf
// end interface    !InTf
// ARGUMENTS
#define S 1.0
#define Z 1.0
void I_DWTHAAR_2D_split_inplace(float *x, int ni, int lni, int nj){    // InTc
//****
  int jeven = (nj+1) >> 1;
  int jodd  = nj >> 1;
  int ieven = (ni+1) >> 1;
  int i, j;
  float *rowe, *rowo, *rowt, *rowx;
  int lni2 = lni+lni;
  float t[(nj*ni)] ;      // temporary storage for odd rows

  // all even rows are in the bottom part of x, and odd rows are in the upper part of x
  // even rows will be unupdated and moved to t, odd rows will be unpredicted in place
  rowe = x + lni * (jeven - 1) ;         // last valid even row in x
  rowo = x + lni * (nj - 1) ;            // last valid odd row in x is at top of x
  rowt = (&(t[0])) + ni * (jeven - 1) ;  // last valid even row
  if(jodd != jeven){                     // nj is odd, unupdate and copy last row
    for(i = 0 ; i < ni ; i++)  rowt[i] = rowe[i] - B * rowo[i] ; // unupdate last even row first
    rowe -= lni ;                        // update rowe and rowt but rowo will be reused
    rowt -= ni  ;
  }
  for(j = 0 ; j < jodd ; j++){            // top to bottom row sweep
    for(i = 0 ; i < ni ; i++) {
      rowt[i] = rowe[i] - B * rowo[i] ;   // unupdate even row and copy it into t
      rowo[i] = rowo[i] + A * rowt[i] ;   // unpredict odd row
    }
    rowe -= lni ;                         // update rowe, rowo, and rowt 
    rowo -= lni ;
    rowt -= ni ;
  }
  // at this point, all even rows are in t, and odd rows are in the upper part of x
  rowt = &(t[0]) ;              // even rows in t
  rowe = x ;                    // even rows into x
  rowx = x + lni*jeven ;        // odd rows from upper part of x
  rowo = x + lni ;              // odd rows into x
  for(j = 0 ; j < jodd ; j++){
    I_DWTHAAR_1D_split(rowe, rowt, rowt+ieven, ni) ;    // inverse transform even row, move it to proper position in x
    if(rowo == rowx) {
      I_DWTHAAR_1D_split_inplace(rowo, ni) ;            // if nj is even, last odd row must be done in place
    }else{
      I_DWTHAAR_1D_split(rowo, rowx, rowx+ieven, ni) ;  // inverse transform odd row, move it to proper position in x
    }
    rowx += lni ;
    rowo += lni2 ;
    rowt += ni ;
    rowe += lni2 ;
  }
  if(jodd != jeven){    // nj is odd, one more even row to unupdate and copy
    I_DWTHAAR_1D_split(rowe, rowt, rowt+ieven, ni) ;  // inverse transform last even row
  }
}

//****f* librkl/F_DWTHAAR_2D_split_inplace_n
// Synopsis
//  
// Forward DWT transform (analysis) (in place, x is overwritten)
// ni     : number of point in rows
// nj     : number of points in columns
// lni    : storage length of rows
// levels : number of recursive transforms
// x[nj][lni]  [C]        array to transform
// x(lni,lnj)  [Fortran]  array to transform
//
// FORTRAN interface
// interface        !InTf
// subroutine F_DWTHAAR_2D_split_inplace_n(x, ni, lni, nj, levels) BIND(C,name='F_DWTHAAR_2D_split_inplace_n') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni, levels    !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine F_DWTHAAR_2D_split_inplace_n              !InTf
// end interface    !InTf
// ARGUMENTS
void F_DWTHAAR_2D_split_inplace_n(float *x, int ni, int lni, int nj, int levels){    // InTc
//****
  F_DWTHAAR_2D_split_inplace(x, ni, lni, nj);
  if(levels > 1) {
    F_DWTHAAR_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels - 1);
  }
}

//****f* librkl/I_DWTHAAR_2D_split_inplace_n
// Synopsis
//  
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// ni     : number of point in rows
// nj     : number of points in columns
// lni    : storage length of rows
// levels : number of recursive transforms
// x[nj][lni]  [C]        array to transform
// x(lni,lnj)  [Fortran]  array to transform
//
// FORTRAN interface
// interface        !InTf
// subroutine I_DWTHAAR_2D_split_inplace_n(x, ni, lni, nj, levels) BIND(C,name='I_DWTHAAR_2D_split_inplace_n') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni, levels    !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine I_DWTHAAR_2D_split_inplace_n                !InTf
// end interface    !InTf
// ARGUMENTS
void I_DWTHAAR_2D_split_inplace_n(float *x, int ni, int lni, int nj, int levels){    // InTc
//****
  if(levels > 1) {
    I_DWTHAAR_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels - 1);
  }
  I_DWTHAAR_2D_split_inplace(x, ni, lni, nj);
}

#if defined(SELF_TEST)
#include <stdlib.h>
#include <math.h>

#if ! defined(NPTS)

#define NPTS 8
#define NPTSM2 3
#endif
#define NPTSM (NPTS-1)
#define NPTS2 (NPTS/2)
#define NPTSM2 ((NPTS-1)/2)

void set_zero(float *a, float *b, int n){
  int i;
  for(i=0 ; i<n ; i++) { a[i] = 0.0f ; b[i] = 0.0f ; }
}

int main() {
  float x[NPTS+1], y[NPTS+1], e[NPTS+1], o[NPTS+1], z[NPTS+1], d[NPTS+1];
  float xy[NPTS][NPTS] ;
  int i, j, k;
  double sum2;
  float quantum;
  int npts2 = (NPTS+1)/2;
  int npts4 = (npts2+1)/2;

  for (i=0;i<NPTS;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  for (i=0;i<NPTS;i++) y[i]=x[i];
  for (i=0;i<NPTS;i++) d[i]=x[i];
  for (i=0;i<NPTS;i++) { e[i] = 0 ; o[i] = 0 ; }
  printf("------------------------------------ in place even ------------------------------------\n");
  printf("O ");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_inplace_N_even(x, NPTS);
  printf("F0");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  I_DWTHAAR_1D_inplace_N_even(x, NPTS);
  printf("I0");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("O ");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_inplace(x, NPTS);
  printf("F1");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  I_DWTHAAR_1D_inplace(x, NPTS);
  printf("I1");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("------------------------------------ split even ------------------------------------\n");
  printf("O ");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_split_N_even(x, e, o, NPTS);
  printf("F2");
  for(i=0 ; i<NPTS2 ; i++) printf(" %10.5g %10.5g",e[i],o[i]); printf("\n");
  I_DWTHAAR_1D_split_N_even(x, e, o, NPTS);
  printf("I2");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  set_zero(e, o, NPTS);
  printf("O ");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_split(x, e, o, NPTS);
  printf("F3");
  for(i=0 ; i<NPTS2 ; i++) printf(" %10.5g %10.5g",e[i],o[i]); printf("\n");
  I_DWTHAAR_1D_split(x, e, o, NPTS);
  printf("I3");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("O ");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_split_inplace(x, NPTS);
  printf("F4");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  I_DWTHAAR_1D_split_inplace(x, NPTS);
  printf("I4");
  for(i=0 ; i<NPTS ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("------------------------------------ in place odd ------------------------------------\n");
  for (i=0;i<NPTSM;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  printf("O ");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_inplace_N_odd(x, NPTSM);
  printf("F5");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  I_DWTHAAR_1D_inplace_N_odd(x, NPTSM);
  printf("I5");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("O ");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_inplace(x, NPTSM);
  printf("F6");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  I_DWTHAAR_1D_inplace(x, NPTSM);
  printf("I6");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("------------------------------------ split odd ------------------------------------\n");
  printf("O ");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_split_N_odd(x, e, o, NPTSM);
  printf("F7");
  for(i=0 ; i<NPTSM2 ; i++) printf(" %10.5g %10.5g",e[i],o[i]); printf(" %10.5g",e[NPTSM2]) ; printf("\n");
  I_DWTHAAR_1D_split_N_odd(x, e, o, NPTSM);
  printf("I7");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  set_zero(e, o, NPTS);
  printf("O ");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_split(x, e, o, NPTSM);
  printf("F8");
  for(i=0 ; i<NPTSM2 ; i++) printf(" %10.5g %10.5g",e[i],o[i]); printf(" %10.5g",e[NPTSM2]) ; printf("\n");
  I_DWTHAAR_1D_split(x, e, o, NPTSM);
  printf("I8");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  printf("O ");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  F_DWTHAAR_1D_split_inplace(x, NPTSM);
  printf("F9");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  I_DWTHAAR_1D_split_inplace(x, NPTSM);
  printf("I9");
  for(i=0 ; i<NPTSM ; i++) printf(" %10.5g",x[i]); printf("\n");
  printf("\n");

  for (j=0;j<NPTS;j++) {  // Make a fancy cubic signal
//     for (i=0;i<NPTS;i++) xy[j][i] = (3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j);
//     for (i=0;i<NPTS;i++) xy[j][i] = 5+j+0.4*j*j-0.02*j*j*j;
//     for (i=0;i<NPTS;i++) xy[j][i] = 5+i+0.4*i*i-0.02*i*i*i;
    for (i=0;i<NPTS;i++) xy[j][i] = sqrt((i-7.45)*(i-7.45) + (j-7.55)*(j-7.55));
  }
  
  printf("Original 2D signal (even dimensions):\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
  F_DWTHAAR_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, 1);
//   F_DWTHAAR_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   F_DWTHAAR_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   F_DWTHAAR_2D_split_inplace((float *)xy, npts4, NPTS, npts4);

  printf("Transformed 2D signal:\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }

  I_DWTHAAR_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, 1);
  printf("Restored 2D signal:\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }

  for (j=0;j<NPTS;j++) {  // Make a fancy cubic signal
//     for (i=0;i<NPTS;i++) xy[j][i] = (3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j);
//     for (i=0;i<NPTS;i++) xy[j][i] = 5+j+0.4*j*j-0.02*j*j*j;
    for (i=0;i<NPTS;i++) xy[j][i] = sqrt((i-7.45)*(i-7.45) + (j-7.55)*(j-7.55));
  }
  printf("Original 2D signal (odd dimensions):\n");
  for (j=NPTSM-1;j>=0;j--) {
    for (i=0;i<NPTSM;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
  F_DWTHAAR_2D_split_inplace_n((float *)xy, NPTSM,  NPTS, NPTSM, 1);
//   F_DWTHAAR_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   F_DWTHAAR_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   F_DWTHAAR_2D_split_inplace((float *)xy, npts4, NPTS, npts4);

  printf("Transformed 2D signal:\n");
  for (j=NPTSM-1;j>=0;j--) {
    for (i=0;i<NPTSM;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }

  I_DWTHAAR_2D_split_inplace_n((float *)xy, NPTSM,  NPTS, NPTSM, 1);
  printf("Restored 2D signal:\n");
  for (j=NPTSM-1;j>=0;j--) {
    for (i=0;i<NPTSM;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
#if 0

  quantum = .025;
  printf("quantum used = %8.3f\n",quantum);
  for (j=0;j<NPTS;j++) {               // quantification pass
    for (i=0;i<NPTS;i++) { k = xy[j][i] / quantum + .5f ; xy[j][i] = k * quantum ; }
  }
  I_DWTHAAR_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, 3);
//   I_DWTHAAR_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
//   I_DWTHAAR_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   I_DWTHAAR_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   printf("Restored 2D signal:\n");
//   for (j=NPTS-1;j>=0;j--) {
//     for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
//     printf("\n");
//   }
  printf("Restored 2D signal error:\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i] - ((3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j)));
    printf("\n");
  }
#endif
}
#endif
