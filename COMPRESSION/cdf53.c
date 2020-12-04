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
// Cohen-Daubechies-Favreau 5/3 wavelets
// https://en.wikipedia.org/wiki/Cohen%E2%80%93Daubechies%E2%80%93Feauveau_wavelet
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
#include <cdf53.h>

#define A      (-0.5f)
#define B      0.25f
#define S      1.414214f
#define Z      0.7071068f

//****f* librkl/F_CDF53_1D_split_N_even
// Synopsis
//
// Forward DWT transform (analysis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_CDF53_1D_split_N_even(x, e, o, n) bind(C,name='F_CDF53_1D_split_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_CDF53_1D_split_N_even                 !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_split_N_even(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = neven;

  for(i = 0 ; i < nodd-1 ; i++) o[i] = x[i+i+1] + A * (x[i+i] + x[i+i+2]);  // predict odd terms
  o[nodd-1] = x[n-1] + 2 * A * x[n-2];  

  e[0 ] = x[0] + 2 * B * o[0];
  for(i = 1; i < neven ; i++) e[i] = x[i+i] + B * (o[i] + o[i-1]);          // update even terms

  e[0] = S * (e[0]);                                                        // scale
  for(i = 1; i < neven ; i++) { e[i] = S * (e[i]) ; o[i-1] *= (-Z); }
  o[nodd-1] *= (-Z);
}

//****f* librkl/F_CDF53_1D_inplace_N_even
// Synopsis
//
// Forward DWT transform (analysis) (in place, x is overwritten)
// n           : number of data points (even)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_CDF53_1D_inplace_N_even(x, e, o, n) bind(C,name='F_CDF53_1D_inplace_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF53_1D_inplace_N_even               !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_inplace_N_even(float *x, int n){    // InTc
//****
  int i;
  
  for (i = 1; i < n - 2; i += 2) x[i] += A * (x[i-1] + x[i+1]);  // predict odd terms
  x[n-1] += 2 * A * x[n-2];                                      // last term is odd
  
  x[0] += 2 * B * x[1];                                          // update even terms
  for (i = 2; i < n; i += 2) x[i] += B * (x[i+1] + x[i-1]);
  
  x[0] = S * (x[0]);                                             // scale
  for (i = 2; i < n; i += 2) { x[i] = S * (x[i]) ; x[i-1] *= (-Z); }
  x[n-1] *= (-Z);                                                // scale last (odd) term
}

//****f* librkl/F_CDF53_1D_split_N_odd
// Synopsis
//
// Forward DWT transform (analysis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_CDF53_1D_split_N_odd(x, e, o, n) bind(C,name='F_CDF53_1D_split_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_CDF53_1D_split_N_odd                  !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_split_N_odd(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;

  for(i = 0 ; i < nodd ; i++) o[i] = x[i+i+1] + A * (x[i+i] + x[i+i+2]);  // predict odd terms

  e[0 ] = x[0] + 2 * B * o[0];
  for(i = 1; i < neven-1 ; i++) e[i] = x[i+i] + B * (o[i] + o[i-1]);      // update even terms
  e[neven-1] = x[n-1] + 2 * B * o[nodd-1];

  e[0] = S * (e[0]);                                                      // scale
  for(i = 1; i < neven-1 ; i++) { e[i] = S * (e[i]) ; o[i-1] *= (-Z); }
  e[neven-1] = S * (e[neven-1]);
  o[nodd-1] *= (-Z);
}

//****f* librkl/F_CDF53_1D_inplace_N_odd
// Synopsis
//
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_CDF53_1D_inplace_N_odd(x, e, o, n) bind(C,name='F_CDF53_1D_inplace_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF53_1D_inplace_N_odd                !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_inplace_N_odd(float *x, int n){    // InTc
//****
  int i;

  for (i = 1; i < n - 1; i += 2) x[i] += A * (x[i-1] + x[i+1]);  // predict odd terms
  
  x[0] += 2 * B * x[1];                                          // update even terms
  for (i = 2; i < n - 2; i += 2) x[i] += B * (x[i+1] + x[i-1]);
  x[n - 1] += 2 * B * x[n - 2];                                  // last term is even
  
  x[0] = S * (x[0]);                                             // scale
  for (i = 2; i < n - 2; i += 2) { x[i] = S * (x[i]);  x[i-1] *= (-Z); }
  x[n - 1] = S * (x[n - 1]);                  // last term is even
  x[n - 2] *= (-Z);                                             // scale last odd term
  }
  
//****f* librkl/F_CDF53_1D_split
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
//   subroutine F_CDF53_1D_split(x, e, o, n) bind(C,name='F_CDF53_1D_split')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_CDF53_1D_split                        !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_split(float *x, float *e, float *o, int n){    // InTc
//****
  if(n < 3) {
    if(n > 0) e[0] = x[0];
    if(n > 1) o[0] = x[1];
    return;
  }
  if(n & 1){
    F_CDF53_1D_split_N_odd(x, e, o, n);
  }else{
    F_CDF53_1D_split_N_even(x, e, o, n);
  }
}
  
//****f* librkl/F_CDF53_1D_inplace
// Synopsis
//
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_CDF53_1D_inplace(x, e, o, n) bind(C,name='F_CDF53_1D_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF53_1D_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_inplace(float *x, int n){    // InTc
//****
  if(n < 3) return;
  if(n & 1){
    F_CDF53_1D_inplace_N_odd(x, n);
  }else{
    F_CDF53_1D_inplace_N_even(x, n);
  }
}
  
//****f* librkl/F_CDF53_1D_split_inplace
// Synopsis
//
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine F_CDF53_1D_split_inplace(x, e, o, n) bind(C,name='F_CDF53_1D_split_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF53_1D_split_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_1D_split_inplace(float *x, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int lasteven = neven - 1;
  int nodd  = n >> 1;
  int lastodd = nodd - 1;
  float odd[n];

  if(n < 3) return;  
  for (i = 0; i < nodd; i ++) { odd[i] = x[i+i+1] ; x[i] = x[i+i] ; }  // copy odd terms to temporary odd array
  x[lasteven] = x[lasteven+lasteven]; x[neven] = x[lasteven] ;         // even terms to beginning of array

  for (i = 0; i < nodd; i ++) odd[i] += A * (x[i] + x[i+1]) ;          // predict odd terms
  odd[nodd] = odd[lastodd] ;                                           // pad end of odd to make it at least neven terms

  x[0] += 2 * B * odd[0] ;                                             // update even terms
  for (i = 1; i < neven; i ++) x[i] += B * (odd[i] + odd[i-1]);
  x[neven] = x[lasteven] ;                                             // pad end of even to make it at least nodd + 1 terms

  x[0] = S * (x[0]) ;                                                  // scale
  for (i = 1; i < neven; i ++) x[i] = S * ( x[i]);
  for (i = 0; i < nodd; i ++) x[neven+i] = (-Z) * odd[i] ;
}

//****f* librkl/I_CDF53_1D_split_N_even
// Synopsis
//
// Inverse DWT transform (synthesis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_CDF53_1D_split_N_even(x, e, o, n) bind(C,name='I_CDF53_1D_split_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_CDF53_1D_split_N_even                 !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_split_N_even(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = neven;

  for(i = 0 ; i < neven ; i++){ x[i+i] = e[i]*Z ; x[i+i+1] = o[i] * (-S) ; }  // unscale and move to x

  x[0] -= 2 * B * x[1];
  for (i = 2; i < n; i += 2) x[i] -= B * (x[i+1] + x[i-1]);         // unupdate even terms

  for (i = 1; i < n - 2; i += 2) x[i] -= A * (x[i-1] + x[i+1]);     // unpredict odd terms
  x[n - 1] -= 2 * A * x[n - 2];
}

//****f* librkl/I_CDF53_1D_inplace_N_even
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_CDF53_1D_inplace_N_even(x, e, o, n) bind(C,name='I_CDF53_1D_inplace_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF53_1D_inplace_N_even               !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_inplace_N_even(float *x, int n){    // InTc
//****
  int i;
  
  for(i = 0 ; i < n-1 ; i+=2 ) { x[i] *= Z ; x[i+1] *= (-S) ; }     // unscale
  
  x[0] -= 2 * B * x[1];                                             // unupdate even terms
  for (i = 2; i < n; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  
  for (i = 1; i < n - 2; i += 2)  x[i] -= A * (x[i-1] + x[i+1]);    // unpredict odd terms
  x[n - 1] -= 2 * A * x[n - 2];
}

//****f* librkl/I_CDF53_1D_split_N_odd
// Synopsis
//
// Inverse DWT transform (synthesis)
// n           : number of data points (odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_CDF53_1D_split_N_odd(x, e, o, n) bind(C,name='I_CDF53_1D_split_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_CDF53_1D_split_N_odd                  !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_split_N_odd(float *x, float *e, float *o, int n){    // InTc
//****
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;

  for(i = 0 ; i < nodd ; i++){ x[i+i] = e[i] * Z ; x[i+i+1] = o[i] * (-S) ; }  // unscale and move to x
  x[n-1] = e[neven-1] * Z;

  x[0] -= 2 * B * x[1];                                             // unupdate even terms
  for (i = 2; i < n - 2; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  x[n - 1] -= 2 * B * x[n - 2];

  for (i = 1; i < n - 1; i += 2) x[i] -= A * (x[i-1] + x[i+1]);     // unpredict odd terms
}

//****f* librkl/I_CDF53_1D_inplace_N_odd
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_CDF53_1D_inplace_N_odd(x, e, o, n) bind(C,name='I_CDF53_1D_inplace_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF53_1D_inplace_N_odd                !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_inplace_N_odd(float *x, int n){    // InTc
//****
  int i;
  
  for(i = 0 ; i < n-2 ; i+=2 ) { x[i] *= Z ; x[i+1] *= (-S) ; }     // unscale odd and even terms
  x[n-1] *= Z;

  x[0] -= 2 * B * x[1];                                             // unupdate even terms
  for (i = 2; i < n - 2; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  x[n - 1] -= 2 * B * x[n - 2];

  for (i = 1; i < n - 1; i += 2) x[i] -= A * (x[i-1] + x[i+1]);     // unpredict odd terms
}

//****f* librkl/I_CDF53_1D_split
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
//   subroutine I_CDF53_1D_split(x, e, o, n) bind(C,name='I_CDF53_1D_split')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_CDF53_1D_split                        !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_split(float *x, float *e, float *o, int n){    // InTc
//****
  if(n < 3) {   // 3 points minimum
    if(n > 0) x[0] = e[0] ;
    if(n > 1) x[1] = o[0] ;
    return;
  }
  if(n & 1){
    I_CDF53_1D_split_N_odd(x, e, o, n);
  }else{
    I_CDF53_1D_split_N_even(x, e, o, n);
  }
}

//****f* librkl/I_CDF53_1D_split_inplace
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_CDF53_1D_split_inplace(x, e, o, n) bind(C,name='I_CDF53_1D_split_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF53_1D_split_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_split_inplace(float *x, int n){    // InTc
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
    I_CDF53_1D_inplace_N_odd(x, n);
  }else{
    I_CDF53_1D_inplace_N_even(x, n);
  }
}

//****f* librkl/I_CDF53_1D_inplace
// Synopsis
//
// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
//
// FORTRAN interface
// interface        !InTf
//   subroutine I_CDF53_1D_inplace(x, e, o, n) bind(C,name='I_CDF53_1D_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF53_1D_inplace                      !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_1D_inplace(float *x, int n){    // InTc
//****
  if(n < 3) return;   // 3 points minimum
  if(n & 1){
    I_CDF53_1D_inplace_N_odd(x, n);
  }else{
    I_CDF53_1D_inplace_N_even(x, n);
  }
}

#define VCONTRIB(DEST,SCALE,SRC1,SRC2,N) { int i; for(i=0 ; i<N ; i++) {DEST[i] += SCALE *(SRC1[i] + SRC2[i]) ; } }
#define VSCALE(WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHAT[i] *= FAC ; } } 
#define VSCALE2(WHERE,WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHERE[i] = WHAT[i] *FAC ; } } 

//****f* librkl/I_CDF53_2D_split_inplace
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
// subroutine I_CDF53_2D_split_inplace(x, ni, lni, nj) BIND(C,name='I_CDF53_2D_split_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine I_CDF53_2D_split_inplace                        !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_2D_split_inplace(float *x, int ni, int lni, int nj){    // InTc
//****
  int neven = (nj+1) >> 1;
  int ieven = (ni+1) >> 1;
  int nodd = nj >> 1;
  int j;
  float *rowd, *row0, *rows, *row1, *row2;
  int lni2 = lni+lni;
  float t[(nj*ni) >> 1];
  float *rowt, *oddrow;;
  // unscale
  rows = x + neven * lni ;                  // first odd row in x
  rowt = t ;                                // first temporary odd row
  row1 = x ; row2 = row1 ;                  // first even row
  VSCALE(row1,(Z),ni) ;                     // unscale first even row
  VSCALE2(rowt,rows,(-S),ni) ;              // unscale first odd row
  for(j = 1 ; j < nodd ; j++){
    row1 = row2 ; row2 = row1 + lni ;       // next even row pair in x
    VSCALE(row2,(Z),ni) ;                   // unscale next even row
    rows += lni ; rowt += ni ;              // next odd row (x and temporary)
    row0 = rowt - ni;                       // previous (unscaled) odd row
    VSCALE2(rowt,rows,(-S),ni) ;            // unscale next odd row
  }
  if(nodd == neven){                         // last row is odd
  }else{                                     // last row is even
    row1 = row2 ; row2 = row1 + lni ;        // next even row pair in x
    VSCALE(row2,(Z),ni) ;                    // unscale last even row
  }
  // combined unupdate, unpredict
  rowt = t ;                                 // first temporary odd row
  row1 = x ; row2 = row1 ;                   // first even row
  VCONTRIB(row1, (-B), rowt, rowt, ni);      // unupdate first even row
  for(j = 1 ; j < nodd ; j++){
    row1 = row2 ; row2 = row1 + lni ;        // next even row pair in x
    rowt += ni ; row0 = rowt - ni;           // current and previous odd rows
    VCONTRIB(row2, (-B), row0, rowt, ni);    // un update even row
    VCONTRIB(row0, (-A), row1, row2, ni);    // unpredict odd row below 
  }
  if(nodd == neven){                         // nj even, last row is odd
    VCONTRIB(rowt, (-A), row2, row2, ni);    // unpredict last odd row
  }else{                                     // last row is even, unupdate it, then unpredict odd row below
    row1 = row2 ; row2 = row1 + lni ;        // next even row pair in x
    VCONTRIB(row2, (-B), rowt, rowt, ni);    // unupdate even row
    VCONTRIB(rowt, (-A), row1, row2, ni);
  }
  // last pass, 1D split inverse transform
  rows = x + (nj -1 ) * lni ;                // last row in x
  rowt = t + (nodd - 1) * ni ;               // last odd row in t
  rowd = x + (neven - 1) * lni ;             // last even row in split x
  for(j = nj - 1 ; j > 0 ; j--){             // last pass for 1D transform
    if(j & 1){                               // odd row, from rowt to rows
      I_CDF53_1D_split(rows, rowt, rowt + ieven, ni) ;
      rowt -= ni ;
    }else{                                   // even row, from rowd to rows
      I_CDF53_1D_split(rows, rowd, rowd + ieven, ni);
      rowd -= lni;
    }
    rows -= lni;
  }
  I_CDF53_1D_split_inplace(x, ni);           // first even row is in place
}

//****f* librkl/F_CDF53_2D_split_inplace
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
// subroutine F_CDF53_2D_split_inplace(x, ni, lni, nj) BIND(C,name='F_CDF53_2D_split_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine F_CDF53_2D_split_inplace                  !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_2D_split_inplace(float *x, int ni, int lni, int nj){    // InTc
//****
  int neven = (nj+1) >> 1;
  int ieven = (ni+1) >> 1;
  int nodd = nj >> 1;
  int i, j;
  float *rowd, *row0, *row1, *row2;
  int lni2 = lni+lni;
  float t[(nj*ni) >> 1];
  float *rowt, *rowx, *rows, *top ;

  top = x + nj * lni ;                // top of array (one row above last, rows must be < top)
  rowx = x;                           // first even row storage
  rowt = t;                           // first odd row temporary storage
  // 1 D transform for the first 4 rows
  rows = x; F_CDF53_1D_split_inplace(rowx, ni); rowx += lni ;  // split, inplace 1D transform (first even row)
  rows += lni ; if(rows < top) F_CDF53_1D_split(rows, rowt, rowt + ieven, ni) ; rowt += ni ;  // first odd row
  rows += lni ; if(rows < top) F_CDF53_1D_split(rows, rowx, rowx + ieven, ni) ; rowx += lni ; // next even row
  rows += lni ; if(rows < top) F_CDF53_1D_split(rows, rowt, rowt + ieven, ni) ; rowt += ni ;  // next odd row
  // combined pass : predict odd rows, update even rows
  // perform the 1D transform on the fly in the first pass before row is needed
  rowd = t ; row1 = x ; row2 = row1 + lni;
  VCONTRIB(rowd, A, row1, row2, ni) ;         // predict first odd row
  VCONTRIB(row1, B, rowd, rowd, ni) ;         // update first even row
  for(j = 1 ; j < neven-1 ; j++){
    rows += lni ; if(rows < top) F_CDF53_1D_split(rows, rowx, rowx + ieven, ni) ; rowx += lni ; // next even row
    rows += lni ; if(rows < top) F_CDF53_1D_split(rows, rowt, rowt + ieven, ni) ; rowt += ni ;  // next odd row
    rowd += ni ;                              // next odd row
    row1 = row2 ; row2 = row1 + lni ;         // next pair of even rows
    VCONTRIB(rowd, A, row1, row2, ni) ;       // predict odd row
    row0 = rowd - ni ; 
    VCONTRIB(row1, B, rowd, row0, ni) ;       // update even row below odd row
  }
  if(nodd == neven){                          // nj even
    rowd += ni ;                              // next (last) odd row
    VCONTRIB(rowd, A, row2, row2, ni);        // last row is odd, predict it using last even row
    row0 = rowd-ni ; 
    VCONTRIB(row2, B, rowd, row0, ni);        // update even row below last odd row
  }else{                                      // nj odd, 
    VCONTRIB(row2, B, rowd, rowd, ni);        // last row is even, update it using last odd row
  }
// combined pass : scaling and move back into original array
  rowd = t ; row1 = x ; row2 = row1 + lni;
  VSCALE(row1,(S),ni);                        // scale first even row
  rows = x + neven * lni ;                    // first odd row in split x array 
  for(j = 1 ; j < neven-1 ; j++){
    rowd += ni ;                              // next odd row
    row1 = row2 ; row2 = row1 + lni ;         // next pair of even rows
    row0 = rowd - ni ; 
    VSCALE2(rows,row0,(-Z),ni);               // scale low odd row and transfer it to upper part of x
    VSCALE(row1,(S),ni) ;                     // scale lower even row
    rows += lni ;
  }
  if(nodd == neven){                          // nj even
    rowd += ni ;                              // next odd row
    row0 = rowd - ni ; 
    VSCALE(row2,(S),ni);                      // scale last even row
    VSCALE2(rows,row0,(-Z),ni);               // scale last 2 odd rows
    rows += lni ;
    VSCALE2(rows,rowd,(-Z),ni);
  }else{                                      // nj odd, 
    VSCALE2(rows,rowd,(-Z),ni);               // scale last odd row
    VSCALE(row2,(S),ni);                      // scale last even row
  }
}

//****f* librkl/I_CDF53_2D_split_inplace_n
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
// subroutine I_CDF53_2D_split_inplace_n(x, ni, lni, nj, levels) BIND(C,name='I_CDF53_2D_split_inplace_n') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni, levels    !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine I_CDF53_2D_split_inplace_n                !InTf
// end interface    !InTf
// ARGUMENTS
void I_CDF53_2D_split_inplace_n(float *x, int ni, int lni, int nj, int levels){    // InTc
//****
  if(levels > 1) {
    I_CDF53_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels - 1);
  }
  I_CDF53_2D_split_inplace(x, ni, lni, nj);
}

//****f* librkl/F_CDF53_2D_split_inplace_n
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
// subroutine F_CDF53_2D_split_inplace_n(x, ni, lni, nj, levels) BIND(C,name='F_CDF53_2D_split_inplace_n') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni, levels    !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine F_CDF53_2D_split_inplace_n                !InTf
// end interface    !InTf
// ARGUMENTS
void F_CDF53_2D_split_inplace_n(float *x, int ni, int lni, int nj, int levels){    // InTc
//****
  F_CDF53_2D_split_inplace(x, ni, lni, nj);
  if(levels > 1) {
    F_CDF53_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels - 1);
  }
}

#if defined(SELF_TEST)
#include <stdlib.h>
#include <math.h>

#if ! defined(NPTS)
#define NPTS 16
#endif

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
  for (j=0;j<NPTS;j++) {
  // Makes a fancy cubic signal
    for (i=0;i<NPTS;i++) xy[j][i] = (3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j);
//     for (i=0;i<NPTS;i++) xy[j][i] = sqrt((i-7.45)*(i-7.45) + (j-7.55)*(j-7.55));
  }
  
  printf("Original 2D signal:\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
  F_CDF53_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, 3);
//   F_CDF53_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   F_CDF53_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   F_CDF53_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
  
  printf("Transformed 2D signal (before quantification):\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
  quantum = .025;
  printf("quantum used = %8.3f\n",quantum);
  for (j=0;j<NPTS;j++) {               // quantification pass
    for (i=0;i<NPTS;i++) { k = xy[j][i] / quantum + .5f ; xy[j][i] = k * quantum ; }
  }
  I_CDF53_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, 3);
//   I_CDF53_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
//   I_CDF53_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   I_CDF53_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
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

}
#endif
