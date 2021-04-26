//  useful routines for C and FORTRAN programming
//  Copyright (C) 2020  Environnement Canada
//
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
//   interface                                                   !InTf!
//     subroutine AnalyzeField(a, ni, lni, nj, small, stats, str) bind(C,name='AnalyzeField')   !InTf!
//       import :: C_INTPTR_T, C_INT, C_FLOAT, C_CHAR            !InTf!
//       integer(C_INTPTR_T), intent(IN), value :: a             !InTf!
//       integer(C_INT), intent(IN), value :: ni, lni, nj        !InTf!
//       real(C_FLOAT), intent(IN), value :: small               !InTf!
//       real(C_FLOAT), intent(OUT), dimension(*) :: stats       !InTf!
//       character(C_CHAR), dimension(*) :: str                  !InTf!
//     end subroutine AnalyzeField                               !InTf!
//   end interface                                               !InTf!
void AnalyzeField(float *fa, int ni, int lni, int nj, float small, float *stats, char *str){
  int i, j, jm, jp ;
//   int im, ip ;
  double ci, cj, grad, gradbar, t, ddi, ddj ; // for gradient computation
  double sum, sum2;
  float fmin, fmax ;
  

// compute gradient grad(fa) = sqrt( [d(fa)/di]**2 +  [d(fa)/dj]**2 )
// compute maximum gradient and average gradient
  grad    = 0.0 ; 
  gradbar = 0.0 ;
  sum     = 0.0 ;
  sum2    = 0.0 ;
  fmin    = fa[0] ;
  fmax    = fa[0] ;
  
  for(j = 0 ; j < nj ; j++){
    jm = (j > 0) ? -lni : 0 ;                // row below current row (or row itself if row 0)
    jp = (j < nj-1) ? lni : 0 ;              // row above current row (or row itself if top row)
    cj = ((jp - jm) > 1 ) ? 0.5f : 1.0f ;    // scaling factor (.5 if centered difference, 1.0 otherwise)
//     im = (i > 0) ? i-1 : i ;
//     ip = (i < ni-1) ? i+1 : i ;
    ddi = fa[1] - fa[0] ; ddj = cj * (fa[0+jp] - fa[0+jm]) ;  // first point, i forward difference
    t = sqrt(ddi*ddi + ddj*ddj) ; grad = (t > grad) ? t : grad ; gradbar += t;
    sum += fa[0] ; sum2 += (fa[0] * fa[0]) ;
    fmin = (fa[0] < fmin) ? fa[0] : fmin ; fmax = (fa[0] > fmax) ? fa[0] : fmax ;
    for(i = 1 ; i < ni-1 ; i++){
      ddi = (fa[i+1] - fa[i-1]) * 0.5 ; ddj = cj * (fa[i+jp] - fa[i+jm]) ; 
      t = sqrt(ddi*ddi + ddj*ddj) ; grad = (t > grad) ? t : grad ; gradbar += t;
      sum += fa[i] ; sum2 += (fa[i] * fa[i]) ;
      fmin = (fa[i] < fmin) ? fa[i] : fmin ; fmax = (fa[i] > fmax) ? fa[i] : fmax ;
    }
    ddi = fa[ni-1] - fa[ni-2] ; ddj = cj * (fa[ni-1+jp] - fa[ni-1+jm]) ;  // last point, i backward difference
    t = sqrt(ddi*ddi + ddj*ddj) ; grad = (t > grad) ? t : grad ; gradbar += t;
    sum += fa[ni-1] ; sum2 += (fa[ni-1] * fa[ni-1]) ;
    fmin = (fa[ni-1] < fmin) ? fa[ni-1] : fmin ; fmax = (fa[ni-1] > fmax) ? fa[ni-1] : fmax ;
    fa += lni ;   // next base row
  }
  i = 0;
  stats[i++] = fmin ;            // min
  stats[i++] = fmax ;            // max
  stats[i++] = sum / (ni*nj) ;   // average
  stats[i++] = sum2 ;            // sum of squares
  stats[i++] = grad ;            // largest gradient
  gradbar /= (ni*nj) ;
  stats[i++] = gradbar ;         // average gradient
  printf("[%4d %4d %6s] min/max/avg/rng = %9.3g %9.3g %9.3g %9.3g, max/avg gradients = %9.3g %9.3g [%6.3f %6.3f]\n",
         ni, nj, str, fmin, fmax, sum/(ni*nj), fmax-fmin, grad, gradbar, grad/(fmax-fmin), gradbar/(fmax-fmin)    );
}

//   interface                                                   !InTf!
//     subroutine AnalyzeFields(a, b, ni, lni, nj, small, stats) bind(C,name='AnalyzeFields')   !InTf!
//       import :: C_INTPTR_T, C_INT, C_FLOAT                    !InTf!
//       integer(C_INTPTR_T), intent(IN), value :: a, b          !InTf!
//       integer(C_INT), intent(IN), value :: ni, lni, nj        !InTf!
//       real(C_FLOAT), intent(IN), value :: small               !InTf!
//       real(C_FLOAT), intent(OUT), dimension(*) :: stats       !InTf!
//     end subroutine AnalyzeFields                              !InTf!
//   end interface                                               !InTf!
void AnalyzeFields(float *fa, float *fb, int ni, int lni, int nj, float small, float *stats){
  int i, j, im, ip, jm, jp ;
  double ci, cj, grada, gradb, gradbara, gradbarb, t, ddi, ddj ;

// compute gradients grad(fa) = sqrt( [d(fa)/di]**2 +  [d(fa)/dj]**2 )
// compute maximum gradient and average gradient for fields fa and fb
  grada = 0.0f ; gradbara = 0.0 ;
  gradb = 0.0f ; gradbarb = 0.0 ;
  for(j = 0 ; j < nj ; j++){
    jm = (j > 0) ? -lni : 0 ;                // row below current row
    jp = (j < nj-1) ? lni : 0 ;              // row above current row
    cj = ((jp - jm) > 1 ) ? 0.5f : 1.0f ;    // scaling factor (.5 if centered difference, 1.0 otherwise)
    im = (i > 0) ? i-1 : i ;
    ip = (i < ni-1) ? i+1 : i ;
    ddi = fa[1] - fa[0] ; ddj = cj * (fa[0+jp] - fa[0+jm]) ;  // first point, i forward difference
    t = sqrt(ddi*ddi + ddj*ddj) ; grada = (t > grada) ? t : grada ; gradbara += t;
    ddi = fb[1] - fb[0] ; ddj = cj * (fb[0+jp] - fb[0+jm]) ;  // first point, i forward difference
    t = sqrt(ddi*ddi + ddj*ddj) ; gradb = (t > gradb) ? t : gradb ; gradbarb += t;
    for(i = 1 ; i < ni-1 ; i++){
      ddi = (fa[i+1] - fa[i-1]) * 0.5 ; ddj = cj * (fa[i+jp] - fa[i+jm]) ; 
      t = sqrt(ddi*ddi + ddj*ddj) ; grada = (t > grada) ? t : grada ; gradbara += t;
      ddi = (fb[i+1] - fb[i-1]) * 0.5 ; ddj = cj * (fb[i+jp] - fb[i+jm]) ; 
      t = sqrt(ddi*ddi + ddj*ddj) ; gradb = (t > gradb) ? t : gradb ; gradbarb += t;
    }
    ddi = fa[ni-1] - fa[ni-2] ; ddj = cj * (fa[ni-1+jp] - fa[ni-1+jm]) ;  // last point, i backward difference
    t = sqrt(ddi*ddi + ddj*ddj) ; grada = (t > grada) ? t : grada ; gradbara += t;
    ddi = fb[ni-1] - fb[ni-2] ; ddj = cj * (fb[ni-1+jp] - fb[ni-1+jm]) ;  // last point, i backward difference
    t = sqrt(ddi*ddi + ddj*ddj) ; gradb = (t > gradb) ? t : gradb ; gradbarb += t;
    fa += lni ;   // next base row
    fb += lni ;   // next base row
  }
  stats[0] = grada ;
  stats[1] = gradb ;
  stats[2] = gradbara/(ni*nj) ;
  stats[3] = gradbarb/(ni*nj) ;
  printf("max/avg gradients a = %8.3g %8.3g, max/avg gradients b = %8.3g %8.3g\n",grada, gradbara/(ni*nj), gradb, gradbarb/(ni*nj));
}

//   interface                                                   !InTf!
//     subroutine CompareFields(a, b, n, small, str) bind(C,name='CompareFields')   !InTf!
//       import :: C_INTPTR_T, C_INT, C_FLOAT, C_CHAR            !InTf!
//       integer(C_INTPTR_T), intent(IN), value :: a, b          !InTf!
//       integer(C_INT), intent(IN), value :: n                  !InTf!
//       real(C_FLOAT), intent(IN), value :: small               !InTf!
//       character(C_CHAR), dimension(*) :: str                  !InTf!
//     end subroutine CompareFields                              !InTf!
//   end interface                                               !InTf!
// quick and dirty error analysis
void CompareFields(float *fa, float *fb, int np, float small, char *str){  // will have to add a few options
  int i;
  float maxval, minval, relerr, rdiff, snr;
  double err, errmax, errsum, errsuma, sum2, acc2, acc0;
  double suma, sumb, suma2, sumb2, sumab;
  double vara, varb, varab, avga, avgb, rab;
  double ssim;  // structural similarity
  uint32_t *ia = (uint32_t *)fa;
  uint32_t *ib = (uint32_t *)fb;
  uint32_t ierr = 0;
  int32_t idiff;
  int indx, iacc, iabs;
  int accuracy, n;
  uint64_t idif64;

  small  = fabsf(small); // in case small is negative
  err    = 0.0f;
  sum2   = 0.0;      // sum of errors**2 (for RMS)
  suma2  = 0.0;      // sum of squares, array fa
  sumb2  = 0.0;      // sum of squares, array fb
  suma   = 0.0;      // sum of terms, array fa
  sumb   = 0.0;      // sum of terms, array fb
  sumab  = 0.0;      // sum of fa*fb products (for covariance)
  indx   = 0;        // position of largest relative error
  iacc   = 0;        // position of largest bit inaccuracy
  iabs   = 0;        // position of largest absolute error
  ierr   = 0;        // largest bit inaccuracy
  idif64 = 0;        // sum of bit inaccuracies
  relerr = 0.0f;     // largest relative error
  errmax = 0.0;     // largest absolute error
  errsum = 0.0;     // sum of errors (for BIAS)
  errsuma = 0.0;    // sum of absolute errors
  maxval = fa[0];    // highest signed value in array fa
  minval = fa[0];    // lowest signed value in array fa
  n      = 0;
  for(i=0 ; i < np ; i++){
    suma  += fa[i] ;
    suma2 += ( fa[i] * fa[i] ) ;
    sumb  += fb[i] ;
    sumb2 += ( fb[i] * fb[i] ) ;
    sumab += ( fb[i] * fa[i] ) ;
    maxval = (fa[i] > maxval ) ? fa[i] : maxval ;
    minval = (fa[i] < minval ) ? fa[i] : minval ;
//     err = 1.0*fb[i] - 1.0*fa[i] ;               // signed error
    err = fb[i] ;
    err -= fa[i] ;
    sum2 = sum2 + err * err ;           // sum of squared error (to compute RMS)
    errsum += err ;                     // sum of signed errors (to compute BIAS)
    err = fabs(err) ;                   // absolute error
    errsuma += err ;                    // sum of absolute errors
    if(err > errmax) {errmax = err ; iabs = i; } ;
    if(fabsf(fa[i]) <= small) continue ;  // ignore absolute values smaller than threshold
    if(fa[i] * fb[i] < 0.0f ) continue ;  // opposite signs, ignore
    n++;
    rdiff = fabsf(fa[i] - fb[i]) ;
    rdiff = rdiff / fa[i] ;              // fa[i] should never be zero at this point
    if(rdiff > relerr) { relerr = rdiff; indx = i ; }   // largest relative error
    idiff = (ia[i] > ib[i]) ? (ia[i] - ib[i]) : (ib[i] - ia[i]) ;
    idif64 += idiff;
    if(idiff > ierr) { ierr = idiff ; iacc = i ; }
  }
  if(n == 0) {
//     printf(">>%g>>>",small);
    n = 1 ;
  }
  avga  = suma/np;
  vara  = suma2/np - avga*avga ; // variance of a
  avgb  = sumb / np;
  varb  = sumb2/np - avgb*avgb ; // variance of b
  varab = sumab/np - avga*avgb ; // covariance of a and b
  ssim  = (2.0 * avga * avgb + .01) * (2.0 * varab + .03) / ((avga*avga + avgb*avgb + .01)*(vara + varb + .03)) ; // structural similarity
  rab   = (np*sumab - suma*sumb) / (sqrt(np*suma2 - suma*suma) * sqrt(np*sumb2 - sumb*sumb)) ;  // Pearson coefficient
  idif64 = idif64/n;                                  // average ULP difference
  acc2 = log2(1.0+idif64);                     // accuracy (in agreed upon bits)
  if(acc2 < 0) acc2 = 0;
  if(acc2 > 24) acc2 = 24;
  acc0 = log2(1.0+ierr);                       // worst accuracy
  if(acc0 < 0) acc0 = 0;
  if(acc0 > 24) acc0 = 24;
  sum2 = sum2 / np;                                    // average quadratic error
  snr  = .25 * (maxval - minval) * (maxval - minval);
  snr  = 10.0 * log10(snr / sum2);                    // Peak Signal / Noise Ratio
  if(snr < 0.0) snr = 0.0 ;
//   if(relerr < .000001f) relerr = .000001f;
  if(relerr == 0.0) relerr = 1.0E-10;
  printf("%s[%6.4f %6d] ",str,small,n);
  printf("max/avg/bias/rms/rel err (%8.6f, %8.6f, %8.6f, %8.6f, 1/%8.2g), range = %10.6f",errmax, errsuma/n, errsum/n, sqrt(sum2),1.0/relerr, maxval-minval);
  printf(", worst/avg accuracy (%6.2f,%6.2f)b, PSNR = %5.0f",24-acc0, 24-acc2, snr);  // probably not relevant for FCST verif
//   printf(", DISSIM = %12.6g, Pearson = %12.6g", 1.0 - ssim,1.0-rab);   // not very useful for packing error analysis, maybe for FCST verif ?
  printf(" [%.0f]\n",(maxval-minval)/errmax);
  accuracy = 0;
  while(ierr >>= 1) accuracy ++;
//   printf("max error  = %8.6f, bias = %10.6f, avg error = %8.6f, rms = %8.6f, npts = %d, min,max = (%10.6f %10.6f), ",
//      errmax, errsum/(n), errsuma/(n), sqrt(sum2), n, minval, maxval);
//   printf("range/errormax = %g\n",(maxval-minval)/errmax);
//   printf("accuracy(bits) = %9d [%8d](%8.3g,%8.3g)(%9.6f)\n",accuracy,iacc, fa[iacc], fb[iacc],fb[iacc] - fa[iacc]);
//   printf("max rel error  = %9.5g [%8d](%8.3g,%8.3g)(%9.6f)\n",relerr,indx, fa[indx], fb[indx],fb[indx] - fa[indx]);
//   printf("max abs error  = %9.6f [%8d](%8.3g,%8.3g)(%9.6f)\n",errmax,iabs, fa[iabs], fb[iabs],fb[iabs] - fa[iabs]);
//   accuracy = 24;
//   while(idif64 >>= 1) accuracy --;
//   printf("average accuracy bits = %6.2f, PSNR = %f\n",acc2, snr);
}
