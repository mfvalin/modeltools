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
/* minimal code example showing how to call the zfp (de)compressor */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "zfp.h"

#if defined(DEBUG)
static int debug = DEBUG;
#else
static int debug = 0;
#endif

//   interface                                                   !InTf!
//     function ZfpCompress3D(z, nx, ny, nz, approx, stream, streamsize) result(addr) bind(C,name='ZfpCompress3D')  !InTf!
//       import :: C_INTPTR_T, C_INT, C_DOUBLE                   !InTf!
//       integer(C_INTPTR_T), intent(IN), value :: z, stream     !InTf!
//       integer(C_INT), intent(IN), value :: nx, ny, nz         !InTf!
//       integer(C_INT), intent(INOUT) :: streamsize             !InTf!
//       real(C_DOUBLE), intent(IN), value :: approx             !InTf!
//       integer(C_INTPTR_T) :: addr                             !InTf!
//     end function ZfpCompress3D                                !InTf!
//     function ZfpCompress3D_debug(flag) result(old) bind(C,name='ZfpCompress3D_debug')   !InTf!
//       import :: C_INT                                         !InTf!
//       integer(C_INT), intent(IN), value :: flag               !InTf!
//       integer(C_INT) :: old                                   !InTf!
//     end function ZfpCompress3D_debug                          !InTf!
//   end interface                                               !InTf!
// quick and dirty error analysis
void AnalyzeCompressionErrors(float *fa, float *fb, int np, float small, char *str){  // will have to add a few options
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
  errmax = 0.0f;     // largest absolute error
  errsum = 0.0f;     // sum of errors (for BIAS)
  errsuma = 0.0f;    // sum of absolute errors
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
    err = fb[i] - fa[i] ;               // signed error
    sum2 = sum2 + err * err ;           // sum of squared error (to compute RMS)
    errsum += err ;                     // sum of signed errors (to compute BIAS)
    err = fabs(err) ;                   // absolute error
    errsuma += err ;                    // sum of absolute errors
    if(err > errmax) {errmax = err ; iabs = i; } ;
    if(fabsf(fa[i]) <= small) continue ;  // ignore absolute values smaller than threshold
    if(fa[i] < 0.0f || fb[i] < 0.0f ) continue ;  // opposite signs, ignore
    n++;
    rdiff = fabsf(fa[i] - fb[i]) ;
    rdiff = rdiff / fa[i] ;              // fa[i] should never be zero at this point
    if(rdiff > relerr) { relerr = rdiff; indx = i ; }   // largest relative error
    idiff = (ia[i] > ib[i]) ? (ia[i] - ib[i]) : (ib[i] - ia[i]) ;
    idif64 += idiff;
    if(idiff > ierr) { ierr = idiff ; iacc = i ; }
  }
//   printf("np, n = %d %d\n",np,n);
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
  acc0 = log2(1.0+ierr);                       // worst accuracy
  if(acc0 < 0) acc0 = 0;
  sum2 = sum2 / np;                                    // average quadratic error
  snr  = .25 * (maxval - minval) * (maxval - minval);
  snr  = 10.0 * log10(snr / sum2);                    // Peak Signal / Noise Ratio
//   if(relerr < .000001f) relerr = .000001f;
  if(relerr == 0.0) relerr = 1.0E-10;
  printf("%s[%6.4f] ",str,small);
  printf("max/avg/bias/rms/rel err = (%8.6f, %8.6f, %8.6f, %8.6f, 1/%8.2g), range = %10.6f",errmax, errsuma/n, errsum/n, sqrt(sum2),1.0/relerr, maxval-minval);
  printf(", worst/avg accuracy = (%6.2f,%6.2f) bits, PSNR = %5.0f",24-acc0, 24-acc2, snr);  // probably not relevant for FCST verif
//   printf(", DISSIM = %12.6g, Pearson = %12.6g", 1.0 - ssim,1.0-rab);   // not very useful for packing error analysis, maybe for FCST verif ?
  printf(" [%.0f]\n",(maxval-minval)/errmax);
  accuracy = 0;
  while(ierr >>= 1) accuracy ++;
//   printf("max error  = %8.6f, bias = %10.6f, avg error = %8.6f, rms = %8.6f, npts = %d, min,max = (%10.6f %10.6f), ",
// 	errmax, errsum/(n), errsuma/(n), sqrt(sum2), n, minval, maxval);
//   printf("range/errormax = %g\n",(maxval-minval)/errmax);
//   printf("accuracy(bits) = %9d [%8d](%8.3g,%8.3g)(%9.6f)\n",accuracy,iacc, fa[iacc], fb[iacc],fb[iacc] - fa[iacc]);
//   printf("max rel error  = %9.5g [%8d](%8.3g,%8.3g)(%9.6f)\n",relerr,indx, fa[indx], fb[indx],fb[indx] - fa[indx]);
//   printf("max abs error  = %9.6f [%8d](%8.3g,%8.3g)(%9.6f)\n",errmax,iabs, fa[iabs], fb[iabs],fb[iabs] - fa[iabs]);
//   accuracy = 24;
//   while(idif64 >>= 1) accuracy --;
//   printf("average accuracy bits = %6.2f, PSNR = %f\n",acc2, snr);
}

int ZfpCompress3D_debug(int flag){ 
  int oldval = debug ;
  debug = flag ; 
  return oldval;
}

/* compress or decompress array */
void *ZfpCompress3D(float* array, int nx, int ny, int nz, double approx, void *Stream, int *Ssize)
{
  int status = 0;    /* return value: 0 = success */
  zfp_type type;     /* array scalar type */
  zfp_field* field;  /* array meta data */
  zfp_stream* zfp;   /* compressed stream */
  void* buffer;      /* storage for compressed stream */
  size_t bufsize;    /* byte size of compressed buffer */
  bitstream* stream; /* bit stream to write to or read from */
  size_t zfpsize;    /* byte size of compressed stream */
  double time, ratio;
  clock_t c;
  int i;
  int decompress;
  size_t BufSiz = *Ssize;
  uint32_t precision = 16;  // default precision
  double tolerance = 0.0;
  double rate = 0.0;

  /* allocate meta data for a compressed stream */
  zfp = zfp_stream_open(NULL);

  if(approx < 0) { precision = -approx; tolerance = 0.0 ; }  // precision mode
  if(approx > 0) { tolerance = approx;  precision = 0 ; }    // error bound mode
  decompress = (Stream != NULL) ;                            // compressed stream is NULL is compressing

  if(decompress){
    buffer = Stream ;
    bufsize = BufSiz ;
    field = zfp_field_alloc() ;
    zfp_field_set_pointer(field, array);
  }else{
    //  allocate meta data for the array
    type = zfp_type_float;                          // set field type and dimensions
    if(nz > 1 && debug)             fprintf(stderr, "3D array, ");
    if(nz <= 1 && ny > 1 && debug)  fprintf(stderr, "2D array, ");
    if(ny <= 1 && nz <= 1 && debug) fprintf(stderr, "1D array, ");
    if(nz > 1)             field = zfp_field_3d(array, type, nx, ny, nz);  // 3D array
    if(nz <= 1 && ny > 1)  field = zfp_field_2d(array, type, nx, ny);      // 2D array
    if(ny <= 1 && nz <= 1) field = zfp_field_1d(array, type, nx);          // 1D array
    field = zfp_field_3d(array, type, nx, ny, nz);  // set field dimensions
    //  set compression mode and parameters via one of three functions */
    if(precision > 0) {         // precision mode
      rate = -approx*.01;
      if(rate > 1.0) ratio = zfp_stream_set_rate(zfp, rate, type, (nz > 0) ? 3 : 2, 0);  // bit rate mode
      else           ratio = zfp_stream_set_precision(zfp, precision);                   // bit precision mode
      if(debug) fprintf(stderr, "precision in/out = %g, %g\n",-approx,ratio);
    }
    if(tolerance > 0) {         // error bound mode
      ratio = zfp_stream_set_accuracy(zfp, tolerance);
      if(debug) fprintf(stderr, "tolerance in/out = %g, %g\n",tolerance,ratio);
    }
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer = malloc(bufsize);                       // allocate buffer for compressed data 
  }

  stream = stream_open(buffer, bufsize);            // open stream
  zfp_stream_set_bit_stream(zfp, stream);           // associate bit stream with allocated buffer
  zfp_stream_rewind(zfp);

  if (decompress) {                                 // decompress array from stream buffer
    zfpsize = BufSiz;
    if(zfp_read_header(zfp, field, ZFP_HEADER_FULL) != 0){
      c = clock();
      if (!zfp_decompress(zfp, field)) {
	if(debug) fprintf(stderr, "decompression failed\n");
	status = 1;
      }
      time = (double)(clock() - c) / CLOCKS_PER_SEC;
      ratio = nx*ny*nz*4 ; ratio /= zfpsize;
      if(debug) fprintf(stderr, "(decompress) compressed = %ld, bufsize = %ld, compression = %4.1f, time = %f, speed = %fMB/s\n",
	      zfpsize, bufsize, ratio, time, (nx*ny*nz*4/1024.0f/1024.0f)/time);
    }
  } else {                                          // compress array into output buffer
    if(zfp_write_header(zfp, field, ZFP_HEADER_FULL) != 0) {   // store compression metadata in stream
      c = clock();
      zfpsize = zfp_compress(zfp, field);
      time = (double)(clock() - c) / CLOCKS_PER_SEC;
      if (!zfpsize) {
	if(debug) fprintf(stderr, "ERROR: compression failed\n");
	status = 1;
      }
      *Ssize = zfpsize ;
      ratio = nx*ny*nz*4 ; ratio /= zfpsize;
//       if(debug) 
        fprintf(stderr, "(compress  ) compressed = %ld, bufsize = %ld, compression = %4.1f, time = %f, speed = %fMB/s\n",
	      zfpsize, bufsize, ratio, time, (nx*ny*nz*4/1024.0f/1024.0f)/time);
    }else{
      if(debug) fprintf(stderr, "ERROR: write_header failed\n");
      status = 1;
    }
  }
  /* clean up */
  zfp_field_free(field);
  zfp_stream_close(zfp);
  stream_close(stream);
  if(debug > 1) {
    for(i = nx*ny/2 ; i < nx*ny*nz ; i += nx*ny*nz/10) {
      fprintf(stderr, "%9.4f ",array[i]); 
    }
  fprintf(stderr, "\n");
  }

  if(decompress) free(buffer);
  return status ? NULL : buffer;    // return buffer address if successful, NULL otherwise
}
#if defined(SELF_TEST)
int main(int argc, char* argv[])
{
  /* allocate array of floats */
  int nx = 1000;
  int ny = 1000;
  int nz = 1;
  float* array = malloc(nx * ny * nz * sizeof(float));
  float* brray = malloc(nx * ny * nz * sizeof(float));
  float maxval, minval;
  double err, errmax, errsum, errsuma;
  int i, j, k;
  void *Stream;
  int Ssize;
  float toler;

  /* initialize array to be compressed */
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
        float x = 2.0 * i / nx;
        float y = 2.0 * j / ny;
        float z = 2.0 * k / nz;
        array[i + nx * (j + ny * k)] = (x * x + y * y + z * z) + .002f;
//         array[i + nx * (j + ny * k)] = exp(-(x * x + y * y + z * z));
      }

  argc--;
  while(argc-- > 0) {
    sscanf(argv[1],"%f",&toler);                                  // get precision control
    Stream = ZfpCompress3D(array, nx, ny, nz, toler, NULL, &Ssize);    // compress
    Stream = ZfpCompress3D(brray, nx, ny, nz,   0.0, Stream, &Ssize);  // decompress
    AnalyzeCompressionErrors(array, brray, nx*ny*nz, 0.0f, "Ctest");           // evaluate errors
    argv++;
  }
}
#endif
