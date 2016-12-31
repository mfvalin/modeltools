#include <immintrin.h>

// enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
// enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
// void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
//                  const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
//                  const int K, const double alpha, const double *A,
//                  const int lda, const double *B, const int ldb,
//                  const double beta, double *C, const int ldc);

void dgemm_shortnk(int ni, int nj, int nk, double *ma, int lda1, double *mb, int ldb1, double *mc, int ldc1) {
//   double alpha = 1.0;
//   double beta = 0.0;
//   if(nk > 100) {   // use regular dgemm if nk > 100
//     cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ni,nj,nk,alpha,ma,lda,mb,ldb,beta,mc,ldc)
//   }
#if defined(UNROLL_8x4)
  __m256d c00, c10, c01, c11, c02, c12, c03, c13 ; // unroll 8 x along i, 4x along j
  __m256d b0, b1, b2, b3 ;
  __m256d a0, a1, t0, t1 ;
  double *ma0, *ma00, *mb0, *mb00, *mc0, *mc00 ;
  int i, i0, j, j0, k, k0 ;
  int ldb2, ldb3, ldb4;

  ldb2 = ldb1 + ldb1 ;
  ldb3 = ldb2 + ldb1 ;

  // basic block, c[i,j] = c[i,j] + a[i,k]*b[k,j]
  for (j=0 ; j<nj-4 ; j+=4) {
    for (i=0 ; i<ni-8 ; i+=8) {
      ma0 = &ma[i];                            // &a[i,0]
      mb0 = &mb[j*ldb1];                       // &b[0,j]
      a0 = _mm256_loadu_pd (&ma0[0]) ;         // a[i:i+3,0]
      b0 = _mm256_broadcast_sd (&mb0[   0]) ;  // b[0,  j]
      a1 = _mm256_loadu_pd (&ma0[4]) ;         // a[i+4:i+7,0]
      b1 = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[0,j+1]
// k = 0, fma is replaced with multiply
      ma0 += lda1 ;                            // next value of k for matrix a
      t0 = _mm256_loadu_pd (&ma0[0]) ;         // a[i  :i+3,k+1]
      t1 = _mm256_loadu_pd (&ma0[4]) ;         // a[i+4:i+7,k+1]
      b2 = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k,j+2]
      c00 =  _mm256_mul_pd (a0,b0) ;
      c10 =  _mm256_mul_pd (a1,b0) ;
      b3 = _mm256_broadcast_sd (&mb0[ldb3]) ;  // b[k,j+3]
      mb0 ++;                                  // next value of k for matrix b
      c01 =  _mm256_mul_pd (a0,b1) ;
      c11 =  _mm256_mul_pd (a1,b1) ;
      b0 = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k+1,  j]
      c02 =  _mm256_mul_pd (a0,b2) ;
      c12 =  _mm256_mul_pd (a1,b2) ;
      b1 = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k+1,j+1]
      c03 =  _mm256_mul_pd (a0,b3) ;
      c13 =  _mm256_mul_pd (a1,b3) ;

      for (k=1 ; k<nk-1 ; k+=2){    // do iterations k and k+1
// k
	ma0 += lda1 ;                            // next value of k for matrix a
	b2 = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k,j+2]
	a0 = _mm256_loadu_pd (&ma0[0]) ;         // a[i  :i+3,k+1]
	a1 = _mm256_loadu_pd (&ma0[4]) ;         // a[i+4:i+7,k+1]
	c00 =  _mm256_fmadd_pd (t0,b0,c00) ;
	c10 =  _mm256_fmadd_pd (t1,b0,c10) ;
	b3 = _mm256_broadcast_sd (&mb0[ldb3]) ;  // b[k,j+3]
	mb0 ++;                                  // next value of k for matrix b
	c01 =  _mm256_fmadd_pd (t0,b1,c01) ;
	c11 =  _mm256_fmadd_pd (t1,b1,c11) ;
	b0 = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k+1,  j]
	c02 =  _mm256_fmadd_pd (t0,b2,c02) ;
	c12 =  _mm256_fmadd_pd (t1,b2,c12) ;
	b1 = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k+1,j+1]
	c03 =  _mm256_fmadd_pd (t0,b3,c03) ;
	c13 =  _mm256_fmadd_pd (t1,b3,c13) ;
// k + 1
	ma0 += lda1 ;                            // next value of k for matrix a
	b2 = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k+1,j+2]
	t0 = _mm256_loadu_pd (&ma0[0]) ;         // a[i  :i+3,k+2]
	t1 = _mm256_loadu_pd (&ma0[4]) ;         // a[i+4:i+7,k+2]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	c10 =  _mm256_fmadd_pd (a1,b0,c10) ;
	b3 = _mm256_broadcast_sd (&mb0[ldb3]) ;  // b[k+1,j+3]
	mb0 ++;                                  // next value of k for matrix b
	c01 =  _mm256_fmadd_pd (a0,b1,c01) ;
	c11 =  _mm256_fmadd_pd (a1,b1,c11) ;
	b0 = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k+2,  j]
	c02 =  _mm256_fmadd_pd (a0,b2,c02) ;
	c12 =  _mm256_fmadd_pd (a1,b2,c12) ;
	b1 = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k+2,j+1]
	c03 =  _mm256_fmadd_pd (a0,b3,c03) ;
	c13 =  _mm256_fmadd_pd (a1,b3,c13) ;
      }  // for k
      if(k == nk-1){   // nk was even, one item left, compute, no prefetch
	b2 = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k,j+2]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	c10 =  _mm256_fmadd_pd (a1,b0,c10) ;
	b3 = _mm256_broadcast_sd (&mb0[ldb3]) ;  // b[k,j+3]
	c01 =  _mm256_fmadd_pd (a0,b1,c01) ;
	c11 =  _mm256_fmadd_pd (a1,b1,c11) ;
	c02 =  _mm256_fmadd_pd (a0,b2,c02) ;
	c12 =  _mm256_fmadd_pd (a1,b2,c12) ;
	c03 =  _mm256_fmadd_pd (a0,b3,c03) ;
	c13 =  _mm256_fmadd_pd (a1,b3,c13) ;
// 	printf(".");
      }
// store c[i:i+7,j:j+4]
      mc0 = &mc[i + j*ldc1] ;
      _mm256_storeu_pd (&mc0[0], c00) ;
      _mm256_storeu_pd (&mc0[4], c10) ;
      mc0 += ldc1 ;
      _mm256_storeu_pd (&mc0[0], c01) ;
      _mm256_storeu_pd (&mc0[4], c11) ;
      mc0 += ldc1 ;
      _mm256_storeu_pd (&mc0[0], c02) ;
      _mm256_storeu_pd (&mc0[4], c12) ;
      mc0 += ldc1 ;
      _mm256_storeu_pd (&mc0[0], c03) ;
      _mm256_storeu_pd (&mc0[4], c13) ;
    }
  }
#else
  __m256d c00, c10, c20, c30, c01, c11, c21, c31, c02, c12, c22, c32 ; // unroll 16 x along i, 3x along j (12 registers)
  __m256d b0, b1, b2 ; // j unroll for B (3 registers)
  __m256d a0 ; // used to fetch consecutive terms from matrix a
  double *ma0, *ma00, *mb0, *mc0, *mc00 ;
  double sc0, sc1, sc2;
  int i, j, k ;
  int ldb2, ldb3, ldc2;

  ldb2 = ldb1 + ldb1 ;
  ldb3 = ldb2 + ldb1 ;
  ldc2 = ldc1 + ldc1 ;

  for (j=0 ; j<nj-2 ; j+=3) {                    // blocks of 3 values of j
    for (i=0 ; i<ni-15 ; i+=16) {                // blocks of 16 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;   // b[0   ,j+1]
      c01 = _mm256_mul_pd (a0,b1) ;
      b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;   // b[0   ,j+2]
      c02 = _mm256_mul_pd (a0,b2) ;

      a0  = _mm256_loadu_pd (&ma0[ 4]) ;         // a[i+ 4,0  ]
      c10 = _mm256_mul_pd (a0,b0) ;
      c11 = _mm256_mul_pd (a0,b1) ;
      c12 = _mm256_mul_pd (a0,b2) ;

      a0  = _mm256_loadu_pd (&ma0[ 8]) ;         // a[i+ 8,0  ]
      c20 = _mm256_mul_pd (a0,b0) ;
      c21 = _mm256_mul_pd (a0,b1) ;
      c22 = _mm256_mul_pd (a0,b2) ;

      a0  = _mm256_loadu_pd (&ma0[12]) ;         // a[i+12,0  ]
      c30 = _mm256_mul_pd (a0,b0) ;
      c31 = _mm256_mul_pd (a0,b1) ;
      c32 = _mm256_mul_pd (a0,b2) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k   ,j+1]
	c01 =  _mm256_fmadd_pd (a0,b1,c01) ;
	b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k   ,j+2]
	c02 =  _mm256_fmadd_pd (a0,b2,c02) ;
	a0  = _mm256_loadu_pd (&ma0[4]) ;         // a[i+ 4,k  ]
	c10 = _mm256_fmadd_pd (a0,b0,c10) ;
	c11 = _mm256_fmadd_pd (a0,b1,c11) ;
	c12 = _mm256_fmadd_pd (a0,b2,c12) ;
	a0  = _mm256_loadu_pd (&ma0[8]) ;         // a[i+ 8,k  ]
	c20 = _mm256_fmadd_pd (a0,b0,c20) ;
	c21 = _mm256_fmadd_pd (a0,b1,c21) ;
	c22 = _mm256_fmadd_pd (a0,b2,c22) ;
	a0  = _mm256_loadu_pd (&ma0[12]) ;        // a[i+12,k  ]
	c30 = _mm256_fmadd_pd (a0,b0,c30) ;
	c31 = _mm256_fmadd_pd (a0,b1,c31) ;
	c32 = _mm256_fmadd_pd (a0,b2,c32) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+15,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      _mm256_storeu_pd (&mc0[ 4], c10) ;
      _mm256_storeu_pd (&mc0[ 8], c20) ;
      _mm256_storeu_pd (&mc0[12], c30) ;
      mc0 += ldc1 ;                               // c[i:i+15,j+1]
      _mm256_storeu_pd (&mc0[ 0], c01) ;
      _mm256_storeu_pd (&mc0[ 4], c11) ;
      _mm256_storeu_pd (&mc0[ 8], c21) ;
      _mm256_storeu_pd (&mc0[12], c31) ;
      mc0 += ldc1 ;                               // c[i:i+15,j+2]
      _mm256_storeu_pd (&mc0[ 0], c02) ;
      _mm256_storeu_pd (&mc0[ 4], c12) ;
      _mm256_storeu_pd (&mc0[ 8], c22) ;
      _mm256_storeu_pd (&mc0[12], c32) ;
    }
    for (i=i ; i<ni-11 ; i+=12) {                   // blocks of 12 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;   // b[0   ,j+1]
      c01 = _mm256_mul_pd (a0,b1) ;
      b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;   // b[0   ,j+2]
      c02 = _mm256_mul_pd (a0,b2) ;

      a0  = _mm256_loadu_pd (&ma0[ 4]) ;         // a[i+ 4,0  ]
      c10 = _mm256_mul_pd (a0,b0) ;
      c11 = _mm256_mul_pd (a0,b1) ;
      c12 = _mm256_mul_pd (a0,b2) ;

      a0  = _mm256_loadu_pd (&ma0[ 8]) ;         // a[i+ 8,0  ]
      c20 = _mm256_mul_pd (a0,b0) ;
      c21 = _mm256_mul_pd (a0,b1) ;
      c22 = _mm256_mul_pd (a0,b2) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k   ,j+1]
	c01 =  _mm256_fmadd_pd (a0,b1,c01) ;
	b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k   ,j+2]
	c02 =  _mm256_fmadd_pd (a0,b2,c02) ;
	a0  = _mm256_loadu_pd (&ma0[4]) ;         // a[i+ 4,k  ]
	c10 = _mm256_fmadd_pd (a0,b0,c10) ;
	c11 = _mm256_fmadd_pd (a0,b1,c11) ;
	c12 = _mm256_fmadd_pd (a0,b2,c12) ;
	a0  = _mm256_loadu_pd (&ma0[8]) ;         // a[i+ 8,k  ]
	c20 = _mm256_fmadd_pd (a0,b0,c20) ;
	c21 = _mm256_fmadd_pd (a0,b1,c21) ;
	c22 = _mm256_fmadd_pd (a0,b2,c22) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+15,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      _mm256_storeu_pd (&mc0[ 4], c10) ;
      _mm256_storeu_pd (&mc0[ 8], c20) ;
      mc0 += ldc1 ;                               // c[i:i+15,j+1]
      _mm256_storeu_pd (&mc0[ 0], c01) ;
      _mm256_storeu_pd (&mc0[ 4], c11) ;
      _mm256_storeu_pd (&mc0[ 8], c21) ;
      mc0 += ldc1 ;                               // c[i:i+15,j+2]
      _mm256_storeu_pd (&mc0[ 0], c02) ;
      _mm256_storeu_pd (&mc0[ 4], c12) ;
      _mm256_storeu_pd (&mc0[ 8], c22) ;
    }
    for (i=i ; i<ni-7 ; i+=8) {                   // blocks of 8 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;   // b[0   ,j+1]
      c01 = _mm256_mul_pd (a0,b1) ;
      b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;   // b[0   ,j+2]
      c02 = _mm256_mul_pd (a0,b2) ;

      a0  = _mm256_loadu_pd (&ma0[ 4]) ;         // a[i+ 4,0  ]
      c10 = _mm256_mul_pd (a0,b0) ;
      c11 = _mm256_mul_pd (a0,b1) ;
      c12 = _mm256_mul_pd (a0,b2) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k   ,j+1]
	c01 =  _mm256_fmadd_pd (a0,b1,c01) ;
	b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k   ,j+2]
	c02 =  _mm256_fmadd_pd (a0,b2,c02) ;
	a0  = _mm256_loadu_pd (&ma0[4]) ;         // a[i+ 4,k  ]
	c10 = _mm256_fmadd_pd (a0,b0,c10) ;
	c11 = _mm256_fmadd_pd (a0,b1,c11) ;
	c12 = _mm256_fmadd_pd (a0,b2,c12) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+11,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      _mm256_storeu_pd (&mc0[ 4], c10) ;
      mc0 += ldc1 ;                               // c[i:i+11,j+1]
      _mm256_storeu_pd (&mc0[ 0], c01) ;
      _mm256_storeu_pd (&mc0[ 4], c11) ;
      mc0 += ldc1 ;                               // c[i:i+11,j+2]
      _mm256_storeu_pd (&mc0[ 0], c02) ;
      _mm256_storeu_pd (&mc0[ 4], c12) ;
    }
    for (i=i ; i<ni-3 ; i+=4) {                   // blocks of 4 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;   // b[0   ,j+1]
      c01 = _mm256_mul_pd (a0,b1) ;
      b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;   // b[0   ,j+2]
      c02 = _mm256_mul_pd (a0,b2) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	b1  = _mm256_broadcast_sd (&mb0[ldb1]) ;  // b[k   ,j+1]
	c01 =  _mm256_fmadd_pd (a0,b1,c01) ;
	b2  = _mm256_broadcast_sd (&mb0[ldb2]) ;  // b[k   ,j+2]
	c02 =  _mm256_fmadd_pd (a0,b2,c02) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+3 ,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      mc0 += ldc1 ;                               // c[i:i+3 ,j+1]
      _mm256_storeu_pd (&mc0[ 0], c01) ;
      mc0 += ldc1 ;                               // c[i:i+3 ,j+2]
      _mm256_storeu_pd (&mc0[ 0], c02) ;
    }
    for (i=i ; i<ni ; i++) {                     // leftovers
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]
      mc0 = &mc[i + j*ldc1] ;                     // c[i ,j  ]
      sc0 = ma0[0]*mb0[0   ] ;
      sc1 = ma0[0]*mb0[ldb1] ;
      sc2 = ma0[0]*mb0[ldb2] ;
      for (k=1 ; k<nk ; k++){
	mb0++ ;
	ma0 += lda1 ;
	sc0 = sc0 + ma0[0]*mb0[0   ] ;
	sc1 = sc1 + ma0[0]*mb0[ldb1] ;
	sc2 = sc2 + ma0[0]*mb0[ldb2] ;
      }
      mc0[0   ] = sc0 ;
      mc0[ldc1] = sc1 ;
      mc0[ldc2] = sc2 ;
    }
//     mb00 += ldb3;
  }
  for (j=j ; j<nj ; j++){   // if 1 or 2 j values left
    for (i=0 ; i<ni-15 ; i+=16) {
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      a0  = _mm256_loadu_pd (&ma0[ 4]) ;         // a[i+ 4,0  ]
      c10 = _mm256_mul_pd (a0,b0) ;
      a0  = _mm256_loadu_pd (&ma0[ 8]) ;         // a[i+ 8,0  ]
      c20 = _mm256_mul_pd (a0,b0) ;
      a0  = _mm256_loadu_pd (&ma0[12]) ;         // a[i+12,0  ]
      c30 = _mm256_mul_pd (a0,b0) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	a0  = _mm256_loadu_pd (&ma0[4]) ;         // a[i+ 4,k  ]
	c10 = _mm256_fmadd_pd (a0,b0,c10) ;
	a0  = _mm256_loadu_pd (&ma0[8]) ;         // a[i+ 8,k  ]
	c20 = _mm256_fmadd_pd (a0,b0,c20) ;
	a0  = _mm256_loadu_pd (&ma0[12]) ;        // a[i+12,k  ]
	c30 = _mm256_fmadd_pd (a0,b0,c30) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+15,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      _mm256_storeu_pd (&mc0[ 4], c10) ;
      _mm256_storeu_pd (&mc0[ 8], c20) ;
      _mm256_storeu_pd (&mc0[12], c30) ;
    }
    for (i=i ; i<ni-11 ; i+=12) {                   // blocks of 12 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      a0  = _mm256_loadu_pd (&ma0[ 4]) ;         // a[i+ 4,0  ]
      c10 = _mm256_mul_pd (a0,b0) ;
      a0  = _mm256_loadu_pd (&ma0[ 8]) ;         // a[i+ 8,0  ]
      c20 = _mm256_mul_pd (a0,b0) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	a0  = _mm256_loadu_pd (&ma0[4]) ;         // a[i+ 4,k  ]
	c10 = _mm256_fmadd_pd (a0,b0,c10) ;
	a0  = _mm256_loadu_pd (&ma0[8]) ;         // a[i+ 8,k  ]
	c20 = _mm256_fmadd_pd (a0,b0,c20) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+11,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      _mm256_storeu_pd (&mc0[ 4], c10) ;
      _mm256_storeu_pd (&mc0[ 8], c20) ;
    }
    for (i=i ; i<ni-7 ; i+=8) {                   // blocks of 8 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      a0  = _mm256_loadu_pd (&ma0[ 4]) ;         // a[i+ 4,0  ]
      c10 = _mm256_mul_pd (a0,b0) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
	a0  = _mm256_loadu_pd (&ma0[4]) ;         // a[i+ 4,k  ]
	c10 = _mm256_fmadd_pd (a0,b0,c10) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+11,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
      _mm256_storeu_pd (&mc0[ 4], c10) ;
    }
    for (i=i ; i<ni-3 ; i+=4) {                   // blocks of 4 values of i
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]

      b0  = _mm256_broadcast_sd (&mb0[   0]) ;   // b[0   ,j  ]
      a0  = _mm256_loadu_pd (&ma0[ 0]) ;         // a[i   ,0  ]
      c00 = _mm256_mul_pd (a0,b0) ;
      for (k=1 ; k<nk ; k++){    // do iterations 1 up to nk-1
	mb0++ ;                                  // &b[k  ,j  ]
	ma0 += lda1 ;                            // &a[i  ,k  ]
	b0  = _mm256_broadcast_sd (&mb0[   0]) ;  // b[k   ,j  ]
	a0  = _mm256_loadu_pd (&ma0[0]) ;         // a[i   ,k  ]
	c00 =  _mm256_fmadd_pd (a0,b0,c00) ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i:i+11,j  ]
      _mm256_storeu_pd (&mc0[ 0], c00) ;
    }
    for (i=i ; i<ni ; i++) {                     // leftovers
      ma0 = &ma[i];                              // &a[i  ,0  ]
      mb0 = &mb[j*ldb1];                         // &b[0  ,j  ]
      sc0 = ma0[0]*mb0[0   ] ;
      for (k=1 ; k<nk ; k++){
	mb0++ ;
	ma0 += lda1 ;
	sc0 = sc0 + ma0[0]*mb0[0   ] ;
      }
      mc0 = &mc[i + j*ldc1] ;                     // c[i ,j  ]
      mc0[0   ] = sc0 ;
    }
  }

#endif
}
#if defined SELF_TEST
#include <stdio.h>
#include <stdlib.h>

#define NI 260
#define NK 32
#define NJ 201
int main(int argc, char **argv){
  double c[NI*NJ];
  double r[NI*NJ];   // reference result
  double a[NI*NK];
  double b[NJ*NK];
  int i,j,k,errors;
  double t0, t1, flops;
  int nrep = atoi(argv[1]);
  double MPI_Wtime();

  for (i=0 ; i<NI*NK ; i++) a[i] =  (1.0*i)/(NI*NK) + 1.0;
  for (i=0 ; i<NK*NJ ; i++) b[i] =  (1.0*i)/(NJ*NK) + 2.0;
  for (i=0 ; i<NJ*NI ; i++) c[i] = 0;
  for (i=0 ; i<NJ*NI ; i++) r[i] = 0;

  for (j=0 ; j<NJ ; j++){
    for (i=0 ; i<NI ; i++){
      for(k=0 ; k<NK ; k++){
	r[i + j*NI] = r[i + j*NI] + a[i + k*NI] * b[k + j*NK];
      }
    }
  }

  t0 = MPI_Wtime();
  for (k=0 ; k<nrep ; k++)
    dgemm_shortnk(NI, NJ, NK, a, NI, b, NK, c, NI);
  t1 = MPI_Wtime();
  errors = 0;
  flops = NI*NJ*NK*2 ;
  flops = flops * nrep ;
  flops = flops / (t1-t0) * 1.0E-9;
  fprintf(stdout,"speed = %6.3f GFlops\n",flops);

  for (i=0 ; i<NJ*NI ; i++) {
//     if(c[i] != NK) errors++;
    if(c[i] != r[i]) errors++;
  }
  fprintf(stdout,"errors = %d\n",errors);
  return(0);

  for(j=NJ-1 ; j>=0 ; j--){
    printf("%3d ",j);
    for (i=0 ; i<NI ; i++){
      printf("%3.0f",c[i+NI*j]);
    }
    printf("\n");
  }
  return(0);
}
#endif
