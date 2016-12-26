#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

// NOTE: Fortran order assumed for matrices
int TransposeBy4bytes(void *a, int la1, void *b, int lb1, int ni, int nj){
#if defined(NEVER_TRUE)
  basic block, 8x8 transpose of 32 bit elements within 256 bit AVX/AVX2 registers
  3 pass shuffle, first 2 operations are within 128 bit lane, last one is across 128 bit lanes

     low part   high part                   low part   high part                   low part   high part
   <---------> <--------->  32 bit op     <---------> <--------->  64 bit op     <---------> <--------->  128 bit op    <---------> <--------->
0  a0 a1 a2 a3 a4 a5 a6 a7  unpacklo(0,1) a0 b0 a1 b1 a4 b4 a5 b5  unpacklo(0,2) a0 b0 c0 d0 a1 b1 c1 d1  perm(0,4)LL   a0 b0 c0 d0 e0 f0 g0 h0
1  b0 b1 b2 b3 b4 b5 b6 b7  unpackhi(0,1) a2 b2 a3 b3 a6 b6 a7 b7  unpackhi(0,2) a4 b4 c4 d4 a5 b5 c5 d5  perm(0,4)HH   a1 b1 c1 d1 e1 f1 g1 h1
2  c0 c1 c2 c3 c4 c5 c6 c7  unpacklo(2,3) c0 d0 c1 d1 c4 d4 c5 d5  unpacklo(1,3) a2 b2 c2 d2 a3 b3 c3 d3  perm(2,6)LL   a2 b2 c2 d2 e2 f2 g2 h2
3  d0 d1 d2 d3 d4 d5 d6 d7  unpackhi(2,3) c2 d2 c3 d3 c6 d6 c7 d7  unpackhi(1,3) a6 b6 c6 d6 a7 b7 c7 d7  perm(2,6)HH   a3 b3 c3 d3 e3 f3 g3 h3
4  e0 e1 e2 e3 e4 e5 e6 e7  unpacklo(4,5) e0 f0 e1 f1 e4 f4 e5 f5  unpacklo(4,6) e0 f0 g0 h0 e1 f1 g1 h1  perm(1,5)LL   a4 b4 c4 d4 e4 f4 g4 h4
5  f0 f1 f2 f3 f4 f5 f6 f7  unpackhi(4,5) e2 f2 e3 f3 e6 f6 e7 f7  unpackhi(4,6) e4 f4 g4 h4 e5 f5 g5 h5  perm(1,5)HH   a5 b5 c5 d5 e5 f5 g5 h5
6  g0 g1 g2 g3 g4 g5 g6 g7  unpacklo(6,7) g0 h0 g1 h1 g4 h4 g5 h5  unpacklo(5,7) e2 f2 g2 h2 e3 f3 g3 h3  perm(3,7)LL   a6 b6 c6 d6 e6 f6 g6 h6
7  h0 h1 h2 h3 h4 h5 h6 h7  unpackhi(6,7) g2 h2 g3 h3 g6 h6 g7 h7  unpackhi(5,7) e6 f6 g6 h6 e7 f7 g7 h7  perm(3,7)HH   a7 b7 c7 d7 e7 f7 g7 h7

permLL uses _mm256_permute2f128_si256 (a,b, 0x20), bot128 = bot128 from a, top128 = bot128 from b
permHH uses _mm256_permute2f128_si256 (a,b, 0x31), bot128 = top128 from a, top128 = top128 from b
#endif
  int *a00, *a01, *b00, *b01;
  int *a000, *b000;
  int i, j, j0, ja0, ja00, jb0, jb00;
  int la2 = la1 + la1;
  int la3 = la2 + la1;
  int la4 = la3 + la1;
  int lb2 = lb1 + lb1;
  int lb3 = lb2 + lb1;
  int lb4 = lb3 + lb1;
  int ni8 ;
  int nj8 ;
  __m256i t0, t1, t2, t3, t4, t5, t6, t7;
  __m256i y0, y1, y2, y3, y4, y5, y6, y7;
  __m256i x0, x1, x2, x3, x4, x5, x6, x7;

  if (ni <= 0) ni = la1 ;
  if (nj <= 0) nj = lb1;
  ni8 = ni & 0x7FFFFFF8 ;   // lower multiple of 8
  nj8 = nj & 0x7FFFFFF8 ;   // lower multiple of 8
  if ( (la1 < ni) || (lb1 < nj) ) {
    fprintf(stderr,"ERROR: TransposeBy4bytes, la1=%d < ni=%d or lb1=%d < nj=%d\n",la1,ni,lb1,nj);
    return(1);
  }
// b[j,i] = a[i,j]
  a000 = (int *)a;
  b000 = (int *)b;
// basic block for transpose is 8x8
  for (j=0 ; j<nj8 ; j+= 8){
    a00 = a000;   // base address for row
    b00 = b000;   // base address for column
    for ( i=0 ; i<ni8 ; i+=8){     // transpose a[i:i+7,J:j+7] into b[j:j+7,i:i+7]
      a01 = a00 + la4;   // base of next block of 4x8
      b01 = b00 + lb4;   // base of next block of 4x8
      y0 = _mm256_loadu_si256 ((__m256i const *) &a00[  0]);
      y1 = _mm256_loadu_si256 ((__m256i const *) &a00[la1]);
      y2 = _mm256_loadu_si256 ((__m256i const *) &a00[la2]);
      y3 = _mm256_loadu_si256 ((__m256i const *) &a00[la3]);
      y4 = _mm256_loadu_si256 ((__m256i const *) &a01[  0]);
      y5 = _mm256_loadu_si256 ((__m256i const *) &a01[la1]);
      y6 = _mm256_loadu_si256 ((__m256i const *) &a01[la2]);
      y7 = _mm256_loadu_si256 ((__m256i const *) &a01[la3]);

      t0 = _mm256_unpacklo_epi32(y0,y1);
      t1 = _mm256_unpackhi_epi32(y0,y1);
      t2 = _mm256_unpacklo_epi32(y2,y3);
      t3 = _mm256_unpackhi_epi32(y2,y3);
      t4 = _mm256_unpacklo_epi32(y4,y5);
      t5 = _mm256_unpackhi_epi32(y4,y5);
      t6 = _mm256_unpacklo_epi32(y6,y7);
      t7 = _mm256_unpackhi_epi32(y6,y7);

      y0 = _mm256_unpacklo_epi64(t0,t2);
      y1 = _mm256_unpackhi_epi64(t0,t2);
      y2 = _mm256_unpacklo_epi64(t1,t3);
      y3 = _mm256_unpackhi_epi64(t1,t3);
      y4 = _mm256_unpacklo_epi64(t4,t6);
      y5 = _mm256_unpackhi_epi64(t4,t6);
      y6 = _mm256_unpacklo_epi64(t5,t7);
      y7 = _mm256_unpackhi_epi64(t5,t7);

      t0 = _mm256_permute2f128_si256 (y0,y4, 0x20);
      t1 = _mm256_permute2f128_si256 (y1,y5, 0x20);
      t2 = _mm256_permute2f128_si256 (y2,y6, 0x20);
      t3 = _mm256_permute2f128_si256 (y3,y7, 0x20);
      t4 = _mm256_permute2f128_si256 (y0,y4, 0x31);
      t6 = _mm256_permute2f128_si256 (y2,y6, 0x31);
      t5 = _mm256_permute2f128_si256 (y1,y5, 0x31);
      t7 = _mm256_permute2f128_si256 (y3,y7, 0x31);

      _mm256_storeu_si256((__m256i *) &b00[  0],t0);
      _mm256_storeu_si256((__m256i *) &b00[lb1],t1);
      _mm256_storeu_si256((__m256i *) &b00[lb2],t2);
      _mm256_storeu_si256((__m256i *) &b00[lb3],t3);
      _mm256_storeu_si256((__m256i *) &b01[  0],t4);
      _mm256_storeu_si256((__m256i *) &b01[lb1],t5);
      _mm256_storeu_si256((__m256i *) &b01[lb2],t6);
      _mm256_storeu_si256((__m256i *) &b01[lb3],t7);
      a00 = a00 + 8;           // bump by 8 cols
      b00 = b00 + lb4 + lb4;   // bump by 8 rows
    } // for i
    for ( ; i<ni ; i++){  // leftovers along i , j -> j+7
      b00[0] = a00[  0] ;  // b[j  ,i] = a[i,j  ]
      b00[1] = a00[la1] ;  // b[j+1,i] = a[i,j+1]
      b00[2] = a00[la2] ;  // b[j+2,i] = a[i,j+2]
      b00[3] = a00[la3] ;  // b[j+3,i] = a[i,j+3]
      b00[4] = a00[la4    ] ;  // b[j+4,i] = a[i,j+4]
      b00[5] = a00[la4+la1] ;  // b[j+5,i] = a[i,j+5]
      b00[6] = a00[la4+la2] ;  // b[j+6,i] = a[i,j+6]
      b00[7] = a00[la4+la3] ;  // b[j+7,i] = a[i,j+7]
      b00 += lb1;
      a00++;
    }
    a000 = a000 + la4 + la4;   // bump by 8 rows
    b000 = b000 + 8;           // bump by 8 cols
  }  // for j
  a000 = (int *)a;
  b000 = (int *)b;
  jb00 = j;
  ja00 = j*la1;
  for (i = 0 ; i<ni ; i++){    // leftovers along j, i = 0 -> nj-1
    ja0 = ja00;
    jb0 = jb00;
    for(j0=j ; j0<nj ; j0++){
      b000[jb0] = a000[ja0];   // b[j,i] = a[i,j]
      ja0 += la1;
      jb0++;
    }
    jb00 += lb1;
    ja00++;
  }
  return(0);
}

int TransposeBy8bytes(void *a, int la1, void *b, int lb1, int ni, int nj){
  // after transpose, b[j,i] = a[i,j]
  int la2 = la1 + la1;
  int la3 = la2 + la1;
  int la4 = la3 + la1;
  int lb2 = lb1 + lb1;
  int lb3 = lb2 + lb1;
  int lb4 = lb3 + lb1;
  int i, j, jj, ni8, nj8;
  long long *a000, *a00, *a01, *a10, *a11, *b000, *b00, *b01, *b10, *b11;
  __m256i y0, y1, y2, y3;
  __m256i t0, t1, t2, t3;
  __m256i x0, x1, x2, x3;

  if (ni <= 0) ni = la1 ;
  if (nj <= 0) nj = lb1;
  if ( (la1 < ni) || (lb1 < nj) ) {
    fprintf(stderr,"ERROR: TransposeBy8bytes, la1=%d < ni=%d or lb1=%d < nj=%d\n",la1,ni,lb1,nj);
    return(1);
  }
//   printf("la1,ni,lb1,nj= %d %d %d %d\n",la1,ni,lb1,nj);
  a000 = (long long *)a;
  b000 = (long long *)b;
  ni8 = ni & 0x7FFFFFF8 ;   // lower multiple of 8
  nj8 = nj & 0x7FFFFFF8 ;   // lower multiple of 8
  for (j=0 ; j<nj8 ; j+= 8){
    for ( i=0 ; i<ni8 ; i+=8){
      a00 = &a000[i + j*la1] ; a01 = a00 + 4 ; a10 = a00 + la4 ; a11 = a10 + 4;
      b00 = &b000[j + i*lb1] ; b01 = b00 + 4 ; b10 = b00 + lb4 ; b11 = b10 + 4;
      // basic 8 x 8 transpose (a00, a01, a10, a11 into b00, b10, b01, b11) ( 4 4x4 sub matrices)
      y0 = _mm256_loadu_si256((__m256i const *)&a00[  0]);   // fetch 4 x 4 from a00
      y1 = _mm256_loadu_si256((__m256i const *)&a00[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a00[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a00[la3]);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a00 -> b00
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      y0 = _mm256_loadu_si256((__m256i const *)&a01[  0]);   // prefetch 4 x 4 from a01
      y1 = _mm256_loadu_si256((__m256i const *)&a01[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a01[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a01[la3]);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a00 -> b00
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b00[  0],x0);   // store 4 x 4 into b00 (from a00)
      _mm256_storeu_si256((__m256i *)&b00[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b00[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b00[lb3],x3);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a01 -> b10
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      y0 = _mm256_loadu_si256((__m256i const *)&a10[  0]);   // prefetch 4 x 4 from a10
      y1 = _mm256_loadu_si256((__m256i const *)&a10[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a10[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a10[la3]);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a01 -> b10
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b10[  0],x0);   // store 4 x 4 into b10 (from a01)
      _mm256_storeu_si256((__m256i *)&b10[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b10[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b10[lb3],x3);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a10 -> b01
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      y0 = _mm256_loadu_si256((__m256i const *)&a11[  0]);   // prefetch 4 x 4 from a11
      y1 = _mm256_loadu_si256((__m256i const *)&a11[la1]);
      y2 = _mm256_loadu_si256((__m256i const *)&a11[la2]);
      y3 = _mm256_loadu_si256((__m256i const *)&a11[la3]);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a10 -> b01
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b01[  0],x0);   // store 4 x 4 into b01 (from a10)
      _mm256_storeu_si256((__m256i *)&b01[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b01[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b01[lb3],x3);

      t0 = _mm256_unpackhi_epi64(y0,y1);    // shuffle pass 1 on a11 -> b11
      t1 = _mm256_unpackhi_epi64(y2,y3);
      t2 = _mm256_unpacklo_epi64(y0,y1);
      t3 = _mm256_unpacklo_epi64(y2,y3);

      x0 = _mm256_permute2f128_si256(t2,t3,0x20);    // shuffle pass 2 on a11 -> b11
      x1 = _mm256_permute2f128_si256(t0,t1,0x20);
      x2 = _mm256_permute2f128_si256(t2,t3,0x31);
      x3 = _mm256_permute2f128_si256(t0,t1,0x31);

      _mm256_storeu_si256((__m256i *)&b11[  0],x0);   // store 4 x 4 into b11 (from a11)
      _mm256_storeu_si256((__m256i *)&b11[lb1],x1);
      _mm256_storeu_si256((__m256i *)&b11[lb2],x2);
      _mm256_storeu_si256((__m256i *)&b11[lb3],x3);
    }
    for ( ; i<ni ; i++){
      for ( jj=j ; jj<j+8 ; jj++){   // b[jj,i] = a[i,jj]
	b000[jj + lb1*i] = a000[i + jj*la1];
      }
    }
  }
  for ( ; j<nj ; j++){
    for (i=0 ; i<ni ; i++){   // b[j,i] = a[i,j]
      b000[j + lb1*i] = a000[i + j*la1];
    }
  }
  return(0);
}

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
main(int argc, char **argv){
  double *mtx1_4x4;
  double *mtx2_4x4;
  int *mat1_4x4;
  int *mat2_4x4;
  int i, j, k;
  double MPI_Wtime();
  double t0, t1;
  double bytes;
  int NI, NJ, ni, nj;
  int nni=atoi(argv[1]);
  int nnj=atoi(argv[2]);
  int nk=atoi(argv[5]);
  int errors;

  NI = atoi(argv[3]);
  NJ = atoi(argv[4]);
  ni = nni;
  nj = nnj;
  if(ni <= 0) ni = NI;
  if(nj <= 0) nj = NJ;
  mtx1_4x4 = (double *) malloc(sizeof(double)*NI*NJ);
  mtx2_4x4 = (double *) malloc(sizeof(double)*NI*NJ);
  mat1_4x4 = (int *) malloc(sizeof(int)*NI*NJ);
  mat2_4x4 = (int *) malloc(sizeof(int)*NI*NJ);
  printf("ni = %d, nj = %d, nk= %d\n",ni,nj,nk);
  for (i=0 ; i<NI*NJ ; i++){
    mtx1_4x4[i] = i;
    mtx2_4x4[i] = 999;
    mat1_4x4[i] = i;
    mat2_4x4[i] = 999;
  }
  printf("\nInput Matrix\n");
  if(NI*NJ < 1000) {
    for(j=NJ-1 ; j>= 0 ; j--){
      printf("row %2d ",j);
      for (i=0 ; i<NI ; i++){
	printf("%4.0f",mtx1_4x4[i+NI*j]);
      }
      printf("\n");
    }
    printf("Col  < ");
    for (i=0 ; i<NI ; i++) printf("%4d",i) ;
    printf(">\n");
  }

  if(TransposeBy8bytes(mtx1_4x4, NI, mtx2_4x4, NJ, nni, nnj)) exit(1);
  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk/10 ; k++) {
    for(j=0 ; j<nj ; j++){
      for(i=0 ; i<ni ; i++){
	mat2_4x4[j + i*NJ] = mat1_4x4[i + NI*j] ;
      }
    }
    bytes = bytes + ni * nj * 16;
  }
  t1 = MPI_Wtime() ;
  printf("\nT8ref = %f, %f GB/s, %f ns/pt, %f us/transpose\n",
	 (t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*10*1000*1000*1000,(t1-t0)/nk*10*1000*1000);
  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk ; k++) {
    if(TransposeBy8bytes(mtx1_4x4, NI, mtx2_4x4, NJ, nni, nnj)) exit(1);
    bytes = bytes + ni * nj * 16;
  }
  t1 = MPI_Wtime() ;
  printf("\nT8 = %f, %f GB/s, %f ns/pt, %f us/transpose\n",
	 (t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*1000*1000*1000,(t1-t0)/nk*1000*1000);
  errors = 0;
  for(j=0 ; j<nj ; j++){
    for(i=0 ; i<ni ; i++){
      if( mtx1_4x4[i + NI*j] != mtx2_4x4[j + i*NJ] ) errors++;
    }
  }
  printf("\nOutput Matrix (double), errors = %d\n\n",errors);
  if(NI*NJ < 1000) {
    for(j=NI-1 ; j>= 0 ; j--){
      printf("row %2d ",j);
      for (i=0 ; i<NJ ; i++){
	printf("%4.0f",mtx2_4x4[i+NJ*j]);
      }
      printf("\n");
    }
    printf("Col  < ");
    for (i=0 ; i<NJ ; i++) printf("%4d",i) ;
    printf(">\n");
  }

  if(TransposeBy4bytes(mat1_4x4, NI, mat2_4x4, NJ, nni, nnj)) exit(1);
  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk/10 ; k++) {
    for(j=0 ; j<nj ; j++){
      for(i=0 ; i<ni ; i++){
	mat2_4x4[j + i*NJ] = mat1_4x4[i + NI*j] ;
      }
    }
    bytes = bytes + ni * nj * 8;
  }
  t1 = MPI_Wtime() ;
  printf("\nT4ref = %f, %f GB/s, %f ns/pt, %f us/transpose\n",
	 (t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*10*1000*1000*1000,(t1-t0)/nk*10*1000*1000);
  t0 = MPI_Wtime();
  bytes = 0;
  for (k=0 ; k<nk ; k++) {
    if(TransposeBy4bytes(mat1_4x4, NI, mat2_4x4, NJ, nni, nnj)) exit(1);
    bytes = bytes + ni * nj * 8;
  }
  t1 = MPI_Wtime() ;
  errors = 0;
  for(j=0 ; j<nj ; j++){
    for(i=0 ; i<ni ; i++){
      if( mat1_4x4[i + NI*j] != mat2_4x4[j + i*NJ] ) errors++;
    }
  }
  printf("\nOutput Matrix (int), errors = %d\n",errors);
  printf("T4 = %f, %f GB/s, %f ns/pt, %f us/transpose\n",(t1-t0),bytes/1000./1000./1000./(t1-t0),(t1-t0)/ni/nj/nk*1000*1000*1000,(t1-t0)/nk*1000*1000);
  if(NI*NJ < 1000) {
    for(j=NI-1 ; j>= 0 ; j--){
      printf("row %2d ",j);
      for (i=0 ; i<NJ ; i++){
	printf("%4d",mat2_4x4[i+NJ*j]);
      }
      printf("\n");
    }
    printf("Col  < ");
    for (i=0 ; i<NJ ; i++) printf("%4d",i) ;
    printf(">\n");
  }
}
#endif