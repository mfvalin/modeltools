#include <immintrin.h>

void transpose_double_by_8(double *a, int la1, double *b, int lb1, int ni, int nj){
  // after transpose, b[j,i] = a[i,j]
  int la2 = la1 + la1;
  int la3 = la2 + la1;
  int la4 = la3 + la1;
  int lb2 = lb1 + lb1;
  int lb3 = lb2 + lb1;
  int lb4 = lb3 + lb1;
  int i, j;
  double *a00, *a01, *a10, *a11, *b00, *b01, *b10, *b11;
  __m256d y0, y1, y2, y3;
  __m256d t0, t1, t2, t3;
  __m256d x0, x1, x2, x3;

  a00 = a ; a01 = a + 4 ; a10 = a00 + la4 ; a11 = a10 + 4;
  b00 = b ; b01 = b + 4 ; b10 = b00 + lb4 ; b11 = b10 + 4;
  for (j=0 ; j<nj ; j+= 8){
    for ( i=0 ; i<ni ; i+=8){
      // basic 8 x 8 transpose (a00, a01, a10, a11 into b00, b10, b01, b11) ( 4 4x4 sub matrices)
      y0 = _mm256_loadu_pd(&a00[  0]);   // fetch 4 x 4 from a00
      y1 = _mm256_loadu_pd(&a00[la1]);
      y2 = _mm256_loadu_pd(&a00[la2]);
      y3 = _mm256_loadu_pd(&a00[la3]);

      t0 = _mm256_unpackhi_pd(y0,y1);    // shuffle pass 1 on a00 -> b00
      t1 = _mm256_unpackhi_pd(y2,y3);
      t2 = _mm256_unpacklo_pd(y0,y1);
      t3 = _mm256_unpacklo_pd(y2,y3);

      y0 = _mm256_loadu_pd(&a01[  0]);   // prefetch 4 x 4 from a01
      y1 = _mm256_loadu_pd(&a01[la1]);
      y2 = _mm256_loadu_pd(&a01[la2]);
      y3 = _mm256_loadu_pd(&a01[la3]);

      x0 = _mm256_permute2f128_pd(t2,t3,0x20);    // shuffle pass 2 on a00 -> b00
      x1 = _mm256_permute2f128_pd(t0,t1,0x20);
      x2 = _mm256_permute2f128_pd(t2,t3,0x31);
      x3 = _mm256_permute2f128_pd(t0,t1,0x31);

      _mm256_storeu_pd(&b00[  0],x0);   // store 4 x 4 into b00 (from a00)
      _mm256_storeu_pd(&b00[lb1],x1);
      _mm256_storeu_pd(&b00[lb2],x2);
      _mm256_storeu_pd(&b00[lb3],x3);

      t0 = _mm256_unpackhi_pd(y0,y1);    // shuffle pass 1 on a01 -> b10
      t1 = _mm256_unpackhi_pd(y2,y3);
      t2 = _mm256_unpacklo_pd(y0,y1);
      t3 = _mm256_unpacklo_pd(y2,y3);

      y0 = _mm256_loadu_pd(&a10[  0]);   // prefetch 4 x 4 from a10
      y1 = _mm256_loadu_pd(&a10[la1]);
      y2 = _mm256_loadu_pd(&a10[la2]);
      y3 = _mm256_loadu_pd(&a10[la3]);

      x0 = _mm256_permute2f128_pd(t2,t3,0x20);    // shuffle pass 2 on a01 -> b10
      x1 = _mm256_permute2f128_pd(t0,t1,0x20);
      x2 = _mm256_permute2f128_pd(t2,t3,0x31);
      x3 = _mm256_permute2f128_pd(t0,t1,0x31);

      _mm256_storeu_pd(&b10[  0],x0);   // store 4 x 4 into b10 (from a01)
      _mm256_storeu_pd(&b10[lb1],x1);
      _mm256_storeu_pd(&b10[lb2],x2);
      _mm256_storeu_pd(&b10[lb3],x3);

      t0 = _mm256_unpackhi_pd(y0,y1);    // shuffle pass 1 on a10 -> b01
      t1 = _mm256_unpackhi_pd(y2,y3);
      t2 = _mm256_unpacklo_pd(y0,y1);
      t3 = _mm256_unpacklo_pd(y2,y3);

      y0 = _mm256_loadu_pd(&a11[  0]);   // prefetch 4 x 4 from a11
      y1 = _mm256_loadu_pd(&a11[la1]);
      y2 = _mm256_loadu_pd(&a11[la2]);
      y3 = _mm256_loadu_pd(&a11[la3]);

      x0 = _mm256_permute2f128_pd(t2,t3,0x20);    // shuffle pass 2 on a10 -> b01
      x1 = _mm256_permute2f128_pd(t0,t1,0x20);
      x2 = _mm256_permute2f128_pd(t2,t3,0x31);
      x3 = _mm256_permute2f128_pd(t0,t1,0x31);

      _mm256_storeu_pd(&b01[  0],x0);   // store 4 x 4 into b01 (from a10)
      _mm256_storeu_pd(&b01[lb1],x1);
      _mm256_storeu_pd(&b01[lb2],x2);
      _mm256_storeu_pd(&b01[lb3],x3);

      t0 = _mm256_unpackhi_pd(y0,y1);    // shuffle pass 1 on a11 -> b11
      t1 = _mm256_unpackhi_pd(y2,y3);
      t2 = _mm256_unpacklo_pd(y0,y1);
      t3 = _mm256_unpacklo_pd(y2,y3);

      x0 = _mm256_permute2f128_pd(t2,t3,0x20);    // shuffle pass 2 on a11 -> b11
      x1 = _mm256_permute2f128_pd(t0,t1,0x20);
      x2 = _mm256_permute2f128_pd(t2,t3,0x31);
      x3 = _mm256_permute2f128_pd(t0,t1,0x31);

      _mm256_storeu_pd(&b11[  0],x0);   // store 4 x 4 into b11 (from a11)
      _mm256_storeu_pd(&b11[lb1],x1);
      _mm256_storeu_pd(&b11[lb2],x2);
      _mm256_storeu_pd(&b11[lb3],x3);
    }
  }
}

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#define NI 17
#define NJ 15
main(){
  double mtx1_4x4[NI*NJ];
  double mtx2_4x4[NI*NJ];
  int i, j;
  for (i=0 ; i<NI*NJ ; i++){
    mtx1_4x4[i] = i;
    mtx2_4x4[i] = 999;
  }
  printf("\nInput Matrix\n");
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
  transpose_double_by_8(mtx1_4x4, NI, mtx2_4x4, NJ, 8, 8);  // 4 x 4 transpose
  printf("\nOutput Matrix\n");
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
#endif