#include <immintrin.h>
#include <stdint.h>

static int I0123[] = {0,1,2,3};

void VGather4x4x4(int *src, int *dst, int ni, int ninj){
  int *base1, *base2;
#if defined(__AVX2__) && defined(__x86_64__)
  __m128i v128a, v128b;
  __m256i vindx;
  __m256i v1, v2, v3, v4, v5, v6, v7, v8;
#endif

#if defined(__AVX2__) && defined(__x86_64__)
  v128a = _mm_loadu_si128((__m128i const*) I0123);
  v128b = _mm_set1_epi32(ni);
  v128b = _mm_add_epi32(v128a,v128b);
  vindx = _mm256_inserti128_si256(vindx,v128a,0);
  vindx = _mm256_inserti128_si256(vindx,v128b,1);

  base1 = src;
  base2 = src + ni + ni;
  v1 = _mm256_i32gather_epi32((int const*) base1, vindx, 4);
  v2 = _mm256_i32gather_epi32((int const*) base2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst +  0), v1);
  _mm256_storeu_si256((__m256i *) (dst +  8), v2);

  base1 += ninj;
  base2 += ninj;
  v3 = _mm256_i32gather_epi32((int const*) base1, vindx, 4);
  v4 = _mm256_i32gather_epi32((int const*) base2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst + 16), v3);
  _mm256_storeu_si256((__m256i *) (dst + 24), v4);

  base1 += ninj;
  base2 += ninj;
  v5 = _mm256_i32gather_epi32((int const*) base1, vindx, 4);
  v6 = _mm256_i32gather_epi32((int const*) base2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst + 32), v5);
  _mm256_storeu_si256((__m256i *) (dst + 40), v6);

  base1 += ninj;
  base2 += ninj;
  v7 = _mm256_i32gather_epi32((int const*) base1, vindx, 4);
  v8 = _mm256_i32gather_epi32((int const*) base2, vindx, 4);
  _mm256_storeu_si256((__m256i *) (dst + 48), v7);
  _mm256_storeu_si256((__m256i *) (dst + 56), v8);

#else
  dst[0] = src[0];
#endif
}

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#define DATA_SIZE 9000
main(){
  int my_data[DATA_SIZE];
  int i, j;
  int ni, ninj, limit;
  int collect[64];

  for (i=0 ; i<DATA_SIZE ; i++) my_data[i] = i;
  ni = 100 ;
  ninj = 1000;
  limit = DATA_SIZE-4*ninj+1;
  printf("limit=%d\n",limit);
  for (j=0 ; j<100000 ; j++) {
    for (i=0 ; i<limit ; i++) VGather4x4x4(my_data+i, collect, ni, ninj);
  }
//   for (i=0 ; i<64 ; i++) collect[i] = -1;
//   VGather4x4x4(my_data, collect, 100, 10000);
//   for (i=0 ; i<64 ; i++) printf("I = %d, G = %d \n",i,collect[i]);
}
#endif
