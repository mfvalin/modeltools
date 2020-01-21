#include <stdint.h>
#include <immintrin.h>

static int64_t __attribute__((aligned(64))) mask64[] = { -1, -1, -1, -1, 0, 0, 0, 0 };
static int32_t *mask32 = (int32_t *) mask64;

double dot_product_multi(int n, float *a, float *b, float *c, float *d, float *e){
  double result = 0;
  float tresult;
  int res = n & 7;

  __m256d d0, d1, d2, d3, td0, td1, td2, td3;
  __m128  f0, f1, f2, f3, tf0, tf1, tf2, tf3;

  tf0 = _mm_load_ps((float *) (mask32 + 4 + 4 - res) );
  f0 = _mm_load_ps(a);
  f1 = _mm_load_ps(b);
  f0 = _mm_mul_ps(f0, f1);
  _mm_store_ss(&tresult, f0);
  result = tresult;

  return result;
}
