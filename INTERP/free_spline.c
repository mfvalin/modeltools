/** Numerical Analysis 9th ed - Burden, Faires (Ch. 3 Natural Cubic Spline, Pg. 149) */
#include <stdio.h>

void free_variablespace_spline(float *x, float *ax, float *bx, float *cx, float *dx, int n){
  float u[n + 1], z[n + 1], A[n], l[n + 1], h[n];
  float third = 1.0 / 3.0;
  int i, j;
  /** Step 1 */
  for (i = 0; i <= n - 1; ++i) h[i] = x[i + 1] - x[i];

  /** Step 2 */
  for (i = 1; i <= n - 1; ++i)
      A[i] = 3 * (ax[i + 1] - ax[i]) / h[i] - 3 * (ax[i] - ax[i - 1]) / h[i - 1];

  /** Step 3 */
  l[0] = 1;
  u[0] = 0;
  z[0] = 0;

  /** Step 4 */
  for (i = 1; i <= n - 1; ++i) {
      l[i] = 2 * (x[i + 1] - x[i-1]) - h[i - 1] * u[i - 1];
      u[i] = h[i] / l[i];
      z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
  }

  /** Step 5 */
  l[n] = 1;
  z[n] = 0;
  cx[n] = 0;

  /** Step 6 */
  for (j = n - 1; j >= 0; --j) {
      cx[j] = z[j] - u[j] * cx[j + 1];
      bx[j] = (ax[j + 1] - ax[j]) / h[j] - h[j] * (cx[j + 1] + 2 * cx[j]) / 3;
      dx[j] = (cx[j + 1] - cx[j]) / (3 * h[j]);
  }
}

void free_constantspace_spline(float *ax, float *bx, float *cx, float *dx, int n){
  float u[n + 1], z[n + 1];
  float third = 1.0 / 3.0;
  int i, j;
  /** Step 3 */
  u[0] = 0;
  z[0] = 0;

  /** Step 4 */
  for (i = 1; i <= n - 1; ++i) {
      u[i] = 1 / ( 4 - u[i - 1]);
      z[i] = ( 3 * (ax[i + 1] - 2 * ax[i] + ax[i - 1]) - z[i - 1] ) * u[i];
  }

  /** Step 5 */
  z[n] = 0;
  cx[n] = 0;

  /** Step 6 */
  for (j = n - 1; j >= 0; --j) {
      cx[j] = z[j] - u[j] * cx[j + 1];
      bx[j] = (ax[j + 1] - ax[j]) -(cx[j + 1] + 2 * cx[j]) * third;
      dx[j] = (cx[j + 1] - cx[j]) * third;
  }
}

int main() {
    /** Step 0 */
  int n, i, j;
  float fx, deltax, fy;
  scanf("%d", &n);
  n--;
  float x[n + 1], ax[n + 1], h[n], A[n], l[n + 1],
    u[n + 1], z[n + 1], cx[n + 1], bx[n], dx[n],
    ay[n + 1], cy[n + 1], by[n], dy[n];
  for (i = 0; i < n + 1; ++i) scanf("%f",  &ax[i]);    // x coordinate
  for (i = 0; i < n + 1; ++i) scanf("%f", &ay[i]);    // y coordinate
    
  for (i = 0; i < n + 1 ; ++i) x[i] = i * 1.0;

  free_variablespace_spline(x, ax, bx, cx, dx, n);
  fprintf(stderr,"%2s %10s %10s %10s %10s\n", "i", "ai", "bi", "ci", "di");
  for (i = 0; i < n; ++i)
      fprintf(stderr,"%2d %10.5f %10.5f %10.5f %10.5f %10.5f\n", i, ax[i], bx[i], cx[i], dx[i], ax[i]+bx[i]+cx[i]+dx[i]);
  fprintf(stderr,"===========================================================\n");
  // explicit monotonic x increasing by 1
  free_constantspace_spline(ax, bx, cx, dx, n);
  /** Step 7 */
  fprintf(stderr,"%2s %10s %10s %10s %10s\n", "i", "ai", "bi", "ci", "di");
  for (i = 0; i < n; ++i)
      fprintf(stderr,"%2d %10.5f %10.5f %10.5f %10.5f %10.5f\n", i, ax[i], bx[i], cx[i], dx[i], ax[i]+bx[i]+cx[i]+dx[i]);
  fprintf(stderr,"===========================================================\n");
  free_constantspace_spline(ay, by, cy, dy, n);
  for (i = 0; i < n; ++i) {
    for (j = 0 ; j<100 ; j++){
      deltax = j * .01;
      fx = ( ( ( (dx[i] * deltax) + cx[i] ) * deltax ) + bx[i] ) * deltax + ax[i] ;
      fy = ( ( ( (dy[i] * deltax) + cy[i] ) * deltax ) + by[i] ) * deltax + ay[i] ;
      printf(" %10.5f %10.5f\n",fx,fy);
    }
  }
  printf(" %10.5f %10.5f\n",ax[n],ay[n]);
  return 0;
}
