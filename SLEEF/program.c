#define _GNU_SOURCE_xx

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include<math.h>

#include<sleef.h>
#include <functions.h>

#define ROW 500
#define COL 1000000

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

int main(){
  int i, j;
  long long t0,t1, t2,t3,t4;

  float array[COL] ;
  t0 = rdtsc();
  for(i=0 ; i<COL ; i++) array[i] = i / (COL*1.0f) * 3.1415926535f * 2.0f;
  t1 = rdtsc();
  Vsl_sin_f(array, array, COL);
  t2 = rdtsc();
  Vm_sin_f(array, array, COL);
  t3 = rdtsc();
  Vsl_sin_f35(array, array, COL);
  t4 = rdtsc();
  t0 /= 3700 ;
  t1 /= 3700 ;
  t2 /= 3700 ;
  t3 /= 3700 ;
  t4 /= 3700 ;
  printf("setup: %Ld, vsin10: %Ld, sinf: %Ld, vsin35: %Ld\n",t1-t0,t2-t1, t3-t2, t4-t3);
}
