#include <stdio.h>
#include <stdint.h>
#include <xxhash.h>

static uint64_t rdtsc(void) {   // version "in order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#define SIZE 1024*1024*4
int main(int argc, char **argv){
  int buffer[SIZE];
  int i, j, k;
  uint64_t t0, t1;
  uint32_t result[23], status, length;
  double t;
  hash_t *h=Hash32_create();

  length = SIZE*4;
  for(k = 0 ; k<23 ; k++) {
    fprintf(stderr,"%8d ",length);
    length /= 2;
  }
  fprintf(stderr,"\n");
  for(i=0 ; i<SIZE ; i++) buffer[i] = i + 133;
  for(j=0 ; j<10; j++){
    length = SIZE*4;
    for(k = 0 ; k<23 ; k++){  // divide length by 4 at each iteration
      t0 = rdtsc();
      if(length > 15){
	result[k] = Hash32(12345, (void *)buffer, length);
      }else{
	result[k] = Hash32_short((const char *)buffer, length);
      }
      t1 = rdtsc();
      t = t1-t0 ; t /= length;
      if(length > 15){
	fprintf(stderr,"%8.3f ",t);
      }else{
	fprintf(stderr,"%8.3f:",t);
      }
      length = length/2;
    }
    fprintf(stderr,"\n");
    if(j==0 || j==9) {
      for(k = 0 ; k<23 ; k++) fprintf(stderr,"%8.8x ",result[k]);
      fprintf(stderr,"\n");
    }
    buffer[0] = j*1234;
  }
  return 0;
}

