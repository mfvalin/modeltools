#include <stdint.h>
#pragma weak rdtscp_=rdtscp
uint64_t rdtscp_(void);
uint64_t rdtscp(void) {   // version "in order" avec "serialization"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
#if defined(TEST_RDTSC)
#include <stdio.h>
main(){   // X86_64 timing demo:  cc -DSELF_TEST rdtscp.c ; ./a.out
uint64_t t1;
uint64_t t2;
int i;
int ta, tb, tc;
float xc;
t1=rdtscp();
sleep(1);
t2=rdtscp();
t2=(t2-t1)/1000000;  // combien de megaticks dans une seconde
tc=t2;
xc=tc*.001;
for(i=0;i<50;i++){
  t1=rdtscp();
  t2=rdtscp();
  ta=t2-t1; ta/=xc;   // conversion en nanosecondes
  t1=rdtscp();
  sleep(1);
  t2=rdtscp();
  t2=(t2-t1)/1000000;
  tc=t2;
  t1=rdtsc();
  t2=rdtsc();
  tb=t2-t1; tb/=xc;   // conversion en nanosecondes
  printf("%d ns, %d ns, %d ticks/us\n",ta,tb,tc);
}
}
#endif

