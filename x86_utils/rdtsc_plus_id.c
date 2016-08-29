#include <stdlib.h>
#include<stdio.h>

unsigned long long rdtscp(int *socket, int *processor)  // get tsc/socket/processor
{
#if defined(__x86_64__) &&  defined(__linux__)
   unsigned int a, d, c;
   // rdtscp instruction
   // EDX:EAX contain TimeStampCounter
   // ECX contains IA32_TSC_AUX[31:0] (MSR_TSC_AUX value set by OS, lower 32 bits contain socket+processor)
   __asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
    *socket = (c & 0xFFF000)>>12;
    *processor = c & 0xFFF;

   return ((unsigned long long)a) | (((unsigned long long)d) << 32);;
#else
   *socket = 0;
   *processor = 0;
   return 0;
#endif
}
#if defined(SELF_TEST)
int main(int argc, char **argv) {
  int socket, processor;
  socket = -1; processor = 01;
  unsigned long long tsc = rdtscp(&socket, &processor);
  printf("socket = %d, processor = %d \n",socket,processor);
  return 0;
}
#endif
