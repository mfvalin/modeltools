#include <stdlib.h>
#include<stdio.h>

unsigned long tacc_rdtscp(int *chip, int *core)
{
   unsigned long int x;
   unsigned a, d, c;

   __asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
    *chip = (c & 0xFFF000)>>12;
    *core = c & 0xFFF;

   return ((unsigned long)a) | (((unsigned long)d) << 32);;
}
#if defined SELF_TEST
int main(int argc, char **argv) {
  int chip, core;
  chip = -1; core = 01;
  unsigned long tsc = tacc_rdtscp(&chip, &core);
  printf("chip = %d, core = %d \n",chip,core);
  return 0;
}
#endif
