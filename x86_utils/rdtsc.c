#include <stdint.h>
#include <stdio.h>

static uint32_t vstring[12];            /* version string */
static unsigned char *cstring = (unsigned char *)vstring;  /* pointer to version string */
static uint64_t hz=0;
static double cycle=0.0;

static uint64_t time0=0;  // pseudo origin for calls to Ix86_tsc, Ix86_tscp, Dx86_tsc, Dx86_tscp

static uint64_t rdtsc(void);
static uint64_t rdtscp(void);

static void X86_cpuid(uint32_t eax, uint32_t ecx, uint32_t* regs)  /* interface to x86 cpuid instruction */
{
#if defined(__x86_64__) || defined( __i386__ )
    uint32_t ebx, edx;
#if defined(DEBUG)
    uint32_t ieax=eax, iecx=ecx;
#endif
# if defined( __i386__ ) && defined ( __PIC__ )
     /* PIC under 32-bit EBX must not be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
   ebx = 0;
    __asm__ ( "cpuid" : "+b" (ebx),
# endif
              "+a" (eax), "+c" (ecx), "=d" (edx) );
    regs[0] = eax; regs[1] = ebx; regs[2] = ecx; regs[3] = edx;
#if defined(DEBUG)
    printf("BEFORE: eax = %8.8x,                 ecx = %8.8x\n",ieax,iecx);
    printf("AFTER : eax = %8.8x, ebx = %8.8x, ecx = %8.8x, edx = %8.8x\n",regs[0],regs[1],regs[2],regs[3]);
#endif
#else
    regs[0] = 0; regs[1] = 0; regs[2] = 0; regs[3] = 0;
#endif
}     
static void init(void){
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t regs[4];
  int j;
  float freq;

  X86_cpuid( 0x80000002, 0, regs );  /* version string (3 calls) */
  for(j=0 ; j<4 ; j++){
    vstring[j] = regs[j];
  }
  X86_cpuid( 0x80000003, 0, regs );
  for(j=0 ; j<4 ; j++){
    vstring[j+4] = regs[j];
  }
  X86_cpuid( 0x80000004, 0, regs );
  for(j=0 ; j<4 ; j++){
    vstring[j+8] = regs[j];
  }
  printf("cstring='%s'\n",cstring);
  j = 0 ;
  while(cstring[j++]) ;
  j = j - 3;
  while(cstring[j] != ' ') j--;
  j++;
  sscanf((const char *__restrict__)(cstring+j),"%f",&freq);
  hz = freq*1000.0 + .5;  // MHz
  hz = hz * 1000000;      // Hz
  cycle = hz;
  cycle = 1.0 / cycle;    // seconds
  time0 = rdtsc();
#endif
}

static uint64_t rdtscp(void) {   // version "in order" avec "serialization"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}

static uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}

uint64_t Ix86_tsc(void) {
  if(time0 == 0) { init ; }
  return (rdtsc()-time0);
}

uint64_t Ix86_tscp(void) {
  if(time0 == 0) { init ; }
  return (rdtscp()-time0);
}

double Dx86_tsc(void) {
  if(time0 == 0) { init ; }
  return  ((int64_t)(rdtsc()-time0)) * cycle ;
}

double Dx86_tscp(void) {
  if(time0 == 0) { init ; }
  return ((int64_t)(rdtscp()-time0)) * cycle ;
}

#if defined(TEST_RDTSC)
#include <stdio.h>
main(){   // X86_64 timing demo:  cc -DSELF_TEST rdtscp.c ; ./a.out
uint64_t t1;
uint64_t t2;
int i;
int ta, tb, tc;
float xc;
double tm1, tm2, tma, tmb ;

init();
t1=rdtscp();
sleep(1);
t2=rdtscp();
t2=(t2-t1)/1000000;  // combien de megaticks dans une seconde
tc=t2;
xc=tc*.001;
printf("hz = %Ld, cycle = %g\n",hz,cycle);
for(i=0;i<50;i++){
  t1=rdtscp();
  t2=rdtscp();
  ta=t2-t1; ta/=xc;   // conversion en nanosecondes
//   t1=rdtscp();
//   sleep(1);
//   t2=rdtscp();
//   t2=(t2-t1)/1000000;
//   tc=t2;
  t1=rdtsc();
  t2=rdtsc();
  tb=t2-t1; tb/=xc;   // conversion en nanosecondes
  printf("%d ns, %d ns, %d ticks/us   -  ",ta,tb,tc);

//   t1=Ix86_tscp();
//   t2=Ix86_tscp();
//   ta=t2-t1; ta/=xc;   // conversion en nanosecondes
  tm1=Dx86_tscp();
  tm2=Dx86_tscp();
  tma=tm2-tm1;          // conversion en nanosecondes
  ta = tma*1.0E+9;
//   t1=Ix86_tscp();
//   sleep(1);
//   t2=Ix86_tscp();
//   t2=(t2-t1)/1000000;
//   tc=t2;
  tm1=Dx86_tsc();
  tm2=Dx86_tsc();
  tmb=tm2-tm1; 
  tb = tmb*1.0E+9;   // conversion en nanosecondes
  printf("%d ns, %d ns, %d ticks/us\n",ta,tb,tc);
}
}
#endif

