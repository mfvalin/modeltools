#include <stdint.h>
// cc -DSELF_TEST_APICID -I. cpu_type.c
#include <stdio.h>

#include <cpu_type.h>

#define HAS_FLAG(flag,flagword) (flag & flagword)
#define X86_FLAG(flag) (flag & ProcessorCapabilities)

static int ProcessorCapabilities = 0 ;  /* by default, no capabilities are indicated as available */
static uint32_t vstring[12];      /* version string */
static unsigned char *cstring = (char *)vstring;
static uint64_t hz=0;

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

static void get_cpu_capabilities()
{
#if defined(__x86_64__)
  uint32_t regs[4], eax1;;
  int j;
  float freq;

  X86_cpuid( 1, 0, regs );  /* get CPU capabilities EAX=1, ECX=0 */
  eax1 = regs[0];
  if((1 <<  0) & regs[2]) ProcessorCapabilities |= FLAG_SSE3 ;  /* SSE3   ECX bit  0 */
  if((1 << 12) & regs[2]) ProcessorCapabilities |= FLAG_FMA ;   /* FMA    ECX bit 12 */
  if((1 << 20) & regs[2]) ProcessorCapabilities |= FLAG_SSE4 ;  /* SSE4.2 ECX bit 20 */
  if((1 << 28) & regs[2]) ProcessorCapabilities |= FLAG_AVX ;   /* AVX    ECX bit 28 */
  if((1 << 25) & regs[3]) ProcessorCapabilities |= FLAG_SSE ;   /* SSE    EDX bit 25 */
  if((1 << 26) & regs[3]) ProcessorCapabilities |= FLAG_SSE2 ;  /* SSE2   EDX bit 26 */

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
  j = 0 ;
  while(cstring[j++]) ;
  j = j - 3;
  while(cstring[j] != ' ') j--;
  j++;
  sscanf(cstring+j,"%f",&freq); hz = freq*1000000.0 + .5;

  if((ProcessorCapabilities & FLAG_FMA) == 0) return ; /* if FMA flag not present, AVX2 will not be */

  X86_cpuid( 7, 0, regs );    /* get more CPU capabilities EAX=7, ECX=0 */
  if((1 << 5) & regs[1]) ProcessorCapabilities |= FLAG_AVX2 ;   /* AVX2  EBX bit 5 */
  if((1 << 8) & regs[1]) ProcessorCapabilities |= FLAG_BMI  ;   /* BMI2  EBX bit 8 needed to set our BMI flag */
#endif
}

#pragma weak Get_cpu_clock__=Get_cpu_clock
#pragma weak Get_cpu_clock_=Get_cpu_clock
uint64_t Get_cpu_clock__();
uint64_t Get_cpu_clock_();
uint64_t Get_cpu_clock()
{
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return( hz );
}

#pragma weak Cpu_has_feature__=Cpu_has_feature
#pragma weak Cpu_has_feature_=Cpu_has_feature
int Cpu_has_feature__(int flag);
int Cpu_has_feature_(int flag);
int Cpu_has_feature(int flag)
{
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  return( X86_FLAG(flag) );
}

#pragma weak rdtscp__=rdtscp
#pragma weak rdtscp_=rdtscp
uint64_t rdtscp__(void);
uint64_t rdtscp_(void);
uint64_t rdtscp(void) {   // version "in order" avec "serialization"
  uint32_t lo, hi;
#if defined(__x86_64__)
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
#else
  lo = 0;
  hi = 0;
#endif
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#pragma weak rdtsc__=rdtsc
#pragma weak rdtsc_=rdtsc
uint64_t rdtsc__(void);
uint64_t rdtsc_(void);
uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
#if defined(__x86_64__)
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
#else
  lo = 0;
  hi = 0;
#endif
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#if defined(TEST_CPUID)
int get_cpu_core_thread()  /* Intel CPUs only and even in this case not always reliable */
{
  uint32_t regs[4];
  if(ProcessorCapabilities == 0) get_cpu_capabilities();
  X86_cpuid( 0x0B, 0, regs );
  return( ((regs[3]>>regs[0]) << 8) | (regs[3]-((regs[3]>>regs[0])<<regs[0])) );  /* core << 8 + thread  */
}
#endif

#pragma weak Get_cpu_apicid__=Get_cpu_apicid
#pragma weak Get_cpu_apicid_=Get_cpu_apicid
int Get_cpu_apicid__();
int Get_cpu_apicid_();
int Get_cpu_apicid()  /* Intel CPUs only */
{
  uint32_t regs[4];
  X86_cpuid( 0x0B, 0, regs );
  return( regs[3] );  /* x2APIC id from EDX  */
}

#if defined(TEST_CPUID)
int main_cpuid(int argc, char** argv)
{
  int core_and_thread;
#if defined(DEBUG)
  uint32_t regs[4];
#endif
  core_and_thread = get_cpu_core_thread();
  printf("core = %d, thread = %d\n",core_and_thread>>8,core_and_thread&0xFF);
#if defined(DEBUG)
  X86_cpuid( 0x0B, 0, regs );
  X86_cpuid( 0x0B, 1, regs );
  X86_cpuid( 0x0B, 2, regs );
#endif
  printf("CPU speed: %lu Hz\n",Get_cpu_clock());
  printf("FLAGS: ");
  if(Cpu_has_feature(FLAG_SSE))  printf(" SSE");
  if(Cpu_has_feature(FLAG_SSE2)) printf(" SSE2");
  if(Cpu_has_feature(FLAG_SSE3)) printf(" SSE3");
  if(Cpu_has_feature(FLAG_SSE4)) printf(" SSE4");
  if(Cpu_has_feature(FLAG_AVX))  printf(" AVX");
  if(Cpu_has_feature(FLAG_FMA))  printf(" FMA");
  if(Cpu_has_feature(FLAG_AVX2)) printf(" AVX2");
  printf("\n");
  printf("CPU Version string: '%s'\n",cstring);
  return (0);
}
#endif

#if defined(TEST_RDTSC)
int main_rdtsc(int argc, char** argv){   // X86_64 timing demo:  cc -DSELF_TEST rdtscp.c ; ./a.out
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
  for(i=0;i<10;i++){
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
  return(0);
}
#endif
#if defined(TEST_CPUID) || defined(TEST_RDTSC)
int main(int argc, char** argv){
  int status = 0;
#if defined(TEST_RDTSC)
  status = main_rdtsc(argc,argv);
  if(status) return(status);
#endif
#if defined(TEST_CPUID)
  status = main_cpuid(argc,argv);
#endif
  return(status);
}
#endif