#define _GNU_SOURCE
#include <fenv.h>

static unsigned short fp_trap_status  = 0xFFFF;

#define FP_TRAP_FLAGS (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)

#pragma weak x86_fp_trap_control__=x86_fp_trap_control
#pragma weak x86_fp_trap_control_=x86_fp_trap_control
int x86_fp_trap_control__ (unsigned int onoff);
int x86_fp_trap_control_ (unsigned int onoff);
int x86_fp_trap_control (unsigned int onoff)
{
  unsigned short int x87_status, sse_status;
  unsigned short int old_status = fp_trap_status ;

  if(onoff > 1) return 0xFFFF;  // ERROR

  if(fp_trap_status == onoff) return fp_trap_status ;  // already set properly

  __asm__ ("fstcw %0" : "=m" (*&x87_status));   // get legacy X87 FPU control word

  if(onoff == 0) {                    // disable
    x87_status |= FP_TRAP_FLAGS;      // set flag bits
  }else{                              // enable
    x87_status &= (~FP_TRAP_FLAGS);   // clear flag bits
  }

  __asm__ ("fldcw %0" : : "m" (*&x87_status));   // set legacy X87 FPU control word

  __asm__ ("stmxcsr %0" : "=m" (*&sse_status));  // get new SSE MXCSR register

  /* The SSE exception masks are shifted left by 7 bits.  */
  if(onoff == 0) {                           // disable
    sse_status |= (FP_TRAP_FLAGS << 7);      // set flag bits
  }else{                                     // enable
    sse_status &= ~(FP_TRAP_FLAGS << 7);     // clear flag bits
  }

  __asm__ ("ldmxcsr %0" : : "m" (*&sse_status));  // set new SSE MXCSR register

  fp_trap_status = onoff;
  return old_status;
}

#if defined(SELF_TEST)
#include <stdio.h>
int main(void) {
  float a = 1., b = 0.;
  float c;
  int old_status;
  old_status = x86_fp_trap_control(0) ;
  c = a/b;
  printf("after first zero divide (SUCCESS, this message must appear), c=%G, old_status=%d\n",c,old_status);
   old_status = x86_fp_trap_control(1) ;
  c = a/b;
  printf("after second zero divide (ERROR, this must not appear), c=%G, old_status=%d\n",c,old_status);
    return 0;
}
#endif
