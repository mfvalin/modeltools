#define _GNU_SOURCE
#include <fenv.h>
#include <stdio.h>

#if defined(__x86_64__)
#include <xmmintrin.h>
#endif

#include <signal.h>
#include <fenv.h>

#define FPU_DIVBYZERO 4
#define FPU_INEXACT 32
#define FPU_INVALID 1
#define FPU_OVERFLOW 8
#define FPU_UNDERFLOW 16
#define FPU_ALL_EXCEPT 61
#define FPU_TONEAREST 0
#define FPU_UPWARD 2048
#define FPU_DOWNWARD 1024
#define FPU_TOWARDZERO 3072
#define FPU_DENORMAL (1<<8)
#define FPU_DAZ ((1<<6)|(1<<8))
#define FPU_FTZ ((1<<11)|(1<<15))

unsigned int disable_fpu_exceptions(unsigned int mode);
unsigned int enable_fpu_exceptions(unsigned int mode);

#if defined(__x86_64__)
static unsigned short translate_x87(unsigned short mode){   // translate "standard" modes into x87 legacy modes
  unsigned short new = 0;

  if(mode & FPU_DIVBYZERO)  new |= FE_DIVBYZERO ;
  if(mode & FPU_INEXACT)    new |= FE_INEXACT ;
  if(mode & FPU_INVALID)    new |= FE_INVALID ;
  if(mode & FPU_OVERFLOW)   new |= FE_OVERFLOW ;
  if(mode & FPU_UNDERFLOW)  new |= FE_UNDERFLOW ;
  if(mode & FPU_TONEAREST)  new |= FE_TONEAREST ;
  if(mode & FPU_UPWARD)     new |= FE_UPWARD ;
  if(mode & FPU_DOWNWARD)   new |= FE_DOWNWARD ;
  if(mode & FPU_TOWARDZERO) new |= FE_TOWARDZERO ;
  return new;
}

static unsigned short translate_xmm(unsigned short mode){   // translate "standard" modes into xmm/sse modes
  unsigned short new = 0;

  if(mode & FPU_DIVBYZERO)  new |= (FE_DIVBYZERO << 7) ;
  if(mode & FPU_INEXACT)    new |= (FE_INEXACT << 7) ;
  if(mode & FPU_INVALID)    new |= (FE_INVALID << 7) ;
  if(mode & FPU_OVERFLOW)   new |= (FE_OVERFLOW << 7) ;
  if(mode & FPU_UNDERFLOW)  new |= (FE_UNDERFLOW << 7) ;
  if(mode & FPU_TONEAREST)  new |= (FE_TONEAREST << 3);
  if(mode & FPU_UPWARD)     new |= (FE_UPWARD << 3);
  if(mode & FPU_DOWNWARD)   new |= (FE_DOWNWARD << 3) ;
  if(mode & FPU_TOWARDZERO) new |= (FE_TOWARDZERO << 3) ;
  if((mode & FPU_DENORMAL))       new |= FPU_DENORMAL ;
  if((mode & FPU_DAZ) == FPU_DAZ) new |= FPU_DAZ ;
  if((mode & FPU_FTZ) == FPU_FTZ) new |= FPU_FTZ ;
  return new;
}

static unsigned int set_xmm_control(unsigned short mode){   // disable interrupts by setting bits
  unsigned int old_mode;
  unsigned int new_mode;

  old_mode = _mm_getcsr();
  old_mode >>= 6;         // get rid of xmm flags
  old_mode <<= 6;
  new_mode = old_mode | mode;
  _mm_setcsr(new_mode);   // set control
  return old_mode;
}
static unsigned int clr_xmm_control(unsigned short mode){   // enable interrupts by clearing bits
  unsigned int old_mode;
  unsigned int new_mode;

  old_mode = _mm_getcsr();
  old_mode >>= 6;         // get rid of xmm flags
  old_mode <<= 6;
  new_mode = old_mode & (~ mode) ;
  _mm_setcsr(new_mode);   // set control
  return old_mode;
}

static unsigned int set_x87_control(unsigned short mode){   // disable interrupts by setting bits
  unsigned short old_mode, new_mode;

  old_mode = 0;
  asm ( "fstcw %w0" : "=m" (old_mode));
  new_mode = old_mode | mode;
  asm ( "fldcw %w0" :: "m" (new_mode));
  return old_mode;
}

static unsigned int clr_x87_control(unsigned short mode){   // enable interrupts by clearing bits
  short old_mode, new_mode;

  old_mode = 0;
  asm ( "fstcw %w0" : "=m" (old_mode));
  new_mode = old_mode & (~ mode) ;
  asm ( "fldcw %w0" :: "m" (new_mode));
  return old_mode;
}
#endif

// Fortran and C callable functions
//    #include <fpu_exceptions.h>       C
//    include "fpu_exceptions.inc"      Fortran
//
//    old_flags = ieee_fp_trap_off()    Fortran and C  (disable floating point traps)
//    old_flags = ieee_fp_trap_on()     Fortran and C  (enable floating point traps)
//    old_flags = ieee_denormal_off()   Fortran and C  (denormals become ZERO)
//    old_flags = ieee_denormal_on()    Fortran and C  (denormals processed again)
//
unsigned int disable_fpu_exceptions(unsigned int modei){   // disable interrupts by setting bits
  unsigned int flags = 0;
  unsigned short mode = modei;
#if defined(__x86_64__)
  flags = set_x87_control(translate_x87(modei));
  flags = (flags << 16) | set_xmm_control(translate_xmm(mode));
#endif
  return flags;
}

unsigned int enable_fpu_exceptions(unsigned int modei){   // enable interrupts by clearing bits
  unsigned int flags = 0;
  unsigned short mode = modei;
#if defined(__x86_64__)
  flags = clr_x87_control(translate_x87(mode));
  flags = (flags << 16) | clr_xmm_control(translate_xmm(mode));
#endif
  return flags;
}

unsigned int ieee_fp_trap_off(){
  unsigned short mode = FPU_ALL_EXCEPT;
  return disable_fpu_exceptions(mode);
}

unsigned int ieee_denormal_off(){
  unsigned short mode = FPU_DAZ+FPU_UNDERFLOW+FPU_FTZ+FPU_INEXACT;
  return disable_fpu_exceptions(mode);
}

unsigned int ieee_fp_trap_on(){
  unsigned short mode = FPU_ALL_EXCEPT;
  return enable_fpu_exceptions(mode);
}

unsigned int ieee_denormal_on(){
  unsigned short mode = (1<<6) | (1<<15);
  unsigned short temp = disable_fpu_exceptions(FPU_UNDERFLOW+FPU_INEXACT);
  return clr_xmm_control(mode);
//   return enable_fpu_exceptions(mode);
}

// to create object file
// cc -c fp_x87_xmm_control
//
// to create include files for Fortran and C
// cc -DPRINT_DEFS fp_x87_xmm_control.c
// ./a.out f >fpu_exceptions.inc   # include file for Fortran (PARAMETER definitions, function prototypes)
// ./a.out c >fpu_exceptions.h     # include file for C (macros + function prototypes)
//
#if defined(PRINT_DEFS)
int main(int argc, char **argv){
  if(argc <2) {
    printf("usage: %s f  (Fortran definitions)\n       %s c  (C definitions)\n",argv[0],argv[0]);
    exit(0);
  }
  if(*argv[1] == 'f'){
    printf("integer(C_INT), parameter :: FPU_DIVBYZERO = %d\n",FPU_DIVBYZERO);
    printf("integer(C_INT), parameter :: FPU_INEXACT = %d\n",FPU_INEXACT);
    printf("integer(C_INT), parameter :: FPU_INVALID = %d\n",FPU_INVALID);
    printf("integer(C_INT), parameter :: FPU_OVERFLOW = %d\n",FPU_OVERFLOW);
    printf("integer(C_INT), parameter :: FPU_UNDERFLOW = %d\n",FPU_UNDERFLOW);
    printf("integer(C_INT), parameter :: FPU_ALL_EXCEPT = %d\n",FPU_ALL_EXCEPT);
    printf("integer(C_INT), parameter :: FPU_TONEAREST = %d\n",FPU_TONEAREST);
    printf("integer(C_INT), parameter :: FPU_UPWARD = %d\n",FPU_UPWARD);
    printf("integer(C_INT), parameter :: FPU_DOWNWARD = %d\n",FPU_DOWNWARD);
    printf("integer(C_INT), parameter :: FPU_TOWARDZERO = %d\n",FPU_TOWARDZERO);
    printf("integer(C_INT), parameter :: FPU_DENORMAL = %d\n",FPU_DENORMAL);
    printf("integer(C_INT), parameter :: FPU_DAZ = %d\n",FPU_DAZ);
    printf("integer(C_INT), parameter :: FPU_FTZ = %d\n",FPU_FTZ);
    printf("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
    "interface",
    "  integer function disable_fpu_exceptions(mode) bind(C,name='disable_fpu_exceptions')",
    "  import :: C_INT",
    "  integer(C_INT), intent(IN), value :: mode",
    "  end function disable_fpu_exceptions",
    "  integer function enable_fpu_exceptions(mode) bind(C,name='enable_fpu_exceptions')",
    "  import :: C_INT",
    "  integer(C_INT), intent(IN), value :: mode",
    "  end function enable_fpu_exceptions",
    "  integer function ieee_fp_trap_off() bind(C,name='ieee_fp_trap_off')",
    "  end function ieee_fp_trap_off",
    "  integer function ieee_fp_trap_on() bind(C,name='ieee_fp_trap_on')",
    "  end function ieee_fp_trap_on",
    "  integer function ieee_denormal_off() bind(C,name='ieee_denormal_off')",
    "  end function ieee_denormal_off",
    "  integer function ieee_denormal_on() bind(C,name='ieee_denormal_on')",
    "  end function ieee_denormal_on",
    "end interface"
    );
  }
  if(*argv[1] == 'c'){
    printf("#define FPU_DIVBYZERO %d\n",FPU_DIVBYZERO);
    printf("#define FPU_INEXACT %d\n",FPU_INEXACT);
    printf("#define FPU_INVALID %d\n",FPU_INVALID);
    printf("#define FPU_OVERFLOW %d\n",FPU_OVERFLOW);
    printf("#define FPU_UNDERFLOW %d\n",FPU_UNDERFLOW);
    printf("#define FPU_ALL_EXCEPT %d\n",FPU_ALL_EXCEPT);
    printf("#define FPU_TONEAREST %d\n",FPU_TONEAREST);
    printf("#define FPU_UPWARD %d\n",FPU_UPWARD);
    printf("#define FPU_DOWNWARD %d\n",FPU_DOWNWARD);
    printf("#define FPU_TOWARDZERO %d\n",FPU_TOWARDZERO);
    printf("#define FPU_DENORMAL %d\n",FPU_DENORMAL);
    printf("#define FPU_DAZ %d\n",FPU_DAZ);
    printf("#define FPU_FTZ %d\n",FPU_FTZ);
    printf("%s\n%s\n%s\n%s\n%s\n%s\n",
    "unsigned int disable_fpu_exceptions(unsigned short mode);",
    "unsigned int enable_fpu_exceptions(unsigned short mode);",
    "unsigned int ieee_fp_trap_off();",
    "unsigned int ieee_fp_trap_on();",
    "unsigned int ieee_denormal_off();",
    "unsigned int ieee_denormal_on();"
    );
  }
  return 0;
}
#endif

#if defined(X86_CODE_FROM_WEB)
#include "xmmintrin.h"
#include "memory.h"

#define X87FLAGBITS         6
#define DAZ_BIT            6
#define FTZ_BIT            15
#define DENORMAL_EXCEPTION_MASK   8
#define UNDERFLOW_EXCEPTION_MASK   11

void set_mxcsr_on(int bit_num)
{
__m128    state[32];
__int32   x;
__asm fxsave   state      
 memcpy( (void*)&x, (char*)state+24, 4);
x |= (1 << bit_num);       
 __asm ldmxcsr   x   
 }

void set_mxcsr_off(int bit_num)
{
__m128    state[32];
__int32   x;
__asm fxsave   state      
 memcpy( (void*)&x, (char*)state+24, 4);
x &= ~(1 << bit_num);
__asm ldmxcsr   x   
 }

void clear_flags()
{
__m128    state[32];
__int32   x;
__asm fxsave   state      
 memcpy( (void*)&x, (char*)state+24, 4);
x = x >> X87FLAGBITS;
x = x << X87FLAGBITS;
__asm ldmxcsr   x   
}

void make_denormal()
{
__m128   denormal;
int      den_vec[4] = {1,1,1,1};
memcpy( &denormal, den_vec, sizeof(int)*4 );
denormal = _mm_add_ps( denormal , denormal );
}

void main()
{
// UNDERFLOWS
set_mxcsr_on(FTZ_BIT);
set_mxcsr_off(UNDERFLOW_EXCEPTION_MASK);
make_denormal();
clear_flags();

// DENORMALS
set_mxcsr_off(DAZ_BIT);
set_mxcsr_on(DENORMAL_EXCEPTION_MASK);
make_denormal();
clear_flags();
}   
#endif