#define _GNU_SOURCE
#include <fenv.h>
#include <stdio.h>

#include <xmmintrin.h>

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
#define FPU_DAZ ((1<<6)|(1<<8))
#define FPU_FTZ ((1<<11)|(1<<15))

static unsigned short translate_x87(unsigned short mode){   // translate "standard" modes into x87 legacy modes
  short new = 0;

  if(mode & FPU_DIVBYZERO)  new |= FE_DIVBYZERO ;
  if(mode & FPU_INEXACT)    new |= FE_INEXACT ;
  if(mode & FPU_INVALID)    new |= FE_INVALID ;
  if(mode & FPU_OVERFLOW)   new |= FE_OVERFLOW ;
  if(mode & FPU_UNDERFLOW)  new |= FE_UNDERFLOW ;
  if((mode & FPU_ALL_EXCEPT) == FPU_ALL_EXCEPT) new |= FE_ALL_EXCEPT ;
  if(mode & FPU_TONEAREST)  new |= FE_TONEAREST ;
  if(mode & FPU_UPWARD)     new |= FE_UPWARD ;
  if(mode & FPU_DOWNWARD)   new |= FE_DOWNWARD ;
  if(mode & FPU_TOWARDZERO) new |= FE_TOWARDZERO ;
  return new;
}

static unsigned short translate_xmm(unsigned short mode){   // translate "standard" modes into xmm/sse modes
  short new = 0;

  if(mode & FPU_DIVBYZERO) new |= (FE_DIVBYZERO << 7) ;
  if(mode & FPU_INEXACT)    new |= (FE_INEXACT << 7) ;
  if(mode & FPU_INVALID)    new |= (FE_INVALID << 7) ;
  if(mode & FPU_OVERFLOW)   new |= (FE_OVERFLOW << 7) ;
  if(mode & FPU_UNDERFLOW)  new |= (FE_UNDERFLOW << 7) ;
  if((mode & FPU_ALL_EXCEPT) == FPU_ALL_EXCEPT) new |= (FE_ALL_EXCEPT << 7) ;
  if(mode & FPU_TONEAREST)  new |= (FE_TONEAREST << 3);
  if(mode & FPU_UPWARD)     new |= (FE_UPWARD << 3);
  if(mode & FPU_DOWNWARD)   new |= (FE_DOWNWARD << 3) ;
  if(mode & FPU_TOWARDZERO) new |= (FE_TOWARDZERO << 3) ;
  if((mode & FPU_DAZ) == FPU_DAZ) new |= FPU_DAZ ;
  if((mode & FPU_FTZ) == FPU_FTZ) new |= FPU_FTZ ;
  return new;
}

static unsigned int set_xmm_control(unsigned short mode){   // disable interrupts by setting bits
  short old_mode, new_mode;

  old_mode = _mm_getcsr();
  old_mode >>= 6;         // get rid of xmm flags
  old_mode <<= 6;
  new_mode = old_mode | mode;
  _mm_setcsr(new_mode);   // set control
  return old_mode;
}
static unsigned int clr_xmm_control(unsigned short mode){   // disable interrupts by setting bits
  short old_mode, new_mode;

  old_mode = _mm_getcsr();
  old_mode >>= 6;         // get rid of xmm flags
  old_mode <<= 6;
  new_mode = old_mode & (~ mode) ;
  _mm_setcsr(new_mode);   // set control
  return old_mode;
}

static unsigned int set_x87_control(unsigned short mode){   // disable interrupts by setting bits
  short old_mode, new_mode;

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

unsigned int set_fpu_control(unsigned short mode){   // disable interrupts by setting bits
  int flags = 0;
  flags = set_x87_control(translate_x87(mode));
  flags = (flags << 16) | set_xmm_control(translate_xmm(mode));
  return flags;
}

unsigned int clr_fpu_control(unsigned short mode){   // disable interrupts by setting bits
  int flags = 0;
  flags = clr_x87_control(translate_x87(mode));
  flags = (flags << 16) | clr_xmm_control(translate_xmm(mode));
  return flags;
}

// unsigned int clr_fpu_flags(){
//   int flags, new;
//   flags = _mm_getcsr();
//   new = flags;
//   new >>= 6;
//   new <<= 6;
//   _mm_setcsr(new);
//   return flags;
// }

#if defined(PRINT_FLAGS)
int main(){
  printf("integer, parameter :: FPU_DIVBYZERO = %d\n",FPU_DIVBYZERO);
  printf("integer, parameter :: FPU_INEXACT = %d\n",FPU_INEXACT);
  printf("integer, parameter :: FPU_INVALID = %d\n",FPU_INVALID);
  printf("integer, parameter :: FPU_OVERFLOW = %d\n",FPU_OVERFLOW);
  printf("integer, parameter :: FPU_UNDERFLOW = %d\n",FPU_UNDERFLOW);
  printf("integer, parameter :: FPU_ALL_EXCEPT = %d\n",FPU_ALL_EXCEPT);
  printf("integer, parameter :: FPU_TONEAREST = %d\n",FPU_TONEAREST);
  printf("integer, parameter :: FPU_UPWARD = %d\n",FPU_UPWARD);
  printf("integer, parameter :: FPU_DOWNWARD = %d\n",FPU_DOWNWARD);
  printf("integer, parameter :: FPU_TOWARDZERO = %d\n",FPU_TOWARDZERO);
  printf("integer, parameter :: FPU_DAZ = %d\n",FPU_DAZ);
  printf("integer, parameter :: FPU_FTZ = %d\n",FPU_FTZ);
  return 0;
}
#endif
#if defined(PRINT_FLAGS_C)
int main(){
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
  printf("#define FPU_DAZ %d\n",FPU_DAZ);
  printf("#define FPU_FTZ %d\n",FPU_FTZ);
  return 0;
}
#endif

#if defined(CODE_FROM_WEB)
// With Visual Studio Proc Pack
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