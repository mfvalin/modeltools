/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2017  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#define _GNU_SOURCE
#include <fenv.h>

static unsigned short fp_trap_status  = 0xFFFF;

#define FP_TRAP_FLAGS (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)

// activate/deactivate the following 3 Floating Point traps: FE_INVALID  FE_DIVBYZERO FE_OVERFLOW
// onoff = 0  : deactivate
// onoff = 1  : activate
// return value : 1 if traps were already active, 0 if traps were not active
//                any other value returned means unknown previous state
int ieee_fp_trap_control (unsigned int onoff)
{
#if defined(__x86_64__)
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
#else
  return 0;
#endif
}

#pragma weak ieee_fp_trap_on__=ieee_fp_trap_on
#pragma weak ieee_fp_trap_on_=ieee_fp_trap_on
int ieee_fp_trap_on__(void);
int ieee_fp_trap_on_(void);
int ieee_fp_trap_on(void)    // activate FE_INVALID  FE_DIVBYZERO FE_OVERFLOW traps
{
  return ieee_fp_trap_control(1);
}

#pragma weak ieee_fp_trap_off__=ieee_fp_trap_off
#pragma weak ieee_fp_trap_off_=ieee_fp_trap_off
int ieee_fp_trap_off__(void);
int ieee_fp_trap_off_(void);
int ieee_fp_trap_off(void)    // deactivate FE_INVALID  FE_DIVBYZERO FE_OVERFLOW traps
{
  return ieee_fp_trap_control(0);
}

#if defined(SELF_TEST)
#include <stdio.h>
int main(void) {
  float a = 1., b = 0.;
  float c;
  int old_status;
  old_status = ieee_fp_trap_control(0) ;
  c = a/b;
  printf("after first zero divide (SUCCESS, this message must appear), c=%G, old_status=%d\n",c,old_status);
   old_status = ieee_fp_trap_control(1) ;
  c = a/b;
  printf("after second zero divide (ERROR, this must not appear), c=%G, old_status=%d\n",c,old_status);
    return 0;
}
#endif
