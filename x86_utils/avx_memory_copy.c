/*
 * Copyright (C) 2001-2004 the xine project
 *
 * This file is part of xine, a free video player.
 *
 * xine is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * xine is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
 *
 * These are the MMX/MMX2/SSE optimized versions of memcpy
 *
 * This code was adapted from Linux Kernel sources by Nick Kurshev to
 * the mplayer program. (http://mplayer.sourceforge.net)
 *
 * Miguel Freitas split the #ifdefs into several specialized functions that
 * are benchmarked at runtime by xine. Some original comments from Nick
 * have been preserved documenting some MMX/SSE oddities.
 * Also added kernel memcpy function that seems faster than libc one.
 *
 */


#include <xmmintrin.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define HAVE_SYS_TIMES_H 1
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#else
#include <time.h>
#endif

#define small_memcpy(to,from,n)\
{\
register uintptr_t dummy;\
__asm__ __volatile__(\
  "rep; movsb"\
  :"=&D"(to), "=&S"(from), "=&c"(dummy)\
  :"0" (to), "1" (from),"2" (n)\
  : "memory");\
}

#define LOG_MODULE "memcpy"
#define LOG_VERBOSE
/*
#define LOG
*/
// static void small_memcpy(char * to, const char * from, size_t n){
//   int i;
//   for(i=0 ; i<n ; i++) to[i] = from[i];
// }
/* linux kernel __memcpy (from: /include/asm/string.h) */

static __inline__ void * linux_kernel_memcpy_impl (
			       void * to,
			       const void * from,
			       size_t n)
{
int d0, d1, d2;

  if( n < 4 ) {
    small_memcpy(to,from,n);
  }
  else
    __asm__ __volatile__(
    "rep ; movsl\n\t"
    "testb $2,%b4\n\t"
    "je 1f\n\t"
    "movsw\n"
    "1:\ttestb $1,%b4\n\t"
    "je 2f\n\t"
    "movsb\n"
    "2:"
    : "=&c" (d0), "=&D" (d1), "=&S" (d2)
    :"0" (n/4), "q" (n),"1" ((uintptr_t) to),"2" ((uintptr_t) from)
    : "memory");

  return (to);
}
#define HAVE_AVX 1

#define CACHE_LINE_SIZE 64
#define AVX_MMREG_SIZE 32
#define SSE_MMREG_SIZE 16
#define MMX_MMREG_SIZE 8

#define MMX1_MIN_LEN 0x800  /* 2K blocks */
#define MIN_LEN 0x40  /* 64-byte blocks */

void * fast_avx_mem_set(unsigned char * to, unsigned char what, size_t len)
{
  void *retval;
  size_t i;
  retval = to;
  __m256i ymm0;
  register uintptr_t delta;
// printf("len = %ld\n",len);

  /* Align destinition to MMREG_SIZE -boundary */
  delta = ((uintptr_t)to)&(CACHE_LINE_SIZE-1);
  if(delta){
    delta=CACHE_LINE_SIZE-delta;
    len -= delta;
    for(i=0 ; i<delta ; i++) to[i] = what;   // less than 32 at this point
//     printf("delta = %ld, len = %ld\n",delta,len);
    to += delta;
  }
  i = len >> 7; // 128 byte blocks
  len &= 127;   // remainder
//   printf("len = %ld , i = %ld\n",len,i);
  /* SRC is presumed misaligned */
//   ymm0 = _mm256_xor_si256(ymm0,ymm0);
  ymm0 = _mm256_set1_epi8(what);
  for(; i>0; i--) {    // 128 byte blocks
//     printf("i = %ld, to = %p\n",i,to);
    _mm256_stream_si256 ((__m256i *)&to[ 0], ymm0);                // non temporal store
    _mm256_stream_si256 ((__m256i *)&to[32], ymm0);                // using write combining buffers
    _mm256_stream_si256 ((__m256i *)&to[64], ymm0);
    _mm256_stream_si256 ((__m256i *)&to[96], ymm0);
    to = ((unsigned char *)to) + 128;
  }
// printf("len = %ld\n",len);
  if(len >= 64){    // 64 byte block
    _mm256_stream_si256 ((__m256i *)&to[ 0], ymm0);
    _mm256_stream_si256 ((__m256i *)&to[32], ymm0);
    to = ((unsigned char *)to) + 64;
    len -= 64;
  }
    _mm_sfence();  // _mm256_stream_si256 is weakly-ordered, "sfence" is needed to become ordered again.
  for(i=0 ; i<len ; i++) to[i] = what;  // less than 64 at this point
  return retval;
}

void * fast_avx_mem_clr(unsigned char * to, size_t len)
{
  return ( fast_avx_mem_set(to, (unsigned char) 0, len) );
}

void * fast_avx_dbl2flt(float * to, double * from, int nval)
{
  void *retval;
  size_t i;
  retval = to;
  __m256d ymm0, ymm1, ymm2, ymm3;
  __m128 xmm0, xmm1, xmm2, xmm3;
  register uintptr_t delta;
// printf("nval = %d\n",nval);

  _mm_prefetch (&from[ 0], _MM_HINT_NTA);  // prefetch 3 cache lines
  _mm_prefetch (&from[ 8], _MM_HINT_NTA);
  _mm_prefetch (&from[16], _MM_HINT_NTA);

  /* Align destinition to MMREG_SIZE -boundary */
// printf("to = %p, from = %p\n",to,from);
  delta = ((uintptr_t)to)&(CACHE_LINE_SIZE-1);
  if(delta){
    delta=CACHE_LINE_SIZE-delta;
    delta >>= 2;                            // byte to float
    nval -= delta;
    for(i=0 ; i<delta ; i++) to[i] = from[i];   // less than 16 at this point
//     printf("delta = %ld, nval = %d\n",delta,nval);
    to += delta;
    from += delta;
  }
  i = nval >> 4; // 128 byte blocks (16 floats -> 16 doubles)
  nval &= 15;    // remainder
//   printf("len = %ld , i = %ld\n",len,i);
// printf("to = %p, from = %p\n",to,from);
  /* SRC is presumed misaligned */
  for(; i>0; i--) {    // 128 byte blocks
    _mm_prefetch (&from[24], _MM_HINT_NTA);  // prefetch 2 more cache lines
    _mm_prefetch (&from[32], _MM_HINT_NTA);
    ymm0 = _mm256_loadu_pd (&from[ 0]);         // unaligned load
    ymm1 = _mm256_loadu_pd (&from[ 4]);
    ymm2 = _mm256_loadu_pd (&from[ 8]);
    ymm3 = _mm256_loadu_pd (&from[12]);
    xmm0 = _mm256_cvtpd_ps (ymm0);
    xmm1 = _mm256_cvtpd_ps (ymm1);
    xmm2 = _mm256_cvtpd_ps (ymm2);
    xmm3 = _mm256_cvtpd_ps (ymm3);
    _mm_stream_ps (&to[ 0], xmm0);                // non temporal store
    _mm_stream_ps (&to[ 4], xmm1);                // using write combining buffers
    _mm_stream_ps (&to[ 8], xmm2);
    _mm_stream_ps (&to[12], xmm3);
    from = from + 16;
    to = to + 16;
  }
    _mm_sfence();  // _mm256_stream_si256 is weakly-ordered, "sfence" is needed to become ordered again.
  for(i=0 ; i<nval ; i++) to[i] = from[i];  // less than 16 at this point
  return retval;
}

void * fast_avx_flt2dbl(double * to, float * from, int nval)
{
  void *retval;
  size_t i;
  retval = to;
  __m256d ymm0, ymm1, ymm2, ymm3;
  __m128 xmm0, xmm1, xmm2, xmm3;
  register uintptr_t delta;
// printf("nval = %d\n",nval);

  _mm_prefetch (&from[ 0], _MM_HINT_NTA);  // prefetch 3 cache lines
  _mm_prefetch (&from[16], _MM_HINT_NTA);
  _mm_prefetch (&from[32], _MM_HINT_NTA);

  /* Align destinition to MMREG_SIZE -boundary */
// printf("to = %p, from = %p\n",to,from);
  delta = ((uintptr_t)to)&(CACHE_LINE_SIZE-1);
  if(delta){
    delta=CACHE_LINE_SIZE-delta;
    delta >>= 3;                            // byte to double
    nval -= delta;
    for(i=0 ; i<delta ; i++) to[i] = from[i];   // less than 16 at this point
//     printf("delta = %ld, nval = %d\n",delta,nval);
    to += delta;
    from += delta;
  }
  i = nval >> 4; // 128 byte blocks (16 floats -> 16 doubles)
  nval &= 15;    // remainder
//   printf("len = %ld , i = %ld\n",len,i);
// printf("to = %p, from = %p\n",to,from);
  /* SRC is presumed misaligned */
  for(; i>0; i--) {    // 128 byte blocks
    _mm_prefetch (&from[48], _MM_HINT_NTA);  // prefetch 1 more cache lines
    xmm0 = _mm_loadu_ps (&from[ 0]);         // unaligned load
    xmm1 = _mm_loadu_ps (&from[ 4]);
    xmm2 = _mm_loadu_ps (&from[ 8]);
    xmm3 = _mm_loadu_ps (&from[12]);
    ymm0 = _mm256_cvtps_pd (xmm0);
    ymm1 = _mm256_cvtps_pd (xmm1);
    ymm2 = _mm256_cvtps_pd (xmm2);
    ymm3 = _mm256_cvtps_pd (xmm3);
    _mm256_stream_pd (&to[ 0], ymm0);                // non temporal store
    _mm256_stream_pd (&to[ 4], ymm1);                // using write combining buffers
    _mm256_stream_pd (&to[ 8], ymm2);
    _mm256_stream_pd (&to[12], ymm3);
    from = from + 16;
    to = to + 16;
  }
    _mm_sfence();  // _mm256_stream_si256 is weakly-ordered, "sfence" is needed to become ordered again.
  for(i=0 ; i<nval ; i++) to[i] = from[i];  // less than 16 at this point
  return retval;
}

void * fast_avx_memcpy(unsigned char * to, const unsigned char * from, size_t len)
{
  void *retval;
  size_t i;
  retval = to;
  __m256i ymm0, ymm1, ymm2, ymm3;
  register uintptr_t delta;
// printf("len = %ld\n",len);

  _mm_prefetch (&from[  0], _MM_HINT_NTA);  // prefetch 5 cache lines
  _mm_prefetch (&from[ 64], _MM_HINT_NTA);
  _mm_prefetch (&from[128], _MM_HINT_NTA);
  _mm_prefetch (&from[192], _MM_HINT_NTA);
  _mm_prefetch (&from[256], _MM_HINT_NTA);

  /* Align destinition to MMREG_SIZE -boundary */
  delta = ((uintptr_t)to)&(CACHE_LINE_SIZE-1);
  if(delta){
    delta=CACHE_LINE_SIZE-delta;
    len -= delta;
//     small_memcpy(to,from,delta);
    for(i=0 ; i<delta ; i++) to[i] = from[i];   // less than 64 at this point
    to += delta;
    from += delta;
//     printf("delta = %ld, len = %ld\n",delta,len);
  }
  i = len >> 7; // 128 byte blocks
  len &= 127;   // remainder
//   printf("len = %ld , i = %ld\n",len,i);
  /* SRC is presumed misaligned */
  for(; i>0; i--) {    // 128 byte blocks
    _mm_prefetch (&from[320], _MM_HINT_NTA);  // prefetch 2 more cache lines
    _mm_prefetch (&from[384], _MM_HINT_NTA);
    ymm0 = _mm256_loadu_si256 ((__m256i const *)&from[ 0]);       // unaligned load
    ymm1 = _mm256_loadu_si256 ((__m256i const *)&from[32]);
    ymm2 = _mm256_loadu_si256 ((__m256i const *)&from[64]);
    ymm3 = _mm256_loadu_si256 ((__m256i const *)&from[96]);
    _mm256_stream_si256 ((__m256i *)&to[ 0], ymm0);                // non temporal store
    _mm256_stream_si256 ((__m256i *)&to[32], ymm1);                // using write combining buffers
    _mm256_stream_si256 ((__m256i *)&to[64], ymm2);
    _mm256_stream_si256 ((__m256i *)&to[96], ymm3);
    from = ((const unsigned char *)from) + 128;
    to = ((unsigned char *)to) + 128;
  }
// printf("len = %ld\n",len);
  if(len >= 64){    // 64 byte block
    ymm0 = _mm256_loadu_si256 ((__m256i const *)&from[ 0]);
    ymm1 = _mm256_loadu_si256 ((__m256i const *)&from[32]);
    _mm256_stream_si256 ((__m256i *)&to[ 0], ymm0);
    _mm256_stream_si256 ((__m256i *)&to[32], ymm1);
    from = ((const unsigned char *)from) + 64;
    to = ((unsigned char *)to) + 64;
    len -= 64;
  }
  _mm_sfence();  // _mm256_stream_si256 is weakly-ordered, "sfence" is needed to become ordered again.
//   small_memcpy(to,from,len);
  for(i=0 ; i<len ; i++) to[i] = from[i];  // less than 64 at this point
  return retval;
}

static void * avx_memcpy_ori(void * to, const void * from, size_t len)
{
  void *retval;
  size_t i;
  retval = to;

  /* PREFETCH has effect even for MOVSB instruction ;) */
  __asm__ __volatile__ (
    "   prefetchnta (%0)\n"
    "   prefetchnta 32(%0)\n"
    "   prefetchnta 64(%0)\n"
    "   prefetchnta 96(%0)\n"
    "   prefetchnta 128(%0)\n"
    "   prefetchnta 160(%0)\n"
    "   prefetchnta 192(%0)\n"
    "   prefetchnta 224(%0)\n"
    "   prefetchnta 256(%0)\n"
    "   prefetchnta 288(%0)\n"
    : : "r" (from) );

  if(len >= MIN_LEN)
  {
    register uintptr_t delta;
    /* Align destinition to MMREG_SIZE -boundary */
    delta = ((uintptr_t)to)&(AVX_MMREG_SIZE-1);
    if(delta)
    {
      delta=AVX_MMREG_SIZE-delta;
      len -= delta;
      small_memcpy(to, from, delta);
    }
    i = len >> 7; /* len/128 */
    len&=127;
    if(((uintptr_t)from) & 31)
      /* if SRC is misaligned */
      for(; i>0; i--)
      {
        __asm__ __volatile__ (
        "prefetchnta 320(%0)\n"
        "prefetchnta 352(%0)\n"
        "prefetchnta 384(%0)\n"
        "prefetchnta 416(%0)\n"
        "vmovups    (%0), %%ymm0\n"
        "vmovups  32(%0), %%ymm1\n"
        "vmovups  64(%0), %%ymm2\n"
        "vmovups  96(%0), %%ymm3\n"
        "vmovntps %%ymm0,   (%1)\n"
        "vmovntps %%ymm1, 32(%1)\n"
        "vmovntps %%ymm2, 64(%1)\n"
        "vmovntps %%ymm3, 96(%1)\n"
        :: "r" (from), "r" (to) : "memory");
        from = ((const unsigned char *)from) + 128;
        to = ((unsigned char *)to) + 128;
      }
     /* since movntq is weakly-ordered, a "sfence"
     * is needed to become ordered again. */
    __asm__ __volatile__ ("sfence":::"memory");
  }
  /*
   *	Now do the tail of the block
   */
  if(len) linux_kernel_memcpy_impl(to, from, len);
  return retval;
}

#if defined(TEST)
static inline uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
static inline uint64_t rdtscp()  // get tsc/socket/processor
{
   unsigned int a, d, c;
   // rdtscp instruction
   // EDX:EAX contain TimeStampCounter
   // ECX contains IA32_TSC_AUX[31:0] (MSR_TSC_AUX value set by OS, lower 32 bits contain socket+processor)
   __asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));

   return ((uint64_t)a) | (((uint64_t)d) << 32);;
}
#define NPTS 1024*1024*128
#define NREP 10

int main(int argc, char **argv){
  int src[NPTS];
  float srcf[NPTS];
  int dst[NPTS];
  double dstd[NPTS];
  unsigned char *srcc, *dstc;
  void *dummy;
  int i;
  unsigned char what = 193;
//   uint64_t tm1[NREP], tm2[NREP], tm3[NREP], tm4[NREP], tm5[NREP], tm6[NREP] ;
  uint64_t t0, t1, t2, t3, t4, t5;
  double f1, f2, f3, f4, f5;
  double clock = 3.7;

  if(argc > 1)  clock = atof(argv[1]);
  srcc = (unsigned char *) &src[0];
  dstc = (unsigned char *) &dst[0];

  for(i=0 ; i<NPTS ; i++) { src[i] = i ; dst[i] = 0 ; }
  src[0] = -1; src[NPTS-1] = (NPTS - 1) | (255 <<24);
  for(i=0 ; i<NPTS ; i++) { srcf[i] = (i+.5) ; dstd[i] = 0 ; }
  fast_avx_flt2dbl(&dstd[0], &srcf[0], NPTS);
  srcf[0] = -1.0 ; srcf[NPTS-1] = -1.0;

  for (i=0 ; i<NREP ; i++){
    t0 = rdtsc();
    dummy = fast_avx_memcpy(&dstc[3], &srcc[3], sizeof(src) - 6);  // MSB of last element expected to be missing
    t1 = rdtsc();
    dummy = fast_avx_flt2dbl(&dstd[1], &srcf[1], NPTS-2);
    t2 = rdtsc();
    dummy = fast_avx_dbl2flt(&srcf[1], &dstd[1], NPTS-2);
    t3 = rdtsc();
    dummy = fast_avx_mem_clr((unsigned char *) &dst[1], sizeof(src) - 6);
    t4 = rdtsc();
    dummy = fast_avx_mem_set((unsigned char *) &dstc[1], what, sizeof(src) - 2);
    t5 = rdtsc();
//     printf("copy = %10ld, flt2dbl=%10ld, dbl2flt=%10ld, clr=%10ld, set=%10ld\n",t1-t0,t2-t1,t3-t2,t4-t3,t5-t4);
//     printf("copy = %10d, flt2dbl=%10d, dbl2flt=%10d, clr=%10d, set=%10d\n",NPTS*8,NPTS*12,NPTS*12,NPTS*4,NPTS*4);
    f1 = NPTS*8 ; f2 = NPTS*12 ; f3 = f2 ; f4 = NPTS*4 ; f5 = f4;
    f1 /= (t1-t0) ; f2 /= (t2-t1) ; f3 /= (t3-t2) ; f4 /= (t4-t3) ; f5 /= (t5-t4) ; 
    f1 *= clock ; f2 *= clock ; f3 *= clock ; f4 *= clock ; f5 *= clock ; 
    printf("copy = %7.3f, flt2dbl=%7.3f, dbl2flt=%7.3f, clr=%7.3f, set=%7.3f  (GBytes/sec at %3.1f GHz)\n",f1,f2,f3,f4,f5,clock);
    if(clock == 1.23456) {   // useless code to get rid of gcc warnings
      dummy = avx_memcpy_ori(&dstc[0], &srcc[0], (size_t) 0);
      printf("%p\n",dummy);
    }
  }

//   for(i=0 ; i<NPTS ; i++) {
//     if(src[i] != dst[i]){
//       printf("i = %d, src[i] = %d, dst[i] = %d\n",i, src[i], dst[i]);
//     }
//   }
  printf(" %15f %15f %15f %15f\n",srcf[0],srcf[1],srcf[NPTS - 2],srcf[NPTS - 1]);
  printf(" %15f %15f %15f %15f\n",dstd[0],dstd[1],dstd[NPTS - 2],dstd[NPTS - 1]);

  for(i=0 ; i<NPTS ; i++) { src[i] = i ; dst[i] = 0 ; }
  src[0] = -1; src[NPTS-1] = (NPTS - 1) | (255 <<24);
  dummy = fast_avx_memcpy(&dstc[3], &srcc[3], sizeof(src) - 5);
  printf("\n");
  printf(" %8.8x %8.8x %8.8x %8.8x\n",src[0],src[1],src[NPTS - 2],src[NPTS - 1]);
  printf(" %8.8x %8.8x %8.8x %8.8x\n",dst[0],dst[1],dst[NPTS - 2],dst[NPTS - 1]);

  dummy = fast_avx_mem_clr((unsigned char *) &dst[0], sizeof(src));
  dummy = fast_avx_mem_set((unsigned char *) &dstc[2], what, sizeof(src) - 3);
  printf("\n");
  printf(" %8.8x %8.8x %8.8x %8.8x\n",src[0],src[1],src[NPTS - 2],src[NPTS - 1]);
  printf(" %8.8x %8.8x %8.8x %8.8x\n",dst[0],dst[1],dst[NPTS - 2],dst[NPTS - 1]);
  return (0) ;
}
#endif

