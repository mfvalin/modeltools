#define _GNU_SOURCE

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include<math.h>

#if defined(WITH_INLINE_FUNCTIONS)

#include <string.h>
#include <x86intrin.h>
#define SLEEF_ALWAYS_INLINE inline
#define SLEEF_INLINE static inline
#define SLEEF_CONST const
#define INT_MIN  -2147483648
#define INT_MAX  2147483647
#define DBL_MIN  2.2250738585072014E-308
#define FLT_MIN  1.175494e-38
#include <sleefinline_avx2.h>

#else

#include <sleef.h>

#endif

#if defined(NO_X86)
#define USE_RDTSCP
#undef __x86_64__
#endif

#if defined(NO_SLEEF)
#undef __SLEEF_H__
#endif

#include <functions.h>
