#include <stdint.h>

/*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */


/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and 
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state; otherwise, we
   rather suggest to use a xoroshiro128+ (for moderately parallel
   computations) or xorshift1024* (for massively parallel computations)
   generator. */

static uint64_t x = 123456789 ; /* The state can be seeded with any value. */

static inline uint64_t rani() {
  uint64_t z = (x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

/*  Added in 2017 by M.Valin */
static const double INVM31   =  4.65661287307739257812e-010 ;        /* 1.0 / 2^31 */
static const double INVM32   =  2.32830643653869628906e-010 ;        /* 1.0 / 2^32 */
#define CVTDBLS_32(i1)      ((int)(i1)  * INVM31 + (INVM32))

// introduce a "fuzz" in some values of a 32 bit floating point array
//
// for the time being only one random generator is used (see above)
// f       : 32bit float array to be "fuzzed"
// n       : size of array f
// start   : starting position (origin 0) in array f
// stride  : increment between points to "fuzz"
// mask    : if non zero, xor mask with selected points of array f to flip bits
// epsilon : if non zero, multiply selected points of array f by (1.0 + epsilon)
// mode    : if non zero, perform a randomized flip or scaling on selected points of array f
//
// in the future different non zero values of mode will be used to select a random generator
#pragma weak float_fuzz__=float_fuzz
#pragma weak float_fuzz_=float_fuzz
void float_fuzz__(float *f, int n, int start, int stride, int mask, float epsilon, int mode);
void float_fuzz_(float *f, int n, int start, int stride, int mask, float epsilon, int mode);
void float_fuzz(float *f, int n, int start, int stride, int mask, float epsilon, int mode)
{
  int i;
  uint32_t j;
  uint32_t *jf = (uint32_t *) f;

  if(mode == 0){                 // flip bits or multiply by 1 + epsilon
    if(mask != 0){                   // flip bits using mask
      for(i=start ; i<n ; i+=stride){
	jf[i] = jf[i] ^ mask;
      }
    }else if(epsilon != 0.0){        // scale by 1 + epsilon
      for(i=start ; i<n ; i+=stride){
	f[i] = f[i] * (1.0 + epsilon);
      }
    }
  }else{                             // randomized behavior
    if(mask != 0){
      for(i=start ; i<n ; i+=stride){  // randomly flip (or not) bits according to mask
	j = rani() ;
	jf[i] = jf[i] ^ (mask & j);
      }
    }else if(epsilon != 0.0){          // random scaling by : 1 - epsilon < random factor < 1 + epsilon
      for(i=start ; i<n ; i+=stride){
	j = rani();
	f[i] = f[i] * (1.0 + epsilon * CVTDBLS_32(j) );   // multiply by 1 + random*epsilon
      }
    }
  }
}
