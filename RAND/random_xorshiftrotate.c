#if defined(ROTATE)
/*  Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
    minor revision Dec 2017 M.Valin (mfvalin@gmail.com)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>

/* This is the successor to xorshift128+. It is the fastest full-period
   generator passing BigCrush without systematic failures, but due to the
   relatively short period it is acceptable only for applications with a
   mild amount of parallelism; otherwise, use a xorshift1024* generator.

   Beside passing BigCrush, this generator passes the PractRand test suite
   up to (and included) 16TB, with the exception of binary rank tests, as
   the lowest bit of this generator is an LFSR of degree 128. The next bit
   can be described by an LFSR of degree 8256, but in the long run it will
   fail linearity tests, too. The other bits needs a much higher degree to
   be represented as LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   Note that the generator uses a simulated rotate operation, which most C
   compilers will turn into a single instruction. In Java, you can use
   Long.rotateLeft(). In languages that do not make low-level rotation
   instructions accessible xorshift128+ could be faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static uint64_t s[2];
#define RES_SIZE 256
#define RES_MASK 255
static uint32_t res[RES_SIZE] ;
static uint32_t res_index=0;

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}
uint32_t next(void) {       // revision Dec 2017 M.Valin (introduce a "reservoir")
  uint64_t s0;
  uint64_t s1;
  uint64_t result;
  uint32_t result32;
  int i;

  if(res_index == 0){  // buffer is empty
    s0 = s[0];
    s1 = s[1];
    for(i=0 ; i<RES_SIZE2 ; i++){
      result = s0 + s1;
      res64[i] = result;
      s1 ^= s0;
      s0 = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
      s1 = rotl(s1, 36); // c
    }
    s[0] = s0;
    s[1] = s1;
  }
  result32 = res[res_index];
  res_index = (res_index+1) & RES_MASK;
  return result32;
}

static uint64_t next_orig(void) {
  const uint64_t s0 = s[0];
  uint64_t s1 = s[1];
  const uint64_t result = s0 + s1;

  s1 ^= s0;
  s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
  s[1] = rotl(s1, 36); // c

  return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jump(void) {
  static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
    for(int b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
        s0 ^= s[0];
        s1 ^= s[1];
      }
      next();
    }

  s[0] = s0;
  s[1] = s1;
}
#else

#if defined(USE1024)
/*  Written in 2017 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <string.h>

/* NOTE: as of 2017-10-08, this generator has a different multiplier (a
   fixed-point representation of the golden ratio), which eliminates
   linear dependencies from one of the lowest bits. The previous
   multiplier was 1181783497276652981 (M_8 in the paper). If you need to
   tell apart the two generators, you can refer to this generator as
   xorshift1024*phi and to the previous one as xorshift1024*M_8.

   This is a fast, high-quality generator. If 1024 bits of state are too
   much, try a xoroshiro128+ generator.

   Note that the two lowest bits of this generator are LFSRs of degree
   1024, and thus will fail binary rank tests. The other bits needs a much
   higher degree to be represented as LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

uint64_t s[16]; 
int p;

uint64_t next(void) {
  const uint64_t s0 = s[p];
  uint64_t s1 = s[p = (p + 1) & 15];
  s1 ^= s1 << 31; // a
  s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
  return s[p] * 0x9e3779b97f4a7c13;
}


/* This is the jump function for the generator. It is equivalent
   to 2^512 calls to next(); it can be used to generate 2^512
   non-overlapping subsequences for parallel computations. */

void jump(void) {
  static const uint64_t JUMP[] = { 0x84242f96eca9c41d,
    0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
    0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70,
    0xaac17d8efa43cab8, 0xc4cb815590989b13, 0x5ee975283d71c93b,
    0x691548c86c1bd540, 0x7910c41d10a1e6a5, 0x0b5fc64563b3e2a8,
    0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3
  };

  uint64_t t[16] = { 0 };
  for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
    for(int b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b)
        for(int j = 0; j < 16; j++)
          t[j] ^= s[(j + p) & 15];
      next();
    }

  for(int j = 0; j < 16; j++)
    s[(j + p) & 15] = t[j];
}

#else

#if defined(USE64)
/*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>

/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and 
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state; otherwise, we
   rather suggest to use a xoroshiro128+ (for moderately parallel
   computations) or xorshift1024* (for massively parallel computations)
   generator. */

uint64_t x; /* The state can be seeded with any value. */

uint64_t next() {
  uint64_t z = (x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

#else

/*  Written in 2014-2016 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>

/* This generator has been replaced by xoroshiro128plus, which is
   significantly faster and has better statistical properties.

   It might be nonetheless useful for languages in which low-level rotate
   instructions are not available. Due to the relatively short period it
   is acceptable only for applications with a mild amount of parallelism;
   otherwise, use a xorshift1024* generator.

   Note that the lowest bit of this generator is an LFSR of degree 128;
   thus, it will fail linearity tests. The next bit can be described by an
   LFSR of degree 8256, but in the long run it will fail linearity tests,
   too. The other bits needs a much higher degree to be represented as
   LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.

   A previous version of this generator was adding the two halves of the
   newly computed state. This version adds the two halves of the *current*
   state (as xoroshiro128plus does), which improves speed due to better
   internal parallelization from the CPU. The resulting streams are off by
   one step. */

uint64_t s[2];

uint64_t next(void) {
  uint64_t s1 = s[0];
  const uint64_t s0 = s[1];
  const uint64_t result = s0 + s1;
  s[0] = s0;
  s1 ^= s1 << 23; // a
  s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
  return result; 
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jump(void) {
  static const uint64_t JUMP[] = { 0x8a5cd789635d2dff, 0x121fd2155c472f96 };

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++){
    for(int b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
        s0 ^= s[0];
        s1 ^= s[1];
      }
      next();
    }
  }

  s[0] = s0;
  s[1] = s1;
}
#endif

#endif

#endif
