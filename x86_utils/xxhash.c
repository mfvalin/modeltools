// MIT License
// 
// Copyright (c) 2018 Stephan Brumme
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the Software
// is furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.#include <stdlib.h>
//
// C++ to C adaptation, added a routine for short lengths  M. Valin may 2019
//
#include <stdlib.h>
#include <stdint.h>
#include <xxhash.h>

#if defined(TIMING)
uint64_t rdtsc(void) {   // version "in order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}
#endif

uint32_t Hash32_h(hash_t *h, uint32_t seed, const void* input, uint64_t length)
{

  XXHash32_init(h, seed);
  XXHash32_add(h, input, length);
  return XXHash32_hash(h->bufferSize, h->totalLength, h->state, h->buffer);
}

hash_t *Hash32_create(void){
  hash_t *h = (hash_t *) malloc(sizeof(hash_t));
  return h;
}

void Hash32_init(hash_t *h, uint32_t seed){
  XXHash32_init(h, seed);
}

uint32_t Hash32_add(hash_t *h, const void* input, uint64_t length){
  return  XXHash32_add(h, input, length);
}

uint32_t Hash32_finalize(hash_t *h){
  return XXHash32_hash(h->bufferSize, h->totalLength, h->state, h->buffer);
}

uint32_t Hash32(uint32_t seed, const void* input, uint32_t length){
  uint32_t result;
  hash_t h;

//   if(length < 16) return xxhash_short((const char *) input, length);

  XXHash32_init(&h, seed);
  XXHash32_add(&h, input, length);
  result = XXHash32_hash(h.bufferSize, h.totalLength, h.state, h.buffer);
  return result;
}

// for the following macros come from uthash, see below
/*
Copyright (c) 2003-2018, Troy D. Hanson     http://troydhanson.github.com/uthash/
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#define HASH_BER(key,keylen,hashv)                                               \
do {                                                                             \
  uint32_t _hb_keylen = (uint32_t)keylen;                                        \
  const unsigned char *_hb_key = (const unsigned char*)(key);                    \
  (hashv) = 0;                                                                   \
  while (_hb_keylen-- != 0U) {                                                   \
    (hashv) = (((hashv) << 5) + (hashv)) + *_hb_key++;                           \
  }                                                                              \
} while (0)

#define HASH_SAX(key,keylen,hashv)                                               \
do {                                                                             \
  unsigned _sx_i;                                                                \
  const unsigned char *_hs_key = (const unsigned char*)(key);                    \
  hashv = 0;                                                                     \
  for (_sx_i=0; _sx_i < keylen; _sx_i++) {                                       \
    hashv ^= (hashv << 5) + (hashv >> 2) + _hs_key[_sx_i];                       \
  }                                                                              \
} while (0)
/* FNV-1a variation */
#define HASH_FNV(key,keylen,hashv)                                               \
do {                                                                             \
  unsigned _fn_i;                                                                \
  const unsigned char *_hf_key = (const unsigned char*)(key);                    \
  (hashv) = 2166136261U;                                                         \
  for (_fn_i=0; _fn_i < keylen; _fn_i++) {                                       \
    hashv = hashv ^ _hf_key[_fn_i];                                              \
    hashv = hashv * 16777619U;                                                   \
  }                                                                              \
} while (0)

#define HASH_OAT(key,keylen,hashv)                                               \
do {                                                                             \
  unsigned _ho_i;                                                                \
  const unsigned char *_ho_key=(const unsigned char*)(key);                      \
  hashv = 0;                                                                     \
  for(_ho_i=0; _ho_i < keylen; _ho_i++) {                                        \
      hashv += _ho_key[_ho_i];                                                   \
      hashv += (hashv << 10);                                                    \
      hashv ^= (hashv >> 6);                                                     \
  }                                                                              \
  hashv += (hashv << 3);                                                         \
  hashv ^= (hashv >> 11);                                                        \
  hashv += (hashv << 15);                                                        \
} while (0)

uint32_t Hash32_ber(const char * data, int len) {
  uint32_t result = 123456;
  HASH_BER(data, len, result);
  return result;
}

uint32_t Hash32_sax(const char * data, int len) {
  uint32_t result = 123456;
  HASH_SAX(data, len, result);
  return result;
}

uint32_t Hash32_fnv(const char * data, int len) {
  uint32_t result = 123456;
  HASH_FNV(data, len, result);
  return result;
}

uint32_t Hash32_oat(const char * data, int len) {
  uint32_t result = 123456;
  HASH_OAT(data, len, result);
  return result;
}


// use the FNV-1a for short byte sequences
#pragma weak Hash32_short = Hash32_fnv
uint32_t Hash32_short(const char * data, int len);
