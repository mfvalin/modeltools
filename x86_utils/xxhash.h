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
// transformed into C from C++  M. Valin may 2019
//
#include <stdlib.h>
#include <stdint.h>
/// magic constants :-)
static const uint32_t Prime1 = 2654435761U;
static const uint32_t Prime2 = 2246822519U;
static const uint32_t Prime3 = 3266489917U;
static const uint32_t Prime4 =  668265263U;
static const uint32_t Prime5 =  374761393U;

/// temporarily store up to 15 bytes between multiple add() calls
#define MaxBufferSize 15+1
#define true 1
#define false 0

typedef struct{
  // internal state and temporary buffer
  // state[2] == seed if totalLength < MaxBufferSize
  uint32_t      state[4];
  unsigned char buffer[MaxBufferSize];
  uint32_t      bufferSize;
  uint32_t      totalLength;
} hash32;

/// rotate bits, should compile to a single CPU instruction (ROL)
static inline uint32_t rotateLeft(uint32_t x, unsigned char bits)
{
  return (x << bits) | (x >> (32 - bits));
}

/// process a block of 4x4 bytes, this is the main part of the XXHash32 algorithm
static inline void process(const void* data, uint32_t *state0, uint32_t *state1, uint32_t *state2, uint32_t *state3)
{
  const uint32_t* block = (const uint32_t*) data;
  *state0 = rotateLeft(*state0 + block[0] * Prime2, 13) * Prime1;
  *state1 = rotateLeft(*state1 + block[1] * Prime2, 13) * Prime1;
  *state2 = rotateLeft(*state2 + block[2] * Prime2, 13) * Prime1;
  *state3 = rotateLeft(*state3 + block[3] * Prime2, 13) * Prime1;
}

static inline void XXHash32_init(uint32_t seed, uint32_t bufferSize, uint32_t totalLength, uint32_t *state)
{
  state[0] = seed + Prime1 + Prime2;
  state[1] = seed + Prime2;
  state[2] = seed;
  state[3] = seed - Prime1;
  bufferSize  = 0;
  totalLength = 0;
}

static inline int XXHash32_add(const void* input, uint64_t length, 
			       uint32_t bufferSize, uint32_t totalLength, uint32_t *state, unsigned char *buffer) {
  unsigned int i;
  // no data ?
  if (input == NULL || length == 0)
    return false;

  totalLength += length;
  // byte-wise access
  const unsigned char* data = (const unsigned char*)input;

  // unprocessed old data plus new data still fit in temporary buffer ?
  if (bufferSize + length < MaxBufferSize)
  {
    // just add new data
    while (length-- > 0)
      buffer[bufferSize++] = *data++;
    return true;
  }

  // point beyond last byte
  const unsigned char* stop      = data + length;
  const unsigned char* stopBlock = stop - MaxBufferSize;

  // some data left from previous update ?
  if (bufferSize > 0)
  {
    // make sure temporary buffer is full (16 bytes)
    while (bufferSize < MaxBufferSize)
      buffer[bufferSize++] = *data++;

    // process these 16 bytes (4x4)
    process(buffer, &state[0], &state[1], &state[2], &state[3]);
  }

  // copying state to local variables helps optimizer A LOT
  uint32_t s0 = state[0], s1 = state[1], s2 = state[2], s3 = state[3];
  // 16 bytes at once
  while (data <= stopBlock)
  {
    // local variables s0..s3 instead of state[0]..state[3] are much faster
    process(data, &s0, &s1, &s2, &s3);
    data += 16;
  }
  // copy back
  state[0] = s0; state[1] = s1; state[2] = s2; state[3] = s3;

  // copy remainder to temporary buffer
  bufferSize = stop - data;
  for (i = 0; i < bufferSize; i++)
    buffer[i] = data[i];

  // done
  return true;
}
  
static inline uint32_t XXHash32_hash(uint32_t bufferSize, uint32_t totalLength, uint32_t *state, unsigned char *buffer) {
  uint32_t result = (uint32_t)totalLength;

  // fold 128 bit state into one single 32 bit value
  if (totalLength >= MaxBufferSize)
    result += rotateLeft(state[0],  1) +
	      rotateLeft(state[1],  7) +
	      rotateLeft(state[2], 12) +
	      rotateLeft(state[3], 18);
  else
    // internal state wasn't set in add(), therefore original seed is still stored in state2
    result += state[2] + Prime5;

  // process remaining bytes in temporary buffer
  const unsigned char* data = buffer;
  // point beyond last byte
  const unsigned char* stop = data + bufferSize;

  // at least 4 bytes left ? => eat 4 bytes per step
  for (; data + 4 <= stop; data += 4)
    result = rotateLeft(result + *(uint32_t*)data * Prime3, 17) * Prime4;

  // take care of remaining 0..3 bytes, eat 1 byte per step
  while (data != stop)
    result = rotateLeft(result +        (*data++) * Prime5, 11) * Prime1;

  // mix bits
  result ^= result >> 15;
  result *= Prime2;
  result ^= result >> 13;
  result *= Prime3;
  result ^= result >> 16;
  return result;
}

uint32_t xxhash32(hash32 *h, const void* input, uint64_t length, uint32_t seed);
hash32 *xxhash32_create(void);
hash32 *xxhash32_init(hash32 *h, uint32_t seed);
int xxhash32_add(hash32 *h, const void* input, uint64_t length);
uint32_t xxhash32_finalize(hash32 *h);
uint32_t xxhash32_auto(const void* input, uint64_t length, uint32_t seed);
