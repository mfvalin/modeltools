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
// C++ to C adaptation:  M. Valin may 2019
//
#include <stdlib.h>
#include <stdint.h>

// magic constants :-)
static const uint32_t Prime1 = 2654435761U;
static const uint32_t Prime2 = 2246822519U;
static const uint32_t Prime3 = 3266489917U;
static const uint32_t Prime4 =  668265263U;
static const uint32_t Prime5 =  374761393U;

/// temporarily store up to MaxBufferSize-1 bytes between multiple add() calls
#define MaxBufferSize 16
#define true 1
#define false 0

typedef struct{
  // internal state and temporary buffer
  // state[2] == seed if totalLength < MaxBufferSize
  uint32_t      state[4];
  unsigned char buffer[MaxBufferSize];
  uint32_t      bufferSize;
  uint32_t      totalLength;
} hash_t;

/// rotate bits, should compile to a single CPU instruction (ROL)
static inline uint32_t rotateLeft(uint32_t x, unsigned char bits)
{
  return (x << bits) | (x >> (32 - bits));
}

/// process a block of 4x4 bytes, this is the main part of the XXHash32 algorithm
// static inline void process(const void* data, uint32_t *state0, uint32_t *state1, uint32_t *state2, uint32_t *state3)
// {
//   const uint32_t* block = (const uint32_t*) data;
//   *state0 = rotateLeft(*state0 + block[0] * Prime2, 13) * Prime1;
//   *state1 = rotateLeft(*state1 + block[1] * Prime2, 13) * Prime1;
//   *state2 = rotateLeft(*state2 + block[2] * Prime2, 13) * Prime1;
//   *state3 = rotateLeft(*state3 + block[3] * Prime2, 13) * Prime1;
// }
static inline void processv4(const void* data, uint32_t *state) // process a 16 byte block
{
  const uint32_t* block = (const uint32_t*) data;
  int i;
  for(i=0 ; i<4 ; i++) state[i] = rotateLeft(state[i] + block[i] * Prime2, 13) * Prime1;
}

static inline void XXHash32_init(hash_t *h, uint32_t seed)
{
  h->state[0] = seed + Prime1 + Prime2;
  h->state[1] = seed + Prime2;
  h->state[2] = seed;
  h->state[3] = seed - Prime1;
  h->bufferSize  = 0;
  h->totalLength = 0;
}

static inline int XXHash32_add(hash_t *h, const void* input, uint64_t length) {
  unsigned int i;
  uint32_t bufferSize = h->bufferSize;
  uint32_t *state = h->state;
  unsigned char *buffer = h->buffer;

  // no data ?
  if (input == NULL || length == 0)
    return false;

  h->totalLength += length;
  // byte-wise access
  const unsigned char* data = (const unsigned char*)input;

  // unprocessed old data plus new data still fit in temporary buffer ?
  if (bufferSize + length < MaxBufferSize) {    // just add new data
    while (length-- > 0)  buffer[bufferSize++] = *data++;
    h->bufferSize = bufferSize;
    return true;
  }
  // point beyond last byte
  const unsigned char* stop      = data + length;
  const unsigned char* stopBlock = stop - MaxBufferSize;

  // some data left from previous update ?
  if (bufferSize > 0) {    // make sure temporary buffer is full (MaxBufferSize bytes)
    while (bufferSize < MaxBufferSize)  buffer[bufferSize++] = *data++;
    processv4(buffer, state);    // process these MaxBufferSize bytes
  }

  while (data <= stopBlock) {   // process MaxBufferSize bytes at a time
    processv4(data, state);
    data += MaxBufferSize;
  }
  // copy remainder to temporary buffer
  bufferSize = stop - data;
  for (i = 0; i < bufferSize; i++)
    buffer[i] = data[i];

  h->bufferSize = bufferSize;
  return true;  // done
}
  
static inline uint32_t XXHash32_hash(uint32_t bufferSize, uint32_t totalLength, uint32_t *state, unsigned char *buffer) {
  uint32_t result = (uint32_t)totalLength;

  // fold 128 bit state into one single 32 bit value
  if (totalLength >= MaxBufferSize)
    result += rotateLeft(state[0],  1) + rotateLeft(state[1],  7) + rotateLeft(state[2], 12) + rotateLeft(state[3], 18);
  else
    result += state[2] + Prime5;  // internal state wasn't set in add(), therefore original seed is still stored in state2

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

  // mix result bits
  result ^= result >> 15;
  result *= Prime2;
  result ^= result >> 13;
  result *= Prime3;
  result ^= result >> 16;
  return result;
}

// exported function definitions
uint32_t Hash32_h(hash_t *h, uint32_t seed, const void* input, uint64_t length);
hash_t   *Hash32_create(void);
void     Hash32_init(hash_t *h, uint32_t seed);
uint32_t Hash32_add(hash_t *h, const void* input, uint64_t length);
uint32_t Hash32_finalize(hash_t *h);
uint32_t Hash32(uint32_t seed, const void* input, uint32_t length);
uint32_t Hash32_short(const char * data, int len);
#if defined(MUST_NEVER_EVER_BE_TRUE)
// Fortran interfaces
interface                                                                 !InTf
  function Hash32(seed, input, length) result(hash) BIND(C,name='Hash32') !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: input                               !InTf
    integer(C_INT), intent(IN), value :: seed, length                     !InTf
    integer(C_INT) :: hash                                                !InTf
  end function Hash32                                                     !InTf
  function Hash32_short(input, length) result(hash) BIND(C,name='Hash32_short') !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: input                               !InTf
    integer(C_INT), intent(IN), value :: length                           !InTf
    integer(C_INT) :: hash                                                !InTf
  end function Hash32_short                                               !InTf
  function Hash32_create() result(hash_struct) BIND(C,name='Hash32_create') !InTf
    import :: C_PTR                                                       !InTf
    type(C_PTR) :: hash_struct                                            !InTf
  end function Hash32_create                                              !InTf
  subroutine Hash32_init(hash_struct, seed) BIND(C,name='Hash32_init')    !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: hash_struct                         !InTf
    integer(C_INT), intent(IN), value :: seed                             !InTf
  end subroutine Hash32_init                                              !InTf
  function Hash32_add(hash_struct, input, length) result(status) BIND(C,name='Hash32_add') !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: hash_struct, input                  !InTf
    integer(C_INT), intent(IN), value :: length                           !InTf
    integer(C_INT) :: status                                              !InTf
  end function Hash32_add                                                 !InTf
  function Hash32_finalize(hash_struct) result(hash) BIND(C,name='Hash32_finalize')   !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: hash_struct                         !InTf
    integer(C_INT) :: hash                                                !InTf
  end function Hash32_finalize                                            !InTf
end interface                                                             !InTf
#endif
