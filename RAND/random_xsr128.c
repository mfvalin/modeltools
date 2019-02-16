
/*  
heavvily modified code from original code written by Sebastiano Vigna (vigna@acm.org)
that had been put in the public domain

"To the extent possible under law, the author has dedicated all copyright
 and related and neighboring rights to this software to the public domain
 worldwide. This software is distributed without any warranty.            "

See <http://creativecommons.org/publicdomain/zero/1.0/>.
*/

// xorshift generator with 128 bit state   // !InTc!
#if defined(NEVER_TO_BE_TRUE)

! void F_Ran_XSR128_new_stream(xsr128_state *clone, int *piSeed, int cSeed)               !InTf!
 interface                                                                                !InTf!
   subroutine Ran_XSR128_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_XSR128_new_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(OUT) :: stream                                             !InTf!
   type(RANDOM_STREAM), intent(IN) :: clone                                               !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine Ran_XSR128_new_stream                                                   !InTf!
 end interface                                                                            !InTf!
#endif
// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

/*==========================================================================
 * XSR128 generators, function names consistent with the other generators
 *==========================================================================*/
/*------------------------ start of XSR128 routines --------------------------*/

typedef struct{
  GENERIC_STATE
  uint32_t part128;
  uint64_t s[2];
  uint64_t res128;
} xsr128_state;                  // XSR128 generator stream control structure

static void jump128(uint64_t *s);

// 4 32 bit quantities are needed for seeding. 
// if cSeed > 4, only the first 4 values will be used
// if cSeed is negative, this is interpreted as a jump ahead in the sequence call
void RanSetSeed_XSR128_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  xsr128_state *state = (xsr128_state *) stream ; 
  if(cSeed < 0) {            // jump call
    jump128(state->s);
    state->part128 = 0 ;     // force "no residual" condition
    state->res128 = 0;
    state->cur = state->top+1 ;
    return;
  }
  if(cSeed >= 4 && piSeed != NULL) {   // enough bits to seed
    state->s[0] = piSeed[0] ;  state->s[0] = (state->s[0] << 32) | piSeed[1] ;
    state->s[1] = piSeed[2] ;  state->s[1] = (state->s[1] << 32) | piSeed[3] ;
  }else{
    state->s[0] = 123456;
    state->s[1] = 456789;
  }
  state->part128 = 0 ;       // force "no residual" condition
  state->res128 = 0;
  state->cur = state->top+1 ;
}

static inline uint64_t advance128(uint64_t *s) {
  uint64_t s1 = s[0];
  const uint64_t s0 = s[1];
  const uint64_t result = s0 + s1;
  s[0] = s0;
  s1 ^= s1 << 23; // a
  s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
  return result; 
}
  
static void FillBuffer_XSR128_stream(generic_state *stream){ 
  xsr128_state *state = (xsr128_state *) stream ; 
  int topp1 = state->top + 1;
  uint32_t *buf = state->buf;
  int i;
  uint64_t t;

  for(i = 0 ; i < topp1 ; i+=2){
    t = advance128(state->s);
    buf[i] = (t >> 32);
    buf[i+1] = t;
  }
  state->cur = 0;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

static void jump128(uint64_t *s) {
  static const uint64_t JUMP[] = { 0x8a5cd789635d2dff, 0x121fd2155c472f96 };

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  int i, b;
  for(i = 0; i < sizeof JUMP / sizeof *JUMP; i++){
    for(b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
        s0 ^= s[0];
        s1 ^= s[1];
      }
      advance128(s);
    }
  }

  s[0] = s0;
  s[1] = s1;
}

unsigned int IRan_XSR128_stream(generic_state *stream)	  // !InTc!	/* returns a random unsigned integer */
{
  xsr128_state *state = (xsr128_state *) stream ;
  uint32_t result;

  state->part128 = state->part128 ^ 1;         // complement previous residual flag
  if(state->part128) {                  // there was no previous residual
    state->res128 = advance128(state->s);  // strore residual
    result = state->res128 >> 32 ;    // return upper part
  }else{                         // there was a previous residual, return lower part
    result = state->res128;
  }
  return result;
}

double DRan_XSR128_stream(generic_state *stream)	  // !InTc!	/* returns a random double (0.0 , 1.0) */
{
  xsr128_state *state = (xsr128_state *) stream ;
  uint32_t result;

  state->part128 = state->part128 ^ 1;         // complement previous residual flag
  if(state->part128) {                  // there was no previous residual
    state->res128 = advance128(state->s);  // strore residual
    result = state->res128 >> 32 ;    // return upper part
  }else{                         // there was a previous residual, return lower part
    result = state->res128;
  }
  return CVTDBL_32(result);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_XSR128_stream(generic_state *stream)	  // !InTc!	/* returns a random double (-1.0 , 1.0) */
{
  xsr128_state *state = (xsr128_state *) stream ;
  uint32_t result;

  state->part128 = state->part128 ^ 1;         // complement previous residual flag
  if(state->part128) {                  // there was no previous residual
    state->res128 = advance128(state->s);  // strore residual
    result = state->res128 >> 32 ;    // return upper part
  }else{                         // there was a previous residual, return lower part
    result = state->res128;
  }
  return CVTDBLS_32(result);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_XSR128_stream(generic_state *stream, unsigned int *ranbuf, int n)  // !InTc!
{
  int i;
  uint64_t t;
  xsr128_state *state = (xsr128_state *) stream ;

  ranbuf[0] = state->res128;               // in case there was a previous residual, store lower part in dest
  for(i = state->part128 ; i < (n-1) ; i+=2){
    t = advance128(state->s);
    ranbuf[i] = (t >> 32);
    ranbuf[i+1] = t;
  }
  state->part128 = (i < n);
  if(state->part128){                   // there will be a residual, save it, store upper part in dest
    t = advance128(state->s);
    state->res128 = t;
    ranbuf[i] = (t >> 32);
  }
}

void VecDRan_XSR128_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!
{
  int i;
  uint64_t t;
  xsr128_state *state = (xsr128_state *) stream ;
  uint32_t tt;

  tt = state->res128;
  ranbuf[0] = CVTDBL_32(tt);               // in case there was a previous residual, store lower part in dest
  for(i = state->part128 ; i < (n-1) ; i+=2){
    t = advance128(state->s);
    tt = t >> 32;
    ranbuf[i] = CVTDBL_32(tt);
    tt = t;
    ranbuf[i+1] = CVTDBL_32(tt);
  }
  state->part128 = (i < n);
  if(state->part128){                   // there will be a residual, save it, store upper part in dest
    t = advance128(state->s);
    state->res128 = t;
    tt = t >> 32;
    ranbuf[i] = CVTDBL_32(tt);
  }
}

void VecDRanS_XSR128_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!
{
  int i;
  uint64_t t;
  xsr128_state *state = (xsr128_state *) stream ;
  uint32_t tt;

  tt = state->res128;
  ranbuf[0] = CVTDBLS_32(tt);               // in case there was a previous residual, store lower part in dest
  for(i = state->part128 ; i < (n-1) ; i+=2){
    t = advance128(state->s);
    tt = t >> 32;
    ranbuf[i] = CVTDBLS_32(tt);
    tt = t;
    ranbuf[i+1] = CVTDBLS_32(tt);
  }
  state->part128 = (i < n);
  if(state->part128){                   // there will be a residual, save it, store upper part in dest
    t = advance128(state->s);
    state->res128 = t;
    tt = t >> 32;
    ranbuf[i] = CVTDBLS_32(tt);
  }
}

static xsr128_state XSR128 = { 
  (REFILLBUFFUN) FillBuffer_XSR128_stream, 
  (RANSETSEEDFUN) RanSetSeed_XSR128_stream, 
  (IRANFUN) IRan_XSR128_stream, 
  (DRANFUN) DRan_XSR128_stream, 
  (DRANSFUN) DRanS_XSR128_stream, 
  (IVECRANFUN) VecIRan_XSR128_stream, 
  (DVECRANFUN) VecDRan_XSR128_stream, 
  (DVECSRANFUN) VecDRanS_XSR128_stream, 
  NULL, 
  -1,
  -1,
  NULL, 
  0,
  0,
  {123456, 456789 },
  0 } ;

void *Ran_XSR128_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   // !InTc! // create and seed a new stream
{
  xsr128_state *source ;
  xsr128_state *clone = (xsr128_state *)clone_in;
  xsr128_state *new_state = (xsr128_state *) memalign(64,sizeof(xsr128_state)) ;
  int i;

  if(cSeed < 0 && piSeed==NULL){  // clone a stream (mostly used for testing)
    source = clone ? clone : &XSR128;         // clone == NULL means clone default stream
    new_state->ngauss = source->ngauss ;
    new_state->gauss = source->gauss ;
    new_state->seed  = source->seed;
    new_state->refill = source->refill ;
    new_state->iran = source->iran ;
    new_state->dran = source->dran ;
    new_state->drans = source->drans ;
    new_state->vec_iran = source->vec_iran ;
    new_state->vec_dran = source->vec_dran ;
    new_state->vec_drans = source->vec_drans ;
    new_state->s[0] = source->s[0];
    new_state->s[1] = source->s[1];
    new_state->res128 = source->res128;
    new_state->part128 = source->part128;
    new_state->buf = malloc(sizeof(uint32_t) * 256) ;
    new_state->cur = source->cur ;
    new_state->top = source->top ;
    for (i=0 ; i<256 ; i++) new_state->buf[i] = source->buf[i] ;
  }else{
    new_state->ngauss = 0;
    new_state->gauss = NULL;
    new_state->seed  = (RANSETSEEDFUN) RanSetSeed_XSR128_stream;
    new_state->refill = (REFILLBUFFUN) FillBuffer_XSR128_stream;
    new_state->iran   = (IRANFUN) IRan_XSR128_stream;
    new_state->dran   = (DRANFUN) DRan_XSR128_stream;
    new_state->drans  = (DRANSFUN) DRanS_XSR128_stream;
    new_state->vec_iran  = (IVECRANFUN) VecIRan_XSR128_stream;
    new_state->vec_dran  = (DVECRANFUN) VecDRan_XSR128_stream;
    new_state->vec_drans = (DVECSRANFUN) VecDRanS_XSR128_stream;
    RanSetSeed_XSR128_stream((generic_state *) new_state, piSeed, cSeed);  // seed the new stream
    new_state->res128 = 0; // no residual
    new_state->part128 = 0;
    new_state->buf = malloc(sizeof(uint32_t) * 256) ;
    new_state->cur = 256 ;
    new_state->top = 255 ;
  }

  return ( (void *) new_state) ;
}
void F_Ran_XSR128_new_stream(statep *s, statep *c, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  s->p = Ran_XSR128_new_stream( c->p, piSeed, cSeed);
}

/*------------------------- end of XSR128 routines ---------------------------*/
