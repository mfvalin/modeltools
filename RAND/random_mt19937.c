#include <stdint.h>
#include <randomgeneric.h>

#define MT_SIZE 624
#define MT_PERIOD 397
// MT_DIFF = 227
#define MT_DIFF (MT_SIZE-MT_PERIOD)

#define BIT31 0x80000000
#define B30_0 0x7FFFFFFF
#define MAGIC 0x9908b0df

static int32_t MT[MT_SIZE];
static int32_t MT2[MT_SIZE];
static int32_t mt_index = 0;

typedef struct{
  REFILLBUFFUN  refill;
  RANSETSEEDFUN seed;
  IRANFUN       iran;
  DRANFUN       dran;
  DRANSFUN      drans;
  IVECRANFUN    vec_iran;
  DVECRANFUN    vec_dran;
  DVECSRANFUN   vec_drans;
  unsigned int *gauss;
  int ngauss;
  int index;
  int bufsz;
  unsigned int *mt;
  unsigned int *mt2;
}mt19937_state ;                  // MT19937 generator stream control structure

static void FillBuffer_MT19937_stream(mt19937_state *MT19937);

static mt19937_state mt19937 = {
  (REFILLBUFFUN) FillBuffer_MT19937_stream, 
  (RANSETSEEDFUN) RanSetSeed_MT19937_stream, 
  (IRANFUN) IRan_MT19937_stream, 
  (DRANFUN) DRan_MT19937_stream, 
  (DRANSFUN) DRanS_MT19937_stream, 
  (IVECRANFUN) VecIRan_MT19937_stream, 
  (DVECRANFUN) VecDRan_MT19937_stream, 
  (DVECSRANFUN) VecDRanS_MT19937_stream, 
  NULL, 
  0,
  0 , 
  MT_SIZE,
  NULL,
  NULL };

generic_state *Ran_MT19937_new_stream(generic_state *clone_in, unsigned int *piSeed, int cSeed)   // !InTc! // create and seed a new stream
{
  mt19937_state *source ;
  mt19937_state *clone = (mt19937_state *)clone_in;
  int i;

  mt19937_state *new_state = (mt19937_state *) memalign(64,sizeof(mt19937_state)) ;

  if(cSeed < 0 && piSeed==NULL){  // clone a stream (mostly used for testing)
    source = clone ? clone : &mt19937;         // clone == NULL means clone default static state
    if(source->mt == NULL){                    // static copy not allocated/initialized
      source->mt = malloc(2 * sizeof(uint32_t) * MT_SIZE);
      source->mt2 = source->mt + MT_SIZE;
      RanSetSeed_MT19937_stream((generic_state *) source, piSeed, cSeed);
       FillBuffer_MT19937_stream(source);
    }
    new_state->index = source->index ;
    new_state->bufsz = source->bufsz ;
    new_state->ngauss = source->ngauss ;
    new_state->gauss = source->gauss ;
    new_state->seed = source->seed ;
    new_state->refill = source->refill ;
    new_state->iran = source->iran ;
    new_state->dran = source->dran ;
    new_state->drans = source->drans ;
    new_state->vec_iran = source->vec_iran ;
    new_state->vec_dran = source->vec_dran ;
    new_state->vec_drans = source->vec_drans ;
    new_state->mt = malloc(2 * sizeof(uint32_t) * MT_SIZE);
    new_state->mt2 = new_state->mt + MT_SIZE;
    for (i=0 ; i<MT_SIZE ; i++) new_state->mt[i] = source->mt[i] ;
    for (i=0 ; i<MT_SIZE ; i++) new_state->mt2[i] = source->mt2[i] ;
  }else{
    new_state->index = 0;
    new_state->bufsz = MT_SIZE;
    new_state->ngauss = 0;
    new_state->gauss = NULL;
    new_state->seed  = (RANSETSEEDFUN) RanSetSeed_MT19937_stream;
    new_state->refill = (REFILLBUFFUN) FillBuffer_MT19937_stream;
    new_state->iran   = (IRANFUN) IRan_MT19937_stream;
    new_state->dran   = (DRANFUN) DRan_MT19937_stream;
    new_state->drans  = (DRANSFUN) DRanS_MT19937_stream;
    new_state->vec_iran  = (IVECRANFUN) VecIRan_MT19937_stream;
    new_state->vec_dran  = (DVECRANFUN) VecDRan_MT19937_stream;
    new_state->vec_drans = (DVECSRANFUN) VecDRanS_MT19937_stream;
    new_state->mt = malloc(2 * sizeof(uint32_t) * MT_SIZE);
    new_state->mt2 = new_state->mt + MT_SIZE;
    RanSetSeed_MT19937_stream((generic_state *) new_state, piSeed, cSeed);  // seed the new stream
    FillBuffer_MT19937_stream(new_state);                                   // fill state buffer
  }

  return ( (void *) new_state) ;
}

void RanSetSeed_MT19937_stream(generic_state *stream, unsigned int *piSeed, int cSeed)    // !InTc!  initial seed
{
  mt19937_state *mt19937 = (mt19937_state *) stream;
  int i;
  uint32_t *mt = mt19937->mt ;
  /*
   * The equation below is a linear congruential generator (LCG),
   * one of the oldest known pseudo-random number generator
   * algorithms, in the form X_(n+1) = = (a*X_n + c) (mod m).
   *
   * We've implicitly got m=32 (mask + word size of 32 bits), so
   * there is no need to explicitly use modulus.
   *
   * What is interesting is the multiplier a.  The one we have
   * below is 0x6c07865 --- 1812433253 in decimal, and is called
   * the Borosh-Niederreiter multiplier for modulus 2^32.
   *
   * It is mentioned in passing in Knuth's THE ART OF COMPUTER
   * PROGRAMMING, Volume 2, page 106, Table 1, line 13.  LCGs are
   * treated in the same book, pp. 10-26
   *
   * You can read the original paper by Borosh and Niederreiter
   * as well.  It's called OPTIMAL MULTIPLIERS FOR PSEUDO-RANDOM
   * NUMBER GENERATION BY THE LINEAR CONGRUENTIAL METHOD (1983) at
   * http://www.springerlink.com/content/n7765ku70w8857l7/
   *
   * You can read about LCGs at:
   * http://en.wikipedia.org/wiki/Linear_congruential_generator
   *
   * From that page, it says:
   * "A common Mersenne twister implementation, interestingly
   * enough, uses an LCG to generate seed data.",
   *
   * Since we're using 32-bits data types for our MT array, we can skip the
   * masking with 0xFFFFFFFF below.
   */
  mt19937->index = 0;
  if((piSeed != NULL) && (cSeed > 0)){
    mt[0] = *piSeed;
  }else{
    mt[0] = 1;   // default seed if piSeed == NULL or cSeed <= 0
  }
  for(i=1; i<MT_SIZE; ++i ) {
    mt[i] = 0x6c078965*(mt[i-1] ^ mt[i-1]>>30) + i;
  }
}

static void FillBuffer_MT19937_stream(mt19937_state *stream)
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  uint32_t i, y, t, v;
  int32_t t2;
  int32_t *MT = MT19937->mt;
  int32_t *MT2 = MT19937->mt2;

  for(i=0 ; i<MT_DIFF ; i++){            // 0 - 396
    y = (MT[i] & BIT31) | (MT[i+1] & B30_0);
    v = MT[i+MT_PERIOD] ^ (y >> 1);
//     t = MATRIX[odd(y)]
    t2 = y ;
    t2 = (t2 << 31) ;   // lower bit into signt bit position
    t2 = t2 >> 31 ;     // propagate sign bit to all bits
    t = t2 ;
    t = t & MAGIC ;     // 0 if even, 0x9908b0df if odd
    MT[i] = v ^ t;
  }

  for(i=MT_DIFF ; i<MT_SIZE-1 ; i++){     // 397 - 622
    y = (MT[i] & BIT31) | (MT[i+1] & B30_0);
    v = MT[i-MT_DIFF] ^ (y >> 1);
//     t = MATRIX[odd(y)]
    t2 = y ;
    t2 = (t2 << 31) ;
    t2 = t2 >> 31 ;
    t = t2 ;
    t = t & MAGIC ;
    MT[i] = v ^ t;
  }
  y = (MT[MT_SIZE-1] & BIT31) | (MT[0] & B30_0);  // 623
  v = MT[MT_SIZE-1-MT_DIFF] ^ (y >> 1);
//     t = MATRIX[odd(y)]
    t2 = y ;
    t2 = (t2 << 31) ;
    t2 = t2 >> 31 ;
    t = t2 ;
  t = t & MAGIC ;
  MT[MT_SIZE-1] = v ^ t;

  for(i=0 ; i<MT_SIZE ; i++){
    y = MT[i];
    y ^= y>>11;
    y ^= y<< 7 & 0x9d2c5680;
    y ^= y<<15 & 0xefc60000;
    y ^= y>>18;
    MT2[i] = y;
  }
  MT19937->index = 0;
}

unsigned int IRan_MT19937_stream(generic_state *stream)       // !InTc!  returns a random unsigned integer
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  register uint32_t index = MT19937->index;
  register uint32_t y ;

  if (index >= MT_SIZE) FillBuffer_MT19937_stream((mt19937_state *)MT19937);

  y = MT19937->mt2[MT19937->index];
  MT19937->index = MT19937->index + 1 ;
  return y;
}

double DRan_MT19937_stream(generic_state *stream)       // !InTc!  returns a random double ( 0.0 , 1.0 )
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  register uint32_t index = MT19937->index;
  register uint32_t y ;

  if (index >= MT_SIZE) FillBuffer_MT19937_stream((mt19937_state *)MT19937);

  y = MT19937->mt2[MT19937->index];
  MT19937->index = MT19937->index + 1 ;
  return CVTDBL_32(y);
}

double DRanS_MT19937_stream(generic_state *stream)     // !InTc!   returns a random double (-1.0 , 1.0)
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  register uint32_t index = MT19937->index;
  register uint32_t y ;

  if (index >= MT_SIZE) FillBuffer_MT19937_stream((mt19937_state *)MT19937);

  y = MT19937->mt2[MT19937->index];
  MT19937->index = MT19937->index + 1 ;
  return CVTDBLS_32(y);
}

void VecIRan_MT19937_stream(generic_state *stream, uint32_t *ranbuf, int n)  // !InTc!  returns a vector of unsigned integers
{
  mt19937_state *MT19937 = (mt19937_state *) stream ; //
  unsigned int *mt2 = MT19937->mt2 ;
  int index = MT19937->index ;
  int i, navail;

  while(n > 0) {
    if (index >= MT_SIZE) FillBuffer_MT19937_stream((mt19937_state *)MT19937);
    index = MT19937->index ;
    navail = MT_SIZE - index;
    if(navail >= n){
      for(i=0 ; i<n ; i++) ranbuf[i] = mt2[index+i] ;
      MT19937->index = index + n;
      return;
    }else{
      for(i=0 ; i<navail ; i++) ranbuf[i] = mt2[index+i] ;
      index = MT_SIZE;
      n -= navail;
    }
  }
}

void VecDRan_MT19937_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!   returns a vector of doubles ( 0.0 , 1.0 )
{
  mt19937_state *MT19937 = (mt19937_state *) stream ; //
  unsigned int *mt2 = MT19937->mt2 ;
  int index = MT19937->index ;
  int i, navail;

  while(n > 0) {
    if (index >= MT_SIZE) FillBuffer_MT19937_stream((mt19937_state *)MT19937);
    index = MT19937->index ;
    navail = MT_SIZE - index;
    if(navail >= n){
      for(i=0 ; i<n ; i++) ranbuf[i] = CVTDBL_32(mt2[index+i]) ;
      MT19937->index = index + n;
      return;
    }else{
      for(i=0 ; i<navail ; i++) ranbuf[i] = CVTDBL_32(mt2[index+i]) ;
      index = MT_SIZE;
      n -= navail;
    }
  }
}

void VecDRanS_MT19937_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!   returns a vector of doubles ( -1.0 , 1.0 )
{
  mt19937_state *MT19937 = (mt19937_state *) stream ; //
  unsigned int *mt2 = MT19937->mt2 ;
  int index = MT19937->index ;
  int i, navail;

  while(n > 0) {
    if (index >= MT_SIZE) FillBuffer_MT19937_stream((mt19937_state *)MT19937);
    index = MT19937->index ;
    navail = MT_SIZE - index;
    if(navail >= n){
      for(i=0 ; i<n ; i++) ranbuf[i] = CVTDBLS_32(mt2[index+i]) ;
      MT19937->index = index + n;
      return;
    }else{
      for(i=0 ; i<navail ; i++) ranbuf[i] = CVTDBLS_32(mt2[index+i]) ;
      index = MT_SIZE;
      n -= navail;
    }
  }
}
