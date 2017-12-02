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
  unsigned int *state;
  unsigned int *state2;
}mt19937_state ;                  // MT19937 generator stream control structure

void *Ran_MT19937_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   ;
void RanSetSeed_MT19937_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  ;
void VecIRan_MT19937_stream(generic_state *stream, unsigned int *ranbuf, int n)  ;
void VecDRan_MT19937_stream(generic_state *stream, double *ranbuf, int n)  ;
void VecDRanS_MT19937_stream(generic_state *stream, double *ranbuf, int n)  ;

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

static void FillBuffer_MT19937_stream(mt19937_state *stream)
{
  uint32_t i, y, t, v;
  int32_t t2;
  int32_t *MT, *MT2;

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
}

unsigned int IRan_MT19937_stream(generic_state *stream)       // !InTc!  returns a random unsigned integer
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  register uint32_t index = MT19937->index;
  register uint32_t y ;

  if (index == 0) FillBuffer_MT19937_stream((mt19937_state *)MT19937);
  y = MT19937->state2[MT19937->index];
  if ( ++index == MT_SIZE ) index = 0;
  MT19937->index = index ;
  return y;
}

double DRan_MT19937_stream(generic_state *stream)       // !InTc!  returns a random double ( 0.0 , 1.0 )
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  register uint32_t index = MT19937->index;
  register uint32_t y ;

  if (index == 0) FillBuffer_MT19937_stream((mt19937_state *)MT19937);
  y = MT19937->state2[MT19937->index];
  if ( ++index == MT_SIZE ) index = 0;
  MT19937->index = index ;
  return CVTDBL_32(y);
}

double DRanS_MT19937_stream(generic_state *stream)     // !InTc!   returns a random double (-1.0 , 1.0)
{
  mt19937_state *MT19937 = (mt19937_state *) stream;
  register uint32_t index = MT19937->index;
  register uint32_t y ;

  if (index == 0) FillBuffer_MT19937_stream((mt19937_state *)MT19937);
  y = MT19937->state2[MT19937->index];
  if ( ++index == MT_SIZE ) index = 0;
  MT19937->index = index ;
  return CVTDBLS_32(y);
}
