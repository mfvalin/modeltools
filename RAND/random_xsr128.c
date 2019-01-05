/* 
 * Copyright (C) 2019 Recherche en Prevision Numerique
 * 
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
// xorshift generator with 128 bit state   // !InTc!
#if defined(NEVER_TO_BE_TRUE)

! void F_Ran_XSR128_new_stream(xsr128_state *clone, int *piSeed, int cSeed)                   !InTf!
 interface                                                                                !InTf!
   subroutine Ran_XSR128_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_XSR128_new_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(OUT) :: stream                                             !InTf!
   type(RANDOM_STREAM), intent(IN) :: clone                                               !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine Ran_XSR128_new_stream                                                     !InTf!
 end interface                                                                            !InTf!
#endif
// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

/*==========================================================================
 * XSR128 generators, naming consistent with the MWC8222 code
 *==========================================================================*/
/*------------------------ start of XSR128 routines --------------------------*/

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
  int32_t ngauss;
  uint32_t part128;
  uint64_t s[2];
  uint64_t res128;
} xsr128_state;                  // XSR128 generator stream control structure
  
static void FillBuffer_XSR128_stream(generic_state *stream){ }

static void jump128(uint64_t *s);
void RanSetSeed_XSR128_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  xsr128_state *state = (xsr128_state *) stream ; 
  if(cSeed == -64) {
    jump128(state->s);
    return;
  }
  if(piSeed == NULL || cSeed == 0) return ; // null call, nothing to do
  if(cSeed == 4) { 
    state->s[0] = piSeed[0] ;  state->s[0] = (state->s[0] << 32) | piSeed[1] ;
    state->s[1] = piSeed[2] ;  state->s[1] = (state->s[1] << 32) | piSeed[3] ;
  }
  state->s[0] = (unsigned int) *piSeed ;
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
  0,
  0,
  {123456, 456789 },
  0 } ;

// #define XSR128 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)

void *Ran_XSR128_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   // !InTc! // create and seed a new stream
{
  xsr128_state *source ;
  xsr128_state *clone = (xsr128_state *)clone_in;
  xsr128_state *new_state = (xsr128_state *) memalign(64,sizeof(xsr128_state)) ;

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
  }

  return ( (void *) new_state) ;
}
void F_Ran_XSR128_new_stream(statep *s, statep *c, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  s->p = Ran_XSR128_new_stream( c->p, piSeed, cSeed);
}

/*------------------------- end of XSR128 routines ---------------------------*/

#if defined(SELF_TEST)
#include <mpi.h>
int main(int argc, char **argv){
  unsigned int lr;
  int i, j;
  double t0, t1, rval;
  double MPI_Wtime() ;
  unsigned int ranbuf[1200000];
  double ranbuf2[1200000];
  int pos, neg, mask, postot, negtot ;
  double dmax, dmin, avg;
  unsigned long long *idmax, *idmin ;
  unsigned int maxpos, maxneg;
  int gaussdist[10];
  int index;
  generic_state *xsr128;
  unsigned int piSeed = 123456789;

  MPI_Init(&argc,&argv);

  xsr128 = Ran_XSR128_new_stream(NULL, &piSeed, 1);

  for(i=0 ; i<1200000 ; i++) ranbuf[i] = 0;
  for(i=0 ; i<1200000 ; i++) ranbuf2[i] = 0.0;
  maxpos = 0x7FFFFFFF ;
  maxneg = 0x80000000 ;
  idmax = (unsigned long long *)&dmax;
  idmin = (unsigned long long *)&dmin;
  dmax = CVTDBL_32(maxpos) ;
  dmin = CVTDBL_32(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBL_32  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = CVTDBLS_32(maxpos) ;
  dmin = CVTDBLS_32(maxneg) ;
  printf("maxpos, maxneg transformed with CVTDBLS_32 : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

  for( i=0 ; i < 1000000 ; i++) lr = IRan_generic_stream(xsr128);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_generic_stream(xsr128);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random XSR128 integer value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_generic_stream(xsr128);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_generic_stream(xsr128);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random XSR128 double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecIRan_generic_stream(xsr128, ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_generic_stream(xsr128, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_generic_stream(xsr128, ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  printf("%5d ",pos-neg) ;
    }
  }
//   printf("%d\n",postot-negtot);
  printf("time for 1E+3 x 1E+6 random XSR128 integer values = %6.3f , pos - neg = %d\n",t1-t0,postot-negtot);

  for( i=0 ; i < 10 ; i++) VecDRan_generic_stream(xsr128, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_generic_stream(xsr128, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random XSR128 double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_generic_stream(xsr128, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_generic_stream(xsr128, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random XSR128 integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_generic_stream(xsr128, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_generic_stream(xsr128, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random XSR128 double values = %6.3f \n",t1-t0);

  MPI_Finalize();
  return(0);
}
#endif
