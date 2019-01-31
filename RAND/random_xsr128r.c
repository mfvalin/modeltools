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
// xorshiftrotate generator with 128 bit state   // !InTc!
#if defined(NEVER_TO_BE_TRUE)

! void F_Ran_XSR128R_new_stream(xsr128r_state *clone, int *piSeed, int cSeed)                   !InTf!
 interface                                                                                !InTf!
   subroutine Ran_XSR128R_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_XSR128R_new_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(OUT) :: stream                                             !InTf!
   type(RANDOM_STREAM), intent(IN) :: clone                                               !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine Ran_XSR128R_new_stream                                                     !InTf!
 end interface                                                                            !InTf!
#endif
// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

/*==========================================================================
 * XSR128R generators, function names consistent with the other generators
 *==========================================================================*/
/*------------------------ start of XSR128R routines --------------------------*/

typedef struct{
  GENERIC_STATE
  uint32_t part128r;
  uint64_t s[2];
  uint64_t res128r;
} xsr128r_state;                  // XSR128R generator stream control structure
  
static void FillBuffer_XSR128R_stream(generic_state *stream){ }

static void jump128r(uint64_t *s);

// 4 32 bit quantities are needed for seeding. 
// if cSeed > 4, only the first 4 values will be used
// if cSeed is negative, this is interpreted as a jump ahead in the sequence call
void RanSetSeed_XSR128R_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  xsr128r_state *state = (xsr128r_state *) stream ; 
  if(cSeed < 0) {            // jump call
    jump128r(state->s);
    state->part128r = 0 ;     // force "no residual" condition
    state->res128r = 0;
    return;
  }
  if(cSeed >= 4 && piSeed != NULL) {   // enough bits to seed
    state->s[0] = piSeed[0] ;  state->s[0] = (state->s[0] << 32) | piSeed[1] ;
    state->s[1] = piSeed[2] ;  state->s[1] = (state->s[1] << 32) | piSeed[3] ;
  }else{
    state->s[0] = 123456;
    state->s[1] = 456789;
  }
  state->part128r = 0 ;     // force "no residual" condition
  state->res128r = 0;
}

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

static inline uint64_t advance128r(uint64_t *s) {
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

static void jump128r(uint64_t *s) {
  static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  int i, b;
  for(i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
    for(b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
        s0 ^= s[0];
        s1 ^= s[1];
      }
      advance128r(s);
    }

  s[0] = s0;
  s[1] = s1;
}

unsigned int IRan_XSR128R_stream(generic_state *stream)	  // !InTc!	/* returns a random unsigned integer */
{
  xsr128r_state *state = (xsr128r_state *) stream ;
  uint32_t result;

  state->part128r = state->part128r ^ 1;         // complement previous residual flag
  if(state->part128r) {                  // there was no previous residual
    state->res128r = advance128r(state->s);  // strore residual
    result = state->res128r >> 32 ;    // return upper part
  }else{                         // there was a previous residual, return lower part
    result = state->res128r;
  }
  return result;
}

double DRan_XSR128R_stream(generic_state *stream)	  // !InTc!	/* returns a random double (0.0 , 1.0) */
{
  xsr128r_state *state = (xsr128r_state *) stream ;
  uint32_t result;

  state->part128r = state->part128r ^ 1;         // complement previous residual flag
  if(state->part128r) {                  // there was no previous residual
    state->res128r = advance128r(state->s);  // strore residual
    result = state->res128r >> 32 ;    // return upper part
  }else{                         // there was a previous residual, return lower part
    result = state->res128r;
  }
  return CVTDBL_32(result);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_XSR128R_stream(generic_state *stream)	  // !InTc!	/* returns a random double (-1.0 , 1.0) */
{
  xsr128r_state *state = (xsr128r_state *) stream ;
  uint32_t result;

  state->part128r = state->part128r ^ 1;         // complement previous residual flag
  if(state->part128r) {                  // there was no previous residual
    state->res128r = advance128r(state->s);  // strore residual
    result = state->res128r >> 32 ;    // return upper part
  }else{                         // there was a previous residual, return lower part
    result = state->res128r;
  }
  return CVTDBLS_32(result);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_XSR128R_stream(generic_state *stream, unsigned int *ranbuf, int n)  // !InTc!
{
  int i;
  uint64_t t;
  xsr128r_state *state = (xsr128r_state *) stream ;

  ranbuf[0] = state->res128r;               // in case there was a previous residual, store lower part in dest
  for(i = state->part128r ; i < (n-1) ; i+=2){
    t = advance128r(state->s);
    ranbuf[i] = (t >> 32);
    ranbuf[i+1] = t;
  }
  state->part128r = (i < n);
  if(state->part128r){                   // there will be a residual, save it, store upper part in dest
    t = advance128r(state->s);
    state->res128r = t;
    ranbuf[i] = (t >> 32);
  }
}

void VecDRan_XSR128R_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!
{
  int i;
  uint64_t t;
  xsr128r_state *state = (xsr128r_state *) stream ;
  uint32_t tt;

  tt = state->res128r;
  ranbuf[0] = CVTDBL_32(tt);               // in case there was a previous residual, store lower part in dest
  for(i = state->part128r ; i < (n-1) ; i+=2){
    t = advance128r(state->s);
    tt = t >> 32;
    ranbuf[i] = CVTDBL_32(tt);
    tt = t;
    ranbuf[i+1] = CVTDBL_32(tt);
  }
  state->part128r = (i < n);
  if(state->part128r){                   // there will be a residual, save it, store upper part in dest
    t = advance128r(state->s);
    state->res128r = t;
    tt = t >> 32;
    ranbuf[i] = CVTDBL_32(tt);
  }
}

void VecDRanS_XSR128R_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!
{
  int i;
  uint64_t t;
  xsr128r_state *state = (xsr128r_state *) stream ;
  uint32_t tt;

  tt = state->res128r;
  ranbuf[0] = CVTDBLS_32(tt);               // in case there was a previous residual, store lower part in dest
  for(i = state->part128r ; i < (n-1) ; i+=2){
    t = advance128r(state->s);
    tt = t >> 32;
    ranbuf[i] = CVTDBLS_32(tt);
    tt = t;
    ranbuf[i+1] = CVTDBLS_32(tt);
  }
  state->part128r = (i < n);
  if(state->part128r){                   // there will be a residual, save it, store upper part in dest
    t = advance128r(state->s);
    state->res128r = t;
    tt = t >> 32;
    ranbuf[i] = CVTDBLS_32(tt);
  }
}

static xsr128r_state XSR128R = { 
  (REFILLBUFFUN) FillBuffer_XSR128R_stream, 
  (RANSETSEEDFUN) RanSetSeed_XSR128R_stream, 
  (IRANFUN) IRan_XSR128R_stream, 
  (DRANFUN) DRan_XSR128R_stream, 
  (DRANSFUN) DRanS_XSR128R_stream, 
  (IVECRANFUN) VecIRan_XSR128R_stream, 
  (DVECRANFUN) VecDRan_XSR128R_stream, 
  (DVECSRANFUN) VecDRanS_XSR128R_stream, 
  NULL, 
  -1,
  -1,
  NULL, 
  0,
  0,
  { 123456, 456789 },
  0 } ;

void *Ran_XSR128R_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   // !InTc! // create and seed a new stream
{
  xsr128r_state *source ;
  xsr128r_state *clone = (xsr128r_state *)clone_in;
  xsr128r_state *new_state = (xsr128r_state *) memalign(64,sizeof(xsr128r_state)) ;

  if(cSeed < 0 && piSeed==NULL){  // clone a stream (mostly used for testing)
    source = clone ? clone : &XSR128R;         // clone == NULL means clone default stream
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
    new_state->res128r = source->res128r;
    new_state->part128r = source->part128r;
  }else{
    new_state->ngauss = 0;
    new_state->gauss = NULL;
    new_state->seed  = (RANSETSEEDFUN) RanSetSeed_XSR128R_stream;
    new_state->refill = (REFILLBUFFUN) FillBuffer_XSR128R_stream;
    new_state->iran   = (IRANFUN) IRan_XSR128R_stream;
    new_state->dran   = (DRANFUN) DRan_XSR128R_stream;
    new_state->drans  = (DRANSFUN) DRanS_XSR128R_stream;
    new_state->vec_iran  = (IVECRANFUN) VecIRan_XSR128R_stream;
    new_state->vec_dran  = (DVECRANFUN) VecDRan_XSR128R_stream;
    new_state->vec_drans = (DVECSRANFUN) VecDRanS_XSR128R_stream;
    RanSetSeed_XSR128R_stream((generic_state *) new_state, piSeed, cSeed);  // seed the new stream
    new_state->res128r = 0; // no residual
    new_state->part128r = 0;
  }

  return ( (void *) new_state) ;
}
void F_Ran_XSR128R_new_stream(statep *s, statep *c, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  s->p = Ran_XSR128R_new_stream( c->p, piSeed, cSeed);
}

/*------------------------- end of XSR128R routines ---------------------------*/

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
  generic_state *XSR128R;
  unsigned int piSeed[4] = {123456, 234567, 345678, 456789 };
  double dcount;

  MPI_Init(&argc,&argv);

  XSR128R = Ran_XSR128R_new_stream(NULL, piSeed, 4);

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

  for( i=0 ; i < 1000000 ; i++) lr = IRan_generic_stream(XSR128R);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) lr = IRan_generic_stream(XSR128R);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random XSR128R integer value = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_generic_stream(XSR128R);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) rval = DRan_generic_stream(XSR128R);
  t1 = MPI_Wtime();
  printf("time for 1E+9 x 1 random XSR128R double value = %6.3f \n",t1-t0);

  for( i=0 ; i < 10 ; i++) VecIRan_generic_stream(XSR128R, ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_generic_stream(XSR128R, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if((int)ranbuf[i] > 0) pos++ ; else neg++ ;
//   pos = 0 ; neg = 0 ; for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & 1) pos++ ; else neg++ ;
  postot = 0 ; negtot = 0;
  for (j=0 ; j<100 ; j++) {
    VecIRan_generic_stream(XSR128R, ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  printf("%5d ",pos-neg) ;
    }
  }
//   printf("%d\n",postot-negtot);
  dcount = 32.0 * 100 * 1000000 ; dcount = sqrt(dcount) ;
  printf("time for 1E+3 x 1E+6 random XSR128R integer values = %6.3f , pos - neg = %d (%9.0f)\n",t1-t0,postot-negtot,dcount);

  for( i=0 ; i < 10 ; i++) VecDRan_generic_stream(XSR128R, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_generic_stream(XSR128R, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+3 x 1E+6 random XSR128R double values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecIRan_generic_stream(XSR128R, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_generic_stream(XSR128R, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random XSR128R integer values = %6.3f \n",t1-t0);

  for( i=0 ; i < 1000 ; i++) VecDRan_generic_stream(XSR128R, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_generic_stream(XSR128R, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  printf("time for 1E+6 x 1E+3 random XSR128R double values = %6.3f \n",t1-t0);

  MPI_Finalize();
  return(0);
}
#endif
