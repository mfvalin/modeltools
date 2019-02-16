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
// SHR3 3-shift-register generator (state = 32 bit integer)   // !InTc!
#if defined(NEVER_TO_BE_TRUE)

! void F_Ran_SHR3_new_stream(shr3_state *clone, int *piSeed, int cSeed)                   !InTf!
 interface                                                                                !InTf!
   subroutine Ran_SHR3_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_SHR3_new_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(OUT) :: stream                                             !InTf!
   type(RANDOM_STREAM), intent(IN) :: clone                                               !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine Ran_SHR3_new_stream                                                     !InTf!
 end interface                                                                            !InTf!
#endif
// macro used to instrument code with counters
#define INSTRUMENT(A) ;

#include <randomgeneric.h>

/*==========================================================================
 * SHR3 generators, function names consistent with the other generators
 *==========================================================================*/
/*------------------------ start of SHR3 routines --------------------------*/

typedef struct{
  GENERIC_STATE
  uint32_t jsr;
} shr3_state;                  // SHR3 generator stream control structure
  
static void FillBuffer_SHR3_stream(generic_state *stream){ 
  shr3_state *state = (shr3_state *) stream ;
  int i;
  int topp1 = state->top + 1;
  uint32_t *buf = state->buf;
  uint32_t jsr, jz;

  jsr = state->jsr ;
  for(i=0 ; i<topp1 ; i++) {
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    buf[i] = (jz+jsr);
  };
  state->jsr=jsr;
  state->cur = 0;
}

void RanSetSeed_SHR3_stream(generic_state *stream, uint32_t *piSeed, int cSeed)  // !InTc!
{
  shr3_state *state = (shr3_state *) stream ; 
  if(piSeed == NULL || cSeed == 0) { // default reseeding call
    state->jsr = 123456;
  }else{
    state->jsr = (uint32_t) *piSeed ;
  }
  state->cur = state->top + 1;   // mark buffer as empty
}

uint32_t IRan_SHR3_stream(generic_state *stream)	  // !InTc!	/* returns a random unsigned integer */
{
  shr3_state *state = (shr3_state *) stream ;
  uint32_t jz;
  uint32_t jsr;

  jsr=state->jsr ;
  jz = jsr;
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  state->jsr=jsr;
  return (jz+jsr) ;
}

double DRan_SHR3_stream(generic_state *stream)	  // !InTc!	/* returns a random double (0.0 , 1.0) */
{
  shr3_state *state = (shr3_state *) stream ;
  uint32_t jz;
  uint32_t jsr;

  jsr=state->jsr ;
  jz = jsr;
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  state->jsr=jsr;
  return CVTDBL_32(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
}

double DRanS_SHR3_stream(generic_state *stream)	  // !InTc!	/* returns a random double (-1.0 , 1.0) */
{
  shr3_state *state = (shr3_state *) stream ;
  uint32_t jz;
  uint32_t jsr;

  jsr=state->jsr ;
  jz = jsr;
  jsr^=(jsr<<13) ;
  jsr^=(jsr>>17) ;
  jsr^=(jsr<<5) ;
  state->jsr=jsr;
  return CVTDBLS_32(jz+jsr);   // convert from 32 bit int to (-1.0 , 1.0)
}

void VecIRan_SHR3_stream(generic_state *stream, uint32_t *ranbuf, int n)  // !InTc!
{
  int i;
  shr3_state *state = (shr3_state *) stream ;
  uint32_t jz;
  uint32_t jsr;

  jsr = state->jsr ;
  for(i=0 ; i<n ; i++) {
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = (jz+jsr);
  };
  state->jsr=jsr;
}

void VecDRan_SHR3_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!
{
  int i;
  shr3_state *state = (shr3_state *) stream ;
  uint32_t jz;
  uint32_t jsr;

  jsr = state->jsr ;
  for(i=0 ; i<n ; i++){
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = CVTDBL_32(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
  }
  state->jsr=jsr;
}

void VecDRanS_SHR3_stream(generic_state *stream, double *ranbuf, int n)  // !InTc!
{
  int i;
  shr3_state *state = (shr3_state *) stream ;
  uint32_t jz;
  uint32_t jsr;

  jsr = state->jsr ;
  for(i=0 ; i<n ; i++){
    jz=jsr ;
    jsr^=(jsr<<13) ;
    jsr^=(jsr>>17) ;
    jsr^=(jsr<<5) ;
    ranbuf[i] = CVTDBLS_32(jz+jsr);   // convert from 32 bit int to (0.0 , 1.0)
  }
  state->jsr=jsr;
}

static shr3_state SHR3 = { 
  (REFILLBUFFUN) FillBuffer_SHR3_stream, 
  (RANSETSEEDFUN) RanSetSeed_SHR3_stream, 
  (IRANFUN) IRan_SHR3_stream, 
  (DRANFUN) DRan_SHR3_stream, 
  (DRANSFUN) DRanS_SHR3_stream, 
  (IVECRANFUN) VecIRan_SHR3_stream, 
  (DVECRANFUN) VecDRan_SHR3_stream, 
  (DVECSRANFUN) VecDRanS_SHR3_stream, 
  NULL, 
  -1,
  -1,
  NULL, 
  0,
  123456789 } ;

void *Ran_SHR3_new_stream(void *clone_in, uint32_t *piSeed, int cSeed)   // !InTc! // create and seed a new stream
{
  shr3_state *source ;
  shr3_state *clone = (shr3_state *)clone_in;
  int i;

  shr3_state *new_state = (shr3_state *) memalign(64,sizeof(shr3_state)) ;

  if(cSeed < 0 && piSeed==NULL){            // clone a stream (mostly used for testing)
    source = clone ? clone : &SHR3;         // clone == NULL means clone default stream
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
    new_state->jsr = source->jsr;
    new_state->buf = malloc(sizeof(uint32_t) * 256) ;
    new_state->cur = source->cur ;
    new_state->top = source->top ;
    for (i=0 ; i<256 ; i++) new_state->buf[i] = source->buf[i] ;
  }else{
    new_state->ngauss = 0;
    new_state->gauss = NULL;
    new_state->seed  = (RANSETSEEDFUN) RanSetSeed_SHR3_stream;
    new_state->refill = (REFILLBUFFUN) FillBuffer_SHR3_stream;
    new_state->iran   = (IRANFUN) IRan_SHR3_stream;
    new_state->dran   = (DRANFUN) DRan_SHR3_stream;
    new_state->drans  = (DRANSFUN) DRanS_SHR3_stream;
    new_state->vec_iran  = (IVECRANFUN) VecIRan_SHR3_stream;
    new_state->vec_dran  = (DVECRANFUN) VecDRan_SHR3_stream;
    new_state->vec_drans = (DVECSRANFUN) VecDRanS_SHR3_stream;
    new_state->buf = malloc(sizeof(uint32_t) * 256) ;
    new_state->cur = 256 ;
    new_state->top = 255 ;
    RanSetSeed_SHR3_stream((generic_state *) new_state, piSeed, cSeed);  // seed the new stream
  }

  return ( (void *) new_state) ;
}
void F_Ran_SHR3_new_stream(statep *s, statep *c, uint32_t *piSeed, int cSeed)  // Fortran interface using derived type
{
  s->p = Ran_SHR3_new_stream( c->p, piSeed, cSeed);
}

/*------------------------- end of SHR3 routines ---------------------------*/
