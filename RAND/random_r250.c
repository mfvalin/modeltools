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

// R250 random number generator, from S.Kirkpatrick and E.  Stoll (state = 250 32 bit integers)  // !InTc!
// extremely fast generator (can use SIMD)

// Fortran interfaces for automated extraction
#if defined(NEVER_TO_BE_TRUE)

! void F_RanSetSeed_R250_stream(statep *s   , int *piSeed, int cSeed)                     !InTf!
 interface                                                                                !InTf!
   subroutine RanSetSeed_R250_stream(stream, piSeed, cSeed) bind(C,name='F_RanSetSeed_R250_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanSetSeed_R250_stream                                                  !InTf!
 end interface                                                                            !InTf!

! void F_Ran_R250_new_stream(r250_state *clone, int *piSeed, int cSeed)                   !InTf!
 interface                                                                                !InTf!
   subroutine Ran_R250_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_R250_new_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(OUT) :: stream                                             !InTf!
   type(RANDOM_STREAM), intent(IN) :: clone                                               !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine Ran_R250_new_stream                                                     !InTf!
 end interface                                                                            !InTf!

! unsigned int F_IRan_R250_stream(statep *s   )                                           !InTf!
 interface                                                                                !InTf!
   function IRan_R250_stream(stream) result(ran) bind(C,name='F_IRan_R250_stream')        !InTf!
   import :: C_INT,RANDOM_STREAM                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT) :: ran                                                                  !InTf!
   end function IRan_R250_stream                                                          !InTf!
 end interface                                                                            !InTf!

! double F_DRan_R250_stream(statep *s   )                                                 !InTf!
 interface                                                                                !InTf!
   function DRan_R250_stream(stream) result(ran) bind(C,name='F_DRan_R250_stream')        !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRan_R250_stream                                                          !InTf!
 end interface                                                                            !InTf!

! double F_DRanS_R250_stream(statep *s   )                                                !InTf!
 interface                                                                                !InTf!
   function DRanS_R250_stream(stream) result(ran) bind(C,name='F_DRanS_R250_stream')      !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRanS_R250_stream                                                         !InTf!
 end interface                                                                            !InTf!

! void F_VecIRan_R250_stream(statep *s   , unsigned int *ranbuf, int n)                   !InTf!
 interface                                                                                !InTf!
   subroutine VecIRan_R250_stream(stream, ranbuf, n) bind(C,name='F_VecIRan_R250_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecIRan_R250_stream                                                     !InTf!
 end interface                                                                            !InTf!

! void F_VecDRanS_R250_stream(statep *s   , double *ranbuf, int n)                        !InTf!
 interface                                                                                !InTf!
   subroutine VecDRanS_R250_stream(stream, ranbuf, n) bind(C,name='F_VecDRanS_R250_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecDRanS_R250_stream                                                    !InTf!
 end interface                                                                            !InTf!

! void F_VecDRan_R250_stream(statep *s   , double *ranbuf, int n)                         !InTf!
 interface                                                                                !InTf!
   subroutine VecDRan_R250_stream(stream, ranbuf, n) bind(C,name='F_VecDRan_R250_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecDRan_R250_stream                                                     !InTf!
 end interface                                                                            !InTf!

#endif

// get generic definitions
#include <randomgeneric.h>

/*==========================================================================
 * R250 generators, function names consistent with the other generators
 *==========================================================================*/

// the first items of the following struct MUST match ALL elements of the generic_state struct
// as defined in file randomgeneric.h
// the extra elements are specific to THIS generator type
typedef struct{
  GENERIC_STATE
//   int index;                   // index into internal state buffer
//   uint32_t *scrap;            // internal state (250 32 bit integers + 1 pad)
}r250_state ;                  // R250 generator stream control structure

// void VecIRan_R250_stream(void *stream, unsigned int *ranbuf, int n);
// static void FillBuffer_R250_stream(r250_state *R250);

// ======================= functions using stream provided state =============================

static void FillBuffer_R250_stream(r250_state *R250){
  int i;
  unsigned int *r250_buffer = R250->buf ;

  R250->cur = 0;
  R250->top = 249;

  for (i=0 ; i< 144 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i + 103 ];
  }
  r250_buffer[144] = r250_buffer[144] ^ r250_buffer[247];
  r250_buffer[145] = r250_buffer[145] ^ r250_buffer[248];
  r250_buffer[146] = r250_buffer[146] ^ r250_buffer[249];
  for (i=147 ; i<251 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i - 147 ];
  }
}

void RanSetSeed_R250_stream(void *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  int i;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer
  unsigned int *r250_buffer ;
  unsigned int seed ;

  r250_buffer = R250->buf ;
  if (cSeed == 250 && piSeed != NULL) {
    for(i=0 ; i<250; i++) r250_buffer[ i ] = piSeed[ i ] ;
  }else{
    seed = ( piSeed && (cSeed > 0) ) ? piSeed[0] : 0 ;
    for (i = 0 ; i < 250 ; ) {   /* see Knuth p.106, Table 1(16) and Numerical Recipes p.284 (ranqd1)*/
      seed = 1664525 * seed + 1013904223 ;
      if (seed <= 0) continue ;
      r250_buffer[i++] = seed ;
    }
  }
  R250->cur = 0;
  R250->top = 249;
}

void F_RanSetSeed_R250_stream(statep *s, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  RanSetSeed_R250_stream( s->p, piSeed, cSeed);
}

unsigned int IRan_R250_stream(void *stream)       // !InTc!     /* returns a single random unsigned integer */
{
  register unsigned int new_rand;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer

  if ( R250->cur > 249 ) FillBuffer_R250_stream(R250);
  new_rand = R250->buf[R250->cur++];
  return new_rand;
}
unsigned int F_IRan_R250_stream(statep *s)  // Fortran interface using derived type
{
  return(IRan_R250_stream(s->p));
}

double DRan_R250_stream(void *stream)	  // !InTc!	/* returns a random double (0.0 , 1.0) */
{
  register unsigned int new_rand;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer

  if ( R250->cur > 249 ) FillBuffer_R250_stream(R250);
  new_rand = R250->buf[R250->cur++] ;
  return CVTDBL_32(new_rand);   // convert from 32 bit int to (0.0 , 1.0)
}
unsigned int F_DRan_R250_stream(statep *s)  // Fortran interface using derived type
{
  return(DRan_R250_stream(s->p));
}

double DRanS_R250_stream(void *stream)	  // !InTc!	/* returns a random double (-1.0 , 1.0) */
{
  register unsigned int new_rand;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer

  if ( R250->cur > 249 ) FillBuffer_R250_stream(R250);
  new_rand = R250->buf[R250->cur++] ;
  return CVTDBLS_32(new_rand);   // convert from 32 bit int to (0.0 , 1.0)
}
unsigned int F_DRanS_R250_stream(statep *s)  // Fortran interface using derived type
{
  return(DRanS_R250_stream(s->p));
}

void VecIRan_R250_stream(void *stream, unsigned int *ranbuf, int n)  // !InTc!
{
  int k = 0;
  int i;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer
  unsigned int *r250_buffer ;
  int r250_index;

  r250_buffer = R250->buf ;
  r250_index = R250->cur;
  while( r250_index < 250 && n > 0 ){
    ranbuf[k++] = r250_buffer[r250_index++] ;
    n-- ;
  }
  R250->cur = r250_index;
  if ( n == 0 ) return;
  FillBuffer_R250_stream(R250);     // we get here if buffer is empty before n is satisfied
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = r250_buffer[i] ;
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250_stream(R250) ;
  }
  r250_index = R250->cur;
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = r250_buffer[r250_index++] ;
    n-- ;
  }
  R250->cur = r250_index;
}
void F_VecIRan_R250_stream(statep *s, unsigned int *ranbuf, int n)  // Fortran interface using derived type
{
  VecIRan_R250_stream(s->p,ranbuf,n);
}

void VecDRan_R250_stream(void *stream, double *ranbuf, int n)  // !InTc!
{
  int k = 0;
  int i;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer
  unsigned int *r250_buffer ;
  int r250_index;

  r250_buffer = R250->buf ;
  r250_index = R250->cur;
  while( r250_index < 250 && n > 0 ){
    ranbuf[k++] = CVTDBL_32(r250_buffer[r250_index++]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n-- ;
  }
  R250->cur = r250_index;
  if ( n == 0 ) return;
  FillBuffer_R250_stream(R250);     // we get here if buffer is empty
  r250_index = R250->cur;
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = CVTDBL_32(r250_buffer[i]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250_stream(R250) ;
  }
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = CVTDBL_32(r250_buffer[r250_index++]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n-- ;
  }
  R250->cur = r250_index;
}
void F_VecDRan_R250_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRan_R250_stream(s->p,ranbuf,n);
}

void VecDRanS_R250_stream(void *stream, double *ranbuf, int n)  // !InTc!
{
  int k = 0;
  int i;
  r250_state *R250 = stream ; //? (r250_state *) stream : &r250 ;   // use default stream if NULL stream pointer
  unsigned int *r250_buffer ;
  int r250_index;

  r250_buffer = R250->buf ;
  r250_index = R250->cur;
  while( r250_index < 250 && n > 0 ){
    ranbuf[k++] = CVTDBLS_32(r250_buffer[r250_index++]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n-- ;
  }
  R250->cur = r250_index;
  if ( n == 0 ) return;
  FillBuffer_R250_stream(R250);     // we get here if buffer is empty
  r250_index = R250->cur;
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = CVTDBLS_32(r250_buffer[i]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250_stream(R250) ;
  }
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = CVTDBLS_32(r250_buffer[r250_index++]) ;   // convert from 32 bit int to (0.0 , 1.0)
    n-- ;
  }
  R250->cur = r250_index;
}
void F_VecDRanS_R250_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRanS_R250_stream(s->p,ranbuf,n);
}

// ======================= internal static state =============================
static uint32_t sbuf[251] =
  {
    0x17617168 ,0x17a0d192 ,0x7cba449a ,0x86d91b38 ,0x7455bfb5 ,0x3bb194f2 ,0xd4cf89e2 ,0xee85f453 ,0x916c5ad3 ,0x8e5c32bb ,
    0x0d035201 ,0xc45d76fd ,0x04d799be ,0x409aa63c ,0x4c56a663 ,0xe0404838 ,0x88141c55 ,0xeb963899 ,0x83f34bc2 ,0xcf784dfe ,
    0x303c4f3d ,0xe2e580de ,0x9213c2a7 ,0x13a77f75 ,0xd0b5ca05 ,0xfc21d156 ,0x1e4c4c31 ,0xe4d70461 ,0xc8d3eb8a ,0x19d5a61d ,
    0x40d8f714 ,0x051457d9 ,0x60dd589d ,0xadcae8f4 ,0xae50f8e4 ,0xc0b1623d ,0x019840e2 ,0x9fd91e54 ,0x205fd793 ,0x9e17be66 ,
    0xc55b8dfa ,0x2fc0cfc5 ,0x7f292738 ,0x826e9e1f ,0x50e160e4 ,0xb3a2c76f ,0x5e03e011 ,0xe05601d4 ,0x5ef981a0 ,0x8676f0f1 ,
    0x6d8690cb ,0x8af13789 ,0xab00a515 ,0x2bcff371 ,0xd4adf100 ,0xdca5d50d ,0xc4ebad34 ,0x0d75d8f3 ,0xf77b9800 ,0x5593a16e ,
    0x1504475c ,0x453356ee ,0x1302dbd3 ,0xad1d21af ,0xcf668d00 ,0x76d68c1c ,0x91b61c7e ,0x831fbcb8 ,0x96258a86 ,0xbd1f4791 ,
    0xdf26ff85 ,0x93f5aa8d ,0xf6f6d072 ,0x55622d0f ,0x64905fb1 ,0x2d4b22c3 ,0x1ca372c3 ,0xda494a98 ,0xcfa8f513 ,0xf4737e2e ,
    0xe77f1a9c ,0x238ff343 ,0xafa9f948 ,0x20823fdf ,0x010c0e27 ,0x0020f353 ,0x31654bfd ,0x9637394d ,0x5ffd52b6 ,0x06db7185 ,
    0x76fadbb5 ,0xc90c3dd1 ,0x71364d2a ,0x49b411a8 ,0xda7556c9 ,0x86a61957 ,0xba798498 ,0x442d2c73 ,0xc89268b0 ,0x77d1be52 ,
    0x326e6b29 ,0x727a3092 ,0xba7b0780 ,0xf751def0 ,0xc1c05141 ,0x774f8101 ,0xde7495b4 ,0x2250658f ,0xd896e6e2 ,0x04a1649d ,
    0xf0dbdc2f ,0xab8a806e ,0xcccaef53 ,0x24a9daf8 ,0x9c472081 ,0xe88534a7 ,0xc6400afc ,0x5db66a6c ,0xf1bcdb68 ,0xe283ac5c ,
    0x71f6bb93 ,0x29bc3b72 ,0x6cc0d398 ,0x7ad7b430 ,0x07f8067d ,0x9bc21534 ,0xcad1a0be ,0x2eb81bf9 ,0x866e596e ,0xe6150f60 ,
    0x4f0cdac2 ,0xf14feeb3 ,0x63559552 ,0xee360d3e ,0xe950794d ,0x5674b5b7 ,0x8866636f ,0xb3bb5604 ,0xf894278a ,0xc4a13b2b ,
    0xc07f38ca ,0x1742de68 ,0xda9d902f ,0x2cb57fbc ,0x050fabca ,0x6b471777 ,0x4bdee5be ,0x60d3d78a ,0x6b31001b ,0x237156a5 ,
    0x617f2b07 ,0x4cbe6264 ,0x05610cd5 ,0x3d222b24 ,0x097ffff2 ,0x38975335 ,0x8682e0dc ,0x7e993479 ,0xe1ecf67b ,0xc619babd ,
    0x0f1ac989 ,0x1eab4e4b ,0x8c3cb3fc ,0x5787c983 ,0x74f19f89 ,0x968f257d ,0x95cc62b2 ,0x11e6bd09 ,0xd1a57e05 ,0x67358a7f ,
    0x95e23779 ,0x30efec41 ,0xe46c4803 ,0x2d2414e1 ,0x352c0fda ,0x1da1a740 ,0x28aea00b ,0xfe1dec28 ,0xae7b6c47 ,0xcfd1de31 ,
    0xa468360d ,0x544fc9e5 ,0xbcd04aa4 ,0xd2cfc38b ,0xb8a48f82 ,0xa8718902 ,0x5bd8a509 ,0x9c40dd86 ,0x6a3dadd0 ,0xd0a0d65f ,
    0xc62298c6 ,0x46393aca ,0x0b7436f2 ,0x99ddd69c ,0x839b79a7 ,0xa155be69 ,0x2e4f0458 ,0x474bd538 ,0x73d65578 ,0xa49ab70f ,
    0xbe2a3c0b ,0x69e550db ,0x9e38abcb ,0x9e483578 ,0xdabc5814 ,0x2e73f8ef ,0x4ed45df8 ,0x05f8d621 ,0x0259c01e ,0xf3927074 ,
    0xfda21b64 ,0x3476f241 ,0x9aa5d95a ,0xef86ea14 ,0x8f3fce06 ,0x8bff6bfa ,0x706ab0a2 ,0x7322f175 ,0x4e8acb27 ,0x336889cc ,
    0x373ea2e0 ,0x0cc5f5ce ,0x35a5cc68 ,0x93169549 ,0xea31a7b1 ,0x6a6569bc ,0xa776f509 ,0x5b0f310e ,0x96322244 ,0x64568c56 ,
    0x08aa6767 ,0x491799f1 ,0x17735c88 ,0x71c32f7e ,0xed0a2ec6 ,0xebd94777 ,0x9b1e1086 ,0xdc740f7a ,0x03c48151 ,0xafcb9f88 ,
    0xd835a40a ,0x21308fc2 ,0x0f459e5e ,0x0358b165 ,0x6422fa89 ,0xdd9cf11b ,0x03daccf5 ,0xec9e2bd9 ,0xe300013e ,0xa97d54e4 , 0
  };
  
static r250_state r250 = {
  (REFILLBUFFUN) FillBuffer_R250_stream, 
  (RANSETSEEDFUN) RanSetSeed_R250_stream, 
  (IRANFUN) IRan_R250_stream, 
  (DRANFUN) DRan_R250_stream, 
  (DRANSFUN) DRanS_R250_stream, 
  (IVECRANFUN) VecIRan_R250_stream, 
  (DVECRANFUN) VecDRan_R250_stream, 
  (DVECSRANFUN) VecDRanS_R250_stream, 
  &sbuf[0],
  -1,
  -1,
  NULL, 
  0
};

void *Ran_R250_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)   // !InTc! // create and seed a new stream
{
  r250_state *source ;
  r250_state *clone = (r250_state *)clone_in;
  int i;

  r250_state *new_state = (r250_state *) memalign(64,sizeof(r250_state)) ;

  if(cSeed < 0 && piSeed==NULL){            // clone a stream (mostly used for testing)
    source = clone ? clone : &r250;         // clone == NULL means clone default internal static stream
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
    new_state->buf = malloc(sizeof(uint32_t) * 251) ;
    new_state->cur = source->cur ;
    new_state->top = source->top ;
    for (i=0 ; i<250 ; i++) new_state->buf[i] = source->buf[i] ;
  }else{
    new_state->ngauss = 0;
    new_state->gauss = NULL;
    new_state->seed  = (RANSETSEEDFUN) RanSetSeed_R250_stream;
    new_state->refill = (REFILLBUFFUN) FillBuffer_R250_stream;
    new_state->iran   = (IRANFUN) IRan_R250_stream;
    new_state->dran   = (DRANFUN) DRan_R250_stream;
    new_state->drans  = (DRANSFUN) DRanS_R250_stream;
    new_state->vec_iran  = (IVECRANFUN) VecIRan_R250_stream;
    new_state->vec_dran  = (DVECRANFUN) VecDRan_R250_stream;
    new_state->vec_drans = (DVECSRANFUN) VecDRanS_R250_stream;
    new_state->buf = malloc(sizeof(uint32_t) * 251) ;
    new_state->cur = 250 ;
    new_state->top = 249 ;
    RanSetSeed_R250_stream(new_state, piSeed, cSeed);  // seed the new stream
  }

  return ( (void *) new_state) ;
}
void F_Ran_R250_new_stream(statep *s, statep *c, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  s->p = Ran_R250_new_stream( c->p, piSeed, cSeed);
}

// ======================= functions using internal static state =============================

static void FillBuffer_R250_static(){
  int i;
  unsigned int *r250_buffer = r250.buf ;
  int r250_index = r250.cur;

//   while(r250_index > 249) r250_index -= 250;
  r250_index = 0;
  for (i=0 ; i< 147 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i + 103 ];
  }
  for (i=147 ; i<250 ; i++) {
    r250_buffer[ i ] = r250_buffer[ i ] ^ r250_buffer[ i - 147 ];
  }
  r250.cur = r250_index;
}

void RanSetSeed_R250_static(unsigned int *piSeed, int cSeed)  // !InTc!
{
  int i;
  unsigned int *r250_buffer ;
  unsigned int seed ;

  r250_buffer = r250.buf ;
  if (cSeed == 250 && piSeed != NULL) {
    for(i=0 ; i<250; i++) r250_buffer[ i ] = piSeed[ i ] ;
  }else{
    seed = piSeed && (cSeed > 0) ? piSeed[0] : 0 ;
    for (i = 0 ; i < 250 ; ) {   /* see Knuth p.106, Table 1(16) and Numerical Recipes p.284 (ranqd1)*/
      seed = 1664525 * seed + 1013904223 ;
      if (seed <= 0) continue ;
      r250_buffer[i++] = seed ;
    }
  }
}

unsigned int IRan_R250_static()	  // !InTc!	/* returns a random unsigned integer */
{
  register unsigned int new_rand;

  if ( r250.cur > 249 ) FillBuffer_R250_static();
  new_rand = r250.buf[r250.cur++];
  return new_rand;
}

void VecIRan_R250_static(unsigned int *ranbuf, int n)  // !InTc!
{
  int k = 0;
  int i;
  unsigned int *r250_buffer ;
  int r250_index;

  r250_buffer = r250.buf ;
  r250_index = r250.cur;
  while( r250_index < 250 && n > 0 ){
    ranbuf[k++] = r250_buffer[r250_index++] ;
    n-- ;
  }
  r250.cur = r250_index;
  if ( n == 0 ) return;
  FillBuffer_R250_static();     // we get here if buffer is empty before n is satisfied
  while(n >= 250){        // chunks of 250 values
    for(i=0 ; i<250 ; i++) ranbuf[k+i] = r250_buffer[i] ;
    n -= 250 ;
    k += 250 ;
    FillBuffer_R250_static() ;
  }
  r250_index = r250.cur;
  while( n > 0 ){  // n < 250 at this point
    ranbuf[k++] = r250_buffer[r250_index++] ;
    n-- ;
  }
  r250.cur = r250_index;
}

/*------------------------- end of R250 routines ---------------------------*/
