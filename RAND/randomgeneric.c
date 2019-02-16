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

// generic interface to random functions using streams  // !InTc!

#include <randomgeneric.h>

// Fortran interfaces for automated extraction
#if defined(NEVER_TO_BE_TRUE)

  type, bind(C) :: RANDOM_STREAM                                                          !InTf!
    type(C_PTR) :: p                                                                      !InTf!
  end type                                                                                !InTf!

! void F_RanSetSeed_generic_stream(statep *s   , int *piSeed, int cSeed)                  !InTf!
 interface                                                                                !InTf!
   subroutine RanSetSeed_generic_stream(stream, piSeed, cSeed) bind(C,name='F_RanSetSeed_generic_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanSetSeed_generic_stream                                               !InTf!
 end interface                                                                            !InTf!

! unsigned int F_IRan_generic_stream(statep *s   )                                        !InTf!
 interface                                                                                !InTf!
   function IRan_generic_stream(stream) result(ran) bind(C,name='F_IRan_generic_stream')  !InTf!
   import :: C_INT,RANDOM_STREAM                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT) :: ran                                                                  !InTf!
   end function IRan_generic_stream                                                       !InTf!
 end interface                                                                            !InTf!

! double F_DRan_generic_stream(statep *s   )                                              !InTf!
 interface                                                                                !InTf!
   function DRan_generic_stream(stream) result(ran) bind(C,name='F_DRan_generic_stream')  !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRan_generic_stream                                                       !InTf!
 end interface                                                                            !InTf!

! double F_DRanS_generic_stream(statep *s   )                                             !InTf!
 interface                                                                                !InTf!
   function DRanS_generic_stream(stream) result(ran) bind(C,name='F_DRanS_generic_stream')  !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRanS_generic_stream                                                      !InTf!
 end interface                                                                            !InTf!

! void F_VecIRan_generic_stream(statep *s   , unsigned int *ranbuf, int n)                !InTf!
 interface                                                                                !InTf!
   subroutine VecIRan_generic_stream(stream, ranbuf, n) bind(C,name='F_VecIRan_generic_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecIRan_generic_stream                                                  !InTf!
 end interface                                                                            !InTf!

! void F_VecDRanS_generic_stream(statep *s   , double *ranbuf, int n)                     !InTf!
 interface                                                                                !InTf!
   subroutine VecDRanS_generic_stream(stream, ranbuf, n) bind(C,name='F_VecDRanS_generic_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecDRanS_generic_stream                                                 !InTf!
 end interface                                                                            !InTf!

! void F_VecDRan_generic_stream(statep *s   , double *ranbuf, int n)                      !InTf!
 interface                                                                                !InTf!
   subroutine VecDRan_generic_stream(stream, ranbuf, n) bind(C,name='F_VecDRan_generic_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecDRan_generic_stream                                                  !InTf!
 end interface                                                                            !InTf!
#endif

void RanSetSeed_generic_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  // !InTc!
{
  generic_state *state = stream ;
  state->seed(stream, piSeed, cSeed);
}
void F_RanSetSeed_generic_stream(statep *s, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  RanSetSeed_generic_stream( s->p, piSeed, cSeed);
}

uint32_t IRan_generic_stream(generic_state *stream)       // !InTc!
{
//   generic_state *state = stream ;
//   return state->iran(stream);
  uint32_t value;
// fprintf(stderr,"IRan_generic_stream\n");
  if(stream->cur > stream->top) stream->refill(stream);
  value = stream->buf[stream->cur] ;
  stream->cur = stream->cur + 1;
  return value;
}
uint32_t F_IRan_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(IRan_generic_stream(s->p));
}

double DRan_generic_stream(generic_state *stream)       // !InTc!
{
//   generic_state *state = stream ;
//   return state->dran(stream);
  uint32_t value;
  if(stream->cur > stream->top) stream->refill(stream);
  value = stream->buf[stream->cur] ;
  stream->cur = stream->cur + 1;
  return CVTDBL_32(value);
}
double F_DRan_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(DRan_generic_stream(s->p));
}

double DRanS_generic_stream(generic_state *stream)       // !InTc!
{
//   generic_state *state = stream ;
//   return state->drans(stream);
  uint32_t value;
  if(stream->cur > stream->top) stream->refill(stream);
  value = stream->buf[stream->cur] ;
  stream->cur = stream->cur + 1;
  return CVTDBLS_32(value);
}
double F_DRanS_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(DRanS_generic_stream(s->p));
}

void VecIRan_generic_stream(generic_state *stream, unsigned int *ranbuf, int n)       // !InTc!
{
//   generic_state *state = stream ;
//   state->vec_iran(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  int *buf = stream->buf;
  int i;
  while(n > 0 && cur < topp1){
    *ranbuf = buf[cur] ; ranbuf++ ; cur++; n-- ;
  }
  if(cur < topp1) {
    stream->cur = cur ;
  }else{
    stream->refill(stream);   // this will set stream->cur
  }
  if(n <= 0) return ;         // done
  while(n > 0){
    if(n >= topp1){             // full buffers
      for(cur=0 ; cur<topp1 ; cur++) {ranbuf[cur] = buf[cur] ; }
      ranbuf += topp1;
      stream->refill(stream);
      n -= topp1;
    }else{                    // less than 1 buffer, no need for refill
      for(cur=0 ; cur<n ; cur++) {ranbuf[cur] = buf[cur] ; }
      stream->cur = cur ;     // store current pointer
      n = 0;
    }
  }
}
void F_VecIRan_generic_stream(statep *s, unsigned int *ranbuf, int n)  // Fortran interface using derived type
{
  VecIRan_generic_stream(s->p,ranbuf,n);
}

void VecDRan_generic_stream(generic_state *stream, double *ranbuf, int n)       // !InTc!
{
//   generic_state *state = stream ;
//   state->vec_dran(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  int *buf = stream->buf;
  int i;
  while(n > 0 && cur < topp1){
    *ranbuf = CVTDBL_32(buf[cur]) ; ranbuf++ ; cur++; n-- ;
  }
  if(cur < topp1) {
    stream->cur = cur ;
  }else{
    stream->refill(stream);   // this will set stream->cur
  }
  if(n <= 0) return ;         // done
  while(n > 0){
    if(n >= topp1){             // full buffers
      for(cur=0 ; cur<topp1 ; cur++) {ranbuf[cur] = CVTDBL_32(buf[cur]) ; }
      ranbuf += topp1;
      stream->refill(stream);
      n -= topp1;
    }else{                    // less than 1 buffer, no need for refill
      for(cur=0 ; cur<n ; cur++) {ranbuf[cur] = CVTDBL_32(buf[cur]) ; }
      stream->cur = cur ;     // store current pointer
      n = 0;
    }
  }
}
void F_VecDRan_generic_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRan_generic_stream(s->p,ranbuf,n);
}

void VecDRanS_generic_stream(generic_state *stream, double *ranbuf, int n)       // !InTc!
{
//   generic_state *state = stream ;
//   state->vec_drans(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  int *buf = stream->buf;
  int i;
  while(n > 0 && cur < topp1){
    *ranbuf = CVTDBLS_32(buf[cur]) ; ranbuf++ ; cur++; n-- ;
  }
  if(cur < topp1) {
    stream->cur = cur ;
  }else{
    stream->refill(stream);   // this will set stream->cur
  }
  if(n <= 0) return ;         // done
  while(n > 0){
    if(n >= topp1){             // full buffers
      for(cur=0 ; cur<topp1 ; cur++) {ranbuf[cur] = CVTDBLS_32(buf[cur]) ; }
      ranbuf += topp1;
      stream->refill(stream);
      n -= topp1;
    }else{                    // less than 1 buffer, no need for refill
      for(cur=0 ; cur<n ; cur++) {ranbuf[cur] = CVTDBLS_32(buf[cur]) ; }
      stream->cur = cur ;     // store current pointer
      n = 0;
    }
  }
}
void F_VecDRanS_generic_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRanS_generic_stream(s->p,ranbuf,n);
}


