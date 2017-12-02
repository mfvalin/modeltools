#include <randomgeneric.h>

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

unsigned int IRan_generic_stream(generic_state *stream)       // !InTc!
{
  generic_state *state = stream ;
  return state->iran(stream);
}
unsigned int F_IRan_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(IRan_generic_stream(s->p));
}

double DRan_generic_stream(generic_state *stream)       // !InTc!
{
  generic_state *state = stream ;
  return state->dran(stream);
}
double F_DRan_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(DRan_generic_stream(s->p));
}

double DRanS_generic_stream(generic_state *stream)       // !InTc!
{
  generic_state *state = stream ;
  return state->drans(stream);
}
double F_DRanS_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(DRanS_generic_stream(s->p));
}

void VecIRan_generic_stream(generic_state *stream, unsigned int *ranbuf, int n)       // !InTc!
{
  generic_state *state = stream ;
  state->vec_iran(stream,ranbuf,n);
}
void F_VecIRan_generic_stream(statep *s, unsigned int *ranbuf, int n)  // Fortran interface using derived type
{
  VecIRan_generic_stream(s->p,ranbuf,n);
}

void VecDRan_generic_stream(generic_state *stream, double *ranbuf, int n)       // !InTc!
{
  generic_state *state = stream ;
  state->vec_dran(stream,ranbuf,n);
}
void F_VecDRan_generic_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRan_generic_stream(s->p,ranbuf,n);
}

void VecDRanS_generic_stream(generic_state *stream, double *ranbuf, int n)       // !InTc!
{
  generic_state *state = stream ;
  state->vec_drans(stream,ranbuf,n);
}
void F_VecDRanS_generic_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRanS_generic_stream(s->p,ranbuf,n);
}


