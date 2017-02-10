  type, bind(C) :: RANDOM_STREAM
    type(C_PTR) :: p
  end type
! void RanSetSeed_R250_stream(void *stream, int *piSeed, int cSeed)
 interface
   subroutine RanSetSeed_R250_stream(stream, piSeed, cSeed) bind(C,name='RanSetSeed_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanSetSeed_R250_stream
 end interface
! r250_state *Ran_R250_new_stream(r250_state *clone, int *piSeed, int cSeed)
 interface
   function Ran_R250_new_stream(clone, piSeed, cSeed) result(stream) bind(C,name='Ran_R250_new_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: clone
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   type(RANDOM_STREAM) :: stream
   end function Ran_R250_new_stream
 end interface
! unsigned int IRan_R250_stream(void *stream)
 interface
   function IRan_R250_stream(stream) result(ran) bind(C,name='IRan_R250_stream')
   import :: C_INT,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT) :: ran
   end function IRan_R250_stream
 end interface
! double DRan_R250_stream(void *stream)
 interface
   function DRan_R250_stream(stream) result(ran) bind(C,name='DRan_R250_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN), value :: stream
   real(C_DOUBLE) :: ran
   end function DRan_R250_stream
 end interface
! double DRanS_R250_stream(void *stream)
 interface
   function DRanS_R250_stream(stream) result(ran) bind(C,name='DRanS_R250_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN), value :: stream
   real(C_DOUBLE) :: ran
   end function DRanS_R250_stream
 end interface
! void VecIRan_R250_stream(void *stream, unsigned int *ranbuf, int n)
 interface
   subroutine VecIRan_R250_stream(stream, ranbuf, n) bind(C,name='VecIRan_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecIRan_R250_stream
 end interface
! void VecDRanS_R250_stream(void *stream, double *ranbuf, int n)
 interface
   subroutine VecDRanS_R250_stream(stream, ranbuf, n) bind(C,name='VecDRanS_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecDRanS_R250_stream
 end interface
! void VecDRan_R250_stream(void *stream, double *ranbuf, int n)
 interface
   subroutine VecDRan_R250_stream(stream, ranbuf, n) bind(C,name='VecDRan_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecDRan_R250_stream
 end interface
! double DRanNormalZigVec(void *stream)
 interface
   function DRanNormalZigVec(stream) result(ran) bind(C,name='DRanNormalZigVec')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN), value :: stream
   real(C_DOUBLE) :: ran
   end function DRanNormalZigVec
 end interface
! void RanNormalSetSeedZig(void *stream, int *piSeed, int cSeed)
 interface
   subroutine RanNormalSetSeedZig(stream, piSeed, cSeed) bind(C,name='RanNormalSetSeedZig')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanNormalSetSeedZig
 end interface
! void RanNormalSetSeedZigVec(void *stream, int *piSeed, int cSeed)
 interface
   subroutine RanNormalSetSeedZigVec(stream, piSeed, cSeed) bind(C,name='RanNormalSetSeedZigVec')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN), value :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanNormalSetSeedZigVec
 end interface
