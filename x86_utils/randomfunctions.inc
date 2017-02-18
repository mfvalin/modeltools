  type, bind(C) :: RANDOM_STREAM
    type(C_PTR) :: p
  end type
! void F_RanSetSeed_R250_stream(statep *s   , int *piSeed, int cSeed)
 interface
   subroutine RanSetSeed_R250_stream(stream, piSeed, cSeed) bind(C,name='F_RanSetSeed_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanSetSeed_R250_stream
 end interface
! void F_Ran_R250_new_stream(r250_state *clone, int *piSeed, int cSeed)
 interface
   subroutine Ran_R250_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_R250_new_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(OUT) :: stream
   type(RANDOM_STREAM), intent(IN) :: clone
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine Ran_R250_new_stream
 end interface
! unsigned int F_IRan_R250_stream(statep *s   )
 interface
   function IRan_R250_stream(stream) result(ran) bind(C,name='F_IRan_R250_stream')
   import :: C_INT,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT) :: ran
   end function IRan_R250_stream
 end interface
! double F_DRan_R250_stream(statep *s   )
 interface
   function DRan_R250_stream(stream) result(ran) bind(C,name='F_DRan_R250_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function DRan_R250_stream
 end interface
! double F_DRanS_R250_stream(statep *s   )
 interface
   function DRanS_R250_stream(stream) result(ran) bind(C,name='F_DRanS_R250_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function DRanS_R250_stream
 end interface
! void F_VecIRan_R250_stream(statep *s   , unsigned int *ranbuf, int n)
 interface
   subroutine VecIRan_R250_stream(stream, ranbuf, n) bind(C,name='F_VecIRan_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecIRan_R250_stream
 end interface
! void F_VecDRanS_R250_stream(statep *s   , double *ranbuf, int n)
 interface
   subroutine VecDRanS_R250_stream(stream, ranbuf, n) bind(C,name='F_VecDRanS_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecDRanS_R250_stream
 end interface
! void F_VecDRan_R250_stream(statep *s   , double *ranbuf, int n)
 interface
   subroutine VecDRan_R250_stream(stream, ranbuf, n) bind(C,name='F_VecDRan_R250_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecDRan_R250_stream
 end interface
! double F_DRanNormalZigVec(statep *s   )
 interface
   function DRanNormalZigVec(stream) result(ran) bind(C,name='F_DRanNormalZigVec')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function DRanNormalZigVec
 end interface
! void F_RanNormalSetSeedZig(statep *s   , int *piSeed, int cSeed)
 interface
   subroutine RanNormalSetSeedZig(stream, piSeed, cSeed) bind(C,name='F_RanNormalSetSeedZig')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanNormalSetSeedZig
 end interface
! void F_RanNormalSetSeedZigVec(statep *s   , int *piSeed, int cSeed)
 interface
   subroutine RanNormalSetSeedZigVec(stream, piSeed, cSeed) bind(C,name='F_RanNormalSetSeedZigVec')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanNormalSetSeedZigVec
 end interface
