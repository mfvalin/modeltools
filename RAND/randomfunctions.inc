  type, bind(C) :: RANDOM_STREAM
    type(C_PTR) :: p
  end type
! void F_RanSetSeed_generic_stream(statep *s   , int *piSeed, int cSeed)
 interface
   subroutine RanSetSeed_generic_stream(stream, piSeed, cSeed) bind(C,name='F_RanSetSeed_generic_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanSetSeed_generic_stream
 end interface
! unsigned int F_IRan_generic_stream(statep *s   )
 interface
   function IRan_generic_stream(stream) result(ran) bind(C,name='F_IRan_generic_stream')
   import :: C_INT,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT) :: ran
   end function IRan_generic_stream
 end interface
! double F_DRan_generic_stream(statep *s   )
 interface
   function DRan_generic_stream(stream) result(ran) bind(C,name='F_DRan_generic_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function DRan_generic_stream
 end interface
! double F_DRanS_generic_stream(statep *s   )
 interface
   function DRanS_generic_stream(stream) result(ran) bind(C,name='F_DRanS_generic_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function DRanS_generic_stream
 end interface
! void F_VecIRan_generic_stream(statep *s   , unsigned int *ranbuf, int n)
 interface
   subroutine VecIRan_generic_stream(stream, ranbuf, n) bind(C,name='F_VecIRan_generic_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecIRan_generic_stream
 end interface
! void F_VecDRanS_generic_stream(statep *s   , double *ranbuf, int n)
 interface
   subroutine VecDRanS_generic_stream(stream, ranbuf, n) bind(C,name='F_VecDRanS_generic_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecDRanS_generic_stream
 end interface
! void F_VecDRan_generic_stream(statep *s   , double *ranbuf, int n)
 interface
   subroutine VecDRan_generic_stream(stream, ranbuf, n) bind(C,name='F_VecDRan_generic_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: n
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf
   end subroutine VecDRan_generic_stream
 end interface
! double F_DRan_NormalZig_stream(statep *s   )
 interface
   function DRan_Normal_stream(stream) result(ran) bind(C,name='F_DRan_NormalZig_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function DRan_Normal_stream
 end interface
! double F_D64Ran_NormalZig_stream(statep *s   )
 interface
   function D64Ran_Normal_stream(stream) result(ran) bind(C,name='F_D64Ran_NormalZig_stream')
   import :: C_DOUBLE,RANDOM_STREAM
   type(RANDOM_STREAM), intent(IN) :: stream
   real(C_DOUBLE) :: ran
   end function D64Ran_Normal_stream
 end interface
! void F_RanNormalZigSetSeed(statep *s   , int *piSeed, int cSeed)
 interface
   subroutine RanNormalZigSetSeed(stream, piSeed, cSeed) bind(C,name='F_RanNormalZigSetSeed')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(IN) :: stream
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine RanNormalZigSetSeed
 end interface
! void F_Ran_MT19937_new_stream(mt19937_state *clone, int *piSeed, int cSeed)
 interface
   subroutine Ran_MT19937_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_MT19937_new_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(OUT) :: stream
   type(RANDOM_STREAM), intent(IN) :: clone
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine Ran_MT19937_new_stream
 end interface
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
! void F_Ran_SHR3_new_stream(shr3_state *clone, int *piSeed, int cSeed)
 interface
   subroutine Ran_SHR3_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_SHR3_new_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(OUT) :: stream
   type(RANDOM_STREAM), intent(IN) :: clone
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine Ran_SHR3_new_stream
 end interface
! void F_Ran_XSR128_new_stream(xsr128_state *clone, int *piSeed, int cSeed)
 interface
   subroutine Ran_XSR128_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_XSR128_new_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(OUT) :: stream
   type(RANDOM_STREAM), intent(IN) :: clone
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine Ran_XSR128_new_stream
 end interface
! void F_Ran_XSR128R_new_stream(xsr128r_state *clone, int *piSeed, int cSeed)
 interface
   subroutine Ran_XSR128R_new_stream(stream, clone, piSeed, cSeed) bind(C,name='F_Ran_XSR128R_new_stream')
   import :: RANDOM_STREAM,C_INT
   type(RANDOM_STREAM), intent(OUT) :: stream
   type(RANDOM_STREAM), intent(IN) :: clone
   integer(C_INT), intent(IN), value :: cSeed
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
   end subroutine Ran_XSR128R_new_stream
 end interface
