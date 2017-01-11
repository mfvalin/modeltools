interface
! void RanSetSeed_R250_stream(r250_state *R250_, int *piSeed, int cSeed)                  !InTf!
   subroutine set_r250_seed(stream, seeds, nseeds) bind(C,name='RanSetSeed_R250_stream')  !InTf!
   import :: C_PTR,C_INT                                                                  !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   integer(C_INT), intent(IN), value :: nseeds                                            !InTf!
   integer(c_INT), dimension(nseeds), intent(IN) :: seeds                                 !InTf!
   end subroutine set_r250_seed                                                           !InTf!
! r250_state *Ran_R250_new_stream(r250_state *clone, int *piSeed, int cSeed)              !InTf!
   function new_r250_stream(clone, seeds, nseeds) result(stream) bind(C,name='Ran_R250_new_stream')  !InTf!
   import :: C_PTR,C_INT                                                                  !InTf!
   type(C_PTR), intent(IN) :: clone                                                       !InTf!
   integer(C_INT), intent(IN), value :: nseeds                                            !InTf!
   integer(c_INT), dimension(nseeds), intent(IN) :: seeds                                 !InTf!
   type(C_PTR) :: stream                                                                  !InTf!
   end function new_r250_stream                                                           !InTf!
! double DRan_R250_stream(r250_state *R250_)                                              !InTf!
   function dran_r250(stream) result(ran) bind(C,name='DRan_R250_stream')                 !InTf!
   import :: C_DOUBLE,C_PTR                                                               !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function dran_r250                                                                 !InTf!
! double DRanS_R250_stream(r250_state *R250_)                                             !InTf!
   function drans_r250(stream) result(ran) bind(C,name='DRanS_R250_stream')               !InTf!
   import :: C_DOUBLE,C_PTR                                                               !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function drans_r250                                                                !InTf!
! unsigned int IRan_R250_stream(r250_state *R250_)                                        !InTf!
   function iran_r250(stream) result(ran) bind(C,name='IRan_R250_stream')                 !InTf!
   import :: C_INT,C_PTR                                                                  !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   integer(C_INT) :: ran                                                                  !InTf!
   end function iran_r250                                                                 !InTf!
! void VecIRan_R250_stream(r250_state *R250_, unsigned int *ranbuf, int n)                !InTf!
   subroutine v_iran_r250(stream, buf, nval) bind(C,name='VecIRan_R250_stream')           !InTf!
   import :: C_PTR,C_INT                                                                  !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   integer(C_INT), intent(IN), value :: nval                                              !InTf!
   integer(c_INT), dimension(nval), intent(OUT) :: buf                                    !InTf!
   end subroutine v_iran_r250                                                             !InTf!
! void VecDRan_R250_stream(r250_state *R250_, unsigned int *ranbuf, int n)                !InTf!
   subroutine v_dran_r250(stream, buf, nval) bind(C,name='VecDRan_R250_stream')           !InTf!
   import :: C_PTR,C_INT                                                                  !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   integer(C_INT), intent(IN), value :: nval                                              !InTf!
   integer(c_INT), dimension(nval), intent(OUT) :: buf                                    !InTf!
   end subroutine v_dran_r250                                                             !InTf!
! void VecDRans_R250_stream(r250_state *R250_, unsigned int *ranbuf, int n)               !InTf!
   subroutine v_drans_r250(stream, buf, nval) bind(C,name='VecDRanS_R250_stream')         !InTf!
   import :: C_PTR,C_INT                                                                  !InTf!
   type(C_PTR), intent(IN) :: stream                                                      !InTf!
   integer(C_INT), intent(IN), value :: nval                                              !InTf!
   integer(c_INT), dimension(nval), intent(OUT) :: buf                                    !InTf!
   end subroutine v_drans_r250                                                            !InTf!
end interface
