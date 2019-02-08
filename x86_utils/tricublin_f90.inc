  interface                                                                      !InTf!
    function vsearch_setup(levels, nk, ni, nj) result (ztab) bind(C,name='Vsearch_setup')     !InTf!
      import :: C_PTR, C_DOUBLE, C_INT                                           !InTf!
      real(C_DOUBLE), dimension(nk), intent(IN) :: levels                        !InTf!
      integer(C_INT), intent(IN), value :: nk, ni, nj                            !InTf!
      type(C_PTR) :: ztab                                                        !InTf!
    end function vsearch_setup                                                   !InTf!
    function vsearch_setup_plus(levels, nk, ni, nj, offseti, offsetj) result (ztab) bind(C,name='Vsearch_setup_plus')     !InTf!
      import :: C_PTR, C_DOUBLE, C_INT                                           !InTf!
      real(C_DOUBLE), dimension(nk), intent(IN) :: levels                        !InTf!
      integer(C_INT), intent(IN), value :: nk, ni, nj, offseti, offsetj          !InTf!
      type(C_PTR) :: ztab                                                        !InTf!
    end function vsearch_setup_plus                                              !InTf!
    subroutine tricublin_zyx1_n(d,f,pxyz,lv,n) bind(C,name='Tricublin_zyx1_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d                                     !InTf!
      real, dimension(*), intent(IN)    :: f                                     !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n                                     !InTf!
    end subroutine tricublin_zyx1_n                                              !InTf!
    subroutine tricublin_mono_zyx_n(d,l,mi,ma,f,pxyz,lv,n) bind(C,name='Tricublin_mono_zyx_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d, l, mi, ma                          !InTf!
      real, dimension(*), intent(IN)    :: f                                     !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n                                     !InTf!
    end subroutine tricublin_mono_zyx_n                                              !InTf!
    subroutine tricublin_zyx3_n(d,f,pxyz,lv,n) bind(C,name='Tricublin_zyx3_n')   !InTf!
      import :: C_PTR                                                            !InTf!
      real, dimension(*), intent(OUT)   :: d                                     !InTf!
      real, dimension(*), intent(IN)    :: f                                     !InTf!
      real, dimension(*), intent(IN)  :: pxyz                                    !InTf!
      type(C_PTR), intent(IN), value    :: lv                                    !InTf!
      integer, intent(IN), value        :: n                                     !InTf!
    end subroutine tricublin_zyx3_n                                              !InTf!
  end interface                                                                  !InTf!
