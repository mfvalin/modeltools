interface
  subroutine tricublin_zyx1_n(d, f, pxpypz, lv, n) bind(C,name='Tricublin_zyx1_n')
    import :: C_DOUBLE, C_FLOAT, C_INT, C_PTR
    real(C_FLOAT), dimension(*), intent(OUT) :: d         ! real dimension is (n)
    real(C_FLOAT), dimension(*), intent(IN) :: f          ! real dimension is (ni,nj,nk) as in lv
    real(C_DOUBLE), dimension(*), intent(IN) :: pxpypz    ! real dimension is (3,n)
    type(C_PTR), intent(IN), value :: lv
    integer(C_INT), intent(IN), value :: n
  end subroutine tricublin_zyx1_n

  function search_list_setup(zaxis, nk, ni, nj) result(lv) bind(C,name='Vsearch_setup')
    import :: C_DOUBLE, C_INT, C_PTR
    integer(C_INT), intent(IN), value :: nk, ni, nj
    real(C_DOUBLE), dimension(nk), intent(IN) :: zaxis
    type(C_PTR) :: lv
  end function search_list_setup

  subroutine tricublin_zyxf_beta(d, f, px, py, pz, ni, ninj) bind(C,name='Tricublin_zyxf_beta')
    import :: C_DOUBLE, C_FLOAT, C_INT
    real(C_FLOAT), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyxf_beta

  subroutine tricubic_coeffs_d(px, py, pz, x, y, z) bind(C,name='Tricubic_coeffs_d')
    import C_DOUBLE
    real(C_DOUBLE), dimension(4), intent(OUT) :: px, py, pz
    real(C_DOUBLE), intent(IN), value :: x, y, z
  end subroutine tricubic_coeffs_d

  subroutine tricublin_zyx3f_beta(d, f, px, py, pz, ni, ninj) bind(C,name='Tricublin_zyx3f_beta')
    import :: C_DOUBLE, C_FLOAT, C_INT
    real(C_FLOAT), dimension(3), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyx3f_beta

  function search_list_inc_n(ix, zaxis, lv, n) result(indx) bind(C,name='Vsearch_list_inc_n')
    import :: C_DOUBLE, C_INT, C_PTR
    integer(C_INT), intent(IN), value :: n
    type(C_PTR), intent(IN), value :: lv
    integer(C_INT), dimension(n), intent(OUT) :: ix
    real(C_DOUBLE), dimension(n), intent(IN) :: zaxis
    integer(C_INT) :: indx
  end function search_list_inc_n

  function search_list_inc_dn(posz, zaxis, lv, n) result(last) bind(C,name='Vsearch_list_inc_dn')
    import :: C_DOUBLE, C_INT, C_PTR
    integer(C_INT), intent(IN), value :: n
    type(C_PTR), intent(IN), value :: lv
    real(C_DOUBLE), dimension(n), intent(IN) :: zaxis
    real(C_DOUBLE), dimension(n), intent(OUT) :: posz
    real(C_DOUBLE) :: last
  end function search_list_inc_dn

  function coef_xyz_incr(ixyz, cxyz, PX, PY, PZ, lv, n) result(indx) bind(C,name='Vcoef_xyz_incr')
    import :: C_DOUBLE, C_INT, C_PTR
    integer(C_INT), intent(IN), value :: n
    type(C_PTR), intent(IN), value :: lv
    real(C_DOUBLE), dimension(24,n), intent(OUT) :: cxyz
    integer(C_INT), dimension(n), intent(OUT) :: ixyz
    real(C_DOUBLE), dimension(n), intent(IN) :: PX, PY, PZ
    integer(C_INT) :: indx
  end function coef_xyz_incr

  function coef_xyz_inc(ixyz, cxyz, PX, PY, PZ, lv, n) result(indx) bind(C,name='Vcoef_xyz_inc')
    import :: C_DOUBLE, C_INT, C_PTR
    integer(C_INT), intent(IN), value :: n
    type(C_PTR), intent(IN), value :: lv
    real(C_DOUBLE), dimension(24,n), intent(OUT) :: cxyz
    integer(C_INT), dimension(n), intent(OUT) :: ixyz
    real(C_DOUBLE), dimension(n), intent(IN) :: PX, PY, PZ
    integer(C_INT) :: indx
  end function coef_xyz_inc

end interface
