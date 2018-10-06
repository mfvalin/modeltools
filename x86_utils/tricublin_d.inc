interface
  subroutine v2cubic_coeffs_d(pxy, xy) bind(C,name='V2cubic_coeffs_d')
    import C_DOUBLE
    real(C_DOUBLE), dimension(8), intent(OUT) :: pxy
    real(C_DOUBLE), dimension(2), intent(IN) :: xy
  end subroutine v2cubic_coeffs_d
  subroutine cubic_coeffs_d(px, x) bind(C,name='Cubic_coeffs_d')
    import C_DOUBLE
    real(C_DOUBLE), dimension(4), intent(OUT) :: px
    real(C_DOUBLE), intent(IN), value :: x
  end subroutine cubic_coeffs_d
  subroutine bicubic_coeffs_d(px, py, x, y) bind(C,name='Bicubic_coeffs_d')
    import C_DOUBLE
    real(C_DOUBLE), dimension(4), intent(OUT) :: px, py
    real(C_DOUBLE), intent(IN), value :: x, y
  end subroutine bicubic_coeffs_d
  subroutine tricubic_coeffs_d(px, py, pz, x, y, z) bind(C,name='Tricubic_coeffs_d')
    import C_DOUBLE
    real(C_DOUBLE), dimension(4), intent(OUT) :: px, py, pz
    real(C_DOUBLE), intent(IN), value :: x, y, z
  end subroutine tricubic_coeffs_d
  subroutine tricublin_zyx3f_d(d, f, px, py, pz, ni, ninj) bind(C,name='Tricublin_zyx3f_d')
    import :: C_DOUBLE, C_FLOAT, C_INT
    real(C_FLOAT), dimension(3), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyx3f_d
  subroutine tricublin_zyxf3_d(d, f1, f2, f3, px, py, pz, ni, ninj) bind(C,name='Tricublin_zyxf3_d')
    import :: C_DOUBLE, C_FLOAT, C_INT
    real(C_FLOAT), dimension(3), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f1, f2, f3
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyxf3_d
  subroutine tricublin_zyxf_d(d, f, px, py, pz, ni, ninj) bind(C,name='Tricublin_zyxf_d')
    import :: C_DOUBLE, C_FLOAT, C_INT
    real(C_FLOAT), intent(OUT) :: d
    real(C_FLOAT), dimension(*), intent(IN) :: f
    real(C_DOUBLE), dimension(4), intent(IN) :: px, py, pz
    integer(C_INT), intent(IN), value :: ni, ninj
  end subroutine tricublin_zyxf_d
end interface
