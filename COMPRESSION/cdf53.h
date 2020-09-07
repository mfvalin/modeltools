#if defined(FORTRAN_SOURCE)
 interface
   subroutine F_CDF53_1D_split_N_even(x, e, o, n) bind(C,name='F_CDF53_1D_split_N_even')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(IN) :: x
     real(C_FLOAT), dimension(*), intent(OUT) :: e, o
   end subroutine F_CDF53_1D_split_N_even
 end interface
 interface
   subroutine F_CDF53_1D_inplace_N_even(x, e, o, n) bind(C,name='F_CDF53_1D_inplace_N_even')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine F_CDF53_1D_inplace_N_even
 end interface
 interface
   subroutine F_CDF53_1D_split_N_odd(x, e, o, n) bind(C,name='F_CDF53_1D_split_N_odd')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(IN) :: x
     real(C_FLOAT), dimension(*), intent(OUT) :: e, o
   end subroutine F_CDF53_1D_split_N_odd
 end interface
 interface
   subroutine F_CDF53_1D_inplace_N_odd(x, e, o, n) bind(C,name='F_CDF53_1D_inplace_N_odd')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine F_CDF53_1D_inplace_N_odd
 end interface
 interface
   subroutine F_CDF53_1D_split(x, e, o, n) bind(C,name='F_CDF53_1D_split')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(IN) :: x
     real(C_FLOAT), dimension(*), intent(OUT) :: e, o
   end subroutine F_CDF53_1D_split
 end interface
 interface
   subroutine F_CDF53_1D_inplace(x, e, o, n) bind(C,name='F_CDF53_1D_inplace')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine F_CDF53_1D_inplace
 end interface
 interface
   subroutine F_CDF53_1D_split_inplace(x, e, o, n) bind(C,name='F_CDF53_1D_split_inplace')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine F_CDF53_1D_split_inplace
 end interface
 interface
   subroutine I_CDF53_1D_split_N_even(x, e, o, n) bind(C,name='I_CDF53_1D_split_N_even')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(OUT) :: x
     real(C_FLOAT), dimension(*), intent(IN) :: e, o
   end subroutine I_CDF53_1D_split_N_even
 end interface
 interface
   subroutine I_CDF53_1D_inplace_N_even(x, e, o, n) bind(C,name='I_CDF53_1D_inplace_N_even')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine I_CDF53_1D_inplace_N_even
 end interface
 interface
   subroutine I_CDF53_1D_split_N_odd(x, e, o, n) bind(C,name='I_CDF53_1D_split_N_odd')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(OUT) :: x
     real(C_FLOAT), dimension(*), intent(IN) :: e, o
   end subroutine I_CDF53_1D_split_N_odd
 end interface
 interface
   subroutine I_CDF53_1D_inplace_N_odd(x, e, o, n) bind(C,name='I_CDF53_1D_inplace_N_odd')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine I_CDF53_1D_inplace_N_odd
 end interface
 interface
   subroutine I_CDF53_1D_split(x, e, o, n) bind(C,name='I_CDF53_1D_split')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(OUT) :: x
     real(C_FLOAT), dimension(*), intent(IN) :: e, o
   end subroutine I_CDF53_1D_split
 end interface
 interface
   subroutine I_CDF53_1D_split_inplace(x, e, o, n) bind(C,name='I_CDF53_1D_split_inplace')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine I_CDF53_1D_split_inplace
 end interface
 interface
   subroutine I_CDF53_1D_inplace(x, e, o, n) bind(C,name='I_CDF53_1D_inplace')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: n
     real(C_FLOAT), dimension(n), intent(INOUT) :: x
   end subroutine I_CDF53_1D_inplace
 end interface
 interface
 subroutine I_CDF53_2D_split_inplace(x, ni, lni, nj) BIND(C,name='I_CDF53_2D_split_inplace')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: ni, nj, lni
     real(C_FLOAT), dimension(*), intent(INOUT) :: x
 end subroutine I_CDF53_2D_split_inplace
 end interface
 interface
 subroutine F_CDF53_2D_split_inplace(x, ni, lni, nj) BIND(C,name='F_CDF53_2D_split_inplace')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: ni, nj, lni
     real(C_FLOAT), dimension(*), intent(INOUT) :: x
 end subroutine F_CDF53_2D_split_inplace
 end interface
 interface
 subroutine I_CDF53_2D_split_inplace_n(x, ni, lni, nj, levels) BIND(C,name='I_CDF53_2D_split_inplace_n')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: ni, nj, lni, levels
     real(C_FLOAT), dimension(*), intent(INOUT) :: x
 end subroutine I_CDF53_2D_split_inplace_n
 end interface
 interface
 subroutine F_CDF53_2D_split_inplace_n(x, ni, lni, nj, levels) BIND(C,name='F_CDF53_2D_split_inplace_n')
     import :: C_FLOAT, C_INT
     integer, intent(IN), value :: ni, nj, lni, levels
     real(C_FLOAT), dimension(*), intent(INOUT) :: x
 end subroutine F_CDF53_2D_split_inplace_n
 end interface
#else
void F_CDF53_1D_split_N_even(float *x, float *e, float *o, int n);
void F_CDF53_1D_inplace_N_even(float *x, int n);
void F_CDF53_1D_split_N_odd(float *x, float *e, float *o, int n);
void F_CDF53_1D_inplace_N_odd(float *x, int n);
void F_CDF53_1D_split(float *x, float *e, float *o, int n);
void F_CDF53_1D_inplace(float *x, int n);
void F_CDF53_1D_split_inplace(float *x, int n);
void I_CDF53_1D_split_N_even(float *x, float *e, float *o, int n);
void I_CDF53_1D_inplace_N_even(float *x, int n);
void I_CDF53_1D_split_N_odd(float *x, float *e, float *o, int n);
void I_CDF53_1D_inplace_N_odd(float *x, int n);
void I_CDF53_1D_split(float *x, float *e, float *o, int n);
void I_CDF53_1D_split_inplace(float *x, int n);
void I_CDF53_1D_inplace(float *x, int n);
void I_CDF53_2D_split_inplace(float *x, int ni, int lni, int nj);
void F_CDF53_2D_split_inplace(float *x, int ni, int lni, int nj);
void I_CDF53_2D_split_inplace_n(float *x, int ni, int lni, int nj, int levels);
void F_CDF53_2D_split_inplace_n(float *x, int ni, int lni, int nj, int levels);
#endif
