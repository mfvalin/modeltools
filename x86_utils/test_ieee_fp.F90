program tst_iee_fp
  use ISO_C_BINDING
  include 'fpu_exceptions.inc'
  real :: a, b, c
  integer old
  interface
    integer(C_LONG) function f_fork() BIND(C,name='fork')
    import :: C_LONG
    end function f_fork
  end interface

!   old = disable_fpu_exceptions(FPU_UNDERFLOW+FPU_FTZ)
  old = ieee_fp_trap_on()
  old = ieee_denormal_on()
  a = 1.0E-22
  b = 1.0E-22
  c = a * b
  print 100,"c = ",c,c
  old = ieee_denormal_off()
  a = 1.0E-22
  b = 1.0E-22
  c = a * b
  print 100,"c = ",c,c
  old = ieee_denormal_on()
  a = 1.0E-22
  b = 1.0E-22
  c = a * b
  print 100,"c = ",c,c
100 format(A,G15.6,2X,Z8.8)
  print *,"end of denormal test"

  a = 1.0
  b = 0.0
  old = ieee_fp_trap_off()
  c = a/b
  print *,"after first zero divide (SUCCESS, this message must appear), c=",c

  if(f_fork() == 0) then
    print *,"restoring FP interrupts zero divide in child"
    old = ieee_fp_trap_on()
    c = a/b
    print 100,"c = ",c,c
    print *,"after second zero divide (ERROR, this must not appear)"
  endif
  old = ieee_denormal_off()
  a = 1.0E-20
  b = 1.0E-20
  c = a * b
  print 100,"denormals OFF, c = ",c,c

  if(f_fork() == 0) then
    print *,"restoring FP interrupts in child"
    old = ieee_denormal_on()
    old = ieee_fp_trap_on()
    a = 1.0E-20
    b = 1.0E-20
    c = a * b
    print 100,"after restoring FP interrupts (ERROR, this must not appear) c = ",c,c
  endif
  print *,"end of test"
  stop
end program
