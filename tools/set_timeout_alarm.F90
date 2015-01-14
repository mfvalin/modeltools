function set_timeout_alarm(seconds) result(seconds_since)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: seconds
  integer :: seconds_since

  interface
  function c_alarm(seconds) result(seconds_since) BIND(C,name='alarm')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), value :: seconds
    integer(C_INT) :: seconds_since
  end function c_alarm
  end interface

  seconds_since = c_alarm(seconds)
  print *,'alarm set to ',seconds,' seconds'
  return
end function set_timeout_alarm

#ifdef TEST
program test
  implicit none
  integer, external :: set_timeout_alarm
  integer :: seconds_since
  integer :: cond

  seconds_since = set_timeout_alarm(10)
  print *,'seconds_since =',seconds_since
  read(5,*) cond
  seconds_since = set_timeout_alarm(10)
  read(5,*) cond
  print *,'seconds_since =',seconds_since
  seconds_since = set_timeout_alarm(10)
  read(5,*) cond
  print *,'seconds_since =',seconds_since
  seconds_since = set_timeout_alarm(10)
  read(5,*) cond
  print *,'seconds_since =',seconds_since
  
end program test
#endif
