module qad_trace_mod
  use ISO_C_BINDING
  integer, save :: curstp = 1
  integer, save :: curtmg = 1
  integer, save :: maxstp = 0
  integer, save :: maxtmg = 0
  integer(C_LONG), save :: t0 = 0
  integer(C_LONG), save, dimension(:,:), pointer :: t => NULL()
interface
  subroutine time_of_day(tv,tz) bind(C,name='gettimeofday')
    import :: C_LONG, C_PTR
    integer(C_LONG), dimension(2), intent(OUT) :: tv
    type(C_PTR), value :: tz
  end subroutine time_of_day
end interface
end module qad_trace_mod

subroutine qad_trace_init(maxi, maxj)
  use qad_trace_mod
  implicit none
  integer, intent(IN) :: maxi          ! max number of timing points
  integer, intent(IN) :: maxj          ! max number of timesteps
  integer(C_LONG), dimension(2) :: tv

  if(associated(t)) then
    print *,'ERROR: timing table already initialized with dimension',maxtmg,' by',maxstp
    return
  endif
  allocate(t(maxi, maxj))
  t = -1
  maxtmg = maxi
  maxstp = maxj
  call time_of_day(tv,C_NULL_PTR)
  t0 = tv(1)      ! base time in seconds that weill be substracted from timings
  return
end subroutine qad_trace_init

subroutine qad_trace_new_step
  use qad_trace_mod
  implicit none
  curstp = min(curstp+1,maxstp)
  curtmg = 1
end subroutine qad_trace_new_step

subroutine qad_trace_point
  use qad_trace_mod
  implicit none
  integer(C_LONG), dimension(2) :: tv

  if(.not. associated(t)) then
    print *,'ERROR: timing table not initialized (qad_trace_init has not been called)'
    return
  endif
  call time_of_day(tv,C_NULL_PTR)
  t(curtmg,curstp) = (( tv(1) - t0)*1000000 + tv(2))
  curtmg = min(curtmg+1,maxtmg)
end subroutine qad_trace_point

subroutine qad_trace_dump(iun,fmt)
  use qad_trace_mod
  implicit none
  integer, intent(IN) :: iun
  character(len=*), intent(IN) :: fmt
  integer :: i

  do i = 1, curstp, 1
    if(fmt(1:1) == '(') then
      write(iun,fmt)i,maxtmg,t(1:maxtmg,i)
    else
      write(iun)    i,maxtmg,t(1:maxtmg,i)
    endif
  enddo
end subroutine qad_trace_dump

subroutine qad_trace_index(tmg, stp)  ! get current timing "point" and "step"
  use qad_trace_mod
  implicit none
  integer, intent(OUT) :: tmg, stp

  tmg = curtmg
  stp = curstp
end subroutine qad_trace_index

#if defined(SELF_TEST)
program test_time
implicit none
integer :: step

call qad_trace_point      ! deliberate error
call qad_trace_init(4,3)
call qad_trace_init(8,5)  ! deliberate error

call qad_trace_point
print *,'hello world 0a', 123.456
call qad_trace_point
call qad_trace_point
print *,'hello world 0b'
call qad_trace_point

do step = 1,2
  call qad_trace_new_step
  call qad_trace_point
  print *,'hello world, step=',step
  call qad_trace_point
  print *,'second print',123
  call qad_trace_point
  call qad_trace_point
enddo

call qad_trace_dump(0,'(I6,(3I8))')
call qad_trace_dump(20,' ')

stop
end
#endif
