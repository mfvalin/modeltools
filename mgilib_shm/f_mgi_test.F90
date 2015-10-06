!*
!* Copyright (C) 2014       ESCER center, UQAM
!*
!* This code is free software; you can redistribute it and/or
!* modify it under the terms of the GNU Lesser General Public
!* License as published by the Free Software Foundation,
!* version 2.1 of the License.
!*
!* This code is distributed in the hope that it will be useful,
!* but WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* Lesser General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with this code; if not, write to the
!* Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!* Boston, MA 02111-1307, USA.
!*
program model_main
  implicit none

#if defined(MODEL)

  include 'mpif.h'
  integer :: ierr, me, world
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,me,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,world,ierr)
  print *,'I am PE',me+1,' of',world
  if(me==0) then
#if defined(ATM)
    call atmospheric_model
#endif
#if defined(OCE)
    call ocean_model
#endif
#if ! defined(ATM) && ! defined(OCE)
    print *,'base test with mpi'
    call main_mgi_test
#endif
    call mgi_perf_print
  endif
  call mpi_finalize(ierr)

#else

  call main_mgi_test
  call mgi_perf_print

#endif
  stop
end program model_main

subroutine atmospheric_model
  implicit none
  print *,'atmospheric_model'
!               IN         OUT
  call model('oce2atm','atm2oce')
end subroutine atmospheric_model

subroutine ocean_model
  implicit none
  print *,'ocean_model'
!               IN         OUT
  call model('atm2oce','oce2atm')
end subroutine ocean_model

subroutine main_mgi_test
! usage: f_mgi_test shared_memory_id R|W scenario_file
  use ISO_C_BINDING
  implicit none
!
  interface
    subroutine sleep_a_bit(duration) bind(C,name='sleep')
    use ISO_C_BINDING
    integer(C_INT), intent(IN), value :: duration
    integer(C_INT) :: status
    end subroutine sleep_a_bit
  end interface
!
  integer status, status_r, status_w
  integer, external :: mgi_open, mgi_init, mgi_clos, mgi_term
  character(len=1) :: testmode_r, testmode_w
  character(len=1024) :: string, testfile
  integer :: iostat, nargs, channel_r, channel_w, i, minargs
  character(len=128)channel_name, channel_name_r, channel_name_W

  nargs = command_argument_count()

  if(nargs >= 1) then
    call getarg(1,string)
    testmode_r = ' '
    testmode_W = ' '
    if(string(1:1) == 'R' .or. string(1:1) == 'r') testmode_r = 'R'
    if(string(1:1) == 'W' .or. string(1:1) == 'w') testmode_w = 'W'
    if(string(2:2) == 'R' .or. string(2:2) == 'r') testmode_r = 'R'
    if(string(2:2) == 'W' .or. string(2:2) == 'w') testmode_w = 'W'
    print *,'INFO: test mode = ','"'//testmode_r//testmode_w//'"'
  else
    print *,'ERROR: test mode (R/W/RW) is mandatory'
    return
  endif
  if(testmode_r .ne. 'R' .and. testmode_w .ne. 'W') then
    print *,'ERROR: one of R or W must be specified for test mode'
    return
  endif
  minargs = 1
  if(testmode_r=='R') minargs = minargs + 1
  if(testmode_w=='W') minargs = minargs + 1

  if(nargs >= minargs) then
    call get_command_argument(2, string)
    if(testmode_r=='R') then
      channel_name_r = trim(string)
    else
      channel_name_w = trim(string)
    endif
    if(minargs == 3) then
      call get_command_argument(3, string)
      channel_name_w=trim(string)
    endif
  else
    print *,'ERROR: insufficient number of channel names'
    return
  endif

  if(testmode_r=='R') print *,'INFO: reading from channel = ',trim(channel_name_r)
  if(testmode_w=='W') print *,'INFO: writing into channel = ',trim(channel_name_w)

  testfile='mgi_test.txt'
  if(nargs >= minargs+1) then
    call getarg(minargs+1,testfile)
  endif
  print *,'INFO: scenario file: '//'"'//trim(testfile)//'"'
  if(nargs > minargs+1) then
    print *,'ERROR: too many arguments'
    return
  endif

  call sleep_a_bit(1)
  call mgi_perf_on
!stop
  channel_r = -1
  if(testmode_r=='R') then
    channel_r = mgi_init(trim(channel_name_r))
    print *,'channel_r=',channel_r
    if(trim(testfile) == 'none') channel_r=channel_r+1000    ! force open
    status_r = mgi_open(channel_r,'R')
    print *,'status_r=',status_r
  endif

  channel_w = -1
  if(testmode_w=='W') then
    channel_w = mgi_init(trim(channel_name_w))
    print *,'channel_w=',channel_w
    if(trim(testfile) == 'none') channel_w=channel_w+1000    ! force open
    status_w = mgi_open(channel_w,'W')
    print *,'status_w=',status_w
  endif

  if(trim(testfile) == 'none') goto 777
  if(status_r < 0) then
    call print_mgi_error(status_r)
    print *,'ERROR: cannot open channel ',trim(channel_name_r)
    stop
  endif
  if(status_w < 0) then
    call print_mgi_error(status_r)
    print *,'ERROR: cannot open channel ',trim(channel_name_w)
    stop
  endif
  call sleep_a_bit(1)

  print *,'===== List of error codes ====='
  do i = -1,-11,-1
    call print_mgi_error(i)
  enddo

  open(unit=10,file=trim(testfile),form='FORMATTED',iostat=iostat)  ! open scenario file

  if(iostat == 0) call mgi_test_body(channel_r,channel_w)

777 continue   ! the end
  call sleep_a_bit(1)
  if(testmode_r=='R') status_r = mgi_clos(channel_r)
  if(testmode_w=='W') status_w = mgi_clos(channel_w)
  call sleep_a_bit(1)
  status = mgi_term()
  return
end subroutine main_mgi_test

subroutine model(channel_name_r,channel_name_w)
! usage: f_mgi_test shared_memory_id R|W scenario_file
  use ISO_C_BINDING
  implicit none
  interface
    subroutine sleep_a_bit(duration) bind(C,name='sleep')
    use ISO_C_BINDING
    integer(C_INT), intent(IN), value :: duration
    integer(C_INT) :: status
    end subroutine sleep_a_bit
  end interface

  integer :: status, status_r, status_w
  integer, external :: mgi_open, mgi_init, mgi_clos, mgi_term
  integer :: iostat, channel_r, channel_w, i
  character(len=*), intent(IN) :: channel_name_r
  character(len=*), intent(IN) :: channel_name_w

  call sleep_a_bit(1)
  call mgi_perf_on
  print *,'reading from '//channel_name_r//', writing to '//channel_name_w
!  return

  channel_r = mgi_init(trim(channel_name_r))
  print *,'channel_r=',channel_r
  status_r = mgi_open(channel_r,'R')
  if(status_r < 0) then
    call print_mgi_error(status_r)
    print *,'ERROR: cannot open channel ',trim(channel_name_r)
    return
  else
    print *,'status_r=',status_r
  endif

  channel_w = mgi_init(trim(channel_name_w))
  print *,'channel_w=',channel_w
  status_w = mgi_open(channel_w,'W')
  if(status_w < 0) then
    call print_mgi_error(status_r)
    print *,'ERROR: cannot open channel ',trim(channel_name_w)
    return
  else
    print *,'status_w=',status_w
  endif


  print *,'===== List of error codes ====='
  do i = -1,-11,-1
    call print_mgi_error(i)
  enddo

  open(unit=10,file='oce_acm.data',form='FORMATTED',iostat=iostat)  ! open scenario file

  if(iostat == 0) then
    call mgi_test_body(channel_r,channel_w)
  else
    print *,'ERROR: cannot open file oce_acm.data'
  endif

777 continue   ! the end
  call sleep_a_bit(1)
  status_r = mgi_clos(channel_r)
  status_w = mgi_clos(channel_w)
  call sleep_a_bit(1)
  status = mgi_term()
  return
end subroutine model

subroutine mgi_test_body(channel_r,channel_w)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: channel_r,channel_w

  interface
    subroutine sleep_a_bit(duration) bind(C,name='sleep')
    use ISO_C_BINDING
    integer(C_INT), intent(IN), value :: duration
    integer(C_INT) :: status
    end subroutine sleep_a_bit
  end interface

  integer, external :: mgi_read, mgi_write, mgi_read_c, mgi_write_c
  integer, parameter :: MAXVAL=1000000
  integer, dimension(MAXVAL) :: is, is2
  real, dimension(MAXVAL) :: fs, fs2
  real*8, dimension(MAXVAL) :: ds, ds2
  character(len=1) :: what
  real :: start, end, delta
  integer :: i, nval, status, iostat
  character(len=1024) :: string, string2

  read(10,*,iostat=iostat)what
  do while(iostat == 0)
    backspace(10)
    call sleep_a_bit(1)
    if(what .ne. 'C') then
      read(10,*)what,start,end,delta
      print *,'"'//what//'",',start,end,delta
      nval = min(MAXVAL,nint(1.0+(end-start)/delta))
      select case(what)
      case('I')
        do i=1,nval
          is(i)=nint(start+(i-1)*delta)
        enddo
        if(nval <= 10) print *,is(1:nval)
        if(channel_w .ne. -1) then
          status = mgi_write(channel_w,is,nval,what)
          if(status /= 0) print *,'ERROR: write error, status=',status
        endif
        if(channel_r .ne. -1) then
          status = mgi_read(channel_r,is2,nval,what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(.not. all(is(1:nval)==is2(1:nval))) print *,'ERROR: did not read what was expected (I)'
          if( all(is(1:nval)==is2(1:nval))) print *,'INFO: read what was expected (I)'
        endif
      case('R')
        do i=1,nval
          fs(i)=start+(i-1)*delta
        enddo
        if(nval <= 10) print *,fs(1:nval)
        if(channel_w .ne. -1) then
          status = mgi_write(channel_w,fs,nval,what)
          if(status /= 0) print *,'ERROR: write error, status=',status
        endif
        if(channel_r .ne. -1) then
          status = mgi_read(channel_r,fs2,nval,what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(.not. all(fs(1:nval)==fs2(1:nval))) print *,'ERROR: did not read what was expected (R)'
          if( all(fs(1:nval)==fs2(1:nval))) print *,'INFO: read what was expected (R)'
        endif
      case('D')
        do i=1,nval
          ds(i)=start+(i-1)*delta
        enddo
        if(nval <= 10) print *,ds(1:nval)
        if(channel_w .ne. -1) then
          status = mgi_write(channel_w,ds,nval,what)
          if(status /= 0) print *,'ERROR: write error, status=',status
        endif
        if(channel_r .ne. -1) then
          status = mgi_read(channel_r,ds2,nval,what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(.not. all(ds(1:nval)==ds2(1:nval))) print *,'ERROR: did not read what was expected (D)'
          if( all(ds(1:nval)==ds2(1:nval))) print *,'INFO: read what was expected (D)'
        endif
      end select
    else   !  what .ne. 'C'
      read(10,*)what,string
      print *,'"'//what//'",','"'//trim(string)//'"'
        if(channel_w .ne. -1) then
          status = mgi_write_c(channel_w,trim(string),len(trim(string)),what)
          if(status /= 0) print *,'ERROR: write error, status=',status
      endif
      if(channel_r .ne. -1) then
          status = mgi_read_c(channel_r,string2,len(trim(string)),what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(trim(string) == trim(string2)) print *,'INFO: read what was expected (C)'
          if(trim(string) /= trim(string2)) print *,'ERROR: expected "'//trim(string)//'"'//' READ back : "'//trim(string2)//'"'
        endif
    endif
    read(10,*,iostat=iostat)what
  enddo
  close(unit=10)

  return
end subroutine mgi_test_body