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

program f_mgi_test
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

  integer, parameter :: MAXVAL=100
  integer status
  integer, external :: mgi_open, mgi_init, mgi_clos, mgi_term, mgi_read, mgi_write
  character(len=1) :: what, testmode
  character(len=1024) :: string, testfile, string2
  real :: start, end, delta
  integer, dimension(MAXVAL) :: is, is2
  real, dimension(MAXVAL) :: fs, fs2
  real*8, dimension(MAXVAL) :: ds, ds2
  integer :: iostat, nval, i, shmemid, nargs, channel
  character(len=128)channel_name

  nargs = command_argument_count()
  if(nargs >= 1) then
!    call getarg(1,string)
    call get_command_argument(1, string)
    channel_name=trim(string)
    print *,'shared memory channel = ',trim(channel_name)
  endif
  if(nargs >= 2) then
    call getarg(2,string)
    testmode=string(1:1)
    print *,' test mode = ','"'//testmode//'"'
  endif
  testfile='mgi_test.txt'
  if(nargs >= 3) then
    call getarg(3,testfile)
  endif
  print *,' scenario file: '//'"'//trim(testfile)//'"'
  if(nargs >= 4) then
  endif

  channel = mgi_init(trim(channel_name))
  print *,'channel=',channel
  call sleep_a_bit(1)
  if(trim(testfile) == 'none') channel=channel+1000    ! force open
  if(testmode=='R') status = mgi_open(channel,'R')
  if(testmode=='W') status = mgi_open(channel,'W')
  print *,'status=',status
  if(trim(testfile) == 'none') goto 777
  if(status < 0) then
    call print_mgi_error(status)
    print *,'ERROR: cannot open channel ',trim(channel_name)
    stop
  endif
  call sleep_a_bit(1)

  print *,'===== List of error codes ====='
  do i = -1,-11,-1
    call print_mgi_error(i)
  enddo
  open(unit=10,file=trim(testfile),form='FORMATTED')
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
        print *,is(1:nval)
        if(testmode == 'R') then
          status = mgi_read(channel,is2,nval,what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(.not. all(is(1:nval)==is2(1:nval))) print *,'ERROR: did not read what was expected (I)'
          if( all(is(1:nval)==is2(1:nval))) print *,'INFO: read what was expected (I)'
        else
          status = mgi_write(channel,is,nval,what)
          if(status /= 0) print *,'ERROR: write error, status=',status
        endif
      case('R')
        do i=1,nval
          fs(i)=start+(i-1)*delta
        enddo
        print *,fs(1:nval)
        if(testmode == 'R') then
          status = mgi_read(channel,fs2,nval,what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(.not. all(fs(1:nval)==fs2(1:nval))) print *,'ERROR: did not read what was expected (R)'
          if( all(fs(1:nval)==fs2(1:nval))) print *,'INFO: read what was expected (R)'
        else
          status = mgi_write(channel,fs,nval,what)
          if(status /= 0) print *,'ERROR: write error, status=',status
        endif
      case('D')
        do i=1,nval
          ds(i)=start+(i-1)*delta
        enddo
        print *,ds(1:nval)
        if(testmode == 'R') then
          status = mgi_read(channel,ds2,nval,what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(.not. all(ds(1:nval)==ds2(1:nval))) print *,'ERROR: did not read what was expected (D)'
          if( all(ds(1:nval)==ds2(1:nval))) print *,'INFO: read what was expected (D)'
        else
          status = mgi_write(channel,ds,nval,what)
          if(status /= 0) print *,'ERROR: write error, status=',status
        endif
      end select
    else
      read(10,*)what,string
      print *,'"'//what//'",','"'//trim(string)//'"'
      if(testmode == 'R') then
          status = mgi_read(channel,string2,len(trim(string)),what)
          if(status <= 0) print *,'ERROR: read error, status=',status
          if(trim(string) == trim(string2)) print *,'INFO: read what was expected (C)'
          if(trim(string) /= trim(string2)) print *,'ERROR: expected "'//trim(string)//'"'//' READ back : "'//trim(string2)//'"'
      else
          status = mgi_write(channel,trim(string),len(trim(string)),what)
          if(status /= 0) print *,'ERROR: write error, status=',status
      endif
    endif
    read(10,*,iostat=iostat)what
  enddo
  close(unit=10)
777 continue   ! the end
  call sleep_a_bit(1)
  status = mgi_clos(channel)
  call sleep_a_bit(1)
  status = mgi_term()
  stop
end