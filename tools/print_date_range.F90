program print_date_range
  use ISO_C_BINDING
  implicit none
  integer, external :: newdate
  integer :: yyyymmdd, hhmmss, stamp, p1, p2, stamp1, stamp2, p3, p4, stamp3
  integer, dimension(2) :: printable1, printable2, printable3
  real *8 :: delta, diff
  integer :: status, argcount
  character(len=128) :: date1, date2, interval, name, options
  character(len=32) :: arg1, arg2
  character(len=4096) :: oldpath, newpath, dirpath
  character(len=4096) :: nest_rept, nest_exp, anal
  character(len=1), dimension(1) :: template
  integer(C_INT) :: mode

  interface
    function f_mkdir(path,mode) result(status) bind(C,name='mkdir')
      import :: C_CHAR, C_INT
      character(C_CHAR), intent(IN) :: path
      integer(C_INT), intent(IN), value :: mode
      integer(C_INT) :: status
    end function f_mkdir
    function f_link(oldpath,newpath) result(status) bind(C,name='link')
      import :: C_CHAR, C_INT
      character(C_CHAR), intent(IN) :: oldpath, newpath
      integer(C_INT) :: status
    end function f_link
    function f_symlink(oldpath,newpath) result(status) bind(C,name='symlink')
      import :: C_CHAR, C_INT
      character(C_CHAR), intent(IN) :: oldpath, newpath
      integer(C_INT) :: status
    end function f_symlink
  end interface

  CALL get_command_argument(0, name)

  argcount = command_argument_count()
  if(argcount < 6 .or. argcount > 7) goto 777
  mode = o'0777'

  CALL get_command_argument(1, date1)
! newpath = './soft_link'
! status = f_symlink( transfer(trim(oldpath)//achar(0),template(1)), transfer(trim(newpath)//achar(0),template(1)) )
! print *,'status =',status
! stop
  date1 = trim(date1)//'000000'
  read(date1,11,err=777)printable1(1),printable1(2)  ! start date in YYYYMMDD format
11 format(I8,1x,I6)

  CALL get_command_argument(2, date2)
  date2 = trim(date2)//'000000'
  read(date2,11,err=777)printable2(1),printable2(2)  ! end date in YYYYMMDD format

  CALL get_command_argument(3, interval)
  read(interval,*,err=777)delta                      ! interval in seconds
  delta = delta/3600.0                               ! interval in hours

  CALL get_command_argument(4, anal)
  CALL get_command_argument(5, nest_rept)
  CALL get_command_argument(6, nest_exp)
!   nest_rept = 'CLIMAT_nest_rept'
!   nest_exp  = 'CLIMAT_nest_exp'
!   anal = 'GEM_anal'

  if(argcount == 7) then                             ! year=option is present
    CALL get_command_argument(4, options)
    call NewDate_Options(trim(options),'set')        ! set calendar option
    write(0,*),'INFO: using calendqar option '//trim(options)
  endif

  open(unit=11,file='content',form='FORMATTED')
  write(11,'(A/A)')'1','GEM_input_file_0001'
  close(unit=11)

  status = newdate(stamp1,printable1(1),printable1(2)*100,3) ! stamp for start date
  status = newdate(stamp2,printable2(1),printable2(2)*100,3) ! stamp for end date
  call difdatr(stamp2,stamp1,diff)                           ! end - start in hours

  do while(diff >= 0)
    status = newdate(stamp1,p1,p2,-3)                 ! convert to printable
    p3 = p1
    if(p2 == 0) then                                  ! hhmmss = 0, get previus day
      call incdatr(stamp3,stamp1,-24.0_8)
      status = newdate(stamp3,p3,p4,-3)
    endif
!     print 102,p1,'.',p2/100,p3                        ! print it
    write(arg1,'(I8.8,A1,I6.6)')p1,'.',p2/100
    write(arg2,'(I8.8)')p3
    print *,trim(arg1),' ',trim(arg2)
    write(dirpath,'(A)')'VALID_'//trim(arg1)
    status = f_mkdir( transfer(trim(dirpath)//achar(0),template(1)), mode )
    oldpath = 'content'
    newpath = trim(dirpath)//'/content'
    status = f_link( transfer(trim(oldpath)//achar(0),template(1)), transfer(trim(newpath)//achar(0),template(1)) )
    oldpath = trim(nest_rept) // '/' // trim(nest_exp) // '_' // arg2(1:6) // '/' // trim(nest_exp) // '_' // arg2(1:8)
    newpath = 'VALID_' // trim(arg1) // '/GEM_input_file_0001'
    status = f_symlink( transfer(trim(oldpath)//achar(0),template(1)), transfer(trim(newpath)//achar(0),template(1)) )
    call incdatr(stamp,stamp1,delta)                  ! increment
    stamp1 = stamp
    call difdatr(stamp2,stamp1,diff)                  ! end - next date
  enddo
  stop
101 format(3X,I3,3x,i10,3x,i10,3x,I8,3x,I6)
102 format(3x,I8.8,A,I6.6,3x,i8.8)
777 continue
  write(0,*),'USAGE: '//trim(name)//' start_date end_date interval [year=gregorian|360_day|365_day]'
  write(0,*),'       start, end : YYYYMMDD.HHMMSS'
  write(0,*),'       interval in seconds'
  stop
end program
