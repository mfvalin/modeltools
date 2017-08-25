program print_date_range
  use ISO_C_BINDING
  implicit none
  integer, external :: newdate
  integer :: stamp, p1, p2, stamp1, stamp2, p3, p4, stamp3
  integer, dimension(2) :: printable1, printable2, printable3
  real *8 :: delta, diff
  integer :: status
  character(len=128) :: date1, date2, interval, name, sym
  character(len=32) :: arg1, arg2
  character(len=4096) :: oldpath, newpath, dirpath, option
  character(C_CHAR), dimension(4096) :: oldp, newp, dirp
  character(len=4096) :: nest_rept, set_name, anal
  integer(C_INT) :: mode
  logical :: use_anal
  integer :: cur_arg, nargs, arg_len, ntimes

  interface
    function f_mkdir(path,mode) result(status) bind(C,name='mkdir')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: path
      integer(C_INT), intent(IN), value :: mode
      integer(C_INT) :: status
    end function f_mkdir
    function f_unlink(path) result(status) bind(C,name='unlink')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: path
      integer(C_INT) :: status
    end function f_unlink
    function f_link(oldpath,newpath) result(status) bind(C,name='link')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: oldpath, newpath
      integer(C_INT) :: status
    end function f_link
    function f_symlink(oldpath,newpath) result(status) bind(C,name='symlink')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: oldpath, newpath
      integer(C_INT) :: status
    end function f_symlink
  end interface

  CALL get_command_argument(0, name)

  mode = o'0777'

  cur_arg = 1
  nargs = command_argument_count()
  printable1 = -1
  printable2 = -1
  printable3 = -1
  delta = -1
  anal = ''
  set_name = ''
  nest_rept = ''
  do while(cur_arg <= nargs)      ! process command line options
    call get_command_argument(cur_arg,option,arg_len,status)
    if(option(1:13)      == '--start_date=' ) then       ! YYYYMMDD.HHMMSS
      date1 = trim(option(14:4096))//'000000'
      read(date1,11,err=777)printable1(1),printable1(2)  ! start date in YYYYMMDD format
    else if(option(1:11) == '--end_date=' ) then         ! YYYYMMDD.HHMMSS
      date2 = trim(option(12:4096))//'000000'
      read(date2,11,err=777)printable2(1),printable2(2)  ! end date in YYYYMMDD format
    else if(option(1:9)  == '--nhours=' ) then           ! hours
      interval = option(10:4096)
      read(interval,*,err=777)delta                      ! interval in hours
    else if(option(1:11) == '--nseconds=' ) then         ! seconds
      interval = option(10:4096)
      read(interval,*,err=777)delta                      ! interval in seconds
      delta = delta/3600.0                               ! interval in hours
    else if(option(1:12) == '--start_sym=' ) then        ! YYYYMMDD.HHMMSS
      sym = trim(option(13:4096))//'000000'
      read(sym,11,err=777)printable3(1),printable3(2)    ! start of simulation in YYYYMMDD format
    else if(option(1:13) == '--start_anal=' ) then       ! initial analysis (only necessary if start_sym == start_date)
      anal = option(14:4096)
    else if(option(1:13) == '--pilot_data=' ) then       ! directory for boundary conditions
      nest_rept = option(14:4096)
    else if(option(1:11) == '--set_name=' ) then         ! experiment name
      set_name = option(12:4096)
    else if(option(1:7)  == '--year=' ) then             ! calendar option (optional)
      call NewDate_Options(trim(option(8:4096)),'set')       ! set calendar option
      write(0,*),'INFO: using calendar option '//trim(option)
    else if(option(1:9)  == '--version' ) then              ! version option
      write(0,*),'version = 1.0.0 2017/08/25'
      goto 777
    else if(option(1:6)  == '--help' ) then              ! help option
      goto 777
    else if(option(1:2)  == '-h' ) then                  ! help option
      goto 777
    else 
      print *,"ERROR: unrecognized option '"//trim(option)//"'"
    endif
    cur_arg = cur_arg + 1
  enddo
  use_anal = (printable3(1) == printable1(1)) .and. (printable3(2) == printable1(2))
  if(printable1(1) == -1 .or. printable2(1) == -1) then
    write(0,*),'ERROR: bad start date'
    goto 777
  endif
  if(printable2(1) == -1 .or. printable2(1) == -1) then
    write(0,*),'ERROR: bad end date'
    goto 777
  endif
  if(use_anal .and. trim(anal) == '' ) then
     write(0,*),'ERROR: initial conditions missing'
    goto 777
  endif
  if(trim(nest_rept) == '' ) then
    write(0,*),'ERROR: boundary conditions directory missing'
    goto 777
  endif
  if(trim(set_name) == '' ) then
    write(0,*),'ERROR: dataset name missing'
    goto 777
  endif
  write(0,*),'INFO: from ',trim(date1),' to ',trim(date2),' every',delta,' hours'
  write(0,*),"INFO: nesting data in '"//trim(nest_rept)//"'"
  if(use_anal) write(0,*),"INFO: initial conditions '"//trim(anal)//"'"

  open(unit=11,file='content',form='FORMATTED')
  write(11,'(A/A)')'1','GEM_input_file_0001'
  close(unit=11)

  status = newdate(stamp1,printable1(1),printable1(2)*100,3) ! stamp for start date
  status = newdate(stamp2,printable2(1),printable2(2)*100,3) ! stamp for end date
  call difdatr(stamp2,stamp1,diff)                           ! end - start in hours

  ntimes = 0
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
!     print *,trim(arg1),' ',trim(arg2)
    write(dirpath,'(A)')'VALID_'//trim(arg1)
    dirp = transfer(trim(dirpath)//achar(0),dirp)
    status = f_mkdir( dirp, mode )
    oldpath = 'content'
    oldp = transfer(trim(oldpath)//achar(0),oldp)
    newpath = trim(dirpath)//'/content'
    newp = transfer(trim(newpath)//achar(0),newp)
    if(use_anal) status = f_unlink( newp )
    status = f_link( oldp, newp )
    oldpath = trim(nest_rept) // '/' // trim(set_name) // '_' // arg2(1:6) // '/' // trim(set_name) // '_' // arg2(1:8)
    if(use_anal) oldpath = trim(anal)
    oldp = transfer(trim(oldpath)//achar(0),oldp)
    newpath = 'VALID_' // trim(arg1) // '/GEM_input_file_0001'
    newp = transfer(trim(newpath)//achar(0),newp)
    if(use_anal) status = f_unlink( newp )
    status = f_symlink( oldp, newp )
    call incdatr(stamp,stamp1,delta)                  ! increment
    stamp1 = stamp
    call difdatr(stamp2,stamp1,diff)                  ! end - next date
    use_anal = .false.
    ntimes = ntimes + 1
  enddo
  write(0,*),"INFO: ",ntimes," directory/link sets created"
  stop
11  format(I8,1x,I6)
777 continue
  write(0,*),'USAGE: '//trim(name)//' [-h|--help] --start_date= --end_date= --nhours= --nseconds= [--start_sym=] \'
  write(0,*),'        [--version] [--start_anal=] --pilot_data= --set_name= [--year=gregorian|360_day|365_day]'
  write(0,*),'       start_date : YYYYMMDD.HHMMSS , start end end of this simulation slice'
  write(0,*),'       end_date   : YYYYMMDD.HHMMSS , start end end of this simulation slice'
  write(0,*),'       nseconds   : interval in seconds between boundary condition files'
  write(0,*),'       nhours     : interval in hours between boundary condition files'
  write(0,*),'       start_sym  : YYYYMMDD.HHMMSS , start of entire simulation'
  write(0,*),'       start_anal : initial conditions (only used if start_date == start_sym)'
  write(0,*),'       pilot_data : directory containing the boundary condition files'
  write(0,*),'       set_name   : dataset name'
  write(0,*),'       year=...   : (optional) argument ,  calendar to be used (gregorian by default)'
  write(0,*),'       arguments between [] are optional'
  stop
end program
