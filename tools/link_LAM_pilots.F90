program print_date_range
  use ISO_C_BINDING
  implicit none
#include "clib_interface.cdk"
#define CLIB_OK 1
  integer, external :: newdate
  integer :: stamp, p1, p2, stamp1, stamp2, p3, p4, stamp3
  integer, dimension(2) :: printable1, printable2, printable3
  real *8 :: delta, diff
  integer :: status
  character(len=128) :: date1, date2, interval, name, sym
  character(len=32) :: arg1, arg2
  character(len=4096) :: oldpath, newpath, dirpath, option, oldmonth, month_name
  character(len=1024) :: set_pattern
  character(C_CHAR), dimension(4096) :: oldp, newp, dirp
  character(len=4096) :: nest_rept, set_name, anal
  integer(C_INT) :: mode
  logical :: use_anal, first_in_month
  integer :: cur_arg, nargs, arg_len, ntimes
  integer :: month_is_file = 0
  character(len=128) :: version = 'version 1.0.6a 2017/11/13'
  integer, parameter :: MAXGLOB=2
  character(len=4096), dimension(MAXGLOB) :: globs
  integer :: nglob, arg2_nc
  character(len=16) :: template

  interface
    function f_mkdir(path,mode) result(status) bind(C,name='mkdir')   ! interface to libc mkdir
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: path
      integer(C_INT), intent(IN), value :: mode
      integer(C_INT) :: status
    end function f_mkdir
    function f_unlink(path) result(status) bind(C,name='unlink')   ! interface to libc unlink (rm)
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: path
      integer(C_INT) :: status
    end function f_unlink
    function f_link(oldpath,newpath) result(status) bind(C,name='link')   ! interface to libc link (hard link)
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: oldpath, newpath
      integer(C_INT) :: status
    end function f_link
    function f_symlink(oldpath,newpath) result(status) bind(C,name='symlink')   ! interface to libc symlink (soft link)
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: oldpath, newpath
      integer(C_INT) :: status
    end function f_symlink
  end interface

  CALL get_command_argument(0, name)   ! program name as seen by OS

  mode = o'0777'  ! to be "anded" with user's umask for mkdir
  oldmonth = ' '
  set_pattern = '*'  ! default filename pattern for set name "globbing"

  cur_arg = 1
  nargs = command_argument_count()
  if(nargs == 0) goto 777
  printable1 = -1
  printable2 = -1
  printable3 = -1
  delta = -1
  anal = ''
  set_name = ''
  nest_rept = ''
  arg2_nc = 8
  first_in_month = .true.
  template = 'YYYYMM????????'

  do while(cur_arg <= nargs)      ! process command line options
    call get_command_argument(cur_arg,option,arg_len,status)
    if(option(1:13)      == '--start_date=' ) then       ! YYYYMMDD.HHMMSS or YYYYMMDDHHMMSS
      date1 = trim(option(14:4096))//'000000'
      if(date1(9:9) == '.') then
        read(date1,11,err=777)printable1(1),printable1(2)  ! start date in YYYYMMDD format
      else
        read(date1,12,err=777)printable1(1),printable1(2)  ! start date in YYYYMMDD format
      endif
    else if(option(1:11) == '--end_date=' ) then         ! YYYYMMDD.HHMMSS or YYYYMMDDHHMMSS
      date2 = trim(option(12:4096))//'000000'
      if(date2(9:9) == '.') then
        read(date2,11,err=777)printable2(1),printable2(2)  ! end date in YYYYMMDD format
      else
        read(date2,12,err=777)printable2(1),printable2(2)  ! end date in YYYYMMDD format
      endif
    else if(option(1:9)  == '--nhours=' ) then           ! hours
      interval = option(10:4096)
      read(interval,*,err=777)delta                      ! interval in hours
    else if(option(1:11) == '--nseconds=' ) then         ! seconds
      interval = option(12:4096)
      read(interval,*,err=777)delta                      ! interval in seconds
      delta = delta/3600.0                               ! interval in hours
    else if(option(1:12) == '--start_sym=' ) then        ! YYYYMMDD.HHMMSS or YYYYMMDDHHMMSS
      sym = trim(option(13:4096))//'000000'
      if(sym(9:9) == '.') then
        read(sym,11,err=777)printable3(1),printable3(2)  ! start of simulation in YYYYMMDD format
      else
        read(sym,12,err=777)printable3(1),printable3(2)  ! start of simulation in YYYYMMDD format
      endif
    else if(option(1:13) == '--start_anal=' ) then       ! initial analysis (only necessary if start_sym == start_date)
      anal = option(14:4096)
    else if(option(1:13) == '--pilot_data=' ) then       ! directory for boundary conditions
      nest_rept = option(14:4096)
    else if(option(1:11) == '--set_name=' ) then         ! experiment name
      set_name = option(12:4096)
    else if(option(1:11) == '--set_pattern=' ) then      ! disambiguation pattern for file "globbing", default is '*'
      set_pattern = option(15:4096)
    else if(option(1:7)  == '--year=' ) then             ! calendar option (optional)
      call NewDate_Options(trim(option(3:4096)),'set')       ! set calendar option
      write(0,*),'INFO: using calendar option '//trim(option)
    else if(option(1:9)  == '--version' ) then              ! version option
      write(0,*),version
      stop
    else if(option(1:6)  == '--help' ) then              ! help option
      goto 777
    else if(option(1:2)  == '-h' ) then                  ! help option
      goto 777
    else 
      print *,"ERROR: unrecognized option '"//trim(option)//"'"
      goto 777
    endif
    cur_arg = cur_arg + 1
  enddo
  use_anal = (printable3(1) == printable1(1)) .and. (printable3(2) == printable1(2))
  if(printable1(1) == -1 .or. printable2(1) == -1) then
    write(0,*),'ERROR: missing start/end date(s)'
    goto 777
  endif
  if(printable1(1) == -1 .or. printable1(1) == -1) then
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
  write(0,*),"INFO: boundary conditions data in '"//trim(nest_rept)//"'"
  if(use_anal) write(0,*),"INFO: using initial conditions file '"//trim(anal)//"'"

  open(unit=11,file='content',form='FORMATTED')
  write(11,'(A/A)')'1','GEM_input_file_0001'   ! there will be 1 file in the directory and it will be called GEM_input_file_0001
  close(unit=11)

  status = newdate(stamp1,printable1(1),printable1(2)*100,3) ! stamp for start date
  status = newdate(stamp2,printable2(1),printable2(2)*100,3) ! stamp for end date
  call difdatr(stamp2,stamp1,diff)                           ! end - start in hours

  ntimes = 0
  do while(diff >= 0)                                 ! end date - next date
    status = newdate(stamp1,p1,p2,-3)                 ! convert to printable
    p3 = p1                                           ! YYYYMMDD
    if(p2 == 0 .and. (.not. use_anal)) then           ! hhmmss = 0, use previous day, except if use_anal is true
      call incdatr(stamp3,stamp1,-24.0_8)
      status = newdate(stamp3,p3,p4,-3)
    endif
!     print 102,p1,'.',p2/100,p3                        ! print itname
    write(arg1,'(I8.8,A1,I6.6)')p1,'.',p2/100          ! YYYYMMDD.hhmmss
    write(arg2,'(2I8.8)')p3,p2                         ! YYYYMMDDhhmmss00
!     print *,trim(arg1),' ',trim(arg2)
    write(dirpath,'(A)')'VALID_'//trim(arg1)           ! VALID_YYYYMMDD.hhmmss
    dirp = transfer(trim(dirpath)//achar(0),dirp)      ! C null terminated string from Fortran string
    status = f_mkdir( dirp, mode )                     ! directory containing boundary files for this time interval

    oldpath = 'content'
    oldp = transfer(trim(oldpath)//achar(0),oldp)      ! C null terminated string from Fortran string
    newpath = trim(dirpath)//'/content'                ! content file in directory created above
    newp = transfer(trim(newpath)//achar(0),newp)      ! C null terminated string from Fortran string
    status = f_unlink( newp )                          ! in case file/link named ....../content already exists
    status = f_link( oldp, newp )                      ! hard link to file named 'content' in upper directory

    month_name = trim(nest_rept) // '/' // trim(set_name) // '_' // arg2(1:6)   ! can be a file or a directory
    if(trim(oldmonth) .ne. month_name) then      ! new month name
      oldmonth = month_name
      first_in_month = .true.
      month_is_file = clib_isfile( month_name )  ! it is a file name
!       write(0,*),month_name,month_is_file
      if(month_is_file == 1) then
        write(0,*),'INFO: using monthly boundary file '//trim(oldmonth)
      else
        if(clib_isdir( month_name ) .ne. 1) then ! is it a directory name
          write(0,*),'ERROR: '//trim(oldmonth)//' is neither a directory nor a file, ABORTING'
          stop
        endif
!         write(0,*),'INFO: using monthly boundary files directory '//trim(oldmonth)
      endif
    endif

    if(use_anal) then        ! use initial conditions file instead of boundary conditions file
      oldpath = trim(anal)
    else
      if(month_is_file == 1) then   ! monthly boundary contitions file, may get linked to multiple times
        oldpath = month_name
      else    ! same month, another day
        if(first_in_month) then  ! see if name ends in YYYYMMDDhh, if so use 10 chars from arg2
          arg2_nc = 8
          do while(arg2_nc <= 14)
            oldpath = trim(month_name) // '/' // trim(set_pattern) // arg2(1:arg2_nc)   ! look for 'pattern'YYYYMMDD[hh][mm][ss] file name
!              write(0,*),'DEBUG: trying ' // trim(oldpath)
            status = clib_glob(globs,nglob,trim(oldpath),MAXGLOB)            ! find file name match(es)
            if(status == CLIB_OK .and. nglob == 1) exit                      ! found unique match , exit loop
            arg2_nc = arg2_nc + 2                                            ! try longer match
          enddo
          write(0,*),'INFO: using boundary files pattern '//trim(oldmonth)//'/'//trim(set_pattern)//arg2(1:6)//template(7:arg2_nc)
          if(arg2_nc > 14) then  ! OOPS
            write(0,*),'ERROR: cannot determine file name pattern for ' // trim(month_name)
            stop
          endif
        endif
        oldpath = trim(month_name) // '/' // trim(set_pattern) // arg2(1:arg2_nc) ! look for 'pattern'YYYYMMDD[hh[mmdd]]  ( default pattern is * )
        globs(1) = 'UnknownFile'
!         write(0,*),'INFO: looking for '//trim(oldpath)
        status = clib_glob(globs,nglob,trim(oldpath),MAXGLOB)            ! find file name match(es)
        if(status .ne. CLIB_OK .or. nglob > 1) then                      ! there must be one and only one match
           write(0,*),'ERROR: '//trim(oldpath)//' is ambiguous or does not exist'
           stop
        endif
        oldpath = globs(1)     ! use file name that matches pattern
      endif
    endif

    oldp = transfer(trim(oldpath)//achar(0),oldp)
    newpath = 'VALID_' // trim(arg1) // '/GEM_input_file_0001'
    newp = transfer(trim(newpath)//achar(0),newp)     ! C null terminated string from Fortran string
    status = f_unlink( newp )                         ! in case ....../GEM_input_file_0001 already exists
    status = f_symlink( oldp, newp )                  ! soft link to initial/boundary conditions file

    call incdatr(stamp,stamp1,delta)                  ! increment next date
    stamp1 = stamp
    call difdatr(stamp2,stamp1,diff)                  ! end date - next date

    first_in_month = .false.
    if(use_anal) first_in_month = .true.
    use_anal = .false.                                ! use_anal can only be true for the first time frame
    ntimes = ntimes + 1                               ! counter for time frames
  enddo
  write(0,*),"INFO: ",ntimes," directory/link sets created"
  stop
11  format(I8,1x,I6)
12  format(I8,I6)
777 continue
  write(0,*),'USAGE: '//trim(name)//' [-h|--help] --start_date= --end_date= --nhours= --nseconds= [--start_sym=] \'
  write(0,*),'        [--version] [--start_anal=] --pilot_data= --set_name= [--set_pattern] [--year=gregorian|360_day|365_day]'
  write(0,*),'       '//version
  write(0,*),''
  write(0,*),'       start_date : YYYYMMDD.HHMMSS or YYYYMMDDHHMMSS , start of this simulation slice'
  write(0,*),'       end_date   : YYYYMMDD.HHMMSS or YYYYMMDDHHMMSS , end of this simulation slice'
  write(0,*),'       nseconds   : interval in seconds between boundary condition files'
  write(0,*),'       nhours     : interval in hours between boundary condition files'
  write(0,*),'       start_sym  : YYYYMMDD.HHMMSS or YYYYMMDDHHMMSS , start of entire simulation'
  write(0,*),'       start_anal : initial conditions (only used if start_date == start_sym)'
  write(0,*),'       pilot_data : directory containing the boundary condition files'
  write(0,*),'       set_name   : dataset/experiment name'
  write(0,*),'       set_pattern: disambiguation pattern for file "globbing"'
  write(0,*),'       year=...   : (optional) argument ,  calendar to be used (gregorian by default)'
  write(0,*),'       arguments between [] are optional'
  write(0,*),'       one of nhours/nseconds is necessary'
  write(0,*),'       for date parameters, the trailing 0s in the HHMMSS part can be omitted'
  stop
end program
