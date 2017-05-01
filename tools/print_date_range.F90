program print_date_range
  implicit none
  integer, external :: newdate
  integer :: yyyymmdd, hhmmss, stamp, p1, p2, stamp1, stamp2
  integer, dimension(2) :: printable1, printable2
  real *8 :: delta, diff
  integer :: status, argcount
  character(len=128) :: date1, date2, interval, name, options

  CALL get_command_argument(0, name)

  argcount = command_argument_count()
  if(argcount < 3 .or. argcount > 4) goto 777

  CALL get_command_argument(1, date1)
  read(date1,11,err=777)printable1(1),printable1(2)  ! start date in YYYYMMDD format
11 format(I8,1x,I6)

  CALL get_command_argument(2, date2)
  read(date2,11,err=777)printable2(1),printable2(2)  ! end date in YYYYMMDD format

  CALL get_command_argument(3, interval)
  read(interval,*,err=777)delta                      ! interval in seconds
  delta = delta/3600.0                               ! interval in hours

  if(argcount == 4) then                             ! year=option is present
    CALL get_command_argument(4, options)
    call NewDate_Options(trim(options),'set')        ! set calendar option
    write(0,*),'INFO: using calendqar option '//trim(options)
  endif

  status = newdate(stamp1,printable1(1),printable1(2)*100,3) ! stamp for start date
  status = newdate(stamp2,printable2(1),printable2(2)*100,3) ! stamp for end date
  call difdatr(stamp2,stamp1,diff)                           ! end - start in hours

  do while(diff >= 0)
    status = newdate(stamp1,p1,p2,-3)                 ! convert to printable
    print 102,p1,'.',p2/100                           ! print it
    call incdatr(stamp,stamp1,delta)                  ! increment
    stamp1 = stamp
    call difdatr(stamp2,stamp1,diff)                  ! end - next date
  enddo
  stop
101 format(3X,I3,3x,i10,3x,i10,3x,I8,3x,I6)
102 format(3x,I8.8,A,I6.6)
777 continue
  write(0,*),'USAGE: '//trim(name)//' start_date end_date interval [year=gregorian|360_day|365_day]'
  write(0,*),'       start, end : YYYYMMDD.HHMMSS'
  write(0,*),'       interval in seconds'
  stop
end program