subroutine flush_listing_stream(iun)
! flush output from Fortran unit iun (normally unit 6 connected to stdout)
! call C routine to flush stdout and stderr
  implicit none
  integer, intent(IN) :: iun
  call flush(iun)
  call flush_listing_c_stream
  return
end