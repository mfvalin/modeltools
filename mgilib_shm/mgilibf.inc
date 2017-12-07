interface
function mgi_init(name) result(status)
  character(len=*), intent(IN) :: name
  integer :: status
end function mgi_init
function mgi_open(channel, mode) result(status)
  integer, intent(IN) :: channel
  character(len=*), intent(IN) :: mode
  integer :: status
end function mgi_open
function mgi_clos(channel) result(status)
  integer, intent(IN) :: channel
  integer :: status
end function mgi_clos
subroutine mgi_term()
end subroutine mgi_term
end interface
integer, external :: mgi_read , mgi_write
