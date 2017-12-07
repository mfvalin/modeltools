interface
function mgi_init(name) result(status)         !InTf!
  character(len=*), intent(IN) :: name         !InTf!
  integer :: status                            !InTf!
end function mgi_init                          !InTf!
function mgi_open(channel, mode) result(status) !InTf!
  integer, intent(IN) :: channel                !InTf!
  character(len=*), intent(IN) :: mode          !InTf!
  integer :: status                             !InTf!
end function mgi_open                           !InTf!
function mgi_clos(channel) result(status)       !InTf!
  integer, intent(IN) :: channel                !InTf!
  integer :: status                             !InTf!
end function mgi_clos                           !InTf!
subroutine mgi_term()     !InTf!
end subroutine mgi_term   !InTf!
end interface
integer, external :: mgi_read , mgi_write
