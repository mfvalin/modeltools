module c_mgi_interfaces
  use ISO_C_BINDING
  interface
    function c_mgi_init(name) result(status) BIND(C,name='MGI_Init')
    import :: C_INT
      character(len=1), dimension(*), intent(IN) :: name
      integer(C_INT) :: status
    end function c_mgi_init
  end interface
end module c_mgi_interfaces
function mgi_init(name) result(status)
  use c_mgi_interfaces
  implicit none
  character(len=*), intent(IN) :: name
  integer :: status

  character(len=1), dimension(128) :: temp
  integer :: l

  l = len(trim(name))
  temp(1:l+1) = transfer( trim(name)//achar(0) , temp)

  status = c_mgi_init(temp)

end function mgi_init
