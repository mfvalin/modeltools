module c_mgi_interfaces
  use ISO_C_BINDING
  interface
    function c_mgi_init(name) result(status) BIND(C,name='MGI_Init')
      import :: C_INT
      character(len=1), dimension(*), intent(IN) :: name
      integer(C_INT) :: status
    end function c_mgi_init
    function c_mgi_open(channel, mode) result(status) BIND(C,name='MGI_Open')
      import :: C_INT
      integer(C_INT), intent(IN) :: channel
      character(len=1), intent(IN) :: mode
      integer(C_INT) :: status
    end function c_mgi_open
    function c_mgi_clos(channel) result(status) BIND(C,name='MGI_Clos')
      import :: C_INT
      integer(C_INT), intent(IN) :: channel
      integer(C_INT) :: status
    end function c_mgi_clos
    function c_mgi_read(channel, data, n, typ) result(status) BIND(C,name='MGI_Read')
      import :: C_INT
      integer(C_INT), intent(IN) :: channel, n
      integer(C_INT), dimension(*), intent(OUT) :: data
      character(len=1), intent(IN) :: typ
      integer(C_INT) :: status
    end function c_mgi_read
    function c_mgi_write(channel, data, n, typ) result(status) BIND(C,name='MGI_Write')
      import :: C_INT
      integer(C_INT), intent(IN) :: channel, n
      integer(C_INT), dimension(*), intent(IN) :: data
      character(len=1), intent(IN) :: typ
      integer(C_INT) :: status
    end function c_mgi_write
    subroutine c_mgi_term() BIND(C,name='MGI_Term')
    end subroutine c_mgi_term
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

function mgi_open(channel, mode) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  character(len=*), intent(IN) :: mode
  integer :: status

  status = c_mgi_open(channel, mode(1:1))

end function mgi_open

function mgi_clos(channel) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  integer :: status

  status = c_mgi_clos(channel)

end function mgi_clos

function mgi_read(channel, data, n, typ) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel, n
  integer, dimension(*), intent(OUT) :: data
  character(len=*), intent(IN) :: typ
  integer :: status

  status = c_mgi_read(channel, data, n, typ(1:1))

end function mgi_read

function mgi_write(channel, data, n, typ) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel, n
  integer, dimension(*), intent(IN) :: data
  character(len=*), intent(IN) :: typ
  integer :: status

  status = c_mgi_write(channel, data, n, typ(1:1))

end function mgi_write

subroutine mgi_term()
  use c_mgi_interfaces
  implicit none

  call c_mgi_term()

end subroutine mgi_term
