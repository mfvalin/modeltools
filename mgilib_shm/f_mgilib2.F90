!
! RMNLIB
! Copyright (C) 1995-2017 Environment Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.
!
module c_mgi_interfaces
 use ISO_C_BINDING
 implicit none
 interface
    function c_mgi_init(alias) result(channel) bind(C,name='MGI_Init')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: alias
      integer(C_INT) :: channel
    end function c_mgi_init
    function c_mgi_open(channel, mode) result(status) bind(C,name='MGI_Open')
      import :: C_CHAR, C_INT
      integer(C_INT), intent(IN), value :: channel
      character(C_CHAR), intent(IN), value :: mode
      integer(C_INT) :: status
    end function c_mgi_open
    function c_mgi_read(channel, data, nelm, dtyp) result(status) bind(C,name='MGI_Read')
      import :: C_CHAR, C_INT, C_PTR
      integer(C_INT), intent(IN), value :: channel, nelm
      type(C_PTR), intent(IN), value :: data
      character(C_CHAR), intent(IN), value :: dtyp
      integer(C_INT) :: status         ! nelm if OK, negetive error code otherwise
    end function c_mgi_read
    function c_mgi_write(channel, data, nelm, dtyp) result(status) bind(C,name='MGI_Write')
      import :: C_CHAR, C_INT, C_PTR
      integer(C_INT), intent(IN), value :: channel, nelm
      type(C_PTR), intent(IN), value :: data
      character(C_CHAR), intent(IN), value :: dtyp
      integer(C_INT) :: status         ! 0 if OK, number of unwritten toekns otherwise
    end function c_mgi_write
    function c_mgi_close(channel) result(status) bind(C,name='MGI_Close')
      import :: C_INT
      integer(C_INT), intent(IN), value :: channel
      integer(C_INT) :: status
    end function c_mgi_close
    function c_mgi_term() result(status) bind(C,name='MGI_Term')
      import :: C_INT
      integer(C_INT) :: status
    end function c_mgi_term
    function c_mgi_set_timeout(channel, timeout) result(status) bind(C,name='MGI_Set_timeout')
      import :: C_INT
      integer(C_INT), intent(IN), value :: channel
      integer(C_INT), intent(IN), value :: timeout
      integer(C_INT) :: status
    end function c_mgi_set_timeout
    function str_loc(what) result(address)
      import C_PTR
      character(len=*), intent(IN) :: what
      type(C_PTR) :: address
    end function str_loc
  end interface
end module

! ftnword f77_name (mgi_init) (char *channel_name, F2Cl lname);
function mgi_init(name) result(channel)
  use c_mgi_interfaces
  implicit none
  character(len=*), intent(IN) :: name
  integer :: channel

  character(C_CHAR), dimension(len(name)+1) :: temp
  temp = transfer(trim(name)//char(0),temp)
  channel = c_mgi_init(temp)
end function mgi_init
! ftnword f77_name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode);
function mgi_open(channel, mode) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  integer :: status
  character(len=1) :: mode

  status = c_mgi_open(channel, mode)
end 
! ftnword f77_name (mgi_read) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
function mgi_read(channel, data, nelm, dtype) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  integer, dimension(*), intent(OUT), target :: data
  integer, intent(IN) :: nelm
  character(len=1), intent(IN) :: dtype
  integer :: status

  status = c_mgi_read(channel, C_LOC(data), nelm, dtype)
end function mgi_read
function mgi_read_c(channel, data, nelm, dtype) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  character(len=*), intent(OUT) :: data
  integer, intent(IN) :: nelm
  character(len=1), intent(IN) :: dtype
  integer :: status

  if(len(data) < nelm) then
    status = -1
    return
  endif

  status = c_mgi_read(channel, STR_LOC(data), nelm, dtype)
end function mgi_read_c
! ftnword f77_name (mgi_write) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
function mgi_write(channel, data, nelm, dtype) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  integer, dimension(*), intent(IN), target :: data
  integer, intent(IN) :: nelm
  character(len=1), intent(IN) :: dtype
  integer :: status

  status = c_mgi_write(channel, C_LOC(data), nelm, dtype)
end function mgi_write
function mgi_write_c(channel, data, nelm, dtype) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  character(len=*), intent(IN) :: data
  integer, intent(IN) :: nelm
  character(len=1), intent(IN) :: dtype
  integer :: status

  if(len(data) < nelm) then
    status = -1
    return
  endif

  status = c_mgi_write(channel, STR_LOC(data), nelm, dtype)
end function mgi_write_c
! ftnword f77_name (mgi_clos) (ftnword *f_chan);
function mgi_clos(channel) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  integer :: status

  status = c_mgi_close(channel)
end function mgi_clos
! ftnword f77_name (mgi_term) ();
function mgi_term() result(status)
  use c_mgi_interfaces
  implicit none
  integer :: status

  status = c_mgi_term()
end function mgi_term
! void f77_name (mgi_set_timeout) (ftnword *chan, ftnword *timeout);
function mgi_set_timeout(channel, timeout) result(status)
  use c_mgi_interfaces
  implicit none
  integer, intent(IN) :: channel
  integer, intent(IN) :: timeout
  integer :: status

  status = c_mgi_set_timeout(channel, timeout)
end function mgi_set_timeout
