	type(UT_SYSTEM_PTR) function f_ut_read_xml(path)
	use ISO_C_BINDING
	implicit none
	include  'f_udunits_2.inc'
	character (len=*), intent(IN) :: path

	character (len=1) :: temp

	interface
	function ut_read_xml(mypath) result(ut_system) bind(C,name='ut_read_xml')
	use ISO_C_BINDING
	implicit none
	type(c_ptr) :: ut_system
	character (len=1) :: mypath
	end function ut_read_xml
	end interface

	f_ut_read_xml%ptr = ut_read_xml(transfer( trim(path)//achar(0) , temp ))

	return
	end function f_ut_read_xml
