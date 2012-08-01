	program test_f_udunits_2
!	use ISO_C_BINDING
	use f_udunits_2
	implicit none
!	include 'f_udunits_2.inc'
	integer, target :: scrap
	type(ut_system_ptr) :: sys

	sys%ptr = c_loc(scrap)
	stop
	end