	program test_f_udunits_2
	use ISO_C_BINDING
	use f_udunits_2
	implicit none

	integer, target :: scrap
	type(UT_SYSTEM_PTR) :: sys
	type(CV_CONVERTER_PTR) :: conv1, conv2, conv3
	type(UT_UNIT_PTR) :: unit1, unit2, unit3
	character (len=*), parameter :: machin1 = "km/hour"
	character (len=*), parameter :: machin2 = "m/second"
	character (len=*), parameter :: machin3 = "knot"
	character (len=128) :: buffer
	real :: A,B
	real *8 :: AA, BB
	integer :: junk

	sys = f_ut_read_xml("")

	unit1 = f_ut_parse(sys,machin1)
	if(.not. c_associated(unit1%ptr)) print *,'ERROR in unit1'
	junk = f_ut_format(unit1,buffer,UT_NAMES)
	print *,"DEBUG: unit ",machin1,"='",trim(buffer),"'"
	junk = f_ut_format(unit1,buffer,UT_DEFINITION)
	print *,"DEBUG: unit ",machin1,"='",trim(buffer),"'"

	unit2 = f_ut_parse(sys,machin2)
	if(.not. c_associated(unit2%ptr)) print *,'ERROR in unit2'
	junk = f_ut_format(unit2,buffer,UT_NAMES)
	print *,"DEBUG: unit ",machin2,"='",trim(buffer),"'"
	junk = f_ut_format(unit2,buffer,UT_DEFINITION)
	print *,"DEBUG: unit ",machin2,"='",trim(buffer),"'"

	conv1 = f_ut_get_converter(unit1,unit2)   ! machin2 to machin1
	if(.not. c_associated(conv1%ptr)) print *,'ERROR in conv1'
	a =1.5
	B = f_cv_convert_float(conv1,A)
	print *,a,machin1,' is the same as',b,machin2

	unit3 = f_ut_parse(sys,machin3)
	if(.not. c_associated(unit3%ptr)) print *,'ERROR in unit3'
	junk = f_ut_format(unit3,buffer,UT_NAMES)
	print *,"DEBUG: unit ",machin3,"='",trim(buffer),"'"
	junk = f_ut_format(unit3,buffer,UT_DEFINITION)
	print *,"DEBUG: unit ",machin3,"='",trim(buffer),"'"

	conv2 = f_ut_get_converter(unit3,unit1)   ! machin3 to machin1
	if(.not. c_associated(conv2%ptr)) print *,'ERROR in conv2'
	a = 1.5
	B = f_cv_convert_float(conv2,A)
	print *,a,machin3,' is the same as',b,machin1

	conv3 = f_ut_get_converter(unit3,unit2)   ! machin3 to machin2
	if(.not. c_associated(conv3%ptr)) print *,'ERROR in conv3'
	a = 1.5
	B = f_cv_convert_float(conv3,A)
	print *,a,machin3,' is the same as',b,machin2

	call f_ut_free(unit1)
	call f_ut_free(unit2)
	call f_ut_free(unit3)
	call f_cv_free(conv1)
	call f_cv_free(conv2)
	call f_cv_free(conv3)
	call f_ut_free_system(sys)
!
	stop
	end