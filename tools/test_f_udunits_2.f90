	program test_f_udunits_2
!	use ISO_C_BINDING
	use f_udunits_2
	implicit none

	integer, target :: scrap
	type(UT_SYSTEM_PTR) :: sys
	type(CV_CONVERTER_PTR) :: conv1, conv2, conv3, time_cvt, time_cvt0
	type(UT_UNIT_PTR) :: unit1, unit2, unit3, unit4, unit5
	type(UT_UNIT_PTR) :: unit0, unit6, unit7, unit8, unit9, unit06
	type(UT_UNIT_PTR) :: erg, watt, watt_s, sec, time_base, sec0, time_base0
	character (len=*), parameter :: machin1 = "km/hour"
	character (len=*), parameter :: machin2 = "m/second"
	character (len=*), parameter :: machin3 = "2 knots"
	character (len=*), parameter :: machin4 = "watt"
	character (len=*), parameter :: machin5 = "radian"
	character (len=128) :: buffer
	real *8, parameter :: ONE = 1.01
	real *8, parameter :: TWO = 2.0
	real *8, parameter :: THREE = 3.03
	real *8, parameter :: FOUR = 4.04
	real *8, parameter :: TEN = 10.0
	real *8, parameter :: ZERO = 0.0
	real :: A,B
	real *8 :: AA, BB
	integer :: junk
	integer :: charset
	integer hour,minute,year,month,day
	real(C_DOUBLE) :: second, time, resolution

	charset = UT_ASCII
	sys = f_ut_read_xml("")

	year=-2000 ; month=8 ; day=8 ; hour=15 ; minute=45 ; second=15.55 ; resolution=0.0
	time = f_ut_encode_time(year,month,day,hour,minute,second)
	print *,"TIME=",time," =<",year,month,day,hour,minute,real(second)
	year=0 ; month=0 ; day=0 ; hour=0 ; minute=0 ; second=0.0 ; resolution=-999.999
	call f_ut_decode_time(time,year,month,day,hour,minute,second,resolution)
	print *,"TIME=",time," =>",year,month,day,hour,minute,real(second)," +-",real(resolution)

	sec0 = f_ut_parse(sys,"second",charset)
	time_base = f_ut_offset_by_time(sec0,f_ut_encode_time(2012,08,08,15,00,ZERO))
	time_base0 = f_ut_offset_by_time(sec0,f_ut_encode_time(2001,01,01,00,00,ZERO))
	time_cvt = f_ut_get_converter(sec0,time_base)
	if(.not. c_associated(time_cvt%ptr)) print *,'ERROR in time_cvt'
	time_cvt0 = f_ut_get_converter(sec0,time_base0)
	if(.not. c_associated(time_cvt0%ptr)) print *,'ERROR in time_cvt0'
	junk = f_ut_format(time_base,buffer,UT_NAMES)
	print *,"DEBUG: new unit of time_base='",trim(buffer),"'"
	print *,"converted time =",f_cv_convert_double(time_cvt,time)
	junk = f_ut_format(time_base0,buffer,UT_NAMES)
	print *,"DEBUG: new unit of time_base='",trim(buffer),"'"
	print *,"converted time =",f_cv_convert_double(time_cvt0,time)

	unit1 = f_ut_parse(sys,machin1,charset)
	if(.not. c_associated(unit1%ptr)) print *,'ERROR in unit1'
	junk = f_ut_format(unit1,buffer,UT_NAMES)
	print *,"DEBUG: unit UT_NAMES of ",machin1,"='",trim(buffer),"'"
	junk = f_ut_format(unit1,buffer,UT_DEFINITION)
	print *,"DEBUG: unit UT_DEFINITION of ",machin1,"='",trim(buffer),"'"

	unit2 = f_ut_parse(sys,machin2,charset)
	if(.not. c_associated(unit2%ptr)) print *,'ERROR in unit2'
	junk = f_ut_format(unit2,buffer,UT_NAMES)
	print *,"DEBUG: unit UT_NAMES of ",machin2,"='",trim(buffer),"'"
	junk = f_ut_format(unit2,buffer,UT_DEFINITION)
	print *,"DEBUG: unit UT_DEFINITION of ",machin2,"='",trim(buffer),"'"

	conv1 = f_ut_get_converter(unit1,unit2)   ! machin2 to machin1
	if(f_ut_are_convertible(unit1,unit2))  print *,machin1," can be converted into ",machin2
	if(.not. c_associated(conv1%ptr)) print *,'ERROR in conv1'
	a =1.5
	B = f_cv_convert_float(conv1,A)
	print *,a,machin1,' is the same as',b,machin2

	unit3 = f_ut_parse(sys,machin3,charset)
	if(.not. c_associated(unit3%ptr)) print *,'ERROR in unit3'
	junk = f_ut_format(unit3,buffer,UT_NAMES)
	print *,"DEBUG: unit UT_NAMES of ",machin3,"='",trim(buffer),"'"
	junk = f_ut_format(unit3,buffer,UT_DEFINITION)
	print *,"DEBUG: unit UT_DEFINITION of ",machin3,"='",trim(buffer),"'"

	conv2 = f_ut_get_converter(unit3,unit1)   ! machin3 to machin1
	if(f_ut_are_convertible(unit3,unit1))  print *,machin3," can be converted into ",machin1
	if(.not. c_associated(conv2%ptr)) print *,'ERROR in conv2'
	a = 1.5
	B = f_cv_convert_float(conv2,A)
	print *,a,machin3,' is the same as',b,machin1

	conv3 = f_ut_get_converter(unit3,unit2)   ! machin3 to machin2
	if(f_ut_are_convertible(unit3,unit2))  print *,machin3," can be converted into ",machin2
	if(.not. c_associated(conv3%ptr)) print *,'ERROR in conv3'
	a = 1.5
	B = f_cv_convert_float(conv3,A)
	print *,a,machin3,' is the same as',b,machin2
	print *,"COMPARE ",machin1," with ",machin1," = ",f_ut_compare(unit1,unit1), &
                ", same system =",f_ut_same_system(unit1,unit1)
	print *,"COMPARE ",machin1," with ",machin2," = ",f_ut_compare(unit1,unit2), &
                ", same system =",f_ut_same_system(unit1,unit2)
	print *,"COMPARE ",machin1," with ",machin3," = ",f_ut_compare(unit1,unit3), &
                ", same system =",f_ut_same_system(unit1,unit3)
	print *,"COMPARE ",machin2," with ",machin3," = ",f_ut_compare(unit2,unit3), &
                ", same system =",f_ut_same_system(unit2,unit3)

	unit4 = f_ut_parse(sys,machin4,charset)
	if(.not. c_associated(unit4%ptr)) print *,'ERROR in unit4'
	if(.not. f_ut_are_convertible(unit3,unit4))  print *,machin3," cannot be converted into ",machin4, &
	                                             ", same system =",f_ut_same_system(unit3,unit4)
	print *,machin4," is dimensionless =",f_ut_is_dimensionless(unit4)
	unit5 = f_ut_parse(sys,machin5,charset)
	if(.not. c_associated(unit5%ptr)) print *,'ERROR in unit5'
	print *,machin5," is dimensionless =",f_ut_is_dimensionless(unit5)

	junk = f_ut_format(unit4,buffer,UT_NAMES)
	print *,"DEBUG: unit ='",trim(buffer),"'"
	junk = f_ut_format(unit5,buffer,UT_NAMES)
	print *,"DEBUG: unit ='",trim(buffer),"'"
	unit6=f_ut_divide(unit4,unit5)
	print *,"DIVIDE ",machin4," by ",machin5
	junk = f_ut_format(unit6,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	print *,"MULTIPLY by ",machin5
	unit06=f_ut_multiply(unit6,unit5)
	junk = f_ut_format(unit06,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	print *,"SCALE by",THREE
	unit7=f_ut_scale(THREE,unit06)
	junk = f_ut_format(unit7,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	print *,"OFFSET by", FOUR
	unit8=f_ut_offset(unit7,FOUR)
	junk = f_ut_format(unit8,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	junk = f_ut_format(unit7,buffer,UT_NAMES)
	print *,"RAISE 3 ",trim(buffer)
	unit9=f_ut_raise(unit7,3)
	junk = f_ut_format(unit9,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	print *,"INVERT"
	unit0=f_ut_invert(unit9)
	junk = f_ut_format(unit0,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	call f_ut_free(unit0)
	print *,"root 2"
	unit0=f_ut_root(unit9,2)
	print *,"LOG",TWO
	unit0=f_ut_log(TWO,unit9)
	junk = f_ut_format(unit0,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	call f_ut_free(unit0)
	print *,"LOG",TEN
	unit0=f_ut_log(TEN,unit9)
	junk = f_ut_format(unit0,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	call f_ut_free(unit0)
	print *,"LOG",FOUR
	unit0=f_ut_log(FOUR,unit9)
	junk = f_ut_format(unit0,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"
	call f_ut_free(unit0)
	print *,"root 3"
	unit0=f_ut_root(unit9,3)
	junk = f_ut_format(unit0,buffer,UT_NAMES)
	print *,"DEBUG: new unit = '",trim(buffer),"'"

	watt = f_ut_parse(sys,"watt",charset)
	sec = f_ut_parse(sys,"second",charset)
	erg = f_ut_parse(sys,"erg",charset)
	watt_s = f_ut_multiply(watt,sec)
	junk = f_ut_format(watt_s,buffer,UT_NAMES)
	print *,"DEBUG: watt_s = '",trim(buffer),"'"
	junk = f_ut_format(erg,buffer,UT_NAMES)
	print *,"DEBUG: erg = '",trim(buffer),"'"
	if(f_ut_are_convertible(watt_s,erg))  print *,"watt.s can be converted into erg"

	call f_ut_free(unit1)
	call f_ut_free(unit2)
	call f_ut_free(unit3)
	call f_ut_free(unit4)
	call f_cv_free(conv1)
	call f_cv_free(conv2)
	call f_cv_free(conv3)
	call f_ut_free_system(sys)
!
	stop
	end