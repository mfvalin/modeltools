	module f_udunits_2
	use ISO_C_BINDING
	implicit none
	include 'f_udunits_2.inc'

	contains
!=============================================================================
	type(UT_SYSTEM_PTR) function f_ut_read_xml(path)
	use ISO_C_BINDING
	implicit none
	character (len=*), intent(IN) :: path

	character (len=1), dimension(len_trim(path)+1), target :: temp

	interface
	function ut_read_xml(mypath) result(ut_system) bind(C,name='ut_read_xml')
	use ISO_C_BINDING
	implicit none
	type(C_PTR) :: ut_system
	type(C_PTR), value :: mypath
	end function ut_read_xml
	end interface

	if(path == "" )then
	  f_ut_read_xml%ptr = ut_read_xml(C_NULL_PTR)
	else
          temp = transfer( trim(path)//achar(0) , temp )
	  f_ut_read_xml%ptr = ut_read_xml(c_loc(temp))
	endif

	return
	end function f_ut_read_xml
!=============================================================================
	subroutine f_ut_free_system(ut_system)
	use ISO_C_BINDING
	implicit none
	type(UT_SYSTEM_PTR), intent(IN), target :: ut_system

	interface
	subroutine ut_free_system(system) bind(C,name='ut_free_system')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: system
        end subroutine ut_free_system
	end interface

	call ut_free_system(ut_system%ptr)
	return
	end subroutine f_ut_free_system
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_get_unit_by_name(ut_system,name)  ! not recommended
	use ISO_C_BINDING
	implicit none
	type(UT_SYSTEM_PTR), intent(IN), target :: ut_system
	character (len=*), intent(IN) :: name

	character (len=1), dimension(len_trim(name)+1), target :: temp

	interface
	type(C_PTR) function ut_get_unit_by_name(ut_system,name) bind(C,name='ut_get_unit_by_name')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: ut_system
	type(C_PTR), value :: name
	end function ut_get_unit_by_name
	end interface

	temp = transfer( trim(name)//achar(0) , temp )
	f_ut_get_unit_by_name%ptr = ut_get_unit_by_name(ut_system%ptr,c_loc(temp))

	end function f_ut_get_unit_by_name
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_get_unit_by_symbol(ut_system,symbol)  ! not recommended
	use ISO_C_BINDING
	implicit none
	type(UT_SYSTEM_PTR), intent(IN), target :: ut_system
	character (len=*), intent(IN) :: symbol

	character (len=1), dimension(len_trim(symbol)+1), target :: temp

	interface
	type(C_PTR) function ut_get_unit_by_symbol(ut_system,symbol) bind(C,name='ut_get_unit_by_symbol')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: ut_system
	type(C_PTR), value :: symbol
	end function ut_get_unit_by_symbol
	end interface

	temp = transfer( trim(symbol)//achar(0) , temp )
	f_ut_get_unit_by_symbol%ptr = ut_get_unit_by_symbol(ut_system%ptr,c_loc(temp))

	end function f_ut_get_unit_by_symbol
!=============================================================================
	integer function f_ut_format(ut_unit,buffer,options)  ! this is the most useful routine
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN), target :: ut_unit
	character (len=*), intent(OUT) :: buffer
	integer, intent(IN) :: options

	integer(C_SIZE_T) :: buflen
	character (len=1), dimension(len(buffer)), target :: temp
	integer(C_INT) :: opt
	integer :: i, blen

	interface
	integer(C_INT) function ut_format(ut_unit,buffer,buflen,options)  bind(C,name='ut_format')
	use ISO_C_BINDING
	type(C_PTR), value :: ut_unit
	type(C_PTR), value :: buffer
	integer(C_SIZE_T), value :: buflen
	integer(C_INT), value :: options
	end function ut_format
	end interface

	buflen=len(buffer)
	opt = options + UT_ASCII
	temp=" "
	blen = ut_format(ut_unit%ptr,c_loc(temp),buflen,opt)
	f_ut_format = blen
	if(blen <= 0) then
	  buffer="ERROR"
	  return
	endif
	buffer = ""
	do i=1,blen
	  buffer(i:i)=temp(i)
	enddo

	end function f_ut_format
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_parse(ut_system,symbol)  ! this is the most useful routine
	use ISO_C_BINDING
	implicit none
	type(UT_SYSTEM_PTR), intent(IN), target :: ut_system
	character (len=*), intent(IN) :: symbol

	integer(C_INT) :: encoding
	character (len=1), dimension(len_trim(symbol)+1), target :: temp

	interface
	type(C_PTR) function ut_parse(ut_system,symbol,encoding) bind(C,name='ut_parse')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: ut_system
	type(C_PTR), value :: symbol
	integer(C_INT), value :: encoding
	end function ut_parse
	end interface

	encoding = UT_ASCII
	temp = transfer( trim(symbol)//achar(0) , temp )
	f_ut_parse%ptr = ut_parse(ut_system%ptr,c_loc(temp),encoding)

	end function f_ut_parse
!=============================================================================
	subroutine f_ut_free(ut_unit)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: ut_unit

	interface
	subroutine ut_free(ut_unit) bind(C,name='ut_free')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: ut_unit
        end subroutine ut_free
	end interface

	call ut_free(ut_unit%ptr)
	return
	end subroutine f_ut_free
!=============================================================================
	integer function f_ut_compare(unit1,unit2)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1,unit2

	interface
	integer(C_INT) function ut_compare(unit1,unit2) bind(C,name='ut_compare')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	type(C_PTR), value :: unit2
	end function ut_compare
	end interface

	f_ut_compare = ut_compare(unit1%ptr,unit2%ptr)

	end function f_ut_compare
!=============================================================================
	logical function f_ut_same_system(unit1,unit2)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1,unit2

	interface
	integer(C_INT) function ut_same_system(unit1,unit2) bind(C,name='ut_same_system')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	type(C_PTR), value :: unit2
	end function ut_same_system
	end interface

	f_ut_same_system = ut_same_system(unit1%ptr,unit2%ptr) .ne. 0

	end function f_ut_same_system
!=============================================================================
	logical function f_ut_is_dimensionless(unit1)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1

	interface
	integer(C_INT) function ut_is_dimensionless(unit1) bind(C,name='ut_is_dimensionless')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	end function ut_is_dimensionless
	end interface

	f_ut_is_dimensionless = ut_is_dimensionless(unit1%ptr) .ne. 0

	end function f_ut_is_dimensionless
!=============================================================================
	logical function f_ut_are_convertible(unit1,unit2)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1,unit2

	interface
	integer(C_INT) function ut_are_convertible(unit1,unit2) bind(C,name='ut_are_convertible')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	type(C_PTR), value :: unit2
	end function ut_are_convertible
	end interface

	f_ut_are_convertible = ut_are_convertible(unit1%ptr,unit2%ptr) .ne. 0

	end function f_ut_are_convertible
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_root(unit1,base)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1
	integer(C_INT), intent(IN) :: base

	interface
	type(C_PTR) function ut_root(unit1,base) bind(C,name='ut_root')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	integer(C_INT), value :: base
	end function ut_root
	end interface

	f_ut_root%ptr = ut_root(unit1%ptr,base)

	end function f_ut_root
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_raise(unit1,base)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1
	integer(C_INT), intent(IN) :: base

	interface
	type(C_PTR) function ut_raise(unit1,base) bind(C,name='ut_raise')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	integer(C_INT), value :: base
	end function ut_raise
	end interface

	f_ut_raise%ptr = ut_raise(unit1%ptr,base)

	end function f_ut_raise
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_offset(unit1,base)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1
	real(C_DOUBLE), intent(IN) :: base

	interface
	type(C_PTR) function ut_offset(unit1,base) bind(C,name='ut_offset')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	real(C_DOUBLE), value :: base
	end function ut_offset
	end interface

	f_ut_offset%ptr = ut_offset(unit1%ptr,base)

	end function f_ut_offset
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_scale(base,unit1)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1
	real(C_DOUBLE), intent(IN) :: base

	interface
	type(C_PTR) function ut_scale(base,unit1) bind(C,name='ut_scale')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	real(C_DOUBLE), value :: base
	end function ut_scale
	end interface

	f_ut_scale%ptr = ut_scale(base,unit1%ptr)

	end function f_ut_scale
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_log(base,unit1)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1
	real(C_DOUBLE), intent(IN) :: base

	interface
	type(C_PTR) function ut_log(base,unit1) bind(C,name='ut_log')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	real(C_DOUBLE), value :: base
	end function ut_log
	end interface

	f_ut_log%ptr = ut_log(base,unit1%ptr)

	end function f_ut_log
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_clone(unit1)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1

	interface
	type(C_PTR) function ut_clone(unit1) bind(C,name='ut_clone')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	end function ut_clone
	end interface

	f_ut_clone%ptr = ut_clone(unit1%ptr)

	end function f_ut_clone
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_invert(unit1)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1

	interface
	type(C_PTR) function ut_invert(unit1) bind(C,name='ut_invert')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	end function ut_invert
	end interface

	f_ut_invert%ptr = ut_invert(unit1%ptr)

	end function f_ut_invert
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_multiply(unit1,unit2)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1,unit2

	interface
	type(C_PTR) function ut_multiply(unit1,unit2) bind(C,name='ut_multiply')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	type(C_PTR), value :: unit2
	end function ut_multiply
	end interface

	f_ut_multiply%ptr = ut_multiply(unit1%ptr,unit2%ptr)

	end function f_ut_multiply
!=============================================================================
	type(UT_UNIT_PTR) function f_ut_divide(unit1,unit2)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: unit1,unit2

	interface
	type(C_PTR) function ut_divide(unit1,unit2) bind(C,name='ut_divide')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: unit1
	type(C_PTR), value :: unit2
	end function ut_divide
	end interface

	f_ut_divide%ptr = ut_divide(unit1%ptr,unit2%ptr)

	end function f_ut_divide
!=============================================================================
	type(CV_CONVERTER_PTR) function f_ut_get_converter(from,to)
	use ISO_C_BINDING
	implicit none
	type(UT_UNIT_PTR), intent(IN) :: from, to

	interface
	type(C_PTR) function ut_get_converter(from,to) bind(C,name='ut_get_converter')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: from
	type(C_PTR), value :: to
	end function ut_get_converter
	end interface

	f_ut_get_converter%ptr = ut_get_converter(from%ptr,to%ptr)

	end function f_ut_get_converter
!=============================================================================
	subroutine f_cv_free(converter)
	use ISO_C_BINDING
	implicit none
	type(CV_CONVERTER_PTR), intent(IN) :: converter

	interface
	subroutine cv_free(converter) bind(C,name='cv_free')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: converter
        end subroutine cv_free
	end interface

	call cv_free(converter%ptr)
	return
	end subroutine f_cv_free
!=============================================================================
	real function f_cv_convert_float(converter,what)
	use ISO_C_BINDING
	implicit none
	type(CV_CONVERTER_PTR), intent(IN) :: converter
	real(C_FLOAT), intent(IN) :: what

	interface
	real(C_FLOAT) function cv_convert_float(converter,what) bind(C,name='cv_convert_float')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: converter
	real(C_FLOAT), value :: what
	end function
	end interface

	f_cv_convert_float = cv_convert_float(converter%ptr,what)

	end function f_cv_convert_float
!=============================================================================
	real function f_cv_convert_double(converter,what)
	use ISO_C_BINDING
	implicit none
	type(CV_CONVERTER_PTR), intent(IN) :: converter
	real(C_DOUBLE), intent(IN) :: what

	interface
	real(C_DOUBLE) function cv_convert_double(converter,what) bind(C,name='cv_convert_double')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: converter
	real(C_DOUBLE), value :: what
	end function
	end interface

	f_cv_convert_double = cv_convert_double(converter%ptr,what)

	end function f_cv_convert_double
!=============================================================================
	subroutine f_cv_convert_floats(converter,what,count,dest)
	use ISO_C_BINDING
	implicit none
	type(CV_CONVERTER_PTR), intent(IN) :: converter
	real(C_FLOAT), intent(IN), dimension(*) :: what
	real(C_FLOAT), intent(OUT), dimension(*) :: dest
	integer, intent(IN) :: count

	type(C_PTR) :: dummy
	integer(C_SIZE_T) :: temp

	interface
	type(C_PTR) function cv_convert_floats(converter,what,count,dest) bind(C,name='cv_convert_floats')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: converter
	real(C_FLOAT), intent(IN), dimension(*) :: what
	real(C_FLOAT), intent(OUT), dimension(*) :: dest
	integer(C_SIZE_T), value :: count
	end function
	end interface

	temp = count
	dummy = cv_convert_floats(converter%ptr,what,temp,dest)

	end subroutine f_cv_convert_floats
!=============================================================================
	subroutine f_cv_convert_doubles(converter,what,count,dest)
	use ISO_C_BINDING
	implicit none
	type(CV_CONVERTER_PTR), intent(IN) :: converter
	real(C_DOUBLE), intent(IN), dimension(*) :: what
	real(C_DOUBLE), intent(OUT), dimension(*) :: dest
	integer, intent(IN) :: count

	type(C_PTR) :: dummy
	integer(C_SIZE_T) :: temp

	interface
	type(C_PTR) function cv_convert_doubles(converter,what,count,dest) bind(C,name='cv_convert_doubles')
	use ISO_C_BINDING
	implicit none
	type(C_PTR), value :: converter
	real(C_DOUBLE), intent(IN), dimension(*) :: what
	real(C_DOUBLE), intent(OUT), dimension(*) :: dest
	integer(C_SIZE_T), value :: count
	end function
	end interface

	temp = count
	dummy = cv_convert_doubles(converter%ptr,what,temp,dest)

	end subroutine f_cv_convert_doubles
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================
	end module f_udunits_2
