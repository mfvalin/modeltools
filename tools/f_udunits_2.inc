	type, bind(C) :: UT_SYSTEM_PTR
	type(C_PTR) :: ptr
	end type

	type, bind(C) :: UT_UNIT_PTR
	type(C_PTR) :: ptr
	end type

	type, bind(C) :: CV_CONVERTER_PTR
	type(C_PTR) :: ptr
	end type

	type, bind(C) :: UT_STATUS
	integer(C_INT) :: value
	end type

	integer, parameter :: UT_ASCII = 0
	integer, parameter :: UT_ISO_8859_1 = 1
	integer, parameter :: UT_LATIN1 = UT_ISO_8859_1
	integer, parameter :: UT_UTF8 = 2

	integer, parameter :: UT_NAMES = 4
	integer, parameter :: UT_DEFINITION = 8
