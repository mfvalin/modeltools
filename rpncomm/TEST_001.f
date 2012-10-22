	program test_001
	implicit none
	external :: RPN_COMM_grid_redist_test
	integer :: RPN_COMM_grid_redist_test
	external RPN_COMM_init, TestUserInit, get_a_free_unit
        integer :: get_a_free_unit
	integer :: RPN_COMM_dist_test
	external RPN_COMM_dist_test
	integer :: Pelocal,Petotal,Pex,Pey,ierr,iun,test_to_perform

	Pex = 0
	Pey = 0
!       UserInit supplied by TEST_helpers.f
	call RPN_COMM_init(TestUserInit,Pelocal,Petotal,Pex,Pey)
!	print *,' Pelocal,Petotal,Pex,Pey =',
!     %          Pelocal,Petotal,Pex,Pey
        iun=get_a_free_unit()
        open(UNIT=iun,FILE='TEST_001.cfg',STATUS='OLD')
        read(UNIT=iun,FMT=*)test_to_perform
        close(UNIT=iun)
        if(IAND(test_to_perform,1)==1)then
          ierr=RPN_COMM_dist_test(Petotal)
        endif
        if(IAND(test_to_perform,2)==2)then
          ierr=RPN_COMM_grid_redist_test()
        endif
        call RPN_COMM_finalize(ierr)
	stop
	end
