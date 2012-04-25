	program test_001
	implicit none
	external RPN_COMM_init, UserInit
	integer RPN_COMM_dist_test
	external RPN_COMM_dist_test
	integer Pelocal,Petotal,Pex,Pey,ierr

	Pex = 0
	Pey = 0
	call RPN_COMM_init(UserInit,Pelocal,Petotal,Pex,Pey)
!	print *,' Pelocal,Petotal,Pex,Pey =',
!     %          Pelocal,Petotal,Pex,Pey
        ierr=RPN_COMM_dist_test(Petotal)
        call RPN_COMM_finalize(ierr)
	stop
	end
	subroutine UserInit(NX,NY)
	return
	end
