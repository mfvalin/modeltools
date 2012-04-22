!
!    Test verifying that rpn_comm_xch_halo and rpn_comm_transpose
!    are giving correct results. Furthermore, it returns a timing
!    for rpn_comm_xch_halo
!
        program test_002
	logical :: asynchronous  ! use asynchronous code
!        print *,'I am version ',VeRsion
        call rpn_comm_test_002
        stop
        end
        subroutine get_n_domains(ndomains, offset, err)
        integer n_domains
        common/ndomains/n_domains
        integer ndomains, offset, err
        character (len=128) SYSTEM_COMMAND
	SYSTEM_COMMAND="1"
        call get_env_var("TEST_DOMAINS",SYSTEM_COMMAND)
        if(SYSTEM_COMMAND == "" )SYSTEM_COMMAND="1"
        read(SYSTEM_COMMAND,*)ndomains
        n_domains=ndomains
        offset = 0
        err = 0
        return
        end
        subroutine rpn_comm_test_002
        implicit none
!
!       dummy test program implemented as a subroutne
!       because some compilers do not like POINTER
!       statements with variable dimensions in a main program
!
        integer  nptsx,nptsy,nptsz,ihalox,ihaloy
        common /the_problem/ nptsx,nptsy,nptsz,ihalox,ihaloy  ! common containing problem dimensions and halo size

        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
        integer, allocatable, dimension(:)  ::iarr,jarr
        integer, allocatable, dimension(:,:,:) ::data,data2,glob
!
        external sss
        integer npex,npey
        integer min3,max3,ierr, istat
        integer mini,maxi,nil,nilmax,ni0,i0
        integer minj,maxj,njl,njlmax,nj0,j0
        integer nzl, nzlmzx, nz0,irep,iter
        real*8 time1,time2,MPI_wtime,time_min,time_max,time_tot
        external MPI_wtime
!
        integer RPN_COMM_topo
        logical RPN_COMM_ngrank, RPN_COMM_grank
        integer RPN_COMM_init_multi_level, mygrid, mygridcomm, mymultigridcomm,myworldcomm,myallgridcomm
        integer RPN_COMM_colors
        integer RPN_COMM_comm
        integer test_grids
	integer RPN_COMM_option_L
        external RPN_COMM_mydomain
        external get_n_domains
        external test_grids
	external RPN_COMM_option_L
        integer domains, mydomain,irank,isize
        integer n_domains
        integer peercomm, npeers
        common/ndomains/n_domains
        character *6 is_async(2)
!
!
!       force a "callback" type initialization
!
        npex=0
        npey=0
        call RPN_COMM_mydomain (get_n_domains, mydomain,ierr)
!        print *,'This is member',mydomain+1,' of',n_domains,' domains'
!
        mygrid = RPN_COMM_init_multi_level(sss,Pelocal,Petotal,npex,npey,n_domains,1)
!
!       ============= determine resolution of MPI_wtine
!
        time1 = MPI_WTIME()
        time_min=1.0
        time_max=0.0
        time_tot = 0.0
        do irep = 1 , 100
          time2 = MPI_WTIME()
          time_min=min(time_min,time2-time1)
          time_max=max(time_max,time2-time1)
          time_tot = time_tot + time2-time1
          time1 = time2
        enddo
        if (Pelocal.eq.0) then
          print *,'MPI_WTIME cost(min,max,avg)=',real(time_min),real(time_max),real(time_tot)*.01
        endif
!
!	============= TEST for WORLD/MULTIGRID/GRID ====================
!
        if(test_grids(mygrid,Pelocal) /= 0) goto 9999
!
!	================= TEST for halo exchanger ===================
!
        call RPN_COMM_BCAST(nptsx,5,'MPI_INTEGER',0,'GRID',ierr)  ! broadcast problem dimensions to GRID
!
!       get data mapping (should really make 2 calls to
!       the 1D function RPN_COMM_topo)
!
        if(.not.RPN_COMM_grank('GRID')) goto 2222    ! if this PE is not a member of 'GRID', bail out
!
        istat= RPN_COMM_topo(nptsx,mini,maxi,nil,nilmax,ihalox,ni0,.TRUE.,.FALSE.)
        if(istat.ne.0) then
           write(*,*) 'Invalid distribution over x, abort',nptsx,npex
           goto 9999
        endif
!
        istat= RPN_COMM_topo(nptsy,minj,maxj,njl,njlmax,ihaloy,nj0,.FALSE.,.FALSE.)
        if(istat.ne.0) then
           write(*,*) 'Invalid distribution over y, abort',nptsy,npey
           goto 9999
        endif
!
        if (Pelocal.eq.0) then  ! print dimensions from PE 0
!          write(*,*) 'I AM ',Pelocal
           print *,'nptsx,nptsy,nptsz,mini,maxi,minj,maxj=',nptsx,nptsy,nptsz,mini,maxi,minj,maxj
           print *,'nil,njl,ihalox,ihaloy=',nil,njl,ihalox,ihaloy
           print *,'nilmax,njlmax=',nilmax,njlmax
        endif
!
!     allocate local arrays
!
        allocate(iarr(mini:maxi))
        allocate(jarr(minj:maxj))
        allocate(data(mini:maxi,minj:maxj,nptsz*2))
        allocate(glob((nptsx+2*ihalox),minj:maxj,nptsz))
!
!     fill arrays with markers
!
        i0=ni0
        j0=nj0
        call set_ijarr(iarr,jarr,mini,maxi,minj,maxj,i0,j0,nptsx,nptsy)
        call setup_arr(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl)
!
!     timing loop for halo exchange
!
	is_async(1) = 'ASYNC'
	is_async(2) = 'SYNC'
	do iter = 1 , 2
	ierr = RPN_COMM_option_L('async_exch',iter == 1)
        time1 = MPI_WTIME()
        time_min=1.0
        time_max=0.0
        time_tot = 0.0
        do irep=1,100 ! 100 repetitions of fully periodic halo exchange
          call RPN_COMM_xch_halo(data,mini,maxi,minj,maxj,nil,njl,nptsz,ihalox,ihaloy,.TRUE.,.TRUE.,nptsx,0)
          time2 = MPI_WTIME()
          time_min=min(time_min,time2-time1)
          time_max=max(time_max,time2-time1)
          time_tot = time_tot + time2-time1
          time1 = time2
        enddo
        if (Pelocal.eq.0) then
        endif
!
!          call affichage(data,mini,maxi,minj,maxj,nptsz) 
!
        if(Pelocal.eq.0 )then
           print *,is_async(iter)//'time (min,max,avg)=',real(time_min),real(time_max),real(time_tot)*.01
           print *,'pe=',Pelocal,                                          &
     &             ' Number of exchanges=', irep-1,                           &
     &             ' Time per exchange=',nint(1000000*time_tot/(irep-1)),            &
     &             ' microseconds'
        endif
	enddo
!
!     verify that there were no errors
!
        call vfy_arr(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
!
!	call RPN_COMM_FINALIZE(ierr)
!	STOP
!
!	================= TEST for transpose ===================
!
           call setup_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl)
           istat =  RPN_COMM_topo(nptsz,min3,max3,nzl,nzlmzx,0,nz0,.true.,.true.)
           allocate(data2((max3-min3+1),(maxj-minj+1),(nptsx+10)*2))
           if(.true.) then
             if(Pelocal.eq.0 )then
              print *,' size of data =',                                      &
     &             (maxi-mini+1),(maxj-minj+1),nptsz,                         &
     &             (maxi-mini+1)*(maxj-minj+1)*nptsz*2
              print *,' size of data2 =',                                     &
     &             (maxj-minj+1),(max3-min3+1),nptsx,                         &
     &             (max3-min3+1)*(maxj-minj+1)*(nptsx+10)*2
              print *,'nptsz,min3,max3,nzl,nzlmzx,nz0',                       &
     &             nptsz,min3,max3,nzl,nzlmzx,nz0
             endif
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,1,2)
!
             call vfy_xpos(data2,jarr,minj,maxj,njl,nz0,nzl,min3,max3,nptsx)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,-1,2)
!
             call vfy_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,1,2)
!
             call vfy_xpos(data2,jarr,minj,maxj,njl,nz0,nzl,min3,max3,nptsx)
!
             call RPN_COMM_transpose(data,mini,maxi,nptsx,(maxj-minj+1),min3,max3,nptsz,data2,-1,2)
!
             call vfy_arr2(data,iarr,jarr,mini,maxi,minj,maxj,nptsz,nil,njl,ihalox,ihaloy)
           endif
!
!     THE END !!
!
!
 2222   Continue
        call RPN_COMM_BARRIER('WORLD',ierr)
 9999	continue
        call RPN_COMM_FINALIZE(ierr)
        stop
        end
!
        subroutine vfy_xpos(z,jtab,minj,maxj,njl,nz0,nzl,min3,max3,nptsx)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers
!
!
        integer minj,maxj,njl,min3,max3,nptsx,nzl,nz0
        integer z(2,minj:maxj,min3:max3,nptsx)
        integer jtab(minj:maxj)
        integer i,j,k
        integer error
        error=0
        do i=1,nptsx
        do k=1,nzl
        do j=1,njl
            if(z(1,j,k,i)/32768 .ne. i+100) error=error+1
            if(iand(z(1,j,k,i),32767) .ne. jtab(j)) error=error+1
            if(z(2,j,k,i).ne.( nz0-1+k )) error = error + 1
        enddo
        enddo
        enddo
        if(error.ne.0) then 
          print *,error,' VERIFICATION ERRORS for Pe ',Pelocal
        endif
      if(Pelocal.eq.0) print *,'vfy_xpos : Verification performed, number of errors was ',error
        return
        end
!
        subroutine vfy_arr2(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers ignoring k
!
        integer minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy
        integer z(2,minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k,k0
        integer error
        error=0
        do k=1,nk
         k0=mod(k,3)
          do j=1,njl
          do i=1,nil
            if(z(1,i,j,k).ne.( itab(i)*32768 + jtab(j))) error = error + 1
            if(z(2,i,j,k).ne.( k ))  error = error + 1
          enddo
          enddo
        enddo
        if(error.ne.0) then 
          print *,error,' VERIFICATION ERRORS for Pe ',Pelocal
        endif
        if(Pelocal.eq.0) print *,'vfy_arr2 : Verification performed, number of errors was ',error
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(1,i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine vfy_arr(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        verify array containing markers
!
        integer minx,maxx,miny,maxy,nk,nil,njl,ihalox,ihaloy
        integer z(minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer error,k0
        integer i,j,k
        error=0
        do k=1,nk
         k0=mod(k,3)
          do j=1-ihaloy,njl+ihaloy
          do i=1-ihalox,nil+ihalox
            if(z(i,j,k).ne.( itab(i)*16384 + jtab(j)*4 )*4 + k0) error = error + 1
          enddo
          enddo
        enddo
        if(error.ne.0) then 
          print *,error,' VERIFICATION ERRORS for Pe ',Pelocal
        endif
        if(Pelocal.eq.0) print *,'vfy_arr : Verification performed, number of errors was ',error
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine set_ijarr(itab,jtab,minx,maxx,miny,maxy,i0,j0,nptsx,nptsy)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        precompute initial position tables
!
        integer minx,maxx,miny,maxy,nptsx,nptsy,i0,j0
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j
        do i=minx,maxx
          itab(i)=i0+i-1 + 100
          if(i+i0-1 .lt. 1)itab(i)=itab(i)+nptsx
          if(i+i0-1 .gt. nptsx) itab(i)=itab(i)-nptsx
        enddo
        do j=miny,maxy
          jtab(j)=j0+j-1 + 100
          if(j+j0-1.lt.1)jtab(j)=jtab(j)+nptsy
          if(j+j0-1.gt.nptsy) jtab(j)=jtab(j)-nptsy
        enddo
!        if(Pelocal.eq.3) print *,'i0,nptsx,itab=',i0,nptsx,itab
!        if(Pelocal.eq.3) print *,'j0,nptsy,jtab=',j0,nptsy,jtab
        return
        end
!
        subroutine setup_arr2(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl)
!!!!        use rpn_comm
        implicit none 
!
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!        fill array with markers ignoring k subscript
!
        integer minx,maxx,miny,maxy,nk
        integer nil,njl,z(2,minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k
        do k=1,nk
          do j=miny,maxy
          do i=minx,maxx
            z(1,i,j,k)=-1
            z(2,i,j,k)=-1
          enddo
          enddo
          do j=1,njl
          do i=1,nil
            z(1,i,j,k)=( itab(i)*32768 + jtab(j) )
            z(2,i,j,k)=k
          enddo
          enddo
        enddo
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(1,i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        subroutine setup_arr(z,itab,jtab,minx,maxx,miny,maxy,nk,nil,njl)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
!        fill array with markers
!
        integer minx,maxx,miny,maxy,nk
        integer nil,njl,z(minx:maxx,miny:maxy,nk)
        integer itab(minx:maxx),jtab(miny:maxy)
        integer i,j,k,k0
        do k=1,nk
         k0=mod(k,3)
          do j=miny,maxy
          do i=minx,maxx
            z(i,j,k)=-1
          enddo
          enddo
          do j=1,njl
          do i=1,nil
            z(i,j,k)=( itab(i)*16384 + jtab(j)*4 )*4 + k0
          enddo
          enddo
        enddo
        if((maxy-miny+1)*(maxx-minx+1) .gt. 25) return
        if(Pelocal.eq.0) then
          do j=maxy,miny,-1
            print 100,(z(i,j,1),i=minx,maxx)
          enddo
100            format(10z10)
        endif
        return
        end
!
        SUBROUTINE sss(nx,ny)
        implicit none 
        integer zouf,nx,ny
!
!        "callback routine" used to get initial topology
!        information
!
        common /the_problem/ nptsx,nptsy,nptsz,ihalox,ihaloy
        integer nptsx,nptsy,nptsz,ihalox,ihaloy
        open(5,file='TEST_data_001',form='FORMATTED')
        print *,'PEs =',nx*ny
        read(5,*)nx,ny,nptsx,nptsy,nptsz,ihalox,ihaloy
        print *, ' problem size is ',nptsx,' by ',nptsy,' by ',nptsz
        print *, ' halo size is ',ihalox,' by ',ihaloy
        print *, ' topology = ',nx,' by ',ny
        return
        end
!
        SUBROUTINE affichage(g,minx,maxx,miny,maxy,nk)
!!!!        use rpn_comm
        implicit none 
        integer Pelocal,Petotal
        common /PEs/ Pelocal,Petotal
!
        integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
        real g(minx:maxx,miny:maxy,nk)
!
        integer i,j,k,m
           if(Pelocal.eq.0) then
              write(*,*) 'matrice',Pelocal
 100       format(30F12.0)
           do j=maxy,miny,-1
              write(*,100) (g(i,j,1),i=minx,maxx)
           enddo
           endif
!
!
        return
        end
!
	integer function test_grids(mygrid,Pelocal)
	implicit none
	integer mygrid, Pelocal

	integer RPN_COMM_comm, RPN_COMM_colors
	external RPN_COMM_comm, RPN_COMM_colors
	integer mygridcomm, mymultigridcomm,myworldcomm,myallgridcomm,peercomm
	integer irank, isize, ierr, npeers

        test_grids = -1
        mymultigridcomm=RPN_COMM_comm('MULTIGRID')
        myallgridcomm=RPN_COMM_comm('ALLGRIDS')
        mygridcomm=RPN_COMM_comm('GRID')
        peercomm=RPN_COMM_comm('GRIDPEERS')
        myworldcomm=RPN_COMM_comm('WORLD')
!
        call RPN_COMM_BARRIER('GRID',ierr)
        call RPN_COMM_rank( 'GRID', irank ,ierr )
        call RPN_COMM_size( 'GRID', isize ,ierr )
        if(irank.eq.0)     print *,'BARRIER on GRID',mygridcomm,isize

        call RPN_COMM_BARRIER('MULTIGRID',ierr)
        call RPN_COMM_rank( 'MULTIGRID', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on MULTIGRID',mymultigridcomm

        call RPN_COMM_BARRIER('ALLGRIDS',ierr)
        call RPN_COMM_rank( 'ALLGRIDS', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on ALLGRIDS',myallgridcomm

        call RPN_COMM_BARRIER('WORLD',ierr)
        call RPN_COMM_rank( 'WORLD', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on WORLD',myworldcomm

        call RPN_COMM_BARRIER('GRIDPEERS',ierr)
        call RPN_COMM_rank( 'GRIDPEERS', irank ,ierr )
        if(irank.eq.0)     print *,'BARRIER on GRIDPEERS',peercomm

        call RPN_COMM_size( 'GRIDPEERS', npeers ,ierr )
        if(Pelocal.eq.0 )then
          print *,'mygrid=',mygrid,'myworld',myworldcomm,'mymultigridcomm=',mymultigridcomm,'mygridcomm=',mygridcomm
          print *,'ID=',RPN_COMM_colors('WORLD'),'/',RPN_COMM_colors('MULTIGRID'),'/',RPN_COMM_colors('GRID')
          print *,'peer to peer comm=',peercomm,' npeers=',npeers
        endif
        call RPN_COMM_BARRIER('WORLD',ierr)
        test_grids = 0
	return
	end