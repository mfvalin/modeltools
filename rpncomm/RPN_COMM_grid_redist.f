*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
      integer function RPN_COMM_grid_redist_test()
      use rpn_comm
      implicit none
      integer :: RPN_COMM_grid_redist
      external :: RPN_COMM_grid_redist
      integer, pointer, dimension(:,:,:) :: localarray
      integer, pointer, dimension(:,:,:) :: globalarray
      integer :: gni, gnj
      integer :: lminx,lmaxx,lminy,lmaxy
      integer, dimension(pe_nx) :: countx,offsetx
      integer, dimension(pe_ny) :: county,offsety
      integer :: i, ii, j, jj, k, ierr, ltok, status
      integer :: ix0, ixn, jy0, jyn
      integer :: i0, in, j0, jn
      integer, parameter :: lni=3
      integer, parameter :: lnj=5
      integer, parameter :: nox=2
      integer, parameter :: noy=2
      integer, dimension(nox) :: outpe_x
      integer, dimension(noy) :: outpe_y
      integer, dimension(2) :: zlist
      integer :: nk, nz
      
      if(pe_me==0) write(rpn_u,*)'grid redistribution test',
     %    pe_tot_grid,pe_nx,pe_ny

      outpe_x(nox) = pe_nx-1
      outpe_y(noy) = pe_ny-1
      outpe_x(1) = 0
      outpe_y(1) = 0
      gni = pe_nx*lni
      gnj = pe_ny*lnj
      nk = 4
      nk = nox*noy
      call RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      call RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      allocate(localarray(lminx-1:lmaxx+2,lminy-2:lmaxy+1,nk))
      localarray = 0
      do k = 1,nk
      do j = lminy,lmaxy
      do i = lminx,lmaxx
        localarray(i,j,k) = k + 10*j + 1000*i
      enddo
      enddo
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(pe_me==0)then
        write(rpn_u,100)'outpe_x =',outpe_x
        write(rpn_u,100)'outpe_y =',outpe_y
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(pe_me==(pe_nx*pe_ny-1)) then
        write(rpn_u,100)'local array allocated',pe_mex,pe_mey,
     %                lminx,lmaxx,lminy,lmaxy,nk
        do j = lmaxy,lminy,-1
          write(rpn_u,100)' ',localarray(lminx:lmaxx,j,max(1,nk-1))
        enddo
100     format(A,10I6)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      nullify(globalarray)
      ix0 = lni+(lni+1)/2-1
      ixn = gni-2*lni+(lni+1)/2+1
      jy0 = (lnj+1)/2 + lnj-1
      jyn = gnj-2*lnj+(lnj+1)/2+1
      i0=1
      if(ix0>lminx) i0 = i0 + (ix0-lminx)
      in=lni
      if(lmaxx>ixn) in = in - (lmaxx-ixn)
      j0=1
      if(jy0>lminy) j0 = j0 + (jy0-lminy)
      jn=lnj
      if(lmaxy>jyn) jn = jn - (lmaxy-jyn)
      do j = 1,noy
      do i = 1,nox
        if(outpe_x(i)==pe_mex .and. outpe_y(j)==pe_mey)then
!          write(rpn_u,*)'global array allocated',ix0, ixn, jy0, jyn
          allocate(globalarray(ix0:ixn,jy0:jyn,2))
        endif
      enddo
      enddo
      if(.not. associated(globalarray))allocate(globalarray(1,1,2))
      globalarray = 0
      do i=1,pe_nx*pe_ny
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(pe_me==i-1)
     %  write(rpn_u,*)'global array dims',
     %                pe_mex,pe_mey,shape(globalarray)
      enddo
      zlist = -1
      nz = 1
      ltok = 1
      if(.true.)then
      status = RPN_COMM_grid_redist(
     % localarray,1-1,lni+2,i0,in,1-2,lnj+1,j0,jn,nk,
     % globalarray,ix0,ixn,jy0,jyn,nz,zlist,
     % gni,gnj,outpe_x,nox,outpe_y,noy,
     % ltok)
      RPN_COMM_grid_redist_test = status
      if(pe_me==0)write(rpn_u,*)'status=',status
      endif
      if(status==-1)goto 777
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!      if(pe_me==outpe_x(nox)+pe_nx*outpe_x(noy)) then
      do i = 1,nox
      do j = 1,noy
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(outpe_x(i)==pe_mex .and. outpe_y(j)==pe_mey)then
          write(rpn_u,101)'global array =',
     %          pe_mex,pe_mey,ix0,ixn,jy0,jyn
101       format(A,20I6.5)
          do jj = jyn,jy0,-1
!            write(rpn_u,*)jj,i0,in,5
            write(rpn_u,101)' ',globalarray(ix0:ixn,jj,1)
          enddo
        endif
      enddo
      enddo
!      endif
777   if(associated(globalarray)) deallocate(globalarray)
      deallocate(localarray)
      return
      end function RPN_COMM_grid_redist_test
*****************************************
*                                       *
*             W A R N I N G             *
*            C O D E    I S             *
*                B E T A                *
*                                       *
*****************************************
      integer function RPN_COMM_grid_redist(
     %           zin,mini,maxi,i0,in,minj,maxj,j0,jn,nk,
     %           zout,ix0,ixn,jy0,jyn,nz,zlist,
     %           gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y,
     %           ltok)
      use rpn_comm
      implicit none
      integer, intent(IN) :: mini,maxi,i0,in,minj,maxj,j0,jn,nk,ltok
      integer, intent(IN) :: ix0,ixn,jy0,jyn,nz
      integer, intent(IN) :: noutpe_x,noutpe_y,gni,gnj
!      integer, dimension(ltok,mini:maxi,minj:maxj,nk), 
      integer, dimension(mini:maxi,minj:maxj,nk), 
     %         intent(IN) :: zin  ! zin(:,io:in,j0:jn,:) is useful
!      integer, dimension(ltok,ix0:ixn,jy0:jyn,nz), intent(OUT) :: zout
      integer, dimension(ix0:ixn,jy0:jyn,nz), intent(OUT) :: zout
      integer, dimension(nz), intent(OUT) :: zlist  ! list of 2D fields returned to this processor
      integer, dimension(noutpe_x) :: outpe_x    ! list of columns where there are PEs doing IO
      integer, dimension(noutpe_y) :: outpe_y    ! list of rows where there are PEs doing IO
!
      integer, dimension(pe_nx) :: istart_g, iend_g, countx, offsetx
      integer, dimension(pe_ny) :: jstart_g, jend_g, county, offsety
      integer, dimension(pe_nx) :: gdispl_x, gcounts_x, gsize_x
      integer, dimension(pe_ny) :: gdispl_y, gcounts_y
      integer, dimension(noutpe_x) :: level_x
      integer :: istart, iend, jstart, jend
      integer, pointer, dimension(:)     :: dest, dest2, ptr1d
      integer, pointer, dimension(:,:)   :: source2
      integer, pointer, dimension(:,:,:) :: source
      integer :: lminx, lmaxx, lminy, lmaxy
      integer :: i, j, k, l, kbot, ktop, n2d, kout
      integer :: lev0, narrays, ierr, temp
      logical :: needed_for_pass2, size_error, no_holes
!
      RPN_COMM_grid_redist = -1  ! return -1 if an error occurred
      if(ltok .ne. 1 ) return   ! for now, ltok MUST be = 1 (feature not implemented yet)
      if(nk > noutpe_x*noutpe_y) return !  OUCH !! , or maybe only do the first noutpe_x*noutpe_y
      nullify(dest,dest2,source,source2,ptr1d)
      if(pe_me==0 .and. .false.)then
        write(rpn_u,*)mini,maxi,i0,in,minj,maxj,j0,jn,nk
        write(rpn_u,*)ix0,ixn,jy0,jyn,nz
        write(rpn_u,*)zlist
        write(rpn_u,*)gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y,ltok
      endif
      narrays = nk
      needed_for_pass2 = .false.
      size_error = .false.
!
      call RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      call RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
!
      level_x = 0
      temp = 0
      do k = 1 , narrays  ! distribute narrays 2D fields over noutpe_x columns
        temp = temp + 1
        if(temp > noutpe_x) temp = 1
        level_x(temp) = level_x(temp) + 1
      enddo
!
      do j = 1 , pe_ny   ! start and end of valid data on PE(any,j) along y in global space
        jstart_g(j) = max(jy0,1+offsety(j))
        jend_g(j)   = min(jyn,offsety(j) + county(j))
      enddo
      jstart = jstart_g(pe_mey+1)  ! same as above but for my row
      jend   = jend_g(pe_mey+1)
      do i = 1 , pe_nx   ! start and end of valid data on PE(i,any) along x in global space
        istart_g(i) = max(ix0,1+offsetx(i))
        iend_g(i)   = min(ixn,offsetx(i) + countx(i))
      enddo
      istart = istart_g(pe_mex+1)  ! same as above but for my column
      iend   = iend_g(pe_mex+1)
      if((in-i0).ne.(iend-istart) .or. (jn-j0).ne.(jend-jstart)) then  ! size consistency problem ?
        if(istart<=iend) size_error = .true.  ! add error message 
        if(istart<=iend)then
          write(rpn_u,100)'error on pe',pe_mex,pe_mey,' =',
     %          i0,in,istart,iend,j0,jn,jstart,jend
100       format(A,2I2,A,10I8)
        endif
      endif
      if(pe_me==0 .and. .false.)then
        write(rpn_u,001)'countx=',countx
        write(rpn_u,001)'offsetx=',offsetx
        write(rpn_u,001)'county=',county
        write(rpn_u,001)'offsety=',offsety
        write(rpn_u,001)'JS limits=',jstart_g
        write(rpn_u,001)'JE limits=',jend_g
        write(rpn_u,001)'IS limits=',istart_g
        write(rpn_u,001)'IE limits=',iend_g
001     format(A,20I3)
      endif
!     check if a size error occurred somewhere in pe_grid, if so return -1
      call MPI_allreduce(MPI_IN_PLACE,size_error,1,MPI_LOGICAL,MPI_LOR,
     %                   pe_grid,ierr)
      if(size_error) return
!
      if(.false.)then
      write(rpn_u,001)'start of pass 1, pe=',
     %      pe_me,pe_mex,pe_mey,istart,iend,jstart,jend,
     %      istart_g,iend_g,jstart_g,jend_g,
     %      ix0,ixn,jy0,jyn
!      return
      endif
      n2d = -1
      if(jstart <= jend) then  ! there is something to do during pass 1 for this PE row
        do i = 1 , pe_nx                                  ! horizontal size of 2D fragments
          gsize_x(i)  = max(0,iend_g(i)-istart_g(i)+1) *  ! gsize_x may be zero for some PEs
     %                  max(0,jend-jstart+1)
        enddo
!
        kbot = 1
        do l = 1 , noutpe_x
          if(level_x(l)==0) cycle  ! no output to do on this column
          ktop = kbot + level_x(l) - 1
!         calculate counts and displacements for gatherv along x
          gcounts_x = 1 + gsize_x * level_x(l)  ! 1 + "horizontal size" * "levels to gather"
          gdispl_x(1) = 0
          do i = 2 , pe_nx
            gdispl_x(i) = gdispl_x(i-1) + gcounts_x(i-1)
          enddo
!
          if(istart<=iend) then ! there is valid data on this PE
            call flush(rpn_u)
            allocate(source(istart:iend,jstart:jend,kbot:ktop))     ! source buffer
            source(istart:iend,jstart:jend,kbot:ktop) = 
     %         zin(i0:in,j0:jn,kbot:ktop)   ! extract subarray from input array
          else                  ! no valid data on this PE
            allocate(source(1,1,1))
            source = 0
          endif
          if(pe_mex == outpe_x(l)) then  ! PEs doing gathering on this column
            n2d = level_x(l)       ! number of 2D fields for this PE column
            needed_for_pass2 = .true.
            allocate( dest( gdispl_x(pe_nx)+gcounts_x(pe_nx)+2 ) )   ! destination buffer 
            allocate( source2(ix0:ixn,jstart:jend) )
          endif
          call MPI_Gatherv(source,gcounts_x(pe_mex+1),MPI_INTEGER,
     %                     dest,gcounts_x,gdispl_x,MPI_INTEGER,
     %                     outpe_x(l),pe_myrow,ierr)
          deallocate(source)
          kbot = kbot + level_x(l)
        enddo
      else  ! no action on this row, allocate dummy source2 used in y direction gatherv
        do l = 1 , noutpe_x
          if(pe_mex == outpe_x(l)) then
            n2d = level_x(l)
            needed_for_pass2 = .true.
            allocate( source2(1,1))  ! dummy array for gatherv pass along y
          endif
        enddo
      endif   ! (jstart <= jend)
!
!      data gathering along X is done, we have the whole X axis in processor now
!
!      write(rpn_u,*)'end pass 1 :',pe_mex,pe_mey,needed_for_pass2,n2d
      if(needed_for_pass2 .and. n2d>0) then  ! nothing to do if n2d=0
        nullify(dest2)
        lev0 = 0
        no_holes=.true.  ! will be true if one or more PE rows contribute no data
                         ! this rigamarole is needed because some MPI implementations
                         ! of gatherv fail with zero data counts. the alternative is 
                         ! to create a new communicator for part of the column.
        do j = 1 , pe_ny ! gcounts_y contains the size of a "full x" "partial y" slab
          gcounts_y(j) = (ixn-ix0+1) * max(0,jend_g(j)-jstart_g(j)+1)
          no_holes = no_holes .and. gcounts_y(j)>0
          if(lev0==0 .and. gcounts_y(j)>0) lev0 = j ! find first non empty slab
          if(gcounts_y(j)==0) gcounts_y(j)=1  ! minimum count set to 1 or gatherv might fail
        enddo
        gdispl_y(1) = 0
        do j = 2 , pe_ny
          gdispl_y(j) = gdispl_y(j-1) + gcounts_y(j-1)
        enddo
        do l = 1 , noutpe_y
          if(pe_mey==outpe_y(l)) then  ! this PE is a data gatherer
            zlist = 0
            if(.not. no_holes)         ! will need an extra copy
     %        allocate( dest2( gdispl_y(pe_ny)+gcounts_y(pe_ny)+1 ) )
          endif
        enddo
        do k = 1 , n2d  ! loop over number of 2D fields to distribute along y
!         reassemble data for this "level" from dest(localx,localy,localk,npex)
!                                          into source2(globalx,localy,localk)
          do i = 1 , pe_nx   ! pe_nx pieces to reassemble
            if(jstart<=jend) then
              ptr1d => dest(gdispl_x(i)+1:gdispl_x(i)+gcounts_x(i))
              call place(
     %           ptr1d,istart_g(i),iend_g(i),jstart,jend,n2d,
     %           source2,ix0,ixn,k)
              endif
          enddo
          kout = 1 + (k-1)/noutpe_y
          l = 1 + mod(k-1,noutpe_y)
          zlist(kout) = k
!          write(rpn_u,002)'mid pass 2, pe=',pe_mex,pe_mey,
!     %      k,n2d,temp,jstart,jend,ix0,ixn,jy0,jyn,gdispl_y,kout
          if(no_holes) then  ! all PEs contribute, can gather directly into zout
002   format(A,20I5)
            call MPI_Gatherv(
     %         source2,gcounts_y(pe_mey+1),MPI_INTEGER,
     %         zout(ix0,jy0,kout),gcounts_y,gdispl_y,MPI_INTEGER,
     %         outpe_y(l),pe_mycol,ierr)
          else  ! not all PEs contribute, nedd an intermediate array
            call MPI_Gatherv(         !  gatherv, then copy into zout
     %         source2,gcounts_y(pe_mey+1),MPI_INTEGER,
     %         dest2,gcounts_y,gdispl_y,MPI_INTEGER,
     %         outpe_y(l),pe_mycol,ierr)
            if(pe_mey == outpe_y(l)) then  ! we are on root PE of gather
!              write(rpn_u,002)'mid pass 2B, pe=',pe_mex,pe_mey,
!     %           k,n2d,temp,jstart,jend,ix0,ixn,jy0,jyn,gdispl_y,kout
              temp = 1 + gdispl_y(lev0) ! point to start of useful data
              do j = jy0 , jyn
              do i = ix0 , ixn
                zout(i,j,kout) = dest2(temp)
                temp = temp + 1
              enddo  ! j
              enddo  ! i
            endif
          endif
        enddo  ! k = 1 , n2d
      endif
7777  continue
      if(associated(dest2)) deallocate(dest2)
      if(associated(dest)) deallocate(dest)
      if(associated(source2)) deallocate(source2)
!      call mpi_barrier(pe_mycol,ierr)
!
!      write(rpn_u,*)'end of pass 1+2, pe=',pe_me
      RPN_COMM_grid_redist = 1  ! number of 2D reassembled fields returned
      return
!
      contains
      subroutine place(src,i0,in,j0,jn,nk,dest,ix0,ixn,k)
      implicit none
      integer :: i0,in,j0,jn,nk,ix0,ixn,k
      integer, dimension(i0:in,j0:jn,nk), intent(IN) :: src
      integer, dimension(ix0:ixn,j0:jn), intent(OUT) :: dest
      dest(i0:in,j0:jn) = src(i0:in,j0:jn,k)
      end subroutine place
!
      end function RPN_COMM_grid_redist
!
      subroutine RPN_COMM_grid_redist_2(zin,lni,lnj,nk,zout,gni,gnj,
     %           outpe_x,noutpe_x,outpe_y,noutpe_y)
!     zout only needs to be defined on PEs doing I/O (pe_mex in outpe_x AND pe_mey in outpe_y)
!
! NOTE:
!     lni should probably really be nilmax
!     lnj should probably really be njlmax
!
!     this MUST be fixed or disaster will strike
!     we may want to use RPN_COMM_topo2 here to compute things properly
!
      use rpn_comm
      implicit none
      integer, intent(IN) :: lni,lnj,nk,gni,gnj,noutpe_x,noutpe_y
      integer, dimension(lni,lnj,nk), intent(IN) :: zin
      integer, dimension(gni,gnj), intent(OUT) :: zout
      integer, dimension(noutpe_x) :: outpe_x    ! list of columns where there are PEs doing IO
      integer, dimension(noutpe_y) :: outpe_y    ! list of rows where there are PEs doing IO
!
!      include 'mpif.h'
!
      integer, allocatable, dimension(:)   :: sendcounts, sdispls
      integer, allocatable, dimension(:)   :: recvcounts, rdispls
      integer, allocatable, dimension(:)   :: recvbuf
      integer, allocatable, dimension(:,:) :: reorder
      integer :: i, j, k, l, i0, ilast, offset, ierr, the_target
      integer :: lminx, lmaxx, mxi
      integer, allocatable, dimension(:)   :: countx, offsetx
      integer :: lminy, lmaxy, mxj, slot
      integer, allocatable, dimension(:)   :: county, offsety
      integer, dimension(noutpe_x) :: level_x
!
      if(nk > noutpe_x*noutpe_y) call abort !  OUCH !! , or maybe only do the first noutpe_x*noutpe_y
!
      allocate(recvcounts(max(pe_nx,pe_ny)))
      allocate(rdispls(max(pe_nx,pe_ny)))
      allocate(sendcounts(max(pe_nx,pe_ny)))
      allocate(sdispls(max(pe_nx,pe_ny)))
!
      allocate(countx(pe_nx),offsetx(pe_nx))
      allocate(county(pe_ny),offsety(pe_ny))
!
      call RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      mxi = countx(1)   ! max tile size along x
      call RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      mxj = county(1)   ! max tile size along y
!
      level_x = 0
      do i = 1 , nk
        slot = 1 + mod(i-1,noutpe_x)
        level_x(slot) = level_x(slot) + 1
      enddo
      k = 1
      do l = 1 , noutpe_x
        if(pe_mex == outpe_x(l)) 
     %    allocate(recvbuf(mxi*lnj*level_x(l)*pe_nx))
        recvcounts(1) = countx(1)*lnj*level_x(l)
        rdispls(1) = 1
        do i = 2 , pe_nx
          recvcounts(i) = countx(i)*lnj*level_x(l)
          rdispls(i) = rdispls(i-1) + recvcounts(i-1)
        enddo
        call MPI_Gatherv(zin(1,1,k),lni*lnj*level_x(l),MPI_INTEGER,
     %                   recvbuf,recvcounts,rdispls,MPI_INTEGER,
     %                   outpe_x(l),pe_myrow,ierr)
        k = k + level_x(l)
      enddo
      if(pe_mex /= outpe_x(l)) goto 111   ! if I am not in an active column, nothing more to do
      allocate(reorder(gni,lnj))  ! reorder buffer
      do l = 1 , noutpe_y
!       fill reorder(gni,lnj) from recvbuf("lni",lnj,"lnk",pe_nx) level l
        recvcounts(1) = gni*county(1)
        rdispls(1) = 1
        do j = 2 , pe_ny
          recvcounts(j) = gni*county(j)
          rdispls(j) = rdispls(j-1) + recvcounts(j-1)
        enddo
        call MPI_Gatherv(reorder,gni*lnj,MPI_INTEGER,
     %                   zout,recvcounts,rdispls,MPI_INTEGER,
     %                   outpe_y(l),pe_mycol,ierr)
      enddo
111   deallocate(recvbuf)
      deallocate(reorder)
!
      return
!
!     ============================== OLD CODE kept for history ==============================
!
!     zin data shape  (lni,lnj,nk) on PE( 0:pe_nx-1  ,0:pe_ny-1 )
!
      if(nk <= pe_nx) then    ! alltoallv along x, then reorder, then gather along y
        sendcounts = 1        ! because zero length can hang some implementations of allgatherv
        sendcounts(1:nk) = lni*lnj   ! sending nk pieces of size lni by lnj, rest with size 1
        if(pe_mex < nk) then                ! real data will only be sent to first nk PEs along X
          recvcounts(1:pe_nx) = countx(1:pe_nx)*lnj  ! received data pieces countx(i) by lnj
          allocate(recvbuf(mxi*lnj*pe_nx))  ! receive buffer for pe_nx (max size along x) data pieces
          allocate(reorder(gni,lnj))        ! reorder buffer, input to gatherv
        else
          recvcounts = 1                     ! dummy collect with no need to reorder
          allocate(recvbuf(pe_nx))           ! will receive pe_nx words
          allocate(reorder(1,1))             ! no reorder buffer
        endif
        do k = 1 , pe_nx
          sdispls(k) = 1 + (min(k,nk)-1)*lni*lnj  ! offset of level k, no greater than nk
          if(pe_mex < nk) then         ! this pe will collect pe_nx pieces of data
            rdispls(k) = 1 + (k-1)*mxi*lnj
          else                         ! this pe only gets one word (because zero length might hang)
            rdispls(k) = k
          endif
        enddo
        call MPI_alltoallv(zin,    sendcounts,sdispls,MPI_INTEGER,
     %                     recvbuf,recvcounts,rdispls,MPI_INTEGER,
     %                     pe_myrow,ierr)
!
!       recvbuf data shape (lni,lnj,pe_nx) on PE( 0:nk-1 , 0:pe_ny-1 )
!
        if(pe_mex < nk) then   ! reorder and collect along Y into zout (first nk PEs only)
          do k = 1 , pe_nx     ! time to get full x axis together
            offset = (k-1) * mxi*lnj  ! offset in recv buffer for i piece no k, row 1
            i0 = offsetx(k) + 1
            ilast = i0 + countx(k) -1
            do j = 1 , lnj      ! copy and insert one slice
              reorder(i0:ilast,j) = recvbuf(offset:offset+ilast-i0)  ! i piece no k, row j
              offset = offset + lnj                                  ! next row
            enddo
          enddo
!
!       reorder data shape (gni,lnj) on PE( 0:nk , 0:pe_ny-1 )
!
          recvcounts(1:pe_ny) = gni*county(1:pe_ny)
          rdispls(1:pe_ny)    = gni*offsety(1:pe_ny) + 1
          the_target = mod(pe_mex,pe_ny)
          call MPI_gatherv(reorder,gni*lnj,            MPI_INTEGER,
     %                     zout,   recvcounts,rdispls, MPI_INTEGER,
     %                     the_target,pe_mycol,ierr)      ! root PEs on diagonal
!
!       zout data shape (gni,gnj) on PE(l,m) (l=0:nk-1 , m=  l mod pe_ny)
!
        endif
        goto 1111        ! clean exit
      endif
!
!     zin data shape  (lni,lnj,nk) on PE( 0:pe_nx-1  ,0:pe_ny-1 )
!
      if(nk <= pe_ny) then    ! alltoallv along y, then gather along x, then reorder
        sendcounts = 1        ! should be 0 but zero length can hang allgatherv
        sendcounts(1:nk) = lni*lnj
        if(pe_mey < nk) then  ! will be collecting and reordering lni by gnj
          recvcounts(1:pe_ny) = county(1:pe_ny)*lni
          allocate(recvbuf(lni*gnj))
          allocate(reorder(mxi*gnj,pe_nx))
        else                  ! dummy collect and no reorder
          recvcounts = 1
          allocate(recvbuf(pe_ny))
          allocate(reorder(1,1))
        endif
        do k = 1 , pe_ny
          sdispls(k) = 1 + (min(k,nk)-1)*lni*lnj  ! offset of level k, no greater than nk
          if(pe_mey < nk) then         ! this pe will collect pe_nx pieces of data
            rdispls(k) = lni*offsety(k)
          else                         ! this pe only gets one word (because zero length might hang)
            rdispls(k) = k
          endif
        enddo
        call MPI_alltoallv(zin,    sendcounts,sdispls,MPI_INTEGER,
     %                     recvbuf,recvcounts,rdispls,MPI_INTEGER,
     %                     pe_mycol,ierr)
!
!       recvbuf data shape (lni,lnj,pe_ny) on PE( 0:nk-1 ,0:pe_ny-1 )
!
        if(pe_mey < nk) then                ! gather and reorder
          recvcounts(1:pe_nx) = countx(1:pe_nx)*gnj              ! lni by gnj slices
          do i = 1 , pe_nx
            rdispls(i) = (i-1)*mxi*gnj + 1
          enddo
          the_target = mod(pe_mey,pe_nx)
          call MPI_gatherv(recvbuf, lni*gnj,            MPI_INTEGER,
     %                     reorder, recvcounts,rdispls, MPI_INTEGER,
     %                     the_target,pe_myrow,ierr)      ! root PEs on diagonal
!
!       reorder data shape (mxi*gnj,pe_nx) on PE(m,l) (l=0:nk-1 , m=  l mod pe_nx)
!
          if(pe_mex == the_target) then
            do i = 1 , pe_nx      ! time to reorder, pe_nx slices along i
              i0 = offsetx(i) + 1
              ilast = i0 + countx(i) - 1
              do j = 1 , gnj
                offset = (j-1) * mxi
                zout(i0:ilast,j) = reorder(offset:offset+ilast-i0,i)
              enddo
            enddo
          endif
        endif
!
!       zout data shape (gni,lnj) on PE(l,m) (l=0:nk-1 , m=  l mod pe_ny)
!

        goto 1111        ! clean exit
      endif
!
!     if we get here OUCH !! (split nk ib slices of max(pe_nx,pe_ny) ?
!
1111  deallocate(sendcounts,recvcounts,sdispls,rdispls)
      deallocate(recvbuf,reorder)
      return
      end
