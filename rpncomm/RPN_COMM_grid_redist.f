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

*****************************************
*                                       *
*             W A R N I N G             *
*            C O D E    I S             *
*            U N T E S T E D            *
*          I N C O M P L E T E          *
*              B U G G E D              *
*                                       *
*****************************************
      integer function RPN_COMM_grid_redist(ltok,
     %           zin,mini,maxi,i0,in,minj,maxj,j0,jn,nk,
     %           zout,ix0,ixn,jy0,jyn,nz,zlist,
     %           gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y)
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
      integer :: temp, narrays, ierr
      logical :: needed_for_pass2, size_error, no_holes
!
      RPN_COMM_grid_redist = -1  ! return -1 if an error occurred
      if(ltok .ne. 1 ) return   ! for now, ltok MUST be = 1 (feature not implemented yet)
      if(nk > noutpe_x*noutpe_y) return !  OUCH !! , or maybe only do the first noutpe_x*noutpe_y
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
        jstart_g(j) = max(jy0,offsety(j))
        jend_g(j)   = min(jyn,offsety(j) + county(j) -1)
      enddo
      jstart = jstart_g(pe_mey)  ! same as above but for my row
      jend   = jend_g(pe_mey)
      do i = 1 , pe_nx   ! start and end of valid data on PE(i,any) along x in global space
        istart_g(i) = max(ix0,offsetx(i))
        iend_g(i)   = min(ixn,offsetx(i) + countx(i) -1)
      enddo
      istart = istart_g(pe_mex)  ! same as above but for my column
      iend   = iend_g(pe_mex)
      if((in-i0).ne.(iend-istart) .or. (jn-j0).ne.(jend-jstart)) then  ! size consistency problem ?
        if(istart<=iend) size_error = .true.  ! add error message 
      endif
!     check if a size error occurred somewhere in pe_grid, if so return -1
      call MPI_allreduce(MPI_IN_PLACE,size_error,1,MPI_LOGICAL,MPI_LOR,
     %                   pe_grid,ierr)
      if(size_error) return
!
      if(jstart <= jend) then  ! there is something to do during pass 1 for this PE row
        do i = 1 , pe_nx                                  ! horizontal size of 2D fragments
          gsize_x(i)  = max(0,iend_g(i)-istart_g(i)+1) *  ! gsize_x may be zero for some PEs
     %                  max(0,jend-jstart+1)
        enddo
!
        kbot = 0
        do l = 1 , noutpe_x
          if(level_x(l)==0) cycle  ! no output to do on this column
          kbot = kbot + level_x(l)
          ktop = kbot + level_x(l) - 1
!         calculate counts and displacements for gatherv along x
          gcounts_x = gsize_x * level_x(l)  ! horizontal size * number of "levels" to gather
          gdispl_x(1) = 1
          do i = 2 , pe_nx
            gdispl_x(i) = gdispl_x(i-1) + gcounts_x(i-1) + 1
          enddo
!
          if(istart<=iend) then ! something to gather from this PE
            allocate(source(istart:iend,jstart:jend,kbot:ktop))     ! source buffer
            source(istart:iend,jstart:jend,kbot:ktop) = 
     %         zin(i0:in,j0:jn,kbot:ktop)
          else
            allocate(source(1,1,1))
            source = 0
          endif
          if(pe_mex == outpe_x(l)) then
            n2d = level_x(l)       ! number of 2D fields for this PE column
            needed_for_pass2 = .true.
            allocate( dest( gdispl_x(pe_nx)+gcounts_x(pe_nx)+1 ) )   ! destination buffer 
            allocate( source2(ix0:ixn,jstart:jend) )
          endif
          call MPI_Gatherv(source,gcounts_x(pe_mex)+1,MPI_INTEGER,
     %                     dest,gcounts_x+1,gdispl_x,MPI_INTEGER,
     %                     outpe_x(l),pe_myrow,ierr)
          deallocate(source)
        enddo
      else  ! no action on this row, allocate dummy source2 used in y direction gatherv
        do l = 1 , noutpe_x
          if(pe_mex == outpe_x(l)) then
            allocate( source2(1,1))  ! dummy array for gatherv pass along y
            source2 = 0
            needed_for_pass2 = .true.
          endif
        enddo
      endif   ! (jstart <= jend)
      if(needed_for_pass2) then
        nullify(dest2)
        temp = 0
        do j = 1 , pe_ny
          gcounts_y(j) = (ixn-ix0+1) * max(0,jend_g(j)-jstart_g(j)+1)
          no_holes = no_holes .and. gcounts_y(j)>0
          if(temp==0 .and. gcounts_y(j)>0) temp = j ! find first non hole
        enddo
        if(.not. no_holes) gcounts_y = gcounts_y + 1  ! add 1 to all counts if there are holes
        gdispl_y(1) = 1
        do j = 2 , pe_ny
          gdispl_y(j) = gdispl_y(j-1) + gcounts_y(j-1)
        enddo
        do l = 1 , noutpe_y
          if(pe_mey==outpe_y(l)) then
            zlist = 0
            if(.not. no_holes)
     %      allocate(dest2(gdispl_y(pe_ny)+gcounts_y(pe_ny)))
          endif
        enddo
        do k = 1 , n2d  ! loop over number of 2D fields to distribute along y
!         reorder data for this "level" from dest into source2
          do i = 1 , pe_nx   ! pe_nx pieces to reassemble
            ptr1d => dest(gdispl_x(i):gdispl_x(i)+gcounts_x(i)-1)
            call place(
     %           ptr1d,istart_g(i),iend_g(i),jstart,jend,n2d,
     %           source2,ix0,ixn,k)
          enddo
          kout = 1 + (k-1)/noutpe_y
          l = 1 + mod(k-1,noutpe_y)
          zlist(kout) = k
          if(no_holes) then
            call MPI_Gatherv(
     %         source2,gcounts_y(pe_mey),MPI_INTEGER,
     %         zout(ix0,jy0,kout),gcounts_y,gdispl_y,MPI_INTEGER,
     %         outpe_y(l),pe_mycol,ierr)
          else
!           gatherv, move into zout
            call MPI_Gatherv(
     %         source2,gcounts_y(pe_mey),MPI_INTEGER,
     %         dest2,gcounts_y,gdispl_y,MPI_INTEGER,
     %         outpe_y(l),pe_mycol,ierr)
            temp = gdispl_y(temp)
            do j = jy0 , jyn
            do i = ix0 , ixn
              zout(i,i,kout) = dest2(temp)
              temp = temp + 1
            enddo
            enddo
          endif
        
        enddo  ! k = 1 , n2d
        if(associated(dest2)) deallocate(dest2)
        deallocate(dest,source2)
      endif
!
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
