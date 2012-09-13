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
     %           zout,ix0,ixn,jy0,jyn,nz,
     %           gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y)
      use rpn_comm
      implicit none
      integer, intent(IN) :: mini,maxi,i0,in,minj,maxj,j0,jn,nk,ltok
      integer, intent(IN) :: ix0,ixn,jy0,jyn,nz
      integer, intent(IN) :: noutpe_x,noutpe_y,gni,gnj
      integer, dimension(ltok,mini:maxi,minj:maxj,nk), 
     %         intent(IN) :: zin  ! zin(:,io:in,j0:jn,:) is useful
      integer, dimension(ix0:ixn,jy0:jyn,nz), intent(OUT) :: zout
      integer, dimension(noutpe_x) :: outpe_x    ! list of columns where there are PEs doing IO
      integer, dimension(noutpe_y) :: outpe_y    ! list of rows where there are PEs doing IO
!
      integer, dimension(pe_nx) :: istart_g, iend_g, countx, offsetx
      integer, dimension(pe_nx) :: gdispl_x, gcounts_x, gsize_x
      integer, dimension(pe_ny) :: jstart_g, jend_g, county, offsety
      integer, dimension(pe_ny) :: gdispl_y, gcounts_y, gsize_y
      integer, dimension(noutpe_x) :: level_x
      integer :: istart, iend, jstart, jend, mxi, mxj
      integer, allocatable, dimension(:)   :: source, dest
      integer, allocatable, dimension(:,:) :: source2
      integer :: lminx, lmaxx, lminy, lmaxy
      integer :: i, j, k, l
      integer :: temp, narrays, ierr
      logical :: needed_for_pass2
!
      needed_for_pass2 = .false.
      RPN_COMM_grid_redist = -1
      if(nk > noutpe_x*noutpe_y) return !  OUCH !! , or maybe only do the first noutpe_x*noutpe_y
      narrays = nk
!
      call RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      mxi = countx(1)   ! max tile size along x
      call RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      mxj = county(1)   ! max tile size along y
!
!
      level_x = 0
      do k = 1 , narrays
        temp = 1 + mod(k-1,noutpe_x)
        level_x(temp) = level_x(temp) + 1
      enddo
!
      do j = 1 , pe_ny   ! start and end of tile(any,j) along y in global space
        jstart_g(j) = max(jy0,offsety(j))
        jend_g(j)   = min(jyn,offsetx(j) + county(j) -1)
      enddo
      jstart = jstart_g(pe_mey)  ! same as above but for this row
      jend = jend_g(pe_mey)
      do i = 1 , pe_nx   ! start and end of tile(i,any) along x in global space
        istart_g(i) = max(ix0,offsetx(i))
        iend_g(i)   = min(ixn,offsetx(i) + countx(i) -1)
        gsize_x(i)  = max(0,iend_g(i)-istart_g(i)+1)  ! horizontal size of tiles on this row
        gsize_x(i)  = gsize_x(i) * max(0,jend-jstart+1) ! gsize_x may be zero for some tiles
      enddo
      istart = istart_g(pe_mex)  ! same as above but for this column
      iend = iend_g(pe_mex)
!
      if( (in-i0-1) * (jn-j0-1) .ne. gsize_x(pe_mex) ) then  ! size consistency problem
        return  ! add error message above
      endif
!
      if(jstart <= jend) then  ! 
        k = 1
        do l = 1 , noutpe_x
!         calculate counts and displacements for gatherv along x
          gcounts_x = 1 + gsize_x * level_x(l)  ! horizontal size * number "levels" to gather
          gdispl_x(1) = 1
          do i = 2 , pe_nx
            gdispl_x(i) = gdispl_x(i-1) + gcounts_x(i-1)
          enddo
!
          allocate(source(gcounts_x(pe_mex)))     ! source buffer
          if(gcounts_x(pe_mex) > 1) then
!           fill source from zin
          endif
          if(pe_mex == outpe_x(l)) then
            needed_for_pass2 = .true.
            allocate( dest( gdispl_x(pe_nx)+gcounts_x(pe_nx)) )   ! destination buffer 
            allocate( source2(ix0:ixn,jstart:jend) )
          endif
          call MPI_Gatherv(source,gcounts_x(pe_mex),MPI_INTEGER,
     %                     dest,gcounts_x,gdispl_x,MPI_INTEGER,
     %                     outpe_x(l),pe_myrow,ierr)
          k = k + level_x(l)
!
          deallocate(source)
        enddo
      else  ! no action on this row, allocate dummy source2 used in y direction gatherv
        do l = 1 , noutpe_x
          if(pe_mex == outpe_x(l)) then
            allocate( source2(1,1))
            needed_for_pass2 = .true.
          endif
        enddo
      endif   ! (jstart <= jend)
!
      if(.not. needed_for_pass2) goto 111
!
      do l = 1 , noutpe_y
        do j = 1 , pe_ny
          gcounts_y(j) = (ixn-ix0+1) * max(0,jend_g(j)-jstart_g(j)+1)
        enddo
        gdispl_y(1) = 1
        do j = 2 , pe_ny
          gdispl_y(j) = gdispl_y(j-1) + gcounts_y(j-1)
        enddo
!       reorder data for this "level" from dest into source2
        call MPI_Gatherv(source2,gcounts_y(pe_mey),MPI_INTEGER,
     %                   zout,gcounts_y,gdispl_y,MPI_INTEGER,
     %                   outpe_y(l),pe_mycol,ierr)
      enddo
      deallocate(dest,source2)
!
111   RPN_COMM_grid_redist = 0
      return
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
      integer :: i, j, k, l, i0, ilast, j0, offset, ierr, the_target
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
