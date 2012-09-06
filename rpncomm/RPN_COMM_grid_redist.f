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
      subroutine RPN_COMM_grid_redist_2(zin,lni,lnj,nk,zout,gni,gnj,outpe_x,noutpe_x,outpe_y,noutpe_y)
!     zout only needs to be defined on PEs doing I/O (pe_mex in outpe_x AND pe_mey in outpe_y)
!
! NOTE:
!     lni should probably really be nilmax
!     lnj should probably really be njlmax
!     ltok necessaire ?
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
      if(nk > noutpe_x*noutpe_y) call abort !  ouch, or maybe only do the  first noutpe_x*noutpe_y
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
        if(pe_mex == outpe_x(l)) allocate(recvbuf(mxi*lnj*level_x(l)*pe_nx))
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
