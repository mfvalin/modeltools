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
      subroutine RPN_COMM_grid_redist(zin,lni,lnj,nk,zout,gni,gnj)
!
! NOTE:
!     lni should really be nilmax
!     lnj should really be njlmax
!
!     this MUST be fixed or disaster will strike
!     we may want to use RPN_COMM_topo2 here to compute things properly
!
      use rpn_comm
      implicit none
      integer, intent(IN) :: lni,lnj,nk,gni,gnj
      integer, dimension(lni,lnj,nk), intent(IN) :: zin
      integer, dimension(gni,gnj), intent(OUT) :: zout
!
      include 'mpif.h'
!
      integer, allocatable, dimension(:)   :: sendcounts, sdispls
      integer, allocatable, dimension(:)   :: recvcounts, rdispls
      integer, allocatable, dimension(:)   :: recvbuf
      integer, allocatable, dimension(:,:) :: reorder
      integer :: i, j, k, i0, ilast, j0, offset, ierr
!
      allocate(recvcounts(max(pe_nx,pe_ny)))
      allocate(rdispls(max(pe_nx,pe_ny)))
      allocate(sendcounts(max(pe_nx,pe_ny)))
      allocate(sdispls(max(pe_nx,pe_ny)))
!
!     zin data shape  (lni,lnj,nk) on PE( 0:pe_nx-1  ,0:pe_ny-1 )
!
      if(nk <= pe_nx) then    ! alltoallv along x, then reorder, then gather along y
        sendcounts = 1        ! should be 0 but zero length can hang allgatherv
        sendcounts(1:nk) = lni*lnj
        if(pe_mex < nk) then    ! real data will only be sent to first nk PEs along X
          recvcounts = lni*lnj  ! data pieces lni by lnj
          allocate(recvbuf(lni*lnj*pe_nx))  ! receive buffer for pe_nx data pieces
          allocate(reorder(gni,lnj))   ! reorder buffer, input to gatherv
        else
          recvcounts = 1               ! dummy collect and no need to reorder
          allocate(recvbuf(pe_nx))
          allocate(reorder(1,1))
        endif
        do i = 1 , pe_nx
          sdispls(i) = 1 + (min(i,nk)-1)*lni*lnj  ! offset of level I, no greater than nk
          if(pe_mex < nk) then
            rdispls(i) = 1 + (i-1)*lni*lnj  ! this pe will collect pe_nx pieces of data
          else
            rdispls(i) = I      ! this pe only gets one word (because zero hangs)
          endif
        enddo
        call MPI_alltoallv(zin,    sendcounts,sdispls,MPI_INTEGER,
     %                     recvbuf,recvcounts,rdispls,MPI_INTEGER,
     %                     pe_myrow,ierr)
!
!       recvbuf data shape (lni,lnj,pe_nx) on PE( 0:nk , 0:pe_ny-1 )
!
        if(pe_mex < nk) then   ! reorder and collect along Y into zout (first nk PEs only)
          i0 = 1
          do k = 1 , pe_nx     ! time to get full x axis together
            offset = (k-1) * lni*lnj  ! offset in recv buffer for i piece no k, row 1
            ilast = min( gni , i0+lni-1 )
            do j = 1 , lnj      ! copy one slice
              reorder(i0:ilast,j) = recvbuf(offset:offset+ilast-i0)  ! i piece no k, row j
              offset = offset + lnj                                  ! next row
            enddo
            i0 = i0 + lni   ! next piece of i values
          enddo
!
!       reorder data shape (gni,lnj) on PE( 0:nk , 0:pe_ny-1 )
!
          recvcounts        = gni * lnj
          recvcounts(pe_ny) = gni * (gnj - lnj*(pe_ny-1))
          do j = 1 , pe_ny
            rdispls(j) = (j-1) * gni*lnj + 1
          enddo
          call MPI_gatherv(reorder,gni*lnj,            MPI_INTEGER,
     %                     zout,   recvcounts,rdispls, MPI_INTEGER,
     %                     mod(pe_mex,pe_ny),pe_mycol,ierr)      ! root PEs on diagonal
!
!       zout data shape (gni,lnj) on PE(l,m) (l=0:nk-1 , m=  l mod pe_ny)
!
        endif
        goto 1111        ! clean exit
      endif

      if(nk <= pe_ny) then    ! alltoallv along y, then gather along x, then reorder
        sendcounts = 1        ! should be 0 but zero length can hang allgatherv
        sendcounts(1:nk) = lni*lnj
        if(pe_mey < nk) then  ! will be collecting and reordering lni by gnj
          recvcounts = lni*lnj
          allocate(recvbuf(lni*lnj*pe_ny))
          allocate(reorder(lni*lnj*pe_ny,pe_nx))
        else                  ! dummy collect and no reorder
          recvcounts = 1
          allocate(recvbuf(pe_ny))
          allocate(reorder(1,1))
        endif
        call MPI_alltoallv(zin,    sendcounts,sdispls,MPI_INTEGER,
     %                     recvbuf,recvcounts,rdispls,MPI_INTEGER,
     %                     pe_mycol,ierr)
!
!       recvbuf data shape (lni,lnj,pe_ny) on PE( 0:nk-1 ,0:pe_ny-1 )
!
        if(pe_mey < nk) then                ! gather and reorder
          recvcounts = lni*gnj              ! lni by gnj slices
          recvcounts(pe_nx) = gnj * (gni - lni*(pe_nx-1))  ! last PEs lni may be shorter
          do i = 1 , pe_nx
            rdispls(i) = (i-1) *lni*gnj + 1
          enddo
          call MPI_gatherv(recvbuf, lni*gnj,            MPI_INTEGER,
     %                     reorder, recvcounts,rdispls, MPI_INTEGER,
     %                     mod(pe_mey,pe_nx),pe_myrow,ierr)      ! root PEs on diagonal
!
!       reorder data shape (lni,gnj,pe_nx) on PE(m,l) (l=0:nk-1 , m=  l mod pe_nx)
!
         i0 = 1
         do k = 1 , pe_nx      ! time to reorder, pe_nx slices along i
           ilast = min( gni , i0+lni-1 )
           do j = 1 , gnj
             offset = (j-1) * lni
             zout(i0:ilast,j) = reorder(offset:offset+ilast-i0,k)
           enddo
           i0 = i0 + lni
         enddo
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