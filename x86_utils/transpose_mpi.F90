!    useful routines for C and FORTRAN programming
!    Copyright (C) 2019  Division de Recherche en Prevision Numerique
!   
!    This is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!   
!    This software is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
subroutine transpose_64_mpi(comm, npe, pe_me, forward, za, zb, mini1, maxi1, n1g, ln1max, ln2, n3g, ln3max) ! MPI transpose
! comm       : MPI communicator
! npe        : number of PEs in communicator
! pe_me      : rank of this PE in communicator
! forward    : if .true. transpose za -> zb, if .false. reverse transpose zb -> za
! za, zb     : source/destination of transpose
!              za(mini1:maxi1, ln2, n3g) NOTE: dimensions 2 and 3 are collapsed in the code
!              zb(ln2, ln3max, n1g)      NOTE: dimensions 1 and 2 are collapsed in the code
! mini1      : lower bound of first subscript of local array za
! maxi1      : upper bound of first subscript of local array za
! n1g        : GLOBAL first dimension of za (all PEs in comm)
! ln1max     : largest value of local useful first dimension of local za ( <= maxi1)
! ln2        : second dimension of za, first dimension of zb
! n3g        : full third dimension of array za
! ln3max     : largest value of local useful second dimension of zb (after distribution over npe PEs)
!
! NOTE: the "forward" transpose from za -> zb is done with subscript rotation. 
!       the "gathered" first subscript of za ends up as the last subscript of zb
!       za(local1 , local2, global3) -> zb(local2, local3, global1)
!       the "inverse" transpose from zb -> za does the opposite.
!       the "redistributed" last subscript of zb ends up as the first subscript of za
!       zb(local2, local3, global1) -> za(local1 , local2, global3)
!
!       it may happen during a "forward" transpose that some of the last PEs may get NOTHING, if the third 
!       subscript of za is not large enough
!       ex.  n3g = 80, npe = 36 : 
!            the first 26 PEs will get 3 levels, one PE will get 2 levels , the last 9 PEs will get NOTHING (26*3+2 = 80)
! (PE is a synonym for process in the comments)
!
  use ISO_FORTRAN_ENV
  implicit none
  include 'mpif.h'
  integer, intent(IN) :: comm, npe, pe_me
  integer, intent(IN) :: mini1, maxi1, n1g, ln1max, ln2, n3g, ln3max
  logical, intent(IN) :: forward
  real(kind=REAL64), dimension(mini1:maxi1, ln2*n3g), intent(INOUT) :: za  ! dimensions 2 and 3 are collapsed
  real(kind=REAL64), dimension(ln2*ln3max, n1g), intent(INOUT)      :: zb  ! dimensions 1 and 2 are collapsed

  real(kind=kind(1.d0)), dimension(ln2*ln3max, ln1max, npe)         :: ta  ! temporary array, ready for alltoallv
  integer, dimension(npe) :: send_cnt, send_dsp, recv_cnt, recv_dsp
  integer :: base3, p, ln3, mybase1, myln1, mybase3, myln3, ln1
  integer :: ierr! , i , j, k
  real(kind=REAL64) :: t0,t1, t2, t3

  send_cnt = 0
  recv_cnt = 0
  send_dsp = 0
  recv_dsp = 0
  ! for THIS PE
  mybase1 = 1 + (pe_me * ln1max)                       ! first point along the first dimension (za array)
  myln1   = min(ln1max, n1g - mybase1 + 1)             ! number of points along first dimension (za array)
  mybase3 = 1 + (pe_me * ln3max)                       ! first point along the third dimension (zb/ta arrays)
  myln3   = max( 0 , min(ln3max, n3g - mybase3 + 1) )  ! number of points along the third dimension (zb/ta arrays)
  if(forward) then         ! forward transpose
    t1 = MPI_wtime()
    do  p = 1, npe
      base3 = (p - 1) * ln3max                         ! first point to send to PE p  along third dimension (za array)
      ln3 = max( 0 , min(ln3max, n3g - base3) )        ! number of points to send to PE p along the third dimension
      if(ln3 > 0) then
        call transpose_64(za(1,1+base3*ln2), maxi1-mini1+1, myln1, ta(1,1,p), ln2*ln3max, ln2*ln3)
        send_cnt(p) = myln1 * ln2 * ln3max             ! number of points to send to PE p
        send_dsp(p) = ln1max * ln2 * (p - 1) * ln3max  ! offset for points to send to PE p
      endif
      ln1 = min(ln1max, n1g - (p - 1) * ln1max)        ! number of points to receive from PE p along the first dimension
      if(myln3 > 0) then
        recv_cnt(p) = ln1 * ln2 * ln3max               ! number of points to receive from PE p
        recv_dsp(p) = (p - 1) * ln1max * ln2 * ln3max  ! offset for points to receive from PE p
      endif
    enddo
    t2 = MPI_wtime()
    call MPI_alltoallv(ta, send_cnt, send_dsp, MPI_REAL8, zb, recv_cnt, recv_dsp, MPI_REAL8, COMM, ierr)
    t3 = MPI_wtime()

  else                     ! inverse transpose
    t0 = MPI_wtime()
    do  p = 1, npe
      base3 = (p - 1) * ln3max
      ln3 = max( 0 , min(ln3max, n3g - base3) )        ! number of points to receive from PE p along the third dimension
      if(ln3 > 0) then
        recv_cnt(p) = myln1 * ln2 * ln3max             ! number of points to receive from PE p
        recv_dsp(p) = ln1max * ln2 * (p - 1) * ln3max  ! offset for points to receive from PE p
      endif
      ln1 = min(ln1max, n1g - (p - 1) * ln1max)        ! number of points to send to PE p along the first dimension
      if(myln3 > 0) then
        send_cnt(p) = ln1 * ln2 * ln3max               ! number of points to send to PE p
        send_dsp(p) = (p - 1) * ln1max * ln2 * ln3max  ! offset for points to send to PE p
      endif
    enddo
    t1 = MPI_wtime()
    call MPI_alltoallv(zb, send_cnt, send_dsp, MPI_REAL8, ta, recv_cnt, recv_dsp, MPI_REAL8, COMM, ierr)
    t2 = MPI_wtime()
    do  p = 1, npe
      base3 = (p - 1) * ln3max 
      ln3 = max( 0 , min(ln3max, n3g - base3) )
      if(ln3 > 0) then
        call transpose_64(ta(1,1,p), ln2*ln3max, ln2*ln3, za(1,1+base3*ln2), maxi1-mini1+1, myln1)
      endif
    enddo
    t3 = MPI_wtime()
  endif
#if defined(SELF_TEST)
! write(6,101)send_cnt, send_dsp, recv_cnt, recv_dsp
  if(pe_me == 0) then
!     if(forward)       print '(A,2G12.3)','forward transpose timing TR AA =',t2-t1,t3-t2
!     if(.not. forward) print '(A,2G12.3)','inverse transpose timing TR AA =',(t3-t2)+(t1-t0),t2-t1
  endif
101 format(20I4)
#endif
  return
contains
  subroutine transpose_64(src,ls,ns,dst,ld,nd)
    use ISO_FORTRAN_ENV
    implicit none
    integer, intent(IN) :: ls, ns, ld, nd
    real(kind=REAL64), dimension(ls,nd), intent(INOUT) :: src
    real(kind=REAL64), dimension(ld,ns), intent(INOUT) :: dst
    integer :: i, j, i0, j0
    integer, parameter :: BLOCK=48   ! heuristically determined

    do i0 = 1, ns, BLOCK
      do j0 = 1, nd, BLOCK
        do i = i0, min(i0+BLOCK-1,ns)  ! BLOCK x BLOCK cache blocking
        do j = j0, min(j0+BLOCK-1,nd)
          dst(j,i) = src(i,j)
        enddo
        enddo
      enddo
    enddo
    return
! 100 format(20F8.0)
  end subroutine transpose_64
end
#if defined(SELF_TEST)
program test    ! calling sequence : ./a.out nrows mni lnj nk (arguments are positional)
! test program for transpose_64_mpi
! the program performs nrows simultaneous transposes along nclos processes
! (ncols x nrows) processes are needed
! (PE is a synonym for process in the comments)
  use ISO_FORTRAN_ENV
  implicit none
  integer :: MNI    = 26    ! maximum size along the first dimension of the local tile
  integer :: LNJ    = 22    ! maximum size along the second dimension of the local tile
  integer :: NK     = 80    ! maximum size along the third dimension of the local tile
  integer :: nrows  = 1     ! npey, number of rows of PEs (3 rows with 36 PEs mean a 8x3 MPI configuration)
  integer :: pernode= 36    ! assuming 36 processes per node
  integer :: blockx = 2
  integer :: blocky = 2
  include 'mpif.h'

  character(len=128) :: arg
  integer, dimension(:,:), allocatable :: map
  real(kind=REAL64), dimension(:,:,:), allocatable  :: za  ! source array
  real(kind=REAL64), dimension(:,:,:), allocatable  :: zb  ! transposed array
  real(kind=REAL64), dimension(:,:,:), allocatable  :: zc  ! back transposed array
  real(kind=REAL64) :: t1, t2, t3
  integer :: i, j, k, ierr, rank, npes, lnkmax, mylnk, lni, gni, gnj, errors, errtot, lnkbase, err2, errtot2, i0, j0
  integer, dimension(288) :: err2ss
  integer, dimension(:), pointer :: err2s  ! this does not work on swan, there is something fishy
  integer :: my_row, pe_mey, nrow, my_col, pe_mex, ncols, npts, irep, nodes
  integer(kind=8) :: netvol

  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,npes,ierr)

  call GET_COMMAND_ARGUMENT(1, arg)
  if(arg .ne. " ") then
    if( arg(1:2) == "-h") then
      if(rank == 0) print *,"program nrows(1) mni(26) lnj(22) nk(80) pernode(36) blockx(2) blocky(2)"
      goto 1
    endif
    read(arg,*)nrows
  endif
  ncols = (npes+nrows-1)/nrows       ! number of columns given number of PEs and number of rows
  if(nrows*ncols .ne. npes) then
    if(rank == 0) print 102,"ERROR: nrows*ncols .ne. npes, nrows, ncols, npes",nrows, ncols, npes
    goto 1
  endif

  allocate(map(ncols,nrows))

  call GET_COMMAND_ARGUMENT(2, arg)
  if(arg .ne. " ") read(arg,*)MNI

  call GET_COMMAND_ARGUMENT(3, arg)
  if(arg .ne. " ") read(arg,*)LNJ

  call GET_COMMAND_ARGUMENT(4, arg)
  if(arg .ne. " ") read(arg,*)NK

  call GET_COMMAND_ARGUMENT(5, arg)
  if(arg .ne. " ") read(arg,*)pernode

  call GET_COMMAND_ARGUMENT(6, arg)
  if(arg .ne. " ") read(arg,*)blockx

  call GET_COMMAND_ARGUMENT(7, arg)
  if(arg .ne. " ") read(arg,*)blocky

  npts = MNI*LNJ*NK*npes
  nodes = (npes + pernode -1) / pernode
  netvol = npts * 2
  netvol = netvol * (nodes-1) / nodes
  if(rank == 0) then
    print 102,' nrows, ncols, nodes, MNI, LNJ, NK =',nrows, ncols, nodes, MNI, LNJ, NK
    print *,'NPTS =',npts
   endif

  k = 0
  do j0 = 1, nrows-blocky+1, blocky
  do i0 = 1, ncols-blockx+1, blockx
  do j = j0, j0+blocky-1
  do i = i0, i0+blockx-1
   if(k == rank) then
     pe_mex = i - 1
     pe_mey = j - 1
   endif
   map(i,j) = k
   k = k + 1
  enddo
  enddo
  enddo
  enddo
  if(mod(ncols,blockx) .ne. 0 .or. mod(nrows,blocky) .ne. 0) then
    if(rank == 0) print 102,"ERROR: bad block distribution, nrows, ncols, npes, blockx, blocky = ",nrows, ncols, npes, blockx, blocky
    goto 1
  endif
  if(rank == 0) then
    do j = nrows, 1, -1
      print 100,map(:,j)
    enddo
  endif
  ! split into rows
!   pe_mey = mod(rank,nrows)
  call MPI_comm_split(MPI_COMM_WORLD, pe_mey, rank, my_row, ierr)
  call MPI_comm_split(MPI_COMM_WORLD, pe_mex, rank, my_col, ierr)
!   call MPI_comm_rank(my_row, pe_mex, ierr)
  call MPI_comm_size(my_row, nrow, ierr)
  allocate(err2s(nrow))
! goto 1
  lnkmax = (NK + nrow -1)/nrow                       ! max number of levels after transpose (all processes)
  mylnk = max( min(lnkmax , NK - lnkmax * pe_mex) , 0) ! number of levels after transpose for THIS process
  lnkbase = lnkmax * pe_mex + 1                        ! first k level for this PE
  gni = MNI*ncols - 2                                ! global first dimension (last tile deliberately 2 points short)
  lni = min(MNI, gni - (MNI*pe_mex))                 ! first dimension for THIS tile
  gnj = LNJ*nrows
  allocate(zb(0:LNJ+1, lnkmax, gni))    ! allocate arrays
  allocate(za(0:MNI+1,0:LNJ+1,NK))
  allocate(zc(0:MNI+1,0:LNJ+1,NK))
  zb = -1.0                         ! initialize arrays
  za = 0.0
  zc = 0.0
! print *,'gni, lni, lnkmax, mylnk =',gni,lni,lnkmax,mylnk

  do k = 1,NK
  do j = 1,LNJ
  do i = 1,lni
    za(i,j,k) = (10000*(i+pe_mex*MNI) + 100*j + k) 
  enddo
  enddo
  enddo
#if defined(DIAG)
  if(pe_mey == 0) then
    print *,'=============== 1 =================='
    do k = NK,1,-1
      print *,'k = ',k
      do j = LNJ,1,-1
        print 101,j*1.0_8,za(:,j,k)
      enddo
    enddo
  endif
#endif
  ! 3 rounds to "prime the pump"
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call transpose_64_mpi(my_row, nrow, pe_mex, .true., za, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call transpose_64_mpi(my_row, nrow, pe_mex, .false., zc, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)

  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call transpose_64_mpi(my_row, nrow, pe_mex, .true., za, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call transpose_64_mpi(my_row, nrow, pe_mex, .false., zc, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)

  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call transpose_64_mpi(my_row, nrow, pe_mex, .true., za, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call transpose_64_mpi(my_row, nrow, pe_mex, .false., zc, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)

  
  do irep = 1, 5
    call MPI_barrier(MPI_COMM_WORLD,ierr)
    if(pe_mey+pe_mex == 0) print *
    call MPI_barrier(my_row,ierr)
    t1 = MPI_wtime()*1000000.
    call transpose_64_mpi(my_row, nrow, pe_mex, .true., za, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)
#if defined(DIAG)
    if(pe_mey == 0 .and. irep == 1) then
      print *,'=============== 2 =================='
      do k = mylnk,1,-1
        print *,'k = ',k
        do i = gni,1,-1
          print 101,i*1.0_8,zb(:,k,i)
        enddo
      enddo
    endif
#endif
    call MPI_barrier(my_row,ierr)
    t2 = MPI_wtime()*1000000.
    call transpose_64_mpi(my_row, nrow, pe_mex, .false., zc, zb, 0, MNI+1, gni, MNI, LNJ+2, NK, lnkmax)
    call MPI_barrier(my_row,ierr)
    t3 = MPI_wtime()*1000000.
    if(pe_mex == 0) print 102,'row, fwd, inv time =',pe_mey,nint(t2-t1), nint(t3-t2)
  enddo

  call MPI_barrier(MPI_COMM_WORLD,ierr)
  if(rank == 0) then
    print 102,'global grid =',gni,LNJ*nrows,NK
    print '(A,F10.0)','MBytes/sec =',8*(npts*4.0/(t3-t1))
    print '(A,2F10.0)','network MBytes/sec =',8*(netvol/(t3-t1)),8*(netvol/(t3-t1))/max(1,nodes-1)
  endif
  call MPI_barrier(MPI_COMM_WORLD,ierr)

#if defined(DIAG)
  if(pe_mey == 0) then
    print *,'=============== 3 =================='
    do k = NK,1,-1
      print *,'k = ',k
      do j = LNJ,1,-1
        print 101,j*1.0_8,zc(:,j,k)
      enddo
    enddo
  endif
#endif
  err2 = 0
  do i = 1,gni
  do k = 1, mylnk
  do j = 1,LNJ
    if(zb(j,k,i) .ne. (10000*i + 100*j + k + lnkbase - 1)) err2 = err2 + 1
  enddo
  enddo
  enddo
  if(err2 .ne. 0) print 102,'pe_mex, pe_mey,mylnk,err2 =',pe_mex, pe_mey,mylnk,err2
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  err2s = 0
  err2ss = 0
!   if(pe_mex == 0) print *,'before gather',my_row,nrow
  call MPI_gather(err2, 1, MPI_INTEGER, err2ss, nrow, MPI_INTEGER, 0, my_row, ierr)
  err2s(1:nrow) = err2ss(1:nrow)
!   if(pe_mex == 0) print *,'after gather',my_row
  call MPI_barrier(MPI_COMM_WORLD,ierr)
!   if(pe_mex == 0) print 102,'pe_mey, pts, err2s =',pe_mey,gni*mylnk*LNJ,err2s(1:nrow)
  call MPI_reduce(err2, errtot2, 1, MPI_INTEGER, MPI_SUM, 0, my_row, ierr)
!   if(pe_mex == 0) print *,'after first reduce',my_row
  errors = 0
  do k = 1,NK
  do j = 1,LNJ
  do i = 1,lni
    if(zc(i,j,k) .ne. za(i,j,k)) errors = errors + 1
  enddo
  enddo
  enddo
!   if(pe_mex == 0) print *,'before second reduce',my_row
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  call MPI_reduce(errors, errtot, 1, MPI_INTEGER, MPI_SUM, 0, my_row,ierr)
!   if(pe_mex == 0) print *,'after second reduce',my_row
  call MPI_barrier(MPI_COMM_WORLD,ierr)
!   if(pe_mex == 0) print 102,'row, fwd, inv time =',pe_mey,nint(t2-t1), nint(t3-t2)h
  call MPI_barrier(MPI_COMM_WORLD,ierr)
  if(pe_mex == 0) print 102,'row, number of errors =',pe_mey,errtot, errtot2

  deallocate(za,zb,zc)
  deallocate(err2s)

1 call MPI_finalize(ierr)
  stop
100 format(20I8)
101 format(20F8.0)
102 format(A,20I8)
end program test
#endif
