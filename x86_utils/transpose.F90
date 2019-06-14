subroutine transpose_m32(src,ls,ns,dst,ld,nd)
  implicit none
  integer, intent(IN) :: ls, ns, ld, nd
  integer, dimension(ls,nd), intent(IN) :: src
  integer, dimension(ld,ns), intent(OUT) :: dst
  integer :: i, j, i0, j0

  do i0 = 1, ns, 8
    do j0 = 1, nd, 8
      do i = i0, min(i0+7,ns)
      do j = j0, min(j0+7,nd)
        dst(j,i) = src(i,j)
      enddo
      enddo
    enddo
  enddo
end subroutine transpose_m32

subroutine transpose_m64(src,ls,ns,dst,ld,nd)
  implicit none
  integer, intent(IN) :: ls, ns, ld, nd
  integer(kind=8), dimension(ls,nd), intent(IN) :: src
  integer(kind=8), dimension(ld,ns), intent(OUT) :: dst
  integer :: i, j, i0, j0

  do i0 = 1, ns, 8
    do j0 = 1, nd, 8
      do i = i0, min(i0+7,ns)
      do j = j0, min(j0+7,nd)
        dst(j,i) = src(i,j)
      enddo
      enddo
    enddo
  enddo
end subroutine transpose_m64
#if defined(SELF_TEST)
#if defined(WITH_MPI)
program test
  implicit none
  include 'mpif.h'
  integer, parameter :: LI = 7
  integer, parameter :: LJ = 5
  integer, parameter :: LK = 4
  integer :: ierr, rank, npes
  integer :: lni, lnj

  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,npes,ierr)
  lnj = LNJ
  lni = (LI+npes-1)/npes
  lni = min(lni,LI-lni*rank)
  print *,'PE',rank+1,' of',npes,' lni =',lni
  call mpi_finalize(ierr)
end program test
#else
program test
  implicit none
  integer, parameter :: LI = 170
  integer, parameter :: LJ = 150
  integer, parameter :: LK = 81
  integer, dimension(LI,LJ) :: src
  integer, dimension(LJ,LI) :: dst
  integer(kind=8), dimension(LI,LJ) :: src8
  integer(kind=8), dimension(LJ,LI) :: dst8
  integer, dimension(-1:LI+1,LJ,LK) :: za, zc
  integer, dimension(LJ,LK,LI) :: zb
  integer :: i, j, k, m, ni, nj
  logical :: ok

  k = 1
  ni = 7
  nj = 5
  src = 0
  dst = -1
  src8 = 0
  dst8 = -1

  do j = 1, nj
  do i = 1, ni
    src(i,j) = k
    src8(i,j) = k
    k = k + 1
  enddo
  enddo

  call transpose_m32(src,LI,ni,dst,LJ,nj)
  do j = 1, LJ
    if(LI*LJ < 100) print 101,src(1:7,j),dst(1:7,j)
  enddo
  print *,'==================================================================='
  call transpose_m64(src8,LI,ni,dst8,LJ,nj)
  do j = 1, LJ
    if(LI*LJ < 100) print 101,src8(1:7,j),dst8(1:7,j)
  enddo
  print *,'==================================================================='

  zb = -1
  za = 0
  m = 1
  do k = 1, LK
  do j = 1, LJ
  do i = 1, LI
    za(i,j,k) = m
    m = m + 1
  enddo
  enddo
  enddo
  do k = LK, 1, -1
    do j = LJ, 1, -1
      if(LI*LJ < 100) print 102,za(:,j,k)
    enddo
    if(LI*LJ < 100) print *,''
  enddo
  print *,'==================================================================='
  call transpose_m32(za(1,1,1),LI+3,LI,zb,LK*LJ,LJ*LK)
  do k = LI, 1, -1
    do j = LK, 1, -1
      if(LI*LJ < 100) print 102,zb(:,j,k)
    enddo
    if(LI*LJ < 100) print *,''
  enddo
  print *,'==================================================================='
  zc = 0
  call transpose_m32(zb,LK*LJ,LJ*LK,zc(1,1,1),LI+3,LI)
  do k = LK, 1, -1
    do j = LJ, 1, -1
      if(LI*LJ < 100) print 102,zc(:,j,k)
    enddo
    if(LI*LJ < 100) print *,''
  enddo
  ok = any(za .ne. zc)
  print *,'errors =',ok

101 format(7I4,4x,7I4)
102 format(30I4)
end program
#endif
#endif
