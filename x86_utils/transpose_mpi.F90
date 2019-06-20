subroutine transpose_64_mpi(comm, npe, pe_me, forward, za, zb, mini1, maxi1, n1g, ln1max, ln2, n3g, ln3max)
  use ISO_FORTRAN_ENV
  implicit none
  include 'mpif.h'
  integer, intent(IN) :: comm, npe, pe_me
  integer, intent(IN) :: mini1, maxi1, n1g, ln1max, ln2, n3g, ln3max
  logical, intent(IN) :: forward
  real(kind=REAL64), dimension(mini1:maxi1,ln2*n3g), intent(INOUT) :: za  ! dimensions 2 and 3 are collapsed
  real(kind=REAL64), dimension(ln2*ln3max,n1g), intent(INOUT) :: zb       ! dimensions 1 and 2 are collapsed

  real(kind=kind(1.d0)), dimension(ln2*ln3max,ln1max,npe) :: ta               ! temporary array, ready for alltoallv
  integer, dimension(npe) :: send_cnt, send_dsp, recv_cnt, recv_dsp
  integer :: base1, base3, p, ln3, mybase1, myln1, mybase3, myln3, ln1
  integer :: ierr !, i, j, k

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
    do  p = 1, npe
      base3 = (p - 1) * ln3max + 1                     ! first point to send to PE p  along third dimension (za array)
      ln3 = max( 0 , min(ln3max, n3g - base3 + 1) )    ! number of points to send to PE p along the third dimension
      if(ln3 > 0) then
        call transpose_64(za(1,1+(base3-1)*ln2), maxi1-mini1+1, myln1, ta(1,1,p), ln2*ln3max, ln2*ln3)
        send_cnt(p) = myln1 * ln2 * ln3                ! number of points to send to PE p
        send_dsp(p) = myln1 * ln2 * (p - 1) * ln3max
      endif
      base1 = (p - 1) * ln1max + 1                     ! first point to receive from PE p  along the first dimension
      ln1 = min(ln1max, n1g - base1 +1 )               ! number of points to receive from PE p along the first dimension
      if(ln1 > 0) then
        recv_cnt(p) = ln1 * ln2 * myln3                ! number of points to receive from PE p
        recv_dsp(p) = (p - 1) * ln1max * ln2 * ln3max
      endif
    enddo
  write(6,*)'------------------------'
  write(6,101)send_cnt,send_dsp,recv_cnt,recv_dsp
  write(6,*)'------------------------'
  call flush(6)
  call MPI_barrier(COMM,ierr)
    call MPI_alltoallv(ta, send_cnt, send_dsp, MPI_REAL8, zb, recv_cnt, recv_dsp, MPI_REAL8, COMM, ierr)
  else                     ! inverse transpose
    do  p = 1, npe
      base3 = (p - 1) * ln3max + 1                     ! first point to send to PE p  along third dimension (za array)
      ln3 = max( 0 , min(ln3max, n3g - base3 + 1) )    ! number of points to send to PE p along the third dimension
      if(ln3 > 0) then
        recv_cnt(p) = myln1 * ln2 * ln3                ! number of points to send to PE p
        recv_dsp(p) = myln1 * ln2 * (p - 1) * ln3max
      endif
      base1 = (p - 1) * ln1max + 1                     ! first point to receive from PE p  along the first dimension
      ln1 = min(ln1max, n1g - base1 +1 )               ! number of points to receive from PE p along the first dimension
      if(ln1 > 0) then
        send_cnt(p) = myln1 * ln2 * ln3                ! number of points to receive from PE p
        send_dsp(p) = (p - 1) * ln1max * ln2 * ln3
      endif
    enddo
    call MPI_alltoallv(zb, send_cnt, send_dsp, MPI_REAL8, ta, recv_cnt, recv_dsp, MPI_REAL8, COMM, ierr)
    do  p = 1, npe
      base3 = (p - 1) * ln3max + 1
      ln3 = max( 0 , min(ln3max, n3g - base3 + 1) )
      if(ln3 > 0) then
        call transpose_64(ta(1,1,p), ln2*ln3max, ln2*ln3, za(1,1+(base3-1)*ln2), maxi1-mini1+1, myln1)
      endif
    enddo
  endif
!   do p = npe,1,-1
!   do i = ln1max,1,-1
!     print 100,1.0_8*i,1.0_8*p,ta(:,i,p)
!   enddo
!   enddo
!   print *,'------------------------'
100 format(20F8.0)
101 format(20I8)
  return
contains
  subroutine transpose_64(src,ls,ns,dst,ld,nd)
    use ISO_FORTRAN_ENV
    implicit none
    integer, intent(IN) :: ls, ns, ld, nd
    real(kind=REAL64), dimension(ls,nd), intent(INOUT) :: src
    real(kind=REAL64), dimension(ld,ns), intent(INOUT) :: dst
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
!     do j = nd,1,-1
!       print 100,1.0_8*j,src(1:ns,j)
!     enddo
!     print *,'+++++++++++++++++++++++'
!     do i = ns,1,-1
!       print 100,1.0_8*i,dst(1:nd,i)
!     enddo
!     print *,'+++++++++++++++++++++++'
    return
100 format(20F8.0)
  end subroutine transpose_64
end
#if defined(SELF_TEST)
program test
  use ISO_FORTRAN_ENV
  implicit none
  integer, parameter :: LNI = 3
  integer, parameter :: LNJ = 4
  integer, parameter :: NK = 3
  include 'mpif.h'

  real(kind=REAL64), dimension(LNI,LNJ,NK)  :: za
  real(kind=REAL64), dimension(LNI,LNJ,NK*2)  :: zc
  real(kind=REAL64), dimension(:,:,:), pointer :: zb
  integer :: i, j, k, ierr, rank, npes, lnkmax

  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_comm_size(MPI_COMM_WORLD,npes,ierr)

  lnkmax = (NK + npes -1)/npes
  allocate(zb(LNJ,lnkmax,LNI*npes))
print *,'size of zb =',LNJ,lnkmax,LNI*npes
  do k = 1,NK
  do j = 1,LNJ
  do i = 1,LNI
    za(i,j,k) = 10000*(i+LNI*rank) + 100*j + k
  enddo
  enddo
  enddo
  do k = NK,1,-1
  do j = LNJ,1,-1
    print 101,j*1.0_8,k*1.0_8,za(:,j,k)
  enddo
  enddo
  print *,'===================================='

  call transpose_64_mpi(MPI_COMM_WORLD, npes, rank, .true., za, zb, 1, LNI, npes*LNI, LNI, LNJ, NK, lnkmax)
  do i = 1,LNI*2
    print 101,i*1.0_8,zb(:,1,i)
  enddo
  print *,'===================================='
goto 1
  call transpose_64_mpi(MPI_COMM_WORLD, npes, rank, .false., zc, zb, 1, LNI, npes*LNI, LNI, LNJ, NK, lnkmax)
  do k = NK,1,-1
  do j = LNJ,1,-1
    print 101,j*1.0_8,k*1.0_8,zc(:,j,k)
  enddo
  enddo
  print *,'===================================='

1 call MPI_finalize(ierr)
  stop
100 format(20I8)
101 format(20F8.0)
end program test
#endif
