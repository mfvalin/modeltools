program test
  use rkl_sync_mod
  use mpi_f08
  implicit none
  integer :: npes, rank, sid, i, errors, status
  integer(C_INT), target :: dummyvar
  type(C_PTR) :: p0, p1, p2
  integer, dimension(:), pointer :: array
  integer, parameter :: MEMSIZE = 4096
  real(kind=8) :: t0, t1
  integer :: t
  type(rkl_lock) :: mylock, dummylock
  type(rkl_barrier) :: mybarrier

  call mpi_init()
  call mpi_comm_size(MPI_COMM_WORLD, npes)
  call mpi_comm_rank(MPI_COMM_WORLD, rank)

!   if(rank == 0) then
    call dummylock%init(dummyvar)
    call dummylock%lock(999, 1)
    if(dummylock%owner() .ne. 999) print *,'dummylock%lock : ERROR, not owner of lock when it should be'
  !   if(dummylock%owner() .eq. 999) print *,'dummylock%lock : SUCCESS, owner of lock when it should be'
!     status = TryReleaseIdLock(dummyvar, 999, 1)
!     if(status .ne. 1) print *,'TryReleaseIdLock 1 : ERROR, expected 1, status =',status
    status = dummylock%tryunlock(999, 1)
    if(status .ne. 1) print *,'dummylock%tryunlock 1 : ERROR, expected 1, status =',status
! 
!     call dummylock%lock(999, 1)
    status = dummylock%trylock(999, 1)    ! trylock should succeed
    if(status .ne. 1) print *,'dummylock%trylock : expected 1, status =',status
    if(dummylock%owner() .ne. 999) print *,'dummylock%lock : ERROR, not owner of lock when it should be'
  !   if(dummylock%owner() .eq. 999) print *,'dummylock%lock : SUCCESS, owner of lock when it should be'

    call dummylock%reset()
    if(dummylock%owner() .ne. -1) print *,'dummylock%lock : ERROR, owner is not -1 when it should be'
  !   if(dummylock%owner() .eq. -1) print *,'dummylock%lock : SUCCESS, owner is -1 when it should be'

    call ReleaseIdLock(dummyvar, 999, 1)  ! lock has already been reset, test that this will not deadlock
    status = dummylock%tryunlock(999, 1)  ! tryunlock is expected to fail
    if(status .ne. 0) print *,'dummylock%tryunlock 2 : ERROR, expected 0, status =',status
!     status = TryReleaseIdLock_(C_LOC(dummyvar), 999, 1)
!     if(status .ne. 0) print *,'TryReleaseIdLock_ 2 : ERROR, expected 0, status =',status
!   endif

  p0 = C_NULL_PTR
  sid = -1
  if(rank == 0) then
    p0 = allocate_safe_shared_memory(sid, MEMSIZE)
  endif
  call mpi_bcast(sid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
  if(rank > 0) then
    p0 = allocate_safe_shared_memory(sid, MEMSIZE)
  endif
  if(.not. C_ASSOCIATED(p0)) write(6,*) "PE",rank+1,' of',npes,', sid =',sid,', attached :',C_ASSOCIATED(p0)
  call C_F_POINTER(p0, array, [MEMSIZE/4])
  if(rank == 0) array = 0

  if(mylock%valid()) write(6,*) 'rank =',rank,', valid =',mylock%valid()
  call mpi_barrier(MPI_COMM_WORLD)
  if(rank == 0) call mylock%init(array(257))
  call mpi_barrier(MPI_COMM_WORLD)
  p1 = C_LOC(array(257))
  call mylock%initp(p1)
  if(.not. mylock%valid()) write(6,*) 'rank =',rank,', valid =',mylock%valid()

  if(rank == 0) array(512) = 0
  call mpi_barrier(MPI_COMM_WORLD)
  t0 = mpi_wtime()
  do i = 1, 100
    call mylock%lock(rank, 1)
    array(512) = array(512) + 1
    call mylock%unlock(rank, 1)
  enddo
  t1 = mpi_wtime() - t0
  t1 = t1 / 100
  t = (t1 * 1.0E+9) + .5
  call mpi_barrier(MPI_COMM_WORLD)
  write(6,*) 'rank =',rank,', value =', array(512),', tlc fenced =',t
  call flush(6)
  call mpi_barrier(MPI_COMM_WORLD)

  if(rank == 0) array(512) = 0
  call mpi_barrier(MPI_COMM_WORLD)
  t0 = mpi_wtime()
  do i = 1, 100
    call mylock%lock(rank, 0)
    array(512) = array(512) + 1
    call mylock%unlock(rank, 0)
  enddo
  t1 = mpi_wtime() - t0
  call mylock%lock(rank, 0)
  if(mylock%owner() .ne. rank) print *,'mylock%owner : ERROR, not owner of lock when it should be'
!   if(mylock%owner() .eq. rank) print *,'mylock%owner : SUCCESS, owner of lock when it should be'
  call mylock%unlock(rank, 0)
  t1 = t1 / 100
  t = (t1 * 1.0E+9) + .5
  call mpi_barrier(MPI_COMM_WORLD)
  write(6,*) 'rank =',rank,', value =', array(512),', tlc =',t
  call mpi_barrier(MPI_COMM_WORLD)

  if(rank == 0) array(512) = 0
  call mpi_barrier(MPI_COMM_WORLD)
  t0 = mpi_wtime()
  do i = 1, 100
    call AcquireIdLock(array(256), rank, 1)
    array(512) = array(512) + 1
    call ReleaseIdLock(array(256), rank, 1)
  enddo
  t1 = mpi_wtime() - t0
  call AcquireIdLock(array(256), rank, 1)
  if(rank .ne. LockOwner(array(256))) print *,'LockOwner : ERROR, not owner of lock when it should be'
!   if(rank .eq. LockOwner(array(256))) print *,'LockOwner : SUCCESS, owner of lock when it should be'
  call ReleaseIdLock(array(256), rank, 1)
  t1 = t1 / 100
  t = (t1 * 1.0E+9) + .5
  call mpi_barrier(MPI_COMM_WORLD)
  write(6,*) 'rank =',rank,', value =', array(512),', tl fenced =',t
  call flush(6)
  call mpi_barrier(MPI_COMM_WORLD)

  array = 0
  if(rank == 0) array(512) = 0
  call mpi_barrier(MPI_COMM_WORLD)
  t0 = mpi_wtime()
  do i = 1, 100
    call AcquireIdLock(array(256), rank, 0)
    array(512) = array(512) + 1
    call ReleaseIdLock(array(256), rank, 0)
  enddo
  t1 = mpi_wtime() - t0
  t1 = t1 / 100
  t = (t1 * 1.0E+9) + .5
  call mpi_barrier(MPI_COMM_WORLD)
  write(6,*) 'rank =',rank,', value =', array(512),', tl =',t
  call flush(6)
  call mpi_barrier(MPI_COMM_WORLD)

  if(rank == 0) array(512) = -1
  call mpi_barrier(MPI_COMM_WORLD)
  errors = 0
  t0 = mpi_wtime()
  do i = 1, 100
    if(rank == 0) array(512) = i
    call BasicNodeBarrier(array(1), array(128), npes)
    if(array(512) .ne. i) errors = errors + 1
    call BasicNodeBarrier(array(1), array(128), npes)
  enddo
  t1 = mpi_wtime() - t0
  t1 = t1 / 200
  t = (t1 * 1.0E+9) + .5
  call mpi_barrier(MPI_COMM_WORLD)
  write(6,*) 'rank =',rank,', errors =',errors,', tbarr1 =',t
  call flush(6)
  call mpi_barrier(MPI_COMM_WORLD)

  call mybarrier%init(array(1), array(128))
  p1 = C_LOC(array(2))
  p2 = C_LOC(array(129))
  call mybarrier%initp( p1,  p2)
  if(rank == 0) array(512) = -1
  call mpi_barrier(MPI_COMM_WORLD)
  errors = 0
  t0 = mpi_wtime()
  do i = 1, 100
    if(rank == 0) array(512) = i
    call mybarrier%set(npes)
    if(array(512) .ne. i) errors = errors + 1
    call mybarrier%set(npes)
  enddo
  t1 = mpi_wtime() - t0
  t1 = t1 / 200
  t = (t1 * 1.0E+9) + .5
  call mpi_barrier(MPI_COMM_WORLD)
  write(6,*) 'rank =',rank,', errors =',errors,', tbarrMPI =',t
  call flush(6)
  
  call mpi_finalize()
end program

