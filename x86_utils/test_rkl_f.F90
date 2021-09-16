program test
  use rkl_sync_mod
  use mpi_f08
  implicit none
  integer :: npes, rank, sid, i, errors
  type(C_PTR) :: p0, p1, p2
  integer, dimension(:), pointer :: array
  integer, parameter :: MEMSIZE = 4096
  real(kind=8) :: t0, t1
  integer :: t
  type(rkl_lock) :: mylock
  type(rkl_barrier) :: mybarrier

  call mpi_init()
  call mpi_comm_size(MPI_COMM_WORLD, npes)
  call mpi_comm_rank(MPI_COMM_WORLD, rank)

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

