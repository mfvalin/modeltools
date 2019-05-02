! Copyright (C) 2019 Recherche en Prevision Numerique
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation, 
! version 2 of the License.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public
! License along with this program; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.
program test_memory_arena
  use ISO_C_BINDING
  implicit none
  integer, parameter :: NSYM=128
  integer, parameter :: DBLK=20
  include 'memory_arena.inc'
  include 'circular_buffer.inc'
  include 'mpif.h'
  interface
    function gethostid() result(id) BIND(C,name='gethostid')
      import :: C_LONG
      integer(C_LONG) :: id
    end function gethostid
  end interface
  integer :: err, rank, isiz, id, myhost, myhost0, myhost1, tempcomm, MY_World
  type(C_PTR) :: shmaddr, p
  integer(C_INT) :: shmsz, shmid
  character(len=128) :: command, myblock
  integer :: bsz, flags, i, errors, j
  integer, dimension(:), pointer :: fp

  call mpi_init(err)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,err)
  call mpi_comm_size(MPI_COMM_WORLD,isiz,err)

  myhost  = gethostid()                    ! host id
  myhost0 = iand(myhost , Z'7FFFFFFF')     ! lower 31 bits
  myhost1 = iand( ishft(myhost, -31) , 1 ) ! upper bit
  call MPI_Comm_split(MPI_COMM_WORLD,myhost0,rank,tempcomm,err)   ! split WORLD using the lower 31 bits of host id , weight=rank in base
  call MPI_Comm_split(tempcomm      ,myhost1,rank,MY_World,err)   ! re split using the upper bit of host id , weight=rank in base
  call MPI_Comm_rank(MY_World,rank,err);                          ! rank of this PE on this SMP node
  call MPI_Comm_size(MY_World,isiz,err);                          ! number of PEs on this SMP node

  id = memory_arena_set_id(rank) ! send my rank as id
  if(rank == 0) then             ! PE0 creates arena and some blocks
    shmsz = 1024 * 1024 * 4
    shmaddr = memory_arena_create_shared(shmid, NSYM, shmsz)
    do id = 1, isiz
      write(myblock,100)'BLOCK',id-1
100   format(A5,I3.3)
      p = memory_block_create(shmaddr, DBLK*id, trim(myblock))
    enddo
    write(command,*)'ipcs -m -u -i ',shmid
    call system(trim(command))       ! list shared memory blocks on system
    call memory_arena_print_status(shmaddr)  ! print arena metadata
  endif

  call MPI_Bcast(shmid, 1, MPI_INTEGER, 0, MY_World, err)          ! broadcast id of shared memory segment
  if(rank .ne. 0) shmaddr = memory_arena_from_id(shmid)            ! everybody but PE0 gets the segment address

  write(myblock,100)'BLOCK',rank
  p = memory_block_find(shmaddr, bsz, flags, trim(myblock))        ! get MY block
  call c_f_pointer(p, fp, [bsz])                                   ! make Fortran pointer
  do i = 1, bsz                    ! fill array
    fp(i) = i + ishft(rank,24)
  enddo
  write(0,*) trim(myblock), ' values =',bsz
  p = memory_block_mark_init(shmaddr, trim(myblock))               ! mark block as initialized
  call MPI_Barrier(MY_World, err)

  if(rank == isiz -1 ) then        ! last PE prints  arena metadata
    write(command,*)'ipcs -m -u -i ',shmid
    call system(command)         ! list shared memory blocks on system
    call memory_arena_print_status(shmaddr)
  endif
  call MPI_Barrier(MY_World, err)

  if(rank == 0) then             ! PE0 checks everything
    do i = 1, isiz
      write(myblock,100)'BLOCK',i-1
      p = memory_block_find(shmaddr, bsz, flags, trim(myblock)) 
      call c_f_pointer(p, fp, [bsz])
      errors = 0
      do j = 1, bsz
        if(fp(j) .ne. j + ishft(i-1,24)) errors = errors + 1
      enddo
      write(0,*) trim(myblock), ' values =',bsz,' errors =',errors
    enddo
  endif

  call MPI_Barrier(MY_World, err)
  write(0, *)'I am process',rank+1,' of',isiz,' on node'
  call mpi_finalize(err)
end program
