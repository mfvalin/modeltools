*/* RPNCOMM - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 2012  Centre ESCER U.Q.A.M
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
      module local_distribution
      implicit none
      SAVE
!
!     dimension of clients is (3,1:nclients)
!         clients(:,N) = exchange number, PE number, number of values (for Nth client)
!
!     with this PE as root of transfer, there will be nclients transfers
!     clients(3,N) values will be sent to / received from  PE clients(2,N),
!     data will be read from / written into the position associated with 
!     exchange numbered clients(1,N) in buffer root_data
!
!     dimension of sources is (3,1:nsources)
!         sources(:,N) = exchange number, PE number, number of values (for Nth source)
!
!     with this PE as client for transfer, there will be nsources transfers
!     sources(3,N) values will be sent to / received from  PE sources(2,N),
!     data will be read from / written into the position associated with 
!     exchange numbered sources(1,N) in buffer client_data
!
!     a "root" PE has multiple clients to send data to or receive data from
!     an "exchange" is a "root" sending to / receiving from one or more "client(s)"
!
!     function RPN_COMM_add_dist_entry receives a multi row array where each row
!     describes such an exchange (number of values, root PE, client PE 1, ...., client PE n)
!     a value of -1 for the client PE means a non existent PE (not all exchanges have the same
!     number of "clients" and all rows mus have the same dimension)
!     this function adds an entry into the local PE's distribution table and returns the
!     table index of this entry to the caller
!
!     this table index, a data buffer, and a metadata buffer will then be passed to 
!     function RPN_COMM_do_dist that will perform the transfer and fill the metadata 
!     buffer with pertinent information
!
      type :: dist_table_entry                        ! exchange table entry 
        integer, pointer, dimension(:,:) :: clients   ! list of exchange clients for this PE
        integer, pointer, dimension(:,:) :: sources   ! list of exchange "roots" for this PE
        integer :: nclients                           ! number of exchanges where I am the root
        integer :: nsources                           ! number  of exchanges where I am a client
        integer :: comm                               ! mpi communicator
      end type dist_table_entry
      type(dist_table_entry), pointer, dimension(:) :: dist_table
      integer :: dist_table_size=-1
      integer :: dist_table_entries=-1

      end module local_distribution
!
! ========================================================================================
!
      integer function allocate_dist_table(nentries)
!     create the distribution table with nentries entries (one needed per exchange pattern)
!     function returns -1 if failure
!     function returns nentries upon success
      use local_distribution
      implicit none
      integer, intent(IN) :: nentries

      allocate_dist_table = -1       ! precondition to FAILED
      if(dist_table_size == -1 .and. nentries > 0)then
        allocate(dist_table(1:nentries))
        dist_table_size = nentries       ! number of possible entries in table
        allocate_dist_table = nentries   ! SUCCESS
        dist_table_entries = 0           ! table is EMPTY
      endif
      return
      end
!
! ========================================================================================
!
      subroutine RPN_COMM_get_dist_dims(pattern,ndata,nmeta)
!
!     get needed sizes for data and metadata buffers when performing
!     exchange with table index = pattern performed by function RPN_COMM_do_dist
!
      use local_distribution
      implicit none
      integer, intent(IN) :: pattern
      integer, intent(OUT) :: ndata, nmeta

      integer nclients, nsources, i

      ndata = -1
      nmeta = -1
      if(pattern > dist_table_entries) return  ! out of table
      if(pattern <= 0) return                  ! out of table

      ndata = 0
      nclients = dist_table(pattern)%nclients
      nsources = dist_table(pattern)%nsources
      do i = 1, nclients
         ndata = ndata + dist_table(pattern)%clients(3,i)  ! add nb of values in this exchange
      enddo
      do i = 1, nsources
         ndata = ndata + dist_table(pattern)%sources(3,i)  ! add nb of values in this exchange
      enddo
      nmeta = nclients + nsources
      end
!
! ========================================================================================
!
      integer function RPN_COMM_add_dist_entry
     %        (table,max_clients,nrows,communicator)
!
!     build a table entry describing a group of data exchanges described in array table
!
!     table(-1,j) = number of values for this exchange (must be > 0)
!     table( 0,j) = ordinal of root in communicator for this exchange (must be >= 0)
!     table( i,j) = ordinal of a PE in communicator for this exchange (-1 means entry is void)
!     table(1:max_clients,j) = ordinals of all client PEs for this exchange
!     information will be sent from root PE to client PEs or gathered on root PE from client PEs
!     max_client: maximum number of clients in exchanges
!     nrows : number of exchanges described in table
!     communicator : RPN COMM communicator for this exchange ('GRID','SUPERGRID',...)
!
!     the function will return a pattern index to be used for the actual exchange 
!     or -1 in case of failure
!
!     this function builds the exchange tables for the current PE from the
!     global table
!
      use local_distribution

      implicit none
      integer, intent(IN) :: max_clients,nrows
      integer, dimension(-1:max_clients,nrows), intent(IN) :: table
      character *(*), intent(IN) :: communicator

      integer RPN_COMM_comm
      external RPN_COMM_comm
      integer the_comm, this_pe
      integer my_clients, my_sources, i, j, ierr
      integer no_client, no_source

      RPN_COMM_add_dist_entry = -1   ! precondition to FAILED
      if(dist_table_entries >= dist_table_size) return    ! OUCH, table is full

      the_comm = RPN_COMM_comm(communicator)
      call RPN_COMM_rank( communicator, this_pe ,ierr )
      my_clients = 0
      my_sources = 0
      dist_table_entries = dist_table_entries + 1           ! new table entry
      do j =  1 , nrows                                     ! find number of clients and sources
         if(table(0,j) == this_pe) then  ! I am the source for this row
            do i = 2, max_clients                           ! count valid clients in row
               if(table(i,j) >= 0) then                   ! I have a client
                  my_clients = my_clients + 1
               endif
            enddo
         else                             ! I am not the root, am i a client ?
            do i = 2, max_clients
               if(table(i,j) == this_pe) then ! I am a client
                  my_sources = my_sources + 1
                  exit       ! I can only be client once in an exchange
               endif
            enddo
         endif
      enddo

      allocate    ! allocate client table
     %  (dist_table(dist_table_entries)%clients(3,1:my_clients))
      allocate    ! allocate sources table
     %  (dist_table(dist_table_entries)%sources(3,1:my_sources))
      no_client = 0
      dist_table(dist_table_entries)%nclients = my_clients  ! number of clients
      dist_table(dist_table_entries)%nsources = my_sources  ! number of root sources
      dist_table(dist_table_entries)%comm = the_comm        ! grid communicator
      no_source = 0
!
!     rows in clients contain the exchange number, client PE ordinal, and data length
!     rows in sources contain the exchange number, root PE ordinal, and data length
!
      do j =  1 , nrows
        if(table(0,j) == this_pe) then  ! I am the root for this exchange
          no_client = no_client + 1
          do i = 2, max_clients                                  ! count valid clients in row
            if(table(i,j) >= 0) then
               no_client = no_client + 1
               dist_table(dist_table_entries)%clients(1,no_client) = j   ! exchange number
               dist_table(dist_table_entries)%clients(2,no_client)
     %             = table(i,j)                                          ! PE ordinal
               dist_table(dist_table_entries)%clients(3,no_client)
     %             = table(-1,j)                                         ! data length
            endif
          enddo
        else                                        ! I am not the root, am i a client ?
           do i = 2, max_clients
             if(table(i,j) == this_pe) then ! I am a client
               no_source = no_source + 1
               dist_table(dist_table_entries)%sources(1,no_source) = j   ! row number
               dist_table(dist_table_entries)%sources(2,no_source)
     %             = table(1,j)                                          ! PE ordinal
               dist_table(dist_table_entries)%sources(3,no_source)
     %             = table(-1,j)                                         ! data length
               exit        ! I can only be client once in an exchange
             endif
           enddo
         endif
      enddo

      RPN_COMM_add_dist_entry = dist_table_entries   !  SUCCESS
      return
      end
!
! ========================================================================================
!
      integer function RPN_COMM_do_dist(pattern,values,nvalues,
     %        metadata,nmeta,nrows,from_root)
!
!     pattern is the pattern number returned by RPN_COMM_add_dist_entry for this data exchange
!
!     values : 1D array containing data to distribute or data to receive
!     nvalues : dimension of array values, max number of values to be distributed/collected
!     metadata : describes what was sent/received
!         metadata(1,N) is offset of data in array values
!         metadata(2,N) is number of data elements to be collected
!         metadata(3,N) is the exchange number
!     nrows : number of exchanges
!     from_root : distribution or collection flag
!          .true. : data flows from root PE to client PEs
!          .false. : data flows from client PEs to root PE
!
      use local_distribution
      implicit none
      integer, intent(IN) :: pattern, nvalues, nrows, nmeta
      integer , dimension(3,nrows), intent(OUT) :: metadata
      integer , dimension(nvalues), intent(INOUT) :: values
      logical, intent(IN) :: from_root

      include 'mpif.h'

      integer row, j, nexchanges, request, offset, length
      integer, pointer, dimension(:,:) :: clients, sources
      integer, allocatable, dimension(:,:) :: statuses
      integer, allocatable, dimension(:) :: requests
      integer nclients, nsources, the_comm, ierr

      RPN_COMM_do_dist = -1   ! precondition to FAILED
      if(pattern > dist_table_entries) return  ! out of table ?
      if(pattern <= 0) return

      clients => dist_table(pattern)%clients  ! array describing my clients
      sources => dist_table(pattern)%sources  ! array describing my root data sources
      nclients = dist_table(pattern)%nclients ! number of clients
      nsources = dist_table(pattern)%nsources ! number of root data sources
      the_comm = dist_table(pattern)%comm     ! MPI communicator
      allocate( statuses(MPI_STATUS_SIZE,nclients+nsources) )
      allocate( requests(nclients+nsources) )

      request = 0
      offset = 1
      do j = 1 , nclients  ! as root send to / receive from clients
         request = request + 1
         row = clients(1,j)
         length = clients(3,j)
         metadata(1,request) = offset
         metadata(2,request) = length
         metadata(3,request) = row
         if(from_root) then ! send
            call MPI_isend(values(offset),length,MPI_INTEGER,
     %           clients(2,j),row,the_comm,requests(request),ierr)
         else
            call MPI_irecv(values(offset),length,MPI_INTEGER,
     %           clients(2,j),row,the_comm,requests(request),ierr)
         endif
         offset = offset + length
      enddo
      do j = 1 , nsources  ! as client receive from / send to roots
         request = request + 1
         row = sources(1,j)
         length = sources(3,j)
         metadata(1,request) = offset
         metadata(2,request) = length
         metadata(3,request) = row
         if(from_root) then ! receive
            call MPI_irecv(values(offset),length,MPI_INTEGER,
     %           sources(2,j),row,the_comm,requests(request),ierr)
         else
            call MPI_isend(values(offset),length,MPI_INTEGER,
     %           sources(2,j),row,the_comm,requests(request),ierr)
         endif
         offset = offset + length
      enddo
!
!     wait for all transfers to complete
!
      call MPI_waitall(request,requests,statuses,ierr)
!
      deallocate(statuses)
      deallocate(requests)
      RPN_COMM_do_dist = 0    ! SUCCESS
      return
      end
!
! ========================================================================================
!
      integer function RPN_COMM_test_dist()
      implicit none
      RPN_COMM_test_dist = 0
      return
      end