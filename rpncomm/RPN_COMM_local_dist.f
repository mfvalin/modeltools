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
      type :: dist_table_entry
        integer, pointer, dimension(:,:) :: clients, sources
      end type dist_table_entry
      type(dist_table_entry), pointer, dimension(:) :: dist_table
      integer :: dist_table_size=-1
      integer :: dist_table_entries=-1

      end module local_distribution

      integer function allocate_dist_table(nentries)

      use local_distribution
      implicit none
      integer nentries

      allocate_dist_table = -1       ! precondition to FAILED
      if(dist_table_size == -1 .and. nentries > 0)then
        allocate(dist_table(1:nentries))
        allocate_dist_table = 0
        dist_table_size = nentries   ! SUCCESS
        dist_table_entries = 0
      endif
      return
      end

      integer function RPN_COMM_add_dist_entry(table,max_clients,nrows)
!
!     table(1,i,j) = X coordinate of a PE (0 means entry is void)
!     table(2,i,j) = Y coordinate of a PE (0 means entry is void)
!     table(:,1,j) = coordinates of root PE for the exchange described by the row
!     table(:,2:max_clients,j) = coordinates of client PEs for the exchange described by the row
!     information will be sent from root PE to client PEs or gathered on root PE from client PEs
!     max_client: max number of clients in rows
!     nrows : number of exchanges described in table
!
!     the rows must be consistent with the data matrix used in RPN_COMM_do_dist
!
!     the function will return an index to be used for the actual exchange or -1
!     in case of failure
!
!     add_dist_table_entry builds the exchange table for the current PE from the
!     global table
!
      use local_distribution
      use rpn_comm

      implicit none
      integer max_clients,nrows
      integer table(2,max_clients+1,nrows)

      integer my_clients, my_sources, i, j
      integer no_client, no_source

      RPN_COMM_add_dist_entry = -1   ! precondition to FAILED
      if(dist_table_entries >= dist_table_size) return    ! OUCH, table is full

      my_clients = 0
      my_sources = 0
      dist_table_entries = dist_table_entries + 1           ! new table entry
      do j =  1 , nrows                                     ! find number of clients and sources
         if(pe_id(table(1,1,j),table(2,1,j)) == pe_me) then  ! I am the source for this row
            do i = 2, max_clients                           ! count valid clients in row
               if(table(1,i,j) <> 0) then                   ! I have a client
                  my_clients = my_clients + 1
               endif
            enddo
         else                                                    ! not the source, am i a client ?
            do i = 2, max_clients
               if(pe_id(table(1,i,j),table(2,i,j)) == pe_me) then ! I am a client
                  my_sources = my_sources + 1
               endif
            enddo
         endif
      enddo

      allocate    ! allocate client table
     %  (dist_table(dist_table_entries)%clients(2,0:my_clients))
      allocate    ! allocate sources table
     %  (dist_table(dist_table_entries)%sources(2,0:my_sources))
      no_client = 0
      no_source = 0
!
!     entry 0 in clients is different, it contains the nb of clients and the communicator
!     the other entries contain the data row number and client PE number
!
!     entry 0 in sources is different, it contains the nb of sources and the communicator
!     the other entries contain the data row number and source PE number
!
      dist_table(dist_table_entries)%clients(1,0) = my_clients  ! number of clients
      dist_table(dist_table_entries)%clients(2,0) = pe_grid     ! grid communicator
      dist_table(dist_table_entries)%sources(1,0) = my_sources  ! number of clients
      dist_table(dist_table_entries)%sources(2,0) = pe_grid     ! grid communicator
      do j =  1 , nrows
         if(pe_id(table(1,1,j),table(2,1,j)) == pe_me) then  ! I am the source for this row
            do i = 2, max_clients                                      ! count valid clients in row
               if(table(1,i,j) <> 0) then
                  no_client = no_client + 1
                  dist_table(dist_table_entries)%clients(1,no_client)
     %             = j     ! data row number
                  dist_table(dist_table_entries)%clients(2,no_client)
     %             = pe_id(table(1,i,j),table(2,i,j))   ! PE number
               endif
            enddo
         else                                                    ! not the source, am i a client ?
            do i = 2, max_clients
               if(pe_id(table(1,i,j),table(2,i,j)) == pe_me) then ! I am a client
                  no_source = no_source + 1
                  dist_table(dist_table_entries)%sources(1,no_source)
     %             = j     ! data row number
                  dist_table(dist_table_entries)%sources(2,no_source)
     %             = pe_id(table(1,1,j),table(2,1,j))   ! PE number
                  exit    ! can only be client once in a row
               endif
            enddo
         endif
      enddo

      RPN_COMM_add_dist_entry = dist_table_entries   !  SUCCESS
      return
      end
      integer function RPN_COMM_do_dist(pattern,values,nvalues,
     %        metadata,nrows,dist_collect)
!
!     pattern is the pattern number returned by RPN_COMM_add_dist_entry for this data exchange
!
!     values : 1D array containing data to distribute or data to receive
!     nvalues : dimension of data array, max number of values to be distributed/collected
!     metadata : describes what is to be sent/received
!         metadata(1,row) is offset of data in array values
!         metadata(2,row) is number of data elements to be collected
!     nrows : number of exchanges
!     dist_collect : distribution or collection flag
!          .true. : data flows from root PE to client PEs
!          .false. : data flows from client PEs to root PE
!
      use local_distribution
      use rpn_comm
      implicit none
      integer :: pattern, nvalues, nrows
      integer , dimension(2,nrows) :: metadata
      integer , dimension(nvalues) :: values
      logical :: dist_collect

      include 'mpif.h'

      integer row, j, nexchanges, request, offset, length
      integer, pointer, dimension(:,:) :: clients, sources
      integer, allocatable, dimension(:,:) :: statuses
      integer, allocatable, dimension(:) :: requests

      RPN_COMM_do_dist = -1   ! precondition to FAILED
      if(pattern > dist_table_entries) return  ! out of table ?
      if(pattern <= 0) return

      clients => dist_table(pattern)%clients
      sources => dist_table(pattern)%sources
      allocate( statuses(MPI_STATUS_SIZE,clients(1,0)+sources(1,0)) )
      allocate( requests(clients(1,0)+sources(1,0)) )

      request = 0
      do j = 1 , clients(1,0)  ! as root send to / receive from clients
         request = request + 1
         row = clients(1,j)
         offset = metadata(1,row)
         length = metadata(2,row)
         if(dist_collect) then ! send
            call MPI_isend(values(offset),length,MPI_INTEGER,
     %           clients(2,j),0,clients(2,0),requests(request))             
         else
            call MPI_irecv(values(offset),length,MPI_INTEGER,
     %           clients(2,j),0,clients(2,0),requests(request))             
         endif
      enddo
      do j = 1 , sources(1,0)  ! as root receive from / send to clients
         request = request + 1
         row = sources(1,j)
         offset = metadata(1,row)
         length = metadata(2,row)
         if(dist_collect) then ! receive
            call MPI_irecv(values(offset),length,MPI_INTEGER,
     %           sources(2,j),0,sources(2,0),requests(request))             
         else
            call MPI_isend(values(offset),length,MPI_INTEGER,
     %           sources(2,j),0,sources(2,0),requests(request))             
         endif
      enddo

      deallocate(statuses)
      RPN_COMM_do_dist = 0    ! SUCCESS
      return
      end
