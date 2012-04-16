      module local_distribution
      implicit none
      SAVE
      type :: dist_table_entry
        integer, pointer, dimension(:,:) :: send_table, recv_table
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

      integer function add_dist_table_entry(table,max_clients,nrows)

      use local_distribution
      use rpn_comm

      implicit none
      integer max_clients,nrows
      integer table(2,max_clients+1,nrows)

      integer my_clients, i, j
      integer no_client

      add_dist_table_entry = -1   ! precondition to FAILED
      if(dist_table_entries >= dist_table_size) return    ! OUCH, table is full

      my_clients = 0
      dist_table_entries = dist_table_entries + 1
      do j =  1 , nrows
         if(pe_mex == table(1,1,j) .and. pe_mey == table(2,1,j)) then  ! I am the source for this row
            do i = 2, max_clients                                      ! count valid clients in row
              if(table(1,i,j) .ne. 0) my_clients = my_clients + 1
            enddo
         endif
      enddo

      allocate    ! allocate client table
     %  (dist_table(dist_table_entries)%send_table(2,0:my_clients))
      allocate    ! allocate supplier table
     %  (dist_table(dist_table_entries)%recv_table(2,0:my_clients))

      no_client = 0
      dist_table(dist_table_entries)%send_table(1,0) = my_clients  ! number of clients
      dist_table(dist_table_entries)%send_table(2,0) = pe_grid     ! grid communicator
      do j =  1 , nrows
         if(pe_mex == table(1,1,j) .and. pe_mey == table(2,1,j)) then  ! I am the source for this row
            do i = 2, max_clients                                      ! count valid clients in row
              if(table(1,i,j) .ne. 0) then
                 no_client = no_client + 1
                 dist_table(dist_table_entries)%send_table(1,no_client)
     %            = pe_id(table(1,i,j),table(2,i,j))
                 dist_table(dist_table_entries)%send_table(2,no_client)
     %            = j
              endif
            enddo
         endif
      enddo

      add_dist_table_entry = 0   !  SUCCESS
      return
      end
