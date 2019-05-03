 interface
   function circular_buffer_init(p, nbytes) result(limit) bind(C,name='circular_buffer_init')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(IN), value :: nbytes
   end function circular_buffer_init
   function circular_buffer_create_shared(shmid, nbytes) result(p) BIND(C,name='circular_buffer_create_shared')
     import :: C_PTR, C_INT
     integer(C_INT), intent(OUT) :: shmid
     integer(C_INT), intent(IN), value :: nbytes
     type(C_PTR) :: p
   end function circular_buffer_create_shared
   function circular_buffer_detach_shared(p) result(status) BIND(C,name='circular_buffer_detach_shared')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT) :: status
   end function circular_buffer_detach_shared
   function circular_buffer_create(nbytes) result(p) BIND(C,name='')
     import :: C_PTR, C_INT
     integer(C_INT), intent(IN), value :: nbytes
     type(C_PTR) :: p
   end function circular_buffer_create
   function circular_buffer_space_available(p) result(n) BIND(C,name='circular_buffer_space_available')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT) :: n
   end function circular_buffer_space_available
   function circular_buffer_wait_space_available(p, na) result(n) BIND(C,name='circular_buffer_wait_space_available')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(IN), value :: na
     integer(C_INT) :: n
   end function circular_buffer_wait_space_available
   function circular_buffer_data_available(p) result(n) BIND(C,name='circular_buffer_data_available')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT) :: n
   end function circular_buffer_data_available
   function circular_buffer_wait_data_available(p, na) result(n) BIND(C,name='circular_buffer_wait_data_available')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(IN), value :: na
     integer(C_INT) :: n
   end function circular_buffer_wait_data_available
   function circular_buffer_start(p) result(start) BIND(C,name='circular_buffer_start')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     type(C_PTR) :: start
   end function circular_buffer_start
   function circular_buffer_data_in(p) result(inp) BIND(C,name='circular_buffer_data_in')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     type(C_PTR) :: inp
   end function circular_buffer_data_in
   function circular_buffer_data_out(p) result(outp) BIND(C,name='circular_buffer_data_out')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     type(C_PTR) :: outp
   end function circular_buffer_data_out
   function circular_buffer_advance_in(p, n1, n2) result(inp) BIND(C,name='circular_buffer_advance_in')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(OUT) :: n1, n2
     type(C_PTR) :: inp
   end function circular_buffer_advance_in
   function circular_buffer_advance_out(p, n1, n2) result(outp) BIND(C,name='circular_buffer_advance_out')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(OUT) :: n1, n2
     type(C_PTR) :: outp
   end function circular_buffer_advance_out
   function circular_buffer_atomic_get(p, dst, ndst) result(n) BIND(C,name='circular_buffer_atomic_get')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(IN), value :: ndst
     integer(C_INT), dimension(*), intent(OUT) :: dst
     integer(C_INT) :: n
   end function circular_buffer_atomic_get
   function circular_buffer_atomic_put(p, src, nsrc) result(n) BIND(C,name='circular_buffer_atomic_put')
     import :: C_PTR, C_INT
     type(C_PTR), intent(IN), value :: p
     integer(C_INT), intent(IN), value :: nsrc
     integer(C_INT), dimension(*), intent(IN) :: src
     integer(C_INT) :: n
   end function circular_buffer_atomic_put
 end interface
