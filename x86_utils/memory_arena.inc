interface
function memory_arena_set_id(id) result(me) BIND(C,name='memory_arena_set_id')
  import :: C_INT
  integer(C_INT), intent(IN), value :: id
  integer(C_INT) :: me
end function memory_arena_set_id
subroutine memory_arena_print_status(mem) BIND(C,name='memory_arena_print_status')
  import :: C_PTR
  type(C_PTR), intent(IN), value :: mem
end subroutine memory_arena_print_status
function memory_arena_init(mem, nsym, size) result(me) BIND(C,name='memory_arena_init')
  import :: C_PTR, C_INT
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(IN), value :: nsym, size
  integer(C_INT) :: me
end function memory_arena_init
function update_local_table(mem) result(status) BIND(C,name='update_local_table')
  import :: C_PTR, C_INT
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT) :: status
end function update_local_table
function master_arena_init(mem, nsym, size) result(id) BIND(C,name='master_arena_init')
  import :: C_PTR, C_INT
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(IN), value :: nsym, size
  integer(C_INT) :: id
end function master_arena_init
function memory_block_find(mem, size, flags, name) result(ptr) BIND(C,name='memory_block_find')
  import :: C_PTR, C_INT, C_CHAR
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(OUT) :: size, flags
  character(C_CHAR), dimension(*), intent(IN) :: name
  type(C_PTR) :: ptr
end function memory_block_find
function memory_block_find_wait(mem, size, flags, name, timeout) result(ptr) BIND(C,name='memory_block_find_wait')
  import :: C_PTR, C_INT, C_CHAR
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(OUT) :: size, flags
  character(C_CHAR), dimension(*), intent(IN) :: name
  integer(C_INT), intent(IN), value :: timeout
  type(C_PTR) :: ptr
end function memory_block_find_wait
function memory_block_mark_init(mem, name) result(ptr) BIND(C,name='memory_block_mark_init')
  import :: C_PTR, C_CHAR
  type(C_PTR), intent(IN), value :: mem
  character(C_CHAR), dimension(*), intent(IN) :: name
  type(C_PTR) :: ptr
end function memory_block_mark_init
function memory_block_create(mem, size, name) result(ptr) BIND(C,name='memory_block_create')
  import :: C_PTR, C_INT, C_CHAR
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(IN), value :: size
  character(C_CHAR), dimension(*), intent(IN) :: name
  type(C_PTR) :: ptr
end function memory_block_create
function memory_allocate_shared(shmid, size) result(ptr) BIND(C,name='memory_allocate_shared')
  import :: C_PTR, C_INT
  integer, intent(OUT) :: shmid
  integer, intent(IN), value :: size
  type(C_PTR) :: ptr
end function memory_allocate_shared
function memory_arena_create_shared(shmid, nsym, size) result(ptr) BIND(C,name='memory_arena_create_shared')
  import :: C_PTR, C_INT
  integer, intent(OUT) :: shmid
  integer, intent(IN), value :: nsym, size
  type(C_PTR) :: ptr
end function memory_arena_create_shared
function master_arena_create_shared(shmid, nsym, size) result(ptr) BIND(C,name='master_arena_create_shared')
  import :: C_PTR, C_INT
  integer, intent(OUT) :: shmid
  integer, intent(IN), value :: nsym, size
  type(C_PTR) :: ptr
end function master_arena_create_shared
function memory_address_from_id(shmid) result(ptr) BIND(C,name='memory_address_from_id')
  import :: C_PTR, C_INT
  integer, intent(IN), value :: shmid
  type(C_PTR) :: ptr
end function memory_address_from_id
function memory_arena_from_master(mem) result(ptr) BIND(C,name='memory_arena_from_master')
  import :: C_PTR
  type(C_PTR), intent(IN), value :: mem
  type(C_PTR) :: ptr
end function memory_arena_from_master
function memory_arena_from_master_id(shmid) result(ptr) BIND(C,name='memory_arena_from_master_id')
  import :: C_PTR, C_INT
  integer, intent(IN), value :: shmid
  type(C_PTR) :: ptr
end function memory_arena_from_master_id
end interface
