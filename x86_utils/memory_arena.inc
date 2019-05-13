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
function memory_arena_init(mem, nsym, size) result(id) BIND(C,name='function memory_arena_init')
  import :: C_PTR, C_INT
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(IN), value :: nsym, size
  integer(C_INT) :: id
end function memory_arena_init
function memory_block_find(mem, size, flags, name) result(ptr) BIND(C,name='memory_block_find')
  import :: C_PTR, C_INT, C_CHAR
  type(C_PTR), intent(IN), value :: mem
  integer(C_INT), intent(OUT) :: size, flags
  character(C_CHAR), dimension(*), intent(IN) :: name
  type(C_PTR) :: ptr
end function memory_block_find
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
function memory_arena_create_shared(id, nsym, size) result(ptr) BIND(C,name='memory_arena_create_shared')
  import :: C_PTR, C_INT
  integer, intent(OUT) :: id
  integer, intent(IN), value :: nsym, size
  type(C_PTR) :: ptr
end function memory_arena_create_shared
function memory_arena_from_id(id) result(ptr) BIND(C,name='memory_arena_from_id')
  import :: C_PTR, C_INT
  integer, intent(IN), value :: id
  type(C_PTR) :: ptr
end function memory_arena_from_id
end interface