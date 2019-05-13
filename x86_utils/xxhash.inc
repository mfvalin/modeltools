interface                                                                 !InTf
  function Hash32(seed, input, length) result(hash) BIND(C,name='Hash32') !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: input                               !InTf
    integer(C_INT), intent(IN), value :: seed, length                     !InTf
    integer(C_INT) :: hash                                                !InTf
  end function Hash32                                                     !InTf
  function Hash32_short(input, length) result(hash) BIND(C,name='Hash32_short') !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: input                               !InTf
    integer(C_INT), intent(IN), value :: length                           !InTf
    integer(C_INT) :: hash                                                !InTf
  end function Hash32_short                                               !InTf
  function Hash32_create() result(hash_struct) BIND(C,name='Hash32_create') !InTf
    import :: C_PTR                                                       !InTf
    type(C_PTR) :: hash_struct                                            !InTf
  end function Hash32_create                                              !InTf
  subroutine Hash32_init(hash_struct, seed) BIND(C,name='Hash32_init')    !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: hash_struct                         !InTf
    integer(C_INT), intent(IN), value :: seed                             !InTf
  end subroutine Hash32_init                                              !InTf
  function Hash32_add(hash_struct, input, length) result(status) BIND(C,name='Hash32_add') !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: hash_struct, input                  !InTf
    integer(C_INT), intent(IN), value :: length                           !InTf
    integer(C_INT) :: status                                              !InTf
  end function Hash32_add                                                 !InTf
  function Hash32_finalize(hash_struct) result(hash) BIND(C,name='Hash32_finalize')   !InTf
    import :: C_PTR, C_INT                                                !InTf
    type(C_PTR), intent(IN), value :: hash_struct                         !InTf
    integer(C_INT) :: hash                                                !InTf
  end function Hash32_finalize                                            !InTf
end interface                                                             !InTf
