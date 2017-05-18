interface
!
  function libc_fopen(path,mode) result(stream) bind(C,name='fopen')
    import :: C_CHAR, C_PTR
    character(C_CHAR), dimension(*), intent(IN) :: path, mode
    type(C_PTR) :: stream
  end function libc_fopen
!
  function libc_fflush(stream) result(status) bind(C,name='fflush')
    import :: C_INT, C_PTR
    type(C_PTR), intent(IN), value :: stream
    integer(C_INT) :: status
  end function libc_fflush
!
  function libc_fclose(stream) result(status) bind(C,name='fclose')
    import :: C_INT, C_PTR
    type(C_PTR), intent(IN), value :: stream
    integer(C_INT) :: status
  end function libc_fclose
!
  function libc_fileno(stream) result(fd) bind(C,name='fileno')
    import :: C_INT, C_PTR
    type(C_PTR), intent(IN), value :: stream
    integer(C_INT) :: fd
  end function libc_fileno
!
  function libc_write(fd,buf,count) result(n) bind(C,name='write')
    import :: C_INT, C_SIZE_T, C_PTR
    type(C_PTR), intent(IN), value :: buf
    integer(C_INT), intent(IN), value :: fd
    integer(C_SIZE_T), intent(IN), value :: count
    integer(C_SIZE_T) :: n
  end function libc_write
!
  function libc_read(fd,buf,count) result(n) bind(C,name='read')
    import :: C_INT, C_SIZE_T, C_PTR
    type(C_PTR), intent(IN), value :: buf
    integer(C_INT), intent(IN), value :: fd
    integer(C_SIZE_T), intent(IN), value :: count
    integer(C_SIZE_T) :: n
  end function libc_read
!
  function libc_close(fd) result(status) bind(C,name='close')
    import :: C_INT, C_SIZE_T
    integer(C_INT), intent(IN), value :: fd
    integer(C_INT) :: status
  end function libc_close
!
  function libc_strlen(str) result(slen) bind(C,name='strlen')
    import :: C_CHAR, C_SIZE_T
    character(C_CHAR), dimension(*), intent(IN) :: str
    integer(C_SIZE_T) :: slen
  end function libc_strlen
!
  function libc_strnlen(str,mlen) result(slen) bind(C,name='strnlen')
    import :: C_CHAR, C_SIZE_T
    character(C_CHAR), dimension(*), intent(IN) :: str
    integer(C_SIZE_T), intent(IN), value :: mlen
    integer(C_SIZE_T) :: slen
  end function libc_strnlen
!
  function libc_alarm(seconds) result(seconds_since) BIND(C,name='alarm')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), value :: seconds
    integer(C_INT) :: seconds_since
  end function libc_alarm
!
  function libc_gethostid() result(id) BIND(C,name='gethostid')  ! interface to libc gethostid
    import :: C_INT
    integer(C_INT) :: id
  end function libc_gethostid
!
  function libc_malloc(siz) result (p) bind(C,name='malloc')
    import C_SIZE_T, C_PTR
    integer(C_SIZE_T), intent(IN), value :: siz
    type(C_PTR) :: p
  end function libc_malloc
!
  subroutine libc_free(p) bind(C,name='free')
    import C_PTR
    type(C_PTR), intent(IN), value :: p
  end subroutine libc_free
!
  function libc_usleep(duration) result(status) bind(C,name='usleep')
    use ISO_C_BINDING
    integer(C_INT), intent(IN), value :: duration
    integer(C_INT) :: status
  end function libc_usleep
!
  subroutine libc_sleep(duration) bind(C,name='sleep')
    use ISO_C_BINDING
    integer(C_INT), intent(IN), value :: duration
  end subroutine libc_sleep
!
end interface

