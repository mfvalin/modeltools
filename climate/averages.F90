  module averages_common
    use iso_c_binding
    implicit none

    interface
      subroutine f_exit(code) BIND(C,name='exit')
        import :: C_INT
        integer(C_INT), intent(IN), value :: code
      end subroutine f_exit
    end interface

    type :: field
      real *8, dimension(:,:,:), pointer :: stats
      integer :: ni, nj
      integer :: nsamples
      integer :: npas_min, npas_max
      integer :: ip1
      character(len=12) :: etiket
      character(len=4)  :: nomvar
      character(len=2)  :: typvar
      character(len=1)  :: grtyp
    end type

    integer, parameter :: PAGE_SIZE  = 256   ! 256 fields per page in table
    integer, parameter :: ENTRY_MASK = 255   ! mask to get slot number from index
    type :: page
      type(field), dimension(:), pointer :: p
    end type

    integer, parameter :: MAX_PAGES = 1024   ! 1024 pages at most
    integer, parameter :: PAGE_SHIFT = -10   ! right shift count to get page number from index
    type(page), dimension(:), pointer, save :: ptab => NULL()

    integer, save :: next = -1
    integer, save :: verbose = 2           ! ERROR + WARNING

  contains

    subroutine create_page(n)  ! create page n that will contain PAGE_SIZE entries
      implicit none
      integer, intent(IN) :: n
      integer :: i, status

      if(n < 0 .or. n >= MAX_PAGES) then
        print *,"FATAL: page number",n,' is invalid. it must be between 0 and',MAX_PAGES-1
        call f_exit(1)
      endif
      if(associated(ptab(n)%p)) then
        print *,"FATAL: page number",n,' is already allocated'
        call f_exit(1)
      endif
      allocate(ptab(n)%p(0:PAGE_SIZE-1),STAT=status)  ! allocate page 1
      if(status .ne. 0) then
        print *,"FATAL: page number",n,' cannot be allocated'
        call f_exit(1)
      endif
      do i = 0, PAGE_SIZE-1
        nullify(ptab(n)%p(i)%stats)         ! nullify entries in page 1
      enddo
    end subroutine

    function create_page_entry(pg,slot,ni,nj) result(ix) ! create entry slot in page pg
      implicit none
      integer, intent(IN) :: pg, slot, ni, nj
      integer ix, status

      ix = -1
      if(slot < 0 .or. slot >= PAGE_SIZE ) then
        if(verbose > 0) print *,"ERROR: slot number",slot," not valid"
        return
      endif
      if(pg < 0 .or. pg >= MAX_PAGES ) then
        if(verbose > 0) print *,"ERROR: page number",pg," not valid"
        return
      endif

      if(.not. associated(ptab(pg)%p)) then
        call create_page(pg)
      endif

      if(associated(ptab(pg)%p(slot)%stats)) return

      ix = slot + PAGE_SIZE * pg
      allocate(ptab(pg)%p(slot)%stats(ni,nj,2),STAT=status)
      if(verbose > 4) print *,'DEBUG: allocated stats - slot, ni,nj =',ix,ni,nj
      if(status .ne. 0) then
        print *,"FATAL: entry number",slot,'in page',pg,' cannot be allocated'
        call f_exit(1)
      endif
      ptab(pg)%p(slot)%stats = 0.0   ! set sum and sum of squares to zero
      ptab(pg)%p(slot)%ni = ni       ! dimensions of field
      ptab(pg)%p(slot)%nj = nj
      ptab(pg)%p(slot)%nsamples = 0  ! no samples
      ptab(pg)%p(slot)%npas_min = 999999999
      ptab(pg)%p(slot)%npas_max = -1
      ptab(pg)%p(slot)%ip1 = -1
      ptab(pg)%p(slot)%etiket = ""
      ptab(pg)%p(slot)%nomvar = ""
      ptab(pg)%p(slot)%typvar = ""
      ptab(pg)%p(slot)%grtyp = ""

      return
    end function

    function new_page_entry(ni,nj) result(ix) ! create a new entry in tables
      implicit none
      integer, intent(IN) :: ni, nj
      integer :: ix
      integer :: pg, slot

      next = next + 1
      pg = next / PAGE_SIZE
      slot = mod(next,PAGE_SIZE)
      ix = create_page_entry(pg,slot,ni,nj)
      if(ix < 0) then
        print *,"FATAL: error creating entry",slot,' in page',pg
        call f_exit(1)
      endif
      return
    end function new_page_entry

    function process_entry(z,ni,nj,ip1,npas,etiket,nomvar,typvar,grtyp) result(ix)
      implicit none
      integer, intent(IN) :: ni, nj, ip1, npas
      character(len=*), intent(IN) :: etiket, nomvar, typvar, grtyp
      real, dimension(ni,nj), intent(IN) :: z
      integer :: ix
      integer :: i, pg, slot
      type(field), pointer :: p

      ix = -1
      do i = 0 , next
        slot = iand(i,ENTRY_MASK)
        pg = ishft(i,PAGE_SHIFT)
        p => ptab(pg)%p(slot)
        if(trim(p%etiket) .ne. trim(etiket)) cycle
        if(trim(p%nomvar) .ne. trim(nomvar)) cycle
        if(trim(p%typvar) .ne. trim(typvar)) cycle
        if(trim(p%grtyp) .ne. trim(grtyp)) cycle
        if(p%ip1 .ne. ip1) cycle
        if(p%ni .ne. ni) cycle
        if(p%nj .ne. nj) cycle
        ix = i                      ! entry has been found
        if(verbose > 4) print *,'DEBUG: found at',ix,p%nsamples
        exit
      enddo
      if(ix == -1) then
        ix = new_page_entry(ni,nj)
        slot = iand(ix,ENTRY_MASK)
        pg = ishft(ix,PAGE_SHIFT)
        p => ptab(pg)%p(slot)
        p%etiket = etiket
        p%nomvar = nomvar
        p%typvar = typvar
        p%grtyp = grtyp
        p%ip1 = ip1
        p%ni = ni
        p%nj = nj
        if(verbose > 4) print *,'DEBUG: created at',ix
      endif
      p%npas_max = max(p%npas_max,npas)
      p%npas_min = min(p%npas_min,npas)
      p%nsamples = p%nsamples + 1
      p%stats(1:ni,1:nj,1) = p%stats(1:ni,1:nj,1) + z(1:ni,1:nj)
      p%stats(1:ni,1:nj,2) = p%stats(1:ni,1:nj,2) + z(1:ni,1:nj)*z(1:ni,1:nj)
    end function process_entry
  end module averages_common

  program averages
    use averages_common
    implicit none
    integer :: curarg
    character (len=8) :: option
    character (len=2048) :: filename, progname
    integer :: arg_count, arg_len, status, i
    integer :: first_file
    logical :: file_exists, strict, test
    real, dimension(:,:), pointer :: z
    integer :: ix

    curarg = 1
    first_file = 0
    strict = .false.                   ! do not abort on error
    test = .false.
    arg_count = command_argument_count()
    call get_command_argument(0,progname,arg_len,status)
    do while(curarg <= arg_count)
      call get_command_argument(curarg,option,arg_len,status)
      if(option(1:1) .ne. '-') exit      ! does not start with -, must be a file name
      curarg = curarg + 1
      if(status .ne. 0) then
        print *,"FATAL: option is too long :'"//trim(option)//"..."//"'"
        call f_exit(1)
      endif
      if( option == '-h' .or. option == '--help' ) then
        print *,'USAGE: '//trim(progname)//' [-h|--help] --strict] [-test] [-q[q]] [-v[v][v]] [--|-f] file_1 ... file_n'
        call f_exit(0)
      else if( option == '-strict' ) then
        strict = .true.                 ! abort on ERROR 
        verbose = 1
      else if( option == '-test' ) then
        test = .true.                   ! test mode
      else if( option == '-q' ) then
        if(.not. strict) verbose = 1    ! ERROR 
      else if( option == '-qq' ) then
        if(.not. strict) verbose = 0    ! NO messages
      else if( option == '-v' ) then
        verbose = 3                     ! + INFO
      else if( option == '-vv' ) then
        verbose = 4                     ! + NOTE
      else if( option == '-vvv' ) then
        verbose = 5                     ! + DEBUG
      else if( option == '--' .or. option == '-f' ) then
        first_file = curarg
        if(verbose > 4) print *,'DEBUG: first file is argument',curarg
        exit
      else
        if(verbose > 1) print *,"WARNING: unrecognized option = '"//trim(option)//"'"
        cycle
      endif
      if(verbose > 3) print *,"NOTE: option = '"//trim(option)//"'"
    enddo

    allocate(ptab(0:MAX_PAGES-1))   ! allocate and initialize page table
    do i = 0 , MAX_PAGES-1
      nullify(ptab(i)%p)
    enddo

    call create_page(0)
    if(test) then
      allocate(z(135,231))
      z = 1.1
      ix = process_entry(z,135,231,12345,1,"etiket","NOM","TP","Z")
      print *,'ix=',ix
      deallocate(z)
      allocate(z(231,135))
      z = 1.1
      ix = process_entry(z,231,135,12345,1,"etiket","NOM","TP","Z")
      print *,'ix=',ix
      deallocate(z)
      allocate(z(135,231))
      z = 1.1
      ix = process_entry(z,135,231,12345,1,"etiket","NOM","TP","Z")
      print *,'ix=',ix
      deallocate(z)
      allocate(z(231,135))
      z = 1.1
      ix = process_entry(z,231,135,12345,1,"etiket","NOM","TP","Z")
      print *,'ix=',ix
      deallocate(z)
      allocate(z(231,135))
      z = 1.1
      ix = process_entry(z,231,135,12345,1,"etiket","NOM2","TP","Z")
      print *,'ix=',ix
      ix = process_entry(z,231,135,12345,1,"etiket","NOM2","TP","Z")
      print *,'ix=',ix
      ix = process_entry(z,231,135,12345,1,"etiket","NOM2","TP","Z")
      print *,'ix=',ix
      ix = process_entry(z,231,135,12345,1,"etiket","NOM2","TP","Z")
      print *,'ix=',ix
      deallocate(z)
    endif

    do i = curarg, arg_count
      call get_command_argument(i,filename,arg_len,status)
      inquire(FILE=trim(filename),EXIST=file_exists)
      if(.not. file_exists)then
        if(verbose > 0) print *,"ERROR: file '"//trim(filename)//"' not found"
        if(strict) call f_exit(1)
      else
        if(verbose > 2) print *,"INFO: processing file '"//trim(filename)//"'"
      endif
    enddo
    
  end program
