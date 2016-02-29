#define VERSION '1.0_rc3'
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
      integer :: ip1, dateo, deet
      integer :: ig1, ig2, ig3, ig4
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
    integer, parameter :: PAGE_SHIFT = -8   ! right shift count to get page number from index
    type(page), dimension(:), pointer, save :: ptab => NULL()

    integer, save :: next = -1
    integer, save :: verbose = 2           ! ERROR + WARNING
    logical, save :: variance = .true.
    logical, save :: firstfile = .true.
    integer, save :: fstdmean = 0
    integer, save :: fstdvar = 0
    logical, save :: std_dev = .false.     ! output standard deviation rather than variance
    logical, save :: newtags = .false.     ! use old ip1/2/3 taggging style by default

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
      if(variance) then ! need sum and sum of squares
        allocate(ptab(pg)%p(slot)%stats(ni,nj,2),STAT=status)
      else              ! only need sum
        allocate(ptab(pg)%p(slot)%stats(ni,nj,1),STAT=status)
      endif

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
      ptab(pg)%p(slot)%dateo = -1
      ptab(pg)%p(slot)%deet = -1
      ptab(pg)%p(slot)%ip1 = -1
      ptab(pg)%p(slot)%ig1 = -1
      ptab(pg)%p(slot)%ig2 = -1
      ptab(pg)%p(slot)%ig3 = -1
      ptab(pg)%p(slot)%ig4 = -1
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

    function process_entry(z,ni,nj,ip1,dateo,deet,npas,etiket,nomvar,typvar,grtyp,ig1,ig2,ig3,ig4) result(ix)
      implicit none
      integer, intent(IN) :: ni, nj, ip1, npas, dateo, deet
      integer, intent(IN) :: ig1, ig2, ig3, ig4
      character(len=*), intent(IN) :: etiket, nomvar, typvar, grtyp
      real, dimension(ni,nj), intent(IN) :: z
      integer :: ix
      integer :: i, j, pg, slot
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
        if(p%dateo .ne. dateo) cycle
        if(p%ip1 .ne. ip1) cycle
        if(p%ni .ne. ni) cycle
        if(p%nj .ne. nj) cycle
        ix = i                      ! entry has been found
        if(verbose > 4) print *,'DEBUG: found entry, previous samples',ix,p%nsamples
        exit
      enddo

      if(ix == -1) then
        ix = new_page_entry(ni,nj)
        slot = iand(ix,ENTRY_MASK)
        pg = ishft(ix,PAGE_SHIFT)
        p => ptab(pg)%p(slot)
        p%etiket = trim(etiket)
        p%nomvar = trim(nomvar)
        p%typvar = trim(typvar)
        p%grtyp = trim(grtyp)
        p%dateo = dateo
        p%deet = deet
        p%ip1 = ip1
        p%ig1 = ig1
        p%ig2 = ig2
        p%ig3 = ig3
        p%ig4 = ig4
        p%ni = ni
        p%nj = nj
      endif
      p%npas_max = max(p%npas_max,npas)
      p%npas_min = min(p%npas_min,npas)
      p%nsamples = p%nsamples + 1
!      p%stats(1:ni,1:nj,1) = p%stats(1:ni,1:nj,1) + z(1:ni,1:nj)
!      if(variance) p%stats(1:ni,1:nj,2) = p%stats(1:ni,1:nj,2) + z(1:ni,1:nj)*z(1:ni,1:nj)
      do j = 1 , nj
      do i = 1 , ni
         p%stats(i,j,1) = p%stats(i,j,1) + z(i,j)
         if(variance) p%stats(i,j,2) = p%stats(i,j,2) + z(i,j)*z(i,j)
      enddo
      enddo
    end function process_entry
  end module averages_common

  subroutine print_usage(name)
    implicit none
    character(len=*) :: name
    print *,'USAGE: '//trim(name)//' [-h|--help] [-newtags] [-strict] [-novar] [-stddev] \'
    print *,'           [-mean mean_out] [-var var_out] [-test] [-q[q]] [-v[v][v]] [--|-f] \'
    print *,'           [mean_out] [var_out] in_1 ... in_n'
    print *,'        var_out  MUST NOT be present if -novar or -var is used'
    print *,'        mean_out MUST NOT be present if -mean is used'
    print *,'        -var -novar are mutually exclusive and may not be used together'
    return
  end
  program averages
    use averages_common
    implicit none
    integer :: curarg
    character (len=8) :: option
    character (len=2048) :: filename, progname, meanfile, varfile
    integer :: arg_count, arg_len, status, i
    integer :: first_file
    logical :: file_exists, strict, test, missing
    real, dimension(:,:), pointer :: z
    type(field), pointer :: p
    integer :: ix, pg, slot, interval, expected
    integer, external :: process_file, write_stats

    print *,"Version "//VERSION
    meanfile = "/dev/null"
    varfile = "/dev/null"
    curarg = 1
    first_file = 0
    strict = .false.                   ! do not abort on error
    test = .false.
    arg_count = command_argument_count()
    call get_command_argument(0,progname,arg_len,status)
    if(arg_count < 1) then
      call print_usage(progname)
      call f_exit(1)
    endif
    do while(curarg <= arg_count)
      call get_command_argument(curarg,option,arg_len,status)
print *,'DEBUG: option=',option
      if(option(1:1) .ne. '-') exit      ! does not start with -, must be a file name
      curarg = curarg + 1
      if(status .ne. 0) then
        print *,"FATAL: option is too long :'"//trim(option)//"..."//"'"
        call f_exit(1)
      endif
      if( option == '-h' .or. option == '--help' ) then
        call print_usage(progname)
!        print *,'USAGE: '//trim(progname)//' [-h|--help] [-newtags] [-strict] [-novar] [-stddev] [-test] [-q[q]] [-v[v][v]] [--|-f] mean_out [var_out] in_1 ... in_n'
        call f_exit(0)
      else if( option == '-strict' ) then
        strict = .true.                 ! abort on ERROR 
        verbose = 1
      else if( option == '-mean' ) then
        if(curarg > arg_count) then
          print *,'FATAL: missing argument after -mean'
          call print_usage(progname)
          call f_exit(1)
        endif
        call get_command_argument(curarg,meanfile,arg_len,status) ! get filename of means file
        curarg = curarg + 1
      else if( option == '-var' ) then
        if(.not. variance) then  ! -var and -novar are mutually exclusive
          call print_usage(progname)
          call f_exit(1)
        endif
        if(curarg > arg_count) then
          print *,'FATAL: missing argument after -var'
          call print_usage(progname)
          call f_exit(1)
        endif
        call get_command_argument(curarg,varfile,arg_len,status) ! get filename of variances/stddevs file
        curarg = curarg + 1
      else if( option == '-novar' ) then
        variance = .false.                   ! test mode
        if(trim(varfile) .ne. "/dev/null") then  ! -var and -novar are mutually exclusive
          call print_usage(progname)
          call f_exit(1)
        endif
      else if( option == '-test' ) then
        test = .true.                   ! test mode
      else if( option == '-newtags' ) then
        newtags = .true.                   ! new ip1/2/3 tagging style
      else if( option == '-stddev' ) then
        std_dev = .true.                   ! test mode
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
      ix = process_entry(z,135,231,12345,0,0,1,"etiket","NOM","TP","Z",0,0,0,0)
      print *,'ix=',ix
      deallocate(z)
      allocate(z(231,135))
      z = 1.1
      ix = process_entry(z,231,135,12345,0,0,1,"etiket","NOM","TP","Z",0,0,0,0)
      print *,'ix=',ix
      deallocate(z)
      allocate(z(135,231))
      z = 1.1
      ix = process_entry(z,135,231,12345,0,0,1,"etiket","NOM","TP","Z",0,0,0,0)
      print *,'ix=',ix
      deallocate(z)
      allocate(z(231,135))
      z = 1.1
      ix = process_entry(z,231,135,12345,0,0,1,"etiket","NOM","TP","Z",0,0,0,0)
      print *,'ix=',ix
      deallocate(z)
      allocate(z(231,135))
      z = 1.1
      ix = process_entry(z,231,135,12345,0,0,1,"etiket","NOM2","TP","Z",0,0,0,0)
      print *,'ix=',ix
      ix = process_entry(z,231,135,12345,0,0,1,"etiket","NOM2","TP","Z",0,0,0,0)
      print *,'ix=',ix
      ix = process_entry(z,231,135,12345,0,0,1,"etiket","NOM2","TP","Z",0,0,0,0)
      print *,'ix=',ix
      ix = process_entry(z,231,135,12345,0,0,1,"etiket","NOM2","TP","Z",0,0,0,0)
      print *,'ix=',ix
      deallocate(z)
    endif

    if(verbose < 4) call fstopi("MSGLVL",4,.false.)
    if(trim(meanfile) == "/dev/null") then 
      call get_command_argument(curarg,meanfile,arg_len,status) ! first file name is mean output file
      curarg = curarg + 1
    endif
    inquire(FILE=trim(meanfile),EXIST=file_exists)
    if(file_exists) then
      if(verbose > 1) print *,"WARNING: mean output file '"//trim(meanfile)//" exists"
      if(strict) call f_exit(1)
    else
      if(verbose > 2) print *,"INFO: mean output file is '"//trim(meanfile)//"'"
    endif

    if(variance .and. (trim(varfile) == "/dev/null")) then
      call get_command_argument(curarg,varfile,arg_len,status) ! second file name is variance output file
      curarg = curarg + 1
    endif
    if(variance .and. (trim(varfile) .ne. "/dev/null")) then
      inquire(FILE=trim(varfile),EXIST=file_exists)
      if(file_exists) then
        if(verbose > 1) print *,"WARNING: variance output file '"//trim(varfile)//" exists"
        if(strict) call f_exit(1)
      else
        if(verbose > 2) print *,"INFO: variance output file is '"//trim(varfile)//"'"
      endif
    endif

    status = write_stats(meanfile,varfile)   ! first call just opens the files

    do i = curarg, arg_count         ! loop over input files
      call get_command_argument(i,filename,arg_len,status)
      inquire(FILE=trim(filename),EXIST=file_exists)
      if(.not. file_exists)then
        if(verbose > 0) print *,"ERROR: file '"//trim(filename)//"' not found"
        if(strict) call f_exit(1)
      else
        if(verbose > 2) print *,"INFO: processing file '"//trim(filename)//"'"
        if(.not. test) status = process_file(trim(filename))
      endif
      firstfile = .false.            ! after first file
    enddo

    missing = .false.
    do i = 0 , next
      slot = iand(i,ENTRY_MASK)
      pg = ishft(i,PAGE_SHIFT)
      p => ptab(pg)%p(slot)
      if(p%nsamples > 1) then
        interval = (p%npas_max - p%npas_min) / (p%nsamples - 1)
        expected = 1 + (p%npas_max - p%npas_min) / interval
!        if(( p%npas_min + (p%nsamples - 1) * interval) .ne. p%npas_max) then
        if((expected - p%nsamples) .ne.0) then
          if(verbose > 1) print *,"WARNING:",expected - p%nsamples ," missing sample(s) for variable '"//p%nomvar//"'"
          missing = .true.
        endif
      endif
    enddo

    if(.not. (missing .and. strict) ) status = write_stats(meanfile,varfile)
    call fstfrm(fstdmean)
    call fclos(fstdmean)
    if(variance) then
      call fstfrm(fstdvar)
      call fclos(fstdvar)
    endif
    if(missing .and. strict) then
       print *,"FATAL: missing sample(s) for some variables"
      call f_exit(1)
    endif
  end program
  function process_file(filename) result(status)
    use averages_common
    implicit none
    character(len=*), intent(IN) :: filename
    integer :: status
    integer, external :: fnom, fstouv, fstinf, fstsui, process_record
    integer :: fstdin
    integer :: ni, nj, nk, key

    fstdin = 0
    status = fnom(fstdin,trim(filename),'STD+RND+OLD+R/O',0)
    if(status <0) return
    status = fstouv(fstdin,'RND')
    if(status <0) return
    if(verbose > 3) print *,"NOTE: opened file '"//trim(filename)//"'"
    key = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","")
    do while(key >= 0)
      if(nk > 1)    cycle
      status = process_record(key,ni,nj,nk)
      key = fstsui(fstdin,ni,nj,nk)
    enddo
    call fstfrm(fstdin)
    call fclos(fstdin)
  end function process_file

  function process_record(key,lni,lnj,lnk) result(status)
    use averages_common
    implicit none
    integer, intent(IN) :: key, lni, lnj, lnk
    integer :: status
    integer, external :: fstluk
    real, dimension(lni,lnj) :: z
    integer :: ni, nj, nk
    integer :: nbits,datyp,ip1,ip2,ip3,dateo,deet,npas
    integer :: datev
    integer :: ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3
    character(len=1) :: grtyp
    character(len=4) :: nomvar
    character(len=2) :: typvar
    character(len=12) :: etiket
    integer :: dtyp
    real, dimension(:,:), pointer :: z8
    real*8 :: hours

    status = 0
    call fstprm(key,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
                typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                swa,lng,dltf,ubc,extra1,extra2,extra3)
    hours = deet
    hours = hours / 3600.0
    hours = hours * npas
    call incdatr(datev,dateo,hours)
    if(nomvar == ">>  " .or. nomvar == "^^  " .or. nomvar == "HY  " .or. nomvar == "!!  ") then  ! special record
      if(.not. firstfile) return
      allocate(z8(ni,nj*2))
      status = fstluk(z8,key,lni,lnj,nk)         ! read record
      if(nbits == 64) call fst_data_length(8)    ! 64 bit data
      if(verbose > 4) print *,'DEBUG: special record - dateo, datev datev', dateo, datev, extra1
      call fstecr(z8,z8,-nbits,fstdmean,dateo,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
      if(variance) then  ! if there is a variance file, write it there too
        if(nbits == 64) call fst_data_length(8)    ! 64 bit data
        call fstecr(z8,z8,-nbits,fstdvar,dateo,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
      endif
      deallocate(z8)
      return
    endif
    if(ni == 1 .or. nj== 1) then
      if(verbose > 4) print *,'DEBUG:skipping record (1D record)'
      return
    endif
!    if(npas == 0 .or. npas == 2205 .or. npas == 2178) then   ! test code for missign samples
    if(npas == 0) then
      if(verbose > 4) print *,'DEBUG:skipping record (npas==0)'
      return
    endif
    dtyp = mod(datyp,16)
    if(dtyp == 0 .or. dtyp == 2 .or. dtyp == 3 .or. &
       dtyp == 4 .or. dtyp == 7 .or. dtyp == 8) then
      if(verbose > 4) print *,'DEBUG:skipping record dtyp =', dtyp
      return  ! datatyp must be one of 1, 5, 6 (float)
    endif
    if(lni .ne. ni .or. lnj .ne. nj .or. lnk .ne. nk) then  ! this should NEVER happen
      if(verbose > 1) print *,'WARNING:skipping record - dimension mismatch (should never happen)'
      status = -1
      return
    endif
    status = fstluk(z,key,lni,lnj,nk)    ! read record
    status = process_entry(z,ni,nj,ip1,dateo,deet,npas,etiket,nomvar,typvar,grtyp,ig1,ig2,ig3,ig4)  ! process record
  end function process_record

  function write_stats(filename,varfile) result(status)
    use averages_common
    implicit none
    character(len=*), intent(IN) :: filename, varfile
    integer :: status
    integer, external :: fnom, fstouv
    integer :: i, ii, jj, slot, pg
    type(field), pointer :: p
    real *8 :: avg, var
    real, dimension(:,:), pointer :: z
    integer :: deet, npas, ip1, ip2, ip3, datev
    character(len=2) :: vartag
    real *8 :: hours, hours_min, hours_max, ov_sample
    integer, dimension(14) :: date_array
    integer :: new_dateo

    if(fstdmean == 0 .and. fstdvar == 0) then  ! first call just opens the files

      status = fnom(fstdmean,trim(filename),'STD+RND',0)
      if(status <0) return
      status = fstouv(fstdmean,'RND')
      if(status <0) return
      if(verbose > 2) print *,"INFO: opened mean output file '"//trim(filename)//"'"

      if(variance) then  ! only open variance file if it is required
        status = fnom(fstdvar,trim(varfile),'STD+RND',0)
        if(status <0) return
        status = fstouv(fstdvar,'RND')
        if(status <0) return
        if(verbose > 2) print *,"INFO: opened variance output file '"//trim(filename)//"'"
      endif

      return
    endif

    status = 0
    if(verbose > 2) print *,"INFO:",next," records will be written into mean/variance files"

    vartag = 'VA'
    if(std_dev) vartag = 'ST'
    do i = 0 , next
      slot = iand(i,ENTRY_MASK)
      pg = ishft(i,PAGE_SHIFT)
      p => ptab(pg)%p(slot)
      if(verbose > 4) then
        print 100,"DEBUG: ",p%nsamples,p%nomvar,p%typvar,p%etiket,p%grtyp,p%ip1,associated(p%stats),p%ni,p%nj,p%npas_min,p%npas_max
100     format(A,I5,A5,A3,A13,A2,I10,L2,2I5,2I8)
      endif
      allocate(z(p%ni,p%nj))
      ov_sample = 1.0
      ov_sample = ov_sample / p%nsamples
      do jj = 1 , p%nj
      do ii = 1 , p%ni
        avg = p%stats(ii,jj,1) * ov_sample
        z(ii,jj) = real(avg)
        if(variance) then
          var = p%stats(ii,jj,2) * ov_sample - avg * avg
          var = max(var,0.0*var)
          if(std_dev) then   ! output standard deviation
            p%stats(ii,jj,2) = sqrt(var)
          else               ! output variance (default)
            p%stats(ii,jj,2) = var
          endif
        endif
      enddo
      enddo

      hours = p%deet
      hours = hours / 3600             ! p%deet in hours
      hours_min = hours * p%npas_min   ! start of period = p%dateo + hours_min
      call incdatr(new_dateo,p%dateo,hours_min) ! start of period timestamp
      hours_max = hours * p%npas_max   ! end   of period = p%dateo + hours_max
      call incdatr(datev,p%dateo,hours_max)  ! datev = p%dateo + hours_max
      hours = (p%npas_max - p%npas_min) * hours   ! span in hours of period
      deet = 3600                      ! deet of written records forced to 1 hour
      npas = nint(hours)
      ip2 = npas                       ! to be adjusted if period does not start at 00Z
      date_array(14) = new_dateo
      call datmgp2(date_array)
      ip2 = ip2 + date_array(5)             ! zulu hour at start of period
      ip2 = ip2 + 24 * (date_array(3)-1)    ! force back to first day of month
      ip3 = p%nsamples
      ip1 = p%ip1
      if(newtags) then ! new tagging style (this code is a placeholder and a NO-OP for now)
        ip1 = ip1      ! beginning time (related to hours_max)
        ip2 = ip2      ! ending time    (related to hours_min)
        ip3 = ip3      ! number of samples
        new_dateo = new_dateo  ! consistent with ip1 and ip2
      endif
!     call fstecr(z8,z8,-nbits,fstdmean, &
!                datev,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket, &
!                grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
      call fstecr(z,z,-32,fstdmean, &
                  new_dateo,deet,npas,p%ni,p%nj,1,ip1,ip2,ip3,"MN",p%nomvar,p%etiket, &
                  p%grtyp,p%ig1,p%ig2,p%ig3,p%ig4,128+5,.false.)
      if(variance) then
        call fst_data_length(8)
        call fstecr(p%stats(1,1,2),p%stats(1,1,2),-64,fstdvar, &
                    new_dateo,deet,npas,p%ni,p%nj,1,ip1,ip2,ip3,vartag,p%nomvar,p%etiket, &
                    p%grtyp,p%ig1,p%ig2,p%ig3,p%ig4,5,.false.)
      endif

      deallocate(z)
    enddo
    return
  end function write_stats
