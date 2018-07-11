#define VERSION '1.0_rc18 2017/10/24'
#define AVG_MARKER '/'
#define VAR_MARKER '%'
  module averages_common   ! tables and table management routines
    use iso_c_binding
    implicit none

    interface
      subroutine f_exit(code) BIND(C,name='exit')
        import :: C_INT
        integer(C_INT), intent(IN), value :: code
      end subroutine f_exit
    end interface

    type :: field
      integer*8 :: date_lo, date_hi                ! earliest date, latest date collected
      real *8, dimension(:,:,:), pointer :: stats  ! (:,:,1) = sum, (:,:,2) = sum of squares
      integer :: ni, nj                            ! field dimensions
      integer :: nsamples, sample                  ! number of samples collected, duration of sample
!       integer :: npas_min, npas_max                ! lowest and highest time step collected
      integer :: ip1, ip2, ip3, dateo, deet        ! from field's standard file parameters
      integer :: ig1, ig2, ig3, ig4                ! grid
      integer :: level2                            ! second level if 2 level data
      character(len=12) :: etiket                  ! from field's standard file parameters
      character(len=4)  :: nomvar                  ! from field's standard file parameters
      character(len=2)  :: typvar                  ! from field's standard file parameters
      character(len=1)  :: grtyp                   ! from field's standard file parameters
      logical :: special                           ! special record, ignore stats
    end type

    integer, parameter :: PAGE_SIZE  = 256   ! 256 fields per page in table (MUST BE A POWER OF 2)
    integer, parameter :: ENTRY_MASK = 255   ! PAGE_SIZE-1, mask to get slot number from index slot=and(ix,ENTRY_MASK)
    type :: page                             ! used to implement an array of pointers
      type(field), dimension(:), pointer :: p
    end type

    integer, parameter :: MAX_PAGES = 128   ! 128 pages at most (32768 records)
    integer, parameter :: PAGE_SHIFT = -8   ! right shift count to get page number from index
    type(page), dimension(:), pointer, save :: ptab => NULL()  ! statistics table, MAX_PAGES pages of size PAGE_SIZE

    integer, save :: next = -1             ! index of las valid record in statistics table
    integer, save :: verbose = 2           ! ERROR + WARNING
    logical, save :: variance = .true.     ! variance or standard deviation required [true unless -novar is used]
    integer, save :: fstdmean = 0          ! unit number for standard file into which averages are written
    integer, save :: fstdvar = 0           ! unit number for standard file into which variances/std deviatoins are written
    logical, save :: std_dev = .false.     ! output standard deviation rather than variance
    logical, save :: newtags = .false.     ! use new ip1/2/3 taggging style (not implemented yet)
    logical, save :: check_dateo = .false. ! check that all samples have the same date of origin (single model run)
    logical, save :: skip_npas0 = .true.   ! skip record if npas == 0 (default)
    logical, save :: weight_ip3 = .false.  ! use ip3 as a weight
    logical, save :: weight_time = .false. ! use time as a weight
    logical, save :: weight_abs = .false.  ! use a specific constant weight
    integer, save :: time_weight = 24      ! weight is in days, set to 1 for weight in hours
    logical, save :: strict = .false.      ! non strict mode by default
    logical, save :: select_etiket = .false. ! etiket 1s a significant item if .true.

    character(len=4), dimension(1024), save :: specials
    integer, save :: nspecials=0
  contains

    subroutine create_page(n)  ! create and initialize page n that will contain PAGE_SIZE entries
      implicit none
      integer, intent(IN) :: n  ! page number
      integer :: i, status

      if(n < 0 .or. n >= MAX_PAGES) then
        print *,"FATAL: page number",n,' is invalid. it must be between 0 and',MAX_PAGES-1
        call f_exit(1)
      endif
      if(associated(ptab(n)%p)) then
        print *,"FATAL: page number",n,' is already allocated'
        call f_exit(1)
      endif
      allocate(ptab(n)%p(0:PAGE_SIZE-1),STAT=status)  ! allocate page
      if(status .ne. 0) then
        print *,"FATAL: page number",n,' cannot be allocated'
        call f_exit(1)
      endif
      do i = 0, PAGE_SIZE-1
        nullify(ptab(n)%p(i)%stats)         ! nullify entries in page
      enddo
    end subroutine

    function create_page_entry(pg,slot,ni,nj) result(ix) ! create and initialize entry slot in page pg
      implicit none
      integer, intent(IN) :: pg, slot    ! page number, slot number
      integer, intent(IN) :: ni, nj      ! field dimensions (neede for stats)
      integer :: ix                      ! return value is index in table (ix < 0 indicates failure)
      integer :: status

      ix = -1
      if(slot < 0 .or. slot >= PAGE_SIZE ) then
        if(verbose > 0) print *,"ERROR: slot number",slot," not valid"
        return
      endif
      if(pg < 0 .or. pg >= MAX_PAGES ) then
        if(verbose > 0) print *,"ERROR: page number",pg," not valid"
        return
      endif

      if(.not. associated(ptab(pg)%p)) then  ! page does not exist, try to create it
        call create_page(pg)
      endif

      if(associated(ptab(pg)%p(slot)%stats)) return ! stats array exists, return with failure code

      ix = slot + PAGE_SIZE * pg                    ! build proper index number
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
!     initialize entry to EMPTY but READY
      ptab(pg)%p(slot)%stats = 0.0   ! set sum and sum of squares to zero
      ptab(pg)%p(slot)%ni = ni       ! dimensions of field
      ptab(pg)%p(slot)%nj = nj
      ptab(pg)%p(slot)%nsamples = 0  ! no samples yet
!       ptab(pg)%p(slot)%npas_min = 999999999
!       ptab(pg)%p(slot)%npas_max = -1
      ptab(pg)%p(slot)%dateo = -1
      ptab(pg)%p(slot)%date_lo = 0
      ptab(pg)%p(slot)%date_hi = 0
      ptab(pg)%p(slot)%sample = 0
      ptab(pg)%p(slot)%deet = -1
      ptab(pg)%p(slot)%ip1 = -1
      ptab(pg)%p(slot)%level2 = -1
      ptab(pg)%p(slot)%ip2 = -1
      ptab(pg)%p(slot)%ip3 = -1
      ptab(pg)%p(slot)%ig1 = -1
      ptab(pg)%p(slot)%ig2 = -1
      ptab(pg)%p(slot)%ig3 = -1
      ptab(pg)%p(slot)%ig4 = -1
      ptab(pg)%p(slot)%etiket = ""
      ptab(pg)%p(slot)%nomvar = ""
      ptab(pg)%p(slot)%typvar = ""
      ptab(pg)%p(slot)%grtyp = ""
      ptab(pg)%p(slot)%special = .false.   ! by default not a special record
      return
    end function create_page_entry

    function new_page_entry(ni,nj) result(ix) ! create a new entry in tables
      implicit none
      integer, intent(IN) :: ni, nj
      integer :: ix
      integer :: pg, slot

      next = next + 1             ! next valid index
      pg = next / PAGE_SIZE       ! page number
      slot = mod(next,PAGE_SIZE)  ! slot number
      ix = create_page_entry(pg,slot,ni,nj)
      if(ix < 0) then
        print *,"FATAL: error creating entry",slot,' in page',pg
        call f_exit(1)
      endif
      return
    end function new_page_entry

    function date_stamp_64(date_stamp) result(date_64)  ! compute date_64 from CMC date_stamp
      implicit none
      integer, intent(IN) :: date_stamp
      integer*8 :: date_64  ! in seconds

      integer :: t1, t2 
      integer year, month, day, hour, minute, second, julian

      call newdate(date_stamp,t1,t2,-3)  ! date stamp to t1(YYYYMMDD), t2(HHMMSShh)

      year = t1 / 10000
      t1 = t1 - (10000*year)
      month = t1 / 100
      day = t1 - (month*100)

      hour = t2 / 1000000
      t2 = t2 - (1000000*hour)
      minute = t2 / 10000
      t2 = t2 - (10000*minute)
      second = t2 / 100

      call jdatec(julian,year,month,day)
      date_64 = julian
!       date_64 = (date_64*86400) + (hour*3600) + (minute*60) + second     ! 86400 seconds in a day
      date_64 = date_64*86400
      date_64 = date_64 + hour*3600
      date_64 = date_64 + minute*60
      date_64 = date_64 + second
      return
    end function date_stamp_64

    function date_stamp_32(date_64) result(date_stamp)  ! convert date_64 to CMC datestamp
      implicit none
      integer*8, intent(IN) :: date_64  ! in seconds
      integer :: date_stamp

      integer :: t1, t2, t0  ! used later when implementing this function
      integer year, month, day, hour, minute, second, hms, ymd
      integer*8 :: ymd8

      ymd8 = date_64 / 86400    ! julian day
      ymd = ymd8
      call datec(ymd,year,month,day)
      t1 = year*10000 + month*100 + day   ! YYYYMMDD
      hms = date_64 -(ymd*86400)
      hour = hms / 3600
      hms = hms - (hour*3600)
      minute = hms / 60
      second = hms - (minute*60)
      t2 = hour*1000000 + minute*10000 + second*100  ! HHMMSS00
      call newdate(t0,t1,t2,3)  ! t1(YYYYMMDD), t2(HHMMSShh) to  date stamp
      date_stamp = t0
      return
    end function date_stamp_32
!
!   process record read from one of the input standard files
!
    function process_entry(z,ni,nj,ip1,ip2,ip3,dateo,deet,npas,etiket,nomvar,typvar,grtyp,ig1,ig2,ig3,ig4) result(ix)
      implicit none
      integer, intent(IN) :: ni, nj, ip1, ip2, ip3, npas, dateo, deet
      integer, intent(IN) :: ig1, ig2, ig3, ig4
      character(len=*), intent(IN) :: etiket, nomvar, typvar, grtyp
      real, dimension(ni,nj), intent(IN) :: z
      integer :: ix

      integer :: i, j, pg, slot, sample, it1, it2, it3, level2
      integer*8 :: date_lo, date_hi, dnp
      type(field), pointer :: p
      real :: weight, r1, r2, r3
      logical :: is_special
      character(len=128) :: string

      ix = -1
      level2 = -1  ! a priori, one level data
      call convip_plus( ip1, r1, it1, -1, string, .false. )    ! convert ip1
      if(ip3 == 0 .and. ip2 < 1000000) then    ! special "old style coding" for time in hours
        it2 = 10
        r2 = ip2
      else
	call convip_plus( ip2, r2, it2, -1, string, .false. )  ! convert ip2
      endif
      call convip_plus( ip3, r3, it3, -1, string, .false. )    ! convert ip3

      if( (it1 == it3) .and. (ip1 > 0) .and. (ip3 > 0) ) level2 = ip3   ! 2 level data
      if( (it1 == it2) .and. (ip1 > 0) .and. (ip2 > 0) ) level2 = ip2   ! 2 level data

      dnp = npas
      dnp = dnp * deet
      is_special = any(nomvar == specials(1:nspecials))
      weight = 1.0
      sample = 0  ! "instantaneous" sample
      if((it2 == 10) .and. (it3 == 10))then  ! both ip2 and ip3 are time tags
        sample = 3600*ABS(r2 - r3)           ! sample interval in seconds
        if(verbose > 4) print *,'DEBUG: sample =',sample
      endif
      if(is_special) then   !  special names
        date_lo = 0
        date_hi = 0
      else                  ! regular data / averages
        if( weight_ip3 .or. typvar(2:2) .eq. AVG_MARKER ) then ! weight is IP3 or number of samples in IP3
          i = ip3
          if(ishft(ip3,-24) == 15) i = iand(ip3,Z'000FFFFF')   ! keep lower 20 bits (type 15)
          weight = max(1,i)     ! number of samples
        endif
        if(weight_time) then                                 ! time weight
          weight = (dnp) / (3600.0 * time_weight)
        endif
        if(weight_abs) then
          weight = time_weight                  ! explicit weight
        endif
        date_lo = date_stamp_64(dateo)      ! compute 64 bit date_lo from dateo
        date_hi = date_lo + dnp             ! date of validity of sample
        if(weight == 1.0) then
          date_lo = date_hi                 ! one date kept (date of validity)
        else
          sample = (date_hi - date_lo)      ! interval in seconds between the 2 dates
          if(weight_ip3) sample = nint( (date_hi - date_lo) / (weight-1.0) )
        endif
      endif
      do i = 0 , next               ! do we have an entry that matches this record's parameters
        slot = iand(i,ENTRY_MASK)   ! slot from index
        pg = ishft(i,PAGE_SHIFT)    ! page from index
        p => ptab(pg)%p(slot)       ! pointer to data
        if(p%ip1 .ne. ip1 .or. p%level2 .ne. level2) cycle    ! not same level(s)
        if(p%ni .ne. ni .or. p%nj .ne. nj) cycle              ! not same dimensions
        if(trim(p%nomvar) .ne. trim(nomvar)) cycle            ! not same name
        if(trim(p%grtyp) .ne. trim(grtyp)) cycle              ! not same grid type
        if((p%ig1 .ne. ig1) .or. (p%ig2 .ne. ig2) .or. (p%ig3 .ne. ig3) .or. (p%ig4 .ne. ig4) ) cycle  ! not same grid
        if(trim(p%etiket) .ne. trim(etiket) .and. select_etiket) cycle            ! not same experiment
        if((p%dateo .ne. dateo) .and. check_dateo) cycle      ! dateo verification is optional
        if(is_special)then
          if(p%ip2 .ne. ip2 .or. p%ip3 .ne. ip3) cycle        ! special records must have same ip1/ip2/ip3
        endif
        if(p%typvar(2:2) .eq. AVG_MARKER) p%typvar = trim(typvar)  ! force typvar into p%typvar if average
        if(p%typvar(1:1) .ne. typvar(1:1)) cycle               ! check first character of typvar
        if(sample .ne. p%sample .and. weight == 1.0) then
           if(verbose > 1) print *,'WARNING: sample interval mismatch, got',sample,' expected',p%sample
           if(strict) call f_exit(1)    ! abort if strict mode
        endif
        ix = i                      ! a matching entry has been found
        if(verbose > 4) print *,'DEBUG: found entry, previous samples',ix,p%nsamples
        exit
      enddo

      if(ix == -1) then              ! not found, make a new entry
        ix = new_page_entry(ni,nj)
        slot = iand(ix,ENTRY_MASK)   ! slot from index
        pg = ishft(ix,PAGE_SHIFT)    ! page from index
        p => ptab(pg)%p(slot)        ! pointer to data
        p%etiket = trim(etiket)      ! initialize entry using parameters from record
        p%nomvar = trim(nomvar)
        p%typvar = trim(typvar)
        p%grtyp = trim(grtyp)
        p%dateo = dateo
        p%date_lo = date_lo
        p%date_hi = date_hi
        p%sample = sample
        p%deet = deet
        p%ip1 = ip1
        p%level2 = level2
        p%ip2 = ip2
        p%ip3 = ip3
        p%ig1 = ig1
        p%ig2 = ig2
        p%ig3 = ig3
        p%ig4 = ig4
        p%ni = ni
        p%nj = nj
        p%nsamples = 0
        p%special = is_special
!         print *,'DEBUG: dateo, deet,npas =',dateo, deet,npas,date_lo
      endif
      p%nsamples = p%nsamples + weight ! add weight to number of samples
!       if(p%nsamples == 48) print *,'DEBUG: dateo, deet,npas =',dateo, deet,npas,p%date_lo,date_hi
      if(is_special) return  ! no stats for special records, just add one to sample count
!       p%npas_max = max(p%npas_max,npas)
!       p%npas_min = min(p%npas_min,npas)
      p%date_lo = min(p%date_lo , date_lo) ! update earliest/latest date
      p%date_hi = max(p%date_hi , date_hi)
      do j = 1 , nj
      do i = 1 , ni
         p%stats(i,j,1) = p%stats(i,j,1) + z(i,j)*weight               ! update sum
         if(variance) p%stats(i,j,2) = p%stats(i,j,2) + z(i,j)*z(i,j)*weight  ! update sum of squares if necessary
      enddo
      enddo
    end function process_entry
  end module averages_common
!================================================================================================================
  subroutine print_usage(name)
    implicit none
    character(len=*) :: name
    print *,'USAGE: '//trim(name)//' [-h|--help] [-version] [-newtags] [-strict] [-novar] [-stddev] [-tag nomvar] \'
    print *,'           [-npas0] [-dateo] [-mean mean_out] [-var var_out] [-weight ip3|time|hours|days|nnn] \'
    print *,'           [-etiket] [-status path/to/status/file] [-test] [-q[q]] [-v[v][v]] [--|-f] \'
    print *,'           [mean_out] [var_out] in_1 ... in_n'
    print *,'        var_out  MUST NOT be present if -novar or -var is used'
    print *,'        mean_out MUST NOT be present if -mean is used'
    print *,'        -var -novar are mutually exclusive and may not be used together'
    print *,'        options are order independent but -- or -f MUST BE THE LAST ONE'
    print *,"        default special tag names = '>>  ', '^^  ', '!!  ', 'HY  '"
    print *,'        the -tag option may be used than once to add to this list'
    print *,'        etiket is ignored except if -etiket used on the command line'
    print *,'example :'
    print *,"  "//trim(name)//" -status stat.dot -vv -mean mean.fst -var var.fst -tag HY -tag '>>' my_dir/dm*"
    return
  end

  subroutine set_status(filename,message)
    implicit none
    character (len=*), intent(IN) :: filename,message
    integer :: iun, status
    integer, external :: fnom

    iun = 0
    status = fnom(iun,trim(filename),'FTN+FMT',0)
    if(iun > 0) then
      write(iun,1) trim(message)
1     format(A)
      call fclos(iun)
    endif
    return
  end subroutine set_status

  program averages
    use averages_common
    implicit none
    integer :: curarg
    character (len=8) :: option
    character (len=2048) :: filename, progname, meanfile, varfile, statusfile
    integer :: arg_count, arg_len, status, i
    integer :: first_file
    logical :: file_exists, test, missing
    real, dimension(:,:), pointer :: z
    type(field), pointer :: p
    integer :: ix, pg, slot, interval, expected
    integer, external :: process_file, write_stats

    print *,"Version "//VERSION       ! print version number
    statusfile = "/dev/null"
    meanfile = "/dev/null"
    varfile = "/dev/null"
    curarg = 1                         ! current argument number
    first_file = 0                     ! argument number of first input file
    strict = .false.                   ! do not abort on error (default)
    test = .false.                     ! internal flag for development purposes
    arg_count = command_argument_count()
    call get_command_argument(0,progname,arg_len,status)  ! get program name

    nspecials = 4
    specials(1) = ">>  "
    specials(2) = "^^  "
    specials(3) = "!!  "
    specials(4) = "HY  "

    if(arg_count < 1) then             ! no arguments, OUCH !!
      call print_usage(progname)
      call f_exit(1)
    endif

    do while(curarg <= arg_count)      ! process options
      call get_command_argument(curarg,option,arg_len,status)
      if(option(1:1) .ne. '-') exit      ! does not start with -, must be a file name (NO MORE OPTIONS)
      if(verbose > 3) print *,"NOTE: option = '"//trim(option)//"'"
      curarg = curarg + 1
      if(status .ne. 0) then
        print *,"FATAL: option is too long :'"//trim(option)//"..."//"'"
        call f_exit(1)
      endif

      if( option == '-h' .or. option == '--help' ) then
        call print_usage(progname)
        call f_exit(0)

      else if( option == '-version' ) then
        call f_exit(0)

      else if( option == '-strict' ) then   ! set strict mode
        strict = .true.                     ! abort on ERROR 

      else if( option == '-etiket' ) then
        select_etiket = .true.

      else if( option == '-weight' ) then     ! -mean file_for_averages
        if(curarg > arg_count) then
          print *,'FATAL: missing argument after -weight'
          call print_usage(progname)
          call f_exit(1)
        endif
        option = "none"
        weight_abs = .false.
        weight_time = .false.
        weight_ip3 = .false.
        call get_command_argument(curarg,option,arg_len,status) ! get weight option (ip3 or time)
        if(trim(option) == 'ip3') then
          weight_ip3 = .true.
        else if(trim(option) == 'hours') then
          weight_time = .true.
          time_weight = 1
        else if(trim(option) == 'days' .or. trim(option) == 'time') then
          weight_time = .true.
          time_weight = 24
        else if(trim(option) .ne. 'none') then
          weight_abs = .true.
          read(option,*) time_weight
        endif
        if(trim(option) .ne. 'none' .and. verbose > 2) print *,"INFO: WEIGHTING active, using '"//trim(option)//"'"
        curarg = curarg + 1

      else if( option == '-mean' ) then     ! -mean file_for_averages
        if(curarg > arg_count) then
          print *,'FATAL: missing argument after -mean'
          call print_usage(progname)
          call f_exit(1)
        endif
        call get_command_argument(curarg,meanfile,arg_len,status) ! get filename of averages file
        curarg = curarg + 1

      else if( option == '-status' ) then     ! -mean path/to/status/file
        if(curarg > arg_count) then
          print *,'FATAL: missing argument after -mean'
          call print_usage(progname)
          call f_exit(1)
        endif
        call get_command_argument(curarg,statusfile,arg_len,status) ! get filename of averages file
        curarg = curarg + 1

      else if( option == '-tag' ) then     ! -tag
        if(curarg > arg_count) then
          print *,'FATAL: missing argument after -tag'
          call print_usage(progname)
          call f_exit(1)
        endif
        call get_command_argument(curarg,option,arg_len,status) ! get filename of averages file
        curarg = curarg + 1
        nspecials = nspecials + 1
        specials(nspecials) = trim(option)
        if(verbose > 2) print *,"INFO: adding '"//specials(nspecials)//"' to special variable list"

      else if( option == '-var' ) then      ! -var file_for_variances
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

      else if( option == '-dateo' ) then
        check_dateo = .true.                ! check date of origin

      else if( option == '-npas0' ) then
        skip_npas0 = .false.                ! do not skip timestep 0

      else if( option == '-test' ) then
        test = .true.                   ! test mode (for development use only)

      else if( option == '-newtags' ) then
        newtags = .true.                   ! use new ip1/2/3 tagging style (not implemented yet)

      else if( option == '-stddev' ) then
        std_dev = .true.                   ! output standard deviation rather than variance

      else if( option == '-q' ) then
        if(.not. strict) verbose = 1    ! ERROR messages only

      else if( option == '-qq' ) then
        if(.not. strict) verbose = 0    ! NO messages at all

      else if( option == '-v' ) then
        verbose = 3                     ! + INFO messages

      else if( option == '-vv' ) then
        verbose = 4                     ! + NOTE messages

      else if( option == '-vvv' ) then
        verbose = 5                     ! + DEBUG messages

      else if( option == '--' .or. option == '-f' ) then  ! end of options, what follows is file names
        first_file = curarg
        if(verbose > 4) print *,'DEBUG: first file is argument',curarg
        exit

      else
        if(verbose > 1) print *,"WARNING: unrecognized option = '"//trim(option)//"'"
        if(strict) call f_exit(1)    ! abort if strict mode
        cycle
      endif
    enddo

    if( trim(statusfile) .ne. '/dev/null' ) call set_status(statusfile,'status="ABORT"')

    if(strict)  verbose = max(2 , verbose)          ! ERROR + WARNING messages at least

    allocate(ptab(0:MAX_PAGES-1))   ! allocate and initialize page table
    do i = 0 , MAX_PAGES-1          ! nullify all page pointers
      nullify(ptab(i)%p)
    enddo

    call create_page(0)

    if(verbose < 4) call fstopi("MSGLVL",4,.false.)  ! unless in debug mode, keep standard file package silent

    if(trim(meanfile) == "/dev/null") then   ! -mean option was not used
      call get_command_argument(curarg,meanfile,arg_len,status) ! get averages output file name
      curarg = curarg + 1
    endif
    inquire(FILE=trim(meanfile),EXIST=file_exists)
    if(file_exists) then
      if(verbose > 1) print *,"WARNING: mean output file '"//trim(meanfile)//" exists"
      if(strict) call f_exit(1)
    else
      if(verbose > 2) print *,"INFO: mean output file is '"//trim(meanfile)//"'"
    endif

    if(variance .and. (trim(varfile) == "/dev/null")) then     ! variance mode on and -var option was not used
      call get_command_argument(curarg,varfile,arg_len,status) ! get variance/std deviations output file name
      curarg = curarg + 1
    endif
    if(variance .and. (trim(varfile) .ne. "/dev/null")) then   ! a file name was given and -novar was not used
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
      call get_command_argument(i,filename,arg_len,status)   ! get file name
      inquire(FILE=trim(filename),EXIST=file_exists)         ! check that it exists
      if(.not. file_exists)then
        if(verbose > 0) print *,"ERROR: file '"//trim(filename)//"' not found"
        if(strict) call f_exit(1)   ! strict mode, error is fatal
      else
        if(verbose > 2) print *,"INFO: processing file '"//trim(filename)//"'"
        status = process_file(trim(filename))
        if(status<0 .and. strict) call f_exit(1)   ! strict mode, error is fatal
      endif
    enddo

    missing = .false.                ! missing samples flag
    do i = 0 , next                  ! check all valid entries in tables
      slot = iand(i,ENTRY_MASK)      ! slot from index
      pg = ishft(i,PAGE_SHIFT)       ! page from index
      p => ptab(pg)%p(slot)          ! pointer to entry
      if(p%special) cycle            ! ignore "special" records
      if(p%nsamples > 1) then        ! what follows is irrelevant if only one sample
!         interval = (p%npas_max - p%npas_min) / (p%nsamples - 1)  ! interval asssuming constant interval between samples 
!         expected = 1 + (p%npas_max - p%npas_min) / interval      ! expected number of samples given lowes/highest timestep numbers
        interval = (p%date_hi - p%date_lo) / (p%nsamples - 1)
        if(p%sample == 0) p%sample = interval
        if(interval == 0) then
	  expected = p%nsamples
        else
	  expected = 1 + (p%date_hi - p%date_lo) / interval
        endif
!        if(( p%npas_min + (p%nsamples - 1) * interval) .ne. p%npas_max) then
        if((expected - p%nsamples) .ne.0) then
          if(verbose > 1) print *,"WARNING:",expected - p%nsamples ," missing sample(s) for variable '"//p%nomvar//"'"
          missing = .true.
        endif
      endif
    enddo

    if(.not. (missing .and. strict) ) status = write_stats(meanfile,varfile)  ! write into statistics file(s)
111 continue
    call fstfrm(fstdmean)    ! close averages file
    call fclos(fstdmean)
    if(variance) then
      call fstfrm(fstdvar)   ! close variance/std deviation file
      call fclos(fstdvar)
    endif
    if(missing .and. strict) then
       print *,"FATAL: missing sample(s) for some variables"
      call f_exit(1)
    endif
    if( trim(statusfile) .ne. '/dev/null' ) call set_status(statusfile,'status="SUCCESS"')
  end program

  function process_file(filename) result(status)   ! process input standard file
    use averages_common
    implicit none
    character(len=*), intent(IN) :: filename
    integer :: status
    integer, external :: fnom, fstouv, fstinf, fstsui, process_record
    integer :: fstdin
    integer :: ni, nj, nk, key

    fstdin = 0
    status = fnom(fstdin,trim(filename),'STD+RND+OLD+R/O',0)   ! connect file
    if(status <0) return
    status = fstouv(fstdin,'RND')                              ! open file
    if(status <0) return
    if(verbose > 3) print *,"NOTE: opened file '"//trim(filename)//"'"

    key = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","")  ! first record
    do while(key >= 0)                                  ! loop while valid records
      if(nk > 1)    cycle                               ! xy 2D records only
      status = process_record(key,ni,nj,nk)             ! process record
      key = fstsui(fstdin,ni,nj,nk)                     ! next record 
    enddo
    call fstfrm(fstdin)                                 ! close file
    call fclos(fstdin)
  end function process_file

  function process_record(key,lni,lnj,lnk) result(status)  ! read record data and process it
    use averages_common
    implicit none
    integer, intent(IN) :: key, lni, lnj, lnk
    integer :: status

    integer, external :: fstluk
    real, dimension(lni,lnj,2) :: z      ! in case data type needs 64 bits
    integer :: ni, nj, nk
    integer :: nbits,datyp,ip1,ip2,ip3,dateo,deet,npas
    integer :: ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3
    character(len=1) :: grtyp
    character(len=4) :: nomvar
    character(len=2) :: typvar
    character(len=12) :: etiket
    integer :: dtyp, ix, slot, pg
    real, dimension(:,:), pointer :: z8
    type(field), pointer :: p

    status = 0
    call fstprm(key,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
                typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                swa,lng,dltf,ubc,extra1,extra2,extra3)

    if(any(nomvar == specials(1:nspecials))) then  ! special record, copy if first time seen, skip otherwise
      ix = process_entry(z,ni,nj,ip1,ip2,ip3,dateo,deet,npas,etiket,nomvar,typvar,grtyp,ig1,ig2,ig3,ig4)  ! process data
      slot = iand(ix,ENTRY_MASK)  ! slot from index
      pg = ishft(ix,PAGE_SHIFT)   ! page from index
      p => ptab(pg)%p(slot)       ! pointer to data
      if(p%nsamples > 1) then     ! we have already seen this record
        if(verbose > 3) print *,"NOTE: skipping non first occurrence of '"//p%nomvar//"'"
        return
      endif
!
!     first time around, copy to mean (and variance) file(s)
!
      status = fstluk(z,key,ni,nj,nk)            ! read record data
      if(status == -1) return                    ! ERROR
      if(nbits == 64) call fst_data_length(8)    ! if 64 bit data
      if(verbose > 4) print *,'DEBUG: special record - dateo', dateo
      call fstecr(z,z,-nbits,fstdmean,dateo,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
      if(variance) then  ! if there is a variance file, write it there too
        if(nbits == 64) call fst_data_length(8)    ! 64 bit data
        call fstecr(z,z,-nbits,fstdvar,dateo,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
      endif
      return
    endif   ! end of if special record

    if(datyp .ne. 5 .and. datyp .ne. 1 .and. datyp .ne. 6 .and. datyp .ne. 133 .and. datyp .ne. 134) then
      if(verbose > 2) print *,"INFO: skipping record '"//nomvar//"' (not floating point)"
      return
    endif

    if(typvar(2:2) == VAR_MARKER) then  ! cannot process variance records
      if(verbose > 2) print *,"INFO: skipping variance record"
      return
    endif

    if(typvar(2:2) .ne. " " .and.typvar(2:2) .ne.AVG_MARKER ) then  ! typvar(2:2) ignored except AVG_MARKER (averages)
      if(verbose > 3) print *,"NOTE: type '"//typvar//"' processed as type '"//typvar(1:1)//"'"
      typvar(2:2) = " "
    endif

    if(nbits > 32) then
      if(verbose > 3) print *,'NOTE: skipping record (nbits > 32)'
      return
    endif

    if(ni == 1 .or. nj== 1) then
      if(verbose > 3) print *,'NOTE: skipping record (1D record)'
      return
    endif

    if((npas == 0) .and. skip_npas0) then
      if(verbose > 3) print *,'NOTE: skipping record (npas==0)'
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

    status = fstluk(z,key,ni,nj,nk)            ! read record data
    if(status == -1) return                    ! ERROR

    status = process_entry(z,ni,nj,ip1,ip2,ip3,dateo,deet,npas,etiket,nomvar,typvar,grtyp,ig1,ig2,ig3,ig4)  ! process data

  end function process_record

  function write_stats(filename,varfile) result(status)   ! write statistics into statistics output file(s)
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
!    real *8 :: hours2
    integer, dimension(14) :: date_array
    integer :: new_dateo
    character(len=16) :: string
    real :: r4
    integer :: my_kind
    integer :: wtype

    if(fstdmean == 0 .and. fstdvar == 0) then  ! first call just opens the files

      status = fnom(fstdmean,trim(filename),'STD+RND',0) ! open averages file
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
    if(verbose > 2) print *,"INFO:",next+1," records will be written into mean/variance files"
    vartag = 'VA'               ! typvar for variance/std deviation file
    if(std_dev) vartag = 'ST'
    do i = 0 , next             ! loop over valid entries
      slot = iand(i,ENTRY_MASK) ! slot from index
      pg = ishft(i,PAGE_SHIFT)  ! page from index
      p => ptab(pg)%p(slot)     ! pointer to data
      if(p%special) then
        if(verbose > 2) print *,"INFO: ignoring special record '"//p%nomvar//"'"
        cycle
      endif
      if(verbose > 4) then
        print 100,"DEBUG: ",p%nsamples,p%nomvar,p%typvar,p%etiket,p%grtyp,p%ip1,associated(p%stats),p%ni,p%nj
100     format(A,I5,A5,A3,A13,A2,I10,L2,2I5,2I8)
      endif
      allocate(z(p%ni,p%nj))   ! allocate space for averages
      ov_sample = 1.0
      ov_sample = ov_sample / p%nsamples
      do jj = 1 , p%nj
      do ii = 1 , p%ni
        avg = p%stats(ii,jj,1) * ov_sample     ! average = sum / number of samples
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

      new_dateo = date_stamp_32(p%date_lo)
      datev = date_stamp_32(p%date_hi)
      hours = (p%date_hi - p%date_lo)
      hours = hours / 3600.0
      deet = 3600                      ! deet of written records forced to 1 hour
      npas = nint(hours)
      ip2 = npas                       ! to be adjusted if period does not start at 00Z
      date_array(14) = new_dateo
      call datmgp2(date_array)
      ip2 = ip2 + date_array(5)             ! zulu hour at start of period
      ip2 = ip2 + 24 * (date_array(3)-1)    ! force back to first day of month
      ip3 = p%nsamples                      ! number of samples
      ip1 = p%ip1
!       if(verbose > 2) print *,'INFO: ',p%nsamples,' '//p%nomvar//' "samples" every',p%sample/3600.0,' hours'
      if(verbose > 2) print *,'INFO: ',p%nsamples,' '//p%nomvar//' "samples"'
      if(newtags) then ! new tagging style (work in progress)
        if(p%sample > 0) then
          deet = p%sample
        else
          deet = (p%date_hi - p%date_lo) / (p%nsamples - 1)
        endif
	npas = nint(hours*3600/deet)
        r4 = ip3       ! force ip3 to new style coding
        call convip_plus( ip3, r4, 15, 2, string, .false. ) ! ip kind 15,  number of samples
        if( p%level2 > 0 ) then   ! 2 level data, put level2 in ip2 for output instead of hours since beginning of integration
          ip2 = p%level2
        else
	  call difdatr(new_dateo,p%dateo,hours)   ! ip2 = hours from start of model run
	  r4 = hours       ! force ip2 to new style coding
	  call convip_plus( ip2, r4, 10, 2, string, .false. ) ! ip kind 10,  hours
        endif
      endif
      p%typvar(2:2) = AVG_MARKER
      wtype = 5   ! E32
      if(p%ni > 16 .and. p%nj > 16) wtype = wtype + 128   ! activate E32 compression for averages
      call fstecr(z,z,-32,fstdmean, &
                  new_dateo,deet,npas,p%ni,p%nj,1,ip1,ip2,ip3,p%typvar,p%nomvar,p%etiket, &
                  p%grtyp,p%ig1,p%ig2,p%ig3,p%ig4,wtype,.false.)
      if(variance) then            ! use IEEE 64 bit format for variances
        call fst_data_length(8)    ! 64 bit format
        p%typvar(2:2) = VAR_MARKER
        call fstecr(p%stats(1,1,2),p%stats(1,1,2),-64,p%typvar, &
                    new_dateo,deet,npas,p%ni,p%nj,1,ip1,ip2,ip3,vartag,p%nomvar,p%etiket, &
                    p%grtyp,p%ig1,p%ig2,p%ig3,p%ig4,5,.false.)
      endif

      deallocate(z)
    enddo
    return
  end function write_stats
