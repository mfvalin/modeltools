program re_tag_scale
  use ISO_C_BINDING
  implicit none
  type :: retags
    character(len=4) :: o_name    ! original variable name
    character(len=4) :: n_name    ! new variable name
    real    :: sf                 ! scale factor
    real *8 :: ts                 ! time shift in hours
    integer :: ip1_old            ! if not -1, retag to ip1_new if ip1 = ip1_old
    integer :: ip1_new
  end type
  integer, parameter :: MAXTABLE=128
  type(retags), dimension(MAXTABLE) :: tags
  integer :: ntags
  integer :: nargs
  character(len=1024) :: old_file, new_file, dir_file
  integer :: i, status
  integer :: fstdin, fstdout
  integer, external :: fnom, fstouv, fstinf, fstsui, fstluk
  real,    dimension(:,:), allocatable :: array
  integer :: ni, nj, nk, nrec, realloc
  integer :: nia, nja
  integer :: nbits,datyp,ip1,ip2,ip3,dateo,deet,npas
  integer :: datyp_m
  integer :: new_dateo, datev
  integer :: new_ip1
  integer :: ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=1) :: grtyp
  character(len=4) :: nomvar, new_nomvar
  character(len=2) :: typvar
  character(len=12) :: etiket
  logical :: no_change
  real *8 :: deetnpas

  nargs = command_argument_count()
  if(nargs /= 3) call print_usage
  call GET_COMMAND_ARGUMENT(1,old_file)
  call GET_COMMAND_ARGUMENT(2,new_file)
  call GET_COMMAND_ARGUMENT(3,dir_file)

  print *,'u.re_tag_scale '//trim(old_file)//' '//trim(new_file)//' '//trim(dir_file)

  open(1,file=trim(dir_file),FORM='FORMATTED',status='OLD')
  print *,'INFO: contents of directive table'
  ntags = 0
  do i = 1 , MAXTABLE
    read(1,*,end=1) tags(i)
    print *,i,tags(i)
    ntags = i
  enddo
1 close(1)

  print *,'INFO: opening input file '//trim(old_file)
  fstdin = 0
  i = fnom(fstdin,trim(old_file),'STD+RND+OLD+R/O',0)
  print *,'DEBUG: fstdin=',fstdin
  i = fstouv(fstdin,'RND')
  print *,'DEBUG: status of fstouv fstdin =',i

  print *,'INFO: opening output file '//trim(new_file)
  fstdout = 0
  i = fnom(fstdout,trim(new_file),'STD+RND',0)
  print *,'DEBUG: fstdout=',fstdout
  i = fstouv(fstdout,'RND')
  print *,'DEBUG: status of fstouv fstdout =',i

  nrec = 0
  realloc = 0
  nia = 1
  nja = 1
  allocate(array(1,1))
  status = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","")
  do while(status >= 0)
    if(ni .ne. nia .or. nj .ne. nja) then
      deallocate(array)
      allocate(array(ni,nj))
      realloc = realloc + 1
      nia = ni
      nja = nj
    endif
    i = fstluk(array,status,ni,nj,nk)
    nrec = nrec + 1
!
!   transform the record as and if needed
!
    call fstprm(status,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
                typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                swa,lng,dltf,ubc,extra1,extra2,extra3)
    if(nbits > 32) goto 2   ! for now encoding with more than 32 bits is not supported
    deetnpas = deet
    deetnpas = deetnpas * npas
    no_change = .true.
    do i = 1 , ntags
      new_nomvar = tags(i)%o_name
      new_dateo = dateo
      new_ip1 = ip1
      if(tags(i)%o_name == nomvar) then
        print *,'DEBUG: transforming '//nomvar//' using ',tags(i)
        new_nomvar = tags(i)%n_name                                       ! new name
        if(tags(i)%ip1_old == -1 .or. tags(i)%ip1_old == ip1) then
          if(tags(i)%ip1_new .ne. -1) new_ip1 = tags(i)%ip1_new           ! new ip1
        endif
        datyp_m = mod(datyp,64)
        if(datyp_m == 1 .or. datyp_m == 5 .or. datyp_m == 6) then   ! can scale only 32 bit floating point data
          if(tags(i)%sf .ne. 1.0) array = array * tags(i)%sf
        else
          print *,'WARNING: cannot scale non 32 bit floating point data'
        endif
        call incdatr(new_dateo,dateo,tags(i)%ts)                          ! time tag shift
        call incdatr(datev,new_dateo,deetnpas)
        call fstecr(array,array,-nbits,fstdout,datev,deet,npas,ni,nj,nk,new_ip1,ip2,ip3,typvar,new_nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
        no_change = .false.
      endif
    enddo
    if(no_change) then
      call incdatr(datev,dateo,deetnpas)
      call fstecr(array,array,-nbits,fstdout,datev,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
    endif
2   continue
    status = fstsui(fstdin,ni,nj,nk)
  enddo
  print *,'INFO: number of records read and written =',nrec
  print *,'INFO: number of reallocations =',realloc

  call fstfrm(fstdin)
  call fstfrm(fstdout)
  stop
end program
subroutine print_usage()
  implicit none
  print *,'USAGE: u.re_tag_scale old_standard_file new_standard_file directives file'
  call qqexit(1)
end subroutine