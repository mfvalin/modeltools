program re_tag_scale
  use ISO_C_BINDING
  implicit none
  integer :: nargs
  character(len=1024) :: old_file, new_file, the_new_date
  integer :: i, status
  integer :: fstdin, fstdout
  integer, external :: fnom, fstouv, fstinf, fstsui, fstluk
  real,    dimension(:,:), allocatable :: array
  integer :: ni, nj, nk, nrec
  integer :: nbits,datyp,ip1,ip2,ip3,dateo,deet,npas
  integer :: datev, new_datev
  integer :: ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=1) :: grtyp
  character(len=4) :: nomvar
  character(len=2) :: typvar
  character(len=12) :: etiket

  nargs = command_argument_count()
  if(nargs /= 3) call print_usage
  call GET_COMMAND_ARGUMENT(1,old_file)
  call GET_COMMAND_ARGUMENT(2,new_file)
  call GET_COMMAND_ARGUMENT(3,the_new_date)

  print *,'u.re_tag_scale '//trim(old_file)//' '//trim(new_file)//' '//trim(the_new_date)
  read(the_new_date,*)new_datev

  print *,'INFO: new date of validity: '//the_new_date,new_datev

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
  allocate(array(1000,1000))
  status = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","")
  do while(status >= 0)
    i = fstluk(array,status,ni,nj,nk)
    nrec = nrec + 1
!
!   transform the record as and if needed
!
    call fstprm(status,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
                typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                swa,lng,dltf,ubc,extra1,extra2,extra3)
    deet = 0
    npas = 0
    datev = new_datev
    call fstecr(array,array,-nbits,fstdout,datev,deet,npas,ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
    status = fstsui(fstdin,ni,nj,nk)
  enddo
  print *,'INFO: number of records read and written =',nrec

  call fstfrm(fstdin)
  call fstfrm(fstdout)
  stop
end program
subroutine print_usage()
  implicit none
  print *,'USAGE: u.re_tag_date old_standard_file new_standard_file new_date'
  call qqexit(1)
end subroutine