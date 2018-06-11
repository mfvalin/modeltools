program re_tag_scale
  use ISO_C_BINDING
  implicit none
  integer, parameter :: MAXTABLE=128
  integer :: nargs
  character(len=1024) :: old_file, new_file, dir_file
  integer :: i, status, j
  integer :: fstdin, fstdout
  integer, external :: fnom, fstouv, fstinf, fstsui, fstluk
  real,    dimension(:), allocatable :: array
  real, dimension(:), allocatable :: latitudes, longitudes
  integer :: ni, nj, nk, nrec, realloc
  integer :: nia, nja
  integer :: nbits,datyp,ip1,ip2,ip3,dateo,deet,npas,datev
  integer :: ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=1) :: grtyp
  character(len=4) :: nomvar
  character(len=2) :: typvar
  character(len=12) :: etiket
  character(len=4), dimension(128) :: names
  integer, dimension(128) :: plist
  integer :: nnames, nlist, nlatitudes, nlongitudes
  real, dimension(128) :: values
  real *8 :: deetnpas
  integer, dimension(:,:), allocatable :: stations, snames

  nargs = command_argument_count()
  if(nargs /= 3) call print_usage
  call GET_COMMAND_ARGUMENT(1,old_file)
  call GET_COMMAND_ARGUMENT(2,new_file)
  call GET_COMMAND_ARGUMENT(3,dir_file)

  print *,'u.clip_grid '//trim(old_file)//' '//trim(new_file)//' '//trim(dir_file)

  plist = -1
  open(1,file=trim(dir_file),FORM='FORMATTED',status='OLD')
  print *,'INFO: contents of directive table'
  read(1,*)nlist,plist(1:nlist)    ! list of points to extract
  write(6,100)'copying points :',plist(1:nlist)
100 format(A,(7I5))
  read(1,*)nnames,names(1:nnames)    ! variable names
  write(6,101)'clipping variabels :',names(1:nnames)
101 format(A,16A5)
1 close(1)

  print *,'INFO: opening input file '//trim(old_file)
  fstdin = 0
  status = fnom(fstdin,trim(old_file),'STD+RND+OLD+R/O',0)
  if(status < 0) stop
  print *,'DEBUG: fstdin=',fstdin
  status = fstouv(fstdin,'RND')
  if(status < 0) stop
  print *,'DEBUG: status of fstouv fstdin =',status

  print *,'INFO: opening output file '//trim(new_file)
  fstdout = 0
  status = fnom(fstdout,trim(new_file),'STD+RND',0)
  if(status < 0) stop
  print *,'DEBUG: fstdout=',fstdout
  status = fstouv(fstdout,'RND')
  if(status < 0) stop
  print *,'DEBUG: status of fstouv fstdout =',status

! clip and copy longitudes, '>>'
  status = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"",">>")
  if(status < 0) stop
  call fstprm(status,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
              typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
              swa,lng,dltf,ubc,extra1,extra2,extra3)
  nlongitudes = ni*nj
  allocate(longitudes(nlongitudes))
  status = fstluk(longitudes,status,ni,nj,nk)  !  read data
  do i=1,nlist
    values(i) = longitudes(plist(i))
  enddo
  deetnpas = deet*npas
  call incdatr(datev,dateo,deetnpas)
  call fstecr(values,values,-nbits,fstdout,datev,deet,npas,nlist,1,1,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)

! clip and copy latitudes, '^^
  status = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","^^")
  if(status < 0) stop
  call fstprm(status,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
              typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
              swa,lng,dltf,ubc,extra1,extra2,extra3)
  nlatitudes = ni*nj
  if(nlatitudes .ne. nlongitudes) stop
  allocate(latitudes(nlatitudes))
  status = fstluk(latitudes,status,ni,nj,nk)  !  read data
  do i=1,nlist
    values(i) = latitudes(plist(i))
  enddo
  deetnpas = deet*npas
  call incdatr(datev,dateo,deetnpas)
  call fstecr(values,values,-nbits,fstdout,datev,deet,npas,nlist,1,1,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)

! clip and copy station names, 'STNS
  status = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","STNS")
  if(status < 0) stop
  allocate(stations(32,nj),snames(32,nlist))
  call fstprm(status,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
              typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
              swa,lng,dltf,ubc,extra1,extra2,extra3)
  status = fstluk(stations,status,ni,nj,nk)  !  read data
  if(status < 0) stop
  do i=1,nlist
    snames(:,i) = stations(:,plist(i))
  enddo
  deetnpas = deet*npas
  call incdatr(datev,dateo,deetnpas)
  call fstecr(snames,snames,-nbits,fstdout,datev,deet,npas,128,nlist,1,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
  

  nrec = 0
  realloc = 0
  nia = 1
  nja = 1
  allocate(array(1))
  status = fstinf(fstdin,ni,nj,nk,-1,"",-1,-1,-1,"","")
  do while(status >= 0)
    if(ni*nj > nia*nja) then
      deallocate(array)
      allocate(array(ni*nj))
      realloc = realloc + 1
      nia = ni
      nja = nj
    endif
    call fstprm(status,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3, &
                typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                swa,lng,dltf,ubc,extra1,extra2,extra3)
!
!   transform the record as and if needed
!
    do i = 1 , nnames
      if(names(i) == nomvar .and. ip3 ==0) then
        nrec = nrec + 1
        status = fstluk(array,status,ni,nj,nk)  !  read data
        print *,'DEBUG: transforming '//nomvar
        do j=1,nlist
          values(j) = array(plist(j))
        enddo
        deetnpas = deet*npas
        call incdatr(datev,dateo,deetnpas)
        call fstecr(values,values,-nbits,fstdout,datev,deet,npas,nlist,1,1,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.false.)
      endif
    enddo
2   continue
    status = fstsui(fstdin,ni,nj,nk)
  enddo   ! while status >= 0
  print *,'INFO: number of records processed =',nrec
  print *,'INFO: number of reallocations =',realloc

  call fstfrm(fstdin)
  call fstfrm(fstdout)
  stop
end program
subroutine print_usage()
  implicit none
  print *,'USAGE: u.clip_grid old_standard_file new_standard_file directives_file'
  call qqexit(1)
end subroutine