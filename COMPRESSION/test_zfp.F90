!   useful routines for C and FORTRAN programming
!   Copyright (C) 2020  Environnement Canada
! 
!   This is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Lesser General Public
!   License as published by the Free Software Foundation,
!   version 2.1 of the License.
! 
!   This software is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Lesser General Public License for more details.

program test_zfp
  use ISO_C_BINDING
  implicit none
  include 'model_tools.inc'

  integer, parameter :: nx=1000
  integer, parameter :: ny=1000
  integer, parameter :: nz=1
  real, dimension(nx,ny,nz) :: x
  integer(C_INTPTR_T) :: buf, dum
  integer(C_INT) :: streamsize
  integer :: i, j, k, l, nargs
  real a, b, c
  character(len=128) :: str, filename
  real(C_DOUBLE) :: toler
  integer :: debug, status, iun, nrec, key, ni, nj, nk, iunout
  integer, external :: fstouv, fnom, fstnbr, fstinf, fstsui
  real, dimension(:,:), pointer :: zi=>NULL()
  real, dimension(:,:), pointer :: zo=>NULL()
!   real, dimension(:,:,:), pointer :: zt=>NULL()
  integer :: date,deet,npas,nbits,datyp,ip1,ip2,ip3,ig1,ig2,ig3,ig4
  integer :: swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=2) :: typvar
  character(len=4) :: nomvar
  character(len=12) :: etiket
  character(len=1) :: grtyp
  real(C_FLOAT) :: small
  real(C_FLOAT), dimension(4) :: stats

  nargs = command_argument_count()
  do k = 1, nz
  do j = 1, ny
  do i = 1, nx
    a = 2.0 * (i-1) / nx
    b = 2.0 * (j-1) / ny
    c = 2.0 * (k-1) / nz
!     x(i,j,k) = exp(-(a * a + b * b + c * c))
    x(i,j,k) = (a * a + b * b + c * c) + .002
  enddo
  enddo
  enddo
  debug = ZfpCompress3D_debug(1)
  call get_command_argument(1, filename)
  iunout = 0
  status = fnom(iunout,'diag_file','RND+STD',0)
  if(status < 0) goto 777
  status = fstouv(iunout,'RND')
  if(status < 0) goto 777

  iun = 0
  status = fnom(iun,trim(filename),'RND+STD+R/O+OLD',0)
  if(status < 0) goto 777
  status = fstouv(iun,'RND')
  if(status < 0) goto 777
  nrec = fstnbr(iun)
  write(6,*) 'FILE='//trim(filename),' nrec =',nrec

  do l=2,nargs
    call get_command_argument(l, str)  ! each following argument is a precision/rate value
    read(str,*) toler                  ! <0 means rate in bits, float > 0 means tolerable error
    key = fstinf(iun,ni,nj,nk,-1,'            ',-1,-1,-1,'  ','   ')
    do while(key >= 0)
      if(ni <10 .or. nj < 10) goto 1
      call fstprm(key,date,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3,  &
                  typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                  swa,lng,dltf,ubc,extra1,extra2,extra3)
!       if(trim(nomvar) .ne. "ES") goto 1
!       if(trim(nomvar) .eq. "VV") goto 1
!       if(trim(nomvar) .eq. "UU") goto 1
!       if(trim(nomvar) .eq. "TT") goto 1
      if(trim(nomvar) .eq. "ES") small = .05
      if(trim(nomvar) .eq. "TT") small = .01
      if(trim(nomvar) .eq. "UU") small = .01
      if(trim(nomvar) .eq. "VV") small = .01
      if(trim(nomvar) .eq. "GZ") small = .1
      if(associated(zi)) then
        deallocate(zi)
        deallocate(zo)
      endif
      allocate(zi(ni,nj), zo(ni,nj))
!       allocate(zt(4,ni/2,nj/2))
      call fstluk(zi,key,ni,nj,nk)   ! read record
!       zt(1,:,:) = zi(1:ni:2,1:nj:2)
!       zt(2,:,:) = zi(2:ni:2,1:nj:2)
!       zt(3,:,:) = zi(1:ni:2,2:nj:2)
!       zt(4,:,:) = zi(2:ni:2,2:nj:2)
      buf = ZfpCompress3D(loc(zi), ni, nj, 1, toler, 0_8, streamsize)  ! compress
      dum = ZfpCompress3D(loc(zo), ni, nj, 1, 0.0_8, buf, streamsize)  ! decompress
      call CompareFields(loc(zi), loc(zo) , ni*nj, small, nomvar//achar(0))          ! evaluate
!       if(ip1 == 500) goto 1
!       call AnalyzeFields(loc(zi), loc(zo) , ni, ni, nj, small, stats)
      zo = zi
      nj = (nj/8) * 8                ! adjust nj to multiple of 8
      call F_CDF97_2D_split_inplace_n(zo, ni, ni, nj, 3)
      call AnalyzeField(loc(zi), ni, ni, nj, small, stats)                       ! original field
      call AnalyzeField(loc(zo), ni, ni, nj, small, stats)                       ! transformed field

      call AnalyzeField(loc(zo(1,1))           , ni/8, ni, nj/8, small, stats)   ! LLLLLL quadrant
      call AnalyzeField(loc(zo(1+ni/8,1))      , ni/8, ni, nj/8, small, stats)   ! LLLLHL quadrant
      call AnalyzeField(loc(zo(1,1+ni/8))      , ni/8, ni, nj/8, small, stats)   ! LLLLLH quadrant
      call AnalyzeField(loc(zo(1+ni/8,1+ni/8)) , ni/8, ni, nj/8, small, stats)   ! LLLLHH quadrant
      print *,""
      call AnalyzeField(loc(zo(1+ni/4,1))      , ni/4, ni, nj/4, small, stats)   ! LLHL quadrant
      call AnalyzeField(loc(zo(1,1+ni/4))      , ni/4, ni, nj/4, small, stats)   ! LLLH quadrant
      call AnalyzeField(loc(zo(1+ni/4,1+ni/4)) , ni/4, ni, nj/4, small, stats)   ! LLHH quadrant
      print *,""
      call AnalyzeField(loc(zo(1+ni/2,1))      , ni/2, ni, nj/2, small, stats)   ! HL quadrant
      call AnalyzeField(loc(zo(1,1+nj/2))      , ni/2, ni, nj/2, small, stats)   ! LH quadrant
      call AnalyzeField(loc(zo(1+ni/2,1+nj/2)) , ni/2, ni, nj/2, small, stats)   ! HH quadrant
      write(6,*)'------------------------------------------------------------'
      if(trim(nomvar) .eq. "UU") nomvar = "VEUU"
      if(trim(nomvar) .eq. "VV") nomvar = "VEVV"
      call fstecr(zo(1:ni/8     ,1:nj/8)     ,zi,-nbits,iunout,date,deet,npas,ni/8,nj/8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLLL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+ni/8:ni/4,1:nj/8)     ,zi,-nbits,iunout,date,deet,npas,ni/8,nj/8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLHL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1:ni/8     ,1+nj/8:nj/4),zi,-nbits,iunout,date,deet,npas,ni/8,nj/8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLLH",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+ni/8:ni/4,1+nj/8:nj/4),zi,-nbits,iunout,date,deet,npas,ni/8,nj/8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLHH",'A',0,0,0,0,datyp,.false.)

      call fstecr(zo(1:ni/4     ,1:nj/4)     ,zi,-nbits,iunout,date,deet,npas,ni/4,nj/4,1,ip1,ip2,ip3,typvar,nomvar,"LLLL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+ni/4:ni/2,1:nj/4)     ,zi,-nbits,iunout,date,deet,npas,ni/4,nj/4,1,ip1,ip2,ip3,typvar,nomvar,"LLHL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1:ni/4     ,1+nj/4:nj/2),zi,-nbits,iunout,date,deet,npas,ni/4,nj/4,1,ip1,ip2,ip3,typvar,nomvar,"LLLH",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+ni/4:ni/2,1+nj/4:nj/2),zi,-nbits,iunout,date,deet,npas,ni/4,nj/4,1,ip1,ip2,ip3,typvar,nomvar,"LLHH",'A',0,0,0,0,datyp,.false.)

      call fstecr(zo(1:ni/2     ,1:nj/2)     ,zi,-nbits,iunout,date,deet,npas,ni/2,nj/2,1,ip1,ip2,ip3,typvar,nomvar,"LL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+ni/2:ni/1,1:nj/2)     ,zi,-nbits,iunout,date,deet,npas,ni/2,nj/2,1,ip1,ip2,ip3,typvar,nomvar,"HL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1:ni/2     ,1+nj/2:nj/1),zi,-nbits,iunout,date,deet,npas,ni/2,nj/2,1,ip1,ip2,ip3,typvar,nomvar,"LH",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+ni/2:ni/1,1+nj/2:nj/1),zi,-nbits,iunout,date,deet,npas,ni/2,nj/2,1,ip1,ip2,ip3,typvar,nomvar,"HH",'A',0,0,0,0,datyp,.false.)

      call fstecr(zi                         ,zi,-nbits,iunout,date,deet,npas,ni/1,nj/1,1,ip1,ip2,ip3,typvar,nomvar,"ORIGINAL",'A',0,0,0,0,datyp,.false.)
1     key = fstsui(iun,ni,nj,nk)
    enddo
  enddo
  goto 888

  stop
777 continue
  stop
888 continue
  call fstfrm(iun)
  call fstfrm(iunout)
  stop
end program

