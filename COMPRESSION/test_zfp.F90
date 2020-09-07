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
  integer :: debug, status, iun, nrec, key, ni, nj, nk
  integer, external :: fstouv, fnom, fstnbr, fstinf, fstsui
  real, dimension(:,:), pointer :: zi=>NULL()
  real, dimension(:,:), pointer :: zo=>NULL()
  real, dimension(:,:,:), pointer :: zt=>NULL()
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
      if(trim(nomvar) .eq. "ES") small = .05
      if(trim(nomvar) .eq. "TT") small = .01
      if(trim(nomvar) .eq. "UU") small = .01
      if(trim(nomvar) .eq. "VV") small = .01
      if(trim(nomvar) .eq. "GZ") small = .1
      if(associated(zi)) then
        deallocate(zi)
        deallocate(zo)
        deallocate(zt)
      endif
      allocate(zi(ni,nj), zo(ni,nj), zt(4,ni/2,nj/2))
      call fstluk(zi,key,ni,nj,nk)   ! read record
      zt(1,:,:) = zi(1:ni:2,1:nj:2)
      zt(2,:,:) = zi(2:ni:2,1:nj:2)
      zt(3,:,:) = zi(1:ni:2,2:nj:2)
      zt(4,:,:) = zi(2:ni:2,2:nj:2)
      buf = ZfpCompress3D(loc(zi), ni, nj, 1, toler, 0_8, streamsize)  ! compress
      dum = ZfpCompress3D(loc(zo), ni, nj, 1, 0.0_8, buf, streamsize)  ! decompress
      call AnalyzeErrors(loc(zi), loc(zo) , ni*nj, small, nomvar//achar(0))          ! evaluate
      call AnalyzeFields(loc(zi), loc(zo) , ni, ni, nj, small, stats)
      zo = zi
      call F_CDF97_2D_split_inplace_n(zo, ni, ni, nj, 3)
      call AnalyzeFields(loc(zi(1,1)), loc(zo(1,1)) , ni/2, ni, nj/2, small, stats)
      call AnalyzeFields(loc(zi(1+ni/2,1)), loc(zo(1+ni/2,1)) , ni/2, ni, nj/2, small, stats)
      call AnalyzeFields(loc(zi(1,1+nj/2)), loc(zo(1,1+nj/2)) , ni/2, ni, nj/2, small, stats)
      call AnalyzeFields(loc(zi(1+ni/2,1+nj/2)), loc(zo(1+ni/2,1+nj/2)) , ni/2, ni, nj/2, small, stats)
      write(6,*)'------------------------------------------------------------'
1     key = fstsui(iun,ni,nj,nk)
    enddo
  enddo
  goto 888

  stop
777 continue
  stop
888 continue
  call fstfrm(iun)
  stop
end program

