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
#define STR(fstr) fstr//achar(0)
program test_zfp
  use ISO_C_BINDING
  implicit none
#include <rpn_pack_header.hf>
  include 'model_tools.inc'

  integer, parameter :: nx=1000
  integer, parameter :: ny=1000
  integer, parameter :: nz=1
!   real, dimension(nx,ny,nz) :: x
  integer(C_INTPTR_T) :: buf, dum
  integer(C_INT) :: streamsize
  integer :: i, j, k, l, nargs
  real a, b, c
  character(len=128) :: str, filename
  real(C_DOUBLE) :: toler
  integer :: debug, status, iun, nrec, key, ni, nj, nk, iunout
  integer :: nie8, nie4, nie2, nje8, nje4, nje2, ni0
  integer :: nio8, nio4, nio2, njo8, njo4, njo2, nj0
  integer, external :: fstouv, fnom, fstnbr, fstinf, fstsui
  integer, dimension(:,:), pointer :: izi=>NULL()
  integer, dimension(:,:), pointer :: izo=>NULL()
  integer*2, dimension(:,:), pointer :: izo16=>NULL()
  real, dimension(:,:), pointer :: zi=>NULL()
  real, dimension(:,:), pointer :: zo=>NULL()
  real, dimension(:,:), pointer :: zt=>NULL()
  real, dimension(:,:), pointer :: pe=>NULL()    ! 2D Lorenzo prediction error
!   real, dimension(:,:,:), pointer :: zt=>NULL()
  integer :: date,deet,npas,nbits,datyp,ip1,ip2,ip3,ig1,ig2,ig3,ig4
  integer :: swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=2) :: typvar
  character(len=4) :: nomvar
  character(len=12) :: etiket
  character(len=1) :: grtyp
  real(C_FLOAT) :: small, small4, small8
  real(C_FLOAT), dimension(4) :: stats
  integer :: total_length
  real :: ratio
  integer(C_SIZE_T) :: nbytes, dstCapacity, srcSize, nbytes16, nbytesf
  interface
  function compress_cycle(z, ni, nj, nbits) result(l)
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN) :: ni, nj
    real(C_FLOAT), intent(INOUT), dimension(ni,nj), target :: z
    integer(C_INT), intent(INOUT) :: nbits
    integer(C_INT) :: l
  end function compress_cycle
!   size_t FSE_compress(void* dst, size_t dstCapacity, const void* src, size_t srcSize);
  function FSE_compress(dst, dstCapacity, src, srcSize) result(nbytes) bind(C,name='FSE_compress')
    import :: C_PTR, C_SIZE_T
    type(C_PTR), intent(IN), value :: src, dst
    integer(C_SIZE_T), intent(IN), value :: dstCapacity, srcSize
    integer(C_SIZE_T) :: nbytes
  end function FSE_compress
  function FSE_compress16(dst, dstCapacity, src, srcSize, maxsym, tbllog) result(nbytes) bind(C,name='FSE_compressU16')
    import :: C_PTR, C_SIZE_T, C_INT
    type(C_PTR), intent(IN), value :: src, dst
    integer(C_SIZE_T), intent(IN), value :: dstCapacity, srcSize
    integer(C_INT), intent(IN), value :: maxsym, tbllog
    integer(C_SIZE_T) :: nbytes
  end function FSE_compress16
  end interface

  nargs = command_argument_count()
!   do k = 1, nz
!   do j = 1, ny
!   do i = 1, nx
!     a = 2.0 * (i-1) / nx
!     b = 2.0 * (j-1) / ny
!     c = 2.0 * (k-1) / nz
! !     x(i,j,k) = exp(-(a * a + b * b + c * c))
!     x(i,j,k) = (a * a + b * b + c * c) + .002
!   enddo
!   enddo
!   enddo
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
!       if(trim(nomvar) .eq. "WW") goto 2
!       if(trim(nomvar) .eq. "TT") goto 2
!       if(trim(nomvar) .eq. "QC") goto 2
      if(trim(nomvar) .eq. "TT") goto 2
!       if(trim(nomvar) .eq. "HU") goto 2
      goto 1
!       if(trim(nomvar) .ne. "ES") goto 1
!       if(trim(nomvar) .eq. "VV") goto 1
!       if(trim(nomvar) .eq. "UU") goto 1
!       if(trim(nomvar) .ne. "TT") goto 1
!       if(trim(nomvar) .ne. "PN") goto 1
2     continue
      small = .005
      if(ip1 == 93423264) goto 1
      if(trim(nomvar) .eq. "ES") small = .004
      if(trim(nomvar) .eq. "TT") small = .01
      if(trim(nomvar) .eq. "UU") small = .004
      if(trim(nomvar) .eq. "VV") small = .004
      if(trim(nomvar) .eq. "GZ") small = .002
      if(trim(nomvar) .eq. "QC") small = 0.000001
      if(associated(zi)) then
        deallocate(zi)
        deallocate(izi)
        deallocate(zo)
        deallocate(izo)
        deallocate(izo16)
        deallocate(zt)
        deallocate(pe)
      endif
      allocate(zi(ni,nj), zo(ni,nj), zt(ni,nj))
      allocate(izi(ni,nj), izo(ni,nj), izo16(ni,nj), pe(ni,nj))
!       allocate(zt(4,ni/2,nj/2))
      call fstluk(zi,key,ni,nj,nk)   ! read record
      call fstecr(zi                         ,zi,-nbits,iunout,date,deet,npas,ni  ,nj  ,1,ip1,ip2,ip3,typvar,nomvar,"OLD",'A',0,0,0,0,datyp,.false.)
      pe(1,1) = 0.0  ! by definition there is no prediction error for lower left corner
      pe(2:ni,1) = zi(2:ni,1) - zi(1:ni-1,1)   ! first row, error is difference with previous row point
      do j = 2, nj
        pe(1,j) = zi(1,j) - zi(1,j-1)          ! first column, error is difference with previous column point
        pe(2:ni,j) = zi(2:ni,j) - zi(1:ni-1,j) + zi(1:ni-1,j-1) - zi(2:ni,j-1)  ! 2D Lorenzo prediction error
      enddo
      call AnalyzeField(loc(zi), ni, ni, nj, small, stats, STR("zi"))                       ! original field
      call AnalyzeField(loc(pe), ni, ni, nj, small, stats, STR("LO"))
!       zt(1,:,:) = zi(1:ni:2,1:nj:2)
!       zt(2,:,:) = zi(2:ni:2,1:nj:2)
!       zt(3,:,:) = zi(1:ni:2,2:nj:2)
!       zt(4,:,:) = zi(2:ni:2,2:nj:2)
      buf = ZfpCompress3D(loc(zi), ni, nj, 1, toler, 0_8, streamsize)  ! compress
      dum = ZfpCompress3D(loc(zo), ni, nj, 1, 0.0_8, buf, streamsize)  ! decompress
      call CompareFields(loc(zi), loc(zo) , ni*nj, small, nomvar//achar(0))          ! evaluate

      izo(1:ni,1:nj) = nint(zi(1:ni,1:nj) / small)
      izo16(1,1) = 0  ! by definition there is no prediction error for lower left corner
      izo16(2:ni,1) = izo(2:ni,1) - izo(1:ni-1,1)   ! first row, error is difference with previous row point
      do j = 2, nj
        izo16(1,j) = izo(1,j) - izo(1,j-1)          ! first column, error is difference with previous column point
        izo16(2:ni,j) = izo(2:ni,j) - izo(1:ni-1,j) + izo(1:ni-1,j-1) - izo(2:ni,j-1)  ! 2D Lorenzo prediction error
      enddo
      izo16(1:ni,1:nj) = izo16(1:ni,1:nj) - minval(izo16(1:ni,1:nj))
!       izo16(1:ni,1:nj) = izo(1:ni,1:nj) - minval(izo(1:ni,1:nj))
      srcSize = ni * nj
      dstCapacity = ni * nj * 4
      nbytes = 0
      nbytes16 = 0
      if(maxval(izo16) > 285) then
        nbytes = nbytes + FSE_compress(C_LOC(zt(1,1)), dstCapacity, C_LOC(izo16(1,1)), srcSize*2)
        nbytes16 = (4.0 * ni * nj)
      else
        nbytes = FSE_compress(C_LOC(zt(1,1)), dstCapacity, C_LOC(izo16(1,1)), srcSize*2)
        nbytesf = FSE_compress16(C_LOC(zt(1,1)), dstCapacity, C_LOC(izo16(1,1)), srcSize, 0, 12)
        do j = 1, nj-31 , nj/32
          nbytes16 = nbytes16 + FSE_compress16(C_LOC(zt(1,1)), dstCapacity, C_LOC(izo16(1,j)), srcSize/32, 0, 12)
        enddo
        print *,'nbytesf, nbytes16 =',nbytesf, nbytes16
        if(nbytes16 < 2) print *,'FSE_compress16 ERROR',nbytes16,maxval(izo16),minval(izo16), nbytesf
      endif

      zt(1:ni,1:nj) = zi(1:ni,1:nj)
      nbits = 16
      total_length = compress_cycle(zt, ni, nj, nbits)
      ratio = (4.0 * ni * nj) / total_length
      print *,'overall compression ratio zt(16), fse, fse16 =',ratio,(4.0 * ni * nj) / nbytes,(4.0 * ni * nj) / nbytes16, minval(izo), maxval(izo), minval(izo16), maxval(izo16)
      call AnalyzeField(loc(zt), ni, ni, nj, small, stats, STR("zt"))                       ! transformed field
      call CompareFields(loc(zi), loc(zt) , ni*nj, small, nomvar//achar(0))          ! evaluate
!       if(ip1 == 500) goto 1
!       call AnalyzeFields(loc(zi), loc(zo) , ni, ni, nj, small, stats)
      zo = zi
!       nj = (nj/8) * 8                ! adjust nj to multiple of 8
      nie2 = (ni+1)/2 
      nio2 = ni/2 
      nje2 = (nj+1)/2
      njo2 = nj/2
      nie4 = (nie2+1)/2
      nio4 = nie2/2
      nje4 = (nje2+1)/2
      njo4 = nje2/2
      nie8 = (nie4+1)/2
      nio8 = nie4/2
      nje8 = (nje4+1)/2
      njo8 = nje4/2
      ni0  = ni
      nj0  = nj
      
      call F_CDF97_2D_split_inplace_n(zo, ni, ni, nj, 3)
!       call AnalyzeField(loc(zo), ni, ni, nj, small, stats, STR("zo"))                       ! transformed field

      call AnalyzeField(loc(zo(1     ,1     )) , nie8, ni, nje8, small, stats, STR("LLLLLL"))   ! LLLLLL quadrant
      call AnalyzeField(loc(zo(1+nie8,1     )) , nio8, ni, nje8, small, stats, STR("LLLLHL"))   ! LLLLHL quadrant
      call AnalyzeField(loc(zo(1     ,1+nje8)) , nie4, ni, njo8, small, stats, STR("LLLL.H"))   ! LLLL.H quadrants
!       call AnalyzeField(loc(zo(1     ,1+nje8)) , nie8, ni, njo8, small, stats, STR("LLLLLH"))   ! LLLLLH quadrant
!       call AnalyzeField(loc(zo(1+nie8,1+nje8)) , nio8, ni, njo8, small, stats, STR("LLLLHH"))   ! LLLLHH quadrant
      
      print *,""
      call AnalyzeField(loc(zo(1+nie4,1     )) , nio4, ni, nje4, small, stats, STR("LLHL"))   ! LLHL quadrant
      call AnalyzeField(loc(zo(1     ,1+nje4)) , nie2, ni, njo4, small, stats, STR("LLLH"))   ! LL.H quadrants
!       call AnalyzeField(loc(zo(1     ,1+nje4)) , nie4, ni, njo4, small, stats, STR("LLLH"))   ! LLLH quadrant
!       call AnalyzeField(loc(zo(1+nie4,1+nje4)) , nio4, ni, njo4, small, stats, STR("LLHH"))   ! LLHH quadrant
      print *,""
      call AnalyzeField(loc(zo(1+nie2,1     )) , nio2, ni, nje2, small, stats, STR("HL"))   ! HL quadrant
      call AnalyzeField(loc(zo(1     ,1+nje2)) , ni  , ni, njo2, small, stats, STR("LH"))   ! .H quadrants
!       call AnalyzeField(loc(zo(1     ,1+nje2)) , nie2, ni, njo2, small, stats, STR("LH"))   ! LH quadrant
!       call AnalyzeField(loc(zo(1+nie2,1+nje2)) , nio2, ni, njo2, small, stats, STR("HH"))   ! HH quadrant
      if(trim(nomvar) .eq. "UU") nomvar = "VEUU"
      if(trim(nomvar) .eq. "VV") nomvar = "VEVV"
      datyp = 134
      nbits = 16
!       izo(1     :ni  ,1+nje2:nj) =     nint(zo(1     :ni  ,1+nje2:nj) / small)       ! .H quadrants
!       zo( 1     :ni  ,1+nje2:nj) = small * izo(1     :ni  ,1+nje2:nj)
!       izo(1+nie2:ni  ,1     :nje2) =     nint(zo(1+nie2:ni  ,1     :nje2) / small)   ! HL quadrant
!       zo( 1+nie2:ni  ,1     :nje2) = small * izo(1+nie2:ni  ,1     :nje2)
! 
!       small4 = small *.5
!       izo(1     :nie2,1+nje4:nje2) =      nint(zo(1     :nie2,1+nje4:nje2) / small4) ! LL.H quadrants
!       zo( 1     :nie2,1+nje4:nje2) = small4 * izo(1     :nie2,1+nje4:nje2)
!       izo(1+nie4:nie2,1     :nje4) =     nint(zo(1+nie4:nie2,1     :nje4) / small4)  ! LLHL quadrant
!       zo( 1+nie4:nie2,1     :nje4) = small4 *izo(1+nie4:nie2,1     :nje4)
! 
!       small8 = small4 *.5
!       izo(1     :nie4,1+nje8:nje4) =      nint(zo(1     :nie4,1+nje8:nje4) / small8) ! LLLL.H quadrants
!       zo( 1     :nie4,1+nje8:nje4) = small8 * izo(1     :nie4,1+nje8:nje4)
!       izo(1+nie8:nie4,1     :nje8) =      nint(zo(1+nie8:nie4,1     :nje8) / small8) ! LLLLHL quadrant
!       zo( 1+nie8:nie4,1     :nje8) = small8 * izo(1+nie8:nie4,1     :nje8)

      call fstecr(zo(1:nie8     ,1:nje8     ),zi,-nbits,iunout,date,deet,npas,nie8,nje8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLLL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1+nie8:nie4,1:nje8     ),zi,-nbits,iunout,date,deet,npas,nio8,nje8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLHL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1:nie4     ,1+nje8:nje4),zi,-nbits,iunout,date,deet,npas,nie4,njo8,1,ip1,ip2,ip3,typvar,nomvar,"LLLL.H",'A',0,0,0,0,datyp,.false.)
      total_length = 4 * nie8 * nje8
      total_length = total_length + compress_cycle(zo(1+nie8:nie4,1:nje8     ), nio8,nje8, nbits)
      total_length = total_length + compress_cycle(zo(1:nie4     ,1+nje8:nje4), nie4,njo8, nbits)
!       call fstecr(zo(1:nie8     ,1+nje8:nje4),zi,-nbits,iunout,date,deet,npas,nie8,njo8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLLH",'A',0,0,0,0,datyp,.false.)
!       call fstecr(zo(1+nie8:nie4,1+nje8:nje4),zi,-nbits,iunout,date,deet,npas,nio8,njo8,1,ip1,ip2,ip3,typvar,nomvar,"LLLLHH",'A',0,0,0,0,datyp,.false.)

!       call fstecr(zo(1:nie4     ,1:nje4     ),zi,-nbits,iunout,date,deet,npas,nie4,nje4,1,ip1,ip2,ip3,typvar,nomvar,"LLLL",'A',0,0,0,0,datyp,.false.)
      nbits = nbits - 3
      call fstecr(zo(1+nie4:nie2,1:nje4     ),zi,-nbits,iunout,date,deet,npas,nio4,nje4,1,ip1,ip2,ip3,typvar,nomvar,"LLHL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1:nie2     ,1+nje4:nje2),zi,-nbits,iunout,date,deet,npas,nie2,njo4,1,ip1,ip2,ip3,typvar,nomvar,"LL.H",'A',0,0,0,0,datyp,.false.)
      total_length = total_length + compress_cycle(zo(1+nie4:nie2,1:nje4     ), nio4, nje4, nbits)
      total_length = total_length + compress_cycle(zo(1:nie2     ,1+nje4:nje2), nie2,njo4, nbits)
!       call fstecr(zo(1:nie4     ,1+nje4:nje2),zi,-nbits,iunout,date,deet,npas,nie4,njo4,1,ip1,ip2,ip3,typvar,nomvar,"LLLH",'A',0,0,0,0,datyp,.false.)
!       call fstecr(zo(1+nie4:nie2,1+nje4:nje2),zi,-nbits,iunout,date,deet,npas,nio4,njo4,1,ip1,ip2,ip3,typvar,nomvar,"LLHH",'A',0,0,0,0,datyp,.false.)

!       call fstecr(zo(1:nie2     ,1:nje2     ),zi,-nbits,iunout,date,deet,npas,nie2,nje2,1,ip1,ip2,ip3,typvar,nomvar,"LL",'A',0,0,0,0,datyp,.false.)
      nbits = nbits - 3
      call fstecr(zo(1+nie2:ni  ,1:nje2     ),zi,-nbits,iunout,date,deet,npas,nio2,nje2,1,ip1,ip2,ip3,typvar,nomvar,"HL",'A',0,0,0,0,datyp,.false.)
      call fstecr(zo(1:ni       ,1+nje2:nj  ),zi,-nbits,iunout,date,deet,npas,ni  ,njo2,1,ip1,ip2,ip3,typvar,nomvar,".H",'A',0,0,0,0,datyp,.false.)
      total_length = total_length + compress_cycle(zo(1+nie2:ni  ,1:nje2     ), nio2,nje2, nbits)
      total_length = total_length + compress_cycle(zo(1:ni       ,1+nje2:nj  ), ni  ,njo2, nbits)
      ratio = (4.0 * ni * nj) / total_length
      print *,'overall compression ratio zo(16/14/12) =',ratio
!       call fstecr(zo(1:nie2     ,1+nje2:nj  ),zi,-nbits,iunout,date,deet,npas,nie2,njo2,1,ip1,ip2,ip3,typvar,nomvar,"LH",'A',0,0,0,0,datyp,.false.)
!       call fstecr(zo(1+nie2:ni  ,1+nje2:nj  ),zi,-nbits,iunout,date,deet,npas,nio2,njo2,1,ip1,ip2,ip3,typvar,nomvar,"HH",'A',0,0,0,0,datyp,.false.)

      call fstecr(pe                         ,zi,-nbits,iunout,date,deet,npas,ni  ,nj  ,1,ip1,ip2,ip3,typvar,nomvar,"LORENZO",'A',0,0,0,0,datyp,.false.)
      call I_CDF97_2D_split_inplace_n(zo, ni, ni, nj, 3)
      call fstecr(zo                         ,zi,-nbits,iunout,date,deet,npas,ni  ,nj  ,1,ip1,ip2,ip3,typvar,nomvar,"NEW",'A',0,0,0,0,datyp,.false.)
      call CompareFields(loc(zi), loc(zo) , ni*nj, small, nomvar//achar(0))          ! evaluate
      call AnalyzeField(loc(zo), ni, ni, nj, small, stats, STR("zo"))                       ! transformed field
      write(6,*)'------------------------------------------------------------'
1     key = fstsui(iun,ni,nj,nk)
    enddo
  enddo
  goto 888

  stop
777 continue
  stop
888 continue
!   print 9,'ni, nie2, nio2, nie4, nio4, nie8, nio8 =',ni0, nie2, nio2, nie4, nio4, nie8, nio8
!   print 9,'nj, nje2, njo2, nje4, njo4, nje8, njo8 =',nj0, nje2, njo2, nje4, njo4, nje8, njo8
  call fstfrm(iun)
  call fstfrm(iunout)
  stop
9 format(A,7I4)
end program
function compress_cycle(z, ni, nj, nbits) result(l)
  use ISO_C_BINDING
  implicit none
  real(C_FLOAT), intent(INOUT), dimension(ni,nj), target :: z
  integer(C_INT), intent(IN) :: ni, nj
  integer(C_INT), intent(INOUT) :: nbits
  integer(C_INT) :: l

  integer(C_INT), dimension(ni*nj*2), target :: stream
  integer(C_INT), dimension(3), target :: header
  integer(C_INT) :: n
  real :: ratio
  interface
    ! INT_32 c_float_packer(float *source, INT_32 nbits, INT_32 *header, INT_32 *stream, INT_32 npts)
    ! INT_32 c_float_unpacker(float *dest, INT_32 *header, INT_32 *stream, INT_32 npts, INT_32 *nbits)
    ! int armn_compress(unsigned char *fld, int ni, int nj, int nk, int nbits, int op_code)
    subroutine float_packer(src, nbits, header, stream, npts) bind(C,name='c_float_packer')
      import :: C_PTR, C_INT
      type(C_PTR), intent(IN), value :: src, stream, header
      integer(C_INT), intent(IN), value :: nbits, npts
    end subroutine float_packer
    subroutine float_unpacker(dst, header, stream, npts, nbits) bind(C,name='c_float_unpacker')
      import :: C_PTR, C_INT
      type(C_PTR), intent(IN), value :: dst, stream, header
      integer(C_INT), intent(IN), value :: npts
      integer(C_INT), intent(IN) :: nbits
    end subroutine float_unpacker
    function armn_compress(fld, ni, nj, nk, nbits, opcode) result(l) bind(C,name='armn_compress')
      import :: C_PTR, C_INT
      type(C_PTR), intent(IN), value :: fld
      integer(C_INT), intent(IN), value :: ni, nj, nk, nbits, opcode
      integer(C_INT) :: l
    end function armn_compress
  end interface
  call float_packer( C_LOC(z(1,1)), nbits, C_LOC(header(1)), C_LOC(stream(1)), ni*nj)
  l = armn_compress(C_LOC(stream(1)), ni, nj, 1, nbits, 1)
  if(l < 0) then
    call float_unpacker(C_LOC(z(1,1)), C_LOC(header(1)), C_LOC(stream(1)), ni*nj, nbits)
    l = ni*ni*nbits/8
    return
  endif
  ratio = ni * nj * 4.0
  ratio = ratio / l
!   print 1, 'ni, nj, nbits, l, ratio = ',ni, nj, nbits, l, ratio
1 format(A,4I8,f6.2)
  n = armn_compress(C_LOC(stream(1)), ni, nj, 1, nbits, 2)
! if(nj > 400) then
!   z = z + .001
!   return
! endif
  call float_unpacker(C_LOC(z(1,1)), C_LOC(header(1)), C_LOC(stream(1)), ni*nj, nbits)
end function compress_cycle

