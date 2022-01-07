! Copyright (C) 2021  Environnement et Changement climatique Canada
!
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
!
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! Author:
!     M. Valin,   Recherche en Prevision Numerique, 2021

#if ! defined(VL)
#define VL 8
#endif

module clip_trajectories_d_mod
  use ISO_C_BINDING
  implicit none

  private
  interface
    subroutine SetClippingLimits_d_avx512(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj) bind(C,name='SetClippingLimitsD_avx')
      import :: C_DOUBLE
      implicit none
      real(C_DOUBLE), intent(IN), value :: cmini, cmaxi, cminj, cmaxj
      real(C_DOUBLE), intent(IN), value :: gmini, gmaxi, gminj, gmaxj
    end subroutine
    subroutine ClipTrajectories_d_avx512(alpha, beta, l1, l2, ni, nj, mini,maxi,minj,maxj, nk, indx1, indx2) bind(C,name='ClipTrajectoriesD_avx512')
      import :: C_DOUBLE, C_INT
      implicit none
      integer(C_INT), intent(IN), value :: ni, nj, mini,maxi,minj,maxj, nk    ! dimensions of alpha and beta
      real(C_DOUBLE), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: alpha   ! targets along x (i)
      real(C_DOUBLE), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: beta    ! targets along y (j)
      integer(C_INT), dimension(*), intent(OUT) :: l1            ! clipping list for tile (cell)
      integer(C_INT), dimension(*), intent(OUT) :: l2            ! clipping list for grid
      integer(C_INT), intent(INOUT) :: indx1                     ! number of points clipped in l1 list
      integer(C_INT), intent(INOUT) :: indx2                     ! number of points clipped in l2 list
    end subroutine
    function ClipTrajectories_avx512_d() result(available) bind(C,name='ClipTrajectories_avx512D')
      import :: C_INT
      implicit none
      integer(C_INT) :: available
    end function
    function ClipTrajectories_avx2_d() result(available) bind(C,name='ClipTrajectories_avx2D')
      import :: C_INT
      implicit none
      integer(C_INT) :: available
    end function
  end interface

  save
  real(C_DOUBLE) :: c_mini = 0.0
  real(C_DOUBLE) :: c_maxi = 0.0
  real(C_DOUBLE) :: c_minj = 0.0
  real(C_DOUBLE) :: c_maxj = 0.0
  real(C_DOUBLE) :: g_mini = 0.0
  real(C_DOUBLE) :: g_maxi = 0.0
  real(C_DOUBLE) :: g_minj = 0.0
  real(C_DOUBLE) :: g_maxj = 0.0
  integer(C_INT) :: avx512 = -1
  integer(C_INT) :: avx2   = -1

  public :: SetClippingLimits_d, ClipTrajectories_d, ClipTrajectories_avx512_d, ClipTrajectories_avx2_d

  contains

  subroutine SetClippingLimits_d(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj) bind(C,name='SetClippingLimitsD')
    implicit none
    real(C_DOUBLE), intent(IN), value :: cmini, cmaxi, cminj, cmaxj
    real(C_DOUBLE), intent(IN), value :: gmini, gmaxi, gminj, gmaxj

    if(avx512 == -1) then
      avx512 = ClipTrajectories_avx512_d()
      if(avx512 == 1) print *,'AVX512 version available'
    endif
    if(avx2 == -1) then
      avx2 = ClipTrajectories_avx2_d()
      if(avx2 == 1) print *,'AVX2 version available'
    endif
    if(avx512 == 1 .or. avx2 == 1) call SetClippingLimits_d_avx512(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj)
    c_mini = cmini
    c_maxi = cmaxi
    c_minj = cminj
    c_maxj = cmaxj
    g_mini = gmini
    g_maxi = gmaxi
    g_minj = gminj
    g_maxj = gmaxj
  end subroutine

  subroutine ClipTrajectories_d(alpha, beta, l1, l2, ni, nj, mini,maxi,minj,maxj, nk, indx1, indx2) bind(C,name='ClipTrajectoriesD')
    implicit none
    integer(C_INT), intent(IN), value :: ni, nj, mini,maxi,minj,maxj, nk    ! dimensions of alpha and beta
    real(C_DOUBLE), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: alpha   ! targets along x (i)
    real(C_DOUBLE), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: beta    ! targets along y (j)
    integer(C_INT), dimension(*), intent(OUT) :: l1            ! clipping list for tile (cell)
    integer(C_INT), dimension(*), intent(OUT) :: l2            ! clipping list for grid
    integer(C_INT), intent(INOUT) :: indx1                     ! number of points clipped in l1 list
    integer(C_INT), intent(INOUT) :: indx2                     ! number of points clipped in l2 list

    integer :: i, j, k, l, nim
    integer :: ix1, ix2, ijk, jk
    logical, dimension(1:ni) :: clip1, clip2

    if(avx512 == 1) then
      call ClipTrajectories_d_avx512(alpha, beta, l1, l2, ni, nj, mini,maxi,minj,maxj, nk, indx1, indx2)
      return
    endif
    ix1 = indx1 + 1  ! next store position
    ix2 = indx2 + 1
    nim = mod(ni,VL)
    if(nim == 0) nim = VL
    do k = 1, nk
      do j = 1, nj
        jk = j * 4096 + k * 4096 * 4096
        if(VL > ni) then
          do i = 1, ni
            clip2(i) = (alpha(i,j,k) < g_mini) .or. (alpha(i,j,k) > g_maxi) .or. (beta(i,j,k) < g_minj) .or. (beta(i,j,k) > g_maxj)
            clip1(i) = (alpha(i,j,k) < c_mini) .or. (alpha(i,j,k) > c_maxi) .or. (beta(i,j,k) < c_minj) .or. (beta(i,j,k) > c_maxj)
            clip1(i) = clip1(i) .and. (.not. clip2(i))
          enddo
        else
          do i = 1, VL     ! force a vector pass instead of a scalar pass for slice 1
            clip2(i) = (alpha(i,j,k) < g_mini) .or. (alpha(i,j,k) > g_maxi) .or. (beta(i,j,k) < g_minj) .or. (beta(i,j,k) > g_maxj)
            clip1(i) = (alpha(i,j,k) < c_mini) .or. (alpha(i,j,k) > c_maxi) .or. (beta(i,j,k) < c_minj) .or. (beta(i,j,k) > c_maxj)
            clip1(i) = clip1(i) .and. (.not. clip2(i))
          enddo
          do i = nim+1, ni, VL
            do l=1,VL
              clip2(i+l) = (alpha(i+l,j,k) < g_mini) .or. (alpha(i+l,j,k) > g_maxi) .or. (beta(i+l,j,k) < g_minj) .or. (beta(i+l,j,k) > g_maxj)
              clip1(i+l) = (alpha(i+l,j,k) < c_mini) .or. (alpha(i+l,j,k) > c_maxi) .or. (beta(i+l,j,k) < c_minj) .or. (beta(i+l,j,k) > c_maxj)
              clip1(i+l) = clip1(i+l) .and. (.not. clip2(i+l))
            enddo
          enddo
        endif
        do i = 1, ni
          ijk = i + jk
          l1(ix1) = ijk
          if(clip1(i)) ix1 = ix1 + 1
          l2(ix2) = ijk
          if(clip2(i)) ix2 = ix2 + 1
        enddo
      enddo
    enddo
    indx1 = ix1 - 1  ! position of last store
    indx2 = ix2 - 1
  end subroutine ClipTrajectories_d

end module

module clip_trajectories_f_mod
  use ISO_C_BINDING
  implicit none

  private
  interface
    subroutine SetClippingLimits_f_avx512(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj) bind(C,name='SetClippingLimitsF_avx512')
      import :: C_FLOAT
      implicit none
      real(C_FLOAT), intent(IN), value :: cmini, cmaxi, cminj, cmaxj
      real(C_FLOAT), intent(IN), value :: gmini, gmaxi, gminj, gmaxj
    end subroutine
    subroutine ClipTrajectories_f_avx512(alpha, beta, l1, l2, ni, nj, mini,maxi,minj,maxj, nk, indx1, indx2) bind(C,name='ClipTrajectoriesF_avx512')
      import :: C_FLOAT, C_INT
      implicit none
      integer(C_INT), intent(IN), value :: ni, nj, mini,maxi,minj,maxj, nk    ! dimensions of alpha and beta
      real(C_FLOAT), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: alpha   ! targets along x (i)
      real(C_FLOAT), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: beta    ! targets along y (j)
      integer(C_INT), dimension(*), intent(OUT) :: l1            ! clipping list for tile (cell)
      integer(C_INT), dimension(*), intent(OUT) :: l2            ! clipping list for grid
      integer(C_INT), intent(INOUT) :: indx1                     ! number of points clipped in l1 list
      integer(C_INT), intent(INOUT) :: indx2                     ! number of points clipped in l2 list
    end subroutine
    function ClipTrajectories_avx512_f() result(available) bind(C,name='ClipTrajectories_avx512F')
      import :: C_INT
      implicit none
      integer(C_INT) :: available
    end function
  end interface
  save
  real(C_FLOAT) :: c_mini = 0.0
  real(C_FLOAT) :: c_maxi = 0.0
  real(C_FLOAT) :: c_minj = 0.0
  real(C_FLOAT) :: c_maxj = 0.0
  real(C_FLOAT) :: g_mini = 0.0
  real(C_FLOAT) :: g_maxi = 0.0
  real(C_FLOAT) :: g_minj = 0.0
  real(C_FLOAT) :: g_maxj = 0.0
  integer(C_INT) :: avx512 = -1

  public :: SetClippingLimits_f, ClipTrajectories_f, ClipTrajectories_avx512_f

  contains

  subroutine SetClippingLimits_f(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj) bind(C,name='SetClippingLimitsF')
    implicit none
    real(C_FLOAT), intent(IN), value :: cmini, cmaxi, cminj, cmaxj
    real(C_FLOAT), intent(IN), value :: gmini, gmaxi, gminj, gmaxj

!     if(avx512 == -1) avx512 = ClipTrajectories_avx512_f()
!     if(avx512 == 1) call SetClippingLimits_f_avx512(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj)
    c_mini = cmini
    c_maxi = cmaxi
    c_minj = cminj
    c_maxj = cmaxj
    g_mini = gmini
    g_maxi = gmaxi
    g_minj = gminj
    g_maxj = gmaxj
  end subroutine

  subroutine ClipTrajectories_f(alpha, beta, l1, l2, ni, nj, mini,maxi,minj,maxj, nk, indx1, indx2) bind(C,name='ClipTrajectoriesF')
    implicit none
    integer(C_INT), intent(IN), value :: ni, nj, mini,maxi,minj,maxj, nk    ! dimensions of alpha and beta
    real(C_FLOAT), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: alpha   ! targets along x (i)
    real(C_FLOAT), dimension(mini:maxi,minj:maxj,nk), intent(IN) :: beta    ! targets along y (j)
    integer(C_INT), dimension(*), intent(OUT) :: l1            ! clipping list for tile (cell)
    integer(C_INT), dimension(*), intent(OUT) :: l2            ! clipping list for grid
    integer(C_INT), intent(INOUT) :: indx1                     ! number of points clipped in l1 list
    integer(C_INT), intent(INOUT) :: indx2                     ! number of points clipped in l2 list

    integer :: i, j, k, l, nim
    integer :: ix1, ix2, ijk, jk
    logical, dimension(1:ni) :: clip1, clip2

    ix1 = indx1 + 1  ! next store position
    ix2 = indx2 + 1
    nim = mod(ni,VL)
    if(nim == 0) nim = VL
    do k = 1, nk
      do j = 1, nj
        jk = j * 4096 + k * 4096 * 4096
        if(VL > ni) then
          do i = 1, ni
            clip2(i) = (alpha(i,j,k) < g_mini) .or. (alpha(i,j,k) > g_maxi) .or. (beta(i,j,k) < g_minj) .or. (beta(i,j,k) > g_maxj)
            clip1(i) = (alpha(i,j,k) < c_mini) .or. (alpha(i,j,k) > c_maxi) .or. (beta(i,j,k) < c_minj) .or. (beta(i,j,k) > c_maxj)
            clip1(i) = clip1(i) .and. (.not. clip2(i))
          enddo
        else
          do i = 1, VL     ! force a vector pass instead of a scalar pass for slice 1
            clip2(i) = (alpha(i,j,k) < g_mini) .or. (alpha(i,j,k) > g_maxi) .or. (beta(i,j,k) < g_minj) .or. (beta(i,j,k) > g_maxj)
            clip1(i) = (alpha(i,j,k) < c_mini) .or. (alpha(i,j,k) > c_maxi) .or. (beta(i,j,k) < c_minj) .or. (beta(i,j,k) > c_maxj)
            clip1(i) = clip1(i) .and. (.not. clip2(i))
          enddo
          do i = nim+1, ni, VL
            do l=1,VL
              clip2(i+l) = (alpha(i+l,j,k) < g_mini) .or. (alpha(i+l,j,k) > g_maxi) .or. (beta(i+l,j,k) < g_minj) .or. (beta(i+l,j,k) > g_maxj)
              clip1(i+l) = (alpha(i+l,j,k) < c_mini) .or. (alpha(i+l,j,k) > c_maxi) .or. (beta(i+l,j,k) < c_minj) .or. (beta(i+l,j,k) > c_maxj)
              clip1(i+l) = clip1(i+l) .and. (.not. clip2(i+l))
            enddo
          enddo
        endif
        do i = 1, ni
          ijk = i + jk
          l1(ix1) = ijk
          if(clip1(i)) ix1 = ix1 + 1
          l2(ix2) = ijk
          if(clip2(i)) ix2 = ix2 + 1
        enddo
      enddo
    enddo
    indx1 = ix1 - 1  ! position of last store
    indx2 = ix2 - 1
  end subroutine ClipTrajectories_f

end module

module clip_trajectories_mod
  use clip_trajectories_d_mod
  use clip_trajectories_f_mod
  implicit none
end module
