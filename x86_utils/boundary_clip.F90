!
! RMNLIB 
! Copyright (C) 2019 Environment Canada
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation, 
! version 2 of the License.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public
! License along with this library; if not, write to the
! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! Boston, MA 02111-1307, USA.
!
function boundary_clip(p0,p1,p2,l) result(clipped)
  use ISO_C_BINDING
  implicit none
  type :: point
    SEQUENCE
    real :: x  ! x coordinate
    real :: y  ! y coordinate
  end type
  type :: limits
    SEQUENCE
    real :: xmin  ! left boundary
    real :: xmax  ! right boundary
    real :: ymin  ! top boundary
    real :: ymax  ! bottom boundary
  end type
  type(point), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  type(point), intent(IN)  :: p1   ! point to check
  type(point), intent(OUT) :: p2   ! 
  type(limits), intent(IN) :: l    ! window used to clip p0 -> p1 segment
  logical :: clipped               ! .true. if clipping occurred

  logical :: north, south, east, west

  north = p1%y > l%ymax     ! point 1 is above upper y limit
  south = p1%y < l%ymin     ! point 1 is below lower y limit
  p2 = p1
  if(north) then            ! clip using upper y limit
    p2%y = l%ymax           ! clipped y coordinate is upper y limit
    p2%x = p0%x + ((p1%x-p0%x)/(p1%y-p0%y)) * (p2%y-p0%y)
  endif
  if(south) then            ! clip using lower y limit
    p2%y = l%ymin           ! clipped y coordinate is lower y limit
    p2%x = p0%x + ((p1%x-p0%x)/(p1%y-p0%y)) * (p2%y-p0%y)
  endif
  east = p2%x > l%xmax     ! point 1 (possibly after clipping) is right of upper x limit
  west = p2%x < l%xmin     ! point 1 (possibly after clipping) is left of lower x limit
  if(east .or. west) then
    if(east) then            ! clip using upper x limit
     p2%x = l%xmax           ! clipped x coordinate is upper x limit
     p2%y = p0%y + ((p1%y-p0%y)/(p1%x-p0%x)) * (p2%x-p0%x)
    endif
    if(west) then            ! clip using lower x limit
     p2%x = l%xmin           ! clipped x coordinate is lower x limit
     p2%y = p0%y + ((p1%y-p0%y)/(p1%x-p0%x)) * (p2%x-p0%x)
    endif
  endif
  clipped = north .or. south .or. east .or. west  ! has clipping occurred
end function boundary_clip
program test
  implicit none
  type :: point
    SEQUENCE
    real :: x  ! x coordinate
    real :: y  ! y coordinate
  end type
  type :: limits
    SEQUENCE
    real :: xmin  ! left boundary
    real :: xmax  ! right boundary
    real :: ymin  ! top boundary
    real :: ymax  ! bottom boundary
  end type
  interface
    function boundary_clip(p0,p1,p2,l) result(clipped)
    import :: point, limits
    type(point), intent(IN)  :: p0, p1
    type(point), intent(OUT) :: p2
    type(limits), intent(IN) :: l
    logical :: clipped
    end function boundary_clip
  end interface
  type(point) :: p2
  logical :: clipped
  clipped = boundary_clip( point(0.0,0.0), point(1.0,1.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip( point(0.0,0.0), point(2.0,2.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip( point(0.0,0.0), point(3.0,2.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip( point(0.0,0.0), point(2.0,3.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
end program
