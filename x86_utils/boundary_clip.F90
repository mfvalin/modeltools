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
!-------------------------------------------------------------------------
!
! un peu de geometrie basee sur les relations entre "triangles semblables"
!
!                     (x1,y1)*---------*(x0,y1)                  *(X1,Y1)
!                             \        |                        /|
!                              \       |                       / |
!                               \      |                      /  |
!             ymax  +------------*-----*--------------------+/   |
!                   |      (x2,y2)\    |(x0,y2)      (X2,Y2)*----*(X1,Y2)
!                   |              \   |                   /|    |
!                   |               \  |                  / |    |
!                   |                \ |                 /  |    |
!                   |                 \|                /   |    |
!                   |                  *(x0,y0)        /    |    |
!                   |                                 /     |    |
!                   |                         (X0,Y0)*------*----*(X1,Y0)
!                   |                                       |(X2,Y0)
!                   |                                       |
!                   |                                       |
!                   |                                       |
!                   |                                       |
!             ymin  +---------------------------------------+
!                   xmin                                 xmax
!
!  pour le segment de droite (x0,y0) ---> (x1,y1)
!
!  point d'intersection (x2,y2): (intersection avec y = constante)
!
!                 y2 = constante y (ymin ou ymax)
!
!                 (x1-x0)   (y1-y0)
!                 ------- = -------
!                 (x2-x0)   (y2-y0)  (y2 = ymax ou ymin selon le cas)
!
!         donc    x2 = x0 + (x1-x0)/(y1-y0) * (y2-y0)
!
!  pour le segment de droite (X0,Y0) ---> (X1,Y1)
!
!  point d'intersection (X2,Y2): (intersection avec x = constante)
!
!                 X2 = constante x (xmin ou xmax)
!
!                 (X1-X0)   (Y1-Y0)
!                 ------- = -------
!                 (X2-X0)   (Y2-Y0)  (X2 = xmax ou xmin selon le cas)
!
!         donc    Y2 = Y0 + (Y1-Y0)/(X1-X0) * (X2-X0)
!
! recette:
!
! (x2,y2) = (x1,y1) , on suppose que (x2,y2) est dans la fenetre
!
! si le point (x1,y1) est au "nord" ou au "sud" de la fenetre
!   on calcule l'intersection (x2,y2) (y2 = ymin ou ymax)
!
! si x2 est a l'exterieur de l'intervalle xmin --- xmax
!   on passe au cas "est" ou "ouest"
!   on calcule l'intersection (X2,Y2)  (X2 = xmin ou xmax)
!
! et dans ce cas, on est sur que Y2 sera dans l'intervalle (ymin --- ymax)
!
function boundary_clip_point(p0,p1,p2,l) result(clipped)
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
end function boundary_clip_point

! for p0, p1, p2 index 1 is x, index 2 is y
! l : [xmin, xmax, ymin, ymax]
function boundary_clip_coord_1(p0,p1,p2,l) result(clipped)
  use ISO_C_BINDING
  implicit none
  real, dimension(2), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  real, dimension(2), intent(IN)  :: p1   ! point to check
  real, dimension(2), intent(OUT) :: p2   ! 
  real, dimension(4), intent(IN) :: l     ! window used to clip p0 -> p1 segment
  logical :: clipped                      ! .true. if clipping occurred

  logical :: north, south, east, west

  north = p1(2) > l(4)     ! point 1 is above upper y limit
  south = p1(2) < l(3)     ! point 1 is below lower y limit
  p2 = p1
  if(north) then            ! clip using upper y limit
    p2(2) = l(4)           ! clipped y coordinate is upper y limit
    p2(1) = p0(1) + ((p1(1)-p0(1))/(p1(2)-p0(2))) * (p2(2)-p0(2))
  endif
  if(south) then            ! clip using lower y limit
    p2(2) = l(3)           ! clipped y coordinate is lower y limit
    p2(1) = p0(1) + ((p1(1)-p0(1))/(p1(2)-p0(2))) * (p2(2)-p0(2))
  endif
  east = p2(1) > l(2)     ! point 1 (possibly after clipping) is right of upper x limit
  west = p2(1) < l(1)     ! point 1 (possibly after clipping) is left of lower x limit
  if(east .or. west) then
    if(east) then            ! clip using upper x limit
     p2(1) = l(2)           ! clipped x coordinate is upper x limit
     p2(2) = p0(2) + ((p1(2)-p0(2))/(p1(1)-p0(1))) * (p2(1)-p0(1))
    endif
    if(west) then            ! clip using lower x limit
     p2(1) = l(1)           ! clipped x coordinate is lower x limit
     p2(2) = p0(2) + ((p1(2)-p0(2))/(p1(1)-p0(1))) * (p2(1)-p0(1))
    endif
  endif
  clipped = north .or. south .or. east .or. west  ! has clipping occurred
end function boundary_clip_coord_1

! for p0, p1, p2 index 1 is x, index 2 is y
! l : [xmin, xmax, ymin, ymax]
subroutine boundary_clip_coord_2n(p0,p1,p2,l,clipped,n)
  use ISO_C_BINDING
  implicit none
  real, dimension(2,n), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  real, dimension(2,n), intent(IN)  :: p1   ! point to check
  real, dimension(2,n), intent(OUT) :: p2   ! 
  real, dimension(4), intent(IN) :: l     ! window used to clip p0 -> p1 segment
  integer, intent(IN) :: n
  logical, dimension(N) :: clipped                      ! .true. if clipping occurred

  logical :: north, south, east, west
  integer :: i

  do i = 1 , n
    north = p1(2,i) > l(4)     ! point 1 is above upper y limit
    south = p1(2,i) < l(3)     ! point 1 is below lower y limit
    p2 = p1
    if(north) then            ! clip using upper y limit
      p2(2,i) = l(4)           ! clipped y coordinate is upper y limit
      p2(1,i) = p0(1,i) + ((p1(1,i)-p0(1,i))/(p1(2,i)-p0(2,i))) * (p2(2,i)-p0(2,i))
    endif
    if(south) then            ! clip using lower y limit
      p2(2,i) = l(3)           ! clipped y coordinate is lower y limit
      p2(1,i) = p0(1,i) + ((p1(1,i)-p0(1,i))/(p1(2,i)-p0(2,i))) * (p2(2,i)-p0(2,i))
    endif
    east = p2(1,i) > l(2)     ! point 1 (possibly after clipping) is right of upper x limit
    west = p2(1,i) < l(1)     ! point 1 (possibly after clipping) is left of lower x limit
    if(east .or. west) then
      if(east) then            ! clip using upper x limit
      p2(1,i) = l(2)           ! clipped x coordinate is upper x limit
      p2(2,i) = p0(2,i) + ((p1(2,i)-p0(2,i))/(p1(1,i)-p0(1,i))) * (p2(1,i)-p0(1,i))
      endif
      if(west) then            ! clip using lower x limit
      p2(1,i) = l(1)           ! clipped x coordinate is lower x limit
      p2(2,i) = p0(2,i) + ((p1(2,i)-p0(2,i))/(p1(1,i)-p0(1,i))) * (p2(1,i)-p0(1,i))
      endif
    endif
    clipped(i) = north .or. south .or. east .or. west  ! has clipping occurred
  enddo
end subroutine boundary_clip_coord_2n

! for p0, p1, p2 index 1 is x, index 2 is y
! l : [xmin, xmax, ymin, ymax]
subroutine boundary_clip_coord_n2(p0,p1,p2,l,clipped,n)
  use ISO_C_BINDING
  implicit none
  real, dimension(n,2), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  real, dimension(n,2), intent(IN)  :: p1   ! point to check
  real, dimension(n,2), intent(OUT) :: p2   ! 
  real, dimension(4), intent(IN) :: l     ! window used to clip p0 -> p1 segment
  integer, intent(IN) :: n
  logical, dimension(N) :: clipped                      ! .true. if clipping occurred

  logical :: north, south, east, west
  integer :: i

  do i = 1 , n
    north = p1(i,2) > l(4)     ! point 1 is above upper y limit
    south = p1(i,2) < l(3)     ! point 1 is below lower y limit
    p2 = p1
    if(north) then            ! clip using upper y limit
      p2(i,2) = l(4)           ! clipped y coordinate is upper y limit
      p2(i,1) = p0(i,1) + ((p1(i,1)-p0(i,1))/(p1(i,2)-p0(i,2))) * (p2(i,2)-p0(i,2))
    endif
    if(south) then            ! clip using lower y limit
      p2(i,2) = l(3)           ! clipped y coordinate is lower y limit
      p2(i,1) = p0(i,1) + ((p1(i,1)-p0(i,1))/(p1(i,2)-p0(i,2))) * (p2(i,2)-p0(i,2))
    endif
    east = p2(i,1) > l(2)     ! point 1 (possibly after clipping) is right of upper x limit
    west = p2(i,1) < l(1)     ! point 1 (possibly after clipping) is left of lower x limit
    if(east .or. west) then
      if(east) then            ! clip using upper x limit
      p2(i,1) = l(2)           ! clipped x coordinate is upper x limit
      p2(i,2) = p0(i,2) + ((p1(i,2)-p0(i,2))/(p1(i,1)-p0(i,1))) * (p2(i,1)-p0(i,1))
      endif
      if(west) then            ! clip using lower x limit
      p2(i,1) = l(1)           ! clipped x coordinate is lower x limit
      p2(i,2) = p0(i,2) + ((p1(i,2)-p0(i,2))/(p1(i,1)-p0(i,1))) * (p2(i,1)-p0(i,1))
      endif
    endif
    clipped(i) = north .or. south .or. east .or. west  ! has clipping occurred
  enddo
end subroutine boundary_clip_coord_n2

#if defined(SELF_TEST)
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
    function boundary_clip_point(p0,p1,p2,l) result(clipped)
    import :: point, limits
    type(point), intent(IN)  :: p0, p1
    type(point), intent(OUT) :: p2
    type(limits), intent(IN) :: l
    logical :: clipped
    end function boundary_clip_point
  end interface
  logical, external :: boundary_clip_coord_1
  type(point) :: p2
  real, dimension(2) :: a2
  logical :: clipped
  logical, dimension(1) :: clippen

  clipped = boundary_clip_point( point(0.0,0.0), point(1.0,1.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip_coord_1( [0.0,0.0], [1.0,1.0], a2, [-1.0,1.0,-1.0,1.0])
  print *,clipped,a2
  call boundary_clip_coord_2n( [0.0,0.0], [1.0,1.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  call boundary_clip_coord_n2( [0.0,0.0], [1.0,1.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  print *,""
  clipped = boundary_clip_point( point(0.0,0.0), point(2.0,2.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip_coord_1( [0.0,0.0], [2.0,2.0], a2, [-1.0,1.0,-1.0,1.0])
  print *,clipped,a2
  call boundary_clip_coord_2n( [0.0,0.0], [2.0,2.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  call boundary_clip_coord_n2( [0.0,0.0], [2.0,2.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  print *,""
  clipped = boundary_clip_point( point(0.0,0.0), point(3.0,2.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip_coord_1( [0.0,0.0], [3.0,2.0], a2, [-1.0,1.0,-1.0,1.0])
  print *,clipped,a2
  call boundary_clip_coord_2n( [0.0,0.0], [3.0,2.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  call boundary_clip_coord_n2( [0.0,0.0], [3.0,2.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  print *,""
  clipped = boundary_clip_point( point(0.0,0.0), point(2.0,3.0), p2, limits(-1.0,1.0,-1.0,1.0))
  print *,clipped,p2
  clipped = boundary_clip_coord_1( [0.0,0.0], [2.0,3.0], a2, [-1.0,1.0,-1.0,1.0])
  print *,clipped,a2
  call boundary_clip_coord_2n( [0.0,0.0], [2.0,3.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  call boundary_clip_coord_n2( [0.0,0.0], [2.0,3.0], a2, [-1.0,1.0,-1.0,1.0],clippen,1)
  print *,clippen,a2
  print *,""
end program
#endif
