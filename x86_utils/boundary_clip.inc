interface
function boundary_clip_coord_1(p0,p1,p2,l) result(clipped) BIND(C,name='BoundaryClipCoord1')
  real, dimension(2), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  real, dimension(2), intent(IN)  :: p1   ! point to check
  real, dimension(2), intent(OUT) :: p2   ! clipped output
  real, dimension(4), intent(IN) :: l     ! window used to clip p0 -> p1 segment
  logical :: clipped                      ! .true. if clipping occurred
end function boundary_clip_coord_1
subroutine boundary_clip_coord_2n(p0,p1,p2,l,clipped,n) BIND(C,name='BoundaryClipCoord2n')
  real, dimension(2,n), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  real, dimension(2,n), intent(IN)  :: p1   ! point to check
  real, dimension(2,n), intent(OUT) :: p2   ! intersection with clipping window
  real, dimension(4), intent(IN) :: l       ! window used to clip p0 -> p1 segment
  integer, intent(IN) :: n                  ! number of points
  logical, dimension(N) :: clipped          ! .true. if clipping occurred
end subroutine boundary_clip_coord_2n
subroutine boundary_clip_coord_n2(p0,p1,p2,l,clipped,n) BIND(C,name='BoundaryClipCoordn2')
  real, dimension(n,2), intent(IN)  :: p0   ! ASSUMED to be inside the window defined by l
  real, dimension(n,2), intent(IN)  :: p1   ! point to check
  real, dimension(n,2), intent(OUT) :: p2   ! intersection with clipping window
  real, dimension(4), intent(IN) :: l       ! window used to clip p0 -> p1 segment
  integer, intent(IN) :: n                  ! number of points
  logical, dimension(N) :: clipped          ! .true. if clipping occurred
end subroutine boundary_clip_coord_n2
end interface
