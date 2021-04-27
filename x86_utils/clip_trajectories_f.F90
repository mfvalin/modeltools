
module clip_trajectories_data
  save
  real(kind=4) :: c_mini = 0.0
  real(kind=4) :: c_maxi = 0.0
  real(kind=4) :: c_minj = 0.0
  real(kind=4) :: c_maxj = 0.0
  real(kind=4) :: g_mini = 0.0
  real(kind=4) :: g_maxi = 0.0
  real(kind=4) :: g_minj = 0.0
  real(kind=4) :: g_maxj = 0.0
end module

subroutine Set_Clipping_Limits(cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj)
  use clip_trajectories_data
  implicit none
  real(kind=4), intent(IN) :: cmini, cmaxi, cminj, cmaxj , gmini, gmaxi, gminj, gmaxj

  c_mini = cmini
  c_maxi = cmaxi
  c_minj = cminj
  c_maxj = cmaxj
  g_mini = gmini
  g_maxi = gmaxi
  g_minj = gminj
  g_maxj = gmaxj
end subroutine

#define VL 16
subroutine clip_trajectories_f(alpha, beta, l1, l2, ni, indx1, indx2)
  use clip_trajectories_data
  implicit none
  integer, intent(IN)                     :: ni
  real(kind=4), dimension(ni), intent(IN) :: alpha, beta
  integer, dimension(*), intent(OUT)      :: l1, l2
  integer, intent(INOUT)                  :: indx1       ! number of points in list l1 (output)
                                                         ! last valid entry in l1 or l2 (input) (0 = none)
  integer, intent(INOUT)                  :: indx2       ! number of points in list l2 (output)

  logical, dimension(0:VL-1) :: clip1, clip2
  integer :: i, k ,vl7, ix1, ix2

  vl7 = mod(ni,VL)
  ix1 = indx1 + 1
  ix2 = indx2 + 1

  do i = 1, ni-VL, VL     !  VL long slices
    do k = 0,VL-1
      clip2(k) = (alpha(i+k) < g_mini) .or. (alpha(i+k) > g_maxi) .or. (beta(i+k) < g_minj) .or. (beta(i+k) > g_maxj)
      clip1(k) = (alpha(i+k) < c_mini) .or. (alpha(i+k) > c_maxi) .or. (beta(i+k) < c_minj) .or. (beta(i+k) > c_maxj)
      clip1(k) = clip1(k) .and. (.not. clip2(k))
    enddo
    do k = 0,VL-1
      l1(ix1) = i+k
      if(clip1(k)) ix1 = ix1 + 1
      l2(ix2) = i+k
      if(clip2(k)) ix2 = ix2 + 1
    enddo
  enddo
  i = ni - vl7 + 1
  do k = 0, vl7-1       ! letfovers after VL slices
    clip2(k) = (alpha(i+k) < g_mini) .or. (alpha(i+k) > g_maxi) .or. (beta(i+k) < g_minj) .or. (beta(i+k) > g_maxj)
    clip1(k) = (alpha(i+k) < c_mini) .or. (alpha(i+k) > c_maxi) .or. (beta(i+k) < c_minj) .or. (beta(i+k) > c_maxj)
    clip1(k) = clip1(k) .and. (.not. clip2(k))
  enddo
  do k = 0,vl7-1
    l1(ix1) = i+k
    if(clip1(k)) ix1 = ix1 + 1
    l2(ix2) = i+k
    if(clip2(k)) ix2 = ix2 + 1
  enddo

  indx1 = ix1-1   ! set to number of points
  indx2 = ix2-1   ! set to number of points
end subroutine

#if defined(SELF_TEST)
#if ! defined(NPTS)
#define NPTS 128007
#endif

function Wall_Time() result(t)
  use ISO_C_BINDING
  implicit none
  real(kind=8) :: t
  interface
    function gettimeofday(tv, tz) result(status) bind(C,name='gettimeofday')
      import :: C_PTR, C_INT
      type(C_PTR), intent(IN), value :: tv, tz
      integer(C_INT) :: status
    end function gettimeofday
  end interface
  type, bind(C) :: timeval
    integer(C_LONG_LONG) :: tv_sec
    integer(C_LONG_LONG) :: tv_usec
  end type
  integer :: status
  type(timeval), target :: tv
  status = gettimeofday(C_LOC(tv), C_NULL_PTR)
  t = tv % tv_sec
  t = t + tv % tv_usec * 1E-6 ;
end function Wall_Time

program test_clip
  implicit none
  real(kind=4), dimension(NPTS) :: alpha, beta
  integer, dimension(NPTS) :: l1, l2
  integer :: i, ii, indx1, indx2
  real(kind=8) :: t0, t1
  real(kind=8), external :: Wall_Time
  integer :: ierr

  call Set_Clipping_Limits(0.06, 0.94, 0.06, 0.94, 0.03, .97, 0.03, .97)
  ii = 1
  do i = 0,NPTS-1
     alpha(ii) = (1.0 * i) / (NPTS - 1)
     beta(ii)  = (1.0 * i) / (NPTS - 1)
     ii = ii + 41
     if(ii > NPTS) ii = ii - NPTS
  enddo
  indx1 = 0 ;
  indx2 = 0 ;
  t0 = Wall_Time()
  call clip_trajectories_f(alpha, beta, l1, l2, NPTS, indx1, indx2)
  t1 = Wall_Time()
  print *, indx1, indx2,(t1 - t0)/NPTS * 1.0E9
  if(indx1 < 10) write(6,'(A,10I8)')'indx1 =',l1(1:indx1)
  if(indx2 < 10) write(6,'(A,10I8)')'indx2 =',l2(1:indx2)
end program
#endif
