  function interp_linz_cubyx(src, z, x, y, ini, ininj) result(val)
  implicit none
  real*4, dimension(*), intent(IN) :: src
  real*4, intent(IN) :: x, y, z
  real*4 :: val
  integer, intent(IN) :: ini, ininj

  real*4, parameter :: cp133 = 1.0 / 6.0
  real*4, parameter :: cm133 = -cp133
  real*4, parameter :: cp5   = .5
  real*4, parameter :: cm5   = -cp5
  real*4, parameter :: one   = 1.0
  real*4, parameter :: two   = 2.0

  integer :: i, ni2, ninj, ni
  real*4, dimension(0:3) :: va4, vb4, vc4, vd4, dst
  real*4 :: x0, x1, x2, x3, y0, y1, y2, y3, z0, z1

  ninj = ininj
  ni = ini
  ni2 = ni + ni

  z0 = (one - z)
  z1 = z

  y0 = cm133*y*(y-one)*(y-two)
  y1 = cp5*(y+one)*(y-one)*(y-two)
  y2 = cm5*y*(y+one)*(y-two)
  y3 = cp133*y*(y+one)*(y-one)

  x0 = cm133*x*(x-one)*(x-two)
  x1 = cp5*(x+one)*(x-one)*(x-two)
  x2 = cm5*x*(x+one)*(x-two)
  x3 = cp133*x*(x+one)*(x-one)

  do i = 0,3
    va4(i) = src(i-ni )*z0 + src(i-ni +ninj)*z1
    vb4(i) = src(i    )*z0 + src(i    +ninj)*z1
    vc4(i) = src(i+ni )*z0 + src(i+ni +ninj)*z1
    vd4(i) = src(i+ni2)*z0 + src(i+ni2+ninj)*z1
    dst(i) = va4(i)*y0 + vb4(i)*y1 + vc4(i)*y2 + vd4(i)*y3
  enddo

  val = dst(0)*x0 + dst(1)*x1 + dst(2)*x2 + dst(3)*x3

  return
  end


