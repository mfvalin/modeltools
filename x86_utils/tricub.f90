  real function tricub_x86_f(src, abcd, x, y, ini, ininj)
  implicit none
  real*4, dimension(0:*), intent(IN) :: src
  real*4, dimension(0:3), intent(IN) :: abcd
  real*4, intent(IN) :: x, y
  integer, intent(IN) :: ini, ininj

  real*4, parameter :: cp133 = 0.166666666667
  real*4, parameter :: cm133 = -0.16666666667
  real*4, parameter :: cp5 = .5
  real*4, parameter :: cm5 = -.5
  real*4, parameter :: one = 1.0
  real*4, parameter :: two = 2.0

  integer :: i, ni2, ni3, ninj2, ninj3, ninj, ni
  real*4, dimension(0:3) :: va4, vb4, vc4, vd4, dst
  real*4 :: x0, x1, x2, x3, y0, y1, y2, y3

  ninj = ininj
  ni = ini
  ninj2 = ninj + ninj
  ninj3 = ninj2 + ninj
  ni2 = ni + ni
  ni3 = ni2 + ni

  y0 = cm133*y*(y-one)*(y-two)
  y1 = cp5*(y+one)*(y-one)*(y-two)
  y2 = cm5*y*(y+one)*(y-two)
  y3 = cp133*y*(y+one)*(y-one)

  do i = 0,3
    va4(i) = src(i    )*abcd(0) + src(i    +ninj)*abcd(1) +  src(i    +ninj2)*abcd(2) + src(i    +ninj3)*abcd(3)
    vb4(i) = src(i+ni )*abcd(0) + src(i+ni +ninj)*abcd(1) +  src(i+ni +ninj2)*abcd(2) + src(i+ni +ninj3)*abcd(3)
    vc4(i) = src(i+ni2)*abcd(0) + src(i+ni2+ninj)*abcd(1) +  src(i+ni2+ninj2)*abcd(2) + src(i+ni2+ninj3)*abcd(3)
    vd4(i) = src(i+ni3)*abcd(0) + src(i+ni3+ninj)*abcd(1) +  src(i+ni3+ninj2)*abcd(2) + src(i+ni3+ninj3)*abcd(3)
    dst(i) = va4(i)*y0 + vb4(i)*y1 + vc4(i)*y2 + vd4(i)*y3
  enddo

  x0 = cm133*x*(x-one)*(x-two)
  x1 = cp5*(x+one)*(x-one)*(x-two)
  x2 = cm5*x*(x+one)*(x-two)
  x3 = cp133*x*(x+one)*(x-one)

  tricub_x86_f = dst(0)*x0 + dst(1)*x1 + dst(2)*x2 + dst(3)*x3

  return
  end


