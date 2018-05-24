program test_semilagx3
  use ISO_C_BINDING
  implicit none
#define NI 300
#define NJ 200
#define NK 85
#define fn(x,y,z) ((x)/128.+1.5) * ((y)/128.+2.5) * ((z)/32.+3.5)

  real, dimension(3,NI,NJ,NK) :: Z
  real, dimension(3) :: expected, base, res
  real, dimension(3) :: offsets
  real, dimension(4) :: abcd
  integer :: i, j, k, l
  real :: dx, dy, dz
  real :: px, py, pz
  integer :: ix, iy, iz
interface
  subroutine tricub_x86_f3(dd, ff, abcd, x, y, ni, nj) bind(C,name='tricub_x86_f3')
    import :: C_FLOAT, C_INT
    real(C_FLOAT),  intent(OUT), dimension(3)   :: dd
    real(C_FLOAT),  intent(IN),  dimension(3,*) :: ff
    real(C_FLOAT),  intent(IN),  dimension(4)   :: abcd
    real(C_FLOAT),  intent(IN),  value :: x, y
    integer(C_INT), intent(IN),  value :: ni, nj
  end subroutine tricub_x86_f3
end interface

  offsets = [0., 10., 100.]
  do k = 1, NK
  do J = 1, NJ
  do i = 1, NI
    Z(1,i,j,k) = fn(i,j,k) + offsets(1)
    Z(2,i,j,k) = Z(1,i,j,k) + offsets(2)
    Z(3,i,j,k) = Z(1,i,j,k) + offsets(3)
  enddo
  enddo
  enddo
  do k = 1, NK-1
    pz = k + .5
    dz = .5
    if(k == 1 .or. k == NK-1) then
      iz = k
      abcd(1) = 1. - dz
      abcd(2) = dz
      abcd(3:4) = 0.
    else
      iz = k - 1
    endif
    do j = 2, NJ - 2
      py = j + .5
      iy = j - 1
      do i = 2, NI - 2
        px = i + .5
        ix = i - 1
        base = fn(px, py, pz)
        expected = base + offsets
        call tricub_x86_f3(res, z(1,ix,iy,iz), abcd, dx, dy, NI, NJ)
      enddo
    enddo
   enddo
end program
