program test_xxhash
  use ISO_C_BINDING
  implicit none
  include 'xxhash.inc'
  interface
    function rdtsc() result(tim) BIND(C,name='rdtsc')
      import :: C_LONG
      integer(C_LONG) :: tim
    end function rdtsc
  end interface
  integer, parameter :: ARRAY_SZ=1024*1024*4
  integer :: i, j, k, ilen
  integer, dimension(23) :: answer
  integer, dimension(0:ARRAY_SZ-1), target :: buffer
  integer *8 :: t0, t1
  real *8, dimension(128) :: t
  integer, dimension(23) :: iln

  do i=0,ARRAY_SZ-1
    buffer(i) = i + 133
  enddo
  ilen = ARRAY_SZ*4
  do i = 1,23
    iln(i) = ilen
    ilen = ilen/2
  enddo
  print 100, iln(1:23)
100 format(23I9)
  do j = 0, 9
    ilen = ARRAY_SZ*4
    do k=1,23
      t0 = rdtsc()
      if(ilen > 15) then
	answer(k) = Hash32(12345, c_loc(buffer(0)), ilen);
      else
	answer(k) = Hash32_short(c_loc(buffer(0)), ilen);
      endif
      t1 = rdtsc()
102 format(23Z9.8)
      t(k) = (t1-t0)
      t(k) = t(k)/ilen
      ilen = ilen / 2
    enddo
    print 101,t(1:23)
    if(j==0 .or. j==9) print 102,answer(1:23)
101 format(23F9.3)
    buffer(0) = j * 1234
  enddo
end program test_xxhash
