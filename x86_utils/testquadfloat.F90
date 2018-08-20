program testquad
  implicit none
  integer, parameter :: ARRSIZE = 400000000
  real *4, dimension(ARRSIZE) :: r4
  real *4  :: s4
  real *8  :: s8, r8, sc
  real *16 :: sf, rf
  integer :: i
  real *8, external :: comp_sum

  do i = 1,ARRSIZE
   r4(i) = (ARRSIZE + 1 - i) / 31415.9
  enddo
  print 101,'min. max',r4(ARRSIZE), r4(1)
  s4 = 0 
  s8 = 0
  sf = 0
  do i = 1,ARRSIZE
    s4 = s4 + r4(i)
    r8 = r4(i)
    s8 = s8 + r8
    rf = r4(i)
    sf = sf + rf
  enddo
  sc = comp_sum(r4,ARRSIZE)
  print 101,'sc sf s8 s4',sc,sf,s8,s4
101 format(A,4G25.16)
  print 101,'s8 -s4', s8 -s4
  print 101,'sc -s4', sc -s4
  print 101,'sc -s8', sc -s8
  print 101,'sc -sf', sc -sf
end 

real*8 function comp_sum(z,n)
  implicit none
  real*4, intent(IN), dimension(n) :: z
  integer, intent(IN) :: n
  real*8 :: s, c, t, y
  integer :: j
  s = 0
  c = 0
  do j=1,n
    y = z(j) - c
    t = s
    s = s + y
    c = (s - t) - y
  enddo
  comp_sum = s
101 format(A,4G25.16)
  print 101,'s, c',s,c
  return
end