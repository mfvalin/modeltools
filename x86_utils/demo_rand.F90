program demo
  use ISO_C_BINDING
  implicit none
  include 'randomfunctions.inc'

  type(RANDOM_STREAM) :: stream
  type(RANDOM_STREAM) :: null_stream
  integer(C_INT) :: piSeed
  integer(C_INT) :: cSeed
  integer :: i, indx
  real(C_DOUBLE) :: gauss
  integer, dimension(-5:5) :: dist

  cSeed = 1
  piSeed = 123456
  dist = 0
  stream = Ran_R250_new_stream(null_stream, [ piSeed ], cSeed)
  call RanNormalSetSeedZigVec(stream, [ piSeed ], cSeed)

  do i = 1,1000000000
    gauss = DRanS_R250_stream(stream) * 5.0
    indx =nint(  max(-5.4,min(5.4,gauss)) )
    dist(indx) = dist(indx)+1
  enddo
  print 101,dist
  dist = 0
  do i = 1,1000000000
    gauss = DRanNormalZigVec(stream)
    indx =nint(  max(-5.4,min(5.4,gauss)) )
    dist(indx) = dist(indx)+1
  enddo
  print 101,dist
101 format(11I10)
  stop
end

  
