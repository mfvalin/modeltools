program demo
  use ISO_C_BINDING
  implicit none
  include 'randomfunctions.inc'

  integer, parameter :: LSTREAMS = 40000000
  integer, parameter :: NSTREAMS = 12
  type(RANDOM_STREAM), dimension(NSTREAMS) :: stream
  real(C_DOUBLE), dimension(LSTREAMS,NSTREAMS) :: g
  type(RANDOM_STREAM) :: null_stream
  integer(C_INT) :: piSeed
  integer(C_INT) :: cSeed
  integer :: i, j, indx
  real(C_DOUBLE) :: gauss
  integer, dimension(-5:5,NSTREAMS) :: dist

  do i = 1, NSTREAMS
    cSeed = 1
    piSeed = 123456 + mod(i,2)
    stream(i) = Ran_R250_new_stream(null_stream, [ piSeed ], cSeed)
    call RanNormalSetSeedZigVec(stream(i), [ piSeed ], cSeed)
  enddo

  dist = 0
!$OMP PARALLEL private(i,j,indx) shared(stream,g,dist)
!$OMP DO
  do j = 1,NSTREAMS
    do i = 1, LSTREAMS
      g(i,j) = DRanS_R250_stream(stream(j)) * 5.0
      indx =nint(  max(-5.4,min(5.4,g(i,j))) )
      dist(indx,j) = dist(indx,j)+1
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
  print 101,(dist(:,j),j=1,NSTREAMS)
  print *,'---------------------------------------'

  dist = 0
!$OMP PARALLEL private(i,j,indx) shared(stream,g,dist)
!$OMP DO
  do j = 1,NSTREAMS
    do i = 1, LSTREAMS
      g(i,j) = DRanNormalZigVec(stream(j))
      indx =nint(  max(-5.4,min(5.4,g(i,j))) )
      dist(indx,j) = dist(indx,j)+1
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
  print 101,(dist(:,j),j=1,NSTREAMS)
101 format(11I10)
  stop
end

  
