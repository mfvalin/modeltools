program demo
  use ISO_C_BINDING
  implicit none
  include 'randomfunctions.inc'

!   interface
!     integer(C_LONG) function time() bind(C,name='time')
!     import :: C_LONG
!     end function time
!   end interface

  integer, parameter :: LSTREAMS = 22000000
  integer, parameter :: NSTREAMS = 8
  integer, parameter :: NREP = 1
  type(RANDOM_STREAM), dimension(NSTREAMS) :: stream
  real(C_DOUBLE), dimension(LSTREAMS,NSTREAMS) :: g
  type(RANDOM_STREAM) :: null_stream
  integer(C_INT) :: piSeed
  integer(C_INT) :: cSeed
  integer :: i, j, k, indx
  real(C_DOUBLE) :: gauss
  integer, dimension(-5:5,NSTREAMS) :: dist
  real*8 :: t0, t1, t2, t3
  real*8 :: omp_get_wtime

  do i = 1, NSTREAMS
    cSeed = 1
    piSeed = 123456 + mod(i,2)
    call Ran_R250_new_stream( stream(i),null_stream, [ piSeed ], cSeed)
    call RanNormalZigSetSeed(stream(i), [ piSeed ], cSeed)
  enddo

  dist = 0
  t0 = omp_get_wtime()
!$OMP PARALLEL private(i,j,indx) shared(stream,g,dist)
!$OMP DO
  do j = 1,NSTREAMS
   do k = 1, NREP
    do i = 1, LSTREAMS
      g(i,j) = DRanS_R250_stream(stream(j)) * 5.5   ! (-5.5 , 5.5) range
      indx =nint( g(i,j) )
      dist(indx,j) = dist(indx,j)+1
    enddo
   enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
  print *,'----------------- 11 bins, values from -5.5 to 5.5 ----------------'
  print *,'-------------------- uniform distribution test  -------------------'
  print *,'expecting',LSTREAMS/11*NREP,' samples per interval'
  print 101,(dist(:,j),j=1,NSTREAMS)
  t1 = omp_get_wtime()
  print *, 'time =',t1-t0
  print *,'-------------------- gaussian distribution test -------------------'

  dist = 0
!$OMP PARALLEL private(i,j,indx) shared(stream,g,dist)
!$OMP DO
  do j = 1,NSTREAMS
   do k = 1, NREP
    do i = 1, LSTREAMS
      g(i,j) = DRan_Normal_stream(stream(j))
      indx =nint(  max(-5.4,min(5.4,g(i,j))) )   ! clamp extreme values to force [-5,+5] range for indx
      dist(indx,j) = dist(indx,j)+1
    enddo
   enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
  print 101,(dist(:,j),j=1,NSTREAMS)
  t2 = omp_get_wtime()
  print *, 'time =',t2-t1
101 format(11I10)
  stop
end

  
