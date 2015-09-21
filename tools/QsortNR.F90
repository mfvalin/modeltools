! An iterative implementation of quick sort 
! A --> Array to be sorted, l  --> Starting index, h  --> Ending index, indx index for sort
subroutine quickSortIterativeI (arr, indx, ll, hh)
  implicit none
  integer, intent(IN) :: ll, hh
  integer, dimension(1:hh), intent(IN) :: arr
  integer, dimension(1:hh), intent(INOUT) :: indx
! Create an auxiliary stack
  integer, dimension( hh - ll + 1 ) ::  stack
  integer :: top, maxdepth, pcalls, p, h, l
!  integer, external :: partition

! initialize top of stack
  top = 0
  maxdepth = 0
  pcalls = 0

! push initial values of l and h to stack
  top = top + 1 ; stack( top ) = ll
  top = top + 1 ; stack( top ) = hh

! Keep popping from stack while is not empty
  do while ( top >= 1 )
    ! Pop h and l
    h = stack( top ); top = top -1
    l = stack( top ); top = top -1

    ! Set pivot element at its correct position in sorted array
    p = partition( arr, indx, l, h );
    pcalls = pcalls + 1

    ! If there are elements on left side of pivot, then push left
    ! side to stack
    if ( p-1 > l ) then
      top = top + 1 ; stack( top ) = l
      top = top + 1 ; stack( top ) = p - 1
    endif

    ! If there are elements on right side of pivot, then push right
    ! side to stack
    if ( p+1 < h ) then
      top = top + 1 ; stack( top ) = p + 1
      top = top + 1 ; stack( top ) = h
    endif
    maxdepth = max(top , maxdepth )
  enddo
!  print *,'maxdepth=',maxdepth,' pcalls=',pcalls
contains
  integer function partition (arr, indx, l, h)
    implicit none
    integer, intent(IN) :: l, h
    integer, dimension(1:hh), intent(IN) :: arr
    integer, dimension(1:hh), intent(INOUT) :: indx

    integer :: x, i, j, temp, temp2

!    x = arr( (l+h)/2 );
    x = arr( indx((l+h)/2) );
    i = l;
    j = h;
    do while(.true.)
!      do while(arr(j) > x)
      do while(arr(indx(j)) > x)
        j = j - 1
      enddo
!      do while(arr(i) < x)
      do while(arr(indx(i)) < x)
        i = i + 1
      enddo
      if(i<j) then
!        if(arr(i) == x .and. arr(j) == x) then
        if(arr(indx(i)) == x .and. arr(indx(j)) == x) then
          i = i + 1
        else 
!          temp  = arr(i)  ; arr(i) = arr(j)   ; arr(j) = temp
          temp2 = indx(i) ; indx(i) = indx(j) ; indx(j) = temp2
        endif
      else
        partition = j
        return
      endif
    enddo
  end function partition
end subroutine quickSortIterativeI
 
subroutine quickSortIterativeR (arr, indx, ll, hh)
  implicit none
  integer, intent(IN) :: ll, hh
  real, dimension(1:hh), intent(IN) :: arr
  integer, dimension(1:hh), intent(INOUT) :: indx
! Create an auxiliary stack
  integer, dimension( hh - ll + 1 ) ::  stack
  integer :: top, maxdepth, pcalls, p, h, l
!  integer, external :: partition

! initialize top of stack
  top = 0
  maxdepth = 0
  pcalls = 0

! push initial values of l and h to stack
  top = top + 1 ; stack( top ) = ll
  top = top + 1 ; stack( top ) = hh

! Keep popping from stack while is not empty
  do while ( top >= 1 )
    ! Pop h and l
    h = stack( top ); top = top -1
    l = stack( top ); top = top -1

    ! Set pivot element at its correct position in sorted array
    p = partition( arr, indx, l, h );
    pcalls = pcalls + 1

    ! If there are elements on left side of pivot, then push left
    ! side to stack
    if ( p-1 > l ) then
      top = top + 1 ; stack( top ) = l
      top = top + 1 ; stack( top ) = p - 1
    endif

    ! If there are elements on right side of pivot, then push right
    ! side to stack
    if ( p+1 < h ) then
      top = top + 1 ; stack( top ) = p + 1
      top = top + 1 ; stack( top ) = h
    endif
    maxdepth = max(top , maxdepth )
  enddo
!  print *,'maxdepth=',maxdepth,' pcalls=',pcalls
contains
  integer function partition (arr, indx, l, h)
    implicit none
    integer, intent(IN) :: l, h
    real, dimension(1:hh), intent(IN) :: arr
    integer, dimension(1:hh), intent(INOUT) :: indx

    integer :: x, i, j, temp, temp2

!    x = arr( (l+h)/2 );
    x = arr( indx((l+h)/2) );
    i = l;
    j = h;
    do while(.true.)
!      do while(arr(j) > x)
      do while(arr(indx(j)) > x)
        j = j - 1
      enddo
!      do while(arr(i) < x)
      do while(arr(indx(i)) < x)
        i = i + 1
      enddo
      if(i<j) then
!        if(arr(i) == x .and. arr(j) == x) then
        if(arr(indx(i)) == x .and. arr(indx(j)) == x) then
          i = i + 1
        else 
!          temp  = arr(i)  ; arr(i) = arr(j)   ; arr(j) = temp
          temp2 = indx(i) ; indx(i) = indx(j) ; indx(j) = temp2
        endif
      else
        partition = j
        return
      endif
    enddo
  end function partition
end subroutine quickSortIterativeR
 
! Driver program to test above functions
program test
  implicit none
  integer, parameter :: N = 133
  integer, dimension(N) :: arr, indx, orig, indx2
  real, dimension(N) :: r
  integer :: i

  arr = 10
  arr(1:33) = [4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 12]
  arr(101:133) = [4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 4, 3, 5, 2, 1, 3, 2, 3, 12]
  orig = arr
  R = arr
  print 100,arr
  do i=1,N
    indx(i) = i
    indx2(i) = i
  enddo
  print *,'============'
  print 100,(arr(indx(i)),i=1,N)
  call  quickSortIterativeI( arr, indx, 1, N )
  print *,'============'
  print 100,(orig(indx(i)),i=1,N)
  call  quickSortIterativeR( R, indx2, 1, N )
  print *,'============'
  print 101,(R(indx2(i))-orig(indx2(i)),i=1,N)
100 format (20I4)
101 format (20F4.0)
  stop
end
