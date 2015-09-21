module qsort_mod
 
implicit none
 
type group
    integer :: order    ! original order of unsorted data
    real :: value       ! values to be sorted by
end type group
integer, save :: depth, maxdepth
contains
 
recursive subroutine QSort(a,indx,na,dir)
 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
type (group), dimension(nA), intent(inout) :: A
integer, dimension(nA), intent(inout) :: indx
character(len=*) :: dir
 
! LOCAL VARIABLES
integer :: left, right
real :: random
real :: pivot
type (group) :: temp
integer :: marker, index, tmp
    index = -1
!    if(na == 1)print *,'depth=',depth,index,nA,dir
    if (nA > 1) then
    depth = depth + 1
    maxdepth=max(maxdepth,depth)
 
        call random_number(random)
        index = int(random*real(nA-1))+1
        index = 1
!        index = na/2
!        print *,'depth=',depth,index,nA,dir
        pivot = A(index)%value   ! random pivot (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right)%value > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%value < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
                tmp = indx(left)
                indx(left) = indx(right)
                indx(right) = tmp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
!        call QSort(A(:marker-1),marker-1,' left')
!        call QSort(A(marker:),nA-marker+1,' right')
        call QSort(A,indx,marker-1,' left')
        call QSort(A(marker),indx(marker),nA-marker+1,' right')
 
    depth = depth - 1
    end if
 
end subroutine QSort
 
recursive subroutine QSortI(A,indx,na,dir)
 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
integer, dimension(nA), intent(in) :: A
integer, dimension(nA), intent(inout) :: indx
character(len=*) :: dir
 
! LOCAL VARIABLES
integer :: left, right
real :: random
real :: pivot
type (group) :: temp
integer :: marker, index, tmp
    if (nA > 1) then
    depth = depth + 1
    maxdepth=max(maxdepth,depth)
 
        call random_number(random)
        index = int(random*real(nA-1))+1
        index = ishft(na,-1)
!        print *,'depth=',depth,index,nA,dir
        pivot = A(indx(index))   ! random pivot (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(indx(right)) > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(indx(left)) < pivot)
                left = left + 1
            end do
            if (left < right) then
                tmp = indx(left)
                indx(left) = indx(right)
                indx(right) = tmp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
!        call QSort(A(:marker-1),marker-1,' left')
!        call QSort(A(marker:),nA-marker+1,' right')
        call QSortI(A,indx,marker-1,' left')
        call QSortI(A,indx(marker),nA-marker+1,' right')
 
    depth = depth - 1
    end if
 
end subroutine QSortI
 
recursive subroutine QSortR(A,indx,na,dir)
 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
real, dimension(nA), intent(in) :: A
integer, dimension(nA), intent(inout) :: indx
character(len=*) :: dir
 
! LOCAL VARIABLES
integer :: left, right
real :: random
real :: pivot
type (group) :: temp
integer :: marker, index, tmp
    index = -1
!    if(na == 1)print *,'depth=',depth,index,nA,dir
    if (nA > 1) then
    depth = depth + 1
    maxdepth=max(maxdepth,depth)
 
        call random_number(random)
        index = int(random*real(nA-1))+1
        index = na
!        print *,'index,na,pivot=',index,na,A(indx(index))
!        index = na/2
!        print *,'depth=',depth,index,nA,dir
        pivot = A(indx(index))   ! random pivot (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(indx(right)) > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(indx(left)) < pivot)
                left = left + 1
            end do
            if (left < right) then
                tmp = indx(left)
                indx(left) = indx(right)
                indx(right) = tmp
            end if
        end do
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
!        call QSort(A(:marker-1),marker-1,' left')
!        call QSort(A(marker:),nA-marker+1,' right')
        marker = min(na,marker)
        call QSortR(A,indx,marker-1,' left')
        call QSortR(A,indx(marker),nA-marker+1,' right')
 
    depth = depth - 1
    end if
 
end subroutine QSortR
 
end module qsort_mod
 
! Test Qsort Module
program qsort_test
use qsort_mod
implicit none
 
integer, parameter :: l = 100000
type (group), dimension(l) :: A
integer, dimension(l) :: indx
real, dimension(l) :: B
integer, dimension(l) :: C
integer, dimension(12) :: seed = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
integer :: i
real :: random
 
    write (*,*) "Unsorted Values:"
    call random_seed(put = seed)
    do i = 1, l
        call random_number(random)
        A(i)%value = random
        if(i == l) a(i)%value = .999999
        B(I) = A(i)%value
        A(i)%order = i
        c(i) = 1000000*B(i)
        indx(i) = i
        if(l>32) cycle
        if (mod(i,4) == 0) write (*,"(4(I5,1X,F8.6))") A(i-3:i)
        if (mod(i,4) == 0) write (*,"(4(6x,F8.6))") B(i-3:i)
        if (mod(i,4) == 0) write (*,"(4(6x,I8))") C(i-3:i)
    end do
    depth = 0
    maxdepth = 0
    call QSortR(B,indx,l,' center')
    write (*,*) "Sorted Values: , depth",maxdepth
    do i = 4, l, 4
        if(l>32) cycle
        if (mod(i,4) == 0) write (*,"(4(6X,F8.6))") B(indx(i-3:i))
    end do
    if(l<=32) print 100,indx
    do i = 1, l
      indx(i) = i
    enddo
    depth = 0
    maxdepth = 0
    call QSortI(C,indx,l,' center')
    write (*,*) "Sorted Values: , depth",maxdepth
    do i = 4, l, 4
        if(l>32) cycle
        if (mod(i,4) == 0) write (*,"(4(6X,I8))") C(indx(i-3:i))
    end do
    if(l<=32) print 100,indx
    do i = 1, l
      indx(i) = i
    enddo
    depth = 0
    depth = 0
    maxdepth = 0
    call QSort(A,indx,l,' center')
!    print *,'depth=',depth
    write (*,*) "Sorted Values: , depth",maxdepth
    do i = 4, l, 4
        if(l>32) cycle
        if (mod(i,4) == 0) write (*,"(4(I5,1X,F8.6))") A(i-3:i)
    end do
    if(l<=32) print 100,indx
100 format(20I5)
 
end program qsort_test
