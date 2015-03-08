! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd
! Modified by Mathias Ljungberg 2008

module qsort_c_module
  use parameters
  implicit none
  public :: qsort
  private :: Partition,QsortC

contains
  
  subroutine qsort(A,B)
    real(kind=wp), intent(inout), dimension(:) :: A
    integer,intent(out), dimension(:) :: B
    
    integer:: i
    
    ! make array containing indicies 
    do i=1,size(A)
       B(i)=i
    end do

    call QsortC(A,B)

  end subroutine qsort
  
 recursive subroutine QsortC(A,B)
  real(kind=wp), intent(inout), dimension(:) :: A
  integer,intent(inout), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq, B)
     call QsortC(A(:iq-1), B(:iq-1))
     call QsortC(A(iq:), B(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker,B)
  real(kind=wp), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer,intent(inout), dimension(:) :: B  
  integer :: i, j, temp_int
  real(kind=wp) :: temp
  real(kind=wp) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ! do the same to B(i) and B(j)
        temp_int = B(i)
        B(i) = B(j)
        B(j) = temp_int
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module

!program sortdriver
!  use qsort_c_module
!  implicit none
!  integer, parameter :: r = 10
!  real, dimension(1:r) :: myarray = &        ! (1:r)
!     (/0, 50, 20, 25, 90, 10, 5, 1, 99, 75/)
!  print *, "myarray is ", myarray
!  call QsortC(myarray)
!  print *, "sorted array is ", myarray
!end program sortdriver
