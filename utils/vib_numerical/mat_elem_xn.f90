subroutine mat_elem_xn(A,n,xn)
  implicit none
  
  !passed variables
  integer, intent(in)::  n, xn
  !real(8), dimension(n,n), intent(in):: A
  real(8), dimension(n,n), intent(out):: A

  !local variables
  integer:: i,j,dim
  real(8), dimension(:,:), allocatable::C1,C2,C3
  
  ! this subroutine calculates the matrix elements of x^n with harmonic 
  ! oscillator basis functions, i.e. < m | x^n | k >, xn is the exponent 
  ! n is the number of basis functions. The returned matrix A should
  ! be multiplied by (hbar/(2*m*omega))^(xn/2) 
  !
  ! Mathias Ljungberg, 2006-10-02


  ! check dimensions, allocate
  if(n>=xn+1) then
     dim=n
  else
     dim=xn+1
  end if

  allocate(C1(dim,dim),C2(dim,dim),C3(dim,dim)) 


  !initialize C3
  C3=0
  do i=1,dim
     C3(i,i) = 1
  end do

  ! do raising and lowering operations
  do i=1,xn
     call raising(C3,C1,dim)
     call lowering(C3,C2,dim)
     C3 = C1 + C2
  end do

  !copy back to A, of dimension n*n
  A=0
  do i=1,n
     do j=1,n
        
        A(i,j) = C3(i,j)
     end do
  end do


deallocate(C1,C2,C3)

end subroutine mat_elem_xn



subroutine raising(A,B,n)
  implicit none
  
  !passed variables
  integer, intent(in)::  n
  real(8), dimension(n,n), intent(in):: A
  real(8), dimension(n,n), intent(out):: B

  !local variables
  integer:: i,j

  B=0
  
  do i=1,n
     do j=1,n-1
        
        B(j+1,i) = sqrt(dfloat(j))*A(j,i) 

     end do
  end do

end subroutine raising



subroutine lowering(A,B,n)
  implicit none
  
  !passed variables
  integer, intent(in)::  n
  real(8), dimension(n,n), intent(in):: A
  real(8), dimension(n,n), intent(out):: B

  !local variables
  integer:: i,j
  
  B=0

  do i=1,n
     do j=2,n
        
        B(j-1,i) = sqrt(dfloat(j-1))*A(j,i) 

     end do
  end do


end subroutine lowering
