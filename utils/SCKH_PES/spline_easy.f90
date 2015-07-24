subroutine spline_easy(x,y,n,x2,y2,n2)
  implicit none
  
  !passed variables
  integer,intent(in):: n,n2
  real(8), intent(in), dimension(n) :: x,y
  real(8), intent(in),dimension(n2) :: x2
  real(8), intent(out),dimension(n2) :: y2
  
  !local variables
  real(8),dimension(n):: yder
  integer::i
  
  ! This routine works like octaves spline function. 
  !
  !   x :  input array of tabulated x-values
  !   y :  input array of tabulated y-values
  !   n :  length of x1,y1 
  !   x2 : input array of new x-values
  !   y2 : input array of new interpolated y-values
  !   n2 : length of x2,y2 
  ! 
  !   The subroutine uses spline.f and slpint.f from numerical recepies
  !   modified for real(8) 
  !   
  !   Mathias Ljungberg, 2006-09-21
  !
  
  ! call spline, natural splines
  call spline(x,y,n,1.0d30,1.0d30,yder)

  ! call splint, for all x2 make y2
  do i=1,n2
     call splint(x,y,yder,n,x2(i),y2(i))
  end do
  
end subroutine spline_easy
