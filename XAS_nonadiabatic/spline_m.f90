module spline_m
  use parameters
  implicit none
  
  interface 
     SUBROUTINE spline(x,y,n,yp1,ypn,y2)
       INTEGER :: n,NMAX
       REAL(8) :: yp1,ypn,x(n),y(n),y2(n)
       PARAMETER (NMAX=50000)
       INTEGER :: i,k
       REAL(8) :: p,qn,sig,un,u(NMAX)
     end SUBROUTINE spline
  end interface
  interface
     SUBROUTINE splint(xa,ya,y2a,n,x,y)
       INTEGER :: n
       REAL(8) :: x,y,xa(n),y2a(n),ya(n)
       INTEGER :: k,khi,klo
       REAL(8) :: a,b,h
     end SUBROUTINE splint
  end interface  

contains

subroutine linspace(x,start,end,npoints)
  !passed variables
  integer, intent(in):: npoints
  real(kind=wp), intent(in):: start,end
  real(kind=wp), intent(out),dimension(npoints)::x

  ! local variables
  integer::i
  real(kind=wp)::dx
  
  dx = (end-start)/dfloat(npoints-1)

  
 ! write(6,*) dx, end,start,end-start,dfloat(npoints-1),npoints

  do i=1,npoints
     x(i)=start +(i-1)*dx
  end do
  
end subroutine linspace

subroutine spline_one(x,y,n, x2,y2)
 
  !passed variables
  integer,intent(in):: n
  real(kind=wp), intent(in), dimension(n):: x,y
  real(kind=wp), intent(in) :: x2
  real(kind=wp), intent(out) :: y2
  
  !local variables
  real(kind=wp),dimension(n):: yder
  integer::i

  ! call spline, natural splines
  call spline(x,y,n,1.0d30,1.0d30,yder)

  ! call splint, for all x2 make y2
  call splint(x,y,yder,n,x2,y2)

end subroutine spline_one

subroutine spline_easy(x,y,n,x2,y2,n2)
  
  !passed variables
  integer,intent(in):: n,n2
  real(kind=wp), intent(in), dimension(n) :: x,y
  real(kind=wp), intent(in),dimension(n2) :: x2
  real(kind=wp), intent(out),dimension(n2) :: y2
  
  !local variables
  real(kind=wp),dimension(n):: yder
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

end module spline_m
