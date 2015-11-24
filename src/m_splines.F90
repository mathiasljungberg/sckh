module m_splines
  use m_precision, only: wp
  implicit none
  
!  interface 
!     SUBROUTINE spline(x,y,n,yp1,ypn,y2)
!       INTEGER :: n,NMAX
!       REAL(8) :: yp1,ypn,x(n),y(n),y2(n)
!       PARAMETER (NMAX=50000)
!       INTEGER :: i,k
!       REAL(8) :: p,qn,sig,un,u(NMAX)
!     end SUBROUTINE spline
!  end interface
!  interface
!     SUBROUTINE splint(xa,ya,y2a,n,x,y)
!       INTEGER :: n
!       REAL(8) :: x,y,xa(n),y2a(n),ya(n)
!       INTEGER :: k,khi,klo
!       REAL(8) :: a,b,h
!     end SUBROUTINE splint
!  end interface  

contains

  
  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
    INTEGER, intent(in):: n !,NMAX
    REAL(8), intent(in):: yp1,ypn,x(n),y(n)
    REAL(8), intent(out):: y2(n)
    !PARAMETER (NMAX=50000)

    INTEGER:: i,k
    REAL(8):: p,qn,sig,un
    real(8), allocatable:: u(:)

    allocate(u(n))

    if (yp1.gt..99d30) then
      y2(1)=0.d0
      u(1)=0.d0
    else
      y2(1)=-0.5d0
      u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    
    do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.d0
      y2(i)=(sig-1.d0)/p
      u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1) &
      -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
      u(i-1))/p
    end do
    
    if (ypn.gt..99d30) then
      qn=0.d0
      un=0.d0
    else
      qn=0.5d0
      un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

    do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
    end do

  end SUBROUTINE spline

  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    INTEGER n
    REAL(8) x,y,xa(n),y2a(n),ya(n)

    INTEGER k,khi,klo
    REAL(8) a,b,h
    klo=1
    khi=n

    do while (khi-klo.gt.1)
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
    end do

    h=xa(khi)-xa(klo)
      
    ! if (abs(h) .le. 1.0d-20) then
    ! if (h .eq. 0.0d0) then
    !  write(6,*) "bad xa input in splint"
    !  stop
    !endif
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0

  END SUBROUTINE splint

  ! this subroutine splines a whole array, last index time dimension
  subroutine splint_array3_d(xa,ya,y2a,x,y)
    real(8), intent(in):: x, xa(:),y2a(:,:,:),ya(:,:,:)
    real(8), intent(out):: y(:,:) 
    integer:: k
    integer:: khi,klo
    real(8):: a,b,h
    klo=1
    khi= size(xa)

    do while (khi-klo.gt.1)
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
    end do

    h=xa(khi)-xa(klo)

    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y(:,:)=a*ya(:,:,klo)+b*ya(:,:,khi)+((a**3-a)*y2a(:,:,klo)+(b**3-b)*y2a(:,:,khi))*(h**2)/6.d0

  end subroutine splint_array3_d



subroutine linspace(x,start,end,npoints)
  !passed variables
  integer, intent(in):: npoints
  real(kind=wp), intent(in):: start,end
  real(kind=wp), intent(out),dimension(npoints)::x

  ! local variables
  integer::i
  real(kind=wp)::dx
  
  if (npoints .ne. 1) then
     dx = (end-start)/dfloat(npoints-1)
     
     ! write(6,*) dx, end,start,end-start,dfloat(npoints-1),npoints
     
     do i=1,npoints
        x(i)=start +(i-1)*dx
     end do
  else
     x(1) = start
  end if

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

end module m_splines
