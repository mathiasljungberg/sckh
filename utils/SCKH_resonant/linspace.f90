subroutine linspace(x,start,end,npoints)
  implicit none

  !passed variables
  integer, intent(in):: npoints
  real(8), intent(in):: start,end
  real(8), intent(out),dimension(npoints)::x

  ! local variables
  integer::i
  real(8)::dx
  
  dx = (end-start)/dfloat(npoints-1)

  
 ! write(6,*) dx, end,start,end-start,dfloat(npoints-1),npoints

  do i=1,npoints
     x(i)=start +(i-1)*dx
  end do
  
end subroutine linspace
