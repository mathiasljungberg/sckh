subroutine vib_finite_diff_subr(x_in, y_in, nin, npoints, nstates, my_in,freq_out )
  use parameters
  implicit none

  ! passed variables
  integer,intent(in)::nin,npoints,nstates
  real(kind=wp), dimension(nin),intent(in):: x_in,y_in
  real(kind=wp),intent(in)::  my_in
  real(kind=wp),intent(out):: freq_out

  ! loop variables
  integer::i,j

  ! other variables
  real(kind=wp), dimension(:),allocatable:: x,y, x_SI,y_SI,x_new_SI,y_new_SI, &
       diag,subdiag
  real(kind=wp)::dx,  my_SI
    
  
  ! lapack
  real(kind=wp):: abstol
  integer::n_eigfound,info
  integer, dimension(:),allocatable::iwork,ifail
  real(kind=wp), dimension(:),allocatable:: eigenval,work
  real(kind=wp), dimension(:,:),allocatable:: eigenvec

  !functions
  real(kind=wp)::dlamch

  ! my_in is the diatomic reduced mass

  ! allocate
  allocate(iwork(5*npoints), work(5*npoints), eigenval(npoints), &
       eigenvec(npoints,nstates),ifail(npoints))
  allocate(x(nin),y(nin),  x_SI(nin),y_SI(nin),x_new_SI(npoints), &
       y_new_SI(npoints),diag(npoints),subdiag(npoints-1))


  my_SI = my_in*amu
  x = x_in
  y = y_in

  ! go to SI units
  y = y -minval(y)
  y_SI=y*hartree
  x_SI=x*1.0d-10
     
  ! spline, make more points
  call linspace(x_new_SI, x_SI(1), x_SI(nin), npoints)
  dx = x_new_SI(2)-x_new_SI(1) 
  call spline_easy(x_SI, y_SI, nin, x_new_SI, y_new_SI, npoints)

  ! s�tt upp hamiltonianen -hbar^2/(2*m) (d/dx)^2 + V(x) 
  
  ! diagonalen
  do i=1,npoints
     diag(i) = (-hbar**2/(dx**2*2*my_SI) )*(-2.0_wp) + y_new_SI(i)
  end do

  ! subdiagonalen
  do i=1,npoints-1
     subdiag(i) = -hbar**2/(dx**2 *2.0_wp *my_SI) 
  end do

  !subdiag(npoints)=0
  
  ! solve eigenvalue problem
  abstol=2d0*dlamch('s')
  call dstevx('v','i',npoints,diag,subdiag, 0.d0,1.d0, 1, nstates, abstol, &
       n_eigfound, eigenval , eigenvec, npoints,work,iwork,ifail,info)


  freq_out = (eigenval(2)-eigenval(1))*cm


! deallocate
  deallocate(iwork, work, eigenval, &
       eigenvec,ifail)
  deallocate(x,y,  x_SI,y_SI,x_new_SI, &
       y_new_SI,diag,subdiag)
  
end subroutine vib_finite_diff_subr
