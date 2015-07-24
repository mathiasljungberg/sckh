program vib_finite_diff
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile
  integer::nin,npoints,nstates,alt3
  real(kind=wp), dimension(:),allocatable:: x,y

  ! loop variables
  integer::i,j

  ! other variables
  real(kind=wp), dimension(:),allocatable:: x_SI,y_SI,x_new_SI,y_new_SI, &
       diag,subdiag
  real(kind=wp)::dx, my, my_SI, my_in
  
  ! lapack
  real(kind=wp):: abstol,ifail,info
  integer::n_eigfound
  integer, dimension(:),allocatable::iwork
  real(kind=wp), dimension(:),allocatable:: eigenval,work
  real(kind=wp), dimension(:,:),allocatable:: eigenvec

  !functions
  real(kind=wp)::dlamch


  ! read from standard input
  read(5,*) inputfile
  read(5,*) outputfile
  read(5,*) nin, npoints, nstates
  read(5,*) alt3, my_in


  ! allocate
  allocate(iwork(5*npoints), work(5*npoints), eigenval(npoints), &
       eigenvec(npoints,nstates))

  allocate( x(nin),y(nin), x_SI(nin),y_SI(nin),x_new_SI(npoints), &
       y_new_SI(npoints),diag(npoints),subdiag(npoints))

  ! read inputfile
  open(10,file=inputfile,status='old')
  
  do i=1,nin
     read(10,*) x(i),y(i)
  enddo

  close(10)

  ! my

  my= 18._wp/19._wp
!  write(6,*) "my", my
  if(alt3.eq.1) my=my_in
!  write(6,*) "my", my
  my_SI = my*amu


  ! go to SI units
  y = y -minval(y)
  y_SI=y*hartree
  x_SI=x*1.0e-10
  
  !dx = x1_SI(2)-x1_SI(1)

  ! spline, make more points
  call linspace(x_new_SI, x_SI(1), x_SI(nin), npoints)

  call spline_easy(x_SI, y_SI, nin, x_new_SI, y_new_SI, npoints)

  ! test old spline routines
  open(10,file="splines_old.txt",status='unknown')
  do i=1, npoints
    write(10,*) i, x_new_SI(i), y_new_SI(i)
  end do
  close(10)

  dx = x_new_SI(2)-x_new_SI(1)

  ! sätt upp hamiltonianen -hbar^2/(2*m) (d/dx)^2 + V(x) 
  
  ! diagonalen
  do i=1,npoints
     diag(i) = (-hbar**2/(dx**2*2*my_SI) )*(-2) + y_new_SI(i)
    ! write(6,*) diag
  end do

  ! subdiagonalen
  do i=1,npoints-1
     subdiag(i) = (-hbar**2/(dx**2*2*my_SI) )*(1)
  end do

  subdiag(npoints)=0
  
  ! solve eigenvalue problem

  abstol=2d0*dlamch('s')

  call dstevx('v','i',npoints,diag,subdiag, 0d0,1d0, 1, nstates, abstol, &
       n_eigfound, eigenval , eigenvec, npoints,work,iwork,ifail,info)


  !open outputfile
  open(10,file=outputfile,status='unknown')
  
  ! write eigenvectors  
  do i=1,npoints
     write(10,*) x_new_SI(i), (eigenvec(i,j), j=1,2)
  end do
  
  close(10)
  
  ! write eigenvalues
  write(6,*) (eigenval(i)*cm, i=1,nstates) 

  ! write fundamental frequency
  write(6,*)
  write(6,*) "Fundamental frequency ", (eigenval(2)-eigenval(1))*cm
  
end program vib_finite_diff
