program vib_finite_diff
  use m_precision, only: wp
  use m_constants, only: const
  use m_splines, only: linspace, spline_easy
  use m_vib_finite_diff, only: solve_finite_diff
  use m_fourier_grid, only: solve_fourier_grid
  implicit none 

  character(75)::inputfile,outputfile
  integer::nin,npoints,nstates

  ! loop variables
  integer::i,j

  real(kind=wp)::mu_in, mu, dx
  real(wp), allocatable:: x(:), y(:), x_new(:), y_new(:)
  real(wp), allocatable:: eigvec(:,:), eigval(:)
  complex(wp), allocatable:: eigvec_z(:,:)

  ! read from standard input
  read(5,*) inputfile
  read(5,*) outputfile
  read(5,*) nin, npoints, nstates
  read(5,*) mu_in
  
  allocate( x(nin),y(nin),  x_new(npoints), &
         y_new(npoints))
  
  ! read inputfile
  open(10,file=inputfile,status='old')
  
  do i=1,nin
     read(10,*) x(i),y(i)
  enddo

  close(10)

  mu = mu_in * const % u
  y = y -minval(y)
  y =y * const % hartree
  x =x * 1.0e-10

  ! spline, make more points
  call linspace(x_new, x(1), x(nin), npoints)
  dx = x_new(2)-x_new(1) 
  call spline_easy(x, y, nin, x_new, y_new, npoints)
  
  ! finite differences
  call solve_finite_diff(dx, y_new, nstates, eigval, eigvec, mu, "SI")  
  
  write(6,*) "solve_finite_diff: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  ! fourier grid
  allocate(eigvec_z(npoints, npoints))
  deallocate(eigval)
  allocate(eigval(npoints))
  call solve_fourier_grid(dx, y_new, eigval, eigvec_z, mu, "SI")  
  
  write(6,*) "solve_finite_diff: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  

end program vib_finite_diff
