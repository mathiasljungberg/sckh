program vib_finite_diff
  use m_precision, only: wp
  use m_constants, only: const
  use m_splines, only: linspace, spline_easy
  use m_vib_finite_diff, only: solve_finite_diff
  use m_fourier_grid, only: solve_fourier_grid
  use m_fourier_grid, only: solve_fourier_grid_real
  use m_fourier_grid, only: solve_hartley_grid
  use m_fourier_grid, only: hartley_grid_d1_fast
  use m_fourier_grid, only: fourier_grid_d1_fast
  use m_fourier_grid, only: fourier_grid_d1_real_fast
  use m_KH_functions, only: solve_sinc_DVR
  use m_KH_functions, only: dsinc_mat_elems
  implicit none 

  character(75)::inputfile,outputfile
  integer::nin,npoints,nstates

  ! loop variables
  integer::i,j

  real(wp)::mu_in, mu, dx
  real(wp), allocatable:: x(:), y(:), x_new(:), y_new(:)
  real(wp), allocatable:: eigvec(:,:), eigval(:)
  complex(wp), allocatable:: eigvec_z(:,:)
  complex(wp), allocatable:: first_der(:,:)
  real(wp), allocatable:: first_der_sinc(:,:)
  real(wp), allocatable:: first_der_real(:,:)

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

  ! open outputfile
  open(11,file=outputfile,status='unknown')

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

  write(6,*) eigval * const % cm 
  write(6,*) "solve_finite_diff: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  do i=1, nstates
    write(11,*) eigval(i) * const % cm
  end do
  
  ! fourier grid, slow but reliable
  allocate(eigvec_z(npoints, npoints))
  deallocate(eigval)
  allocate(eigval(npoints))
  call solve_fourier_grid(dx, y_new, eigval, eigvec_z, mu, "SI", "fast")  
  
  write(6,*) eigval(1:nstates) * const % cm 
  write(6,*) "solve_fourier_grid: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  do i=1, nstates
    write(11,*) eigval(i) * const % cm
  end do
  
  ! fourier grid with sin cos, fast
  eigval=0.0_wp
  eigvec=0.0_wp
  call solve_fourier_grid_real(dx, y_new, eigval, eigvec, mu, "SI", "fast")  
  
  write(6,*) eigval(1:nstates) * const % cm 
  write(6,*) "solve_fourier_grid_real fast: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  do i=1, nstates
    write(11,*) eigval(i) * const % cm
  end do
  
  ! hartely grid with cos + sin, fast
  eigval=0.0_wp
  eigvec=0.0_wp
  call solve_hartley_grid(dx, y_new, eigval, eigvec, mu, "SI", "fast")  
  
  write(6,*) eigval(1:nstates) * const % cm 
  write(6,*) "solve_hatrley_grid fast: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  do i=1, nstates
    write(11,*) eigval(i) * const % cm
  end do
  
  ! sinc
  eigval=0.0_wp
  eigvec=0.0_wp
  call solve_sinc_DVR(dx, mu, y_new, eigvec, eigval)

  write(6,*) eigval(1:nstates) * const % cm 
  write(6,*) "solve_sinc_DVR: fundamental frequency", (eigval(2)-eigval(1)) * const % cm

  do i=1, nstates
    write(11,*) eigval(i) * const % cm
  end do

  ! look at first derivatives
  allocate(first_der(npoints, npoints))
  call fourier_grid_d1_fast(npoints, dx, first_der)

  write(6,*) "first derivative matrix", abs(first_der * const % cm * const % hbar * dcmplx(0.0_wp, -1.0_wp))

  ! look at first derivatives
  allocate(first_der_real(npoints, npoints))
  call fourier_grid_d1_real_fast(npoints, dx, first_der_real)

  write(6,*) "first derivative fourier real", first_der_real * const % cm * const % hbar

  ! look at first derivatives
  !allocate(first_der_real(npoints, npoints))
  first_der_real = 0.0_wp
  call hartley_grid_d1_fast(npoints, dx, first_der_real)

  write(6,*) "first derivative Hartley", first_der_real * const % cm * const % hbar 


  allocate(first_der_sinc(npoints, npoints))
  call dsinc_mat_elems(first_der_sinc,x_new)

  write(6,*) "first derivaive sinc", first_der_sinc * const % cm * const % hbar 
  
  write(6,*) "differnce abs first derivaives", (abs(first_der_sinc) - abs(first_der)) * const % cm * const % hbar 
  
  
end program vib_finite_diff
