program wpd
  use m_precision, only: wp
  use m_constants, only: const
  use m_splines, only: linspace, spline_easy
  use m_vib_finite_diff, only: solve_finite_diff
  use m_fourier_grid, only: solve_fourier_grid
  use m_fourier_grid, only: fourier_grid_d1_fast
  use m_KH_functions, only: solve_sinc_DVR
  use m_KH_functions, only: dsinc_mat_elems
  use m_wave_packet_dynamics, only: wpd_eigenfun_z
  use m_algebra, only: matmul_Adagger_x_z

  implicit none 

  character(75)::inputfile_i, inputfile_n, outputfile
  integer::nin,npoints,nstates

  ! loop variables
  integer::i,j

  real(kind=wp)::mu_in, mu, dx
    real(wp), allocatable:: x_i(:), y_i(:), x_n(:), y_n(:)
  real(wp), allocatable:: x_new(:), y_i_new(:), y_n_new(:)
  real(wp), allocatable:: eigval_i(:), eigval_n(:)
  complex(wp), allocatable:: eigvec_z_i(:,:), eigvec_z_n(:,:)
  complex(wp), allocatable:: S_ni(:)
  complex(wp), allocatable:: c_l_t(:,:)
  
  integer:: t, l, n, ntimes
  real(wp), allocatable:: times(:)
  character(80):: file
  real(wp), allocatable:: x_mean(:)

  ! read from standard input
  read(5,*) inputfile_i
  read(5,*) inputfile_n
  read(5,*) outputfile
  read(5,*) nin, npoints, nstates
  read(5,*) mu_in
  
  allocate( x_i(nin), y_i(nin), x_n(nin), y_n(nin), x_new(npoints), &
         y_i_new(npoints), y_n_new(npoints))
  allocate(S_ni(npoints))
  
  ! read inputfile_gs
  open(10,file=inputfile_i,status='old')
  do i=1,nin
    read(10,*) x_i(i),y_i(i)
  enddo
  close(10)

  open(10,file=inputfile_n,status='old')
  do i=1,nin
    read(10,*) x_n(i),y_n(i)
  enddo
  close(10)

  mu = mu_in * const % u
  y_i = y_i -minval(y_i)
  y_i =y_i * const % hartree
  x_i =x_i * 1.0e-10
  y_n = y_n -minval(y_n)
  y_n =y_n * const % hartree
  x_n =x_n * 1.0e-10
  
  ! spline, make more points
  call linspace(x_new, x_i(1), x_i(nin), npoints)
  dx = x_new(2)-x_new(1) 

  call spline_easy(x_i, y_i, nin, x_new, y_i_new, npoints)
  call spline_easy(x_n, y_n, nin, x_new, y_n_new, npoints)
  
  ! fourier grid, state i
  allocate(eigvec_z_i(npoints, npoints))
  allocate(eigval_i(npoints))
  call solve_fourier_grid(dx, y_i_new, eigval_i, eigvec_z_i, mu, "SI", "fast")  
  
  write(6,*) eigval_i(1:nstates) * const % cm 
  write(6,*) "solve_fourier_grid: fundamental frequency", (eigval_i(2)-eigval_i(1)) * const % cm

  ! fourier grid, state n
  allocate(eigvec_z_n(npoints, npoints))
  allocate(eigval_n(npoints))
  call solve_fourier_grid(dx, y_n_new, eigval_n, eigvec_z_n, mu, "SI", "fast")  
  
  write(6,*) eigval_n(1:nstates) * const % cm 
  write(6,*) "solve_fourier_grid: fundamental frequency", (eigval_n(2)-eigval_n(1)) * const % cm

  ! create intial wave packet, S_ni = sum_l <n | l> <l | i> = sum_l c^*_ln c_li 
  ! slow, shoudl use matmul
  !S_ni = 0.0d0
  !do l=1, size(eigvec_z_n, 1)
  !  do n=1, size(eigvec_z_n, 2)  
  !    S_ni(n) =  S_ni(n) + conjg(eigvec_z_n(l,n)) * eigvec_z_i(l,1)  
  !  end do
  !end do

  ! this is faster
  call matmul_Adagger_x_z(eigvec_z_n, eigvec_z_i(:,1), S_ni)
    
  ! now do wpd
  ntimes=20
  allocate(times(ntimes),  c_l_t(npoints, ntimes))
  call linspace(times, 0.0d0, 20.0d-15, ntimes)
  
  do t = 1, size(times)
    call wpd_eigenfun_z(S_ni, eigval_n / const % hbar, eigvec_z_n, times(t), c_l_t(:,t))
  end do

  ! compute the mean position of the squared wave packet as a function of time
  allocate(x_mean(ntimes))
  do t = 1, size(times)
    x_mean(t) = sum(x_new(:) * abs(c_l_t(:,t))**2) / sum(abs(c_l_t(:,t))**2)
  end do
  
  open(11,file="x_mean_t.txt",status='unknown')

  do t = 1, size(times)
    write(11,*) times(t), x_mean(t) 
  end do
  
  close(11)
   
   
  ! write wave packet file
  do t=1, ntimes
    write(file, *) t-1 
    file = "wpd_t_" //  trim(adjustl(file))  // ".txt"
    open(11,file=file,status='unknown')

    do l=1, size(x_i)
      write(11,*) x_i(l), real(c_l_t(l,t),8), aimag(c_l_t(l,t)), abs(c_l_t(l,t))
    end do
    close(11)
  end do
  
  
  
end program wpd
