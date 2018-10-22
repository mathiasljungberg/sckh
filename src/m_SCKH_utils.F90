module m_SCKH_utils

  implicit none
  
contains

  subroutine sample_x_mom(x, funct, x_sampl, mom_sampl, x_mom_sampl, mode)
    use m_precision, only: wp
    use m_constants, only: const
    use m_FFT, only: SFFTEU, next_power_of_2, reorder_sigma1

    real(wp), dimension(:), intent(in):: x, funct
    real(wp), dimension(:), intent(out):: x_sampl, mom_sampl
    real(wp), dimension(:,:), intent(out)::x_mom_sampl
    integer, intent(in):: mode
    !local variables
    integer::  npoints_x_sampl, npoints, npoints_mom_sampl,npoints_x_mom_sampl, &
         npoints_pad, npoints_pad_pow
    real(wp),dimension(:),allocatable:: funct_real, funct_imag, x_mom
    integer::i,j,k
    real(wp):: dx, pes_l

    npoints = size(x) 
    npoints_x_sampl = size(x_sampl) 
    npoints_mom_sampl = size(mom_sampl) 
    npoints_x_mom_sampl =  size(x_mom_sampl,1)

    dx=x(2)-x(1)

    call next_power_of_2(npoints, npoints_pad, npoints_pad_pow)

    pes_l = (npoints_pad -1) * dx

    allocate( funct_real(npoints_pad), funct_imag(npoints_pad), &
         x_mom(npoints_pad))
    
    ! alt 1: do an 'even' sampling of N points
    ! 1) calcualte integral = I(x) of |c_i|^2
    ! 2) use inverse function x = I^-1(y) spline points where I = (n+0.5) / N, n=0,N-1
    
    ! space: sample from funct ** 2
    
    if(mode .eq. 1) then
       call sample_even(x, funct ** 2, x_sampl)
    else if(mode .eq. 2) then
       call sample_random(x, funct ** 2, x_sampl)
    else
       write(6,*) "error in sample_x_mom, mode must be either 1 or 2"
       stop
    end if

    do i=1,npoints
       write(16,*) x(i), funct(i) ** 2
    end do
       
    ! momentum: fourier transform funct and sample F(q) ** 2
    funct_real = 0
    funct_real(1:npoints) =funct(1:npoints)
    funct_imag = 0
    
    call SFFTEU( funct_real, funct_imag, npoints_pad, npoints_pad_pow, 1 )   
    
    j=1
    do i=npoints_pad/2, 1, -1
       x_mom(j) = -2.0_wp * const % pi * i * const % hbar / pes_l 
       j=j+1
    end do

    do i=0, npoints_pad/2-1
       x_mom(j) =  2.0_wp * const % pi * i * const % hbar / pes_l 
       j=j+1
    end do

    call reorder_sigma1(funct_real)
    call reorder_sigma1(funct_imag)

    do i=1,npoints_pad
       write(19,'(5F16.6)') x_mom(i), funct_real(i), funct_imag(i), &
            funct_real(i) ** 2 + funct_imag(i) ** 2, &
            (funct_real(i) ** 2 + funct_imag(i) ** 2) * 4.0_wp * const % pi * x_mom(i)**2
    end do

    if(mode .eq. 1) then
       call sample_even(x_mom, funct_real ** 2 &
            + funct_imag ** 2, mom_sampl)
       
       k=1
       do i=1,npoints_x_sampl
          do j=1, npoints_mom_sampl
             x_mom_sampl(k,1) = x_sampl(i)
             x_mom_sampl(k,2) = mom_sampl(j)
             k=k+1
          end do
       end do
    
    else if(mode .eq. 2) then
       call sample_random(x_mom, funct_real ** 2 &
            + funct_imag ** 2, mom_sampl)
       
       do i=1,npoints_x_mom_sampl
          x_mom_sampl(i,1) = x_sampl(i)
          x_mom_sampl(i,2) = mom_sampl(i)
       end do

    else
       write(6,*) "error in sample_x_mom, mode must be either 1 or 2"
       stop   
    end if
    
   deallocate(funct_real, funct_imag, &
        x_mom)
    
  end subroutine sample_x_mom

  subroutine sample_x_mom2(x, funct, x_sampl, x_mom_sampl, E_i_inp, E_gs, my_SI)
    use m_precision, only: wp
    use m_splines, only: spline_easy
    !use m_constants, only: const

    real(wp), dimension(:), intent(in):: x, funct, E_i_inp
    real(wp), dimension(:), intent(out):: x_sampl
    real(wp), dimension(:,:), intent(out):: x_mom_sampl
    real(wp):: E_gs, my_SI
    !local variables
    integer:: npoints, npoints_x_sampl, npoints_x_mom_sampl
    real(wp),dimension(:),allocatable:: V_sampl
    integer::i

    ! sample space and assign momentum according to conservation of energy
    ! does not seem to work very well! especially since the quantum space probablility distribution
    ! lies outside of the classically allowed region
    
    npoints = size(x)
    npoints_x_sampl = size(x_sampl) 
    npoints_x_mom_sampl = size(x_mom_sampl,1)

    if(npoints_x_mom_sampl .ne. 2*npoints_x_sampl) then
       write(6,*) "error in dimension in sample_x_mom2", npoints_x_mom_sampl, 2*npoints_x_sampl
    end if
    
    allocate(V_sampl(npoints_x_sampl))
    
    ! alt 1: do an 'even' sampling of N points
    ! 1) calcualte integral = I(x) of |c_i|^2
    ! 2) use inverse function x = I^-1(y) spline points where I = (n+0.5) / N, n=0,N-1
    
    ! space: sample from funct ** 2
    call sample_even(x, funct ** 2, x_sampl)
     
    do i=1,npoints
       write(16,*) x(i), funct(i) ** 2
    end do
     
    ! Conservation of energy: E = p^2 / 2m +V => p = +-sqrt(2*m*(E -V) )        
    call spline_easy(x, E_i_inp, npoints, x_sampl, V_sampl, npoints_x_sampl)

    do i=1,npoints_x_sampl
       write(6,*) x_sampl(i), V_sampl(i), E_gs
    end do

    !p = sqrt( 2*my_SI*(E_gs - V_sampl))
    do i=1,npoints_x_sampl
       x_mom_sampl(i,1) = x_sampl(i)
       if( E_gs - V_sampl(i) .gt. 0d0) then
          x_mom_sampl(i,2) = sqrt( 2.0_wp * my_SI * (E_gs - V_sampl(i)) )
       else
          x_mom_sampl(i,2) = 0
       end if
    end do
    do i=1,npoints_x_sampl !+1,npoints_x_sampl*2
       x_mom_sampl(i + npoints_x_sampl,1) = x_sampl(i)
       if( E_gs - V_sampl(i) .gt. 0d0) then
          x_mom_sampl(i + npoints_x_sampl,2) = -sqrt( 2.0_wp * my_SI * (E_gs - V_sampl(i) ))  
       else
          x_mom_sampl(i + npoints_x_sampl,2) = 0
       end if
    end do

    deallocate(V_sampl)
  end subroutine sample_x_mom2

  
  subroutine sample_even(x, funct, x_sampl)
    use m_precision, only: wp
    use m_splines, only: spline_easy, linspace

    !use m_constants, only: const

    ! passed variables
    real(wp), dimension(:), intent(in)::x, funct
    real(wp), dimension(:), intent(out)::x_sampl
    ! local variables
    integer,parameter:: npoints_new=10000
    integer:: npoints, npoints_sampl
    real(wp), dimension(npoints_new)::funct_spl, I_x, x_new
    real(wp), dimension(:),allocatable::points_sampl
    integer:: i,j
    
    npoints =size(x) 
    npoints_sampl = size(x_sampl) 
   
    allocate(points_sampl(npoints_sampl))
    
    do i=1,npoints
       if(funct(i) .lt. 0) then
          write(6,*) "function has a negative value in sample_even!", i, funct(i)
          stop
       end if
     end do
    

    ! first spline function to npoints_new points, then integrate up to
    call linspace(x_new, x(1), x(npoints), npoints_new)
    call spline_easy(x, funct, npoints, x_new, funct_spl, npoints_new)

    funct_spl =abs(funct_spl) ! safeguard against bad splining

    I_x=0
    I_x(1) = funct_spl(1) 
    do i=2,npoints_new
       I_x(i) = I_x(i-1) + funct_spl(i)
    end do
    
    do i=1,npoints_sampl
       points_sampl(i) = ( (i-0.5D0) * I_x(npoints_new)) / npoints_sampl
    end do

    do i=1,npoints_sampl
       do j=1,npoints_new
          if(I_x(j) .le. points_sampl(i)) then
             x_sampl(i) =x_new(j)
          else
             exit
          end if
       end do
    end do

    deallocate(points_sampl)
    
  end subroutine sample_even

  subroutine sample_random(x, funct, x_sampl)
    use m_precision, only: wp
    use m_splines, only: spline, splint
    !use m_constants, only: const

    ! passed variables
    real(wp), dimension(:), intent(in)::x, funct
    real(wp), dimension(:), intent(out)::x_sampl
    ! local variables
    integer,parameter:: npoints_new=10000
    integer:: npoints, npoints_sampl
    real(wp), dimension(:),allocatable::yder
    real(wp):: funct_x_tmp, x_tmp, y_tmp, tol, x_min, x_max,y_min, y_max, x_l, y_l,ran1,ran2
    integer:: i,j
    
    npoints = size(x)
    npoints_sampl = size(x_sampl) 
   
    allocate(yder(npoints))
    
    do i=1,npoints
       if(funct(i) .lt. 0) then
          write(6,*) "function has a negative value in sample_random!", i, funct(i)
          stop
       end if
    end do
    
    tol = maxval(funct) * 5.0d-3

    !find the interval of sampling
    y_min = 0
    y_max = maxval(funct)
    y_l = y_max - y_min

    do i=1,npoints
       if(funct(i) .le. tol) then
          x_min = x(i)
       else
          exit
       end if
    end do

    do i=npoints, 1, -1
       if(funct(i) .le. tol) then
          x_max = x(i)
       else
          exit
       end if
    end do

    x_l = x_max - x_min

    write(6,*) "x_min, x_max, x_l", x_min, x_max, x_l

    ! sample!

    j=1
    call random_seed

    call spline(x,funct,npoints,1.0d30,1.0d30,yder) 
   
    do while(j .le. npoints_sampl )
       call random_number(ran1)
       call random_number(ran2)
       
       x_tmp = x_l * ran1 + x_min
       y_tmp = y_l * ran2 + y_min

       !call spline_one(x, funct, npoints, x_tmp, funct_x_tmp)
       call splint(x,funct,yder,npoints,x_tmp, funct_x_tmp)

       if(y_tmp .le. funct_x_tmp) then
          x_sampl(j) = x_tmp
          j=j+1
       end if
    end do

    deallocate(yder)

  end subroutine sample_random

  ! wrapper routine for sample_x_mom
  subroutine sample_x_mom_modes(npoints_x_sampl, npoints_mom_sampl, samplemode, X_dvr, c_i, x_mom_sampl)
    use m_precision, only: wp
    use m_io, only: get_free_handle
    !use m_sckh_params_t, only: sckh_params_t 

    
    integer, intent(in):: npoints_x_sampl
    integer, intent(inout):: npoints_mom_sampl
    integer, intent(in):: samplemode
    real(wp), intent(in):: X_dvr(:)
    real(wp), intent(in):: c_i(:)
    real(wp), intent(out), allocatable:: x_mom_sampl(:,:)
    
    real(wp), allocatable:: x_sampl(:)
    real(wp), allocatable:: mom_sampl(:)
    integer:: npoints_x_mom_sampl
    integer:: i, ifile
    
    if(samplemode .eq. 1) then
      npoints_x_mom_sampl = npoints_x_sampl * npoints_mom_sampl
    else if(samplemode .eq. 2) then
      npoints_x_mom_sampl = npoints_x_sampl !* 2 !npoints_x_sampl * npoints_mom_sampl
      npoints_mom_sampl =  npoints_x_sampl
    else
      write(6,*) "samplemode muste be either 1 or 2!"
     stop
   end if
   
   allocate(x_sampl(npoints_x_sampl))
   allocate(mom_sampl(npoints_mom_sampl))
   allocate(x_mom_sampl(npoints_x_mom_sampl, 2))
   
   if(samplemode .eq. 1) then  
     call sample_x_mom(X_dvr, c_i(:), x_sampl, mom_sampl, x_mom_sampl, 1)
   else if(samplemode .eq. 2) then  
     call sample_x_mom(X_dvr, c_i(:), x_sampl, mom_sampl, x_mom_sampl, 2)
   end if
   
   write(6,*) "Sampling done"
   
   !check
   ifile = get_free_handle()
   open(ifile, file="sampled_points.txt", action='write')
   do i=1,npoints_x_sampl
     write(ifile,*) x_sampl(i)
   end do
   close(ifile)
   
   ifile = get_free_handle()
   open(ifile, file="x_mom_sampl.txt", action='write')
   do i=1,npoints_x_mom_sampl
     ! modified by O.Takahashi 2018/06/29
!     write(22,*) x_mom_sampl(i,1), x_mom_sampl(i,2)
     write(ifile,*) x_mom_sampl(i,1), x_mom_sampl(i,2)
   end do
   close(ifile)
   
 end subroutine sample_x_mom_modes


!  subroutine verlet_trajectory(x_in, v_in, V_x, V_y, dt, my_SI, x_out)
!    use m_precision, only: wp
!    !use m_constants, only: const
!
!    implicit none
!    ! passed variables
!    real(wp), dimension(:), intent(in):: V_x, V_y
!    real(wp), intent(in):: x_in,v_in, dt, my_SI    
!    real(wp), dimension(:), intent(out)::x_out
!        
!    ! local variables
!    integer:: ntsteps, V_len
!    integer:: i,j
!    real(wp):: a, x, v, x_new, v_new, a_new,t
!    
!    V_len =size(V_x)     
!    ntsteps = size(x_out)
!    
!    ! first force evaluation on PES
!    x = x_in
!    v = v_in
!    a = force(x, V_x, V_y ) / my_SI
!
!    ! write first (zeroth) step to trajectory
!    t=0
!    !call trajectory_add(traj,x,v,a,t)
!    
!    i=0
!    !write(14,'(I4,4ES16.6)') i, t, x, v, a
!    ! run n_steps steps
!    j=0
!    do i =1,ntsteps
!       call verlet_step(x,v,a, V_x, V_y, V_len, x_new, v_new, a_new, dt, my_SI)
!       
!       t=t+dt
!       x = x_new
!       v= v_new
!       a = a_new
!       
!       j=j+1
!       x_out(j) = x 
!
!       !write(14,'(I4,4ES16.6)') i, t, x, v, a
!       !call trajectory_add(traj,x,v,a,t)
!       
!    end do
!    
!  end subroutine verlet_trajectory

  subroutine verlet_trajectory(x_in, v_in, V_x, V_y, dt, my_SI, x_out)
    use m_precision, only: wp
    !use m_constants, only: const

    implicit none
    ! passed variables
    real(wp), dimension(:), intent(in):: V_x, V_y
    real(wp), intent(in):: x_in,v_in, dt, my_SI    
    real(wp), dimension(:), intent(out)::x_out
        
    ! local variables
    integer:: ntsteps, V_len
    integer:: i
    real(wp):: a, x, v, x_new, v_new, a_new,t
    
    V_len =size(V_x)     
    ntsteps = size(x_out)
    
    ! first force evaluation on PES
    x = x_in
    v = v_in
    a = force(x, V_x, V_y ) / my_SI

    t = 0.0_wp

    x_out(1) = x

    do i =2,ntsteps
       call verlet_step(x,v,a, V_x, V_y, V_len, x_new, v_new, a_new, dt, my_SI)
       
       t=t+dt
       x = x_new
       v= v_new
       a = a_new

       x_out(i) = x 
       
    end do
    
  end subroutine verlet_trajectory

  
  subroutine verlet_trajectory_xva(x_in, v_in, V_x, V_y, dt, my_SI, x_out, v_out, a_out, &
       use_abs_bc_in, abs_bc_max_x_in)
    use m_precision, only: wp
    !use m_constants, only: const

    implicit none
    ! passed variables
    real(wp), intent(in):: V_x(:), V_y(:)
    real(wp), intent(in):: x_in,v_in, dt, my_SI    
    real(wp), intent(out):: x_out(:)
    real(wp), intent(out):: v_out(:)
    real(wp), intent(out):: a_out(:)
    logical, intent(in), optional::  use_abs_bc_in
    real(wp), intent(in), optional::  abs_bc_max_x_in 
    
    ! local variables
    integer:: ntsteps, V_len
    integer:: i
    real(wp):: a, x, v, x_new, v_new, a_new,t

    logical::  use_abs_bc
    real(wp)::  abs_bc_max_x

    if(present(use_abs_bc_in)) then
      use_abs_bc =use_abs_bc_in
    else
      use_abs_bc = .false.
    end if

    if(present(abs_bc_max_x_in)) then
      abs_bc_max_x = abs_bc_max_x_in
    else
      abs_bc_max_x = 1.0d100
    end if
    
    V_len =size(V_x)     
    ntsteps = size(x_out)
    
    ! first force evaluation on PES
    x = x_in
    v = v_in
    a = force(x, V_x, V_y ) / my_SI

    t = 0.0_wp

    x_out(1) = x
    v_out(1) = v
    a_out(1) = a 
    
    do i =2,ntsteps

      ! absorbing boundary conditions
      if(use_abs_bc) then
        if(x .gt. abs_bc_max_x) then
          write(6,*) "absorbing boundary at x", abs_bc_max_x, "reached, exit loop"
          x_out(i:ntsteps) = x
          v_out(i:ntsteps) = 0.0_wp
          a_out(i:ntsteps) = 0.0_wp
          exit
        end if
      end if
      
      call verlet_step(x,v,a, V_x, V_y, V_len, x_new, v_new, a_new, dt, my_SI)
       
       t=t+dt
       x = x_new
       v= v_new
       a = a_new
       
       !j=j+1
       !x_out(j) = x
       !v_out(j) = v
       !a_out(j) = a
       x_out(i) = x
       v_out(i) = v
       a_out(i) = a 
    end do
    
  end subroutine verlet_trajectory_xva

  
  
  subroutine verlet_step(x,v,a, V_x, V_y, V_len, x_new, v_new, a_new, dt, my_SI)
    use m_precision, only: wp
    !use m_constants, only: const

    ! passed variables
    real(wp), intent(in):: x,v,a,my_SI, dt
    real(wp), intent(out):: x_new,v_new,a_new
    integer, intent(in)::V_len
    real(wp), dimension(V_len), intent(in):: V_x,V_y

    ! velocity verlet
    ! 1) a(t) = F/m  
    ! 2) calculate x(t+dt) 
    !    x(t+dt)  = x(t) + v(t)*dt + 0.5*a(t)*dt**2  
    ! 3) again evaluate force, now at time t+dt
    ! 4) calcualte v(t+dt)
    !    v(t+dt) = v(t) +0.5*(a(t) +a(t+dt))*dt
    
    x_new = x + v * dt + 0.5d0 * a * dt ** 2
    a_new  = force(x_new, V_x, V_y) / my_SI
    v_new = v +0.5d0 * (a + a_new) * dt  
    
    !write(6,*) x_new, a_new, v_new
    
  end subroutine verlet_step
  
  
  real(wp) function force(x, V_x, V_y, dx_in)
    use m_precision, only: wp
    use m_splines, only: spline_easy
    
    !passed
    real(wp), dimension(:), intent(in):: V_x,V_y
    real(wp), intent(in) :: x
    real(wp), intent(in), optional::dx_in
    !local
    integer::V_len
    real(wp), dimension(2):: x2, V2 
    real(wp)::dx    

    if (present(dx_in)) then
       dx=dx_in
    else
       dx=1.0d-13
    end if

    V_len = size(V_x) 
    
    ! spline two nearby points on he PES
    ! calculmate force by central finite difference formula:
    ! F = - d/dx V(x(t))  \approx - (V(x(t)+h/2) - V(x(t)-h/2) ) / h   
    
    x2(1) = x - dx/2
    x2(2) = x + dx/2
    
    call spline_easy(V_x, V_y, V_len, x2, V2, 2)

    force = - (V2(2) -V2(1)) / dx
  
  end function force
  
  
  
  subroutine compute_SCKH(E_n, E_f, E_fn_mean, D_fn, time, sigma_m, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_FFT, only: FFT_complex 

    ! passed variables
    real(wp), dimension(:), intent(in):: E_n,time
    real(wp), dimension(:,:), intent(in):: E_f
    real(wp),dimension(:,:,:),intent(in):: D_fn 
    complex(wp), dimension(:,:,:),intent(out) ::  sigma_m
    real(wp), intent(in):: E_fn_mean, gamma
    ! local variables
    integer:: ntsteps, ntsteps_pad, nfinal
    complex(wp), dimension(:),allocatable:: funct
    complex(wp), dimension(:,:),allocatable::  e_factor1
    real(wp), dimension(:,:),allocatable:: int_W_I
    real(wp),dimension(:),allocatable::  omega_out 
    real(wp):: delta_t
    integer:: i,k,m 

    ntsteps = size(time) 
    nfinal = size(E_f,1) 
    ntsteps_pad = size(sigma_m,2) 
    
    allocate(funct(ntsteps), int_W_I(nfinal, ntsteps), &
         e_factor1(nfinal,ntsteps), omega_out(ntsteps_pad) )
    
   
    delta_t = time(2)-time(1)
        
    int_W_I(:,1) = E_f(:,1) - E_n(1)  - E_fn_mean  
    
    do i = 2, ntsteps
       int_W_I(:,i) = int_W_I(:,i-1) + ( E_f(:,i) - E_n(i)) - E_fn_mean  
    end do
    int_W_I =  int_W_I * delta_t
    
    do i = 1,nfinal 
       e_factor1(i,:) = exp(dcmplx(0, -(const % eV  / const % hbar) *int_W_I(i,:)  ))
    end do
    
    ! compute A_{fm}(t) = D^m_{fn}(t) * e_factor_f(t) * exp(-gamma * eV * time(:) / hbar)
    ! and fourer transform
    
    sigma_m = 0.0_wp
    
    !do a FFT
    do k =1,nfinal! final state 
       do m=1,3 ! polarization     
          
          funct = D_fn(k,:,m) * e_factor1(k,:) * exp(-gamma * const % eV * time(:) / const % hbar)

          call FFT_complex(time,funct, sigma_m(k,:,m), omega_out)
          
       end do ! m
    end do ! k
    
    deallocate(funct, int_W_I)
    
  end subroutine compute_SCKH

!  ! same as above, but use FFTW to make things more standard  
!  subroutine compute_F_if_omp_m(E_n, E_f, E_nf_mean, D_fn, time, F_if_omp_m, gamma)
!    use m_precision, only: wp
!    use m_constants, only: const
!    !use m_fftw3, only: fft_c2c_1d_forward
!    use m_fftw3, only: fft_c2c_1d_backward
!    use m_fftw3, only: reorder_sigma_fftw_z
!    
!    real(wp), dimension(:), intent(in):: E_n,time
!    real(wp), dimension(:,:), intent(in):: E_f
!    real(wp),dimension(:,:,:),intent(in):: D_fn 
!    complex(wp), dimension(:,:,:),intent(out) ::  F_if_omp_m
!    real(wp), intent(in):: E_nf_mean, gamma
!
!    integer:: ntsteps, nfinal
!    complex(wp), dimension(:),allocatable:: funct
!    complex(wp), dimension(:,:),allocatable::  e_factor1
!    integer:: f_e, m 
!
!    ntsteps = size(time) 
!    nfinal = size(E_f,1) 
!    
!    allocate(funct(ntsteps), &
!         e_factor1(nfinal,ntsteps))
!   
!    do f_e = 1,nfinal 
!      call compute_efactor(E_n, E_f(f_e,:), E_nf_mean, time, e_factor1(f_e,:), .true.)
!    end do
!  
!    F_if_omp_m = 0.0_wp
!    
!    do f_e =1,nfinal! final state 
!      do m=1,3 ! polarization     
!        
!        funct = D_fn(f_e,:,m) * e_factor1(f_e,:) * exp(-gamma * const % eV * time(:) / const % hbar)
!        
!        !call fft_c2c_1d_forward(funct, F_if_omp_m(f_e,:,m))
!        call fft_c2c_1d_backward(funct, F_if_omp_m(f_e,:,m))
!        call reorder_sigma_fftw_z(F_if_omp_m(f_e,:,m))
!
!      end do ! m
!    end do ! f_e
!    
!  end subroutine compute_F_if_omp_m

  ! mimics old nonresonant routine, only one intermediate state
  subroutine compute_F_if_omp(E_n, E_f, E_nf_mean, D_fn, D_ni, time, F_if_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), dimension(:), intent(in):: E_n,time
    real(wp), dimension(:,:), intent(in):: E_f
    real(wp),dimension(:,:,:),intent(in):: D_fn 
    real(wp),dimension(:),intent(in):: D_ni 
    complex(wp), dimension(:,:,:,:),intent(out) ::  F_if_omp
    real(wp), intent(in):: E_nf_mean, gamma
    
    complex(wp), allocatable::  e_factor(:)
    integer:: nfinal, f_e, ntsteps

    ntsteps = size(time)     
    
    allocate(e_factor(ntsteps))

    nfinal = size(E_f,1) 
    
    do f_e =1,nfinal! final state 
      !call compute_F_ifn_omp(E_n, E_f(f_e,:), E_nf_mean, D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
      call compute_efactor(E_n, E_f(f_e,:), E_nf_mean, time, e_factor, .true.)
      call compute_F_ifn_omp(e_factor, D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
    end do ! f_e
    
  end subroutine compute_F_if_omp

  ! with matrix solution
  subroutine compute_F_if_omp_matrix(E_n_inp, E_n, E_f_inp, E_n_mean, E_f_mean, D_fn, D_ni, &
       time_A, time, F_if_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    use m_crank_nicolson, only: solve_crank_nicolson_matrix_z
    use m_crank_nicolson, only: solve_crank_nicolson_matrix_z_spline

    real(wp), intent(in):: E_n_inp(:)
    real(wp), intent(in):: E_n(:)
    real(wp),  intent(in):: E_f_inp(:,:)
    real(wp), intent(in):: E_n_mean
    real(wp), intent(in):: E_f_mean
    real(wp), intent(in):: D_fn(:,:,:) 
    real(wp), intent(in):: D_ni(:)
    real(wp), intent(in):: time_A(:)
    real(wp), intent(in):: time(:) 
    complex(wp), intent(out) ::  F_if_omp(:,:,:,:)
    real(wp), intent(in):: gamma
    
    complex(wp), allocatable::  e_factor(:)
    complex(wp), allocatable::  A(:,:,:)
    complex(wp), allocatable::  y(:,:,:), y_0(:,:)
    complex(wp), allocatable::  D_fn_new(:,:,:)
    integer:: nfinal, f_e, ntsteps, ntsteps_A, m1

    ntsteps_A = size(time_A)
    ntsteps = size(time)     
    
    allocate(e_factor(ntsteps))

    nfinal = size(E_f_inp,1) 

    allocate(A(nfinal,nfinal,ntsteps_A))
    allocate(y(nfinal,nfinal,ntsteps))
    allocate(D_fn_new(nfinal,ntsteps,3))
    allocate(y_0(nfinal,nfinal))
        
    ! set up diagonal A matrix for final states only
    y_0 = 0.0_wp
    A = 0.0_wp
    do f_e=1, nfinal
      y_0(f_e,f_e) = 1.0_wp
      !A(f_e,f_e,:) = (0.0_wp, -1.0_wp) * (const % eV  / const % hbar) *((E_n(:) - E_f(f_e,:)) &
      !     - (E_n_mean - E_f_mean))
      A(f_e,f_e,:) = (0.0_wp, -1.0_wp) * (const % eV  / const % hbar) *( E_f_inp(f_e,:) - E_f_mean)
    end do

    ! put some stuff into the martix

    !do f_e=1, nfinal
    !  A(f_e,2,:) = A(f_e,2,:) + (0.0_wp, -1.0_wp)*(const % eV  / const % hbar) * 5.e-2*f_e
    !  A(2,f_e,:) = A(2,f_e,:) + (0.0_wp, 1.0_wp)*(const % eV  / const % hbar) * 5.e-2*f_e
    !
    !  A(f_e,4,:) = A(f_e,4,:) + (const % eV  / const % hbar) * 0.95e-2*f_e
    !  A(4,f_e,:) = A(4,f_e,:) + (const % eV  / const % hbar) * 0.95e-2*f_e
    !end do
    
    !write(6,*) "solving matrix"
    !call solve_crank_nicolson_matrix_z(y_0,  time, A, y)

    write(6,*) "solving matrix with splining"
    call solve_crank_nicolson_matrix_z_spline(y_0,  time, time_A, A, y)

    ! compute efactor for E_n
    call compute_efactor_one(E_n, E_n_mean, time, e_factor, .true.)
    
    ! compute \sum_{f'} U^{\dagger}_{ff'}(t) D_{f'n}(t) = \sum_{f'} U^*_{f'f}(t) D_{f'n}(t)
    do f_e =1,nfinal
      do m1=1,3
        D_fn_new(f_e,:,m1) = sum( D_fn(:,:,m1) * conjg(y(:,f_e,:)),1) * e_factor(:)
      end do
    end do
    
    do f_e =1,nfinal! final state 
      !call compute_F_ifn_omp(E_n, E_f(f_e,:), E_nf_mean, D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
      !call compute_efactor(E_n, E_f(f_e,:), E_nf_mean, time, e_factor, .true.)
      !e_factor = y(f_e,f_e, :)
      !call compute_F_ifn_omp(e_factor*conjg(y(f_e,f_e, :)), D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
      call compute_F_ifn_omp_new(D_fn_new(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
    end do ! f_e
    
  end subroutine compute_F_if_omp_matrix

!  ! with matrix solution
!  subroutine compute_F_if_omp_matrix(E_n, E_f, E_nf_mean, D_fn, D_ni, &
!       time_A, time, F_if_omp, gamma)
!    use m_precision, only: wp
!    use m_constants, only: const
!    use m_fftw3, only: fft_c2c_1d_backward
!    use m_fftw3, only: reorder_sigma_fftw_z
!    use m_crank_nicolson, only: solve_crank_nicolson_matrix_z
!    use m_crank_nicolson, only: solve_crank_nicolson_matrix_z_spline
!
!    real(wp), intent(in):: E_n(:)
!    real(wp),  intent(in):: E_f(:,:)
!    real(wp), intent(in):: E_nf_mean
!    real(wp), intent(in):: D_fn(:,:,:) 
!    real(wp), intent(in):: D_ni(:)
!    real(wp), intent(in):: time_A(:)
!    real(wp), intent(in):: time(:) 
!    complex(wp), intent(out) ::  F_if_omp(:,:,:,:)
!    real(wp), intent(in):: gamma
!    
!    complex(wp), allocatable::  e_factor(:)
!    complex(wp), allocatable::  A(:,:,:)
!    complex(wp), allocatable::  y(:,:,:), y_0(:,:)
!    integer:: nfinal, f_e, ntsteps, ntsteps_A
!
!    ntsteps_A = size(time_A)
!    ntsteps = size(time)     
!    
!    allocate(e_factor(ntsteps))
!
!    nfinal = size(E_f,1) 
!
!    allocate(A(nfinal,nfinal,ntsteps_A))
!    allocate(y(nfinal,nfinal,ntsteps))
!    allocate(y_0(nfinal,nfinal))
!        
!    ! set up diagonal A matrix for final states only
!    y_0 = 0.0_wp
!    A = 0.0_wp
!    do f_e=1, nfinal
!      y_0(f_e,f_e) = 1.0_wp
!      !A(f_e,f_e,:) = (0.0_wp, -1.0_wp) * (const % eV  / const % hbar) *((E_n(:) - E_f(f_e,:))  - E_nf_mean)
!      A(f_e,f_e,:) = (0.0_wp, -1.0_wp) * (const % eV  / const % hbar) *((E_n(:) - E_f(f_e,:))  - E_nf_mean)
!    end do
!
!    !call compute_efactor(E_n, E_f(f_e,:), E_nf_mean, time, e_factor, .true.)
!    
!    !write(6,*) "solving matrix"
!    !call solve_crank_nicolson_matrix_z(y_0,  time, A, y)
!
!    write(6,*) "solving matrix with splining"
!    call solve_crank_nicolson_matrix_z_spline(y_0,  time, time_A, A, y)
!    
!    do f_e =1,nfinal! final state 
!      !call compute_F_ifn_omp(E_n, E_f(f_e,:), E_nf_mean, D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
!      !call compute_efactor(E_n, E_f(f_e,:), E_nf_mean, time, e_factor, .true.)
!      e_factor = y(f_e,f_e, :)
!      call compute_F_ifn_omp(e_factor, D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
!    end do ! f_e
!    
!  end subroutine compute_F_if_omp_matrix



  
  ! for a single combination of i,f,n, now including the D_ni matrix
  subroutine compute_F_ifn_omp(e_factor, D_fn, D_ni, time, F_ifn_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    !real(wp), intent(in):: E_n(:)
    !real(wp), intent(in):: E_f(:)
    complex(wp), intent(in)::  e_factor(:)
    real(wp), intent(in):: D_fn(:,:)
    real(wp), intent(in):: D_ni(:) 
    real(wp), intent(in):: time(:)
    complex(wp), intent(out) ::  F_ifn_omp(:,:,:)
    !real(wp), intent(in):: E_nf_mean, gamma
    real(wp), intent(in):: gamma

    integer:: ntsteps, nfinal
    complex(wp), allocatable:: funct(:)
    complex(wp), allocatable:: funct_fft(:)
    !complex(wp), allocatable::  e_factor1(:)
    integer:: m1, m2 

    ntsteps = size(time) 
    
    allocate(funct(ntsteps), &
         funct_fft(ntsteps))
        
    !call compute_efactor(E_n, E_f, E_nf_mean, time, e_factor1, .true.)
      
    F_ifn_omp = 0.0_wp
    
    do m1=1,3 ! polarization
      
      funct = D_fn(:,m1) * e_factor(:) * exp(-gamma * const % eV * time(:) / const % hbar)
        
      call fft_c2c_1d_backward(funct, funct_fft)
      call reorder_sigma_fftw_z(funct_fft)

      do m2=1,3 ! polarization     
        F_ifn_omp(:,m1,m2) = D_ni(m2) * funct_fft(:)
      end do ! m2 
        
    end do ! m1

      
  end subroutine compute_F_ifn_omp


  ! for a single combination of i,f,n, now including the D_ni matrix and complex D_fn matrix
  subroutine compute_F_ifn_omp_new(D_fn, D_ni, time, F_ifn_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    !real(wp), intent(in):: E_n(:)
    !real(wp), intent(in):: E_f(:)
    !complex(wp), intent(in)::  e_factor(:)
    complex(wp), intent(in):: D_fn(:,:)
    real(wp), intent(in):: D_ni(:) 
    real(wp), intent(in):: time(:)
    complex(wp), intent(out) ::  F_ifn_omp(:,:,:)
    !real(wp), intent(in):: E_nf_mean, gamma
    real(wp), intent(in):: gamma

    integer:: ntsteps, nfinal
    complex(wp), allocatable:: funct(:)
    complex(wp), allocatable:: funct_fft(:)
    !complex(wp), allocatable::  e_factor1(:)
    integer:: m1, m2 

    ntsteps = size(time) 
    
    allocate(funct(ntsteps), &
         funct_fft(ntsteps))
        
    F_ifn_omp = 0.0_wp
    
    do m1=1,3 ! polarization
      
      funct = D_fn(:,m1) * exp(-gamma * const % eV * time(:) / const % hbar)
        
      call fft_c2c_1d_backward(funct, funct_fft)
      call reorder_sigma_fftw_z(funct_fft)

      do m2=1,3 ! polarization     
        F_ifn_omp(:,m1,m2) = D_ni(m2) * funct_fft(:)
      end do ! m2 
        
    end do ! m1

      
  end subroutine compute_F_ifn_omp_new



  
  ! several intermediate states
  subroutine compute_F_if_omp_many_n(E_n, E_f, E_nf_mean, D_fn, D_ni, time, F_if_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_n(:,:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_f(:,:)
    real(wp), intent(in):: D_fn(:,:,:) 
    real(wp), intent(in):: D_ni(:,:) 
    complex(wp), dimension(:,:,:,:),intent(out) ::  F_if_omp
    real(wp), intent(in):: E_nf_mean, gamma
    
    integer:: nfinal, f_e, n_e, ntsteps
    complex(wp), allocatable:: F_tmp(:,:,:)

    nfinal = size(E_f,1)

    F_if_omp = 0.0_wp

    ! version with internal sum over n
    do f_e =1,nfinal
      call compute_F_if_omp_sum_n(E_n, E_f(f_e,:), E_nf_mean, D_fn(f_e,:,:), D_ni, time, F_if_omp(f_e,:,:,:), gamma)
    end do 

    ! other version (not implemented yet), compute each n separately and add
    
  end subroutine compute_F_if_omp_many_n




  !  ! for a single combination of i,f,n, now including the D_ni matrix
!  subroutine compute_F_ifn_omp(E_n, E_f, E_nf_mean, D_fn, D_ni, time, F_ifn_omp, gamma)
!    use m_precision, only: wp
!    use m_constants, only: const
!    use m_fftw3, only: fft_c2c_1d_backward
!    use m_fftw3, only: reorder_sigma_fftw_z
!    
!    real(wp), intent(in):: E_n(:)
!    real(wp), intent(in):: time(:)
!    real(wp), intent(in):: E_f(:)
!    real(wp), intent(in):: D_fn(:,:)
!    real(wp), intent(in):: D_ni(:) 
!    complex(wp), intent(out) ::  F_ifn_omp(:,:,:)
!    real(wp), intent(in):: E_nf_mean, gamma
!
!    integer:: ntsteps, nfinal
!    complex(wp), allocatable:: funct(:)
!    complex(wp), allocatable:: funct_fft(:)
!    complex(wp), allocatable::  e_factor1(:)
!    integer:: m1, m2 
!
!    ntsteps = size(time) 
!    
!    allocate(funct(ntsteps), &
!         funct_fft(ntsteps), &
!         e_factor1(ntsteps))
!    
!    call compute_efactor(E_n, E_f, E_nf_mean, time, e_factor1, .true.)
!      
!    F_ifn_omp = 0.0_wp
!    
!    do m1=1,3 ! polarization
!      
!      funct = D_fn(:,m1) * e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar)
!        
!      call fft_c2c_1d_backward(funct, funct_fft)
!      call reorder_sigma_fftw_z(funct_fft)
!
!      do m2=1,3 ! polarization     
!        F_ifn_omp(:,m1,m2) = D_ni(m2) * funct_fft(:)
!      end do ! m2 
!        
!    end do ! m1
!
!      
!  end subroutine compute_F_ifn_omp

  
  ! summing internally over several intermediate states to avoid fft:s, only savings for ninter > 3
  subroutine compute_F_if_omp_sum_n(E_n, E_f, E_nf_mean, D_fn, D_ni, time, F_if_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_n(:,:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_f(:)
    real(wp), intent(in):: D_fn(:,:)
    real(wp), intent(in):: D_ni(:,:) 
    complex(wp), intent(out) ::  F_if_omp(:,:,:)
    real(wp), intent(in):: E_nf_mean, gamma

    integer:: ntsteps, ninter
    complex(wp), allocatable:: funct(:,:,:)
    complex(wp), allocatable::  e_factor1(:)
    integer:: m1, m2, n_e 

    ntsteps = size(time) 
    ninter = size(E_n, 1)
    
    allocate(funct(ntsteps,3,3), &
         e_factor1(ntsteps))

    funct = 0.0_wp
    
    do n_e =1, ninter
      call compute_efactor(E_n(n_e,:), E_f, E_nf_mean, time, e_factor1, .true.)
      
      do m1=1,3 
        do m2=1,3 
          funct(:,m1, m2) = funct(:,m1,m2) + D_ni (n_e, m2) * D_fn(:, m1) * &
               e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar)
        end do
      end do
    end do
    
    F_if_omp = 0.0_wp
    
    do m1=1,3 
      do m2= 1,3 
        
        call fft_c2c_1d_backward(funct(:,m1, m2), F_if_omp(:,m1,m2))
        call reorder_sigma_fftw_z(F_if_omp(:,m1,m2))
        
      end do ! m2 
    end do ! m1
    
  end subroutine compute_F_if_omp_sum_n

  subroutine compute_F_if_omp_one_n(E_n, E_f, E_nf_mean, D_fn, D_ni, time, F_if_omp, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_n(:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_f(:)
    real(wp), intent(in):: D_fn(:,:)
    real(wp), intent(in):: D_ni(:) 
    complex(wp), intent(out) ::  F_if_omp(:,:,:)
    real(wp), intent(in):: E_nf_mean, gamma

    integer:: ntsteps, ninter
    complex(wp), allocatable:: funct(:,:,:)
    complex(wp), allocatable::  e_factor1(:)
    integer:: m1, m2, n_e 

    ntsteps = size(time) 
    ninter = size(E_n, 1)
    
    allocate(funct(ntsteps,3,3), &
         e_factor1(ntsteps))

    funct = 0.0_wp
    
    call compute_efactor(E_n(:), E_f(:), E_nf_mean, time(:), e_factor1(:), .true.)
    
    do m1=1,3 
      do m2=1,3 
        funct(:,m1, m2) = funct(:,m1,m2) + D_ni (m2) * D_fn(:, m1) * &
             e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar)
      end do
    end do
    
    F_if_omp = 0.0_wp
    
    do m1=1,3 
      do m2= 1,3 
        
        call fft_c2c_1d_backward(funct(:,m1, m2), F_if_omp(:,m1,m2))
        call reorder_sigma_fftw_z(F_if_omp(:,m1,m2))
        
      end do ! m2 
    end do ! m1
    
  end subroutine compute_F_if_omp_one_n


  
!  ! summing internally over several intermediate states to avoid fft:s, only savings for ninter > 3
!  subroutine compute_F_ifc_omp_sum_n(E_n, E_fc, E_nf_mean, D_fn, D_ni, time, F_if_omp, gamma)
!    use m_precision, only: wp
!    use m_constants, only: const
!    use m_fftw3, only: fft_c2c_1d_backward
!    use m_fftw3, only: reorder_sigma_fftw_z
!    
!    real(wp), intent(in):: E_n(:,:)
!    real(wp), intent(in):: time(:)
!    real(wp), intent(in):: E_fc(:,:)
!    real(wp), intent(in):: D_fn(:,:,:)
!    real(wp), intent(in):: D_ni(:,:) 
!    complex(wp), intent(out) ::  F_if_omp(:,:,:)
!    real(wp), intent(in):: E_nf_mean, gamma
!
!    integer:: ntsteps, ninter
!    complex(wp), allocatable:: funct(:,:,:)
!    complex(wp), allocatable::  e_factor1(:)
!    integer:: m1, m2, n_e 
!
!    ntsteps = size(time) 
!    ninter = size(E_n, 1)
!    
!    allocate(funct(ntsteps,3,3), &
!         e_factor1(ntsteps))
!
!    funct = 0.0_wp
!    
!    do n_e =1, ninter
!      call compute_efactor(E_n(n_e,:), E_fc(n_e,:), E_nf_mean, time, e_factor1, .true.)
!      
!      do m1=1,3 
!        do m2=1,3 
!          funct(:,m1, m2) = funct(:,m1,m2) + D_ni (n_e, m2) * D_fn(n_e, :, m1) * &
!               e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar)
!        end do
!      end do
!    end do
!    
!    F_if_omp = 0.0_wp
!    
!    do m1=1,3 
!      do m2= 1,3 
!        
!        call fft_c2c_1d_backward(funct(:,m1, m2), F_if_omp(:,m1,m2))
!        call reorder_sigma_fftw_z(F_if_omp(:,m1,m2))
!        
!      end do ! m2 
!    end do ! m1
!    
!  end subroutine compute_F_ifc_omp_sum_n


  

  subroutine compute_F_if_om_omp(F_if_t_omp, E_f, E_fi_mean, time, &
       E_i, gamma_inc, omega_out, E_nf_mean,F_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_omp(:,:,:,:,:)
    real(wp), intent(in):: E_f(:,:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: gamma_inc
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: E_nf_mean
    complex(wp), intent(out) ::  F_if_om_omp(:,:,:,:,:)

    integer:: nfinal, n_omega_in, n_omega_out, f_e, om_out, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    
    nfinal = size(F_if_t_omp,1)
    n_omega_in = size(F_if_t_omp,2)
    n_omega_out = size(F_if_t_omp,3)
    
    allocate(e_factor1(n_omega_in))
    
    do f_e= 1, nfinal
      
      call compute_efactor(E_f(f_e,:), E_i, E_fi_mean, time, e_factor1(:), .false.)
      
      do om_out= 1, n_omega_out
        do m1 =1, 3
          do m2 =1, 3

            ! no broadening here
            call fft_c2c_1d_forward( e_factor1(:) * &
                 !exp(-gamma_inc * const % eV * time(:) / const % hbar) * &
                 !exp(-(gamma_inc * const % eV) ** 2 / (4.0_wp * log(2.0_wp)) * (time(:) / const % hbar)**2) * &
                 exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean)* const % eV * time(:) / const % hbar )) * &
                 F_if_t_omp(f_e, :,om_out, m1,m2), &
                 F_if_om_omp(f_e, :,om_out, m1,m2))
            call reorder_sigma_fftw_z(F_if_om_omp(f_e,:, om_out, m1,m2))
            
          end do
        end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp

  subroutine compute_F_if_om_omp_one_f(F_if_t_omp, E_f, E_fi_mean, time, &
       E_i, gamma_inc, omega_out, E_nf_mean,F_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_omp(:,:,:,:)
    real(wp), intent(in):: E_f(:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: gamma_inc
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: E_nf_mean
    complex(wp), intent(out) ::  F_if_om_omp(:,:,:,:)

    integer:: n_omega_in, n_omega_out, om_out, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    
    n_omega_in = size(F_if_t_omp,1)
    n_omega_out = size(F_if_t_omp,2)
    
    allocate(e_factor1(n_omega_in))
      
    call compute_efactor(E_f(:), E_i, E_fi_mean, time, e_factor1(:), .false.)
      
    do om_out= 1, n_omega_out
      do m1 =1, 3
        do m2 =1, 3
          
          call fft_c2c_1d_forward( e_factor1(:) * &
               exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean)* const % eV * time(:) / const % hbar )) * &
               F_if_t_omp(:,om_out, m1,m2), &
               F_if_om_omp(:,om_out, m1,m2))
          call reorder_sigma_fftw_z(F_if_om_omp(:, om_out, m1,m2))
          
        end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp_one_f

  subroutine compute_F_if_om_omp_one_f_separate(F_if_t_omp, E_f1, E_f2, E_fi_mean, time, &
       E_i1, E_n1, E_n2,  gamma_inc, omega_out, E_nf_mean,F_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_omp(:,:,:,:)
    real(wp), intent(in):: E_f1(:)
    real(wp), intent(in):: E_f2(:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i1(:)
    real(wp), intent(in):: E_n1(:)
    real(wp), intent(in):: E_n2(:)
    real(wp), intent(in):: gamma_inc
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: E_nf_mean
    complex(wp), intent(out) ::  F_if_om_omp(:,:,:,:)

    integer:: n_omega_in, n_omega_out, om_out, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    real(wp):: mean
    
    n_omega_in = size(F_if_t_omp,1)
    n_omega_out = size(F_if_t_omp,2)
    
    allocate(e_factor1(n_omega_in))

    mean = sum(E_n2-E_n1)/size(E_n1) - sum(E_f2-E_f1)/size(E_n1)
    !mean =  sum(E_f2-E_f1)/size(E_n1)
    
    write(6,*) "mean", mean !* const % eV

    call compute_efactor(E_f2(:) -E_n2(:)+mean, E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f2(:) -E_n2(:), E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f2(:) -E_n2(:)+mean, E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .true.)
    !call compute_efactor(E_f2(:), E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f2(:)-mean, E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f1(:), E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f1(:), E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    
    do om_out= 1, n_omega_out
      do m1 =1, 3
        do m2 =1, 3
          
          call fft_c2c_1d_forward( e_factor1(:) * &
               exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean)* const % eV * time(:) / const % hbar )) * &
               F_if_t_omp(:,om_out, m1,m2), &
               F_if_om_omp(:,om_out, m1,m2))
          call reorder_sigma_fftw_z(F_if_om_omp(:, om_out, m1,m2))
          
        end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp_one_f_separate
  
  subroutine compute_F_if_om_omp_one_f_factor_XAS(F_if_t_omp, E_f1, E_f2, E_fi_mean, time, &
       E_i1, E_n1, E_n2,  gamma_inc, omega_out, omega_in, E_nf_mean,F_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_omp(:,:,:,:)
    real(wp), intent(in):: E_f1(:)
    real(wp), intent(in):: E_f2(:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i1(:)
    real(wp), intent(in):: E_n1(:)
    real(wp), intent(in):: E_n2(:)
    real(wp), intent(in):: gamma_inc
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: E_nf_mean
    complex(wp), intent(out) ::  F_if_om_omp(:,:,:,:)

    integer:: n_omega_in, n_omega_out, om_in, om_out, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    real(wp):: mean
    
    n_omega_in = size(F_if_t_omp,1)
    n_omega_out = size(F_if_t_omp,2)
    
    allocate(e_factor1(n_omega_in))

    mean = sum(E_n2-E_n1)/size(E_n1) - sum(E_f2-E_f1)/size(E_n1)
    
    write(6,*) "mean", mean !* const % eV

    call compute_efactor(E_f2(:) -E_n2(:)+mean, E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f2(:) -E_n2(:), E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f2(:) -E_n2(:)+mean, E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .true.)
    !call compute_efactor(E_f2(:), E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_n1(:), E_i1(:), E_nf_mean-E_fi_mean, time, e_factor1(:), .false.)

    
    do om_out= 1, n_omega_out
      do m1 =1, 3
        do m2 =1, 3

          call fft_c2c_1d_forward( e_factor1(:) * &
               exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean)* const % eV * time(:) / const % hbar )) * &
               F_if_t_omp(:,om_out, m1,m2), &
               F_if_om_omp(:,om_out, m1,m2))
          call reorder_sigma_fftw_z(F_if_om_omp(:, om_out, m1,m2))

         end do
      end do
    end do

    ! cut spectrum at absorption edge
    do om_in= 1, n_omega_in
      if (omega_in(om_in) .lt. 200) then
        F_if_om_omp(om_in,:,:,:) = 0.0d0
      end if
    end do
    
  end subroutine compute_F_if_om_omp_one_f_factor_XAS

  

  subroutine compute_F_if_om_omp_one_f_separate_ingoing(F_if_t_omp, E_f1, E_f2, E_fi_mean, time, &
       E_i1, E_i2, E_n1, E_n2,  gamma_inc, omega_in, E_nf_mean, E_ni_mean, F_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_omp(:,:,:,:)
    real(wp), intent(in):: E_f1(:)
    real(wp), intent(in):: E_f2(:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i1(:)
    real(wp), intent(in):: E_i2(:)
    real(wp), intent(in):: E_n1(:)
    real(wp), intent(in):: E_n2(:)
    real(wp), intent(in):: gamma_inc
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: E_nf_mean
    real(wp), intent(in):: E_ni_mean
    complex(wp), intent(out) ::  F_if_om_omp(:,:,:,:)

    integer:: n_omega_in, om_in, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    real(wp):: mean
    
    n_omega_in = size(F_if_t_omp,2)
    !n_omega_out = size(F_if_t_omp,2)
    
    allocate(e_factor1(n_omega_in))

    mean = sum(E_n2-E_n1)/size(E_n1) - sum(E_f2-E_f1)/size(E_n1)
    
    
    write(6,*) "mean", mean !* const % eV

    !call compute_efactor(E_f2(:) -E_n2(:)+mean, E_i1(:)-E_n1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f1(:) -E_n1(:), E_i2(:)-E_n2(:), E_fi_mean, time, e_factor1(:), .false.)
    call compute_efactor(E_f2(:), E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f2(:), E_i1(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f1(:), E_i2(:), E_fi_mean, time, e_factor1(:), .false.)
    !call compute_efactor(E_f1(:), E_i1(:), E_fi_mean, time, e_factor1(:), .true.)
      
    do om_in= 1, n_omega_in
      do m1 =1, 3
        do m2 =1, 3
          
          call fft_c2c_1d_backward( e_factor1(:) * &
               exp(dcmplx(0.0_wp,  (-omega_in(om_in) + E_ni_mean)* const % eV * time(:) / const % hbar )) * &
               F_if_t_omp(:,om_in, m1,m2), &
               F_if_om_omp(om_in,:, m2,m1)) ! just switch om_in and om_out and m1 and m2
          call reorder_sigma_fftw_z(F_if_om_omp(om_in, :, m2,m1))
          
        end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp_one_f_separate_ingoing

  

  

  
  ! here use factorization
  subroutine compute_F_if_om_omp_no_F(E_f, E_fi_mean, time, &
       E_i, gamma_R, omega_out, E_nf_mean, R_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    !complex(wp), intent(in) ::  F_if_t_omp(:,:,:,:,:)
    real(wp), intent(in):: E_f(:,:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: gamma_R
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: E_nf_mean
    complex(wp), intent(out) ::  R_if_om_omp(:,:,:)

    integer:: nfinal, n_omega_in, n_omega_out, f_e, om_out
    complex(wp), allocatable ::  e_factor1(:)
    
    nfinal = size(E_f,1)
    n_omega_in = size(E_f,2)
    n_omega_out = size(omega_out)
    
    allocate(e_factor1(n_omega_in))
    
    do f_e= 1, nfinal
      
      call compute_efactor(E_f(f_e,:), E_i, E_fi_mean, time, e_factor1(:), .false.)
      
      do om_out= 1, n_omega_out
        !do m1 =1, 3
        !  do m2 =1, 3
            
        call fft_c2c_1d_forward( e_factor1(:) * &
             !exp(-gamma_R * const % eV * time(:) / const % hbar) * &
             !exp(-(gamma_inc * const % eV) ** 2 / (4.0_wp * log(2.0_wp)) * (time(:) / const % hbar)**2) * &
             exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean)* const % eV * time(:) / const % hbar )), &
             !* &
             !F_if_t_omp(f_e, :,om_out, m1,m2), &
             R_if_om_omp(f_e, :,om_out))
        call reorder_sigma_fftw_z(R_if_om_omp(f_e,:, om_out))
        
        ! end do
        !end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp_no_F


  ! here use factorization
  subroutine compute_F_if_om_omp_no_F_one(E_f, E_fi_mean, time, &
       E_i, gamma_R, omega_out, E_nf_mean, R_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_f(:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: gamma_R
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: E_nf_mean
    complex(wp), intent(out) ::  R_if_om_omp(:,:)

    integer:: nfinal, n_omega_in, n_omega_out, f_e, om_out
    complex(wp), allocatable ::  e_factor1(:)
    
    n_omega_in = size(E_f,1)
    n_omega_out = size(omega_out)
    
    allocate(e_factor1(n_omega_in))
    
    call compute_efactor(E_f(:), E_i, E_fi_mean, time, e_factor1(:), .false.)
    
    do om_out= 1, n_omega_out
      
      call fft_c2c_1d_forward( e_factor1(:) * &
           exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean)* const % eV * time(:) / const % hbar )), &
           R_if_om_omp(:,om_out))
      call reorder_sigma_fftw_z(R_if_om_omp(:, om_out))
      
      
    end do
    
      
  end subroutine compute_F_if_om_omp_no_F_one

  

  
  !
  ! using F(omega) instead of F(omega') 
  !
  
  ! several intermediate states
  subroutine compute_F_if_om_many_n(E_n, E_i, E_ni_mean, D_fn, D_ni, time, F_if_om, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_n(:,:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: D_fn(:,:,:) 
    real(wp), intent(in):: D_ni(:,:,:) 
    complex(wp), dimension(:,:,:,:),intent(out) ::  F_if_om
    real(wp), intent(in):: E_ni_mean, gamma
    
    integer:: nfinal, f_e, n_e, ntsteps
    complex(wp), allocatable:: F_tmp(:,:,:)

    nfinal = size(D_fn,1)

    F_if_om = 0.0_wp

    ! version with internal sum over n
    do f_e =1,nfinal
      !call compute_F_if_om_sum_n(E_n, E_f(f_e,:), E_nf_mean, D_fn(f_e,:,:,:), D_ni, time, F_if_om(f_e,:,:,:), gamma)
      call compute_F_if_om_sum_n(E_n(:,:), E_i(:), E_ni_mean, D_fn(f_e,:,:), D_ni(:,:,:), time, F_if_om(f_e,:,:,:), gamma)
    end do 

    ! other version (not implemented yet), compute each n separately and add
    
  end subroutine compute_F_if_om_many_n

  
  ! summing internally over several intermediate states to avoid fft:s, only savings for ninter > 3
  subroutine compute_F_if_om_sum_n(E_n, E_i, E_ni_mean, D_fn, D_ni, time, F_if_om, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_n(:,:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: D_fn(:,:)
    real(wp), intent(in):: D_ni(:,:,:) 
    complex(wp), intent(out) ::  F_if_om(:,:,:)
    real(wp), intent(in):: E_ni_mean, gamma

    integer:: ntsteps, ninter
    complex(wp), allocatable:: funct(:,:,:)
    complex(wp), allocatable::  e_factor1(:)
    integer:: m1, m2, n_e 

    ntsteps = size(time) 
    ninter = size(E_n, 1)
    
    allocate(funct(ntsteps,3,3), &
         e_factor1(ntsteps))

    funct = 0.0_wp
    
    do n_e =1, ninter
      call compute_efactor(E_n(ninter,:), E_i(:), E_ni_mean, time, e_factor1, .true.) ! .false. ? back in time
      
      do m1=1,3 
        do m2=1,3 
          funct(:,m1, m2) = funct(:,m1,m2) + D_fn(n_e, m1) * D_ni (n_e, :, m2) * &
               e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar)
        end do
      end do
    end do
    
    F_if_om = 0.0_wp
    
    do m1=1,3 
      do m2= 1,3 
        
        call fft_c2c_1d_backward(funct(:,m1, m2), F_if_om(:,m1,m2))
        call reorder_sigma_fftw_z(F_if_om(:,m1,m2))
        
      end do ! m2 
    end do ! m1
    
  end subroutine compute_F_if_om_sum_n
  
  subroutine compute_F_if_om_omp_ingoing(F_if_t_om, E_f, E_fi_mean, time, &
       E_i, gamma_instr, omega_in, E_ni_mean, F_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_om(:,:,:,:,:)
    real(wp), intent(in):: E_f(:,:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: gamma_instr
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: E_ni_mean
    complex(wp), intent(out) ::  F_if_om_omp(:,:,:,:,:)

    integer:: nfinal, n_omega_in, n_omega_out, f_e, om_in, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    
    nfinal = size(F_if_t_om,1)
    n_omega_in = size(F_if_t_om,2)
    n_omega_out = size(F_if_t_om,3)
    
    allocate(e_factor1(n_omega_out))
    
    do f_e= 1, nfinal
      
      call compute_efactor(E_f(f_e,:), E_i, E_fi_mean, time, e_factor1(:), .true.) ! before: false
      
      do om_in= 1, n_omega_in
        do m1 =1, 3
          do m2 =1, 3
            
            call fft_c2c_1d_forward( e_factor1(:) * &
                 exp(-gamma_instr * const % eV * time(:) / const % hbar) * &
                 exp(dcmplx(0.0_wp,  (omega_in(om_in) - E_ni_mean)* const % eV * time(:) / const % hbar )) * &
                 F_if_t_om(f_e, :, om_in, m1,m2), &
                 F_if_om_omp(f_e, om_in,:, m1,m2)) ! switch places of indices
            call reorder_sigma_fftw_z(F_if_om_omp(f_e, om_in,:, m1,m2))
            
          end do
        end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp_ingoing

  subroutine compute_F_if_om_omp_ingoing_no_F(E_f, E_fi_mean, time, &
       E_i, gamma_instr, omega_in, E_ni_mean, R_if_om_omp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    !complex(wp), intent(in) ::  F_if_t_om(:,:,:,:,:)
    real(wp), intent(in):: E_f(:,:)
    real(wp), intent(in):: E_fi_mean
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: gamma_instr
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: E_ni_mean
    complex(wp), intent(out) ::  R_if_om_omp(:,:,:)

    integer:: nfinal, n_omega_in, n_omega_out, f_e, om_in, m1, m2
    complex(wp), allocatable ::  e_factor1(:)
    
    !nfinal = size(F_if_t_om,1)
    !n_omega_in = size(F_if_t_om,2)
    !n_omega_out = size(F_if_t_om,3)
    nfinal = size(E_f,1)
    n_omega_in = size(omega_in)
    n_omega_out = size(E_f,2)
    
    allocate(e_factor1(n_omega_out))
    
    do f_e= 1, nfinal
      
      call compute_efactor(E_f(f_e,:), E_i, E_fi_mean, time, e_factor1(:), .true.) ! before: false
      
      do om_in= 1, n_omega_in
        !do m1 =1, 3
        !  do m2 =1, 3
        
        call fft_c2c_1d_forward( e_factor1(:) * &
             !exp(-gamma_instr * const % eV * time(:) / const % hbar) * &
             exp(dcmplx(0.0_wp,  (omega_in(om_in) - E_ni_mean)* const % eV * time(:) / const % hbar )), &
             !F_if_t_om(f_e, :, om_in, m1,m2), &
             R_if_om_omp(f_e, om_in,:)) ! switch places of indices
        call reorder_sigma_fftw_z(R_if_om_omp(f_e, om_in,:))
        
        !          end do
        !        end do
      end do
    end do
    
  end subroutine compute_F_if_om_omp_ingoing_no_F

  !
  ! New SCXAS routines, as similar as possible to SCKH
  !

  subroutine compute_SCXAS_sigma_f_om_mm(E_i, E_f, E_fi_mean, D_fi, time, sigma_f_om_mm, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), dimension(:), intent(in):: E_i, time
    real(wp), dimension(:,:), intent(in):: E_f
    real(wp),dimension(:,:,:),intent(in):: D_fi 
    real(wp), dimension(:,:,:,:),intent(out) ::  sigma_f_om_mm
    real(wp), intent(in):: E_fi_mean, gamma
    
    integer:: nfinal, ntsteps, f_e
    complex(wp), allocatable:: funct(:)
    complex(wp), allocatable:: funct_fft(:)
    complex(wp), allocatable:: F_m(:,:)
    complex(wp), allocatable::  e_factor1(:)
    integer:: m1, m2 
    
    nfinal = size(E_f,1) 
    ntsteps = size(time) 
    
    allocate(funct(ntsteps), &
         funct_fft(ntsteps), &
         e_factor1(ntsteps), &
         F_m(ntsteps,3))
    
    do f_e =1,nfinal! final state 

      call compute_efactor(E_f(f_e,:), E_i, E_fi_mean, time, e_factor1, .true.)
      
      do m1=1,3 ! polarization

        funct = D_fi(f_e, :,m1) * e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar)
        
        call fft_c2c_1d_backward(funct, funct_fft)
        call reorder_sigma_fftw_z(funct_fft)

        F_m(:,m1) = funct_fft(:)

      end do ! m1
      
      ! compute full sigma tensor
      do m1=1,3
        do m2=1,3
          sigma_f_om_mm(f_e,:,m1,m2) = real(conjg(F_m(:,m1))* F_m(:,m2))
        end do
      end do
      
    end do ! f_e
    
  end subroutine compute_SCXAS_sigma_f_om_mm

  ! XAS with FC, only vertical transition

  subroutine compute_SCXAS_sigma_FC_f_om_mm(E_i, E_f, D_fi, omega_in, sigma_f_om_mm, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_i
    real(wp), intent(in):: E_f(:)
    real(wp), intent(in):: D_fi(:,:) 
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(out) ::  sigma_f_om_mm(:,:,:,:)
    real(wp), intent(in):: gamma
    
    integer:: nfinal, n_omega_in, f_e, om_in
    complex(wp), allocatable:: F_m(:,:)
    integer:: m1, m2 
    integer:: i_low
    real(wp):: i_low_weight
    
    nfinal = size(E_f,1) 
    n_omega_in = size(omega_in)
    
    allocate( F_m(n_omega_in,3))
    
    do f_e =1,nfinal! final state 

      call put_on_grid(omega_in, E_f(f_e)-E_i, i_low, i_low_weight )

      F_m = 0.0_wp
      do om_in= 1, n_omega_in
        do m1=1,3 ! polarization
          
          F_m(om_in, m1) = F_m(om_in, m1) + D_fi(f_e,m1) * i_low_weight / dcmplx(omega_in(om_in) &
               - omega_in(i_low), gamma)
          F_m(om_in, m1) = F_m(om_in, m1) + D_fi(f_e,m1) * (1.0_wp -i_low_weight) / &
               dcmplx(omega_in(om_in) - omega_in(i_low+1), gamma)
          
        end do ! m1
      end do ! om_in
      
      ! compute full sigma tensor
      do m1=1,3
        do m2=1,3
          sigma_f_om_mm(f_e,:,m1,m2) = real(conjg(F_m(:,m1))* F_m(:,m2))
        end do
      end do
      
    end do ! f_e
    
  end subroutine compute_SCXAS_sigma_FC_f_om_mm


  ! given a equally spaced range x, and a point x_in
  ! find i_low, the index of the point below x_in
  ! and i_low_weigth in that point, assuming "tents"
  ! the weight on i_low +1 is 1-i_low_weigth
  subroutine put_on_grid(x, x_in, i_low, i_low_weight )
    
    use m_precision, only: wp
    
    real(wp), intent(in):: x(:)
    real(wp), intent(in):: x_in
    integer, intent(out):: i_low
    real(wp), intent(out):: i_low_weight
    
    integer:: nx
    real(wp):: dx
    
    nx = size(x)
    dx = x(2)-x(1)
    
    i_low = floor((x_in-x(1))/dx) +1
    
    if(i_low .lt. 1) then
      i_low = 1
      i_low_weight =1.0_wp
    else if(i_low .gt. nx-1) then
      i_low = nx-1
      i_low_weight =0.0_wp
    else
      i_low_weight = -x_in/dx +(1 + x(i_low)/dx)
    end if
    
    !   write(6,*) "x_in, i_low, i_low_weight, x(i_low)", x_in, i_low, i_low_weight, x(i_low)
   
   
  end subroutine put_on_grid

  
  

  !
  ! Old SCXAS routines
  !


  
  subroutine compute_SCXAS(E_i, E_f, E_fi_mean, D_fi, time, sigma_m, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_FFT, only: FFT_complex 

    ! passed variables
    real(wp), dimension(:), intent(in):: E_i, time
    real(wp), dimension(:,:), intent(in):: E_f
    real(wp),dimension(:,:,:),intent(in):: D_fi 
    complex(wp), dimension(:,:,:),intent(out) ::  sigma_m
    real(wp), intent(in):: E_fi_mean, gamma
    ! local variables
    integer:: ntsteps, ntsteps_pad, nfinal
    complex(wp), dimension(:),allocatable:: funct
    complex(wp), dimension(:,:),allocatable::  e_factor1
    real(wp), dimension(:,:),allocatable:: int_W_I
    real(wp),dimension(:),allocatable::  omega_out 
    real(wp):: delta_t
    integer:: i,k,m 

    ntsteps = size(time) 
    nfinal = size(E_f,1) 
    ntsteps_pad = size(sigma_m,2) 
    
    allocate(funct(ntsteps), int_W_I(nfinal, ntsteps), &
         e_factor1(nfinal,ntsteps), omega_out(ntsteps_pad) )
       
    delta_t = time(2)-time(1)
        
    int_W_I(:,1) = E_f(:,1) - E_i(1)  - E_fi_mean  
    
    do i = 2, ntsteps
       int_W_I(:,i) = int_W_I(:,i-1) + ( E_f(:,i) - E_i(i)) - E_fi_mean  
    end do
    int_W_I =  -int_W_I * delta_t
    
    do i = 1,nfinal 
       e_factor1(i,:) = exp(dcmplx(0, -(const % eV  / const % hbar) *int_W_I(i,:)  ))
    end do
    
    ! compute A_{fm}(t) = D^m_{fn}(0) * D^m_{fn}(t) * e_factor_f(t) * exp(-gamma * eV * time(:) / hbar)
    ! and fourer transform
    
    sigma_m = 0.0_wp
    !alpha=4.0_wp*log(2.0_wp)/((gamma*2)**2.0_wp)    

    !do a FFT
    do k =1,nfinal! final state 
       do m=1,3 ! polarization    
          
          funct = D_fi(k,1,m) * D_fi(k,:,m) * e_factor1(k,:) * exp(-gamma * const % eV * time(:) / const % hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !funct = (cos(100 * eV * time(:) / hbar  )  - cos(110 * eV * time(:) / hbar  )) * exp(-gamma * eV * time(:) / hbar) 

          call FFT_complex(time, funct, sigma_m(k,:,m), omega_out) 

       end do ! m
    end do ! k
    
    deallocate(funct, int_W_I)
    
  end subroutine compute_SCXAS

  subroutine compute_SCXAS_new(corr_fn, time, sigma_m, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_FFT, only: FFT_complex 

    ! passed variables
    real(wp), dimension(:), intent(in):: time
    complex(wp), dimension(:,:,:),intent(in) ::  corr_fn
    complex(wp), dimension(:,:,:),intent(out) ::  sigma_m
    real(wp), intent(in):: gamma
    ! local variables
    integer:: ntsteps, ntsteps_pad, nfinal
    complex(wp), dimension(:),allocatable:: funct
    real(wp),dimension(:),allocatable::  omega_out 
    integer:: k,m 

    ntsteps = size(time) 
    nfinal = size(corr_fn,1) 
    ntsteps_pad = size(sigma_m,2) 
    
    allocate(funct(ntsteps),  omega_out(ntsteps_pad) )
       
    ! compute A_{fm}(t) = D^m_{fn}(0) * D^m_{fn}(t) * e_factor_f(t) * exp(-gamma * eV * time(:) / hbar)
    ! and fourer transform
    
    sigma_m = 0.0_wp

    !do a FFT
    do k =1,nfinal! final state 
       do m=1,3 ! polarization    
          
          funct = corr_fn(k,:,m) *exp(-gamma * const % eV * time(:) / const % hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !funct = (cos(5 * eV * time(:) / hbar  )  - cos(11 * eV * time(:) / hbar  )) * exp(-gamma * eV * time(:) / hbar) 

          call FFT_complex(time, funct, sigma_m(k,:,m), omega_out) 
       end do ! m
    end do ! k
    
    deallocate(funct, omega_out)
    
  end subroutine compute_SCXAS_new

  subroutine collect_corr_fun(E_i, E_f, E_fi_mean, D_fi, time, corr_fn)
    use m_precision, only: wp
    use m_constants, only: const

    ! passed variables
    real(wp), dimension(:), intent(in):: E_i, time
    real(wp), dimension(:,:), intent(in):: E_f
    real(wp),dimension(:,:,:),intent(in):: D_fi 
    complex(wp), dimension(:,:,:),intent(inout) ::  corr_fn
    real(wp), intent(in):: E_fi_mean
    ! local variables
    integer:: ntsteps,  nfinal
    !complex(wp), dimension(:),allocatable:: funct
    complex(wp), dimension(:,:),allocatable::  e_factor1
    real(wp), dimension(:,:),allocatable:: int_W_I
    !real(wp),dimension(:),allocatable::  omega_out 
    real(wp):: delta_t
    integer:: i,k,m

    ntsteps = size(time) 
    nfinal = size(E_f,1) 
    !ntsteps_pad = size(sigma_m,2) 
    
    allocate( int_W_I(nfinal, ntsteps), &
         e_factor1(nfinal,ntsteps) )
       
    delta_t = time(2)-time(1)
        
    int_W_I(:,1) = E_f(:,1) - E_i(1)  - E_fi_mean  
    
    do i = 2, ntsteps
       int_W_I(:,i) = int_W_I(:,i-1) + ( E_f(:,i) - E_i(i)) - E_fi_mean  
       !int_W_I(:,i) =  ( E_f(:,1) - E_i(1)) - E_fi_mean  
    end do
    int_W_I =  int_W_I * delta_t
       

    do i = 1,nfinal 
       !e_factor1(i,:) = cos( -(eV  / hbar) *int_W_I(i,:) ) !exp(dcmplx(0, -(eV  / hbar) *int_W_I(i,:)  ))
       e_factor1(i,:) = exp(dcmplx(0, -(const % eV  / const % hbar) *int_W_I(i,:)  ))
    end do
    
    !do i = 1, ntsteps
    !   write(24,'(13ES25.16)') time(i), e_factor1(:,i)
    !end do

    do k =1,nfinal! final state 
       do m=1,3 ! polarization    

          !corr_fn(k,:,m) = corr_fn(k,:,m) +  e_factor1(k,:) ! exp(-gamma * eV * time(:) / hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !corr_fn(k,:,m) = corr_fn(k,:,m) + D_fi(k,1,m) * D_fi(k,:,m)   !* exp(-gamma * eV * time(:) / hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          corr_fn(k,:,m) = corr_fn(k,:,m) + D_fi(k,1,m) * D_fi(k,:,m) * e_factor1(k,:)  !* exp(-gamma * eV * time(:) / hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !funct = (cos(100 * eV * time(:) / hbar  )  - cos(110 * eV * time(:) / hbar  )) * exp(-gamma * eV * time(:) / hbar) 

          !call FFT_complex(time, funct, sigma_m(k,:,m), omega_out) 

       end do ! m
    end do ! k
    
    deallocate(int_W_I, e_factor1)
    
  end subroutine collect_corr_fun


  subroutine compute_corr_fun_FT(E_i, E_f, E_fi_mean, D_fi, time, gamma, corr_fn_FT_tmp)
    use m_precision, only: wp
    use m_constants, only: const
    use m_FFT, only: FFT_complex 

    ! passed variables
    real(wp), dimension(:), intent(in):: E_i, time
    real(wp), dimension(:,:), intent(in):: E_f
    real(wp),dimension(:,:,:),intent(in):: D_fi 
    real(wp),intent(in):: gamma
    complex(wp), dimension(:,:,:,:),intent(out) ::  corr_fn_FT_tmp
    real(wp), intent(in):: E_fi_mean
    ! local variables
    integer:: ntsteps,  nfinal, ntsteps_pad
    complex(wp), dimension(:,:),allocatable::  e_factor1
    real(wp), dimension(:,:),allocatable:: int_W_I
    complex(wp), allocatable::  funct(:), funct_out(:,:)    
    real(wp), allocatable:: omega_out(:)
    real(wp):: delta_t
    integer:: i,k,m, m1, m2

    ntsteps = size(time) 
    nfinal = size(E_f,1) 
    ntsteps_pad = size(corr_fn_FT_tmp,2) 

    allocate(funct_out(ntsteps, 3))    
    allocate(funct(ntsteps),  omega_out(ntsteps_pad) )    
    allocate( int_W_I(nfinal, ntsteps), &
         e_factor1(nfinal,ntsteps) )
       
    delta_t = time(2)-time(1)
        
    int_W_I(:,1) = E_f(:,1) - E_i(1)  - E_fi_mean  
    
    do i = 2, ntsteps
       int_W_I(:,i) = int_W_I(:,i-1) + ( E_f(:,i) - E_i(i)) - E_fi_mean  
       !int_W_I(:,i) =  ( E_f(:,1) - E_i(1)) - E_fi_mean  
    end do
    int_W_I =  int_W_I * delta_t
       
    do i = 1,nfinal 
       e_factor1(i,:) = exp(dcmplx(0, -(const % eV  / const % hbar) *int_W_I(i,:)  ))
    end do
    
    do k =1,nfinal! final state 
       do m=1,3 ! polarization    

          !corr_fn(k,:,m) = corr_fn(k,:,m) +  e_factor1(k,:) ! exp(-gamma * eV * time(:) / hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !corr_fn(k,:,m) = corr_fn(k,:,m) + D_fi(k,1,m) * D_fi(k,:,m)   !* exp(-gamma * eV * time(:) / hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !corr_fn(k,:,m) = corr_fn(k,:,m) + D_fi(k,1,m) * D_fi(k,:,m) * e_factor1(k,:)  !* exp(-gamma * eV * time(:) / hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 
          !funct = (cos(100 * eV * time(:) / hbar  )  - cos(110 * eV * time(:) / hbar  )) * exp(-gamma * eV * time(:) / hbar) 

          funct = D_fi(k,:,m) * e_factor1(k,:) * exp(-gamma * const % eV * time(:) / const % hbar) !*  exp(-gamma * (eV * time(:) / hbar)**2 ) 

          !call FFT_complex(time, funct, sigma_m(k,:,m), omega_out) 
          call FFT_complex(time, funct,  funct_out(:,m), omega_out) 

       end do ! m

       do m1=1,3
         do m2=1,3
           corr_fn_FT_tmp(k,:,m1,m2) = funct_out(:,m1) * conjg(funct_out(:,m2)) 
         end do
       end do

    end do ! k
    
    deallocate(funct_out)
    deallocate(funct, omega_out)
    deallocate(int_W_I, e_factor1)
    
  end subroutine compute_corr_fun_FT


  
  subroutine reorder_sigma(sigma)
    use m_precision, only: wp
    
    complex(wp), dimension(:,:,:), intent(inout)::sigma
    !local variables    
    integer::nfreq
    integer:: i,j
    complex(wp), dimension(:,:,:), allocatable::sigma_tmp
    
    nfreq = size(sigma,2) 

    allocate(sigma_tmp(size(sigma,1),size(sigma,2),size(sigma,3)))
    
    j=1
    do i=nfreq/2, 1, -1
       sigma_tmp(:,j,:) = sigma(:,nfreq-i+1,:) 
       j=j+1
    end do
    
    do i=0, nfreq/2-1
       sigma_tmp(:,j,:) = sigma(:,i+1,:) 
       j=j+1     
    end do
    
  end subroutine reorder_sigma


  subroutine convolute_gaussian(x,y, fwhm)
    use m_precision, only: wp
    use m_constants, only: const
   
    real(wp), dimension(:), intent(in)::x
    real(wp), dimension(:), intent(inout)::y
    real(wp):: fwhm

    real(wp), dimension(size(y)):: y2
    real(wp):: alpha
    integer::i,j

    alpha=4.0_wp*log(2.0_wp)/(fwhm**2.0_wp)

    write(6,*) fwhm, x(1), x(2)
    write(6,*) size(x), size(y), size(y2)
 
    y2 = 0.0_wp
    do i=1, size(x)
       do j=1, size(x)
          y2(i) = y2(i) + y(j) * (alpha / const % pi)**0.5_wp * exp(-alpha*(x(i)-x(j))**2.0_wp)  
       end do
    end do

    y2 = y2 * (x(2)-x(1))

    !do i=1, size(x)
    !   write(45,*) x(i), y(i), y2(i)
    !end do

    y = y2

  end subroutine convolute_gaussian


  subroutine lowpass_filter(x,dt, cutoff)
    use m_precision, only: wp
    use m_constants, only: const

    real(wp), dimension(:), intent(inout)::x
    real(wp), intent(in):: dt, cutoff

    real(wp), dimension(size(x)):: y
    real(wp):: RC, alpha
    integer::i


    ! omega = 1 / RC  -> RC = 1 / omega = hbar / (hbar*omega*eV) =  hbar / (cutoff * eV)  

    RC = const % hbar / (cutoff * const % eV)

    alpha = dt / (dt+ RC)

    y(1) = x(1)
    
    do i=2,size(x)
       y(i) = y(i-1) + alpha * (x(i) -y(i-1))
    end do

    x = y

  end subroutine lowpass_filter


  subroutine sinc_filter(x,y, cutoff)
    use m_precision, only: wp
    use m_constants, only: const
    use m_KH_functions, only: sinc

    real(wp), dimension(:), intent(in)::x
    real(wp), dimension(:), intent(inout)::y
    real(wp), intent(in):: cutoff

    real(wp), dimension(size(x)):: y2
    real(wp):: B, dx
    integer::i,j

    B = cutoff * const % ev / (2.0_wp * const %  pi * const % hbar)
    dx = x(2) -x(1)

    write(6,*) "cutoff, B", cutoff, B

    y2= 0.0_wp
    do i=1,size(x)
       do j=1,size(x)
          
          y2(i) = y2(i) + y(j) * 2 * B * sinc(2 * B *(x(i) -x(j)), 0.0_wp, 1.0_wp )
       
          !write(6,*) y2(i)
       end do
    end do

    !write(6,*) "integrals", sum(y), sum(y2*dx)

    y = y2  * dx

 
  end subroutine sinc_filter

  subroutine compute_efactor_one(E_n, E_n_mean, time, efactor, negative)
    use m_precision,only:wp     
    use m_constants, only: const
    
    real(wp), intent(in):: E_n(:), time(:), E_n_mean
    complex(wp), intent(out):: efactor(:)
    logical:: negative
    
    real(wp):: delta_t, factor
    real(wp), allocatable:: int_W_I(:)
    integer:: i, ntsteps
    
    ntsteps = size(time,1)
    
    allocate(int_W_I(ntsteps) )
    
    if (negative) then
      factor = -1.0_wp
    else
      factor = 1.0_wp
    end if
    
    delta_t = time(2)-time(1)
    
    int_W_I(1) = E_n(1) - E_n_mean  
    do i = 2, ntsteps
      int_W_I(i) = int_W_I(i-1) + ( E_n(i) - E_n_mean)  
    end do
    int_W_I =  int_W_I * delta_t
    efactor = exp(dcmplx(0, factor * (const % eV  / const % hbar) *int_W_I ))    
    
    deallocate(int_W_I)
    
  end subroutine compute_efactor_one

  subroutine compute_efactor(E_n, E_f, E_nf_mean, time, efactor, negative)
    use m_precision,only:wp
    use m_constants, only: const
    
    real(wp), intent(in):: E_f(:), E_n(:), time(:), E_nf_mean
    complex(wp), intent(out):: efactor(:)
    logical:: negative
    
    real(wp):: delta_t, factor
    real(wp), allocatable:: int_W_I(:)
    integer:: i, ntsteps
    
    ntsteps = size(time,1)
    
    allocate(int_W_I(ntsteps) )

    if (negative) then
      factor = -1.0_wp
    else
      factor = 1.0_wp
    end if
    
    delta_t = time(2)-time(1)
    
    int_W_I(1) = E_n(1) - E_f(1)  - E_nf_mean
    do i = 2, ntsteps
      int_W_I(i) = int_W_I(i-1) + ( E_n(i) - E_f(i) ) - E_nf_mean
    end do
    int_W_I =  int_W_I * delta_t
    efactor = exp(dcmplx(0, factor * (const % eV  / const % hbar) *int_W_I ))

    deallocate(int_W_I)
    
  end subroutine compute_efactor

  ! computes the integral from t to inf, complex version
  subroutine compute_int_t_inf_z(integrand, time, integral)
    use m_precision,only:wp
    use m_constants, only: const
    
    complex(wp), intent(in):: integrand(:)
    real(wp), intent(in):: time(:)
    complex(wp), intent(out):: integral(:)
    
    real(wp):: delta_t
    integer:: i, ntsteps
    
    ntsteps = size(time,1)
    
    delta_t = time(2)-time(1)

    integral(ntsteps) = integrand(ntsteps)
    do i = ntsteps-1, 1, -1
      integral(i) = integral(i+1) + integrand(i)
    end do

    integral = integral * delta_t
    
  end subroutine compute_int_t_inf_z

!  ! computes the integral from 0 to t
!  subroutine compute_int_0_t(integrand, time, integral)
!    use m_precision,only:wp
!    use m_constants, only: const
!    
!    real(wp), intent(in):: E_f(:), E_n(:), time(:), E_nf_mean
!    complex(wp), intent(out):: efactor(:)
!    logical:: negative
!    
!    real(wp):: delta_t, factor
!    real(wp), allocatable:: int_W_I(:)
!    integer:: i, ntsteps
!    
!    ntsteps = size(time,1)
!    
!    delta_t = time(2)-time(1)
!
!    integral(nsteps) = integrand(nsteps)
!    do i = nsteps-1, 1, -1
!      integral(i) = int_W_I(i+1) + integrand(i)
!    end do
!
!    integral = integral * delta_t
!    
!  end subroutine compute_int_t_inf


  
  
  subroutine read_projections(p)
    use m_precision, only: wp
    use m_io, only: get_free_handle
    use m_sckh_params_t, only: sckh_params_t

    type(sckh_params_t), intent(inout):: p   

    integer:: ifile, i
    real(kind=wp):: dnrm2

    if(p % use_proj) then
      ifile = get_free_handle()
      open(ifile, file= p % proj_file, action='read')    
      read(ifile,*) p % nproj

      allocate(p % projvec(p % nproj,3))

      do i=1, p % nproj                                         
        read(ifile,*) p % projvec(i,1), p % projvec(i,2), p % projvec(i,3) 

        !normalize projvec                               
        p % projvec(i,:) = p % projvec(i,:) / dnrm2(3,p % projvec(i,:),1)
        write(6,*) "projvector", i,  p % projvec(i,:)            
      end do

      close(ifile)
    else
      ! the three cartesian directions
      p % nproj =3 
      allocate(p % projvec(p % nproj,3))

      p % projvec =0.0_wp
      p % projvec(1,1) =1.0_wp
      p % projvec(2,2) =1.0_wp
      p % projvec(3,3) =1.0_wp

    end if

  end subroutine read_projections

  subroutine read_one_sckh_traj(ntsteps_inp, nfinal, traj_file, time_inp, &
       time_inp2, E_gs_inp, E_n_inp, E_IP1s, E_trans, E_f_inp, D_fn_inp, check_times_in)
    use m_precision, only: wp
    use m_constants, only: const
    use m_io, only: get_free_handle

    integer, intent(in):: ntsteps_inp, nfinal
    character(*), intent(in):: traj_file
    real(kind=wp), intent(in):: time_inp(:)
    real(kind=wp), intent(out)::  time_inp2(:), E_gs_inp(:),  E_n_inp(:),&
         E_IP1s(:), E_trans(:,:), E_f_inp(:,:), D_fn_inp(:,:,:)
    logical, intent(in), optional:: check_times_in

    integer:: ifile, i, j, ntrans
    character(80):: dummy
    logical:: check_times

    if(present(check_times_in)) then
      check_times =check_times_in
    else
      check_times =.true.
    end if
    
    ifile = get_free_handle()
    open(ifile,file=traj_file,status='old')  

    do i=1,ntsteps_inp

      read(ifile,*) time_inp2(i)
      read(ifile,*) E_gs_inp(i)
      read(ifile,*) E_n_inp(i)
      read(ifile,*) E_IP1s(i)  

      read(ifile,*) dummy, ntrans

      ! check that
      if ( ntrans .ne. nfinal ) then
        write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
      end if

      do j =1,nfinal
        read(ifile,*) E_trans(j,i), D_fn_inp(j,i,1), D_fn_inp(j,i,2), D_fn_inp(j,i,3)
      end do

      !compute E_f_inp
      E_f_inp(:,i) = E_gs_inp(i) - E_trans(:,i) + E_IP1s(i) * (const % eV / const % Hartree)

      !check that time_inp(i) = time_inp2(i) 
      if(check_times) then
        if ( abs(time_inp(i) - time_inp2(i)*1.d-15 ) .gt. 1.d-30) then
          !write(6,*) "Error in time too big", i, abs(time_inp(i) - time_inp2(i)*1.d-15 )
        end if
      end if
      
    end do !i

    ! convert to eV units
    E_n_inp = E_n_inp * const % hartree / const % eV
    do j=1,nfinal 
      E_f_inp(j,:) = E_f_inp(j,:) * const % hartree / const % eV 
    end do

    close(ifile)
    
  end subroutine read_one_sckh_traj


  subroutine read_one_sckh_res_traj(ntsteps_inp, nfinal, ninter, traj_file, time_inp, &
       time_inp2, E_gs_inp, E_n_inp, &
       E_IP1s, E_trans, &
       E_f_inp, D_fn_inp, &
       !E_XAS_inp, E_IP1s_XAS, E_trans_XAS, &
       D_in_inp, &
       E_n0, &
       E_fn_corr, &
       !norbs_gs, nocc_gs, eps_gs, norbs_exc, nocc_exc, eps_exc,&
       check_times_in)
    use m_precision, only: wp
    use m_constants, only: const
    use m_io, only: get_free_handle

    integer, intent(in):: ntsteps_inp, nfinal, ninter
    character(*), intent(in):: traj_file
    real(kind=wp), intent(in):: time_inp(:)
    real(kind=wp), intent(out):: time_inp2(:), E_gs_inp(:),  E_n_inp(:,:)
    real(kind=wp), intent(out):: E_f_inp(:,:), D_fn_inp(:,:,:)
    real(kind=wp), intent(out)::  D_in_inp(:,:)
    real(kind=wp), intent(out)::  E_n0(:)
    real(kind=wp), intent(out):: E_fn_corr(:,:)
    real(kind=wp), intent(out):: E_IP1s(:), E_trans(:,:)
    logical, intent(in), optional:: check_times_in

    real(kind=wp)::  E_XAS_inp, E_IP1s_XAS
!    real(wp), allocatable:: E_trans_XAS(:)
    integer:: norbs_gs, nocc_gs, norbs_exc, nocc_exc
    real(kind=wp), allocatable::  eps_gs(:,:), eps_exc(:,:)
    integer:: ifile, i, j, ntrans, jj
    character(80):: dummy
    logical:: check_times

    if(present(check_times_in)) then
      check_times =check_times_in
    else
      check_times =.true.
    end if
    
!    allocate(E_trans_XAS(ninter))
!    allocate(E_IP1s(ntsteps_inp))
!    allocate(E_trans(nfinal,ntsteps_inp))

    ifile = get_free_handle()
    open(ifile,file=traj_file,status='old')  

    ! XAS for first time step
    read(ifile,*) dummy         ! E_XAS_inp  ! total energy
    read(ifile,*) dummy         ! E_IP1s_XAS   ! 1s orbital energy
    read(ifile,*) dummy, ntrans ! number of x-ray transitions, should be the same number as the number of unocc states used

    !! check 
    !if ( ntrans .ne. ninter ) then
    !  write(6,*) "Error, ntrans != nfinal", ntrans, ninter
    !end if
    
    do j = 1 , ninter
!      read(ifile,*) E_trans_XAS(j), D_in_inp(j,1), D_in_inp(j,2), D_in_inp(j,3)
      read(ifile,*) dummy, D_in_inp(j,1), D_in_inp(j,2), D_in_inp(j,3)
    end do

    ! now read trajectory with XES and more stuff
    do i = 1 , ntsteps_inp

      read(ifile,*) time_inp2(i)
      read(ifile,*) E_gs_inp(i)
      read(ifile,*) E_n0(i)
      read(ifile,*) E_IP1s(i)  
      ! XES

      read(ifile,*) dummy, ntrans

      !! check that
      !if ( ntrans .ne. nfinal ) then
      !  write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
      !end if

      do j = 1 , nfinal
        jj = j + ntrans-nfinal 
        read(ifile,*) E_trans(jj,i), D_fn_inp(jj,i,1), D_fn_inp(jj,i,2), D_fn_inp(jj,i,3)
      end do

      ! all orbital energies for the ground state
      read(ifile,*) norbs_gs, nocc_gs
      ! add by O.Takahashi 2018/07/04
      if ( i .eq. 1 ) then
        allocate(eps_gs(norbs_gs,ntsteps_inp))
        ! check
        if ( norbs_gs .lt. nocc_gs ) then
          write(6,*) "Error, norbs_gs should be larger than nocc_gs", &
               norbs_gs, nocc_gs
          stop
        end if
        if ( norbs_gs .lt. nocc_gs+ninter-1 ) then
          write(6,*) "Error, norbs_gs should be larger than nocc_gs+ninter-1", &
               norbs_gs, nocc_gs, ninter 
          stop
        end if
      end if

      do j =1, norbs_gs
        read(ifile,*) eps_gs(j,i)
      end do

      ! all orbital energies for the excited state
      read(ifile,*) norbs_exc, nocc_exc
      ! add by O.Takahashi 2018/07/04
      if ( i .eq. 1 ) then
        allocate(eps_exc(norbs_exc,ntsteps_inp))
        ! check
        if ( norbs_exc .lt. nocc_exc ) then
          write(6,*) "Error, norbs_exc should be larger than nocc_exc", &
               norbs_exc, nocc_exc
          stop
        end if
        if ( norbs_exc .lt. nocc_exc+ninter-1 ) then
          write(6,*) "Error, norbs_exc should be larger than nocc_exc+ninter-1", &
               norbs_exc, nocc_exc, ninter 
          stop
        end if
      end if

      do j =1, norbs_exc
        read(ifile,*) eps_exc(j,i)
      end do
      
      ! compute E_f_inp (will later be corrected through orbital energies)
      E_f_inp(:,i) = E_gs_inp(i) - E_trans(:,i) + E_IP1s(i) * (const % eV / const % Hartree)

      ! corrections to the intermediate state E_n0 (note that we start at HOMO!)
      ! modified by O.Takahashi 2018/07/10
      !E_n_inp(:,i) = eps_exc(nocc_exc:nocc_exc+ninter-1 ,i) -eps_exc(nocc_exc,i)
      E_n_inp(:,i) = E_n0(i) + eps_exc(nocc_exc:nocc_exc+ninter-1,i) - eps_exc(nocc_exc,i)
      !E_n_inp(:,i) = E_n0(i) + (E_n0(i) - E_n0(1)) + E_trans_XAS(:) - E_trans_XAS(1)

      ! corrections to final state depending on the interemediate state
      E_fn_corr(:,i) = eps_gs(nocc_gs+1:nocc_gs+ninter,i)
      ! modified by O.Takahashi 2018/07/10
      !E_fn_corr(:,i) = (E_n0(i) - E_n0(1)) + E_trans_XAS(:) - E_trans_XAS(1)
      
      !check that time_inp(i) = time_inp2(i) 
      if(check_times) then
        if ( abs(time_inp(i) - time_inp2(i)*1.d-15 ) .gt. 1.d-30) then
          !write(6,*) "Error in time too big", i, abs(time_inp(i) - time_inp2(i)*1.d-15 )
        end if
      end if
      
    end do !i

    ! convert to eV units
    ! modified by O.Takahashi 2018/07/09
    E_gs_inp = E_gs_inp * const % hartree / const % eV
    E_n0 = E_n0 * const % hartree / const % eV

    E_n_inp = E_n_inp * const % hartree / const % eV
    do j = 1 , nfinal 
      E_f_inp(j,:) = E_f_inp(j,:) * const % hartree / const % eV 
!      E_fn_corr(j,:) = E_fn_corr(j,:) * const % hartree / const % eV 
    end do
    do j = 1 , ninter
      E_fn_corr(j,:) = E_fn_corr(j,:) * const % hartree / const % eV 
    end do

    close(ifile)
    
    ! add by O.Takahashi 2018/10/20
!    deallocate(E_trans_XAS)
    deallocate(eps_gs)
    deallocate(eps_exc)


  end subroutine read_one_sckh_res_traj


  
end module m_SCKH_utils

