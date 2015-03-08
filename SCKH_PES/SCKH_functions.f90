module SCKH_functions
  use parameters
  use spline_m
  use FFT_m
  !use trajectory_class
  implicit none

contains

  subroutine sample_x_mom(x, funct, x_sampl, mom_sampl, x_mom_sampl, mode)
    real(kind=wp), dimension(:), intent(in):: x, funct
    real(kind=wp), dimension(:), intent(out):: x_sampl, mom_sampl
    real(kind=wp), dimension(:,:), intent(out)::x_mom_sampl
    integer, intent(in):: mode
    !local variables
    integer:: npoints, npoints_x_sampl, npoints_mom_sampl, npoints_pad, npoints_pad_pow, &
         npoints_x_mom_sampl
    real(kind=wp):: dx, pes_l
    real(kind=wp),dimension(:),allocatable:: funct_real, funct_imag, x_mom
    integer::i,j,k
    integer,dimension(1)::dime
    integer,dimension(2)::dime2

    dime = shape(x)
    npoints = dime(1)
    dime = shape(x_sampl)
    npoints_x_sampl = dime(1)
    dime = shape(mom_sampl)
    npoints_mom_sampl = dime(1)
    dime2= shape(x_mom_sampl)
    npoints_x_mom_sampl = dime2(1)


    write(6,*) npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl

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
       x_mom(j) = -2 * pi * i * hbar / pes_l 
       j=j+1
    end do

    do i=0, npoints_pad/2-1
       x_mom(j) =  2 * pi * i * hbar / pes_l 
       j=j+1
    end do

    call reorder_sigma1(funct_real)
    call reorder_sigma1(funct_imag)

    do i=1,npoints_pad
       write(19,'(5ES16.6)') x_mom(i), funct_real(i), funct_imag(i), &
            funct_real(i) ** 2 + funct_imag(i) ** 2, &
            (funct_real(i) ** 2 + funct_imag(i) ** 2) * 4 * pi * x_mom(i)**2
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
    real(kind=wp), dimension(:), intent(in):: x, funct, E_i_inp
    real(kind=wp), dimension(:), intent(out):: x_sampl
    real(kind=wp), dimension(:,:), intent(out):: x_mom_sampl
    real(kind=wp):: E_gs, my_SI
    !local variables
    integer:: npoints, npoints_x_sampl, npoints_mom_sampl, npoints_pad, npoints_pad_pow,&
         npoints_x_mom_sampl
    real(kind=wp):: dx, pes_l
    real(kind=wp),dimension(:),allocatable:: funct_real, funct_imag, x_mom, V_sampl
    integer::i,j
    integer,dimension(1)::dime
    integer,dimension(2)::dime2

    ! sample space and assign momentum according to conservation of energy
    ! does not seem to work very well! especially since the quantum space probablility distribution
    ! lies outside of the classically allowed region
    
    dime = shape(x)
    npoints = dime(1)
    dime = shape(x_sampl)
    npoints_x_sampl = dime(1)
    dime2 = shape(x_mom_sampl)
    npoints_x_mom_sampl = dime2(1)

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
          x_mom_sampl(i,2) = sqrt( 2 * my_SI * (E_gs - V_sampl(i)) )
       else
          x_mom_sampl(i,2) = 0
       end if
    end do
    do i=1,npoints_x_sampl !+1,npoints_x_sampl*2
       x_mom_sampl(i + npoints_x_sampl,1) = x_sampl(i)
       if( E_gs - V_sampl(i) .gt. 0d0) then
          x_mom_sampl(i + npoints_x_sampl,2) = -sqrt( 2 * my_SI * (E_gs - V_sampl(i) ))  
       else
          x_mom_sampl(i + npoints_x_sampl,2) = 0
       end if
    end do

    deallocate(V_sampl)
  end subroutine sample_x_mom2

  
  subroutine sample_even(x, funct, x_sampl)
    ! passed variables
    real(kind=wp), dimension(:), intent(in)::x, funct
    real(kind=wp), dimension(:), intent(out)::x_sampl
    ! local variables
    integer,parameter:: npoints_new=10000
    integer:: npoints, npoints_sampl
    real(kind=wp), dimension(npoints_new)::funct_spl, I_x, x_new
    real(kind=wp), dimension(:),allocatable::points_sampl
    integer,dimension(1):: dime
    integer:: i,j
    

    dime = shape(x)
    npoints =dime(1)
    dime = shape(x_sampl)    
    npoints_sampl = dime(1)
   
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
    ! passed variables
    real(kind=wp), dimension(:), intent(in)::x, funct
    real(kind=wp), dimension(:), intent(out)::x_sampl
    ! local variables
    integer,parameter:: npoints_new=10000
    integer:: npoints, npoints_sampl
    real(kind=wp), dimension(npoints_new)::funct_spl, I_x, x_new
    real(kind=wp), dimension(:),allocatable::points_sampl,yder
    real(kind=wp):: funct_x_tmp, x_tmp, y_tmp, tol, x_min, x_max,y_min, y_max, x_l, y_l,ran1,ran2
    integer,dimension(1):: dime
    integer:: i,j
    

    dime = shape(x)
    npoints =dime(1)
    dime = shape(x_sampl)    
    npoints_sampl = dime(1)
   
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
  
  subroutine verlet_trajectory(x_in, v_in, V_x, V_y, dt, my_SI, x_out)
    implicit none
    ! passed variables
    real(kind=wp), dimension(:), intent(in):: V_x, V_y
    real(kind=wp), intent(in):: x_in,v_in, dt, my_SI    
    real(kind=wp), dimension(:), intent(out)::x_out
        
    ! local variables
    integer:: npoints, ntsteps, V_len
    integer:: i,j
    integer,dimension(1):: dime
    real(kind=wp):: a, x, v, x_new, v_new, a_new,t
    
    dime = shape(V_x)
    V_len =dime(1)
    
    dime = shape(x_out)
    ntsteps = dime(1)
    
    ! first force evaluation on PES
    x = x_in
    v = v_in
    a = force(x, V_x, V_y ) / my_SI

    ! write first (zeroth) step to trajectory
    t=0
    !call trajectory_add(traj,x,v,a,t)
    
    i=0
    !write(14,'(I4,4ES16.6)') i, t, x, v, a
    ! run n_steps steps
    j=0
    do i =1,ntsteps
       call verlet_step(x,v,a, V_x, V_y, V_len, x_new, v_new, a_new, dt, my_SI)
       
       t=t+dt
       x = x_new
       v= v_new
       a = a_new
       
       j=j+1
       x_out(j) = x 

       !write(14,'(I4,4ES16.6)') i, t, x, v, a
       !call trajectory_add(traj,x,v,a,t)
       
    end do
    
  end subroutine verlet_trajectory
  
  
  subroutine verlet_step(x,v,a, V_x, V_y, V_len, x_new, v_new, a_new, dt, my_SI)
    implicit none 
    ! passed variables
    real(wp), intent(in):: x,v,a,my_SI, dt
    real(wp), intent(out):: x_new,v_new,a_new
    integer, intent(in)::V_len
    real(wp), dimension(V_len), intent(in):: V_x,V_y

    ! velocity verlet
    ! 1) a(t) = F/m  
    ! 2) calcualte x(t+dt) 
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
    implicit none
    !passed
    real(wp), dimension(:), intent(in):: V_x,V_y
    real(wp), intent(in) :: x
    real(wp), intent(in), optional::dx_in
    !local
    integer::V_len,i
    real(wp), dimension(2):: x2, V2 
    real(wp)::dx    
    integer, dimension(1)::dime

    if (present(dx_in)) then
       dx=dx_in
    else
       dx=1.0d-13
    end if

    dime = shape(V_x)
    V_len = dime(1)
    
    ! spline two nearby points on he PES
    ! calculate force by central finite difference formula:
    ! F = - d/dx V(x(t))  \approx - (V(x(t)+h/2) - V(x(t)-h/2) ) / h   
    
    x2(1) = x - dx/2
    x2(2) = x + dx/2
    
    call spline_easy(V_x, V_y, V_len, x2, V2, 2)

    force = - (V2(2) -V2(1)) / dx
  
  end function force
  
  
  
  subroutine compute_SCKH(E_n, E_f, E_fn_mean, D_fn, time, sigma_m, gamma)
    ! passed variables
    real(kind=wp), dimension(:), intent(in):: E_n,time
    real(kind=wp), dimension(:,:), intent(in):: E_f
    real(kind=wp),dimension(:,:,:),intent(in):: D_fn 
    complex(kind=wp), dimension(:,:,:),intent(out) ::  sigma_m
    real(kind=wp), intent(in):: E_fn_mean, gamma
    ! local variables
    integer:: ntsteps, ntsteps_pad, nfinal
    complex(kind=wp), dimension(:),allocatable:: funct
    complex(kind=wp), dimension(:,:),allocatable::  e_factor1
    real(kind=wp), dimension(:,:),allocatable:: int_W_I
    real(kind=wp),dimension(:),allocatable:: funct_real, funct_imag, omega_out 
    real(kind=wp):: delta_t
    integer:: i,k,m 
    integer, dimension(1)::dime
    integer, dimension(2)::dime2
    integer, dimension(3)::dime3

    dime = shape(time)
    ntsteps = dime(1)
    dime2 = shape(E_f)
    nfinal = dime2(1)
    dime3 = shape(sigma_m)
    ntsteps_pad = dime3(2)
    
    allocate(funct(ntsteps), int_W_I(nfinal, ntsteps), &
         e_factor1(nfinal,ntsteps), omega_out(ntsteps_pad) )
    
   
    delta_t = time(2)-time(1)
        
    int_W_I(:,1) = E_f(:,1) - E_n(1)  
    
    do i = 2, ntsteps
       int_W_I(:,i) = int_W_I(:,i-1) + ( E_f(:,i) - E_n(i)) - E_fn_mean  
    end do
    int_W_I =  int_W_I * delta_t
    
    do i = 1,nfinal 
       e_factor1(i,:) = exp(dcmplx(0, -(eV  / hbar) *int_W_I(i,:)  ))
    end do
    
    ! compute A_{fm}(t) = D^m_{fn}(t) * e_factor_f(t) * exp(-gamma * eV * time(:) / hbar)
    ! and fourer transform
    
    sigma_m = 0
    
    !do a FFT
    do k =1,nfinal! final state 
       do m=1,3 ! polarization     
          
          funct = D_fn(k,:,m) * e_factor1(k,:) * exp(-gamma * eV * time(:) / hbar)

          call FFT_complex(time,funct, sigma_m(k,:,m), omega_out)

       end do ! m
    end do ! k
    
    deallocate(funct, int_W_I)
    
  end subroutine compute_SCKH
  
  subroutine reorder_sigma(sigma)
    complex(kind=wp), dimension(:,:,:), intent(inout)::sigma
    !local variables    
    integer::nfreq
    integer:: i,j
    integer, dimension(3)::dime3
    complex(kind=wp), dimension(:,:,:), allocatable::sigma_tmp
    
    dime3= shape(sigma)
    nfreq = dime3(2)

    allocate(sigma_tmp(dime3(1),dime3(2),dime3(3)))
    
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
  
end module SCKH_functions

