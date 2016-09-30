module m_SCKH_resonant_PES
  implicit none

contains

  subroutine calculate_SCKH_res_PES(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    !use m_SCKH_utils, only: verlet_trajectory
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_SCKH_utils, only: compute_F_if_omp_many_n
    use m_SCKH_utils, only: compute_F_if_om_many_n
    use m_SCKH_utils, only: compute_F_if_om_omp
    use m_SCKH_utils, only: compute_F_if_om_omp_ingoing
    use m_sckh_params_t, only: sckh_params_t 
    !use hist_class, only: hist, hist_init, hist_add
    !use hist_class, only: hist_broadening, hist_write
    use m_io, only: get_free_handle
    use m_splines, only: spline_easy
    use m_PES_io, only: read_dipole_file
    !use m_PES_io, only: read_nac_file
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: get_projections
    use m_PES_io, only: read_file_list
    use m_KH_functions, only: solve_vib_problem
    use m_fftw3, only: get_omega_reordered_fftw
    use m_upper, only : upper
    use m_KH_utils, only: convolute_incoming
    use m_KH_utils, only: convolute_instrumental
    
    type(sckh_params_t), intent(inout):: p 

    integer:: ninter, nfinal, ntsteps  
    integer:: npoints_in

    real(wp), allocatable::  time(:), E_dyn_inp(:),  E_dyn2_inp(:)  
    real(wp), allocatable::  E_lp_corr(:), shift(:)
    real(wp), allocatable::  E_i_inp(:), E_i1(:), E_i2(:)
    real(wp), allocatable::  E_n_inp(:,:), E_n1(:,:), E_n2(:,:)
    real(wp), allocatable::  E_f_inp(:,:), E_f1(:,:), E_f2(:,:)

    real(wp), allocatable:: D_fn_inp(:,:,:,:), D_fn1(:,:,:,:), D_fn2(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:), D_ni1(:,:,:), D_ni2(:,:,:) 
    real(wp):: E_nf_mean, E_fi_mean, E_ni_mean
    
    character(80)::  string, file
    real(wp), allocatable:: x_sampl(:), mom_sampl(:), x_mom_sampl(:,:)
    real(wp), allocatable:: x_new(:), v_new(:), a_new(:)
    real(wp), allocatable:: x_new2(:), v_new2(:), a_new2(:)
    integer:: npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl

    real(wp), allocatable:: omega_in(:), omega_out(:)
    integer:: n_omega_in, n_omega_out

    real(wp), allocatable:: c_i(:,:)
    real(wp), allocatable:: eig_i(:)
    
    real(wp):: gamma, gamma_inc, gamma_instr
    real(wp):: time_l, delta_t, norm, mu_SI, dx
    real(wp):: fac_t
    
    !type(hist), dimension(:), allocatable:: time_h
    !type(hist):: time_h_0, time_h_0_mom
    integer, dimension(1)::ind  
    integer:: ifile
    
    real(wp), allocatable::  X_r(:)
    real(wp) :: dvr_start 
    integer:: npoints, ii
    real(8):: dnrm2
    integer::i,j,m,m1,m2, traj, traj2

    complex(wp), allocatable::  F_if_t_omp(:,:,:,:,:)    
    complex(wp), allocatable::  F_if_om_omp(:,:,:,:,:)    
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)
    
    ! set some local variables
    ntsteps = p % ntsteps
    npoints_x_sampl = p % npoints_x_sampl
    npoints_mom_sampl =  p % npoints_mom_sampl
    npoints_in = p % npoints_in ! this should be renamed, corresponds to p % nstates
    nfinal = p % npesfile_f
    ninter = p % npesfile_n
    mu_SI = p % mu * const % u

    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    delta_t = p % delta_t

    ! use HWHM internally
    gamma = p % gamma_FWHM / 2
    gamma_inc = p % gamma_inc_FWHM / 2
    gamma_instr = p % gamma_instr_FWHM / 2 

    ! projections
    call  get_projections(p)

    allocate( E_i_inp(npoints_in),&
         E_n_inp(ninter,npoints_in), &
         E_dyn_inp(npoints_in),&
         E_dyn2_inp(npoints_in),&
         E_f_inp(nfinal,npoints_in), &
         D_ni_inp(ninter, npoints_in,3), &
         D_fn_inp(nfinal, ninter, npoints_in,3), &
         time(ntsteps),&
         E_i1(ntsteps), &
         E_i2(ntsteps), &
         E_n1(ninter, ntsteps), &
         E_n2(ninter, ntsteps), &
         E_f1(nfinal,ntsteps), &
         E_f2(nfinal,ntsteps), &
         D_ni1(ninter, ntsteps,3), &
         D_ni2(ninter, ntsteps,3), &
         D_fn2(nfinal, ninter, ntsteps,3), &
         D_fn1(nfinal, ninter, ntsteps,3), &
         E_lp_corr(npoints_in),&
         shift(npoints_in), &
         c_i(npoints_in,npoints_in), &
         eig_i(npoints_in), &
         x_new(ntsteps),&
         v_new(ntsteps),&
         a_new(ntsteps),&
         x_new2(ntsteps),&
         v_new2(ntsteps),&
         a_new2(ntsteps),&
         !time_h(ntsteps)&
         )
    allocate(X_r(npoints_in))

    ! set up grid points
    do i = 1, npoints_in
      X_r(i) = (i-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_r, E_i_inp)
    !call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n_inp)

    ! read list of intermediate state pes_files and dipole_files
    call read_file_list(p % pes_file_list_n, p % npesfile_n, p % pes_files_n)
    call read_file_list(p % dipole_file_list_n, p % npesfile_n, p % dipolefile_n)

    do j=1,p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, p % npoints_in, X_r, E_n_inp(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % npoints_in, X_r, D_ni_inp(j,:,:))
    end do
    
    ! read PES file where the first and second dynamics are run 
    if (p % use_dynamics_file) then
      call read_PES_file(p % pes_file_dyn, p % npoints_in, p % npoints_in, X_r, E_dyn_inp)
      call read_PES_file(p % pes_file_dyn2, p % npoints_in, p % npoints_in, X_r, E_dyn2_inp)
    else
      E_dyn_inp = E_i_inp
      E_dyn2_inp = E_n_inp(1,:)
    end if

    ! read list of final state pes_files and dipole_files
    call read_file_list(p % pes_file_list_f, p % npesfile_f, p % pes_files_f)
    call read_file_list(p % dipole_file_list_f, p % npesfile_f, p % dipolefile_f)

    do j=1,p % npesfile_f
      call read_PES_file(p % pes_files_f(j), p % npoints_in, p % npoints_in, X_r, E_f_inp(j,:))
      ! temporary hack: only one intermediate state
      call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % npoints_in, X_r, D_fn_inp(j,1,:,:))
    end do

    ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved
    if( p % shift_PES .eq.  1) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % npoints_in, X_r, E_lp_corr)

      shift = E_lp_corr -E_f_inp(1,:) 

      do j=1,p % npesfile_f
        E_f_inp(j,:) = E_f_inp(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if

    ! Solve the vibrational problem for initial state to be able to sample initial distribution
    call solve_vib_problem(dx, E_i_inp, eig_i, c_i, mu_SI, p % vib_solver)
    write(6,*) "Calculated initial state eigenfunctions"
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm

    ! convert to eV units
    E_i_inp = E_i_inp  / const % eV
    E_n_inp = E_n_inp  / const % eV
    E_dyn_inp = E_dyn_inp  / const % eV
    E_dyn2_inp = E_dyn2_inp  / const % eV
    E_f_inp = E_f_inp / const % eV

    ifile = get_free_handle()
    open(ifile, file="inital_state_eigvec.txt", action='write')
    do i=1,npoints_in
      write(ifile,'(3ES16.6)') X_r(i), c_i(i,1), c_i(i,1) ** 2
    end do
    close(ifile)

    ! sample the positions and momenta
    call sample_x_mom_modes(npoints_x_sampl, npoints_mom_sampl, &
         p % samplemode, X_r, c_i(:,1), x_mom_sampl)
    npoints_x_mom_sampl = size(x_mom_sampl,1)

    delta_t = delta_t * 1.d-15 ! femtoseconds
    time_l = (ntsteps-1) * delta_t

    write(6,*) "outfile", p % outfile
    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
    write(6,*) "mu_SI", mu_SI 
    write(6,*) "time_l", time_l
    write(6,*)
    write(6,*) "Fundamental frequency resolution", 2.0_wp * const % pi * const % hbar /( time_l * const % eV)
    write(6,*) "delta t", delta_t
    write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV), "eV"

    n_omega_in = ntsteps
    n_omega_out = ntsteps
    allocate(omega_in(n_omega_in))
    allocate(omega_out(n_omega_out))
    allocate(F_if_t_omp(nfinal, n_omega_in, n_omega_out,3,3))
    allocate(F_if_om_omp(nfinal, n_omega_in, n_omega_out,3,3))
    allocate(lambda_F(nfinal, n_omega_in, n_omega_out))
    allocate(lambda_G(nfinal, n_omega_in, n_omega_out))
    allocate(lambda_H(nfinal, n_omega_in, n_omega_out))
    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))

    !
    ! Loop over trajectories
    !

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do

    ! introduce factor for possible time development backwards in time
    write(6,*) "upper(p % KH_amplitude_mode)", upper(p % KH_amplitude_mode)
    if(upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
      fac_t = 1.0_wp
    else if(upper(p % KH_amplitude_mode) .eq. "INGOING") then
      fac_t = 1.0_wp
    else
      write(6,*) "p % KH_amplitude_mode must be either OUTGOING or INGOING"
    end if
      
    !call hist_init(time_h_0, 1000, X_r(1), X_r(npoints_in) ) 
    !call hist_init(time_h_0_mom, 1000, minval(x_mom_sampl(:,2)), maxval(x_mom_sampl(:,2)) ) 
    !do i=1, ntsteps
    !   call hist_init(time_h(i), 1000, X_r(1), X_r(npoints_in) ) 
    ! end do

    lambda_F = 0.0_wp
    lambda_G = 0.0_wp
    lambda_H = 0.0_wp

    do traj=1, npoints_x_mom_sampl

      write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl
      
      call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
           E_dyn_inp * const % eV, fac_t * delta_t, mu_SI, x_new, v_new, a_new )
      
      call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
      do i=1,ninter
        call spline_easy(X_r, E_n_inp(i,:), npoints_in, x_new, E_n1(i,:), ntsteps)  
      end do

      do i=1,ninter
        do m=1,3
          call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new, D_ni1(i,:,m) , ntsteps)  
        end do
      end do
      
      do i=1,nfinal
        do m=1,3
          call spline_easy(X_r, D_fn_inp(i,1,:,m) , npoints_in, x_new, D_fn1(i,1,:,m) , ntsteps)  
        end do
      end do
      
      do i=1,nfinal
        call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_f1(i,:), ntsteps)  
      end do

      ! first time, compute the mean transition energy, and frequency
      if (traj .eq. 1) then
        ind = minloc(E_i_inp)

        E_ni_mean =  E_n1(nfinal, ind(1)) - E_i1(ind(1))
        write(6,*) "E_ni_mean", E_ni_mean
        
        E_fi_mean =  E_f1(nfinal, ind(1)) - E_i1(ind(1))
        write(6,*) "E_fi_mean", E_fi_mean
        
        E_nf_mean =  E_n1(ninter, ind(1)) - E_f1(nfinal,ind(1))
        write(6,*) "E_nf_mean", E_nf_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
        omega_in = omega_in + E_ni_mean !E_nf_mean + E_fi_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
        omega_out = omega_out + E_nf_mean
        
      end if

      ! loop over each starting point of first trajectory
      do traj2 =1, ntsteps ! possibly with a stride
        
        call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
             E_dyn2_inp * const % eV, fac_t * delta_t, mu_SI, x_new2, v_new2, a_new2 )

        !call verlet_trajectory_xva(x_new(1), v_new(1), X_r, &
        !     E_dyn2_inp * const % eV, fac_t * delta_t, mu_SI, x_new2, v_new2, a_new2 )
        
        call spline_easy(X_r, E_n_inp, npoints_in, x_new2, E_n2, ntsteps)
        call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
        
        do i=1,nfinal
          call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new2, E_f2(i,:), ntsteps)  
        end do
        
        do i=1,nfinal
          do m=1,3
            call spline_easy(X_r, D_fn_inp(i,1,:,m) , npoints_in, x_new2, D_fn2(i,1,:,m) , ntsteps)  
          end do
        end do

        do i=1,ninter
          do m=1,3
            call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new2, D_ni2(i,:,m) , ntsteps)  
          end do
        end do

        
        if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

          call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,:,:,:), &
               D_ni2(:, 1, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)

        else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
          ! F_if_t_omp will now store F_if_t_om temporarily
          call compute_F_if_om_many_n(E_n2(:,:), E_i2(:), E_ni_mean, D_fn2(:,1,:,:), &
               D_ni2(:, :, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)
          
        end if
        
      end do ! do traj2 =1, ntsteps2


      if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

        ! for test purposes, write F_if_t_omp to file
        call write_F_if_t_omp(p, F_if_t_omp, omega_out, time)
        
        call compute_F_if_om_omp(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
             E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp)

      else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
        call compute_F_if_om_omp_ingoing(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
             E_i1, gamma_instr, omega_in, E_ni_mean, F_if_om_omp)
        
      end if
        
        
      ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
      do m1=1,3
        do m2=1,3
          lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m1)) * F_if_om_omp(:,:,:,m2,m2))
          lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m1,m2))
          lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m2,m1))
        end do
      end do

    end do ! end traj
    
    ! averages according to J. Phys. B. 27, 4169 (1994) 
    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
    lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
    lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
    lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

      !if(gamma_inc .gt. 1d-5) then
      !  sigma_tmp = lambda_lp
      !  call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
      !  sigma_tmp = lambda_ln
      !  call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
      !  sigma_tmp = lambda_cp
      !  call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
      !end if
      !
      !if(gamma_instr .gt. 1d-5) then
      !  sigma_tmp = lambda_lp
      !  call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_lp)
      !  sigma_tmp = lambda_ln
      !  call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_ln)
      !  sigma_tmp = lambda_cp
      !  call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_cp)
      !end if
      
    else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then      
      allocate(sigma_tmp(n_omega_in, n_omega_out))

      if(gamma_inc .gt. 1d-5) then
        sigma_tmp = lambda_lp
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
        sigma_tmp = lambda_ln
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
        sigma_tmp = lambda_cp
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
      end if
      
    end if
    
    ! write spectra to file
    do j=1, n_omega_in
      
      file="_sigma_"
      write(string,'(F6.2)') omega_in(j)   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1, n_omega_out
        
        write(ifile,'(4ES18.10)') omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do
      
      close(ifile) 
      
   end do

  end subroutine calculate_SCKH_res_PES

  subroutine calculate_SCKH_res_PES_factor(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    !use m_SCKH_utils, only: verlet_trajectory
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_SCKH_utils, only: compute_F_if_omp_many_n
    !use m_SCKH_utils, only: compute_F_if_om_many_n
    !use m_SCKH_utils, only: compute_F_if_om_omp
    use m_SCKH_utils, only: compute_F_if_om_omp_no_F
    !use m_SCKH_utils, only: compute_F_if_om_omp_ingoing
    use m_sckh_params_t, only: sckh_params_t 
    !use hist_class, only: hist, hist_init, hist_add
    !use hist_class, only: hist_broadening, hist_write
    use m_io, only: get_free_handle
    use m_splines, only: spline_easy
    use m_PES_io, only: read_dipole_file
    !use m_PES_io, only: read_nac_file
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: get_projections
    use m_PES_io, only: read_file_list
    use m_KH_functions, only: solve_vib_problem
    use m_fftw3, only: get_omega_reordered_fftw
    use m_upper, only : upper
    use m_KH_utils, only: convolute_incoming
    !use m_KH_utils, only: convolute_incoming_fft
    use m_KH_utils, only: convolute_instrumental
    use m_spectrum_utils, only: convolution_lorentzian_grid_fft_many_freq
    
    type(sckh_params_t), intent(inout):: p 

    integer:: ninter, nfinal, ntsteps  
    integer:: npoints_in

    real(wp), allocatable::  time(:), E_dyn_inp(:),  E_dyn2_inp(:)  
    real(wp), allocatable::  E_lp_corr(:), shift(:)
    real(wp), allocatable::  E_i_inp(:), E_i1(:), E_i2(:)
    real(wp), allocatable::  E_n_inp(:,:), E_n1(:,:), E_n2(:,:)
    real(wp), allocatable::  E_f_inp(:,:), E_f1(:,:), E_f2(:,:)

    real(wp), allocatable:: D_fn_inp(:,:,:,:), D_fn1(:,:,:,:), D_fn2(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:), D_ni1(:,:,:), D_ni2(:,:,:) 
    real(wp):: E_nf_mean, E_fi_mean, E_ni_mean
    
    character(80)::  string, file
    real(wp), allocatable:: x_sampl(:), mom_sampl(:), x_mom_sampl(:,:)
    real(wp), allocatable:: x_new(:), v_new(:), a_new(:)
    real(wp), allocatable:: x_new2(:), v_new2(:), a_new2(:)
    integer:: npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl

    real(wp), allocatable:: omega_in(:), omega_out(:)
    integer:: n_omega_in, n_omega_out

    real(wp), allocatable:: c_i(:,:)
    real(wp), allocatable:: eig_i(:)
    
    real(wp):: gamma, gamma_inc, gamma_instr
    real(wp):: time_l, delta_t, norm, mu_SI, dx
    real(wp):: fac_t
    
    !type(hist), dimension(:), allocatable:: time_h
    !type(hist):: time_h_0, time_h_0_mom
    integer, dimension(1)::ind  
    integer:: ifile
    
    real(wp), allocatable::  X_r(:)
    real(wp) :: dvr_start 
    integer:: npoints, ii
    real(8):: dnrm2
    integer::i,j,m,m1,m2, traj, traj2, om_in, f_e

    complex(wp), allocatable::  F_if_t_omp(:,:,:,:,:)    
    complex(wp), allocatable::  F_if_om_omp(:,:,:,:,:)
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)


    complex(wp), allocatable::  R_if_om_omp(:,:,:)    
    real(wp), allocatable:: R_factor(:,:,:)
    
    ! set some local variables
    ntsteps = p % ntsteps
    npoints_x_sampl = p % npoints_x_sampl
    npoints_mom_sampl =  p % npoints_mom_sampl
    npoints_in = p % npoints_in ! this should be renamed, corresponds to p % nstates
    nfinal = p % npesfile_f
    ninter = p % npesfile_n
    mu_SI = p % mu * const % u

    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    delta_t = p % delta_t

    ! use HWHM internally
    gamma = p % gamma_FWHM / 2
    gamma_inc = p % gamma_inc_FWHM / 2
    gamma_instr = p % gamma_instr_FWHM / 2 

    ! projections
    call  get_projections(p)

    allocate( E_i_inp(npoints_in),&
         E_n_inp(ninter,npoints_in), &
         E_dyn_inp(npoints_in),&
         E_dyn2_inp(npoints_in),&
         E_f_inp(nfinal,npoints_in), &
         D_ni_inp(ninter, npoints_in,3), &
         D_fn_inp(nfinal, ninter, npoints_in,3), &
         time(ntsteps),&
         E_i1(ntsteps), &
         E_i2(ntsteps), &
         E_n1(ninter, ntsteps), &
         E_n2(ninter, ntsteps), &
         E_f1(nfinal,ntsteps), &
         E_f2(nfinal,ntsteps), &
         D_ni1(ninter, ntsteps,3), &
         D_ni2(ninter, ntsteps,3), &
         D_fn2(nfinal, ninter, ntsteps,3), &
         D_fn1(nfinal, ninter, ntsteps,3), &
         E_lp_corr(npoints_in),&
         shift(npoints_in), &
         c_i(npoints_in,npoints_in), &
         eig_i(npoints_in), &
         x_new(ntsteps),&
         v_new(ntsteps),&
         a_new(ntsteps),&
         x_new2(ntsteps),&
         v_new2(ntsteps),&
         a_new2(ntsteps),&
         !time_h(ntsteps)&
         )
    allocate(X_r(npoints_in))

    ! set up grid points
    do i = 1, npoints_in
      X_r(i) = (i-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_r, E_i_inp)
    !call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n_inp)

    ! read list of intermediate state pes_files and dipole_files
    call read_file_list(p % pes_file_list_n, p % npesfile_n, p % pes_files_n)
    call read_file_list(p % dipole_file_list_n, p % npesfile_n, p % dipolefile_n)

    do j=1,p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, p % npoints_in, X_r, E_n_inp(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % npoints_in, X_r, D_ni_inp(j,:,:))
    end do
    
    ! read PES file where the first and second dynamics are run 
    if (p % use_dynamics_file) then
      call read_PES_file(p % pes_file_dyn, p % npoints_in, p % npoints_in, X_r, E_dyn_inp)
      call read_PES_file(p % pes_file_dyn2, p % npoints_in, p % npoints_in, X_r, E_dyn2_inp)
    else
      E_dyn_inp = E_i_inp
      E_dyn2_inp = E_n_inp(1,:)
    end if

    ! read list of final state pes_files and dipole_files
    call read_file_list(p % pes_file_list_f, p % npesfile_f, p % pes_files_f)
    call read_file_list(p % dipole_file_list_f, p % npesfile_f, p % dipolefile_f)

    do j=1,p % npesfile_f
      call read_PES_file(p % pes_files_f(j), p % npoints_in, p % npoints_in, X_r, E_f_inp(j,:))
      ! temporary hack: only one intermediate state
      call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % npoints_in, X_r, D_fn_inp(j,1,:,:))
    end do

    ! XXX dangerousn hack!
    ! E_dyn_inp = (E_i_inp + E_f_inp(1,:)) /2.0_wp
    


    
    ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved
    if( p % shift_PES .eq.  1) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % npoints_in, X_r, E_lp_corr)

      shift = E_lp_corr -E_f_inp(1,:) 

      do j=1,p % npesfile_f
        E_f_inp(j,:) = E_f_inp(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if

    ! Solve the vibrational problem for initial state to be able to sample initial distribution
    call solve_vib_problem(dx, E_i_inp, eig_i, c_i, mu_SI, p % vib_solver)
    write(6,*) "Calculated initial state eigenfunctions"
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm

    ! convert to eV units
    E_i_inp = E_i_inp  / const % eV
    E_n_inp = E_n_inp  / const % eV
    E_dyn_inp = E_dyn_inp  / const % eV
    E_dyn2_inp = E_dyn2_inp  / const % eV
    E_f_inp = E_f_inp / const % eV

    ifile = get_free_handle()
    open(ifile, file="inital_state_eigvec.txt", action='write')
    do i=1,npoints_in
      write(ifile,'(3ES16.6)') X_r(i), c_i(i,1), c_i(i,1) ** 2
    end do
    close(ifile)

    ! sample the positions and momenta
    call sample_x_mom_modes(npoints_x_sampl, npoints_mom_sampl, &
         p % samplemode, X_r, c_i(:,1), x_mom_sampl)
    npoints_x_mom_sampl = size(x_mom_sampl,1)

    delta_t = delta_t * 1.d-15 ! femtoseconds
    time_l = (ntsteps-1) * delta_t

    write(6,*) "outfile", p % outfile
    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
    write(6,*) "mu_SI", mu_SI 
    write(6,*) "time_l", time_l
    write(6,*)
    write(6,*) "Fundamental frequency resolution", 2.0_wp * const % pi * const % hbar /( time_l * const % eV)
    write(6,*) "delta t", delta_t
    write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV), "eV"

    n_omega_in = ntsteps
    n_omega_out = ntsteps
    allocate(omega_in(n_omega_in))
    allocate(omega_out(n_omega_out))
    allocate(F_if_t_omp(nfinal, 1, n_omega_out,3,3))
    allocate(F_if_om_omp(nfinal, 1, n_omega_out,3,3))

    allocate(R_if_om_omp(nfinal, n_omega_in, n_omega_out))
    allocate(R_factor(nfinal, n_omega_in, n_omega_out))

    allocate(lambda_F(nfinal, 1, n_omega_out))
    allocate(lambda_G(nfinal, 1, n_omega_out))
    allocate(lambda_H(nfinal, 1, n_omega_out))
    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))

    !
    ! Loop over trajectories
    !

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do

    !! introduce factor for possible time development backwards in time
    write(6,*) "upper(p % KH_amplitude_mode)", upper(p % KH_amplitude_mode)
    if(upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
      fac_t = 1.0_wp
    else if(upper(p % KH_amplitude_mode) .eq. "INGOING") then
      fac_t = 1.0_wp
    else
      write(6,*) "p % KH_amplitude_mode must be either OUTGOING or INGOING"
    end if
      
    !call hist_init(time_h_0, 1000, X_r(1), X_r(npoints_in) ) 
    !call hist_init(time_h_0_mom, 1000, minval(x_mom_sampl(:,2)), maxval(x_mom_sampl(:,2)) ) 
    !do i=1, ntsteps
    !   call hist_init(time_h(i), 1000, X_r(1), X_r(npoints_in) ) 
    ! end do

    lambda_F = 0.0_wp
    lambda_G = 0.0_wp
    lambda_H = 0.0_wp
    R_factor = 0.0_wp
    
    do traj=1, npoints_x_mom_sampl

      write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl
      
      call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
           E_dyn_inp * const % eV, fac_t * delta_t, mu_SI, x_new, v_new, a_new )
      
      call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
      do i=1,ninter
        call spline_easy(X_r, E_n_inp(i,:), npoints_in, x_new, E_n1(i,:), ntsteps)  
      end do
      
      do i=1,ninter
        do m=1,3
          call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new, D_ni1(i,:,m) , ntsteps)  
        end do
      end do
      
      do i=1,nfinal
        do m=1,3
          call spline_easy(X_r, D_fn_inp(i,1,:,m) , npoints_in, x_new, D_fn1(i,1,:,m) , ntsteps)  
        end do
      end do
      
      do i=1,nfinal
        call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_f1(i,:), ntsteps)  
      end do

      ! first time, compute the mean transition energy, and frequency
      if (traj .eq. 1) then
        ind = minloc(E_i_inp)

        E_ni_mean =  E_n1(nfinal, ind(1)) - E_i1(ind(1))
        write(6,*) "E_ni_mean", E_ni_mean
        
        E_fi_mean =  E_f1(nfinal, ind(1)) - E_i1(ind(1))
        write(6,*) "E_fi_mean", E_fi_mean
        
        E_nf_mean =  E_n1(ninter, ind(1)) - E_f1(nfinal,ind(1))
        write(6,*) "E_nf_mean", E_nf_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
        omega_in = omega_in + E_ni_mean !E_nf_mean + E_fi_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
        omega_out = omega_out + E_nf_mean
        
      end if

      ! compute the resonance factor
      call compute_F_if_om_omp_no_F(E_f1(:,:), E_fi_mean, time, &
           E_i1, gamma_inc, omega_out, E_nf_mean, R_if_om_omp)
      R_factor = R_factor + abs(R_if_om_omp)**2  
      
    end do! do traj=1, npoints_x_mom_sampl
    ! here ends the ground state propagation

    ! loop over each sampling point again, to compute the non-resonant contribution
    do traj2 =1, npoints_x_mom_sampl !ntsteps ! possibly with a stride
        
      call verlet_trajectory_xva(x_mom_sampl(traj2,1), x_mom_sampl(traj2,2)/mu_SI, X_r, &
             E_dyn2_inp * const % eV, fac_t * delta_t, mu_SI, x_new2, v_new2, a_new2 )

        call spline_easy(X_r, E_n_inp, npoints_in, x_new2, E_n2, ntsteps)
        call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
        
        do i=1,nfinal
          call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new2, E_f2(i,:), ntsteps)  
        end do
        
        do i=1,nfinal
          do m=1,3
            call spline_easy(X_r, D_fn_inp(i,1,:,m) , npoints_in, x_new2, D_fn2(i,1,:,m) , ntsteps)  
          end do
        end do

        do i=1,ninter
          do m=1,3
            call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new2, D_ni2(i,:,m) , ntsteps)  
          end do
        end do

        ! non-resonant spectrum, now only one time point needed, hence the one in F_if_t_omp(:,1,:,:,:)
        call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,:,:,:), &
             D_ni2(:, 1, :), time,  F_if_t_omp(:,1,:,:,:), gamma)
        
        !if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
        !
        !  call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,:,:,:), &
        !       D_ni2(:, 1, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)
        !
        !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
        !  ! F_if_t_omp will now store F_if_t_om temporarily
        !  call compute_F_if_om_many_n(E_n2(:,:), E_i2(:), E_ni_mean, D_fn2(:,1,:,:), &
        !       D_ni2(:, :, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)
        !  
        !end if

        ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
        do m1=1,3
          do m2=1,3
            lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m1)) * F_if_t_omp(:,:,:,m2,m2))
            lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m1,m2))
            lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m2,m1))
          end do
        end do
        
      end do ! do traj2 =1, ntsteps2

      ! here ends second trajectory loop
      
     !if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
     !
     !  ! for test purposes, write F_if_t_omp to file
     !  ! call write_F_if_t_omp(p, F_if_t_omp, omega_out, time)
     !  
     !  call compute_F_if_om_omp(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
     !       E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp)
     !
     !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
     !  call compute_F_if_om_omp_ingoing(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
     !       E_i1, gamma_instr, omega_in, E_ni_mean, F_if_om_omp)
        
    !end if
        
        
      !! perform spherical average according to J. Phys. B. 27, 4169 (1994)
      !do m1=1,3
      !  do m2=1,3
      !    lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m1)) * F_if_om_omp(:,:,:,m2,m2))
      !    lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m1,m2))
      !    lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m2,m1))
      !  end do
      !end do

    !end do ! end traj
    
    ! averages according to J. Phys. B. 27, 4169 (1994) 
    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
      !lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
      !lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
      !lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)
      
      lambda_lp =0.0_wp
      lambda_ln =0.0_wp
      lambda_cp =0.0_wp
      
      do om_in =1, n_omega_in
        do f_e =1, nfinal
          lambda_lp(om_in,:) = lambda_lp(om_in,:) + (2.0_wp * lambda_F(f_e,1,:) + 2.0_wp * &
               lambda_G(f_e,1,:)  + 2.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
          lambda_ln(om_in,:) = lambda_ln(om_in,:) + (-1.0_wp * lambda_F(f_e,1,:) + 4.0_wp  &
               * lambda_G(f_e,1,:) -1.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
          lambda_cp(om_in,:) = lambda_cp(om_in,:) + (-2.0_wp * lambda_F(f_e,1,:) + 3.0_wp * &
               lambda_G(f_e,1,:) + 3.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
        end do
      end do

    if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

      write(6,*) "Entering convolute_incoming"
      
      if(gamma_inc .gt. 1d-5) then
        sigma_tmp = lambda_lp
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_lp, 1, "GAUSSIAN")
        sigma_tmp = lambda_ln
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_ln, 1,"GAUSSIAN")
        sigma_tmp = lambda_cp
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_cp, 1,"GAUSSIAN")
      end if

      write(6,*) "Done"
      write(6,*) "Entering convolute_instrumental"
      
      if(gamma_instr .gt. 1d-5) then
        sigma_tmp = lambda_lp
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_lp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_lp, 2,"GAUSSIAN")
        sigma_tmp = lambda_ln
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_ln)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_ln, 2,"GAUSSIAN")
        sigma_tmp = lambda_cp
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_cp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_cp, 2,"GAUSSIAN")
      end if

      write(6,*) "Done"
      
    else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then      
      allocate(sigma_tmp(n_omega_in, n_omega_out))

      if(gamma_inc .gt. 1d-5) then
        sigma_tmp = lambda_lp
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
        sigma_tmp = lambda_ln
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
        sigma_tmp = lambda_cp
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
      end if
      
    end if
    
    ! write spectra to file
    do j=1, n_omega_in
      
      file="_sigma_"
      write(string,'(F6.2)') omega_in(j)   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1, n_omega_out
        
        write(ifile,'(4ES18.10)') omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do
      
      close(ifile) 
      
   end do

 end subroutine calculate_SCKH_res_PES_factor


  
  
  subroutine write_F_if_t_omp(p, F_if_t_omp, omega_out, time)
    use m_precision, only: wp
    use m_constants, only: const
    use m_sckh_params_t, only: sckh_params_t
    use m_io, only: get_free_handle
    
    type(sckh_params_t), intent(in):: p 
    complex(wp), intent(in)::  F_if_t_omp(:,:,:,:,:)
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: time(:)
    
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)

    integer:: nfinal, n_omega_in, n_omega_out, ifile, i, j, m1,m2
    real(wp):: delta_t
    character(80)::  string, file
    
    nfinal = size(F_if_t_omp,1)
    n_omega_in = size(F_if_t_omp,2)
    n_omega_out = size(F_if_t_omp,3)
    delta_t = time(2)-time(1)
    
    allocate(lambda_F(nfinal, n_omega_in, n_omega_out))
    allocate(lambda_G(nfinal, n_omega_in, n_omega_out))
    allocate(lambda_H(nfinal, n_omega_in, n_omega_out))
    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))
    
    lambda_F = 0.0d0
    lambda_G = 0.0d0
    lambda_H = 0.0d0
    
    ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
    do m1=1,3
      do m2=1,3
        lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m1)) * F_if_t_omp(:,:,:,m2,m2))
        lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m1,m2))
        lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m2,m1))
      end do
    end do
    
    ! averages according to J. Phys. B. 27, 4169 (1994) 
    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
    lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
    lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
    lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    ! now write to file
    do j=1, n_omega_in
      
      file="_F_if_t_omp_"
      write(string,'(F6.2)') time(j)/delta_t   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1, n_omega_out
        
        write(ifile,'(4ES18.10)') omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
        
      end do
      
      close(ifile) 
      
    end do

  end subroutine write_F_if_t_omp
  
end module m_SCKH_resonant_PES
