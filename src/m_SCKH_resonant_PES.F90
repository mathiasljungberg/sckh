module m_SCKH_resonant_PES
  implicit none

contains

  subroutine calculate_SCKH_res_PES(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_SCKH_utils, only: compute_F_if_omp_many_n
    use m_SCKH_utils, only: compute_F_if_omp_sum_n
    use m_SCKH_utils, only: compute_F_if_omp_one_n
    use m_SCKH_utils, only: compute_F_if_om_omp
    use m_SCKH_utils, only: compute_F_if_om_omp_one_f
    use m_SCKH_utils, only: compute_F_if_om_omp_no_F
    use m_SCKH_utils, only: compute_F_if_om_omp_no_F_one
    use m_sckh_params_t, only: sckh_params_t 
    use m_io, only: get_free_handle
    use m_splines, only: spline_easy
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: get_projections
    use m_PES_io, only: read_file_list
    use m_KH_functions, only: solve_vib_problem
    use m_fftw3, only: get_omega_reordered_fftw
    use m_upper, only : upper
    use m_KH_utils, only: convolute_incoming
    use m_KH_utils, only: convolute_instrumental
    use m_spectrum_utils, only: convolution_lorentzian_grid_fft_many_freq
    
    type(sckh_params_t), intent(inout):: p 

    real(wp), allocatable::  time(:)
    real(wp), allocatable:: E_dyn_inp(:)
    real(wp), allocatable:: E_dyn2_inp(:)  
    real(wp), allocatable:: E_n0(:)
    real(wp), allocatable:: E_fn_corr(:,:)
    real(wp), allocatable:: E_lp_corr(:)
    real(wp), allocatable:: shift(:)
    real(wp), allocatable:: E_i_inp(:)
    real(wp), allocatable:: E_i1(:)
    real(wp), allocatable:: E_i2(:)
    real(wp), allocatable:: E_n_inp(:,:)
    real(wp), allocatable:: E_n1(:,:)
    real(wp), allocatable:: E_n2(:,:)
    real(wp), allocatable:: E_f_inp(:,:)
    real(wp), allocatable:: E_f1(:,:)
    real(wp), allocatable:: E_f2(:,:)
    real(wp), allocatable:: E_fc1(:,:,:)
    real(wp), allocatable:: E_fc2(:,:,:)
    real(wp), allocatable:: D_fn_inp(:,:,:,:)
    real(wp), allocatable:: D_fn1(:,:,:,:)
    real(wp), allocatable:: D_fn2(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:)
    real(wp), allocatable:: D_ni1(:,:,:)
    real(wp), allocatable:: D_ni2(:,:,:) 
    real(wp):: E_nf_mean
    real(wp):: E_fi_mean
    real(wp):: E_ni_mean
    
    character(80)::  string
    character(80)::  file
    real(wp), allocatable:: x_sampl(:)
    real(wp), allocatable:: mom_sampl(:)
    real(wp), allocatable:: x_mom_sampl(:,:)
    real(wp), allocatable:: x_new(:)
    real(wp), allocatable:: v_new(:)
    real(wp), allocatable:: a_new(:)
    
    real(wp), allocatable:: x_new2(:)
    real(wp), allocatable:: v_new2(:)
    real(wp), allocatable:: a_new2(:)

    integer:: npoints_x_sampl
    integer:: npoints_mom_sampl
    integer:: npoints_x_mom_sampl
    integer:: ninter
    integer:: nfinal
    integer:: ntsteps  
    integer:: npoints_in
    integer:: nfinal_tot

    real(wp), allocatable:: omega_in(:)
    real(wp), allocatable:: omega_out(:)
    integer:: n_omega_in
    integer:: n_omega_out

    real(wp), allocatable:: c_i(:,:)
    real(wp), allocatable:: eig_i(:)
    
    real(wp):: gamma
    real(wp):: gamma_inc
    real(wp):: gamma_instr
    real(wp):: time_l
    real(wp):: delta_t
    real(wp):: norm
    real(wp):: mu_SI
    real(wp):: dx
    real(wp):: fac_t
    
    integer, dimension(1)::ind  
    integer:: ifile
    
    real(wp), allocatable::  X_r(:)
    real(wp):: dvr_start 
    integer:: npoints, ii
    real(8):: dnrm2
    integer:: i
    integer:: j
    integer:: m
    integer:: m1
    integer:: m2
    integer:: m3
    integer:: m4
    integer:: traj
    integer:: traj2
    integer:: n_e
    integer:: f_e
    integer:: fn_e
    integer:: om_in

    complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
    complex(wp), allocatable::  F_if_om_omp(:,:,:,:)

    complex(wp), allocatable::  F_if_omp(:,:,:,:)    
    complex(wp), allocatable::  F_if_omp_tmp(:,:,:,:)    

    complex(wp), allocatable::  R_if_om_omp(:,:)    
    real(wp), allocatable:: R_factor(:,:)

    !complex(wp), allocatable::  R_if_om_omp(:,:,:)    
    !real(wp), allocatable:: R_factor(:,:,:)

    real(wp), allocatable:: lambda_F(:,:)
    real(wp), allocatable:: lambda_G(:,:)
    real(wp), allocatable:: lambda_H(:,:)

    real(wp), allocatable:: lambda_F_tmp(:,:)
    real(wp), allocatable:: lambda_G_tmp(:,:)
    real(wp), allocatable:: lambda_H_tmp(:,:)

    real(wp), allocatable:: lambda_lp(:,:)
    real(wp), allocatable:: lambda_ln(:,:)
    real(wp), allocatable:: lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)
    
    real(wp), allocatable::  F2_tensor(:,:,:,:,:,:)
    
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

    allocate( E_i_inp(npoints_in))
    allocate(E_n_inp(ninter,npoints_in))
    allocate(E_dyn_inp(npoints_in))
    allocate(E_dyn2_inp(npoints_in))
    allocate(E_f_inp(nfinal,npoints_in))
    allocate(E_n0(npoints_in))
    allocate(E_lp_corr(npoints_in))
    allocate(D_ni_inp(ninter, npoints_in,3))
    allocate(D_fn_inp(nfinal, ninter, npoints_in,3))
    allocate(time(ntsteps))
    allocate(E_i1(ntsteps))
    allocate(E_i2(ntsteps))
    allocate(E_n1(ninter, ntsteps))
    allocate(E_n2(ninter, ntsteps))
    allocate(D_ni1(ninter, ntsteps,3))
    allocate(D_ni2(ninter, ntsteps,3))
    allocate(D_fn2(nfinal, ninter, ntsteps,3))
    allocate(D_fn1(nfinal, ninter, ntsteps,3))
    !allocate(E_i_mean(ntsteps))
    allocate(shift(npoints_in))
    allocate(c_i(npoints_in,npoints_in))
    allocate(eig_i(npoints_in))
    allocate(x_new(ntsteps))
    allocate(v_new(ntsteps))
    allocate(a_new(ntsteps))
    allocate(x_new2(ntsteps))
    allocate(v_new2(ntsteps))
    allocate(a_new2(ntsteps))
    allocate(X_r(npoints_in))
    
    ! in ORB mode there are nfinal * ninter final states
    write(6,*) "p % KH_states_mode = ", p % KH_states_mode
    if (upper(p % KH_states_mode) .eq. "STATES") then
      nfinal_tot = nfinal

      allocate(E_fc1(nfinal,1, ntsteps), &
           E_fc2(nfinal,1,ntsteps))
      
    else if (upper(p % KH_states_mode) .eq. "ORBS") then
      nfinal_tot = nfinal * ninter

      allocate(E_fc1(nfinal,ninter, ntsteps), &
           E_fc2(nfinal,ninter, ntsteps))
    else
      write(6,*) "p % KH_states_mode should be either 'STATES' or 'ORBS'"
    end if

    
    ! set up grid points
    do i = 1, npoints_in
      X_r(i) = (i-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_r, E_i_inp)

    ! intermediate state reference energy (lowest state) 
    if(p % use_n0_state) then 
      call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n0)
    else
      E_n0 =0.0d0
    end if

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

    ! ugly hack again to play
    !E_dyn_inp = (E_i_inp +  E_f_inp(1,:))/2.0d0
    !E_dyn2_inp = (E_i_inp + E_n_inp(1,:) +  E_f_inp(1,:))/3.0d0

    ! read list of corrections to the final state files coming from the excited electron
    if (upper(p % KH_states_mode) .eq. "ORBS") then    
      allocate( E_fn_corr(p % npesfile_n, p % npoints_in))

      call read_file_list(p % pes_file_list_fn_corr, p % npesfile_n, p % pes_files_fn_corr)
      
      do j=1,p % npesfile_n
        call read_PES_file(p % pes_files_fn_corr(j), p % npoints_in, &
             p % npoints_in, X_r, E_fn_corr(j,:))
      end do
    else
      
    end if

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
    E_n0 = E_n0  / const % eV
    E_dyn_inp = E_dyn_inp  / const % eV
    E_dyn2_inp = E_dyn2_inp  / const % eV
    E_f_inp = E_f_inp / const % eV
    if(allocated( E_fn_corr))  E_fn_corr =  E_fn_corr / const % eV
    
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

!    allocate(F_if_t_omp(nfinal_tot, n_omega_in, n_omega_out,3,3),&
!         F_if_om_omp(nfinal_tot, n_omega_in, n_omega_out,3,3),&
!         lambda_F(nfinal_tot, n_omega_in, n_omega_out),&
!         lambda_G(nfinal_tot, n_omega_in, n_omega_out),&
!         lambda_H(nfinal_tot, n_omega_in, n_omega_out),&
!         lambda_lp(n_omega_in, n_omega_out),&
!         lambda_ln(n_omega_in, n_omega_out),&
!         lambda_cp(n_omega_in, n_omega_out),&
!         sigma_tmp(n_omega_in, n_omega_out))

    if (upper(p % runmode_sckh_res) .eq. "FACTOR") then 
        allocate(F_if_omp(nfinal_tot, n_omega_out,3,3))
        allocate(F_if_omp_tmp(nfinal_tot, n_omega_out,3,3))
        
        allocate(lambda_F_tmp(nfinal_tot, n_omega_out))
        allocate(lambda_G_tmp(nfinal_tot, n_omega_out))
        allocate(lambda_H_tmp(nfinal_tot, n_omega_out))

      else if (upper(p % runmode_sckh_res) .eq. "FACTOR_TRAJ") then
      write(6,*) "p % runmode_sckh_res = FACTOR_TRAJ"
      allocate(F2_tensor(n_omega_in, n_omega_out,3,3,3,3))        
      F2_tensor = 0.0_wp
    else
      allocate(F_if_t_omp(n_omega_in, n_omega_out,3,3))
      allocate(F_if_om_omp(n_omega_in, n_omega_out,3,3))
      allocate(F2_tensor(n_omega_in, n_omega_out,3,3,3,3))        
      F2_tensor = 0.0_wp
    end if

    allocate(lambda_F(n_omega_in, n_omega_out))
    allocate(lambda_G(n_omega_in, n_omega_out))
    allocate(lambda_H(n_omega_in, n_omega_out))

    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))
    allocate(sigma_tmp(n_omega_in, n_omega_out))

    
    !
    ! Loop over trajectories
    !

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do

    lambda_F = 0.0_wp
    lambda_G = 0.0_wp
    lambda_H = 0.0_wp
    
    if (upper(p % runmode_sckh_res) .eq. "FACTOR") then

      F_if_omp = 0.0_wp
      lambda_F_tmp = 0.0_wp
      lambda_G_tmp = 0.0_wp
      lambda_H_tmp = 0.0_wp
      
      allocate(R_if_om_omp(ntsteps, ntsteps))
      allocate(R_factor(ntsteps, ntsteps))

      ! compute mean energies on E_dyn_inp that will match with the other routines
      call verlet_trajectory_xva(x_mom_sampl(1,1), x_mom_sampl(1,2)/mu_SI, X_r, &
           E_dyn_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )
      
      call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      call spline_easy(X_r, E_n_inp(1,:) + E_n0, npoints_in, x_new, E_n1(1,:), ntsteps)
      call spline_easy(X_r, E_f_inp(1,:) + E_fn_corr(1,:), npoints_in, x_new, E_fc1(1,1,:), ntsteps)  
        
      
      if (upper(p % KH_states_mode) .eq. "STATES") then
      
                ! run nonresoannt spectrum on E_dyn2_inp
        do traj=1, npoints_x_mom_sampl
          
          write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl
          call compute_traj_and_spline(p, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, &
               X_r, E_dyn2_inp, delta_t, mu_SI, &
               x_new, v_new, a_new, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
               E_i1, E_n1, E_fc1, D_ni1, D_fn1 )

          !if (traj .eq. 1) then
          !  call compute_E_means_and_omegas(E_i_inp, E_i1, E_n1, E_fc1, &
          !       E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)
          !end if
          
          call compute_F_if_omp_many_n(E_n1(:,:), E_fc1(:,1,:), E_nf_mean, D_fn1(:,1,:,:), &
               D_ni1(:, 1, :), time, F_if_omp(:,:,:,:), gamma)

        ! sum over trajectories
        !F_if_omp(:,:,:,:) = F_if_omp(:,:,:,:) + F_if_omp_tmp(:,:,:,:)

        ! temporary lambdas (f_e, omega_out), sum goes over square 
        do m1=1,3
          do m2=1,3
              lambda_F_tmp(:, :) =  lambda_F_tmp(:, :) + real(conjg(F_if_omp(:,:,m1,m1)) * &
                   F_if_omp(:,:,m2,m2)) 
              lambda_G_tmp(:, :) =  lambda_G_tmp(:, :) + real(conjg(F_if_omp(:,:,m1,m2)) * &
                   F_if_omp(:,:,m1,m2)) 
              lambda_H_tmp(:, :) =  lambda_H_tmp(:, :) + real(conjg(F_if_omp(:,:,m1,m2)) * &
                   F_if_omp(:,:,m2,m1)) 
          end do
        end do

      end do !do traj=1, npoints_x_mom_sampl

      ! compute second set of trajectories on E_dyn_inp
      do f_e =1, nfinal
        do traj=1, npoints_x_mom_sampl

        write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl
        call compute_traj_and_spline(p, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, &
             X_r, E_dyn_inp, delta_t, mu_SI, &
             x_new, v_new, a_new, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
             E_i1, E_n1, E_fc1, D_ni1, D_fn1 )
        
        ! compute R-factor
        call compute_F_if_om_omp_no_F_one(E_fc1(f_e,1,:), E_fi_mean, time, &
             E_i1, gamma_inc, omega_out, E_nf_mean, R_if_om_omp)
        !call compute_F_if_om_omp_no_F(E_fc1(:,1,:), E_fi_mean, time, &
        !     E_i1, gamma_inc, omega_out, E_nf_mean, R_if_om_omp)
        
        R_factor = (abs(R_if_om_omp))**2
        
        
        
        ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
        do om_in =1, ntsteps
          lambda_F(om_in, :) =  lambda_F(om_in, :) + lambda_F_tmp(f_e,:) * R_factor( om_in, :)
          lambda_G(om_in, :) =  lambda_G(om_in, :) + lambda_G_tmp(f_e,:) * R_factor( om_in, :)
          lambda_H(om_in, :) =  lambda_H(om_in, :) + lambda_H_tmp(f_e,:) * R_factor( om_in, :)
        end do

      end do

    end do !do traj=1, npoints_x_mom_sampl

  else if (upper(p % KH_states_mode) .eq. "ORBS") then

    call compute_E_means_and_omegas_one(E_i_inp, E_i1, E_n1(1,:), E_fc1(1,1,:), &
         E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)
   
    ! run nonresonant spectrum on E_dyn2_inp
    do traj=1, npoints_x_mom_sampl

        !write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl
        call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
             E_dyn2_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )

        do f_e =1, nfinal
          do n_e =1, ninter

            fn_e = ninter * (f_e -1) + n_e  
            
            call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
            call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new, E_n1(n_e,:), ntsteps)
            call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new, E_fc1(f_e,n_e,:), ntsteps)  
            
            do m=1,3
              call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new, D_fn1(f_e,1,:,m) , ntsteps)  
              call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new, D_ni1(n_e,:,m) , ntsteps)  
            end do
            
            call compute_F_if_omp_one_n(E_n1(n_e,:), E_fc1(f_e,n_e,:), E_nf_mean, D_fn1(f_e,1,:,:), &
                 D_ni1(n_e, 1, :), time, F_if_omp(1,:,:,:), gamma)

            ! temporary lambdas (f_e, omega_out), sum goes over square 
            do m1=1,3
              do m2=1,3
                lambda_F_tmp(fn_e, :) =  lambda_F_tmp(fn_e, :) + real(conjg(F_if_omp(1,:,m1,m1)) * &
                     F_if_omp(1,:,m2,m2)) 
                lambda_G_tmp(fn_e, :) =  lambda_G_tmp(fn_e, :) + real(conjg(F_if_omp(1,:,m1,m2)) * &
                     F_if_omp(1,:,m1,m2)) 
                lambda_H_tmp(fn_e, :) =  lambda_H_tmp(fn_e, :) + real(conjg(F_if_omp(1,:,m1,m2)) * &
                     F_if_omp(1,:,m2,m1)) 
              end do
            end do

          end do ! do n_e =1, ninter
        end do ! do f_e =1, nfinal
        
      end do !do traj=1, npoints_x_mom_sampl

      ! compute second set of trajectories on E_dyn_inp
      do f_e =1, nfinal
        do n_e =1, ninter
          fn_e = ninter * (f_e -1) + n_e  

          do traj=1, npoints_x_mom_sampl
            
            !write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl
            !call compute_traj_and_spline(p, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, &
            !     X_r, E_dyn_inp, delta_t, mu_SI, &
            !     x_new, v_new, a_new, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
            !     E_i1, E_n1, E_fc1, D_ni1, D_fn1 )
            
            call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
                 E_dyn2_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )
            
            call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
            call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new, E_n1(n_e,:), ntsteps)
            call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new, E_fc1(f_e,n_e,:), ntsteps)  
            
            do m=1,3
              call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new, D_fn1(f_e,1,:,m) , ntsteps)  
              call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new, D_ni1(n_e,:,m) , ntsteps)  
            end do
            
            ! compute R-factor
            call compute_F_if_om_omp_no_F_one(E_fc1(f_e,n_e,:), E_fi_mean, time, &
                 E_i1, gamma_inc, omega_out, E_nf_mean, R_if_om_omp)
            !call compute_F_if_om_omp_no_F(E_fc1(:,1,:), E_fi_mean, time, &
            !     E_i1, gamma_inc, omega_out, E_nf_mean, R_if_om_omp)
            
            R_factor = (abs(R_if_om_omp))**2
            
            
            ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
            do om_in =1, ntsteps
              lambda_F(om_in, :) =  lambda_F(om_in, :) + lambda_F_tmp(fn_e,:) * R_factor( om_in, :)
              lambda_G(om_in, :) =  lambda_G(om_in, :) + lambda_G_tmp(fn_e,:) * R_factor( om_in, :)
              lambda_H(om_in, :) =  lambda_H(om_in, :) + lambda_H_tmp(fn_e,:) * R_factor( om_in, :)
            end do
            
          end do !do traj=1, npoints_x_mom_sampl
          
        end do ! do n_e =1, ninter
      end do! do f_e =1, nfinal
      
    end if
    
  else
    
    do traj=1, npoints_x_mom_sampl
      
      write(6,*) "Computing trajectory ", traj, "out of", npoints_x_mom_sampl

      
      ! option to do separate dynamcis of D_nf and D_in
      if(.true.) then
        call compute_traj_and_spline(p, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, &
             X_r, E_dyn_inp, delta_t, mu_SI, &
             x_new, v_new, a_new, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
             E_i1, E_n1, E_fc1, D_ni1, D_fn1 )
        
      else
        ! this is dynamics on E_dyn_inp
        call compute_traj_and_spline(p, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, &
             X_r, E_dyn_inp, delta_t, mu_SI, &
             x_new, v_new, a_new, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
             E_i1, E_n1, E_fc1, D_ni1, D_fn1 )

        ! this is dynamics on E_dyn2_inp
        call compute_traj_and_spline(p, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, &
             X_r, E_dyn2_inp, delta_t, mu_SI, &
             x_new2, v_new2, a_new2, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
             E_i2, E_n2, E_fc2, D_ni2, D_fn2 )


      end if

        
      ! compute mean energies on E_dyn_inp that will match with the other routines
      if (traj .eq. 1) then
        call compute_E_means_and_omegas(E_i_inp, E_i1, E_n1, E_fc1, &
             E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)
      end if

      if(upper(p % runmode_sckh_res) .eq. "FULL") then
        
      if(.true.) then
        ! FULL case, run trajectories with all points in x_new as starting points 
        call compute_F_SCKH_res_PES_full(p, x_new, v_new, X_r, E_dyn2_inp, delta_t, mu_SI, &
             E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
             E_i1, E_fc1, gamma_inc, omega_out, F2_tensor) 
      else if (.true.) then
        call compute_F_SCKH_res_PES_full_separate(p, x_new2, v_new2, X_r, E_dyn2_inp, delta_t, mu_SI, &
             E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
             E_i1, E_n1, E_fc1, D_ni1, gamma_inc, omega_out, omega_in, F2_tensor) 
      else
        call compute_F_SCKH_res_PES_full_separate_ingoing(p, x_new2, v_new2, X_r, E_dyn2_inp, delta_t, mu_SI, &
             E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean, E_ni_mean,time, gamma, &
             E_i1, E_n1, E_fc1, D_ni1, D_fn1, gamma_inc, omega_out, omega_in, F2_tensor) 
      end if


    else if(upper(p % runmode_sckh_res) .eq. "FACTOR_TRAJ") then

        ! use same initial conditions x_mom_sampl to run another trajectory to compute R(omega -omega')
        call compute_F_SCKH_res_PES_factor_each_traj(p, x_mom_sampl, traj, X_r, E_dyn2_inp, delta_t, mu_SI, &
             E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
             E_i1, E_fc1, gamma_inc, omega_out, F2_tensor)
        
      end if
        
        ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
        do m1=1,3
          do m2=1,3
            lambda_F(:, :) =  lambda_F(:, :) + F2_tensor(:,:,m1,m1,m2,m2) !  lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
            lambda_G(:, :) =  lambda_G(:, :) + F2_tensor(:,:,m1,m2,m1,m2)  !  lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
            lambda_H(:, :) = lambda_H(:, :) + F2_tensor(:,:,m1,m2,m2,m1)  !lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
          end do
        end do
      
    end do ! end traj      

  end if
    
    ! averages according to J. Phys. B. 27, 4169 (1994) 
    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
    lambda_lp = 2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H
    lambda_ln = -1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H
    lambda_cp = -2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H

    write(6,*) "Entering convolute_incoming, broadening ",  upper(p % broadening_func_inc)
    
    if(gamma_inc .gt. 1d-5) then
      sigma_tmp = lambda_lp
      !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_lp, 1,  upper(p % broadening_func_inc))
      sigma_tmp = lambda_ln
      !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(p % broadening_func_inc))
      sigma_tmp = lambda_cp
      !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_cp, 1, upper(p % broadening_func_inc))
    end if
    
    write(6,*) "Done"
    !write(6,*) "Entering convolute_instrumental"
    write(6,*) "Entering convolute_instrumental, broadening  ",  upper(p % broadening_func_instr) 
    
    if(gamma_instr .gt. 1d-5) then
      sigma_tmp = lambda_lp
      !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_lp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_ln
      !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_ln)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_cp
      !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_cp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(p % broadening_func_instr))
    end if
    
    write(6,*) "Done"
    
    ! write spectra to individual files
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

    ! write spectra to file
   file="_sigma_all"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in
     do i=1, n_omega_out
       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
   end do
   
   close(ifile) 


   ! write spectra to file
   file="_sigma_all_nogrid"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in
     do i=1, n_omega_out
       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
     write(ifile, *) 
   end do
   
   close(ifile) 
   
 end subroutine calculate_SCKH_res_PES


  
  subroutine calculate_SCKH_res_PES_factor_each_traj(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    !use m_SCKH_utils, only: verlet_trajectory
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_SCKH_utils, only: compute_F_if_omp_many_n
    use m_SCKH_utils, only: compute_F_if_om_many_n
    !use m_SCKH_utils, only: compute_F_if_om_omp
    use m_SCKH_utils, only: compute_F_if_om_omp_no_F
    use m_SCKH_utils, only: compute_F_if_om_omp_ingoing_no_F
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
    
    real(wp):: gamma, gamma_inc, gamma_instr, gamma_R
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
    integer::i,j,m,m1,m2, traj, traj2, om_in, om_out, f_e

    complex(wp), allocatable::  F_if_t_omp(:,:,:,:,:)    
    complex(wp), allocatable::  F_if_om_omp(:,:,:,:,:)
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)


    complex(wp), allocatable::  R_if_om_omp(:,:,:)    
    real(wp), allocatable:: R_factor(:,:,:)
    real(wp):: threshold
    
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
    gamma_R = p % gamma_R_FWHM / 2 

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
         a_new2(ntsteps)&
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

    ! XXX dangerous hack!
    !E_dyn_inp = (E_i_inp + E_f_inp(1,:)) /2.0_wp
    
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

    allocate(sigma_tmp(n_omega_in, n_omega_out))

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

!    lambda_F = 0.0_wp
!    lambda_G = 0.0_wp
!    lambda_H = 0.0_wp
    R_factor = 0.0_wp

    lambda_lp =0.0_wp
    lambda_ln =0.0_wp
    lambda_cp =0.0_wp


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
      if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
        call compute_F_if_om_omp_no_F(E_f1(:,:), E_fi_mean, time, &
             E_i1, gamma_R, omega_out, E_nf_mean, R_if_om_omp)

      else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then

        call compute_F_if_om_omp_ingoing_no_F(E_f1(:,:), E_fi_mean, time, &
             E_i1, gamma_R, omega_in, E_ni_mean, R_if_om_omp)
      end if
      
      !R_factor = R_factor + (abs(R_if_om_omp))**2
      !R_factor = R_factor + (aimag(R_if_om_omp))**2
      R_factor = (abs(R_if_om_omp))**2

      
!!!    end do! do traj=1, npoints_x_mom_sampl

    !! broaden R_factor
    !write(6,*) "broadening R-factor"
    !do f_e =1, nfinal
    !  sigma_tmp = R_factor(f_E,:,:)
    !  call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
    !       2.0_wp * gamma_R, omega_out, R_factor(f_e,:,:), 1, "LORENTZIAN", -10.0_wp)
    !end do
    !write(6,*) "Done"
    
    !! hack to restrict R_factor
    !do f_e =1, nfinal
    !  do om_in =1, n_omega_in
    !    threshold = maxval(R_factor(f_e, om_in,:)) * 0.001
    !    do om_out=1, n_omega_out
    !      if(R_factor(f_e, om_in, om_out) .lt. threshold) then
    !        R_factor(f_e, om_in, om_out) =0.0_wp
    !      end if
    !    end do
    !  end do
    !end do
          
    ! here ends the ground state propagation

    ! loop over each sampling point again, to compute the non-resonant contribution
!!!    do traj2 =1, npoints_x_mom_sampl !ntsteps ! possibly with a stride
        
      call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
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
        if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
          call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,1,:,:), &
               D_ni2(:, 1, :), time,  F_if_t_omp(:,1,:,:,:), gamma)
          
        else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
          ! F_if_t_omp will now store F_if_t_om temporarily
          write(6,*) "Ingoing, traj", traj
          call compute_F_if_om_many_n(E_n2(:,:), E_i2(:), E_ni_mean, D_fn2(:,1,:,:), &
               D_ni2(:, :, :), time,  F_if_t_omp(:,1,:,:,:), gamma)
          
        end if



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
        lambda_F = 0.0_wp
        lambda_G = 0.0_wp
        lambda_H = 0.0_wp

        do m1=1,3
          do m2=1,3
!            lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m1)) * F_if_t_omp(:,:,:,m2,m2))
!            lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m1,m2))
!            lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m2,m1))
            lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m1)) * F_if_t_omp(:,:,:,m2,m2))
            lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m1,m2))
            lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_t_omp(:,:,:,m1,m2)) * F_if_t_omp(:,:,:,m2,m1))
          end do
        end do
        
!      end do ! do traj2 =1, ntsteps2

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
      
!      lambda_lp =0.0_wp
!      lambda_ln =0.0_wp
!      lambda_cp =0.0_wp

      if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

        do om_in =1, n_omega_in
        do f_e =1, nfinal
          lambda_lp(om_in,:) = lambda_lp(om_in,:) + (2.0_wp * lambda_F(f_e,1,:) + 2.0_wp * &
               lambda_G(f_e,1,:)  + 2.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
               !(0.1_wp * maxval(R_factor(f_e,om_in, :)) + (1.0_wp -0.1_wp)*R_factor(f_e,om_in, :)) ! superdangerous hack!!
          lambda_ln(om_in,:) = lambda_ln(om_in,:) + (-1.0_wp * lambda_F(f_e,1,:) + 4.0_wp  &
               * lambda_G(f_e,1,:) -1.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
          lambda_cp(om_in,:) = lambda_cp(om_in,:) + (-2.0_wp * lambda_F(f_e,1,:) + 3.0_wp * &
               lambda_G(f_e,1,:) + 3.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
        end do
      end do

    else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then

      do om_in =1, n_omega_in
        do f_e =1, nfinal
          lambda_lp(om_in,:) = lambda_lp(om_in,:) + (2.0_wp * lambda_F(f_e,1,om_in) + 2.0_wp * &
               lambda_G(f_e,1,om_in)  + 2.0_wp * lambda_H(f_e,1,om_in)) * R_factor(f_e,om_in, :)
          lambda_ln(om_in,:) = lambda_ln(om_in,:) + (-1.0_wp * lambda_F(f_e,1,om_in) + 4.0_wp  &
               * lambda_G(f_e,1,om_in) -1.0_wp * lambda_H(f_e,1,om_in)) * R_factor(f_e,om_in, :)
          lambda_cp(om_in,:) = lambda_cp(om_in,:) + (-2.0_wp * lambda_F(f_e,1,om_in) + 3.0_wp * &
               lambda_G(f_e,1,om_in) + 3.0_wp * lambda_H(f_e,1,om_in)) * R_factor(f_e,om_in, :)
        end do
      end do

      
    end if

  end do ! do traj
    
      write(6,*) "Entering convolute_incoming, broadening ",  upper(p % broadening_func_inc) 
      
      if(gamma_inc .gt. 1d-5) then
        sigma_tmp = lambda_lp
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_lp, 1, upper(p % broadening_func_inc))
        sigma_tmp = lambda_ln
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(p % broadening_func_inc))
        sigma_tmp = lambda_cp
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_cp, 1, upper(p % broadening_func_inc))
      end if

      write(6,*) "Done"
      write(6,*) "Entering convolute_instrumental, broadening  ",  upper(p % broadening_func_instr) 
      
      if(gamma_instr .gt. 1d-5) then
        sigma_tmp = lambda_lp
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_lp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(p % broadening_func_instr))
        sigma_tmp = lambda_ln
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_ln)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(p % broadening_func_instr))
        sigma_tmp = lambda_cp
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_cp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(p % broadening_func_instr))
      end if

      write(6,*) "Done"
      
    !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then      
    !  allocate(sigma_tmp(n_omega_in, n_omega_out))
    !
    !  if(gamma_inc .gt. 1d-5) then
    !    sigma_tmp = lambda_lp
    !    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
    !    sigma_tmp = lambda_ln
    !    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
    !    sigma_tmp = lambda_cp
    !    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
    !  end if
    !  
    !end if
    
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

 end subroutine calculate_SCKH_res_PES_factor_each_traj


 ! complete factorization, trajs separately
  subroutine calculate_SCKH_res_PES_factor(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    !use m_SCKH_utils, only: verlet_trajectory
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_SCKH_utils, only: compute_F_if_omp_many_n
    use m_SCKH_utils, only: compute_F_if_om_many_n
    !use m_SCKH_utils, only: compute_F_if_om_omp
    use m_SCKH_utils, only: compute_F_if_om_omp_no_F
    use m_SCKH_utils, only: compute_F_if_om_omp_ingoing_no_F
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
    
    real(wp):: gamma, gamma_inc, gamma_instr, gamma_R
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
    integer::i,j,m,m1,m2, traj, traj2, om_in, om_out, f_e

    complex(wp), allocatable::  F_if_t_omp(:,:,:,:,:)    
    complex(wp), allocatable::  F_if_om_omp(:,:,:,:,:)
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)


    complex(wp), allocatable::  R_if_om_omp(:,:,:)    
    real(wp), allocatable:: R_factor(:,:,:)
    real(wp):: threshold
    
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
    gamma_R = p % gamma_R_FWHM / 2 

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
         a_new2(ntsteps)&
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

    ! XXX dangerous hack!
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

    allocate(sigma_tmp(n_omega_in, n_omega_out))

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
      if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
        call compute_F_if_om_omp_no_F(E_f1(:,:), E_fi_mean, time, &
             E_i1, gamma_R, omega_out, E_nf_mean, R_if_om_omp)

      else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then

        call compute_F_if_om_omp_ingoing_no_F(E_f1(:,:), E_fi_mean, time, &
             E_i1, gamma_R, omega_in, E_ni_mean, R_if_om_omp)
      end if
      
      R_factor = R_factor + (abs(R_if_om_omp))**2
      !R_factor = R_factor + (aimag(R_if_om_omp))**2
      
    end do! do traj=1, npoints_x_mom_sampl

    !! broaden R_factor
    !write(6,*) "broadening R-factor"
    !do f_e =1, nfinal
    !  sigma_tmp = R_factor(f_E,:,:)
    !  call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
    !       2.0_wp * gamma_R, omega_out, R_factor(f_e,:,:), 1, "LORENTZIAN", -10.0_wp)
    !end do
    !write(6,*) "Done"
    
    !! hack to restrict R_factor
    !do f_e =1, nfinal
    !  do om_in =1, n_omega_in
    !    threshold = maxval(R_factor(f_e, om_in,:)) * 0.001
    !    do om_out=1, n_omega_out
    !      if(R_factor(f_e, om_in, om_out) .lt. threshold) then
    !        R_factor(f_e, om_in, om_out) =0.0_wp
    !      end if
    !    end do
    !  end do
    !end do
          
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
        if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
          call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,1,:,:), &
               D_ni2(:, 1, :), time,  F_if_t_omp(:,1,:,:,:), gamma)
          
        else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
          ! F_if_t_omp will now store F_if_t_om temporarily
          write(6,*) "Ingoing, traj2", traj2
          call compute_F_if_om_many_n(E_n2(:,:), E_i2(:), E_ni_mean, D_fn2(:,1,:,:), &
               D_ni2(:, :, :), time,  F_if_t_omp(:,1,:,:,:), gamma)
          
        end if



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

      if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

        do om_in =1, n_omega_in
        do f_e =1, nfinal
          lambda_lp(om_in,:) = lambda_lp(om_in,:) + (2.0_wp * lambda_F(f_e,1,:) + 2.0_wp * &
               lambda_G(f_e,1,:)  + 2.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
               !(0.1_wp * maxval(R_factor(f_e,om_in, :)) + (1.0_wp -0.1_wp)*R_factor(f_e,om_in, :)) ! superdangerous hack!!
          lambda_ln(om_in,:) = lambda_ln(om_in,:) + (-1.0_wp * lambda_F(f_e,1,:) + 4.0_wp  &
               * lambda_G(f_e,1,:) -1.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
          lambda_cp(om_in,:) = lambda_cp(om_in,:) + (-2.0_wp * lambda_F(f_e,1,:) + 3.0_wp * &
               lambda_G(f_e,1,:) + 3.0_wp * lambda_H(f_e,1,:)) * R_factor(f_e,om_in, :)
        end do
      end do

    else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then

      do om_in =1, n_omega_in
        do f_e =1, nfinal
          lambda_lp(om_in,:) = lambda_lp(om_in,:) + (2.0_wp * lambda_F(f_e,1,om_in) + 2.0_wp * &
               lambda_G(f_e,1,om_in)  + 2.0_wp * lambda_H(f_e,1,om_in)) * R_factor(f_e,om_in, :)
          lambda_ln(om_in,:) = lambda_ln(om_in,:) + (-1.0_wp * lambda_F(f_e,1,om_in) + 4.0_wp  &
               * lambda_G(f_e,1,om_in) -1.0_wp * lambda_H(f_e,1,om_in)) * R_factor(f_e,om_in, :)
          lambda_cp(om_in,:) = lambda_cp(om_in,:) + (-2.0_wp * lambda_F(f_e,1,om_in) + 3.0_wp * &
               lambda_G(f_e,1,om_in) + 3.0_wp * lambda_H(f_e,1,om_in)) * R_factor(f_e,om_in, :)
        end do
      end do

      
    end if
      
      write(6,*) "Entering convolute_incoming"
      
      if(gamma_inc .gt. 1d-5) then
        sigma_tmp = lambda_lp
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_lp, 1, upper(p % broadening_func_inc))
        sigma_tmp = lambda_ln
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(p % broadening_func_inc))
        sigma_tmp = lambda_cp
        !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_inc, omega_out, lambda_cp, 1,upper(p % broadening_func_inc))
      end if

      write(6,*) "Done"
      write(6,*) "Entering convolute_instrumental"
      
      if(gamma_instr .gt. 1d-5) then
        sigma_tmp = lambda_lp
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_lp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(p % broadening_func_instr))
        sigma_tmp = lambda_ln
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_ln)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(p % broadening_func_instr))
        sigma_tmp = lambda_cp
        !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_cp)
        call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
             2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(p % broadening_func_instr))
      end if

      write(6,*) "Done"
      
    !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then      
    !  allocate(sigma_tmp(n_omega_in, n_omega_out))
    !
    !  if(gamma_inc .gt. 1d-5) then
    !    sigma_tmp = lambda_lp
    !    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
    !    sigma_tmp = lambda_ln
    !    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
    !    sigma_tmp = lambda_cp
    !    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
    !  end if
    !  
    !end if
    
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

  subroutine calculate_SCKH_res_PES_old(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    !use m_SCKH_utils, only: verlet_trajectory
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_SCKH_utils, only: compute_F_if_omp_many_n
    use m_SCKH_utils, only: compute_F_if_omp_sum_n
    use m_SCKH_utils, only: compute_F_if_omp_one_n
    !use m_SCKH_utils, only: compute_F_if_om_many_n
    use m_SCKH_utils, only: compute_F_if_om_omp
    use m_SCKH_utils, only: compute_F_if_om_omp_one_f
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
    use m_KH_utils, only: convolute_instrumental
    use m_spectrum_utils, only: convolution_lorentzian_grid_fft_many_freq
    
    type(sckh_params_t), intent(inout):: p 

    real(wp), allocatable::  time(:), E_dyn_inp(:),  E_dyn2_inp(:)  
    real(wp), allocatable:: E_n0(:), E_fn_corr(:,:)
    real(wp), allocatable::  E_lp_corr(:), shift(:)
    real(wp), allocatable::  E_i_inp(:), E_i1(:), E_i2(:)
    real(wp), allocatable::  E_n_inp(:,:), E_n1(:,:), E_n2(:,:)
    real(wp), allocatable::  E_f_inp(:,:)
    real(wp), allocatable::  E_f1(:,:), E_f2(:,:)
    real(wp), allocatable::  E_fc1(:,:,:), E_fc2(:,:,:)

    real(wp), allocatable:: D_fn_inp(:,:,:,:), D_fn1(:,:,:,:), D_fn2(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:), D_ni1(:,:,:), D_ni2(:,:,:) 
    real(wp):: E_nf_mean, E_fi_mean, E_ni_mean
    real(wp), allocatable:: E_i_mean(:)
    
    character(80)::  string, file
    real(wp), allocatable:: x_sampl(:), mom_sampl(:), x_mom_sampl(:,:)
    real(wp), allocatable:: x_new(:), v_new(:), a_new(:)
    real(wp), allocatable:: x_new2(:), v_new2(:), a_new2(:)

    integer:: npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl
    integer:: ninter, nfinal, ntsteps  
    integer:: npoints_in
    integer:: nfinal_tot

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
    integer::i,j,m,m1,m2, traj, traj2, n_e, f_e, fn_e

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
         E_n0(npoints_in),&
         E_lp_corr(npoints_in),&
         D_ni_inp(ninter, npoints_in,3), &
         D_fn_inp(nfinal, ninter, npoints_in,3), &
         time(ntsteps),&
         E_i1(ntsteps), &
         E_i2(ntsteps), &
         E_n1(ninter, ntsteps), &
         E_n2(ninter, ntsteps), &
         !E_f1(nfinal,ntsteps), &
         !E_f2(nfinal,ntsteps), &
         D_ni1(ninter, ntsteps,3), &
         D_ni2(ninter, ntsteps,3), &
         D_fn2(nfinal, ninter, ntsteps,3), &
         D_fn1(nfinal, ninter, ntsteps,3), &
         E_i_mean(ntsteps), &
         shift(npoints_in), &
         c_i(npoints_in,npoints_in), &
         eig_i(npoints_in), &
         x_new(ntsteps),&
         v_new(ntsteps),&
         a_new(ntsteps),&
         x_new2(ntsteps),&
         v_new2(ntsteps),&
         a_new2(ntsteps), &
         X_r(npoints_in))
    
    ! in ORB mode there are nfinal * ninter final states
    write(6,*) "p % KH_states_mode = ", p % KH_states_mode
    if (upper(p % KH_states_mode) .eq. "STATES") then
      nfinal_tot = nfinal

      allocate(E_fc1(nfinal,1, ntsteps), &
           E_fc2(nfinal,1,ntsteps))
      
    else if (upper(p % KH_states_mode) .eq. "ORBS") then
      nfinal_tot = nfinal * ninter

      allocate(E_fc1(nfinal,ninter, ntsteps), &
           E_fc2(nfinal,ninter, ntsteps))
    else
      write(6,*) "p % KH_states_mode should be either 'STATES' or 'ORBS'"
    end if

    
    ! set up grid points
    do i = 1, npoints_in
      X_r(i) = (i-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_r, E_i_inp)
    !call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n_inp)

    ! intermediate state reference energy (lowest state) 
    if(p % use_n0_state) then 
      call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n0)
    else
      E_n0 =0.0d0
    end if

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

    ! read list of corrections to the final state files coming from the excited electron
    if (upper(p % KH_states_mode) .eq. "ORBS") then    
      allocate( E_fn_corr(p % npesfile_n, p % npoints_in))

      call read_file_list(p % pes_file_list_fn_corr, p % npesfile_n, p % pes_files_fn_corr)
      
      do j=1,p % npesfile_n
        call read_PES_file(p % pes_files_fn_corr(j), p % npoints_in, &
             p % npoints_in, X_r, E_fn_corr(j,:))
      end do
    else
      
    end if


    ! XXX dangerous hack!
    !E_dyn_inp = (E_i_inp + E_f_inp(1,:)) /2.0_wp
    !E_dyn_inp = (E_f_inp(1,:) + E_n_inp(1,:) + E_i_inp) /3.0_wp
    !E_dyn_inp = (E_n_inp(1,:) + E_f_inp(1,:)) /2.0_wp
    !E_dyn_inp = (E_i_inp + E_f_inp(1,:)) /2.0_wp
    !E_dyn2_inp = (E_n_inp(1,:) + E_f_inp(1,:)) /2.0_wp
    
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
    E_n0 = E_n0  / const % eV
    E_dyn_inp = E_dyn_inp  / const % eV
    E_dyn2_inp = E_dyn2_inp  / const % eV
    E_f_inp = E_f_inp / const % eV
    if(allocated( E_fn_corr))  E_fn_corr =  E_fn_corr / const % eV
    
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

    allocate(F_if_t_omp(nfinal_tot, n_omega_in, n_omega_out,3,3),&
         F_if_om_omp(nfinal_tot, n_omega_in, n_omega_out,3,3),&
         lambda_F(nfinal_tot, n_omega_in, n_omega_out),&
         lambda_G(nfinal_tot, n_omega_in, n_omega_out),&
         lambda_H(nfinal_tot, n_omega_in, n_omega_out),&
         lambda_lp(n_omega_in, n_omega_out),&
         lambda_ln(n_omega_in, n_omega_out),&
         lambda_cp(n_omega_in, n_omega_out),&
         sigma_tmp(n_omega_in, n_omega_out))

!    allocate(F_if_t_omp_tmp(n_omega_in, n_omega_out,3,3),&
!         F_if_t_omp(n_omega_in, n_omega_out,3,3),&
!         F_if_om_omp(n_omega_in, n_omega_out,3,3),&
!         lambda_F(n_omega_in, n_omega_out),&
!         lambda_G(n_omega_in, n_omega_out),&
!         lambda_H(n_omega_in, n_omega_out),&
!         lambda_lp(n_omega_in, n_omega_out),&
!         lambda_ln(n_omega_in, n_omega_out),&
!         lambda_cp(n_omega_in, n_omega_out),&
!         sigma_tmp(n_omega_in, n_omega_out))
    
    !
    ! Loop over trajectories
    !

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do

    ! introduce factor for possible time development backwards in time
    !write(6,*) "upper(p % KH_amplitude_mode)", upper(p % KH_amplitude_mode)
    !if(upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
    !  fac_t = 1.0_wp
    !else if(upper(p % KH_amplitude_mode) .eq. "INGOING") then
    !  fac_t = 1.0_wp
    !else
    !  write(6,*) "p % KH_amplitude_mode must be either OUTGOING or INGOING"
    !end if
      
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
           E_dyn_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )
      
      call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
      do i=1,ninter
        call spline_easy(X_r, E_n_inp(i,:) + E_n0, npoints_in, x_new, E_n1(i,:), ntsteps)  
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
      
      !do i=1,nfinal
      !  call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_f1(i,:), ntsteps)  
      !end do

      if (upper(p % KH_states_mode) .eq. "STATES") then
        do i=1,nfinal
          call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_fc1(i,1,:), ntsteps)  
        end do
            
      else if (upper(p % KH_states_mode) .eq. "ORBS") then
        do i=1,nfinal
          do j=1,ninter
            call spline_easy(X_r, E_f_inp(i,:) + E_fn_corr(j,:), npoints_in, x_new, E_fc1(i,j,:), ntsteps)  
          end do
        end do
        
      end if
      
      ! first time, compute the mean transition energy, and frequency
      if (traj .eq. 1) then
        ind = minloc(E_i_inp)

        !E_ni_mean =  E_n1(nfinal, ind(1)) - E_i1(ind(1))
        E_ni_mean =  E_n1(1, ind(1)) - E_i1(ind(1))
        write(6,*) "E_ni_mean", E_ni_mean
        
        !E_fi_mean =  E_f1(nfinal, ind(1)) - E_i1(ind(1))
        E_fi_mean =  E_fc1(1, 1, ind(1)) - E_i1(ind(1))
        write(6,*) "E_fi_mean", E_fi_mean

        !E_nf_mean =  E_n1(ninter, ind(1)) - E_f1(nfinal,ind(1))
        E_nf_mean =  E_n1(1, ind(1)) - E_fc1(1,1,ind(1))
        write(6,*) "E_nf_mean", E_nf_mean

        ! this is to take away the initial state vib effects
        E_i_mean =  E_i1(ind(1))
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
        omega_in = omega_in + E_ni_mean !E_nf_mean + E_fi_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
        omega_out = omega_out + E_nf_mean
        
      end if

      ! loop over each starting point of first trajectory
      do traj2 =1, ntsteps ! possibly with a stride
        
        !if(upper(p % runmode_sckh_res) .eq. "FULL") then

        call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
             E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )

        !else if(upper(p % runmode_sckh_res) .eq. "APPROX") then
        !  ! run from same starting point, then add the energies of the runs on intermediate and initial states
        !  call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
        !       E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
        !else
        !  write(6,*) "upper(p % sckh_res_runmode) must be either 'FULL' or 'APPROX' "
        !end if

        
        !call verlet_trajectory_xva(x_new(1), v_new(1), X_r, &
        !     E_dyn2_inp * const % eV, fac_t * delta_t, mu_SI, x_new2, v_new2, a_new2 )
        
        call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)

        do i=1,ninter
          call spline_easy(X_r, E_n_inp(i,:) + E_n0, npoints_in, x_new2, E_n2(i,:), ntsteps)
        end do
        
        !do i=1,nfinal
        !  call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new2, E_f2(i,:), ntsteps)  
        !end do
        
        if (upper(p % KH_states_mode) .eq. "STATES") then
          do i=1,nfinal
            call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new2, E_fc2(i,1,:), ntsteps)  
          end do
          
        else if (upper(p % KH_states_mode) .eq. "ORBS") then
          do i=1,nfinal
            do j=1,ninter
              call spline_easy(X_r, E_f_inp(i,:) + E_fn_corr(j,:), npoints_in, x_new2, E_fc2(i,j,:), ntsteps)  
            end do
          end do

        end if

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

        !if(upper(p % runmode_sckh_res) .eq. "APPROX") then
        !  E_i2 =  E_i2 +  (E_i1(traj2) - E_i1(1))
        !
        !  do i=1,ninter
        !    E_n2(i,:) =  E_n2(i,:) +  (E_n1(i,traj2) - E_n1(i,1))
        !  end do
        !
        !  do i=1,nfinal
        !    E_f2(i,:) =  E_f2(i,:) +  (E_f1(i,traj2) - E_f1(i,1))
        !  end do
        !  
        !end if
        
        
        !if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then
          
        !call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,:,:,:), &
        !     D_ni2(:, 1, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)

        if (upper(p % KH_states_mode) .eq. "STATES") then
          do f_e = 1, nfinal 
            call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
                 D_ni2(:,1,:), time, F_if_t_omp(f_e,traj2,:,:,:), gamma)
            !call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
            !     D_ni2(:,1,:), time, F_if_t_omp_tmp(:,:,:), gamma)
            !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
          end do

          !call compute_F_if_omp_many_n(E_n2(:,:), E_fc2(:,1,:), E_nf_mean, D_fn2(:,1,:,:), &
          !     D_ni2(:, 1, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)

        else if (upper(p % KH_states_mode) .eq. "ORBS") then
          do f_e = 1, nfinal
            do n_e = 1, ninter
              fn_e = ninter * (f_e -1) + n_e  
              call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
                   D_ni2(n_e,1,:), time, F_if_t_omp(fn_e,traj2,:,:,:), gamma)
              !call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
              !     D_ni2(n_e,1,:), time, F_if_t_omp_tmp(:,:,:), gamma)
              !
              !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
            end do
          end do
          
        end if
        
        !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
        !  ! F_if_t_omp will now store F_if_t_om temporarily
        !  call compute_F_if_om_many_n(E_n2(:,:), E_i2(:), E_ni_mean, D_fn2(:,1,:,:), &
        !       D_ni2(:, :, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)
        !  
        !end if
        
      end do ! do traj2 =1, ntsteps2


      !if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

      ! for test purposes, write F_if_t_omp to file
      !call write_F_if_t_omp(p, F_if_t_omp, omega_out, time)
     ! F_if_om_omp = 0.0d0
      if (upper(p % KH_states_mode) .eq. "STATES") then

        do f_e = 1, nfinal         
          call compute_F_if_om_omp_one_f(F_if_t_omp(f_e, :,:,:,:), E_fc1(f_e,1,:), E_fi_mean, time, &
               E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(f_e,:,:,:,:)) !F_if_om_omp)
        end do
        
        !call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fi_mean, time, &
        !     E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)
   
        !call compute_F_if_om_omp(F_if_t_omp(:, :,:,:,:), E_fc1(:,1,:), E_fi_mean, time, &
        !     E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:,:)) !F_if_om_omp)
        
      else if (upper(p % KH_states_mode) .eq. "ORBS") then

        do f_e = 1, nfinal
          do n_e = 1, ninter
            fn_e = ninter * (f_e -1) + n_e  
            call compute_F_if_om_omp_one_f(F_if_t_omp(fn_e, :,:,:,:), E_fc1(f_e,n_e,:), E_fi_mean, time, &
                 E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(fn_e,:,:,:,:)) !F_if_om_omp)
          end do
        end do
        
      end if

      !! below to remove dependence on E_i1, so to quantize that energy
        !call compute_F_if_om_omp(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
        !     E_i_mean, gamma_inc, omega_out, E_nf_mean, F_if_om_omp)
        
      !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then
      !  call compute_F_if_om_omp_ingoing(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
      !       E_i1, gamma_instr, omega_in, E_ni_mean, F_if_om_omp)
      !  
      !end if
        
        
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

!if (.false.) then
    !write(6,*) "Entering convolute_incoming"
    write(6,*) "Entering convolute_incoming, broadening ",  upper(p % broadening_func_inc)
    
    if(gamma_inc .gt. 1d-5) then
      sigma_tmp = lambda_lp
      !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_lp, 1,  upper(p % broadening_func_inc))
      sigma_tmp = lambda_ln
      !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(p % broadening_func_inc))
      sigma_tmp = lambda_cp
      !call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_cp, 1, upper(p % broadening_func_inc))
    end if
    
    write(6,*) "Done"
    !write(6,*) "Entering convolute_instrumental"
    write(6,*) "Entering convolute_instrumental, broadening  ",  upper(p % broadening_func_instr) 
    
    if(gamma_instr .gt. 1d-5) then
      sigma_tmp = lambda_lp
      !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_lp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_ln
      !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_ln)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_cp
      !call convolute_instrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_cp)
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(p % broadening_func_instr))
    end if
    
    write(6,*) "Done"
!  end if
    
    
    !if (upper(p % KH_amplitude_mode) .eq. "OUTGOING") then

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
      
    !else if (upper(p % KH_amplitude_mode) .eq. "INGOING") then      
    !  allocate(sigma_tmp(n_omega_in, n_omega_out))
!
!      if(gamma_inc .gt. 1d-5) then
!        sigma_tmp = lambda_lp
!        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_lp)
!        sigma_tmp = lambda_ln
!        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_ln)
!        sigma_tmp = lambda_cp
!        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_cp)
!      end if
      
    !  end if
    
    ! write spectra to individual files
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

    ! write spectra to file
   file="_sigma_all"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in
     do i=1, n_omega_out
       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
   end do
   
   close(ifile) 


   ! write spectra to file
   file="_sigma_all_nogrid"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in
     do i=1, n_omega_out
       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
     write(ifile, *) 
   end do
   
   close(ifile) 
   
 end subroutine calculate_SCKH_res_PES_old



 ! lifted out routine to compute full resonant sckh pes
 subroutine compute_F_SCKH_res_PES_full(p, x_new, v_new, X_r, E_dyn2_inp, delta_t, mu_SI, &
      E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
      E_i1, E_fc1, gamma_inc, omega_out, F2_tensor) 

   use m_precision, only: wp
   use m_constants, only: const
   use m_sckh_params_t, only: sckh_params_t 
   use m_upper, only : upper
   use m_SCKH_utils, only: compute_F_if_omp_sum_n
   use m_SCKH_utils, only: compute_F_if_omp_one_n
   use m_SCKH_utils, only: compute_F_if_om_omp
   use m_SCKH_utils, only: compute_F_if_om_omp_one_f
   use m_splines, only: spline_easy
   use m_SCKH_utils, only: verlet_trajectory_xva
   
   type(sckh_params_t), intent(in):: p
   real(wp), intent(in):: x_new(:)
   real(wp), intent(in):: v_new(:)
   real(wp), intent(in)::  X_r(:)
   real(wp), intent(in):: E_dyn2_inp(:)
   real(wp), intent(in):: delta_t
   real(wp), intent(in):: mu_SI
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_n_inp(:,:)
   real(wp), intent(in):: E_n0(:)
   real(wp), intent(in):: E_f_inp(:,:)
   real(wp), intent(in):: E_fn_corr(:,:)
   real(wp), intent(in):: D_fn_inp(:,:,:,:)
   real(wp), intent(in):: D_ni_inp(:,:,:)
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in)::  time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: E_i1(:)
   real(wp), intent(in):: E_fc1(:,:,:)   
   real(wp), intent(in):: gamma_inc
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(out)::  F2_tensor(:,:,:,:,:,:)

   real(wp), allocatable:: x_new2(:)
   real(wp), allocatable:: v_new2(:)
   real(wp), allocatable:: a_new2(:)
   real(wp), allocatable:: E_i2(:)
   real(wp), allocatable:: E_n2(:,:)
   real(wp), allocatable:: E_fc2(:,:,:)
   real(wp), allocatable:: D_fn2(:,:,:,:)
   real(wp), allocatable:: D_ni2(:,:,:)

   complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
   complex(wp), allocatable::  F_if_om_omp(:,:,:,:)    

   integer:: f_e
   integer:: n_e
   integer:: nfinal
   integer:: ninter
   integer:: traj2
   integer:: ntsteps
   integer:: npoints_in
   integer:: m
   integer:: m1
   integer:: m2
   integer:: m3
   integer:: m4

   ntsteps = size(time)
   ninter = size(E_n_inp,1)
   nfinal = size(E_f_inp,1)
   npoints_in = size(E_n_inp, 2)
   
   allocate(x_new2(ntsteps))
   allocate(v_new2(ntsteps))
   allocate(a_new2(ntsteps))

   allocate(E_i2(ntsteps))
   allocate(E_n2(ninter, ntsteps))
   allocate(D_ni2(ninter, ntsteps,3))
   allocate(D_fn2(nfinal, ninter, ntsteps,3))

   allocate(F_if_t_omp(ntsteps, ntsteps,3,3))
   allocate(F_if_om_omp(ntsteps, ntsteps,3,3))

   F2_tensor = 0.0_wp

   ! here two options depending on p % KH_states_mode
   if (upper(p % KH_states_mode) .eq. "STATES") then

     allocate(E_fc2(nfinal,1, ntsteps))
     
     ! now the f_e loop is outside of traj2 to save memory
     do f_e = 1,nfinal
       
       do traj2 =1, ntsteps ! possibly with a stride
         
         ! stupid to recalcualte trajectory all the time... could be lifted out of f_e loop by saving x_new2 for each traj2
         call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
              E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )


         ! also E_i2 and E_n2 could be lifted out if they take a lot of time
         call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)

         do n_e=1,ninter
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
         end do
         
         call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new2, E_fc2(f_e,1,:), ntsteps)  
         
         do m=1,3
           call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
         end do

         do n_e=1,ninter
           do m=1,3
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do
         end do

         call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
              D_ni2(:,1,:), time, F_if_t_omp(traj2,:,:,:), gamma)

         
         !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
         
       end do !do traj2 =1, ntsteps

       call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fi_mean, time, &
            E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)



       
       !! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       !do m1=1,3
       !  do m2=1,3
       !    lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
       !    lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
       !    lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
       !  end do
       !end do

       ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       do m1=1,3
         do m2=1,3
           do m3=1,3
             do m4=1,3
               !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
               !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
               !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
               F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
                    real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
             end do
           end do
         end do
       end do
       write(6,*) "here 3"
       
     end do !do f_e = 1,nfinal
     
   else if (upper(p % KH_states_mode) .eq. "ORBS") then

     allocate(E_fc2(nfinal,ninter, ntsteps))

     do f_e = 1,nfinal
       do n_e = 1,ninter
         
         !fn_e = ninter * (f_e -1) + n_e  
         
         do traj2 =1, ntsteps ! possibly with a stride              
           
           call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
                E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
           
           call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
           
           ! this one really belongs here now
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
           
           call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new2, E_fc2(f_e,n_e,:), ntsteps)  
           
           do m=1,3
             call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do
           
           call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
                D_ni2(n_e,1,:), time, F_if_t_omp(traj2,:,:,:), gamma)
           
         end do ! do traj2 =1, ntsteps
         
         call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,n_e,:), E_fi_mean, time, &
              E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) 
         
         ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
         do m1=1,3
           do m2=1,3
             do m3=1,3
               do m4=1,3
                 !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
                 !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
                 !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
                 F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
                      real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
               end do
             end do
           end do
         end do

         
       end do! do n_e = 1,ninter
     end do! do f_e = 1,nfinal
     
   end if !  if (upper(p % KH_states_mode) .eq. "STATES") then
   
 end subroutine compute_F_SCKH_res_PES_full

 ! lifted out routine to compute full resonant sckh pes
 subroutine compute_F_SCKH_res_PES_full_separate(p, x_new, v_new, X_r, E_dyn2_inp, delta_t, mu_SI, &
      E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
      E_i1, E_n1, E_fc1, D_ni1, gamma_inc, omega_out, omega_in, F2_tensor) 

   use m_precision, only: wp
   use m_constants, only: const
   use m_sckh_params_t, only: sckh_params_t 
   use m_upper, only : upper
   use m_SCKH_utils, only: compute_F_if_omp_sum_n
   use m_SCKH_utils, only: compute_F_if_omp_one_n
   use m_SCKH_utils, only: compute_F_if_om_omp
   use m_SCKH_utils, only: compute_F_if_om_omp_one_f
   use m_SCKH_utils, only: compute_F_if_om_omp_one_f_separate
   use m_SCKH_utils, only: compute_F_if_om_omp_one_f_factor_XAS
   use m_splines, only: spline_easy
   use m_SCKH_utils, only: verlet_trajectory_xva
   
   type(sckh_params_t), intent(in):: p
   real(wp), intent(in):: x_new(:)
   real(wp), intent(in):: v_new(:)
   real(wp), intent(in)::  X_r(:)
   real(wp), intent(in):: E_dyn2_inp(:)
   real(wp), intent(in):: delta_t
   real(wp), intent(in):: mu_SI
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_n_inp(:,:)
   real(wp), intent(in):: E_n0(:)
   real(wp), intent(in):: E_f_inp(:,:)
   real(wp), intent(in):: E_fn_corr(:,:)
   real(wp), intent(in):: D_fn_inp(:,:,:,:)
   real(wp), intent(in):: D_ni_inp(:,:,:)
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in)::  time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: E_i1(:)
   real(wp), intent(in):: E_n1(:,:)   
   real(wp), intent(in):: E_fc1(:,:,:)   
   real(wp), intent(in):: D_ni1(:,:,:)   
   real(wp), intent(in):: gamma_inc
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(out)::  F2_tensor(:,:,:,:,:,:)

   real(wp), allocatable:: x_new2(:)
   real(wp), allocatable:: v_new2(:)
   real(wp), allocatable:: a_new2(:)
   real(wp), allocatable:: E_i2(:)
   real(wp), allocatable:: E_n2(:,:)
   real(wp), allocatable:: E_fc2(:,:,:)
   real(wp), allocatable:: D_fn2(:,:,:,:)
   real(wp), allocatable:: D_ni2(:,:,:)

   complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
   complex(wp), allocatable::  F_if_om_omp(:,:,:,:)    

   integer:: f_e
   integer:: n_e
   integer:: nfinal
   integer:: ninter
   integer:: traj2
   integer:: ntsteps
   integer:: npoints_in
   integer:: m
   integer:: m1
   integer:: m2
   integer:: m3
   integer:: m4

   ntsteps = size(time)
   ninter = size(E_n_inp,1)
   nfinal = size(E_f_inp,1)
   npoints_in = size(E_n_inp, 2)
   
   allocate(x_new2(ntsteps))
   allocate(v_new2(ntsteps))
   allocate(a_new2(ntsteps))

   allocate(E_i2(ntsteps))
   allocate(E_n2(ninter, ntsteps))
   allocate(D_ni2(ninter, ntsteps,3))
   allocate(D_fn2(nfinal, ninter, ntsteps,3))

   allocate(F_if_t_omp(ntsteps, ntsteps,3,3))
   allocate(F_if_om_omp(ntsteps, ntsteps,3,3))

   F2_tensor = 0.0_wp

   ! here two options depending on p % KH_states_mode
   if (upper(p % KH_states_mode) .eq. "STATES") then

     allocate(E_fc2(nfinal,1, ntsteps))
     
     ! now the f_e loop is outside of traj2 to save memory
     do f_e = 1,nfinal
       
       do traj2 =1, ntsteps ! possibly with a stride
         
         ! stupid to recalcualte trajectory all the time... could be lifted out of f_e loop by saving x_new2 for each traj2
         call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
              E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )


         ! also E_i2 and E_n2 could be lifted out if they take a lot of time
         call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)

         do n_e=1,ninter
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
         end do
         
         call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new2, E_fc2(f_e,1,:), ntsteps)  
         
         do m=1,3
           call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
         end do

         do n_e=1,ninter
           do m=1,3
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do
         end do

         ! here changed D_ni2 to D_ni1
         call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
              D_ni1(:,traj2,:), time, F_if_t_omp(traj2,:,:,:), gamma)

         !call compute_F_if_omp_sum_n(E_n1(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
         !     D_ni1(:,traj2,:), time, F_if_t_omp(traj2,:,:,:), gamma)

         
         !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
         
       end do !do traj2 =1, ntsteps

       ! here we should have exp(-\int E_in(t) -E_nf(t) - E_fi_mean) where E_in comes from grouind state dynamcis
       ! and E_nf comes from excited state dynamics
       ! since n is present we need to also include this before the sum over n, but for one n this should be ok

       call compute_F_if_om_omp_one_f_separate(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fc2(f_e,1,:), E_fi_mean, time, &
            E_i1, E_n1(1,:), E_n2(1,:), gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)

       !call compute_F_if_om_omp_one_f_factor_XAS(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fc2(f_e,1,:), E_fi_mean, time, &
       !     E_i1, E_n1(1,:), E_n2(1,:), gamma_inc, omega_out, omega_in, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)



       
       !! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       !do m1=1,3
       !  do m2=1,3
       !    lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
       !    lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
       !    lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
       !  end do
       !end do

       ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       do m1=1,3
         do m2=1,3
           do m3=1,3
             do m4=1,3
               !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
               !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
               !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
               F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
                    real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
             end do
           end do
         end do
       end do
       write(6,*) "here 3"
       
     end do !do f_e = 1,nfinal
     
   else if (upper(p % KH_states_mode) .eq. "ORBS") then

     allocate(E_fc2(nfinal,ninter, ntsteps))

     do f_e = 1,nfinal
       do n_e = 1,ninter
         
         !fn_e = ninter * (f_e -1) + n_e  
         
         do traj2 =1, ntsteps ! possibly with a stride              
           
           call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
                E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
           
           call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
           
           ! this one really belongs here now
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
           
           call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new2, E_fc2(f_e,n_e,:), ntsteps)  
           
           do m=1,3
             call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do
           
           
           call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
                D_ni2(n_e,1,:), time, F_if_t_omp(traj2,:,:,:), gamma)
           
         end do ! do traj2 =1, ntsteps
         
         call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,n_e,:), E_fi_mean, time, &
              E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) 
         
         ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
         do m1=1,3
           do m2=1,3
             do m3=1,3
               do m4=1,3
                 !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
                 !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
                 !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
                 F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
                      real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
               end do
             end do
           end do
         end do

         
       end do! do n_e = 1,ninter
     end do! do f_e = 1,nfinal
     
   end if !  if (upper(p % KH_states_mode) .eq. "STATES") then
   
 end subroutine compute_F_SCKH_res_PES_full_separate


 ! lifted out routine to compute full resonant sckh pes
 subroutine compute_F_SCKH_res_PES_full_separate_ingoing(p, x_new, v_new, X_r, E_dyn2_inp, delta_t, mu_SI, &
      E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean, E_ni_mean, time, gamma, &
      E_i1, E_n1, E_fc1, D_ni1, D_fn1, gamma_inc, omega_out, omega_in, F2_tensor) 

   use m_precision, only: wp
   use m_constants, only: const
   use m_sckh_params_t, only: sckh_params_t 
   use m_upper, only : upper
   use m_SCKH_utils, only: compute_F_if_omp_sum_n
   use m_SCKH_utils, only: compute_F_if_om_sum_n
   use m_SCKH_utils, only: compute_F_if_omp_one_n
   use m_SCKH_utils, only: compute_F_if_om_omp
   use m_SCKH_utils, only: compute_F_if_om_omp_one_f
   use m_SCKH_utils, only: compute_F_if_om_omp_one_f_separate
      use m_SCKH_utils, only: compute_F_if_om_omp_one_f_separate_ingoing
   use m_splines, only: spline_easy
   use m_SCKH_utils, only: verlet_trajectory_xva
   
   type(sckh_params_t), intent(in):: p
   real(wp), intent(in):: x_new(:)
   real(wp), intent(in):: v_new(:)
   real(wp), intent(in)::  X_r(:)
   real(wp), intent(in):: E_dyn2_inp(:)
   real(wp), intent(in):: delta_t
   real(wp), intent(in):: mu_SI
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_n_inp(:,:)
   real(wp), intent(in):: E_n0(:)
   real(wp), intent(in):: E_f_inp(:,:)
   real(wp), intent(in):: E_fn_corr(:,:)
   real(wp), intent(in):: D_fn_inp(:,:,:,:)
   real(wp), intent(in):: D_ni_inp(:,:,:)
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in):: E_ni_mean
   real(wp), intent(in)::  time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: E_i1(:)
   real(wp), intent(in):: E_n1(:,:)   
   real(wp), intent(in):: E_fc1(:,:,:)   
   real(wp), intent(in):: D_ni1(:,:,:)
   real(wp), intent(in):: D_fn1(:,:,:,:)   
   real(wp), intent(in):: gamma_inc
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(out)::  F2_tensor(:,:,:,:,:,:)

   real(wp), allocatable:: x_new2(:)
   real(wp), allocatable:: v_new2(:)
   real(wp), allocatable:: a_new2(:)
   real(wp), allocatable:: E_i2(:)
   real(wp), allocatable:: E_n2(:,:)
   real(wp), allocatable:: E_fc2(:,:,:)
   real(wp), allocatable:: D_fn2(:,:,:,:)
   real(wp), allocatable:: D_ni2(:,:,:)

   complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
   complex(wp), allocatable::  F_if_om_omp(:,:,:,:)    

   integer:: f_e
   integer:: n_e
   integer:: nfinal
   integer:: ninter
   integer:: traj2
   integer:: ntsteps
   integer:: npoints_in
   integer:: m
   integer:: m1
   integer:: m2
   integer:: m3
   integer:: m4

   ntsteps = size(time)
   ninter = size(E_n_inp,1)
   nfinal = size(E_f_inp,1)
   npoints_in = size(E_n_inp, 2)
   
   allocate(x_new2(ntsteps))
   allocate(v_new2(ntsteps))
   allocate(a_new2(ntsteps))

   allocate(E_i2(ntsteps))
   allocate(E_n2(ninter, ntsteps))
   allocate(D_ni2(ninter, ntsteps,3))
   allocate(D_fn2(nfinal, ninter, ntsteps,3))

   allocate(F_if_t_omp(ntsteps, ntsteps,3,3))
   allocate(F_if_om_omp(ntsteps, ntsteps,3,3))

   F2_tensor = 0.0_wp

   ! here two options depending on p % KH_states_mode
   if (upper(p % KH_states_mode) .eq. "STATES") then

     allocate(E_fc2(nfinal,1, ntsteps))
     
     ! now the f_e loop is outside of traj2 to save memory
     do f_e = 1,nfinal
       
       do traj2 =1, ntsteps ! possibly with a stride
         
         ! stupid to recalcualte trajectory all the time... could be lifted out of f_e loop by saving x_new2 for each traj2
         call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
              E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )


         ! also E_i2 and E_n2 could be lifted out if they take a lot of time
         call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)

         do n_e=1,ninter
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
         end do
         
         call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new2, E_fc2(f_e,1,:), ntsteps)  
         
         do m=1,3
           call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
         end do

         do n_e=1,ninter
           do m=1,3
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do
         end do

!         ! here changed D_ni2 to D_ni1
!         call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
!              D_ni1(:,traj2,:), time, F_if_t_omp(traj2,:,:,:), gamma)

         ! ingoing now
         call compute_F_if_om_sum_n(E_n2(:,:), E_i2(:), E_ni_mean, D_fn1(f_e,:,traj2,:),  D_ni2(:,:,:),&
              time, F_if_t_omp(traj2,:,:,:), gamma)

         
         !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
         
       end do !do traj2 =1, ntsteps

       ! here we should have exp(-\int E_in(t) -E_nf(t) - E_fi_mean) where E_in comes from grouind state dynamcis
       ! and E_nf comes from excited state dynamics
       ! since n is present we need to also include this before the sum over n, but for one n this should be ok


       !call compute_F_if_om_omp_one_f_separate(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fc2(f_e,1,:), E_fi_mean, time, &
       !     E_i1, E_n1(1,:), E_n2(1,:), gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)

       call compute_F_if_om_omp_one_f_separate_ingoing(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fc2(f_e,1,:), E_fi_mean, time, &
            E_i1, E_i2, E_n1(1,:), E_n2(1,:), gamma_inc, omega_in, E_ni_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)



       
       !! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       !do m1=1,3
       !  do m2=1,3
       !    lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
       !    lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
       !    lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
       !  end do
       !end do

       ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       do m1=1,3
         do m2=1,3
           do m3=1,3
             do m4=1,3
               !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
               !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
               !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
               F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
                    real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
             end do
           end do
         end do
       end do
       write(6,*) "here 3"
       
     end do !do f_e = 1,nfinal
     
   else if (upper(p % KH_states_mode) .eq. "ORBS") then

     allocate(E_fc2(nfinal,ninter, ntsteps))

     do f_e = 1,nfinal
       do n_e = 1,ninter
         
         !fn_e = ninter * (f_e -1) + n_e  
         
         do traj2 =1, ntsteps ! possibly with a stride              
           
           call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
                E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
           
           call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
           
           ! this one really belongs here now
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
           
           call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new2, E_fc2(f_e,n_e,:), ntsteps)  
           
           do m=1,3
             call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do
           
           
           call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
                D_ni2(n_e,1,:), time, F_if_t_omp(traj2,:,:,:), gamma)
           
         end do ! do traj2 =1, ntsteps
         
         call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,n_e,:), E_fi_mean, time, &
              E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) 
         
         ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
         do m1=1,3
           do m2=1,3
             do m3=1,3
               do m4=1,3
                 !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
                 !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
                 !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
                 F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
                      real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
               end do
             end do
           end do
         end do

         
       end do! do n_e = 1,ninter
     end do! do f_e = 1,nfinal
     
   end if !  if (upper(p % KH_states_mode) .eq. "STATES") then
   
 end subroutine compute_F_SCKH_res_PES_full_separate_ingoing

 

 
 ! lifted out routine to compute sckh where each traj is factorized 
 subroutine compute_F_SCKH_res_PES_factor_each_traj(p, x_mom_sampl, traj, X_r, E_dyn2_inp, delta_t, mu_SI, &
      E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
      E_i1, E_fc1, gamma_R, omega_out, F2_tensor) 

   use m_precision, only: wp
   use m_constants, only: const
   use m_sckh_params_t, only: sckh_params_t 
   use m_upper, only : upper
   use m_SCKH_utils, only: compute_F_if_omp_sum_n
   use m_SCKH_utils, only: compute_F_if_omp_one_n
   use m_SCKH_utils, only:compute_F_if_om_omp_no_F_one

    !use m_SCKH_utils, only: compute_F_if_omp_one_n
   !use m_SCKH_utils, only: compute_F_if_om_omp
   !use m_SCKH_utils, only: compute_F_if_om_omp_one_f
   use m_splines, only: spline_easy
   use m_SCKH_utils, only: verlet_trajectory_xva
   use m_SCKH_utils, only: compute_F_if_om_omp_no_F
   
   type(sckh_params_t), intent(in):: p
   real(wp), intent(in):: x_mom_sampl(:,:)
   integer, intent(in):: traj
   !real(wp), intent(in):: x_new(:)
   !real(wp), intent(in):: v_new(:)
   real(wp), intent(in)::  X_r(:)
   real(wp), intent(in):: E_dyn2_inp(:)
   real(wp), intent(in):: delta_t
   real(wp), intent(in):: mu_SI
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_n_inp(:,:)
   real(wp), intent(in):: E_n0(:)
   real(wp), intent(in):: E_f_inp(:,:)
   real(wp), intent(in):: E_fn_corr(:,:)
   real(wp), intent(in):: D_fn_inp(:,:,:,:)
   real(wp), intent(in):: D_ni_inp(:,:,:)
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in)::  time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: E_i1(:)
   real(wp), intent(in):: E_fc1(:,:,:)   
   real(wp), intent(in):: gamma_R
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(out)::  F2_tensor(:,:,:,:,:,:)

   real(wp), allocatable:: x_new2(:)
   real(wp), allocatable:: v_new2(:)
   real(wp), allocatable:: a_new2(:)
   real(wp), allocatable:: E_i2(:)
   real(wp), allocatable:: E_n2(:,:)
   real(wp), allocatable:: E_fc2(:,:,:)
   real(wp), allocatable:: D_fn2(:,:,:,:)
   real(wp), allocatable:: D_ni2(:,:,:)

   complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
   complex(wp), allocatable::  F_if_om_omp(:,:,:,:)    

   complex(wp), allocatable::  R_if_om_omp(:,:,:)    
   real(wp), allocatable:: R_factor(:,:,:)
   
   integer:: f_e
   integer:: n_e
   integer:: nfinal
   integer:: ninter
   integer:: traj2
   integer:: ntsteps
   integer:: npoints_in
   integer:: m
   integer:: m1
   integer:: m2
   integer:: m3
   integer:: m4
   integer:: om_in

   ntsteps = size(time)
   ninter = size(E_n_inp,1)
   nfinal = size(E_f_inp,1)
   npoints_in = size(E_n_inp, 2)
   
   allocate(x_new2(ntsteps))
   allocate(v_new2(ntsteps))
   allocate(a_new2(ntsteps))

   allocate(E_i2(ntsteps))
   allocate(E_n2(ninter, ntsteps))
   allocate(D_ni2(ninter, ntsteps,3))
   allocate(D_fn2(nfinal, ninter, ntsteps,3))

   allocate(F_if_t_omp(1, ntsteps,3,3))
   !allocate(F_if_om_omp(ntsteps, ntsteps,3,3))

   F2_tensor = 0.0_wp

   ! here two options depending on p % KH_states_mode
   if (upper(p % KH_states_mode) .eq. "STATES") then

     allocate(R_if_om_omp(nfinal, ntsteps, ntsteps))
     allocate(R_factor(nfinal, ntsteps, ntsteps))
   
     ! compute R-factor
     call compute_F_if_om_omp_no_F(E_fc1(:,1,:), E_fi_mean, time, &
          E_i1, gamma_R, omega_out, E_nf_mean, R_if_om_omp)
     
     R_factor = (abs(R_if_om_omp))**2
     
     allocate(E_fc2(nfinal,1, ntsteps))
     
     !       do traj2 =1, ntsteps ! possibly with a stride
         
     !    ! stupid to recalcualte trajectory all the time... could be lifted out of f_e loop by saving x_new2 for each traj2
     !call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
     !         E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
     
     call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
          E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
     
     ! also E_i2 and E_n2 could be lifted out if they take a lot of time
     call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
         
     do n_e=1,ninter
       call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
     end do

     do n_e=1,ninter
       do m=1,3
         call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
       end do
     end do
         
     do f_e = 1,nfinal
           
       call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new2, E_fc2(f_e,1,:), ntsteps)  
           
       do m=1,3
         call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
       end do
           
       call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
            D_ni2(:,1,:), time, F_if_t_omp(1,:,:,:), gamma)

       !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
       
       !end do !do traj2 =1, ntsteps
       
       !call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fi_mean, time, &
       !     E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)
       

       !write(6,*) "here before big sum"
       ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
       do m1=1,3
         do m2=1,3
           do m3=1,3
             do m4=1,3
               do om_in =1, ntsteps
                 !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
                 !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
                 !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
                 F2_tensor(om_in, :, m1,m2,m3,m4) =  F2_tensor(om_in, :, m1,m2,m3,m4) + &
                      real(conjg(F_if_t_omp(1,:,m1,m2)) * F_if_t_omp(1,:,m3,m4)) * R_factor(f_e, om_in, :) 
               end do
             end do
           end do
         end do
       end do
       !write(6,*) "here after big sum"

     end do !do f_e = 1,nfinal
            
     
   else if (upper(p % KH_states_mode) .eq. "ORBS") then

     allocate(R_if_om_omp(1, ntsteps, ntsteps))
     allocate(R_factor(1, ntsteps, ntsteps))
     allocate(E_fc2(nfinal,ninter, ntsteps))
     
     do f_e = 1,nfinal
       do n_e = 1,ninter
         
         write(6,*) "f_e, n_e", f_e, n_e

         !call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
         !     E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )

         call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
              E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
         
         call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)

           ! this one really belongs here now
           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)

           call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new2, E_fc2(f_e,n_e,:), ntsteps)  
           
           do m=1,3
             call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
           end do

           call compute_F_if_om_omp_no_F_one(E_fc1(f_e,n_e,:), E_fi_mean, time, &
                E_i1, gamma_R, omega_out, E_nf_mean, R_if_om_omp(1,:,:))

           R_factor = (abs(R_if_om_omp))**2
           
           call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
                D_ni2(n_e,1,:), time, F_if_t_omp(1,:,:,:), gamma)

           
         !! perform spherical average according to J. Phys. B. 27, 4169 (1994)
         !do m1=1,3
         !  do m2=1,3
         !    do m3=1,3
         !      do m4=1,3
         !        !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
         !        !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
         !        !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
         !        F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
         !             real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
         !      end do
         !    end do
         !  end do
         !end do

           !write(6,*) "here before big sum"
           ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
           do m1=1,3
             do m2=1,3
               do m3=1,3
                 do m4=1,3
                   do om_in =1, ntsteps
                     !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
                     !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
                     !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
                     F2_tensor(om_in, :, m1,m2,m3,m4) =  F2_tensor(om_in, :, m1,m2,m3,m4) + &
                          real(conjg(F_if_t_omp(1,:,m1,m2)) * F_if_t_omp(1,:,m3,m4)) * R_factor(1, om_in, :) 
                   end do
                 end do
               end do
             end do
           end do
           !write(6,*) "here after big sum"
           
         end do! do n_e = 1,ninter
     end do! do f_e = 1,nfinal
     
   end if !  if (upper(p % KH_states_mode) .eq. "STATES") then
   
 end subroutine compute_F_SCKH_res_PES_factor_each_traj

! ! lifted out routine to compute sckh where |F(omega')|^2 is completely separated from |R(omega-omega')|^2 
! subroutine compute_F_SCKH_res_PES_factor(p, x_mom_sampl, traj, X_r, E_dyn2_inp, delta_t, mu_SI, &
!      E_i_inp, E_n_inp, E_n0, E_fn_corr, E_f_inp, D_fn_inp, D_ni_inp, E_nf_mean, E_fi_mean,time, gamma, &
!      E_i1, E_fc1, gamma_R, omega_out, F2_omp_tensor) 
!
!   use m_precision, only: wp
!   use m_constants, only: const
!   use m_sckh_params_t, only: sckh_params_t 
!   use m_upper, only : upper
!   use m_SCKH_utils, only: compute_F_if_omp_sum_n
!   !use m_SCKH_utils, only: compute_F_if_omp_one_n
!   !use m_SCKH_utils, only: compute_F_if_om_omp
!   !use m_SCKH_utils, only: compute_F_if_om_omp_one_f
!   use m_splines, only: spline_easy
!   use m_SCKH_utils, only: verlet_trajectory_xva
!   use m_SCKH_utils, only: compute_F_if_om_omp_no_F
!   
!   type(sckh_params_t), intent(in):: p
!   real(wp), intent(in):: x_mom_sampl(:,:)
!   integer, intent(in):: traj
!   !real(wp), intent(in):: x_new(:)
!   !real(wp), intent(in):: v_new(:)
!   real(wp), intent(in)::  X_r(:)
!   real(wp), intent(in):: E_dyn2_inp(:)
!   real(wp), intent(in):: delta_t
!   real(wp), intent(in):: mu_SI
!   real(wp), intent(in):: E_i_inp(:)
!   real(wp), intent(in):: E_n_inp(:,:)
!   real(wp), intent(in):: E_n0(:)
!   real(wp), intent(in):: E_f_inp(:,:)
!   real(wp), intent(in):: E_fn_corr(:,:)
!   real(wp), intent(in):: D_fn_inp(:,:,:,:)
!   real(wp), intent(in):: D_ni_inp(:,:,:)
!   real(wp), intent(in):: E_nf_mean
!   real(wp), intent(in):: E_fi_mean
!   real(wp), intent(in)::  time(:)
!   real(wp), intent(in):: gamma
!   real(wp), intent(in):: E_i1(:)
!   real(wp), intent(in):: E_fc1(:,:,:)   
!   real(wp), intent(in):: gamma_R
!   real(wp), intent(in):: omega_out(:)
!   real(wp), intent(out)::  F2_tensor(:,:,:,:,:,:)
!
!   real(wp), allocatable:: x_new2(:)
!   real(wp), allocatable:: v_new2(:)
!   real(wp), allocatable:: a_new2(:)
!   real(wp), allocatable:: E_i2(:)
!   real(wp), allocatable:: E_n2(:,:)
!   real(wp), allocatable:: E_fc2(:,:,:)
!   real(wp), allocatable:: D_fn2(:,:,:,:)
!   real(wp), allocatable:: D_ni2(:,:,:)
!
!   complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
!   complex(wp), allocatable::  F_if_om_omp(:,:,:,:)    
!
!   complex(wp), allocatable::  R_if_om_omp(:,:,:)    
!   real(wp), allocatable:: R_factor(:,:,:)
!   
!   integer:: f_e
!   integer:: n_e
!   integer:: nfinal
!   integer:: ninter
!   integer:: traj2
!   integer:: ntsteps
!   integer:: npoints_in
!   integer:: m
!   integer:: m1
!   integer:: m2
!   integer:: m3
!   integer:: m4
!   integer:: om_in
!
!   ntsteps = size(time)
!   ninter = size(E_n_inp,1)
!   nfinal = size(E_f_inp,1)
!   npoints_in = size(E_n_inp, 2)
!   
!   allocate(x_new2(ntsteps))
!   allocate(v_new2(ntsteps))
!   allocate(a_new2(ntsteps))
!
!   allocate(E_i2(ntsteps))
!   allocate(E_n2(ninter, ntsteps))
!   allocate(D_ni2(ninter, ntsteps,3))
!   allocate(D_fn2(nfinal, ninter, ntsteps,3))
!
!   allocate(F_if_t_omp(1, ntsteps,3,3))
!   !allocate(F_if_om_omp(ntsteps, ntsteps,3,3))
!
!   allocate(R_if_om_omp(nfinal, ntsteps, ntsteps))
!   allocate(R_factor(nfinal, ntsteps, ntsteps))
!   
!   F2_tensor = 0.0_wp
!
!   ! here two options depending on p % KH_states_mode
!   if (upper(p % KH_states_mode) .eq. "STATES") then
!
!     ! compute R-factor
!     !call compute_F_if_om_omp_no_F(E_fc1(:,1,:), E_fi_mean, time, &
!     !     E_i1, gamma_R, omega_out, E_nf_mean, R_if_om_omp)
!     !
!     !R_factor = (abs(R_if_om_omp))**2
!     
!     allocate(E_fc2(nfinal,1, ntsteps))
!     
!     !       do traj2 =1, ntsteps ! possibly with a stride
!         
!     !    ! stupid to recalcualte trajectory all the time... could be lifted out of f_e loop by saving x_new2 for each traj2
!     !call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
!     !         E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
!     
!     call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
!          E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
!     
!     ! also E_i2 and E_n2 could be lifted out if they take a lot of time
!     call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
!         
!     do n_e=1,ninter
!       call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
!     end do
!
!     do n_e=1,ninter
!       do m=1,3
!         call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
!       end do
!     end do
!         
!     do f_e = 1,nfinal
!           
!       call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new2, E_fc2(f_e,1,:), ntsteps)  
!           
!       do m=1,3
!         call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
!       end do
!           
!       call compute_F_if_omp_sum_n(E_n2(:,:), E_fc2(f_e,1,:), E_nf_mean, D_fn2(f_e,1,:,:), &
!            D_ni2(:,1,:), time, F_if_t_omp(1,:,:,:), gamma)
!
!       
!       !F_if_t_omp(traj2,:,:,:) = F_if_t_omp(traj2,:,:,:) + F_if_t_omp_tmp
!       
!       !end do !do traj2 =1, ntsteps
!       
!       !call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,1,:), E_fi_mean, time, &
!       !     E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) !F_if_om_omp)
!       
!
!       !write(6,*) "here before big sum"
!       ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
!       do m1=1,3
!         do m2=1,3
!           do m3=1,3
!             do m4=1,3
!               do om_in =1, ntsteps
!                 !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
!                 !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
!                 !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
!                 F2_tensor(om_in, :, m1,m2,m3,m4) =  F2_tensor(om_in, :, m1,m2,m3,m4) + &
!                      real(conjg(F_if_t_omp(1,:,m1,m2)) * F_if_t_omp(1,:,m3,m4)) * R_factor(f_e, om_in, :) 
!               end do
!             end do
!           end do
!         end do
!       end do
!       !write(6,*) "here after big sum"
!
!     end do !do f_e = 1,nfinal
!            
!     
!   else if (upper(p % KH_states_mode) .eq. "ORBS") then
!!
!!     allocate(E_fc2(nfinal,ninter, ntsteps))
!!
!!     do f_e = 1,nfinal
!!       do n_e = 1,ninter
!!         
!!         !fn_e = ninter * (f_e -1) + n_e  
!!         
!!!         do traj2 =1, ntsteps ! possibly with a stride              
!!           
!!           call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
!!                E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
!!           
!!           call spline_easy(X_r, E_i_inp, npoints_in, x_new2, E_i2, ntsteps)
!!           
!!           ! this one really belongs here now
!!           call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new2, E_n2(n_e,:), ntsteps)
!!           
!!           call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new2, E_fc2(f_e,n_e,:), ntsteps)  
!!           
!!           do m=1,3
!!             call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new2, D_fn2(f_e,1,:,m) , ntsteps)  
!!             call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new2, D_ni2(n_e,:,m) , ntsteps)  
!!           end do
!!           
!!           call compute_F_if_omp_one_n(E_n2(n_e,:), E_fc2(f_e,n_e,:), E_nf_mean, D_fn2(f_e,1,:,:), &
!!                D_ni2(n_e,1,:), time, F_if_t_omp(traj2,:,:,:), gamma)
!!           
!! !        end do ! do traj2 =1, ntsteps
!!         
!!         call compute_F_if_om_omp_one_f(F_if_t_omp(:,:,:,:), E_fc1(f_e,n_e,:), E_fi_mean, time, &
!!              E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp(:,:,:,:)) 
!!         
!!         ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
!!         do m1=1,3
!!           do m2=1,3
!!             do m3=1,3
!!               do m4=1,3
!!                 !lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m1)) * F_if_om_omp(:,:,m2,m2))
!!                 !lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m1,m2))
!!                 !lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m2,m1))
!!                 F2_tensor(:, :, m1,m2,m3,m4) =  F2_tensor(:, :, m1,m2,m3,m4) + &
!!                      real(conjg(F_if_om_omp(:,:,m1,m2)) * F_if_om_omp(:,:,m3,m4))
!!               end do
!!             end do
!!           end do
!!         end do
!!
!!         
!!       end do! do n_e = 1,ninter
!!     end do! do f_e = 1,nfinal
!     
!   end if !  if (upper(p % KH_states_mode) .eq. "STATES") then
!   
! end subroutine compute_F_SCKH_res_PES_factor



 subroutine compute_traj_and_spline(p, x_0, v_0, X_r, E_dyn_inp, delta_t, mu_SI, &
      x_new, v_new, a_new, E_i_inp, E_n_inp, E_n0, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp, &
      E_i1, E_n1, E_fc1, D_ni1, D_fn1 )

   use m_precision, only: wp
   use m_constants, only: const
   use m_sckh_params_t, only: sckh_params_t 
   use m_upper, only : upper
   use m_splines, only: spline_easy
   use m_SCKH_utils, only: verlet_trajectory_xva
   
   type(sckh_params_t), intent(in):: p
   real(wp), intent(in):: x_0
   real(wp), intent(in):: v_0
   real(wp), intent(in)::  X_r(:)
   real(wp), intent(in):: E_dyn_inp(:)
   real(wp), intent(in):: delta_t
   real(wp), intent(in):: mu_SI
   real(wp), intent(out):: x_new(:)
   real(wp), intent(out):: v_new(:)
   real(wp), intent(out):: a_new(:)
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_n_inp(:,:)
   real(wp), intent(in):: E_n0(:)
   real(wp), intent(in):: E_f_inp(:,:)
   real(wp), intent(in):: E_fn_corr(:,:)
   real(wp), intent(in):: D_fn_inp(:,:,:,:)
   real(wp), intent(in):: D_ni_inp(:,:,:)
   real(wp), intent(out):: E_i1(:)
   real(wp), intent(out):: E_n1(:,:)
   real(wp), intent(out):: E_fc1(:,:,:)
   real(wp), intent(out):: D_ni1(:,:,:)
   real(wp), intent(out):: D_fn1(:,:,:,:)

   integer:: f_e
   integer:: n_e
   integer:: m
   integer:: ninter
   integer:: nfinal
   integer:: npoints_in
   integer:: ntsteps
   
   ntsteps = size(E_i1)
   npoints_in = size(X_r)
   ninter = size(E_n_inp,1)
   nfinal = size(E_f_inp,1)
   
   !call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
   !     E_dyn_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )

   call verlet_trajectory_xva(x_0, v_0, X_r, &
        E_dyn_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )
   
   call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
   do n_e=1,ninter
     call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new, E_n1(n_e,:), ntsteps)  
   end do
   
   do n_e=1,ninter
     do m=1,3
       call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new, D_ni1(n_e,:,m) , ntsteps)  
     end do
   end do
   
   do f_e=1,nfinal
     do m=1,3
       call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new, D_fn1(f_e,1,:,m) , ntsteps)  
     end do
   end do
   
   !do i=1,nfinal
   !  call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_f1(i,:), ntsteps)  
   !end do
   
   if (upper(p % KH_states_mode) .eq. "STATES") then
     do f_e=1,nfinal
       call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new, E_fc1(f_e,1,:), ntsteps)  
     end do
     
   else if (upper(p % KH_states_mode) .eq. "ORBS") then
     do f_e=1,nfinal
       do n_e=1,ninter
         call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new, E_fc1(f_e,n_e,:), ntsteps)  
       end do
     end do
     
   end if
   
   
 end subroutine compute_traj_and_spline

 subroutine compute_E_means_and_omegas(E_i_inp, E_i1, E_n1, E_fc1, &
      E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: get_omega_reordered_fftw
   
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_i1(:)
   real(wp), intent(in):: E_n1(:,:)
   real(wp), intent(in):: E_fc1(:,:,:)
   real(wp), intent(out):: E_ni_mean
   real(wp), intent(out):: E_fi_mean
   real(wp), intent(out):: E_nf_mean
   real(wp), intent(in):: time_l
   real(wp), intent(out):: omega_in(:)
   real(wp), intent(out):: omega_out(:)

   integer, dimension(1)::ind  
   
!   ! first time, compute the mean transition energy, and frequency
!   !if (traj .eq. 1) then
!   ind = minloc(E_i_inp)
!   
!   E_ni_mean =  E_n1(1, ind(1)) - E_i1(ind(1))
!   !E_ni_mean =  E_n1(size(E_fc1,1), ind(1)) - E_i1(ind(1))
!   write(6,*) "E_ni_mean", E_ni_mean
!   
!   E_fi_mean =  E_fc1(1, 1, ind(1)) - E_i1(ind(1))
!   !E_fi_mean =  E_fc1(size(E_fc1,1), 1, ind(1)) - E_i1(ind(1))
!   write(6,*) "E_fi_mean", E_fi_mean
!   
!   E_nf_mean =  E_n1(1, ind(1)) - E_fc1(1,1,ind(1))
!   !E_nf_mean =  E_n1(size(E_n1,1), ind(1)) - E_fc1(size(E_fc1,1),1,ind(1))
!   write(6,*) "E_nf_mean", E_nf_mean
!   
!   call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
!   omega_in = omega_in + E_ni_mean !E_nf_mean + E_fi_mean
!   
!   call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
!   omega_out = omega_out + E_nf_mean
   
   call compute_E_means_and_omegas_one(E_i_inp, E_i1, E_n1(1,:), E_fc1(1,1,:), &
        E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)

  end subroutine compute_E_means_and_omegas

 subroutine compute_E_means_and_omegas_one(E_i_inp, E_i1, E_n1, E_fc1, &
      E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: get_omega_reordered_fftw
   
   real(wp), intent(in):: E_i_inp(:)
   real(wp), intent(in):: E_i1(:)
   real(wp), intent(in):: E_n1(:)
   real(wp), intent(in):: E_fc1(:)
   real(wp), intent(out):: E_ni_mean
   real(wp), intent(out):: E_fi_mean
   real(wp), intent(out):: E_nf_mean
   real(wp), intent(in):: time_l
   real(wp), intent(out):: omega_in(:)
   real(wp), intent(out):: omega_out(:)

   integer, dimension(1)::ind  
   
   ! first time, compute the mean transition energy, and frequency
   !if (traj .eq. 1) then
   ind = minloc(E_i_inp)

   write(6,*) "ind", ind
   
   E_ni_mean =  E_n1(ind(1)) - E_i1(ind(1))
   !E_ni_mean =  E_n1(size(E_fc1,1), ind(1)) - E_i1(ind(1))
   write(6,*) "E_ni_mean", E_ni_mean
   
   E_fi_mean =  E_fc1(ind(1)) - E_i1(ind(1))
   !E_fi_mean =  E_fc1(size(E_fc1,1), 1, ind(1)) - E_i1(ind(1))
   write(6,*) "E_fi_mean", E_fi_mean
   
   E_nf_mean =  E_n1(ind(1)) - E_fc1(ind(1))
   !E_nf_mean =  E_n1(size(E_n1,1), ind(1)) - E_fc1(size(E_fc1,1),1,ind(1))
   write(6,*) "E_nf_mean", E_nf_mean
   
   call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
   omega_in = omega_in + E_ni_mean !E_nf_mean + E_fi_mean
   
   call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
   omega_out = omega_out + E_nf_mean
   
   
 end subroutine compute_E_means_and_omegas_one

 
 

end module m_SCKH_resonant_PES
