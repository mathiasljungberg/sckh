module m_SCKH_resonant_PES_FC
  implicit none

contains

  subroutine calculate_SCKH_res_PES_FC(p)
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
    use m_SCKH_resonant_PES, only: compute_E_means_and_omegas_one
    use m_rixs_io, only: write_RIXS_spectra
    use m_rixs_io, only: write_RIXS_map
    
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
    real(wp):: gamma_R
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

    !complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
    !complex(wp), allocatable::  F_if_om_omp(:,:,:,:,:)

    complex(wp), allocatable::  F_tmp(:,:,:,:)
    complex(wp), allocatable::  F_tmp2(:,:,:,:)

    !complex(wp), allocatable::  F_if_omp(:,:,:,:)    
    !complex(wp), allocatable::  F_if_omp_tmp(:,:,:,:)    

    !complex(wp), allocatable::  R_if_om_omp(:,:)    
    !real(wp), allocatable:: R_factor(:,:)

    !complex(wp), allocatable::  R_if_om_omp(:,:,:)    
    !real(wp), allocatable:: R_factor(:,:,:)

    !real(wp), allocatable:: lambda_F(:,:,:)
    !real(wp), allocatable:: lambda_G(:,:,:)
    !real(wp), allocatable:: lambda_H(:,:,:)

    real(wp), allocatable:: lambda_F(:,:)
    real(wp), allocatable:: lambda_G(:,:)
    real(wp), allocatable:: lambda_H(:,:)

    !real(wp), allocatable:: lambda_F_tmp(:,:)
    !real(wp), allocatable:: lambda_G_tmp(:,:)
    !real(wp), allocatable:: lambda_H_tmp(:,:)

    real(wp), allocatable:: lambda_lp(:,:)
    real(wp), allocatable:: lambda_ln(:,:)
    real(wp), allocatable:: lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)
    
    real(wp):: E_ni
    real(wp):: E_nf
    real(wp):: E_i_min, E_i_av
    
    !real(wp), allocatable::  F2_tensor(:,:,:,:,:,:)


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
    write(6,*) "p % sckh_pes_dyn_mode  = ", p % sckh_pes_dyn_mode
    write(6,*) "p % sckh_alt_mode  = ", p % sckh_alt_mode
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
      E_dyn2_inp = E_n_inp(1,:) + E_n0 ! add reference state to dynamics 
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

    ind = minloc(E_i_inp)
    E_i_min = E_i_inp(ind(1))
    E_i_av = sum(E_i_inp * c_i(:,1)**2) - E_i_min! ZPE

    write(6,*) "Average potenital energy, from minimum", E_i_av
    
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

    !allocate(lambda_F(nfinal_tot, n_omega_in, n_omega_out))
    !allocate(lambda_G(nfinal_tot, n_omega_in, n_omega_out))
    !allocate(lambda_H(nfinal_tot, n_omega_in, n_omega_out))

    allocate(lambda_F(n_omega_in, n_omega_out))
    allocate(lambda_G(n_omega_in, n_omega_out))
    allocate(lambda_H(n_omega_in, n_omega_out))

    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))
    allocate(sigma_tmp(n_omega_in, n_omega_out))

    !allocate(F_if_om_omp(nfinal_tot, n_omega_in, n_omega_out,3,3))
    allocate(F_tmp(n_omega_in, n_omega_out,3,3))
    allocate(F_tmp2(n_omega_in, n_omega_out,3,3))

    !
    ! Here starts the program proper
    ! 

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do
    
    ! get mean energies and omegas

    !call compute_E_means_and_omegas_one(E_i_inp, E_i_inp, E_n_inp(1,:)+ E_n0, &
    !     E_f_inp(1,:), &
    !     E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)    

    if (p % use_E_mean) then
      write(6,*) "p % use_E_mean", p % use_E_mean
      E_nf_mean = p % E_nf_mean
      E_ni_mean = p % E_ni_mean
    else
      write(6,*) "p % use_E_mean", p % use_E_mean
      E_nf_mean =   E_n_inp(1,ind(1))+ E_n0(ind(1)) - E_f_inp(1,ind(1)) !E_n1(1,ind(1)) - E_fc1(1,1,ind(1))
      E_ni_mean =   E_n_inp(1,ind(1))+ E_n0(ind(1)) - E_i_inp(ind(1))   !E_n1(1,ind(1)) - E_i1(ind(1))
    end if
    E_fi_mean = E_ni_mean -E_nf_mean

    write(6,*) "E_nf_mean", E_nf_mean
    write(6,*) "E_ni_mean", E_ni_mean
    write(6,*) "E_fi_mean", E_fi_mean
    
    call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
    omega_in = omega_in + E_ni_mean 
    
    call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
    omega_out = omega_out + E_nf_mean

    
    !write(6,*) "omega_in", omega_in
    !write(6,*) "omega_out", omega_out
    
    !F_if_om_omp = 0.0_wp
    
    do traj=1, npoints_x_mom_sampl

      write(6,*) "traj ", traj, "out of ", npoints_x_mom_sampl
      
      ! new loop here over intermediate states i order to do separate dynamics
      do n_e=1,ninter

        if(upper(p % sckh_pes_dyn_mode) .eq. "SEPARATE") then

          write(6,*) "dynamics mode: SEPARATE"
            ! compute mean energies on E_dyn_inp that will match with the other routines
            call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
                 (E_n_inp(n_e,:) + E_n0) * const % eV, delta_t, mu_SI, x_new, v_new, a_new, &
                 p % use_abs_bc, p % abs_bc_max_x * 1.0d-10)
          
        else if (upper(p % sckh_pes_dyn_mode) .eq. "SINGLE") then

          ! only do dynamics at first n
          if(n_e .eq. 1) then

            write(6,*) "dynamics mode: SINGLE"
            
            ! compute mean energies on E_dyn_inp that will match with the other routines
            call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
                 E_dyn2_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new, &
                 p % use_abs_bc, p % abs_bc_max_x * 1.0d-10)
          end if

        else
          write(6,*) "p % sckh_pes_dyn_mode must be 'SINGLE' or 'SEPARATE'"
          
        end if ! if(upper(p % dynamics_mode .eq. "SEPARATE")) then

        ! a bit unnecessary splines but well...
        call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
        !do n_e=1,ninter
        call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new, E_n1(n_e,:), ntsteps)  
        !  end do
        
        ! options for dipole moments in XAS
        if (upper(p % dipole_mode) .eq. "DIPOLE") then
          !do n_e=1,ninter
          do m=1,3
            call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new, D_ni1(n_e,:,m) , ntsteps)  
            ! end do
        end do
      else if(upper(p % dipole_mode) .eq. "FC") then
        D_ni1 = 1.0_wp
      else if(upper(p % dipole_mode) .eq. "DIPOLE_X0") then
        do i=1, npoints_in !p % nstates
          D_ni1(:,i,:) = D_ni_inp(:,ind(1),:) 
        end do
      else
        write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
        stop
      end if


      ! options for dipole moments in XES
      if (upper(p % dipole_mode) .eq. "DIPOLE") then        
        do f_e=1,nfinal
          do m=1,3
            call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new, D_fn1(f_e,1,:,m) , ntsteps)  
          end do
        end do
      else if(upper(p % dipole_mode) .eq. "FC") then
        D_fn1 = 1.0_wp
      else if(upper(p % dipole_mode) .eq. "DIPOLE_X0") then
        do i=1, npoints_in !p % nstates
          D_fn1(:,1,i,:) = D_fn_inp(:,1,ind(1),:) 
        end do
      else
        write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
        stop
      end if

      !F_if_om_omp = 0.0_wp
      
      if (upper(p % KH_states_mode) .eq. "STATES") then
        
        do f_e=1,nfinal
          call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new, E_fc1(f_e,1,:), ntsteps)  
        end do
        
        do f_e=1,nfinal

          F_tmp =0.0_wp
          
        !  do n_e=1,ninter
            
            !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV

            if(p % include_ZPE) then
              E_ni = E_n1(n_e,1)- (E_i_min + E_i_av)
              !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            else
              E_ni = E_n1(n_e,1)-E_i1(1) 
              !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            end if
            
            E_nf = E_n1(n_e,1)- E_fc1(f_e,1,1)
            

            ! change the sckh_alt_mode NORMAL to what before was ALT4, the old normal will be ALT5 
            if(upper(p % sckh_alt_mode) .eq. "NORMAL") then
              call compute_F_FC_if_om_omp_alt4(E_ni, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT") then
              call compute_F_FC_if_om_omp_alt(E_ni, E_nf, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, omega_out, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT2") then
              call compute_F_FC_if_om_omp_alt2(E_ni, E_nf, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, omega_out, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT3") then
              call compute_F_FC_if_om_omp_alt3(E_ni, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, omega_out, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT4") then
              call compute_F_FC_if_om_omp_alt4(E_ni, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT5") then
              call compute_F_FC_if_om_omp(E_ni, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            end if
              
            !write(6,*) "Here... 3"
            !F_if_om_omp(f_e,:,:,:,:) = F_if_om_omp(f_e,:,:,:,:) + F_tmp(:,:,:,:)  
            
            ! internal sum only over intermediate states
           ! F_tmp = F_tmp + F_tmp2
            
          !end do
          
            ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
            do m1=1,3
              do m2=1,3
                lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_tmp(:,:,m1,m1)) * F_tmp(:,:,m2,m2))
                lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m1,m2))
                lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m2,m1))
              end do
            end do

          end do
        
      else if (upper(p % KH_states_mode) .eq. "ORBS") then

        do f_e=1,nfinal
          !do n_e=1,ninter

            fn_e = ninter * (f_e -1) + n_e  

            write(6,*) "f_e, n_e ", f_e, n_e, " out of ", nfinal, ninter
            
            call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new, E_fc1(f_e,n_e,:), ntsteps)  
            
            if(p % include_ZPE) then
              E_ni = E_n1(n_e,1)- (E_i_min + E_i_av)
            else
              E_ni = E_n1(n_e,1)-E_i1(1) !+ (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            end if

            E_nf = E_n1(n_e,1)- E_fc1(f_e,n_e,1)
            
            ! change the sckh_alt_mode NORMAL to what before was ALT4, the old normal will be ALT5 
            if(upper(p % sckh_alt_mode) .eq. "NORMAL") then
              call compute_F_FC_if_om_omp_alt4(E_ni, E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT") then
              call compute_F_FC_if_om_omp_alt(E_ni, E_nf,E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, omega_out, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT2") then
              call compute_F_FC_if_om_omp_alt2(E_ni, E_nf,E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, omega_out, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT3") then
              call compute_F_FC_if_om_omp_alt3(E_ni, E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, omega_out, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT4") then
              call compute_F_FC_if_om_omp_alt4(E_ni, E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            else if (upper(p % sckh_alt_mode) .eq. "ALT5") then
              call compute_F_FC_if_om_omp(E_ni, E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                   D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            end if
              
              !write(6,*) "Here... 3"
            !F_if_om_omp(fn_e,:,:,:,:) = F_if_om_omp(fn_e,:,:,:,:) + F_tmp(:,:,:,:)  

            ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
            do m1=1,3
              do m2=1,3
                lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_tmp(:,:,m1,m1)) * F_tmp(:,:,m2,m2))
                lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m1,m2))
                lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m2,m1))
              end do
            end do
            
          !end do
          end do ! do f_e=1,nfinal
        
        end if ! if (upper(p % KH_states_mode) .eq. "STATES") then
      
!      ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
!      do m1=1,3
!        do m2=1,3
!          lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m1)) * F_if_om_omp(:,:,:,m2,m2))
!          lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m1,m2))
!          lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m2,m1))
!        end do
!      end do
      
      end do ! do n_e=1,ninter
    end do ! do traj=1, npoints_x_mom_sampl

!      ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
!      do m1=1,3
!        do m2=1,3
!          lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m1)) * F_if_om_omp(:,:,:,m2,m2))
!          lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m1,m2))
!          lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m2,m1))
!        end do
!      end do


    
!    ! averages according to J. Phys. B. 27, 4169 (1994) 
!    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
!    lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
!    lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
!    lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    ! averages according to J. Phys. B. 27, 4169 (1994) 
    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
    lambda_lp = 2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H
    lambda_ln = -1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H
    lambda_cp = -2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H

    ! write unbroadened spectrum (well, only lifetime and extra lifetime on omega' here)
    call write_RIXS_map(omega_in, omega_out, lambda_lp, lambda_ln, lambda_cp, &
         p % outfile, p % kh_print_stride)

    write(6,*) "Entering convolute_incoming, broadening ",  upper(p % broadening_func_inc)
    
    if(gamma_inc .gt. 1d-5) then
      sigma_tmp = lambda_lp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_lp, 1,  upper(p % broadening_func_inc))
      sigma_tmp = lambda_ln
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(p % broadening_func_inc))
      sigma_tmp = lambda_cp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_cp, 1, upper(p % broadening_func_inc))
    end if
    
    write(6,*) "Done"
    write(6,*) "Entering convolute_instrumental, broadening  ",  upper(p % broadening_func_instr) 
    
    if(gamma_instr .gt. 1d-5) then
      sigma_tmp = lambda_lp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_ln
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_cp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(p % broadening_func_instr))
    end if
    
    write(6,*) "Done"

    call write_RIXS_spectra(omega_in, omega_out, lambda_lp, lambda_ln, lambda_cp, p % outfile, p %kh_print_stride)

!    ! write spectra to individual files
!    do j=1, n_omega_in, p % kh_print_stride
!      
!      file="_sigma_"
!      write(string,'(F6.2)') omega_in(j)   
!
!      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!      
!      ifile = get_free_handle()
!      open(ifile,file=file,status='unknown')
!      
!      do i=1, n_omega_out, p % kh_print_stride
!        
!        write(ifile,'(4ES18.10)') omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
!      end do
!      
!      close(ifile) 
!      
!   end do
!
!   if(.true.) then
!   
!    ! write spectra to file
!   file="_sigma_all"
!   !write(string,'(F6.2)') omega_in(j)   
!   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
!   
!   ifile = get_free_handle()
!   open(ifile,file=file,status='unknown')
!
!   do j=1, n_omega_in, p % kh_print_stride
!     do i=1, n_omega_out, p % kh_print_stride
!       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
!     end do
!     write(ifile, *) 
!   end do
!   
!   close(ifile) 
!
!
!   ! write spectra to file
!   file="_sigma_all_nogrid"
!   !write(string,'(F6.2)') omega_in(j)   
!   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
!   
!   ifile = get_free_handle()
!   open(ifile,file=file,status='unknown')
!
!   do j=1, n_omega_in, p % kh_print_stride
!     do i=1, n_omega_out, p % kh_print_stride
!       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
!     end do
!     write(ifile, *) 
!     write(ifile, *) 
!   end do
!   
!   close(ifile) 
!
!   ! write spectra summed over incoming frequencies
!   file="_nonres"
!   !write(string,'(F6.2)') omega_in(j)   
!   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
!   
!   ifile = get_free_handle()
!   open(ifile,file=file,status='unknown')
!
!   do i=1, n_omega_out, p % kh_print_stride
!     write(ifile,'(5ES18.10)') omega_out(i), sum(lambda_lp(:,i)), sum(lambda_ln(:,i)), sum(lambda_cp(:,i))
!   end do
!   
!   close(ifile) 
!
!   ! write spectra summed over outgoing frequencies
!   file="_xas"
!   !write(string,'(F6.2)') omega_in(j)   
!   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
!   
!   ifile = get_free_handle()
!   open(ifile,file=file,status='unknown')
!   
!   do i=1, n_omega_in, p % kh_print_stride
!     write(ifile,'(5ES18.10)') omega_in(i), sum(lambda_lp(i,:)), sum(lambda_ln(i,:)), sum(lambda_cp(i,:))
!   end do
!   
!   close(ifile) 
!
!   
! end if! if(.false.) then
   
    
 end subroutine calculate_SCKH_res_PES_FC


 subroutine calculate_SCKH_res_PES_FC_old(p)
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
    use m_SCKH_resonant_PES, only: compute_E_means_and_omegas_one
    
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
    real(wp):: gamma_R
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

    !complex(wp), allocatable::  F_if_t_omp(:,:,:,:)    
    !complex(wp), allocatable::  F_if_om_omp(:,:,:,:,:)

    complex(wp), allocatable::  F_tmp(:,:,:,:)
    complex(wp), allocatable::  F_tmp2(:,:,:,:)

    !complex(wp), allocatable::  F_if_omp(:,:,:,:)    
    !complex(wp), allocatable::  F_if_omp_tmp(:,:,:,:)    

    !complex(wp), allocatable::  R_if_om_omp(:,:)    
    !real(wp), allocatable:: R_factor(:,:)

    !complex(wp), allocatable::  R_if_om_omp(:,:,:)    
    !real(wp), allocatable:: R_factor(:,:,:)

    !real(wp), allocatable:: lambda_F(:,:,:)
    !real(wp), allocatable:: lambda_G(:,:,:)
    !real(wp), allocatable:: lambda_H(:,:,:)

    real(wp), allocatable:: lambda_F(:,:)
    real(wp), allocatable:: lambda_G(:,:)
    real(wp), allocatable:: lambda_H(:,:)

    !real(wp), allocatable:: lambda_F_tmp(:,:)
    !real(wp), allocatable:: lambda_G_tmp(:,:)
    !real(wp), allocatable:: lambda_H_tmp(:,:)

    real(wp), allocatable:: lambda_lp(:,:)
    real(wp), allocatable:: lambda_ln(:,:)
    real(wp), allocatable:: lambda_cp(:,:)
    real(wp), allocatable:: sigma_tmp(:,:)
    
    real(wp):: E_ni
    real(wp):: E_i_min, E_i_av
    
    !real(wp), allocatable::  F2_tensor(:,:,:,:,:,:)


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
      E_dyn2_inp = E_n_inp(1,:) + E_n0 ! add reference state to dynamics 
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

    ind = minloc(E_i_inp)
    E_i_min = E_i_inp(ind(1))
    E_i_av = sum(E_i_inp * c_i(:,1)**2) - E_i_min! ZPE

    write(6,*) "Average potenital energy, from minimum", E_i_av
    
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

    !allocate(lambda_F(nfinal_tot, n_omega_in, n_omega_out))
    !allocate(lambda_G(nfinal_tot, n_omega_in, n_omega_out))
    !allocate(lambda_H(nfinal_tot, n_omega_in, n_omega_out))

    allocate(lambda_F(n_omega_in, n_omega_out))
    allocate(lambda_G(n_omega_in, n_omega_out))
    allocate(lambda_H(n_omega_in, n_omega_out))

    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))
    allocate(sigma_tmp(n_omega_in, n_omega_out))

    !allocate(F_if_om_omp(nfinal_tot, n_omega_in, n_omega_out,3,3))
    allocate(F_tmp(n_omega_in, n_omega_out,3,3))
    allocate(F_tmp2(n_omega_in, n_omega_out,3,3))

    !
    ! Here starts the program proper
    ! 

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do
    
    call compute_E_means_and_omegas_one(E_i_inp, E_i_inp, E_n_inp(1,:)+ E_n0, &
         E_f_inp(1,:), &
         E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)    

    !write(6,*) "omega_in", omega_in
    !write(6,*) "omega_out", omega_out
    
    !F_if_om_omp = 0.0_wp
    
    do traj=1, npoints_x_mom_sampl

      write(6,*) "traj ", traj, "out of ", npoints_x_mom_sampl
      
      ! compute mean energies on E_dyn_inp that will match with the other routines
      call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
           E_dyn2_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )
      
      call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
      do n_e=1,ninter
        call spline_easy(X_r, E_n_inp(n_e,:) + E_n0, npoints_in, x_new, E_n1(n_e,:), ntsteps)  
      end do
                
      ! options for dipole moments in XAS
      if (upper(p % dipole_mode) .eq. "DIPOLE") then
        do n_e=1,ninter
          do m=1,3
            call spline_easy(X_r, D_ni_inp(n_e,:,m) , npoints_in, x_new, D_ni1(n_e,:,m) , ntsteps)  
          end do
        end do
      else if(upper(p % dipole_mode) .eq. "FC") then
        D_ni1 = 1.0_wp
      else if(upper(p % dipole_mode) .eq. "DIPOLE_X0") then
        do i=1, npoints_in !p % nstates
          D_ni1(:,i,:) = D_ni_inp(:,ind(1),:) 
        end do
      else
        write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
        stop
      end if


      ! options for dipole moments in XES
      if (upper(p % dipole_mode) .eq. "DIPOLE") then        
        do f_e=1,nfinal
          do m=1,3
            call spline_easy(X_r, D_fn_inp(f_e,1,:,m) , npoints_in, x_new, D_fn1(f_e,1,:,m) , ntsteps)  
          end do
        end do
      else if(upper(p % dipole_mode) .eq. "FC") then
        D_fn1 = 1.0_wp
      else if(upper(p % dipole_mode) .eq. "DIPOLE_X0") then
        do i=1, npoints_in !p % nstates
          D_fn1(:,1,i,:) = D_fn_inp(:,1,ind(1),:) 
        end do
      else
        write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
        stop
      end if


        
      !F_if_om_omp = 0.0_wp
      
      if (upper(p % KH_states_mode) .eq. "STATES") then

        do f_e=1,nfinal
          call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new, E_fc1(f_e,1,:), ntsteps)  
        end do
        
        do f_e=1,nfinal

          F_tmp =0.0_wp
          
          do n_e=1,ninter
            
            !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV

            if(p % include_ZPE) then
              E_ni = E_n1(n_e,1)- (E_i_min + E_i_av)
            else
              E_ni = E_n1(n_e,1)-E_i1(1) 
              !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            end if
            
            
            call compute_F_FC_if_om_omp(E_ni, E_n1(n_e,:), E_fc1(f_e,1,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                 D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp2) !F_if_om_omp(:,om_in,:,:,:), gamma)
            
            !write(6,*) "Here... 3"
            !F_if_om_omp(f_e,:,:,:,:) = F_if_om_omp(f_e,:,:,:,:) + F_tmp(:,:,:,:)  
            
            ! internal sum only over intermediate states
            F_tmp = F_tmp + F_tmp2
            
          end do
          
            ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
            do m1=1,3
              do m2=1,3
                lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_tmp(:,:,m1,m1)) * F_tmp(:,:,m2,m2))
                lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m1,m2))
                lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m2,m1))
              end do
            end do

        end do
        
      else if (upper(p % KH_states_mode) .eq. "ORBS") then

        do f_e=1,nfinal
          do n_e=1,ninter

            fn_e = ninter * (f_e -1) + n_e  

            write(6,*) "f_e, n_e ", f_e, n_e, " out of ", nfinal, ninter
            
            call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new, E_fc1(f_e,n_e,:), ntsteps)  
            
            if(p % include_ZPE) then
              E_ni = E_n1(n_e,1)- (E_i_min + E_i_av)
            else
              E_ni = E_n1(n_e,1)-E_i1(1) !+ (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            end if
            
            call compute_F_FC_if_om_omp(E_ni, E_n1(n_e,:), E_fc1(f_e,n_e,:), E_ni_mean, E_nf_mean, E_fi_mean, &
                 D_fn1(f_e,1,:,:), D_ni1(n_e,1,:), omega_in, time, gamma, gamma_R, F_tmp) !F_if_om_omp(:,om_in,:,:,:), gamma)
            
            !write(6,*) "Here... 3"
            !F_if_om_omp(fn_e,:,:,:,:) = F_if_om_omp(fn_e,:,:,:,:) + F_tmp(:,:,:,:)  

            ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
            do m1=1,3
              do m2=1,3
                lambda_F(:, :) = lambda_F(:, :) +  real(conjg(F_tmp(:,:,m1,m1)) * F_tmp(:,:,m2,m2))
                lambda_G(:, :) = lambda_G(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m1,m2))
                lambda_H(:, :) = lambda_H(:, :) +  real(conjg(F_tmp(:,:,m1,m2)) * F_tmp(:,:,m2,m1))
              end do
            end do
            
          end do
        end do
        
      end if
      
!      ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
!      do m1=1,3
!        do m2=1,3
!          lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m1)) * F_if_om_omp(:,:,:,m2,m2))
!          lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m1,m2))
!          lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m2,m1))
!        end do
!      end do

      
    end do ! do traj=1, npoints_x_mom_sampl

!      ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
!      do m1=1,3
!        do m2=1,3
!          lambda_F(:, :, :) = lambda_F(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m1)) * F_if_om_omp(:,:,:,m2,m2))
!          lambda_G(:, :, :) = lambda_G(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m1,m2))
!          lambda_H(:, :, :) = lambda_H(:, :, :) +  real(conjg(F_if_om_omp(:,:,:,m1,m2)) * F_if_om_omp(:,:,:,m2,m1))
!        end do
!      end do


    
!    ! averages according to J. Phys. B. 27, 4169 (1994) 
!    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
!    lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
!    lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
!    lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    ! averages according to J. Phys. B. 27, 4169 (1994) 
    ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
    lambda_lp = 2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H
    lambda_ln = -1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H
    lambda_cp = -2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H


    write(6,*) "Entering convolute_incoming, broadening ",  upper(p % broadening_func_inc)
    
    if(gamma_inc .gt. 1d-5) then
      sigma_tmp = lambda_lp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_lp, 1,  upper(p % broadening_func_inc))
      sigma_tmp = lambda_ln
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(p % broadening_func_inc))
      sigma_tmp = lambda_cp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_cp, 1, upper(p % broadening_func_inc))
    end if
    
    write(6,*) "Done"
    write(6,*) "Entering convolute_instrumental, broadening  ",  upper(p % broadening_func_instr) 
    
    if(gamma_instr .gt. 1d-5) then
      sigma_tmp = lambda_lp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_ln
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(p % broadening_func_instr))
      sigma_tmp = lambda_cp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(p % broadening_func_instr))
    end if
    
    write(6,*) "Done"

    ! write spectra to individual files
    do j=1, n_omega_in, p % kh_print_stride
      
      file="_sigma_"
      write(string,'(F6.2)') omega_in(j)   

      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1, n_omega_out, p % kh_print_stride
        
        write(ifile,'(4ES18.10)') omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do
      
      close(ifile) 
      
   end do

   if(.true.) then
   
    ! write spectra to file
   file="_sigma_all"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in, p % kh_print_stride
     do i=1, n_omega_out, p % kh_print_stride
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

   do j=1, n_omega_in, p % kh_print_stride
     do i=1, n_omega_out, p % kh_print_stride
       write(ifile,'(5ES18.10)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
     write(ifile, *) 
   end do
   
   close(ifile) 

   ! write spectra summed over incoming frequencies
   file="_nonres"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do i=1, n_omega_out, p % kh_print_stride
     write(ifile,'(5ES18.10)') omega_out(i), sum(lambda_lp(:,i)), sum(lambda_ln(:,i)), sum(lambda_cp(:,i))
   end do
   
   close(ifile) 

   ! write spectra summed over outgoing frequencies
   file="_xas"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')
   
   do i=1, n_omega_in, p % kh_print_stride
     write(ifile,'(5ES18.10)') omega_in(i), sum(lambda_lp(i,:)), sum(lambda_ln(i,:)), sum(lambda_cp(i,:))
   end do
   
   close(ifile) 

   
 end if! if(.false.) then
   
    
end subroutine calculate_SCKH_res_PES_FC_old


 
 subroutine compute_F_FC_if_om_omp(E_ni, E_n1, E_fc1, E_ni_mean, E_nf_mean, E_fi_mean, &
      D_fn1, D_ni1, omega_in, time, gamma, gamma_F, F_if_om_omp)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: fft_c2c_1d_backward
   use m_fftw3, only: fft_c2c_1d_forward
   use m_fftw3, only: reorder_sigma_fftw_z
   use m_SCKH_utils, only:compute_efactor_one
   use m_SCKH_utils, only:compute_efactor
   use m_SCKH_utils, only: put_on_grid
   
   real(wp), intent(in):: E_ni
   real(wp), intent(in):: E_n1(:)
   real(wp), intent(in):: E_fc1(:)
   real(wp), intent(in):: E_ni_mean
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in):: D_fn1(:,:)
   real(wp), intent(in):: D_ni1(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(in):: time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: gamma_F
   complex(wp), intent(out):: F_if_om_omp(:,:,:,:)
   
   integer:: ntsteps, ninter
   ! complex(wp), allocatable:: funct(:,:,:)
   complex(wp), allocatable::  e_factor1(:)
   complex(wp), allocatable::  F_tmp(:)
   integer:: m1, m2, om_in, n_omega_in, n_omega_out

   integer:: i_low
   real(wp):: i_low_weight
   
   n_omega_in = size(F_if_om_omp,1)
   n_omega_out = size(F_if_om_omp,2)
   
   allocate(e_factor1(n_omega_in))
   allocate(F_tmp(n_omega_out))

   !write(6,*) "Here... 1...", n_omega_in, n_omega_out
   !write(6,*) "E_n1(:)", E_n1(:)
   !write(6,*) "E_fc1(f_e,1,:)",E_fc1(:)
   !write(6,*) "E_fi_mean", E_nf_mean
   !write(6,*) "time", time
   !call compute_efactor_one(E_fc1(:)- E_n1(:), E_fi_mean, time, e_factor1(:), .true.)

!   write(6,*) "E_ni", E_ni
   
   call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .true.)
   !call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .false.)
   
   F_if_om_omp = 0.0_wp

   call put_on_grid(omega_in, E_ni, i_low, i_low_weight )

   !write(6,*) "Here... 2"

   !write(6,*) e_factor1
   
   do om_in= 1, n_omega_in
     do m1 =1, 3
       ! the factor of exp(-gamma * const % eV * time(:) / const % hbar) seems to belong here, but I must figure out why.
       call fft_c2c_1d_backward( e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar) *&
            !exp(dcmplx(0.0_wp,  (-omega_in(om_in) + E_ni )* const % eV * time(:) / const % hbar )) * &
            !( i_low_weight * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low))* const % eV * time(:) / const % hbar ))  &
            !+ (1.0_wp -i_low_weight) * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low+1))* const % eV *&
            !time(:) / const % hbar )) ) * &
            D_fn1(:,m1), &
            F_tmp(:))
       !call fft_c2c_1d_forward( e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar) *&
       !     exp(dcmplx(0.0_wp,  (omega_in(om_in) - E_ni )* const % eV * time(:) / const % hbar )) * &
        !    D_fn1(:,m1), &
        !    F_tmp(:))
       call reorder_sigma_fftw_z(F_tmp(:))
       
       do m2 =1, 3
         !F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) / dcmplx(omega_in(om_in) - E_ni, gamma)

         F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) * &
              i_low_weight / dcmplx(omega_in(om_in) - omega_in(i_low), gamma)
         F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) * &
              (1.0_wp -i_low_weight) / dcmplx(omega_in(om_in) - omega_in(i_low+1), gamma) 
       end do
       
     end do
   end do
   
 end subroutine compute_F_FC_if_om_omp

! various failed attempts to approximate things
 subroutine compute_F_FC_if_om_omp_alt(E_ni, E_nf, E_n1, E_fc1, E_ni_mean, E_nf_mean, E_fi_mean, &
      D_fn1, D_ni1, omega_in, omega_out, time, gamma, gamma_F, F_if_om_omp)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: fft_c2c_1d_backward
   use m_fftw3, only: fft_c2c_1d_forward
   use m_fftw3, only: reorder_sigma_fftw_z
   use m_SCKH_utils, only:compute_efactor_one
   use m_SCKH_utils, only:compute_efactor
   use m_SCKH_utils, only: put_on_grid
   
   real(wp), intent(in):: E_ni
   real(wp), intent(in):: E_nf
   real(wp), intent(in):: E_n1(:)
   real(wp), intent(in):: E_fc1(:)
   real(wp), intent(in):: E_ni_mean
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in):: D_fn1(:,:)
   real(wp), intent(in):: D_ni1(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(in):: time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: gamma_F
   complex(wp), intent(out):: F_if_om_omp(:,:,:,:)
   
   integer:: ntsteps, ninter
   ! complex(wp), allocatable:: funct(:,:,:)
   complex(wp), allocatable::  e_factor1(:)
   complex(wp), allocatable::  F_tmp(:)
   integer:: m1, m2, om_in, om_out, n_omega_in, n_omega_out

   integer:: i_low
   real(wp):: i_low_weight
   
   n_omega_in = size(F_if_om_omp,1)
   n_omega_out = size(F_if_om_omp,2)
   
   allocate(e_factor1(n_omega_in))
   allocate(F_tmp(n_omega_out))

   !write(6,*) "Here... 1...", n_omega_in, n_omega_out
   !write(6,*) "E_n1(:)", E_n1(:)
   !write(6,*) "E_fc1(f_e,1,:)",E_fc1(:)
   !write(6,*) "E_fi_mean", E_nf_mean
   !write(6,*) "time", time
   !call compute_efactor_one(E_fc1(:)- E_n1(:), E_fi_mean, time, e_factor1(:), .true.)

!   write(6,*) "E_ni", E_ni
   
   call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .true.)
   !call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .false.)
   
   F_if_om_omp = 0.0_wp

   !call put_on_grid(omega_in, E_ni, i_low, i_low_weight )
   !call put_on_grid(omega_out, E_nf, i_low2, i_low_weight2 )

   !write(6,*) "Here... 2"

   !write(6,*) e_factor1
   
   do m1 =1, 3
     ! the factor of exp(-gamma * const % eV * time(:) / const % hbar) seems to belong here, but I must figure out why.
     call fft_c2c_1d_backward( e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar) *&
          !exp(dcmplx(0.0_wp,  (-omega_in(om_in) + E_ni )* const % eV * time(:) / const % hbar )) * &
          !( i_low_weight * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low))* const % eV * time(:) / const % hbar ))  &
          !+ (1.0_wp -i_low_weight) * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low+1))* const % eV *&
          !time(:) / const % hbar )) ) * &
          D_fn1(:,m1), &
          F_tmp(:))
     !call fft_c2c_1d_forward( e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar) *&
     !     exp(dcmplx(0.0_wp,  (omega_in(om_in) - E_ni )* const % eV * time(:) / const % hbar )) * &
     !    D_fn1(:,m1), &
     !    F_tmp(:))
     call reorder_sigma_fftw_z(F_tmp(:))
     
     do om_in= 1, n_omega_in
       do om_out= 1, n_omega_out

         do m2 =1, 3
         !F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) / dcmplx(omega_in(om_in) - E_ni, gamma)

         !F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) * &
         !     i_low_weight / dcmplx(omega_in(om_in) - omega_in(i_low), gamma)
         !F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) * &
           !     (1.0_wp -i_low_weight) / dcmplx(omega_in(om_in) - omega_in(i_low+1), gamma)

           ! Dyson type separation  
           !F_if_om_omp(om_in,om_out,m1,m2) = F_if_om_omp(om_in,om_out,m1,m2) +  F_tmp(om_out) * D_ni1(m2)  &
           !     / dcmplx(omega_in(om_in) - E_ni, gamma)

           !F_if_om_omp(om_in,om_out,m1,m2) = F_if_om_omp(om_in,om_out,m1,m2) +  F_tmp(om_out) * D_ni1(m2)  &
           !     / (dcmplx(omega_in(om_in) - E_ni, gamma) * dcmplx(omega_in(om_in) - E_ni -(omega_out(om_out) - E_nf), gamma) )

           
           !F_if_om_omp(om_in,om_out,m1,m2) = F_if_om_omp(om_in,om_out,m1,m2) +  F_tmp(om_out) * D_ni1(m2) * &
           !     (omega_out(om_out) - E_nf)**2 / (dcmplx(omega_in(om_in) - E_ni,gamma)* &
           !     dcmplx(omega_in(om_in) - E_ni -(omega_out(om_out) - E_nf), gamma))

           !F_if_om_omp(om_in,om_out,m1,m2) = F_if_om_omp(om_in,om_out,m1,m2) +  F_tmp(om_out) * D_ni1(m2) * &
           !     (omega_out(om_out) - E_nf)/dcmplx(omega_in(om_in) - E_ni,gamma)**2 

           ! with coupling
           F_if_om_omp(om_in,om_out,m1,m2) = F_if_om_omp(om_in,om_out,m1,m2) +  F_tmp(om_out) * D_ni1(m2) / &
                !dcmplx(omega_in(om_in) - E_ni -(omega_out(om_out) - E_nf), gamma)
                dcmplx((omega_out(om_out) - E_nf) - (omega_in(om_in) - E_ni), gamma)


           
         end do
       
       end do
     end do

   end do! do m1 =1, 3

 
 end subroutine compute_F_FC_if_om_omp_alt

 ! Lars' suggestion
 subroutine compute_F_FC_if_om_omp_alt2(E_ni, E_nf, E_n1, E_fc1, E_ni_mean, E_nf_mean, E_fi_mean, &
      D_fn1, D_ni1, omega_in, omega_out, time, gamma, gamma_F, F_if_om_omp)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: fft_c2c_1d_backward
   use m_fftw3, only: fft_c2c_1d_forward
   use m_fftw3, only: reorder_sigma_fftw_z
   use m_SCKH_utils, only:compute_efactor_one
   use m_SCKH_utils, only:compute_efactor
   use m_SCKH_utils, only: put_on_grid
   use m_SCKH_utils, only: compute_int_t_inf_z
   
   real(wp), intent(in):: E_ni
   real(wp), intent(in):: E_nf
   real(wp), intent(in):: E_n1(:)
   real(wp), intent(in):: E_fc1(:)
   real(wp), intent(in):: E_ni_mean
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in):: D_fn1(:,:)
   real(wp), intent(in):: D_ni1(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(in):: time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: gamma_F
   complex(wp), intent(out):: F_if_om_omp(:,:,:,:)
   
   integer:: ntsteps, ninter
   ! complex(wp), allocatable:: funct(:,:,:)
   complex(wp), allocatable::  e_factor1(:)
   complex(wp), allocatable::  F_tmp(:)
   complex(wp), allocatable::  integral(:)
   integer:: m1, m2, om_in, om_out, n_omega_in, n_omega_out

   integer:: i_low
   real(wp):: i_low_weight
   
   n_omega_in = size(F_if_om_omp,1)
   n_omega_out = size(F_if_om_omp,2)
   
   allocate(e_factor1(n_omega_out))
   allocate(F_tmp(n_omega_in))
   allocate(integral(n_omega_out))

   !write(6,*) "Here... 1...", n_omega_in, n_omega_out
   !write(6,*) "E_n1(:)", E_n1(:)
   !write(6,*) "E_fc1(f_e,1,:)",E_fc1(:)
   !write(6,*) "E_fi_mean", E_nf_mean
   !write(6,*) "time", time
   !call compute_efactor_one(E_fc1(:)- E_n1(:), E_fi_mean, time, e_factor1(:), .true.)

!   write(6,*) "E_ni", E_ni
   
   call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .true.)
   !call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .false.)
   
   F_if_om_omp = 0.0_wp

   !call put_on_grid(omega_in, E_ni, i_low, i_low_weight )
   !call put_on_grid(omega_out, E_nf, i_low2, i_low_weight2 )

   !write(6,*) "Here... 2"

   !write(6,*) e_factor1
   
   do m2 =1, 3

     do om_out= 1, n_omega_out
       
       ! \int_{t}^{\infty} ds D(s) e^{-i \int_0^{s} d\tau E_{nf}(\tau)} e^{i(\omega' +i\Gamma)s}
       call compute_int_t_inf_z(e_factor1(:) * D_fn1(:,m2) * exp(dcmplx(-gamma, omega_out(om_out) - E_nf_mean) *&
            const % eV * time(:) / const % hbar), &
            time(:), integral(:))

         call fft_c2c_1d_forward(  exp(dcmplx(+gamma-gamma_F, (E_ni-E_ni_mean)) * const % eV * time(:) / const % hbar ) *&
              integral(:), &
              F_tmp(:))
         call reorder_sigma_fftw_z(F_tmp(:))

         do m1=1,3
           F_if_om_omp(:, om_out, m1,m2) = D_ni1(m1) * F_tmp(:)
         end do 

     end do ! do om_out= 1, n_omega_out
   end do !do m2 =1, 3
   
 end subroutine compute_F_FC_if_om_omp_alt2
 

 ! include a t-dependence  e^{-i (\int^t_ d\tau E_nf{\tau} -\omega' t) in the XES spectrum to get the connection  
 subroutine compute_F_FC_if_om_omp_alt3(E_ni, E_n1, E_fc1, E_ni_mean, E_nf_mean, E_fi_mean, &
      D_fn1, D_ni1, omega_in, omega_out, time, gamma, gamma_F, F_if_om_omp)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: fft_c2c_1d_backward
   use m_fftw3, only: fft_c2c_1d_forward
   use m_fftw3, only: reorder_sigma_fftw_z
   use m_SCKH_utils, only:compute_efactor_one
   use m_SCKH_utils, only:compute_efactor
   use m_SCKH_utils, only: put_on_grid
   
   real(wp), intent(in):: E_ni
   real(wp), intent(in):: E_n1(:)
   real(wp), intent(in):: E_fc1(:)
   real(wp), intent(in):: E_ni_mean
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in):: D_fn1(:,:)
   real(wp), intent(in):: D_ni1(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(in):: omega_out(:)
   real(wp), intent(in):: time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: gamma_F
   complex(wp), intent(out):: F_if_om_omp(:,:,:,:)
   
   integer:: ntsteps, ninter
   ! complex(wp), allocatable:: funct(:,:,:)
   complex(wp), allocatable::  e_factor1(:)
   complex(wp), allocatable::  e_factor2(:)
   complex(wp), allocatable::  F_tmp(:)
   complex(wp), allocatable::  F_tmp_om_out(:,:)
   integer:: m1, m2, om_in, om_out, n_omega_in, n_omega_out

   real(wp), allocatable:: E_n1_tmp(:)
   
   integer:: i_low
   real(wp):: i_low_weight
   
   n_omega_in = size(F_if_om_omp,1)
   n_omega_out = size(F_if_om_omp,2)
   
   allocate(e_factor1(n_omega_in))
   allocate(e_factor2(n_omega_in))
   allocate(F_tmp_om_out(n_omega_out, 3))
   allocate(F_tmp(n_omega_in))
   allocate(E_n1_tmp(n_omega_in))
   
   !write(6,*) "Here... 1...", n_omega_in, n_omega_out
   !write(6,*) "E_n1(:)", E_n1(:)
   !write(6,*) "E_fc1(f_e,1,:)",E_fc1(:)
   !write(6,*) "E_fi_mean", E_nf_mean
   !write(6,*) "time", time
   !call compute_efactor_one(E_fc1(:)- E_n1(:), E_fi_mean, time, e_factor1(:), .true.)

   !   write(6,*) "E_ni", E_ni
   
   call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .true.)

   E_n1_tmp = E_ni
   call compute_efactor( E_n1_tmp(:), E_fc1(:), E_nf_mean, time, e_factor2(:), .true.)
   !call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .false.)
   
   F_if_om_omp = 0.0_wp

   call put_on_grid(omega_in, E_ni, i_low, i_low_weight )

   !write(6,*) "Here... 2"

   !write(6,*) e_factor1
   

   ! XES spectrum
   do m1 =1, 3
     ! the factor of exp(-gamma * const % eV * time(:) / const % hbar) seems to belong here, but I must figure out why.
     call fft_c2c_1d_backward( e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar) *&
          !exp(dcmplx(0.0_wp,  (-omega_in(om_in) + E_ni )* const % eV * time(:) / const % hbar )) * &
          !( i_low_weight * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low))* const % eV * time(:) / const % hbar ))  &
          !+ (1.0_wp -i_low_weight) * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low+1))* const % eV *&
          !time(:) / const % hbar )) ) * &
          D_fn1(:,m1), &
          F_tmp_om_out(:,m1))
     call reorder_sigma_fftw_z(F_tmp_om_out(:, m1))
     
   end do

   ! XAS and coupling to XES
   do om_out= 1, n_omega_out
     do m2 =1, 3
       call fft_c2c_1d_forward(  exp(dcmplx(-gamma_F, (E_ni-E_ni_mean)) * const % eV * time(:) / const % hbar ) * &
            e_factor2(:) * exp(dcmplx(0.0_wp,  (omega_out(om_out) - E_nf_mean +E_ni)* const % eV * time(:) / const % hbar )) * &
            D_ni1(m2), &
            F_tmp(:))
       call reorder_sigma_fftw_z(F_tmp(:))
       
       do m1=1,3
         F_if_om_omp(:,om_out,m1,m2) = F_if_om_omp(:,om_out,m1,m2) +  F_tmp(:) *  F_tmp_om_out(om_out, m1)
       end do

     end do
   end do
   
 end subroutine compute_F_FC_if_om_omp_alt3

 ! this is the method I tried implementing first, that actually seemed to work more or less, but without an additional broadening
 ! it worked less well. It inlcudes the omega -E_in in the emission
 subroutine compute_F_FC_if_om_omp_alt4(E_ni, E_n1, E_fc1, E_ni_mean, E_nf_mean, E_fi_mean, &
      D_fn1, D_ni1, omega_in, time, gamma, gamma_F, F_if_om_omp)

   use m_precision, only: wp
   use m_constants, only: const
   use m_fftw3, only: fft_c2c_1d_backward
   use m_fftw3, only: fft_c2c_1d_forward
   use m_fftw3, only: reorder_sigma_fftw_z
   use m_SCKH_utils, only:compute_efactor_one
   use m_SCKH_utils, only:compute_efactor
   use m_SCKH_utils, only: put_on_grid
   
   real(wp), intent(in):: E_ni
   real(wp), intent(in):: E_n1(:)
   real(wp), intent(in):: E_fc1(:)
   real(wp), intent(in):: E_ni_mean
   real(wp), intent(in):: E_nf_mean
   real(wp), intent(in):: E_fi_mean
   real(wp), intent(in):: D_fn1(:,:)
   real(wp), intent(in):: D_ni1(:)
   real(wp), intent(in):: omega_in(:)
   real(wp), intent(in):: time(:)
   real(wp), intent(in):: gamma
   real(wp), intent(in):: gamma_F
   complex(wp), intent(out):: F_if_om_omp(:,:,:,:)
   
   integer:: ntsteps, ninter
   ! complex(wp), allocatable:: funct(:,:,:)
   complex(wp), allocatable::  e_factor1(:)
   complex(wp), allocatable::  F_tmp(:)
   integer:: m1, m2, om_in, n_omega_in, n_omega_out

   integer:: i_low
   real(wp):: i_low_weight
   
   n_omega_in = size(F_if_om_omp,1)
   n_omega_out = size(F_if_om_omp,2)
   
   allocate(e_factor1(n_omega_in))
   allocate(F_tmp(n_omega_out))

   !write(6,*) "Here... 1...", n_omega_in, n_omega_out
   !write(6,*) "E_n1(:)", E_n1(:)
   !write(6,*) "E_fc1(f_e,1,:)",E_fc1(:)
   !write(6,*) "E_fi_mean", E_nf_mean
   !write(6,*) "time", time
   !call compute_efactor_one(E_fc1(:)- E_n1(:), E_fi_mean, time, e_factor1(:), .true.)

!   write(6,*) "E_ni", E_ni
   
   call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .true.)
   !call compute_efactor( E_n1(:), E_fc1(:), E_nf_mean, time, e_factor1(:), .false.)
   
   F_if_om_omp = 0.0_wp

   call put_on_grid(omega_in, E_ni, i_low, i_low_weight )

   !write(6,*) "Here... 2"

   !write(6,*) e_factor1
   
   do om_in= 1, n_omega_in
     do m1 =1, 3
       ! the factor of exp(-gamma * const % eV * time(:) / const % hbar) seems to belong here, but I must figure out why.
       call fft_c2c_1d_backward( e_factor1(:) * exp(-gamma_F * const % eV * time(:) / const % hbar) *&
            exp(dcmplx(0.0_wp,  (-omega_in(om_in) + E_ni )* const % eV * time(:) / const % hbar )) * &
            !( i_low_weight * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low))* const % eV * time(:) / const % hbar ))  &
            !+ (1.0_wp -i_low_weight) * exp(dcmplx(0.0_wp,  (-omega_in(om_in) + omega_in(i_low+1))* const % eV *&
            !time(:) / const % hbar )) ) * &
            D_fn1(:,m1), &
            F_tmp(:))
       !call fft_c2c_1d_forward( e_factor1(:) * exp(-gamma * const % eV * time(:) / const % hbar) *&
       !     exp(dcmplx(0.0_wp,  (omega_in(om_in) - E_ni )* const % eV * time(:) / const % hbar )) * &
        !    D_fn1(:,m1), &
        !    F_tmp(:))
       call reorder_sigma_fftw_z(F_tmp(:))
       
       do m2 =1, 3
         !F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) / dcmplx(omega_in(om_in) - E_ni, gamma)

         F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) * &
              i_low_weight / dcmplx(omega_in(om_in) - omega_in(i_low), gamma)
         F_if_om_omp(om_in,:,m1,m2) = F_if_om_omp(om_in,:,m1,m2) +  F_tmp(:) * D_ni1(m2) * &
              (1.0_wp -i_low_weight) / dcmplx(omega_in(om_in) - omega_in(i_low+1), gamma) 
       end do
       
     end do
   end do
   
 end subroutine compute_F_FC_if_om_omp_alt4

 

! ! given a equally spaced range x, and a point x_in
! ! find i_low, the index of the point below x_in
! ! and i_low_weigth in that point, assuming "tents"
! ! the weight on i_low +1 is 1-i_low_weigth
! subroutine put_on_grid(x, x_in, i_low, i_low_weight )
!
!   use m_precision, only: wp
!
!   real(wp), intent(in):: x(:)
!   real(wp), intent(in):: x_in
!   integer, intent(out):: i_low
!   real(wp), intent(out):: i_low_weight
!
!   integer:: nx
!   real(wp):: dx
!   
!   nx = size(x)
!   dx = x(2)-x(1)
!
!   i_low = floor((x_in-x(1))/dx) +1
!      
!   if(i_low .lt. 1) then
!     i_low = 1
!     i_low_weight =1.0_wp
!   else if(i_low .gt. nx-1) then
!     i_low = nx-1
!     i_low_weight =0.0_wp
!   else
!     i_low_weight = -x_in/dx +(1 + x(i_low)/dx)
!   end if
!
!!   write(6,*) "x_in, i_low, i_low_weight, x(i_low)", x_in, i_low, i_low_weight, x(i_low)
!   
!   
! end subroutine put_on_grid
 
end module m_SCKH_resonant_PES_FC



