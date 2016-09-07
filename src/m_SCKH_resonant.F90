module m_SCKH_resonant
  implicit none

contains

  subroutine calculate_SCKH_res_PES(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    use m_SCKH_utils, only: compute_SCKH
    use m_SCKH_utils, only: compute_SCKH_FFTW
    use m_SCKH_utils, only: verlet_trajectory
    use m_sckh_params_t, only: sckh_params_t 
    use hist_class, only: hist, hist_init, hist_add
    use hist_class, only: hist_broadening, hist_write
    use m_io, only: get_free_handle
    use m_splines, only: spline_easy
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_nac_file
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: get_projections
    use m_PES_io, only: read_file_list
    use m_KH_functions, only: solve_vib_problem
    use m_fftw3, only: get_omega_reordered_fftw

    type(sckh_params_t), intent(inout):: p 

    integer:: nfinal, ntsteps  
    integer:: npoints_in  
    real(kind=wp),dimension(:),allocatable::  E_n_inp, time, E_n, E_dyn_inp 
    real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f
    real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn 
    
    character(80)::  string, file
    real(kind=wp),dimension(:),allocatable:: sigma_tot
    real(kind=wp),allocatable:: sigma_mm(:,:,:,:)
    real(kind=wp),allocatable:: sigma_f(:,:)
    real(kind=wp),dimension(:),allocatable:: x_sampl, mom_sampl, x_new, omega
    real(kind=wp),dimension(:,:),allocatable:: sigma_proj, c_i, x_mom_sampl
    real(kind=wp):: gamma, time_l, delta_t, norm, E_nf_mean, mu_SI, dx
    integer:: n_omega, npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl
    complex(kind=wp), dimension(:,:,:),allocatable::  F_if_omp_m
    real(kind=wp),dimension(:),allocatable:: eig_i, E_i_inp, E_lp_corr, shift 
    type(hist), dimension(:), allocatable:: time_h
    type(hist):: time_h_0, time_h_0_mom
    integer, dimension(1)::ind  
    integer:: ifile
    
    real(kind=wp), dimension(:),allocatable::  X_r
    real(kind=wp) :: dvr_start 
    integer:: npoints, ii
    real(8):: dnrm2
    integer::i,j,m,m1,m2, traj

    ! set some local variables
    ntsteps = p % ntsteps
    npoints_x_sampl = p % npoints_x_sampl
    npoints_mom_sampl =  p % npoints_mom_sampl
    npoints_in = p % npoints_in ! this should be renamed, corresponds to p % nstates
    nfinal = p % npesfile_f
    mu_SI = p % mu * const % u

    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    delta_t = p % delta_t

    ! use HWHM internally
    gamma = p % gamma_FWHM / 2 

    ! projections
    call  get_projections(p)

    allocate( E_i_inp(npoints_in),&
         E_n_inp(npoints_in), &
         E_dyn_inp(npoints_in),&
         E_f_inp(nfinal,npoints_in), &
         D_fn_inp(nfinal,npoints_in,3), &
         time(ntsteps),&
         E_n(ntsteps), &
         E_f(nfinal,ntsteps), &
         D_fn(nfinal,ntsteps,3), &
         E_lp_corr(npoints_in),&
         shift(npoints_in), &
         c_i(npoints_in,npoints_in), &
         eig_i(npoints_in), &
         x_new(ntsteps),&
         time_h(ntsteps))
    allocate(X_r(npoints_in))

    ! set up grid points
    do i = 1, npoints_in
      X_r(i) = (i-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_r, E_i_inp)
    call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n_inp)

    ! read PES file where the dynamcis is run 
    if (p % use_dynamics_file) then
      call read_PES_file(p % pes_file_dyn, p % npoints_in, p % npoints_in, X_r, E_dyn_inp)
    else
      E_dyn_inp = E_n_inp
    end if

    ! read list of final state pes_files and dipole_files
    call read_file_list(p % pes_file_list_f, p % npesfile_f, p % pes_files_f)
    call read_file_list(p % dipole_file_list_f, p % npesfile_f, p % dipolefile_f)

    do j=1,p % npesfile_f
      call read_PES_file(p % pes_files_f(j), p % npoints_in, p % npoints_in, X_r, E_f_inp(j,:))
      call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % npoints_in, X_r, D_fn_inp(j,:,:))
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

    do j=1,nfinal
      E_f_inp(j,:) = E_f_inp(j,:) / const % eV
    end do

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
    !time_l2 = (ntsteps-1) * delta_t 

    write(6,*) "outfile", p % outfile
    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
    write(6,*) "mu_SI", mu_SI 
    write(6,*) "time_l", time_l
    write(6,*)
    write(6,*) "Fundamental frequency resolution", 2.0_wp * const % pi * const % hbar /( time_l * const % eV)
    write(6,*) "delta t", delta_t
    write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV), "eV"

    n_omega = ntsteps
    allocate(F_if_omp_m(nfinal,n_omega,3), sigma_f(nfinal,n_omega), sigma_tot(n_omega), &
         sigma_proj(p % nproj,n_omega), sigma_mm(nfinal,n_omega,3,3), omega(n_omega))

    !
    ! Loop over trajectories
    !

    do i=1, ntsteps
      time(i)= (i-1)*delta_t
    end do

    call hist_init(time_h_0, 1000, X_r(1), X_r(npoints_in) ) 
    call hist_init(time_h_0_mom, 1000, minval(x_mom_sampl(:,2)), maxval(x_mom_sampl(:,2)) ) 
    do i=1, ntsteps
      call hist_init(time_h(i), 1000, X_r(1), X_r(npoints_in) ) 
    end do

    sigma_mm=0.0_wp

    !
    !
    !
    
    do traj=1, npoints_x_mom_sampl

      call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
           E_dyn_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new )
      
      call spline_easy(X_r, E_n_inp, npoints_in, x_new, E_n, ntsteps)
      do i=1,nfinal
        call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_f(i,:), ntsteps)  
      end do

      ! here, compute efactor for E_fi
      do f_e = 1,nfinal 
        call compute_efactor(E_f(f_e,:), E_i, E_fi_mean, time, e_factor1(f_e,:), .true.)
      end do
      
      do traj2 =1, ntsteps2 
        traj_pos = (traj2-1) * traj_stride +1
        
        call verlet_trajectory_xva(x_new(traj_pos), v_new(traj_pos), X_r, &
             E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
        
        call spline_easy(X_r, E_n_inp, npoints_in, x_new2, E_n, ntsteps)

        do i=1,nfinal
          call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new2, E_f(i,:), ntsteps)  
        end do
        
        do i=1,nfinal
          do m=1,3
            call spline_easy(X_r, D_fn_inp(i,:,m) , npoints_in, x_new2, D_fn(i,:,m) , ntsteps)  
          end do
        end do
        
        ! first time, compute the mean transition energy, and frequency
        if (traj .eq. 1) then
          ind = minloc(E_i_inp)
          E_nf_mean =  E_n(ninter, ind(1)) - E_f(nfinal,ind(1))
          write(6,*) "E_nf_mean", E_nf_mean
          
          call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega)
          omega = omega + E_nf_mean
        end if
        
        call compute_F_if_omp_many_n(E_n(:,:), E_f(:,:), E_nf_mean, D_fn(:,:,:,:), D_ni(:, traj_pos, :), time,  F_if_omp_t(:,:,traj2,:,:), gamma)
        
      end do ! do traj2 =1, ntsteps2


      call compute_F_if_omp_om(F_if_omp_t, e_factor1, F_if_omp_om)
      ! should contain the things below
      
!      !fourier transform F_if_omp_t
!      do f_e= 1, nfinal
!        do om_out= 1, n_omega_out
!          do m1 =1, 3
!            do m2 =1, 3
!               call fft_c2c_1d_backward( e_factor1(f_e,:) * F_if_omp_t(f_e, om_out, :,m1,m2), &
!                                         F_if_omp_om(f_e, om_out, :,m1,m2))
!              call reorder_sigma_fftw_z(F_if_omp_om(f_e,om_out, :, m1,m2))
!          end do
!        end do
!      end do

      ! here collect spectrum parts 
      
    end do ! end traj
    
    ! compute other sigmas
    sigma_f(:,:) = 0.0_wp
    do m1=1,3
      sigma_f(:,:) = sigma_f(:,:) + sigma_mm(:,:,m1,m1)
    end do
    
    sigma_tot = sum(sigma_f(:,:),1)
    
    sigma_proj =0.0_wp
      do i=1,p % nproj
        do m1=1,3
          do m2=1,3
            sigma_proj(i,:) = sigma_proj(i,:) + p % projvec(i,m1) * sum( sigma_mm(:,:,m1,m2),1) * p % projvec(i,m2)  
          end do
        end do
      end do
      
      write(6,*)
      write(6,*) "Averaged over", npoints_x_mom_sampl, "trajectories"
      
    
    !
    ! Write output
    !


    !write distribution
    file="_time_h_0"
    file = trim(adjustl(p % outfile)) //  trim(adjustl(file))  // ".dat"
    call hist_broadening(time_h_0, 0.01d-10)
    call hist_write(time_h_0, file)

    file="_time_h_0_mom"
    file = trim(adjustl(p % outfile)) //  trim(adjustl(file))  // ".dat"
    !call hist_broadening(time_h_0_mom, 0.01d-10)
    call hist_write(time_h_0_mom, file)

    do i=1,10 !ntsteps
      file="_time_h_"
      write(string,*) i
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      call hist_broadening(time_h(i), 0.01d-10)
      call hist_write(time_h(i), file)
    end do

    !normalize 
    norm = sum(sigma_tot) *  2.0_wp * const % pi  * const % hbar / (time_l * const % eV ) 
    sigma_tot = sigma_tot / norm
    sigma_f = sigma_f / norm
    sigma_proj = sigma_proj / norm

    ! write sigma to file
    file="_sigma"
    file = trim(adjustl(p % outfile)) // trim(adjustl(file)) // ".dat"
    ifile = get_free_handle()
    open(ifile,file=file,status='unknown')

    do i=1,n_omega
      write(ifile,*)  omega(i), sigma_tot(i)
    end do

    close(ifile)

    ! write spectrum from different final states
    do j=1, nfinal

      file="_sigma_final_"
      write(string,*) j
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1,n_omega
        write(ifile,*)  omega(i), sigma_f(j,i)
      end do

      close(ifile)
    end do ! j

    ! write spectrum from projections  
    do j=1, p % nproj

      file="_sigma_proj_"
      write(string,*) j
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1,n_omega
        write(ifile,*)  omega(i), sigma_proj(j,i)
      end do

      close(ifile)
    end do !j

  end subroutine calculate_SCKH_res_PES

  
  
 ! a simplified routine, doing the simplest possible thing, also using FFTW
 subroutine SC_Kramers_heisenberg_resonant(pes_num, ntsteps, ntsteps2, delta_t)

    complex(kind=wp), allocatable:: efactor_fi(:,:) 

    complex(kind=wp), allocatable::  F_if_omp_t(:,:,:,:,:) 
    complex(kind=wp), allocatable::  F_if_omp_om(:,:,:,:,:)
    complex(kind=wp), allocatable::  F_if_omp_om_tot(:,:,:,:,:)
    complex(kind=wp), allocatable::  F_if_omp_0_tot(:,:,:,:)
    real(kind=wp):: E_fi_mean, E_ni_mean, E_nf_mean
    
    real(kind=wp), allocatable:: sigma_trace(:,:), sigma_trace_0(:)
    real(kind=wp), allocatable:: omega(:), omega_p(:) ! , omega_tmp(:)

    integer:: omp, om, traj, traj2

    delta_t = delta_t * 1.d-15

    call init_traj(trajectory1, ntsteps, ninterm, nfinal, delta_t)
    call init_traj(trajectory2, ntsteps2, ninterm, nfinal, delta_t)

    allocate(efactor_fi(nfinal, ntsteps))
    
    allocate(F_if_omp_t(nfinal,3,3, ntsteps2, ntsteps)) 
    allocate(F_if_omp_om(nfinal,3,3, ntsteps2, ntsteps))
    allocate(F_if_omp_om_tot(nfinal, 3,3, ntsteps2, ntsteps))
    allocate(F_if_omp_0_tot(nfinal, 3,3, ntsteps)
    
    allocate(sigma_trace(ntsteps2, ntsteps))
    allocate(sigma_trace_0(ntsteps2))
    
    allocate(omega(ntsteps))
    allocate(omega_p(ntsteps2))
    
    ! initialize variables
    F_f_omp_t = 0.0_wp
    F_f_omp_om = 0.0_wp
    F_f_omp_om_tot = 0.0_wp
    F_f_omp_0_tot = 0.0_wp

    write(6,*) "SCKH resonant"
    
    ! eV units (should probably be converted already when reading the files)
    !E_i_inp = E_i_inp / eV
    !E_n_inp = E_n_inp / eV
    !E_f_inp = E_f_inp / eV
    !E_dyn_inp = E_dyn_inp / eV

    ! alt 1: do an 'even' sampling of N points
    if(samplemode .eq. 1) then  
      call sample_x_mom(x, c_i(:,1), x_sampl, mom_sampl, x_mom_sampl, 1)
    else if(samplemode .eq. 2) then  
      call sample_x_mom(x, c_i(:,1), x_sampl, mom_sampl, x_mom_sampl, 2)
    end if
    
    !
    !write(6,*) "Sampling done"
    !
    !! give a kick in momentum!
    !!x_mom_sampl = x_mom_sampl * 1.0e3_wp
    !
    !!check
    !do i=1,npoints_x_sampl
    !  write(17,*) x_sampl(i)
    !end do
    !
    !do i=1,npoints_x_mom_sampl
    !  write(22,*) x_mom_sampl(i,1), x_mom_sampl(i,2)
    !end do
    !
    !do i=1,npoints_pes
    !  write(20,'(3ES16.6)') x(i), c_i(i,1), c_i(i,1) ** 2
    !end do

    !call init_traj_stuff(E_i_inp, E_n_inp, E_f_inp, E_fi_mean, E_ni_mean, E_nf_mean,  omega, omega_p) 
    call init_traj_stuff(pes % E_i, pes % E_n, pes % E_f, E_fi_mean, E_ni_mean, E_nf_mean,  omega, omega_p) 


    ! test
    !do traj=1, npoints_x_mom_sampl
    !  
    !  if(pes_type .eq. 'polynomial') then
    !
    !    call comp_SCKH_traj_pes_poly()
    !
    !  else if (pes_type .eq. "numerical") then
    !  
    !    call comp_SCKH_traj_pes_num()
    !    
    !  end if
    !
    !  F_if_omp_om_tot =  F_if_omp_om_tot + F_if_omp_om
    !
    !end do


    !initialize fourier transforms, I think I must pass the in and out arrays to my functions
    call dfftw_plan_dft_1d ( F_if_omp_plan, ntsteps, fftw_in1, fftw_out1, &
       FFTW_BACKWARD, FFTW_ESTIMATE )
    call dfftw_plan_dft_1d ( F_if_omp_om_plan, ntsteps2, fftw_in2, fftw_out2, &
         FFTW_FORWARD, FFTW_ESTIMATE )


    !
    ! Loop over trajectories
    !

    do traj=1, npoints_x_mom_sampl
      
      write(6,*) "running traj", traj
      
      call verlet_trajectory_traj(trajectory1, x_mom_sampl(traj,1), x_mom_sampl(traj,2)/my_SI, pes % x, pes % E_i * eV, delta_t, my_SI)
      
      do i=1,nfinal
        call compute_efactor(E_f_t(i,:), E_i_t, E_fi_mean, time, efactor_fi(i,:), .false.)      
      end do
      
      do traj2=1, ntsteps
        
        write(6,*) "running traj2", traj2
        
        call verlet_trajectory_traj(trajectory2, trajectory1 % x(traj2), trajectory1 % v(traj2), pes % E_dyn * eV, delta_t, my_SI)
        
        call compute_F_if_omp_t(trajectory2, E_nf_mean, F_if_omp_plan, fftw_in1, fftw_out1, F_if_omp_t(:,:,:,:,traj2))

      end do !do traj2=1, ntsteps
      
    end do ! traj2

    call compute_F_if_omp_om(F_if_omp_t, E_nf_mean, F_if_omp_om_plan, fftw_in2, fftw_out2, F_if_omp_om)

    F_if_omp_om_tot = F_if_omp_om_tot + F_if_omp_om

    write(6,*) "Computed trajectory", traj, x_mom_sampl(traj,1), x_mom_sampl(traj,2)
    
  end do ! end traj






  
  sigma_trace = 0.0_wp
  sigma_trace_0 = 0.0_wp
  do fin =1,nfinal
    do m1=1,3
      sigma_trace = sigma_trace + dreal(F_f_omp_om_tot(fin, m1,m1,:,:) *conjg( F_f_omp_om_tot(fin, m1,m1,:,:)))
      sigma_trace_0 = sigma_trace_0 + dreal(F_f_omp_0_tot(fin, m1,m1,:) *conjg( F_f_omp_0_tot(fin, m1,m1,:)))
    end do
  end do
  
  !sigma_ion_sc = sum(sigma_ion_sc_f,1)
  !sigma_unpol_sc_tmp = sum(real(sigma_unpol_sc_f),1)
  
  !! convolute
  !write(6,*) "convoluting instruemntal broadening"
  !sigma_unpol_sc = sigma_unpol_sc_tmp 
  !do om_out= 1, n_omega_out_sc
  !   do om_in= 1, n_omega_in_sc
  !      do om_in2= 1, n_omega_in_sc
  !
  !         sigma_unpol_sc(om_out,om_in2) = sigma_unpol_sc(om_out,om_in2) + sigma_unpol_sc_tmp(om_out,om_in) * & !(omega_out_sc(om_out)/omega_in_sc(om_in2)) 
  !           gaussian(0.0d0, omega_in_sc(om_in) -omega_in_sc(om_in2),  instrument_FWHM) 
  !   
  !      end do
  !   end do
  !end do ! om_out
  
  !write(6,*) "convoluting detector broadening"
  !sigma_unpol_sc_tmp = sigma_unpol_sc 
    !sigma_unpol_sc = 0
    !do om_in= 1, n_omega_in_sc        
    !   do om_out= 1, n_omega_out_sc
    !      do om_out2= 1, n_omega_out_sc
    !
    !         sigma_unpol_sc(om_out2,om_in) = sigma_unpol_sc(om_out2,om_in) + sigma_unpol_sc_tmp(om_out,om_in) * & !(omega_out_sc(om_out)/omega_in_sc(om_in2)) 
    !           gaussian(0.0d0, omega_out_sc(om_out) -omega_out_sc(om_out2),  detector_FWHM) 
    !   
    !      end do
    !   end do
    !end do ! om_out


    write(6,*)
    write(6,*) "Averaged over", npoints_x_mom_sampl, "trajectories"

    !!    

! ! normalize
! norm=sum(sigma_unpol_sc(:,1)) *(omega_out_sc(2) -omega_out_sc(1)) 
! write(6,*) "norm",norm
! do j=1,n_omega_in_sc
!   !   sigma_unpol_sc(:,j) = sigma_unpol_sc(:,j)/norm
! end do
!
! norm=sum(sigma_ion_sc) *(omega_out_sc(2) -omega_out_sc(1)) 
! write(6,*) "norm",norm
! !sigma_ion_sc = sigma_ion_sc/norm
!
!
! !write distributions
! file="_time_h_0"
! file = trim(adjustl(outfile)) //  trim(adjustl(file))  // ".dat"
! call hist_broadening(time_h_0, 0.01d-10)
! call hist_write(time_h_0, file)
!
! file="_time_h_0_mom"
! file = trim(adjustl(outfile)) //  trim(adjustl(file))  // ".dat"
! !call hist_broadening(time_h_0_mom, 0.01d-10)
! call hist_write(time_h_0_mom, file)
!
! do i=1,10 !ntsteps
!   file="_time_h_"
!   write(string,*) i
!   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!   call hist_broadening(time_h(i), 0.01d-10)
!   call hist_write(time_h(i), file)
! end do
!
!
!
!    file=outfile
!    !write(string,'(F6.2)') omega_in(j)   
!    file = trim(adjustl(file)) // "_sigma_ion.dat"
!
!    open(10,file=file,status='unknown')
!
!    do i=1, n_omega_out_sc
!      write(10,*) omega_out_sc(i), sigma_ion_sc(i)
!    end do
!
!    close(10) 
!
!    ! write sigma to file
!    do j=1,n_omega_in_sc
!
!      file=outfile
!      write(string,'(F6.2)') omega_in_sc(j)   
!      write(6,'(F6.2)') omega_in_sc(j)   
!
!      file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!
!      write(31,*) file
!
!      open(10,file=file,status='unknown')
!
!      do i=1,n_omega_out_sc
!        write(10,*) omega_out_sc(i), sigma_unpol_sc(i,j)
!      end do
!
!      close(10) 
!
!    end do
!

!    ! write sigma to file
!    do j=1, n_omega
!      
!      file=outfile
!      write(string,'(F7.2)') omega(j) + E_ni_mean  
!      write(6,*) omega(j) + E_ni_mean, string  
!      
!      file = trim(adjustl(file)) // "_sc_" // trim(adjustl(string)) // ".dat"
!
!      write(31,*) file
!
!      open(10,file=file,status='unknown')
!      
!      do i=1,n_omega_p
!        !write(10,*) omega(j)- omega_p(i) + E_in_mean - E_if_mean , sigma_trace(j,i)
!        write(10,*) omega(j)+ E_ni_mean - omega_p(i) - E_fi_mean , sigma_trace(j,i)
!      end do
!      
!      close(10) 
!
!    end do

    ! write sigma to file
    do j=1, n_omega
      
      file=outfile
      write(string,'(F7.2)') omega(j) + E_ni_mean  !omega_p(j) + E_ni_mean  
      !write(6,*) omega(j) + E_ni_mean, string  
      
      file = trim(adjustl(file)) // "_sc_" // trim(adjustl(string)) // ".dat"

      write(31,*) file

      open(10,file=file,status='unknown')
 
      do i=1, n_omega_p
        !write(10,*) omega(j)- omega_p(i) + E_in_mean - E_if_mean , sigma_trace(j,i)
        !write(10,*) omega(j)+ E_ni_mean - omega_p(i) - E_fi_mean , sigma_trace(j,i)
        !write(10,*) omega(i) - E_nf_mean, sigma_trace(i, j -i)
        write(10,*) omega_p(i) + E_nf_mean, sigma_trace(i, j)
      end do
      
      close(10) 

    end do

    ! write xes spectrum at time t to file
    !do j=1, n_omega
    !  
    !  file=outfile
    !  write(string,'(F7.2)') omega(j) + E_ni_mean  !omega_p(j) + E_ni_mean  
    !  !write(6,*) omega(j) + E_ni_mean, string  
    !  
    !  file = trim(adjustl(file)) // "_omp_t_" // trim(adjustl(string)) // ".dat"
    !
    !  write(31,*) file
    !
    !  open(10,file=file,status='unknown')
    !
    !  do i=1, n_omega_p
    !    write(10,*) omega_p(i) + E_nf_mean, F_f_omp_t(fin, m1,m2,omp,:)
    !  end do
    !  
    !  close(10) 
    !
    !end do



    ! write the non-resonant XES cross section |\tilde F_f_{omp -Ef}(0)|^2
    file=outfile
    file = trim(adjustl(file)) // "_sc_XES_0.dat"
    open(10,file=file,status='unknown')

    do i=1,n_omega
      write(10,*) -omega(i) + E_ni_mean - E_fi_mean , sigma_trace_0(i)
      !write(10,*) i , sigma_trace_0(i)
    end do

    close(10) 


  end subroutine SC_Kramers_heisenberg_resonant



  
end module m_SCKH_resonant
