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
    use m_SCKH_utils, only: compute_F_if_om_omp
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

    type(sckh_params_t), intent(inout):: p 

    integer:: ninter, nfinal, ntsteps  
    integer:: npoints_in

    real(wp), allocatable::  time(:), E_dyn_inp(:),  E_dyn2_inp(:)  
    real(wp), allocatable::  E_i_inp(:), E_i1(:), E_lp_corr(:), shift(:)
    real(wp), allocatable::  E_n_inp(:,:), E_n1(:,:), E_n2(:,:)
    real(wp), allocatable::  E_f_inp(:,:), E_f1(:,:), E_f2(:,:)

    real(wp), allocatable:: D_fn_inp(:,:,:,:), D_fn2(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:), D_ni1(:,:,:) 
    real(wp):: E_nf_mean, E_fi_mean
    
    character(80)::  string, file
    real(wp), allocatable:: x_sampl(:), mom_sampl(:), x_mom_sampl(:,:)
    real(wp), allocatable:: x_new(:), v_new(:), a_new(:)
    real(wp), allocatable:: x_new2(:), v_new2(:), a_new2(:)
    integer:: npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl

    real(wp), allocatable:: omega_in(:), omega_out(:)
    integer:: n_omega_in, n_omega_out

    real(wp), allocatable:: c_i(:,:)
    real(wp), allocatable:: eig_i(:)
    
    real(wp):: gamma, gamma_inc
    real(wp):: time_l, delta_t, norm, mu_SI, dx
    
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
         E_n1(ninter, ntsteps), &
         E_n2(ninter, ntsteps), &
         E_f1(nfinal,ntsteps), &
         E_f2(nfinal,ntsteps), &
         D_ni1(ninter, ntsteps,3), &
         D_fn2(nfinal, ninter, ntsteps,3), &
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
        call spline_easy(X_r, E_n_inp(i,:), npoints_in, x_new, E_n1(i,:), ntsteps)  
      end do

      do i=1,ninter
        do m=1,3
          call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new, D_ni1(i,:,m) , ntsteps)  
        end do
      end do

      do i=1,nfinal
        call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new, E_f1(i,:), ntsteps)  
      end do

      ! first time, compute the mean transition energy, and frequency
      if (traj .eq. 1) then
        ind = minloc(E_i_inp)
        E_fi_mean =  E_f1(nfinal, ind(1)) - E_i1(ind(1))
        write(6,*) "E_fi_mean", E_fi_mean
        
        E_nf_mean =  E_n1(ninter, ind(1)) - E_f1(nfinal,ind(1))
        write(6,*) "E_nf_mean", E_nf_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_in)
        omega_in = omega_in + E_nf_mean + E_fi_mean
        
        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
        omega_out = omega_out + E_nf_mean
        
      end if

      ! loop over each starting point of first trajectory
      do traj2 =1, ntsteps ! possibly with a stride
        
        call verlet_trajectory_xva(x_new(traj2), v_new(traj2), X_r, &
             E_dyn2_inp * const % eV, delta_t, mu_SI, x_new2, v_new2, a_new2 )
        
        call spline_easy(X_r, E_n_inp, npoints_in, x_new2, E_n2, ntsteps)
        
        do i=1,nfinal
          call spline_easy(X_r, E_f_inp(i,:), npoints_in, x_new2, E_f2(i,:), ntsteps)  
        end do
        
        do i=1,nfinal
          do m=1,3
            call spline_easy(X_r, D_fn_inp(i,1,:,m) , npoints_in, x_new2, D_fn2(i,1,:,m) , ntsteps)  
          end do
        end do
        
!        ! first time, compute the mean transition energy, and frequency
!        if (traj .eq. 1 .and. traj2 .eq. 1) then
!          ind = minloc(E_i_inp)
!          E_nf_mean =  E_n2(ninter, ind(1)) - E_f2(nfinal,ind(1))
!          write(6,*) "E_nf_mean", E_nf_mean
!          
!          call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega_out)
!          omega_out = omega_out + E_nf_mean
!
!          write(6,*) "here 1"
!        end if
        
        call compute_F_if_omp_many_n(E_n2(:,:), E_f2(:,:), E_nf_mean, D_fn2(:,:,:,:), &
             D_ni1(:, traj2, :), time,  F_if_t_omp(:,traj2,:,:,:), gamma)

        !write(6,*) "here 2", traj2
        
      end do ! do traj2 =1, ntsteps2

      call compute_F_if_om_omp(F_if_t_omp, E_f1(:,:), E_fi_mean, time, &
           E_i1, gamma_inc, omega_out, E_nf_mean, F_if_om_omp)
      
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
  
end module m_SCKH_resonant_PES
