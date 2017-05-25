module m_SCXAS_PES

  implicit none

contains

  subroutine calculate_SCXAS_PES(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    !use m_SCKH_utils, only: compute_F_if_omp
    use m_SCKH_utils, only: compute_SCXAS_sigma_f_om_mm
    use m_SCKH_utils, only: compute_SCXAS_sigma_FC_f_om_mm
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
    use m_upper, only: upper
    
    type(sckh_params_t), intent(inout):: p 

    integer:: nfinal, ntsteps  
    integer:: npoints_in  
    real(kind=wp),dimension(:),allocatable::  time, E_dyn_inp 
    real(kind=wp),dimension(:,:),allocatable:: E_n_inp, E_n
    real(kind=wp),dimension(:,:,:),allocatable:: D_ni_inp, D_ni
    
    character(80)::  string, file
    real(kind=wp),dimension(:),allocatable:: sigma_tot
    real(kind=wp),allocatable:: sigma_n_om_mm(:,:,:,:)
    real(kind=wp),allocatable:: sigma_tmp(:,:,:,:)
    real(kind=wp),allocatable:: sigma_n(:,:)
    real(kind=wp),dimension(:),allocatable:: x_sampl, mom_sampl, x_new, omega
    real(kind=wp),dimension(:,:),allocatable:: sigma_proj, c_i, x_mom_sampl
    real(kind=wp):: gamma, time_l, delta_t, norm, E_ni_mean, mu_SI, dx

    integer:: n_omega, npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl
    !complex(kind=wp), allocatable::  F_if_omp(:,:,:,:)
    real(kind=wp),dimension(:),allocatable:: eig_i, E_i_inp, E_i, E_lp_corr, shift 
    type(hist), dimension(:), allocatable:: time_h
    type(hist):: time_h_0, time_h_0_mom
    integer, dimension(1)::ind  
    integer:: ifile
    
    real(kind=wp), dimension(:),allocatable::  X_r
    real(kind=wp) :: dvr_start 
    integer:: npoints, ii
    real(8):: dnrm2
    integer::i,j,m,m1,m2, traj
    !real(kind=wp),dimension(:),allocatable:: D_ni(:)

    real(wp):: E_i_mean
    real(wp):: E_ni
    real(wp):: E_i_min, E_i_av
    real(wp), allocatable:: E_n0(:)
    
    ! set some local variables
    ntsteps = p % ntsteps
    npoints_x_sampl = p % npoints_x_sampl
    npoints_mom_sampl =  p % npoints_mom_sampl
    npoints_in = p % npoints_in ! this should be renamed, corresponds to p % nstates
    !nfinal = p % npesfile_n
    mu_SI = p % mu * const % u

    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    delta_t = p % delta_t

    ! use HWHM internally
    gamma = p % gamma_FWHM / 2 

    ! projections
    call  get_projections(p)

    allocate( E_i_inp(npoints_in),&
         E_i(ntsteps), &
         !E_n_inp(npoints_in), &
         E_dyn_inp(npoints_in),&
         E_n_inp(p % npesfile_n,npoints_in), &
         D_ni_inp(p % npesfile_n,npoints_in,3), &
         time(ntsteps),&
         !E_n(ntsteps), &
         E_n(p % npesfile_n,ntsteps), &
         D_ni(p % npesfile_n,ntsteps,3), &
         E_lp_corr(npoints_in),&
         shift(npoints_in), &
         c_i(npoints_in,npoints_in), &
         eig_i(npoints_in), &
         x_new(ntsteps),&
         time_h(ntsteps))
    allocate(X_r(npoints_in))
    allocate(E_n0(npoints_in))
    
    !allocate(D_ni(3))

    !D_ni = 1.0_wp
    
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
    
    ! read PES file where the dynamcis is run 
    if (p % use_dynamics_file) then
      call read_PES_file(p % pes_file_dyn, p % npoints_in, p % npoints_in, X_r, E_dyn_inp)
    else
      E_dyn_inp =  E_i_inp !E_n_inp
    end if

    ! read list of final state pes_files and dipole_files
    call read_file_list(p % pes_file_list_n, p % npesfile_n, p % pes_files_n)
    call read_file_list(p % dipole_file_list_n, p % npesfile_n, p % dipolefile_n)

    do j=1,p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, p % npoints_in, X_r, E_n_inp(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % npoints_in, X_r, D_ni_inp(j,:,:))
    end do

    ! Shift orbital energies so that E_n(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved
    if( p % shift_PES .eq.  1) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % npoints_in, X_r, E_lp_corr)

      shift = E_lp_corr -E_n_inp(1,:) 

      do j=1,p % npesfile_n
        E_n_inp(j,:) = E_n_inp(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if

    ! Solve the vibrational problem for initial state to be able to sample initial distribution
    call solve_vib_problem(dx, E_i_inp, eig_i, c_i, mu_SI, p % vib_solver)
    write(6,*) "Calculated initial state eigenfunctions"
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm

    
    ! convert to eV units
    E_i_inp = E_i_inp  / const % eV
    E_n0 = E_n0  / const % eV
    !E_n_inp = E_n_inp  / const % eV
    E_dyn_inp = E_dyn_inp  / const % eV
    E_n_inp = E_n_inp / const % eV

    ind = minloc(E_i_inp)
    E_i_min = E_i_inp(ind(1))
    E_i_av = sum(E_i_inp * c_i(:,1)**2) - E_i_min

    write(6,*) "Average potenital energy, from minimum", E_i_av
    
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

    n_omega = ntsteps
    allocate(sigma_n_om_mm(p % npesfile_n,n_omega,3,3), &
         sigma_tmp(p % npesfile_n,n_omega,3,3), &
         sigma_n(p % npesfile_n,n_omega), &
         sigma_tot(n_omega), &
         sigma_proj(p % nproj,n_omega), &
         omega(n_omega))

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

    sigma_n_om_mm = 0.0_wp

    do traj=1, npoints_x_mom_sampl

      call hist_add(time_h_0, x_mom_sampl(traj,1), 1.0d0)   
      call hist_add(time_h_0_mom, x_mom_sampl(traj,2), 1.0d0)    

      ! FC code 
      if (upper(p % runmode_scxas) .eq. "FC" ) then
        !write(6,*) "rumode_scxas FC"
        
        if (traj .eq. 1) then
          ind = minloc(E_i_inp)
          E_ni_mean =  E_n_inp(max(1,p % npesfile_n/2),ind(1)) + E_n0(ind(1))- E_i_inp(ind(1))
          write(6,*) "E_ni_mean", E_ni_mean
          E_i_mean =  E_i_inp(ind(1))
          write(6,*) "E_i_mean", E_i_mean
          
          call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega)
          omega = omega + E_ni_mean
        end if
        
        call spline_easy(X_r, E_i_inp, npoints_in, x_mom_sampl(traj,1), E_i(1), 1)

        do i=1,p % npesfile_n
          call spline_easy(X_r, E_n_inp(i,:) + E_n0, npoints_in, x_mom_sampl(traj,1), E_n(i,1), 1)  
        end do
        

        ! options for dipole moments in XAS
        if (upper(p % dipole_mode) .eq. "DIPOLE") then
          do i=1,p % npesfile_n
            do m=1,3
              call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new, D_ni(i,1,m) , 1)  
            end do
          end do
        else if(upper(p % dipole_mode) .eq. "FC") then
          D_ni = 1.0_wp
        else if(upper(p % dipole_mode) .eq. "DIPOLE_X0") then
          D_ni(:,1,:) = D_ni_inp(:,ind(1),:) 
        else
          write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
          stop
        end if
        
        if(p % include_ZPE) then
          call compute_SCXAS_sigma_FC_f_om_mm(E_i_min + E_i_av, E_n(:,1), D_ni(:,1,:), omega, sigma_tmp, gamma)
        else
          call compute_SCXAS_sigma_FC_f_om_mm(E_i(1), E_n(:,1), D_ni(:,1,:), omega, sigma_tmp, gamma)
        end if

      ! SCXAS code
      else
        
      call verlet_trajectory(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, E_dyn_inp * const % eV, delta_t, mu_SI, x_new )

      do i=1, ntsteps
        call hist_add(time_h(i), x_new(i), 1.0d0)
      end do

      ! look up energies as a function of distance (which in turn is a function of time)
      !call spline_easy(X_r, E_n_inp, npoints_in, x_new, E_n, ntsteps)
      call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i, ntsteps)

      do i=1,p % npesfile_n
        call spline_easy(X_r, E_n_inp(i,:) + E_n0, npoints_in, x_new, E_n(i,:), ntsteps)  
      end do

      do i=1,p % npesfile_n
        do m=1,3
          call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new, D_ni(i,:,m) , ntsteps)  
        end do
      end do


      ! options for dipole moments in XAS
      if (upper(p % dipole_mode) .eq. "DIPOLE") then
        do i=1,p % npesfile_n
          do m=1,3
            call spline_easy(X_r, D_ni_inp(i,:,m) , npoints_in, x_new, D_ni(i,:,m) , 1)  
          end do
        end do
      else if(upper(p % dipole_mode) .eq. "FC") then
        D_ni = 1.0_wp
      else if(upper(p % dipole_mode) .eq. "DIPOLE_X0") then
        do i=1, p % nstates
          D_ni(:,i,:) = D_ni_inp(:,ind(1),:) 
        end do
      else
        write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
        stop
      end if
      
        
      ! first time, compute the mean transition energy, and frequency
      if (traj .eq. 1) then
        ind = minloc(E_i_inp)
        E_ni_mean =  E_n_inp(max(1,p % npesfile_n/2),ind(1)) + E_n0(ind(1)) - E_i(ind(1))
        write(6,*) "E_ni_mean", E_ni_mean
        E_i_mean =  E_i_inp(ind(1))
        write(6,*) "E_i_mean", E_i_mean
        !write(6,*) "E_i_mean", E_i_mean

        call get_omega_reordered_fftw(time_l * const % eV /  const % hbar, omega)
        omega = omega + E_ni_mean
      end if

      if(p % include_ZPE) then
        E_i = E_i_min + E_i_av !E_i_mean
        call compute_SCXAS_sigma_f_om_mm(E_i, E_n, E_ni_mean, D_ni, time, sigma_tmp, gamma)
      else
        call compute_SCXAS_sigma_f_om_mm(E_i, E_n, E_ni_mean, D_ni, time, sigma_tmp, gamma)
      end if
      
    end if !  if (upper(p % runmode_scxas) .eq. "FC" ) then
    
    sigma_n_om_mm = sigma_n_om_mm + sigma_tmp
      
  end do ! end traj
    
    
    ! compute other sigmas
    sigma_n(:,:) = 0.0_wp
    do m1=1,3
      sigma_n(:,:) = sigma_n(:,:) + sigma_n_om_mm(:,:,m1,m1)
    end do

    sigma_tot = sum(sigma_n(:,:),1)

    sigma_proj =0.0_wp
    do i=1,p % nproj
      do m1=1,3
        do m2=1,3
          sigma_proj(i,:) = sigma_proj(i,:) + p % projvec(i,m1) * &
               sum( sigma_n_om_mm(:,:,m1,m2),1) * p % projvec(i,m2)  
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
    sigma_n = sigma_n / norm
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
    do j=1, p % npesfile_n

      file="_sigma_final_"
      write(string,*) j
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1,n_omega
        write(ifile,*)  omega(i), sigma_n(j,i)
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

  end subroutine calculate_SCXAS_PES

end module m_SCXAS_PES
