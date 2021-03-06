module m_XAS_eigenstates
!  use parameters
!  use m_XAS_io
!  use KH_functions
!  use spline_m
!  use m_XAS_functions
  implicit none

contains
  
  subroutine calculate_XAS(p)
    use m_precision, only: wp
    use m_sckh_params_t, only: sckh_params_t 
    use m_constants, only: const
    use m_splines, only: spline_easy
    use m_splines, only: linspace
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_file_list
    use m_io, only: get_free_handle
    use m_KH_functions, only: solve_vib_problem
    use m_upper, only : upper
    use m_XAS_functions, only : calculate_dipoles_XAS
    use m_XAS_functions, only : compute_XAS_spectrum
    
    type(sckh_params_t), intent(inout):: p
    
    character(80):: file,string
    real(wp):: mu_SI, dx,dvr_start, gamma, E_n_mean 
    integer:: npoints
    real(wp), dimension(:),allocatable:: X_r,E_i, eig_i, eig_na
    real(wp), dimension(:),allocatable::  sigma, omega_in
    real(wp), dimension(:,:),allocatable::  c_i, &
         eig_n, E_n, sigma_states
    real(wp), allocatable:: E_n0(:)
    real(wp), dimension(:,:,:),allocatable:: c_n, dipole_n, D_ni, c_na
    real(wp), dimension(:,:,:,:),allocatable::  nac
    real(wp), dimension(:,:),allocatable::  M_dsinc
    real(wp), allocatable:: sigma_final(:,:,:,:)
    real(wp), allocatable:: sigma_tensor(:,:,:)
    real(wp):: gamma_instr, gamma_inc 
    integer::i,j,ii,jj,k,l,m, ifile
    real(wp):: norm
    integer:: ind(1)
    
    !
    ! This progam calculates the XAS cross section including vibrational effects
    !
    
    gamma = p % gamma_FWHM /  2.0_wp 
    gamma_instr = p % gamma_instr_FWHM /  2.0_wp
    gamma_inc = p % gamma_inc_FWHM /  2.0_wp
    
    !nstates_na = p % npesfile_n *p % nstates
    
    ! allocate everything
    allocate( X_r(p % nstates), E_i(p % nstates), &
         E_n(p % npesfile_n,p % nstates), eig_i(p % nstates), eig_n(p % npesfile_n,p % nstates))
    allocate( omega_in(p % n_omega_in), sigma(p % n_omega_in) )
    allocate(c_i(p % nstates,p % nstates),c_n(p % npesfile_n,p % nstates,p % nstates) )
    allocate(dipole_n(p % npesfile_n, p % nstates, 3), sigma_states(p % npesfile_n,p % n_omega_in))
    allocate(D_ni(p % npesfile_n,p % nstates,3) )
    allocate(sigma_final(p % npesfile_n, p % n_omega_in, 3, 3))
    allocate(sigma_tensor(p % n_omega_in, 3, 3))
    allocate(E_n0(p % nstates))
    
    npoints = (p % nstates-1)/2
    mu_SI = p % mu * const % u
    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    
    ! set up grid points
    do i = -npoints,npoints
      ii = i + npoints +1
      X_r(ii) = (ii-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % nstates, X_r, E_i)

    ind = minloc(E_i)

    ! intermediate state reference energy (lowest state) 
    if(p % use_n0_state) then 
      call read_PES_file(p % pes_file_n, p % npoints_in, p % nstates, X_r, E_n0)
    else
      E_n0 =0.0d0
    end if

    ! read list of final state pes_files and dipole_files
    call read_file_list(p % pes_file_list_n, p % npesfile_n, p % pes_files_n)
    call read_file_list(p % dipole_file_list_n, p % npesfile_n, p % dipolefile_n)
    
    do j=1,p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, p % nstates, X_r, E_n(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % nstates, X_r, dipole_n(j,:,:))
    end do
    
    !create omega_in
    call linspace(omega_in, p % omega_in_start,p % omega_in_end, p % n_omega_in ) 
    
    !
    ! Solve the vibrational problem for all eigenfunctions
    !
    
    ! initial state
    call solve_vib_problem(dx, E_i, eig_i, c_i, mu_SI, p % vib_solver)
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1)) * const % cm
    
    ! final states
    do j=1,p % npesfile_n
      call solve_vib_problem(dx, E_n0 + E_n(j,:), eig_n(j,:), c_n(j,:,:), mu_SI, p % vib_solver)
      write(6,*) "Final state fundamental",j, (eig_n(j,2) -eig_n(j,1)) * const % cm
    end do
    
    ! calculate transition dipoles
    call calculate_dipoles_XAS( c_i, c_n, dipole_n, D_ni, p % dipole_mode, ind(1))
    
    ! convert eigenvalues to eV units
    eig_i =eig_i / const % eV
    eig_n =eig_n / const % eV
    
    write(6,*) eig_n(1,1) - eig_i(1)
    
    do j=1,p % npesfile_n    
      call compute_XAS_spectrum(eig_i(1), eig_n(j,:), D_ni(j,:,:), &
           omega_in, gamma, gamma_inc, sigma_final(j,:,:,:))
    end do

    write(6,*) "Calculated XAS spectrum"
    
    sigma_states(:,:) = 0.0_wp
    do m=1,3
      sigma_states(:,:) = sigma_states(:,:) + sigma_final(:,:,m,m)
    end do

    sigma_tensor(:,:,:) = sum(sigma_final(:,:,:,:),1)
    
    sigma = sum(sigma_states,1)
    
    norm = sum(sigma) *(omega_in(2)-omega_in(1))
    sigma = sigma / norm
    sigma_states = sigma_states / norm
    sigma_tensor = sigma_tensor / norm
    
    !
    ! Write to files
    ! 
    
    ! write spectra to file
    file="_sigma"
    file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
    
    ifile = get_free_handle()
    open(ifile,file=file,status='unknown')
    
    do i=1, p % n_omega_in
      write(ifile,'(2ES18.10)') omega_in(i), sigma(i) !, sigma_tensor(i,:,:)
    end do
    
    close(ifile) 
    
    
    do j= 1,p % npesfile_n
      file="vib_eigenvalues_final_"
      write(string,*) j
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1,p % nstates
        write(ifile,*) eig_n(j,i), D_ni(j,i,1), D_ni(j,i,2), D_ni(j,i,3)
      end do
      
      write(ifile,*) 
      
      close(ifile)
    end do
  
  end subroutine calculate_XAS

!  subroutine calculate_XAS_nonadiabatic(p)
!    use m_precision, only: wp
!    use m_sckh_params_t, only: sckh_params_t 
!    use m_constants, only: const
!    use m_splines, only: spline_easy
!    use m_splines, only: linspace
!    use m_PES_io, only: read_PES_file
!    use m_PES_io, only: read_dipole_file
!    use m_PES_io, only: read_file_list
!    use m_io, only: get_free_handle
!    use m_KH_functions, only: solve_vib_problem
!    use m_upper, only : upper
!    use m_XAS_functions, only : calculate_dipoles_XAS
!    use m_XAS_functions, only : compute_XAS_spectrum
!    
!    type(sckh_params_t), intent(inout):: p
!    
!    character(80):: file,string
!    real(wp):: mu_SI, dx,dvr_start, gamma, E_n_mean 
!    integer:: npoints, nstates_na 
!    real(wp), dimension(:),allocatable:: X_r,E_i, eig_i, eig_na
!    real(wp), dimension(:),allocatable::  sigma, omega_in
!    real(wp), dimension(:,:),allocatable::  c_i, &
!         eig_n, E_n, sigma_states
!    real(wp), dimension(:,:,:),allocatable:: c_n, dipole_n, D_ni, c_na
!    real(wp), dimension(:,:,:,:),allocatable::  nac
!    real(wp), dimension(:,:),allocatable::  M_dsinc
!    real(wp), allocatable:: sigma_final(:,:,:,:)
!    real(wp), allocatable:: sigma_tensor(:,:,:)
!    real(wp):: gamma_instr, gamma_inc 
!    integer::i,j,ii,jj,k,l,m, ifile
!    real(wp):: norm
!    
!    !
!    ! This progam calculates the XAS cross section including non-adiabatic couplings
!    !
!    
!    gamma = p % gamma_FWHM /  2.0_wp 
!    gamma_instr = p % gamma_instr_FWHM /  2.0_wp
!    gamma_inc = p % gamma_inc_FWHM /  2.0_wp
!    
!    nstates_na = p % npesfile_n *p % nstates
!    
!    ! allocate everything
!    allocate( X_r(p % nstates), E_i(p % nstates), &
!         E_n(p % npesfile_n,p % nstates), eig_i(p % nstates), eig_n(p % npesfile_n,p % nstates))
!    allocate( omega_in(p % n_omega_in), sigma(p % n_omega_in) )
!    allocate(c_i(p % nstates,p % nstates),c_n(p % npesfile_n, p % nstates,p % nstates) )
!    allocate(dipole_n(p % npesfile_n, p % nstates, 3), sigma_states(p % npesfile_n,p % n_omega_in))
!    allocate(D_ni(p % npesfile_n,p % nstates,3) )
!    allocate(sigma_final(p % npesfile_n, p % n_omega_in, 3, 3))
!    allocate(sigma_tensor(p % n_omega_in, 3, 3))
!
!    !if (p % nonadiabatic .eq. 1) then
!    allocate(nac(p % npesfile_n, p % npesfile_n, p % nstates, 2) )
!    allocate(eig_na(nstates_na), c_na(nstates_na, nstates_na) )
!    allocate(M_dsinc(p % nstates, p % nstates) )
!    ! end if
!    
!    npoints = (p % nstates-1)/2
!    mu_SI = p % mu * const % u
!    dvr_start = p % dvr_start_in * 1.0d-10
!    dx = p % dx_in * 1.0d-10
!    
!    !
!    ! set up DVR points
!    !
!    
!    do i = -npoints,npoints
!      ii = i + npoints +1
!      X_r(ii) = (ii-1)*dx + dvr_start
!    end do
!    
!    ! read PES files
!    call read_PES_file(p % pes_file_i, p % npoints_in, p % nstates, X_r, E_i)
!
!    ! read list of final state pes_files and dipole_files
!    call read_file_list(p % pes_file_list_n, p % npesfile_n, p % pes_files_n)
!    call read_file_list(p % dipole_file_list_n, p % npesfile_n, p % dipolefile_n)
!    
!    do j=1,p % npesfile_n
!      call read_PES_file(p % pes_files_n(j), p % npoints_in, p % nstates, X_r, E_n(j,:))
!      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % nstates, X_r, dipole_n(j,:,:))
!    end do
!    
!    !  if (p % nonadiabatic .eq. 1) then
!    !     !call read_nac_file(p % nac_file, p % npoints_in, p %nstates, X_r, p % npesfile_f, nac)
!    !  end if
!    
!    !create omega_in
!    call linspace(omega_in, p % omega_in_start,p % omega_in_end, p % n_omega_in ) 
!    
!    !
!    ! Solve the vibrational problem for all eigenfunctions
!    !
!    
!    ! initial state
!    call solve_vib_problem(dx, E_i, eig_i, c_i, mu_SI, p % vib_solver)
!    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1)) * const % cm
!    
!    ! final states
!    do j=1,p % npesfile_n
!      call solve_vib_problem(dx, E_n(j,:), eig_n(j,:), c_n(j,:,:), mu_SI, p % vib_solver)
!      write(6,*) "Final state fundamental",j, (eig_n(j,2) -eig_n(j,1)) * const % cm
!    end do
!   
!    ! calculate transition dipoles
!    call calculate_dipoles_XAS( c_i, c_n, dipole_n, D_ni)
!
!    ! solve non-adiabatic problem
!    !if (p % nonadiabatic .eq. 1) then  
!    nac = 0.0_wp
!    !nac(1,2,:,2) = 1d-19     
!    !nac(2,1,:,2) = 1d-19     
!    
!    ! compute matrix elements <x_i | d/dx | x_j >
!    !XXX not implemented
!    call compute_d1_mat_elems(X_r, M_d1, p % vib_solver)
!
!    !implemented ?
!    call solve_non_adiabatic(eig_n, c_n, nac, M_d1, eig_na, c_na)
!    
!    ! convert eigenvalues to eV units
!    eig_i =eig_i / const % eV
!    eig_n =eig_n / const % eV
!    
!    write(6,*) eig_n(1,1) - eig_i(1)
!    
!    ! calculate transition dipoles for non-adiabatic transitions
!    !XXX not implemented
!    call calculate_dipoles_nac( D_ni, c_na, D_na_i)
!    
!    do j=1,p % npesfile_n    
!      call compute_XAS_spectrum(eig_i(1), eig_na(j,:), D_na_i(j,:,:), &
!           omega_in, gamma, gamma_inc, sigma_final(j,:,:,:))
!    end do
!
!    
!    write(6,*) "Calculated non-adiabatic XAS spectrum"
!    
!    sigma_states(:,:) = 0.0_wp
!    do m=1,3
!      sigma_states(:,:) = sigma_states(:,:) + sigma_final(:,:,m,m)
!    end do
!
!    sigma_tensor(:,:,:) = sum(sigma_final(:,:,:,:),1)
!    
!    sigma = sum(sigma_states,1)
!    
!    norm = sum(sigma) *(omega_in(2)-omega_in(1))
!    sigma = sigma / norm
!    sigma_states = sigma_states / norm
!    sigma_tensor = sigma_tensor / norm
!    
!    !
!    ! Write to files
!    ! 
!    
!    ! write spectra to file
!    file="_sigma"
!    file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
!    
!    ifile = get_free_handle()
!    open(ifile,file=file,status='unknown')
!    
!    do i=1, p % n_omega_in
!      write(ifile,'(2ES18.10)') omega_in(i), sigma(i) !, sigma_tensor(i,:,:)
!    end do
!    
!    close(ifile) 
!    
!    
!    do j= 1,p % npesfile_n
!      file="vib_eigenvalues_final_"
!      write(string,*) j
!      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
!      
!      ifile = get_free_handle()
!      open(ifile,file=file,status='unknown')
!      
!      do i=1,p % nstates
!        write(ifile,*) eig_n(j,i)
!      end do
!      
!      write(ifile,*) 
!      
!      close(ifile)
!    end do
!    
!  end subroutine calculate_XAS_nonadiabatic


  
end module m_XAS_eigenstates


