module m_XAS_eigenstates
!  use parameters
!  use m_XAS_io
!  use KH_functions
!  use spline_m
!  use m_XAS_functions
  implicit none

contains
  
  subroutine calculate_XAS_nonadiabatic(p)
    use m_precision, only: wp
    use m_sckh_params_t, only: sckh_params_t 
    use m_constants, only: const
    use m_splines, only: spline_easy
    use m_splines, only: linspace
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: read_dipole_file
    !use m_PES_io, only: read_nac_file
    !use m_KH_utils, only: calculate_dipoles_KH_res
    !use m_KH_utils, only: compute_XES_res
    !use m_KH_utils, only: compute_XES_res_alt
    use m_io, only: get_free_handle
    use m_KH_functions, only: solve_vib_problem
    use m_upper, only : upper
    use m_XAS_functions, only : calculate_dipoles_XAS
    use m_XAS_functions, only : spectrum_XAS
    
    type(sckh_params_t), intent(inout):: p
    
    character(80):: file,string
    real(kind=wp):: mu_SI, dx,dvr_start, gamma, E_n_mean 
    integer:: npoints, nstates_na 
    real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, eig_i, eig_na
    real(kind=wp), dimension(:),allocatable::  sigma, omega_in
    real(kind=wp), dimension(:,:),allocatable::  c_i, &
         eig_f, E_f, sigma_states
    real(kind=wp), dimension(:,:,:),allocatable:: c_f, dipole_f, D_fi, c_na
    real(kind=wp), dimension(:,:,:,:),allocatable::  nac
    real(kind=wp), dimension(:,:),allocatable::  M_dsinc
    real(kind=wp):: gamma_instr, gamma_inc 
    integer::i,j,ii,jj,k,l,m, ifile
    
    !
    ! This progam calculates the XAS cross section including non-adiabatic couplings
    !
    
    gamma = p % gamma_FWHM /  2.0_wp 
    gamma_instr = p % gamma_instr_FWHM /  2.0_wp
    gamma_inc = p % gamma_inc_FWHM /  2.0_wp
    
    nstates_na = p % npesfile_f *p % nstates
    
    ! allocate everything
    allocate( X_dvr(p % nstates), E_i(p % nstates), &
         E_f(p % npesfile_f,p % nstates), eig_i(p % nstates), eig_f(p % npesfile_f,p % nstates))
    allocate( omega_in(p % n_omega_in), sigma(p % n_omega_in) )
    allocate(c_i(p % nstates,p % nstates),c_f(p % npesfile_f,p % nstates,p % nstates) )
    allocate(dipole_f(p % npesfile_f, p % nstates, 3), sigma_states(p % npesfile_f,p % n_omega_in))
    allocate(D_fi(p % npesfile_f,p % nstates,3) )
    
    !if (p % nonadiabatic .eq. 1) then
    allocate(nac(p % npesfile_f, p % npesfile_f, p % nstates, 2) )
    allocate(eig_na(nstates_na), c_na(p % npesfile_f, p % nstates, nstates_na) )
    allocate( M_dsinc(p % nstates, p % nstates) )
    ! end if
    
    npoints = (p % nstates-1)/2
    mu_SI = p % mu * const % u
    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    
    !
    ! set up DVR points
    !
    
    do i = -npoints,npoints
      ii = i + npoints +1
      X_dvr(ii) = (ii-1)*dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % nstates, X_dvr, E_i)

    ! final state pes_files and dipole_files
    ifile = get_free_handle()
    open(ifile, file= p % pes_file_list_f, action='read')
    
    allocate(p % pes_files_f(p % npesfile_f))
    do i=1, p % npesfile_f
      read(ifile,*) p % pes_files_f(i)
      write(6,*) p % pes_files_f(i)
    end do
    
    close(ifile)
    
    ifile = get_free_handle()
    open(ifile, file= p % dipole_file_list_f, action='read')
    
    write(6,*) "p % dipole_file_list_f", p % dipole_file_list_f
    allocate(p % dipolefile_f(p % npesfile_f))
    do i=1, p % npesfile_f
      read(ifile,*) p % dipolefile_f(i)
    end do

    close(ifile)
    
    do j=1,p % npesfile_f
      call read_PES_file(p % pes_files_f(j), p % npoints_in, p % nstates, X_dvr, E_f(j,:))
      call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % nstates, X_dvr, dipole_f(j,:,:))
    end do
    
    !  if (p % nonadiabatic .eq. 1) then
    !     !call read_nac_file(p % nac_file, p % npoints_in, p %nstates, X_dvr, p % npesfile_f, nac)
    !  end if
    
    !create omega_in
    call linspace(omega_in, p % omega_in_start,p % omega_in_end, p % n_omega_in ) 
    
    !
    ! Solve the vibrational problem for all eigenfunctions
    !
    
    
    ! initial state
    !call solve_sinc_DVR(dx,mu_SI, E_i, c_i, eig_i)
    call solve_vib_problem(dx, E_i, eig_i, c_i, mu_SI, p % vib_solver)
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1)) * const % cm
    
    ! final states
    do j=1,p % npesfile_f
      !call solve_sinc_DVR(dx,mu_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))
      call solve_vib_problem(dx, E_f(j,:), eig_f(j,:), c_f(j,:,:), mu_SI, p % vib_solver)
      write(6,*) "Final state fundamental",j, (eig_f(j,2) -eig_f(j,1)) * const % cm
    end do
    
    !  ! solve non-adiabatic problem
    !  if (p % nonadiabatic .eq. 1) then  
    !     nac = 0.0_wp
    !     nac(1,2,:,2) = 1d-19     
    !     nac(2,1,:,2) = 1d-19     
    !
    !     ! compute matrix elements numerically, could take time... 
    !     call dsinc_mat_elems(M_dsinc, X_dvr)
    !
    !     call solve_non_adiabatic(eig_f, c_f, nac, eig_na, c_na, M_dsinc)
    !  end if
    
    ! calculate transition dipoles
    call calculate_dipoles_XAS( c_i, c_f, dipole_f, D_fi)
    
    ! convert eigenvalues to eV units
    eig_i =eig_i / const % eV
    eig_f =eig_f / const % eV
    eig_na =eig_na / const % eV
    
    write(6,*) eig_i(1) -eig_f(1,1)
    
    ! calculate spectrum
    if (p % nonadiabatic .eq. 1) then
      !   call spectrum_XAS_nonadiabatic(eig_na, eig_i, c_na, D_fi, omega_in, sigma, sigma_states, gamma)
    else
      call spectrum_XAS(eig_f, eig_i, D_fi, omega_in, sigma, sigma_states, gamma)
    end if
    
    !
    ! Write to files
    ! 
    
    ! write spectra to file
    file="_sigma_"
    file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
    
    ifile = get_free_handle()
    open(ifile,file=file,status='unknown')
    
    do i=1, p % n_omega_in
      write(ifile,'(2ES18.10)') omega_in(i), sigma(i) !, sigma_tensor(i,:,:)
    end do
    
    close(ifile) 
    
    
    do j= 1,p % npesfile_f
      file="vib_eigenvalues_final_"
      write(string,*) j
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1,p % nstates
        write(ifile,*) eig_f(j,i)
      end do
      
      write(ifile,*) 
      
      close(ifile)
    end do
  
end subroutine calculate_XAS_nonadiabatic

end module m_XAS_eigenstates


