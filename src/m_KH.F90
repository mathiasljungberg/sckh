

module m_KH
  implicit none

#include "m_define_macro.F90"
contains

  subroutine calculate_KH_nonres(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_splines, only: spline_easy
    use m_splines, only: linspace
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_nac_file
    use m_sckh_params_t, only: sckh_params_t 
    use m_KH_functions, only: solve_sinc_DVR
    use m_KH_utils, only: calculate_dipoles_KH_nonres
    use m_KH_utils, only: spectrum_XES
    use m_KH_utils, only: compute_XES_nonres
    use m_io, only: get_free_handle

    type(sckh_params_t), intent(inout):: p 

    ! loop variables
    integer::i,j,ii,jj,k,l,m,t 

    ! other variables
    character(80):: file,string
    real(kind=wp):: mu_SI, dx,dvr_start, gamma, E_n_mean 
    integer:: npoints 
    real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, E_n, E_lp_corr, eig_i, eig_n, shift
    real(kind=wp), dimension(:),allocatable::  sigma, omega
    real(kind=wp), dimension(:,:),allocatable::  c_i, c_n, &
         eig_f, E_f, sigma_states, D_ni
    real(kind=wp), dimension(:,:,:),allocatable:: c_f, dipole, D_fi
    real(kind=wp), dimension(:,:,:,:),allocatable:: D_fn, nac
    real(wp), allocatable:: sigma_final(:,:,:,:)
    integer:: ifile
    real(wp):: norm


    !
    ! This progam calculates the KH emission spectrum using the eigenstate basis
    !

    gamma = p % gamma_FWHM /  2 

    ! allocate everything
    allocate( X_dvr(p % nstates), E_i(p % nstates), E_n(p % nstates), E_lp_corr(p % nstates), &
         E_f(p % npesfile_f,p % nstates), eig_i(p % nstates),eig_n(p % nstates),eig_f(p % npesfile_f,p % nstates), &
         shift(p % nstates))
    allocate( omega(p % n_omega), sigma(p % n_omega), D_ni(p % nstates, 3) )
    allocate(c_i(p % nstates,p % nstates),c_f(p % npesfile_f,p % nstates,p % nstates),c_n(p % nstates,p % nstates), &
         D_fn(p % npesfile_f,p % nstates,p % nstates,3))
    allocate(dipole(p % npesfile_f, p % nstates, 3), sigma_states(p % npesfile_f,p % n_omega))
    allocate(D_fi(p % npesfile_f,p % nstates,3) )
    allocate(sigma_final(p % npesfile_f,p % n_omega,3,3))

    if (p % nonadiabatic .eq. 1) then
      allocate(nac(p % npesfile_f, p % npesfile_f, p % nstates, 2) )
    end if

    if (mod(p % nstates,2).ne.1 ) then
      write(6,*) "nstates must be an odd number"
      stop
    end if

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

    write(6,*) "X_dvr", X_dvr

    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % nstates, X_dvr, E_i)
    call read_PES_file(p % pes_file_n, p % npoints_in, p % nstates, X_dvr, E_n)

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
      call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % nstates, X_dvr, dipole(j,:,:))
    end do

    !  if (p % nonadiabatic .eq. 1) then
    !     call read_nac_file(p % nac_file, p % npoints_in, p %nstates, X_dvr, p % npesfile_f, nac)
    !  end if

    ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved

    if( p % shift_PES .eq.  1) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % nstates, X_dvr, E_lp_corr)

      shift = E_lp_corr -E_f(1,:) 

      do j=1,p % npesfile_f
        E_f(j,:) = E_f(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if

    !create omega
    call linspace(omega, p % omega_start,p % omega_end, p % n_omega ) 

    !
    ! Solve the vibrational problem for all eigenfunctions
    !


    ! initial state
    call solve_sinc_DVR(dx,mu_SI, E_i, c_i, eig_i)
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm

    ! intermediate state
    call solve_sinc_DVR(dx,mu_SI, E_n, c_n, eig_n)
    write(6,*) "Intermediate state fundamental", (eig_n(2) -eig_n(1))*const % cm

    ! final states
    do j=1,p % npesfile_f
      call solve_sinc_DVR(dx,mu_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))
      write(6,*) "Final state fundamental",j, (eig_f(j,2) -eig_f(j,1))*const % cm
    end do

    ! solve non-adiabatic problem
    !call solve_non_adiabatic(eig_f, c_f, eig_na, c_na)

    ! calculate transition dipoles
    call calculate_dipoles_KH_nonres( c_i, c_n, c_f, dipole, D_ni, D_fn, D_fi)

    ! convert eigenvalues to eV units
    eig_i =eig_i / const % eV
    eig_n =eig_n / const % eV
    eig_f =eig_f / const % eV

    write(6,*) eig_n(1) -eig_f(1,1)

    ! calculate "mean energy" for intermediate state
    E_n_mean = sum(eig_n(:) * D_ni(:,1) ** 2) 

    write(6,*) "sum (D_ni(:,1)**2)", sum( D_ni(:,1) ** 2)
    write(6,*) "E_n_mean", E_n_mean

    !
    ! calculate spectrum
    !

    if (p % nonadiabatic .eq. 1) then
      !call spectrum_XES_nonadiabatic(eig_na, eig_n, c_na, D_ni, D_fn, D_fi, omega, sigma, sigma_states)
      stop
    else
      !call spectrum_XES(eig_f, eig_n, D_ni, D_fn, D_fi, omega, sigma, sigma_states, gamma)

      do j=1,p % npesfile_f    
        call compute_XES_nonres(eig_i(1), eig_n, eig_f(j,:), D_ni(:,:), D_fn(j,:,:,:), &
             omega, gamma, sigma_final(j,:,:,:))
      end do

      sigma_states(:,:) = 0.0_wp
      do m=1,3
        sigma_states(:,:) = sigma_states(:,:) + sigma_final(:,:,m,m)
      end do

      sigma = sum(sigma_states,1)

      norm = sum(sigma) *(omega(2)-omega(1))
      sigma = sigma / norm
      sigma_states = sigma_states / norm

    end if

    !
    ! Write to files
    ! 

    ! write sigma to file
    file="_sigma"
    file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"

    open(10,file=file,status='unknown')

    do i=1,p % n_omega
      write(10,*) omega(i), sigma(i)
    end do

    close(10) 


    !! write sigma_dir to file
    !
    !file="_sigma_direct"
    !file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"
    !
    !!write(string,*) i
    !!file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
    !
    !open(10,file=file,status='unknown')
    !
    !do i=1,p % n_omega
    !   write(10,*) omega(i), sigma_dir(i)
    !end do
    !
    !close(10)
    !
    !! write sigma_max_int to file
    !file="sigma_max_int"
    !file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"
    !
    !open(10,file=file,status='unknown')
    !
    !do i=1,p % n_omega
    !   write(10,*) omega(i), sigma_max_int(i)
    !end do
    !
    !close(10)
    !
    !! write sigma_states
    !do j= 1,p % npesfile_f
    !   file="_sigma_states_"
    !   write(string,*) j
    !   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
    !
    !   open(10,file=file,status='unknown')
    !
    !   do i=1,p % n_omega
    !      write(10,*) omega(i), sigma_states(j,i)
    !   end do
    !
    !   close(10)
    !end do
    !
    !! write sigma_dir_states
    !do j= 1,p % npesfile_f
    !   file="_sigma_dir_states_"
    !   write(string,*) j
    !   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
    !
    !   open(10,file=file,status='unknown')
    !
    !   do i=1,p % n_omega
    !      write(10,*) omega(i), sigma_dir_states(j,i)
    !   end do
    !
    !   close(10)
    !end do
    !
    !! write sigma_max_int_states
    !do j= 1,p % npesfile_f
    !   file="_sigma_max_int_states_"
    !   write(string,*) j
    !   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
    !
    !   open(10,file=file,status='unknown')
    !
    !   do i=1,p % n_omega
    !      write(10,*) omega(i), sigma_max_int_states(j,i)
    !   end do
    !
    !   close(10)
    !end do
    !
    !! write vibrational eigenvalues
    !file="_vib_eigenvalues_initial"
    !file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"
    !
    !open(10,file=file,status='unknown')
    !do i=1,p % nstates
    !   write(10,*) eig_i(i), 1 ,10
    !end do
    !close(10)
    !
    !file="_vib_eigenvalues_intermediate"
    !file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"
    !
    !open(10,file=file,status='unknown')
    !do i=1,p % nstates
    !   write(10,*) eig_n(i), 1 , D_ni(i)
    !end do
    !close(10)


    !! a lot of output disabled
    !if (0) then
    !   
    !   do j= 1,p % npesfile_f
    !      file="_vib_eigenvalues_final_"
    !      write(string,*) j
    !      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
    !      
    !      open(10,file=file,status='unknown')
    !      
    !      do i=1,p % nstates
    !         write(10,*) eig_f(j,i), 1, 10
    !      end do
    !      
    !      close(10)
    !   end do
    !   
    !
    !
    !   ! write transition moments
    !   file="_transition_moments_i_to_n"
    !   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"
    !   
    !   open(10,file=file,status='unknown')
    !   do i=1,p % nstates
    !      write(10,*) D_ni(i)
    !   end do
    !   close(10)
    !   
    !   !file="transition_moments_n_to_f"
    !   !do i=1,p % npesfile_f
    !   !   do j=1,p % nstates ! final
    !   !      do k=1,p % nstates ! intermediate
    !   !         do l=1,3
    !   !            D_fn(i,j,k,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
    !   !         end do
    !   !      end do
    !   !   end do
    !   !end  do
    !        
    !
    !
    !   ! write transition energies for the 10 first intermediate states
    !   do j= 1,10
    !      file="_delta_e_lp_"
    !      write(string,*) j
    !      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
    !      
    !      open(10,file=file,status='unknown')
    !      
    !      do i=1,p % nstates
    !         write(10,*) delta_e_lp(j,i)
    !      end do
    !      close(10)
    !   end do
    !   
    !end if
    !
    !! write vibrational eigenfucntions, normally disabled
    !if(0) then 
    !
    !   file="_vib_eigenfunctions_initial"
    !   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"     
    !
    !   open(10,file=file,status='unknown')
    !   do i=1,p % nstates
    !      do j=1,p % nstates
    !         write(10,*) X_dvr(j), c_i(j,i)
    !      end do
    !      write(10,*) 
    !      write(10,*) 
    !   end do
    !   close(10)
    !
    !   file="_vib_eigenfunctions_intermediate"
    !   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) //  ".dat"
    !
    !  open(10,file=file,status='unknown')
    !   do i=1,p % nstates
    !      do j=1,p % nstates
    !         write(10,*) X_dvr(j), c_n(j,i)
    !      end do
    !      write(10,*) 
    !      write(10,*) 
    !   end do
    !   close(10)
    !
    !   ! write 10 first eigenfunctions to seprate files
    !   do i=1,10
    !      file="_vib_eigenfunctions_intermediate_"
    !      write(string,*) i
    !      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
    !      
    !      open(10,file=file,status='unknown')
    !
    !      do j=1,p % nstates
    !         write(10,*) X_dvr(j), c_n(j,i)
    !      end do
    !      write(10,*)
    !      write(10,*)
    !      close(10)
    !   end do! i
    !
    !   ! write 10 first eigenfunctions of first final state to separate files
    !   do i=1,10
    !      file="_vib_eigenfunctions_final_1_"
    !      write(string,*) i
    !      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
    !
    !      open(10,file=file,status='unknown')
    !
    !      do j=1,p % nstates
    !         write(10,*) X_dvr(j), c_f(1,j,i)
    !      end do
    !      write(10,*)
    !      write(10,*)
    !      close(10)
    !   end do
    !
    !   do k= 1,p % npesfile_f
    !      file="_vib_eigenfunctions_final_"
    !      write(string,*) k
    !      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
    !
    !      open(10,file=file,status='unknown')
    !
    !      do i=1,p % nstates
    !         do j=1,p % nstates
    !            write(10,*) X_dvr(j), c_f(k,j,i)
    !         end do
    !
    !         write(10,*) 
    !         write(10,*) 
    !
    !      end do
    !      close(10)
    !   end do
    !
    !
    !end if ! if(1) 

    do j= 1,p % npesfile_f
      file="vib_eigenvalues_final_"
      write(string,*) j
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

      open(10,file=file,status='unknown')

      do i=1,p % nstates
        write(10,*) eig_f(j,i)
      end do

      write(10,*) 
      write(10,*) 

      close(10)
    end do

  end subroutine calculate_KH_nonres

  subroutine calculate_KH_res(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_splines, only: spline_easy
    use m_splines, only: linspace
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_nac_file
    use m_sckh_params_t, only: sckh_params_t 
    !use m_KH_functions, only: solve_sinc_DVR
    use m_KH_utils, only: calculate_dipoles_KH_res
    use m_KH_utils, only: compute_XES_res
    use m_KH_utils, only: compute_XES_res_alt
    use m_io, only: get_free_handle
    !use m_fourier_grid, only: solve_fourier_grid_real
    !use m_upper, only : upper
    use m_KH_functions, only: solve_vib_problem
    type(sckh_params_t), intent(inout):: p 

    ! loop variables
    integer::i,j,ii,jj,k,l,m,t 

    ! other variables
    character(80):: file,string
    real(kind=wp):: mu_SI, dx,dvr_start, gamma, gamma_instr, gamma_inc, E_n_mean 
    integer:: npoints 
    real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, E_lp_corr, eig_i, shift
    real(kind=wp), dimension(:),allocatable::  sigma, omega, omega_in
    real(kind=wp), dimension(:,:),allocatable::  c_i, &
         eig_f, E_f, sigma_states, E_n, eig_n
    real(kind=wp), dimension(:,:,:),allocatable:: c_n, c_f, dipole_f, dipole_n, D_fi, D_ni
    !real(kind=wp), dimension(:,:,:,:),allocatable:: D_fn, nac
    real(kind=wp), allocatable::D_fn(:,:,:,:,:)
    real(wp), allocatable:: sigma_final(:,:,:,:,:)
    integer:: ifile
    real(wp):: norm
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    real(wp), allocatable:: eig_n_cmp(:), D_ni_cmp(:,:), D_fn_cmp(:,:,:,:)
    integer:: n_e, n_v, n_ev, f_e
    
    !
    ! This progam calculates the KH emission spectrum using the eigenstate basis
    !

    gamma = p % gamma_FWHM /  2
    gamma_instr = p % gamma_instr_FWHM /  2
    gamma_inc = p % gamma_inc_FWHM /  2 

    ! allocate everything
    allocate( X_dvr(p % nstates), E_i(p % nstates), E_n(p % npesfile_n, p % nstates), E_lp_corr(p % nstates), &
         E_f(p % npesfile_f,p % nstates), eig_i(p % nstates),eig_n(p % npesfile_n, p % nstates),&
         eig_f(p % npesfile_f,p % nstates), &
         shift(p % nstates))
    allocate( omega(p % n_omega), sigma(p % n_omega), D_ni(p % npesfile_n, p % nstates, 3) )
    allocate( omega_in(p % n_omega_in))
    allocate(c_i(p % nstates,p % nstates),c_f(p % npesfile_f,p % nstates,p % nstates), &
         c_n(p % npesfile_n, p % nstates,p % nstates), &
         D_fn(p % npesfile_f,p % nstates,p % npesfile_n, p % nstates,3))
    allocate(dipole_f(p % npesfile_f, p % nstates, 3), sigma_states(p % npesfile_f,p % n_omega))
    allocate(dipole_n(p % npesfile_n, p % nstates, 3))
    allocate(D_fi(p % npesfile_f,p % nstates,3) )
    write(6,*) "", p % npesfile_f, p % n_omega_in,  p % n_omega 
    allocate(sigma_final(p % npesfile_f, p % n_omega_in, p % n_omega,3,3))
    allocate(eig_n_cmp(p % npesfile_n * p % nstates))
    allocate(D_ni_cmp(p % npesfile_n * p % nstates, 3))
    allocate(D_fn_cmp(p % npesfile_f,p % nstates, p % npesfile_n * p % nstates, 3))
    allocate(lambda_F(p % npesfile_f, p % n_omega_in, p % n_omega))
    allocate(lambda_G(p % npesfile_f, p % n_omega_in, p % n_omega))
    allocate(lambda_H(p % npesfile_f, p % n_omega_in, p % n_omega))

    if (p % nonadiabatic .eq. 1) then
    !  allocate(nac(p % npesfile_f, p % npesfile_f, p % nstates, 2) )
    end if

    if (mod(p % nstates,2).ne.1 ) then
      write(6,*) "nstates must be an odd number"
      stop
    end if

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

    write(6,*) "X_dvr", X_dvr

    !
    ! read PES files
    !

    ! iniital state
    call read_PES_file(p % pes_file_i, p % npoints_in, p % nstates, X_dvr, E_i)

    ! intermediate states
    ifile = get_free_handle()
    open(ifile, file= p % pes_file_list_n, action='read')

    allocate(p % pes_files_n(p % npesfile_n))
    do i=1, p % npesfile_n
      read(ifile,*) p % pes_files_n(i)
      write(6,*) p % pes_files_n(i)
    end do
    
    close(ifile)

    ifile = get_free_handle()
    open(ifile, file= p % dipole_file_list_n, action='read')

    write(6,*) "p % dipole_file_list_n  ", p % dipole_file_list_n
    allocate(p % dipolefile_n(p % npesfile_n))
    do i=1, p % npesfile_n
      read(ifile,*) p % dipolefile_n(i)
    end do

    close(ifile)
    
    do j=1,p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, &
           p % nstates, X_dvr, E_n(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % nstates, X_dvr, dipole_n(j,:,:))
    end do

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

    write(6,*) "p % dipole_file_list_f  ", p % dipole_file_list_f
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
    !     call read_nac_file(p % nac_file, p % npoints_in, p %nstates, X_dvr, p % npesfile_f, nac)
    !  end if

    ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved

    if( p % shift_PES .eq.  1) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % nstates, X_dvr, E_lp_corr)

      shift = E_lp_corr -E_f(1,:) 

      do j=1,p % npesfile_f
        E_f(j,:) = E_f(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if


    !create omega_in
    call linspace(omega_in, p % omega_in_start,p % omega_in_end, p % n_omega_in )

    !create omega_out
    call linspace(omega, p % omega_start,p % omega_end, p % n_omega ) 


    !
    ! Solve the vibrational problem for all eigenfunctions
    !

    ! initial state
    call solve_vib_problem(dx, E_i, eig_i, c_i, mu_SI, p % vib_solver)
    
    !if (upper(p % vib_solver) .eq. "FOURIER_REAL") then
    !  write(6,*) "FOURIER_REAL vibrational solver"
    !  call solve_fourier_grid_real(dx, E_i, eig_i, c_i, mu_SI, "SI", "fast")
    !else 
    !  write(6,*) "DVR vibrational solver, ", upper(p % vib_solver)
    !  call solve_sinc_DVR(dx,mu_SI, E_i, c_i, eig_i)
    !end if
    
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm

    ! intermediate states
    do j=1,p % npesfile_n

      call solve_vib_problem(dx, E_n(j,:), eig_n(j,:), c_n(j,:,:), mu_SI, p % vib_solver)
      
      !if (upper(p % vib_solver) .eq. "FOURIER_REAL") then
      !  call solve_fourier_grid_real(dx, E_n(j,:), eig_n(j,:), c_n(j,:,:), mu_SI, "SI", "fast")
      !else 
      !  call solve_sinc_DVR(dx,mu_SI, E_n(j,:), c_n(j,:,:), eig_n(j,:))
      !end if

      !call solve_sinc_DVR(dx,mu_SI, E_n(j,:), c_n(j,:,:), eig_n(j,:))
      write(6,*) "Intermediate state fundamental", j, (eig_n(j,2) -eig_n(j,1))*const % cm
    end do

    ! final states
    do j=1,p % npesfile_f

      call solve_vib_problem(dx, E_f(j,:), eig_f(j,:), c_f(j,:,:), mu_SI, p % vib_solver)
      
      !if (upper(p % vib_solver) .eq. "FOURIER_REAL") then
      !  call solve_fourier_grid_real(dx, E_f(j,:), eig_f(j,:), c_f(j,:,:), mu_SI, "SI", "fast")
      !else
      !  call solve_sinc_DVR(dx,mu_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))
      !end if
      
      !call solve_sinc_DVR(dx,mu_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))
      write(6,*) "Final state fundamental",j, (eig_f(j,2) -eig_f(j,1))*const % cm
    end do

    ! solve non-adiabatic problem
    !call solve_non_adiabatic(eig_f, c_f, eig_na, c_na)

    write(6,*) "here 1"
    call calculate_dipoles_KH_res(c_i, c_n, c_f, dipole_f, D_ni, D_fn, D_fi) 
    write(6,*) "here 2"
        
    ! convert eigenvalues to eV units
    eig_i =eig_i / const % eV
    eig_n =eig_n / const % eV
    eig_f =eig_f / const % eV

    !write(6,*) eig_n(1) -eig_f(1,1)

    ! calculate "mean energy" for intermediate state
    E_n_mean = sum(eig_n(:,1) * D_ni(:, 1, 1) ** 2) 

    write(6,*) "sum (D_ni(:,1,1)**2)", sum( D_ni(:,1,1) ** 2)
    write(6,*) "E_n_mean", E_n_mean

    !
    ! calculate spectrum
    !
      ! use composite indexing for intermediate states
      do n_e = 1, p % npesfile_n    
        do n_v = 1, p % nstates
          n_ev = (n_e -1)* p % nstates + n_v
          eig_n_cmp(n_ev) = eig_n(n_e, n_v) 
          D_ni_cmp(n_ev,:) = D_ni(n_e, n_v, :)
          D_fn_cmp(:, :, n_ev, :) = D_fn(:,:, n_e, n_v, :)
        end do
      end do

      do f_e = 1,p % npesfile_f
        call compute_XES_res(eig_i(1), eig_n_cmp(:), eig_f(f_e,:), &
             D_ni_cmp(:,:), D_fn_cmp(f_e,:,:,:), &
             omega_in, omega, gamma, gamma_inc, gamma_instr, .true., .true., &
             sigma_final(f_e,:,:,:,:), &
             lambda_F(f_e,:,:), lambda_G(f_e,:,:), lambda_H(f_e,:,:))
      end do


      ! averages according to J. Phys. B. 27, 4169 (1994) 
      ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
      lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
      lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
      lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    
    !
    ! Write to files
    ! 

    ! write spectra to file
    do j=1, p % n_omega_in

      file="_sigma_"
      write(string,'(F6.2)') omega_in(j)   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1, p % n_omega 
        write(ifile,'(4ES18.10)') omega(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do

      close(ifile) 

    end do

    !
    ! alternative formula with F(\omega) and instrumental broadening instead of incoming broadening
    ! 

      do f_e = 1,p % npesfile_f
        call compute_XES_res_alt(eig_i(1), eig_n_cmp(:), eig_f(f_e,:), &
             D_ni_cmp(:,:), D_fn_cmp(f_e,:,:,:), &
             omega_in, omega, gamma, gamma_inc, gamma_instr, .true., .true., &
             sigma_final(f_e,:,:,:,:), &
             lambda_F(f_e,:,:), lambda_G(f_e,:,:), lambda_H(f_e,:,:))
      end do


      ! averages according to J. Phys. B. 27, 4169 (1994) 
      ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
      lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
      lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
      lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    
    !
    ! Write to files
    ! 

    ! write spectra to file
    do j=1, p % n_omega_in

      file="_sigma_omega_"
      write(string,'(F6.2)') omega_in(j)   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1, p % n_omega 
        write(ifile,'(4ES18.10)') omega(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do

      close(ifile) 

    end do



    
!    !
!    ! only resonant contribution 
!    !
!    
!    do f_e = 1,p % npesfile_f
!      call compute_XES_res(eig_i(1), eig_n_cmp(:), eig_f(f_e,:), D_ni_cmp(:,:), D_fn_cmp(f_e,:,:,:), &
!           omega_in, omega, gamma, gamma_inc, .true., .false., &
!           sigma_final(f_e,:,:,:,:), &
!           lambda_F(f_e,:,:), lambda_G(f_e,:,:), lambda_H(f_e,:,:))
!    end do
!
!
!      ! averages according to J. Phys. B. 27, 4169 (1994) 
!      ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
!      lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
!      lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
!      lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)
!
!    
!    !
!    ! Write to files
!    ! 
!
!    ! write spectra to file
!    do j=1, p % n_omega_in
!
!      file="_sigma_res_"
!      write(string,'(F6.2)') omega_in(j)   
!      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!
!      ifile = get_free_handle()
!      open(ifile,file=file,status='unknown')
!
!      do i=1, p % n_omega 
!        write(ifile,'(4ES18.10)') omega(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
!      end do
!
!      close(ifile) 
!
!    end do
!
!
!
!    !
!    ! only nonresonant contribution 
!    !
!    
!    do f_e = 1,p % npesfile_f
!      call compute_XES_res(eig_i(1), eig_n_cmp(:), eig_f(f_e,:), D_ni_cmp(:,:), D_fn_cmp(f_e,:,:,:), &
!           omega_in, omega, gamma, gamma_inc, .false., .true., &
!           sigma_final(f_e,:,:,:,:), &
!           lambda_F(f_e,:,:), lambda_G(f_e,:,:), lambda_H(f_e,:,:))
!    end do
!
!
!      ! averages according to J. Phys. B. 27, 4169 (1994) 
!      ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
!      lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
!      lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
!      lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)
!
!    
!    !
!    ! Write to files
!    ! 
!
!    ! write spectra to file
!    do j=1, p % n_omega_in
!
!      file="_sigma_nonres_"
!      write(string,'(F6.2)') omega_in(j)   
!      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!
!      ifile = get_free_handle()
!      open(ifile,file=file,status='unknown')
!
!      do i=1, p % n_omega 
!        write(ifile,'(4ES18.10)') omega(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
!      end do
!
!      close(ifile) 
!
!    end do
!
!    
!
!
!
!    
!    
!    ! write sigma to file
!    file="_lambda_ln"
!    file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
!
!    ifile = get_free_handle()
!    open(ifile,file=file,status='unknown')
!
!    do i=1,p % n_omega
!      write(ifile,*) omega(i), sigma(i)
!    end do
!
!    close(ifile) 

!    do j= 1,p % npesfile_f
!      file="vib_eigenvalues_final_"
!      write(string,*) j
!      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
!
!      open(10,file=file,status='unknown')
!
!      do i=1,p % nstates
!        write(10,*) eig_f(j,i)
!      end do
!
!      write(10,*) 
!      write(10,*) 
!
!      close(10)
!    end do

  end subroutine calculate_KH_res

  subroutine calculate_KH_res_el(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_splines, only: spline_easy
    use m_splines, only: linspace
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_nac_file
    use m_sckh_params_t, only: sckh_params_t 
    use m_KH_functions, only: solve_sinc_DVR
    use m_KH_utils, only: calculate_dipoles_KH_res
    use m_KH_utils, only: compute_XES_res
    use m_KH_utils, only: compute_XES_res_alt
    use m_io, only: get_free_handle

    type(sckh_params_t), intent(inout):: p 

    ! loop variables
    integer::i,j,ii,jj,k,l,m,t 

    ! other variables
    character(80):: file,string
    real(kind=wp):: mu_SI, dx,dvr_start, gamma, gamma_instr, gamma_inc, E_n_mean 
    integer:: npoints 
    real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, E_lp_corr, eig_i, shift
    real(kind=wp), dimension(:),allocatable::  sigma, omega, omega_in
    real(kind=wp), dimension(:,:),allocatable::  c_i, &
         eig_f, E_f, sigma_states, E_n, eig_n
    real(kind=wp), dimension(:,:,:),allocatable:: c_n, c_f, dipole_f, dipole_n, D_fi, D_ni
    !real(kind=wp), dimension(:,:,:,:),allocatable:: D_fn, nac
    real(kind=wp), allocatable::D_fn(:,:,:,:,:)
    real(wp), allocatable:: sigma_final(:,:,:,:,:)
    integer:: ifile
    real(wp):: norm
    real(wp), allocatable:: lambda_F(:,:,:), lambda_G(:,:,:), lambda_H(:,:,:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    real(wp), allocatable:: eig_n_cmp(:), D_ni_cmp(:,:), D_fn_cmp(:,:,:,:)
    integer:: n_e, n_v, n_ev, f_e, i_E_eq, i_E_tmp(1), nstates 

    
    !
    ! This progam calculates the KH emission spectrum using the eigenstate basis
    !

    gamma = p % gamma_FWHM /  2
    gamma_instr = p % gamma_instr_FWHM /  2
    gamma_inc = p % gamma_inc_FWHM /  2 

    nstates = 1
    
    ! allocate everything
    allocate( X_dvr(p % nstates), E_i(p % nstates), E_n(p % npesfile_n, p % nstates), E_lp_corr(p % nstates), &
         E_f(p % npesfile_f,p % nstates), eig_i(nstates),eig_n(p % npesfile_n, nstates),&
         eig_f(p % npesfile_f,nstates), &
         shift(p % nstates))
    allocate( omega(p % n_omega), sigma(p % n_omega), D_ni(p % npesfile_n, nstates, 3) )
    allocate( omega_in(p % n_omega_in))
    allocate(c_i(nstates,nstates),c_f(p % npesfile_f,nstates,nstates), &
         c_n(p % npesfile_n, nstates,nstates), &
         D_fn(p % npesfile_f,nstates,p % npesfile_n, nstates,3))
    allocate(dipole_f(p % npesfile_f, p % nstates, 3), sigma_states(p % npesfile_f,p % n_omega))
    allocate(dipole_n(p % npesfile_n, p % nstates, 3))
    allocate(D_fi(p % npesfile_f,nstates,3) )

    write(6,*) "", p % npesfile_f, p % n_omega_in,  p % n_omega 

    allocate(sigma_final(p % npesfile_f, p % n_omega_in, p % n_omega,3,3))
    allocate(eig_n_cmp(p % npesfile_n * nstates))
    allocate(D_ni_cmp(p % npesfile_n * nstates, 3))
    allocate(D_fn_cmp(p % npesfile_f,nstates, p % npesfile_n * nstates, 3))
    allocate(lambda_F(p % npesfile_f, p % n_omega_in, p % n_omega))
    allocate(lambda_G(p % npesfile_f, p % n_omega_in, p % n_omega))
    allocate(lambda_H(p % npesfile_f, p % n_omega_in, p % n_omega))

    write(6,*) "so far... 1"
    
    !if (p % nonadiabatic .eq. 1) then
    !  allocate(nac(p % npesfile_f, p % npesfile_f, p % nstates, 2) )
    !end if

    !if (mod(p % nstates,2).ne.1 ) then
    !  write(6,*) "nstates must be an odd number"
    !  stop
    !end if

    npoints = (p % nstates-1)/2
    mu_SI = p % mu * const % u
    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10

    write(6,*) "so far... 2"
    
    !
    ! set up DVR points
    !

    do i = -npoints,npoints
      ii = i + npoints +1
      X_dvr(ii) = (ii-1)*dx + dvr_start
    end do

    write(6,*) "so far... 3"
    
    write(6,*) "X_dvr", X_dvr

    !
    ! read PES files
    !

    ! iniital state
    call read_PES_file(p % pes_file_i, p % npoints_in, p % nstates, X_dvr, E_i)

    write(6,*) "so far... 4"
    
    ! intermediate states
    ifile = get_free_handle()
    open(ifile, file= p % pes_file_list_n, action='read')

    allocate(p % pes_files_n(p % npesfile_n))
    do i=1, p % npesfile_n
      read(ifile,*) p % pes_files_n(i)
      write(6,*) p % pes_files_n(i)
    end do
    
    close(ifile)

    ifile = get_free_handle()
    open(ifile, file= p % dipole_file_list_n, action='read')

    write(6,*) "p % dipole_file_list_n  ", p % dipole_file_list_n
    allocate(p % dipolefile_n(p % npesfile_n))
    do i=1, p % npesfile_n
      read(ifile,*) p % dipolefile_n(i)
    end do

    close(ifile)
    
    do j=1,p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, &
           p % nstates, X_dvr, E_n(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % nstates, X_dvr, dipole_n(j,:,:))
    end do

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

    write(6,*) "p % dipole_file_list_f  ", p % dipole_file_list_f
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
    !     call read_nac_file(p % nac_file, p % npoints_in, p %nstates, X_dvr, p % npesfile_f, nac)
    !  end if

    ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved

    if( p % shift_PES .eq.  1) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % nstates, X_dvr, E_lp_corr)

      shift = E_lp_corr -E_f(1,:) 

      do j=1,p % npesfile_f
        E_f(j,:) = E_f(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if


    !create omega_in
    call linspace(omega_in, p % omega_in_start,p % omega_in_end, p % n_omega_in )

    !create omega_out
    call linspace(omega, p % omega_start,p % omega_end, p % n_omega ) 


    write(6,*) "so far... 1"
    
    !
    ! locate minimum of initial state PES and use only that point
    !

    i_E_tmp = minloc(E_i) !? 
    i_E_eq = i_E_tmp(1)

    write(6,*) "i_E_eq", i_E_eq
    
    eig_i(1) = E_i(i_E_eq) / const % eV
    eig_n(:,1) = E_n(:,i_E_eq) / const % eV   
    eig_f(:,1) = E_f(:,i_E_eq) / const % eV   

    !write(6,*) "eig_i", eig_i
    !write(6,*) "eig_n", eig_n
    !write(6,*) "eig_f", eig_f

    write(6,*) "eig_n -eig_i"
    do j=1, nstates 
      write(6,*) "n =", j, eig_n(j,:) - eig_i
    end do
      
    write(6,*) "eig_f -eig_n"
    do j=1, nstates 
      do i=1, nstates
        write(6,*) "n, f", j,i,  eig_n(j,:) - eig_f(i,:)
      end do
    end do

    ! dipoles
    D_ni = 1
    do f_e = 1, p %  npesfile_f
      do n_e = 1, p %  npesfile_f
        D_fn(f_e,1,n_e,1,:) = dipole_f(f_e, i_E_eq, :)
        D_fi(f_e,1,:) = dipole_f(f_e, i_E_eq, :) ! ?
      end do
    end do
    
    !
    ! calculate spectrum
    !
      ! use composite indexing for intermediate states
      do n_e = 1, p % npesfile_n    
        do n_v = 1, nstates !p % nstates
          n_ev = (n_e -1)* p % nstates + n_v
          eig_n_cmp(n_ev) = eig_n(n_e, n_v) 
          D_ni_cmp(n_ev,:) = D_ni(n_e, n_v, :)
          D_fn_cmp(:, :, n_ev, :) = D_fn(:,:, n_e, n_v, :)
        end do
      end do

      do f_e = 1,p % npesfile_f
        call compute_XES_res(eig_i(1), eig_n_cmp(:), eig_f(f_e,:), D_ni_cmp(:,:), D_fn_cmp(f_e,:,:,:), &
             omega_in, omega, gamma, gamma_inc, gamma_instr, .true., .true., &
             sigma_final(f_e,:,:,:,:), &
             lambda_F(f_e,:,:), lambda_G(f_e,:,:), lambda_H(f_e,:,:))
      end do


      ! averages according to J. Phys. B. 27, 4169 (1994) 
      ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
      lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
      lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
      lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    
    !
    ! Write to files
    ! 

    ! write spectra to file
    do j=1, p % n_omega_in

      file="_sigma_"
      write(string,'(F6.2)') omega_in(j)   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1, p % n_omega 
        write(ifile,'(4ES18.10)') omega(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do

      close(ifile) 

    end do

    !
    ! alternative formula with F(\omega) and instrumental broadening instead of incoming broadening
    ! 

      do f_e = 1,p % npesfile_f
        call compute_XES_res_alt(eig_i(1), eig_n_cmp(:), eig_f(f_e,:), D_ni_cmp(:,:), D_fn_cmp(f_e,:,:,:), &
             omega_in, omega, gamma, gamma_inc, gamma_instr, .true., .true., &
             sigma_final(f_e,:,:,:,:), &
             lambda_F(f_e,:,:), lambda_G(f_e,:,:), lambda_H(f_e,:,:))
      end do


      ! averages according to J. Phys. B. 27, 4169 (1994) 
      ! lambda_lp: parallel linear, lambda_ln: perpendicular linear, lambda_cp: circularly polarized
      lambda_lp = sum(2.0_wp * lambda_F + 2.0_wp * lambda_G  + 2.0_wp * lambda_H, 1)
      lambda_ln = sum(-1.0_wp * lambda_F + 4.0_wp  * lambda_G -1.0_wp * lambda_H, 1)
      lambda_cp = sum(-2.0_wp * lambda_F + 3.0_wp * lambda_G + 3.0_wp * lambda_H, 1)

    
    !
    ! Write to files
    ! 

    ! write spectra to file
    do j=1, p % n_omega_in

      file="_sigma_omega_"
      write(string,'(F6.2)') omega_in(j)   
      file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i=1, p % n_omega 
        write(ifile,'(4ES18.10)') omega(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do

      close(ifile) 

    end do

  end subroutine calculate_KH_res_el


  


end module m_KH


