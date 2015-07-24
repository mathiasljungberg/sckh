module m_XAS_nonadiabatic
  use parameters
  use KH_functions
  use spline_m

  implicit none

  type input_params
     ! input/output
     character(80):: runmode
     character(80)::pes_file_i,pes_file_n, pes_file_lp_corr, &
          pes_file_ch_corr, outfile, nac_file
     character(80), dimension(:),allocatable:: pes_file_f,dipolefile_f
     integer:: nstates, n_omega, npesfile_f, shift_PES, nonadiabatic, npoints_in
     real(kind=wp):: my, dvr_start_in, dx_in
     real(kind=wp):: omega_start, omega_end, gamma_FWHM 
  end type input_params

contains


subroutine calculate_XAS_nonadiabatic(inp)
  type(input_params), intent(in):: inp

  ! loop variables
  integer::i,j,ii,jj,k,l,m

  ! other variables
  character(80):: file,string
  real(kind=wp):: my_SI, dx,dvr_start, gamma, E_n_mean 
  integer:: npoints 
  real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, eig_i
  real(kind=wp), dimension(:),allocatable::  sigma, omega
  real(kind=wp), dimension(:,:),allocatable::  c_i, &
       eig_f, E_f, sigma_states
  real(kind=wp), dimension(:,:,:),allocatable:: c_f, dipole, D_fi
  real(kind=wp), dimension(:,:,:,:),allocatable::  nac
  !complex(kind=wp), dimension(:,:),allocatable:: 

  !
  ! This progam calculates the XAS cross section including non-adiabatic couplings
  !

  gamma = inp % gamma_FWHM /  2 

  ! allocate everything
  allocate( X_dvr(inp % nstates), E_i(inp % nstates), &
       E_f(inp % npesfile_f,inp % nstates), eig_i(inp % nstates), eig_f(inp % npesfile_f,inp % nstates))
  allocate( omega(inp % n_omega), sigma(inp % n_omega) )
  allocate(c_i(inp % nstates,inp % nstates),c_f(inp % npesfile_f,inp % nstates,inp % nstates) )
  allocate(dipole(inp % npesfile_f, inp % nstates, 3), sigma_states(inp % npesfile_f,inp % n_omega))
  allocate(D_fi(inp % npesfile_f,inp % nstates,3) )
  
  if (inp % nonadiabatic .eq. 1) then
     allocate(nac(inp % npesfile_f, inp % npesfile_f, inp % nstates, 2) )
  end if

  if (mod(inp % nstates,2).ne.1 ) then
     write(6,*) "nstates must be an odd number"
     stop
  end if

  npoints = (inp % nstates-1)/2
  my_SI = inp % my * amu
  dvr_start = inp % dvr_start_in * 1.0d-10
  dx = inp % dx_in * 1.0d-10

  !
  ! set up DVR points
  !

  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx + dvr_start
  end do

  ! read PES files
  call read_PES_file(inp % pes_file_i, inp % npoints_in, inp % nstates, X_dvr, E_i)
  
  do j=1,inp % npesfile_f
     call read_PES_file(inp % pes_file_f(j), inp % npoints_in, inp % nstates, X_dvr, E_f(j,:))
     call read_dipole_file(inp % dipolefile_f(j), inp % npoints_in, inp % nstates, X_dvr, dipole(j,:,:))
  end do

  if (inp % nonadiabatic .eq. 1) then
     call read_nac_file(inp % nac_file, inp % npoints_in, inp %nstates, X_dvr, inp % npesfile_f, nac)
  end if

  !create omega
  call linspace(omega, inp % omega_start,inp % omega_end, inp % n_omega ) 

  !
  ! Solve the vibrational problem for all eigenfunctions
  !


  ! initial state
  call solve_sinc_DVR(dx,my_SI, E_i, c_i, eig_i)
  write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*cm

  ! final states
  do j=1,inp % npesfile_f
     call solve_sinc_DVR(dx,my_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))
     write(6,*) "Final state fundamental",j, (eig_f(j,2) -eig_f(j,1))*cm
  end do

  ! solve non-adiabatic problem
  !call solve_non_adiabatic(eig_f, c_f, eig_na, c_na)

  ! calculate transition dipoles
  call calculate_dipoles_XAS( c_i, c_f, dipole, D_fi)

  ! convert eigenvalues to eV units
  eig_i =eig_i / eV
  eig_f =eig_f / eV

  write(6,*) eig_i(1) -eig_f(1,1)

  ! calculate spectrum
  if (inp % nonadiabatic .eq. 1) then
     !call spectrum_XES_nonadiabatic(eig_na, eig_n, c_na, D_ni, D_fn, D_fi, omega, sigma, sigma_states)
     stop
  else
     call spectrum_XAS(eig_f, eig_i, D_fi, omega, sigma, sigma_states, gamma)
  end if

  !
  ! Write to files
  ! 

  ! write sigma to file
  file="_sigma"
  file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // ".dat"
  
  open(10,file=file,status='unknown')
  
  do i=1,inp % n_omega
     write(10,*) omega(i), sigma(i)
  end do
  
  close(10) 

  do j= 1,inp % npesfile_f
     file="vib_eigenvalues_final_"
     write(string,*) j
     file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

     open(10,file=file,status='unknown')

     do i=1,inp % nstates
        write(10,*) eig_f(j,i)
     end do

     write(10,*) 
     write(10,*) 

     close(10)
  end do
 
end subroutine calculate_XAS_nonadiabatic

subroutine calculate_XES_nonadiabatic(inp)
  type(input_params), intent(in):: inp 

  ! loop variables
  integer::i,j,ii,jj,k,l,m,t 

  ! other variables
  character(80):: file,string
  real(kind=wp):: my_SI, dx,dvr_start, gamma, E_n_mean 
  integer:: npoints 
  real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, E_n, E_lp_corr, eig_i, eig_n, shift
  real(kind=wp), dimension(:),allocatable::  sigma, omega, D_ni
  real(kind=wp), dimension(:,:),allocatable::  c_i, c_n, &
       eig_f, E_f, sigma_states
  real(kind=wp), dimension(:,:,:),allocatable:: c_f, dipole, D_fi
  real(kind=wp), dimension(:,:,:,:),allocatable:: D_fn, nac
  !complex(kind=wp), dimension(:,:),allocatable:: 

  !
  ! This progam calculates the XES cross section including non-adiabatic couplings
  !

  gamma = inp % gamma_FWHM /  2 

  ! allocate everything
  allocate( X_dvr(inp % nstates), E_i(inp % nstates), E_n(inp % nstates), E_lp_corr(inp % nstates), &
       E_f(inp % npesfile_f,inp % nstates), eig_i(inp % nstates),eig_n(inp % nstates),eig_f(inp % npesfile_f,inp % nstates), &
       shift(inp % nstates))
  allocate( omega(inp % n_omega), sigma(inp % n_omega), D_ni(inp % nstates) )
  allocate(c_i(inp % nstates,inp % nstates),c_f(inp % npesfile_f,inp % nstates,inp % nstates),c_n(inp % nstates,inp % nstates), &
       D_fn(inp % npesfile_f,inp % nstates,inp % nstates,3))
  allocate(dipole(inp % npesfile_f, inp % nstates, 3), sigma_states(inp % npesfile_f,inp % n_omega))
  allocate(D_fi(inp % npesfile_f,inp % nstates,3) )
  
  if (inp % nonadiabatic .eq. 1) then
     allocate(nac(inp % npesfile_f, inp % npesfile_f, inp % nstates, 2) )
  end if

  if (mod(inp % nstates,2).ne.1 ) then
     write(6,*) "nstates must be an odd number"
     stop
  end if

  npoints = (inp % nstates-1)/2
  my_SI = inp % my * amu
  dvr_start = inp % dvr_start_in * 1.0d-10
  dx = inp % dx_in * 1.0d-10

  !
  ! set up DVR points
  !

  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx + dvr_start
  end do

  ! read PES files
  call read_PES_file(inp % pes_file_i, inp % npoints_in, inp % nstates, X_dvr, E_i)
  call read_PES_file(inp % pes_file_n, inp % npoints_in, inp % nstates, X_dvr, E_n)

  do j=1,inp % npesfile_f
     call read_PES_file(inp % pes_file_f(j), inp % npoints_in, inp % nstates, X_dvr, E_f(j,:))
     call read_dipole_file(inp % dipolefile_f(j), inp % npoints_in, inp % nstates, X_dvr, dipole(j,:,:))
  end do

  if (inp % nonadiabatic .eq. 1) then
     call read_nac_file(inp % nac_file, inp % npoints_in, inp %nstates, X_dvr, inp % npesfile_f, nac)
  end if

  ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
  ! and the spacing between the intermediate and final states are preserved

  if( inp % shift_PES .eq.  1) then
     call read_PES_file(inp % pes_file_lp_corr, inp % npoints_in, inp % nstates, X_dvr, E_lp_corr)

     shift = E_lp_corr -E_f(1,:) 

     do j=1,inp % npesfile_f
        E_f(j,:) = E_f(j,:) + shift
     end do
     write(6,*) "Shifted PES:s"
  end if

  !create omega
  call linspace(omega, inp % omega_start,inp % omega_end, inp % n_omega ) 

  !
  ! Solve the vibrational problem for all eigenfunctions
  !


  ! initial state
  call solve_sinc_DVR(dx,my_SI, E_i, c_i, eig_i)
  write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*cm

  ! intermediate state
  call solve_sinc_DVR(dx,my_SI, E_n, c_n, eig_n)
  write(6,*) "Intermediate state fundamental", (eig_n(2) -eig_n(1))*cm

  ! final states
  do j=1,inp % npesfile_f
     call solve_sinc_DVR(dx,my_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))
     write(6,*) "Final state fundamental",j, (eig_f(j,2) -eig_f(j,1))*cm
  end do

  ! solve non-adiabatic problem
  !call solve_non_adiabatic(eig_f, c_f, eig_na, c_na)

  ! calculate transition dipoles
  call calculate_dipoles_XES( c_i, c_n, c_f, dipole, D_ni, D_fn, D_fi)

  ! convert eigenvalues to eV units
  eig_i =eig_i / eV
  eig_n =eig_n / eV
  eig_f =eig_f / eV

  write(6,*) eig_n(1) -eig_f(1,1)

  ! calculate "mean energy" for intermediate state
  E_n_mean = sum(eig_n(:) * D_ni(:) ** 2) 

  write(6,*) "sum (D_ni(:)**2)", sum( D_ni(:) ** 2)
  write(6,*) "E_n_mean", E_n_mean

  ! calculate spectrum
  if (inp % nonadiabatic .eq. 1) then
     !call spectrum_XES_nonadiabatic(eig_na, eig_n, c_na, D_ni, D_fn, D_fi, omega, sigma, sigma_states)
     stop
  else
     call spectrum_XES(eig_f, eig_n, D_ni, D_fn, D_fi, omega, sigma, sigma_states, gamma)
  end if

  !
  ! Write to files
  ! 

  ! write sigma to file
  file="_sigma"
  file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // ".dat"
  
  open(10,file=file,status='unknown')
  
  do i=1,inp % n_omega
     write(10,*) omega(i), sigma(i)
  end do
  
  close(10) 


  !! write sigma_dir to file
  !
  !file="_sigma_direct"
  !file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"
  !
  !!write(string,*) i
  !!file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
  !
  !open(10,file=file,status='unknown')
  !
  !do i=1,inp % n_omega
  !   write(10,*) omega(i), sigma_dir(i)
  !end do
  !
  !close(10)
  !
  !! write sigma_max_int to file
  !file="sigma_max_int"
  !file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"
  !
  !open(10,file=file,status='unknown')
  !
  !do i=1,inp % n_omega
  !   write(10,*) omega(i), sigma_max_int(i)
  !end do
  !
  !close(10)
  !
  !! write sigma_states
  !do j= 1,inp % npesfile_f
  !   file="_sigma_states_"
  !   write(string,*) j
  !   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
  !
  !   open(10,file=file,status='unknown')
  !
  !   do i=1,inp % n_omega
  !      write(10,*) omega(i), sigma_states(j,i)
  !   end do
  !
  !   close(10)
  !end do
  !
  !! write sigma_dir_states
  !do j= 1,inp % npesfile_f
  !   file="_sigma_dir_states_"
  !   write(string,*) j
  !   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
  !
  !   open(10,file=file,status='unknown')
  !
  !   do i=1,inp % n_omega
  !      write(10,*) omega(i), sigma_dir_states(j,i)
  !   end do
  !
  !   close(10)
  !end do
  !
  !! write sigma_max_int_states
  !do j= 1,inp % npesfile_f
  !   file="_sigma_max_int_states_"
  !   write(string,*) j
  !   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
  !
  !   open(10,file=file,status='unknown')
  !
  !   do i=1,inp % n_omega
  !      write(10,*) omega(i), sigma_max_int_states(j,i)
  !   end do
  !
  !   close(10)
  !end do
  !
  !! write vibrational eigenvalues
  !file="_vib_eigenvalues_initial"
  !file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"
  !
  !open(10,file=file,status='unknown')
  !do i=1,inp % nstates
  !   write(10,*) eig_i(i), 1 ,10
  !end do
  !close(10)
  !
  !file="_vib_eigenvalues_intermediate"
  !file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"
  !
  !open(10,file=file,status='unknown')
  !do i=1,inp % nstates
  !   write(10,*) eig_n(i), 1 , D_ni(i)
  !end do
  !close(10)


  !! a lot of output disabled
  !if (0) then
  !   
  !   do j= 1,inp % npesfile_f
  !      file="_vib_eigenvalues_final_"
  !      write(string,*) j
  !      file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
  !      
  !      open(10,file=file,status='unknown')
  !      
  !      do i=1,inp % nstates
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
  !   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"
  !   
  !   open(10,file=file,status='unknown')
  !   do i=1,inp % nstates
  !      write(10,*) D_ni(i)
  !   end do
  !   close(10)
  !   
  !   !file="transition_moments_n_to_f"
  !   !do i=1,inp % npesfile_f
  !   !   do j=1,inp % nstates ! final
  !   !      do k=1,inp % nstates ! intermediate
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
  !      file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
  !      
  !      open(10,file=file,status='unknown')
  !      
  !      do i=1,inp % nstates
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
  !   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"     
  !
  !   open(10,file=file,status='unknown')
  !   do i=1,inp % nstates
  !      do j=1,inp % nstates
  !         write(10,*) X_dvr(j), c_i(j,i)
  !      end do
  !      write(10,*) 
  !      write(10,*) 
  !   end do
  !   close(10)
  !
  !   file="_vib_eigenfunctions_intermediate"
  !   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) //  ".dat"
  !
  !  open(10,file=file,status='unknown')
  !   do i=1,inp % nstates
  !      do j=1,inp % nstates
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
  !      file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
  !      
  !      open(10,file=file,status='unknown')
  !
  !      do j=1,inp % nstates
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
  !      file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
  !
  !      open(10,file=file,status='unknown')
  !
  !      do j=1,inp % nstates
  !         write(10,*) X_dvr(j), c_f(1,j,i)
  !      end do
  !      write(10,*)
  !      write(10,*)
  !      close(10)
  !   end do
  !
  !   do k= 1,inp % npesfile_f
  !      file="_vib_eigenfunctions_final_"
  !      write(string,*) k
  !      file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
  !
  !      open(10,file=file,status='unknown')
  !
  !      do i=1,inp % nstates
  !         do j=1,inp % nstates
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

  do j= 1,inp % npesfile_f
     file="vib_eigenvalues_final_"
     write(string,*) j
     file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

     open(10,file=file,status='unknown')

     do i=1,inp % nstates
        write(10,*) eig_f(j,i)
     end do

     write(10,*) 
     write(10,*) 

     close(10)
  end do

end subroutine calculate_XES_nonadiabatic


subroutine read_input(inp)
  type(input_params), intent(out):: inp 

  integer:: i
  
  read(5,*) inp % runmode      ! either XAS or XAS

  if (inp % runmode .ne. "XAS" .and. inp % runmode .ne. "XES") then
     write(6,*) "runmode must be either XAS or XES"
     stop
  end if

  read(5,*) inp % npoints_in  
  read(5,*) inp % pes_file_i
  read(5,*) inp % pes_file_n   ! core ionized state, not used for XAS
  read(5,*) inp % shift_PES, inp % pes_file_lp_corr        

  read(5,*) inp % npesfile_f
  allocate(inp % pes_file_f(inp % npesfile_f), inp % dipolefile_f(inp % npesfile_f))

  do i=1,inp % npesfile_f
     read(5,*) inp % pes_file_f(i), inp % dipolefile_f(i)
  end do

  read(5,*) inp % nonadiabatic, inp % nac_file
  read(5,*) inp % outfile
  read(5,*) inp % my
  read(5,*) inp % nstates, inp % dvr_start_in, inp % dx_in
  read(5,*) inp % omega_start, inp % omega_end, inp % n_omega
  read(5,*) inp % gamma_FWHM
  
  ! runmode: either "XAS" or "XES"
  ! npoints_in: the number of points in the supplied PES files
  ! pes_file_i: PES file for the initial state
  ! pes_file_n: PES file for the intermediate state
  ! shift_PES: =1, then use pes_file_lp_correct to correct the first final state, =0, do nothing
  ! pes_file_lp_corr: PES for first final state, computed more accurately
  ! npesfile_f: the number of final state PES files
  ! pesfile_f(i): PES files for the final states
  ! dipolefile_f(i): dipole transition matrix elements between intermediate state and final state i
  ! nonadiabatic: if 1 yes, if 0 no
  ! nac_file: the file containing the non-adiabatic couplings
  ! outfile: the base name of the output files
  ! nstates: the number of dvr points, equally the number of vibrational eigenstates
  ! dvr_start_in: x coordinate of the starting point of the dvr points [Angstrom]
  ! dx_in: spacing between dvr points [Angstroms]
  ! omega_start: starting emission frequency [eV]
  ! omega_end: ending emission frequency [eV]
  ! n_omega: number of emission frequencies
  ! gamma_FWHM: lifetime broadening [eV]

end subroutine read_input

subroutine read_PES_file(filename, npoints_in, nstates, X_dvr, E)
  character(80), intent(in):: filename
  integer, intent(in):: npoints_in, nstates
  real(kind=wp), intent(inout):: X_dvr(nstates)
  real(kind=wp), intent(out):: E(nstates)

  real(kind=wp):: x_in(npoints_in), E_in(npoints_in)
  integer:: i

  !
  ! read pes from files (allow for PES:s with arbitrary points)
  !

  open(10,file= filename,status='unknown')

  do i=1,npoints_in
     read(10,*) x_in(i), E_in(i)
  end do

  close(10) 

  ! use SI units
  x_in = x_in*1d-10
  E_in = E_in *hartree
  
  call check_dvr_bounds(X_dvr, x_in)
  call spline_easy(x_in, E_in, npoints_in, X_dvr, E, nstates)

end subroutine read_PES_file


subroutine read_dipole_file(filename, npoints_in, nstates, X_dvr, dipole)
  character(80), intent(in):: filename
  integer, intent(in):: npoints_in, nstates
  real(kind=wp), intent(inout):: X_dvr(nstates)
  real(kind=wp), intent(out):: dipole(nstates,3)
  
  real(kind=wp):: x_in(npoints_in), dipole_in(npoints_in,3)
  integer:: i 

  open(10,file=filename ,status='unknown')
  
  do i=1,npoints_in
     read(10,*) x_in(i), dipole_in(i,1),dipole_in(i,2),dipole_in(i,3)
  end do
  
  close(10) 

  ! use SI units
  x_in = x_in*1d-10
     
  call check_dvr_bounds(X_dvr, x_in)
  call spline_easy(x_in, dipole_in(:,1), npoints_in, X_dvr, dipole(:,1), nstates)
  call spline_easy(x_in, dipole_in(:,2), npoints_in, X_dvr, dipole(:,2), nstates)
  call spline_easy(x_in, dipole_in(:,3), npoints_in, X_dvr, dipole(:,3), nstates)

end subroutine read_dipole_file


subroutine read_nac_file(filename, npoints_in, nstates, X_dvr, npesfiles, nac)
  character(80), intent(in):: filename
  integer, intent(in):: npoints_in, nstates
  real(kind=wp), intent(inout):: X_dvr(nstates)
  integer, intent(in):: npesfiles
  real(kind=wp), intent(out):: nac(npesfiles,npesfiles,nstates, 2)  

  real(kind=wp):: nac_in(npesfiles,npesfiles,npoints_in, 2), x_in(npoints_in)  
  integer:: j,k,l 

  open(10,file=filename ,status='unknown')

  do j=1,npesfiles
     do k=1,npesfiles
        do l=1,npoints_in
           ! read < j(R) | d/dR | k(R) > and < j(R) | d2/dR2 | k(R) >    
           read(10,*) x_in(l), nac_in(j,k,l,1), nac_in(j,k,l,2)
        end do

        call spline_easy(x_in, nac_in(j,k,:,1), npoints_in, X_dvr, nac(j,k,:,1), nstates)
        call spline_easy(x_in, nac_in(j,k,:,2), npoints_in, X_dvr, nac(j,k,:,2), nstates)
     end do
  end do

end subroutine read_nac_file

subroutine solve_non_adiabatic(eig_f, c_f, nac, eig_na, c_na) )
  real(kind=wp), intent(in):: eig_f(:,:), c_f(:,:,:), nac(:,:,:)
  real(kind=wp), intent(in):: eig_na(:), c_na(:,:,:)

  !! local variables
  integer:: i,j,ii,jj,npoints,nstates
  real(kind=wp), dimension(:,:), allocatable::Mat_tmp, H_kin
  real(kind=wp), dimension(:),allocatable:: W
  real(kind=wp), dimension(:),allocatable:: WORK
  integer:: INFO, LWORK

  npesfiles = size(eig_f,1)
  nstates = size(eig_f,2)
  nstates_na = size(eig_na,1)

  !write(6,*) "solve_sinc_dvr: nstates=", nstates
 
  LWORK = 3*nstates_na
  
  !
  allocate(Mat_tmp(nstates_na,nstates_na), H_kin(nstates_na,nstates_na), W(nstates_na),&
       WORK(LWORK))
  
  ! set up hamiltonan matrix
 H_tot = 0.0_wp
  do i1=1, npesfiles_f
     do j1=1, nstates
        
        ii = (i1 -1) * nstates + j1  

        do i2=1, npesfiles_f
           do j2=1, nstates

              jj = (i2 -1) * nstates + j2  

              Mat_tmp(ii,jj) = eig_f(i1, j1) * delta(i1,i2) * delta(j1,j2) + &
                   sum(c_f(i1, :, j1) * c_f(i2, : , j2) * nac(i1,i2,:,2) )

              ! here should be a term + 2 * <g_G|<G| d/dR |F> d/dR |f_F> (cannot yet compute it) 

           end do ! j2
        end do ! i2
     end do ! j1
  end do ! i1
  
  ! solve eigenvalue problem 
  call dsyev( "V", "U", nstates_na, Mat_tmp, nstates, W, WORK, LWORK, INFO)
  
  ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
  ! put c_na in a good order

  do i1=1, npesfiles_f
     do j1=1, nstates
        ii = (i1 -1) * nstates + j1  

        c_na(i1, j1, :) = Mat_tmp(ii, :)
     end do
  end do

  eig_na = W


end subroutine solve_non_adiabatic


subroutine calculate_dipoles_XAS(c_i, c_f, dipole, D_fi)
  real(kind=wp), intent(in):: c_i(:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out)::  D_fi(:,:,:)

  integer:: nstates, npesfile_f,i,j,l

  nstates = size(c_i,2)
  npesfile_f = size(c_f,1)

  !
  ! calculate dipole matrix elements between states 
  !

  ! transitions from ground state to final states
  do i=1,npesfile_f
     do j=1,nstates ! final
        do l=1,3
           D_fi(i,j,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        end do
     end do
  end do
  
  write(6,*) "Calculated dipole matrix elements"

end subroutine calculate_dipoles_XAS


subroutine spectrum_XAS(eig_f, eig_i, D_fi, omega, sigma, sigma_states, gamma)
  real(kind=wp), intent(in):: eig_f(:,:), eig_i(:), D_fi(:,:,:), omega(:)
  real(kind=wp), intent(out):: sigma(:), sigma_states(:,:)
  real(kind=wp),intent(in):: gamma

  integer:: nstates, npesfile_f, n_omega, i,j,k,l,m
  complex(kind=wp),allocatable:: F(:,:)
  real(kind=wp):: norm

  nstates = size(eig_f,2)
  npesfile_f = size(eig_f,1)
  n_omega = size(omega)

  allocate(F(n_omega, 3))

  !
  ! Fermi golden rule 
  !

  sigma = 0.0_wp
  sigma_states = 0.0_wp

  do j= 1,npesfile_f ! F
     do l=1,nstates ! f_F
        
        do i =1, n_omega 
           do m=1,3 ! polarization
              
              F(i,m) = F(i,m) + D_fi(j,l ,m)   / ( omega(i) - &
                   (eig_f(j,l) -eig_i(1)) + dcmplx(0,gamma) )
              
           end do
        end do
        
        !square F
        sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
             + real(conjg(F(:,3))*F(:,3)) 
        
     end do ! l
     
     sigma = sigma  + sigma_states(j,:)
     
  end do !j
  
  write(6,*) "Calculated XAS spectrum"

  ! normalize spectrum
  norm=sum(sigma) *(omega(2) -omega(1)) 
  sigma = sigma/norm
  
  do j= 1, npesfile_f
     sigma_states(j,:) = sigma_states(j,:)/norm
  end do

  deallocate(F)

end subroutine spectrum_XAS

subroutine spectrum_XAS_nonadiabatic(eig_na, eig_i, c_na, D_fi, omega, sigma, sigma_states, gamma)
  real(kind=wp), intent(in):: eig_na(:,:), eig_i(:), c_na(:,:), D_fi(:), omega(:)
  real(kind=wp), intent(out):: sigma(:), sigma_states(:,:)
  real(kind=wp),intent(in):: gamma

  !
  ! Fermi golden rule 
  !

!  sigma=0
!  
!  do j= 1,inp % npesfile_f ! F'
!     do j= 1,inp % npesfile_f ! F
!        !sigma_states(j,:) = 0
!        do l=1,nstates ! f_F
!
!           do i =1, inp % n_omega 
!              do m=1,3 ! polarization
!                 
!                 F(i,m) = F(i,m) +    / ( omega(i) - &
!                      (eig_na(j) -eig_i(1)) + dcmplx(0,gamma) )
!                 
!              end do
!           end do
!        end do
!     end do
!
!     !square F
!     !sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
!     !     + real(conjg(F(:,3))*F(:,3)) 
!     
!     sigma = sigma  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
!          + real(conjg(F(:,3))*F(:,3)) 
!
!  end do !j
!  
!
!
!  write(6,*) "Calculated spectrum"
!
!  ! normalize spectrum
!  norm=sum(sigma) *(omega(2) -omega(1)) 
!  sigma = sigma/norm
!  sigma_dir = sigma_dir/norm
!  sigma_max_int = sigma_max_int/norm
!
!  sigma_lp = sigma_lp/norm
!
!  do j= 1,inp % npesfile_f
!     sigma_states(j,:) = sigma_states(j,:)/norm
!     sigma_dir_states(j,:) = sigma_dir_states(j,:)/norm
!     sigma_max_int_states(j,:) = sigma_max_int_states(j,:)/norm
!  end do


end subroutine spectrum_XAS_nonadiabatic

subroutine calculate_dipoles_XES(c_i, c_n, c_f, dipole, D_ni, D_fn, D_fi)
  real(kind=wp), intent(in):: c_i(:,:), c_n(:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out):: D_ni(:), D_fn(:,:,:,:), D_fi(:,:,:)

  integer:: nstates, npesfile_f,i,j,k,l

  nstates = size(c_i,2)
  npesfile_f = size(c_f,1)

  !
  ! calculate dipole matrix elements between states 
  !

  do i=1,nstates 
     D_ni(i) = sum(c_n(:,i)*c_i(:,1))
  end do
     
  do i=1, npesfile_f
     do j=1, nstates ! final
        do k=1, nstates ! intermediate
           do l=1,3
              D_fn(i,j,k,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
              !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
              !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))                  ! no dipole moment
           end do
        end do
     end do
  end  do

  ! transitions from ground state directly to final states
  do i=1, npesfile_f
     do j=1, nstates ! final
        do l=1,3
           D_fi(i,j,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        end do
     end do
  end do

  write(6,*) "Calculated XES dipole matrix elements"

end subroutine calculate_dipoles_XES



subroutine spectrum_XES(eig_f, eig_n,  D_ni, D_fn, D_fi, omega, sigma, sigma_states, gamma)
  real(kind=wp), intent(in):: eig_f(:,:), eig_n(:), D_ni(:), D_fn(:,:,:,:), D_fi(:,:,:), omega(:)
  real(kind=wp), intent(out):: sigma(:), sigma_states(:,:)
  real(kind=wp), intent(in):: gamma

  complex(kind=wp), dimension(:,:),allocatable:: F
  integer:: npesfile_f, nstates, n_omega
  integer:: i,j,k,l,m
  real(kind=wp):: norm


  npesfile_f = size(eig_f,1)
  nstates = size(eig_f,2)
  n_omega = size(omega)


  allocate(F(n_omega,3))


  !
  ! Full Kramers-Heisenberg    
  !


  write(6,*) "Full Kramers-Heisenberg"

  sigma=0
  !sigma_dir=0
  !sigma_max_int = 0

  do j= 1,npesfile_f ! final el
     sigma_states(j,:) = 0
     !sigma_dir_states(j,:) = 0
     !sigma_max_int_states(j,:) = 0
     do l=1,nstates ! final vib

        F=0
        !F2_dir=0
        !F2_max_int =0

        do i =1,  n_omega 
           do k=1,nstates ! intermediate vib
              do m=1,3 ! polarization

                 F(i,m) = F(i,m) + D_fn(j,l,k,m) * D_ni(k) / ( omega(i) - &
                      (eig_n(k) - eig_f(j,l)) + dcmplx(0,gamma) )

                 ! direct contribution (no interference)
                 !F2_dir(i,m) = F2_dir(i,m) + D_fn(j,l,k,m) ** 2 * D_ni(k) ** 2 / (( omega(i) - &
                 !     (eig_n(k) - eig_f(j,l)))**2  + gamma**2 ) 

              end do
           end do

           ! maximal interference, direct transition from initial to final states
           !do m=1,3 ! polarization 
!
!              F2_max_int(i,m) =  F2_max_int(i,m) + D_fi(j,l,m) ** 2 / (( omega(i) - &
!                   (E_n_mean - eig_f(j,l)) )**2  + gamma**2 )
!           end do

        end do

        !square F
        sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
             + real(conjg(F(:,3))*F(:,3)) 

        !sigma_dir_states(j,:) = sigma_dir_states(j,:) + F2_dir(:,1) +F2_dir(:,2) +F2_dir(:,3)
        !sigma_max_int_states(j,:) = sigma_max_int_states(j,:) + F2_max_int(:,1) + F2_max_int(:,2) + F2_max_int(:,3)

     end do !l

     sigma = sigma + sigma_states(j,:)
     !sigma_dir = sigma_dir + sigma_dir_states(j,:)
     !sigma_max_int = sigma_max_int + sigma_max_int_states(j,:)

  end do ! j

  write(6,*) "Calculated XES spectrum"

  ! normalize spectrum
  norm=sum(sigma) *(omega(2) -omega(1)) 
  sigma = sigma/norm
  !sigma_dir = sigma_dir/norm
  !sigma_max_int = sigma_max_int/norm

  !sigma_lp = sigma_lp/norm

  do j= 1,npesfile_f
     sigma_states(j,:) = sigma_states(j,:)/norm
     !sigma_dir_states(j,:) = sigma_dir_states(j,:)/norm
     !sigma_max_int_states(j,:) = sigma_max_int_states(j,:)/norm
  end do

  deallocate(F)

end subroutine spectrum_XES

!subroutine print_output()
!end subroutine print_output

end module m_XAS_nonadiabatic


