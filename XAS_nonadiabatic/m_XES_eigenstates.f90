module m_XES_eigenstates
  use parameters
  use m_XAS_io
  use KH_functions
  use spline_m
  use m_XAS_functions

  implicit none

contains

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

end module m_XES_eigenstates


