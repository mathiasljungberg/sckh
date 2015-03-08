module m_test_XAS_nonadiabatic
  use parameters
  use m_XAS_io
  use KH_functions
  use spline_m

  implicit none

contains

subroutine test_XAS_nonadiabatic(inp)
  type(input_params), intent(in):: inp

  ! run a few unit tests

  !call test_solve_nonadiabatic

end subroutine test_XAS_nonadiabatic


subroutine test_solve_nonadiabatic(inp)
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

  


  stop

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
 
end subroutine test_solve_nonadiabatic

end module m_test_XAS_nonadiabatic


