program KH
  use parameters
  use KH_functions
  use spline_m
  implicit none

  ! input/output
  character(80)::pes_file_i,pes_file_n, pes_file_lp_corr, pes_file_ch_corr, outfile
  character(80), dimension(:),allocatable:: pes_file_f,dipolefile_f
  integer:: nstates,ntsteps,n_omega, npesfile_f, shift_PES, runmode, chdyn, broad, nproj
  real(kind=wp):: my, dvr_start_in,dx_in,tstep
  real(kind=wp):: omega_start, omega_end, gamma_FWHM 
  real(kind=wp),allocatable:: projvec(:,:)

  ! loop variables
  integer::i,j,ii,jj,k,l,m,t, proj 

  ! other variables
  character(80):: file,string
  real(kind=wp):: my_SI, dx,dvr_start, E_el_i,E_el_f, E_el_n,  gamma,time,alpha, norm, E_n_mean 
  integer:: INFO,npoints, k_eq, npoints_in
  real(kind=wp), dimension(:),allocatable:: X_dvr,E_i, E_n, E_lp_corr, E_ch_corr, eig_i,eig_n,e_ch,shift, E_in
  real(kind=wp), dimension(:),allocatable:: x_in, sigma, omega, D_ni, sigma_scl1, sigma_scl2, sigma_dir, sigma_max_int
  real(kind=wp), dimension(:,:),allocatable:: Mat_tmp, H_kin,E_dvr, c_i, c_n, c_ch &
       ,eig_f, E_f,sigma_states, sigma_scl2_m, F2_dir, sigma_time, F2_max_int, sigma_lp &
       , delta_e_lp, sigma_dir_states, sigma_max_int_states, dipole_in, sigma_proj
  real(kind=wp), dimension(:,:,:),allocatable:: c_f,dipole, D_fi, F2_dir_lp
  real(kind=wp), dimension(:,:,:,:),allocatable:: D_fn, sigma_scl2_m_states
  complex(kind=wp), dimension(:,:),allocatable:: F, sigma_scl1_m
   complex(kind=wp), dimension(:),allocatable:: F_proj

  !
  ! This progam calcualtes the XES cross section using the Kramers-Heisenberg formula (non-resonant transitions) 
  !

  !
  ! read input
  !

  read(5,*) npoints_in
  read(5,*) pes_file_i
  read(5,*) pes_file_n   ! core ionized state
  read(5,*) shift_PES, pes_file_lp_corr        

  read(5,*) npesfile_f
  allocate(pes_file_f(npesfile_f), dipolefile_f(npesfile_f))

  do i=1,npesfile_f
     read(5,*) pes_file_f(i), dipolefile_f(i)
     write(6,*) i, pes_file_f(i), dipolefile_f(i)
  end do

  read(5,*) outfile
  read(5,*) my
  read(5,*) nstates, dvr_start_in, dx_in
  read(5,*) omega_start, omega_end, n_omega
  read(5,*) gamma_FWHM
  read(5,*) nproj
  
  allocate(projvec(nproj,3))
  
  do i=1, nproj
     read(5,*) projvec(i,1),projvec(i,2),projvec(i,3)
     
     !normalize projvec
     projvec(i,:) = projvec(i,:) / norm2(projvec(i,:))
     
     write(6,*) "projvector", i,  projvec(i,:) 
  end do

  
  ! npoints_in: the number of points in the supplied PES files
  ! pes_file_i: PES file for the initial state
  ! pes_file_n: PES file for the intermediate state
  ! shift_PES: =1, then use pes_file_lp_correct to correct the first final state, =0, do nothing
  ! pes_file_lp_correct: PES for first final stae, computed more accurately
  ! npesfile_f: the number of final state PES files
  ! pesfile_f(i): PES files for the final states
  ! dipolefile_f(i): dipole transition matrix elements between intermediate state and final state i
  ! outfile: the base name of the output files
  ! nstates: the number of dvr points, equally the number of vibrational eigenstates
  ! dvr_start_in: x coordinate of the starting point of the dvr points [Angstrom]
  ! dx_in: spacing between dvr points [Angstroms]
  ! omega_start: starting emission frequency [eV]
  ! omega_end: ending emission frequency [eV]
  ! n_omega: number of emission frequencies
  ! gamma_FWHM: lifetime broadening [eV]
  ! nproj: number of projection vectors

  gamma = gamma_FWHM /  2 

  ! allocate everything
  allocate( X_dvr(nstates), E_in(npoints_in),E_i(nstates), E_n(nstates), E_lp_corr(nstates), E_ch_corr(nstates), &
       E_f(npesfile_f,nstates), eig_i(nstates),eig_n(nstates),eig_f(npesfile_f,nstates), e_ch(nstates), &
       x_in(npoints_in), shift(nstates))
  allocate( omega(n_omega), sigma(n_omega), F(n_omega,3), D_ni(nstates) )
  allocate( sigma_scl1(n_omega), sigma_scl1_m(n_omega,3),sigma_scl2(n_omega), sigma_scl2_m(n_omega,3), &
       sigma_scl2_m_states(npesfile_f, ntsteps,n_omega,3), sigma_time(ntsteps,n_omega))
  allocate(c_i(nstates,nstates),c_f(npesfile_f,nstates,nstates),c_n(nstates,nstates), &
       c_ch(nstates,nstates),D_fn(npesfile_f,nstates,nstates,3))
  allocate(dipole(npesfile_f, nstates, 3), sigma_states(npesfile_f,n_omega))
  allocate(F2_dir(n_omega,3), sigma_dir(n_omega), F2_max_int(n_omega,3), sigma_max_int(n_omega), &
       D_fi(npesfile_f,nstates,3) )
  allocate(F2_dir_lp(10,n_omega,3), sigma_lp(10,n_omega), delta_e_lp(10,nstates))
  allocate(sigma_dir_states(npesfile_f,n_omega), sigma_max_int_states(npesfile_f,n_omega))
  allocate(dipole_in(npoints_in,3))
  allocate(F_proj(n_omega), sigma_proj(n_omega, nproj))


  if (mod(nstates,2).ne.1 ) then
     write(6,*) "nstates must be an odd number"
     stop
  end if

  npoints = (nstates-1)/2
  my_SI = my * amu
  dvr_start = dvr_start_in * 1.0d-10
  dx = dx_in * 1.0d-10

  !
  ! set up DVR points
  !

  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx + dvr_start
  end do


  !
  ! read pes from files (allow for PES:s with arbitrary points)
  !

  !initial state
  open(10,file=pes_file_i,status='unknown')

  do i=1,npoints_in
     read(10,*) x_in(i), E_in(i)
  end do

  close(10) 

  ! use SI units
  x_in = x_in*1d-10
  E_in = E_in *hartree
  
  call check_dvr_bounds(X_dvr, x_in)
  call spline_easy(x_in, E_in, npoints_in, X_dvr, E_i, nstates)


  ! intermediate state, koopmans
  open(10,file=pes_file_n,status='unknown')

  do i=1,npoints_in
     read(10,*) x_in(i), E_in(i)
  end do

  close(10) 

  ! use SI units
  x_in = x_in*1d-10
  E_in = E_in *hartree
  
  call check_dvr_bounds(X_dvr, x_in)
  call spline_easy(x_in, E_in, npoints_in, X_dvr, E_n, nstates)

  ! correct core hole lone pair  PES
  if( shift_PES .eq.  1 ) then

     open(10,file=pes_file_lp_corr,status='unknown')

     do i=1,npoints_in
        read(10,*) x_in(i), E_in(i)
     end do

     close(10)

     ! use SI units
     x_in = x_in*1d-10
     E_in = E_in *hartree
     
     call check_dvr_bounds(X_dvr, x_in)
     call spline_easy(x_in, E_in, npoints_in, X_dvr, E_lp_corr, nstates)
     
  end if ! if shift_PES

  ! final states

  do j=1,npesfile_f

     !pes file
     open(10,file=pes_file_f(j),status='unknown')

     do i=1,npoints_in
        read(10,*) x_in(i), E_in(i)
     end do

     close(10) 


     ! use SI units
     x_in = x_in*1d-10
     E_in = E_in *hartree
     
     call check_dvr_bounds(X_dvr, x_in)
     call spline_easy(x_in, E_in, npoints_in, X_dvr, E_f(j,:), nstates)

     !dipole file
     open(10,file=dipolefile_f(j),status='unknown')
     
     do i=1,npoints_in
        read(10,*) x_in(i), dipole_in(i,1),dipole_in(i,2),dipole_in(i,3)
     end do

     close(10) 

     ! use SI units
     x_in = x_in*1d-10
     
     call check_dvr_bounds(X_dvr, x_in)
     call spline_easy(x_in, dipole_in(:,1), npoints_in, X_dvr, dipole(j,:,1), nstates)
     call spline_easy(x_in, dipole_in(:,2), npoints_in, X_dvr, dipole(j,:,2), nstates)
     call spline_easy(x_in, dipole_in(:,3), npoints_in, X_dvr, dipole(j,:,3), nstates)

  end do ! j

  write(6,*) "Done reading"


  ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
  ! and the spacing between the intermediate and final states are preserved

  if( shift_PES .eq.  1) then

     shift = E_lp_corr -E_f(1,:) 

     do j=1,npesfile_f
        E_f(j,:) = E_f(j,:) + shift
     end do

     write(6,*) "Shifted PES:s"

  end if


  !create omega
  call linspace(omega, omega_start,omega_end, n_omega ) 


  !
  ! Solve the vibrational problem for all eigenfunctions
  !


  ! initial state

  call solve_sinc_DVR(dx,my_SI, E_i, c_i, eig_i)

 ! write(6,*) "Calculated initial state eigenfunctions"

  write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*cm

  ! intermediate state

  call solve_sinc_DVR(dx,my_SI, E_n, c_n, eig_n)

 ! write(6,*) "Calculated intermediate state eigenfunctions"

  write(6,*) "Intermediate state fundamental", (eig_n(2) -eig_n(1))*cm

  ! final states

  do j=1,npesfile_f

     call solve_sinc_DVR(dx,my_SI, E_f(j,:), c_f(j,:,:), eig_f(j,:))

     !write(6,*) "Calculated final state", j

     write(6,*) "Final state fundamental",j, (eig_f(j,2) -eig_f(j,1))*cm

  end do

  !
  ! calculate dipole matrix elements between states 
  !

  do i=1,nstates 
     D_ni(i) = sum(c_n(:,i)*c_i(:,1))
  end do

  do i=1,npesfile_f
     do j=1,nstates ! final
        do k=1,nstates ! intermediate
           do l=1,3
              D_fn(i,j,k,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
              !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
              !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))                  ! no dipole moment
           end do
        end do
     end do
  end  do

  ! transitions from ground state directly to final states
  do i=1,npesfile_f
     do j=1,nstates ! final
        do l=1,3
           D_fi(i,j,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        end do
     end do
  end  do

  write(6,*) "Calculated dipole matrix elements"

  ! convert eigenvalues to eV units
  eig_i =eig_i / eV
  eig_n =eig_n / eV
  eig_f =eig_f / eV

  write(6,*) eig_n(1) -eig_f(1,1)

  ! calculate "mean energy" for intermediate state
  E_n_mean = sum(eig_n(:) * D_ni(:) ** 2) 

  write(6,*) "sum (D_ni(:)**2)", sum( D_ni(:) ** 2)
  write(6,*) "E_n_mean", E_n_mean

  !
  ! Full Kramers-Heisenberg    
  !


  write(6,*) "Full Kramers-Heisenberg"

  sigma=0
  sigma_dir=0
  sigma_max_int = 0
  sigma_proj = 0

  do j= 1,npesfile_f ! final el
     sigma_states(j,:) = 0
     sigma_dir_states(j,:) = 0
     sigma_max_int_states(j,:) = 0
     do l=1,nstates ! final vib

        F=0
        F2_dir=0
        F2_max_int =0

        do i =1, n_omega 
           do k=1,nstates ! intermediate vib
              do m=1,3 ! polarization

                 F(i,m) = F(i,m) + D_fn(j,l,k,m) * D_ni(k) / ( omega(i) - &
                      (eig_n(k) - eig_f(j,l)) + dcmplx(0,gamma) )

                 ! direct contribution (no interference)
                 F2_dir(i,m) = F2_dir(i,m) + D_fn(j,l,k,m) ** 2 * D_ni(k) ** 2 / (( omega(i) - &
                      (eig_n(k) - eig_f(j,l)))**2  + gamma**2 ) 

              end do
           end do

           ! maximal interference, direct transition from initial to final states
           do m=1,3 ! polarization 

              F2_max_int(i,m) =  F2_max_int(i,m) + D_fi(j,l,m) ** 2 / (( omega(i) - &
                   (E_n_mean - eig_f(j,l)) )**2  + gamma**2 )

           end do

        end do

        !square F
        sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
             + real(conjg(F(:,3))*F(:,3)) 

        sigma_dir_states(j,:) = sigma_dir_states(j,:) + F2_dir(:,1) +F2_dir(:,2) +F2_dir(:,3)
        sigma_max_int_states(j,:) = sigma_max_int_states(j,:) + F2_max_int(:,1) + F2_max_int(:,2) + F2_max_int(:,3)
        
        ! projections
        do proj =1, nproj
           F_proj = matmul(F, projvec(proj,:))
           sigma_proj(:,proj) =  sigma_proj(:,proj) + real(conjg(F_proj)*F_proj)
        end do
        
     end do !l

     sigma = sigma + sigma_states(j,:)
     sigma_dir = sigma_dir + sigma_dir_states(j,:)
     sigma_max_int = sigma_max_int + sigma_max_int_states(j,:)

  end do ! j

  write(6,*) "Calculated spectrum"

  ! normalize spectrum
  norm=sum(sigma) *(omega(2) -omega(1)) 
  sigma = sigma/norm
  sigma_dir = sigma_dir/norm
  sigma_max_int = sigma_max_int/norm

  sigma_lp = sigma_lp/norm

  do j= 1,npesfile_f
     sigma_states(j,:) = sigma_states(j,:)/norm
     sigma_dir_states(j,:) = sigma_dir_states(j,:)/norm
     sigma_max_int_states(j,:) = sigma_max_int_states(j,:)/norm
  end do

  !
  ! Write to files
  ! 

  ! write sigma to file
  file="_sigma"
  file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
  
  open(10,file=file,status='unknown')
  
  do i=1,n_omega
     write(10,*) omega(i), sigma(i)
  end do
  
  close(10) 


  ! write sigma_dir to file

  file="_sigma_direct"
  file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"

  !write(string,*) i
  !file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

  open(10,file=file,status='unknown')

  do i=1,n_omega
     write(10,*) omega(i), sigma_dir(i)
  end do

  close(10)

  ! write sigma_max_int to file
  file="sigma_max_int"
  file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"

  open(10,file=file,status='unknown')

  do i=1,n_omega
     write(10,*) omega(i), sigma_max_int(i)
  end do

  close(10)

  ! write sigma_states
  do j= 1,npesfile_f
     file="_sigma_states_"
     write(string,*) j
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

     open(10,file=file,status='unknown')

     do i=1,n_omega
        write(10,*) omega(i), sigma_states(j,i)
     end do

     close(10)
  end do

  ! write sigma_dir_states
  do j= 1,npesfile_f
     file="_sigma_dir_states_"
     write(string,*) j
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

     open(10,file=file,status='unknown')

     do i=1,n_omega
        write(10,*) omega(i), sigma_dir_states(j,i)
     end do

     close(10)
  end do

  ! write sigma_max_int_states
  do j= 1,npesfile_f
     file="_sigma_max_int_states_"
     write(string,*) j
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

     open(10,file=file,status='unknown')

     do i=1,n_omega
        write(10,*) omega(i), sigma_max_int_states(j,i)
     end do

     close(10)
  end do

  ! write sigma_proj
  do j= 1, nproj
     file="_sigma_proj_"
     write(string,*) j
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

     open(10,file=file,status='unknown')

     do i=1,n_omega
        write(10,*) omega(i), sigma_proj(i,j)
     end do

     close(10)
  end do

  ! write vibrational eigenvalues
  file="_vib_eigenvalues_initial"
  file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"
  
  open(10,file=file,status='unknown')
  do i=1,nstates
     write(10,*) eig_i(i), 1 ,10
  end do
  close(10)

  file="_vib_eigenvalues_intermediate"
  file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"

  open(10,file=file,status='unknown')
  do i=1,nstates
     write(10,*) eig_n(i), 1 , D_ni(i)
  end do
  close(10)


  ! write transition moments
  file="_transition_moments_i_to_n"
  file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"
  
  open(10,file=file,status='unknown')
  do i=1,nstates
     write(10,*) D_ni(i)
  end do
  close(10)

  ! Transition dipoles between all intermediate and final states
  if ( .true.) then
     do i=1,npesfile_f
        file="_transition_moments_n_to_f_"
        write(string,*) i
        file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

        open(10,file=file,status='unknown')

        do k=1,nstates ! intermediate
           write(10,*) "Intermediate vibrational state", k
           write(10,*) 
           
           do j=1,nstates ! final
              write(10,'(I5,5ES18.10)') j, D_fn(i,j,k,:),  sum(D_fn(i,j,k,:)**2), eig_n(k) - eig_f(i,j)
           end do
        end do

        close(10)

     end  do

  end if ! if ( .false.) then

  ! a lot of output disabled
  if ( .false.) then
     
     do j= 1,npesfile_f
        file="_vib_eigenvalues_final_"
        write(string,*) j
        file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
        
        open(10,file=file,status='unknown')
        
        do i=1,nstates
           write(10,*) eig_f(j,i), 1, 10
        end do
        
        close(10)
     end do
     
     ! write transition energies for the 10 first intermediate states
     do j= 1,10
        file="_delta_e_lp_"
        write(string,*) j
        file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
        
        open(10,file=file,status='unknown')
        
        do i=1,nstates
           write(10,*) delta_e_lp(j,i)
        end do
        close(10)
     end do
     
  end if ! if ( .false.) then
  
  ! write vibrational eigenfucntions, normally disabled
  if(.false.) then 

     file="_vib_eigenfunctions_initial"
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"     

     open(10,file=file,status='unknown')
     do i=1,nstates
        do j=1,nstates
           write(10,*) X_dvr(j), c_i(j,i)
        end do
        write(10,*) 
        write(10,*) 
     end do
     close(10)

     file="_vib_eigenfunctions_intermediate"
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) //  ".dat"
 
    open(10,file=file,status='unknown')
     do i=1,nstates
        do j=1,nstates
           write(10,*) X_dvr(j), c_n(j,i)
        end do
        write(10,*) 
        write(10,*) 
     end do
     close(10)

     ! write 10 first eigenfunctions to seprate files
     do i=1,10
        file="_vib_eigenfunctions_intermediate_"
        write(string,*) i
        file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     
        
        open(10,file=file,status='unknown')

        do j=1,nstates
           write(10,*) X_dvr(j), c_n(j,i)
        end do
        write(10,*)
        write(10,*)
        close(10)
     end do! i

     ! write 10 first eigenfunctions of first final state to separate files
     do i=1,10
        file="_vib_eigenfunctions_final_1_"
        write(string,*) i
        file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

        open(10,file=file,status='unknown')

        do j=1,nstates
           write(10,*) X_dvr(j), c_f(1,j,i)
        end do
        write(10,*)
        write(10,*)
        close(10)
     end do

     do k= 1,npesfile_f
        file="_vib_eigenfunctions_final_"
        write(string,*) k
        file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

        open(10,file=file,status='unknown')

        do i=1,nstates
           do j=1,nstates
              write(10,*) X_dvr(j), c_f(k,j,i)
           end do

           write(10,*) 
           write(10,*) 

        end do
        close(10)
     end do


  end if ! if(1) 

  do j= 1,npesfile_f
     file="vib_eigenvalues_final_"
     write(string,*) j
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"     

     open(10,file=file,status='unknown')

     do i=1,nstates
        write(10,*) eig_f(j,i)
     end do

     write(10,*) 
     write(10,*) 

     close(10)
  end do

contains

function norm2(a)
  real(kind=wp):: norm2
  real(kind=wp), intent(in):: a(:)
  norm2 = sqrt(sum(a**2))

end function norm2

end program KH


