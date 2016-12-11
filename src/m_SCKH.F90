module m_SCKH
  implicit none

contains

  subroutine  calculate_SCKH(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: compute_F_if_omp
    use m_SCKH_utils, only: read_projections ! move this?
    use m_SCKH_utils, only: read_one_sckh_traj
    use m_fftw3, only: get_omega_reordered_fftw
    use m_sckh_params_t, only: sckh_params_t
    use m_io, only: get_free_handle
    use m_splines, only: spline_easy
    use m_splines, only: linspace
    use m_FFT, only: next_power_of_2
    use m_KH_utils, only: compute_XES_nonres_elec 

    type(sckh_params_t), intent(inout):: p 

    character(80):: dummy
    character(80), dimension(:),allocatable:: traj_files
    integer:: ntraj, ntsteps_inp, nfinal, ntsteps , n_omega 
    integer:: nproj 
    real(kind=wp),dimension(:),allocatable::  time_inp2, E_n_inp, time, E_n, traj_weights 
    real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f, E_trans, projvec
    real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn 
    real(kind=wp):: gamma_FWHM 

    character(80)::  string, file
    real(kind=wp),dimension(:),allocatable:: sigma_tot, time_inp, freq, funct_real, funct_imag
    real(kind=wp),dimension(:,:),allocatable:: sigma_f, sigma_proj
    real(kind=wp):: gamma,  time_l, E, delta_t, time_l2, norm, E_nf_mean, D_proj
    integer:: ntrans
    integer, dimension(:), allocatable:: freq_ind
    complex(kind=wp), allocatable::  F_if_omp(:,:,:,:)
    real(kind=wp),allocatable:: sigma_mm(:,:,:,:)
    real(kind=wp),dimension(:),allocatable:: E_IP1s, E_gs_inp 
    real(kind=wp),dimension(:),allocatable:: omega
    real(kind=wp),dimension(:),allocatable:: D_ni(:)
    integer:: ifile
    integer::i,j,k,l,m,m1,m2, traj, ll,t

    !functions
    real(kind=wp):: dnrm2

    ! set some local variables
    ntraj = p % ntraj
    ntsteps_inp = p % ntsteps
    ntsteps = p % ntsteps2
    nfinal = p % nfinal 
    delta_t = p % delta_t

    ! use HWHM internally
    gamma = p % gamma_FWHM / 2 

    ! projections
    call read_projections(p)

    allocate(traj_files(ntraj), traj_weights(ntraj))

    ! read names of trajectory files
    ifile = get_free_handle()
    open(ifile,file= p % traj_file_list,status='unknown')

    do i=1,ntraj
      read(ifile,*) traj_files(i), traj_weights(i)
    end do

    close(ifile)

    allocate(time_inp(ntsteps_inp),time_inp2(ntsteps_inp), E_gs_inp(ntsteps_inp),E_n_inp(ntsteps_inp), &
         E_f_inp(nfinal,ntsteps_inp), D_fn_inp(nfinal,ntsteps_inp,3), time(ntsteps),&
         E_n(ntsteps), E_f(nfinal,ntsteps), D_fn(nfinal,ntsteps,3),&
         E_IP1s(ntsteps_inp), E_trans(nfinal,ntsteps_inp))
    allocate(D_ni(3))

    D_ni = 1.0_wp

    ! create time_inp, in s
    time_inp = 0.0_wp  
    do i=1,ntsteps_inp
      time_inp(i) = (i-1) * delta_t * 1.d-15
    end do

    time_l = (ntsteps_inp-1) * delta_t *1.d-15

    call linspace(time, time_inp(1), time_inp(ntsteps_inp), ntsteps)

    delta_t = time(2)-time(1)
    time_l2 = (ntsteps-1) * delta_t 

    write(6,*) "gamma_FWHM (fwhm of lorentzian broadening)", p % gamma_FWHM
    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
    write(6,*) "ntsteps", ntsteps
    write(6,*) "delta_t", delta_t
    write(6,*) "time_l", time_l
    write(6,*)
    write(6,*) "Fundamental frequency resolution", 2 * const % pi * const % hbar /( time_l * const % eV)
    write(6,*) "new delta t", delta_t
    write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV)

    n_omega = ntsteps !_pad 
    allocate(F_if_omp(nfinal,n_omega,3,3), sigma_f(nfinal,n_omega), sigma_tot(n_omega), &
         sigma_proj(p % nproj,n_omega), sigma_mm(nfinal,n_omega,3,3), omega(n_omega))
    !
    ! Loop over trajectories
    !

    sigma_mm=0.0_wp

    do traj = 1, ntraj

      !
      ! read trajectory file (this version reads all quantities every time step and interpolates the quantities of interest)
      !  

      call read_one_sckh_traj(ntsteps_inp, nfinal, traj_files(traj), time_inp, &
           time_inp2, E_gs_inp, E_n_inp, E_IP1s, E_trans, E_f_inp, D_fn_inp)

      !
      ! Compute all relevant properties
      !
      !if(p % runmode .eq. "static") then
      !
      !  ! static spectrum
      !  call compute_XES_nonres_elec(E_gs_inp(1), E_n_inp(1), E_f_inp(:,1), D_ni, D_fn_inp(:,1,:), omega, gamma, sigma_m)
      !
      !else

      call spline_easy(time_inp, E_n_inp, ntsteps_inp, time, E_n, ntsteps)

      do i=1,nfinal
        call spline_easy( time_inp, E_f_inp(i,:), ntsteps_inp, time, E_f(i,:), ntsteps)  
      end do

      do i=1,nfinal
        do m=1,3
          call spline_easy(time_inp, D_fn_inp(i,:,m), ntsteps_inp, time, D_fn(i,:,m) , ntsteps)  
        end do
      end do

      ! first time, compute the mean transition energy
      if (traj .eq. 1) then
        E_nf_mean = sum(E_n(:)-  E_f(nfinal,:)) / max(1,size(E_n(:))) 
        
        call get_omega_reordered_fftw(time_l2 * const % eV /  const % hbar, omega)
        omega =  omega + E_nf_mean !E_nf_mean - omega !omega - E_nf_mean

      end if

      call compute_F_if_omp(E_n, E_f, E_nf_mean, D_fn, D_ni, time,  F_if_omp, gamma)

      sigma_mm = sigma_mm + abs(F_if_omp)**2
      
      write(6,*) "Computed trajectory", traj

    end do !end traj

    close(ifile)

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
    write(6,*) "Averaged over", ntraj, "trajectories"

    !normalize 
    norm = sum(sigma_tot) *  2.0_wp * const % pi  * const % hbar / (time_l2 * const % eV ) 
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
    
  end subroutine calculate_SCKH

!  subroutine read_projections(p)
!    use m_precision, only: wp
!    use m_io, only: get_free_handle
!    use m_sckh_params_t, only: sckh_params_t
!
!    type(sckh_params_t), intent(inout):: p   
!
!    integer:: ifile, i
!    real(kind=wp):: dnrm2
!
!    if(p % use_proj) then
!      ifile = get_free_handle()
!      open(ifile, file= p % proj_file, action='read')    
!      read(ifile,*) p % nproj
!
!      allocate(p % projvec(p % nproj,3))
!
!      do i=1, p % nproj                                         
!        read(ifile,*) p % projvec(i,1), p % projvec(i,2), p % projvec(i,3) 
!
!        !normalize projvec                               
!        p % projvec(i,:) = p % projvec(i,:) / dnrm2(3,p % projvec(i,:),1)
!        write(6,*) "projvector", i,  p % projvec(i,:)            
!      end do
!
!      close(ifile)
!    else
!      ! the three cartesian directions
!      p % nproj =3 
!      allocate(p % projvec(p % nproj,3))
!
!      p % projvec =0.0_wp
!      p % projvec(1,1) =1.0_wp
!      p % projvec(2,2) =1.0_wp
!      p % projvec(3,3) =1.0_wp
!
!    end if
!
!  end subroutine read_projections
!
!  subroutine read_one_sckh_traj(ntsteps_inp, nfinal, traj_file, time_inp, &
!       time_inp2, E_gs_inp, E_n_inp, E_IP1s, E_trans, E_f_inp, D_fn_inp)
!    use m_precision, only: wp
!    use m_constants, only: const
!    use m_io, only: get_free_handle
!
!    integer, intent(in):: ntsteps_inp, nfinal
!    character(*), intent(in):: traj_file
!    real(kind=wp), intent(in):: time_inp(:)
!    real(kind=wp), intent(out)::  time_inp2(:), E_gs_inp(:),  E_n_inp(:),&
!         E_IP1s(:), E_trans(:,:), E_f_inp(:,:), D_fn_inp(:,:,:)
!
!    integer:: ifile, i, j, ntrans
!    character(80):: dummy
!
!    ifile = get_free_handle()
!    open(ifile,file=traj_file,status='old')  
!
!    do i=1,ntsteps_inp
!
!      read(ifile,*) time_inp2(i)
!      read(ifile,*) E_gs_inp(i)
!      read(ifile,*) E_n_inp(i)
!      read(ifile,*) E_IP1s(i)  
!
!      read(ifile,*) dummy, ntrans
!
!      ! check that
!      if ( ntrans .ne. nfinal ) then
!        write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
!      end if
!
!      do j =1,nfinal
!        read(ifile,*) E_trans(j,i), D_fn_inp(j,i,1), D_fn_inp(j,i,2), D_fn_inp(j,i,3)
!      end do
!
!      !compute E_f_inp
!      E_f_inp(:,i) = E_gs_inp(i) - E_trans(:,i) + E_IP1s(i) * (const % eV / const % Hartree)
!
!      !check that time_inp(i) = time_inp2(i) 
!      if ( abs(time_inp(i) - time_inp2(i)*1.d-15 ) .gt. 1.d-30) then
!        !write(6,*) "Error in time too big", i, abs(time_inp(i) - time_inp2(i)*1.d-15 )
!      end if
!
!    end do !i
!
!    ! convert to eV units
!    E_n_inp = E_n_inp * const % hartree / const % eV
!    do j=1,nfinal 
!      E_f_inp(j,:) = E_f_inp(j,:) * const % hartree / const % eV 
!    end do
!
!  end subroutine read_one_sckh_traj
!


!  subroutine compute_sckh_diagonal_nonresonant(p)
!    use m_precision,only:wp
!    use m_func,only:read_overlap,get_diagonal_hamiltonian,funct_complex
!    use m_splines,only: spline_easy,linspace
!    use m_constants, only: const
!    use m_FFT, only:  next_power_of_2
!    use m_sckh_params_t, only: sckh_params_t
!    use m_io, only: get_free_handle
!    use m_SCKH_utils, only: ODE_solver
!
!    type(sckh_params_t),intent(inout)::p
!
!    ! input/output
!    character(80):: dummy
!    ! overlap_files - array which contains the names of overlap t12 files
!    character(80), dimension(:),allocatable:: traj_files
!    integer:: ntraj, ntsteps_inp, nfinal, ntsteps , n_omega, ntsteps_pad
!    integer:: ntsteps_pad_pow, nfreq_inp,nproj
!    real(kind=wp),dimension(:),allocatable::  time_inp2, E_n_inp, time, E_n, traj_weight
!    real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f, E_trans, projvec
!    real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn
!    real(kind=wp):: freq_min,  freq_max
!    !loop variables
!    integer::i,j,k,l,m,traj
!
!    ! other variables
!    character(80)::  string, file
!    real(kind=wp),dimension(:),allocatable:: sigma_tot, time_inp, freq, funct_real, funct_imag
!    real(kind=wp),dimension(:,:),allocatable:: sigma, sigma_proj
!    real(kind=wp):: gamma, gamma2, time_l, delta_t, time_l2, norm, E_fn_mean, D_proj
!    integer:: freq_min_i, freq_max_i,nfreq, ntrans
!    integer, dimension(:), allocatable:: freq_ind
!    complex(kind=wp), dimension(:),allocatable::  funct
!    complex(kind=wp), dimension(:,:),allocatable::   sigma_tmp
!    complex(kind=wp), dimension(:,:,:),allocatable::  sigma_m
!    real(kind=wp),dimension(:),allocatable:: E_IP1s, E_gs, traj_weights
!    ! H_mat contains the Hamiltonian matrix elements for the different time steps
!    complex(kind=wp),dimension(:,:,:),allocatable::A_mat
!    complex(kind=wp),dimension(:,:,:),allocatable::Final_state_sum
!    integer::ifile
!    !functions
!    real(kind=wp):: dnrm2
!    ! set up local variables
!    ntraj=p%ntraj
!    ntsteps_inp=p%ntsteps
!    ntsteps=p%ntsteps2
!    nfinal=p%nfinal
!    delta_t=p%delta_t
!
!
!    ! use HWHM internally
!    gamma = p % gamma_FWHM / 2
!
!    ! projections
!    call read_projections(p)
!
!    allocate(traj_files(ntraj), traj_weights(ntraj))
!
!    ! read names of trajectory files
!    ifile = get_free_handle()
!    open(ifile,file= p % traj_file_list,status='unknown')
!
!    do i=1,ntraj
!      read(ifile,*) traj_files(i), traj_weights(i)
!    end do
!
!    close(ifile)
!
!    call next_power_of_2(ntsteps, ntsteps_pad, ntsteps_pad_pow)! next_power_of_2
!
!
!    ! allocate const%Everything (changed dimension of E_n_inp etc)
!    allocate(time_inp(ntsteps_inp),time_inp2(ntsteps_inp), E_gs(ntsteps_inp),E_n_inp(ntsteps_inp), &
!         E_f_inp(nfinal,ntsteps_inp), D_fn_inp(nfinal,ntsteps_inp,3), time(ntsteps),&
!         E_n(ntsteps), E_f(nfinal,ntsteps), D_fn(nfinal,ntsteps,3), &
!         funct(ntsteps), funct_real(ntsteps_pad), funct_imag(ntsteps_pad),&
!         E_IP1s(ntsteps_inp), E_trans(nfinal,ntsteps_inp))
!    allocate(A_mat(nfinal,nfinal,ntsteps))
!    allocate(Final_state_sum(nfinal,ntsteps,3))
!
!    ! create time_inp, in s
!    time_inp = 0
!    do i=1,ntsteps_inp
!      time_inp(i) = (i-1) * delta_t * 1.d-15
!    end do
!
!    time_l = (ntsteps_inp-1) * delta_t *1.d-15
!    call linspace(time, time_inp(1), time_inp(ntsteps_inp), ntsteps) ! linspace(x,start,end,npoints)
!    delta_t = time(2)-time(1)
!    time_l2 = (ntsteps_pad-1) * delta_t
!
!
!    ! some output
!    write(6,*) "gamma_FWHM (fwhm of lorentzian broadening)", p % gamma_FWHM
!    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
!    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma2
!    write(6,*) "ntsteps", ntsteps
!    write(6,*) "delta_t", delta_t
!    write(6,*) "time_l", time_l
!    write(6,*) "ntsteps_pad", ntsteps_pad, "= 2 **", ntsteps_pad_pow 
!    write(6,*) "time_l2 (padded time)", time_l2
!    write(6,*)
!    write(6,*) "Fundamental frequency resolution", 2 * const % pi * const % hbar /( time_l * const % eV)
!    write(6,*) "Padded frequency resolution", 2 * const % pi * const % hbar /( time_l2 * const % eV)
!    write(6,*) "new delta t", delta_t
!    write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV)
!
!
!    nfreq=ntsteps_pad
!    n_omega = ntsteps_pad ! n_omeg=n_freq
!    allocate(sigma_m(nfinal,n_omega,3), sigma(nfinal,n_omega), sigma_tot(n_omega), &
!         sigma_proj(p % nproj,n_omega), sigma_tmp(nfinal,n_omega))
!
!
!    sigma=0.0_wp
!    sigma_tot=0.0_wp
!    sigma_proj=0.0_wp
!
!    !
!    ! Loop over trajectories
!    !
!    do traj = 1, ntraj
!      ! read trajectory file (this version reads all quantities const%Every time step and interpolates the quantities of interest)
!
!
!      call read_one_sckh_traj(ntsteps_inp, nfinal, traj_files(traj), time_inp, &
!           time_inp2, E_gs, E_n_inp, E_IP1s, E_trans, E_f_inp, D_fn_inp)
!
!      ! Transform energy back to the atomic unit
!      E_f_inp=E_f_inp/const%Hartree2eV
!      E_n_inp=E_n_inp/const%Hartree2eV
!
!
!      !
!      ! Compute all relevant properties
!      !
!
!
!      call spline_easy(time_inp, E_n_inp, ntsteps_inp, time, E_n, ntsteps)
!
!      do i=1,nfinal
!        call spline_easy( time_inp, E_f_inp(i,:), ntsteps_inp, time, E_f(i,:), ntsteps)
!      end do
!
!
!      do i=1,nfinal
!        do m=1,3
!          call spline_easy(time_inp, D_fn_inp(i,:,m), ntsteps_inp, time, D_fn(i,:,m) , ntsteps)
!        end do
!      end do
!
!      A_mat=0.0
!
!      if (traj .eq. 1) then
!        E_fn_mean = sum(E_f(nfinal,:) - E_n(:)) / max(1,size(E_n(:)))
!      end if
!
!      ! Solve coupling matrix equation
!
!      call  get_diagonal_hamiltonian(E_f_inp,E_n_inp,nfinal,ntsteps_inp,(time_inp(2)-time_inp(1))/const%autime, E_fn_mean)
!      call ODE_solver(A_mat,time/const%autime,nfinal,ntsteps)
!
!      ! Write the solution to the file
!
!
!      write(6,*) "Semi-Classical Kramers-Heisenberg"
!
!      ! Spline A matrix
!      ! compute \sum_f(Aff_\prime* Df\primen
!
!      sigma_m = 0
!
!      !do a FFT
!      do k =1,nfinal! final state
!        do m=1,3 ! polarization
!          Final_state_sum(k,:,m)=A_mat(k,k,:)*cmplx(D_fn(k,:,m)*(-E_f(k,:)+E_n(:)))
!          funct=Final_state_sum(k,:,m)*exp(-gamma*const%Ev*time(:)/const%hbar)
!          funct_real = 0
!          funct_real(1:ntsteps) = real(funct)
!          funct_imag = 0
!          funct_imag(1:ntsteps) = aimag(funct)
!          call SFFTEU( funct_real, funct_imag, ntsteps_pad, ntsteps_pad_pow, -1 ) 
!          sigma_m(k,:,m) = dcmplx(funct_real, funct_imag)
!        end do ! m
!      end do ! k
!      ! square sigma, check this expression
!      sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3))
!      sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)
!      ! compute projections
!      do i=1,p % nproj
!        sigma_tmp = sigma_m(:,:,1) * p % projvec(i,1) + sigma_m(:,:,2) * &
!             p % projvec(i,2) + sigma_m(:,:,3) * p % projvec(i,3)
!        sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
!      end do
!
!      write(6,*) "Computed trajectory", traj
!
!    end do !end traj
!    close(10)
!    write(6,*)
!    write(6,*) "Averaged over", ntraj, "trajectories"
!    !normalize
!
!    norm = sum(sigma_tot) *  2 * const%pi  * const%hbar / (time_l2 * const%eV ) 
!    sigma_tot = sigma_tot / norm
!    sigma = sigma / norm
!    sigma_proj = sigma_proj / norm
!    E_fn_mean=E_fn_mean*const%Hartree2eV
!    ! write sigma to file
!    file="_sigma"
!    file = trim(adjustl(p%outfile)) // trim(adjustl(file)) // ".dat"
!    open(ifile,file=file,status='unknown')
!
!    do i=nfreq/2, 1, -1 
!      write(ifile,*)  -2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean , sigma_tot(nfreq -i +1)
!    end do
!
!    do i=0,nfreq/2
!      write(ifile,*)  2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma_tot(i+1)
!    end do
!
!    close(ifile)
!
!    ! write spectrum from different final states
!    do j=1, nfinal
!
!      file="_sigma_final_"
!      write(string,*) j
!      file = trim(adjustl(p%outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!      open(ifile,file=file,status='unknown')
!
!      do i=nfreq/2, 1, -1
!        write(ifile,*)  -2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma(j,nfreq-i+1)
!      end do
!
!      do i=0, nfreq/2
!        write(ifile,*)  2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma(j,i+1)
!      end do
!
!      close(ifile)
!    end do ! j
!
!    ! write spectrum from projections  
!    do j=1, nproj
!
!      file="_sigma_proj_"
!      write(string,*) j
!      file = trim(adjustl(p%outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!      open(ifile,file=file,status='unknown')
!
!      do i=nfreq/2, 1, -1
!        write(ifile,*)  -2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma_proj(j,nfreq-i+1)
!      end do
!
!      do i=0, nfreq/2
!        write(ifile,*)  2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma_proj(j,i+1)
!      end do
!
!      close(ifile)
!    end do !j
!  end subroutine  compute_sckh_diagonal_nonresonant
!
!
!
!  subroutine compute_sckh_offdiagonal(p)
!
!    use m_precision,only:wp
!    use m_func,only:read_overlap,get_hamiltonian_offdiagonal,funct_complex,&
!         transform_dipole_operator
!    use m_splines,only: spline_easy,linspace
!    use m_constants, only: const
!    use m_FFT, only:  next_power_of_2
!    use m_sckh_params_t, only: sckh_params_t
!    use m_io, only: get_free_handle
!    use m_SCKH_utils, only: ODE_solver
!
!    type(sckh_params_t),intent(inout)::p
!
!    ! input/output
!    character(80):: dummy
!    ! overlap_files - array which contains the names of overlap t12 files
!    character(80), dimension(:),allocatable:: traj_files,overlap_files
!    integer:: ntraj, ntsteps_inp, nfinal, ntsteps , n_omega, ntsteps_pad
!    integer:: ntsteps_pad_pow, nfreq_inp,nproj
!    real(kind=wp),dimension(:),allocatable::  time_inp2, E_n_inp, time, E_n, traj_weight
!    real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f, E_trans, projvec
!    real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn
!    real(kind=wp):: freq_min,  freq_max
!    !loop variables
!    integer::i,j,k,l,m,traj
!
!    ! other variables
!    character(80)::  string, file
!    real(kind=wp),dimension(:),allocatable:: sigma_tot, time_inp, freq, funct_real, funct_imag
!    real(kind=wp),dimension(:,:),allocatable:: sigma, sigma_proj
!    real(kind=wp):: gamma, gamma2, time_l, delta_t, time_l2, norm, E_fn_mean, D_proj
!    integer:: freq_min_i, freq_max_i,nfreq, ntrans
!    integer, dimension(:), allocatable:: freq_ind
!    complex(kind=wp), dimension(:),allocatable::  funct
!    complex(kind=wp), dimension(:,:),allocatable::   sigma_tmp
!    complex(kind=wp), dimension(:,:,:),allocatable::  sigma_m
!    real(kind=wp),dimension(:),allocatable:: E_IP1s, E_gs, traj_weights
!    ! H_mat contains the Hamiltonian matrix elements for the different time steps
!    complex(kind=wp),dimension(:,:,:),allocatable::A_mat
!    complex(kind=wp),dimension(:,:,:),allocatable::Final_state_sum
!    integer::ifile
!    !functions
!    real(kind=wp):: dnrm2
!    ! set up local variables
!    ntraj=p%ntraj
!    ntsteps_inp=p%ntsteps
!    ntsteps=p%ntsteps2
!    nfinal=p%nfinal
!    delta_t=p%delta_t
!
!    ! use HWHM internally
!    gamma = p % gamma_FWHM / 2
!
!    ! projections
!    call read_projections(p)
!
!    allocate(traj_files(ntraj), traj_weights(ntraj),overlap_files(ntraj))
!
!    ! read names of trajectory files
!    ifile = get_free_handle()
!    open(ifile,file= p % traj_file_list,status='unknown')
!
!    ! by convention overlap file name is 'traj_file_name'+'.t12' 
!    do i=1,ntraj
!      read(ifile,*) traj_files(i), traj_weights(i)
!      overlap_files(i)=trim(traj_files(i))//'.t12'    
!    end do
!
!    close(ifile)
!
!    call next_power_of_2(ntsteps, ntsteps_pad, ntsteps_pad_pow)! next_power_of_2
!
!
!    ! allocate const%Everything (changed dimension of E_n_inp etc)
!    allocate(time_inp(ntsteps_inp),time_inp2(ntsteps_inp), E_gs(ntsteps_inp),E_n_inp(ntsteps_inp), &
!         E_f_inp(nfinal,ntsteps_inp), D_fn_inp(nfinal,ntsteps_inp,3), time(ntsteps),&
!         E_n(ntsteps), E_f(nfinal,ntsteps), D_fn(nfinal,ntsteps,3), &
!         funct(ntsteps), funct_real(ntsteps_pad), funct_imag(ntsteps_pad),&
!         E_IP1s(ntsteps_inp), E_trans(nfinal,ntsteps_inp))
!    allocate(A_mat(nfinal,nfinal,ntsteps))
!    allocate(Final_state_sum(nfinal,ntsteps,3))
!
!    ! create time_inp, in s
!    time_inp = 0
!    do i=1,ntsteps_inp
!      time_inp(i) = (i-1) * delta_t * 1.d-15
!    end do
!
!    time_l = (ntsteps_inp-1) * delta_t *1.d-15
!    call linspace(time, time_inp(1), time_inp(ntsteps_inp), ntsteps) ! linspace(x,start,end,npoints)
!    delta_t = time(2)-time(1)
!    time_l2 = (ntsteps_pad-1) * delta_t
!
!
!    ! some output
!    write(6,*) "gamma_FWHM (fwhm of lorentzian broadening)", p % gamma_FWHM
!    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
!    write(6,*) "gamma (hwhm of lorentzian broadening)", gamma2
!    write(6,*) "ntsteps", ntsteps
!    write(6,*) "delta_t", delta_t
!    write(6,*) "time_l", time_l
!    write(6,*) "ntsteps_pad", ntsteps_pad, "= 2 **", ntsteps_pad_pow 
!    write(6,*) "time_l2 (padded time)", time_l2
!    write(6,*)
!    write(6,*) "Fundamental frequency resolution", 2 * const % pi * const % hbar /( time_l * const % eV)
!    write(6,*) "Padded frequency resolution", 2 * const % pi * const % hbar /( time_l2 * const % eV)
!    write(6,*) "new delta t", delta_t
!    write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV)
!
!
!    nfreq=ntsteps_pad
!    n_omega = ntsteps_pad ! n_omeg=n_freq
!    allocate(sigma_m(nfinal,n_omega,3), sigma(nfinal,n_omega), sigma_tot(n_omega), &
!         sigma_proj(p % nproj,n_omega), sigma_tmp(nfinal,n_omega))
!
!
!    sigma=0.0_wp
!    sigma_tot=0.0_wp
!    sigma_proj=0.0_wp
!
!    !
!    ! Loop over trajectories
!    !
!    do traj = 1, ntraj
!      ! read trajectory file (this version reads all quantities const%Every time step and interpolates the quantities of interest)
!
!
!      call read_one_sckh_traj(ntsteps_inp, nfinal, traj_files(traj), time_inp, &
!           time_inp2, E_gs, E_n_inp, E_IP1s, E_trans, E_f_inp, D_fn_inp)
!
!      !  Read ovrlap matrix
!      call read_overlap(nfinal+1,ntsteps,overlap_files(traj))
!
!      ! Transform energy back to the atomic unit
!      E_f_inp=E_f_inp/const%Hartree2eV
!      E_n_inp=E_n_inp/const%Hartree2eV
!
!
!      !
!      ! Compute all relevant properties
!      !
!
!
!      call spline_easy(time_inp, E_n_inp, ntsteps_inp, time, E_n, ntsteps)
!
!      do i=1,nfinal
!        call spline_easy( time_inp, E_f_inp(i,:), ntsteps_inp, time, E_f(i,:), ntsteps)
!      end do
!
!
!      do i=1,nfinal
!        do m=1,3
!          call spline_easy(time_inp, D_fn_inp(i,:,m), ntsteps_inp, time, D_fn(i,:,m) , ntsteps)
!        end do
!      end do
!
!      A_mat=0.0
!
!      if (traj .eq. 1) then
!        E_fn_mean = sum(E_f(nfinal,:) - E_n(:)) / max(1,size(E_n(:)))
!      end if
!
!      ! Solve coupling matrix equation
!      call  get_hamiltonian_offdiagonal(E_f_inp,E_n_inp,nfinal,ntsteps_inp,&
!           (time_inp(2)-time_inp(1))/const%autime,E_fn_mean)
!      call ODE_solver(A_mat,time/const%autime,nfinal,ntsteps)
!
!      ! Write the solution to the file
!
!
!      write(6,*) "Semi-Classical Kramers-Heisenberg"
!
!      ! Spline A matrix
!
!      call transform_dipole_operator(D_fn,E_f,E_n,nfinal,ntsteps) 
!      sigma_m = 0
!
!      !do a FFT
!      do k =1,nfinal! final state
!        do m=1,3 ! polarization
!          ! compute Final_state_sum=\sum{f_\prime}A_{ff_\prime}*D_{f_\primen}
!          do j=1,nfinal
!            Final_state_sum(k,:,m)=Final_state_sum(k,:,m)+A_mat(k,j,:)*cmplx(D_fn(nfinal+1-j,:,m))
!          enddo
!          funct=Final_state_sum(k,:,m)*exp(-gamma*const%Ev*time(:)/const%hbar)
!          funct_real = 0
!          funct_real(1:ntsteps) = real(funct)
!          funct_imag = 0
!          funct_imag(1:ntsteps) = aimag(funct)
!          call SFFTEU( funct_real, funct_imag, ntsteps_pad, ntsteps_pad_pow, -1 ) 
!          sigma_m(k,:,m) = dcmplx(funct_real, funct_imag)
!        end do ! m
!      end do ! k
!      ! square sigma, check this expression
!      sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3))
!      sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)
!      ! compute projections
!      do i=1,p % nproj
!        sigma_tmp = sigma_m(:,:,1) * p % projvec(i,1) + sigma_m(:,:,2) * &
!             p % projvec(i,2) + sigma_m(:,:,3) * p % projvec(i,3)
!        sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
!      end do
!
!      write(6,*) "Computed trajectory", traj
!
!    end do !end traj
!    close(10)
!    write(6,*)
!    write(6,*) "Averaged over", ntraj, "trajectories"
!    !normalize
!
!    norm = sum(sigma_tot) *  2 * const%pi  * const%hbar / (time_l2 * const%eV ) 
!    sigma_tot = sigma_tot / norm
!    sigma = sigma / norm
!    sigma_proj = sigma_proj / norm
!    E_fn_mean=E_fn_mean*const%Hartree2eV
!    ! write sigma to file
!    file="_sigma"
!    file = trim(adjustl(p%outfile)) // trim(adjustl(file)) // ".dat"
!    open(ifile,file=file,status='unknown')
!
!    do i=nfreq/2, 1, -1 
!      write(ifile,*)  -2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean , sigma_tot(nfreq -i +1)
!    end do
!
!    do i=0,nfreq/2
!      write(ifile,*)  2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma_tot(i+1)
!    end do
!
!    close(ifile)
!
!    ! write spectrum from different final states
!    do j=1, nfinal
!
!      file="_sigma_final_"
!      write(string,*) j
!      file = trim(adjustl(p%outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!      open(ifile,file=file,status='unknown')
!
!      do i=nfreq/2, 1, -1
!        write(ifile,*)  -2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma(j,nfreq-i+1)
!      end do
!
!      do i=0, nfreq/2
!        write(ifile,*)  2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma(j,i+1)
!      end do
!
!      close(ifile)
!    end do ! j
!
!    ! write spectrum from projections  
!    do j=1, nproj
!
!      file="_sigma_proj_"
!      write(string,*) j
!      file = trim(adjustl(p%outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
!      open(ifile,file=file,status='unknown')
!
!      do i=nfreq/2, 1, -1
!        write(ifile,*)  -2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma_proj(j,nfreq-i+1)
!      end do
!
!      do i=0, nfreq/2
!        write(ifile,*)  2 * const%pi * (i) * const%hbar / (time_l2 * const%eV) - E_fn_mean, sigma_proj(j,i+1)
!      end do
!
!      close(ifile)
!    end do !j
!  end subroutine compute_sckh_offdiagonal


end module m_SCKH
