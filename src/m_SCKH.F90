module m_SCKH
  !use parameters
  implicit none

contains

subroutine  calculate_SCKH(p)
  use m_precision, only: wp
  use m_constants, only: const
  use m_SCKH_utils, only: compute_SCKH
  use m_sckh_params_t, only: sckh_params_t
  use m_io, only: get_free_handle
  use m_splines, only: spline_easy
  use m_splines, only: linspace
  use m_FFT, only: next_power_of_2

  type(sckh_params_t), intent(inout):: p 

  ! input/output                                                                                                                                                                                                  
  character(80):: dummy
  character(80), dimension(:),allocatable:: traj_files
  integer:: ntraj, ntsteps_inp, nfinal, ntsteps , n_omega, ntsteps_pad
  integer:: ntsteps_pad_pow, nproj 
  real(kind=wp),dimension(:),allocatable::  time_inp2, E_n_inp, time, E_n, traj_weights 
  real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f, E_trans, projvec
  real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn 
  real(kind=wp):: gamma_FWHM 

  !loop variables
  integer::i,j,k,l,m,traj, ll,t

  ! other variables
  character(80)::  string, file
  real(kind=wp),dimension(:),allocatable:: sigma_tot, time_inp, freq, funct_real, funct_imag
  real(kind=wp),dimension(:,:),allocatable:: int_W_I, sigma, sigma_proj
  real(kind=wp):: gamma, gamma2, time_l, E, delta_t, time_l2, norm, E_fn_mean, D_proj
  integer:: ntrans
  integer, dimension(:), allocatable:: freq_ind
  complex(kind=wp), dimension(:,:),allocatable::  sigma_tmp
  complex(kind=wp), dimension(:,:,:),allocatable::  sigma_m
  real(kind=wp),dimension(:),allocatable:: E_IP1s, E_gs 
  real(kind=wp),dimension(:),allocatable:: omega
  integer:: ifile

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
  if(p % use_proj) then
    ifile = get_free_handle()
    open(ifile, file= p % proj_file, action='read')    
    read(ifile,*) p % nproj
    
    allocate(p % projvec(p % nproj,3))
    
    do i=1, p % nproj                                         
      read(ifile,*) p % projvec(i,1), p % projvec(i,2), p % projvec(i,3) 
      
      !normalize projvec                               
      p % projvec(i,:) = p % projvec(i,:) / dnrm2(3,p % projvec(i,:),1)
      write(6,*) "projvector", i,  p % projvec(i,:)            
    end do

    close(ifile)
  else
    ! the three cartesian directions
    p % nproj =3 
    allocate(p % projvec(p % nproj,3))
    
    p % projvec =0.0_wp
    p % projvec(1,1) =1.0_wp
    p % projvec(2,2) =1.0_wp
    p % projvec(3,3) =1.0_wp

  end if

  ! 

  allocate(traj_files(ntraj), traj_weights(ntraj))

  ! read names of trajectory files
  ifile = get_free_handle()
  open(ifile,file= p % traj_file_list,status='unknown')

  do i=1,ntraj
     read(ifile,*) traj_files(i), traj_weights(i)
  end do

  close(ifile)
  
  call next_power_of_2(ntsteps, ntsteps_pad, ntsteps_pad_pow)

  ! allocate everything (changed dimension of E_n_inp etc) 
  allocate(time_inp(ntsteps_inp),time_inp2(ntsteps_inp), E_gs(ntsteps_inp),E_n_inp(ntsteps_inp), &
       E_f_inp(nfinal,ntsteps_inp), D_fn_inp(nfinal,ntsteps_inp,3), time(ntsteps),&
       E_n(ntsteps), E_f(nfinal,ntsteps), D_fn(nfinal,ntsteps,3), int_W_I(nfinal,ntsteps),&
       E_IP1s(ntsteps_inp), E_trans(nfinal,ntsteps_inp))
  
  ! create time_inp, in s
  time_inp = 0.0_wp  
  do i=1,ntsteps_inp
     time_inp(i) = (i-1) * delta_t * 1.d-15
  end do

  time_l = (ntsteps_inp-1) * delta_t *1.d-15
  
  call linspace(time, time_inp(1), time_inp(ntsteps_inp), ntsteps)

  delta_t = time(2)-time(1)
  time_l2 = (ntsteps_pad-1) * delta_t 

  ! some output
  write(6,*) "gamma_FWHM (fwhm of lorentzian broadening)", p % gamma_FWHM
  write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
  write(6,*) "gamma (hwhm of lorentzian broadening)", gamma2
  write(6,*) "ntsteps", ntsteps
  write(6,*) "delta_t", delta_t
  write(6,*) "time_l", time_l
  write(6,*) "ntsteps_pad", ntsteps_pad, "= 2 **", ntsteps_pad_pow 
  write(6,*) "time_l2 (padded time)", time_l2
  write(6,*)
  write(6,*) "Fundamental frequency resolution", 2 * const % pi * const % hbar /( time_l * const % eV)
  write(6,*) "Padded frequency resolution", 2 * const % pi * const % hbar /( time_l2 * const % eV)
  write(6,*) "new delta t", delta_t
  write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV)

  n_omega = ntsteps_pad 
  allocate(sigma_m(nfinal,n_omega,3), sigma(nfinal,n_omega), sigma_tot(n_omega), &
       sigma_proj(p % nproj,n_omega), sigma_tmp(nfinal,n_omega), omega(n_omega))

  !
  ! Loop over trajectories
  !

  sigma=0.0_wp
  sigma_tot=0.0_wp
  sigma_proj=0.0_wp

  do traj = 1, ntraj
     
  !
  ! read trajectory file (this version reads all quantities every time step and interpolates the quantities of interest)
  !  

    ifile = get_free_handle()
     open(ifile,file=traj_files(traj),status='old')  
         
     do i=1,ntsteps_inp
        
        read(ifile,*) time_inp2(i)
        read(ifile,*) E_gs(i)
        read(ifile,*) E_n_inp(i)
        read(ifile,*) E_IP1s(i)  

        read(ifile,*) dummy, ntrans

        ! check that
        if ( ntrans .ne. nfinal ) then
           write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
        end if
        
        do j =1,nfinal
           read(ifile,*) E_trans(j,i), D_fn_inp(j,i,1), D_fn_inp(j,i,2), D_fn_inp(j,i,3)
        end do

        !compute E_f_inp
        E_f_inp(:,i) = E_gs(i) - E_trans(:,i) + E_IP1s(i) * (const % eV / const % Hartree)

        !check that time_inp(i) = time_inp2(i) 
        if ( abs(time_inp(i) - time_inp2(i)*1.d-15 ) .gt. 1.d-30) then
           !write(6,*) "Error in time too big", i, abs(time_inp(i) - time_inp2(i)*1.d-15 )
         end if

     end do !i

     ! convert to eV units
     E_n_inp = E_n_inp * const % hartree / const % eV
     do j=1,nfinal 
        E_f_inp(j,:) = E_f_inp(j,:) * const % hartree / const % eV 
     end do

  !
  ! Compute all relevant properties
  !

  call spline_easy(time_inp, E_n_inp, ntsteps_inp, time, E_n, ntsteps)

  do i=1,nfinal
     call spline_easy( time_inp, E_f_inp(i,:), ntsteps_inp, time, E_f(i,:), ntsteps)  
  end do
  
  do i=1,nfinal
     do m=1,3
        call spline_easy(time_inp, D_fn_inp(i,:,m), ntsteps_inp, time, D_fn(i,:,m) , ntsteps)  
     end do
  end do


  ! now everything should be splined and ready. 

  ! first time, compute the mean transition energy
  if (traj .eq. 1) then
    E_fn_mean = sum(E_f(nfinal,:) - E_n(:)) / max(1,size(E_n(:))) 

    ! omega centered around E_fn_mean
    j=1
    do i=n_omega/2, 1, -1
      omega(j) = -2.0_wp * const % pi * (i) * const % hbar / (time_l2 * const % eV) - E_fn_mean
      j=j+1
    end do
    do i=0, n_omega/2-1
      omega(j) =  2.0_wp * const % pi * (i) * const % hbar / (time_l2 * const % eV) - E_fn_mean
      j=j+1
    end do
    
  end if
  
  call compute_SCKH(E_n, E_f, E_fn_mean, D_fn, time,  sigma_m, gamma)
   
  ! square sigma, check this expression
  sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3)) 
  sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)
  
  ! compute projections
  do i=1,p % nproj
    sigma_tmp = sigma_m(:,:,1) * p % projvec(i,1) + sigma_m(:,:,2) * &
         p % projvec(i,2) + sigma_m(:,:,3) * p % projvec(i,3)
    sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
  end do
  
  write(6,*) "Computed trajectory", traj
  
end do !end traj

close(ifile)

write(6,*)
write(6,*) "Averaged over", ntraj, "trajectories"

!normalize 
norm = sum(sigma_tot) *  2.0_wp * const % pi  * const % hbar / (time_l2 * const % eV ) 
sigma_tot = sigma_tot / norm
sigma = sigma / norm
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
     write(ifile,*)  omega(i), sigma(j,i)
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


end module m_SCKH
