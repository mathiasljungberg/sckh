program SCKH
  use parameters
  implicit none

  ! input/output                                                                                                                                                                                                  
  character(80):: files_file, outfile,dummy
  character(80), dimension(:),allocatable:: trajfile
  integer:: ntraj, ntsteps_inp, nfinal, ntsteps , n_omega, ntsteps_pad
  integer:: ntsteps_pad_pow, runmode, nfreq_inp,nproj 
  real(kind=wp),dimension(:),allocatable::  time_inp2, E_n_inp, time, E_n, traj_weight 
  real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f, E_trans, projvec
  real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn 
  real(kind=wp):: freq_min,  freq_max, gamma_FWHM, gamma_FWHM_2, tstep

  !loop variables
  integer::i,j,k,l,m,traj, ll,t

  ! other variables
  character(80)::  string, file
  real(kind=wp),dimension(:),allocatable:: sigma_tot, time_inp, freq, funct_real, funct_imag
  real(kind=wp),dimension(:,:),allocatable:: int_W_I, sigma, sigma_proj
  real(kind=wp):: gamma, gamma2, time_l, E, delta_t, time_l2, norm, E_fn_mean, D_proj
  integer:: freq_min_i, freq_max_i,nfreq, ntrans
  integer, dimension(:), allocatable:: freq_ind
  complex(kind=wp), dimension(:),allocatable::  funct
  complex(kind=wp), dimension(:,:),allocatable::  e_factor1, sigma_tmp
  complex(kind=wp), dimension(:,:,:),allocatable::  sigma_m
  real(kind=wp),dimension(:),allocatable:: E_IP1s, E_gs 

  !functions
  real(kind=wp):: dnrm2

  ! Read input 
  read(5,*) nfinal
  read(5,*) ntraj, tstep, ntsteps_inp, ntsteps
  read(5,*) files_file
  read(5,*) outfile
  read(5,*) freq_min, freq_max, nfreq_inp
  read(5,*) gamma_FWHM, gamma_FWHM_2     
  read(5,*) runmode
  read(5,*) nproj
 
  allocate(projvec(nproj,3))

  do i=1, nproj
     read(5,*) projvec(i,1),projvec(i,2),projvec(i,3)
     
     !normalize projvec
     projvec(i,:) = projvec(i,:) / dnrm2(3,projvec(i,:),1)

     write(6,*) "projvector", i,  projvec(i,:) 
  end do


  allocate(trajfile(ntraj), traj_weight(ntraj))

  ! read names of trajectory files
  open(10,file=files_file,status='unknown')

  do i=1,ntraj
     read(10,*) trajfile(i), traj_weight(i)
  end do

  close(10)
  

  !
  ! Input parameters
  !
  ! nfinal  :  number of final states 
  ! ntraj : number of trajectories
  ! tstep : time step in the trajectories 
  ! ntsteps_inp : the number of timesteps in the input trajectories 
  ! ntsteps : number of timesteps, interpolated from the ones in trajectory. For runmode=1 this can to be like 500, 
  !           since the number of time steps determines the maximum frequency. for runmode=2 don't use more than makes
  !           a smooth spectrum, often nsteps = nsteps_inp will do
  ! files_file  :  a file containing the names of the trajectroy files and their weights 
  ! outfile : the file where the spectrum is written
  ! freq_min : minimum frequency (only for runmode = 2 ) in eV
  ! freq_max : maximum frequency (only for runmode = 2 ) in eV 
  ! nfreq_inp : number of frequency bins (only for runmode = 2 ) 
  ! gamma_FWHM : FWHM of Loretzian in eV
  ! gamma_FWHM_2 : empirical broadening parameter (only for runmode 2). For rigorous results gamma_FWHM_2 = gamma_FWHM in eV
  ! runmode :  1: Semi-classical Kramers-Heisenberg, 2: method of Odelius / Takahashi
  ! nproj: number of projection vectors

  ! The filenames of the trajectories
  !allocate( trajfile(ntraj) )


  ! pad ntsteps to be  apower of 2 for fourier transform
  do i=1,15
     if( 2 ** i .gt. ntsteps) then
        ntsteps_pad_pow = i        
        exit
     end if
  end do

  ntsteps_pad = 2 ** ntsteps_pad_pow ! 2**12
  
  ! allocate everything (changed dimension of E_n_inp etc) 
  allocate(time_inp(ntsteps_inp),time_inp2(ntsteps_inp), E_gs(ntsteps_inp),E_n_inp(ntsteps_inp), &
       E_f_inp(nfinal,ntsteps_inp), D_fn_inp(nfinal,ntsteps_inp,3), time(ntsteps),&
       E_n(ntsteps), E_f(nfinal,ntsteps), D_fn(nfinal,ntsteps,3), int_W_I(nfinal,ntsteps),&
       e_factor1(nfinal,ntsteps), &
       funct(ntsteps), funct_real(ntsteps_pad), funct_imag(ntsteps_pad),&
       E_IP1s(ntsteps_inp), E_trans(nfinal,ntsteps_inp))
  
  ! use HWHM internally
  gamma = gamma_FWHM / 2 

  ! gamma2 is an empirical broadening parameter different from gamma for the purpose of 
  ! computing spectra in the old (erronous) way
  gamma2 = gamma_FWHM_2 / 2 

  ! create time_inp, in s
  time_inp = 0  
  do i=1,ntsteps_inp
     time_inp(i) = (i-1) * tstep * 1.d-15
  end do

  time_l = (ntsteps_inp-1) * tstep *1.d-15
  
  call linspace(time, time_inp(1), time_inp(ntsteps_inp), ntsteps)

  delta_t = time(2)-time(1)
  time_l2 = (ntsteps_pad-1) * delta_t 

  ! some output
  write(6,*) "gamma_FWHM (fwhm of lorentzian broadening)", gamma_FWHM
  write(6,*) "gamma_FWHM_2 (fwhm of lorentzian broadening)", gamma_FWHM_2
  write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
  write(6,*) "gamma (hwhm of lorentzian broadening)", gamma2
  write(6,*) "ntsteps", ntsteps
  write(6,*) "tstep", tstep
  write(6,*) "time_l", time_l
  write(6,*) "ntsteps_pad", ntsteps_pad, "= 2 **", ntsteps_pad_pow 
  write(6,*) "time_l2 (padded time)", time_l2
  write(6,*)
  write(6,*) "Fundamental frequency resolution", 2 * pi * hbar /( time_l * eV)
  write(6,*) "Padded frequency resolution", 2 * pi * hbar /( time_l2 * eV)
  write(6,*) "new delta t", delta_t
  write(6,*) "max freq",  pi * hbar /( delta_t  * eV)


  if (runmode .eq. 1 ) then
     nfreq = ntsteps_pad 
  end if

if (runmode .eq. 2) then

   nfreq = nfreq_inp
   write(6,*) "nfreq", nfreq
   allocate( freq(nfreq), freq_ind(nfreq))  

   write(6,*) "allocated freq", size(freq),nfreq
   freq = 0
   freq_ind =0 
   
   call linspace(freq, freq_min, freq_max, nfreq)
end if



allocate(sigma_m(nfinal,nfreq,3), sigma(nfinal,nfreq), sigma_tot(nfreq),  sigma_proj(nproj,nfreq), sigma_tmp(nfinal,nfreq))

  !
  ! Loop over trajectories
  !

  sigma=0
  sigma_tot=0
  sigma_proj=0

  do traj = 1, ntraj
     
  !
  ! read trajectory file (this version reads all quantities every time step and interpolates the quantities of interest)
  !  

     open(10,file=trajfile(traj),status='old')  
         
     do i=1,ntsteps_inp
        
        read(10,*) time_inp2(i)
        read(10,*) E_gs(i)
        read(10,*) E_n_inp(i)
        read(10,*) E_IP1s(i)  

        read(10,*) dummy, ntrans

        ! check that
        if ( ntrans .ne. nfinal ) then
           write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
        end if
        
        do j =1,nfinal
           read(10,*) E_trans(j,i), D_fn_inp(j,i,1), D_fn_inp(j,i,2), D_fn_inp(j,i,3)
        end do

        !compute E_f_inp
        E_f_inp(:,i) = E_gs(i) - E_trans(:,i) + E_IP1s(i) * (eV / Hartree)

        !check that time_inp(i) = time_inp2(i) 
        if ( abs(time_inp(i) - time_inp2(i)*1.d-15 ) .gt. 1.d-30) then
           !write(6,*) "Error in time too big", i, abs(time_inp(i) - time_inp2(i)*1.d-15 )
         end if

     end do !i

     ! convert to eV units
     E_n_inp = E_n_inp * hartree / eV
     do j=1,nfinal 
        E_f_inp(j,:) = E_f_inp(j,:) * hartree / eV 
     end do

  !
  ! Compute all relevant properties
  !

  call spline_easy(time_inp, E_n_inp, ntsteps_inp, time, E_n, ntsteps)

  do i=1,nfinal
     !call spline_easy( x, E_f_inp(i,:), npoints_pes, dist_OH, E_f(i,:), ntsteps)  
     call spline_easy( time_inp, E_f_inp(i,:), ntsteps_inp, time, E_f(i,:), ntsteps)  
  end do
  
  do i=1,nfinal
     do m=1,3
        call spline_easy(time_inp, D_fn_inp(i,:,m), ntsteps_inp, time, D_fn(i,:,m) , ntsteps)  
     end do
  end do


  ! now everything should be splined and ready. 

  if (runmode .eq. 1) then
     write(6,*) "Semi-Classical Kramers-Heisenberg" 
    
  ! compute  e^{-i \int_0^t W_I(t') dt' }
  
     ! first time, compute the mean transition energy
     if (traj .eq. 1) then
        E_fn_mean = sum(E_f(nfinal,:) - E_n(:)) / max(1,size(E_n(:))) 
     end if
     
  int_W_I(:,1) = E_f(:,1) - E_n(1)  
  
  do i = 2, ntsteps
     int_W_I(:,i) = int_W_I(:,i-1) + ( E_f(:,i) - E_n(i)) - E_fn_mean  
  end do
   int_W_I =  int_W_I * delta_t

   write(6,*) "E_fn_mean", E_fn_mean

  do i = 1,nfinal 
     e_factor1(i,:) = exp(dcmplx(0, -(eV  / hbar) *int_W_I(i,:)  ))
  end do

  ! compute A_{fm}(t) = D^m_{fn}(t) * e_factor_f(t) * exp(-gamma * eV * time(:) / hbar)
  ! and fourer transform
  
  sigma_m = 0
  
  !do a FFT
  do k =1,nfinal! final state 
     do m=1,3 ! polarization     

        !funct = D_fn(k,:,m) * e_factor1(k,:) * exp(dcmplx(0, - E_fn_mean(k) * eV / hbar * time(:))) * &
        !            exp(-time(:) / ( 2 * lifetime * 1.d-15))

        funct = D_fn(k,:,m) * e_factor1(k,:) * exp(-gamma * eV * time(:) / hbar)
        
        funct_real = 0
        funct_real(1:ntsteps) = real(funct)
        funct_imag = 0
        funct_imag(1:ntsteps) = aimag(funct)

        call SFFTEU( funct_real, funct_imag, ntsteps_pad, ntsteps_pad_pow, 1 ) ! 4096 = 2**11

        sigma_m(k,:,m) = dcmplx(funct_real, funct_imag)
        
     end do ! m
  end do ! k

 ! square sigma, check this expression
 sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3)) 
 sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)

 ! compute projections
 do i=1,nproj
    sigma_tmp = sigma_m(:,:,1) * projvec(i,1) + sigma_m(:,:,2) * projvec(i,2) + sigma_m(:,:,3) * projvec(i,3)
    sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
 end do
 
  
else  if(runmode .eq. 2) then
   write(6,*) "Method of Odelius / Takashashi et al. "
   
       sigma_m = 0

       ! compute sigma_tot  
       do k =1,nfinal! final state 
          do m=1,3 ! polarization     
             do l = 1, nfreq
                do t = 1, ntsteps              
                   sigma_m(k,l,m) =  sigma_m(k,l,m) + D_fn(k,t,m) ** 2 *  exp(-(gamma * eV / hbar) * time(t)) &
                        / (  (freq(l) - (E_n(t) - E_f(k,t)) ) ** 2 + gamma2 ** 2) 

                end do !t
             end do ! l
          end do ! m
       end do ! k

       ! compute projections
       do i=1,nproj
         sigma_tmp = sigma_m(:,:,1) * projvec(i,1) + sigma_m(:,:,2) * projvec(i,2) + sigma_m(:,:,3) * projvec(i,3)
         sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
       end do

!       do k =1,nfinal! final state 
!          do m=1,nproj
!             do l = 1, nfreq
!                do t = 1, ntsteps              
!                   D_proj = D_fn(k,t,1) * projvec(m,1) + D_fn(k,t,2) * projvec(m,2) + D_fn(k,t,3) * projvec(m,3)
!                   sigma_proj(k,l,m) =  sigma_proj(k,l,m) + D_proj ** 2 * exp(-(gamma * eV / hbar) * time(t)) &
!                        / (  (freq(l) - (E_n(t) - E_f(k,t)) ) ** 2 + gamma2 ** 2)  !dcmplx(funct_real, funct_imag)
!                end do !t
!             end do ! l
!          end do ! m
!       end do ! k

    !sigma = sigma + real( sum( sigma_m, 3)) 
       !    sigma_tot = sigma_tot +  sum( real( sum( sigma_m , 3)),1)

       sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3)) 
       sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)

    
    !! compute projections (sigma_m already squared!)
    !do i=1,nproj
    !   sigma_proj(i,:) = sigma_proj(i,:) + real( sum( sigma_m(:,:,i) ,1 ) )
    !end do

end if

write(6,*) "Computed trajectory", traj
  

end do !end traj



close(10)

write(6,*)
write(6,*) "Averaged over", ntraj, "trajectories"


if (runmode .eq. 1) then 

!normalize 
norm = sum(sigma_tot) *  2 * pi  * hbar / (time_l2 * eV ) 
sigma_tot = sigma_tot / norm
sigma = sigma / norm
sigma_proj = sigma_proj / norm

! write sigma to file
file="_sigma"
file = trim(adjustl(outfile)) // trim(adjustl(file)) // ".dat"
open(10,file=file,status='unknown')
 
do i=nfreq/2, 1, -1 
   write(10,*)  -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean , sigma_tot(nfreq -i +1)
end do

do i=0,nfreq/2
   write(10,*)  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma_tot(i+1)
end do

close(10)

! write spectrum from different final states
do j=1, nfinal

   file="_sigma_final_"
   write(string,*) j
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')
   
   do i=nfreq/2, 1, -1
      write(10,*)  -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma(j,nfreq-i+1)
   end do
 
   do i=0, nfreq/2
      write(10,*)  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma(j,i+1)
   end do
     
   close(10)
end do ! j

! write spectrum from projections  
do j=1, nproj
  
   file="_sigma_proj_"
   write(string,*) j
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do i=nfreq/2, 1, -1
      write(10,*)  -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma_proj(j,nfreq-i+1)
   end do
     
   do i=0, nfreq/2
      write(10,*)  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma_proj(j,i+1)
   end do
     
   close(10)
end do !j
  

else if (runmode .eq. 2) then
   
!normalize 
norm = sum(sigma_tot) * (freq(2)-freq(1))

sigma_tot = sigma_tot / norm
sigma = sigma / norm
sigma_proj = sigma_proj / norm


! write sigma to file
file="_sigma"
file = trim(adjustl(outfile)) // trim(adjustl(file)) // ".dat"
open(10,file=file,status='unknown')

do i=1,nfreq
   write(10,*)  freq(i), sigma_tot(i)
end do

close(10)


! write spectrum from different final states
do j=1, nfinal

   file="_sigma_final_"
   write(string,*) j
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do i=1,nfreq
      write(10,*)  freq(i), sigma(j,i)
   end do

   close(10)
end do

! write spectrum from projections
do j=1, nproj

   file="_sigma_proj_"
   write(string,*) j
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do i=1,nfreq
      write(10,*)  freq(i), sigma_proj(j,i)
   end do
end do

end if
         
        
end program SCKH
