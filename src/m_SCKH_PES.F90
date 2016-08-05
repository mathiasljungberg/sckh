module m_SCKH_PES

  implicit none

contains

subroutine calculate_SCKH_PES(p)
  use m_precision, only: wp
  use m_constants, only: const
  use m_KH_functions, only: solve_sinc_DVR
  use m_SCKH_utils, only: sample_x_mom
  use m_SCKH_utils, only: compute_SCKH
  use m_SCKH_utils, only: verlet_trajectory
  use m_sckh_params_t, only: sckh_params_t 
  use hist_class, only: hist, hist_init, hist_add
  use hist_class, only: hist_broadening, hist_write
  use m_io, only: get_free_handle
  use m_splines, only: spline_easy
  use m_FFT, only: next_power_of_2
  use m_PES_io, only: read_dipole_file
  use m_PES_io, only: read_nac_file
  use m_PES_io, only: read_PES_file

  type(sckh_params_t), intent(inout):: p 

  ! input/output
  integer:: nfinal, ntsteps,  ntsteps_pad
  integer:: ntsteps_pad_pow, npoints_in  
  real(kind=wp),dimension(:),allocatable::  E_n_inp, time, E_n 
  real(kind=wp),dimension(:,:),allocatable:: E_f_inp, E_f, E_trans
  real(kind=wp),dimension(:,:,:),allocatable:: D_fn_inp, D_fn 

  !loop variables
  integer::i,j,m,traj

  ! other variables
  character(80)::  string, file
  real(kind=wp),dimension(:),allocatable:: sigma_tot
  real(kind=wp),dimension(:),allocatable:: x_sampl, mom_sampl, x_new, omega
  real(kind=wp),dimension(:,:),allocatable:: int_W_I, sigma, sigma_proj, c_i, x_mom_sampl
  real(kind=wp):: gamma, time_l, delta_t, time_l2, norm, E_fn_mean, mu_SI, dx
  integer:: n_omega, npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl
  complex(kind=wp), dimension(:,:),allocatable::   sigma_tmp
  complex(kind=wp), dimension(:,:,:),allocatable::  sigma_m
  real(kind=wp),dimension(:),allocatable:: eig_i, E_IP1s, E_i_inp, E_lp_corr, shift 
  type(hist), dimension(:), allocatable:: time_h
  type(hist):: time_h_0, time_h_0_mom
  integer, dimension(1)::ind  
  integer:: ifile

  real(kind=wp), dimension(:,:,:,:),allocatable::  nac
  real(kind=wp), dimension(:),allocatable::  X_dvr
  real(kind=wp) :: dvr_start 
  integer:: npoints, ii
  real(8):: dnrm2

  ! set some local variables
  ntsteps = p % ntsteps
  npoints_x_sampl = p % npoints_x_sampl
  npoints_mom_sampl =  p % npoints_mom_sampl
  npoints_in = p % npoints_in
  nfinal = p % npesfile_f
  mu_SI = p % mu * const % u

  npoints = (p % npoints_in -1)/2
  dvr_start = p % dvr_start_in * 1.0d-10
  dx = p % dx_in * 1.0d-10
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

  !outfile = p % outfile

  call next_power_of_2(ntsteps, ntsteps_pad, ntsteps_pad_pow)

  write(6,*) "ntsteps_pad_pow", ntsteps_pad, ntsteps_pad_pow

  if(p % samplemode .eq. 1) then
     npoints_x_mom_sampl = npoints_x_sampl * npoints_mom_sampl
  else if(p % samplemode .eq. 2) then
     npoints_x_mom_sampl = npoints_x_sampl !* 2 !npoints_x_sampl * npoints_mom_sampl
     npoints_mom_sampl =  npoints_x_sampl
  else
     write(6,*) "samplemode muste be either 1 or 2!"
     stop
  end if

  ! allocate everything (changed dimension of E_n_inp etc) 
  allocate( E_i_inp(npoints_in),E_n_inp(npoints_in), &
       E_f_inp(nfinal,npoints_in), D_fn_inp(nfinal,npoints_in,3), time(ntsteps),&
       E_n(ntsteps), E_f(nfinal,ntsteps), D_fn(nfinal,ntsteps,3), int_W_I(nfinal,ntsteps),&
       E_IP1s(npoints_in), E_trans(nfinal,npoints_in),E_lp_corr(npoints_in),&
       shift(npoints_in), &
       c_i(npoints_in,npoints_in), eig_i(npoints_in), &
       x_sampl(npoints_x_sampl), mom_sampl(npoints_mom_sampl), x_new(ntsteps),&
       time_h(ntsteps), &
       x_mom_sampl(npoints_x_mom_sampl, 2))

  allocate(X_dvr(npoints_in))
  !allocate(E_f_elastic(1,ntsteps)
  !allocate(D_fn_elastic(1,ntsteps,3)
  
  if (p % nonadiabatic .eq. 1) then
     allocate(nac(p % npesfile_f, p % npesfile_f, p % npoints_in, 2) )
  end if
  
  ! set up DVR points
  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx + dvr_start
  end do

  ! read PES files
  call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_dvr, E_i_inp)
  call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_dvr, E_n_inp)

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
     call read_PES_file(p % pes_files_f(j), p % npoints_in, p % npoints_in, X_dvr, E_f_inp(j,:))
     call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % npoints_in, X_dvr, D_fn_inp(j,:,:))
  end do

  if (p % nonadiabatic .eq. 1) then
     call read_nac_file(p % nac_file, p % npoints_in, p %npoints_in, X_dvr, p % npesfile_f, nac)
  end if

  ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
  ! and the spacing between the intermediate and final states are preserved

  if( p % shift_PES .eq.  1) then
     call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % npoints_in, X_dvr, E_lp_corr)

     shift = E_lp_corr -E_f_inp(1,:) 

     do j=1,p % npesfile_f
        E_f_inp(j,:) = E_f_inp(j,:) + shift
     end do
     write(6,*) "Shifted PES:s"
  end if


  ! convert to eV units
  E_i_inp = E_i_inp  / const % eV
  E_n_inp = E_n_inp  / const % eV
     
  do j=1,nfinal
     E_f_inp(j,:) = E_f_inp(j,:) / const % eV
  end do
  
  !write(6,*) E_n_inp(30) -E_f_inp(1,30)


  !
  ! Solve the vibrational problem for initial state to be able to sample initial distribution
  !  

  call solve_sinc_DVR(dx, mu_SI, E_i_inp * const % eV, c_i, eig_i)
  write(6,*) "Calculated initial state eigenfunctions"
  write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm


  if(p % samplemode .eq. 1) then  
     call sample_x_mom(X_dvr, c_i(:,1), x_sampl, mom_sampl, x_mom_sampl, 1)
  else if(p % samplemode .eq. 2) then  
     call sample_x_mom(X_dvr, c_i(:,1), x_sampl, mom_sampl, x_mom_sampl, 2)
  end if

  write(6,*) "Sampling done"

  !check
  ifile = get_free_handle()
  open(ifile, file="sampled_points.txt", action='write')
  do i=1,npoints_x_sampl
     write(ifile,*) x_sampl(i)
  end do
  close(ifile)

  ifile = get_free_handle()
  open(ifile, file="x_mom_sampl.txt", action='write')
  do i=1,npoints_x_mom_sampl
     write(22,*) x_mom_sampl(i,1), x_mom_sampl(i,2)
  end do
  close(ifile)

  ifile = get_free_handle()
  open(ifile, file="inital_state_eigvec.txt", action='write')
  do i=1,npoints_in
     write(ifile,'(3ES16.6)') X_dvr(i), c_i(i,1), c_i(i,1) ** 2
  end do
  close(ifile)

  delta_t = delta_t * 1.d-15 ! femtoseconds
  time_l = (ntsteps-1) * delta_t
  time_l2 = (ntsteps_pad-1) * delta_t 

  write(6,*) "outfile", p % outfile
  write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
  write(6,*) "mu_SI", mu_SI 
  write(6,*) "time_l", time_l
  write(6,*) "ntsteps_pad", ntsteps_pad, "= 2 **", ntsteps_pad_pow 
  write(6,*) "time_l2 (padded time)", time_l2
  write(6,*)
  write(6,*) "Fundamental frequency resolution", 2.0_wp * const % pi * const % hbar /( time_l * const % eV)
  write(6,*) "Padded frequency resolution", 2.0_wp * const % pi * const % hbar /( time_l2 * const % eV)
  write(6,*) "delta t", delta_t
  write(6,*) "max freq",  const % pi * const % hbar /( delta_t  * const % eV), "eV"

  n_omega = ntsteps_pad 
  allocate(sigma_m(nfinal,n_omega,3), sigma(nfinal,n_omega), sigma_tot(n_omega), &
       sigma_proj(p % nproj,n_omega), sigma_tmp(nfinal,n_omega), omega(n_omega))

  !
  ! Loop over trajectories
  !

  
  do i=1, ntsteps
     time(i)= (i-1)*delta_t
  end do
           
  call hist_init(time_h_0, 1000, X_dvr(1), X_dvr(npoints_in) ) 
  call hist_init(time_h_0_mom, 1000, minval(x_mom_sampl(:,2)), maxval(x_mom_sampl(:,2)) ) 
  do i=1, ntsteps
     call hist_init(time_h(i), 1000, X_dvr(1), X_dvr(npoints_in) ) 
  end do


  sigma=0.0_wp
  sigma_tot=0.0_wp
  sigma_proj=0.0_wp
!  sigma_m =0.0_wp
  
  !sigma_elastic=0.0_wp
  !sigma_tot_elastic=0.0_wp
  !sigma_proj_elastic=0.0_wp

  do traj=1, npoints_x_mom_sampl
     
     !
     ! Compute trajectory
     !
             
     call hist_add(time_h_0, x_mom_sampl(traj,1), 1.0d0)   
     call hist_add(time_h_0_mom, x_mom_sampl(traj,2), 1.0d0)    
     
     call verlet_trajectory(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_dvr, E_n_inp * const % eV, delta_t, mu_SI, x_new )
     
     do i=1, ntsteps
        call hist_add(time_h(i), x_new(i), 1.0d0)
     end do
     
     ! look up energies as a function of distance (which in turn is a function of time)
     call spline_easy(X_dvr, E_n_inp, npoints_in, x_new, E_n, ntsteps)

     !call spline_easy(X_dvr, E_i_inp, npoints_in, x_new, E_f_elastic(1,:), ntsteps)
     
     do i=1,nfinal
        call spline_easy(X_dvr, E_f_inp(i,:), npoints_in, x_new, E_f(i,:), ntsteps)  
     end do
     
     do i=1,nfinal
        do m=1,3
           call spline_easy(X_dvr, D_fn_inp(i,:,m) , npoints_in, x_new, D_fn(i,:,m) , ntsteps)  
        end do
     end do
     
     
     ! now everything should be splined and ready. 
     ! compute  e^{-i \int_0^t W_I(t') dt' }
     ! first time, compute the mean transition energy, and frequency
     if (traj .eq. 1) then
        
        ind = minloc(E_i_inp)
        !                 E_fn_mean = sum(  E_f(nfinal,:) - E_n(:) ) / max(1,size(E_n(:))) 
        E_fn_mean =  E_f(nfinal,ind(1)) - E_n(ind(1))
        
        write(6,*) "E_fn_mean", E_fn_mean
        
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
        
        !do i=0,n_omega-1
        !   omega(i+1) = 2 * pi * i * hbar / (time_l2 * eV) - E_fn_mean
        !end do
        
     end if

     call compute_SCKH(E_n, E_f, E_fn_mean, D_fn, time,  sigma_m, gamma)

     ! static spectrum
     ! call compute_XES_spectrum_novib(E_n(1), E_f(:,1), E_fn_mean, D_fn(1), sigma_m_static, gamma)

     ! square sigma, must be done before summing over trajectories
     sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3)) 
     sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)
     
     ! compute projections
     do i=1,p % nproj
        sigma_tmp = sigma_m(:,:,1) * p % projvec(i,1) + sigma_m(:,:,2) * &
             p % projvec(i,2) + sigma_m(:,:,3) * p % projvec(i,3)
        sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
     end do
     
  end do ! end traj
  
  
  write(6,*)
  write(6,*) "Averaged over", npoints_x_mom_sampl, "trajectories"
  
!write distribution
file="_time_h_0"
file = trim(adjustl(p % outfile)) //  trim(adjustl(file))  // ".dat"
call hist_broadening(time_h_0, 0.01d-10)
call hist_write(time_h_0, file)

file="_time_h_0_mom"
file = trim(adjustl(p % outfile)) //  trim(adjustl(file))  // ".dat"
!call hist_broadening(time_h_0_mom, 0.01d-10)
call hist_write(time_h_0_mom, file)

do i=1,10 !ntsteps
   file="_time_h_"
   write(string,*) i
   file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   call hist_broadening(time_h(i), 0.01d-10)
   call hist_write(time_h(i), file)
end do

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

end subroutine calculate_SCKH_PES

end module m_SCKH_PES
