module m_SCKH_PES
  use KH_functions  
  use parameters  
  use SCKH_functions
  use m_XAS_io
  use spline_m
  use hist_class
  use FFT_m

  implicit none

contains

subroutine calculate_SCKH_PES(inp)
  type(input_params), intent(in):: inp 

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
  real(kind=wp):: gamma, time_l, delta_t, time_l2, norm, E_fn_mean, my_SI, dx
  integer:: n_omega, npoints_x_sampl, npoints_mom_sampl, npoints_x_mom_sampl
  complex(kind=wp), dimension(:,:),allocatable::   sigma_tmp
  complex(kind=wp), dimension(:,:,:),allocatable::  sigma_m
  real(kind=wp),dimension(:),allocatable:: eig_i, E_IP1s, E_i_inp, E_lp_corr, shift 
  type(hist), dimension(:), allocatable:: time_h
  type(hist):: time_h_0, time_h_0_mom
  integer, dimension(1)::ind  

  real(kind=wp), dimension(:,:,:,:),allocatable::  nac
  real(kind=wp), dimension(:),allocatable::  X_dvr
  real(kind=wp) :: dvr_start 
  integer:: npoints, ii

  ! set some local variables
  ntsteps = inp % ntsteps
  npoints_x_sampl = inp % npoints_x_sampl
  npoints_mom_sampl =  inp % npoints_mom_sampl
  npoints_in = inp % npoints_in
  nfinal = inp % npesfile_f
  my_SI = inp % my * amu

  npoints = (inp % npoints_in -1)/2
  dvr_start = inp % dvr_start_in * 1.0d-10
  dx = inp % dx_in * 1.0d-10

  ntsteps = inp % ntsteps
  delta_t = inp % delta_t

  !outfile = inp % outfile

  call next_power_of_2(ntsteps, ntsteps_pad, ntsteps_pad_pow)

  write(6,*) "ntsteps_pad_pow", ntsteps_pad, ntsteps_pad_pow

  if(inp % samplemode .eq. 1) then
     npoints_x_mom_sampl = npoints_x_sampl * npoints_mom_sampl
  else if(inp % samplemode .eq. 2) then
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

  if (inp % nonadiabatic .eq. 1) then
     allocate(nac(inp % npesfile_f, inp % npesfile_f, inp % npoints_in, 2) )
  end if
  
  ! use HWHM internally
  gamma = inp % gamma_FWHM / 2 


  ! set up DVR points
  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx + dvr_start
  end do

  
  ! read PES files
  call read_PES_file(inp % pes_file_i, inp % npoints_in, inp % npoints_in, X_dvr, E_i_inp)
  call read_PES_file(inp % pes_file_n, inp % npoints_in, inp % npoints_in, X_dvr, E_n_inp)

  do j=1,inp % npesfile_f
     call read_PES_file(inp % pes_file_f(j), inp % npoints_in, inp % npoints_in, X_dvr, E_f_inp(j,:))
     call read_dipole_file(inp % dipolefile_f(j), inp % npoints_in, inp % npoints_in, X_dvr, D_fn_inp(j,:,:))
  end do

  if (inp % nonadiabatic .eq. 1) then
     call read_nac_file(inp % nac_file, inp % npoints_in, inp %npoints_in, X_dvr, inp % npesfile_f, nac)
  end if

  ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
  ! and the spacing between the intermediate and final states are preserved

  if( inp % shift_PES .eq.  1) then
     call read_PES_file(inp % pes_file_lp_corr, inp % npoints_in, inp % npoints_in, X_dvr, E_lp_corr)

     shift = E_lp_corr -E_f_inp(1,:) 

     do j=1,inp % npesfile_f
        E_f_inp(j,:) = E_f_inp(j,:) + shift
     end do
     write(6,*) "Shifted PES:s"
  end if


  ! convert to eV units
  E_i_inp = E_i_inp  / eV
  E_n_inp = E_n_inp  / eV
     
  do j=1,nfinal
     E_f_inp(j,:) = E_f_inp(j,:) / eV
  end do
  
  write(6,*) E_n_inp(30) -E_f_inp(1,30)


  !
  ! Solve the vibrational problem for initial state to be able to sample initial distribution
  !  

  call solve_sinc_DVR(dx, my_SI, E_i_inp*eV, c_i, eig_i)
  write(6,*) "Calculated initial state eigenfunctions"
  write(6,*) "First vibrational transition", (eig_i(2)-eig_i(1)) * cm 


  if(inp % samplemode .eq. 1) then  
     call sample_x_mom(X_dvr, c_i(:,1), x_sampl, mom_sampl, x_mom_sampl, 1)
  else if(inp % samplemode .eq. 2) then  
     call sample_x_mom(X_dvr, c_i(:,1), x_sampl, mom_sampl, x_mom_sampl, 2)
  end if

  write(6,*) "Sampling done"

  !check
  do i=1,npoints_x_sampl
     write(17,*) x_sampl(i)
  end do

  do i=1,npoints_x_mom_sampl
     write(22,*) x_mom_sampl(i,1), x_mom_sampl(i,2)
  end do

  do i=1,npoints_in
     write(20,'(3ES16.6)') X_dvr(i), c_i(i,1), c_i(i,1) ** 2
  end do
  
  delta_t = delta_t * 1.d-15
  time_l = (ntsteps-1) * delta_t
  time_l2 = (ntsteps_pad-1) * delta_t 

  write(6,*) "outfile", inp % outfile
  ! some output
  write(6,*) "gamma (hwhm of lorentzian broadening)", gamma
  !write(6,*) "ntsteps", ntsteps
  !write(6,*) "tstep", tstep
  write(6,*) "my_SI", my_SI 
  write(6,*) "time_l", time_l
  write(6,*) "ntsteps_pad", ntsteps_pad, "= 2 **", ntsteps_pad_pow 
  write(6,*) "time_l2 (padded time)", time_l2
  write(6,*)
  write(6,*) "Fundamental frequency resolution", 2 * pi * hbar /( time_l * eV)
  write(6,*) "Padded frequency resolution", 2 * pi * hbar /( time_l2 * eV)
  write(6,*) "delta t", delta_t
  write(6,*) "max freq",  pi * hbar /( delta_t  * eV), "eV"


  n_omega = ntsteps_pad 
  allocate(sigma_m(nfinal,n_omega,3),&
       sigma(nfinal,n_omega), &
       sigma_tot(n_omega),  &
       sigma_proj(inp % nproj,n_omega), &
       sigma_tmp(nfinal,n_omega), &
       omega(n_omega))


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


  sigma=0
  sigma_tot=0
  sigma_proj=0

  do traj=1, npoints_x_mom_sampl
     
     !
     ! Compute trajectory
     !
             
     call hist_add(time_h_0, x_mom_sampl(traj,1), 1.0d0)   
     call hist_add(time_h_0_mom, x_mom_sampl(traj,2), 1.0d0)    
     
     call verlet_trajectory(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/my_SI, X_dvr, E_n_inp * eV, delta_t, my_SI, x_new )
     
     do i=1, ntsteps
        call hist_add(time_h(i), x_new(i), 1.0d0)
     end do
     
     ! look up energies as a function of distance (which in turn is a function of time)
     call spline_easy(X_dvr, E_n_inp, npoints_in, x_new, E_n, ntsteps)
     
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
           omega(j) = -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean
           j=j+1
        end do
        do i=0, n_omega/2-1
           omega(j) =  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean
           j=j+1
        end do
        
        !do i=0,n_omega-1
        !   omega(i+1) = 2 * pi * i * hbar / (time_l2 * eV) - E_fn_mean
        !end do
        
     end if
     
     call compute_SCKH(E_n, E_f, E_fn_mean, D_fn, time,  sigma_m, gamma)
     
     ! square sigma
     sigma = sigma + real( sum( sigma_m * conjg(sigma_m), 3)) 
     sigma_tot = sigma_tot +  sum( real( sum( sigma_m * conjg(sigma_m), 3)),1)
           
     ! compute projections
     do i=1,inp % nproj
        sigma_tmp = sigma_m(:,:,1) * inp % projvec(i,1) + sigma_m(:,:,2) * inp % projvec(i,2) + sigma_m(:,:,3) * inp % projvec(i,3)
        sigma_proj(i,:) = sigma_proj(i,:) + sum( real( sigma_tmp * conjg(sigma_tmp)),1)
     end do
     
  end do ! end traj
  
  
  write(6,*)
  write(6,*) "Averaged over", npoints_x_mom_sampl, "trajectories"
  
!write distribution
file="_time_h_0"
file = trim(adjustl(inp % outfile)) //  trim(adjustl(file))  // ".dat"
call hist_broadening(time_h_0, 0.01d-10)
call hist_write(time_h_0, file)

file="_time_h_0_mom"
file = trim(adjustl(inp % outfile)) //  trim(adjustl(file))  // ".dat"
!call hist_broadening(time_h_0_mom, 0.01d-10)
call hist_write(time_h_0_mom, file)

do i=1,10 !ntsteps
   file="_time_h_"
   write(string,*) i
   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   call hist_broadening(time_h(i), 0.01d-10)
   call hist_write(time_h(i), file)
end do

!normalize 
norm = sum(sigma_tot) *  2 * pi  * hbar / (time_l2 * eV ) 
sigma_tot = sigma_tot / norm
sigma = sigma / norm
sigma_proj = sigma_proj / norm

! write sigma to file
file="_sigma"
file = trim(adjustl(inp % outfile)) // trim(adjustl(file)) // ".dat"
open(10,file=file,status='unknown')
 

  !do i=n_omega/2, 1, -1
  !    write(10,*)  -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma(j,n_omega-i+1)
  ! end do
 !
  ! do i=0, n_omega/2
  !    write(10,*)  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma(j,i+1)
  ! end do

do i=1,n_omega
   write(10,*)  omega(i), sigma_tot(i)
!   write(10,*)  i, sigma_tot(i)
end do

close(10)

! write spectrum from different final states
do j=1, nfinal

   file="_sigma_final_"
   write(string,*) j
   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')
   
   do i=n_omega/2, 1, -1
      write(10,*)  -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma(j,n_omega-i+1)
   end do
 
   do i=0, n_omega/2
      write(10,*)  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma(j,i+1)
   end do
     
   close(10)
end do ! j

! write spectrum from projections  
do j=1, inp % nproj
  
   file="_sigma_proj_"
   write(string,*) j
   file = trim(adjustl(inp % outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do i=n_omega/2, 1, -1
      write(10,*)  -2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma_proj(j,n_omega-i+1)
   end do
     
   do i=0, n_omega/2
      write(10,*)  2 * pi * (i) * hbar / (time_l2 * eV) - E_fn_mean, sigma_proj(j,i+1)
   end do
     
   close(10)
end do !j
  
end subroutine calculate_SCKH_PES

end module m_SCKH_PES
