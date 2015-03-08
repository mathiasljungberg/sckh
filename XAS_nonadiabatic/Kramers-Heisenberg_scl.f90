program DVR_harm
  use parameters
  implicit none

  ! input/output
  character(80)::pes_file_i,pes_file_n, pes_file_lp_corr, pes_file_ch_corr, outfile
  character(80), dimension(:),allocatable:: pes_file_f,dipolefile_f
  integer:: nstates,ntsteps,n_omega, npesfile_f, shift_PES, runmode, chdyn, broad
  real(kind=wp):: my, dvr_start_in,dx_in,tstep
  real(kind=wp):: omega_start, omega_end, lifetime, lifetime2 
  !real(kind=wp), dimension(:),allocatable:: x,y

  ! loop variables
  integer::i,j,ii,jj,k,l,m,t 

  ! other variables
  character(80):: file,string
  real(kind=wp):: my_SI, dx,dvr_start, E_el_i,E_el_f, E_el_n, gamma, gamma2,time,alpha, norm, E_n_mean 
  integer:: INFO,npoints, k_eq
  !real(kind=wp):: 
  real(kind=wp), dimension(:),allocatable:: W, WORK, X_dvr,V_i, V_n, V_lp_corr, V_ch_corr, e_i,e_n,e_ch,shift
  real(kind=wp), dimension(:),allocatable:: x, sigma, omega, D_ni, sigma_scl1, sigma_scl2, sigma_dir, sigma_max_int
  real(kind=wp), dimension(:,:),allocatable:: Mat_tmp, H_kin,V_dvr, c_i, c_n, c_ch &
       ,e_f, V_f,sigma_states, sigma_scl2_m, F2_dir, sigma_time, F2_max_int, sigma_lp &
       , delta_e_lp, sigma_dir_states, sigma_max_int_states
  real(kind=wp), dimension(:,:,:),allocatable:: c_f,dipole, D_fi, F2_dir_lp
  real(kind=wp), dimension(:,:,:,:),allocatable:: D_fn, sigma_scl2_m_states
! real(kind=wp), dimension(:,:),allocatable:: H_exc
  !complex(kind=16), dimension(:),allocatable:: wfn_t, wfn_old, d_wfn2, wfn2_t
  complex(kind=wp), dimension(:,:),allocatable:: F, sigma_scl1_m
  complex(kind=wp), dimension(:),allocatable:: wfn_t, wfn_old, d_wfn2, wfn2_t, wfn_w, wfn_w_d
  complex(kind=wp), dimension(:,:),allocatable:: wfn_all_t

 ! real(kind=wp), dimension(:,:),allocatable::  

  !functions
  !real(kind=wp):: harm_eigenfun
 ! real(kind=wp)::dlamch


!
! This progam propagates a one dimensional wave packet
!

!
! read input
!*******************************************************************
!
  read(5,*) pes_file_i
  read(5,*) pes_file_n   ! core ionized state
  read(5,*) chdyn, pes_file_ch_corr   ! ev. dynamics on core ionized state     
  read(5,*) shift_PES, pes_file_lp_corr        
  ! this is the the correct lone pair hole PES, if shift_PES=1 all PES:s will be shifted acoordingly 
  ! observe that the first pes_file_f must be the LUMO PES!
  ! if chdyn=1 pes_file_ch_corr will be used for the dynamics

  read(5,*) npesfile_f
  allocate(pes_file_f(npesfile_f), dipolefile_f(npesfile_f))
  
  do i=1,npesfile_f
     read(5,*) pes_file_f(i), dipolefile_f(i)
  end do

  read(5,*) outfile
  read(5,*) my
  read(5,*) nstates, dvr_start_in, dx_in
  read(5,*) omega_start, omega_end, n_omega
  read(5,*) lifetime, lifetime2, broad
  read(5,*) runmode, tstep,ntsteps

  ! runmode 1: K-H, 
  !  2: K-H approximate, final state as delta functions   
  !  3: K-H approximate, intermediate state as delta functions   
  !  4: semiclassical 1 
  ! lifetime2 is the broadening used for scl simlations (should be equal to gamma)
  ! broad, 1 = lorentzian, 2= gaussian: for seimclassical simulation


! 6.626e-34/(2*pi*3.6e-15*1.60217e-19)
! ans =  0.18284 eV   for 3.6 fs FWHM lifetime

! in this program gamma is the HWHM!! but the lifetime is defined as hbar / (gamma') = hbar / (2*gamma) 
! with gamma' being the FWHM. so we have gamma = hbar /  (2 * lifetime)
! use lifetime = 3.6 fs for correct results! for comparison to old (erronous) spectra
! use lifetime = 1.8 and lifetime2 = 3.6. Never use lifetime =7.2 (definitions have changed)

gamma = h / (2 * 2 * pi * lifetime * 1d-15 * ev)
gamma2 = h / (2 * 2 * pi * lifetime2 * 1d-15 * ev)

write(6,*) "gamma HWHM =",gamma 
write(6,*) "gamma2 HWHM =",gamma2 


! allocate everything
allocate( X_dvr(nstates), V_i(nstates), V_n(nstates), V_lp_corr(nstates), V_ch_corr(nstates), &
     V_f(npesfile_f,nstates), e_i(nstates),e_n(nstates),e_f(npesfile_f,nstates), e_ch(nstates), &
     x(nstates), shift(nstates))

allocate( omega(n_omega), sigma(n_omega), F(n_omega,3), D_ni(nstates) )

allocate( sigma_scl1(n_omega), sigma_scl1_m(n_omega,3),sigma_scl2(n_omega), sigma_scl2_m(n_omega,3), &
     sigma_scl2_m_states(npesfile_f, ntsteps,n_omega,3), sigma_time(ntsteps,n_omega))

allocate(Mat_tmp(nstates,nstates),H_kin(nstates,nstates),V_dvr(nstates,nstates),c_i(nstates,nstates),&
     c_f(npesfile_f,nstates,nstates),c_n(nstates,nstates), c_ch(nstates,nstates),D_fn(npesfile_f,nstates,nstates,3))

allocate(dipole(npesfile_f, nstates, 3), sigma_states(npesfile_f,n_omega))

allocate(wfn_all_t(ntsteps,nstates))

allocate(F2_dir(n_omega,3), sigma_dir(n_omega), F2_max_int(n_omega,3), sigma_max_int(n_omega), &
     D_fi(npesfile_f,nstates,3) )

allocate(F2_dir_lp(10,n_omega,3), sigma_lp(10,n_omega), delta_e_lp(10,nstates))

allocate(sigma_dir_states(npesfile_f,n_omega), sigma_max_int_states(npesfile_f,n_omega))


! allocate workvectors for lapack
allocate( W(nstates), WORK(3*nstates))

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
! read pes from files
!

!initial state
open(10,file=pes_file_i,status='unknown')

  do i=1,nstates
     read(10,*) x(i), V_i(i)
     
     ! check that the points match the DVR points
     if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
        write(6,*) "i point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i) 
     end if
  end do

close(10) 


! intermediate state, koopmans
open(10,file=pes_file_n,status='unknown')

  do i=1,nstates
     read(10,*) x(i), V_n(i)
     
     ! check that the points match the DVR points
     if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
        write(6,*) "n point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i) 
     end if
  end do

close(10) 


! intermediate state, correct core ionized  PES for dynamics
if( chdyn .eq.  1 ) then

open(10,file=pes_file_ch_corr,status='unknown')

  do i=1,nstates
     read(10,*) x(i), V_ch_corr(i)

     ! check that the points match the DVR points                                                                                                                                                                 
     if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
        write(6,*) "corr point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i)
     end if
  end do

end if

! correct core hole lone pair  PES
if( shift_PES .eq.  1 ) then

open(10,file=pes_file_lp_corr,status='unknown')

  do i=1,nstates
     read(10,*) x(i), V_lp_corr(i)

     ! check that the points match the DVR points
     if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
        write(6,*) "corr point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i)
     end if
  end do

close(10)

end if

! final states

do j=1,npesfile_f

   !pes file
   open(10,file=pes_file_f(j),status='unknown')

   do i=1,nstates
      read(10,*) x(i), V_f(j,i)
      
      ! check that the points match the DVR points
      if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
         write(6,*) "f",j, " point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i) 
      end if
   end do

   close(10) 

   !dipole file
   
   open(10,file=dipolefile_f(j),status='unknown')
 
   do i=1,nstates
      read(10,*) x(i), dipole(j,i,1),dipole(j,i,2),dipole(j,i,3)
      
      ! check that the points match the DVR points
      if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
         write(6,*) "d",j, "point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i) 
      end if
   end do

   close(10) 
   
end do

write(6,*) "Done reading"


! Shift orbital energies so that V_n have energies V_lp_corr
! and the spacing between the intermediate and final states are preserved

if( shift_PES .eq.  1) then
   
   shift = V_lp_corr -V_f(1,:) 
   
   do j=1,npesfile_f
      V_f(j,:) = V_f(j,:) + shift
   end do

   write(6,*) "Shifted PES:s"

end if


!create omega
call linspace(omega, omega_start,omega_end, n_omega ) 



!
! Solve the vibrational problem for all eigenfunctions
!


! initial state

call solve_sinc_DVR(nstates,dx,my_SI, V_i, c_i, e_i)

write(6,*) "Calculated initial state eigenfunctions"

! intermediate state

call solve_sinc_DVR(nstates,dx,my_SI, V_n, c_n, e_n)

write(6,*) "Calculated intermediate state eigenfunctions"


! core hole state for dynamics

if( chdyn .eq.  1) then

   call solve_sinc_DVR(nstates,dx,my_SI, V_ch_corr, c_ch, e_ch)
   
   write(6,*) "Calculated intermediate state eigenfunctions"

end if


if(runmode .eq. 1 .or. runmode .eq. 2 .or. runmode .eq. 3) then  ! K-H

write(6,*) "Kramers-Heisenberg"

! final states

do j=1,npesfile_f

call solve_sinc_DVR(nstates,dx,my_SI, V_f(j,:), c_f(j,:,:), e_f(j,:))

write(6,*) "Calculated final state", j

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
            !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))  !sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))
            D_fn(i,j,k,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
            !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
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

! calculate "mean energy" for intermediate state
E_n_mean = sum(e_n(:) * D_ni(:) ** 2) 

write(6,*) "sum (D_ni(:)**2)", sum( D_ni(:) ** 2)
write(6,*) "E_n_mean", E_n_mean


!
! Full Kramers-Heisenberg    
!

if (runmode .eq. 1) then

write(6,*) "Full Kramers-Heisenberg"

!for each energy, omega, calculate F(omega), the cross section


! sigma(omega) = sum_f [F_f]**2 
!
! F_f(omega) = sum_n <i|D|n><n|D|f> / (omega - (e_n - e_f) + i*gamma)  
!                     

sigma=0
sigma_dir=0
sigma_max_int = 0

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
                    (e_n(k) - e_f(j,l)) + dcmplx(0,gamma) )
               
               ! direct contribution (no interference)
               F2_dir(i,m) = F2_dir(i,m) + D_fn(j,l,k,m) ** 2 * D_ni(k) ** 2 / (( omega(i) - &
                    (e_n(k) - e_f(j,l)))**2  + gamma**2 ) 

            end do
         end do
      
         ! maximal interference, direct transition from initial to final states
         do m=1,3 ! polarization 
            
            F2_max_int(i,m) =  F2_max_int(i,m) + D_fi(j,l,m) ** 2 / (( omega(i) - &
                    (E_n_mean - e_f(j,l)) )**2  + gamma**2 )
         end do

      end do

      !square F
      
      sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
           + real(conjg(F(:,3))*F(:,3)) 
      
      !sigma_dir = sigma_dir + F2_dir(:,1) +F2_dir(:,2) +F2_dir(:,3)

      sigma_dir_states(j,:) = sigma_dir_states(j,:) + F2_dir(:,1) +F2_dir(:,2) +F2_dir(:,3)

      !sigma_max_int = sigma_max_int + F2_max_int(:,1) + F2_max_int(:,2) + F2_max_int(:,3)
      
      sigma_max_int_states(j,:) = sigma_max_int_states(j,:) + F2_max_int(:,1) + F2_max_int(:,2) + F2_max_int(:,3)

      !sigma = sigma +real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
      !     + real(conjg(F(:,3))*F(:,3)) 
      
      !sigma = sigma + sigma_states(j,:)
      
   end do !l

   sigma = sigma + sigma_states(j,:)
   sigma_dir = sigma_dir + sigma_dir_states(j,:)
   sigma_max_int = sigma_max_int + sigma_max_int_states(j,:)

end do ! j

! for lone pair state, compute the individual vibrational transitions
! in the case of no interference

j=1
sigma_lp=0
do l=1,nstates ! final vib
         
   F2_dir_lp=0

   do i =1, n_omega 
      do k=1,10 ! intermediate vib, only 10 first states
         do m=1,3 ! polarization
            ! direct contribution (no interference)
            F2_dir_lp(k,i,m) = F2_dir_lp(k,i,m) + D_fn(j,l,k,m) ** 2 * D_ni(k) ** 2 / (( omega(i) - &
                 (e_n(k) - e_f(j,l)))**2  + gamma**2 ) 
         end do
      end do
   end do

   !save transition energies
   do k=1,10 ! intermediate vib, only 10 first states
      delta_e_lp(k,l) =  e_n(k) - e_f(j,l)
   end do


   sigma_lp = sigma_lp + F2_dir_lp(:,:,1) +F2_dir_lp(:,:,2) +F2_dir_lp(:,:,3)
   
end do !l



else if(runmode .eq. 2) then  

write(6,*) "Approxiamte Kramers-Heisenberg: final vibrational states as delta functions"

!
! Approximative K-H, approximate final state vibrational wave functions as delta functions
!

sigma=0
sigma_dir=0

do j= 1,npesfile_f ! final el
   sigma_states(j,:) =0
   do l=1,nstates  ! sum over R here
      

      F=0
      do i =1, n_omega 
         do k=1,nstates ! intermediate vib
            do m=1,3 ! polarization
               
               !F(i,m) = F(i,m) + D_fn(j,l,k,m) * D_ni(k) / ( omega(i) - &
               !     (e_n(k) - e_f(j,l)) + dcmplx(0,gamma) )

               F(i,m) = F(i,m) + dipole(j,l,m) * D_ni(k) * c_n(l,k)  / ( omega(i) - &
                    (e_n(k) - V_f(j,l)*hartree/eV) + dcmplx(0,gamma) )
               
            end do
         end do
      end do

      !square F
      
      sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
           + real(conjg(F(:,3))*F(:,3)) 
      
   end do !l

   sigma = sigma + sigma_states(j,:)

end do ! j

else if(runmode .eq. 3) then  

write(6,*) "Approxiamte Kramers-Heisenberg: final vibrational states as delta functions"


!
! Approximative K-H, approximate intermediate state vibrational wave functions as delta functions
!

sigma=0
sigma_dir=0

do j= 1,npesfile_f ! final el
   sigma_states(j,:) =0
   do l=1,nstates  ! final vib
      

      F=0
      do i =1, n_omega 
         do k=1,nstates !  sum over R 
            do m=1,3 ! polarization
               
               !F(i,m) = F(i,m) + D_fn(j,l,k,m) * D_ni(k) / ( omega(i) - &
               !     (e_n(k) - e_f(j,l)) + dcmplx(0,gamma) )

               !F(i,m) = F(i,m) + dipole(j,l,m) * D_ni(k) * c_n(l,k)  / ( omega(i) - &
               !     (e_n(k) - V_f(j,l)*hartree/eV) + dcmplx(0,gamma) )

               !write(6,*) minloc(V_i), V_n(k)*hartree/eV, V_n(minloc(V_i))*hartree/eV
               !write(6,*) ( omega(i) - &
               !     ( V_n(k)*hartree/eV - e_f(j,l)) + dcmplx(0,gamma) )

               !write(6,*) ( omega(i) - &
               !     ( V_n(minloc(V_i))*hartree/eV - e_f(j,l)) + dcmplx(0,gamma) )

               !F(i,m) = F(i,m) + dipole(j,k,m) * c_f(j,k,l) * c_i(k,1)  / ( omega(i) - &
               !     ( V_n(k)*hartree/eV - e_f(j,l)) + dcmplx(0,gamma) )
             
               !k_eq = minloc(V_i,1)

               F(i,m) = F(i,m) + dipole(j,k,m) * c_f(j,k,l) * c_i(k,1)  / ( omega(i) - &
                    ( V_n(k) * hartree / eV - e_f(j,l)) + dcmplx(0,gamma) )

  
            end do
         end do
      end do

      !square F
      
      sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
           + real(conjg(F(:,3))*F(:,3)) 
      
   end do !l

   sigma = sigma + sigma_states(j,:)

end do ! j



end if  ! runmode 1 or 2 or 3




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


! write sigma to file
open(10,file=outfile,status='unknown')

  do i=1,n_omega
     write(10,*) omega(i), sigma(i)
  end do

close(10) 


! write sigma_dir to file
file="sigma_direct.dat"

open(10,file=file,status='unknown')

do i=1,n_omega
   write(10,*) omega(i), sigma_dir(i)
end do

close(10)

! write sigma_max_int to file
file="sigma_max_int.dat"

open(10,file=file,status='unknown')

do i=1,n_omega
   write(10,*) omega(i), sigma_max_int(i)
end do

close(10)

! write sigma_states
do j= 1,npesfile_f
   file="sigma_states_"

   write(string,*) j

   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

   open(10,file=file,status='unknown')
   
   do i=1,n_omega
      write(10,*) omega(i), sigma_states(j,i)
   end do

   close(10)
end do

! write sigma_dir_states
do j= 1,npesfile_f
   file="sigma_dir_states_"

   write(string,*) j

   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

   open(10,file=file,status='unknown')
   
   do i=1,n_omega
      write(10,*) omega(i), sigma_dir_states(j,i)
   end do

   close(10)
end do

! write sigma_max_int_states
do j= 1,npesfile_f
   file="sigma_max_int_states_"

   write(string,*) j

   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

   open(10,file=file,status='unknown')
   
   do i=1,n_omega
      write(10,*) omega(i), sigma_max_int_states(j,i)
   end do

   close(10)
end do


! write sigma_lp to file
do j= 1,10
   file="sigma_lp_"

   write(string,*) j
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

   open(10,file=file,status='unknown')
   do i=1,n_omega
      write(10,*) omega(i), sigma_lp(j,i)
   end do

   close(10)
end do

! write vibrational eigenvalues
file="vib_eigenvalues_initial.dat"
open(10,file=file,status='unknown')
do i=1,nstates
   write(10,*) e_i(i), 1 ,10
end do
close(10)

file="vib_eigenvalues_intermediate.dat"
open(10,file=file,status='unknown')
do i=1,nstates
   write(10,*) e_n(i), 1 , D_ni(i)
end do
close(10)

do j= 1,npesfile_f
   file="vib_eigenvalues_final_"
   
   write(string,*) j
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')
   
   do i=1,nstates
      write(10,*) e_f(j,i), 1, 10
   end do
   
   close(10)
end do

! write transition moments

file="transition_moments_i_to_n.dat"
open(10,file=file,status='unknown')
do i=1,nstates
   write(10,*) D_ni(i)
end do
close(10)

!file="transition_moments_n_to_f"
!do i=1,npesfile_f
!   do j=1,nstates ! final
!      do k=1,nstates ! intermediate
!         do l=1,3
!            D_fn(i,j,k,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
!         end do
!      end do
!   end do
!end  do


! write transition energies for the 10 first intermediate states
do j= 1,10
   file="delta_e_lp_"
   write(string,*) j
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

   open(10,file=file,status='unknown')

   do i=1,nstates
      write(10,*) delta_e_lp(j,i)
   end do
   close(10)
end do


! write vibrational eigenfucntions
if(1) then 

file="vib_eigenfunctions_initial.dat"
open(10,file=file,status='unknown')
do i=1,nstates
   do j=1,nstates
      write(10,*) x(j), c_i(j,i)
   end do
   write(10,*) 
   write(10,*) 
end do
close(10)

file="vib_eigenfunctions_intermediate.dat"
open(10,file=file,status='unknown')
do i=1,nstates
   do j=1,nstates
      write(10,*) x(j), c_n(j,i)
   end do
   write(10,*) 
   write(10,*) 
end do
close(10)

! write 10 first eigenfunctions to seprate files
do i=1,10
   file="vib_eigenfunctions_intermediate_"
   write(string,*) i
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do j=1,nstates
      write(10,*) x(j), c_n(j,i)
   end do
   write(10,*)
   write(10,*)
   close(10)
end do




! write 10 first eigenfunctions of first final state to separate files
do i=1,10
   file="vib_eigenfunctions_final_1_"
   write(string,*) i
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do j=1,nstates
      write(10,*) x(j), c_f(1,j,i)
   end do
   write(10,*)
   write(10,*)
   close(10)
end do



do k= 1,npesfile_f
   file="vib_eigenfunctions_final_"
   write(string,*) k
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do i=1,nstates
      do j=1,nstates
         write(10,*) x(j), c_f(k,j,i)
      end do
      
      write(10,*) 
      write(10,*) 
      
   end do
   close(10)
end do


end if

do j= 1,npesfile_f
   file="vib_eigenvalues_final_"

   write(string,*) j
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
   open(10,file=file,status='unknown')

   do i=1,nstates
      write(10,*) e_f(j,i)
   end do

   write(10,*) 
   write(10,*) 

   close(10)
end do




!end if ! runmode =1 , else 2, 3



else if(runmode .eq. 4) then  ! semiclassical 1

write(6,*) "Semiclassical approximation 2"






!
! Semiclassical approximations
!


! Wave packet dynamics on intermediate state

if( chdyn .eq.  1) then

   call wpd(ntsteps,tstep,nstates, c_i, c_ch, e_ch, wfn_all_t )

else

   call wpd(ntsteps,tstep,nstates, c_i, c_n, e_n, wfn_all_t )

end if



write(6,*) "Wpd done"



if(0) then


! here, do semiclassical calculation according to
! \sigma(\omega) = sum_R sum_t |psi(R,t)|^2 |D_{NF}(R)|^2 e^{i(omega -(E_F(R) - E_N(R))) t} \frac{i}{omega -(E_F(R) - E_N(R)) +i\Gamma}
! this assumes only that |n_N>  and  |f_F> = |R>, 

!dipole(j,i,1),dipole(j,i,2),dipole(j,i,3) , first index final state PES, second index space coordinate

sigma_scl1_m=0
do j=1,npesfile_f 
   do k=1,nstates ! R
      do t=1,ntsteps !t
         time = (t-1) * tstep *1.d-15
         do i=1,n_omega
            do m=1,3 ! polarization
               !sigma_scl1(i,m) = sigma_scl1(i,m) + wfn_all_t(t,k) ** 2 * dipole(j,k,m) ** 2 * &
               !     exp(dcmplx(0, omega(i) - (V_n(k) -V_f(j,k)))) * &
               !     (dcmplex(0,1)/( omega(i) - (V_n(k) -V_f(j,k)) + dcmplx(0, Gamma))
               sigma_scl1_m(i,m) = sigma_scl1_m(i,m) + wfn_all_t(t,k) ** 2 * dipole(j,k,m) ** 2 * &
                    exp(dcmplx(-gamma, omega(i) - (V_n(k) -V_f(j,k)))* time) * &
                    (dcmplx(0,1)/( omega(i) - (V_n(k) -V_f(j,k)) + dcmplx(0, gamma)))
            end do
         end do
      end do
   end do
end do

end if ! if (0)

! sigma_scl1 =  real(sigma_scl1_m(:,1) +sigma_scl1_m(:,2) + sigma_scl1_m(:,3)) 



! here, do semiclassical calculation according to
! \sigma(\omega) = sum_R sum_t |psi(R,t)|^2 |D_{NF}(R)|^2 e^{-\Gamma t} \frac{1}{ (E_F(R) -E_N(R) - \hbar omega)^2 + \Gamma^2}
! or  \delta(E_F(R) -E_N(R) - \hbar omega)
! this is equal to the (uncontrolled) approximation we use 


sigma_scl2_m=0

if(broad .eq. 1) then  ! lorentzian broadening

do j=1,npesfile_f 
   
   do t=1,ntsteps !t
      time = (t-1) * tstep  !*1.d-15

      sigma_scl2_m_states(j,t,:,:) = 0

      do k=1,nstates ! R 
         do i=1,n_omega
            do m=1,3 ! polarization
               ! write(6,*) j,k,t,i,m
               
               !sigma_scl2_m(i,m) = sigma_scl2_m(i,m) + dble(conjg(wfn_all_t(t,k)) * wfn_all_t(t,k)) * dipole(j,k,m) ** 2 * exp(-time / lifetime) * &
               !     1.0_wp/(  (omega(i) - (V_n(k) -V_f(j,k))*hartree/eV ) ** 2 + gamma2 ** 2)
               
               sigma_scl2_m_states(j,t,i,m) = sigma_scl2_m_states(j,t,i,m) + dble(conjg(wfn_all_t(t,k)) * wfn_all_t(t,k)) * dipole(j,k,m) ** 2 * &
                    1.0_wp/(  (omega(i) - (V_n(k) -V_f(j,k))*hartree/eV ) ** 2 + gamma2 ** 2)
               
               !sigma_scl2_m(i,m) = sigma_scl2_m(i,m) + sigma_scl2_m_states(j,t,i,m) * exp(-time / lifetime) 
        
            end do !m
         end do ! i
      end do ! k   
    end do  !t
 end do ! j

else if(broad .eq. 2) then   ! gaussian broadening, WATCH OUT, NOT REALLY THE SAME AS ABOVE!

alpha=4.0_wp*log(2.0_wp)/((2.0_wp * gamma)**2.0_wp)          !alpha=4.0_wp*log(2.0_wp)/(fwhm**2.0_wp)

do j=1,npesfile_f 
   do t=1,ntsteps !t
      time = (t-1) * tstep  !*1.d-15
      do k=1,nstates ! R
         do i=1,n_omega
            do m=1,3 ! polarization
!               write(6,*) j,k,t,i,m
               !sigma_scl2(i,m) = sigma_scl2(i,m) + wfn_all_t(t,k) ** 2 * dipole(j,k,m) ** 2 * &
               !     1.0_wp/(  (omega(i) - (V_n(k) -V_f(j,k))) ** 2 + Gamma ** 2)
               sigma_scl2_m(i,m) = sigma_scl2_m(i,m) + dble(conjg(wfn_all_t(t,k)) * wfn_all_t(t,k)) * dipole(j,k,m) ** 2 * exp(-time / (2* lifetime)) * &
                   (alpha / pi)**0.5_wp * exp(-alpha*(omega(i) - (V_n(k) -V_f(j,k))*hartree/eV )**2.0_wp)  ! 1.0_wp/(  (omega(i) - (V_n(k) -V_f(j,k))*hartree/eV ) ** 2 + gamma2 ** 2)

               

            end do
         end do
      end do
   end do
end do


end if

! compute sigma_states and sigma_time

do j=1,npesfile_f
   sigma_states(j,:)=0
   do t=1,ntsteps   
      time = (t-1) * tstep
      sigma_states(j,:) = sigma_states(j,:) + &
           (sigma_scl2_m_states(j,t,:,1) +sigma_scl2_m_states(j,t,:,2) +sigma_scl2_m_states(j,t,:,3)) *  exp(-time / (2 * lifetime))
   end do
end do

sigma_scl2 =0
do j=1,npesfile_f
   sigma_scl2 = sigma_scl2 +  sigma_states(j,:)
end do

do t=1,ntsteps   
   sigma_time(t,:) =0
   do j=1,npesfile_f  

      sigma_time(t,:) = sigma_time(t,:) + &
           sigma_scl2_m_states(j,t,:,1) +sigma_scl2_m_states(j,t,:,2) +sigma_scl2_m_states(j,t,:,3)
   end do
end do


write(6,*) "Calulcated spectrum"

!sigma_scl2 =  sigma_scl2_m(:,1) +sigma_scl2_m(:,2) + sigma_scl2_m(:,3)

norm = sum(sigma_scl2) *(omega(2) -omega(1)) 
! normalize spectrum
sigma_scl2 = sigma_scl2/norm

do j= 1,npesfile_f   
   sigma_states(j,:) = sigma_states(j,:)/norm
end do


! write sigma to file
open(10,file=outfile,status='unknown')

  do i=1,n_omega
     write(10,*) omega(i), sigma_scl2(i)
  end do

close(10) 

! write wavepacket squared (weights) to trajectory file
do t=1,ntsteps 
   write(11,*) "t=",(t-1)*.1d-15  !XXXX why .1 here?
   do i = 1,nstates
      write(11,*) x(i), dble(conjg(wfn_all_t(t,i)) * wfn_all_t(t,i))
   end do
   write(11,*)
end do

! write wavepacket to trajectory file
do t=1,ntsteps 
   write(13,*) "t=",(t-1)*.1d-15
   write(14,*) "t=",(t-1)*.1d-15
   do i = 1,nstates
      write(13,*) x(i), dble(wfn_all_t(t,i))
      write(14,*) x(i), aimag(wfn_all_t(t,i))
   end do
   write(13,*)
   write(14,*)
end do

! write sigma states to file

do j= 1,npesfile_f
  do i=1,n_omega
     write(20+j,*) omega(i), sigma_states(j,i)
  end do
end do


! write sigma time to file

do t= 1,ntsteps
  do i=1,n_omega
     write(30+t,*) omega(i), sigma_time(t,i)
  end do
end do


end if ! runmode 

end program DVR_harm

!contains

subroutine wpd(ntsteps,tstep,nstates, c_gs, c_exc, E_exc, wfn_all_t )
use parameters 
implicit none

!passed variables
integer, intent(in):: ntsteps, nstates
real(kind=wp), intent(in):: tstep
real(kind=wp), dimension(nstates,nstates), intent(in)::c_gs,c_exc
real(kind=wp), dimension(nstates), intent(in):: E_exc
complex(kind=wp), dimension(ntsteps,nstates),intent(out):: wfn_all_t
!real(kind=wp), dimension(ntsteps,nstates), intent(out)::wnf_all_t
 
! local variables
integer:: t,i,j,k
real(kind=wp):: time
!real(kind=wp), dimension(nstates)::wnf_t,wnf_old
complex(kind=wp), dimension(nstates):: wfn_t, wfn_old


write(6,*) ntsteps,tstep,nstates
! calculate wfn at time t, femtoseconds
do t = 1, ntsteps
   time = (t-1) * tstep *1.d-15
   write(6,*) time
   
   if (t .eq. 1) then
      wfn_old = c_gs(:,1)
   else 
      wfn_old = wfn_t
   end if

   wfn_t = 0 
   do i = 1,nstates
      do j=1,nstates 
         do k=1,nstates
            wfn_t(k) = wfn_t(k) + c_gs(i,1) * c_exc(i,j) * c_exc(k,j) * exp(- dcmplx(0, E_exc(j) * ev * time / hbar) ) 
         end do
      end do
   end do

   wfn_all_t(t,:) = wfn_t

end do

end subroutine wpd


subroutine solve_sinc_DVR(nstates,dx,my_SI, V_i, c_i, e_i)
use parameters

implicit none

!passed variables
integer, intent(in):: nstates
real(kind=wp), intent(in):: dx, my_SI
real(kind=wp), dimension(nstates), intent(in):: V_i
real(kind=wp), dimension(nstates,nstates), intent(out):: c_i
real(kind=wp), dimension(nstates), intent(out):: e_i
 
! local variables
integer:: i,j,ii,jj,npoints
real(kind=wp), dimension(nstates,nstates)::Mat_tmp, H_kin
!real(kind=wp):: 
real(kind=wp), dimension(nstates)::wnf_t,wnf_old
! lapack
real(kind=wp), dimension(nstates):: W
real(kind=wp), dimension(3*nstates):: WORK
integer:: INFO


npoints = (nstates-1)/2

!set up kinetic energy for sinc DVR                                                                                                                                                                                

H_kin=0 
do i = -npoints,npoints
   ii = i + npoints +1
   H_kin(ii,ii) = (hbar**2  / (2 * my_SI * dx **2)) * (pi ** 2) /3.0d0
end do

do i = -npoints,npoints
   ii = i + npoints +1
   do j = i +1, npoints
      jj = j + npoints +1
     H_kin(ii,jj) =  (hbar**2 * (-1)**(i-j) / ( my_SI * dx **2)) / (i-j) **2
     H_kin(jj,ii) = H_kin(ii,jj)
  end do
end do



! potential term                                                                                                                                                                                                    
Mat_tmp = H_kin 
!V_dvr = 0
do i=1,nstates
   !V_dvr(i,i) = V_i(i)*hartree
   Mat_tmp(i,i) = Mat_tmp(i,i) + V_i(i)*hartree
end do



! solve eigenvalue problem

call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, 3*nstates, INFO)

! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)

c_i = Mat_tmp
e_i = W / ev

end subroutine solve_sinc_DVR


