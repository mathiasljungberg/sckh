program DVR_harm
  use parameters
  use KH_functions
  implicit none

  ! input/output
  character(80)::pes_file_i, outfile
  character(80), dimension(:),allocatable:: pes_file_f,pes_file_n, dipolefile_f,dipolefile_n
  integer:: nstates,n_omega_in, n_omega_out, n_theta, npesfile_f, npesfile_n
  real(kind=wp):: my, dvr_start_in, dx_in, gamma_FWHM, instrument_FWHM
  real(kind=wp):: omega_in_start, omega_in_end, omega_out_start, omega_out_end, theta_start, theta_end

  ! loop variables
  integer::i,ii,j, n1,n2,f1,f2,m1,m2,om_in,om_out,i_theta 

  ! other variables
  character(80):: file,string
  real(kind=wp):: my_SI, dx,dvr_start, E_el_i, E_el_f, E_el_n, gamma, alpha, norm 
  real(kind=wp):: F_2_unpol, F_2_theta
  integer:: INFO, npoints

  real(kind=wp), dimension(:),allocatable:: W, WORK, X_dvr,V_i,  e_i
  real(kind=wp), dimension(:),allocatable:: x, omega_in, omega_out
  real(kind=wp), dimension(:,:),allocatable:: Mat_tmp, H_kin, V_dvr, c_i &
       ,e_n, e_f, V_n, V_f, sigma_unpol, lambda_F, lambda_G, lambda_H
  real(kind=wp), dimension(:,:,:),allocatable:: c_n, c_f, dipole_ni, dipole_nf, D_fi, D_ni, sigma_theta
  real(kind=wp), dimension(:,:,:,:,:),allocatable:: D_fn 
  complex(kind=wp), dimension(:,:),allocatable:: F

  !functions
  !real(kind=wp):: gaussian

!
! This progam calculates the cross section with the resonant Kramers-Heisenberg formula
!

!
! read input
!*******************************************************************
!
  read(5,*) pes_file_i

  write(6,*) pes_file_i

  ! intermediate states
  read(5,*) npesfile_n  
  
  write(6,*) npesfile_n  

  allocate(pes_file_n(npesfile_n), dipolefile_n(npesfile_n))
  
  do i=1,npesfile_n
     read(5,*) pes_file_n(i), dipolefile_n(i)

     write(6,*) pes_file_n(i), dipolefile_n(i)
  end do
  
  ! final states
  read(5,*) npesfile_f

  write(6,*) npesfile_f

  allocate(pes_file_f(npesfile_f), dipolefile_f(npesfile_f))
  
  do i=1,npesfile_f
     read(5,*) pes_file_f(i), dipolefile_f(i)

     write(6,*) pes_file_f(i), dipolefile_f(i)
  end do

  read(5,*) outfile

  write(6,*) outfile

  read(5,*) my
  read(5,*) nstates, dvr_start_in, dx_in
  read(5,*) omega_in_start, omega_in_end, n_omega_in 
  read(5,*) omega_out_start, omega_out_end, n_omega_out 
  read(5,*) theta_start, theta_end, n_theta
  read(5,*) gamma_FWHM, instrument_FWHM
  
gamma = gamma_FWHM / 2

write(6,*) "nstates", nstates
write(6,*) "gamma FWHM =",gamma_FWHM 
write(6,*) "gamma HWHM =",gamma 

! allocate everything
allocate( X_dvr(nstates), V_i(nstates), V_n(npesfile_n, nstates), &
     V_f(npesfile_f,nstates), e_i(nstates),e_n(npesfile_n,nstates),e_f(npesfile_f,nstates), &
     x(nstates) )

!allocate(Mat_tmp(nstates,nstates),H_kin(nstates,nstates))

allocate(V_dvr(nstates,nstates),c_i(nstates,nstates),&
     c_f(npesfile_f,nstates,nstates),c_n(npesfile_n,nstates,nstates) )

allocate(dipole_ni(npesfile_n, nstates, 3),dipole_nf(npesfile_f, nstates, 3), D_fn(npesfile_f,nstates,npesfile_n,nstates,3), &
     D_ni(npesfile_n,nstates,3))

allocate( omega_in(n_omega_in), omega_out(n_omega_out),  F(3,3), sigma_unpol(n_omega_out, n_omega_in),&
     sigma_theta(n_omega_out, n_omega_in, n_theta), lambda_F(npesfile_f,nstates), lambda_G(npesfile_f,nstates),lambda_H(npesfile_f,nstates))

write(6,*) shape( X_dvr)
write(6,*) shape(V_i)
write(6,*) shape(V_n)
write(6,*) shape(V_f) 
write(6,*) shape(e_i)
write(6,*) shape(e_n)
write(6,*) shape(e_f)
write(6,*) shape(x)

!allocate(Mat_tmp(nstates,nstates),H_kin(nstates,nstates))

write(6,*) shape(V_dvr)
write(6,*) shape(c_i)
write(6,*) shape(c_f)
write(6,*) shape(c_n)
write(6,*) shape(dipole_ni)
write(6,*) shape(dipole_nf)
write(6,*) shape(D_fn)
write(6,*) shape(D_ni)
write(6,*) shape( omega_in)
write(6,*) shape(omega_out)
write(6,*) shape(F)
write(6,*) shape(sigma_unpol)
write(6,*) shape(sigma_theta)
write(6,*) shape(lambda_F)
write(6,*) shape(lambda_G)
write(6,*) shape(lambda_H)


! allocate workvectors for lapack
!allocate( W(nstates), WORK(3*nstates))

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


! intermediate states
do j=1,npesfile_n

   !pes file
   open(10,file=pes_file_n(j),status='unknown')

   do i=1,nstates
      read(10,*) x(i), V_n(j,i)
      
      ! check that the points match the DVR points
      if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
         write(6,*) "f",j, " point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i) 
      end if
   end do

   close(10) 

   !dipole file
   
   open(10,file=dipolefile_n(j),status='unknown')
 
   do i=1,nstates
      read(10,*) x(i), dipole_ni(j,i,1),dipole_ni(j,i,2),dipole_ni(j,i,3)
      
      ! check that the points match the DVR points
      if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
         write(6,*) "d",j, "point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i) 
      end if
   end do

   close(10) 
   
end do


! correct core hole lone pair  PES
!if( shift_PES .eq.  1 ) then
!
!open(10,file=pes_file_lp_corr,status='unknown')
!
!  do i=1,nstates
!     read(10,*) x(i), V_lp_corr(i)
!
!     ! check that the points match the DVR points
!     if (abs(x(i)*1.0d-10-X_dvr(i)) .gt. 1.d-13) then
!        write(6,*) "corr point x(", i, ")=",x(i),"does not match the DVR point X_dvr(",i,")=",X_dvr(i)
!     end if
!  end do
!
!close(10)
!
!end if

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

   !dipole file, need N * F transitions
   
   open(10,file=dipolefile_f(j),status='unknown')
   do i=1,nstates
      read(10,*) x(i), dipole_nf(j,i,1),dipole_nf(j,i,2),dipole_nf(j,i,3)
      
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

!if( shift_PES .eq.  1) then
!   
!   shift = V_lp_corr -V_f(1,:) 
!   
!   do j=1,npesfile_f
!      V_f(j,:) = V_f(j,:) + shift
!   end do
!
!   write(6,*) "Shifted PES:s"
!
!end if


!create omega_in
call linspace(omega_in, omega_in_start,omega_in_end, n_omega_in ) 


!create omega_out
call linspace(omega_out, omega_out_start,omega_out_end, n_omega_out ) 

!write(6,*) "omega_in", omega_in
!write(6,*) "omega_out", omega_out

!do i=1,1000
   
   !write(6,*) gaussian(0.15d0, 0.0d0, 0.15d0 )

   !write(7,*) (i-1)*0.1d0, gaussian((i-1)*0.1d0, 0.9d0, 0.15d0  )

!end do

!stop

!
! Solve the vibrational problem for all eigenfunctions
!


! initial state
write(6,*) "so far..."

call solve_sinc_DVR(dx,my_SI, V_i, c_i, e_i)

write(6,*) "Calculated initial state eigenfunctions"

! intermediate states

do j=1,npesfile_n

write(6,*) "so why crash now?"
call solve_sinc_DVR(dx,my_SI, V_n(j,:), c_n(j,:,:), e_n(j,:))

write(6,*) "Calculated intermediate state eigenfunctions"

end do

! final states

do j=1,npesfile_f

call solve_sinc_DVR(dx,my_SI, V_f(j,:), c_f(j,:,:), e_f(j,:))

write(6,*) "Calculated final state", j

end do



!
! calculate dipole matrix elements between states, for now assume that there is no dependence of the int. el. state on the emission transition dipoles  
!

! absorption dipole matrix elements
do n1=1,npesfile_n
   do n2=1,nstates ! intermediate vib
      do m1=1,3
         !D_ni(n1,n2,m1) = sum(dipole_ni(n1,:,m1) * c_n(n1,:,n2) * c_i(:,1))   ! true dipole moment  
         D_ni(n1,n2,m1) = sum( c_n(n1,:,n2) * c_i(:,1))   ! FC
         !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
      end do
   end do
end  do


!do i=1,nstates 
!   D_ni(i) = sum(c_n(:,i)*c_i(:,1))
!end do

! emission dipole matrix elements
do f1=1,npesfile_f
   do f2=1,nstates ! final
      do n1=1,npesfile_n
         do n2=1,nstates ! intermediate vib
            do m1=1,3
               !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))  !sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))
               !D_fn(f1,f2,n1,n2,m1) =  sum(dipole_nf(f1,:,m1) * c_f(f1,:,f2) * c_n(n1,:,n2))   ! true dipole moment
               D_fn(f1,f2,n1,n2,m1) =  sum(c_f(f1,:,f2) * c_n(n1,:,n2))   ! Frank-Condon
               !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
            end do
         end do
      end do
   end  do
end do

! transitions from ground state directly to final states
!do f1=1,npesfile_f
!   do f2=1,nstates ! final
!      do m1=1,3
!         D_fi(f1,f2,l) = sum(dipole(f1,:,m1) * c_f(f1,:,f2) * c_i(:,1))   ! true dipole moment
!      end do
!   end do
!end  do

write(6,*) "Calculated dipole matrix elements"

!! calculate "mean energy" for intermediate state
!E_n_mean = sum(e_n(:) * D_ni(:) ** 2) 
!
!write(6,*) "sum (D_ni(:)**2)", sum( D_ni(:) ** 2)
!write(6,*) "E_n_mean", E_n_mean


!
! Full Kramers-Heisenberg    
!

write(6,*) "Full Kramers-Heisenberg"

!
! a and b components of the dipole matrix elements (polarization) 
! F_f^{a,b}(omega) = alpha * E_ni E_nf  sum_n D_ni^a D_fn^b / (omega - (e_n - e_i) + i*gamma)     
!                     

!  indices
!  f1 - electronic final state, f2 - vibratioanl final state
!  n1 - electronic intermediate state, n2 - vibratioanl intermediate state
!  m1 - polarization direction 1, m2 - polarization direction 2
!  om_in - index for omega, om_out - index for omega'



sigma_unpol=0
sigma_theta=0

do om_in =1, n_omega_in

   lambda_F =0 
   lambda_G =0
   lambda_H =0  
   do f1= 1,npesfile_f ! final el
      do f2=1,nstates ! final vib
         
         F=0

         do m1=1,3 ! polarization 1
            do m2=1,3 ! polarization 2
               
               do n1 = 1,npesfile_n ! intermediate electronic state 
                  do n2= 1,nstates ! intermediate vib
                     
                     F(m1,m2) = F(m1,m2) + (e_n(n1,n2) -e_i(1)) * (e_n(n1,n2) -e_f(f1,f2))  &
                          * D_fn(f1,f2,n1,n2,m1) * D_ni(n1,n2,m2) / ( omega_in(om_in) - &
                          (e_n(n1,n2) - e_i(1)) + dcmplx(0,gamma) )
                     
                  end do !n2
               end do !n1 
               
               ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
               lambda_F(f1,f2) = lambda_F(f1,f2) +  real(F(m1,m1) * conjg(F(m2,m2)) )
               lambda_G(f1,f2) = lambda_G(f1,f2) +  real(F(m1,m2) * conjg(F(m1,m2)) )
               lambda_H(f1,f2) = lambda_H(f1,f2) +  real(F(m1,m2) * conjg(F(m2,m1)) )
               
            end do !m2
         end do !m1
         

         ! <|F^2|> unpolarized
         F_2_unpol = (10.0d0/3.0d0 ) * lambda_G(f1,f2)
         
         ! convolute
         do om_out= 1, n_omega_out
            
            sigma_unpol(om_out,om_in) = sigma_unpol(om_out,om_in) + (omega_out(om_out)/omega_in(om_in))* F_2_unpol * &
                 gaussian(0.0d0, omega_in(om_in) -omega_out(om_out) - (e_f(f1,f2) -e_i(1)),  instrument_FWHM) 
            
            !write(7,*) gaussian(omega_out(om_out), omega_in(om_in) -omega_out(om_out) - (e_f(f1,f2) -e_i(1)),  instrument_FWHM), &
            !      omega_in(om_in) -omega_out(om_out) - (e_f(f1,f2) -e_i(1)), instrument_FWHM
            
         end do ! om_out

         !! angle dependence
         !do i_theta=1,n_theta
         !   
         !   ! <|F^2|> for angle theta
         !   F_2_theta = (-lambda_F(f1,f2) +4d0 *lambda_G(f1,f2) -lambda_H(f1,f2) ) &
         !        +(3d0 * lambda_F(f1,f2) -2d0 * lambda_G(f1,f2) +3d0 * lambda_H(f1,f2))*(cos(theta(i_theta)))**2
         !   
         !   ! construct cross section by convolution
         !   do om_out= 1, n_omega_out
         !      
         !      sigma_theta(om_out,om_in, i_theta) = sigma_theta(om_out,om_in, i_theta) + F_2_theta * &
         !           gaussian(omega_out(om_out), omega_in(om_in) -omega_out(om_out) - (e_(f1,f2) -e_i(1)),  FWHM_instr) 
         !      
         !   end do ! om_out
         !   
         !end do ! i_theta

      end do !f2
   end do ! f1

end do ! om_in 


write(6,*) "Calculated spectrum"

! normalize spectrum
do j=1,n_omega_in
   norm=sum(sigma_unpol(:,j)) *(omega_out(2) -omega_out(1)) 
   sigma_unpol(:,j) = sigma_unpol(:,j)/norm
end do

! write sigma to file
do j=1,n_omega_in

   file=outfile
   write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(file)) // trim(adjustl(string)) // ".dat"

   open(10,file=file,status='unknown')

   do i=1,n_omega_out
      write(10,*) omega_out(i), sigma_unpol(i,j)
   end do

   close(10) 
   
end do



file="vib_eigenvalues_intermediate.dat"
open(10,file=file,status='unknown')
do i=1,nstates
   write(10,*) e_n(1,i), 1 , D_ni(1,i,1)
end do
close(10)





! write sigma_theta to file...
!open(10,file=outfile,status='unknown')
!
!  do i=1,n_omega
!     write(10,*) omega(i), sigma(i)
!  end do
!
!close(10) 


!contains
!
!real(8) function gaussian(x, x0, FWHM)
!  use parameters
!  implicit none
!  
!  real(wp), intent(in):: x, x0, FWHM
!  !real(wp) :: gaussian
!  real(wp) :: alpha
!  
!  alpha=4.0_wp*log(2.0_wp)/(FWHM**2.0_wp)
!  gaussian = (alpha / pi)**0.5_wp * exp(-alpha*(x-x0)**2.0_wp)
!  
!end function gaussian
!
!
!subroutine solve_sinc_DVR(nstates,dx,my_SI, V_i, c_i, e_i)
!use parameters
!
!implicit none
!
!!passed variables
!
!integer, intent(in):: nstates
!real(kind=wp), intent(in):: dx, my_SI
!real(kind=wp), dimension(:), intent(in):: V_i
!real(kind=wp), dimension(nstates,nstates), intent(out):: c_i
!real(kind=wp), dimension(nstates), intent(out):: e_i
! 
!! local variables
!integer:: i,j,ii,jj,npoints
!!real(kind=wp), dimension(:,:), allocatable::Mat_tmp, H_kin
!real(kind=wp), dimension(nstates,nstates)::Mat_tmp, H_kin
!!!real(kind=wp):: 
!!real(kind=wp), dimension(nstates)::wnf_t,wnf_old
!! lapack
!!real(kind=wp), dimension(:), allocatable:: W
!real(kind=wp), dimension(nstates):: W
!!real(kind=wp), dimension(:), allocatable:: WORK
!real(kind=wp), dimension(nstates):: WORK
!integer:: INFO
!
!!allocate(Mat_tmp(nstates,nstates), H_kin(nstates,nstates), W(nstates),&
!!     WORK(3*nstates))
!
!!write(6,*) nstates,dx,my_SI, V_i, c_i, e_i
!
!
!npoints = (nstates-1)/2
!
!!write(6,*) shape(H_kin), shape(c_i), shape(e_i)
!write(6,*) "lalalalalalal"
!!write(6,*) shape(c_i), shape(e_i)
!
!!write(6,*) shape(Mat_tmp)
!
!!set up kinetic energy for sinc DVR                                                                                                                                                                                
!
!H_kin=0 
!
!do i = -npoints,npoints
!   ii = i + npoints +1
!   !write(6,*) "loop", i, ii
!   H_kin(ii,ii) = (hbar**2  / (2 * my_SI * dx **2)) * (pi ** 2) /3.0d0
!end do
!
!do i = -npoints,npoints
!   ii = i + npoints +1
!   !write(6,*) "second loop", i, ii, j, jj
!   do j = i +1, npoints
!      jj = j + npoints +1
!     H_kin(ii,jj) =  (hbar**2 * (-1)**(i-j) / ( my_SI * dx **2)) / (i-j) **2
!     H_kin(jj,ii) = H_kin(ii,jj)
!  end do
!end do
!
!write(6,*) "really..."
!
!! potential term                                                                                                                                                                                                    
!Mat_tmp = H_kin 
!!V_dvr = 0
!do i=1,nstates
!   !V_dvr(i,i) = V_i(i)*hartree
!   Mat_tmp(i,i) = Mat_tmp(i,i) + V_i(i)*hartree
!end do
!
!write(6,*) "well well.."
!
!! solve eigenvalue problem
!
!call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, 3*nstates, INFO)
!
!write(6,*) "solved it..", INFO
!
!! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
!
!
!c_i = Mat_tmp
!e_i = W / ev
!
!!deallocate(Mat_tmp, H_kin, W, WORK)
!
!
!end subroutine solve_sinc_DVR

end program DVR_harm
