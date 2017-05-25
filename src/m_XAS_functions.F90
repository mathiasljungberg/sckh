module m_XAS_functions

  implicit none

contains

subroutine solve_non_adiabatic(eig_f, c_f, nac, eig_na, c_na, M_dsinc ) 
  use m_precision, only: wp
  use m_KH_functions, only: delta
  
  real(kind=wp), intent(in):: eig_f(:,:), c_f(:,:,:), nac(:,:,:,:), M_dsinc(:,:)
  real(kind=wp), intent(out):: eig_na(:), c_na(:,:,:)

  !! local variables
  integer:: i,j,i1,j1,i2,j2, ii,jj
  integer:: npesfiles_f,npoints,nstates, nstates_na
  real(kind=wp), dimension(:,:), allocatable::Mat_tmp
  real(kind=wp), dimension(:),allocatable:: W
  real(kind=wp), dimension(:),allocatable:: WORK
  integer:: INFO, LWORK

  npesfiles_f = size(eig_f,1)
  nstates = size(eig_f,2)
  nstates_na = size(eig_na,1)

  write(6,*) "solving nonadiabatic matrix elements", nstates_na

  ! M_dsinc are the matrix elements <i| d/dx |j> of sinc functions

  !write(6,*) "solve_sinc_dvr: nstates=", nstates
 
  LWORK = 3*nstates_na
  
  !
  allocate(Mat_tmp(nstates_na,nstates_na), W(nstates_na),&
       WORK(LWORK))
  
  ! compute 


  ! set up hamiltonan matrix
  Mat_tmp = 0.0_wp
  do i1=1, npesfiles_f
     do j1=1, nstates
        
        ii = (i1 -1) * nstates + j1  

        do i2=1, npesfiles_f
           do j2=1, nstates

              jj = (i2 -1) * nstates + j2  


              ! here should be a term + 2 * <g_G|<G| d/dR |F> d/dR |f_F> =  2 * sum_l <g_G|<G| d/dR |F> |l><l|  d/dR |f_F>
              !  2 * sum_k sum_m c^{g_G}_k c^{f_F}_m <k| <G| d/dR |F> |k><k|  d/dR |m> 
              !m_elem2=0.0_wp
              !do k=1,nstates
              !   do m=1, nstates
              !   m_elem2 = m_elem2 +  c_f(i1, k, j1) * c_f(i2, m , j2) * nac(i1,i2,k,1) * M_dsinc(k,m) 
              !  end do
              !end do



              Mat_tmp(ii,jj) = eig_f(i1, j1) * delta(i1,i2) * delta(j1,j2) + &
                   sum(c_f(i1, :, j1) * c_f(i2, : , j2) * nac(i1,i2,:,2) ) + &
                   sum(c_f(i1, :, j1) * nac(i1,i2,:,1) *  matmul( M_dsinc, c_f(i2, : , j2)) )


           end do ! j2
        end do ! i2
     end do ! j1
  end do ! i1
  
  write(6,*) "done setting up Hamiltonian"

  ! solve eigenvalue problem 
  call dsyev( "V", "U", nstates_na, Mat_tmp, nstates_na, W, WORK, LWORK, INFO)
  
  ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
  ! put c_na in a good order

  do i1=1, npesfiles_f
     do j1=1, nstates
        ii = (i1 -1) * nstates + j1  

        c_na(i1, j1, :) = Mat_tmp(ii, :)
     end do
  end do

  eig_na = W

end subroutine solve_non_adiabatic


subroutine calculate_dipoles_XAS(c_i, c_f, dipole, D_fi, mode, ind_in)
  use m_precision, only: wp

  real(kind=wp), intent(in):: c_i(:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out)::  D_fi(:,:,:)
  character(*), intent(in):: mode
  integer, intent(in), optional :: ind_in

  real(wp), allocatable:: factor(:,:,:)
  integer:: nstates, npesfile_f,i,j,l
  integer:: ind

  write(6,*) "calculate_dipoles_XAS: mode ", mode
  
  if (present(ind_in)) then
    ind = ind_in
    write(6,*) "ind", ind
  else
    ind  = 1
  end if
  
  nstates = size(c_i,2)
  npesfile_f = size(c_f,1)

  allocate(factor(npesfile_f, nstates, 3))

  !
  ! calculate dipole matrix elements between states 
  !

  if (mode .eq. "DIPOLE") then
    factor = dipole
  else if(mode .eq. "FC") then
    factor = 1.0_wp
  else if(mode .eq. "DIPOLE_X0") then
    
    write(6,*) "calculate_dipoles_XAS: DIPOLE_X0"

    do i=1, nstates
      factor(:,i,:) = dipole(:, ind, :)
    end do
    
  else
    
    write(6,*) "calculate_dipoles_XAS: mode must be either DIPOLE, FC, DIPOLE_X0"
    
  end if
    
  ! transitions from ground state to final states
  do i=1,npesfile_f
    do j=1,nstates ! final
      do l=1,3
        D_fi(i,j,l) = sum(factor(i,:,l) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        !D_fi(i,j,l) = sum(abs(dipole(i,:,l)) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        !D_fi(i,j,l) = sum(c_f(i,:,j) * c_i(:,1))   ! FC profile 
      end do
    end do
  end do
  
    
  write(6,*) "Calculated dipole matrix elements"

end subroutine calculate_dipoles_XAS


subroutine spectrum_XAS(eig_f, eig_i, D_fi, omega, sigma, sigma_states, gamma)
  use m_precision, only: wp
  
  real(kind=wp), intent(in):: eig_f(:,:), eig_i(:), D_fi(:,:,:), omega(:)
  real(kind=wp), intent(out):: sigma(:), sigma_states(:,:)
  real(kind=wp),intent(in):: gamma

  integer:: nstates, npesfile_f, n_omega, i,j,k,l,m
  complex(kind=wp),allocatable:: F(:,:)
  real(kind=wp):: norm

  nstates = size(eig_f,2)
  npesfile_f = size(eig_f,1)
  n_omega = size(omega)

  write(6,*) "n_omega", n_omega

  allocate(F(n_omega, 3))

  !
  ! Fermi golden rule 
  !


  sigma = 0.0_wp
  sigma_states = 0.0_wp

  do j= 1,npesfile_f ! F
     do l=1,nstates ! f_F

        F=0.0_wp        
        do i =1, n_omega 
           do m=1,3 ! polarization
              
              F(i,m) = F(i,m) + D_fi(j,l ,m)   / ( omega(i) - &
                   (eig_f(j,l) -eig_i(1)) + dcmplx(0,gamma) )

           end do
        end do
        
        !square F
        sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
             + real(conjg(F(:,3))*F(:,3)) 

     end do ! l
     
     sigma = sigma  + sigma_states(j,:)
     
  end do !j
  
  write(6,*) "Calculated XAS spectrum"

  ! normalize spectrum
  norm=sum(sigma) *(omega(2) -omega(1)) 
  sigma = sigma/norm
  
  do j= 1, npesfile_f
     sigma_states(j,:) = sigma_states(j,:)/norm
  end do

  deallocate(F)

end subroutine spectrum_XAS

subroutine compute_XAS_spectrum(eig_i, eig_n, D_ni, omega_in, gamma, gamma_inc, sigma_final)
  use m_precision, only: wp
  
  real(wp), intent(in):: eig_i, eig_n(:), D_ni(:,:), omega_in(:)
  real(wp),intent(in):: gamma, gamma_inc
  real(wp), intent(out):: sigma_final(:,:,:)

  complex(wp),allocatable:: F(:,:)
  real(wp):: norm
  integer:: nstates, n_omega_in, i,j,k,l,m1, m2
  real(wp), allocatable:: sigma_tmp(:)
  
  nstates = size(eig_n)
  n_omega_in = size(omega_in)

  allocate(F(n_omega_in, 3))

  sigma_final= 0.0_wp
  F=0.0_wp

  do l=1,nstates ! f_F

    F=0.0_wp
    do m1=1,3
      F(:,m1) = D_ni(l ,m1)   / ( omega_in(:) - &
           (eig_n(l) -eig_i) + dcmplx(0.0_wp,gamma) )
    end do

    do m1=1,3 
      do m2=1,3 
        sigma_final(:,m1,m2) = sigma_final(:,m1,m2) + real(conjg(F(:,m1))*F(:,m2))
      end do
    end do
    
  end do ! l

  ! convolute
  if (gamma_inc .gt. 0.001_wp) then
    allocate(sigma_tmp(n_omega_in))

    write(6,*) "convoluting XAS spectrum with HWHM", gamma_inc, "eV"
    
    do m1=1,3 
      do m2=1,3 
        sigma_tmp = sigma_final(:,m1,m2) 
        call convolute_XAS_incoming(sigma_tmp, omega_in, gamma_inc, sigma_final(:,m1,m2))
        !sigma_final(:,m1,m2) = sigma_tmp
      end do
    end do

  end if
  
end subroutine compute_XAS_spectrum


subroutine spectrum_XAS_nonadiabatic(eig_na, eig_i, c_na, D_fi, omega, sigma, sigma_states, gamma)
  use m_precision, only: wp
  
  real(kind=wp), intent(in):: eig_na(:), eig_i(:), c_na(:,:,:), D_fi(:,:,:), omega(:)
  real(kind=wp), intent(out):: sigma(:), sigma_states(:,:)
  real(kind=wp),intent(in):: gamma


  integer:: nstates, nstates_na, npesfile_f, n_omega, i,j,k,l,m
  complex(kind=wp),allocatable:: F(:,:)
  real(kind=wp):: norm

  nstates = size(eig_i)
  nstates_na = size(eig_na)
  npesfile_f = size(D_fi,1)
  n_omega = size(omega)

  allocate(F(n_omega, 3))


  !
  ! Fermi golden rule 
  !

  write(6,*) "Computing nonadiabatic XAS spectrum", npesfile_f, nstates, nstates_na

  sigma=0.0_wp
  
  do j= 1,nstates_na ! F'
     
      do k= 1, npesfile_f ! F
        do l=1, nstates ! f_F

           F=0.0_wp         
           do i =1, n_omega 
              do m=1,3 ! polarization
                 
                 F(i,m) = F(i,m) + c_na(k,l,j) *  D_fi(k, l ,m)  / ( omega(i) - &
                      (eig_na(j) -eig_i(1)) + dcmplx(0,gamma) )
                 

              end do !i
           end do ! m

           !square F
           ! the summation is put here to avoid possible interference effects, not sure if they should be there or not
           sigma = sigma  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
                + real(conjg(F(:,3))*F(:,3)) 

           
        end do ! l
     end do ! k

  end do !j
  
  write(6,*) "Calculated spectrum"

  ! normalize spectrum
  norm=sum(sigma) *(omega(2) -omega(1)) 
  sigma = sigma/norm

  deallocate(F)
 
end subroutine spectrum_XAS_nonadiabatic

! convolute with respect to omega_out
subroutine convolute_XAS_incoming(sigma, omega_in, gamma_inc, sigma_out)
  use m_precision, only: wp
  use m_KH_functions, only: gaussian
  
  real(wp), intent(in):: sigma(:)
  real(wp), intent(in):: omega_in(:)
  real(wp), intent(in):: gamma_inc
  real(wp), intent(out):: sigma_out(:)

  integer:: om_in1, om_in2
  
  sigma_out =0.0_wp
  
  do om_in1 = 1, size(omega_in)
    do om_in2 = 1, size(omega_in)
      sigma_out(om_in1) = sigma_out(om_in1) + &
           sigma(om_in2) * gaussian(omega_in(om_in1)- omega_in(om_in2), &
           0.0_wp, 2.0_wp * gamma_inc)
    end do
  end do

  sigma_out = sigma_out * (omega_in(2)- omega_in(1))
  
end subroutine convolute_XAS_incoming


end module m_XAS_functions


