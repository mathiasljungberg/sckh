module m_XAS_functions
  use parameters
  use m_XAS_io
  use KH_functions
  use spline_m
  !use qsort_c_module

  implicit none

contains

subroutine solve_non_adiabatic(eig_f, c_f, nac, eig_na, c_na, M_dsinc ) 
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


subroutine calculate_dipoles_XAS(c_i, c_f, dipole, D_fi)
  real(kind=wp), intent(in):: c_i(:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out)::  D_fi(:,:,:)

  integer:: nstates, npesfile_f,i,j,l

  nstates = size(c_i,2)
  npesfile_f = size(c_f,1)

  !
  ! calculate dipole matrix elements between states 
  !

  ! transitions from ground state to final states
  do i=1,npesfile_f
     do j=1,nstates ! final
        do l=1,3
           D_fi(i,j,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        end do
     end do
  end do
  
  write(6,*) "Calculated dipole matrix elements"

end subroutine calculate_dipoles_XAS


subroutine spectrum_XAS(eig_f, eig_i, D_fi, omega, sigma, sigma_states, gamma)
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

subroutine spectrum_XAS_nonadiabatic(eig_na, eig_i, c_na, D_fi, omega, sigma, sigma_states, gamma)
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

subroutine calculate_dipoles_XES(c_i, c_n, c_f, dipole, D_ni, D_fn, D_fi)
  real(kind=wp), intent(in):: c_i(:,:), c_n(:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out):: D_ni(:), D_fn(:,:,:,:), D_fi(:,:,:)

  integer:: nstates, npesfile_f,i,j,k,l

  nstates = size(c_i,2)
  npesfile_f = size(c_f,1)

  !
  ! calculate dipole matrix elements between states 
  !

  do i=1,nstates 
     D_ni(i) = sum(c_n(:,i)*c_i(:,1))
  end do
     
  do i=1, npesfile_f
     do j=1, nstates ! final
        do k=1, nstates ! intermediate
           do l=1,3
              D_fn(i,j,k,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
              !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
              !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))                  ! no dipole moment
           end do
        end do
     end do
  end  do

  ! transitions from ground state directly to final states
  do i=1, npesfile_f
     do j=1, nstates ! final
        do l=1,3
           D_fi(i,j,l) = sum(dipole(i,:,l) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        end do
     end do
  end do

  write(6,*) "Calculated XES dipole matrix elements"

end subroutine calculate_dipoles_XES


subroutine spectrum_XES(eig_f, eig_n,  D_ni, D_fn, D_fi, omega, sigma, sigma_states, gamma)
  real(kind=wp), intent(in):: eig_f(:,:), eig_n(:), D_ni(:), D_fn(:,:,:,:), D_fi(:,:,:), omega(:)
  real(kind=wp), intent(out):: sigma(:), sigma_states(:,:)
  real(kind=wp), intent(in):: gamma

  complex(kind=wp), dimension(:,:),allocatable:: F
  integer:: npesfile_f, nstates, n_omega
  integer:: i,j,k,l,m
  real(kind=wp):: norm


  npesfile_f = size(eig_f,1)
  nstates = size(eig_f,2)
  n_omega = size(omega)


  allocate(F(n_omega,3))


  !
  ! Full Kramers-Heisenberg    
  !


  write(6,*) "Full Kramers-Heisenberg"

  sigma=0
  !sigma_dir=0
  !sigma_max_int = 0

  do j= 1,npesfile_f ! final el
     sigma_states(j,:) = 0
     !sigma_dir_states(j,:) = 0
     !sigma_max_int_states(j,:) = 0
     do l=1,nstates ! final vib

        F=0.0_wp
        !F2_dir=0
        !F2_max_int =0

        do i =1,  n_omega 
           do k=1,nstates ! intermediate vib
              do m=1,3 ! polarization

                 F(i,m) = F(i,m) + D_fn(j,l,k,m) * D_ni(k) / ( omega(i) - &
                      (eig_n(k) - eig_f(j,l)) + dcmplx(0,gamma) )

                 ! direct contribution (no interference)
                 !F2_dir(i,m) = F2_dir(i,m) + D_fn(j,l,k,m) ** 2 * D_ni(k) ** 2 / (( omega(i) - &
                 !     (eig_n(k) - eig_f(j,l)))**2  + gamma**2 ) 

              end do
           end do

           ! maximal interference, direct transition from initial to final states
           !do m=1,3 ! polarization 
!
!              F2_max_int(i,m) =  F2_max_int(i,m) + D_fi(j,l,m) ** 2 / (( omega(i) - &
!                   (E_n_mean - eig_f(j,l)) )**2  + gamma**2 )
!           end do

        end do

        !square F
        sigma_states(j,:) = sigma_states(j,:)  + real(conjg(F(:,1))*F(:,1)) + real(conjg(F(:,2))*F(:,2)) &
             + real(conjg(F(:,3))*F(:,3)) 

        !sigma_dir_states(j,:) = sigma_dir_states(j,:) + F2_dir(:,1) +F2_dir(:,2) +F2_dir(:,3)
        !sigma_max_int_states(j,:) = sigma_max_int_states(j,:) + F2_max_int(:,1) + F2_max_int(:,2) + F2_max_int(:,3)

     end do !l

     sigma = sigma + sigma_states(j,:)
     !sigma_dir = sigma_dir + sigma_dir_states(j,:)
     !sigma_max_int = sigma_max_int + sigma_max_int_states(j,:)

  end do ! j

  write(6,*) "Calculated XES spectrum"

  ! normalize spectrum
  norm=sum(sigma) *(omega(2) -omega(1)) 
  sigma = sigma/norm
  !sigma_dir = sigma_dir/norm
  !sigma_max_int = sigma_max_int/norm

  !sigma_lp = sigma_lp/norm

  do j= 1,npesfile_f
     sigma_states(j,:) = sigma_states(j,:)/norm
     !sigma_dir_states(j,:) = sigma_dir_states(j,:)/norm
     !sigma_max_int_states(j,:) = sigma_max_int_states(j,:)/norm
  end do

  deallocate(F)

end subroutine spectrum_XES

!subroutine print_output()
!end subroutine print_output

end module m_XAS_functions


