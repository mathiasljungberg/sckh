module m_KH_utils
  implicit none

#include "m_define_macro.F90"
contains

subroutine solve_non_adiabatic(eig_f, c_f, nac, eig_na, c_na) 
  use m_precision, only: wp
  use m_KH_functions, only: delta

  real(kind=wp), intent(in):: eig_f(:,:), c_f(:,:,:), nac(:,:,:,:)
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


  !write(6,*) "solve_sinc_dvr: nstates=", nstates
 
  LWORK = 3*nstates_na
  
  !
  allocate(Mat_tmp(nstates_na,nstates_na), W(nstates_na),&
       WORK(LWORK))
  
  ! set up hamiltonan matrix
  Mat_tmp = 0.0_wp
  do i1=1, npesfiles_f
     do j1=1, nstates
        
        ii = (i1 -1) * nstates + j1  

        do i2=1, npesfiles_f
           do j2=1, nstates

              jj = (i2 -1) * nstates + j2  

              Mat_tmp(ii,jj) = eig_f(i1, j1) * delta(i1,i2) * delta(j1,j2) + &
                   sum(c_f(i1, :, j1) * c_f(i2, : , j2) * nac(i1,i2,:,2) )

              ! here should be a term + 2 * <g_G|<G| d/dR |F> d/dR |f_F> (cannot yet compute it) 

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
  use m_precision, only: wp

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

subroutine calculate_dipoles_KH_nonres(c_i, c_n, c_f, dipole, D_ni, D_fn, D_fi)
  use m_precision, only: wp

  real(kind=wp), intent(in):: c_i(:,:), c_n(:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out):: D_ni(:,:), D_fn(:,:,:,:), D_fi(:,:,:)

  integer:: nstates, npesfile_f,i,j,k,m

  nstates = size(c_i,2)
  npesfile_f = size(c_f,1)

  !
  ! calculate dipole matrix elements between states 
  !

  do i=1,nstates
    do m=1,3
      D_ni(i,m) = sum(c_n(:,i)*c_i(:,1))
    end do
  end do
     
  do i=1, npesfile_f
     do j=1, nstates ! final
        do k=1, nstates ! intermediate
           do m=1,3
              D_fn(i,j,k,m) = sum(dipole(i,:,m) * c_f(i,:,j) * c_n(:,k))   ! true dipole moment
              !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
              !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))                  ! no dipole moment
           end do
        end do
     end do
  end  do

  ! transitions from ground state directly to final states
  do i=1, npesfile_f
     do j=1, nstates ! final
        do m=1,3
           D_fi(i,j,m) = sum(dipole(i,:,m) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
        end do
     end do
  end do

  write(6,*) "Calculated XES dipole matrix elements"

end subroutine calculate_dipoles_KH_nonres

subroutine calculate_dipoles_KH_res(c_i, c_n, c_f, dipole, D_ni, D_fn, D_fi)
  use m_precision, only: wp

  real(kind=wp), intent(in):: c_i(:,:), c_n(:,:,:), c_f(:,:,:), dipole(:,:,:)
  real(kind=wp), intent(out):: D_ni(:,:,:), D_fn(:,:,:,:,:), D_fi(:,:,:)

  integer:: nstates, npesfile_f, npesfile_n, f_e, f_v, n_e, n_v, m

  nstates = size(c_i,2)
  npesfile_n = size(c_n,1)
  npesfile_f = size(c_f,1)
  

  !
  ! calculate dipole matrix elements between states 
  !

  do n_e =1, npesfile_n
    do n_v =1,nstates 
      do m=1,3
        D_ni(n_e, n_v, m) = sum(c_n(n_e,:,n_v)*c_i(:,1))
      end do
    end do
  end do

  
  do f_e = 1, npesfile_f
    do f_v = 1, nstates ! final
      
      do n_e = 1, npesfile_n ! intermediate
        do n_v = 1, nstates ! intermediate
          
          do m=1,3
            !D_fn(f_e, f_v, n_e, n_v, m) = sum(dipole(f_e, n_e, :, m) * c_f(f_e, :, f_v) * c_n(n_e, :, n_v))   ! true dipole moment
            D_fn(f_e, f_v, n_e, n_v, m) = sum(c_f(f_e, :, f_v) * c_n(n_e, :, n_v))  ! FC, dipole moment at eq geom
            !D_fn(i,j,k,l) = dipole(i,21,l) * sum(c_f(i,:,j) * c_n(:,k)) ! FC, dipole moment at eq geom
            !D_fn(i,j,k,l) = sum(c_f(i,:,j) * c_n(:,k))                  ! no dipole moment
          end do
          
        end do
      end do
      
    end  do
  end do
  
  ! transitions from ground state directly to final states
  do f_e = 1, npesfile_f
     do f_v = 1, nstates ! final
       do m=1,3
         !D_fi(f_e, f_v ,m) = sum(dipole(f_e,:,m) * c_f(i,:,j) * c_i(:,1))   ! true dipole moment
         D_fi(f_e, f_v ,m) = sum(c_f(f_e, :, f_v) * c_i(:,1))   ! FC
       end do
     end do
   end do

  write(6,*) "Calculated XES dipole matrix elements"

end subroutine calculate_dipoles_KH_res




subroutine spectrum_XES(eig_f, eig_n,  D_ni, D_fn, D_fi, omega, sigma, sigma_states, gamma)
  use m_precision, only: wp

  real(kind=wp), intent(in):: eig_f(:,:), eig_n(:), D_ni(:,:), D_fn(:,:,:,:), D_fi(:,:,:), omega(:)
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

                 F(i,m) = F(i,m) + D_fn(j,l,k,m) * D_ni(k,1) / ( omega(i) - &
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



!
! These are new routines for general XES spectrum calcultions with resonant and non-resonant on the same footing
!

!
! n_vib_f: is in the case of the final state index being a composite index of electronic and vibrational states
!          f_tot = (f_el-1) * n_vib_f + f_vib
!          the vibrational subindices will be summed.
!          In case n_vib_f = size(E_f) then all states will be summed, if n_vib_f =1, no states will be summed
!

subroutine compute_XES_res(E_i, E_n, E_f, D_ni, D_fn, omega_in, omega_out, &
     gamma, gamma_inc, gamma_instr, flag_res, flag_nonres,&
     sigma_final, lambda_F, lambda_G, lambda_H)
  use m_precision, only: wp
  use m_KH_functions, only: gaussian
  
  real(kind=wp), intent(in):: E_i, E_n(:), E_f(:)
  !integer, intent(in):: n_el_f   
  real(kind=wp), intent(in):: D_ni(:,:), D_fn(:,:,:)
  real(kind=wp), intent(in):: omega_in(:)
  real(kind=wp), intent(in):: omega_out(:)
  real(kind=wp), intent(in):: gamma
  real(kind=wp), intent(in):: gamma_inc
  real(kind=wp), intent(in):: gamma_instr
  logical,intent(in):: flag_res, flag_nonres
  real(kind=wp), intent(out):: sigma_final(:,:,:,:)
  real(kind=wp), intent(out):: lambda_F(:,:)
  real(kind=wp), intent(out):: lambda_G(:,:)
  real(kind=wp), intent(out):: lambda_H(:,:)

  integer:: om_in, om_out, i_f, m1, m2 !, i_vib_f, i_el_f, n_vib_f
  complex(wp):: F(3,3)
  real(wp):: prefac, broadening
  real(wp), allocatable:: sigma_tmp(:,:)
  
  write(6,*) "compute_XES_res 0"

  sigma_final = 0.0_wp
  lambda_F =0 
  lambda_G =0
  lambda_H =0
  
  do om_out = 1, size(omega_out)
      write(6,*) "om_out", om_out

      do i_f = 1, size(E_f)

        call compute_amplitude_F(E_i, E_n, E_f(i_f), &
             D_ni(:,:), D_fn(i_f,:,:), omega_out(om_out), gamma, F(:,:), &
             flag_res, flag_nonres)

        do om_in = 1, size(omega_in)
        
          ! broadening, intitial distribution, including prefactor
          ! broadening = (omega_in(om_in) / omega_out(om_out)) *  &
          ! gamma_instr / ((omega_in(om_in) - omega_out(om_out) - (E_f(i_f)- E_i))**2 + gamma_instr ** 2 )

          ! alt broadening, instrumental, including prefactor
          ! prefac = omega_out(om_iout) / omega_in(om_in)
          prefac = omega_out(om_out) / (omega_out(om_out) + (E_f(i_f)- E_i))
          
          broadening = prefac * gaussian(omega_in(om_in) - omega_out(om_out) - (E_f(i_f)- E_i), 0.0_wp, 2.0_wp *  gamma_inc)
          !gamma_instr / ((omega_in(om_in) - omega_out(om_out) - (E_f(i_f)- E_i))**2 + gamma_instr ** 2 )
          
          sigma_final(om_in, om_out,:,:) =  sigma_final(om_in, om_out,:,:) + &
               abs(F)**2 * broadening
          
          ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
          do m1=1,3
            do m2=1,3
              lambda_F(om_in, om_out) = lambda_F(om_in, om_out) +  real(F(m1,m1) * conjg(F(m2,m2)) ) * broadening
              lambda_G(om_in, om_out) = lambda_G(om_in, om_out) +  real(F(m1,m2) * conjg(F(m1,m2)) ) * broadening
              lambda_H(om_in, om_out) = lambda_H(om_in, om_out) +  real(F(m1,m2) * conjg(F(m2,m1)) ) * broadening
            end do
          end do
          
        end do
        
      end do
    end do

    ! convolute with instrumental broadening
    allocate(sigma_tmp(size(omega_in), size(omega_out)))

    do m1=1,3
      do m2=1,3
        sigma_tmp = sigma_final(:,:,m1,m2)
        call convolute_intrumental(sigma_tmp, omega_in, omega_out, gamma_instr, sigma_final(:,:,m1,m2))
      end do
    end do

    sigma_tmp = lambda_F
    call convolute_intrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_F)
    sigma_tmp = lambda_G
    call convolute_intrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_G)
    sigma_tmp = lambda_H
    call convolute_intrumental(sigma_tmp, omega_in, omega_out, gamma_instr, lambda_H)

    
  end subroutine compute_XES_res

  ! here use F(\omega) instead of F(\omega'). Should be exactly the same, but we can save time if N_{\omega} < N_{\omega'}
  subroutine compute_XES_res_alt(E_i, E_n, E_f, D_ni, D_fn, omega_in, omega_out, &
       gamma, gamma_inc, gamma_instr, flag_res, flag_nonres,&
       sigma_final, lambda_F, lambda_G, lambda_H)
    use m_precision, only: wp
    use m_KH_functions, only: gaussian
    
    real(kind=wp), intent(in):: E_i, E_n(:), E_f(:)
    !integer, intent(in):: n_el_f   
    real(kind=wp), intent(in):: D_ni(:,:), D_fn(:,:,:)
    real(kind=wp), intent(in):: omega_in(:)
    real(kind=wp), intent(in):: omega_out(:)
    real(kind=wp), intent(in):: gamma
    real(kind=wp), intent(in):: gamma_inc
    real(kind=wp), intent(in):: gamma_instr
    logical,intent(in):: flag_res, flag_nonres
    real(kind=wp), intent(out):: sigma_final(:,:,:,:)
    real(kind=wp), intent(out):: lambda_F(:,:)
    real(kind=wp), intent(out):: lambda_G(:,:)
    real(kind=wp), intent(out):: lambda_H(:,:)

    integer:: om_in, om_out, i_f, m1, m2 !, i_vib_f, i_el_f, n_vib_f
    complex(wp):: F(3,3)
    real(wp):: prefac, broadening
    real(wp), allocatable:: sigma_tmp(:,:)
    
    write(6,*) "compute_XES_res 0"
    
    sigma_final = 0.0_wp
    lambda_F =0 
    lambda_G =0
    lambda_H =0
    
    do om_in = 1, size(omega_in)
      write(6,*) "om_in", om_in
      
      do i_f = 1, size(E_f)
        
        ! switcha E_f och E_i, bor bli ekvivalent
        call compute_amplitude_F(E_f(i_f), E_n, E_i, &
             D_ni(:,:), D_fn(i_f,:,:), omega_in(om_in), gamma, F(:,:), &
             flag_res, flag_nonres)

        do om_out = 1, size(omega_out)
        
          ! broadening, intitial distribution, including prefactor
          ! broadening = (omega_in(om_in) / omega_out(om_out)) *  &
          ! gamma_instr / ((omega_in(om_in) - omega_out(om_out) - (E_f(i_f)- E_i))**2 + gamma_instr ** 2 )

          ! alt broadening, instrumental, including prefactor
          ! prefac = omega_out(om_iout) / omega_in(om_in)
          !prefac = omega_out(om_out) / (omega_out(om_out) + (E_f(i_f)- E_i))
          prefac = (omega_in(om_in) - (E_f(i_f)- E_i) ) /  omega_in(om_in) 
          
          broadening = prefac * gaussian(omega_out(om_out) - omega_in(om_in) + (E_f(i_f)- E_i), 0.0_wp, 2.0_wp *  gamma_instr)
          !gamma_instr / ((omega_in(om_in) - omega_out(om_out) - (E_f(i_f)- E_i))**2 + gamma_instr ** 2 )
          
          sigma_final(om_in, om_out,:,:) =  sigma_final(om_in, om_out,:,:) + &
               abs(F)**2 * broadening
          
          ! perform spherical average according to J. Phys. B. 27, 4169 (1994)
          do m1=1,3
            do m2=1,3
              lambda_F(om_in, om_out) = lambda_F(om_in, om_out) +  real(F(m1,m1) * conjg(F(m2,m2)) ) * broadening
              lambda_G(om_in, om_out) = lambda_G(om_in, om_out) +  real(F(m1,m2) * conjg(F(m1,m2)) ) * broadening
              lambda_H(om_in, om_out) = lambda_H(om_in, om_out) +  real(F(m1,m2) * conjg(F(m2,m1)) ) * broadening
            end do
          end do
          
        end do
        
      end do
    end do

    ! convolute with incoming broadening
    allocate(sigma_tmp(size(omega_in), size(omega_out)))

    do m1=1,3
      do m2=1,3
        sigma_tmp = sigma_final(:,:,m1,m2)
        call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, sigma_final(:,:,m1,m2))
      end do
    end do

    sigma_tmp = lambda_F
    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_F)
    sigma_tmp = lambda_G
    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_G)
    sigma_tmp = lambda_H
    call convolute_incoming(sigma_tmp, omega_in, omega_out, gamma_inc, lambda_H)

    
  end subroutine compute_XES_res_alt

  

subroutine compute_XES_nonres_elec(E_i, E_n, E_f,  D_ni, D_fn, omega_out, gamma, sigma_m)
  use m_precision, only: wp
  
  real(kind=wp), intent(in):: E_i, E_n, E_f(:)
  real(kind=wp), intent(in):: D_ni(:), D_fn(:,:)
  real(kind=wp), intent(in):: omega_out(:)
  real(kind=wp), intent(in):: gamma
  complex(kind=wp), intent(out):: sigma_m(:,:,:)

  real(kind=wp), allocatable:: E_n_tmp(:), E_f_tmp(:)
  real(kind=wp), allocatable:: D_ni_tmp(:,:), D_fn_tmp(:,:,:)
  real(kind=wp), allocatable:: sigma_tmp(:,:,:)
  integer:: m, i_f
  
  allocate(E_n_tmp(1))
  allocate(E_f_tmp(1))
  allocate(D_ni_tmp(1,size(D_ni)))
  allocate(D_fn_tmp(1,size(D_fn,1), size(D_fn,2)))
  allocate(sigma_tmp(size(omega_out), 3, 3))
  
  E_n_tmp(1) = E_n
  !E_f_tmp(:) = E_f(:)
  D_ni_tmp(1,:) = D_ni(:)
  D_fn_tmp(1,:,:) = D_fn(:,:)
  
  do i_f =1, size(E_f_tmp)
    E_f_tmp(1) = E_f(i_f)
    call compute_XES_nonres(E_i, E_n_tmp, E_f_tmp, D_ni_tmp, D_fn_tmp, omega_out, gamma, sigma_tmp)

    do m=1,3
      sigma_m(i_f,:,m) = sigma_tmp(:,m,m)
    end do
  end do
end subroutine compute_XES_nonres_elec

subroutine compute_XES_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_out, gamma, sigma_pol)
  use m_precision, only: wp

  real(kind=wp), intent(in):: E_i
  real(kind=wp), intent(in):: E_n(:), E_f(:)
 ! integer, intent(in):: n_el_f   
  real(kind=wp), intent(in):: D_ni(:,:), D_fn(:,:,:)
  real(kind=wp), intent(in):: omega_out(:)
  real(kind=wp), intent(in):: gamma
  real(kind=wp), intent(out):: sigma_pol(:,:,:)

  integer:: i_f, om_out !, i_vib_f, i_el_f, n_vib_f
  complex(wp):: F(3,3)

!  n_vib_f = size(E_f) / n_el_f

  sigma_pol = 0.0_wp
  do om_out=1, size(omega_out)
    
    !do i_el_f = 1, n_el_f !size(E_f)
    !  sigma_final(i_el_f, om_out,:,:) = 0.0_wp
    !  do i_vib_f = 1, n_vib_f
    !    i_f = (i_el_f-1)*n_vib_f + i_vib_f
    
    do i_f =1, size(E_f)
      call compute_amplitude_F(E_i, E_n, E_f(i_f),  D_ni, D_fn(i_f,:,:), &
           omega_out(om_out), gamma, F(:,:), .true., .true.)      
      
      sigma_pol(om_out, :,:) = sigma_pol(om_out, :,:) + abs(F(:,:))**2
      
    end do
    
  end do
    
end subroutine compute_XES_nonres

! scattering amplitude, for each frequency, and final state
! F_f(\omega') = \sum_{n} \frac{D_{in} D_{nf} } {\omega' -(E_n -E_f) + i\gamma }
subroutine compute_amplitude_F(E_i, E_n, E_f, D_ni, D_fn, omega_out, gamma, F, flag_res, flag_nonres)
  use m_precision, only: wp

  real(kind=wp), intent(in):: E_i
  real(kind=wp), intent(in):: E_n(:), E_f
  real(kind=wp), intent(in):: D_ni(:,:), D_fn(:,:)
  real(kind=wp), intent(in):: omega_out
  real(kind=wp), intent(in):: gamma
  logical, intent(in):: flag_res, flag_nonres
  complex(kind=wp), intent(out):: F(:,:)

  integer:: m1, m2
  complex(wp):: F_tmp(3,3)

  F=0.0
  
  if (flag_res) then
    call compute_amplitude_F_res(E_i, E_n, E_f, D_ni, D_fn, omega_out, gamma, F_tmp)    
    F = F + F_tmp
  end if

  if (flag_nonres) then
    call compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_out, gamma, F_tmp)    
    F = F + F_tmp
  end if
  
end subroutine compute_amplitude_F

! scattering amplitude, for each frequency, and final state
! F_f(\omega') = \sum_{n} \frac{D_{in} D_{nf} } {\omega' -(E_n -E_f) + i\gamma }
subroutine compute_amplitude_F_res(E_i, E_n, E_f, D_ni, D_fn, omega_out, gamma, F)
  use m_precision, only: wp
  
  real(kind=wp), intent(in):: E_i
  real(kind=wp), intent(in):: E_n(:), E_f
  real(kind=wp), intent(in):: D_ni(:,:), D_fn(:,:)
  real(kind=wp), intent(in):: omega_out
  real(kind=wp), intent(in):: gamma
  complex(kind=wp), intent(out):: F(:,:)

  integer:: m1, m2

  F =0.0_wp
  
  do m1=1,3
    do m2=1,3
      F(m1,m2) = sum(D_ni(:,m2)*D_fn(:,m1) / &
           (omega_out -(E_n(:)-E_f) + cmplx(0, gamma,8)))
    end do
  end do
  
end subroutine compute_amplitude_F_res


! scattering amplitude, for each frequency, and final state
! F_f(\omega') = \sum_{n} \frac{D_{in} D_{nf} } {\omega' -(E_n -E_f) + i\gamma }
subroutine compute_amplitude_F_nonres(E_i, E_n, E_f, D_ni, D_fn, omega_out, gamma, F)
  use m_precision, only: wp

  real(kind=wp), intent(in):: E_i
  real(kind=wp), intent(in):: E_n(:), E_f
  real(kind=wp), intent(in):: D_ni(:,:), D_fn(:,:)
  real(kind=wp), intent(in):: omega_out
  real(kind=wp), intent(in):: gamma
  complex(kind=wp), intent(out):: F(:,:)

  integer:: m1, m2

  F =0.0_wp

  do m1=1,3
    do m2=1,3
      F(m1,m2) = -sum(D_ni(:,m2)*D_fn(:,m1) /&
           (omega_out + (E_n(:)-E_i)))
    end do
  end do

end subroutine compute_amplitude_F_nonres

! sigma(om_in, om_out)
! convolute with respect to omega_out
subroutine convolute_intrumental(sigma, omega_in, omega_out, gamma_instr, sigma_out)
  use m_precision, only: wp
  use m_KH_functions, only: gaussian
  
  real(wp), intent(in):: sigma(:,:)
  real(wp), intent(in):: omega_in(:)
  real(wp), intent(in):: omega_out(:)
  real(wp), intent(in):: gamma_instr
  real(wp), intent(out):: sigma_out(:,:)

  integer:: om_in1, om_out1, om_out2
  
  sigma_out =0.0_wp
  
  !do om_in1 = 1, size(omega_in)
  do om_out1 = 1, size(omega_out)
    do om_out2 = 1, size(omega_out) 
      sigma_out(:, om_out1) = sigma_out(:, om_out1) + sigma(:, om_out2) * &
           gaussian(omega_out(om_out1)- omega_out(om_out2), &
           0.0_wp, 2.0_wp * gamma_instr)
    end do
  end do
  !end do
  
  sigma_out = sigma_out * (omega_out(2)- omega_out(1))
  
end subroutine convolute_intrumental

! sigma(om_in, om_out)
! convolute with respect to omega_out
subroutine convolute_incoming(sigma, omega_in, omega_out, gamma_inc, sigma_out)
  use m_precision, only: wp
  use m_KH_functions, only: gaussian
  
  real(wp), intent(in):: sigma(:,:)
  real(wp), intent(in):: omega_in(:)
  real(wp), intent(in):: omega_out(:)
  real(wp), intent(in):: gamma_inc
  real(wp), intent(out):: sigma_out(:,:)

  integer:: om_in1, om_in2, om_out1
  
  sigma_out =0.0_wp
  
  !do om_out1 = 1, size(omega_out)
  do om_in1 = 1, size(omega_in)
    do om_in2 = 1, size(omega_in)
      sigma_out(om_in1, :) = sigma_out(om_in1, :) + &
           sigma(om_in2, :) * gaussian(omega_in(om_in1)- omega_in(om_in2), &
           0.0_wp, 2.0_wp * gamma_inc)
    end do
  end do
  !end do

  sigma_out = sigma_out * (omega_in(2)- omega_in(1))
  
end subroutine convolute_incoming

  
end module m_KH_utils


