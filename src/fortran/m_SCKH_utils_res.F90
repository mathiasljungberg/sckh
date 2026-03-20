module m_SCKH_utils_res

  implicit none
  
contains

  ! basic routines for resonant SCKH

  ! computes
  ! F_om(om, m1, m2) = \int_0^{t_max} dt D1(m2) D2(m1,t) e_factor(t) exp(-gamma t) exp(i om t)   
  subroutine compute_F_om(e_factor, D1, D2, time, gamma, F_om)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z

    complex(wp), intent(in)::  e_factor(:)
    real(wp), intent(in):: D1(:) 
    real(wp), intent(in):: D2(:,:)
    real(wp), intent(in):: time(:)
    complex(wp), intent(out) ::  F_ifn_omp(:,:,:)
    real(wp), intent(in):: gamma

    integer:: ntsteps, nfinal
    complex(wp), allocatable:: funct(:)
    complex(wp), allocatable:: funct_fft(:)
    integer:: m1, m2 

    ntsteps = size(time) 
    
    allocate(funct(ntsteps), &
         funct_fft(ntsteps))
        
    F_om = 0.0_wp
    
    do m1=1,3 ! polarization
      
      funct = D2(:,m1) * e_factor(:) * exp(-gamma * const % eV * time(:) / const % hbar)
        
      call fft_c2c_1d_backward(funct, funct_fft)
      call reorder_sigma_fftw_z(funct_fft)

      do m2=1,3 ! polarization     
        F_om(:,m1,m2) = D1(m2) * funct_fft(:)
      end do ! m2 
        
    end do ! m1

  end subroutine compute_F_om


  ! computes
  ! F_om2_om1(om2,om1) = int_0^{t_max} e_factor(t) F_om2_om1(t, om1) exp(-gamma t ) exp(-i (om2-omega_mean)) exp(-i om2 t) d t
  subroutine compute_F_om2_om1(F_t_om1, e_factor, time, &
         gamma, omega1, omega1_mean, F_om2_om1)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_forward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    complex(wp), intent(in) ::  F_if_t_om1(:,:,:,:)
    complex(wp), intent(in) ::  e_factor(:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: gamma 
    real(wp), intent(in):: omega1(:)
    real(wp), intent(in):: omega1_mean
    complex(wp), intent(out) ::  F_om2_om1(:,:,:,:)

    integer:: n_omega1, om1, m1, m2
    
    n_omega1 = size(F_t_om1,1)
    
    do om1= 1, n_omega1
      do m1 =1, 3
        do m2 =1, 3
          
          call fft_c2c_1d_forward( e_factor1(:) * &
               exp(-gamma_instr * const % eV * time(:) / const % hbar) * &
               exp(dcmplx(0.0_wp,  (omega1(om1) - omega1_mean)* const % eV * time(:) / const % hbar )) * &
               F_t_om1(f_e, :, om1, m1,m2), &
               F_om2_om1(f_e, :, om1, m1,m2)) 
          call reorder_sigma_fftw_z(F_if_om2_om1(f_e, :, om1, m1,m2))
          
        end do
      end do
    end do
    
  end subroutine compute_F_om2_om1
  
  subroutine switch_place_om1_om2(F_om2_om1)

    complex(wp), intent(in) ::  F_om2_om1(:,:,:,:)

    complex(wp), allocatable:: F_tmp(:,:,:,:)
    integer:: m1,m2,om1,om2,nomega1

    n_omega1 = size(F_om2_om1,2)

    allocate(F_tmp(n_omega1, n_omega1,3,3))
    
    do m1 =1,3
      do m2 =1,3
        do om1 =1, n_omega1
          do om2 =1, n_omega1
            F_tmp(om1, om2, m1, m2) = F_om2_om1(om2, om1, m2, m1)
          end do
        end do
      end do
    end do

    F_om2_om1 = f_tmp
    
  end subroutine switch_place_om1_om2

  ! several intermediate states
  subroutine compute_F_if_om_many_n(E_n, E_i, E_ni_mean, D_fn, D_ni, time, F_if_om, gamma)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: fft_c2c_1d_backward
    use m_fftw3, only: reorder_sigma_fftw_z
    
    real(wp), intent(in):: E_n(:,:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: E_i(:)
    real(wp), intent(in):: D_fn(:,:,:) 
    real(wp), intent(in):: D_ni(:,:,:) 
    complex(wp), dimension(:,:,:,:),intent(out) ::  F_if_om
    real(wp), intent(in):: E_ni_mean, gamma
    
    integer:: nfinal, f_e, n_e, ntsteps
    complex(wp), allocatable:: F_tmp(:,:,:)

    nfinal = size(D_fn,1)

    F_if_om = 0.0_wp

            
    do n_e =1,ninter
      
      call compute_efactor(E_n(n_e, :), E_i(:), E_ni_mean, time, e_factor(:), .true.)
      
      do f_e =1,nfinal

        call compute_F_om(e_factor, D_ni(n_e,:,:), D_fn(f_e,:,:), time, gamma, F_tmp)
        
        F_if_om(f_e,:) = F_if_om(f_e,:) + F_tmp 

      end do

    end do
    
  end subroutine compute_F_if_om_many_n

  
  
  
end module m_SCKH_utils_res
