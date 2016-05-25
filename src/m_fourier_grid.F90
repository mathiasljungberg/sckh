module m_fourier_grid
  implicit none
contains

  subroutine solve_fourier_grid(dx, V_i, eigval, eigvec, mu, units_in, mode_kin_energy_in)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: get_order_fftw
    
    real(kind=wp), intent(in):: dx
    real(kind=wp), intent(in):: V_i(:)
    real(kind=wp), intent(out):: eigval(:)
    complex(kind=wp), intent(out):: eigvec(:,:)
    real(kind=wp), intent(in):: mu
    character(*), intent(in), optional:: units_in
    character(*), intent(in), optional:: mode_kin_energy_in
    
    complex(wp), allocatable:: H_mat(:,:) 
    
    integer:: nstates, i
    character(80):: units
    character(80):: mode_kin_energy
    real(wp):: hbar
    real(wp):: prefac_me
    complex(wp), allocatable::  H_kin(:,:) 

    if(present(units_in)) then
      units = units_in
    else
      units = "SI"
    end if

    if(present(mode_kin_energy_in)) then
      mode_kin_energy = mode_kin_energy_in      
    else
      mode_kin_energy = "fast"
    end if

    if (units .eq. "SI") then
      hbar = const % hbar
    else if (units .eq. "AU") then
      hbar =1.0_wp
    end if

    nstates = size(V_i)

    !write(6,*) "hello from solve_fourier_grid"

    allocate(H_mat(nstates,nstates))
    allocate(H_kin(nstates,nstates))
    
    H_mat =0.0_wp

    ! potential energy in real space
    do i=1, nstates
      H_mat(i,i) = V_i(i)
    end do
    

    ! kinetic energy, in fourier space
    prefac_me= -hbar**2/ (2.0_wp * mu)

    H_kin = 0.0_wp
    if(mode_kin_energy .eq. "slow" ) then
      ! slow mode for testing purposes, N^3 scaling
      !write(6,*) "slow mode"
      call fourier_grid_d2_slow(nstates, dx, H_kin)
    else
      ! fast mode, default, N^2 scaling by precomputing some factors, noting that only difference i-j contributes
      call fourier_grid_d2_fast(nstates, dx, H_kin)            
    end if

    H_kin = prefac_me * H_kin 

    H_mat = H_mat + H_kin
    
    eigval =0.0_wp
    call zHermitianEigen(H_mat, eigval, eigvec)
    
  end subroutine solve_fourier_grid

  subroutine zHermitianEigen(A, E, X)
    !use m_precision, only : blas_int
    
    implicit none
    complex(8), intent(in)  :: A(:,:)
    real(8), intent(out)    :: E(:)
    complex(8), intent(out) :: X(:,:)
    
    !! internal
    integer :: ndim
    complex(8), allocatable :: A_copy(:,:)
    complex(8), allocatable :: CWORK(:)
    real(8), allocatable    :: RWORK(:)
    !integer(blas_int) :: info, lwork
    integer :: info, lwork
    complex(8) :: CWORK_TMP(1)
    
    ndim = size(A,1)
    allocate(A_copy(ndim, ndim))
    allocate(RWORK(3*ndim-2))
    call ZHEEV('V', 'U', ndim, A_copy, ndim, E, CWORK_TMP, -1, RWORK, info)
    
    if(info/=0) then
      write(0,*)'zSymmetricEigen: info/=0: calculate lwork:', info
      stop 
    endif
    
    lwork = int(CWORK_TMP(1));
    
    allocate(CWORK(lwork))
    A_copy = A(1:ndim, 1:ndim)
    call ZHEEV('V', 'U', ndim, A_copy, ndim, E, CWORK, lwork, RWORK, info)
    
    if(info/=0) then 
      write(0,*)'zSymmetricEigen: info/=0: diagonalize:', info 
      stop 
    endif
    
    X = A_copy
  end subroutine zHermitianEigen


  subroutine fourier_grid_d2_slow(nstates, dx, H_kin)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: get_order_fftw
    
    integer, intent(in):: nstates
    real(wp), intent(in):: dx
    complex(wp), intent(out):: H_kin(:,:)
    
    integer:: i,j,k, k_fftw
    real(wp):: recip_vec, prefac_d2, prefac
    
    recip_vec = 2.0_wp * const % pi  / (dx * nstates) 
    
    H_kin = 0.0_wp

    do i=1,nstates
      do j=1,nstates
        
        do k =1,nstates
          
          k_fftw = get_order_fftw(k, nstates)
          prefac_d2 = (-1.0_wp)* (k_fftw * recip_vec)**2
          !prefac = prefac_me * prefac_d2 / nstates
          prefac = prefac_d2 / nstates
          !          r_diff = (i-j) !dx * (i-j)
          
          !H_mat(i,j) =  H_mat(i,j) + exp( dcmplx(0,2.0_wp * const % pi * (i-j) * k / nstates)) * &
          !     (k_fftw ** 2) * (-hbar**2/ (2.0_wp * mu)) * (-1.0_wp) * (2.0_wp * const % pi)**2 / ((dx * nstates)**2 &
          !     *  nstates)
          H_kin(i,j) =  H_kin(i,j) + prefac * exp( dcmplx(0,2.0_wp * const % pi * (i-j) * k / nstates))
          
        end do
      end do
    end do
  end subroutine fourier_grid_d2_slow
  
  ! fast version of the routine above
  subroutine fourier_grid_d2_fast(nstates, dx, H_kin)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: get_order_fftw
    
    integer, intent(in):: nstates
    real(wp), intent(in):: dx
    complex(wp), intent(out):: H_kin(:,:)
    
    integer:: i,j,k, k_fftw
    real(wp):: recip_vec, prefac_d2, prefac
    complex(wp), allocatable::  phase_and_prefac_neg(:)
    
    recip_vec = 2.0_wp * const % pi  / (dx * nstates) 
    
    H_kin = 0.0_wp

    allocate(phase_and_prefac_neg(nstates))
    
    ! precompute prefactors
    phase_and_prefac_neg =0.0_wp
    do  i =1, nstates
      do k=1, nstates
        k_fftw = get_order_fftw(k, nstates)
        prefac_d2 = (-1.0_wp)* (k_fftw * recip_vec)**2
        prefac = prefac_d2 / nstates
        
        phase_and_prefac_neg(i) = phase_and_prefac_neg(i) +  prefac * exp( dcmplx(0,-2.0_wp * const % pi * (i-1) * k / nstates))
        
      end do
    end do
    
    ! compute matrix element
    do i=1,nstates
      H_kin(i,i) =  real(phase_and_prefac_neg(1), kind=8) ! hermitian matrix must have real diagonal element
      do j=i+1,nstates
        
        H_kin(i,j) =  phase_and_prefac_neg(j-i+1)       
        H_kin(j,i) =  conjg(phase_and_prefac_neg(j-i+1))       
        
      end do
    end do
    
  end subroutine fourier_grid_d2_fast
  
  ! first derivative
  subroutine fourier_grid_d1_fast(nstates, dx, H_kin)
    use m_precision, only: wp
    use m_constants, only: const
    use m_fftw3, only: get_order_fftw
    
    integer, intent(in):: nstates
    real(wp), intent(in):: dx
    complex(wp), intent(out):: H_kin(:,:)
    
    integer:: i,j,k, k_fftw
    real(wp):: recip_vec
    complex(wp):: prefac_d1, prefac
    complex(wp), allocatable::  phase_and_prefac_neg(:)
    
    recip_vec = 2.0_wp * const % pi  / (dx * nstates) 
    
    H_kin = 0.0_wp

    allocate(phase_and_prefac_neg(nstates))
    
    ! precompute prefactors
    phase_and_prefac_neg =0.0_wp
    do  i =1, nstates
      do k=1, nstates
        k_fftw = get_order_fftw(k, nstates)
        prefac_d1 = dcmplx(0.0_wp, -1.0_wp) * (k_fftw * recip_vec)
        prefac = prefac_d1 / nstates
        
        phase_and_prefac_neg(i) = phase_and_prefac_neg(i) +  prefac * exp( dcmplx(0,-2.0_wp * const % pi * (i-1) * k / nstates))

      end do
    end do
    
    ! compute matrix element
    do i=1,nstates
      H_kin(i,i) =  real(phase_and_prefac_neg(1), kind=8) ! hermitian matrix must have real diagonal element
      do j=i+1,nstates
        
        H_kin(i,j) =  phase_and_prefac_neg(j-i+1)       
        H_kin(j,i) =   conjg(phase_and_prefac_neg(j-i+1))       
        
      end do
    end do
    
  end subroutine fourier_grid_d1_fast


  
end module m_fourier_grid
