module m_fourier_grid
  implicit none
contains

  subroutine solve_fourier_grid(dx, V_i, eigval, eigvec, mu, units_in)
    use m_precision, only: wp
    use m_constants, only: const

    real(kind=wp), intent(in):: dx
    real(kind=wp), intent(in):: V_i(:)
    real(kind=wp), intent(out):: eigval(:)
    complex(kind=wp), intent(out):: eigvec(:,:)
    real(kind=wp), intent(in):: mu
    character(*), intent(in), optional:: units_in

    complex(wp), allocatable:: H_mat(:,:) 
    
    integer:: nstates, i,j,k
    character(80):: units    
    real(wp):: hbar, p, r_diff

    if(present(units_in)) then
      units = units_in
    else
      units = "SI"
    end if

    if (units .eq. "SI") then
      hbar = const % hbar
    else if (units .eq. "AU") then
      hbar =1.0_wp
    end if

    nstates = size(V_i)

    write(6,*) "hello form solve_fourier_grid"

    allocate(H_mat(nstates,nstates))
    
    H_mat =0.0_wp

    ! potential energy in real space
    do i=1, nstates
      H_mat(i,i) = V_i(i)
    end do

    ! kinetic energy, in fourier space
    do i=1,nstates
      do j=1,nstates
     
        do k =1,nstates
          !p = 2.0_wp * const % pi * (k-1) / (dx*nstates) 
          !r_diff = (i-j) !dx * (i-j)
          H_mat(i,j) =  H_mat(i,j) + exp( dcmplx(0,2.0_wp * const % pi * (i-j) * (k-1) / nstates)) * &
               ((k-1)**2) * (-hbar**2/ (2.0_wp * mu)) * (-1.0_wp) * (2.0_wp * const % pi)**2 / (dx * nstates)**2 
          !write(6,*) exp( dcmplx(0,r_diff * p)) * &
          !     (-p**2) * (-hbar**2)/ (2.0_wp * mu) / nstates

        end do
      end do
    end do

    eigval =0.0_wp
    call zHermitianEigen(H_mat, eigval, eigvec)
    
    write(6,*) eigval * const % cm 

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
  
  
end module m_fourier_grid
