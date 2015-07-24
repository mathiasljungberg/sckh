module m_fourier_grid
  use parameters
  implicit none
contains

  subroutine fourier_grid(dx, my_SI, V_i, c_i, e_i)
    real(kind=wp), intent(in):: dx, my_SI
    real(kind=wp), dimension(:), intent(in):: V_i
    complex(kind=wp), dimension(:,:), intent(out):: c_i
    real(kind=wp), dimension(:), intent(out):: e_i

    complex(wp), allocatable:: H_mat(:,:), 
    
    integer:: nstates, i,j,k

    nstates = size(V_i)

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
          p = 2.0_wp * pi * (k-1) / nstates 
          r_diff = dx * (i-j)
          H_mat(i,j) =  H_mat(i,j) + exp( dcmplx(0,-rdiff * p)) * (-p**2 / (2.0_wp * my_SI)) / nstates
        end do

      end do
    end do

    call zHermitianEigen(H_mat, c_i, e_i)
    
    
end subroutine fourier_grid

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
  integer(blas_int) :: info, lwork
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
