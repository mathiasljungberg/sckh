module m_fftw3
  implicit none
#include <fftw3.f>

#include "m_define_macro.F90"
  
  contains

!
! Simple 1d fft complex to complex
!

subroutine fft_c2c_1d_forward(vector_in, vector_out) 
  use m_precision, only : fftw_real

  complex(fftw_real), intent(in):: vector_in(:) 
  complex(fftw_real), intent(out):: vector_out(:) 

  integer(8):: plan
  integer:: rank

  rank = size(vector_in)
  call dfftw_plan_dft_1d(plan, rank, &
       vector_in, vector_out, FFTW_FORWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, vector_in, vector_out)
  call dfftw_destroy_plan(plan)

end subroutine fft_c2c_1d_forward

subroutine fft_c2c_1d_backward(vector_in, vector_out) 
  use m_precision, only : fftw_real

  complex(fftw_real), intent(in):: vector_in(:) 
  complex(fftw_real), intent(out):: vector_out(:) 

  integer(8):: plan
  integer:: rank

  rank = size(vector_in)
  call dfftw_plan_dft_1d(plan, rank, &
       vector_in, vector_out, FFTW_BACKWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, vector_in, vector_out)
  call dfftw_destroy_plan(plan)

end subroutine fft_c2c_1d_backward

!
! Simple 2d fftw complex to complex. Observe that no prefactor is included, it should be 1/size(vector_in)  
!
subroutine fft_c2c_2d_forward(vector_in, vector_out) 
  use m_precision, only : fftw_real
  
  complex(fftw_real), intent(in):: vector_in(:,:) 
  complex(fftw_real), intent(out):: vector_out(:,:) 

  integer(8):: plan
  integer:: ranks(2)
  
  ranks = shape(vector_in)
  call dfftw_plan_dft_2d(plan, ranks(1), ranks(2), &
       vector_in, vector_out, FFTW_FORWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, vector_in, vector_out)
  call dfftw_destroy_plan(plan)
  
end subroutine fft_c2c_2d_forward

!
!
!
subroutine fft_c2c_2d_backward(vector_in, vector_out) 
  use m_precision, only : fftw_real  
  complex(fftw_real), intent(in):: vector_in(:,:) 
  complex(fftw_real), intent(out):: vector_out(:,:) 
  
  integer(8):: plan
  integer:: ranks(2)
  
  ranks = shape(vector_in)

  call dfftw_plan_dft_2d(plan, ranks(1), ranks(2), &
       vector_in, vector_out, &
       FFTW_BACKWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, vector_in, vector_out)
  call dfftw_destroy_plan(plan)

end subroutine fft_c2c_2d_backward


  
!
! Simple 3d fftw complex to complex. Observe that no prefactor is included, it should be 1/size(vector_in)  
!
subroutine fft_c2c_3d_forward(vector_in, vector_out) 
  use m_precision, only : fftw_real

  complex(fftw_real), intent(in):: vector_in(:,:,:) 
  complex(fftw_real), intent(out):: vector_out(:,:,:) 

  integer(8):: plan
  integer:: ranks(3)

  ranks = shape(vector_in)
  call dfftw_plan_dft_3d(plan, ranks(1), ranks(2), ranks(3), &
       vector_in, vector_out, FFTW_FORWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, vector_in, vector_out)
  call dfftw_destroy_plan(plan)

end subroutine fft_c2c_3d_forward

!
!
!
subroutine fft_c2c_3d_backward(vector_in, vector_out) 
  use m_precision, only : fftw_real  
  complex(fftw_real), intent(in):: vector_in(:,:,:) 
  complex(fftw_real), intent(out):: vector_out(:,:,:) 
  
  integer(8):: plan
  integer:: ranks(3)
  
  ranks = shape(vector_in)

  call dfftw_plan_dft_3d(plan, ranks(1), ranks(2), ranks(3), &
       vector_in, vector_out, &
       FFTW_BACKWARD, FFTW_ESTIMATE)
  call dfftw_execute_dft(plan, vector_in, vector_out)
  call dfftw_destroy_plan(plan)

end subroutine fft_c2c_3d_backward

!
! for j in [1, N], gives the order [0,1,..., N//2, -(N-N//2), ...,-1]   
! that corresponds to how fftw calculates things from Fortran point of view.
! http://fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029
!
function get_order_fftw(j, n) result(j_fftw)
  implicit none
  !! external
  integer, intent(in):: j, n
  integer:: j_fftw
  !! internal

  if(j-1 .le. n/2) then
    j_fftw = (j-1)
  else
    j_fftw = - ( n - (j-1) )
  end if
  
end function !get_order_fftw

! reorders function as [-(N-N//2), ...,-1, 0, 1, ..., N/2]
subroutine reorder_sigma_fftw_z(fun)
  use m_precision, only : fftw_real
  
  complex(fftw_real), intent(inout):: fun(:)

  complex(fftw_real), allocatable :: fun_tmp(:)
  integer:: n, add, i,j
  
  
  n = size(fun)
  allocate(fun_tmp(n))

  !if( (n/2) * 2 .eq. n) then
  !  ! even function, for 6 tex [0,1,2,-3,-2,-1]
  !  add = 1
  !else
  !  ! odd function, for 5 tex [0,1,2,-2,-1]
  !  add=1
  !end if

  j=1
  do i=n/2 + 1, n
    fun_tmp(j) = fun(i)
    j=j+1
  end do

  do i=1, n/2 
    fun_tmp(j) = fun(i)
    j=j+1
  end do

  fun = fun_tmp
  
end subroutine reorder_sigma_fftw_z


subroutine get_omega_reordered_fftw(x_l, omega)
  use m_precision, only : wp
  use m_constants, only: const
  
  real(wp), intent(in):: x_l
  real(wp), intent(out):: omega(:)

  integer:: n, i,ii, j, add
  
  n = size(omega)
  
  !if( (n/2) * 2 .eq. n) then
  !  ! even function, for 6 tex [0,1,2,-3,-2,-1]
  !  add = 1
  !else
  !  ! odd function, for 5 tex [0,1,2,-2,-1]
  !  add=1
  !end if

  j=1
  do i= n/2 + 1, n
    ii = - ( n - (i-1) )
    omega(j) = 2.0_wp * const % pi * ii  / x_l 
    j=j+1
  end do
  do i=1, n/2  !0, n/2 -1
    ii = (i-1)
    omega(j) =  2.0_wp * const % pi * ii / x_l 
    j=j+1
  end do

  !write(6,*) "j",j
  !write(6,*) "omega", omega * const % cm
end subroutine get_omega_reordered_fftw

end module m_fftw3
