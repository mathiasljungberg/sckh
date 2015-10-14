module m_func_time_ordered_exp
  use m_precision, only: wp
  implicit none

  real(kind=wp), private, allocatable:: A_a(:), omega_a(:), &
       phi_a(:), A_l(:), omega_l(:), phi_l(:)

contains
  
  subroutine func_c ( t, y, yp )
    real(kind=wp), intent(in):: t
    complex(kind=wp), intent(in) :: y(:,:)
    complex(kind=wp), intent(out) :: yp(:,:)
    
    integer:: i
    
    yp=0.0d0
    do i=1, size(yp,1)
      yp(i,i) = dcmplx(0.0d0,1.0d0) * y(i,i) * (A_a(i)*cos(omega_a(i)*t + phi_a(i)) &
           - A_l(i)*cos(omega_l(i)*t + phi_l(i))) 
    end do
    
  end subroutine func_c
  
  subroutine set_A_omega_phi_a(A_in, omega_in, phi_in)
    real(kind=wp), intent(in):: A_in(:), omega_in(:), &
         phi_in(:)
    
    if ( allocated(A_a) ) deallocate(A_a)
    allocate(A_a(size(A_in)))
    
    if ( allocated(omega_a) ) deallocate(omega_a)
    allocate(omega_a(size(omega_in)))
    
    if ( allocated(phi_a) ) deallocate(phi_a)
    allocate(phi_a(size(phi_in)))
    
    A_a     =  A_in    
    omega_a =  omega_in
    phi_a   =  phi_in  
    
  end subroutine set_A_omega_phi_a
  
  subroutine set_A_omega_phi_l(A_in, omega_in, phi_in)
    real(kind=wp), intent(in):: A_in(:), omega_in(:), &
         phi_in(:)
    
    if ( allocated(A_l) ) deallocate(A_l)
    allocate(A_l(size(A_in)))
    
    if ( allocated(omega_l) ) deallocate(omega_l)
    allocate(omega_l(size(omega_in)))
    
    if ( allocated(phi_l) ) deallocate(phi_l)
    allocate(phi_l(size(phi_in)))
    
    A_l     =  A_in    
    omega_l =  omega_in
    phi_l   =  phi_in  
    
  end subroutine set_A_omega_phi_l
  
  
end module m_func_time_ordered_exp
