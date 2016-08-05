module m_wave_packet_dynamics
  implicit none
contains

  ! S_ni(:), complex(8) : the overlap of the initial state i, to eigenstates of the PES n
  !             = <n |i> =  <n |\psi(0)> at t=0
  ! omega_n(:), real(8)     : eigenvalues of the intermediate state, divided by \hbar
  ! c_ln(:,:), complex(8) : representation of the intermediate eigenstates in the basis l, < l | n>
  ! t, real(8)            : time        
  ! c_l_t(:), complex(8)  : the wave packet at time t, output. < l | \psi(t)> 
  
  subroutine wpd_eigenfun_z(S_ni, omega_n, c_ln, time, c_l_t)
    use m_precision, only: wp
    
    complex(wp), intent(in):: S_ni(:)
    real(wp), intent(in):: omega_n(:)
    complex(wp), intent(in):: c_ln(:,:)
    real(wp), intent(in):: time
    complex(wp), intent(out):: c_l_t(:)

    integer:: l
    
    do l=1, size(c_ln,2)
     c_l_t(l) = sum( exp(dcmplx(0,-omega_n(:) * time )) * c_ln(l,:) * S_ni(:)) 
    end do
        
  end subroutine wpd_eigenfun_z

  

  
end module m_wave_packet_dynamics
