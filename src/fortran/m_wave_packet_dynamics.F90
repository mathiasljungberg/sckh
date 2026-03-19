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

  ! n and m ar eigenfunctions in orthogonal basis l (grid points)
  ! psi_n = <n | psi>
  ! psi_m = <m | psi>
  ! c_ln = < l | n >
  ! c_lm = < l | m >
  ! psi_m =  sum_l c_lm^*  sum_l c_ln * psi_n
  subroutine basis_transf_wfn_z(psi_n, c_ln, c_lm, psi_m)
    use m_precision, only: wp
    use m_algebra, only: matmul_Ax_z
    use m_algebra, only: matmul_Adagger_x_z
    
    complex(wp), intent(in):: psi_n(:)
    complex(wp), intent(in):: c_ln(:,:)
    complex(wp), intent(in):: c_lm(:,:)
    complex(wp), intent(out):: psi_m(:)

    complex(wp), allocatable:: psi_l(:)

    allocate(psi_l(size(psi_n)))
    
    call matmul_Ax_z(c_ln, psi_n, psi_l)
    call matmul_Adagger_x_z(c_lm, psi_l, psi_m)    
    
  end subroutine basis_transf_wfn_z
  
end module m_wave_packet_dynamics
