module m_crank_nicolson
  use m_precision, only: wp
  implicit none

contains

  ! this solves the equation y'(t) = A(t) *y by the implicit trapezoid rule, aka Crank-Nicolson
  ! at times given by the array time
  ! A is given at times time_A, and this will be splined to time
  subroutine solve_crank_nicolson_matrix_z_spline(y_0, time, time_A, A, y)
    use m_splines, only: spline_z, splint_array3_z
    
    complex(wp), intent(in):: y_0(:,:)
    real(wp), intent(in):: time(:)
    real(wp), intent(in):: time_A(:)
    complex(wp), intent(in):: A(:,:,:)
    complex(wp), intent(out):: y(:,:,:)

    integer:: ntsteps, i, ntsteps_A, j 
    real(wp):: dt
    complex(wp), allocatable:: A_der(:,:,:)
    complex(wp), allocatable:: A_t(:,:)
    complex(wp), allocatable:: A_tdt(:,:)

    
    ntsteps = size(time)
    ntsteps_A = size(time_A)

    allocate(A_der(size(A,1),size(A,2),size(A,3)))
    allocate(A_t(size(A,1),size(A,2)))
    allocate(A_tdt(size(A,1),size(A,2)))
    
    ! spline all matrix elements
    do i=1,size(A,1)
      do j=1,size(A,2)
        call spline_z(time_A, A(i,j,:), ntsteps_A, 1.0d30,1.0d30, A_der(i,j,:))
      enddo
    enddo

    call splint_array3_z(time_A, A, A_der,time(1), A_t)

    y(:,:,1) = y_0
    
    do i=2, ntsteps

      dt = time(i) - time(i-1)

      call splint_array3_z(time_A, A, A_der, time(i), A_tdt)
      
      call  cn_matrix_z_step(y(:,:,i-1), A_t, A_tdt, dt, y(:,:,i))

      A_t = A_tdt
      
    end do
    
  end subroutine solve_crank_nicolson_matrix_z_spline





  
  ! this solves the equation y'(t) = A(t) *y by the implicit trapezoid rule, aka Crank-Nicolson
  ! at times given by the array time
  subroutine solve_crank_nicolson_matrix_z(y_0, time, A, y)
    complex(wp), intent(in):: y_0(:,:)
    complex(wp), intent(in):: A(:,:,:)
    real(wp), intent(in):: time(:)
    complex(wp), intent(out):: y(:,:,:)

    integer:: ntsteps, i
    real(wp):: dt
    
    ntsteps = size(time)

    y(:,:,1) = y_0
    !y_old = y_0
    
    do i=2, ntsteps

      dt = time(i) - time(i-1)
      
      call  cn_matrix_z_step(y(:,:,i-1), A(:,:,i-1), A(:,:,i), dt, y(:,:,i))
      
      !y(:,:,i) = y_new
      !y_old = y_new
      
    end do
    
  end subroutine solve_crank_nicolson_matrix_z
    

  ! assume y'(t) = A(t)*y(t)
  ! A_t= A(t)
  ! A_tdt = A(t+dt)
  !
  subroutine forward_euler_matrix_z_step(y_t, A_t, dt, y_tdt)
    use m_algebra, only: matmul_AB_z, matinv_z

    complex(wp), intent(in):: y_t(:,:)
    complex(wp), intent(in):: A_t(:,:)
    real(wp), intent(in):: dt
    complex(wp), intent(out):: y_tdt(:,:)

    integer:: n, i
    complex(wp), allocatable:: B1(:,:)
    
    n = size(y_t, 1)
    
    allocate(B1(n,n))
    
    B1 = 0.0d0
    do i=1, n
      B1(i,i) = 1.0d0
    end do
    B1 = B1 + A_t * dt
    
    !call matmul_Ax_z(B, y_t, y2)

    !B2 = 0.0d0
    !do i=1, n
    !  B2(i,i) = 1.0d0
    !end do
    !B2 = B2 - 0.5_wp * A_tdt * dt
    
    !call matinv_z(B2)

    !call matmul_AB_z(B2, B1, C)

    !call matmul_Ax_z(B, y2, y_tdt)

    call matmul_AB_z(B1, y_t, y_tdt)
    
  end subroutine forward_euler_matrix_z_step

  ! assume y'(t) = A(t)*y(t)
  ! A_t= A(t)
  ! A_tdt = A(t+dt)
  !
  subroutine backward_euler_matrix_z_step(y_t, A_tdt, dt, y_tdt)
    use m_algebra, only: matmul_AB_z, matinv_z

    complex(wp), intent(in):: y_t(:,:)
    complex(wp), intent(in):: A_tdt(:,:)
    real(wp), intent(in):: dt
    complex(wp), intent(out):: y_tdt(:,:)

    integer:: n, i
    complex(wp), allocatable:: B1(:,:), B2(:,:), C(:,:)
    
    ! y(t +dt) = (I - A(t +dt)*dt)^{-1} y(t) 
    n = size(y_t, 1)
    !m = size(y_t, 2)
    
    allocate(B2(n,n))
    
    B2 = 0.0d0
    do i=1, n
      B2(i,i) = 1.0d0
    end do
    B2 = B2 - A_tdt * dt
    
    call matinv_z(B2)

    !call matmul_AB_z(B2, B1, C)

    !call matmul_Ax_z(B, y2, y_tdt)

    call matmul_AB_z(B2, y_t, y_tdt)
    
  end subroutine backward_euler_matrix_z_step


  
  ! assume y'(t) = A(t)*y(t)
  ! A_t= A(t)
  ! A_tdt = A(t+dt)
  !
  subroutine cn_matrix_z_step(y_t, A_t, A_tdt, dt, y_tdt)
    use m_algebra, only: matmul_AB_z, matinv_z

    complex(wp), intent(in):: y_t(:,:)
    complex(wp), intent(in):: A_t(:,:)
    complex(wp), intent(in):: A_tdt(:,:)
    real(wp), intent(in):: dt
    complex(wp), intent(out):: y_tdt(:,:)

    integer:: n, i
    complex(wp), allocatable:: B1(:,:), B2(:,:), C(:,:)
    
    ! y(t +dt) = (I - 0.5*A(t +dt)*dt)^{-1} (I + 0.5*A(t)*dt ) y(t) 
    n = size(y_t, 1)
    !m = size(y_t, 2)
    
    allocate(B1(n,n), B2(n,n), C(n,n))
    
    B1 = 0.0d0
    do i=1, n
      B1(i,i) = 1.0d0
    end do
    B1 = B1 + 0.5_wp * A_t * dt
    
    !call matmul_Ax_z(B, y_t, y2)

    B2 = 0.0d0
    do i=1, n
      B2(i,i) = 1.0d0
    end do
    B2 = B2 - 0.5_wp * A_tdt * dt
    
    call matinv_z(B2)

    call matmul_AB_z(B2, B1, C)

    !call matmul_Ax_z(B, y2, y_tdt)

    call matmul_AB_z(C, y_t, y_tdt)
    
  end subroutine cn_matrix_z_step

!
!     TEST
!  

subroutine test1
  use m_precision, only: wp
  use m_func, only: set_H, func, y_exact_1
  use m_rkf45_matrix

!
! Test two uncoupled equations
! 
!*******************************************************************************
!

  implicit none

  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  real(kind=wp)::  dt
  real(kind=wp), allocatable:: time(:)

  complex(kind=wp), allocatable:: y(:,:), y_old(:,:)
  complex(kind=wp), allocatable:: H_input(:,:)
  integer:: i,j, i_step, ntsteps
  real(kind=wp):: h

  allocate(y(2,2), H_input(2,2), y_old(2,2))

  t_start = 0.0e+00_wp
  !t_stop = 2.0e+00_wp * 3.14159265e+00_wp
  t_stop = 10.0e+00_wp * 3.14159265e+00_wp
  ntsteps = 80 !12

  allocate(time(ntsteps))
  
  dt = (t_stop -t_start) / real(ntsteps -1,8)

  do i=1, ntsteps
    time(i) = (i-1)*dt + t_start
  end do

  y=0.0_wp
  y(1,1)=1.0_wp
  y(2,2)=1.0_wp

  ! this will give uncoupled equations
  H_input = 0.0_wp
  H_input(1,1) = (-1.0_wp, 1.0_wp)
  H_input(2,2) = (1.0_wp, -1.0_wp)

  t_out = time(1)

  write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) '       T          Y(1)           Y(2)        Y_exact(1)           Y_exact(2)       Y_error(1)        Y_error(2) '
  write ( *, '(a)' ) ' '

  write ( *, '(13g14.6)' ) t_out, y(1,1), y(2,2), exp(H_input(1,1) * t_out), exp(H_input(2,2) * t_out), &
       y(1,1) -exp(H_input(1,1) * t_out), y(2,2)- exp(H_input(2,2) * t_out)
  
  y_old = y
  
  do i_step = 2, ntsteps

    t_out = time(i_step)

    !call forward_euler_matrix_z_step(y_old, H_input, dt, y)
    call backward_euler_matrix_z_step(y_old, H_input, dt, y)
    !call cn_matrix_z_step(y_old, H_input, H_input, dt, y)
    
    write ( *, '(13g14.6)' ) t_out, y(1,1), y(2,2), exp(H_input(1,1) * t_out), exp(H_input(2,2) * t_out), &
         y(1,1) -exp(H_input(1,1) * t_out), y(2,2)- exp(H_input(2,2) * t_out)
    
    y_old = y
    
  end do

  return
end subroutine test1

subroutine test2
  use m_precision, only: wp
  use m_func, only: set_H, func, y_exact_1
  use m_rkf45_matrix

!
! Test two coupled equations : y'_11 = y_21, y'_21 = - y_11 and y'_12 = y_22, y'_22 = - y_12
! or H_input = [0 1; -1 0]
! eigenvalues +i, -i, period of 2*pi
!*******************************************************************************
!

  implicit none

  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  real(kind=wp)::  dt
  real(kind=wp), allocatable:: time(:)

  complex(kind=wp), allocatable:: y(:,:), y_old(:,:)
  complex(kind=wp), allocatable:: H_input(:,:)
  integer:: i,j, i_step, ntsteps
  real(kind=wp):: h

  allocate(y(2,2), H_input(2,2), y_old(2,2))

  t_start = 0.0e+00_wp
  !t_stop = 2.0e+00_wp * 3.14159265e+00_wp
  t_stop = 10.0e+00_wp * 3.14159265e+00_wp
  ntsteps = 100 !12

  allocate(time(ntsteps))
  
  dt = (t_stop -t_start) / real(ntsteps -1,8)

  do i=1, ntsteps
    time(i) = (i-1)*dt + t_start
  end do

  y=0.0_wp
  y(1,1)=1.0_wp
  y(2,2)=1.0_wp

  ! this will give coupled oscillatory eqs
  H_input = 0.0_wp
  H_input(1,2) = 1.0_wp !(-1.0_wp, 1.0_wp)
  H_input(2,1) = -1.0_wp !(1.0_wp, -1.0_wp)

  t_out = time(1)

  write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) '       T          Y(1)           Y(2)        Y_exact(1)           Y_exact(2)       Y_error(1)        Y_error(2) '
  write ( *, '(a)' ) ' '
  write ( *, '(9g14.6)' ) t_out, y(1,1), y(1,2), y(2,1), y(2,2)
  
  y_old = y
  
  do i_step = 2, ntsteps

    t_out = time(i_step)

    !call forward_euler_matrix_z_step(y_old, H_input, dt, y)
    !call backward_euler_matrix_z_step(y_old, H_input, dt, y)
    call cn_matrix_z_step(y_old, H_input, H_input, dt, y)
    
    write ( *, '(9g14.6)' ) t_out, y(1,1), y(1,2), y(2,1), y(2,2)
    
    y_old = y
    
  end do

  return
end subroutine test2

subroutine test3
  use m_precision, only: wp
  use m_func, only: set_H, func, y_exact_1
  use m_rkf45_matrix

!! Solve the diagonal matrix equation
!! d/dt A_{ab}(t) = i \delta_{ab} (E_a(t) -E_l(t)) A_{aa}(t)
!! A_{ab}(0) = \delta_{ab} 
!! this should be equal to A_{ab}(t) = \delta_{ab} exp(i int_0^t (E_a(\tau) -E_l(\tau)) d\tau )  
!!
!! take E_a(t) as C_a*cos(omega_a t +phi_a) and  E_l(t) as C_l*cos(omega_l t + phi_l)

  implicit none

  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  real(kind=wp)::  dt
  real(kind=wp), allocatable:: time(:)

  complex(kind=wp), allocatable:: y(:,:), y_old(:,:)
  complex(kind=wp), allocatable:: H_input1(:,:)
  complex(kind=wp), allocatable:: H_input2(:,:)
  integer:: i,j, i_step, ntsteps
  real(kind=wp):: h

  real(kind=wp), allocatable:: A_a(:), omega_a(:), &
       phi_a(:), A_l(:), omega_l(:), phi_l(:)
  real(kind=wp):: int_Ea(2)
  real(kind=wp):: int_El(2)
  complex(kind=wp):: y_exact(2,2)

  allocate(y(2,2), &
       H_input1(2,2), &
       H_input2(2,2), &
       y_old(2,2))

  t_start = 0.0e+00_wp
  !t_stop = 2.0e+00_wp * 3.14159265e+00_wp
  t_stop = 10.0e+00_wp * 3.14159265e+00_wp
  ntsteps = 100 !12

  allocate(time(ntsteps))
  
  dt = (t_stop -t_start) / real(ntsteps -1,8)

  do i=1, ntsteps
    time(i) = (i-1)*dt + t_start
  end do

  y=0.0_wp
  y(1,1)=1.0_wp
  y(2,2)=1.0_wp

  allocate(A_a(2))
  allocate(omega_a(2))
  allocate(phi_a(2))
  allocate(A_l(2))
  allocate(omega_l(2))
  allocate(phi_l(2))

  A_a =0.0d0
  omega_a =0.0d0
  phi_a =0.0d0
  A_l =0.0d0
  omega_l =0.0d0
  phi_l =0.0d0
  !h=1.0e-9_wp
  
  do i=1,2
    A_a(i) = i
    omega_a(i) = 0.3d0 / i * 3.14159265e+00_wp
    phi_a(i) = 0.14d0 * i *  3.14159265e+00_wp
    A_l(i) = 0.7d0 / i
    omega_l(i) = 0.8d0 / i * 3.14159265e+00_wp
    phi_l(i) = 0.8d0 * i * 3.14159265e+00_wp
  end do

  !call set_A_omega_phi_a( A_a, omega_a, phi_a)
  !call set_A_omega_phi_l( A_l, omega_l, phi_l)

  t_out = time(1)
  
  write ( *, '(a)' ) 'test_time_ordered_exp'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
  write ( *, '(a)' ) ' '
  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  y_old = y
  
  do i_step = 2, ntsteps

    H_input1 = 0.0_wp
    H_input2 = 0.0_wp
    do j=1,2
      H_input1(j,j) = (0.0_wp, 1.0_wp) * (A_a(j) * cos(omega_a(j) * time(i_step-1) + phi_a(j)) &
           - A_l(j) * cos(omega_l(j) * time(i_step-1) + phi_l(j)))
      H_input2(j,j) = (0.0_wp, 1.0_wp) * (A_a(j) * cos(omega_a(j) * time(i_step) + phi_a(j)) &
           - A_l(j) * cos(omega_l(j) * time(i_step) + phi_l(j)))
    end do
    
    !call forward_euler_matrix_z_step(y_old, H_input1, dt, y)
    !call backward_euler_matrix_z_step(y_old, H_input2, dt, y)    
    call cn_matrix_z_step(y_old, H_input1, H_input2, dt, y)
    
    y_old = y
    
    ! exact solution for t_out
    y_exact = 0.0d0
    t_out = time(i_step)
    do i=1,2
      int_Ea(i) = (A_a(i) / omega_a(i)) * &
           (cos(phi_a(i)) * sin(omega_a(i) * t_out) + sin(phi_a(i)) * (cos(omega_a(i) * t_out) -1.0d0)) 
      int_El(i) = (A_l(i) / omega_l(i)) * &
           (cos(phi_l(i)) * sin(omega_l(i) * t_out) + sin(phi_l(i)) * (cos(omega_l(i) * t_out) -1.0d0)) 
      y_exact(i,i) = exp(dcmplx(0.0d0, (int_Ea(i)- int_El(i))))
    end do
 
    write ( *, '(17g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2), &
         y_exact(1,1), y_exact(2,2), & 
         y(1,1)- y_exact(1,1), y(2,2)- y_exact(2,2) 

  end do

end subroutine test3





!subroutine test2
!  use m_precision, only: wp
!  use m_func
!  use m_rkf45_matrix
!!
!! Test two coupled equations : y'_11 = y_21, y'_21 = - y_11 and y'_12 = y_22, y'_22 = - y_12
!!*******************************************************************************
!!
!
!  implicit none
!
!  real(kind=wp) :: t_out
!  real(kind=wp) :: t_start
!  real(kind=wp)::  t_stop
!  !real(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
!  !     f4(neqn1,neqn2), f5(neqn1,neqn2)
!  !real(kind=wp)::  savre, savae
!
!  real(kind=wp):: y(neqn1, neqn2)
!  real(kind=wp):: H_input(neqn1, neqn2,1)
!  integer:: i,j
!  real(kind=wp):: h
!
!!
!!
!!  abserr = 0.000000001e+00_wp
!!  relerr = 0.000000001e+00_wp
!!
!!  iflag = 1
!
!  t_start = 0.0e+00_wp
!  t_stop = 2.0e+00_wp * 3.14159265e+00_wp
!  n_step = 12
!
!  dt = (t_stop -t_start) / (n_step -1)
!
!  do i=1, ntsteps
!    time(i) = (i-1)*dt + t_start
!  end do
!  
!  y=0.0_wp
!  y(1,1)=1.0_wp
!  y(1,2)=1.0_wp
!
!  H_input = 0.0_wp
!  H_input(1,2,1) = 1.0_wp
!  H_input(2,1,1) = -1.0_wp
!  
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
!  write ( *, '(a)' ) ' '
!  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  do i_step = 1, n_step
!
!    t = ( ( n_step - i_step + 1 ) * t_start &
!            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
!    t_out = ( ( n_step - i_step ) * t_start &
!            + ( i_step ) * t_stop ) / dble ( n_step )
!
!    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )
!
!    !call rkfs_matrix ( func, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
!    !     f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
!    
!  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  end do
!
!  return
!end subroutine test2
!
!
!
!subroutine test3
!  use m_precision, only: wp
!  use m_func
!  use m_rkf45_matrix
!  
!!
!! Test a matrix equation with constant corefficients, Y' = A Y solution:  Y = exp(A * t) * Y(0)
!! we set A to [1,1; 1,-1] because it is easy to solve, use the following octave program to test it 
!!
!!n_step = 12
!!t_start = 0.0
!!t_stop = 2.0 * 3.14159265
!!
!!Y0 = [1.0,0.0; 0.0, 1.0];
!!A = [1.0,1.0; 1.0, -1.0];
!!
!!for i_step = 2: n_step+1
!!  t = ( ( n_step - i_step + 1 ) * t_start + \
!!       ( i_step - 1 ) * t_stop ) / n_step
!!  Y0 * expm(A*t)
!!end
!!
!! 
!!*******************************************************************************
!!
! 
!  implicit none
!!
!  integer, parameter :: neqn1 = 2
!  integer, parameter :: neqn2 = 2
!!
!  real(kind=wp) ::abserr
!  !external f_04
!  integer:: i_step
!  integer:: iflag
!  integer:: iwork(5)
!  integer:: n_step
!  real(kind=wp):: relerr
!  real(kind=wp):: t
!  real(kind=wp) :: t_out
!  real(kind=wp) :: t_start
!  real(kind=wp)::  t_stop
!  real(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
!       f4(neqn1,neqn2), f5(neqn1,neqn2)
!  real(kind=wp)::  savre, savae
!
!  real(kind=wp):: y(neqn1, neqn2)
!  real(kind=wp):: H_input(neqn1, neqn2,1)
!  integer:: i,j
!  real(kind=wp):: h
!
!!
!
!  abserr = 0.000000001e+00_wp
!  relerr = 0.000000001e+00_wp
!
!  iflag = 1
!
!  t_start = 0.0e+00_wp
!  t_stop = 2.0e+00_wp * 3.14159265e+00_wp
!
!  n_step = 12
!
!  t_out = 0.0e+00_wp
!
!  y=0.0_wp
!  y(1,1)=1.0_wp
!  y(2,2)=1.0_wp
!
!  H_input = 1.0_wp
!  H_input(2,2,1) = -1.0_wp
!
!  call set_H( H_input )
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
!  write ( *, '(a)' ) ' '
!  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  do i_step = 1, n_step
!
!    t = ( ( n_step - i_step + 1 ) * t_start &
!            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
!    t_out = ( ( n_step - i_step ) * t_start &
!            + ( i_step ) * t_stop ) / dble ( n_step )
!
!    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )
!
!    call rkfs_matrix ( func, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
!         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
!    
!  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  end do
!
!  return
!end subroutine test3
!
!subroutine test4
!  use m_precision, only: wp
!  use m_func
!  use m_rkf45_matrix
!  
!!
!! Test a cmplex matrix equation with constant corefficients, Y' = iA Y solution:  Y = Y(0)*exp(iA * t)     
!! we set A to [1,1; 1,-1] because it is easy to solve, use the following octave program to test it 
!!
!!n_step = 12
!!t_start = 0.0
!!t_stop = 2.0 * 3.14159265
!!
!!Y0 = [1.0,0.0; 0.0, 1.0];
!!A = [1.0,1.0; 1.0, -1.0];
!!
!!for i_step = 2: n_step+1
!!  t = ( ( n_step - i_step + 1 ) * t_start + \
!!       ( i_step - 1 ) * t_stop ) / n_step
!!  Y0 * expm(j*A*t)
!!end
!!
!! 
!!*******************************************************************************
!!
!
!  implicit none
!!
!  integer, parameter :: neqn1 = 2
!  integer, parameter :: neqn2 = 2
!!
!  real(kind=wp) ::abserr
!  !external f_04
!  integer:: i_step
!  integer:: iflag
!  integer:: iwork(5)
!  integer:: n_step
!  real(kind=wp):: relerr
!  real(kind=wp):: t
!  real(kind=wp) :: t_out
!  real(kind=wp) :: t_start
!  real(kind=wp)::  t_stop
!  complex(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
!       f4(neqn1,neqn2), f5(neqn1,neqn2)
!  real(kind=wp)::  savre, savae
!
!  complex(kind=wp):: y(neqn1, neqn2)
!  real(kind=wp):: H_input(neqn1, neqn2,1)
!  integer:: i,j
!  real(kind=wp):: h
!
!!
!
!  abserr = 0.000000001e+00_wp
!  relerr = 0.000000001e+00_wp
!
!  iflag = 1
!
!  t_start = 0.0e+00_wp
!  t_stop = 2.0e+00_wp * 3.14159265e+00_wp
!
!  n_step = 12
!
!  t_out = 0.0e+00_wp
!
!  y=0.0_wp
!  y(1,1)=1.0_wp
!  y(2,2)=1.0_wp
!
!  H_input = 1.0_wp
!  H_input(2,2,1) = -1.0_wp
!
!  call set_H( H_input )
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
!  write ( *, '(a)' ) ' '
!  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  do i_step = 1, n_step
!
!    t = ( ( n_step - i_step + 1 ) * t_start &
!            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
!    t_out = ( ( n_step - i_step ) * t_start &
!            + ( i_step ) * t_stop ) / dble ( n_step )
!
!    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )
!
!    call rkfs_matrix_c ( func_c, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
!         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
!    
!  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  end do
!
!  return
!end subroutine test4
!
!
!subroutine test_time_ordered_exp
!  use m_precision, only: wp
!  use m_func_time_ordered_exp
!  use m_rkf45_matrix
!  
!!
!! Solve the diagonal matrix equation
!! d/dt A_{ab}(t) = i \delta_{ab} A_{aa}(t)(E_a(t) -E_l(t)) 
!! A_{ab}(0) = \delta_{ab} 
!! this should be equal to A_{ab}(t) = \delta_{ab} exp(i int_0^t (E_a(\tau) -E_l(\tau)) d\tau )  
!!
!! take E_a(t) as C_a*cos(omega_a t +phi_a) and  E_l(t) as C_l*cos(omega_l t + phi_l)
!!
!!
!!*******************************************************************************
!!
!
!  implicit none
!!
!  integer, parameter :: neqn1 = 2
!  integer, parameter :: neqn2 = 2
!!
!  real(kind=wp) ::abserr
!  integer:: i_step
!  integer:: iflag
!  integer:: iwork(5)
!  integer:: n_step
!  real(kind=wp):: relerr
!  real(kind=wp):: t
!  real(kind=wp) :: t_out
!  real(kind=wp) :: t_start
!  real(kind=wp)::  t_stop
!  complex(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
!       f4(neqn1,neqn2), f5(neqn1,neqn2)
!  real(kind=wp)::  savre, savae
!
!  complex(kind=wp):: y(neqn1, neqn2)
!  real(kind=wp):: H_input(neqn1, neqn2,1)
!  integer:: i,j
!  real(kind=wp):: h
!
!  real(kind=wp), allocatable:: A_a(:), omega_a(:), &
!       phi_a(:), A_l(:), omega_l(:), phi_l(:)
!  real(kind=wp):: int_Ea(2)
!  real(kind=wp):: int_El(2)
!  complex(kind=wp):: y_exact(2,2)
!
!!
!
!  abserr = 1.0e-9_wp
!  relerr = 1.0e-9_wp
!
!  iflag = 1
!
!  t_start = 0.0e+00_wp
!  t_stop = 10.0e+00_wp * 3.14159265e+00_wp
!
!  n_step = 12
!
!  t_out = 0.0e+00_wp
!
!  y=0.0_wp
!  y(1,1)=1.0_wp
!  y(2,2)=1.0_wp
!  
!  allocate(A_a(2))
!  allocate(omega_a(2))
!  allocate(phi_a(2))
!  allocate(A_l(2))
!  allocate(omega_l(2))
!  allocate(phi_l(2))
!
!  A_a =0.0d0
!  omega_a =0.0d0
!  phi_a =0.0d0
!  A_l =0.0d0
!  omega_l =0.0d0
!  phi_l =0.0d0
!  h=1.0e-9_wp
!  
!  do i=1,2
!    A_a(i) = i
!    omega_a(i) = 0.3d0 / i * 3.14159265e+00_wp
!    phi_a(i) = 0.14d0 * i *  3.14159265e+00_wp
!    A_l(i) = 0.7d0 / i
!    omega_l(i) = 0.8d0 / i * 3.14159265e+00_wp
!    phi_l(i) = 0.8d0 * i * 3.14159265e+00_wp
!  end do
!
!  call set_A_omega_phi_a( A_a, omega_a, phi_a)
!  call set_A_omega_phi_l( A_l, omega_l, phi_l)
!
!  write ( *, '(a)' ) 'test_time_ordered_exp'
!
!  write ( *, '(a)' ) ' '
!  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
!  write ( *, '(a)' ) ' '
!  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)
!
!  do i_step = 1, n_step
!
!    t = ( ( n_step - i_step + 1 ) * t_start &
!            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
!    t_out = ( ( n_step - i_step ) * t_start &
!            + ( i_step ) * t_stop ) / dble ( n_step )
!
!    call rkfs_matrix_c ( func_c, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
!         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
!
!    !write(6,'(A5,ES18.10)') "h", h
!    !write(6,*) "NFE", iwork(1)
!    
!    ! exact solution for t_out
!    y_exact = 0.0d0
!    do i=1,2
!      int_Ea(i) = (A_a(i) / omega_a(i)) * &
!           (cos(phi_a(i)) * sin(omega_a(i) * t_out) + sin(phi_a(i)) * (cos(omega_a(i) * t_out) -1.0d0)) 
!      int_El(i) = (A_l(i) / omega_l(i)) * &
!           (cos(phi_l(i)) * sin(omega_l(i) * t_out) + sin(phi_l(i)) * (cos(omega_l(i) * t_out) -1.0d0)) 
!      y_exact(i,i) = exp(dcmplx(0.0d0, (int_Ea(i)- int_El(i))))
!    end do
! 
!    write ( *, '(13g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2), &
!         y(1,1)- y_exact(1,1), y(2,2)- y_exact(2,2) 
!
!  end do
!
!  return
!end subroutine test_time_ordered_exp

  
end module m_crank_nicolson
