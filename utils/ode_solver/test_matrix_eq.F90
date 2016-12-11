program test_matrix_eq
  use m_rkf45_matrix, only: timestamp

  implicit none

  call timestamp ( )

  call test1
  call test2
  call test3
  call test4
  call test_time_ordered_exp
  stop
  
contains

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
!
  integer, parameter :: neqn1 = 2
  integer, parameter :: neqn2 = 2
!
  real(kind=wp) ::abserr
  !external f_04
  integer:: i_step
  integer:: iflag
  integer:: iwork(5)
  integer:: n_step
  real(kind=wp):: relerr
  real(kind=wp):: t
  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  real(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
       f4(neqn1,neqn2), f5(neqn1,neqn2)
  real(kind=wp)::  savre, savae

  real(kind=wp):: y(neqn1, neqn2)
  real(kind=wp):: H_input(neqn1, neqn2,1)
  integer:: i,j
  real(kind=wp):: h

!

  abserr = 0.000000001e+00_wp
  relerr = 0.000000001e+00_wp

  iflag = 1

  t_start = 0.0e+00_wp
  t_stop = 2.0e+00_wp * 3.14159265e+00_wp

  n_step = 12


  t_out = 0.0e+00_wp

  y=1.0_wp

  ! this will give uncoupled equations
  H_input = 0.0_wp

  H_input(1,1,1) = 1.0_wp
  H_input(2,2,1) = -1.0_wp

  call set_H( H_input )

  write ( *, '(a)' ) ' '
  !write ( *, '(a)' ) '       T          Y(1)           Y(2)        Y_exact(1)           Y_exact(2)       Y_error(1)        Y_error(2) '
  write ( *, '(a)' ) ' '
  write ( *, '(7g14.6)' ) t_out, y(1,1), y(2,2), y_exact_1(t_out, 1.0_wp), y_exact_1(t_out, 33.0_wp), &
          y_exact_1(t_out, 1.0_wp) -y(1,1), y_exact_1(t_out, 33.0_wp) -y(2,2)

  do i_step = 1, n_step

    t = ( ( n_step - i_step + 1 ) * t_start &
            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
    t_out = ( ( n_step - i_step ) * t_start &
            + ( i_step ) * t_stop ) / dble ( n_step )

    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )

    call rkfs_matrix ( func, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
    

    write ( *, '(7g14.6)' ) t_out, y(1,1), y(2,2), y_exact_1(t_out, 1.0_wp), y_exact_1(t_out, -1.0_wp), &
          y_exact_1(t_out, 1.0_wp) -y(1,1), y_exact_1(t_out, -1.0_wp) -y(2,2)

  end do

  return
end subroutine test1



subroutine test2
  use m_precision, only: wp
  use m_func
  use m_rkf45_matrix
!
! Test two coupled equations : y'_11 = y_21, y'_21 = - y_11 and y'_12 = y_22, y'_21 = - y_12
! They should each match test 4 in rkf45_prb 
!*******************************************************************************
!

  implicit none
!
  integer, parameter :: neqn1 = 2
  integer, parameter :: neqn2 = 2
!
  real(kind=wp) ::abserr
  !external f_04
  integer:: i_step
  integer:: iflag
  integer:: iwork(5)
  integer:: n_step
  real(kind=wp):: relerr
  real(kind=wp):: t
  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  real(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
       f4(neqn1,neqn2), f5(neqn1,neqn2)
  real(kind=wp)::  savre, savae

  real(kind=wp):: y(neqn1, neqn2)
  real(kind=wp):: H_input(neqn1, neqn2,1)
  integer:: i,j
  real(kind=wp):: h

!

  abserr = 0.000000001e+00_wp
  relerr = 0.000000001e+00_wp

  iflag = 1

  t_start = 0.0e+00_wp
  t_stop = 2.0e+00_wp * 3.14159265e+00_wp

  n_step = 12

  t_out = 0.0e+00_wp

  y=0.0_wp
  y(1,1)=1.0_wp
  y(1,2)=1.0_wp

  H_input = 0.0_wp
  H_input(1,2,1) = 1.0_wp
  H_input(2,1,1) = -1.0_wp

  call set_H( H_input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  do i_step = 1, n_step

    t = ( ( n_step - i_step + 1 ) * t_start &
            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
    t_out = ( ( n_step - i_step ) * t_start &
            + ( i_step ) * t_stop ) / dble ( n_step )

    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )

    call rkfs_matrix ( func, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
    
  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  end do

  return
end subroutine test2



subroutine test3
  use m_precision, only: wp
  use m_func
  use m_rkf45_matrix
  
!
! Test a matrix equation with constant corefficients, Y' = A Y solution:  Y = Y(0)*exp(A * t)     
! we set A to [1,1; 1,-1] because it is easy to solve, use the following octave program to test it 
!
!n_step = 12
!t_start = 0.0
!t_stop = 2.0 * 3.14159265
!
!Y0 = [1.0,0.0; 0.0, 1.0];
!A = [1.0,1.0; 1.0, -1.0];
!
!for i_step = 2: n_step+1
!  t = ( ( n_step - i_step + 1 ) * t_start + \
!       ( i_step - 1 ) * t_stop ) / n_step
!  Y0 * expm(A*t)
!end
!
! 
!*******************************************************************************
!

  implicit none
!
  integer, parameter :: neqn1 = 2
  integer, parameter :: neqn2 = 2
!
  real(kind=wp) ::abserr
  !external f_04
  integer:: i_step
  integer:: iflag
  integer:: iwork(5)
  integer:: n_step
  real(kind=wp):: relerr
  real(kind=wp):: t
  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  real(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
       f4(neqn1,neqn2), f5(neqn1,neqn2)
  real(kind=wp)::  savre, savae

  real(kind=wp):: y(neqn1, neqn2)
  real(kind=wp):: H_input(neqn1, neqn2,1)
  integer:: i,j
  real(kind=wp):: h

!

  abserr = 0.000000001e+00_wp
  relerr = 0.000000001e+00_wp

  iflag = 1

  t_start = 0.0e+00_wp
  t_stop = 2.0e+00_wp * 3.14159265e+00_wp

  n_step = 12

  t_out = 0.0e+00_wp

  y=0.0_wp
  y(1,1)=1.0_wp
  y(2,2)=1.0_wp

  H_input = 1.0_wp
  H_input(2,2,1) = -1.0_wp

  call set_H( H_input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
  write ( *, '(a)' ) ' '
  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  do i_step = 1, n_step

    t = ( ( n_step - i_step + 1 ) * t_start &
            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
    t_out = ( ( n_step - i_step ) * t_start &
            + ( i_step ) * t_stop ) / dble ( n_step )

    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )

    call rkfs_matrix ( func, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
    
  write ( *, '(5g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  end do

  return
end subroutine test3

subroutine test4
  use m_precision, only: wp
  use m_func
  use m_rkf45_matrix
  
!
! Test a cmplex matrix equation with constant corefficients, Y' = iA Y solution:  Y = Y(0)*exp(iA * t)     
! we set A to [1,1; 1,-1] because it is easy to solve, use the following octave program to test it 
!
!n_step = 12
!t_start = 0.0
!t_stop = 2.0 * 3.14159265
!
!Y0 = [1.0,0.0; 0.0, 1.0];
!A = [1.0,1.0; 1.0, -1.0];
!
!for i_step = 2: n_step+1
!  t = ( ( n_step - i_step + 1 ) * t_start + \
!       ( i_step - 1 ) * t_stop ) / n_step
!  Y0 * expm(j*A*t)
!end
!
! 
!*******************************************************************************
!

  implicit none
!
  integer, parameter :: neqn1 = 2
  integer, parameter :: neqn2 = 2
!
  real(kind=wp) ::abserr
  !external f_04
  integer:: i_step
  integer:: iflag
  integer:: iwork(5)
  integer:: n_step
  real(kind=wp):: relerr
  real(kind=wp):: t
  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  complex(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
       f4(neqn1,neqn2), f5(neqn1,neqn2)
  real(kind=wp)::  savre, savae

  complex(kind=wp):: y(neqn1, neqn2)
  real(kind=wp):: H_input(neqn1, neqn2,1)
  integer:: i,j
  real(kind=wp):: h

!

  abserr = 0.000000001e+00_wp
  relerr = 0.000000001e+00_wp

  iflag = 1

  t_start = 0.0e+00_wp
  t_stop = 2.0e+00_wp * 3.14159265e+00_wp

  n_step = 12

  t_out = 0.0e+00_wp

  y=0.0_wp
  y(1,1)=1.0_wp
  y(2,2)=1.0_wp

  H_input = 1.0_wp
  H_input(2,2,1) = -1.0_wp

  call set_H( H_input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
  write ( *, '(a)' ) ' '
  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  do i_step = 1, n_step

    t = ( ( n_step - i_step + 1 ) * t_start &
            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
    t_out = ( ( n_step - i_step ) * t_start &
            + ( i_step ) * t_stop ) / dble ( n_step )

    !call rkf45_matrix ( func, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )

    call rkfs_matrix_c ( func_c, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )
    
  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  end do

  return
end subroutine test4


subroutine test_time_ordered_exp
  use m_precision, only: wp
  use m_func_time_ordered_exp
  use m_rkf45_matrix
  
!
! Solve the diagonal matrix equation
! d/dt A_{ab}(t) = i \delta_{ab} A_{aa}(t)(E_a(t) -E_l(t)) 
! A_{ab}(0) = \delta_{ab} 
! this should be equal to A_{ab}(t) = \delta_{ab} exp(i int_0^t (E_a(\tau) -E_l(\tau)) d\tau )  
!
! take E_a(t) as C_a*cos(omega_a t +phi_a) and  E_l(t) as C_l*cos(omega_l t + phi_l)
!
!
!*******************************************************************************
!

  implicit none
!
  integer, parameter :: neqn1 = 2
  integer, parameter :: neqn2 = 2
!
  real(kind=wp) ::abserr
  integer:: i_step
  integer:: iflag
  integer:: iwork(5)
  integer:: n_step
  real(kind=wp):: relerr
  real(kind=wp):: t
  real(kind=wp) :: t_out
  real(kind=wp) :: t_start
  real(kind=wp)::  t_stop
  complex(kind=wp):: yp(neqn1,neqn2) , f1(neqn1,neqn2), f2(neqn1,neqn2), f3(neqn1,neqn2),&
       f4(neqn1,neqn2), f5(neqn1,neqn2)
  real(kind=wp)::  savre, savae

  complex(kind=wp):: y(neqn1, neqn2)
  real(kind=wp):: H_input(neqn1, neqn2,1)
  integer:: i,j
  real(kind=wp):: h

  real(kind=wp), allocatable:: A_a(:), omega_a(:), &
       phi_a(:), A_l(:), omega_l(:), phi_l(:)
  real(kind=wp):: int_Ea(2)
  real(kind=wp):: int_El(2)
  complex(kind=wp):: y_exact(2,2)

!

  abserr = 1.0e-9_wp
  relerr = 1.0e-9_wp

  iflag = 1

  t_start = 0.0e+00_wp
  t_stop = 10.0e+00_wp * 3.14159265e+00_wp

  n_step = 12

  t_out = 0.0e+00_wp

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
  h=1.0e-9_wp
  
  do i=1,2
    A_a(i) = i
    omega_a(i) = 0.3d0 / i * 3.14159265e+00_wp
    phi_a(i) = 0.14d0 * i *  3.14159265e+00_wp
    A_l(i) = 0.7d0 / i
    omega_l(i) = 0.8d0 / i * 3.14159265e+00_wp
    phi_l(i) = 0.8d0 * i * 3.14159265e+00_wp
  end do

  call set_A_omega_phi_a( A_a, omega_a, phi_a)
  call set_A_omega_phi_l( A_l, omega_l, phi_l)

  write ( *, '(a)' ) 'test_time_ordered_exp'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T          Y(1,1)         Y(2,1)       Y(1,2)      Y(2,2)'
  write ( *, '(a)' ) ' '
  write ( *, '(9g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2)

  do i_step = 1, n_step

    t = ( ( n_step - i_step + 1 ) * t_start &
            + ( i_step - 1 ) * t_stop ) / dble ( n_step )
    t_out = ( ( n_step - i_step ) * t_start &
            + ( i_step ) * t_stop ) / dble ( n_step )

    call rkfs_matrix_c ( func_c, neqn1, neqn2, y, t, t_out, relerr, abserr, iflag, yp, h, &
         f1, f2, f3, f4, f5, savre, savae, iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )

    !write(6,'(A5,ES18.10)') "h", h
    !write(6,*) "NFE", iwork(1)
    
    ! exact solution for t_out
    y_exact = 0.0d0
    do i=1,2
      int_Ea(i) = (A_a(i) / omega_a(i)) * &
           (cos(phi_a(i)) * sin(omega_a(i) * t_out) + sin(phi_a(i)) * (cos(omega_a(i) * t_out) -1.0d0)) 
      int_El(i) = (A_l(i) / omega_l(i)) * &
           (cos(phi_l(i)) * sin(omega_l(i) * t_out) + sin(phi_l(i)) * (cos(omega_l(i) * t_out) -1.0d0)) 
      y_exact(i,i) = exp(dcmplx(0.0d0, (int_Ea(i)- int_El(i))))
    end do
 
    write ( *, '(13g14.6)' ) t_out, y(1,1), y(2,1), y(1,2), y(2,2), &
         y(1,1)- y_exact(1,1), y(2,2)- y_exact(2,2) 

  end do

  return
end subroutine test_time_ordered_exp



end program test_matrix_eq



!subroutine f_01 ( t, y, yp )
!!
!!*******************************************************************************
!!
!!! F_01 evaluates the derivative for the ODE.
!!
!  implicit none
!!
!  real t
!  real y(1)
!  real yp(1)
!!
!  yp(1) = 0.25E+00 * y(1) * ( 1.0E+00 - y(1) / 20.0E+00 )
!
!  return
!end subroutine f_01
!
!function y_exact_01 ( t )
!!
!!*******************************************************************************
!!
!!! Y_EXACT_01 evaluates the exact solution of the ODE.
!!
!  implicit none
!!
!  real t
!  real y_exact_01
!!
!  y_exact_01 = 20.0E+00 / ( 1.0E+00 + 19.0E+00 * exp ( - 0.25E+00 * t ) )
!
!  return
!end function y_exact_01
!
!subroutine f_02 ( t, y, yp )
!!
!!*******************************************************************************
!!
!!! F_02 evaluates the derivative for the ODE.
!!
!  implicit none
!!
!  real t
!  real y(2)
!  real yp(2)
!!
!  yp(1) = y(2)
!  yp(2) = - y(1)
!
!  return
!end subroutine f_02
!
!subroutine f_03 ( t, y, yp )
!!
!!*******************************************************************************
!!
!!! F_03 evaluates the derivative for the ODE.
!!
!  implicit none
!!
!  double precision t
!  double precision y(1)
!  double precision yp(1)
!!
!  yp(1) = 0.25D+00 * y(1) * ( 1.0D+00 - y(1) / 20.0D+00 )
!
!  return
!end subroutine f_03
!
!function y_exact_03 ( t )
!!
!!*******************************************************************************
!!
!!! Y_EXACT_03 evaluates the exact solution of the ODE.
!!
!  implicit none
!!
!  double precision t
!  double precision y_exact_03
!!
!  y_exact_03 = 20.0D+00 / ( 1.0D+00 + 19.0D+00 * exp ( - 0.25D+00 * t ) )
!
!  return
!end function y_exact_03
!
!
!subroutine f_04 ( t, y, yp )
!!
!!*******************************************************************************
!!
!!! F_04 evaluates the derivative for the ODE.
!!
!  implicit none
!!
!  double precision t
!  double precision y(2)
!  double precision yp(2)
!!
!  yp(1) = y(2)
!  yp(2) = - y(1)
!
!  return
!end subroutine f_04







