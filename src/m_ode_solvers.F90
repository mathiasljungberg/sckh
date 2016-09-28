module m_ode_solvers

  implicit none
  
contains
  
subroutine ODE_solver(A_matrix,times,nstates,ntsteps)
    use m_precision,only:wp 
    use m_func, only: funct_complex
    use m_rkf45_matrix, only: rkfs_matrix_c
    implicit none
    integer,intent(in)::nstates,ntsteps
    complex(wp),intent(out),dimension(nstates,nstates,ntsteps)::A_matrix
    real(wp),dimension(:,:,:),allocatable::H_matrix
    complex(wp),allocatable,dimension(:,:)::y_value
    complex(wp),allocatable,dimension(:,:)::yp_value
    real(wp),dimension(ntsteps)::times
    real(wp) ::abserr
    !external f_04
    integer:: i_step
    integer:: iflag
    integer:: iwork(5)
    real(wp):: relerr
    real(wp):: t
    real(wp) :: t_out
    complex(wp),allocatable::  f1(:,:), f2(:,:), f3(:,:),f4(:,:), f5(:,:)
    real(wp)::  savre, savae
    integer:: i,j
    real(wp):: h,t_start,t_end
    ! set initial conditions for A matrix at t=0
    A_matrix(:,:,1)=0.0
    do i=1,ntsteps
       A_matrix(i,i,1)=(1.0_wp,0.0_wp)
    enddo
    allocate(H_matrix(nstates,nstates,ntsteps),y_value(nstates,nstates),yp_value(nstates,nstates))
    allocate(f1(nstates,nstates),f2(nstates,nstates),f3(nstates,nstates))
    allocate(f4(nstates,nstates),f5(nstates,nstates))
      abserr = 0.000000001e+00_wp
      relerr = 0.000000001e+00_wp
      iflag = 1
      y_value=A_matrix(:,:,1)
      h=times(2)-times(1)
      t_start=times(1)
      t_end=times(ntsteps)
      write(*,*) 'ODE_solver ', ntsteps 
      do i=1,ntsteps-1
        t = ( ( ntsteps - i + 1 ) * t_start+ ( i - 1 ) * t_end ) / dble( ntsteps )
        t_out = (( ntsteps - i ) * t_start+ ( i) * t_end ) / dble ( ntsteps )
        write(*,*) 't ',times(i),' t_out ',times(i+1),'step ',i
        call rkfs_matrix_c(nstates,nstates,funct_complex,y_value,times(i),times(i+1),relerr,abserr,&
             iflag,yp_value,h,f1,f2,f3,f4,f5,savre,savae,iwork(1),iwork(2),iwork(3),iwork(4), iwork(5))
         A_matrix(:,:,i+1)=y_value
      enddo
    deallocate(H_matrix,y_value,yp_value)
    deallocate(f1,f2,f3,f4,f5)
end subroutine ODE_solver

  
  subroutine ODE_solver_offdiagonal(A_matrix,times,nstates,ntsteps, ntsteps_inp)
    use m_precision, only: wp
    use m_func, only: funct_complex
    use m_rkf45_matrix, only: rkfs_matrix_c
    implicit none
    integer,intent(in)::nstates,ntsteps,ntsteps_inp
    complex(kind=wp),intent(out),dimension(nstates,nstates,ntsteps)::A_matrix
    complex(kind=wp),allocatable,dimension(:,:)::y_value
    complex(kind=wp),allocatable,dimension(:,:)::yp_value
    real(kind=wp),dimension(ntsteps)::times
    real(kind=wp) ::abserr
    !external f_04
    integer:: i_step
    integer:: iflag
    integer:: iwork(5)
    real(kind=wp):: relerr
    real(kind=wp):: t
    real(kind=wp) :: t_out
    complex(kind=wp),allocatable::  f1(:,:), f2(:,:), f3(:,:),f4(:,:), f5(:,:)
    real(kind=wp)::  savre, savae
    integer:: i,j
    real(kind=wp):: h,t_start,t_end
    ! set initial conditions for A matrix at t=0
    write(*,*) 'Time step is ',times(2)-times(1)
    A_matrix(:,:,1)=(0.0_wp,0.0_wp)
    do i=1,ntsteps
       A_matrix(i,i,1)=(1.0_wp,0.0_wp)
    enddo
    allocate(y_value(nstates,nstates),yp_value(nstates,nstates))
    allocate(f1(nstates,nstates),f2(nstates,nstates),f3(nstates,nstates))
    allocate(f4(nstates,nstates),f5(nstates,nstates))
      abserr = 0.000000001e+00_wp
      relerr = 0.000000001e+00_wp
      iflag = 1
      y_value=A_matrix(:,:,1)
      h=times(2)-times(1)
      t_start=times(1)
      t_end=times(ntsteps)
      !write(*,*) 'ODE_solver ', ntsteps
	 ! write(*,*) ' f1 ',shape(f1) 
	  !write(*,*) ' f2 ',shape(f2) 
	  !write(*,*) ' f3 ',shape(f3) 
	  !write(*,*) ' f4 ',shape(f4) 
	  !write(*,*) ' f5 ',shape(f5) 
	  !write(*,*) ' yp_value ',shape(yp_value)
	  !write(*,*) ' y_value ',shape(y_value)	  
	  !write(*,*) ' h ',shape(h)
      do i=1,ntsteps-1
        write(*,*) 't ',times(i),' t_out ',times(i+1),'step ',i,' nstates ',nstates
        call rkfs_matrix_c(nstates,nstates,funct_complex,y_value,times(i),times(i+1),relerr,abserr,&
          iflag,yp_value,h,f1,f2,f3,f4,f5,savre,savae,iwork(1),iwork(2),iwork(3),iwork(4), iwork(5))
		!write(*,*) 'Initialization A_matrix '
         A_matrix(:,:,i+1)=y_value
     !    write(*,*) abs(y_value)
      enddo
    deallocate(y_value,yp_value)
    deallocate(f1,f2,f3,f4,f5)
end subroutine ODE_solver_offdiagonal

 subroutine ODE_solver_diagonal(A_matrix,times,nstates,ntsteps, ntsteps_inp)
    use m_precision, only: wp
    use m_func, only: funct_complex_diagonal
    use m_rkf45_matrix, only: rkfs_matrix_c_diagonal
    implicit none
    integer,intent(in)::nstates,ntsteps,ntsteps_inp
    complex(kind=wp),intent(out),dimension(nstates,nstates,ntsteps)::A_matrix
    complex(kind=wp),allocatable,dimension(:)::y_value,yp_value
    real(kind=wp),dimension(ntsteps)::times
    real(kind=wp) ::abserr
    !external f_04
    integer:: i_step
    integer:: iflag
    integer:: iwork(5)
    real(kind=wp):: relerr
    real(kind=wp):: t
    real(kind=wp) :: t_out
    complex(kind=wp),allocatable,dimension(:)::  f1, f2, f3,f4,f5
    real(kind=wp)::  savre, savae
    integer:: i,j
    real(kind=wp):: h,t_start,t_end
    ! set initial conditions for A matrix at t=0
    allocate(f1(nstates),f2(nstates),f3(nstates))
    allocate(f4(nstates),f5(nstates))
    allocate(y_value(nstates),yp_value(nstates))
    write(*,*) 'Time step is ',times(2)-times(1)
    A_matrix=(0.0_wp,0.0_wp)
    do i=1,ntsteps
       A_matrix(i,i,1)=(1.0_wp,0.0_wp)
    enddo
      abserr = 0.000000001e+00_wp
      relerr = 0.000000001e+00_wp
      iflag = 1
      y_value=(1.0_wp,0.0_wp) ! put the initial conditions to the diagonal case
      h=times(2)-times(1)
      t_start=times(1)
      t_end=times(ntsteps)
      !write(*,*) 'ODE_solver ', ntsteps
	 ! write(*,*) ' f1 ',shape(f1) 
	  !write(*,*) ' f2 ',shape(f2) 
	  !write(*,*) ' f3 ',shape(f3) 
	  !write(*,*) ' f4 ',shape(f4) 
	  !write(*,*) ' f5 ',shape(f5) 
	  !write(*,*) ' yp_value ',shape(yp_value)
	  !write(*,*) ' y_value ',shape(y_value)	  
	  !write(*,*) ' h ',shape(h)
      do i=1,ntsteps-1
        !write(*,*) 't ',times(i),' t_out ',times(i+1),'step ',i,' nstates ',nstates
        call rkfs_matrix_c_diagonal(nstates,funct_complex_diagonal,y_value,times(i),times(i+1),relerr,abserr,&
          iflag,yp_value,h,f1,f2,f3,f4,f5,savre,savae,iwork(1),iwork(2),iwork(3),iwork(4), iwork(5))
		!write(*,*) 'Initialization A_matrix '
        !write(*,*) 'Time step is ', i,' ',abs(y_value)
        do j=1,nstates
         A_matrix(j,j,i+1)=y_value(j)
        enddo
     !    write(*,*) abs(y_value)
      enddo
    deallocate(y_value,yp_value)
    deallocate(f1,f2,f3,f4,f5)
end subroutine ODE_solver_diagonal


end module m_ode_solvers

