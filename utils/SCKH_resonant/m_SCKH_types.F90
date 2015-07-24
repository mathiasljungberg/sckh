module m_SCKH_types
  implicit none

  type pes_num_t
     real(kind=wp), allocatable:: x(:)
     real(kind=wp), allocatable:: E_i(:)
     real(kind=wp), allocatable:: E_n(:,:) ! ninterm, npoints
     real(kind=wp), allocatable:: E_f(:,:) ! nfinal, npoints
     real(kind=wp), allocatable:: D_ni(:,:,:,:) ! ninterm, m1, m2, npoints
     real(kind=wp), allocatable:: D_nf(:,:,:,:,:) ! ninterm, nfinal, m1, m2, npoints
  end type pes_num_t


!  type pes_poly_t
!  end type pes_poly_t


  type traj_t
     !     integer:: ntsteps
     real(kind=wp), allocatable:: x(:)
     real(kind=wp), allocatable:: v(:)
     real(kind=wp), allocatable:: a(:)
     real(kind=wp), allocatable:: time(:)
     real(kind=wp), allocatable:: E_i(:) 
     real(kind=wp), allocatable:: E_n(:,:) ! ninterm, ntsteps
     real(kind=wp), allocatable:: E_f(:,:) ! nfinal, ntsteps
     real(kind=wp), allocatable:: D_ni(:,:,:,:) ! ninterm, m1, m2, ntsteps
     real(kind=wp), allocatable:: D_nf(:,:,:,:,:) ! ninterm, nfinal, m1, m2, ntsteps
  end type traj_t

contains

subroutine init_traj(traj, ntsteps, ninter, nfinal, delta_t)

  integer:: i

  allocate(x(ntsteps))   
  allocate(v(ntsteps))   
  allocate(a(ntsteps))   
  allocate(time(ntsteps))
  allocate(E_i(ntsteps) 
  allocate(E_n(ninter, ntsteps) )
  allocate(E_f(nfinal, ntsteps) )
  allocate(D_ni(ninterm, 3,3, ntsteps) )
  allocate(D_nf(ninterm, nfinal, 3,3, ntsteps) )

  x=0.0_wp
  v=0.0_wp
  a=0.0_wp
  
  do i=1, ntsteps                                                                                                                                                                                                
    time(i)= (i-1)*delta_t                                                                                                                                                                                       
  end do

  E_i = 0.0_wp
  E_n = 0.0_wp
  E_f = 0.0_wp
  D_ni =0.0_wp
  D_nf =0.0_wp

end subroutine init_traj

end module m_SCKH_types
