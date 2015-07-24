program test
  use spline_m
  use KH_functions
  implicit none

  ! loop variables
  integer::i,j,ii,jj,k, j1,j2

  ! other variables
  real(kind=wp):: dx,d_sinc,dx_sinc, prod1, my_SI
  real(kind=wp),allocatable:: x_sinc(:), X_dvr(:)
  integer:: npoints,nstates,nx
  real(kind=wp), dimension(:,:),allocatable::   H_kin, H_kin2, M
  


  dx=0.075e-10_wp
  npoints = 20
  nstates = npoints*2 +1
  nx=100000
  my_SI= amu

allocate(x_sinc(nx), X_dvr(nstates), H_kin(nstates, nstates), H_kin2(nstates, nstates), M(nstates, nstates) )

  ! set up DVR points
  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx 
  end do
  
!write(6,*) X_dvr
!set up kinetic energy for sinc DVR                                                                                                                                                                                 
do i = -npoints,npoints
   ii = i + npoints +1
   H_kin(ii,ii) = (hbar**2  / (2 * my_SI * dx **2)) * (pi ** 2) /3d0
end do

do i = -npoints,npoints
   ii = i + npoints +1
   do j = i +1, npoints
      jj = j + npoints +1
     H_kin(ii,jj) =  (hbar**2 * (-1)**(i-j) / ( my_SI * dx **2)) / (i-j) **2
     H_kin(jj,ii) = H_kin(ii,jj)
  end do
end do



  call linspace(x_sinc,X_dvr(1) -10000*dx,X_dvr(nstates) + 10000*dx, nx)
  
  dx_sinc = x_sinc(2) - x_sinc(1)


! test kinetic energy matrix elements with DVR formula  H_kin{i,i'} = -hbar^2/2m * sum_{n=0}^N-1 chi(x_i) * chi''(x_{i'}) 
! THe above doesn't work at all: first, N -> inf, second: 

  !!do j=1,nstates
  !j=1
  !write(6,*) X_dvr(j)     
  !
  !   do i=2,nx-1
  !      d_sinc =  (1.d0/(dx_sinc * 2))*( &
  !           sinc(x_sinc(i+1), X_dvr(j),dx) &
  !           -  sinc(x_sinc(i-1),X_dvr(j),dx))  
  !      
  !      write(10,*) x_sinc(i), sinc(x_sinc(i), X_dvr(j),dx), &
  !           dsinc(x_sinc(i), X_dvr(j),dx)
  !      
  !   
  !   !d2_sinc = (-hbar**2/(dx_sinc**2*2*my_SI) )*( &
  !   !     sinc(x_sinc(i-1), X_dvr(j),dx) &
  !         !     + (-2) * sinc(x_sinc(i), X_dvr(j),dx)&
  !         !     +  sinc(x_sinc(i+1),X_dvr(j),dx))  
  !         
  !      end do
  !      write(10,*) 
  !      write(10,*) 

  !end do


 ! test scalar product 

  call dsinc_mat_elems(M, X_dvr)

!j2=1
do j1=1,3
   do j2=1,3
      prod1 = 0.0_wp
      !prod2 = 0.0_wp
      do i=2,nx-1
         
         !prod1 =  prod1 + sinc(x_sinc(i), X_dvr(j1),dx) * sinc(x_sinc(i), X_dvr(j2),dx)
         !prod1 =  prod1 + d2sinc(x_sinc(i), X_dvr(j1),dx) * sinc(x_sinc(i), X_dvr(j2),dx)
         prod1 =  prod1 +  sinc(x_sinc(i), X_dvr(j1),dx) * dsinc(x_sinc(i), X_dvr(j2),dx) 
        
      end do

      !prod1 = -(hbar**2  / (2 * my_SI)) * prod1 * dx_sinc
      prod1 = prod1 * dx_sinc

      write(6,*) dx_sinc, j1,j2, X_dvr(j1), X_dvr(j2), prod1, M(j1,j2) !H_kin(j1,j2)  !/ (X_dvr(2)-X_dvr(1)) 
   end do
end do




end program test
