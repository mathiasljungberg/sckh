program DVR_harm
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile,basisfile,potentialfile,dummy
  integer::nin,npoints, nstates, alt1, alt2, alt3,alt4,alt5,alt6
  integer::nstates_in, npoints_in, npoints_pot
  real(kind=wp)::freq, my, dist_eq, dmy_x, dmy_y, dmy_z, dalpha
  !real(kind=wp), dimension(:),allocatable:: x,y

  ! loop variables
  integer::i,j

  ! other variables
  real(kind=wp):: my_SI
  
  integer:: INFO
  real(kind=wp):: omega_SI ,factor,d, alpha, mfreq, dx_in,sigma_my,sigma_alpha,mat_elem_q
  real(kind=wp), dimension(:),allocatable:: x_SI, V_in, W, WORK, X_dvr,V_dvr_1d,E_in
  real(kind=wp), dimension(:),allocatable:: x_pot, V_pot, x_tmp, V_tmp
  real(kind=wp), dimension(:,:),allocatable::  A, X_mat, Mat_tmp, T, H_harm,V_dvr,H_harm_dvr
  real(kind=wp), dimension(:,:),allocatable::  eigenvec, eigenvec_dvr

  !functions
  real(kind=wp):: harm_eigenfun
 ! real(kind=wp)::dlamch


!
! This progam is a potential optimized DVR (PODVR)
!

!
! read files
!*******************************************************************
!
  read(5,*) basisfile
  read(5,*) potentialfile
  read(5,*) my
  read(5,*) nstates
  read(5,*) alt1, npoints_pot, dist_eq
  read(5,*) dmy_x, dmy_y, dmy_z
  read(5,*) dalpha

  ! alt1 = 1 then the energies are exactly at the dvr points,
  ! alt1 = 2 then spline the energies onto the dvr points, taking
  ! equilibrium distance to be dist_eq
  ! my is reduced mass
  ! dmy_x, dmy_y, dmy_z are the dipole derivatives in x,y,z
  ! dalpha is the isotropic polarizability derivative
  !
  ! it computes IR intesities by  (d my_x / dq ) < n | q | m > 
  ! and computes Raman intensities as  (d alpha / dq ) < n | q | m >
  

  allocate(x_pot(nstates), V_pot(nstates))
  if (alt1.eq.1) then
     ! read file with ab initio potential on DVR points
     open(10,file=potentialfile,status='unknown')
     do i=1,nstates
        read(10,*) x_pot(i), V_pot(i)
     end do
     close(10)
  else if (alt1.eq.2) then
     allocate(x_tmp(npoints_pot),V_tmp(npoints_pot))
     open(10,file=potentialfile,status='unknown')
     do i=1,npoints_pot
        read(10,*) x_tmp(i), V_tmp(i)
     end do
     close(10)
  else 
     write(6,*) "alt1 shoud be either 1 or 2!"
     stop
  end if

  ! read eigenfunctions on grid points
  open(10,file=basisfile,status='unknown')
  
  read(10,*) dummy,  nstates_in
  read(10,*) dummy,  npoints_in
  read(10,*) 
  read(10,*) 

  allocate(x_SI(npoints_in), V_in(npoints_in),eigenvec(npoints_in, nstates_in), E_in(nstates_in))

  do i=1,nstates_in
     read(10,*) dummy, E_in(i)
  end do

  read(10,*) 

  do i=1,npoints_in
     read(10,*) x_SI(i), V_in(i)
  end do

  read(10,*)
  do i=1,npoints_in
     read(10,*) x_SI(i), eigenvec(i,1)
  end do
  
  do j=2, nstates_in
     read(10,*)
     do i=1,npoints_in
        read(10,*) dummy, eigenvec(i,j)
     end do
  end do
  
  close(10)

! Construct DVR point and DVR basis
!*******************************************************************


my_SI = my*amu
dx_in = x_SI(2)-x_SI(1)
!V_min = V_in(minloc(y_SI, DIM = 1))

!translate x to minimum
!write(6,*) "minloc",minloc(V_in, DIM = 1)
x_SI = x_SI - x_SI(minloc(V_in, DIM = 1))
 
allocate(A(nstates,nstates), X_mat(nstates,nstates), X_dvr(nstates),V_dvr_1d(nstates),Mat_tmp(nstates,nstates)&
     ,T(nstates,nstates), eigenvec_dvr(nstates,nstates))

! set up the matrix Xnm = < n | x | m > in basis functions that are read from file
do i=1,nstates
   do j=1,nstates
      X_mat(i,j) =  sum(eigenvec(:,i) * x_SI * eigenvec(:,j)) * dx_in
   end do
end do

!write(6,*) X_mat
Mat_tmp =X_mat


! find it's eigenvalues and eigenvectors  TXT_t = X_dvr, X_dvr is the diagonal matrix of quadrature points
! note that we get T_t from lapack!

! allocate workvectors for lapack
allocate( W(nstates), WORK(3*nstates))

! call lapack routine
call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, 3*nstates, INFO)
     
T = Mat_tmp
X_dvr = W

!write(6,*) "so far 7", INFO

! print eigenvalues and first frequency
write(6,*) "dvr points", X_dvr*1.d10


! Set up Hamiltonian and diagonalize
!*******************************************************************


! get potential in DVR points
call spline_easy(x_SI, V_in, npoints_in, X_dvr, V_dvr_1d, nstates)
     

! if alt1=2 get V_pot in dvr points
if(alt1.eq.2) then
   x_pot = X_dvr*1.d10 + dist_eq
   if(x_pot(1).lt.x_tmp(1)) then
      write(6,*) "dvr point lower than supplied potential range!"
      stop
   end if
   if(x_pot(nstates).gt.x_tmp(npoints_pot)) then
      write(6,*) "dvr point higher than supplied potential range!"
      stop
   end if
   call spline_easy(x_tmp, V_tmp, npoints_pot, x_pot, V_pot, nstates)
end if

! check if the DVR points are the same as the points read from potentialfile
!do i=1,nstates
!   if(abs(X_dvr(i)-x_pot(i)).gt.1.d-14) then
!      write(6,*) "error, the dvr point", i, "differs by", abs(X_dvr(i)-x_pot(i))
!      stop
!   end if
!end do

! set up the hamiltonian  
allocate( H_harm(nstates,nstates),H_harm_dvr(nstates,nstates),V_dvr(nstates,nstates) )

H_harm = 0
do i=1,nstates
   H_harm(i,i) = E_in(i) 
end do

! transform hamiltonian to dvr basis
H_harm_dvr = matmul(transpose(T),matmul(H_harm,T))

! add residual potential term at DVR points
V_dvr = 0
do i=1,nstates
   V_dvr(i,i) = V_pot(i)*hartree  -V_dvr_1d(i) 
end do

!test
!do i=1,npoints_in
!   write(7,*) x_SI(i), V_in(i)
!end do
!do i=1,nstates
!   write(8,*) X_dvr(i), V_dvr_1d(i)
!end do
!do i=1,nstates
!   write(9,*) X_dvr(i), d * (1 -exp(-alpha * (X_dvr(i)*1.0d10 -0.1_wp) )) ** 2 * hartree  -V_dvr_1d(i)
!end do
!do i=1,npoints_in
!   write(11,*) x_SI(i), d * (1 -exp(-alpha * (X_SI(i)*1.0d10 -0.1_wp) )) ** 2 * hartree 
!end do
!do i=1,npoints_in
!   write(12,*) x_SI(i), d * (1 -exp(-alpha * (X_SI(i)*1.0d10 -0.1_wp) )) ** 2 * hartree - V_in(i)
!end do


!write(6,*) V_dvr

Mat_tmp = H_harm_dvr + V_dvr


! solve eigenvalue problem
call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, 3*nstates, INFO)

eigenvec_dvr = Mat_tmp

!write(6,*) W*cm
write(6,*) "Fundamental frequency" , (W(2)-W(1))*cm
write(6,*) "Overtones" , (W(3)-W(2))*cm,  (W(4)-W(3))*cm, (W(5)-W(4))*cm

!
! Calculate matrix elements and intensities
!

! calculate matrix element < 0 | q | 1 > = sum_i c^0_i c^1_i x_i
mat_elem_q =  sum(eigenvec_dvr(:,1) * eigenvec_dvr(:,2) * X_dvr)

! calculate cross sections
sigma_my = (dmy_x**2 + dmy_y**2 +dmy_z**2) * mat_elem_q**2  * (W(2)-W(1))
sigma_alpha = dalpha**2 * mat_elem_q**2 * (W(2)-W(1))

write(6,*) "IR intensity", sigma_my
write(6,*) "Raman intensity", sigma_alpha

end program DVR_harm

