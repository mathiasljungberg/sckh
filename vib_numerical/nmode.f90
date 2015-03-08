program nmode
  use parameters
  implicit none
  
  integer:: flag1,flag2,flag3
  real(kind=wp)::r_OH1,r_OH2, theta, deg_to_rad,step, r_out1,r_out2,theta_out
  real(kind=wp):: f0, f2h, fm2h, fhk, fhmk, fmhk, fmhmk, fh,fmh, test, my1, my2, my3
  real(kind=wp):: scan1,scan2,scan3,scan_V, scan_start, scan_nsteps, scan_step
  real(kind=wp):: scan_start1, scan_start2, scan_start3, scan_step1, scan_step2, scan_step3
  real(kind=wp):: toang
  real(kind=wp), dimension(3)::coord_O, coord_H1, coord_H2, tmp
  real(kind=wp), dimension(9):: x, x2,x3 ,mass, nmode1_mw, nmode2_mw, nmode3_mw, nmode1,nmode2,nmode3
  real(kind=wp), dimension(9,9):: Hmat, Hmw
  real(kind=wp):: min_V,min_theta,R_OH1_2
  real(kind=wp), dimension(:), allocatable::scan_additional

  ! lapack variables
  integer:: INFO, LWORK
  real(kind=wp), dimension(:), allocatable:: W, WORK

  ! loop variables
  integer:: i,j,k,l

  ! functions
  real(kind=wp):: dnrm2

  ! the program calculates the normal modes of Tennyssons potential  


  ! put water molecule in coordinate system, z-plane
  
  deg_to_rad = pi/180.0_wp
  
  r_OH1 = 0.95792059
  r_OH2 = r_OH1
  theta = 104.4996469200000
  toang=0.5291772


  coord_O= (/0,0,0/)
  coord_H1= (/sin(theta/2*deg_to_rad)*R_OH1, cos(theta/2*deg_to_rad)*R_OH1, 0/)
  coord_H2= (/ - sin(theta/2*deg_to_rad)*R_OH2, cos(theta/2*deg_to_rad)*R_OH2, 0/)
  
!  mass = (/16,16,16,1,1,1,1,1,1/)
mass(1) = 15.99491463
mass(2) = mass(1)
mass(3) = mass(1)
mass(4) = 1.007825032   
mass(5) = mass(4)
mass(6) = mass(4)
mass(7) = mass(4) !2.014101778
mass(8) = mass(7)
mass(9) = mass(7)

  !molden file
  !write(6,*) 3
  !write(6,*)
  !write(6,*) "O", coord_O
  !write(6,*) "H", coord_H1
  !write(6,*) "H", coord_H2
  !write(6,*)



  ! calculate Hessian matrix in cartesian coordinates
  !***************************************************
  
  x = (/coord_O, coord_H1, coord_H2/) 
  
  ! very small step indeed
  step = 0.0001
  
  !write(6,*) x
  write(6,*) 
  ! diagonal elements, centered difference formula
  !SUBROUTINE POTS(V,Q1,Q2,THETA)

  call coords(x,x2,1,0.0, r_out1, r_out2, theta_out)
  
  call POTS( f0 ,r_out1,r_out2 ,theta_out*deg_to_rad)
    
  do i=1,9
     
     call coords(x,x2,i, step*2, r_out1, r_out2, theta_out)
     call POTS( f2h ,r_out1,r_out2 ,theta_out*deg_to_rad)
     
     call coords(x,x2,i, -step*2, r_out1, r_out2, theta_out)
     call POTS( fm2h ,r_out1,r_out2 ,theta_out*deg_to_rad)
     
     
     Hmat(i,i)= (f2h -2*f0 + fm2h)/(4*step**2)
     
  end do
  
  ! non diagonal elements
  
  do i=1,9
     do j=i+1,9
        
        call coords(x,x2,i, step, r_out1, r_out2, theta_out)
        call coords(x2,x3,j, step, r_out1, r_out2, theta_out)
        call POTS( fhk ,r_out1,r_out2 ,theta_out*deg_to_rad)
        
        call coords(x,x2,i, step, r_out1, r_out2, theta_out)
        call coords(x2,x3,j, -step, r_out1, r_out2, theta_out)
        call POTS( fhmk ,r_out1,r_out2 ,theta_out*deg_to_rad)
        
        call coords(x,x2,i, -step, r_out1, r_out2, theta_out)
        call coords(x2,x3,j, step, r_out1, r_out2, theta_out)
        call POTS( fmhk ,r_out1,r_out2 ,theta_out*deg_to_rad)
        
        call coords(x,x2,i, -step, r_out1, r_out2, theta_out)
        call coords(x2,x3,j, -step, r_out1, r_out2, theta_out)
        call POTS( fmhmk ,r_out1,r_out2 ,theta_out*deg_to_rad)
        
        Hmat(i,j)= (fhk - fhmk - fmhk + fmhmk )/(4*step**2)
        
     end do
  end do
  
  
  ! mass weigh coordinate matrix
  do i=1,9
     do j=i,9
        
        Hmw(i,j) = Hmat(i,j)/sqrt(mass(i)*mass(j)) 
        
     end do
  end do
  
  !  do i=1,9
  !     write(6,'(9F8.4)') (Hmat(i,j), j=1,9)
  !  end do
  
  
  ! diagonalize Hmw
  
  ! allocate workvectors for lapack
  !nstates3 = ind1
  LWORK=3*9 +10
  allocate( W(9), WORK(LWORK) )!, Hmatrixtmp(nstates3,nstates3))
  
  ! call lapack routine
  call dsyev( "V", "U", 9, Hmw, 9, W, WORK, LWORK, INFO)
  
  write(6,*) W 
  write(6,*)
  
  write(6,*) hbar*sqrt(W*hartree*1.d20/amu)*cm, amu, hartree
  
  !do i=1,9
  !   write(12,*)
  !   do j=1,9
  !      write(12,*) Hmw(j,i)
  !   end do
  !end do


! now, normalize normal modes, extract reduced masses, make 3d-scan

! x=M^(-1/2)*C*Q
nmode1_mw = Hmw(:,7)/sqrt(mass)
nmode2_mw = Hmw(:,8)/sqrt(mass)
nmode3_mw = Hmw(:,9)/sqrt(mass)

! normalize normal modes
nmode1 = nmode1_mw/(dnrm2(9, nmode1_mw,1))
nmode2 = nmode2_mw/(dnrm2(9, nmode2_mw,1))
nmode3 = nmode3_mw/(dnrm2(9, nmode3_mw,1))

! definition of reduced mass as in Gaussian 03
my1 = 1./(dnrm2(9, nmode1_mw,1)**2)
my2 = 1./(dnrm2(9, nmode2_mw,1)**2)
my3 = 1./(dnrm2(9, nmode3_mw,1)**2)

write(6,*) "my ",my1, my2,my3
write(6,*) "norm " , dnrm2(9, nmode1,1),dnrm2(9, nmode2,1),dnrm2(9, nmode3,1)
write(6,*) "nmode1", nmode1
write(6,*) "nmode2", nmode2
write(6,*) "nmode3", nmode3


! now perform potential energy surface scan
! *****************************************************'
scan_step1 = 0.04
scan_start1 =  -0.32

scan_step2 = 0.04
scan_start2 =  -0.32

scan_step3 = 0.04
scan_start3 =  -0.32

scan_nsteps = 16

x2=x

 
do i=1, scan_nsteps
   do j=1, scan_nsteps
      do k=1, scan_nsteps
         
         scan1 = scan_start1 +(i-1)*scan_step1 
         scan2 = scan_start2 +(j-1)*scan_step2 
         scan3 = scan_start3 +(k-1)*scan_step3 

         x2 = x + scan1*nmode1 + scan2*nmode2 +scan3*nmode3 

         call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
         call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)

         write(13,'(3F8.4,F18.12)') scan1, scan2, scan3, scan_V

!   write(14,*) "O", x2(1:3)
!   write(14,*) "H", x2(4:6)
!   write(14,*) "H", x2(7:9)

      end do
   end do
end do

! add additional points to fix asymptotic behaviour

if(0) then

allocate(scan_additional(8))

scan_additional = (/ -0.5, -0.4 ,-0.3 , 0.3, 0.4, 0.5 /)

scan_step1 = 0.1
scan_start1 =  -0.6

scan_step2 = 0.1
scan_start2 =  -0.6

scan_step3 = 0.1
scan_start3 =  -0.6

scan_nsteps = 13

x2=x



do i=1,scan_nsteps
   do j=1, scan_nsteps
      do k=1, scan_nsteps
         
         !scan1 = scan_additional(i)
         !scan2 = scan_additional(j)
         !scan3 = scan_additional(k)

         scan1 = scan_start1 +(i-1)*scan_step1 
         scan2 = scan_start2 +(j-1)*scan_step2 
         scan3 = scan_start3 +(k-1)*scan_step3 

         flag1=0
         flag2=0
         flag3=0

         do l=1,6
            if( abs(scan1-scan_additional(l)) < 1d-6 ) flag1=1
            if( abs(scan2-scan_additional(l)) < 1d-6 ) flag2=1
            if( abs(scan3-scan_additional(l)) < 1d-6 ) flag3=1
         end do


         if(flag1.eq.1.or.flag2.eq.1.or.flag3.eq.1) then

            x2 = x + scan1*nmode1 + scan2*nmode2 +scan3*nmode3 
            
            call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
            call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
            
            write(13,'(3F8.4,F18.12)') scan1, scan2, scan3, scan_V
         else
            write(6,*) "skipped", scan1, scan2, scan3

         end if

      end do
   end do
end do

do i=1,8
   scan1 = scan_additional(i)
   
   x2 = x + scan1*nmode1 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
   write(26,*) scan1, scan_V
   
end do
do i=1,8
   scan1 = scan_additional(i)
   
   x2 = x + scan1*nmode2 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
   write(27,*) scan1, scan_V
   
end do
do i=1,8
   scan1 = scan_additional(i)
   
   x2 = x + scan1*nmode3 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
   write(28,*) scan1, scan_V
   
end do


end if






!testing
!**********************

! write molden file
write(14,*) "6"
write(14,*)

scan_start = -0.2
do i=1, 11,10
   do j=1, 1
      do k=1,1

         scan1 = scan_start +(i-1)*scan_step1 
         scan2 = scan_start +(j-1)*scan_step2 
         scan3 = scan_start +(k-1)*scan_step3 
         x2 = x + scan1*nmode1 !+ scan2*nmode2 +scan3*nmode3 
   
!         call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
!         call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
         
         write(14,*) "O", x2(1:3)
         write(14,*) "H", x2(4:6)
         write(14,*) "H", x2(7:9)
         
         
      end do
   end do
end do

write(14,*)

! write PES cuts through normal modes
scan_start = -0.6
scan_step = 0.01
do i=1, 141
   
   scan1 = scan_start +(i-1)*scan_step 
   scan2 = scan_start +(i-1)*scan_step
   scan3 = scan_start +(i-1)*scan_step
   
   x2 = x + scan1*nmode1 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
   write(24,*) r_out1*toang, theta_out

   tmp(1)=scan_V
   
   x2 = x + scan2*nmode2 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
   tmp(2)=scan_V
   
   x2 = x + scan3*nmode3 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
   tmp(3)=scan_V
   
   
   !write(15,'(F8.4,3F14.8)') scan1, tmp

   write(16,*) scan1, tmp(1)
   write(17,*) scan2, tmp(2)
   write(18,*) scan3, tmp(3)

  end do


! write PES cut for bend mode
scan_start = -2.0
scan_step = 0.01
do i=1, 401
   
   scan1 = scan_start +(i-1)*scan_step 
!   scan2 = scan_start +(i-1)*scan_step
!   scan3 = scan_start +(i-1)*scan_step
   
   x2 = x + scan1*nmode1 
   
   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
!   write(24,*) r_out1*toang, theta_out

   tmp(1)=scan_V
   
!   x2 = x + scan2*nmode2 
   
!   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
!   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
!   tmp(2)=scan_V
   
!   x2 = x + scan3*nmode3 
   
!   call coords(x2,x3,1,0.0, r_out1, r_out2, theta_out)
!   call POTS( scan_V ,r_out1,r_out2 ,theta_out*deg_to_rad)
   
!   tmp(3)=scan_V
   
   
   !write(15,'(F8.4,3F14.8)') scan1, tmp

   write(25,*) scan1, tmp(1)
!   write(17,*) scan2, tmp(2)
!   write(18,*) scan3, tmp(3)

  end do








! write one dimensional scan in stretch coordinate
 r_OH1 = 0.95792059
 theta = 104.4996469200000
do i=1,2500
   r_OH2 = r_OH1 -0.7d0 +(i-1)*0.001_wp 

   call POTS( scan_V ,r_OH1/toang,r_OH2/toang ,theta*deg_to_rad)
   
   write(19,*) r_OH2, scan_V

end do


! write one dimensional scan in bend coordinate
  r_OH1 = 0.95792059_wp
  r_OH2 = r_OH1

do i=1,1500
   theta_out = theta  -15  +(i-1)*0.02 
   !r_OH2 = r_OH1 -0.5d0 +(i-1)*0.001 

   call POTS( scan_V ,r_OH1/toang,r_OH2/toang ,theta_out*deg_to_rad)
   
   write(20,*) theta_out, scan_V

end do



! write two dimensional scan in bend and sym stretch coordinate
  !r_OH1 = 0.95792059_wp
  !r_OH2 = r_OH1

do i=1,100

   r_OH2 = r_OH1 -0.3d0 +(i-1)*0.01 

min_V =100000
do j=1,100

   theta_out = theta  -20  +(j-1)*0.4

   call POTS( scan_V ,r_OH2/toang,r_OH2/toang ,theta_out*deg_to_rad)
   
!   write(21,*) theta_out, R_OH2, scan_V

   if (scan_V<min_V) then
      min_V = scan_V
      min_theta =theta_out
   end if

end do

write(21,*) R_OH2, min_theta

end do


! write two dimensional scan in bend and asym stretch coordinate
  !r_OH1 = 0.95792059_wp
  !r_OH2 = r_OH1

do i=1,100

   r_OH2 = r_OH1 -0.3d0 +(i-1)*0.01 
   r_OH1_2 = r_OH1 +0.3d0 -(i-1)*0.01 

min_V =100000
do j=1,100

   theta_out = theta  -20  +(j-1)*0.4

   call POTS( scan_V ,r_OH2/toang,r_OH1_2/toang ,theta_out*deg_to_rad)
   
!   write(21,*) theta_out, R_OH2, scan_V

   if (scan_V<min_V) then
      min_V = scan_V
      min_theta =theta_out
   end if

end do

write(22,*) R_OH2, min_theta

end do



! write two dimensional scan in bend and single stretch coordinate
  !r_OH1 = 0.95792059_wp
  !r_OH2 = r_OH1

do i=1,100

   r_OH2 = r_OH1 -0.3d0 +(i-1)*0.01 
!   r_OH1_2 = r_OH1 +0.3d0 -(i-1)*0.01 

min_V =100000
do j=1,100

   theta_out = theta  -20  +(j-1)*0.4

   call POTS( scan_V ,r_OH2/toang,r_OH1/toang ,theta_out*deg_to_rad)
   
!   write(21,*) theta_out, R_OH2, scan_V

   if (scan_V<min_V) then
      min_V = scan_V
      min_theta =theta_out
   end if

end do

write(23,*) R_OH2, min_theta

end do










end program nmode




subroutine coords(coords_in, coords_out, ncoordmove, step, r_out1, r_out2, theta_out )
  use parameters
  implicit none
  !passed variables
  real(kind=wp), dimension(9), intent(in)::coords_in
  real(kind=wp), dimension(9), intent(out)::coords_out
  integer, intent(in)::ncoordmove
  real(kind=wp), intent(in)::step
  real(kind=wp), intent(out):: r_out1, r_out2, theta_out
  !local variables
  real(kind=wp):: rad_to_deg, toang
  real(kind=wp), dimension(3):: r_tmp1, r_tmp2, r_tmp1_norm, r_tmp2_norm, r_tmp_sum, r_tmp_sum_norm, cross
  !functions
  real(kind=wp):: dnrm2

  ! the subroutine takes a vector of dim 9, moves one coordinate, transforms to bond angle coordinates
  ! and returns those coordinates in bohr and degrees

  rad_to_deg = 180.0_wp/pi
  toang=0.5291772_wp

  ! move atom in cartesian coordiantes
  coords_out = coords_in
  coords_out(ncoordmove) = coords_out(ncoordmove) + step
  

  r_tmp1 = coords_out(1:3)-coords_out(4:6)
  r_tmp2 = coords_out(1:3)-coords_out(7:9)
   
  r_out1 = dnrm2( 3,r_tmp1,1)/toang
  r_out2 = dnrm2( 3,r_tmp2,1)/toang

  r_tmp1_norm = r_tmp1/dnrm2( 3,r_tmp1,1)
  r_tmp2_norm = r_tmp2/dnrm2( 3,r_tmp2,1)
  r_tmp_sum = r_tmp1_norm + r_tmp2_norm
  r_tmp_sum_norm =r_tmp_sum/dnrm2( 3,r_tmp_sum,1)
  
  call crossprod(cross, r_tmp1_norm, r_tmp_sum_norm)

  ! arcsin defined to be -90 to 90 degrees
!  theta_out = 180.0_wp - rad_to_deg*asin( dnrm2( 3,cross,1)/(dnrm2( 3,r_tmp1,1)*dnrm2( 3,r_tmp2,1)))

  theta_out = 2*abs(rad_to_deg*asin( dnrm2( 3,cross,1)))

 end subroutine coords
