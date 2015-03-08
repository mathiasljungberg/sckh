 program vib_normal3
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile
  real(kind=wp), dimension(3,9):: normal
  real(kind=wp), dimension(:,:),allocatable:: q
  real(kind=wp), dimension(9):: x
  real(kind=wp), dimension(:),allocatable:: value
  integer:: nlines
  
  !other
  real(kind=wp), dimension(3):: r_o, r_h1, r_h2, r_tmp1, r_tmp2
  real(kind=wp), dimension(3):: r_tmp1_norm,r_tmp2_norm,r_tmp_sum,r_tmp_sum_norm
  real(kind=wp), dimension(3):: cross
  real(kind=wp), dimension(:),allocatable:: theta_out,stretch1,stretch2
  real(kind=wp):: rad_to_deg, r_OH1, r_OH2, theta, deg_to_rad
  real(kind=wp), dimension(3):: coord_O, coord_H1, coord_H2

  !loop variables
  integer:: i

  !functions
  real(kind=wp):: dnrm2


  r_OH1 = 0.95792059
  r_OH2 = r_OH1
  theta = 104.4996469200000
  deg_to_rad = pi/180.0_wp
  rad_to_deg = 180.0_wp/pi

  !equilibrium geometry of water molecule
  coord_O= (/0,0,0/)
  coord_H1= (/sin(theta/2*deg_to_rad)*R_OH1, cos(theta/2*deg_to_rad)*R_OH1, 0/)
  coord_H2= (/ - sin(theta/2*deg_to_rad)*R_OH2, cos(theta/2*deg_to_rad)*R_OH2, 0/)
  


  ! read from standard input
  read(5,*) inputfile
  read(5,*) outputfile
  read(5,*) nlines
  read(5,*) (normal(1,i), i=1,3)
  read(5,*) (normal(1,i), i=4,6)
  read(5,*) (normal(1,i), i=7,9)
  read(5,*) (normal(2,i), i=1,3)
  read(5,*) (normal(2,i), i=4,6)
  read(5,*) (normal(2,i), i=7,9)
  read(5,*) (normal(3,i), i=1,3)
  read(5,*) (normal(3,i), i=4,6)
  read(5,*) (normal(3,i), i=7,9)

  allocate(q(3,nlines), value(nlines), stretch1(nlines), stretch2(nlines),&
       theta_out(nlines))

  open(10,file=inputfile,status='old')
  open(11,file=outputfile,status='unknown')

  do i=1,nlines
     
     read(10,*) q(1,i),q(2,i),q(3,i), value(i)
     
  end do

  !stop

  do i=1,nlines
 
 ! convert normal coordinates to xyz
  x = q(1,i)*normal(1,:) +q(2,i)*normal(2,:) +q(3,i)*normal(3,:)

 ! write(6,*) i

  ! convert to internal coordinates
  r_o = x(1:3) + coord_O
  r_h1 = x(4:6) + coord_H1
  r_h2 = x(7:9) + coord_H2

  stretch1(i) = dnrm2(3,r_h1 - r_o,1)
  stretch2(i) = dnrm2(3,r_h2 - r_o,1)


  r_tmp1 = r_h1 - r_o 
  r_tmp2 = r_h2 - r_o 
  r_tmp1_norm = r_tmp1/dnrm2( 3,r_tmp1,1)
  r_tmp2_norm = r_tmp2/dnrm2( 3,r_tmp2,1)
  r_tmp_sum = r_tmp1_norm + r_tmp2_norm
  r_tmp_sum_norm =r_tmp_sum/dnrm2( 3,r_tmp_sum,1)
  
  call crossprod(cross, r_tmp1_norm, r_tmp_sum_norm)

  ! arcsin defined to be -90 to 90 degrees

  theta_out(i) = 2*abs(rad_to_deg*asin( dnrm2( 3,cross,1)))


  write (11,'(4F16.8)') stretch1(i), stretch2(i), theta_out(i), value(i)

  !write cut of probability density with stretch1 = r_OH1
  if(abs(stretch1(i)-r_OH1).lt.0.0005 ) then
     write(12,'(4F16.8)') stretch1(i), stretch2(i), theta_out(i), value(i)
  end if


 end do


!compute cross


end program vib_normal3
