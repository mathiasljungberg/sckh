program sinc_DVR
  use parameters
  use KH_functions
  use spline_m
  implicit none

  ! input/output
  character(80)::pes_file, outfile
  real(kind=wp):: my, dvr_start, dx
  integer::npoints_in, nstates

  ! loop variables
  integer::i,ii

  ! other variables
  integer:: npoints
  real(kind=wp):: my_SI
  real(kind=wp), dimension(:), allocatable:: X_dvr, x_in, E_in,eig_i,E_dvr
  real(kind=wp), dimension(:,:), allocatable:: c_i

  !
  ! This progam calculates the frequencies and eigenfunctions of a sinc DVR
  !

  !
  ! read input
  !

  read(5,*) pes_file
  read(5,*) npoints_in
  read(5,*) nstates, dvr_start, dx 
  read(5,*) my

  ! pes_file: the file with the potential energy surface, in Angstroms and Hartree units
  ! npoints_in; the number of points in pes_file
  ! nstates: the number of dvr points, and the number of eigenstates
  ! dvr_start: the atarting dvr point, in Angstroms
  ! dx: the spacing betwen dvr points, in Angstroms
  ! my: reduced mass, in atomic units

  allocate(X_dvr(nstates), x_in(npoints_in), E_in(npoints_in), eig_i(nstates), c_i(nstates,nstates))
  allocate(E_dvr(nstates))

  !
  ! set up DVR points
  !

  if (mod(nstates,2).ne.1 ) then
     write(6,*) "nstates must be an odd number"
     stop
  end if

  npoints = (nstates-1)/2
  my_SI = my * amu
  dvr_start = dvr_start * 1.0d-10
  dx = dx * 1.0d-10
  
  do i = -npoints,npoints
     ii = i + npoints +1
     X_dvr(ii) = (ii-1)*dx + dvr_start
  end do

  !
  ! read pes from file
  !

  open(10,file=pes_file,status='unknown')

  do i=1,npoints_in
     read(10,*) x_in(i), E_in(i)
  end do
  
  close(10) 
  
  ! use SI units
  x_in = x_in*1d-10
  E_in = E_in *hartree
  
  if(X_dvr(1).lt.x_in(1)) then
     if(abs(X_dvr(1)-x_in(1)) .lt. 1d-20 ) then
        X_dvr(1)=x_in(1)
     else
        write(6,*) "dvr point lower than supplied potential range!", X_dvr(1), x_in(1)
        stop
     end if
  end if
  if(X_dvr(nstates).gt.x_in(npoints_in)) then
     if(abs(X_dvr(nstates)-x_in(npoints_in)) .lt. 1d-20 ) then
        X_dvr(nstates)=x_in(npoints_in)
     else
        write(6,*) "dvr point higher than supplied potential range!", X_dvr(nstates), x_in(npoints_in), X_dvr(nstates)- x_in(npoints_in)
        stop
     end if
  end if
  call spline_easy(x_in, E_in, npoints_in, X_dvr, E_dvr, nstates)
    
  !
  ! Solve the vibrational problem for all eigenfunctions
  !
  
  call solve_sinc_DVR(dx,my_SI, E_dvr, c_i, eig_i)
  
  ! eV units
  eig_i =eig_i /eV

  !
  ! write output
  !

  write(6,*) "Inital state fundamental", (eig_i(2)-eig_i(1))*eV*cm
  write(6,*) "1:st overtone", (eig_i(3)-eig_i(2))*eV*cm
  write(6,*) "2:nd overtone", (eig_i(4)-eig_i(3))*eV*cm
  write(6,*) "3:rd overtone", (eig_i(5)-eig_i(4))*eV*cm
  
  
end program sinc_DVR
