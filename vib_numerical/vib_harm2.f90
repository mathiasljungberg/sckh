program vib_harm
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile
  integer::nin,npoints,nstates,alt1, alt2, alt3,alt4,alt5,alt6
  real(kind=wp)::my_in
  real(kind=wp), dimension(:),allocatable:: x,y

  ! loop variables
  integer::i,j

  ! other variables
  real(kind=wp), dimension(:),allocatable:: x_SI,y_SI,x_new,y_new, &
       tmp
  real(kind=wp)::dx, my, my_SI, x_new_min, x_new_max, coeff_tmp
  integer:: ymin_index, nfit_dim, nfit_p,nfit_s,nfit_int, INFO, poly_deg
  real(kind=wp)::omega_h, freq,x0,y0,corr2,corr3,corr4, freq_per
  real(kind=wp), dimension(:),allocatable:: x_inp, y_inp, W,WORK, eigenfun1, coeff,&
       start
  !real(kind=wp), dimension(3):: start
  !real(kind=wp), dimension(5):: start4, coeff4
  !real(kind=wp), dimension(7):: start6, coeff6
  real(kind=wp), dimension(:,:),allocatable:: eigenfun, Hmatrix, Htmp

  !functions
  real(kind=wp):: harm_eigenfun
 ! real(kind=wp)::dlamch


  ! read from standard input
  read(5,*) inputfile
  read(5,*) outputfile
  read(5,*) alt1, poly_deg
  read(5,*) alt3, my_in
  read(5,*) alt6, npoints
  read(5,*) nin, nstates 
  read(5,*) alt4, nfit_p, nfit_s
  read(5,*) alt5, x_new_min, x_new_max


  ! alt1=1, do numerical integration, alt1=2, do variational calculation on polynomial
  !         potential energy surface
  ! alt2   degree of fitted polynomial
  ! alt3=0, use my=18/19, alt3=1, use my_in 
  ! alt4=0, use whole potential surface in polynomial fitting, alt=1, use nfit_p and nfit_s
  ! nfit_p : number of fitting points on each side of middle
  ! nfit_s : points "skipped", =2 means every other point etc. 
  ! alt5=0, use the whole interval, alt5=1, use x_new_min to x_new_max
  ! alt6=1, spline, alt6=0 don't spline

  ! allocate
  allocate( x(nin),y(nin), x_SI(nin),y_SI(nin))

  if (alt6.eq.0) then
     npoints=nin
  end if
  
  allocate (x_new(npoints), &
          y_new(npoints))


  ! read inputfile
  open(10,file=inputfile,status='old')
  
  do i=1,nin
     read(10,*) x(i),y(i)
 enddo

  close(10)

  ! my

  my= 18._wp/19._wp
  if(alt3.eq.1) my=my_in
  my_SI = my*amu

  
  ! go to SI units
  y_SI=y*hartree
  x_SI=x*1.0e-10
  
  if(alt5.eq.0) then
     x_new_min = x_SI(1)
     x_new_max = x_SI(nin)
   else
     x_new_min = x_new_min*1.0d-10
     x_new_max = x_new_max*1.0d-10
  end if
 
  ! if alt6=1 spline
  if(alt6.eq.1) then

     call linspace(x_new, x_new_min, x_new_max, npoints)
     call spline_easy(x_SI, y_SI, nin, x_new, y_new, npoints)
  else  
     x_new = x_SI
     y_new = y_SI
  end if
  
  
  ! fit a polynomial around minimum
  !*********************************************************************  
  
  ! find minimum
  ymin_index = minloc(y_new, DIM = 1)  
  
  write(6,*) "ymin_index", ymin_index, y_new(ymin_index)
  
  nfit_dim = size(y_new)
  nfit_int = 1
  if (alt4.eq.1) then
     nfit_dim= nfit_p*2 + 1
     nfit_int = nfit_p*nfit_s
  end if
  
  write(6,*) "nfit_P", nfit_p, nfit_dim
  
  ! allocate input vectors
  allocate(y_inp(nfit_dim),x_inp(nfit_dim))
   y_inp=0
  x_inp=0
  
  if (alt4.eq.1) then
     y_inp = y_new(ymin_index - nfit_int : ymin_index + nfit_int: nfit_s )
     x_inp = x_new(ymin_index - nfit_int : ymin_index + nfit_int: nfit_s )
   else
     y_inp =y_new
     x_inp =x_new
  end if
  
  allocate(coeff(poly_deg+1), start(poly_deg+1))
  
  ! startvalues for fitting
  start = 0.1d20*hartree
  start(1) = x_new(ymin_index)
  start(poly_deg +1) = y_new(ymin_index)

  !write(6,*) "x_inp",x_inp
  !write(6,*) "y_inp",y_inp
  call poly_n_fit(poly_deg, x_inp, y_inp, nfit_dim, start, coeff)
  
  ! frequency in cm-1 and omega in hartrees/hbar, x0,y0
  freq =  hbar*sqrt(2*coeff(2)/my_SI)*cm
  omega_h = sqrt(2*coeff(2)/my_SI)/hartree
  x0 = coeff(1)*1.0d10
  y0 = coeff(poly_deg)/hartree
 

  !write (6,*) "harmonic freq", freq
  !write (6,*) "coeff2", coeff(2)
  !write (6,*) "omega_h", omega_h

  
  ! test perturbation theory
  call perturb(coeff(2)*1.d-20/hartree, coeff(3)*1.d-30/hartree, &
       coeff(4)*1.d-40/hartree, my, freq_per, corr2, corr3, corr4)

  write(6,*) "freq perturb", freq_per
  write(6,*) "harmonic", corr2
  write(6,*) "third ",corr3
  write(6,*) "fourth",corr4


  ! set up momentum part of the hamiltonian
  !*********************************************************************  

  !allocate
  allocate( Hmatrix(nstates, nstates) , tmp(npoints))

  Hmatrix=0

  ! momentum part, < m |p^2/(2m)  | n >, watch out for the fence post!

  ! diagonal part
  do i=1,nstates
     Hmatrix(i,i) = (2*i-1)*omega_h*hbar/4.
  end do

  ! non diagonal part
  do i=1,nstates-2 
     Hmatrix(i,i+2) = -sqrt(dfloat(i+1))*sqrt(dfloat(i))*omega_h*hbar/4.
  end do
  
  
  ! alt1=1, do numerical integration
  !*********************************************************************  
  
  if(alt1.eq.1) then
     
     ! calculate nstate eigenfunctions on the interval [x_new_min, x_new_max] with npoints 
     
     !allocate
     allocate( eigenfun(npoints,nstates) )
     write(6,*) "so far 2"
     
     dx = (x_new(2)-x_new(1))*1.0e10
    
     ! make eigenfunctions
     do i=1,nstates
        do j=1,npoints
           eigenfun(j,i) = harm_eigenfun(i-1, freq, my, x_new(j)*1.0e10 - x0)
        end do
     end do
     
     
     ! calculate matrix elements Hmatrix =  p^2/(2m) + V(x) 
     ! Hmatrix_mn = < m | Hmatrix | n >
     
     ! potential energy part, trapezoidal rule. Only upper triangular part of matrix
     do i=1,nstates
        do j= i,nstates
           
           tmp = eigenfun(:,i)*eigenfun(:,j)*y_new/hartree
           Hmatrix(i,j) =  Hmatrix(i,j) + ( sum( tmp(2:npoints-1) ) + 0.5*(tmp(1) &
                + tmp(npoints)) )*dx
           
        end do
     end do
     
  end if
  
  ! alt1=2, do variational calculation on polynomial surface
  !*********************************************************************  
  
  if(alt1.eq.2) then
     
     allocate( Htmp(nstates,nstates) )
     
     do i=2,poly_deg
        
        ! matrix elements
        call mat_elem_xn(Htmp,nstates,i)
        
        !allt i SI
        coeff_tmp = coeff(i)/hartree
        
        Hmatrix = Hmatrix + coeff_tmp*(hbar/(2*my_SI*omega_h*hartree))**(dfloat(i)/2.d0)*Htmp

     end do
  end if



     ! solve eigenvalue problem
    !*********************************************************************    
 
     ! allocate workvectors for lapack
     allocate( W(nstates), WORK(3*nstates))
     
     ! call lapack routine
     call dsyev( "V", "U", nstates, Hmatrix, nstates, W, WORK, 3*nstates, INFO)
     
     ! print eigenvalues and first frequency
     write(6,*) "fundamental frequency", (W(2)-W(1))*cm*hartree
     
     ! print the first eigenfunction(s)
     
     !  allocate( eigenfun1(npoints) )
     
     !  do j=1,10
     !     
     !     eigenfun1=0
     !     do i=1,nstates
     !        eigenfun1 = eigenfun1 + Hmatrix(j,i)*eigenfun(:,i)
     !     end do
     !     
     !     do i=1,npoints
     !        write(10+j,*) x_new(i), eigenfun1(i), eigenfun(i,j)
     !     end do
     !     
     !  end do


end program vib_harm
