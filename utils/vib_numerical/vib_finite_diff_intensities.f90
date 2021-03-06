program vib_finite_diff
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile, dipolefile, polfile
  integer::nin,npoints,nstates,alt1, alt2, alt3
  real(kind=wp), dimension(:),allocatable:: x,y, dip_e, dip_e_new,&
       pol_e, pol_e_new, pol, pol_new
  real(kind=wp), dimension(:,:),allocatable:: dip, dip_new
  real(kind=wp)::  my_in, dmy_x, dmy_y, dmy_z, dalpha
  
  ! loop variables
  integer::i,j

  ! other variables
  real(kind=wp), dimension(:),allocatable:: x_SI,y_SI,x_new_SI,y_new_SI, &
       dip_e_SI,pol_e_SI,diag,subdiag
  real(kind=wp)::dx, ddip_e, dpol_e, my, my_SI
  real(kind=wp):: mat_elem_q, mat_elem_pol, sigma_my, sigma_alpha
  real(kind=wp),dimension(3):: mat_elem_dip
  
  
  ! lapack
  real(kind=wp):: abstol
  integer::n_eigfound,info
  integer, dimension(:),allocatable::iwork,ifail
  real(kind=wp), dimension(:),allocatable:: eigenval,work
  real(kind=wp), dimension(:,:),allocatable:: eigenvec

  !functions
  real(kind=wp)::dlamch

  ! this version calculates IR and Raman intensities
  !
  ! al1=1 it computes matrix elements < n | my_x | m > for x,y,z from dipole 
  ! surface file, as well as IR intensities 
  ! al1=2 computes IR intesities by  (d my_x / dq ) < n | q | m > 
  !
  ! alt2=1 it computes matrix elements < n | alpha | m > where alpha is the isotropic
  ! polarizability read from file, as well as Raman intensities
  ! alt2=2 computes Raman intensities as  (d alpha / dq ) < n | q | m >
  
  ! my_in is the diatomic reduced mass


  ! read from standard input
  read(5,*) inputfile
  read(5,*) outputfile
  read(5,*) alt1, dmy_x, dmy_y, dmy_z, dipolefile
  read(5,*) alt2, dalpha, polfile 
  read(5,*) nin, npoints, nstates
  read(5,*) alt3, my_in


  ! allocate
  allocate(iwork(5*npoints), work(5*npoints), eigenval(npoints), &
       eigenvec(npoints,nstates),ifail(npoints))
  allocate( x(nin),y(nin), x_SI(nin),y_SI(nin),x_new_SI(npoints), &
       y_new_SI(npoints),diag(npoints),subdiag(npoints-1))

  ! read inputfile
  open(10,file=inputfile,status='old')
  do i=1,nin
     read(10,*) x(i),y(i)
  enddo
  close(10)

  ! read dipole file
  if(alt1.eq.1) then
     allocate(dip_e(nin), dip(nin,3), dip_e_new(npoints), dip_new(npoints,3),&
          dip_e_SI(nin))
     open(10,file=dipolefile,status='old')
     do i=1,nin
        read(10,*) dip_e(i), dip(i,1), dip(i,2), dip(i,3)
     enddo
     dip_e_SI = dip_e * 1.d-10
     close(10)
  end if

  ! read polariability file
  if(alt2.eq.1) then
     allocate(pol_e(nin), pol(nin), pol_e_new(npoints), pol_new(npoints),&
          pol_e_SI(nin))
     open(10,file=polfile,status='old')
     do i=1,nin
        read(10,*) pol_e(i), pol(i)
     enddo
     pol_e_SI = pol_e * 1.d-10
     close(10)
  end if

  ! my
  my= 18._wp/19._wp
  if(alt3.eq.1) my=my_in
  my_SI = my*amu


  ! go to SI units
  y = y -minval(y)
  y_SI=y*hartree
  x_SI=x*1.0d-10
     
  ! spline, make more points
  call linspace(x_new_SI, x_SI(1), x_SI(nin), npoints)
  dx = x_new_SI(2)-x_new_SI(1) 
  call spline_easy(x_SI, y_SI, nin, x_new_SI, y_new_SI, npoints)

  ! s�tt upp hamiltonianen -hbar^2/(2*m) (d/dx)^2 + V(x) 
  
  ! diagonalen
  do i=1,npoints
     diag(i) = (-hbar**2/(dx**2*2*my_SI) )*(-2.0_wp) + y_new_SI(i)
  end do

  ! subdiagonalen
  do i=1,npoints-1
     subdiag(i) = -hbar**2/(dx**2 *2.0_wp *my_SI) 
  end do

  !subdiag(npoints)=0
  
  ! solve eigenvalue problem
  abstol=2d0*dlamch('s')
  call dstevx('v','i',npoints,diag,subdiag, 0.d0,1.d0, 1, nstates, abstol, &
       n_eigfound, eigenval , eigenvec, npoints,work,iwork,ifail,info)


!
! Calculate matrix elements and intensities
!

  !normalize eigenvectors
  do i=1,nstates
     eigenvec(:,i) = eigenvec(:,i)/(sqrt(sum(eigenvec(:,i)**2)*dx)) 
  end do

  write(6,*) "sum of eigenvector 1,2", sqrt(sum(eigenvec(:,1)**2)*dx),&
       sqrt(sum(eigenvec(:,2)**2)*dx)


  ! calculate matrix element < 0 | q | 1 >
  mat_elem_q =  sum(eigenvec(:,1) * x_new_SI * eigenvec(:,2)) * dx
 

  ! dipole from file
  if(alt1.eq.1) then
     ! spline dipole surface in the same way as potential energy surface
     call linspace(dip_e_new, dip_e_SI(1) , dip_e_SI(nin), npoints)
     ddip_e = dip_e_new(2) - dip_e_new(1)
     call spline_easy(dip_e_SI, dip(:,1), nin, dip_e_new, dip_new(:,1), npoints)
     call spline_easy(dip_e_SI, dip(:,2), nin, dip_e_new, dip_new(:,2), npoints)
     call spline_easy(dip_e_SI, dip(:,3), nin, dip_e_new, dip_new(:,3), npoints)
     
     mat_elem_dip(1) =  sum(eigenvec(:,1) * dip_new(:,1) * eigenvec(:,2)) * ddip_e
     mat_elem_dip(2) =  sum(eigenvec(:,1) * dip_new(:,2) * eigenvec(:,2)) * ddip_e
     mat_elem_dip(3) =  sum(eigenvec(:,1) * dip_new(:,3) * eigenvec(:,2)) * ddip_e
     
  else 
     ! matrix elements (d my_x / dq ) < n | q | m > 
     mat_elem_dip(1) = dmy_x * 1.d10 * mat_elem_q
     mat_elem_dip(2) = dmy_y * 1.d10 * mat_elem_q
     mat_elem_dip(3) = dmy_z * 1.d10 * mat_elem_q
  end if

  ! polarizability from file
  if(alt2.eq.1) then
     ! spline polarizability surface in the same way as potential energy surface
     call linspace(pol_e_new, pol_e_SI(1), pol_e_SI(nin), npoints)
     dpol_e = pol_e_new(2) - pol_e_new(1)
     call spline_easy(pol_e_SI, pol, nin, pol_e_new, pol_new, npoints)
     
     mat_elem_pol =  sum(eigenvec(:,1) * pol_new * eigenvec(:,2)) * dpol_e
     
  else
     ! matrix elements (d alpha / dq ) < n | q | m > 
     mat_elem_pol = dalpha * 1.d10 * mat_elem_q
  end if

  ! cross sections
  sigma_my = sum(mat_elem_dip**2) * (eigenval(2)-eigenval(1))
  sigma_alpha = mat_elem_pol**2 * (eigenval(2)-eigenval(1))


  !
  ! write output
  !


  !open outputfile
  open(10,file=outputfile,status='unknown')
  
  ! write potential and eigenvectors  
   write(10,*) "nstates",  nstates
   write(10,*) "npoints",  npoints
   write(10,*) 
   write(10,*) "energies [J]"
    do i=1,nstates
      write(10,*) i, eigenval(i)
   end do
   write(10,*) "potential"
   do i=1,npoints
      write(10,*) x_new_SI(i), y_new_SI(i)
   end do

   do j=1, nstates
      write(10,*) "state",  j
      do i=1,npoints
         write(10,*) x_new_SI(i), eigenvec(i,j)
      end do
   end do

  close(10)
  
  ! write eigenvalues
  write(6,*) (eigenval(i)*cm, i=1,nstates) 

  ! write fundamental frequency
  write(6,*)
  write(6,*) "Fundamental frequency ", (eigenval(2)-eigenval(1))*cm
  write(6,*) "IR intensity", sigma_my
  write(6,*) "Raman intensity", sigma_alpha
  
end program vib_finite_diff
