  program vib_normal3
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile
  integer::npoints,alt1,alt2, pdeg,pstart, maxdeg_inp,alt_state, pstartdip, alt_dip
  integer::pdegdip
  real(kind=wp)::c0
  real(kind=wp), dimension(3)::dip0
  real(kind=wp), dimension(:,:),allocatable:: q,dip
  real(kind=wp), dimension(:),allocatable:: V
  real(kind=wp), dimension(3)::hfreq,my
  integer, dimension(3)::nstates
  
  ! loop variables
  integer::i,j,k,l,m ,n, deg, point, ind1, ind2, count, ii, deg1, deg2, ind, i1,i2,i3
  integer, dimension(2):: ind_1, ind_2, ind_3
  
  ! other variables
  integer::nstates3, pot1,pot2,pot3,ncols, pot_tot,maxdeg
  integer::ncolsdip
  !integer, dimension(:),allocatable:: q1_vec, q2_vec, q3_vec
  !integer, dimension(:),allocatable:: q1_dip, q2_dip, q3_dip
  integer, dimension(:,:),allocatable:: map, map_dip1, map_dip2, map_dip3
  real(kind=wp), dimension(:),allocatable:: V_SI, sol, rhs, sol_dip1, sol_dip2, sol_dip3 
  real(kind=wp), dimension(:),allocatable:: intensity, intensity_1, intensity_2,intensity_3
  real(kind=wp), dimension(:,:),allocatable:: p_elem1, p_elem2, p_elem3
  real(kind=wp), dimension(:,:),allocatable:: Hmatrix, Htmp,q_SI
  integer, dimension(:,:),allocatable::ind_to_ijk
  real(kind=wp), dimension(:,:,:),allocatable:: Mx1, Mx2, Mx3
  real(kind=wp), dimension(3):: my_SI, coeff2_h, coeff2_SI, omega_h, omega_SI 
  real(kind=wp):: Hupdate, coeff_first, coeff_exc

  ! test varibles
  real(kind=wp), dimension(:,:),allocatable:: tmp1,tmp2,tmp3
  real(kind=wp), dimension(:),allocatable:: poly1,poly2,poly3
  real(kind=wp):: tol, poly1_V, poly2_V, poly3_V, poly_x,tmp
  real(kind=wp), dimension(:,:),allocatable:: Pmatrix, Vmatrix,Hmatrix_harm,&
       Hmatrix_2, Hmatrix_3,Hmatrix_4 ,Hmatrix_numerical


  ! LAPACK variables
  integer::INFO, LWORK
  real(kind=wp), dimension(:),allocatable::W,WORK 
  real(kind=wp), dimension(:,:),allocatable::Hmatrix_lapack

  !SVDCPM variables
  real(kind=wp), dimension(:,:),allocatable::svd_A, svd_V, svd_W_mat, Hmatrix_svd 
  real(kind=wp), dimension(:),allocatable::svd_W
  real(kind=wp):: svd_tol

  !write wavefunction variables
  integer:: ngridpoints
  real(kind=wp):: maxq,minq,dx,psi2,coeff,fun1,fun2,fun3,q1_tmp,q2_tmp,q3_tmp

  !functions
  integer::delta
  real(kind=wp):: harm_eigenfun

  ! also implements intensities

  ! read from standard input
  read(5,*) inputfile
  read(5,*) npoints, pdeg, alt1, alt2, alt_dip
  read(5,*)  pdegdip
  read(5,*) nstates(1),nstates(2), nstates(3), maxdeg_inp, alt_state
  read(5,*) hfreq(1), hfreq(2), hfreq(3)
  read(5,*) my(1), my(2), my(3)
  read(5,*) c0, dip0(1), dip0(2), dip0(3)

  ! alt1=0, use analytical omegas, al1=1, use omegas from polynomial fitting
  ! alt2=1 also fit x,y,x terms
  ! alt_state =0, do nothing, =1,2,3 only use mode 1,2, or three, set all other non-diagonal matrix elements=0
  ! maxdeg= max number of exitations
  ! alt_dip=0 don't read dipole surface , no intensities

  write(6,*) hfreq
  
  ! potential energy surface 
  if(alt1.eq.0) then
   
     ! use analytical omegas, fit from 3 degree
     ! these numbers are actually sum(multinomial(3,r), r=3,pdeg)
     ! maybe this function should be defined?
     pstart=3
     !if(pdeg.eq.2) ncols=0
     if(pdeg.eq.3) ncols=10
     if(pdeg.eq.4) ncols=25
     if(pdeg.eq.5) ncols=46
     if(pdeg.eq.6) ncols=74
     if(pdeg.eq.7) ncols=110
     if(pdeg.eq.8) ncols=155
     if(pdeg.eq.10) ncols=276
     if(pdeg.eq.12) ncols=445
     if(pdeg.eq.13) ncols=550
     if(pdeg.eq.14) ncols=670

  else if(alt1.eq.1) then
     
     ! fit second degree terms
     if(pdeg.eq.2) ncols=6
     if(pdeg.eq.3) ncols=16
     if(pdeg.eq.4) ncols=31
     if(pdeg.eq.5) ncols=52
     if(pdeg.eq.6) ncols=80
     if(pdeg.eq.7) ncols=116
     if(pdeg.eq.8) ncols=161
     if(pdeg.eq.10) ncols=282
     if(pdeg.eq.12) ncols=451
     if(pdeg.eq.13) ncols=556
     if(pdeg.eq.14) ncols=676
 
     if(alt2.eq.1) then

        ! also fit x,y,z terms
        ncols = ncols +3
        pstart=1

     else if(alt2.eq.0) then
        ! dont fit x,y,z terms
        pstart=2
     end if

  end if


! dipole surface
if(alt_dip.eq.1) then
  if(pdegdip.eq.1) ncolsdip=3
  if(pdegdip.eq.2) ncolsdip=6 + 3
  if(pdegdip.eq.3) ncolsdip=16 + 3
  if(pdegdip.eq.4) ncolsdip=31 + 3
  if(pdegdip.eq.5) ncolsdip=52 + 3
  if(pdegdip.eq.6) ncolsdip=80 + 3
  
  ! allocate
  allocate( q(3,npoints), q_SI(3,npoints), V(npoints), V_SI(npoints), dip(3,npoints))
  allocate( rhs(npoints))


  ! read file with values of q1, q2, q3 [Å] and V [hartree], dip
  open(10,file=inputfile,status='old')
  
  do i=1,npoints
     read(10,*) q(1,i), q(2,i), q(3,i), V(i), dip(1,i), dip(2,i), dip(3,i)
  enddo

  close(10)

else if(alt_dip.eq.0) then
   ! allocate
   allocate( q(3,npoints), q_SI(3,npoints), V(npoints), V_SI(npoints))
   allocate( rhs(npoints))
   
   ! read file with values of q1, q2, q3 [Å] and V [hartree], dip
   open(10,file=inputfile,status='old')
   
  do i=1,npoints
     read(10,*) q(1,i), q(2,i), q(3,i), V(i)
  enddo

  close(10)

end if


  ! my
  my_SI = my*amu

  ! maxdeg
  maxdeg =  maxdeg_inp+3


  ! calculate nstates3
  ind=0
  do i=1,nstates(1)
     do j=1,nstates(2)
        do k=1,nstates(3)
           if(i+j+k <= maxdeg) then
              ind=ind+1
           end if
        end do
     end do
  end do

  nstates3= ind

  ! go to SI units
  V_SI=V*hartree
  q_SI=q*1.0e-10
  coeff2_SI = ((hfreq/(hbar*cm))**2)*0.5*my_SI
  coeff2_h = coeff2_SI*1.0d-20/hartree  !((hfreq/(hbar))**2)*0.5*my*1.0d20 
  !coeff2_SI*1.0d-20/amu ! helt fel, kolla snarast!
  omega_SI = hfreq/(hbar*cm)
  omega_h = hfreq/(hbar*cm*hartree)
  
  write(6,*) coeff2_h
  write(6,*) omega_h
  
  ! fit a polynomial arbitrary degree to potential energy surface
  ! linear least squares
  !*********************************************************************  

  ! sätter upp right hand side, hartree
  
  if(alt1.eq.0) then
     
     do point = 1,npoints
        rhs(point) = V(point) - c0 - sum(coeff2_h*(q(:,point)**2))  
     end do
     
  else if(alt1.eq.1) then
     
     do point = 1,npoints
        rhs(point) = V(point) - c0
     end do
     
  end if
  

  allocate(sol(ncols),map(3,ncols))

  !write(6,*) "1", npoints,ncols
  write(6,*) npoints,pstart, pdeg, ncols, size(q), size(rhs), size(map), size(sol)
  call fit_poly_surface(npoints, pstart, pdeg, ncols, q, rhs, map, sol)

  !write(6,*) sol
  !  stop


  ! fit a polynomial arbitrary degree to dipole surface
  ! linear least squares
  !*********************************************************************  
 
  if(alt_dip.eq.1) then

  ! x direction

  ! sätter upp right hand side, dipole
  do point = 1,npoints
     rhs(point) = dip(1,point) - dip0(1)  
  end do

  allocate(sol_dip1(ncolsdip), map_dip1(3,ncolsdip))

  pstartdip=1  
  call fit_poly_surface(npoints, pstartdip, pdegdip, ncolsdip, q, rhs, map_dip1, sol_dip1)

  write(6,*) "dipole, x", sol_dip1
  
  ! y direction

  ! sätter upp right hand side, dipole
  do point = 1,npoints
     rhs(point) = dip(2,point) - dip0(2)  
  end do

  allocate(sol_dip2(ncolsdip), map_dip2(3,ncolsdip))

  pstartdip=1  
  call fit_poly_surface(npoints, pstartdip, pdegdip, ncolsdip, q, rhs, map_dip2, sol_dip2)

  write(6,*) "dipole, y", sol_dip2

  ! z direction

  ! sätter upp right hand side, dipole
  do point = 1,npoints
     rhs(point) = dip(3,point) - dip0(3)  
  end do
  
  allocate(sol_dip3(ncolsdip), map_dip3(3,ncolsdip))
  
  pstartdip=1  
  call fit_poly_surface(npoints, pstartdip, pdegdip, ncolsdip, q, rhs, map_dip3, sol_dip3)

  write(6,*) "dipole, z", sol_dip3

end if

  ! set up Hamiltonian
  !*********************************************************************  

  !allocate
  allocate( Hmatrix(nstates3, nstates3) , p_elem1(nstates(1),nstates(1)), p_elem2(nstates(2),nstates(2)) &
       , p_elem3(nstates(3),nstates(3)), ind_to_ijk(nstates3,3))


!write(6,*) "test************"
  Hmatrix=0


  ! basis functions | i >| j >| k > , nstates(1)*nstates(3)*nstates(3) number of them
  ! no need to precompute the matrix elements < i | p^2/2m | j >, or?

  ! Fitted omegas if alt1=1
  !****************************
  ! momentum part, < m |p^2/(2m)  | n >, watch out for the fence post!
  !precompute it 

if(alt1.eq.1) then

   ! set to zero
   p_elem1=0
   p_elem2=0
   p_elem3=0

  ! diagonal part
  do i=1,nstates(1)
     p_elem1(i,i) = (2*i-1)*hbar*omega_h(1)/4.
  end do

  do i=1,nstates(2)
     p_elem2(i,i) = (2*i-1)*hbar*omega_h(2)/4.
  end do

  do i=1,nstates(3)
     p_elem3(i,i) = (2*i-1)*hbar*omega_h(3)/4.
  end do


  !non diagonal part
  do i=1,nstates(1)-2 
     p_elem1(i,i+2) = -sqrt(dfloat(i+1))*sqrt(dfloat(i))*hbar*omega_h(1)/4.
     p_elem1(i+2,i) = p_elem1(i,i+2)
  end do

  do i=1,nstates(2)-2 
     p_elem2(i,i+2) = -sqrt(dfloat(i+1))*sqrt(dfloat(i))*hbar*omega_h(2)/4.
     p_elem2(i+2,i) = p_elem2(i,i+2)
  end do
  
  do i=1,nstates(3)-2 
     p_elem3(i,i+2) = -sqrt(dfloat(i+1))*sqrt(dfloat(i))*hbar*omega_h(3)/4.
     p_elem3(i+2,i) = p_elem3(i,i+2)
  end do

end if


! Potential energy part, both fitted and analytical omegas
!********************************

allocate( Htmp(nstates(1),nstates(1)),  Mx1(pdeg+1,nstates(1),nstates(1))  )

do i=0,pdeg
   
   call mat_elem_xn(Htmp,nstates(1),i)
   
   !Mx1(i+1,:,:) = (1.0d10)**i*(hbar/(2*my_SI(1)*omega_SI(1)))**(dfloat(i)/2.d0)*Htmp
   Mx1(i+1,:,:) = (hbar/(2*my_SI(1)*omega_SI(1)) )**(dfloat(i)/2.d0)*Htmp

end do
 

  deallocate(Htmp)
  allocate( Htmp(nstates(2),nstates(2)), Mx2(pdeg+1,nstates(2),nstates(2)))
  
  do i=0,pdeg
     
     call mat_elem_xn(Htmp,nstates(2),i)
   
        Mx2(i+1,:,:) = (hbar/(2*my_SI(2)*omega_SI(2) ))**(dfloat(i)/2.d0)*Htmp
     
  end do
  
  deallocate(Htmp)
  allocate( Htmp(nstates(3),nstates(3)) , Mx3(pdeg+1,nstates(3),nstates(3)))
  
  do i=0,pdeg
     
     call mat_elem_xn(Htmp,nstates(3),i)
     
        Mx3(i+1,:,:) = (hbar/(2*my_SI(3)*omega_SI(3)) )**(dfloat(i)/2.d0)*Htmp
     
  end do

deallocate(Htmp)



! reset matrix elements if alt_state=1,2,3
!  
if(alt_state.eq.1) then
   
   Mx2(2:pdeg+1,:,:)=0
   Mx3(2:pdeg+1,:,:)=0
   
else if(alt_state.eq.2) then
   
   Mx1(2:pdeg+1,:,:)=0
   Mx3(2:pdeg+1,:,:)=0
   
else if(alt_state.eq.3) then
   
   Mx1(2:pdeg+1,:,:)=0
   Mx2(2:pdeg+1,:,:)=0
   
end if


! Set up Hamiltonian Matrix
!******************************  

!test
allocate( Pmatrix(nstates3, nstates3),Vmatrix(nstates3, nstates3), Hmatrix_harm(nstates3, nstates3),&
     Hmatrix_2(nstates3, nstates3),Hmatrix_3(nstates3, nstates3), Hmatrix_4(nstates3, nstates3))

 Hmatrix_2=0
 Hmatrix_3=0
 Hmatrix_4=0
 Pmatrix=0

 ind1=0
 !ind2=0
 do i=1,nstates(1)
    do j=1,nstates(2)
       do k=1,nstates(3)
          
          !           ind1 = k + (j-1)*nstates(3) + (i-1)*nstates(2)*nstates(3)
          
          if(i+j+k <= maxdeg) then
             
             ind1 = ind1 + 1
             
             ind2=0
             do l=1,nstates(1)
                do m=1,nstates(2)
                   do n=1,nstates(3)
                      
                      !ind2 = n + (m-1)*nstates(3) + (l-1)*nstates(2)*nstates(3)
                      
                      if(l+m+n <= maxdeg) then
                         
                         ind2 = ind2 + 1
                         
                         ! only upper triangular part of matrix 
                         if(ind1 <= ind2) then 
                            
                            if(alt1.eq.0) then
                               Hmatrix(ind1,ind2) =  ( (0.5 + i-1)*omega_h(1) +(0.5 + j-1)*omega_h(2) &
                                    + (0.5 + k-1)*omega_h(3) )*hbar*delta(i,l)*delta(j,m)*delta(k,n) 
                               
                               
                               
                            else if(alt1.eq.1) then
                               Hmatrix(ind1,ind2) = p_elem1(i,l)*delta(j,m)*delta(k,n) + p_elem2(j,m)*delta(i,l)*delta(k,n) &
                                    + p_elem3(k,n)*delta(i,l)*delta(j,m) 
                               
                               
                            end if
                            
                            do ii=1,ncols
                               pot1 = map(1,ii)+1
                               pot2 = map(2,ii)+1
                               pot3 = map(3,ii)+1
                               pot_tot = map(1,ii) + map(2,ii) + map(3,ii)
                               
                               Hupdate=  ((1.0d10)**pot_tot)*sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n) 
                               
                               Hmatrix(ind1,ind2)= Hmatrix(ind1,ind2) + Hupdate
                               
                            end do
                            
                         end if ! if(ind1 <= ind2) then 
                         
                      end if !l+m+n <= maxdeg
                      
                   end do
                end do
             end do
             
          end if !i+j+k <= maxdeg
          
       end do
    end do
 end do
 
  
 ! write(6,*) "ind1", ind1, "ind2",ind2
  
  ! mapping from ind to i,j,k
  ind=0
  do i=1,nstates(1)
     do j=1,nstates(2)
        do k=1,nstates(3)
           
           if(i+j+k <= maxdeg) then
              
              ind=ind+1
              
              ind_to_ijk(ind,1) = i
              ind_to_ijk(ind,2) = j
              ind_to_ijk(ind,3) = k
           end if
           
        end do
     end do
  end do
  
  
  ! solve eigenvalue problem
  !*********************************************************************    
   
  ! allocate workvectors for lapack
  !nstates3 = ind1
  LWORK=3*nstates3
  allocate( W(nstates3), WORK(LWORK), Hmatrix_lapack(nstates3,nstates3) )!, Hmatrixtmp(nstates3,nstates3))
  
 ! hmatrixtmp=Hmatrix(1:ind1,1:ind1)


  Hmatrix_lapack = Hmatrix

  ! call lapack routine
   call dsyev( "V", "U", nstates3, Hmatrix_lapack, nstates3, W, WORK, LWORK, INFO)
 

   ! calculate transition dipoles
   if(alt_dip.eq.1) then
      
   allocate(intensity(nstates3),intensity_1(nstates3),intensity_2(nstates3), intensity_3(nstates3))
   
   intensity(1) = 0
   ! excited state
   do i=2,nstates3
      intensity(i) =0
      
      ! expansion coefficient ith respect to basis functions
      do i2 = 1,nstates3
         ind_1(1) = ind_to_ijk(i2,1)
         ind_2(1) = ind_to_ijk(i2,2)
         ind_3(1) = ind_to_ijk(i2,3)
         
         ! excited state coefficient 
         coeff_exc = Hmatrix_lapack(i2,i)

         ! expansion coefficient ith respect to basis functions
         do i3 = 1,nstates3
            ind_1(2) = ind_to_ijk(i3,1)
            ind_2(2) = ind_to_ijk(i3,2)
            ind_3(2) = ind_to_ijk(i3,3)
            
            ! ground state coefficient
            coeff_first = Hmatrix_lapack(i3,1)
            
            do j=1,ncolsdip
               pot1 = map_dip1(1,j) +1
               pot2 = map_dip1(2,j) +1
               pot3 = map_dip1(3,j) +1
               pot_tot = map_dip1(1,j) + map_dip1(2,j) + map_dip1(3,j)

               !tmp = sol_dip1(j)*coeff_first*coeff_exc*Mx1(pot1,ind_1(1),ind_1(2))*Mx2(pot2,ind_2(1),ind_2(2))*Mx3(pot3,ind_3(1),ind_3(2))
               
               if(0) then
               if(abs(tmp)>1.0d-50) then
                  write(6,*) "i, i2,i3 ",i, i2,i3
                  write(6,*) "ind exc  ", ind_1(1),ind_2(1),ind_3(1)
                  write(6,*) "ind first",ind_1(2),ind_2(2),ind_3(2)
                  write(6,*) "pot      ", pot1, pot2, pot3
                  write(6,*) tmp
                  write(6,*)
               end if
            end if
               !write(6,*) soldip(j)*Mx1(pot1,ind_1,1)*Mx2(pot2,ind_2,1)*Mx3(pot3,ind_3,1) 
               !write(6,*) soldip(j),Mx1(pot1,ind_1,1),Mx2(pot2,ind_2,1),Mx3(pot3,ind_3,1)
               intensity_1(i) = intensity_1(i)+ &
                    ((1.0d10)**pot_tot)*sol_dip1(j)*coeff_first*coeff_exc*Mx1(pot1,ind_1(1),ind_1(2))*Mx2(pot2,ind_2(1),&
                    ind_2(2))*Mx3(pot3,ind_3(1),ind_3(2))
               
               intensity_2(i) = intensity_2(i)+ &
                    ((1.0d10)**pot_tot)*sol_dip2(j)*coeff_first*coeff_exc*Mx1(pot1,ind_1(1),ind_1(2))*Mx2(pot2,ind_2(1),&
                    ind_2(2))*Mx3(pot3,ind_3(1),ind_3(2))
               
               intensity_3(i) = intensity_3(i)+ &
                    ((1.0d10)**pot_tot)*sol_dip3(j)*coeff_first*coeff_exc*Mx1(pot1,ind_1(1),ind_1(2))*Mx2(pot2,ind_2(1),&
                    ind_2(2))*Mx3(pot3,ind_3(1),ind_3(2))
               
            end do
         end do
      end do
   end do
   
   intensity = intensity_1**2 +  intensity_2**2 + intensity_3**2 
   
end if


   !print eigenvectors and intensities
   write(6,*) "fundamental frequency, dseyv"

   do i=2,20
      write(6,*) (W(i)-W(1))*cm*hartree !, intensity(i)
      write(19,*) (W(i)-W(1))*cm*hartree !, intensity(i)
      !, intensity_1(i), intensity_2(i), intensity_3(i)
   end do




   ! write eigenvectors
   do i=1,nstates3
      write(11,*) "eigenvector ", i
      
      do j=1,nstates3
         write(11,*) ind_to_ijk(j,1), ind_to_ijk(j,2), ind_to_ijk(j,3), Hmatrix_lapack(j,i)
      end do
      
   end do



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
  

   ! print ground state wave function on grid, to be implemented soon
   ! not yet finished!
   ngridpoints = 101
   maxq = 0.3_wp
   minq = -0.3_wp
   dx= (maxq - minq)/real(ngridpoints-1)
   
   do i1 = 1, ngridpoints - 1
      do i2 = 1, ngridpoints - 1
         do i3 = 1, ngridpoints - 1
            q1_tmp = minq + (i1-1)*dx  
            q2_tmp = minq + (i2-1)*dx  
            q3_tmp = minq + (i3-1)*dx  
            
            
            psi2 = 0
            do j=1,nstates3
               coeff = Hmatrix_lapack(j,1)
               fun1 = harm_eigenfun(ind_to_ijk(j,1)-1, hfreq(1), my(1), q1_tmp)
               fun2 = harm_eigenfun(ind_to_ijk(j,2)-1, hfreq(2), my(2), q2_tmp)
               fun3 = harm_eigenfun(ind_to_ijk(j,3)-1, hfreq(3), my(3), q3_tmp)
               psi2 = psi2 + coeff*fun1*fun2*fun3
            end do

            psi2 = psi2**2
            write(20,'(4F16.8)') q1_tmp, q2_tmp, q3_tmp, psi2
      
         end do
      end do
   end do

   !write cuts of wavefunction squared
   ! coord 1
   do i1 = 1, ngridpoints - 1
       q1_tmp = minq + (i1-1)*dx 
       q2_tmp = 0
       q3_tmp = 0

       psi2 = 0
       do j=1,nstates3
          coeff = Hmatrix_lapack(j,1)**2
          fun1 = harm_eigenfun(ind_to_ijk(j,1)-1, hfreq(1), my(1), q1_tmp)**2
          fun2 = harm_eigenfun(ind_to_ijk(j,2)-1, hfreq(2), my(2), q2_tmp)**2
          fun3 = harm_eigenfun(ind_to_ijk(j,3)-1, hfreq(3), my(3), q3_tmp)**2
          psi2 = psi2 + coeff*fun1*fun2*fun3
       end do

       write(21,'(4F16.8)') q1_tmp, q2_tmp, q3_tmp, psi2
    end do

    ! coord 2
    do i1 = 1, ngridpoints - 1
       q1_tmp = 0
       q2_tmp =  minq + (i1-1)*dx 
       q3_tmp = 0

       psi2 = 0
       do j=1,nstates3
          coeff = Hmatrix_lapack(j,1)**2
          fun1 = harm_eigenfun(ind_to_ijk(j,1)-1, hfreq(1), my(1), q1_tmp)**2
          fun2 = harm_eigenfun(ind_to_ijk(j,2)-1, hfreq(2), my(2), q2_tmp)**2
          fun3 = harm_eigenfun(ind_to_ijk(j,3)-1, hfreq(3), my(3), q3_tmp)**2
          psi2 = psi2 + coeff*fun1*fun2*fun3
       end do

       write(22,'(4F16.8)') q1_tmp, q2_tmp, q3_tmp, psi2
    end do

    ! coord 3
    do i1 = 1, ngridpoints - 1
       q1_tmp = 0
       q2_tmp = 0 
       q3_tmp = minq + (i1-1)*dx 

       psi2 = 0
       do j=1,nstates3
          coeff = Hmatrix_lapack(j,1)**2
          fun1 = harm_eigenfun(ind_to_ijk(j,1)-1, hfreq(1), my(1), q1_tmp)**2
          fun2 = harm_eigenfun(ind_to_ijk(j,2)-1, hfreq(2), my(2), q2_tmp)**2
          fun3 = harm_eigenfun(ind_to_ijk(j,3)-1, hfreq(3), my(3), q3_tmp)**2
          psi2 = psi2 + coeff*fun1*fun2*fun3
       end do

       write(23,'(4F16.8)') q1_tmp, q2_tmp, q3_tmp, psi2
    end do


! print intensities along normal modes
  if(alt_dip.eq.1) then
     
     tol=1.d-10
     allocate(tmp1(22,4),tmp2(22,4),tmp3(22,4) )
     
     j=1
     k=1
     l=1

     !write(6,*) "test 1"
     do i=1,npoints
        
        if(abs(q(2,i))<tol.and.abs(q(3,i))<tol) then
           tmp1(j,1) = q(1,i)
           tmp1(j,2)= dip(1,i)
           tmp1(j,3)= dip(2,i)
           tmp1(j,4)= dip(3,i)
           j=j+1
           !write(6,*) q(2,i),q(3,i)
        end if
        if(abs(q(1,i))<tol.and.abs(q(3,i))<tol) then
           tmp2(k,1) = q(2,i)
           tmp2(k,2)= dip(1,i)
           tmp2(k,3)= dip(2,i)
           tmp2(k,4)= dip(3,i)
           k=k+1
        end if
        if(abs(q(1,i))<tol.and.abs(q(2,i))<tol) then
           tmp3(l,1) = q(3,i)
           tmp3(l,2)= dip(1,i)
           tmp3(l,3)= dip(2,i)
           tmp3(l,4)= dip(3,i)
           l=l+1
        end if
        
     enddo
     
     do i=1,j
        write(16,'(4F16.8)')  tmp1(i,:)
        write(17,'(4F16.8)')  tmp2(i,:)
        write(18,'(4F16.8)')  tmp3(i,:)
     end do
     
  end if

end program vib_normal3
