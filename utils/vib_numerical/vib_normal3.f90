program vib_normal3
  use parameters
  implicit none

  ! input/output
  character(75)::inputfile,outputfile
  integer::npoints,alt1,alt2, pdeg,pstart, maxdeg_inp,alt_state
  real(kind=wp)::c0
  real(kind=wp), dimension(:,:),allocatable:: q
  real(kind=wp), dimension(:),allocatable:: V
  real(kind=wp), dimension(3)::hfreq,my
  integer, dimension(3)::nstates
  
  ! loop variables
  integer::i,j,k,l,m ,n, deg, point, ind1, ind2, count, ii, deg1, deg2, ind
  
  ! other variables
  integer::nstates3, pot1,pot2,pot3,ncols, pot_tot,maxdeg
  integer, dimension(:),allocatable:: q1_vec, q2_vec, q3_vec
  real(kind=wp), dimension(:),allocatable:: V_SI, sol, rhs
  real(kind=wp), dimension(:,:),allocatable:: p_elem1, p_elem2, p_elem3
  real(kind=wp), dimension(:,:),allocatable:: Hmatrix, Htmp,q_SI, matrix
  integer, dimension(:,:),allocatable::ind_to_ijk
  real(kind=wp), dimension(:,:,:),allocatable:: Mx1, Mx2, Mx3
  real(kind=wp), dimension(3):: my_SI, coeff2_h, coeff2_SI, omega_h, omega_SI 
  real(kind=wp):: Hupdate

  ! test varibles
  real(kind=wp), dimension(:,:),allocatable:: tmp1,tmp2,tmp3
  real(kind=wp), dimension(:),allocatable:: poly1,poly2,poly3
  real(kind=wp):: tol, poly1_V, poly2_V, poly3_V, poly_x
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

  !functions
  integer::delta

     !         write(11,*) i,j,k
  !real(kind=wp):: harm_eigenfun
  ! real(kind=wp)::dlamch

  ! read from standard input
  read(5,*) inputfile
  !read(5,*) outputfile
  !read(5,*) alt1, alt2
  !read(5,*) nin, npoints, nstates, nfit_p, nfit_s
  !read(5,*) x_new_min, x_new_max
  read(5,*) npoints, pdeg, alt1, alt2
  read(5,*) nstates(1),nstates(2), nstates(3), maxdeg_inp, alt_state
  read(5,*) hfreq(1), hfreq(2), hfreq(3)
  read(5,*) my(1), my(2), my(3)
  read(5,*) c0
  !read(5,*) svd_tol

  ! alt1=0, use analytical omegas, al1=1, use omegas from polynomial fitting
  ! alt2=1 also fit x,y,x terms
  ! alt_state =0, do nothing, =1,2,3 only use mode 1,2, or three, set all other non-diagonal matrix elements=0
  ! maxdeg= max number of exitations


write(6,*) hfreq
  
  if(alt1.eq.0) then
   
     ! use analytical omegas, fit from 3 degree
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


!ncols=25 !6!52



  ! allocate
  !allocate( x(nin),y(nin), x_SI(npoints),y_SI(npoints),x_new(npoints), &
  !     y_new(npoints))

!  write(6,*) "test************"

  ! allocate
  allocate( q(3,npoints), q_SI(3,npoints), V(npoints), V_SI(npoints))
!write(6,*) "test************"
  allocate(q1_vec(ncols), q2_vec(ncols), q3_vec(ncols), matrix(npoints,ncols) )
!write(6,*) "test************"
  allocate( rhs(npoints), sol(npoints))
!write(6,*) "test************"
  ! read file with values of q1, q2, q3 [�] and V, hartree
  open(10,file=inputfile,status='old')
  
  
  do i=1,npoints
     read(10,*) q(1,i), q(2,i), q(3,i), V(i)
  enddo

  close(10)

 !!write(6,*) "1"

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

  ! nstates3
  !nstates3 = nstates(1)* nstates(2)* nstates(3)
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
  
  ! fit a sixth degree polynomial, linear least squares
  !*********************************************************************  
  
  ! s�tter upp matrisen med coefficienter, q [�]
   
  do point=1,npoints
     count=0
     
     do deg=pstart,pdeg
        ! loopar igenom alla n,m,l med summa deg
        do i=0,deg
           do j=0,deg-i

              k=deg-(i+j)
              count = count +1
              
              if(point.eq.1) then
                 q1_vec(count) = i
                 q2_vec(count) = j
                 q3_vec(count) = k                 
              end if


              matrix(point,count)= (q(1,point)**i)*(q(2,point)**j)*(q(3,point)**k)

     !         write(11,*) i,j,k
     !         write(11,*) q(1,point),q(2,point),q(3,point)
     !         write(11,*) (q(1,point)**i)*(q(2,point)**j)*(q(3,point)**k) 
              !write(6,*) "count", count
           end do
        end do
     end do
  end do

write(6,*) "count", count

!do i=1,npoints
!   write(11,*) "row number", i
!   write(11,*) matrix(i,:)
!end do

  ! s�tter upp right hand side, hartree

if(alt1.eq.0) then
  
   do point = 1,npoints
     rhs(point) = V(point) - c0 - sum(coeff2_h*(q(:,point)**2))  
  end do

else if(alt1.eq.1) then
   
   do point = 1,npoints
      rhs(point) = V(point) - c0
   end do
   
end if

!do i=1,count
!write(6,*)
!     write(6,*) q1_vec(i),q2_vec(i),q3_vec(i)
!     write(6,*) matrix(:,i)
!end do
  
 !write(6,*) coeff2_h
 !write(6,*) q(:,1)**2


!write(6,*) q1_vec(1),q2_vec(1),q3_vec(1)
!write(6,*) q(1,1), q(2,1),q(3,1)
!write(6,*) (q(1,1)**0)*(q(2,1)**0)*(q(3,1)**3)
!write(6,*) V
!write(6,*)
!write(6,*) c0
!write(6,*) sum(coeff2_h*(q(:,1)**2))

!write(6,*) 
!write(6,*)

!write(6,*) rhs



  ! solve the  lsq problem 
  !*************************************
 
  sol =rhs

  LWORK = npoints*ncols +100
  allocate( WORK(LWORK))
!write(6,*) "test************"
  ! DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
  call  DGELS( 'N', npoints, ncols, 1, matrix, npoints, sol, npoints, WORK, LWORK, INFO )

  deallocate( WORK )

  !allocate( svd_W(ncols), svd_V(ncols,ncols))
  !CALL SVDCMP(matrix,npoints,ncols,npoints,ncols,svd_W,svd_V)

  !write(6,*) "svd_W"
  !write(6,*) "*************************"
  !write(6,*) svd_W
  !write(6,*) "*************************"

  !deallocate( svd_W, svd_V)
 
  

write(12,*) INFO
write(12,*) sol(1:ncols)
write(12,*)

!write(6,*) "sol" , hbar*cm*sqrt(hartree*1.0d20*2*sol(1)/my_SI(3)),&
!     hbar*cm*sqrt(hartree*1.0d20*2*sol(3)/my_SI(2)),&
!     hbar*cm*sqrt(hartree*1.0d20*2*sol(6)/my_SI(1))


! coeff2_SI = ((hfreq/(hbar*cm))**2)*0.5*my_SI
!  coeff2_h = coeff2_SI*1.0d-20/hartree

! hfreq = hbar*cm*sqrt(hartree*1.0d20*2*coeff2_h/(my_SI))


! write cuts of PES to file
tol=1.d-10
allocate(tmp1(22,2),tmp2(22,2),tmp3(22,2) )

j=1
k=1
l=1

!write(6,*) "test 1"
do i=1,npoints

   if(abs(q(2,i))<tol.and.abs(q(3,i))<tol) then
      tmp1(j,1) = q(1,i)
      tmp1(j,2)= V(i)
      j=j+1
      !write(6,*) q(2,i),q(3,i)
   end if
   if(abs(q(1,i))<tol.and.abs(q(3,i))<tol) then
      tmp2(k,1) = q(2,i)
      tmp2(k,2)= V(i)
      k=k+1
   end if
   if(abs(q(1,i))<tol.and.abs(q(2,i))<tol) then
      tmp3(l,1) = q(3,i)
      tmp3(l,2)= V(i)
      l=l+1
   end if
   
enddo

!write(6,*) i,j,k,l
!write(6,*) "test 2"
!do i=1,17
!   write(16,*)  tmp1(i,:)
!   write(17,*)  tmp2(i,:)
!   write(18,*)  tmp3(i,:)
!end do
!write(6,*) "test 3"


deallocate(tmp1,tmp2,tmp3 )

! write corresponding polynomial fits to file
allocate(poly1(pdeg),poly2(pdeg),poly3(pdeg))
!write(6,*) "test************"

poly1=0
poly2=0
poly3=0

  do ii=1,ncols
     if(q2_vec(ii).eq.0.and.q3_vec(ii).eq.0) then
        poly1(q1_vec(ii)) = sol(ii)
     end if
     if(q1_vec(ii).eq.0.and.q3_vec(ii).eq.0) then
        poly2(q2_vec(ii)) = sol(ii)
     end if
     if(q1_vec(ii).eq.0.and.q2_vec(ii).eq.0) then
        poly3(q3_vec(ii)) = sol(ii)
     end if
  end do

!write(6,*) "poly1"
!write(6,*) poly1
!write(6,*) "test************2"

  do i=1,301

     poly1_V=0
     poly2_V=0
     poly3_V=0

!write(6,*) "test************2.1"
     poly_x = -1.5_wp + (i-1)*0.01_wp

     do j=pstart, pdeg

        poly1_V = poly1_V + poly1(j)*(poly_x)**j
        poly2_V = poly2_V + poly2(j)*(poly_x)**j
        poly3_V = poly3_V + poly3(j)*(poly_x)**j

     end do
!write(6,*) "test************2.2"
!      write(6,*)  poly_x,poly1_V + c0 + coeff2_h(1)*(poly_x)**2
!      write(6,*) poly1_V + c0 + coeff2_h(2)*(poly_x)**2
!      write(6,*) poly1_V + c0 + coeff2_h(1)*(poly_x)**2
!write(6,*) "test************2.2.1"
     write(19,*) poly_x, poly1_V + c0 + coeff2_h(1)*(poly_x)**2

     write(20,*) poly_x, poly2_V + c0 + coeff2_h(2)*(poly_x)**2
     write(21,*) poly_x, poly3_V + c0 + coeff2_h(3)*(poly_x)**2
!write(6,*) "test************2.3"
  end do

!write(6,*) "test************3"
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


!   
!   do i=1,nstates(2)
!      Mx2(1,i,i)=1
!   end do
!   
!   do i=1,nstates(3)
!      Mx3(1,i,i)=1
!   end do
!
!
!else if(alt_state.eq.2) then
!
!   Mx1=0
!   Mx3=0
!   
!   do i=1,nstates(1)
!      Mx1(1,i,i)=1
!   end do
!   
!   do i=1,nstates(3)
!      Mx3(1,i,i)=1
!   end do
!
!else if(alt_state.eq.3) then
! 
!   Mx1=0
!   Mx2=0
!   
!   do i=1,nstates(1)
!      Mx1(1,i,i)=1
!   end do
!   
!   do i=1,nstates(2)
!      Mx2(1,i,i)=1
!   end do!!
!!!
!!
!
!end if





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
                     
                          !test
                          Hmatrix_harm(ind1,ind2) = ( (0.5 + i-1)*omega_h(1) +(0.5 + j-1)*omega_h(2) &
                            + (0.5 + k-1)*omega_h(3) )*hbar*delta(i,l)*delta(j,m)*delta(k,n) 


                       else if(alt1.eq.1) then
                          Hmatrix(ind1,ind2) = p_elem1(i,l)*delta(j,m)*delta(k,n) + p_elem2(j,m)*delta(i,l)*delta(k,n) &
                               + p_elem3(k,n)*delta(i,l)*delta(j,m) 

                          !test
                          Pmatrix(ind1,ind2) = p_elem1(i,l)*delta(j,m)*delta(k,n) + p_elem2(j,m)*delta(i,l)*delta(k,n) &
                               + p_elem3(k,n)*delta(i,l)*delta(j,m) 
                          
                       end if

                       
                       do ii=1,ncols
                          pot1 = q1_vec(ii)+1
                          pot2 = q2_vec(ii)+1
                          pot3 = q3_vec(ii)+1
                          
                          pot_tot = q1_vec(ii) + q2_vec(ii) + q3_vec(ii)
                          !write(6,*) "ijklmn", i,j,k,l,m,n
                          !write(6,*) "pot",pot1,pot2,pot3
                          !write(6,*) Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)
                          !write(6,*) sol(ii)
                          !write(6,*) sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)
                  
                          !Hmatrix(ind1,ind2)= Hmatrix(ind1,ind2) + &
                          !     ((1.0d10)**pot_tot)*sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)
                          

                          !if alt_state=1,2,3 force matrix elements to be 0

                          Hupdate=  ((1.0d10)**pot_tot)*sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)                          


                          !test
                          if(pot_tot.eq.2) then
                             Hmatrix_2(ind1,ind2) = Hmatrix_2(ind1,ind2) + ((1.0d10)**pot_tot)*sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n) 
                             
                            ! write(6,*)
                            ! write(6,*) 
                            ! write(6,*) Hmatrix_2(ind1,ind2), ind1,ind2
                            ! write(6,*) i,l,j,m,k,n
                            ! write(6,*) sol(ii), pot1,pot2,pot3  
                            ! write(6,*) Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)
                          end if
                        
                          !
                          if(pot_tot.eq.3) then
                             Hmatrix_3(ind1,ind2) = Hmatrix_3(ind1,ind2) + ((1.0d10)**pot_tot)*sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)                          
                          end if
                          !
                          if(pot_tot.eq.4) then
                             Hmatrix_4(ind1,ind2) = Hmatrix_4(ind1,ind2) + ((1.0d10)**pot_tot)*sol(ii)*Mx1(pot1,i,l)*Mx2(pot2,j,m)*Mx3(pot3,k,n)                          
                          end if
                          


                          !write(13,*) pot1-1, pot2-1, pot3-1, Hupdate

                          !if(alt_state.eq.1) then
                          !   
                          !   if(pot2.ne.1.or.pot3.ne.1)then
                          !      Hupdate=0
                          !   end if

!                          else if(alt_state.eq.2) then
!
!                             if(pot1.ne.1.or.pot3.ne.1)then
!                                Hupdate=0
!                             end if
!!
!
 !                         else if(alt_state.eq.3) then
  !                         
   !                          if(pot1.ne.1.or.pot2.ne.1)then
    !                            
   !                             write(13,*) pot1, pot2, pot3, Hupdate!
!
!                                Hupdate=0
!                             end if
!                             
!                          end if
                          
                          
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

if(0) then

if (alt1.eq.1) then

 Vmatrix = Hmatrix - Pmatrix

  do i=1,nstates3
    write(23,*) "Pmatrix + Hmatrix_2"
    write(23,*) (Pmatrix(i,j) +Hmatrix_2(i,j), j=1, nstates3)
 end do
 do i=1,nstates3
    write(24,*) "Pmatrix "
    write(24,*) (Pmatrix(i,j), j=1, nstates3)
 end do



end if


do i=1,nstates3
    write(22,*) "Hmatrix"
    write(22,*) (Hmatrix(i,j), j=1, nstates3)
 end do

if (alt1.eq.0) then
do i=1,nstates3
    write(23,*) "Hmatrix_harm"
    write(23,*) (Hmatrix_harm(i,j), j=1, nstates3)
 end do
 
end if


 do i=1,nstates3
    write(25,*) "Hmatrix_2"
    write(25,*) (Hmatrix_2(i,j), j=1, nstates3)
  end do


 do i=1,nstates3
    write(26,*) "Hmatrix_3"
    write(26,*) (Hmatrix_3(i,j), j=1, nstates3)
  end do

 do i=1,nstates3
    write(27,*) "Hmatrix_4"
    write(27,*) (Hmatrix_4(i,j), j=1, nstates3)
  end do

end if


write(6,*) "ind1", ind1, "ind2",ind2

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

! mapping from ind to i,j,k
!  do i=1,nstates(1)
!     do j=1,nstates(2)
!        do k=1,nstates(3)
!           
!           ind1 = k + (j-1)*nstates(3) + (i-1)*nstates(2)*nstates(3)
!           
!           ind_to_ijk(ind1,1) = i
!           ind_to_ijk(ind1,2) = j
!           ind_to_ijk(ind1,3) = k
!           
!        end do
!     end do
!  end do
  

!  write(6,*) "Hmatrix after"
!  do i=1,nstates3
!     write(6,*) (Hmatrix(i,j), j=1, nstates3)
!  end do
!if(alt1.eq.1) then
! do i=1,nstates3
!    write(23,*) "Pmatrix"
!     write(23,*) (Pmatrix(i,j), j=1, nstates3)
!  end do
! do i=1,nstates3
!    write(24,*) "Vmatrix"
!     write(24,*) (Vmatrix(i,j), j=1, nstates3)
!  end do
!
!end if

!  write(6,*) "second degree"
!  write(6,*) "Mx1(1,:,:)"
!  write(6,*) Mx1(1,:,:)
!  write(6,*) "Mx1(2,:,:)"
!  write(6,*) Mx1(2,:,:)
!  write(6,*) "Mx1(3,:,:)"
!  write(6,*) Mx1(3,:,:)
!  write(6,*) "Mx1(4,:,:)"
!  write(6,*) Mx1(4,:,:)
!  write(6,*) "Mx1(5,:,:)"
!  write(6,*) Mx1(5,:,:)
!  write(6,*) "Mx1(6,:,:)"
!  write(6,*) Mx1(6,:,:)
!  write(6,*) "Mx1(7,:,:)"
!  write(6,*) Mx1(7,:,:)

!write(6,*) "second degree"
!write(6,*) "Mx1"
!  write(6,*) Mx1(3,:,:)
!  write(6,*) "Mx2"
!  write(6,*) Mx2(3,:,:)
!  write(6,*) "Mx3"
!  write(6,*) Mx3(3,:,:)


!write(6,*) "third degree"
!write(6,*) "Mx1"
!  write(6,*) Mx1(4,:,:)
!  write(6,*) "Mx2"
!  write(6,*) Mx2(4,:,:)
!  write(6,*) "Mx3"
!  write(6,*) Mx3(4,:,:)

!write(6,*) "fourth degree"
!write(6,*) "Mx1"
!  write(6,*) Mx1(5,:,:)
!  write(6,*) "Mx2"
!  write(6,*) Mx2(5,:,:)
!  write(6,*) "Mx3"
!  write(6,*) Mx3(5,:,:)
!write(6,*) 
 ! write(6,*) "Mx1 (1)"
 ! write(6,*) Mx1(1,:,:)
 ! write(6,*) "Mx2 (1)"
 ! write(6,*) Mx2(1,:,:)
 ! write(6,*) "Mx3 (1)"
 ! write(6,*) Mx3(1,:,:)

!write(6,*) "fifth degree"
!write(6,*) "Mx1"
!  write(6,*) Mx1(6,:,:) 
! write(6,*) "Mx2"
!  write(6,*) Mx2(6,:,:)
!  write(6,*) "Mx3"
!  write(6,*) Mx3(6,:,:)
!write(6,*) "sixth degree"
! write(6,*) "Mx1"
!  write(6,*) Mx1(7,:,:)
!  write(6,*) "Mx2 "
!  write(6,*) Mx2(7,:,:)
!  write(6,*) "Mx3 "
!  write(6,*) Mx3(7,:,:)



  
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
 
   write(6,*) "fundamental frequency, dseyv", (W-W(1))*cm*hartree

  ! call singular value decomposition
  !SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
  !allocate( svd_A(nstates3,nstates3), svd_W(nstates3), svd_V(nstates3,nstates3))
  
  !set up svd_A
  !svd_A = Hmatrix + transpose(Hmatrix)  
  !do i=1,nstates3
  !   svd_A(i,i) =  svd_A(i,i) -Hmatrix(i,i)
  !end do


  !CALL SVDCMP(svd_A,nstates3,nstates3,nstates3,nstates3,svd_W,svd_V)
  
  !write(6,*) "svd_W"
  !write(6,*) "*************************"
  !write(6,*) svd_W
  !write(6,*) "*************************"
  
  ! plocka bort sm� v�rden i svd_W

  !svd_tol=1.d-2
  !do i=1,nstates3
  !   if(svd_W(i)< svd_tol) then
  !      write(6,*) "tog bort en med svd_W=", svd_W(i)
  !      svd_W(i) = 0
  !   end if
  !end do
  
  ! svd_A=U, svd_V=V, svd_W=W -> A = U*W*V^T = svd_A * svd_W *

  ! allokera och fixa en matris med diagonalelementen svd_W 
  !allocate(svd_W_mat(nstates3, nstates3), Hmatrix_svd(nstates3,nstates3))
  !svd_W_mat=0
  !do i=1,nstates3
  !   svd_W_mat(i,i) = svd_W(i)
  !end do


   !test
   allocate(Hmatrix_numerical(nstates3,nstates3))
   Hmatrix_numerical = Pmatrix +Hmatrix_2

  ! call lapack routine
   call dsyev( "V", "U", nstates3, Hmatrix_numerical, nstates3, W, WORK, LWORK, INFO)

  ! print eigenvalues
if(0)  then

  write(6,*) "eigenvalues, Pmatrix + Hmatrix2", W
  write(6,*)
 write(6,*)  "fundamental frequency, dseyv", (W-W(1))*cm*hartree

  do i=1,nstates3
     write(11,*) "eigenvector ", i

     do j=1,nstates3
        write(11,*) ind_to_ijk(j,1), ind_to_ijk(j,2), ind_to_ijk(j,3), Hmatrix_lapack(j,i)

     end do

  end do

end if
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
  
  
  
  
end program vib_normal3
