module m_func
  use m_precision, only: wp
  use m_splines
  implicit none

  complex(kind=wp), allocatable, private:: H_in(:,:,:)
  complex(kind=wp),allocatable,private::Ht_diagonal(:,:)
  real(kind=wp),allocatable,private::Overlap_matrix(:,:,:)
  real(kind=wp),allocatable,public::time(:)
  integer, private :: H_in_shape(3)
  real(kind=wp),allocatable,private::yder(:,:,:),H_effective(:,:,:),time_effective(:)
  integer,private:: effective_ntsteps=800
  real(kind=wp),allocatable::yderr_real(:,:,:),yderr_im(:,:,:)
  real(kind=wp),allocatable::yderr_diagonal(:,:)

contains
  
  subroutine func ( t, y, yp )
    !
    !*******************************************************************************
    !
    !! func should evaluate yp_{ji}(t) = sum_k H_{jk}(t) * y_{kj}(t)
    !! that is, we need   H_{jk} and y_{kj} for an arbitrary time
    !! we can call another subroutine to compute H_jk(t) uisng splines
    !! 
    
    implicit none
    !
    real(kind=wp), intent(in):: t
    real(kind=wp), intent(in) :: y(:,:)
    real(kind=wp), intent(out) :: yp(:,:)
    real(kind=wp), dimension( H_in_shape(1) , H_in_shape(2)) :: Ham

    !
   
    if (allocated(H_in)) then
       call get_H(Ham, t)
       yp = matmul(Ham,y)
    else
       write(6,*) "H_in not allocated!"
       stop
    end if

  end subroutine func

  subroutine func_c ( t, y, yp )
    !
    !*******************************************************************************
    !
    !! func should evaluate yp_{ji}(t) = sum_k H_{jk}(t) * y_{kj}(t)
    !! that is, we need   H_{jk} and y_{kj} for an arbitrary time
    !! we can call another subroutine to compute H_jk(t) uisng splines
    !! 
    
    implicit none
    !
    real(kind=wp), intent(in):: t
    complex(kind=wp), intent(in) :: y(:,:)
    complex(kind=wp), intent(out) :: yp(:,:)
    real(kind=wp), dimension( H_in_shape(1) , H_in_shape(2)) :: Ham

    !
    
    if (allocated(H_in)) then
       call get_H(Ham, t)
       yp = dcmplx(0,1) * matmul(Ham,y)
    else
       write(6,*) "H_in not allocated!"
       stop
    end if

  end subroutine func_c

  subroutine get_H(Ham, t)
    real(kind=wp), intent(out):: Ham(:,:)
    real(kind=wp), intent(in):: t

    ! for now, just return a constant H
    Ham = H_in(:,:,1)

  end subroutine get_H
  

  subroutine set_H(H_input)
    real(kind=wp), intent(in):: H_input(:,:,:)
    
    if ( .not. allocated(H_in) ) then
       H_in_shape = shape(H_input) 
       allocate( H_in(H_in_shape(1),H_in_shape(2),H_in_shape(3) ))
    end if
    
    H_in = H_input

  end subroutine set_H
  

function y_exact_1 ( t, k )  
  implicit none
  real(kind=wp):: y_exact_1
  real(kind=wp), intent(in):: t, k

  y_exact_1 = exp(k*t)
  
end function y_exact_1


subroutine read_overlap(nstates,ntsteps,filename)
   implicit none
   character(len=80),intent(in):: filename
   integer,intent(in):: nstates,ntsteps
   integer::i,j,k,l,m,traj,shape_ov(3)
   integer(kind=2)::n
   real(kind=wp),allocatable,dimension(:,:,:)::ov_temp
   if( allocated(Overlap_matrix)) then
       Overlap_matrix=0.0_wp
   else
       allocate(Overlap_matrix(nstates,nstates,ntsteps-1))
       Overlap_matrix=0.0_wp
   endif
   traj=15
   write(*,*) 'Overlap file reading...', filename
   open(traj,file=filename,form='unformatted',status='old',access='sequential')
   ! Download the overlap elements from the file
   rewind(traj)
   read(traj) n
   write(*,*)'The number of states is read overlap routine***********',n
   if(n .ne. (nstates)) then 
     write(*,*) 'Error in the GetOverlap, n_states from  t12 files and in program are not equal.'
     write(*,*) 'DSTO from t12 files: ',n,' number of states from program: ',nstates
     stop
   endif 
   rewind(traj)
   do k=1,ntsteps-1
   read(traj) n
   read(traj) ((Overlap_matrix(i,j,k),i=1,nstates), j=1,nstates)
    if(n .ne. (nstates)) then 
     write(*,*) 'Error in the GetOverlap, n_states from  t12 files and in program are not equal.'
     write(*,*) 'DSTO from t12 files: ',n,' number of states from program: ',nstates
     stop
  endif 
  enddo
  close(traj) ! close overlap file, reading has been done
  open(11,file='overlap_marix_module.txt',status='unknown')
  ! Print the overlap matrix
   write(11,*) 'Print the overlap matrix  '
   do k=1,ntsteps-1 !nsteps-1
    write(11,*) 'Time step is ',k,nstates
    do j=1,nstates
      write(11,'(9F14.8)') (Overlap_matrix(j,i,k),i=1,nstates)
    enddo
   enddo
  close(11)
  ! Arrange overlap matrix in appropriate way 
   shape_ov=shape(Overlap_matrix)
   allocate(ov_temp(shape_ov(1),shape_ov(2),shape_ov(3)))
   ov_temp=Overlap_matrix
   do i=1,shape_ov(1)-1
     do j=1,shape_ov(2)-1
      Overlap_matrix(j,i,:)=ov_temp(j+1,i+1,:)
     enddo
   enddo
   Overlap_matrix(shape_ov(1),shape_ov(1),:)=ov_temp(1,1,:)
   do i=1,shape_ov(1)-1
     Overlap_matrix(shape_ov(1),i,:)=ov_temp(1,i+1,:)
     Overlap_matrix(i,shape_ov(1),:)=ov_temp(i+1,1,:)
   enddo
   open(11,file='overlap_marix_after_reordering.txt',status='unknown')
 !   Print the overlap matrix after reordering
   write(11,*) 'Print the overlap matrix read from t12 file  ', shape(Overlap_matrix)
   do k=1,shape_ov(3) !nsteps-1
    write(11,*) 'Time step is ',k
    do j=1,shape_ov(1)
      write(11,'(9F14.8)') (Overlap_matrix(j,i,k),i=1,nstates)
    enddo
   enddo
   close(11)
   deallocate(ov_temp)
end subroutine read_overlap



subroutine get_hamiltonian(E_f,E_int,nstates,ntsteps,delta_t,E_f_mean)
   implicit none
   integer,intent(in)::nstates,ntsteps
   real(kind=wp),intent(in),dimension(nstates,ntsteps)::E_f
   real(kind=wp),intent(in),dimension(ntsteps)::E_int
   real(kind=wp),intent(in)::delta_t
   integer::i,j,k,l,m ! loop variables
   real(kind=wp), intent(in)::E_f_mean
   real(kind=wp),allocatable::Unitary_matrix(:,:)
   real(kind=wp):: dummy
   real(kind=wp),allocatable::dummy_array(:)

   if(.not. allocated(H_in)) then
       allocate(H_in(nstates,nstates,ntsteps))
   else
       write(*,*) 'Error in m_func module'
       write(*,*) 'H_in has been already allocated'
       stop
   endif
   ! define unitary matrix
   allocate(Unitary_matrix(nstates,nstates))
   Unitary_matrix=0.0_wp
   do i=1,nstates
      Unitary_matrix(i,i)=1.0_wp
   enddo
   H_in=0.0_wp
   
   ! Put another approximation formula for the first time step
    H_in(:,:,1)=(0.0_wp,1.0_wp)/delta_t*(Overlap_matrix(:,:,1)-Unitary_matrix)

   do k=2,ntsteps
     do i=1,nstates
        do j=1,nstates
             H_in(j,i,k)=(0.0d0,1.0d0)/2*delta_t*( Overlap_matrix(j,i,k)+Overlap_matrix(i,j,k-1) )
        enddo ! j
     enddo ! i
   enddo ! k
   ! Use another formula for the last time step
   H_in(:,:,ntsteps)=(0.0_wp,1.0_wp)/delta_t*(Unitary_matrix-transpose(Overlap_matrix(:,:,ntsteps-1)))
   ! take care about diagonal elements

   do i=1,nstates
     H_in(i,i,:)=(E_f(i,:)-E_int(:)-E_f_mean)
   enddo
   allocate(dummy_array(nstates))
   open(11,file='Eab_module.txt',status='unknown')
   ! save to file the Eab matrix
   do k=1,ntsteps 
    write(11,*) 'Time step is ',k
    do j=1,nstates
      dummy_array(:)=aimag(H_in(j,:,k))+real(H_in(j,:,k))
      write(11,'(9F14.8)') (dummy_array(i),i=1,nstates)
    enddo!j
   enddo!k
   close(11)
   allocate(time(ntsteps))
   do i=1,ntsteps
     time(i)=(i-1)*delta_t
   enddo
   write(*,*) 'delta_t ',delta_t
   deallocate(Unitary_matrix)
end subroutine get_hamiltonian




subroutine get_H_neighboring(E_f,E_int,nstates,ntsteps,delta_t, E_fn_mean)
   implicit none
   integer,intent(in)::nstates,ntsteps
   real(kind=wp),intent(in),dimension(nstates,ntsteps)::E_f
   real(kind=wp),intent(in),dimension(ntsteps)::E_int
   real(kind=wp),intent(in)::delta_t
   integer::i,j,k,l,m ! loop variables
   real(kind=wp),allocatable::H_half(:,:,:),time_half(:),time_n(:),&
    Unitary_matrix(:,:),H_in_temp(:,:,:)
   real(kind=wp),intent(in):: E_fn_mean
   real(kind=wp),allocatable::dummy_array(:)
   real(kind=wp)::dummy


   if(.not. allocated(H_in)) then
       allocate(H_in(nstates,nstates,ntsteps))
       H_in=0.0_wp
   else
       H_in=0.0_wp
   endif
   allocate(H_half(nstates,nstates,ntsteps-1))
   allocate(time_half(ntsteps-1))
   allocate(time_n(ntsteps-2))
   do i=1,ntsteps-2
     time_n(i)=(i)*delta_t
   enddo
! time_n contains all time points except the first and last one
   do i=1,ntsteps-1
     time_half(i)=(i-0.5)*delta_t
   enddo
! time_half contains the time coordinates of the middle time interval for which the 
! derivatives will be computed   
   do i=1,ntsteps-1
    do j=1,nstates
      do k=1,nstates
      ! compute matrix elements according to the formula 27 in Hammers-Schiffer and J.C. Tully: Proton transfer in solution
      H_half(j,k,i)=1/(2*delta_t)*(Overlap_matrix(j,k,i)-Overlap_matrix(k,j,i))
      ! This matrix should be arranged to the middle of the time interval
      enddo
    enddo     
   enddo
   open(11,file='Eab_module_before_symm.txt',status='unknown')
   do k=1,ntsteps-1 
    write(11,*) 'Time step is ',k
    do j=1,nstates
        write(11,'(9F14.8)') (H_half(j,i,k),i=1,nstates)
    enddo!j
   enddo!k

   ! Now the each point in H_half matrix corresponds to the time moments t_i=(i+0.5)*delta_t
   allocate(H_in_temp(nstates,nstates,ntsteps))
   H_in_temp=0.0_wp
   ! Now we may spline the Hamiltonian matrix elements defined at in the middle of time interval to the ends of each intervals
   write(*,*) 'Shape ',shape(time_n),shape(time_half),shape(H_half(1,1,:)),shape(H_in_temp(j,i,2:ntsteps-1))
   do i=1,nstates
     do j=1,nstates
       call spline_easy(time_half,H_half(j,i,:),ntsteps-1,time_n,H_in_temp(j,i,2:ntsteps-1),ntsteps-2)
     enddo
   enddo
   deallocate(time_n,time_half)
   ! Explicitly compute the H elements for the last and first time steps 
   allocate(Unitary_matrix(nstates,nstates))
   Unitary_matrix=0.0_wp
   do i=1,nstates
      Unitary_matrix(i,i)=1.0_wp
   enddo

   ! Put another approximation formula for the first time step
    H_in(:,:,1)=(0.0_wp,1.0_wp)/delta_t*(Overlap_matrix(:,:,1)-Unitary_matrix)
    H_in(:,:,ntsteps)=(0.0_wp,1.0_wp)/delta_t*(Unitary_matrix-transpose(Overlap_matrix(:,:,ntsteps-1)))
   write(11,*) 'Raw data for the last and first point of the H_in array'
   write(11,*) H_in(:,:,1), '****'
   write(11,*) H_in(:,:,1) ,'****'
   ! copy previously splined matrix
   H_in(:,:,2:ntsteps-1)=(0.0_wp,1.0_wp)*H_in_temp(:,:,2:ntsteps-1)
   !  store matrix before symmetrization

   do i=1,nstates
     H_in(i,i,:)=(E_f(i,:)-E_int(:)-E_fn_mean)
   enddo
    allocate(dummy_array(nstates))
   ! save to file the Eab matrix
   do k=1,ntsteps 
    write(11,*) 'Time step is ',k
    do j=1,nstates
      dummy_array(:)=aimag(H_in(j,:,k))+real(H_in(j,:,k))
      write(11,'(9F14.8)') (dummy_array(i),i=1,nstates)
    enddo!j
   enddo!k

   close(11)
   ! symmetrize the offdiagonal elements because otherwise 
   ! matrix coefficients become too large
   do k=1,ntsteps
     do i=1,nstates
       do j=1,nstates
         if(i.ne.j) then

           dummy=( abs(H_in(i,j,k))+abs(H_in(j,i,k)) )/2.0_wp
           if(aimag(H_in(i,j,k)).gt.0) then
             H_in(i,j,k)=(0.0_wp,1.0_wp)*dummy
           else
             H_in(i,j,k)=-(0.0_wp,1.0_wp)*dummy
           endif

           if(aimag(H_in(j,i,k)).gt.0) then
             H_in(j,i,k)=(0.0_wp,1.0_wp)*dummy
           else
             H_in(j,i,k)=-(0.0_wp,1.0_wp)*dummy
           endif            
          endif
       enddo
     enddo
   enddo
   ! For the test cases - modify the Hab elements 3 and 4
   !H_in(8,4,:)=10000*H_in(8,7,:)
   !H_in(4,8,:)=10000*H_in(7,8,:)

   ! write the H matrix to the file
   open(11,file='Eab_module.txt',status='unknown')
   ! save to file the Eab matrix
   do k=1,ntsteps 
    write(11,*) 'Time step is ',k
    do j=1,nstates
      dummy_array(:)=aimag(H_in(j,:,k))+real(H_in(j,:,k))
      write(11,'(9F14.8)') (dummy_array(i),i=1,nstates)
    enddo!j
   enddo!k
   close(11)
   deallocate(Unitary_matrix,dummy_array,H_in_temp,H_half)
   if (.not.allocated(time)) then
     allocate(time(ntsteps))
     do i=1,ntsteps
       time(i)=(i-1)*delta_t
     enddo
   endif

end subroutine get_H_neighboring


subroutine funct_complex(t,y,yp) 
   implicit none
   real(kind=wp), intent(in):: t
   complex(kind=wp),dimension(:,:),intent(in):: y
   complex(kind=wp),dimension(:,:),intent(out):: yp
   complex(kind=wp),allocatable::Ht(:,:)
   real(kind=wp), allocatable::Ht_real(:,:),Ht_im(:,:)
   integer::shape_H(3)
   integer::i,j,k,l ! loop variables

   !write(6,*) 'Inside routine funct_complex'
   if(allocated(H_in)) then
     shape_H=shape(H_in)
     !end if
     ! write(*,*) 'Shape H_in is ',shape_H(1),shape_H(2),shape_H(3)
     if(.not. allocated(Ht)) then
       allocate(Ht(shape_H(1),shape_H(2)))
       allocate(Ht_real(shape_H(1),shape_H(2)),Ht_im(shape_H(1),shape_H(2)))
     endif
     
     if(allocated(yderr_real) .and. allocated(yderr_im)) then
       !write(6,*) 'Before splint array'
       call splint_array3_d(time,real(H_in),yderr_real,t,Ht_real)
       call splint_array3_d(time,aimag(H_in),yderr_im,t,Ht_im)
       !write(6,*) 'After splint array'
       Ht=cmplx(Ht_real,Ht_im)
     else
       !write(*,*) 'Before allocation y_derr array ',shape_H
       allocate( yderr_real(shape_H(1),shape_H(2),shape_H(3)) )
       allocate(yderr_im(shape_H(1),shape_H(2),shape_H(3)))
       !write(*,*) 'After allocation y_derr array '
       do i=1,shape_H(1)
         do j=1,shape_H(2)
           !write(*,*) 'spline real '
           call spline(time,real(H_in(i,j,:)),shape_H(3),1.0d30,1.0d30,yderr_real(i,j,:))
           !write(*,*) 'spline imagine '
           call spline(time,aimag(H_in(i,j,:)),shape_H(3),1.0d30,1.0d30,yderr_im(i,j,:))
         enddo
       enddo
       !write(*,*)'time t',t
       !write(*,*) 'time ',shape(time),'yderr_real ',shape(yderr_real), 'Ht_real ',shape(Ht_real)
       ! write(*,*) ' splint_array3_d real before calling  ','H_in ',shape(H_in)
       !write(*,*) 'Shape time', shape(time)
       !write(*,*) 'Shape yderr_real',shape(yderr_real)
       !write(*,*)'Shape H_in', shape(H_in)
       !write(*,*) 'Shape Ht_real',shape(Ht_real)
       ! Test allocation for different arrays which are calling from the splint_array3_d routine
       !write(*,*) 'time',allocated(time)
       !write(*,*) 'yderr_real',allocated(yderr_real)
       !write(*,*) 'H_in',allocated(H_in)
       !write(*,*) 'Ht_real',allocated(Ht_real)
       !write(*,*) 't ',t
       call splint_array3_d(time,real(H_in,8),yderr_real,t,Ht_real)
       !write(*,*) ' splint_array3_d imagine ', shape(H_in)
       
       call splint_array3_d(time,aimag(H_in),yderr_im,t,Ht_im)
       Ht=cmplx(Ht_real,Ht_im)
     endif
     !write(*,*) 'Before library zgemm routine'
     !call zgemm('N','N',shape_H(1),shape_H(1),shape_H(1),(0.0d0,1.0d0),y,shape_H(1),Ht,shape_H(1),(0.0d0,0.0d0),yp,shape_H(1))
     yp=(0.0_wp,1.0_wp)*yp
     yp=(0.0_wp,1.0_wp)*matmul(y,Ht)
   else
     write(*,*) 'Error in m_func module'
     write(*,*) 'H_in is not allocated'
     stop
   endif
   !write(*,*) 'end routine funct_complex'
 end subroutine funct_complex





subroutine transform_dipole_operator(D_nf,E_final,E_n,nstates,ntimesteps)
   implicit none
   integer,intent(in)::nstates,ntimesteps
   real(kind=wp),dimension(nstates,ntimesteps),intent(in)::E_final
   real(kind=wp),dimension(ntimesteps),intent(in)::E_n
   real(kind=wp),dimension(nstates,ntimesteps,3),intent(out)::D_nf
   real(kind=wp),dimension(:,:,:),allocatable::D_nf_trans
   integer::i,j,k,l,m ! loop variables
   allocate(D_nf_trans(nstates,ntimesteps,3))
   D_nf_trans=0.0d0
   ! Build D_ab^{~}=D_ab(t)*(E_b(t)-E_a(t)) - instanteneous basis
   do m=1,3
       do i=1,nstates
            D_nf_trans(i,:,m)=D_nf(i,:,m)*(E_final(i,:)-E_n(:))
       enddo
   enddo 
   D_nf=0.0d0  
   D_nf(:,:,1)=D_nf_trans(:,:,1) ! keep the dipole operator for step k=1 unchanged
   do m=1,3 ! dipole polarization
     do k=2,ntimesteps
       do i=1,nstates
         do  j=1,nstates
          D_nf(i,k,m)=D_nf(i,k,m)+Overlap_matrix(nstates+1,nstates+1,k-1)*D_nf_trans(j,k,m)*Overlap_matrix(i,nstates+1-j,k-1)
         enddo!j
       enddo!i
     enddo!k
   enddo !m
   deallocate(D_nf_trans)

end subroutine transform_dipole_operator


subroutine ReorderStates(dipole,final_energy,nstates,ntsteps)
   implicit none
   integer::nstates,ntsteps
   real(kind=wp),intent(inout)::dipole(nstates,ntsteps,3)
   real(kind=wp),intent(inout)::final_energy(nstates,ntsteps)
   real(kind=wp),allocatable::dipole_module(:,:)
   real(kind=wp)::dipole_module_tol, energy_tol, ov_tol
   real(kind=wp)::dip_module_small,dip_module_large 
   real(kind=wp)::cut_off_small,cut_off_strict,max_offdiagonal
   integer::shape_ov(3)
   integer::i,j,k,l,m
   integer,allocatable::control_state(:)
   real(kind=wp)::dummy
   integer::a_i,b_i



   dipole_module_tol=0.0025_wp
   energy_tol=0.01_wp
   ov_tol=0.7_wp
   if(allocated(Overlap_matrix)) then
   shape_ov=shape(Overlap_matrix)
   !write(*,*) 'Shape overlap matrix ',shape_ov
   !write(*,*) 'ReorderState routine nstates and ntsteps ',nstates,ntsteps
   allocate(dipole_module(nstates,ntsteps))
   ! compute the dipole module
   dipole_module=0.0_wp
   do i=1,nstates
     do j=1,3
       dipole_module(i,:)=dipole_module(i,:)+dipole(i,:,j)*dipole(i,:,j)
     enddo
   enddo
   dipole_module=dipole_module*1000
   !write(*,*) dipole_module(3,:)
   ! start to check the trajectory data
   allocate(control_state(nstates))
   !do k=1,ntsteps   
   !control_state=0
    !  do j=1,nstates
     !  dummy=abs(dipole_module(j,k+1)-dipole_module(j,k))
      ! if(dummy > dipole_module_tol)then
       !  control_state(j)=1
         !write(*,*) 'Step is ',k,' state is ',j
       !endif
      !enddo
      !do j=1,nstates
       !  if(control_state(j).ne. 0)then 
        !   do l=j+1,nstates
         !    if(control_state(l).ne. 0) then
          !     dummy=abs(final_energy(j,k)-final_energy(l,k))
           !    if(dummy < energy_tol)then
            !     if((abs(Overlap_matrix(j,l,k))> ov_tol) .and.(abs(Overlap_matrix(l,j,k))> ov_tol)) then
                  ! change the state ordering 
                  !write(*,*)'timestep is ',k,'states are ',l,j
                  !write(*,*)'Overlap off-diagonal',Overlap_matrix(j,l,k),' ',Overlap_matrix(l,j,k)
                  !write(*,*)'Module of dipole ',dipole_module(j,k),' ',dipole_module(j,k+1)
                  !write(*,*)'Module of dipole ',dipole_module(l,k),' ',dipole_module(l,k+1)
                  !write(*,*)'Final energy ',final_energy(j,k),' ',final_energy(j,k+1)
                  !write(*,*)'Final energy ',final_energy(l,k),' ',final_energy(l,k+1)
             !     call ExchangeStates(dipole,dipole_module,final_energy,nstates,ntsteps,j,l,k+1)
              !    control_state(j)=0
               !  endif
              ! endif
            ! endif
           !enddo!l
         !endif
      !enddo!j
   !enddo
   ! Try to write some code which should test the ofdiagonal elements
   ! If some zero offdiagonal elements are found this means that such states should be reordered
   !write(*,*)'Find small diagonal elements *****************************************************'
   ! new variables
   ov_tol=0.9
   max_offdiagonal=0
   cut_off_strict=0.9
   cut_off_small=0.4
   dip_module_small=0.5
   dip_module_large=2
   do k=1,ntsteps-1
     do i=1,shape_ov(1)
   if(abs(Overlap_matrix(i,i,k))< ov_tol) then ! large if
        a_i=i
	b_i=0
        max_offdiagonal=0
        ! try to find the maximum offdiaginal element
        do l=1,shape_ov(1)
          !write(*,*)'Overlap_matrix(a_i,l,k) ',abs(Overlap_matrix(a_i,l,k)),' ',a_i,' ',l
           dummy=abs( Overlap_matrix(a_i,l,k))
          if( ( dummy .gt. max_offdiagonal ) .and. (l .ne. a_i) ) then
            max_offdiagonal=dummy
            b_i=l
          endif
        enddo!l
       if ( max_offdiagonal .gt. abs(Overlap_matrix(a_i,a_i,k)) ) then
         if(   (max_offdiagonal .gt. cut_off_strict )  ) then !.and. ( abs( Overlap_matrix(b_i,a_i,k) ) .gt. cut_off_strict)  ) then
           if (    abs( Overlap_matrix(a_i,a_i,k) ) .lt. cut_off_small ) then
             dummy=abs(final_energy(a_i,k)/final_energy(b_i,k))
            if ( abs(1-dummy) .lt. energy_tol) then
           write(*,*)'Found small diagonal element large energy criterion *************************************'
           write(*,*) 'Timestep ',k,' state is',i
           write(*,*) 'Overlap_matrix diagonal value', Overlap_matrix(i,i,k)
           write(*,*) 'b_i ',b_i,'a_i ',a_i
           write(*,*)'Final state energy  a_i',final_energy(a_i,k),'final state energy b_i',final_energy(b_i,k)
           write(*,*)'Max off-diagonal element, position  ',max_offdiagonal,' ',b_i
          ! write(*,'(A20,F15.9,F15.9)') 'Offdiagonal overlap',Overlap_matrix(a_i,b_i,k),Overlap_matrix(a_i,b_i,k)
            dummy=abs(final_energy(a_i,k)-final_energy(b_i,k))
            write(*,*) 'Diagonal energy ',final_energy(a_i,k),'Energy offdiagonal state ',final_energy(b_i,k),&
              '  ',final_energy(a_i,k)/final_energy(b_i,k)
            write(*,*)'Difference in energy ',dummy,energy_tol
            write(*,*)'Diagonal dipole',dipole_module(a_i,k),' ',dipole_module(a_i,k+1)
            write(*,*)'Offdiagonal dipole ', dipole_module(b_i,k),'  ',dipole_module(b_i,k+1)
            write(*,*) 'Dipole dividing ', dipole_module(a_i,k)/dipole_module(a_i,k+1),dipole_module(b_i,k)/dipole_module(b_i,k+1)
        !   if(dummy<energy_tol) then
            call ExchangeStates(dipole,dipole_module,final_energy,nstates,ntsteps,a_i,b_i,k+1)
            write(*,*)'states ',a_i,' and ',b_i,' were exchanged due large difference between diagonal and offdiagonal elements'
           endif
          endif
        else 
          if ( (dipole_module(a_i,k)/dipole_module(a_i,k+1) .lt. dip_module_small ) .or. &
             (dipole_module(a_i,k)/dipole_module(a_i,k+1) .gt. dip_module_large) ) then
              dummy=abs(final_energy(a_i,k)/final_energy(b_i,k))
             if ( abs(1-dummy) .lt. energy_tol ) then
               dummy=abs(final_energy(a_i,k)-final_energy(b_i,k))
                write(*,*)'Found small diagonal element dipole criterion--------------------------------------'
               write(*,*) 'Timestep ',k,' state is',i
               write(*,*) 'Overlap_matrix diagonal value', Overlap_matrix(i,i,k)
               write(*,*) 'b_i ',b_i,'a_i ',a_i
               write(*,*)'Final state energy  a_i',final_energy(a_i,k),'final state energy b_i',final_energy(b_i,k)
                write(*,*)'Max off-diagonal element, position  ',max_offdiagonal,' ',b_i
               write(*,*) 'Diagonal energy ',final_energy(a_i,k),'Energy offdiagonal state ',final_energy(b_i,k),&
                  '  ',final_energy(a_i,k)/final_energy(b_i,k)
               write(*,*)'Difference in energy ',dummy,energy_tol
               write(*,*)'Diagonal dipole',dipole_module(a_i,k),'Offdiagonal dipole ', dipole_module(b_i,k)
               write(*,*) 'Dipole dividing ', dipole_module(a_i,k)/dipole_module(b_i,k)
        !   if(dummy<energy_tol) then
               call ExchangeStates(dipole,dipole_module,final_energy,nstates,ntsteps,a_i,b_i,k+1)
               write(*,*)'states ',a_i,' and ',b_i,' were exchanged using continious dipole module criterion'
            endif
          endif
        endif 
       
       ! Now if maximum off_diagonal element large than diagonal and there is exist jump in dipole moment module - such states should be exchanged
       endif
      endif ! end large if
     enddo!i
   enddo! k
   
   open(11,file='overlap_marix_after_sorting.txt',status='unknown')
  ! Print the overlap matrix
   write(11,*) 'Print the overlap matrix  '
   do k=1,ntsteps-1 !nsteps-1
    write(11,*) 'Time step is ',k
    do j=1,nstates+1
      write(11,'(9F14.8)') (Overlap_matrix(j,i,k),i=1,nstates+1)
    enddo
   enddo
  close(11)
  open(11,file='dipole_module_after_sorting.txt',status='unknown')
  ! Print the overlap matrix
   write(11,*) 'Print the dipole moment  '
  do l=1,nstates
   write(11,*) 'state ',l
   do k=1,ntsteps !nsteps-1
    write(11,'(I5,E15.8)') k,dipole_module(l,k)
    enddo
   enddo
  close(11)
  open(11,file='final_energy_after_sorting.txt',status='unknown')
  ! Print the overlap matrix
   write(11,*) 'Print the final energies '
  do l=1,nstates
   write(11,*) 'state ',l
   do k=1,ntsteps !nsteps-1
    write(11,'(I5,F18.9)') k,final_energy(l,k)
    enddo
   enddo
  deallocate(control_state,dipole_module)
  close(11)
   else
     write(*,*) 'The Overlap matrix is not allocated. Program stop'
     stop
   endif

end subroutine ReorderStates

subroutine ExchangeStates(D_fn,dipole_module,E_f,nstates,ntsteps,i,j,tstep)
   implicit none
   integer,intent(in)::nstates,ntsteps,i,j,tstep
   real(kind=wp),intent(out)::D_fn(nstates,ntsteps,3)
   real(kind=wp),intent(out)::E_f(nstates,ntsteps)
   real(kind=wp),intent(out)::dipole_module(nstates,ntsteps)
   integer::k,l,m,n ! loop variables
   real,allocatable::dummy1d(:)
   integer::shape1d(1),shape3d(3)
   
   ! exchange final state energy
   shape1d=shape(E_f(1,tstep:ntsteps))
   allocate(dummy1d(shape1d(1)))
   dummy1d(:)=E_f(i,tstep:ntsteps)
   E_f(i,tstep:ntsteps)=E_f(j,tstep:ntsteps)
   E_f(j,tstep:ntsteps)=dummy1d(:)
   ! exchange dipole module elements
   dummy1d(:)=dipole_module(i,tstep:ntsteps)
   dipole_module(i,tstep:ntsteps)=dipole_module(j,tstep:ntsteps)
   dipole_module(j,tstep:ntsteps)=dummy1d(:)
   ! exchange dipole components
   do m=1,3 ! dipole polarization
     dummy1d(:)=D_fn(i,tstep:ntsteps,m)
     D_fn(i,tstep:ntsteps,m)=D_fn(j,tstep:ntsteps,m)
     D_fn(j,tstep:ntsteps,m)=dummy1d(:)
   enddo !m
   ! exchange the overlap matrix
   shape3d=shape(Overlap_matrix)
   deallocate(dummy1d)
   !write(*,*) 'Shape overlap matrix ',shape(Overlap_matrix)
   allocate(dummy1d(nstates+1))
   !write(*,*)'Input variables ',nstates,ntsteps
   ! exchange the column for the timestep with large  off diagonal elements
   dummy1d(:)=Overlap_matrix(:,i,tstep-1)
   Overlap_matrix(:,i,tstep-1)=Overlap_matrix(:,j,tstep-1)
   Overlap_matrix(:,j,tstep-1)=dummy1d(:)
   do k=tstep, shape3d(3)
     ! exchange the columns
     dummy1d(:)=Overlap_matrix(:,i,k)
     Overlap_matrix(:,i,k)=Overlap_matrix(:,j,k)
     Overlap_matrix(:,j,k)=dummy1d(:)
     ! exchange the lines
     dummy1d(:)=Overlap_matrix(i,:,k)
     Overlap_matrix(i,:,k)=Overlap_matrix(j,:,k)
     Overlap_matrix(j,:,k)=dummy1d(:)
   enddo
   deallocate(dummy1d)
    

end subroutine ExchangeStates


subroutine get_H_diagonal(E_f,E_int,nstates,ntsteps,delta_t, E_fn_mean)
   implicit none
   integer,intent(in)::nstates,ntsteps
   real(kind=wp),intent(in),dimension(nstates,ntsteps)::E_f
   real(kind=wp),intent(in),dimension(ntsteps)::E_int
   real(kind=wp),intent(in)::delta_t
   integer::i,j,k,l,m ! loop variables
   real(kind=wp),intent(in):: E_fn_mean
   real(kind=wp),allocatable::dummy_array(:)
   real(kind=wp)::dummy


   if(.not. allocated(H_in)) then
       allocate(H_in(nstates,nstates,ntsteps))
       H_in=0.0_wp
   else
       H_in=0.0_wp
   endif
   do i=1,nstates
     H_in(i,i,:)=(E_f(i,:)-E_int(:)-E_fn_mean)
   enddo
   allocate(dummy_array(nstates))
   open(11,file='Eab_module.txt',status='unknown')
   ! save to file the Eab matrix
   do k=1,ntsteps 
    write(11,*) 'Time step is ',k
    do j=1,nstates
      write(11,'(F14.8)') aimag(H_in(j,j,k))+real(H_in(j,j,k))
    enddo!j
   enddo!k
   close(11)
   if (.not.allocated(time)) then
     allocate(time(ntsteps))
     do i=1,ntsteps
       time(i)=(i-1)*delta_t
     enddo
   endif
   deallocate(dummy_array)

end subroutine get_H_diagonal


subroutine funct_complex_diagonal(t,y,yp) 
   implicit none
   real(kind=wp), intent(in):: t
   complex(kind=wp),dimension(:),intent(in):: y
   complex(kind=wp),dimension(:),intent(out):: yp
   real(kind=wp), allocatable::Ht_real(:)
   integer::shape_H(3)
   integer::i,j,k,l ! loop variables
   !write(*,*) 'Inside routine funct_complex'
   if(allocated(H_in)) then
     shape_H=shape(H_in)
     !write(*,*) 'Shape H_in is ',shape_H(1),shape_H(2),shape_H(3)
   if(.not. allocated(Ht_real)) then
        allocate( Ht_real( shape_H(1) ) )
   endif
   Ht_real=0.0_wp
   if (.not. allocated(Ht_diagonal)) then
        allocate(Ht_diagonal(shape_H(1),shape_H(3)))
        Ht_diagonal=0.0_wp
       do i=1,shape_H(1)
         Ht_diagonal(i,:)=H_in(i,i,:)
       enddo
   end if

    if( allocated(yderr_diagonal) ) then
      !write(*,*) 'Before splint array'
      call splint_array2_d(time,real(Ht_diagonal),yderr_diagonal,t,Ht_real)
     else
       !write(*,*) 'Before allocation y_derr_diagonal array ',shape_H
       allocate( yderr_diagonal(shape_H(1),shape_H(3)) ) 
      ! write(*,*) 'After allocation y_derr_diagonal array '
       do i=1,shape_H(1)
	  !write(*,*) shape(time)
          !write(*,*) shape(Ht_diagonal(i,:))
          !write(*,*) shape(yderr_diagonal(i,:)),shape_H(3)
          call spline(time,real(H_in(i,i,:)),shape_H(3),1.0d30,1.0d30,yderr_diagonal(i,:))
          !write(*,*) 'after spline '
       enddo
      call splint_array2_d(time,real(Ht_diagonal),yderr_diagonal,t,Ht_real)
      !Ht_diagonal=cmplx(Ht_diagonal)
     endif
     !write(*,*) 'Before library zgemm routine'
     !call zgemm('N','N',shape_H(1),shape_H(1),shape_H(1),(0.0d0,1.0d0),y,shape_H(1),Ht,shape_H(1),(0.0d0,0.0d0),yp,shape_H(1))
     yp=(0.0_wp,1.0_wp)*y*dcmplx(Ht_real)
   else
      write(*,*) 'Error in m_func module'
      write(*,*) 'H_in is not allocated'
      stop
   endif
  !write(*,*) 'end routine funct_complex'
end subroutine funct_complex_diagonal

subroutine get_diagonal_hamiltonian(E_f,E_int,nstates,ntsteps,delta_t, E_fn_mean)
   implicit none
   integer,intent(in)::nstates,ntsteps
   real(kind=wp),intent(in),dimension(nstates,ntsteps)::E_f
   real(kind=wp),intent(in),dimension(ntsteps)::E_int
   real(kind=wp),intent(in)::delta_t
   integer::i,j,k,l,m ! loop variables
   real(kind=wp),intent(in):: E_fn_mean


   if(allocated(H_in)) deallocate(H_in)
   allocate(H_in(nstates,nstates,ntsteps))
     
   H_in=0.0_wp
   do i=1,nstates
     H_in(i,i,:)=E_f(i,:)-E_int(:)-E_fn_mean
   enddo ! i
   
   ! first time, allocate module variable time(:)
   if(.not. allocated(time)) then
     allocate(time(ntsteps))
     do i=1,ntsteps
       time(i)=(i-1)*delta_t
     enddo
   end if
   
end subroutine get_diagonal_hamiltonian


subroutine clean_variable()

  if (allocated(yderr_diagonal)) then
    deallocate(yderr_diagonal)
  endif
  if (allocated(H_in)) then
    deallocate(H_in)
  endif
  if (allocated(yderr_real)) then
    deallocate(yderr_real)
  endif
  if (allocated(yderr_im)) then
    deallocate(yderr_im)
  endif
  if (allocated(Overlap_matrix)) then
    deallocate(Overlap_matrix)
  endif
  if (allocated(Ht_diagonal)) then
    deallocate(Ht_diagonal)
  endif
  if (allocated(time)) then
    deallocate(time)
  endif
end subroutine clean_variable





end module m_func
