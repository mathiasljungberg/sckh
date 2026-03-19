module m_overlaps
  implicit none
  
contains
  
subroutine ReorderStates(dipole,final_energy,nstates,ntsteps, Overlap_matrix)
   implicit none
   real(kind=wp),intent(inout)::dipole(nstates,ntsteps,3)
   integer, intent(in)::nstates,ntsteps
   real(kind=wp),intent(inout)::final_energy(nstates,ntsteps)
   real(kind=wp),intent(inout)::Overlap_matrix(:,:,:)

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

subroutine ExchangeStates(D_fn,dipole_module,E_f,nstates,ntsteps,i,j,tstep, Overlap_matrix)
   implicit none
   integer,intent(in)::nstates,ntsteps,i,j,tstep
   real(kind=wp),intent(out)::D_fn(nstates,ntsteps,3)
   real(kind=wp),intent(out)::E_f(nstates,ntsteps)
   real(kind=wp),intent(out)::dipole_module(nstates,ntsteps)
   real(kind=wp),intent(inout)::Overlap_matrix(:,:,:)

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


end module m_overlaps

