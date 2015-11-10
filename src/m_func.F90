module m_func
  use m_precision, only: wp
  use m_splines
  implicit none

  real(kind=wp), allocatable, private:: H_in(:,:,:)
  real(kind=wp),allocatable,private::Overlap_matrix(:,:,:)
  real(kind=wp),allocatable,public::time(:)
  integer, private :: H_in_shape(3)

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
   integer::i,j,k,l,m,n,traj


   if( allocated(Overlap_matrix)) then
       write(*,*)'Error in m_func module, routine get_overlap'
       write(*,*)'The overlap matrix has been already allocated'
       stop
   else
       allocate(Overlap_matrix(nstates,nstates,ntsteps))
   endif


   Overlap_matrix(:,:,1)=0.0d0
   do i=1,nstates
     Overlap_matrix(i,i,1)=1.0d0
   enddo


   traj=15
   write(*,*) 'Overlap file reading...', filename
   open(traj,file=filename,form='unformatted',status='old',access='sequential')
   ! Download the overlap elements from the file
   read(traj) n
   if(n .ne. (nstates)) then 
     write(*,*) 'Error in the GetOverlap, n_states from  t12 files and in program are not equal.'
     write(*,*) 'DSTO from t12 files: ',n,' number of states from program: ',nstates
     stop
   endif 
   rewind(traj)


   do k=2,ntsteps
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
   do k=1,ntsteps !nsteps-1
    write(11,*) 'Time step is ',k
    do j=1,nstates
      write(11,'(9F14.8)') (Overlap_matrix(j,i,k),i=1,nstates)
    enddo
   enddo
   close(11)
end subroutine read_overlap



subroutine get_diagonal_hamiltonian(E_f,E_int,nstates,ntsteps,delta_t, E_fn_mean)
   implicit none
   integer,intent(in)::nstates,ntsteps
   real(kind=wp),intent(in),dimension(nstates,ntsteps)::E_f
   real(kind=wp),intent(in),dimension(ntsteps)::E_int
   real(kind=wp),intent(in)::delta_t
   integer::i,j,k,l,m ! loop variables
   real(kind=wp),intent(in):: E_fn_mean


   if(.not. allocated(H_in)) then
       allocate(H_in(nstates,nstates,ntsteps))
   else
       write(*,*) 'Error in m_func module'
       write(*,*) 'H_in has been already allocated'
       stop
   endif
   H_in=0.0_wp
     do i=1,nstates
         H_in(i,i,:)=E_f(i,:)-E_int(:)-E_fn_mean
     enddo ! i

   allocate(time(ntsteps))
   do i=1,ntsteps
     time(i)=(i-1)*delta_t
   enddo

end subroutine get_diagonal_hamiltonian



subroutine funct_complex(t,y,yp)
   implicit none
   real(kind=wp), intent(in):: t
   complex(kind=wp),dimension(:,:),intent(in):: y
   complex(kind=wp),dimension(:,:),intent(out):: yp
   real(kind=wp),allocatable::Ht(:,:)
   integer::shape_H(3)
   integer::i,j,k,l ! loop variables


   if(allocated(H_in)) then
     shape_H=shape(H_in)
     !write(*,*) 'Shape H_in is ',shape_H(1),shape_H(2),shape_H(3)
     allocate(Ht(shape_H(1),shape_H(2)))
     Ht=0.0d0
     do i=1,shape_H(1)
         call spline_one(time,H_in(i,i,:),shape_H(3),t,Ht(i,i)) !spline_one(x,y,n, x2,y2)
     enddo!i
     !yp=(0.0_wp,1.0_wp)*matmul(y,dcmplx(Ht))
     call zgemm('N','N',shape_H(1),shape_H(1),shape_H(1),(1.0d0,0.0d0),y,shape_H(1),&
         dcmplx(Ht),shape_H(1),(0.0d0,0.0d0),yp,shape_H(1))
     yp=(0.0_wp,1.0_wp)*yp
     deallocate(Ht)
   else
      write(*,*) 'Error in m_func module'
      write(*,*) 'H_in is not allocated'
      stop
   endif

end subroutine funct_complex


subroutine clean_memory()

   deallocate(H_in)
   deallocate(time)
   deallocate(Overlap_matrix)

end subroutine clean_memory 


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



subroutine get_hamiltonian_offdiagonal(E_f,E_int,nstates,ntsteps,delta_t,E_f_mean)
   implicit none
   integer,intent(in)::nstates,ntsteps
   real(kind=wp),intent(in),dimension(nstates,ntsteps)::E_f
   real(kind=wp),intent(in),dimension(ntsteps)::E_int
   real(kind=wp),intent(in)::delta_t
   integer::i,j,k,l,m ! loop variables
   real(kind=wp), intent(in)::E_f_mean

   if(.not. allocated(H_in)) then
       allocate(H_in(nstates,nstates,ntsteps))
   else
       write(*,*) 'Error in m_func module'
       write(*,*) 'H_in has been already allocated'
       stop
   endif
   H_in=0.0_wp
   ! for the time=0 there is no transformation
   ! Also we assume that hamiltonian matrix is diagonal at time t=0
     do i=1,nstates
      H_in(i,i,1)=E_f(i,1)
     enddo
   do k=2,ntsteps
     do i=1,nstates
        do j=1,nstates
          do l=1,nstates
             H_in(i,j,k)=H_in(i,j,k)+Overlap_matrix(i,l,k-1)*E_f(nstates+1-l,k)*Overlap_matrix(j,l,k-1)
          enddo ! l
        enddo ! j
     enddo ! i
   enddo ! k
   do i=1,nstates
     H_in(i,i,:)=H_in(i,i,:)-E_int(:)-E_f_mean
   enddo
   open(11,file='Eab_module.txt',status='unknown')
   ! save to file the Eab matrix
   do k=1,ntsteps 
    write(11,*) 'Time step is ',k
    do j=1,nstates
      write(11,'(9F14.8)') (H_in(j,i,k),i=1,nstates)
    enddo!j
   enddo!k
   close(11)
   allocate(time(ntsteps))
   do i=1,ntsteps
     time(i)=(i-1)*delta_t
   enddo
   write(*,*) 'delta_t ',delta_t
end subroutine get_hamiltonian_offdiagonal










end module m_func
