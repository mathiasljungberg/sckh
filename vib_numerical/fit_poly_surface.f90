subroutine fit_poly_surface(npoints, pstart, pdeg, ncols, q, rhs, map, sol)
  use parameters
  implicit none

  ! passed variables
  integer, intent(in) :: npoints, pdeg, pstart, ncols
  real(kind=wp), intent(in), dimension(3,npoints) :: q
  real(kind=wp), intent(in), dimension(npoints) :: rhs
  integer, intent(out), dimension(3,ncols) :: map
  real(kind=wp), intent(out), dimension(ncols) :: sol

  ! local variables
  real(kind=wp), dimension(:,:),allocatable::matrix
!dimension(npoints,ncols) :: matrix
  real(kind=wp), dimension(npoints):: sol_tmp
  integer:: LWORK,INFO
  real(kind=wp), dimension(:),allocatable:: WORK

  ! loop vaiables
  integer:: point,count,deg,i,j,k 

  ! sätter upp matrisen med coefficienter, q [Å]
  allocate(matrix(npoints,ncols))

   do point=1,npoints
     count=0
     
     do deg=pstart,pdeg
        ! loopar igenom alla n,m,l med summa deg
        do i=0,deg
           do j=0,deg-i

              k=deg-(i+j)
              count = count +1
              
              if(point.eq.1) then
                 map(1,count) = i
                 map(2,count) = j
                 map(3,count) = k                 
              end if


              matrix(point,count)= (q(1,point)**i)*(q(2,point)**j)*(q(3,point)**k)

           end do
        end do
     end do
  end do


  ! solve the  lsq problem 
  !*************************************
  sol_tmp =rhs
  LWORK = 2*npoints*ncols +100
  allocate( WORK(LWORK))

  ! DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
  call  DGELS( 'N', npoints, ncols, 1, matrix, npoints, sol_tmp, npoints, WORK, LWORK, INFO )

  deallocate( WORK )
  deallocate( matrix )


  !write(12,*) INFO
  !write(12,*) sol(1:ncols)
  !write(12,*)

 sol = sol_tmp(1:ncols)




end subroutine fit_poly_surface
