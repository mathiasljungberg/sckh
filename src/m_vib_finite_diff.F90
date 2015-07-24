module m_vib_finite_diff
  implicit none
contains
  
  subroutine solve_finite_diff_spline(x_in, y_in, npoints, nstates, mu, eigval, eigvec, units_in)
    use m_precision, only: wp
    use m_constants, only: const
    use m_splines, only: linspace, spline_easy
    implicit none
    
    integer,intent(in)::npoints,nstates
    real(kind=wp), intent(in):: x_in(:),y_in(:)
    real(kind=wp),intent(in)::  mu
    real(wp), intent(out), allocatable:: eigval(:), eigvec(:,:)
    character(*), intent(in), optional:: units_in
    
    ! loop variables
    integer::i,j
    
    ! other variables
    real(kind=wp),allocatable:: x(:),y(:), diag(:),subdiag(:)
    real(kind=wp),allocatable:: x_new(:),y_new(:)
    real(kind=wp)::dx !, hbar
    character(80):: units    
    integer:: nin

    !real(kind=wp)::dlamch
    
    if(present(units_in)) then
      units = units_in
    else
      units = "SI"
    end if
    
    nin = size(x_in)
    
    allocate(x(nin),y(nin), x_new(npoints), &
         y_new(npoints))
    
    if(allocated(eigval)) deallocate(eigval)
    allocate(eigval(nstates))

    if(allocated(eigvec)) deallocate(eigvec)
    allocate(eigvec(npoints,nstates))

    x = x_in
    y = y_in
    
    ! spline, make more points
    call linspace(x_new, x(1), x(nin), npoints)
    dx = x_new(2)-x_new(1) 
    call spline_easy(x, y, nin, x_new, y_new, npoints)

    call solve_finite_diff(dx, y_new, nstates, eigval, eigvec, mu, units_in)

  
  end subroutine solve_finite_diff_spline

  subroutine solve_finite_diff(dx, y_new, nstates, eigval, eigvec, mu, units_in)
    use m_precision, only: wp
    use m_constants, only: const
    !use m_splines, only: linspace, spline_easy
    implicit none
    
    real(kind=wp),intent(in)::  dx
    real(kind=wp), intent(in):: y_new(:)
    integer, intent(in):: nstates
    real(kind=wp),intent(in)::  mu
    real(wp), intent(out), allocatable:: eigval(:), eigvec(:,:)
    character(*), intent(in), optional:: units_in
    
    ! loop variables
    integer::i,j
    
    ! other variables
    real(kind=wp),allocatable:: diag(:),subdiag(:)
    !real(kind=wp),allocatable:: x_new(:),y_new(:)
    real(kind=wp):: hbar
    character(80):: units    
    integer:: npoints

    ! lapack
    real(kind=wp):: abstol
    integer::n_eigfound,info
    integer, allocatable::iwork(:),ifail(:)
    real(kind=wp), allocatable:: eigenval(:), work(:)
    real(kind=wp), allocatable:: eigenvec(:,:)
    
    !functions
    real(kind=wp)::dlamch
    
    if(present(units_in)) then
      units = units_in
    else
      units = "SI"
    end if

    if (units .eq. "SI") then
      hbar = const % hbar
    else if (units .eq. "AU") then
      hbar =1.0_wp
    end if

    npoints = size(y_new)
    
    allocate(iwork(5*npoints), work(5*npoints), eigenval(npoints), &
         eigenvec(npoints,nstates),ifail(npoints))

    allocate(diag(npoints),subdiag(npoints-1))
    
    if(allocated(eigval)) deallocate(eigval)
    allocate(eigval(nstates))

    if(allocated(eigvec)) deallocate(eigvec)
    allocate(eigvec(npoints,nstates))

    !x = x_in
    !y = y_in
    

    ! sätt upp hamiltonianen -hbar^2/(2*m) (d/dx)^2 + V(x) 
    ! diagonalen
    do i=1,npoints
      diag(i) = (-hbar**2/(dx**2 * 2.0_wp * mu) )*(-2.0_wp) + y_new(i)
    end do
    
    ! subdiagonalen
    do i=1,npoints-1
      subdiag(i) = -hbar**2/(dx**2 *2.0_wp *mu) 
    end do
    
    ! solve eigenvalue problem
    abstol=2d0*dlamch('s')
    call dstevx('v','i',npoints,diag,subdiag, 0.d0,1.d0, 1, nstates, abstol, &
         n_eigfound, eigenval , eigenvec, npoints,work,iwork,ifail,info)
    
    eigvec(:,1:nstates) = eigenvec(:,1:nstates)
    eigval(1:nstates) = eigenval(1:nstates)

  end subroutine solve_finite_diff


end module m_vib_finite_diff
