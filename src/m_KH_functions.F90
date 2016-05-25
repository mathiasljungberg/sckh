module m_KH_functions
  !use parameters
  !use spline_m
  implicit none

contains

  function delta(i,j)
    implicit none
    
    integer:: delta
    integer,intent(in):: i,j
    
    if(i.eq.j) then
       delta = 1
    else
       delta = 0
    end if
    
  end function delta

  real(kind=wp) function gaussian(x, x0, FWHM)
    use m_precision, only: wp
    use m_constants, only: const

    real(wp), intent(in):: x, x0, FWHM
    !real(wp) :: gaussian
    real(wp) :: alpha

    alpha=4.0_wp*log(2.0_wp)/(FWHM**2.0_wp)
    gaussian = (alpha / const % pi)**0.5_wp * exp(-alpha*(x-x0)**2.0_wp)

  end function gaussian

  function  sinc(x, x_n, dx) 
    use m_precision, only: wp
    use m_constants, only: const

    real(kind=wp):: sinc
    real(kind=wp), intent(in):: x,x_n,dx
    
    ! normalized sinc function  =   sqrt(dx)*sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))
    ! J. Chem. Phys. 116, 7350 (2002); DOI:10.1063/1.1467055 
    
    ! case when x close to X_n  (expand in taylor series)
    if ( abs((x-x_n)/dx) .lt. 1.0d-13) then
       sinc = sqrt(dx)* ( 1.0d0/dx  - (const % pi**2/(6*dx**3))*(x-x_n)**2 )
       
       !write(7,*) (x-x_n), sinc, sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))
       
       !otherwise 
    else
       sinc = sqrt(dx)*sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))
    end if
    
  end function sinc
  
  function  dsinc(x, x_n, dx) 
    use m_precision, only: wp
    use m_constants, only: const

    real(kind=wp):: dsinc
    real(kind=wp), intent(in):: x,x_n,dx
    
    ! derivative of normalized sinc function  = D/Dx  sqrt(dx)*sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))
    ! =  const % pi * ((1/sqrt(dx)) * cos(const % pi * (x-x_n) / dx)/ (const % pi * (x-x_n))  - sqrt(dx) * sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))**2 ) 
    
    ! case when x close to X_n  (expand in taylor series)
    if ( abs(x-x_n) .lt. 1.0d-13) then
       dsinc = sqrt(dx)* (  - (const % pi**2/(3*dx**3))*(x-x_n) )
       
       !otherwise 
    else
       dsinc = ( (1.0d0 / sqrt(dx)) *  cos(const % pi * (x-x_n) / dx) / ( (x-x_n)) - &
            sqrt(dx) * sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n) **2))
    end if
    
  end function dsinc
  
  function  d2sinc(x, x_n, dx) 
    use m_precision, only: wp
    use m_constants, only: const

    real(kind=wp):: d2sinc
    real(kind=wp), intent(in):: x,x_n,dx
    
    ! second derivative of normalized sinc function  = d/dx  sqrt(dx)*sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))
    ! =  const % pi * ((1/sqrt(dx)) * cos(const % pi * (x-x_n) / dx)/ (const % pi * (x-x_n))  - sqrt(dx) * sin(const % pi * (x-x_n) / dx) / (const % pi * (x-x_n))**2 ) 
    
    ! case when x close to X_n  (expand in taylor series)
    if ( abs(x-x_n) .lt. 1.0d-13) then
       d2sinc = sqrt(dx)* const % pi**2 * (  - (1.0d0 / (3*dx**3))+ &
            (1.0d0 / dx**5 ) * (const % pi* (x-x_n))  **2 )
       
       !otherwise 
    else
       d2sinc =  -sqrt(dx) * const % pi**2 * ( (1.0d0 / dx**2) *  sin(const % pi * (x-x_n) / dx) / ( const % pi*(x-x_n) ) + &
            (2.0d0 / dx) * cos(const % pi * (x-x_n) / dx) / ( (const % pi * (x-x_n)) **2) - &
            2.0d0 * sin(const % pi * (x-x_n) / dx) / ( (const % pi * (x-x_n)) **3) )
    end if
    
  end function d2sinc

  subroutine dsinc_mat_elems(M, X_dvr)
    use m_precision, only: wp
    !use m_constants, only: const
    use m_splines, only: linspace

    real(kind=wp), dimension(:,:),intent(out):: M
    real(kind=wp), dimension(:),intent(in):: X_dvr

    integer:: nstates, nx, i,j1,j2
    real(kind=wp):: dx, dx_sinc, prod1
    real(kind=wp), dimension(:), allocatable:: x_sinc

    ! computes matrix elements of dsinc and sinc numerically

    nstates = size(X_dvr)
    
    nx = 10000
    dx = X_dvr(2) -X_dvr(1)

    allocate(x_sinc(nx))

    call linspace(x_sinc,X_dvr(1) -1000*dx,X_dvr(nstates) + 1000*dx, nx)
    dx_sinc = x_sinc(2) - x_sinc(1)

    do j1 = 1, nstates
       do j2 = 1, nstates
          prod1 = 0.0_wp
          do i=1,nx
             
             prod1 =  prod1 + sinc(x_sinc(i), X_dvr(j1),dx) * dsinc(x_sinc(i), X_dvr(j2),dx)
             
          end do
          
!          write(6,*) j1,j2
          
          M(j1,j2) = prod1 * dx_sinc
          
       end do
    end do

    deallocate(x_sinc)

  end subroutine dsinc_mat_elems

  subroutine solve_sinc_DVR(dx,my_SI, V_i, c_i, e_i)
    use m_precision, only: wp
    use m_constants, only: const

    real(kind=wp), intent(in):: dx, my_SI
    real(kind=wp), dimension(:), intent(in):: V_i
    real(kind=wp), dimension(:,:), intent(out):: c_i
    real(kind=wp), dimension(:), intent(out):: e_i

    !! local variables
    integer:: i,j,ii,jj,npoints,nstates
    real(kind=wp), dimension(:,:), allocatable::Mat_tmp, H_kin
    real(kind=wp), dimension(:),allocatable:: W
    real(kind=wp), dimension(:),allocatable:: WORK
    integer:: INFO, LWORK

    ! assume SI units!

    nstates = size(V_i)
    !write(6,*) "solve_sinc_dvr: nstates=", nstates

    LWORK = 3*nstates

    !
    allocate(Mat_tmp(nstates,nstates), H_kin(nstates,nstates), W(nstates),&
         WORK(LWORK))

    npoints = (nstates-1)/2

    H_kin=0 

    do i = -npoints,npoints
       ii = i + npoints +1
       H_kin(ii,ii) = (const % hbar**2  / (2.0_wp * my_SI * dx **2)) * (const % pi ** 2) / 3.0_wp
    end do

    do i = -npoints,npoints
       ii = i + npoints +1
       do j = i +1, npoints
          jj = j + npoints +1
          H_kin(ii,jj) =  (const % hbar**2 * (-1.0_wp)**(i-j) / ( my_SI * dx **2)) / (i-j) **2
          H_kin(jj,ii) = H_kin(ii,jj)
       end do
    end do


    ! potential term                                                                                                                                                                                                    
    Mat_tmp = H_kin 
    !V_dvr = 0
    do i=1,nstates
       !V_dvr(i,i) = V_i(i)*hartree
       Mat_tmp(i,i) = Mat_tmp(i,i) + V_i(i) !*hartree
    end do


    ! solve eigenvalue problem
    call dsyev( "V", "U", nstates, Mat_tmp, nstates, W, WORK, LWORK, INFO)

    ! obs c_gs = U^T, i.e. the transpose unitary matrix. n:th eigenvector: c_gs(:,n)
    c_i = Mat_tmp

    e_i = W  !/ ev


    deallocate(Mat_tmp, H_kin, W, WORK)

  end subroutine solve_sinc_DVR

  subroutine check_dvr_bounds(X_dvr,x_in)
    !use parameters
    use m_precision, only: wp
    implicit none
    real(kind=wp), dimension(:),intent(inout):: X_dvr
    real(kind=wp), dimension(:),intent(in):: x_in

    integer:: npoints_in, nstates

    nstates =size(X_dvr)
    npoints_in = size(x_in)

    if(X_dvr(1).lt.x_in(1)) then
       if(abs(X_dvr(1)-x_in(1)) .lt. 1d-20 ) then
          X_dvr(1)=x_in(1)
          write(6,*) "changed first dvr point!"
       else
          write(6,*) "dvr point lower than supplied potential range!", X_dvr(1), x_in(1)
          stop
       end if
    end if
    if(X_dvr(nstates).gt.x_in(npoints_in)) then
       if(abs(X_dvr(nstates)-x_in(npoints_in)) .lt. 1d-20 ) then
          X_dvr(nstates)=x_in(npoints_in)
          write(6,*) "changed last dvr point!"
       else
          write(6,*) "dvr point higher than supplied potential range!", X_dvr(nstates), x_in(npoints_in), &
               X_dvr(nstates)- x_in(npoints_in)
          stop
       end if
    end if

  end subroutine check_dvr_bounds


end module M_KH_functions

