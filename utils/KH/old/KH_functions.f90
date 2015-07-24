module KH_functions
  use parameters
  implicit none

contains

  real(kind=wp) function gaussian(x, x0, FWHM)
    real(wp), intent(in):: x, x0, FWHM
    !real(wp) :: gaussian
    real(wp) :: alpha

    alpha=4.0_wp*log(2.0_wp)/(FWHM**2.0_wp)
    gaussian = (alpha / pi)**0.5_wp * exp(-alpha*(x-x0)**2.0_wp)

  end function gaussian

  subroutine solve_sinc_DVR(dx,my_SI, V_i, c_i, e_i)
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
       H_kin(ii,ii) = (hbar**2  / (2 * my_SI * dx **2)) * (pi ** 2) /3.0d0
    end do

    do i = -npoints,npoints
       ii = i + npoints +1
       do j = i +1, npoints
          jj = j + npoints +1
          H_kin(ii,jj) =  (hbar**2 * (-1)**(i-j) / ( my_SI * dx **2)) / (i-j) **2
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
    use parameters
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
          write(6,*) "dvr point higher than supplied potential range!", X_dvr(nstates), x_in(npoints_in), X_dvr(nstates)- x_in(npoints_in)
          stop
       end if
    end if

  end subroutine check_dvr_bounds


end module KH_functions

