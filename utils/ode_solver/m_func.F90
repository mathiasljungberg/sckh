module m_func
  use m_precision, only: wp
  implicit none

  real(kind=wp), allocatable, private:: H_in(:,:,:)
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

end module m_func
