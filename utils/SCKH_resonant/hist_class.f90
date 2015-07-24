! last changed 2007-10-10

module hist_class
  use parameters
  implicit none
  type hist
     integer:: nbins
     real(kind=wp),dimension(:),allocatable :: x,y     
     real(kind=wp):: start,end, interval, dx
  end type hist

contains

!constructor
  subroutine hist_init(a, nb, start, end)
    !passed variables
    type (hist),intent(out):: a
    integer,intent(in):: nb
    real(kind=wp):: start, end

    !local variables
    integer:: i

    !set class variables
    a%nbins=nb
    
    allocate( a%x(nb), a%y(nb) )
  
    a%y=0
    a%start = start
    a%end = end
    a%interval= end-start
    a%dx = (end - start)/dfloat(nb-1)
  
    ! points defined at left of interval
    do i=1,nb
       a%x(i)= start + (i-1)*a%dx
    end do

  end subroutine hist_init

!destructor
  subroutine hist_final(a)
    !passed variables
    type (hist),intent(out):: a
 
    deallocate(a%x)
    deallocate(a%y)

  end subroutine hist_final

! add point to histogram
  subroutine hist_add(a, xin, yin)
    !passed variables
    type (hist),intent(inout):: a
    real(kind=wp):: xin,yin
    
    !local variables
    integer:: rx
    
    rx =int( (xin -a%start )/ a%dx) + 1 
    
    ! force them inte a bin
    if(rx.lt.1) rx=1
    if(rx.ge.a%nbins) rx=a%nbins
    
    !add intensities
    a%y(rx) = a%y(rx)+ yin 
   
  end subroutine hist_add
  
  ! write histogram to file
  subroutine hist_write(a, outfile)
    !passed variables
    type (hist),intent(inout):: a
    character(LEN=*):: outfile
    
    !local variables
    integer:: i
    
    open(10,file=outfile,status='unknown')
    do i=1,a%nbins
       write(10,'((ES16.6E3) (ES16.6E3))') a%x(i), a%y(i)
    end do
    close(10)

  end subroutine hist_write

  subroutine hist_broadening(a,fwhm)
    !passed variables
    type (hist),intent(inout):: a
    real(kind=wp)::fwhm
    !local variables
    real(kind=wp), dimension(:), allocatable::intenshi
    real(kind=wp):: alpha, sum1, sum2
    integer::i,j

    ! calculate sum of initial histogram
    sum1=0
    do i=1,a%nbins
       sum1 = sum1 + a%y(i) 
    end do

    allocate( intenshi(a%nbins) )
    intenshi=0
    alpha=4.0_wp*log(2.0_wp)/(fwhm**2.0_wp)

    !broadening with gaussian
    do i=1,a%nbins
       do j=1,a%nbins
          intenshi(i) = intenshi(i) &
               + a%y(j) * (alpha / pi)**0.5_wp * exp(-alpha*(a%x(i)-a%x(j))**2.0_wp)
       end do
    end do
    
    ! calculate sum of broadened histogram
    sum2=0
    do i=1,a%nbins
       sum2 = sum2 + intenshi(i) 
    end do

    ! set histogram to normalized broadened histogram
    do i=1,a%nbins
       a%y(i) = intenshi(i) *(sum1/sum2)
    end do
    deallocate(intenshi)


  end subroutine hist_broadening


    subroutine normalize_integral(a, norm)
      !passed variables
      type (hist),intent(inout):: a
      real(kind=wp):: norm
 
      !local variables
      integer:: i
      real(kind=wp):: sum2

      !calculate integral      
      sum2=0
      do i=1,a%nbins
         sum2 = sum2 + a%y(i) 
      end do
      sum2=sum2*a%dx

      !normalize
      do i=1,a%nbins
         a%y(i) = a%y(i)*(norm/sum2) 
      end do

    end subroutine normalize_integral

    subroutine normalize_sum(a, norm)
      !passed variables
      type (hist),intent(inout):: a
      real(kind=wp):: norm
 
      !local variables
      integer:: i
      real(kind=wp):: sum2

      !calculate integral      
      sum2=0
      do i=1,a%nbins
         sum2 = sum2 + a%y(i) 
      end do

      !normalize
      do i=1,a%nbins
         a%y(i) = a%y(i)*(norm/sum2) 
      end do

    end subroutine normalize_sum

     real(kind=wp) function get_sum(a)
      !passed variables
      type (hist),intent(inout):: a

       !local variables
      integer:: i
      real(kind=wp):: sum

      !calculate integral      
      sum=0
      do i=1,a%nbins
         sum = sum + a%y(i) 
      end do
      get_sum=sum

    end function get_sum




end module hist_class
