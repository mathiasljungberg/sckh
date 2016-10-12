module m_spectrum_utils
  implicit none
contains

  ! this subroutine takes a stick spectrum and convolutes it with a Lorentzian
 subroutine convolution_lorentzian(energies_in, intensities_in, fwhm, energies_out, intensities_out)

   use m_constants, only: const

   real(8), intent(in):: energies_in(:)
   real(8), intent(in):: intensities_in(:)
   real(8), intent(in):: fwhm
   real(8), intent(in):: energies_out(:)
   real(8), intent(out):: intensities_out(:)
   
   integer:: i
!   real(8):: dx_in, dx_out
   
   !convolution
   !dx_in = energies_in(2)-energies_in(1)
   !dx_out = energies_out(2)-energies_out(1)
   do i=1,size(energies_out)
     intensities_out(i) = (0.5d0 * fwhm / const %pi) * &
          sum( intensities_in(:) / ((energies_in(:)-energies_out(i))**2 +(0.5d0 * fwhm)**2 )) 
   end do
   
   !intensities_out = intensities_out * (dx_out / dx_in)

 end subroutine convolution_lorentzian

 ! this subroutine takes a spectrum on a grid and convolutes it with a Lorentzian
 subroutine convolution_lorentzian_grid(energies_in, intensities_in, fwhm, energies_out, intensities_out)

   use m_constants, only: const

   real(8), intent(in):: energies_in(:)
   real(8), intent(in):: intensities_in(:)
   real(8), intent(in):: fwhm
   real(8), intent(in):: energies_out(:)
   real(8), intent(out):: intensities_out(:)
   
   integer:: i
   real(8):: dx_in, dx_out
   
   !convolution
   dx_in = energies_in(2)-energies_in(1)
   dx_out = energies_out(2)-energies_out(1)
   do i=1,size(energies_out)
     intensities_out(i) = (0.5d0 * fwhm / const %pi) * &
          sum( intensities_in(:) / ((energies_in(:)-energies_out(i))**2 +(0.5d0 * fwhm)**2 )) 
   end do
   
   !intensities_out = intensities_out * (dx_out / dx_in)
   intensities_out = intensities_out * dx_out

 end subroutine convolution_lorentzian_grid

 ! this subroutine takes a spectrum on a grid and convolutes it with a Lorentzian, with a linear brodeneing starting
 ! at linbroad_start and finishing at linbroad_end. 
 ! care is taken to normalize the numerical integral of the lorentzian since if it is too narrow the numerical integral
 ! will differ a lot to the analytic one
 ! loop order changed so that the linear broadening is in the original spectrum
 subroutine convolution_grid_linear(energies_in, intensities_in, fwhm1, fwhm2, linbroad_start, linbroad_end,&
      energies_out, funct_type, intensities_out)
   
   use m_constants, only: const
   
   real(8), intent(in):: energies_in(:)
   real(8), intent(in):: intensities_in(:)
   real(8), intent(in):: fwhm1
   real(8), intent(in):: fwhm2
   real(8), intent(in):: linbroad_start
   real(8), intent(in):: linbroad_end
   real(8), intent(in):: energies_out(:)
   character(*), intent(in):: funct_type
   real(8), intent(out):: intensities_out(:)
   
   integer:: i
   real(8):: dx_in, dx_out, fwhm, norm, fwhm_tol
   real(8), allocatable:: funct(:)

   !convolution
   dx_in = energies_in(2)-energies_in(1)
   dx_out = energies_out(2)-energies_out(1)

   fwhm_tol = dx_in * 1d-3
   !do i=1,size(energies_out)
   !  intensities_out(i) = (0.5d0 * fwhm / const %pi) * &
   !       sum( intensities_in(:) / ((energies_in(:)-energies_out(i))**2 +(0.5d0 * fwhm)**2 )) 
   !end do

   allocate(funct(size(energies_out)))

   intensities_out = 0.0d0 
   do i=1,size(energies_in)
     if(energies_in(i) .lt. linbroad_start) then
       !  intensities_out(i) = intensities_in(i)
       fwhm = fwhm1
     else if (energies_in(i) .ge. linbroad_start .and. &
          energies_in(i) .le. linbroad_end) then
       !fwhm = fwhm1 * linbroad_start + (energies_in(i) - linbroad_end - linbroad_start)* (fwhm2-fwhm1)
       fwhm = ((fwhm2-fwhm1)/(linbroad_end -linbroad_start))*(energies_in(i) -linbroad_start) + fwhm1
     else
       fwhm = fwhm2
     end if
     
     !write(6,*) "fwhm", fwhm

     if (fwhm .lt. fwhm_tol) then
       intensities_out(i) = intensities_out(i) + intensities_in(i)
     else
       
       if( funct_type .eq. "lorentzian") then
         funct = 1.0d0 / ( (energies_out(:)-energies_in(i))**2 +(0.5d0 * fwhm)**2)
       else if ( funct_type .eq. "gaussian") then
         funct = exp( -(4.0d0 * log(2.0d0) / fwhm**2) * (energies_out(:)-energies_in(i))**2)
       else
         write(6,*) "Error: convololution_grid_linear: funct_type must be either 'lorentzian' or 'gaussian' "
         stop
       end if

       !write(6,*) "funct", funct
       norm = sum(funct) !sum(1.0d0 / (energies_out(:)-energies_in(i))**2 +(0.5d0 * fwhm)**2 ))
       !write(6,*) "norm", norm
       !stop

       intensities_out(:) = intensities_out(:) + (1.0d0 / norm) * &
            intensities_in(i) * funct(:) !/ ((energies_out(:)-energies_in(i))**2 +(0.5d0 * fwhm)**2 ) 
     end if
   
   end do
   
   !intensities_out = intensities_out * (dx_out / dx_in)
   intensities_out = intensities_out !x* dx_out
   
 end subroutine convolution_grid_linear
 


 ! this subroutine takes a spectrum on a grid and convolutes it with a Lorentzian
 subroutine convolution_lorentzian_grid_fft(energies_in, intensities_in, fwhm, energies_out, intensities_out)

   use m_constants, only: const

   real(8), intent(in):: energies_in(:)
   real(8), intent(in):: intensities_in(:)
   real(8), intent(in):: fwhm
   real(8), intent(in):: energies_out(:)
   real(8), intent(out):: intensities_out(:)
   
   complex(8), allocatable:: intensities_in_cmplx(:)
   complex(8), allocatable:: intensities_out_cmplx(:)
   complex(8), allocatable:: conv_func(:) 

   integer:: i,s1,s2,s3,f1,f2,f3
   real(8):: dx 
   
   !convolution
   s1 = lbound(intensities_in,1)
   f1 = ubound(intensities_in,1)
   s3 = lbound(intensities_out,1)
   f3 = ubound(intensities_out,1)

   ! assume same dx in input and output (for fft)!
   dx = energies_in(s1+1)-energies_in(s1)
   !dx = energies_out(s3+1)-energies_out(s3)
 
   s2 = min(s1-s3, s1-f3, f1-s3, f1-f3)
   f2 = max(s1-s3, s1-f3, f1-s3, f1-f3)

   allocate(conv_func(s2:f2))
   allocate(intensities_in_cmplx(s1:f1))
   allocate(intensities_out_cmplx(s3:f3))

   intensities_in_cmplx = intensities_in
   intensities_out_cmplx = 0.0d0

   do i=s2,f2 
     conv_func(i) = 1.0d0 / cmplx(i*dx, 0.5d0*fwhm, 8 )
   end do
    
   call convolution_fft_zzz(intensities_in_cmplx, conv_func, intensities_out_cmplx)
   
   intensities_out(:) = -dx *aimag(intensities_out_cmplx) / const % pi 

   !do i=1,size(energies_out)
   !  intensities_out(i) = (0.5d0 * fwhm / const %pi) * &
   !       sum( intensities_in(:) / ((energies_in(:)-energies_out(i))**2 +(0.5d0 * fwhm)**2 )) 
   !end do
   
   !intensities_out = intensities_out * (dx_out / dx_in)
   !intensities_out = intensities_out * dx

 end subroutine convolution_lorentzian_grid_fft

 ! this subroutine takes a spectrum on a grid and convolutes it with a Lorentzian
 ! if dim =1, then the first dimension will be convoluted, if dim=2 the second dimension will be
 subroutine convolution_lorentzian_grid_fft_many_freq(energies_in, intensities_in, fwhm, &
      energies_out, intensities_out, dim, funct_type_in, asym_in)

   use m_precision, only: wp
   use m_constants, only: const
   use m_KH_functions, only: gaussian
   use m_upper, only : upper
   
   real(wp), intent(in):: energies_in(:)
   real(wp), intent(in):: intensities_in(:,:)
   real(wp), intent(in):: fwhm
   real(wp), intent(in):: energies_out(:)
   real(wp), intent(out):: intensities_out(:,:)
   integer, intent(in):: dim
   character(*), intent(in), optional:: funct_type_in
   real(wp), intent(in), optional:: asym_in
   
   complex(wp), allocatable:: intensities_in_cmplx(:,:)
   complex(wp), allocatable:: intensities_out_cmplx(:,:)
   complex(wp), allocatable:: conv_func(:) 

   integer:: i,s1,s2,s3,f1,f2,f3
   integer:: dim2
   real(wp):: dx, asym, gamma_FWHM, integral 
   character(80):: funct_type
   
   if(present(funct_type_in)) then
     funct_type = funct_type_in
   else
     funct_type = "LORENTZIAN"
   end if

   if(present(asym_in)) then
     asym = asym_in
   else
     asym = 1.0_wp
   end if
   
   if(dim .eq. 1) then
     dim2 = 2 
   else if(dim .eq. 2) then
     dim2=1
   else
     write(6,*) "convolution_lorentzian_grid_fft_many_freq: dim must be either 1 or 2"
     stop
   end if
      
   !convolution
   s1 = lbound(intensities_in,dim)
   f1 = ubound(intensities_in,dim)
   s3 = lbound(intensities_out,dim)
   f3 = ubound(intensities_out,dim)
   
   ! assume same dx in input and output (for fft)!
   if(dim .eq. 1) then
     dx = energies_in(s1+1)-energies_in(s1)
   else if(dim .eq. 2) then
     dx = energies_out(s1+1)-energies_out(s1)
   end if
     
   s2 = min(s1-s3, s1-f3, f1-s3, f1-f3)
   f2 = max(s1-s3, s1-f3, f1-s3, f1-f3)
   
   allocate(conv_func(s2:f2))
   allocate(intensities_in_cmplx(size(intensities_in, dim2), s1:f1))
   allocate(intensities_out_cmplx(size(intensities_in,dim2), s3:f3))

   if (dim .eq. 1) then
     intensities_in_cmplx = transpose(intensities_in)
   else if(dim .eq. 2) then
     intensities_in_cmplx = intensities_in
   end if
     
   intensities_out_cmplx = 0.0d0

   if(upper(funct_type) .eq. "LORENTZIAN") then 
     do i=s2,f2 
       !conv_func(i) = 1.0d0 / cmplx(i*dx, 0.5d0*fwhm, 8 )
       conv_func(i) = (0.5_wp *  fwhm / const % pi)  / ((i*dx)**2+ (0.5_wp * fwhm)**2)
       !conv_func(i) = (0.5_wp *  fwhm )  / ((i*dx)**2+ (0.5_wp * fwhm)**2)
     end do
   elseif(upper(funct_type) .eq. "ASYM_LORENTZIAN") then 
     do i=s2,f2 
       !conv_func(i) = 1.0d0 / cmplx(i*dx, 0.5d0*fwhm, 8 )
       gamma_FWHM = (2.0_wp * fwhm) / (1 + exp(asym * i * dx)) 
       conv_func(i) = (0.5_wp *  gamma_FWHM / const % pi)  / ((i*dx)**2+ (0.5_wp * gamma_FWHM)**2)
       !conv_func(i) = (0.5_wp *  gamma_FWHM )  / ((i*dx)**2+ (0.5_wp * gamma_FWHM)**2)
     end do
   else if(upper(funct_type) .eq. "GAUSSIAN") then 
     do i=s2,f2 
       conv_func(i) = gaussian(i*dx, 0.0_wp, fwhm)
     end do
   else if(upper(funct_type) .eq. "ASYM_GAUSSIAN") then 
     do i=s2,f2 
       gamma_FWHM = (2.0_wp * fwhm) / (1 + exp(asym * i * dx)) 
       if(gamma_FWHM .lt. 1d-10) then
         conv_func(i) =0.0_wp
       else
         conv_func(i) = gaussian(i*dx, 0.0_wp, gamma_FWHM)
         write(6,*) gamma_FWHM, conv_func(i)
       end if
     end do

   else if(upper(funct_type) .eq. "ASYM_GAUSSIAN") then 
     do i=s2,f2 
       gamma_FWHM = (2.0_wp * fwhm) / (1 + exp(asym * i * dx)) 
       if(gamma_FWHM .lt. 1d-10) then
         conv_func(i) =0.0_wp
       else
         conv_func(i) = gaussian(i*dx, 0.0_wp, gamma_FWHM)
         write(6,*) gamma_FWHM, conv_func(i)
       end if
     end do
   else
     write(6,*) "convolution_lorentzian_grid_fft_many_freq: mode must be either 'LORENTZIAN' or 'GAUSSIAN''"
   end if

   ! normalize conv_func
   integral = real(sum(conv_func),8)*dx
   conv_func = conv_func / integral
   
   call convolution_a_fft_zzz(conv_func, intensities_in_cmplx, intensities_out_cmplx)
   
   if (dim .eq. 1) then
     intensities_out = dx * transpose(intensities_out_cmplx)
   else if(dim .eq. 2) then
     intensities_out = dx * intensities_out_cmplx
   end if
   
 end subroutine convolution_lorentzian_grid_fft_many_freq
 
 
 

! ! this subroutine takes a spectrum on a grid and convolutes it with a Lorentzian
! subroutine convolution_lorentzian_grid_fft_dd(energies_in, intensities_in, fwhm, energies_out, intensities_out)
!
!   use m_constants, only: const
!
!   real(8), intent(in):: energies_in(:)
!   real(8), intent(in):: intensities_in(:)
!   real(8), intent(in):: fwhm
!   real(8), intent(in):: energies_out(:)
!   real(8), intent(out):: intensities_out(:)
!   
!   !complex(8), allocatable:: intensities_in(:)
!   !complex(8), allocatable:: intensities_out(:)
!   !complex(8), allocatable:: conv_func(:) 
!
!   integer:: i,s1,s2,s3,f1,f2,f3
!   real(8):: dx 
!   
!   !convolution
!   s1 = lbound(intensities_in,1)
!   f1 = ubound(intensities_in,1)
!   s3 = lbound(intensities_out,1)
!   f3 = ubound(intensities_out,1)
!
!   ! assume same dx in input and output (for fft)!
!   dx = energies_in(s1+1)-energies_in(s1)
!   !dx = energies_out(s3+1)-energies_out(s3)
! 
!   s2 = min(s1-s3, s1-f3, f1-s3, f1-f3)
!   f2 = max(s1-s3, s1-f3, f1-s3, f1-f3)
!
!   allocate(conv_func(s2:f2))
!   !allocate(intensities_in_cmplx(s1:f1))
!   !allocate(intensities_out_cmplx(s3:f3))
!
!   !intensities_in_cmplx = intensities_in
!   !intensities_out_cmplx = 0.0d0
!
!   do i=s2,f2 
!     conv_func(i) = 0.5d0*fwhm / ((i*dx)**2 + (0.5d0*fwhm)**2) !cmplx(i*dx, 0.5d0*fwhm, 8 )
!   end do
!    
!   call convolution_fft_zzz(intensities_in, conv_func, intensities_out)
!   
!   !intensities_out(:) = -dx *aimag(intensities_out_cmplx) / const % pi
!   intensities_out(:) = dx * intensities_out(:) / const % pi 
!
!   !do i=1,size(energies_out)
!   !  intensities_out(i) = (0.5d0 * fwhm / const %pi) * &
!   !       sum( intensities_in(:) / ((energies_in(:)-energies_out(i))**2 +(0.5d0 * fwhm)**2 )) 
!   !end do
!   
!   !intensities_out = intensities_out * (dx_out / dx_in)
!   !intensities_out = intensities_out * dx
!
! end subroutine convolution_lorentzian_grid_fft_dd
 
!subroutine convolution_fft_ddd(a,b,c)
!  use m_precision, only : fftw_real
!  use m_convolution, only: Convolution_limits
!  use m_convolution, only: Convolution_calc_dimensions
!  use m_convolution, only: Convolution_create_c2c_plans
!  use m_convolution, only: Convolution_apply_fft_8
!  use m_convolution, only: Convolution_multiply_and_apply_ifft_8
!  use m_convolution, only: Convolution_slow_zzz
!
!  implicit none
!
!  real(8), intent(in), allocatable ::a(:), b(:)
!  real(8), intent(out) ::c(:)
!
!  real(8), allocatable :: cf(:)
!  integer :: s1,f1,n1, s2,f2,n2,s3,f3,ncaux,saux
!  integer :: salloc,falloc
!
!!  integer:: n!,i1,i2,i3
!
!  !! For fast convolution
!  integer    :: nff_conv, nff_conv_r2c
!  integer(8) :: fft_plan, ifft_plan
!  complex(fftw_real), allocatable :: aux_c2c_1(:)
!  complex(fftw_real), allocatable :: aux_c2c_2(:)
!
!  complex(fftw_real), allocatable :: a_r2c(:)
!  complex(fftw_real), allocatable :: b_r2c(:)
!  !! END of For fast convolution
!
!  ! s1: starting index for a, f1: final index for a
!  ! s2: starting index for b, f2: final index for b 
!  ! s3: starting index for c, f3: final index for c 
!
!  ! remember that the dimensions of a,b,c must conform! Look at the integral
!  ! c(\omega) = \int a(\omega-\omega') b(\omega') d\omega
!  ! and from the limits of for example c and b, the limits of a are deduced
!  ! these limits can be negative.
!
!  s1 = lbound(a,1)
!  f1 = ubound(a,1)
!  s2 = lbound(b,1)
!  f2 = ubound(b,1)
!  s3 = lbound(c,1)
!  f3 = ubound(c,1)
!    
!  !! Fast way of doind the same calculation
!  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
!  allocate(cf(salloc:falloc))
!  cf = 0.0d0
!  
!  !! Take care of complex to complex transform(s)
!  call Convolution_calc_r2c_dimensions(n1,n2,nff_conv,nff_conv_r2c)
!  allocate(aux_c2c_1(nff_conv))
!  allocate(aux_c2c_2(nff_conv))
!  allocate(a_c2c(nff_conv))
!  allocate(b_c2c(nff_conv))
!  call Convolution_create_r2c_plans(nff_conv, ifft_plan, fft_plan,aux_c2c_1,aux_c2c_2)
!  !! END of Take care of complex to complex transform(s)
!  
!  call Convolution_apply_r2c_8(n1, a(s1:), a_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)
!  call Convolution_apply_r2c_8(n2, b(s2:), b_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)
!  call Convolution_multiply_and_apply_c2r_8(nff_conv, a_c2c, b_c2c, ncaux, cf(saux:), &
!    ifft_plan, aux_c2c_1, aux_c2c_2);
!
!  c(s3:f3) = cf(s3:f3) 
!
!  deallocate(cf)
!  deallocate(aux_c2c_1)
!  deallocate(aux_c2c_2)
!  deallocate(a_c2c)
!  deallocate(b_c2c)
!
!end subroutine convolution_fft_zzz
 


subroutine convolution_fft_zzz(a,b,c)
  use m_precision, only : fftw_real
  use m_convolution, only: Convolution_limits
  use m_convolution, only: Convolution_calc_dimensions
  use m_convolution, only: Convolution_create_c2c_plans
  use m_convolution, only: Convolution_apply_fft_8
  use m_convolution, only: Convolution_multiply_and_apply_ifft_8
  use m_convolution, only: Convolution_slow_zzz

  implicit none

  complex(8), intent(in), allocatable ::a(:), b(:)
  complex(8), intent(out) ::c(:)

  complex(8), allocatable :: cf(:)
  integer :: s1,f1,n1, s2,f2,n2,s3,f3,ncaux,saux
  integer :: salloc,falloc

!  integer:: n!,i1,i2,i3

  !! For fast convolution
  integer    :: nff_conv, nff_conv_r2c
  integer(8) :: fft_plan, ifft_plan
  complex(fftw_real), allocatable :: aux_c2c_1(:)
  complex(fftw_real), allocatable :: aux_c2c_2(:)

  complex(fftw_real), allocatable :: a_c2c(:)
  complex(fftw_real), allocatable :: b_c2c(:)
  !! END of For fast convolution

  ! s1: starting index for a, f1: final index for a
  ! s2: starting index for b, f2: final index for b 
  ! s3: starting index for c, f3: final index for c 

  ! remember that the dimensions of a,b,c must conform! Look at the integral
  ! c(\omega) = \int a(\omega-\omega') b(\omega') d\omega
  ! and from the limits of for example c and b, the limits of a are deduced
  ! these limits can be negative.

  s1 = lbound(a,1)
  f1 = ubound(a,1)
  s2 = lbound(b,1)
  f2 = ubound(b,1)
  s3 = lbound(c,1)
  f3 = ubound(c,1)

  !write(6,*) "size(a)", size(a)
  !write(6,*) "s1,f2,s2,f2,s3,f3",s1,f1,s2,f2,s3,f3

 !!! Doing the same calculation but using "standart" convolution
 !call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
 !
 !write(6,*) "n1, n2, ncaux,saux,salloc,falloc",n1, n2, ncaux,saux,salloc,falloc
 !
 !allocate(cf(salloc:falloc))
 !cf = 0
 !call Convolution_slow_zzz(n1,a(s1), n2,b(s2),ncaux,cf(saux))

  !c=0.0d0
  !n = size(a)
  !do i1 = 1, n!size(a)
  !  do i2 = 1, n!size(a)
  !    i3 = mod((n-1)/2+ (i1+i2) -1,n) +1 
  !    !write(6,*) i1,i2, i3
  !    c(i3) = c(i3) + a(i1)*b(i2)
  !    !i3 = mod( i1 - i2, n) +1
  !    !c(i1) = c(i1) + a(i2)*b(i3)
  !  end do
  !end do

  !if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  !if(iv>0)write(6,*)  salloc, falloc, saux, ncaux
  !if(iv>0)write(6,*) 'slow'
  !write(6,'(i5,4g25.16)') &
  !  (k, real(c(k)), real(caux(k)), aimag(c(k)), aimag(caux(k)),k=s3,f3);
  !!! END of Doing the same calculation but using "standart" convolution
    
  !! Fast way of doind the same calculation
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(cf(salloc:falloc))
  cf = 0.0d0
  
  !write(6,*) "salloc, falloc", salloc, falloc
  
  !! Take care of complex to complex transform(s)
  call Convolution_calc_dimensions(n1,n2,nff_conv,nff_conv_r2c)
  allocate(aux_c2c_1(nff_conv))
  allocate(aux_c2c_2(nff_conv))
  allocate(a_c2c(nff_conv))
  allocate(b_c2c(nff_conv))
  call Convolution_create_c2c_plans(nff_conv, ifft_plan, fft_plan,aux_c2c_1,aux_c2c_2)
  !! END of Take care of complex to complex transform(s)
  
  call Convolution_apply_fft_8(n1, a(s1:), a_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)
  call Convolution_apply_fft_8(n2, b(s2:), b_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)
  call Convolution_multiply_and_apply_ifft_8(nff_conv, a_c2c, b_c2c, ncaux, cf(saux:), &
    ifft_plan, aux_c2c_1, aux_c2c_2);

  c(s3:f3) = cf(s3:f3) 

  deallocate(cf)
  deallocate(aux_c2c_1)
  deallocate(aux_c2c_2)
  deallocate(a_c2c)
  deallocate(b_c2c)

end subroutine convolution_fft_zzz

subroutine convolution_a_fft_zzz(a,b,c)
  use m_precision, only : fftw_real
  use m_convolution, only: Convolution_limits
  use m_convolution, only: Convolution_calc_dimensions
  use m_convolution, only: Convolution_create_c2c_plans
  use m_convolution, only: Convolution_apply_fft_8
  use m_convolution, only: Convolution_multiply_and_apply_ifft_8
  use m_convolution, only: Convolution_slow_zzz

  implicit none

  complex(8), intent(in), allocatable ::a(:), b(:,:)
  complex(8), intent(out) ::c(:,:)

  complex(8), allocatable :: cf(:)
  integer :: i,s1,f1,n1, s2,f2,n2,s3,f3,ncaux,saux
  integer :: salloc,falloc

  !! For fast convolution
  integer    :: nff_conv, nff_conv_r2c
  integer(8) :: fft_plan, ifft_plan
  complex(fftw_real), allocatable :: aux_c2c_1(:)
  complex(fftw_real), allocatable :: aux_c2c_2(:)

  complex(fftw_real), allocatable :: a_c2c(:)
  complex(fftw_real), allocatable :: b_c2c(:)
  complex(fftw_real), allocatable :: b_tmp(:)
  !! END of For fast convolution

  ! s1: starting index for a, f1: final index for a
  ! s2: starting index for b, f2: final index for b 
  ! s3: starting index for c, f3: final index for c 

  ! remember that the dimensions of a,b,c must conform! Look at the integral
  ! c(\omega) = \int a(\omega-\omega') b(\omega') d\omega
  ! and from the limits of for example c and b, the limits of a are deduced
  ! these limits can be negative.

  s1 = lbound(a,1)
  f1 = ubound(a,1)
  s2 = lbound(b,2)
  f2 = ubound(b,2)
  s3 = lbound(c,2)
  f3 = ubound(c,2)

  !! Fast way of doind the same calculation
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(cf(salloc:falloc))
  cf = 0.0d0
  c = 0.0d0
  
  !write(6,*) "salloc, falloc", salloc, falloc
  
  !! Take care of complex to complex transform(s)
  call Convolution_calc_dimensions(n1,n2,nff_conv,nff_conv_r2c)
  allocate(aux_c2c_1(nff_conv))
  allocate(aux_c2c_2(nff_conv))
  allocate(a_c2c(nff_conv))
  allocate(b_c2c(nff_conv))
  allocate(b_tmp(nff_conv))
  call Convolution_create_c2c_plans(nff_conv, ifft_plan, fft_plan,aux_c2c_1,aux_c2c_2)
  !! END of Take care of complex to complex transform(s)

  call Convolution_apply_fft_8(n1, a(s1:), a_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)

  do i=1,size(b,1)
    !do j=1,size(b,2)
      
      b_tmp(s2:f2) = b(i,s2:f2)
      call Convolution_apply_fft_8(n2, b_tmp(s2:), b_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)
      call Convolution_multiply_and_apply_ifft_8(nff_conv, a_c2c, b_c2c, ncaux, cf(saux:), &
           ifft_plan, aux_c2c_1, aux_c2c_2);
      
      c(i,s3:f3) = cf(s3:f3) 

    !end do
  end do

  call dfftw_destroy_plan(fft_plan)
  call dfftw_destroy_plan(ifft_plan)
  
  deallocate(cf)
  deallocate(aux_c2c_1)
  deallocate(aux_c2c_2)
  deallocate(a_c2c)
  deallocate(b_c2c)

end subroutine convolution_a_fft_zzz


subroutine convolution_ab_fft_zzz(a,b,c)
  use m_precision, only : fftw_real
  use m_convolution, only: Convolution_limits
  use m_convolution, only: Convolution_calc_dimensions
  use m_convolution, only: Convolution_create_c2c_plans
  use m_convolution, only: Convolution_apply_fft_8
  use m_convolution, only: Convolution_multiply_and_apply_ifft_8
  use m_convolution, only: Convolution_slow_zzz

  implicit none

  complex(8), intent(in), allocatable ::a(:), b(:,:,:)
  complex(8), intent(out) ::c(:,:,:)

  complex(8), allocatable :: cf(:)
  integer :: i,j,s1,f1,n1, s2,f2,n2,s3,f3,ncaux,saux
  integer :: salloc,falloc

  !! For fast convolution
  integer    :: nff_conv, nff_conv_r2c
  integer(8) :: fft_plan, ifft_plan
  complex(fftw_real), allocatable :: aux_c2c_1(:)
  complex(fftw_real), allocatable :: aux_c2c_2(:)

  complex(fftw_real), allocatable :: a_c2c(:)
  complex(fftw_real), allocatable :: b_c2c(:)
  complex(fftw_real), allocatable :: b_tmp(:)
  !! END of For fast convolution

  ! s1: starting index for a, f1: final index for a
  ! s2: starting index for b, f2: final index for b 
  ! s3: starting index for c, f3: final index for c 

  ! remember that the dimensions of a,b,c must conform! Look at the integral
  ! c(\omega) = \int a(\omega-\omega') b(\omega') d\omega
  ! and from the limits of for example c and b, the limits of a are deduced
  ! these limits can be negative.

  s1 = lbound(a,1)
  f1 = ubound(a,1)
  s2 = lbound(b,3)
  f2 = ubound(b,3)
  s3 = lbound(c,3)
  f3 = ubound(c,3)

  !! Fast way of doind the same calculation
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(cf(salloc:falloc))
  cf = 0.0d0
  c = 0.0d0
  
  !write(6,*) "salloc, falloc", salloc, falloc
  
  !! Take care of complex to complex transform(s)
  call Convolution_calc_dimensions(n1,n2,nff_conv,nff_conv_r2c)
  allocate(aux_c2c_1(nff_conv))
  allocate(aux_c2c_2(nff_conv))
  allocate(a_c2c(nff_conv))
  allocate(b_c2c(nff_conv))
  allocate(b_tmp(nff_conv))
  call Convolution_create_c2c_plans(nff_conv, ifft_plan, fft_plan,aux_c2c_1,aux_c2c_2)
  !! END of Take care of complex to complex transform(s)

  call Convolution_apply_fft_8(n1, a(s1:), a_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)

  do i=1,size(b,1)
    do j=1,size(b,2)
      
      b_tmp(s2:f2) = b(i,j,s2:f2)
      call Convolution_apply_fft_8(n2, b_tmp(s2:), b_c2c, nff_conv,fft_plan,aux_c2c_1,aux_c2c_2)
      call Convolution_multiply_and_apply_ifft_8(nff_conv, a_c2c, b_c2c, ncaux, cf(saux:), &
           ifft_plan, aux_c2c_1, aux_c2c_2);
      
      c(i,j,s3:f3) = cf(s3:f3) 

    end do
  end do

  call dfftw_destroy_plan(fft_plan)
  call dfftw_destroy_plan(ifft_plan)
  
  deallocate(cf)
  deallocate(aux_c2c_1)
  deallocate(aux_c2c_2)
  deallocate(a_c2c)
  deallocate(b_c2c)

end subroutine convolution_ab_fft_zzz




end module m_spectrum_utils
