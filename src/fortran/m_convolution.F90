module m_convolution

  use m_precision, only : fftw_real

  !public Cauchy_driver
  !public Convolution_driver
  !
  !public Convolution_example_cmplx_cmplx
  !public Convolution_example_cmplx_real
  !public Convolution_example_real_real
  !
  !public Convolution_calc_dimensions
  !public Convolution_limits
  !
  !public Convolution_calc_r2c_dimensions
  !public Convolution_create_r2c_plans
  !public Convolution_apply_r2c_4
  !public Convolution_apply_r2c_8
  !
  !public Convolution_apply_c2r
  !public Convolution_apply_c2r_4
  !public Convolution_apply_c2r_8
  !public Convolution_multiply_and_apply_c2r
  !public Convolution_multiply_and_apply_c2r_8
  !
  !public Convolution_calc_c2c_dimensions
  !public Convolution_create_c2c_plans
  !public Convolution_apply_fft_4
  !public Convolution_apply_fft_8
  !public Convolution_multiply_and_apply_ifft_4
  !public Convolution_multiply_and_apply_ifft_8
  !public Convolution_transform_r2c_to_c2c_4
  !public Convolution_transform_r2c_to_c2c_8
  !
  !public Convolution_slow
  !public Convolution_slow_8
  !public Convolution_slow_zdz
  
  interface Convolution_apply_c2r
    module procedure Convolution_apply_c2r_4
    module procedure Convolution_apply_c2r_8
  end interface ! Convolution_apply_c2r
  
  contains

!
!
!
subroutine demonstrate_hartree_pot_1d(iv, ilog)
  implicit none  
  !! external
  integer, intent(in) :: iv, ilog
  !! internal
  integer :: x1,x2,y1,y2,ic1,ic2,ix,iy,ic
  real(8) :: dr = 0.1D0
  real(8), allocatable :: dens(:), vh1(:), vh2(:), vh3(:)
  real(8), allocatable :: coul(:)

  !! For fast convolution
  integer :: n1, n2, ncaux,saux
  integer :: salloc,falloc
  integer :: nff_complex, nff_real
  integer(8) ::  c2r_plan, r2c_plan
  real(fftw_real), allocatable :: inout_real(:)
  complex(fftw_real), allocatable :: inout_complex(:)
  complex(fftw_real), allocatable :: a_r2c(:)
  complex(fftw_real), allocatable :: b_r2c(:)
  !! END of For fast convolution

  x1 = 1; x2 = 10
  y1 = 1; y2 = 10

  allocate(dens(x1:x2))
  allocate(vh1(y1:y2))
  call random_number(dens)
  
  do ix=x1,x2; write(ilog,*) ix, dens(ix); enddo;

  vh1 = 0
  do iy=y1,y2
    do ix=x1,x2
      vh1(iy) = vh1(iy) + dens(ix) / abs(dr*( iy - ix + 0.5D0 ))
    enddo
  enddo

  do iy=y1,y2; write(ilog,*) iy, vh1(iy); enddo;

  !! Construct Coulomb denominator
  ic1 = min(y1-x1, y1-x2, y2-x1, y2-x2)
  ic2 = max(y1-x1, y1-x2, y2-x1, y2-x2)
  allocate(coul(ic1:ic2))
  do ic=ic1,ic2; coul(ic) = 1D0 / abs(dr*( ic + 0.5D0 )); enddo
  !! END of Construct Coulomb denominator

  allocate(vh2(y1:y2))
  !! Double loop which must be accelerated with FFT(W)
  vh2 = 0
  do ix=x1,x2
  do ic=ic1,ic2
    iy = ix+ic
    if(iy<y1 .or. iy>y2) cycle
    vh2(iy) = vh2(iy) + dens(ix)*coul(ic)
  enddo
  enddo
  !! END of Double loop which must be accelerated with FFT(W)
  write(ilog,*) 'vhartree computed with double loop: second column must be zero'
  do iy=y1,y2; write(ilog,*) iy, vh1(iy)-vh2(iy); enddo;

  !! Fast way of doind the same calculation
  call Convolution_limits(x1,x2,ic1,ic2,y1,y2, n1,n2,ncaux,saux,salloc,falloc);
  allocate(vh3(salloc:falloc))
  vh3 = 0
  call Convolution_calc_r2c_dimensions(n1,n2,nff_complex, nff_real)
  allocate(inout_complex(nff_complex))
  allocate(inout_real(nff_real))
  allocate(a_r2c(nff_complex))
  allocate(b_r2c(nff_complex))

  call Convolution_create_r2c_plans(nff_real, c2r_plan, r2c_plan, inout_complex, inout_real);
  call Convolution_apply_r2c_8(n1, dens(x1:), a_r2c, nff_real, r2c_plan, inout_complex, inout_real);
  call Convolution_apply_r2c_8(n2, coul(ic1:), b_r2c, nff_real, r2c_plan, inout_complex, inout_real);
  call Convolution_multiply_and_apply_c2r_8(nff_complex, a_r2c, b_r2c, ncaux, vh3(saux:), &
    c2r_plan, inout_complex, inout_real);

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux
  write(6,*) 'fast'
  do iy=y1,y2
    write(6,*) iy, vh1(iy)-vh3(iy)
  enddo
  !! END of Fast way of doind the same calculation


end subroutine ! demonstrate_hartree_pot_1d


!
! Shows how to use convolutions in order to speedup Cauchy transformation
!
subroutine Cauchy_driver(iv)
  implicit none
  integer, intent(in) :: iv
  ! internal
  real(4), allocatable :: Spectral_funct_of_Response_pos(:)
  complex(4), allocatable :: Response_pos(:), Response_pos_fast(:), Response_pos_slow(:), Response_pos_fast_cmplx(:)
  integer :: f,s
  real(4), parameter :: d_omega=0.1E0, eps=0.1E0
  integer, parameter :: nff_chi=8
  integer, parameter :: nff_a  =8
  integer            :: nff_cauchy
  !! Auxiliary parameters to rewrite Cauchy as convolution
  real(4), allocatable    :: Cauchy_real(:), Cauchy_imag(:), Response_real(:), Response_imag(:)

  !! Auxiliary parameters for fast convolution 
  integer :: nff_complex, nff_real
  integer(8) :: r2c_plan, c2r_plan
  complex(fftw_real), allocatable :: inout_complex(:)
  real(fftw_real), allocatable :: inout_real(:)
  complex(4), allocatable :: Spectral_funct_of_Response_pos_r2c(:), cauchy_real_r2c(:), cauchy_imag_r2c(:)

  !! Auxiliary variables for fast complex convolution 
  integer(8) :: fft_plan, ifft_plan
  complex(fftw_real), allocatable :: inout_complex_1(:), inout_complex_2(:)
  complex(4), allocatable :: Spectral_funct_of_Response_pos_fft(:), cauchy_complex_fft(:)
  complex(4), allocatable :: Cauchy_complex(:), Response_complex(:)

  if(iv>0) write(6,*)'Cauchy_driver: enter ';

  allocate(Spectral_funct_of_Response_pos(0:nff_a-1));
  allocate(Response_pos(-nff_chi+1:nff_chi-1));
  allocate(Response_pos_fast(-nff_chi+1:nff_chi-1));
  allocate(Response_pos_fast_cmplx(-nff_chi+1:nff_chi-1));
  allocate(Response_pos_slow(-nff_chi+1:nff_chi-1));

  call random_number(Spectral_funct_of_Response_pos)

  !! Slow way to compute the Cauchy transform 
  !! (for positive and negative freq with a positive-support spectral function)
  Response_pos = 0
  do f=-nff_chi+1,nff_chi-1
    do s=0,nff_a-1;
      Response_pos(f) = Response_pos(f) + Spectral_funct_of_Response_pos(s)/cmplx((f-s)*d_omega, eps); 
    enddo;
  enddo;

  !! The Cauchy transform done with convolution
  allocate(Cauchy_real(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(Cauchy_imag(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(Response_real(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(Response_imag(-nff_chi+1-nff_a+1:nff_chi-1));
  cauchy_imag = 0; cauchy_real = 0;
  do s=-nff_chi+1-nff_a+1,nff_chi-1
    cauchy_imag(s) = -eps/((s*d_omega)**2+eps**2);
    cauchy_real(s) = (s*d_omega)/((s*d_omega)**2+eps**2);
  end do
  nff_cauchy = 2*nff_chi+nff_a
  call Convolution_slow(nff_a, Spectral_funct_of_Response_pos, size(cauchy_imag), cauchy_imag, size(Response_imag), Response_imag)
  call Convolution_slow(nff_a, Spectral_funct_of_Response_pos, size(cauchy_real), cauchy_real, size(Response_real), Response_real)
  Response_pos_slow(-nff_chi+1:nff_chi-1) = cmplx(Response_real(-nff_chi+1:nff_chi-1), Response_imag(-nff_chi+1:nff_chi-1))
  deallocate(Cauchy_real)
  deallocate(Cauchy_imag)
  deallocate(Response_real)
  deallocate(Response_imag)
  

  !! Fast way to compute the Cauchy transform : real 2 complex FFTW
  allocate(Cauchy_real(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(Cauchy_imag(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(Response_real(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(Response_imag(-nff_chi+1-nff_a+1:nff_chi-1));
  cauchy_imag = 0; cauchy_real = 0;
  do s=-nff_chi+1-nff_a+1,nff_chi-1
    cauchy_imag(s) = -eps/((s*d_omega)**2+eps**2);
    cauchy_real(s) = (s*d_omega)/((s*d_omega)**2+eps**2);
  end do
  nff_cauchy = 2*nff_chi+nff_a

  call Convolution_calc_r2c_dimensions(nff_a, nff_cauchy, nff_complex, nff_real)
  write(6,*)'nff_a, nff_cauchy: ', nff_a, nff_cauchy
  write(6,*)'nff_complex, nff_real: ', nff_complex, nff_real

  allocate(inout_complex(nff_complex))
  allocate(inout_real(nff_real))
  allocate(Spectral_funct_of_Response_pos_r2c(nff_complex))
  allocate(cauchy_real_r2c(nff_complex))
  allocate(cauchy_imag_r2c(nff_complex))

  call Convolution_create_r2c_plans(nff_real, c2r_plan, r2c_plan, inout_complex, inout_real);

  call Convolution_apply_r2c_4(size(cauchy_real), cauchy_real, cauchy_real_r2c, &
    nff_real, r2c_plan, inout_complex, inout_real);

  call Convolution_apply_r2c_4(size(cauchy_imag), cauchy_imag, cauchy_imag_r2c, &
    nff_real, r2c_plan, inout_complex, inout_real);

  call Convolution_apply_r2c_4(size(Spectral_funct_of_Response_pos), Spectral_funct_of_Response_pos, &
    Spectral_funct_of_Response_pos_r2c, nff_real, r2c_plan, inout_complex, inout_real);

  call Convolution_multiply_and_apply_c2r(nff_complex, Spectral_funct_of_Response_pos_r2c, cauchy_real_r2c, &
    size(Response_real), Response_real, c2r_plan, inout_complex, inout_real);

  call Convolution_multiply_and_apply_c2r(nff_complex, Spectral_funct_of_Response_pos_r2c, cauchy_imag_r2c, &
    size(Response_imag), Response_imag, c2r_plan, inout_complex, inout_real);

  Response_pos_fast(-nff_chi+1:nff_chi-1) = cmplx(Response_real(-nff_chi+1:nff_chi-1), Response_imag(-nff_chi+1:nff_chi-1))

  deallocate(Cauchy_real);
  deallocate(Cauchy_imag);
  deallocate(cauchy_real_r2c)
  deallocate(cauchy_imag_r2c)
  deallocate(Response_real);
  deallocate(Response_imag);
  deallocate(inout_complex)
  deallocate(inout_real)
  deallocate(Spectral_funct_of_Response_pos_r2c)

  !! 
  !! Fast way to compute the Cauchy transform real 2 complex for spectral function and complex 2 complex for "denominator"
  allocate(cauchy_complex(-nff_chi+1-nff_a+1:nff_chi-1));
  allocate(response_complex(-nff_chi+1-nff_a+1:nff_chi-1));
  cauchy_complex = 0;
  do f=-nff_chi+1-nff_a+1,nff_chi-1
    cauchy_complex(f) = 1.0/cmplx(f*d_omega, eps);
  end do
  nff_cauchy = 2*nff_chi+nff_a

  call Convolution_calc_c2c_dimensions(nff_a, nff_cauchy, nff_complex)
  nff_real = nff_complex
  allocate(inout_complex_1(nff_complex))
  allocate(inout_complex_2(nff_complex))
  allocate(inout_complex(nff_complex/2+1))
  allocate(inout_real(nff_real))
  allocate(cauchy_complex_fft(nff_complex))
  allocate(Spectral_funct_of_Response_pos_r2c(nff_complex/2+1))
  allocate(Spectral_funct_of_Response_pos_fft(nff_complex))

  call Convolution_create_r2c_plans(nff_complex, c2r_plan, r2c_plan, inout_complex, inout_real);
  call Convolution_create_c2c_plans(nff_complex, ifft_plan, fft_plan, inout_complex_1, inout_complex_2);

  call Convolution_apply_fft_4(size(cauchy_complex), cauchy_complex, cauchy_complex_fft, &
    nff_complex, fft_plan, inout_complex_1, inout_complex_2);

  call Convolution_apply_r2c_4(size(Spectral_funct_of_Response_pos), Spectral_funct_of_Response_pos(0:), &
    Spectral_funct_of_Response_pos_r2c, nff_complex, r2c_plan, inout_complex, inout_real);
  call Convolution_transform_r2c_to_c2c_4(Spectral_funct_of_Response_pos_r2c, nff_complex, Spectral_funct_of_Response_pos_fft);

  call Convolution_multiply_and_apply_ifft_4(nff_complex, Spectral_funct_of_Response_pos_fft, cauchy_complex_fft, &
    size(response_complex), response_complex, ifft_plan, inout_complex_1, inout_complex_2);

  Response_pos_fast_cmplx(-nff_chi+1:nff_chi-1) = response_complex(-nff_chi+1:nff_chi-1)


  !! Compare results  
  write(6,*) 'real'
  do f=-nff_chi+1,nff_chi-1
    write(6,*) f, real(Response_pos(f)), real(Response_pos_slow(f)), real(Response_pos_fast(f)), real(Response_pos_fast_cmplx(f))
  end do

  write(6,*) 'imag'
  do f=-nff_chi+1,nff_chi-1
    write(6,*) f, aimag(Response_pos(f)), aimag(Response_pos_slow(f)), aimag(Response_pos_fast(f)),aimag(Response_pos_fast_cmplx(f))
  end do

  deallocate(Spectral_funct_of_Response_pos);
  deallocate(Response_pos);
  deallocate(Response_pos_fast);
  deallocate(Response_pos_slow);

  if(iv>0) write(6,*)'Cauchy_driver: exit ';

end subroutine ! Cauchy_driver

!
! Show how to use Convolutions
!
subroutine Convolution_driver(iv)
  implicit none
  integer, intent(in) :: iv
  ! internal
  real(4), allocatable :: input_1(:), input_2(:), output_slow(:), output_fast(:)
  integer, parameter :: start=0
  integer, parameter :: nff_in=8
  integer, parameter :: nff_out=16;
  integer :: f

  !! Auxiliary variables for fast convolutions
  integer :: nff_complex, nff_real
  integer(8) :: r2c_plan, c2r_plan
  complex(fftw_real), allocatable :: inout_complex(:)
  real(fftw_real), allocatable :: inout_real(:)
  complex(4), allocatable :: input_1_r2c(:), input_2_r2c(:)
  

  if(iv>0) write(6,*)'Convolution_driver: enter ';
  allocate(input_1(start:start+nff_in-1));
  allocate(input_2(start:start+nff_in-1));
  allocate(output_slow(start:start+nff_out-1));
  allocate(output_fast(start:start+nff_out-1));
  
  call random_number(input_1)
  call random_number(input_2)

  !! Slow way
  call Convolution_slow(nff_in, input_1(start:), nff_in, input_2(start:),  nff_out, output_slow(start:))

  !! Fast way
  call Convolution_calc_r2c_dimensions(nff_in, nff_in, nff_complex, nff_real)
  allocate(inout_complex(nff_complex))
  allocate(inout_real(nff_real))
  allocate(input_1_r2c(nff_complex))
  allocate(input_2_r2c(nff_complex))


  call Convolution_create_r2c_plans(nff_real, c2r_plan, r2c_plan, inout_complex, inout_real);
  call Convolution_apply_r2c_4(nff_in, input_1(start:), input_1_r2c, nff_real, r2c_plan, inout_complex, inout_real);
  call Convolution_apply_r2c_4(nff_in, input_2(start:), input_2_r2c, nff_real, r2c_plan, inout_complex, inout_real);
  call Convolution_multiply_and_apply_c2r(nff_complex, input_1_r2c, input_2_r2c, nff_out, output_fast(start:), &
    c2r_plan, inout_complex, inout_real);

  deallocate(inout_complex)
  deallocate(inout_real)
  deallocate(input_1_r2c)
  deallocate(input_2_r2c)

  !! Compare results  
  do f=start,start+nff_out-1
    write(6,*) f, output_slow(f), output_fast(f)
  end do


  deallocate(input_1)
  deallocate(input_2)
  deallocate(output_slow)
  deallocate(output_fast)
  if(iv>0) write(6,*)'Convolution_driver: exit ';
  stop 'Convolution_driver: exit'
end subroutine !Convolution_driver

!
!
!
subroutine Convolution_slow(nff_a,a, nff_b,b,  nff_c,c)
! 
  integer, intent(in) :: nff_a, nff_b, nff_c
  real(4), intent(in) :: a(:), b(:)
  real(4), intent(out) :: c(:)
  
  ! internal
  integer :: s1, s2

  c(1:nff_c) = 0
  do s1=1,nff_a; 
    do s2=1,nff_b;
      if(s1+s2-1<1 .or. s1+s2-1>nff_c) cycle
      c(s1+s2-1) = c(s1+s2-1) + a(s1) * b(s2); 
    enddo;
  enddo;
 
end subroutine ! Convolution_slow

!
!
!
subroutine Convolution_slow_8(nff_a,a, nff_b,b,  nff_c,c)
! 
  integer, intent(in) :: nff_a, nff_b, nff_c
  real(8), intent(in) :: a(:), b(:)
  real(8), intent(out) :: c(:)
  
  ! internal
  integer :: s1, s2

  c(1:nff_c) = 0
  do s1=1,nff_a; 
    do s2=1,nff_b;
      if(s1+s2-1<1 .or. s1+s2-1>nff_c) cycle
      c(s1+s2-1) = c(s1+s2-1) + a(s1) * b(s2); 
    enddo;
  enddo;
 
end subroutine ! Convolution_slow_8

!
!
!
subroutine Convolution_slow_zdz(nff_a,a, nff_b,b,  nff_c,c)
! 
  integer, intent(in) :: nff_a, nff_b, nff_c
  complex(8),  intent(in) :: a(:)
  real(8),     intent(in) :: b(:)
  complex(8), intent(out) :: c(:)
  
  ! internal
  integer :: s1, s2

  c(1:nff_c) = 0
  do s1=1,nff_a; 
    do s2=1,nff_b;
      if(s1+s2-1<1 .or. s1+s2-1>nff_c) cycle
      c(s1+s2-1) = c(s1+s2-1) + a(s1) * b(s2); 
    enddo;
  enddo;
 
end subroutine ! Convolution_slow_zdz

!
!
!
subroutine Convolution_slow_zzz(nff_a,a, nff_b,b,  nff_c,c)
! 
  integer, intent(in) :: nff_a, nff_b, nff_c
  complex(8),  intent(in) :: a(:)
  complex(8),     intent(in) :: b(:)
  complex(8), intent(out) :: c(:)
  
  ! internal
  integer :: s1, s2

  c(1:nff_c) = 0
  do s1=1,nff_a; 
    do s2=1,nff_b;
      if(s1+s2-1<1 .or. s1+s2-1>nff_c) cycle
      c(s1+s2-1) = c(s1+s2-1) + a(s1) * b(s2); 
    enddo;
  enddo;
 
end subroutine ! Convolution_slow_zdz

!
!
!
subroutine Convolution_calc_r2c_dimensions(nff_in_1, nff_in_2, nff_complex, nff_real)
  implicit none
  integer, intent(in)  :: nff_in_1, nff_in_2
  integer, intent(out) :: nff_complex, nff_real
  nff_real = (nff_in_1+nff_in_2)
  nff_complex  = nff_real/2+1;
end subroutine ! nff_in_out_2_nff_c2r_conv


!
!
!
subroutine Convolution_calc_c2c_dimensions(nff_in_1, nff_in_2, nfft_complex)
  implicit none
  integer, intent(in)  :: nff_in_1, nff_in_2
  integer, intent(out) :: nfft_complex
  nfft_complex = (nff_in_1+nff_in_2)
end subroutine ! Convolution_calc_c2c_dimensions

!
!
!
subroutine Convolution_calc_dimensions(nff_in_1, nff_in_2, nfft_c2c, nff_r2c)
  implicit none
  integer, intent(in)  :: nff_in_1, nff_in_2
  integer, intent(out) :: nfft_c2c
  integer, intent(out), optional :: nff_r2c
  
  nfft_c2c = (nff_in_1+nff_in_2)
  if(present(nff_r2c)) nff_r2c = nfft_c2c/2+1;
end subroutine ! Convolution_calc_c2c_dimensions

!
! Computes 6 numbers out of limits 
!
subroutine Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc_aux,falloc_aux)
  implicit none
  ! external
  integer, intent(in) :: s1, f1, s2, f2, s3, f3
  integer, intent(out) :: n1,n2,ncaux,saux,salloc_aux,falloc_aux
  
  n1 = f1-s1+1
  n2 = f2-s2+1
  saux = s1+s2
  salloc_aux = min(saux,s3)
  falloc_aux = f3
  ncaux = min(falloc_aux-saux+1, f1+f2-saux+1)
  
end subroutine ! Convolution_limits

!!
!! Example of fast loops with complex input arrays and c2c fftw
!!
subroutine Convolution_example_cmplx_cmplx(iv)
  use m_precision, only : fftw_real

  implicit none
  integer, intent(in) :: iv

  !! internal
  real(8), allocatable :: ar(:), ai(:), br(:), bi(:)
  complex(8), allocatable :: a(:), b(:), c(:), caux(:), cf(:)
  integer :: i,j,k,s1,f1,n1, s2,f2,n2,s3,f3,ncaux,saux
  integer :: salloc,falloc
  
  !! For fast convolution
  integer    :: nff_conv, nff_conv_r2c
  integer(8) :: fft_plan, ifft_plan
  complex(fftw_real), allocatable :: aux_c2c_1(:)
  complex(fftw_real), allocatable :: aux_c2c_2(:)

  complex(fftw_real), allocatable :: a_c2c(:)
  complex(fftw_real), allocatable :: b_c2c(:)
  !! END of For fast convolution
  
  s1 = -10;  f1 = -8;
  s2 = -4;   f2 = 10;
  s3 = -6;   f3 = 2;

  allocate(a(s1:f1))
  allocate(b(s2:f2))
  allocate(c(s3:f3))

  allocate(ai(s1:f1))
  allocate(ar(s1:f1))
  call random_number(ai)
  call random_number(ar)
  a = cmplx(ar,ai,8)

  allocate(bi(s2:f2))
  allocate(br(s2:f2))
  call random_number(bi)
  call random_number(br)
  b = cmplx(br,bi,8)

  !! Double loop which must be accelerated with FFT(W)
  c = 0
  do i=s1,f1
  do j=s2,f2
    k = i+j
    if(k<s3 .or. k>f3) cycle
    c(k) = c(k) + a(i)*b(j)
  enddo
  enddo
  !! END of Double loop which must be accelerated with FFT(W)

  
  !! Doing the same calculation but using "standart" convolution
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(caux(salloc:falloc))
  caux = 0
  call Convolution_slow_zzz(n1,a(s1:), n2,b(s2:),ncaux,caux(saux:))

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux
  if(iv>0)write(6,*) 'slow'
  write(6,'(i5,4g25.16)') &
    (k, real(c(k)), real(caux(k)), aimag(c(k)), aimag(caux(k)),k=s3,f3);
  !! END of Doing the same calculation but using "standart" convolution
    
  !! Fast way of doind the same calculation
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(cf(salloc:falloc))
  cf = 0
  
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

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux, nff_conv'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux, nff_conv
  write(6,*) 'fast'
  write(6,'(i5,4g25.16)') &
    (k, real(c(k)), real(cf(k)), aimag(c(k)), aimag(cf(k)), k=s3,f3);
  !! END of Fast way of doind the same calculation

  read(5,*)
  
end subroutine ! example_cmplx_cmplx


!!
!! Example of fast loops with real input arrays and r2c fftw
!!
subroutine Convolution_example_cmplx_real(iv)
  use m_precision, only : fftw_real

  implicit none
  integer, intent(in) :: iv

  !! internal
  real(8), allocatable :: b(:), ar(:), ai(:)
  complex(8), allocatable :: a(:), c(:), caux(:), cf(:)
  integer :: i,j,k,s1,f1,n1, s2,f2,n2,s3,f3,ncaux,saux
  integer :: salloc,falloc
  
  !! For fast convolution
  integer :: nff_conv_r2c, nff_conv
  integer(8) ::  c2r_plan, r2c_plan, fft_plan, ifft_plan
  complex(fftw_real), allocatable :: aux_c2c_1(:)
  complex(fftw_real), allocatable :: aux_c2c_2(:)
  complex(fftw_real), allocatable :: aux_r2c_1(:)
  real(fftw_real),    allocatable :: aux_r2c_2(:)

  complex(fftw_real), allocatable :: a_c2c(:)
  complex(fftw_real), allocatable :: b_c2c(:)
  complex(fftw_real), allocatable :: b_r2c(:)

  !! END of For fast convolution
  
  s1 = -10;   f1 = -8
  s2 = -4;   f2 = 10;
  s3 = -2;   f3 = 2;
  
  allocate(a(s1:f1))
  allocate(ai(s1:f1))
  allocate(ar(s1:f1))

  allocate(b(s2:f2))
  allocate(c(s3:f3))
  call random_number(ai)
  call random_number(ar)
  a = cmplx(ar,ai,8)
  call random_number(b)

  !! Double loop which must be accelerated with FFT(W)
  c = 0
  do i=s1,f1
  do j=s2,f2
    k = i+j
    if(k<s3 .or. k>f3) cycle
    c(k) = c(k) + a(i)*b(j)
  enddo
  enddo
  !! END of Double loop which must be accelerated with FFT(W)

  
  !! Doing the same calculation but using "standart" convolution
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(caux(salloc:falloc))
  caux = 0
  call Convolution_slow_zdz(n1,a(s1:), n2,b(s2:),ncaux,caux(saux:))

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux
  if(iv>0)write(6,*) 'slow'
  write(6,'(i5,4g25.16)') (k, real(c(k)), real(caux(k)), aimag(c(k)), aimag(caux(k)), k=s3,f3)
  !! END of Doing the same calculation but using "standart" convolution
    
  !! Fast way of doind the same calculation
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(cf(salloc:falloc))
  cf = 0
  
  !! Take care of real to complex transform(s)
  call Convolution_calc_dimensions(n1,n2,nff_conv,nff_conv_r2c)
  allocate(aux_r2c_1(nff_conv_r2c))
  allocate(aux_r2c_2(nff_conv))
  allocate(b_r2c(nff_conv_r2c))
  allocate(b_c2c(nff_conv))
  call Convolution_create_r2c_plans(nff_conv, c2r_plan, r2c_plan, aux_r2c_1, aux_r2c_2);
  !! END of Take care of real to complex transform

  !! Take care of complex to complex transform(s)
  allocate(aux_c2c_1(nff_conv))
  allocate(aux_c2c_2(nff_conv))
  allocate(a_c2c(nff_conv))
  call Convolution_create_c2c_plans(nff_conv, ifft_plan, fft_plan, aux_c2c_1, aux_c2c_2);
  !! END of Take care of complex to complex transform(s)

  call Convolution_apply_fft_8(n1, a(s1:), a_c2c, nff_conv, fft_plan, aux_c2c_1, aux_c2c_2);
  call Convolution_apply_r2c_8(n2, b(s2:), b_r2c, nff_conv, r2c_plan, aux_r2c_1, aux_r2c_2);
  call Convolution_transform_r2c_to_c2c_8(b_r2c, nff_conv, b_c2c);
  call Convolution_multiply_and_apply_ifft_8(nff_conv, a_c2c, b_c2c, ncaux, cf(saux:), &
    ifft_plan, aux_c2c_1, aux_c2c_2);

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux, nff_conv, nff_conv_r2c
  write(6,*) 'fast'
  write(6,'(i5,4g25.16)') (k, real(c(k)), real(cf(k)), aimag(c(k)), aimag(cf(k)), k=s3,f3)
  !! END of Fast way of doind the same calculation

end subroutine ! example_cmplx_real


!!
!! Example of fast loops with real input arrays and r2c fftw
!!
subroutine Convolution_example_real_real(iv)
  use m_precision, only : fftw_real

  implicit none
  integer, intent(in) :: iv

  !! internal
  real(8), allocatable :: a(:), b(:), c(:), caux(:), cf(:)
  integer :: i,j,k,s1,f1,n1, s2,f2,n2, s3,f3, ncaux,saux
  integer :: salloc,falloc
  
  !! For fast convolution
  integer :: nff_complex, nff_real
  integer(8) ::  c2r_plan, r2c_plan
  real(fftw_real), allocatable :: inout_real(:)
  complex(fftw_real), allocatable :: inout_complex(:)
  complex(fftw_real), allocatable :: a_r2c(:)
  complex(fftw_real), allocatable :: b_r2c(:)
  !! END of For fast convolution
  
  s1 = 7;   f1 = 8
  s2 = 0;   f2 = 10;
  s3 = 0;   f3 = 20;
  
  allocate(a(s1:f1))
  allocate(b(s2:f2))
  allocate(c(s3:f3))
  call random_number(a)
  call random_number(b)

  !! Double loop which must be accelerated with FFT(W)
  c = 0
  do i=s1,f1
  do j=s2,f2
    k = i+j
    if(k<s3 .or. k>f3) cycle
    c(k) = c(k) + a(i)*b(j)
  enddo
  enddo
  !! END of Double loop which must be accelerated with FFT(W)

  
  !! Doing the same calculation but using "standart" convolution
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(caux(salloc:falloc))
  caux = 0
  call Convolution_slow_8(n1,a(s1:), n2,b(s2:),ncaux,caux(saux:))

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux
  write(6,*) 'slow'
  do k=s3,f3
    write(6,*) k, c(k), caux(k)
  enddo
  !! END of Doing the same calculation but using "standart" convolution
  
  !! Fast way of doind the same calculation
  call Convolution_limits(s1,f1,s2,f2,s3,f3, n1,n2,ncaux,saux,salloc,falloc);
  allocate(cf(salloc:falloc))
  cf = 0
  call Convolution_calc_r2c_dimensions(n1,n2,nff_complex, nff_real)
  allocate(inout_complex(nff_complex))
  allocate(inout_real(nff_real))
  allocate(a_r2c(nff_complex))
  allocate(b_r2c(nff_complex))

  call Convolution_create_r2c_plans(nff_real, c2r_plan, r2c_plan, inout_complex, inout_real);
  call Convolution_apply_r2c_8(n1, a(s1:), a_r2c, nff_real, r2c_plan, inout_complex, inout_real);
  call Convolution_apply_r2c_8(n2, b(s2:), b_r2c, nff_real, r2c_plan, inout_complex, inout_real);
  call Convolution_multiply_and_apply_c2r_8(nff_complex, a_r2c, b_r2c, ncaux, cf(saux:), &
    c2r_plan, inout_complex, inout_real);

  if(iv>0)write(6,*) 'salloc, falloc, saux, ncaux'
  if(iv>0)write(6,*)  salloc, falloc, saux, ncaux
  write(6,*) 'fast'
  do k=s3,f3
    write(6,*) k, c(k), cf(k)
  enddo
  !! END of Fast way of doind the same calculation

end subroutine ! example_real_real

!
!
!
subroutine Convolution_create_r2c_plans(nff_real, c2r_plan, r2c_plan, input_or_output_complex, input_or_output_real)
  implicit none
#include <fftw3.f>
  integer, intent(in) :: nff_real
  integer(8), intent(inout) :: c2r_plan, r2c_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)

  call dfftw_plan_dft_c2r_1d(c2r_plan, nff_real, input_or_output_complex, input_or_output_real,    FFTW_ESTIMATE);
  call dfftw_plan_dft_r2c_1d(r2c_plan, nff_real, input_or_output_real,    input_or_output_complex, FFTW_ESTIMATE);

end subroutine ! 


!
!
!
subroutine Convolution_create_c2c_plans(nfft, c2c_plan_bwd, c2c_plan_fwd, input_or_output_1, input_or_output_2)
  implicit none
#include <fftw3.f>
  integer, intent(in) :: nfft
  integer(8), intent(inout) :: c2c_plan_bwd, c2c_plan_fwd
  complex(fftw_real), intent(inout) :: input_or_output_1(:)
  complex(fftw_real), intent(inout) :: input_or_output_2(:)

  call dfftw_plan_dft_1d(c2c_plan_bwd, nfft, input_or_output_1, input_or_output_2, FFTW_BACKWARD, FFTW_ESTIMATE);
  call dfftw_plan_dft_1d(c2c_plan_fwd, nfft, input_or_output_2, input_or_output_1, FFTW_FORWARD,  FFTW_ESTIMATE);

end subroutine ! 

!
!
!
subroutine Convolution_apply_r2c_8(nff_a, a, a_r2c, nff_real, r2c_plan, input_or_output_complex, input_or_output_real)
  integer, intent(in) :: nff_a
  real(8), intent(in) :: a(:)
  complex(8), intent(inout) :: a_r2c(:)
  ! auxiliary
  integer, intent(in)       :: nff_real
  integer(8), intent(inout) :: r2c_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)
    
  input_or_output_real(1:nff_a) = a(1:nff_a); 
  input_or_output_real(nff_a+1:nff_real) = 0;
  call dfftw_execute(r2c_plan);
  a_r2c(1:nff_real/2+1) = input_or_output_complex(1:nff_real/2+1)/sqrt(real(1,fftw_real)*nff_real)

end subroutine !Convolution_apply_r2c_8

!
!
!
subroutine Convolution_apply_r2c_4(nff_a, a, a_r2c, nff_real, r2c_plan, input_or_output_complex, input_or_output_real)
  integer, intent(in) :: nff_a
  real(4), intent(in) :: a(:)
  complex(4), intent(inout) :: a_r2c(:)
  ! auxiliary
  integer, intent(in)       :: nff_real
  integer(8), intent(inout) :: r2c_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)
    
  input_or_output_real(1:nff_a) = a(1:nff_a); input_or_output_real(nff_a+1:nff_real) = 0;
  call dfftw_execute(r2c_plan);
  a_r2c(1:nff_real/2+1) = input_or_output_complex(1:nff_real/2+1)/sqrt(real(1,fftw_real)*nff_real)

end subroutine !Convolution_apply_r2c_4

!
!
!
subroutine Convolution_apply_fft_8(nff_a, a, a_fft, nfft, c2c_plan_fwd, input_or_output_1, input_or_output_2)
  integer, intent(in) :: nff_a
  complex(8), intent(in) :: a(:)
  complex(8), intent(inout) :: a_fft(:)
  ! auxiliary
  integer, intent(in)       :: nfft
  integer(8), intent(inout) :: c2c_plan_fwd
  complex(fftw_real), intent(inout) :: input_or_output_1(:), input_or_output_2(:)
    
  input_or_output_2(1:nff_a) = a(1:nff_a); input_or_output_2(nff_a+1:nfft) = 0;
  call dfftw_execute(c2c_plan_fwd);
  a_fft(1:nfft) = input_or_output_1(1:nfft)/sqrt(real(1,fftw_real)*nfft)

end subroutine ! Convolution_apply_fft_8

!
!
!
subroutine Convolution_apply_fft_4(nff_a, a, a_fft, nfft, c2c_plan_fwd, input_or_output_1, input_or_output_2)
  integer, intent(in) :: nff_a
  complex(4), intent(in) :: a(:)
  complex(4), intent(inout) :: a_fft(:)
  ! auxiliary
  integer, intent(in)       :: nfft
  integer(8), intent(inout) :: c2c_plan_fwd
  complex(fftw_real), intent(inout) :: input_or_output_1(:), input_or_output_2(:)
    
  input_or_output_2(1:nff_a) = a(1:nff_a); input_or_output_2(nff_a+1:nfft) = 0;
  call dfftw_execute(c2c_plan_fwd);
  a_fft(1:nfft) = input_or_output_1(1:nfft)/sqrt(real(1,fftw_real)*nfft)

end subroutine ! Convolution_apply_fft_4

!
!
!
subroutine Convolution_multiply_and_apply_c2r(nff_complex, a_r2c, b_r2c, nff_c, c, &
  c2r_plan, input_or_output_complex, input_or_output_real)
!
  implicit none
  integer, intent(in) :: nff_complex, nff_c
  complex(4), intent(in) :: a_r2c(:), b_r2c(:)
  real(4), intent(out) :: c(:)
  ! auxiliary
  integer(8), intent(inout) :: c2r_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)

  input_or_output_complex(1:nff_complex) = a_r2c(1:nff_complex) * b_r2c(1:nff_complex);
  call dfftw_execute(c2r_plan);
  c(1:nff_c) = input_or_output_real(1:nff_c)

end subroutine ! s_Convolution_ss

!
!
!
subroutine Convolution_multiply_and_apply_c2r_8(nff_complex, a_r2c, b_r2c, nff_c, c, &
  c2r_plan, input_or_output_complex, input_or_output_real)
!
  implicit none
  integer, intent(in) :: nff_complex, nff_c
  complex(8), intent(in) :: a_r2c(:), b_r2c(:)
  real(8), intent(out) :: c(:)
  ! auxiliary
  integer(8), intent(inout) :: c2r_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)

  input_or_output_complex(1:nff_complex) = a_r2c(1:nff_complex) * b_r2c(1:nff_complex);
  call dfftw_execute(c2r_plan);
  c(1:nff_c) = input_or_output_real(1:nff_c)

end subroutine ! Convolution_multiply_and_apply_c2r_8


!
!
!
subroutine Convolution_apply_c2r_4(nff_complex, ab_r2c, nff_c, c, &
  c2r_plan, input_or_output_complex, input_or_output_real)
!
  implicit none
  integer, intent(in) :: nff_complex, nff_c
  complex(4), intent(in) :: ab_r2c(:)
  real(4), intent(out) :: c(:)
  ! auxiliary
  integer(8), intent(inout) :: c2r_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)

  input_or_output_complex(1:nff_complex) = ab_r2c(1:nff_complex);
  call dfftw_execute(c2r_plan);
  c(1:nff_c) = input_or_output_real(1:nff_c)

end subroutine ! s_Convolution_ss

!
!
!
subroutine Convolution_apply_c2r_8(nff_complex, ab_r2c, nff_c, c, &
  c2r_plan, input_or_output_complex, input_or_output_real)
!
  implicit none
  integer, intent(in) :: nff_complex, nff_c
  complex(8), intent(in) :: ab_r2c(:)
  real(8), intent(out) :: c(:)
  ! auxiliary
  integer(8), intent(inout) :: c2r_plan
  complex(fftw_real), intent(inout) :: input_or_output_complex(:)
  real(fftw_real), intent(inout)    :: input_or_output_real(:)

  input_or_output_complex(1:nff_complex) = ab_r2c(1:nff_complex);
  call dfftw_execute(c2r_plan);
  c(1:nff_c) = input_or_output_real(1:nff_c)

end subroutine ! s_Convolution_ss

!
!
!
subroutine Convolution_multiply_and_apply_ifft_8(nfft, a_fft, b_fft, nff_c, c, &
  c2c_plan_bwd, input_or_output_1, input_or_output_2)
!
  implicit none
  integer, intent(in) :: nfft, nff_c
  complex(8), intent(in) :: a_fft(:), b_fft(:)
  complex(8), intent(out) :: c(:)
  ! auxiliary
  integer(8), intent(inout) :: c2c_plan_bwd
  complex(fftw_real), intent(inout) :: input_or_output_1(:), input_or_output_2(:)

  input_or_output_1(1:nfft) = a_fft(1:nfft) * b_fft(1:nfft);
  call dfftw_execute(c2c_plan_bwd);
  c(1:nff_c) = input_or_output_2(1:nff_c)

end subroutine ! Convolution_multiply_and_apply_ifft_8

!
!
!
subroutine Convolution_multiply_and_apply_ifft_4(nfft, a_fft, b_fft, nff_c, c, &
  c2c_plan_bwd, input_or_output_1, input_or_output_2)
!
  implicit none
  integer, intent(in) :: nfft, nff_c
  complex(4), intent(in) :: a_fft(:), b_fft(:)
  complex(4), intent(out) :: c(:)
  ! auxiliary
  integer(8), intent(inout) :: c2c_plan_bwd
  complex(fftw_real), intent(inout) :: input_or_output_1(:), input_or_output_2(:)

  input_or_output_1(1:nfft) = a_fft(1:nfft) * b_fft(1:nfft);
  call dfftw_execute(c2c_plan_bwd);
  c(1:nff_c) = input_or_output_2(1:nff_c)

end subroutine ! Convolution_multiply_and_apply_ifft_4


!
! Transforms the half-arrays from a real-to-complex FFT to full arrays
! which we would normally get from a complex-to-complex FFT.
!
subroutine Convolution_transform_r2c_to_c2c_4(r2c_output, nfft, c2c_output)
  implicit none
  complex(4), intent(in) :: r2c_output(:)
  integer, intent(in) :: nfft
  complex(4), intent(out) :: c2c_output(:)

  !! internal
  integer :: f

  c2c_output(1:nfft/2+1) = r2c_output(1:nfft/2+1);

  do f=nfft/2+2,nfft; 
    c2c_output(f) = conjg(c2c_output(nfft-f+2));
  enddo

end subroutine !! Convolution_transform_r2c_to_c2c_4

!
! Transforms the half-arrays from a real-to-complex FFT to full arrays
! which we would normally get from a complex-to-complex FFT.
!
subroutine Convolution_transform_r2c_to_c2c_8(r2c_output, nfft, c2c_output)
  implicit none
  complex(8), intent(in) :: r2c_output(:)
  integer, intent(in) :: nfft
  complex(8), intent(out) :: c2c_output(:)

  !! internal
  integer :: f

  c2c_output(1:nfft/2+1) = r2c_output(1:nfft/2+1);

  do f=nfft/2+2,nfft; 
    c2c_output(f) = conjg(c2c_output(nfft-f+2));
  enddo

end subroutine !! Convolution_transform_r2c_to_c2c_8


end module ! m_convolution
