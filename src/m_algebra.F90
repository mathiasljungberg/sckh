module m_algebra

  use m_precision, only : blas_int
  use m_log, only : die, warn, log_memory_note
#include "m_define_macro.F90"

  implicit none

  interface matinv
    module procedure matinv_c_simple
    module procedure matinv_c_blas
    module procedure matinv_d
    module procedure matinv_s
    module procedure matinv_z
  end interface !matinv 

  interface matinv_c
    module procedure matinv_c_simple
    module procedure matinv_c_blas
  end interface
    
  interface complex_eigen;       module procedure complex_eigen;       end interface
  interface complex_inverse;     module procedure complex_inverse;     end interface
  interface generaleigen;        module procedure generaleigen;        end interface 
  interface matinv_csy;          module procedure matinv_csy;          end interface
  interface matinv_s;            module procedure matinv_s;            end interface
  interface matinv_z;            module procedure matinv_z;            end interface
  interface real_inverse;        module procedure real_inverse;        end interface
  interface simplediagonalize;   module procedure simplediagonalize;   end interface
  interface simpleeigen;         module procedure simpleeigen;         end interface
  interface t11_inverse_sp;      module procedure t11_inverse_sp;      end interface

  interface SymmetricEigen;
    module procedure sSymmetricEigen;
    module procedure dSymmetricEigen;
    module procedure zHermitianEigen;
  end interface

  interface Diagonalize
    module procedure sSymmetricEigen;
    module procedure dSymmetricEigen;
    module procedure sSymmGeneralEigen !(H,S,E,X)
    module procedure dSymmGeneralEigen !(H,S,E,X)
    module procedure dSymmGeneralEigen_blas_like!(n,H,ldh,S,lds,E,X,ldx)
    module procedure zHermitianEigen; !   
    module procedure zHermitGenDiag !(A,B, -> E,X)
    module procedure zHermitGenDiag_evals !(A,B, ->E)
  end interface ! Diagonalize

  interface init_unit_matrix
    module procedure alloc_init_d_unit_matrix
  end interface ! alloc_init_unit_matrix

  interface matmul_AtBC_d
     module procedure matmul_AtBC_d_work
     module procedure matmul_AtBC_d_nowork
  end interface

  !! Variables for t11_inverse
  integer(blas_int) :: ndim_t11_inverse_sp=-999 
  complex(4), allocatable, private :: d(:), aa(:), b(:), c(:);
  !$OMP THREADPRIVATE(aa,b,c,d,ndim_t11_inverse_sp)

  !! Varables for matinv_z
  integer(blas_int) :: ndim_matinv_z=-999
  integer(blas_int), allocatable :: ipiv_matinv_z(:)
  complex(8), allocatable :: work_matinv_z(:)
  !$OMP THREADPRIVATE(ndim_matinv_z, ipiv_matinv_z, work_matinv_z)

  !! Varables for matinv_c
  integer(blas_int) :: ndim_matinv_c=-999
  integer(blas_int), allocatable :: ipiv_matinv_c(:)
  complex(4), allocatable :: work_matinv_c(:)
  !$OMP THREADPRIVATE(ndim_matinv_c, ipiv_matinv_c, work_matinv_c)

  !! Varables for matinv_d
  integer(blas_int) :: ndim_matinv_d=-999
  integer(blas_int), allocatable :: ipiv_matinv_d(:)
  real(8), allocatable :: work_matinv_d(:)
  !$OMP THREADPRIVATE(ndim_matinv_d, ipiv_matinv_d, work_matinv_d)

  !! Varables for matinv_s
  integer(blas_int) :: ndim_matinv_s=-999
  integer(blas_int), allocatable :: ipiv_matinv_s(:)
  real(4), allocatable :: work_matinv_s(:)
  !$OMP THREADPRIVATE(ndim_matinv_s, ipiv_matinv_s, work_matinv_s)

  !! Varables for matinv_csy
  integer(blas_int) :: ndim_matinv_csy=-999
  integer(blas_int), allocatable :: ipiv_matinv_csy(:)
  complex(4), allocatable :: work_matinv_csy(:)
  integer(blas_int) :: lwork_matinv_csy
  complex(4), allocatable :: worki_matinv_csy(:)
  !$OMP THREADPRIVATE(ndim_matinv_csy, ipiv_matinv_csy, work_matinv_csy)
  !$OMP THREADPRIVATE(lwork_matinv_csy, worki_matinv_csy)

  !! Varables for matinv_zsy
  integer(blas_int) :: ndim_matinv_zsy=-999
  integer(blas_int), allocatable :: ipiv_matinv_zsy(:)
  complex(8), allocatable :: work_matinv_zsy(:)
  integer(blas_int) :: lwork_matinv_zsy
  complex(8), allocatable :: worki_matinv_zsy(:)
  !$OMP THREADPRIVATE(ndim_matinv_zsy, ipiv_matinv_zsy, work_matinv_zsy)
  !$OMP THREADPRIVATE(lwork_matinv_zsy, worki_matinv_zsy)

#ifdef USE_OWN_CZDOT
  interface zdotc
    module procedure zdotc_z_z
    module procedure zdotc_d2_z
  end interface ! zdotc

  public cdotc, cdot, sdot, scnrm2, snrm2, zdotc, zdot
#endif

  contains

!
!
!
subroutine alloc_init_d_unit_matrix(n, m)
  implicit none
  integer, intent(in) :: n
  real(8), allocatable, intent(inout) :: m(:,:)
  !
  integer :: i
  _dealloc(m)
  allocate(m(n,n))
  m = 0
  do i=1,n; m(i,i) = 1; enddo 
  
end subroutine ! alloc_init_d_unit_matrix

!
!
!
function cross_product(a,b) result(c)
  implicit none
  !! external
  real(8), intent(in) :: a(3), b(3)
  real(8) :: c(3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
    
end function ! cross_product


#ifdef USE_OWN_CZDOT
!
!
!
complex(4) function cdotc(n,a,inca, b,incb)
  integer, intent(in) :: n, inca, incb
  complex(4), intent(in) :: a(*),b(*)
 
  integer :: ia, ib, i

  if(inca==1 .and. incb==1) then
    cdotc = sum(conjg(a(1:n))*b(1:n))
  else
    cdotc = 0
    ia = 1
    ib = 1
    do i=1,n
      cdotc = cdotc + conjg(a(ia))*b(ib)
      ia = ia + inca
      ib = ib + incb
    enddo
  endif
end function ! cdotc

!
!
!
complex(4) function cdot(n,a,inca, b,incb)
  integer, intent(in) :: n, inca, incb
  complex(4), intent(in) :: a(*),b(*)
 
  integer :: ia, ib, i

  if(inca==1 .and. incb==1) then
    cdot = sum(a(1:n)*b(1:n))
  else
    cdot = 0
    ia = 1
    ib = 1
    do i=1,n
      cdot = cdot + a(ia)*b(ib)
      ia = ia + inca
      ib = ib + incb
    enddo
  endif
end function ! cdot


!
!
!
real(4) function sdot(n,a,inca, b,incb)
  integer, intent(in) :: n, inca, incb
  real(4), intent(in) :: a(*),b(*)

  integer :: ia, ib, i

  if(inca==1 .and. incb==1) then
    sdot = sum(a(1:n)*b(1:n))
  else
    sdot = 0
    ia = 1
    ib = 1
    do i=1,n
      sdot = sdot + a(ia)*b(ib)
      ia = ia + inca
      ib = ib + incb
    enddo
  endif
end function ! sdot

 
!
!
!
REAL(4) FUNCTION SNRM2(N,X,INCX)
  INTEGER, intent(in) :: INCX, N
  REAL(4), intent(in) :: X(*)
!*  Purpose
!*  =======
!*  SNRM2 returns the euclidean norm of a vector via the function
!*  name, so that
!*     SNRM2 := sqrt( x'*x ).
!*  Further Details
!*  ===============
!*  -- This version written on 25-October-1982.
!*     Modified on 14-October-1993 to inline the call to SLASSQ.
!*     Sven Hammarling, Nag Ltd.

  REAL(4) :: ONE,ZERO
  PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
  REAL(4) :: NORM,SCALE,SSQ
  INTRINSIC ABS,SQRT
  IF (N.LT.1 .OR. INCX.LT.1) THEN
    NORM = ZERO
  ELSE IF (N.EQ.1) THEN
    NORM = ABS(X(1))
  ELSE
    SCALE = ZERO
    SSQ = ONE
    CALL SLASSQ( N, X, INCX, SCALE, SSQ )
    NORM = SCALE*SQRT(SSQ)
  END IF

  SNRM2 = NORM

END FUNCTION ! SNRM2

!
!
!
REAL(4) FUNCTION SCNRM2(N,X,INCX)
  INTEGER, intent(in) :: INCX,N
  COMPLEX(4), intent(in) :: X(*)
!*  Purpose
!*  SCNRM2 returns the euclidean norm of a vector via the function
!*  name, so that
!*     SCNRM2 := sqrt( x**H*x )
!*  Further Details
!*  -- This version written on 25-October-1982.
!*     Modified on 14-October-1993 to inline the call to CLASSQ.
!*     Sven Hammarling, Nag Ltd.
  REAL(4) :: ONE,ZERO
  PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
  REAL(4) :: NORM,SCALE,SSQ
  INTRINSIC ABS,AIMAG,REAL,SQRT
  IF (N.LT.1 .OR. INCX.LT.1) THEN
    NORM = ZERO
  ELSE
    SCALE = ZERO
    SSQ = ONE
!*        The following loop is equivalent to this call to the LAPACK
!*        auxiliary routine:
    CALL CLASSQ( N, X, INCX, SCALE, SSQ )
    NORM = SCALE*SQRT(SSQ)
  END IF
  SCNRM2 = NORM
end function ! scnrm2

!
!
!
complex(8) function zdotc_z_z(n,a,inca, b,incb)
  integer, intent(in) :: n, inca, incb
  complex(8), intent(in) :: a(*),b(*)
 
  integer :: ia, ib, i

  if(inca==1 .and. incb==1) then
    zdotc_z_z = sum(conjg(a(1:n))*b(1:n))
  else
    zdotc_z_z = 0
    ia = 1
    ib = 1
    do i=1,n
      zdotc_z_z = zdotc_z_z + conjg(a(ia))*b(ib)
      ia = ia + inca
      ib = ib + incb
    enddo
  endif
end function ! zdotc

!
!
!
complex(8) function zdotc_d2_z(n,a,inca, b,incb)
  integer, intent(in) :: n, inca, incb
  real(8), intent(in) :: a(:,:)
  complex(8), intent(in) :: b(*)

  integer :: ia, ib, i

  if(inca==1 .and. incb==1) then
    zdotc_d2_z = sum(conjg(cmplx(a(1,1:n), a(2,1:n),8))*b(1:n))
  else
    zdotc_d2_z = 0
    ia = 1
    ib = 1
    do i=1,n
      zdotc_d2_z = zdotc_d2_z + conjg(cmplx(a(1,ia), a(2,ia),8))*b(ib)
      ia = ia + inca
      ib = ib + incb
    enddo
  endif
end function ! zdotc_d2_z

!
!
!
complex(8) function zdot(n,a,inca, b,incb)
  integer, intent(in) :: n, inca, incb
  complex(8), intent(in) :: a(*),b(*)
 
  integer :: ia, ib, i

  if(inca==1 .and. incb==1) then
    zdot = sum(a(1:n)*b(1:n))
  else
    zdot = 0
    ia = 1
    ib = 1
    do i=1,n
      zdot = zdot + a(ia)*b(ib)
      ia = ia + inca
      ib = ib + incb
    enddo
  endif
end function ! zdot

#endif

!
!
!
subroutine print_matrix(ifile, name, matrix)
  integer, intent(in) :: ifile
  character(*), intent(in) :: name
  real(8), intent(in) :: matrix(:,:)

  !! internal
  integer :: j, n, m
  n = size(matrix,1)
  m = size(matrix,2)

  write(ifile,*) trim(name)
  do j=1,m
    write(ifile,'(100000g10.3)')matrix(1:n,j)
  enddo
end subroutine ! print_matrix

!
! Inverse of a square complex matrix (wrapper around LAPACK)
!
subroutine matinv_d(A)
  implicit none
  real(8), intent(inout) :: A(:,:)
  ! internal
  integer(blas_int)   :: ndim
  integer(blas_int)   :: info
  ndim = size(A,1)

  if(ndim<1) then; write(0,*)'matinv_d: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_matinv_d) then;
    if(allocated(ipiv_matinv_d)) deallocate(ipiv_matinv_d);
    allocate(ipiv_matinv_d(ndim));
    if(allocated(work_matinv_d)) deallocate(work_matinv_d);
    allocate(work_matinv_d(ndim));
    ndim_matinv_d = ndim;
  end if
  
  call DGETRF( ndim, ndim, A, ndim, IPIV_matinv_d, INFO )
  if (info/=0) write(0,*) 'matinv_d: DGETRF: INFO', INFO;
  call DGETRI( ndim, A, ndim, IPIV_matinv_d, WORK_matinv_d, ndim, INFO );
  if (info/=0) write(0,*) 'matinv_d: DGETRI: INFO', INFO;

end subroutine !

!
! Inverse of a square complex matrix (wrapper around LAPACK)
!
subroutine matinv_s(A)
  implicit none
  real(4), intent(inout) :: A(:,:)
  ! internal
  integer(blas_int)   :: ndim
  integer(blas_int)   :: info
  ndim = size(A,1)

  if(ndim<1) then; write(0,*)'matinv_s: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_matinv_s) then;
    if(allocated(ipiv_matinv_s)) deallocate(ipiv_matinv_s);
    allocate(ipiv_matinv_s(ndim));
    if(allocated(work_matinv_s)) deallocate(work_matinv_s);
    allocate(work_matinv_s(ndim));
    ndim_matinv_s = ndim;
  end if
  
  call SGETRF( ndim, ndim, A, ndim, IPIV_matinv_s, INFO )
  if (info/=0) write(0,*) 'matinv_s: SGETRF: INFO', INFO;
  call SGETRI( ndim, A, ndim, IPIV_matinv_s, WORK_matinv_s, ndim, INFO );
  if (info/=0) write(0,*) 'matinv_s: SGETRI: INFO', INFO;

end subroutine !

!
! Inverse of a square symmetric complex matrix (wrapper around LAPACK)
!
subroutine matinv_csy(A)
  implicit none
  complex(4), intent(inout) :: A(:,:)
  ! internal
  integer(blas_int)   :: ndim
  integer(blas_int)   :: info
  ndim = size(A,1)

  if(ndim<1) then; write(0,*)'matinv_csy: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_matinv_csy) then;
    if(allocated(ipiv_matinv_csy)) deallocate(ipiv_matinv_csy);
    allocate(ipiv_matinv_csy(ndim));
    if(allocated(work_matinv_csy)) deallocate(work_matinv_csy);
    allocate(WORK_matinv_csy(2))
    call CSYTRF( 'U', ndim, A, ndim, IPIV_matinv_csy, WORK_matinv_csy, -1, INFO )
    LWORK_matinv_csy = int(WORK_matinv_csy(1))
    deallocate(WORK_matinv_csy)
    allocate(work_matinv_csy(LWORK_matinv_csy));
    if(allocated(WORKI_matinv_csy)) deallocate(WORKI_matinv_csy);
    allocate(WORKI_matinv_csy(ndim*2));
    ndim_matinv_csy = ndim;
  end if
  
  call CSYTRF( 'U', ndim, A, ndim, IPIV_matinv_csy, WORK_matinv_csy, LWORK_matinv_csy, INFO );
  if (info/=0) write(0,*) 'matinv_csy: CSYTRF: INFO', INFO;
  call CSYTRI( 'U', ndim, A, ndim, IPIV_matinv_csy, WORKI_matinv_csy, INFO );
  if (info/=0) write(0,*) 'matinv_csy: CSYTRI: INFO', INFO;

end subroutine !matinv_csy

!
! Inverse of a square symmetric complex matrix (wrapper around LAPACK)
!
subroutine matinv_zsy(A)
  implicit none
  complex(8), intent(inout) :: A(:,:)
  ! internal
  integer(blas_int)   :: ndim
  integer(blas_int)   :: info
  ndim = size(A,1)

  if(ndim<1) then; write(0,*)'matinv_zsy: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_matinv_zsy) then;
    if(allocated(ipiv_matinv_zsy)) deallocate(ipiv_matinv_zsy);
    allocate(ipiv_matinv_zsy(ndim));
    if(allocated(work_matinv_zsy)) deallocate(work_matinv_zsy);
    allocate(WORK_matinv_zsy(2))
    call ZSYTRF( 'U', ndim, A, ndim, IPIV_matinv_zsy, WORK_matinv_zsy, -1, INFO )
    LWORK_matinv_zsy = int(WORK_matinv_zsy(1))
    deallocate(WORK_matinv_zsy)
    allocate(work_matinv_zsy(LWORK_matinv_zsy));
    if(allocated(WORKI_matinv_zsy)) deallocate(WORKI_matinv_zsy);
    allocate(WORKI_matinv_zsy(ndim*2));
    ndim_matinv_zsy = ndim;
  end if
  
  call ZSYTRF( 'U', ndim, A, ndim, IPIV_matinv_zsy, WORK_matinv_zsy, LWORK_matinv_zsy, INFO );
  if (info/=0) write(0,*) 'matinv_zsy: ZSYTRF: INFO', INFO;
  call ZSYTRI( 'U', ndim, A, ndim, IPIV_matinv_zsy, WORKI_matinv_zsy, INFO );
  if (info/=0) write(0,*) 'matinv_zsy: ZSYTRI: INFO', INFO;

end subroutine !matinv_zsy


!
! Inverse of a square complex matrix (wrapper around LAPACK)
!
subroutine matinv_c_simple(A)
  use m_precision, only : blas_int
  implicit none
  complex(4), intent(inout) :: A(:,:)
  ! internal
  integer(blas_int)   :: ndim
  integer(blas_int)   :: info
  ndim = size(A,1)

  if(ndim<1) then; write(0,*)'matinv_c_simple: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_matinv_c) then;
    if(allocated(ipiv_matinv_c)) deallocate(ipiv_matinv_c);
    allocate(ipiv_matinv_c(ndim));
    if(allocated(work_matinv_c)) deallocate(work_matinv_c);
    allocate(work_matinv_c(ndim));
    ndim_matinv_c = ndim;
  end if
  
  call CGETRF( ndim, ndim, A, ndim, IPIV_matinv_c, INFO )
  if (info/=0) write(0,*) 'matinv_c_simple: CGETRF: INFO', INFO;
  call CGETRI( ndim, A, ndim, IPIV_matinv_c, WORK_matinv_c, ndim, INFO );
  if (info/=0) write(0,*) 'matinv_c_simple: CGETRI: INFO', INFO;

end subroutine !matinv_c_simple

!
! Inverse of a square complex matrix (wrapper around LAPACK)
!
subroutine matinv_c_blas(n,A,lda_in)
  use m_precision, only : blas_int
  implicit none
  integer, intent(in) :: n, lda_in
  complex(4), intent(inout) :: A(lda_in,*)
  ! internal
  integer(blas_int)   :: ndim, lda, info
  ndim = n
  lda  = lda_in

  if(n<1) then; write(0,*)'matinv_c_blas: n<1 ', n; stop; endif
  if(n>ndim_matinv_c) then;
    if(allocated(ipiv_matinv_c)) deallocate(ipiv_matinv_c); allocate(ipiv_matinv_c(n))
    if(allocated(work_matinv_c)) deallocate(work_matinv_c); allocate(work_matinv_c(n))
    ndim_matinv_c = n;
  end if
  
  call CGETRF( ndim, ndim, A, lda, IPIV_matinv_c, INFO )
  if (info/=0) write(0,*) 'matinv_c_blas: CGETRF: INFO', INFO;
  call CGETRI( ndim, A, lda, IPIV_matinv_c, WORK_matinv_c, ndim, INFO );
  if (info/=0) write(0,*) 'matinv_c_blas: CGETRI: INFO', INFO;

end subroutine !matinv_c_blas


!
! Inverse of a square complex matrix (wrapper around LAPACK)
!
subroutine matinv_z(A)
  use m_precision, only : blas_int
  implicit none
  complex(8), intent(inout) :: A(:,:)
  ! internal
  integer(blas_int)   :: ndim
  integer(blas_int)   :: info
  ndim = size(A,1)


  if(ndim<1) then; write(0,*)'matinv_z: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_matinv_z) then;
    if(allocated(ipiv_matinv_z)) deallocate(ipiv_matinv_z);
    allocate(ipiv_matinv_z(ndim));
    if(allocated(work_matinv_z)) deallocate(work_matinv_z);
    allocate(work_matinv_z(ndim));
    ndim_matinv_z = ndim;
  end if
  
  call ZGETRF( ndim, ndim, A, ndim, IPIV_matinv_z, INFO )
  if (info/=0) write(0,*) 'matinv_z: ZGETRF: INFO', INFO;
  call ZGETRI( ndim, A, ndim, IPIV_matinv_z, WORK_matinv_z, ndim, INFO );
  if (info/=0) write(0,*) 'matinv_z: ZGETRI: INFO', INFO;

end subroutine !matinv_z



subroutine real_inverse(dim,A,inverse_A,success)
  use m_precision, only : blas_int
  implicit None
  ! extern:
  integer::dim
  logical::success
  double precision::A(dim,dim),inverse_A(dim,dim)
  ! intern:
  integer::IPIV(dim)
  integer(blas_int) :: info
  double precision::copy_of_A(dim,dim),work(dim)
  !write(*,*) 'enter real_inverse'
   copy_of_A=A
  !write(*,*) 'call DGETRF'
  Call DGETRF( dim, dim, copy_of_A, dim, IPIV, INFO )
  !write(*,*) 'call DGETRI'
  call DGETRI(dim,copy_of_A,dim,IPIV,WORK,dim, INFO )
  !write(*,*) 'exit DGETRI,info=  ',info
  inverse_A=copy_of_A
  if (info ==0 ) then; success=.true. ; endif
  if (info /=0 ) then; success=.false. ; endif
  !write(*,*) 'exit real_inverse'
end subroutine real_inverse

!
!
!
subroutine complex_inverse(dim,A,inverse_A,success)
  use m_precision, only : blas_int
  implicit none
  !! external
  integer :: dim
  complex(8), intent(in)  :: A(dim,dim)
  complex(8), intent(out) :: inverse_A(dim,dim)
  logical, intent(out) :: success

  !! internal
  integer(blas_int) :: dim4
  integer(blas_int) :: info
  integer(blas_int), allocatable :: IPIV(:)
  complex(8), allocatable :: copy_of_A(:,:), work(:)

  dim4 = dim
  allocate(copy_of_A(dim,dim))
  allocate(work(dim))
  allocate(ipiv(dim))

  !write(*,*) 'enter complex_inverse'
  copy_of_A=A
  !write(*,*) 'call ZGETRF'
  Call ZGETRF( dim4, dim4, copy_of_A, dim4, IPIV, INFO )
  !write(*,*) 'call ZGETRF'
  call ZGETRI( dim4, copy_of_A, dim4, IPIV, WORK, dim4, INFO )
  !write(*,*) 'exit ZGETRI,info=  ',info
  inverse_A=copy_of_A
  if (info ==0 ) then; success=.true. ; endif
  if (info /=0 ) then; success=.false. ; endif
  !write(*,*) 'exit complex_inverse'
end subroutine !complex_inverse

!
!
!
subroutine t11_inverse_sp(ndim, A, inverse_A11)
! The tridiagonal matrix algorithm (TDMA), also known as the Thomas algorithm
! only (1,1) element is computed
  implicit none
  integer, intent(in) :: ndim
  complex(4), intent(in)  :: A(:,:)
  complex(4), intent(out) :: inverse_A11

  !! internal
  integer  :: k
  complex(4)  :: m

  if(ndim<1) then; write(0,*)'t11_inverse_sp: ndim<1 ', ndim; stop; endif
  if(ndim/=ndim_t11_inverse_sp) then;
    if(allocated(d)) deallocate(d); allocate(d(ndim));
    if(allocated(aa)) deallocate(aa); allocate(aa(ndim));
    if(allocated(b)) deallocate(b); allocate(b(ndim));
    if(allocated(c)) deallocate(c); allocate(c(ndim));
    ndim_t11_inverse_sp = ndim;
  end if

  d=0;
  d(1)  = 1;
  b(1)  = A(1,1); 
  aa(1) = 0; 
  c(1)  = A(1,2);

  do k = 2, ndim-1
     b(k)  = a(k,k)
     aa(k) = a(k,k-1)
     c(k)  = a(k,k+1)
  end do
  b(ndim) = A(ndim,ndim) ; aa(ndim) = A(ndim,ndim-1) ; c(ndim) = 0; 

  do k = 2, ndim
     m    = aa(k) / b(k-1)
     b(k) = b(k) - m*c(k-1)
     d(k) = d(k) - m*d(k-1)
  end do

  d(ndim) = d(ndim) / b(ndim)
  do k = ndim-1, 1, -1
     d(k) =  ( d(k) - c(k) *d(k+1) ) / b(k) 
  end do
  inverse_A11 = d(1)

end subroutine !t11_inverse_sp


!
!
!
subroutine t11_inverse(dim,A,inverse_A)
! The tridiagonal matrix algorithm (TDMA), also known as the Thomas algorithm
  implicit none
  integer, intent(in) :: dim
  complex(8), intent(in) :: A(dim,dim)
  complex(8), intent(out) :: inverse_A(dim,dim)

  integer :: k
  complex(8) :: m,d(dim),aa(dim),b(dim),c(dim)
  !write(*,*) 'enter complex_inverse'
   d = (0., 0.) ; d(1) = (1.0,0.0)
   b(1) = a(1,1) ; aa(1) = (0., 0.)  ; c(1) = a(1,2)
  do k = 2, dim-1
     b(k) = a(k,k)
     aa(k) = a(k,k-1)
     c(k) = a(k,k+1)
  end do
  b(dim) = a(dim,dim) ; aa(dim) = a(dim,dim-1) ; c(dim) =(0., 0.) 

  do k = 2, dim
!!$     m =  copy_of_a(k,k-1) /  copy_of_a(k-1,k-1) 
!!$     copy_of_a(k,k) =  copy_of_a(k,k) - m* copy_of_a(k-1,k)
!!$     d(k)             =   d(k) - m *d(k-1)
     m    = aa(k)/b(k-1)
     b(k) = b(k) -m*c(k-1)
     d(k) = d(k) - m *d(k-1)
  end do
!  d(dim) = d(dim)/ copy_of_a(dim,dim)
  d(dim) = d(dim)/ b(dim)
  do k = dim-1,1,-1
!     d(k)             =  ( d(k) - copy_of_a(k,k+1) *d(k+1))/ copy_of_a(k,k) 
     d(k)             =  ( d(k) - c(k) *d(k+1))/ b(k) 
  end do
  inverse_A(1,1) = d(1)

end subroutine t11_inverse

!
!
!
subroutine SimpleEigen(dim,H,E,success)
  Implicit None
  integer, intent(in)  :: dim
  logical, intent(out) :: success
  real(8), intent(in)  :: H(dim,dim)
  real(8), intent(out) :: E(dim)

!  internal
  integer :: i
  real(8), allocatable :: S(:,:),X(:,:)

! executable

  allocate(S(dim,dim), X(dim,dim));
  S=0; do i=1,dim; S(i,i)=1; enddo

  call GeneralEigen(dim,H,S,X,E,success)

  deallocate(S,X);
end subroutine SimpleEigen


subroutine GeneralEigen(dim,H,S,X,E,success)
  use m_precision, only : blas_int
  ! H X=E S X verallgemeinertes Eigenwertproblem
  implicit none
  integer, intent(in)   :: dim
  real(8), intent(in)  :: H(:,:),S(:,:)
  real(8), intent(out) :: X(:,:),E(:)
  logical,  intent(out) :: success

  !internal
  real(8), allocatable :: A(:,:), B(:,:), W(:), WORK(:)
  character :: JOBZ, UPLO
  integer(blas_int) :: ITYPE, N, LDA, LDB, LWORK, INFO 

  allocate(A(dim,dim), B(dim,dim), W(dim), WORK(1));

  ! erster Aufruf um workspace zu finden
  ITYPE=1; JOBZ='V'; UPLO='U'; N=int(dim,4); LDA=int(dim,4); LDB=int(dim,4); LWORK=int(-1,4); WORK=int(0,4);
  A=H; B=S;
  call DSYGV(ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,LWORK,INFO);
  ! zweiter Aufruf mit korrekter Dimension von workspace
  LWORK=int(WORK(1),4);
  deallocate(WORK);
  allocate(WORK(LWORK));
  WORK = 0
  call DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,LWORK,INFO)
  X=transpose(A)
  E=W
  success = (INFO==0)
  deallocate(A,B,W,WORK);

end subroutine GeneralEigen

!
!
!
subroutine sSymmGeneralEigen(H,S,E,X)
  use m_precision, only : blas_int
  ! H X=E S X verallgemeinertes Eigenwertproblem
  implicit none
  real(4), intent(in)  :: H(:,:), S(:,:)
  real(4), intent(out) :: X(:,:), E(:)

  !internal
  integer(blas_int)  :: n, lwork, info
  real(4), allocatable :: WORK(:), B(:,:)
  
  n = size(H,1)
  LWORK=3*n;
  allocate(WORK(LWORK));
  allocate(B(n,n))
  X = H
  B = S
  call SSYGV( 1, 'V', 'U', n, X, n, B, n, E, WORK,LWORK,INFO)
  if(info/=0) then; write(0,*) 'sSymmGeneralEigen: (info/=0) ==> stop'; stop; endif

end subroutine !sSymmGeneralEigen

!
!
!
subroutine dSymmGeneralEigen(H,S,E,X)
  use m_precision, only : blas_int
  ! H X=E S X verallgemeinertes Eigenwertproblem
  implicit none
  real(8), intent(in)  :: H(:,:), S(:,:)
  real(8), intent(out) :: X(:,:), E(:)

  !internal
  integer(blas_int)  :: n, lwork, info
  real(8), allocatable :: WORK(:), B(:,:)
  
  n = size(H,1)
  LWORK=3*n;
  allocate(WORK(LWORK));
  allocate(B(n,n))
  X = H
  B = S
  call DSYGV( 1, 'V', 'U', n, X, n, B, n, E, WORK,LWORK,INFO)
  if(info/=0) then; write(0,*) 'dSymmGeneralEigen: (info/=0) ==> stop'; stop; endif

end subroutine !dSymmGeneralEigen

!
!
!
subroutine dSymmGeneralEigen_blas_like(n,H,ldh,S,lds,E,X,ldx)
  use m_precision, only : blas_int
  ! H X=E S X verallgemeinertes Eigenwertproblem
  implicit none
  integer, intent(in) :: n, ldh, lds, ldx
  real(8), intent(in)  :: H(ldh,*), S(lds,*)
  real(8), intent(out) :: X(ldx,*), E(*)

  !internal
  integer(blas_int)  :: nb, ldxb, lwork, info
  real(8), allocatable :: WORK(:), B(:,:)
  
  nb = n
  ldxb = ldx
  LWORK=3*nb;
  allocate(WORK(LWORK));
  allocate(B(nb,nb))
  X(1:nb,1:nb) = H(1:nb,1:nb)
  B = S(1:nb,1:nb)
  call DSYGV( 1, 'V', 'U', nb, X, ldxb, B, nb, E, WORK,LWORK,INFO)

  if(info/=0) then; write(0,*) 'dSymmGeneralEigen_blas_like: (info/=0) ==> stop'; stop; endif

end subroutine !dSymmGeneralEigen_blas_like


!
!
!
subroutine sSymmetricEigen(A, E, X)
  use m_precision, only : blas_int
!! Debug/test
  implicit none
  real(4), intent(in)  :: A(:,:)
  real(4), intent(out) :: E(:), X(:,:)

  !! internal
  integer(blas_int) :: ndim
  real(4), allocatable :: A_copy(:,:), WORK(:)
  integer(blas_int) :: info, lwork
  real(4) :: WORK_TMP(1)

  ndim = size(A,1)
  allocate(A_copy(ndim, ndim));
  call SSYEV('V', 'U', ndim, A_copy, ndim, E, WORK_TMP, -1, info)
  if(info/=0) then; write(0,*)'sSymmetricEigen: info/=0: calculate lwork:', info; stop; endif;

  lwork = int(WORK_TMP(1));
  allocate(WORK(lwork));
  A_copy = A(1:ndim, 1:ndim);
  call SSYEV('V', 'U', ndim, A_copy, ndim, E, WORK, lwork, info);
  if(info/=0) then; write(0,*)'sSymmetricEigen: info/=0: diagonalize:', info; stop; endif;

  X = A_copy
  _dealloc(A_copy)
  _dealloc(WORK)
end subroutine ! DiagonalizeSymmetricMatrix

!
!
!
subroutine dSymmetricEigen(A, E, X)
  use m_precision, only : blas_int
!! Debug/test
  implicit none
  real(8), intent(in)  :: A(:,:)
  real(8), intent(out) :: E(:), X(:,:)

  !! internal
  integer(blas_int) :: ndim
  real(8), allocatable :: A_copy(:,:), WORK(:)
  integer(blas_int) :: info, lwork
  real(8) :: WORK_TMP(1)

  ndim = size(A,1)
  allocate(A_copy(ndim, ndim));
  call DSYEV('V', 'U', ndim, A_copy, ndim, E, WORK_TMP, -1, info)
  if(info/=0) then; write(0,*)'dSymmetricEigen: info/=0: calculate lwork:', info; stop; endif;

  lwork = int(WORK_TMP(1));
  allocate(WORK(lwork));
  A_copy = A(1:ndim, 1:ndim);
  call DSYEV('V', 'U', ndim, A_copy, ndim, E, WORK, lwork, info);
  if(info/=0) then; write(0,*)'dSymmetricEigen: info/=0: diagonalize:', info; stop; endif;

  X = A_copy
  _dealloc(A_copy)
  _dealloc(WORK)

end subroutine ! DiagonalizeSymmetricMatrix

!
! CHECK thisss
!
subroutine zHermitianEigen(A, E, X)
  use m_precision, only : blas_int
!! Debug/test
  implicit none
  !integer, intent(in)     :: ndim
  complex(8), intent(in)  :: A(:,:)
  real(8), intent(out)    :: E(:)
  complex(8), intent(out) :: X(:,:)

  !! internal
  integer(blas_int) :: ndim
  complex(8), allocatable :: A_copy(:,:)
  complex(8), allocatable :: CWORK(:)
  real(8), allocatable    :: RWORK(:)
  integer(blas_int) :: info, lwork
  complex(8) :: CWORK_TMP(1)

  ndim = size(A,1)
  allocate(A_copy(ndim, ndim));
  allocate(RWORK(3*ndim-2))
  call ZHEEV('V', 'U', ndim, A_copy, ndim, E, CWORK_TMP, -1, RWORK, info)
  if(info/=0) then; write(0,*)'zSymmetricEigen: info/=0: calculate lwork:', info; stop; endif;
  lwork = int(CWORK_TMP(1));

  allocate(CWORK(lwork));
  A_copy = A(1:ndim, 1:ndim);
  call ZHEEV('V', 'U', ndim, A_copy, ndim, E, CWORK, lwork, RWORK, info);
  if(info/=0) then; write(0,*)'zSymmetricEigen: info/=0: diagonalize:', info; stop; endif;

  X = A_copy
end subroutine ! zHermitianEigen

!
!
!
subroutine zHermitGenDiag(A,B,  E,X)
  use m_precision, only : blas_int
  implicit none
  complex(8), intent(in)  :: A(:,:), B(:,:)
  real(8), intent(out)    :: E(:)
  complex(8), intent(out) :: X(:,:)

  !! internal
  integer(blas_int) :: ndim
  complex(8), allocatable :: A_copy(:,:), B_copy(:,:)
  complex(8), allocatable :: ZWORK(:)
  complex(8) :: ZWORK_TMP(1)
  real(8), allocatable    :: RWORK(:)
  integer(blas_int) :: info, lwork

  ndim = size(A,1)
  allocate(A_copy(ndim, ndim));
  allocate(B_copy(ndim, ndim));
  allocate(RWORK(3*ndim-2));

  LWORK=-1
  call ZHEGV( 1, 'V', 'U', ndim, A_copy, ndim, B_copy, ndim, E, ZWORK_TMP, LWORK, RWORK, INFO);
  if(info/=0) then;
    write(0,*)'zHermitGenDiag: info/=0: calculate lwork:', info;
    _die('info/=0')
  endif;

  lwork = int(ZWORK_TMP(1));
  allocate(ZWORK(lwork));
  A_copy = A(1:ndim, 1:ndim);
  B_copy = B(1:ndim, 1:ndim);

  call ZHEGV( 1, 'V', 'U', ndim, A_copy, ndim, B_copy, ndim, E, ZWORK, LWORK, RWORK, INFO);
  if(info/=0) then;
    write(0,*) 'ndim', ndim
    write(0,*) 'sum(abs(A))', sum(abs(A))
    write(0,*) 'sum(abs(B))', sum(abs(B))
    !write(0,*) 'B', B
    write(0,*)'zHermitGenDiag: info/=0: diagonalize:', info; 
    _die('info/=0')
  endif;
  X = A_copy

  _dealloc(A_copy);
  _dealloc(B_copy);
  _dealloc(RWORK);
  _dealloc(ZWORK);

end subroutine ! zHermitGenDiag

!
!
!
subroutine zHermitGenDiag_evals(A,B,E)
  use m_precision, only : blas_int
  implicit none
  complex(8), intent(in)  :: A(:,:), B(:,:)
  real(8), intent(inout)    :: E(:)
  !! internal
  complex(8) :: ZTMP(1)
  complex(8), allocatable :: ZWORK(:)
  integer(blas_int) :: n, info, lwork
  real(8), allocatable    :: RWORK(:)

  n = size(A,1)
  if(n<1) _die('n<1')
  if(n>size(A,2)) _die('n>size(A,2)')
  if(any(ubound(B)/=n)) _die('any(ubound(B)<n)')
  if(size(E)<n) _die('size(E)<n')
  
  allocate(RWORK(3*n-2));

  LWORK=-1
  call ZHEGV( 1, 'N', 'U', n, A,n, B,n, E, ZTMP, LWORK, RWORK, INFO);
  if(info/=0) then;
    write(0,*)'zHermitGenDiag_evals: info/=0: calculate lwork:', info;
    _die('info/=0')
  endif;

  lwork = int(ZTMP(1));
  allocate(ZWORK(lwork));

  call ZHEGV( 1, 'N', 'U', n, A,n, B,n, E, ZWORK, LWORK, RWORK, INFO);
  if(info/=0) then;
    write(0,*) 'n', n
    write(0,*) 'sum(abs(A))', sum(abs(A))
    write(0,*) 'sum(abs(B))', sum(abs(B))
    !write(0,*) 'B', B
    write(0,*)'zHermitGenDiag_evals: info/=0: diagonalize:', info; 
    _die('info/=0?')
  endif;

  _dealloc(RWORK);
  _dealloc(ZWORK);

end subroutine ! zHermitGenDiag_ev

!
!
!
subroutine SimpleDiagonalize(ndim,H,X,E,success)
  use m_sort, only : qsort 

  ! Eigenvalues in descending order, largest first as E(1)
  implicit none
  integer, intent(in)  :: ndim
  real(8), intent(in)  :: H(:,:)
  real(8), intent(out) :: X(:,:), E(:)
  logical, intent(out) :: success

  !! internal
  integer  :: i
  integer, allocatable :: t(:)
  real(8), allocatable :: minusE(:), E_disordered(:), S(:,:), X_disordered(:,:)

  allocate(minusE(ndim))
  allocate(E_disordered(ndim))
  allocate(S(ndim,ndim))
  allocate(X_disordered(ndim,ndim))
  allocate(t(ndim))

  S=0;
  if (ndim ==1) then 
      E(1)=H(1,1)
      X(1,1)=1.d0  
  else if (ndim >1) then 
    do i=1,ndim; S(i,i)=dble(1); t(i)=i; enddo;
    call GeneralEigen(ndim,H,S,X_disordered,E_disordered,success);
    minusE=-E_disordered
    call qsort(minusE,ndim,t)
    do i=1,ndim
        E(i)  = E_disordered(t(i)  )
        X(i,:)= X_disordered(t(i),:)
    enddo
  endif ! ndim >1

  _dealloc(minusE)
  _dealloc(E_disordered)
  _dealloc(S)
  _dealloc(X_disordered)
  _dealloc(t)

end subroutine !SimpleDiagonalize

!
!
!
subroutine complex_eigen(dim,A,X,eigenvalues,success)
  use m_precision, only : blas_int
  implicit None
  ! Extern:
  integer::dim
  complex(8), intent(in)    :: A(dim,dim)
  complex(8), intent(inout) :: X(dim,dim)
  real(8), intent(inout)    :: eigenvalues(dim)
  logical::success
  ! Intern
  complex(8)::copy_of_A(dim,dim),WORK(6*dim)
  real(8)::RWORK(1:3*dim)
  integer(blas_int)::INFO,dim4
  copy_of_A=A
  dim4 = dim
  call ZHEEV('V','U',dim4, copy_of_A, dim, eigenvalues, WORK, 6*dim, RWORK, INFO )
  if (info ==0 ) then; success=.true. ; endif
  if (info /=0 ) then; success=.false. ; endif
  X=transpose(copy_of_A)
end subroutine !complex_eigen

!
!
!
subroutine nonhermitian_eigen(A,X_L,X_R,eig,success)
  use m_precision, only : blas_int
  implicit None
  ! Extern:
  real(8), intent(in)    :: A(:,:)
  complex(8), intent(out) :: X_L(:,:), X_R(:,:)
  complex(8), intent(out)    :: eig(:)
  logical::success
  ! Intern
  integer:: N, LDA, LDVL, LDVR, LWORK, INFO
  real(8), allocatable:: WR(:), WI(:), VL(:,:), VR(:,:), WORK(:)
 
  !SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
  !     $                  LDVR, WORK, LWORK, INFO )
  N = size(A,1)
  LDA = N
  allocate(WR(N)) 
  allocate(WI(N)) 
  allocate(VL(N,N)) 
  LDVL = N
  allocate(VR(N,N))
  LDVR = N
  LWORK = 6*N
  allocate(WORK(LWORK))
  
  call DGEEV( 'V', 'V', N, A, LDA, WR, WI, VL, LDVL, VR,LDVR, WORK, LWORK, INFO )
  
  !call ZHEEV('V','U',dim4, copy_of_A, dim, eigenvalues, WORK, 6*dim, RWORK, INFO )
  if (info ==0 ) then
    success=.true. 
  else
    success=.false. 
  endif

  write(6,*) "nonhermitian_eigen, success", success
  
  !X=transpose(copy_of_A)
  
  ! OBS assume real eigenvalues for the eigenvectors 
  if(sum(WI**2) .lt. 1d-8) then    
  else
    write(6,*) "nonhermitian_eigen: complex eigenvalues, the eigenvectors cannot be trusted"
  end if
  
  !write(6,*) "nonhermitian_eigen, sum(WR)", sum(WR)
  !write(6,*) "nonhermitian_eigen, sum(WR**2)", sum(WR**2)
  !write(6,*) "nonhermitian_eigen, sum(WI)", sum(WI)
  !write(6,*) "nonhermitian_eigen, sum(WI**2)", sum(WI**2)
  !write(6,*) "nonhermitian_eigen, sum(VL)", sum(VL)
  !write(6,*) "nonhermitian_eigen, sum(VR)", sum(VR)
  
  !write(6,*) "WR"
  !do i=1, N
  !  write(6,*) WR()
  !end do

  X_L = cmplx(VL,0.0d0,8)
  X_R = cmplx(VR,0.0d0,8)
  eig = cmplx(WR, WI,8)
  
  

end subroutine nonhermitian_eigen

subroutine nonhermitian_eigen_new(A,X_L,X_R,eig,success)
  use m_precision, only : blas_int
  implicit None
  ! Extern:
  real(8), intent(in)    :: A(:,:)
  complex(8), intent(out) :: X_L(:,:), X_R(:,:)
  complex(8), intent(out)    :: eig(:)
  logical::success
  ! Intern
  integer:: N, LDA, LDVL, LDVR, LWORK, INFO
  real(8), allocatable:: WR(:), WI(:), VL(:,:), VR(:,:), WORK(:)
  real(8), allocatable:: S(:,:), S_inv(:,:)

  integer:: i, j
  
  !SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
  !     $                  LDVR, WORK, LWORK, INFO )
  N = size(A,1)
  LDA = N
  allocate(WR(N)) 
  allocate(WI(N)) 
  allocate(VL(N,N)) 
  LDVL = N
  allocate(VR(N,N))
  LDVR = N
  LWORK = 6*N
  allocate(WORK(LWORK))

  allocate(S(N,N))
  allocate(S_inv(N,N))
  
  call DGEEV( 'N', 'V', N, A, LDA, WR, WI, VL, LDVL, VR,LDVR, WORK, LWORK, INFO )
  
  !call ZHEEV('V','U',dim4, copy_of_A, dim, eigenvalues, WORK, 6*dim, RWORK, INFO )
  if (info ==0 ) then
    success=.true. 
  else
    success=.false. 
  endif

  write(6,*) "nonhermitian_eigen, success", success
  
  !X=transpose(copy_of_A)
  
  ! OBS assume real eigenvalues for the eigenvectors 
  if(sum(WI**2) .lt. 1d-8) then    
  else
    write(6,*) "nonhermitian_eigen: complex eigenvalues, the eigenvectors cannot be trusted"
  end if
  
  !write(6,*) "nonhermitian_eigen, sum(WR)", sum(WR)
  !write(6,*) "nonhermitian_eigen, sum(WR**2)", sum(WR**2)
  !write(6,*) "nonhermitian_eigen, sum(WI)", sum(WI)
  !write(6,*) "nonhermitian_eigen, sum(WI**2)", sum(WI**2)
  !write(6,*) "nonhermitian_eigen, sum(VL)", sum(VL)
  !write(6,*) "nonhermitian_eigen, sum(VR)", sum(VR)
  
  !write(6,*) "WR"
  !do i=1, N
  !  write(6,*) WR()
  !end do

  ! reconstruct left eigevectors from the right ones
  S=0.0d0 
  do i=1,N
    do j=1, N
      S(i,j) = sum(VR(:,i)* VR(:,j)) 
    end do
  end do

  ! invert
  S_inv = S
  call matinv_d(S_inv)

  call matmul_AB_d(VR, S_inv, VL)

  !!check orthogonality
  !S=0.0d0
  !call matmul_AtB_d(V_L, V_L, S)
  !do i=1, N
  !  do j=1,N
  !  end do
  !end do

  X_L = cmplx(VL,0.0d0,8)
  X_R = cmplx(VR,0.0d0,8)
  eig = cmplx(WR, WI,8)
  
end subroutine nonhermitian_eigen_new


!! FH |i> = E_i |i>  --> HFH |i> = E_i H |i> -->  HFH = |i> H^{-1} E_j <j| -- FH =|i> H^{-2} E_j <j|
!! <i | FH  =  <i |H-1 HFH  = <i| E_i  
!subroutine pseudohermitian_eigen()
!
!end subroutine pseudohermitian_eigen

subroutine transpose_d(A, A_t)
  implicit none
  real(8), intent(in)::A(:,:)
  real(8), intent(out)::A_t(:,:)

  integer:: i,j 
  
  if(size(A,1) .ne. size(A_t,2)) then
    write(6,*) "transpose_d: the number of rows of A does not match the number of columns in A_t!", &
         size(A,1),size(A_t,2)
  end if
  if(size(A_t,1) .ne. size(A,2)) then
    write(6,*) "transpose_d: the number of rows of A_t does not match the number of columns in A!", &
         size(A_t,1),size(A,2)
  end if
  
  do i=1,size(A,1)
    do j =1, size(A,2)
      A_t(i,j) = A(j,i)
    end do
  end do
  
end subroutine transpose_d

subroutine transpose_A_d(A)
  implicit none
  real(8), intent(inout)::A(:,:)
  !real(8), intent(out)::A_t(:,:)
  
  real(8), allocatable:: A_tmp(:,:)
  integer:: i,j 
  
  if(size(A,1) .ne. size(A,2)) then
    write(6,*) "transpose_A_d: A is not square!!", &
         size(A,1),size(A,2)
  end if
  
  allocate(A_tmp(size(A,1), size(A,1)))
  
  do i=1,size(A,1)
    do j =1, size(A,1)
      A_tmp(i,j) = A(j,i)
    end do
  end do
  
  A =A_tmp
  
end subroutine transpose_A_d

subroutine transpose_dc(A, A_t)
  implicit none
  complex(8), intent(in)::A(:,:)
  complex(8), intent(out)::A_t(:,:)

  integer:: i,j 
  
  if(size(A,1) .ne. size(A_t,2)) then
    write(6,*) "transpose_d: the number of rows of A does not match the number of columns in A_t!", &
         size(A,1),size(A_t,2)
  end if
  if(size(A_t,1) .ne. size(A,2)) then
    write(6,*) "transpose_d: the number of rows of A_t does not match the number of columns in A!", &
         size(A_t,1),size(A,2)
  end if
  
  do i=1,size(A,1)
    do j =1, size(A,2)
      A_t(i,j) = A(j,i)
    end do
  end do
  
end subroutine transpose_dc

subroutine transpose_A_dc(A)
  implicit none
  complex(8), intent(inout)::A(:,:)
  !real(8), intent(out)::A_t(:,:)

  complex(8), allocatable:: A_tmp(:,:)
  integer:: i,j 
  
  if(size(A,1) .ne. size(A,2)) then
    write(6,*) "transpose_A_d: A is not square!!", &
         size(A,1),size(A,2)
  end if
  
  allocate(A_tmp(size(A,1), size(A,1)))

  do i=1,size(A,1)
    do j =1, size(A,1)
      A_tmp(i,j) = A(j,i)
    end do
  end do
  
  A =A_tmp

end subroutine transpose_A_dc



subroutine matmul_Ax_d(A, x, y)
  implicit none
  real(8), intent(in):: A(:,:), x(:)
  real(8), intent(inout):: y(:)
  
  !write(6,*) "in matmul_AB_d"
  ! SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X,
  ! INCX, BETA, Y, INCY)
  
  call DGEMV('N', size(A,1), size(A,2), 1.0d0, A, size(A,1), x, &
       1, 0.0d0, Y, 1)
  
end subroutine matmul_Ax_d

subroutine matmul_Ax_z(A, x, y)
  implicit none
  complex(8), intent(in):: A(:,:), x(:)
  complex(8), intent(inout):: y(:)
  
  !write(6,*) "in matmul_AB_d"
  ! SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X,
  ! INCX, BETA, Y, INCY)
  
  call ZGEMV('N', size(A,1), size(A,2), (1.0d0, 0.0d0), A, size(A,1), x, &
       1, (0.0d0,0.0d0), Y, 1)
  
end subroutine matmul_Ax_z

subroutine matmul_Atx_d(A, x, y)
  implicit none
  real(8), intent(in):: A(:,:), x(:)
  real(8), intent(inout):: y(:)
  
  !write(6,*) "in matmul_AB_d"
  ! SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X,
  ! INCX, BETA, Y, INCY)
  
  call DGEMV('T', size(A,1), size(A,2), 1.0d0, A, size(A,1), x, &
       1, 0.0d0, Y, 1)
  
end subroutine 


!
!
!
subroutine matmul_AB_d(A, B, D)
  implicit none
  real(8), intent(in):: A(:,:), B(:,:)
  real(8), intent(inout):: D(:,:)
  
  !write(6,*) "in matmul_AB_d"
  ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  call DGEMM('N', 'N', size(A,1), size(B,2), size(A,2), 1.0D0, &
       A, size(A,1), B, size(B,1), 0.0D0, D, size(D,1))
  !write(6,*) "in matmul_AB_d"
  
end subroutine matmul_AB_d


subroutine matmul_AB_z(A, B, D)
  implicit none
  complex(8), intent(in):: A(:,:), B(:,:)
  complex(8), intent(inout):: D(:,:)
  
  ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  call ZGEMM('N', 'N', size(A,1), size(B,2), size(A,2), (1.0D0, 0.0D0), &
       A, size(A,1), B, size(B,1), (0.0D0, 0.0d0), D, size(D,1))
  !write(6,*) "in matmul_AB_d"
  
end subroutine matmul_AB_z


subroutine matmul_AtB_d(A, B, D)
  implicit none
  real(8), intent(in):: A(:,:), B(:,:)
  real(8), intent(inout):: D(:,:)

  !write(6,*) "in matmul_AtB_d"
  ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  call DGEMM('T', 'N', size(A,2), size(B,2), size(A,1), 1.0D0, &
       A, size(A,1), B, size(B,1), 0.0D0, D, size(D,1))
  !write(6,*) "in matmul_AtB_d"

end subroutine matmul_AtB_d


subroutine matmul_ABt_d(A, B, D)
  implicit none
  real(8), intent(in):: A(:,:), B(:,:)
  real(8), intent(inout):: D(:,:)

  !write(6,*) "in matmul_ABt_d"
  ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  call DGEMM('N', 'T', size(A,1), size(B,1), size(A,2), 1.0D0, &
       A, size(A,1), B, size(B,1), 0.0D0, D, size(D,1))
  !write(6,*) "in matmul_ABt_d*"

end subroutine matmul_ABt_d

! function for double matrix multiplication 
! D =  matmul(transpose(A), matmul(B, C))
subroutine matmul_AtBC_d_nowork(A, B, C, D)
  implicit none
  real(8), intent(in):: A(:,:), B(:,:), C(:,:)
  real(8), intent(inout):: D(:,:)
  
  !real(8):: BC_tmp(size(B,1), size(C,2))
  real(8), allocatable:: BC_tmp(:,:)

  allocate(BC_tmp(size(B,1), size(C,2)))

  !!write(6,*) "in matmul_AtBC_d"
  !! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  !call DGEMM('N', 'N', size(B,1), size(C,2), size(B,2), 1.0D0, &
  !     B, size(B,1), C, size(C,1), 0.0D0, BC_tmp, size(BC_tmp,1))
  !
  !call DGEMM('T', 'N', size(A,2), size(BC_tmp,2), size(A,1), 1.0D0, &
  !     A, size(A,1), BC_tmp, size(BC_tmp,1), 0.0D0, D, size(D,1))

  call matmul_AtBC_d_work(A, B, C, D, BC_tmp)

  deallocate(BC_tmp)

end subroutine matmul_AtBC_d_nowork

! function for double matrix multiplication 
! D =  matmul(transpose(A), matmul(B, C))
subroutine matmul_AtBC_d_work(A, B, C, D, BC_tmp)
  implicit none
  real(8), intent(in):: A(:,:), B(:,:), C(:,:)
  real(8), intent(inout):: D(:,:), BC_tmp(:,:)
  
  !real(8):: BC_tmp(size(B,1), size(C,2))
  !real(8), allocatable:: BC_tmp(:,:)

  !allocate(BC_tmp(size(B,1), size(C,2)))

  !write(6,*) "in matmul_AtBC_d"
  ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  call DGEMM('N', 'N', size(B,1), size(C,2), size(B,2), 1.0D0, &
       B, size(B,1), C, size(C,1), 0.0D0, BC_tmp, size(BC_tmp,1))

  call DGEMM('T', 'N', size(A,2), size(BC_tmp,2), size(A,1), 1.0D0, &
       A, size(A,1), BC_tmp, size(BC_tmp,1), 0.0D0, D, size(D,1))

end subroutine matmul_AtBC_d_work

! function for double matrix multiplication 
! D =  matmul(transpose(A), matmul(B, C))
subroutine matmul_AtBC_z_nowork(A, B, C, D)
  implicit none
  complex(8), intent(in):: A(:,:), B(:,:), C(:,:)
  complex(8), intent(inout):: D(:,:)
  
  !real(8):: BC_tmp(size(B,1), size(C,2))
  complex(8), allocatable:: BC_tmp(:,:)

  allocate(BC_tmp(size(B,1), size(C,2)))

  !!write(6,*) "in matmul_AtBC_d"
  !! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  !call DGEMM('N', 'N', size(B,1), size(C,2), size(B,2), 1.0D0, &
  !     B, size(B,1), C, size(C,1), 0.0D0, BC_tmp, size(BC_tmp,1))
  !
  !call DGEMM('T', 'N', size(A,2), size(BC_tmp,2), size(A,1), 1.0D0, &
  !     A, size(A,1), BC_tmp, size(BC_tmp,1), 0.0D0, D, size(D,1))

  call matmul_AtBC_z_work(A, B, C, D, BC_tmp)

  deallocate(BC_tmp)

end subroutine matmul_AtBC_z_nowork

! function for double matrix multiplication 
! D =  matmul(transpose(A), matmul(B, C))
subroutine matmul_AtBC_z_work(A, B, C, D, BC_tmp)
  implicit none
  complex(8), intent(in):: A(:,:), B(:,:), C(:,:)
  complex(8), intent(inout):: D(:,:), BC_tmp(:,:)
  
  !real(8):: BC_tmp(size(B,1), size(C,2))
  !real(8), allocatable:: BC_tmp(:,:)

  !allocate(BC_tmp(size(B,1), size(C,2)))

  !write(6,*) "in matmul_AtBC_d"
  ! SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  call ZGEMM('N', 'N', size(B,1), size(C,2), size(B,2), (1.0D0,0.0D0), &
       B, size(B,1), C, size(C,1), (0.0D0,0.0D0), BC_tmp, size(BC_tmp,1))

  call ZGEMM('T', 'N', size(A,2), size(BC_tmp,2), size(A,1), (1.0D0,0.0D0), &
       A, size(A,1), BC_tmp, size(BC_tmp,1), (0.0D0, 0.0D0), D, size(D,1))

end subroutine matmul_AtBC_z_work



! this routine projects transforms A with X (X^T A X), takes the diagonal elements and transforms back 
subroutine make_A_diag_c(A, X, S, A_diag)
  complex(8), intent(in):: A(:,:)
  real(8), intent(in):: X(:,:)
  real(8), intent(in):: S(:,:)
  complex(8), intent(out):: A_diag(:,:)

  complex(8), allocatable:: A_ij(:,:), A_ii(:,:), SXt(:,:)
  real(8), allocatable:: A_ij_real(:,:), a_ij_imag(:,:) 
  integer::i, norbs

  norbs = size(A,1)

  !write(6,*) "shape(A), shape(X), shape(S), shape(A_diag)",shape(A), shape(X), shape(S), shape(A_diag)
  allocate(A_ij(norbs,norbs))
  allocate(A_ij_real(norbs,norbs))
  allocate(A_ij_imag(norbs,norbs))
  allocate(A_ii(norbs,norbs))
  allocate(SXt(norbs,norbs))

  !A_ij = matmul( transpose(X), matmul(A, X))  

  call matmul_AtBC_d(X,real(A),X, A_ij_real)
  call matmul_AtBC_d(X,aimag(A),X, A_ij_imag)

  A_ij = cmplx(A_ij_real, A_ij_imag, 8)

  A_ii = 0.0d0
  do i=1,norbs
    A_ii(i,i)  =  A_ij(i,i)
  end do
 
  ! reexpress self_energy in atomic basis
  SXt = matmul(S, transpose(X))
  A_diag = matmul( transpose(SXt), matmul(A_ii, SXt))  

end subroutine !make_A_diag_c

!
!
!
subroutine dealloc_algebra()
  implicit none

  !$OMP PARALLEL
  _dealloc(d)
  _dealloc(aa)
  _dealloc(b)
  _dealloc(c)
  _dealloc(ipiv_matinv_z)
  _dealloc(work_matinv_z)
  _dealloc(ipiv_matinv_c)
  _dealloc(work_matinv_c)
  _dealloc(ipiv_matinv_d)
  _dealloc(work_matinv_d)
  _dealloc(ipiv_matinv_s)
  _dealloc(work_matinv_s)
  _dealloc(ipiv_matinv_csy)
  _dealloc(work_matinv_csy)
  _dealloc(worki_matinv_csy)
  _dealloc(ipiv_matinv_zsy)
  _dealloc(work_matinv_zsy)
  _dealloc(worki_matinv_zsy)
  !$OMP END PARALLEL

end subroutine !dealloc_algebra

end module !m_algebra


