

real(8) function morse_eigenfun(n, D, a, my, x)
  use parameters
  implicit none

  ! passed variables
  integer, intent(in)::n
  real(8), intent(in):: D, a, my,x

  ! local variables
  real(8):: nn,k,b,z,my_SI,a_SI, D_SI, dy, dy2,ass_laguerre

  !functions
  real(8) ::factrl,gamma

  ! x [Å], a [Å-1], D [Hartree], my reduced mass in [amu]
  ! egenfunktioner som i Gallas, 1979, Physical reviev A
  ! V = D(1-exp(-ax))**2
  ! the subroutine uses valapo from splib.f, and fac and gamma from
  ! numerical recepies

  my_SI = my*amu;
  a_SI = a*1.0d10;
  D_SI = D*hartree;

  k = 2.*sqrt(2.0_wp*my_SI*D_SI)/(a_SI*hbar)
  z = k*exp(-a*x)
  b = k-2.0_wp*n -1.0_wp
  nn = sqrt(a*b*factrl(n)/gamma(k-n)) ! normalized to Å
  
  call valapo(n,b,z,ass_laguerre,dy,dy2)

  morse_eigenfun = nn*exp(-z/2.0_wp)*(z**(b/2.0_wp))*ass_laguerre
  
end function morse_eigenfun

real(8) function harm_eigenfun(n,freq, my, x)
  use parameters
  implicit none
  
  ! passed variables
  integer::n
  real(8):: freq, my,x

  ! local variables
  real(8) :: nn,s,sx,omega_SI,my_SI,y,dy,d2y,x_SI,norm
  
  !functions
  real(8) ::factrl

  ! freq [cm-1], my [amu],  [Å]
  ! egenfunktioner från wikipedia
  ! the routine uses vahepo from splib.f and fac from numerical recepies
  ! latest change 2006-10-01 to make possible higher eigenfunctions
  ! latest change 2009-02-23, just cosmetics
  ! freq = hbar*omega

  x_SI = x*1.0d-10
  my_SI = my*amu
  !omega_SI= omega*c*2.0*pi*100.0
  omega_SI= freq/(cm*hbar)
  s=my_SI*omega_SI/hbar
  sx=sqrt(s)*x_SI
  call vahepo(n,sx,y,dy,d2y) ! call fortran function
  
  nn=1./sqrt((2.0_wp**dfloat(n))*factrl(n))
  
  norm=sqrt(1.d-10)   !normaliserat till Å
  harm_eigenfun = norm*nn*((s/pi)**(1.0_wp/4.0_wp))*exp(-(s/2.0_wp)*(x_SI**2.0_wp))*y;



end function harm_eigenfun
