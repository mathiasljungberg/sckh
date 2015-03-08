module parameters
  implicit none

  !parameters
  integer,parameter::wp=8
  integer,parameter::seed=323087
  real(kind=wp),parameter:: bohr=0.52917725d0
  real(kind=wp), parameter:: pi=3.14159265d0
  real(kind=wp),parameter:: hbar=1.05457168d-34
  real(kind=wp),parameter:: amu=1.66053873d-27
  real(kind=wp),parameter:: c=2.99792458d8
  real(kind=wp),parameter:: hartree=4.35974418d-18

  real(kind=wp),parameter:: h=6.6260688d-34
  real(kind=wp),parameter:: cm = 1.0d0/(100.0d0*h*c)
  real(kind=wp),parameter:: ev =1.60217646d-19
  
end module parameters







