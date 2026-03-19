!> Mathematical and physical constants
module m_constants
  
  implicit none
  
  type constants_t
     !> Bohrradius (in A)
     real(8):: bohr !=0.52917725d0
     !> \f$ \pi \f$
     real(8):: pi !=4.0d0*atan(1.0d0) ! 3.14159265358979d0
     !> \f$ \hbar \f$
     real(8):: hbar !=1.05457168d-34
     !> Atomic mass unit (in kg)
     real(8):: u !=1.66053873d-27
     !> \em c (in m/s)
     real(8):: c !=2.99792458d8
     !> Hartree (in J)
     real(8):: hartree !=4.35974418d-18
     !> Plank constant (in Js)
     real(8):: h !=6.6260688d-34
     !> cm^{-1) in J
     real(8):: cm != 5.03411759319722d22 =1.0d0 / (100.0d0 * h * c)
     !> eV (in J)
     real(8):: eV ! =1.60217646d-19
     ! > Hartree to eV
     real(8):: Hartree2eV ! = 27.211396132d0  
     ! Boltzmann's constant (SI)
     real(8):: k_b !=1.3806503d-23
     ! Au time unit
     real(8)::autime!=2.418884326505d−17

  end type constants_t

  type(constants_t), parameter:: const = constants_t( &
       !> Bohrradius (in A) 
       0.52917725d0, &
       !> \f$ \pi \f$
       3.14159265358979d0, &
       !> \f$ \hbar \f$
       1.05457168d-34, &
       !> Atomic mass unit (in kg)
       1.66053873d-27, &
       !> \em c (in m/s)
       2.99792458d8, &
       !> Hartree (in J)
       4.35974418d-18, &
       !> Plank constant (in Js)
       6.6260688d-34, &
       !> cm^{-1) in J
       5.03411759319722d22,&
       !> eV (in J)
       1.60217646d-19, &
       !> Hartree to eV
       27.211396132d0, &
       !> Boltzmann's constant
       1.3806503d-23, &   
       !> Au time s
       2.418884326505d-17  )

end module m_constants







