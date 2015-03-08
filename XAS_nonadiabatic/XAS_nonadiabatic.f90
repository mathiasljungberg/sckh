program XAS_nonadiabatic
  !use parameters
  !use KH_functions
  !use spline_m
  !use m_XAS_nonadiabatic
  use m_XAS_functions
  use m_test_XAS_nonadiabatic
  use m_XAS_eigenstates
  use m_XES_eigenstates
  use m_SCKH_PES
  use m_SCXAS_PES
  implicit none


  type(input_params):: input_par 

  !
  ! This progam calclates the XAS or XES cross section including non-adiabatic couplings
  !

  call read_input(input_par)

  if(input_par % runmode .eq. "XAS") then
     call calculate_XAS_nonadiabatic(input_par)
  else if (input_par % runmode .eq. "XES") then
     call calculate_XES_nonadiabatic(input_par)
  else if (input_par % runmode .eq. "SCKH_PES") then
     call calculate_SCKH_PES(input_par)
  else if (input_par % runmode .eq. "SCXAS_PES") then
     call calculate_SCXAS_PES(input_par)
!  else if (input_par % runmode .eq. "SCXAS_TRAJ") then
!     call calculate_SCXAS_TRAJ((input_par)
!  else if (input_par % runmode .eq. "SCKH_TRAJ") then
!     call calculate_SCKH_TRAJ((input_par)
!  else if (input_par % runmode .eq. "TEST") then
!     call test_XAS_nonadiabatic(input_par)
  end if

end program XAS_nonadiabatic

