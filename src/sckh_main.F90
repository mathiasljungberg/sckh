program sckh_main
  use m_sckh_params_t, only: sckh_params_t
  use m_sckh_params_t, only: init_sckh_params_t
  use m_upper, only : upper
  use m_KH, only: calculate_XES_nonadiabatic

  implicit none

  type(sckh_params_t):: p

  !
  ! This progam calclates the XAS or XES cross section including non-adiabatic couplings
  !

  !call read_input(input_par)

  call init_sckh_params_t(p)

  if(upper(p % runmode) .eq. "KH") then
    call calculate_XES_nonadiabatic(p)
  else if (upper(p % runmode) .eq. "SCKH") then
    !call calculate_SCKH(p)
  else 
    write(6,*) "runmode must be either 'KH' or 'SCKH' ", upper(p % runmode)
    stop
  end if

end program sckh_main

