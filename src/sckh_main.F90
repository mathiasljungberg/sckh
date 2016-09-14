program sckh_main
  use m_sckh_params_t, only: sckh_params_t
  use m_sckh_params_t, only: init_sckh_params_t
  use m_upper, only : upper
  use m_KH, only: calculate_KH_nonres
  use m_KH, only: calculate_KH_res
  use m_KH, only: calculate_KH_res_el
  use m_log, only : init_logs, ilog, log_timing_note
  use m_SCKH_PES, only : calculate_SCKH_PES
  use m_SCKH, only : calculate_SCKH
  use m_SCKH_nonadiabatic, only : compute_sckh_diagonal_nonresonant
  use m_SCKH_nonadiabatic, only : compute_sckh_offdiagonal
  use m_XAS_eigenstates, only: calculate_XAS
  use m_SCKH_resonant_PES, only: calculate_SCKH_res_PES
  
  implicit none

  type(sckh_params_t):: p
  integer :: iv
  iv = 1

  !
  ! This progam calclates the XAS or XES cross section including non-adiabatic couplings
  !

  call init_logs(-1, '.SCKH.txt', 1, iv, ilog)

  call init_sckh_params_t(p)

  ! routines working in the eigestate basis, 1d vibrational PES
  if(upper(p % runmode) .eq. "XAS") then
    call calculate_XAS(p)
    
  else if(upper(p % runmode) .eq. "KH") then
    call calculate_KH_nonres(p)

  else  if(upper(p % runmode) .eq. "KH_RESONANT") then
    call calculate_KH_res(p)
    
  else  if(upper(p % runmode) .eq. "KH_RESONANT_EL") then
    call calculate_KH_res_el(p)

  ! routines working with 1d  PES, but use classical trajectories on that
  else if (upper(p % runmode) .eq. "SCKH_PES") then
    call calculate_SCKH_PES(p)

  else if (upper(p % runmode) .eq. "SCKH_RESONANT_PES") then
    call calculate_SCKH_res_PES(p)
    
  ! rotines working with general used-supplied trajectories
  else if (upper(p % runmode) .eq. "SCKH") then
    if  (upper(p%runmode_sckh) .eq. "NONRESONANT") then
      call calculate_SCKH(p)
    endif
    if (upper(p%runmode_sckh) .eq. "NONRESONANT_DIAGONAL") then
      call compute_sckh_diagonal_nonresonant(p)
    endif
    if (upper(p%runmode_sckh) .eq. "NONRESONANT_OFFDIAGONAL") then
      call compute_sckh_offdiagonal(p)
    endif
    
  else 
    write(6,*) "runmode must be one of 'XAS', 'KH', 'KH_RESONANT', 'KH_RESONANT_EL',&
         'SCKH_PES', 'SCKH' ", upper(p % runmode)
    stop
  end if
  
end program sckh_main

