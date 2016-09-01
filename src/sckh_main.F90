program sckh_main
  use m_sckh_params_t, only: sckh_params_t
  use m_sckh_params_t, only: init_sckh_params_t
  use m_upper, only : upper
  use m_KH, only: calculate_KH_nonres
  use m_KH, only: calculate_KH_res
  use m_KH, only: calculate_KH_res_el
  use m_log, only : init_logs, ilog, log_timing_note
  use m_SCKH_PES, only : calculate_SCKH_PES
  use m_SCKH, only : calculate_SCKH,compute_sckh_diagonal_nonresonant,&
       compute_sckh_offdiagonal
  !use m_fact, only : init_fact

  implicit none

  type(sckh_params_t):: p
  integer :: iv
  iv = 1


  !
  ! This progam calclates the XAS or XES cross section including non-adiabatic couplings
  !

  !call read_input(input_par)
  call init_logs(-1, '.SBE_multband.txt', 1, iv, ilog)
  !call init_fact(iv, ilog)

  call init_sckh_params_t(p)
  !write(*,*) ' ',p%runmode_sckh

  if(upper(p % runmode) .eq. "KH") then
    call calculate_KH_nonres(p)
  else  if(upper(p % runmode) .eq. "KH_RESONANT") then
    call calculate_KH_res(p)
  else  if(upper(p % runmode) .eq. "KH_RESONANT_EL") then
    call calculate_KH_res_el(p)
  else if (upper(p % runmode) .eq. "SCKH_PES") then
    call calculate_SCKH_PES(p)
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
    write(6,*) "runmode must be either 'KH','SCKH_PES', or 'SCKH' ", upper(p % runmode)
    stop
  end if

end program sckh_main

