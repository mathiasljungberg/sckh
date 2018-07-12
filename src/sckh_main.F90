program sckh_main
  use m_sckh_params_t, only: sckh_params_t
  use m_sckh_params_t, only: init_sckh_params_t
  use m_upper, only : upper
  use m_KH, only: calculate_KH_nonres
  use m_KH, only: calculate_KH_res
  use m_KH, only: calculate_KH_res_el
  use m_log, only : init_logs, ilog, log_timing_note
  use m_SCXAS_PES, only : calculate_SCXAS_PES
  use m_SCKH_PES, only : calculate_SCKH_PES
  use m_SCKH, only : calculate_SCKH
  use m_SCKH_nonadiabatic, only : compute_sckh_diagonal_nonresonant
  use m_SCKH_nonadiabatic, only : calculate_sckh_matrix_diagonal
  use m_SCKH_nonadiabatic, only : compute_sckh_offdiagonal
  use m_XAS_eigenstates, only: calculate_XAS
  use m_SCKH_resonant_PES, only: calculate_SCKH_res_PES
  use m_SCKH_resonant_PES, only: calculate_SCKH_res_PES_factor
  use m_SCKH_resonant_PES, only: calculate_SCKH_res_PES_factor_each_traj
  use m_SCKH_resonant_PES_FC, only: calculate_SCKH_res_PES_FC
  use m_SCKH_resonant_PES_traj, only: calculate_SCKH_res_PES_traj
  use m_SCKH_resonant, only: calculate_SCKH_res_FC
  use m_broaden_rixs, only: read_broaden_write_rixs
  
  implicit none

  type(sckh_params_t):: p
  integer :: iv
  iv = 1

  !
  ! This progam calclates the XAS or XES cross section including non-adiabatic couplings
  !

  call init_logs(-1, '.SCKH.txt', 1, iv, ilog)

  call init_sckh_params_t(p)

  select case(upper(p % runmode))
    
  ! routines working in the eigestate basis, 1d vibrational PES
  case("XAS")
    call calculate_XAS(p)
  case("KH")
    call calculate_KH_nonres(p)
  case( "KH_RESONANT")
    call calculate_KH_res(p)
  case("KH_RESONANT_EL")
    call calculate_KH_res_el(p)
    
  ! routines working with 1d  PES, but use classical trajectories on that
  case("SCXAS_PES")
    call calculate_SCXAS_PES(p)
  case("SCKH_PES")
    call calculate_SCKH_PES(p)
  case("SCKH_RESONANT_PES")
    call calculate_SCKH_res_PES(p)
  case("SCKH_RESONANT_PES_FACTOR")
    call calculate_SCKH_res_PES_factor(p)
  case("SCKH_RESONANT_PES_FACTOR_EACH_TRAJ")
    call calculate_SCKH_res_PES_factor_each_traj(p)
  case("SCKH_RESONANT_PES_FC")
    call calculate_SCKH_res_PES_FC(p)
  ! add by O.Takahashi 2018/06/29
  case("SCKH_RESONANT_PES_TRAJ")
    call calculate_SCKH_res_PES_traj(p)

  ! rotines working with general used-supplied trajectories
  case("SCKH")

    select case(upper(p%runmode_sckh))

    case( "NONRESONANT")
      call calculate_SCKH(p)
    case("MATRIX_DIAGONAL")
      write(6,*) "MATRIX_DIAGONAL"
      call calculate_sckh_matrix_diagonal(p)
    case("NONRESONANT_DIAGONAL")
      write(6,*) "NONRESONANT_DIAGONAL"
      call compute_sckh_diagonal_nonresonant(p)
    case("NONRESONANT_OFFDIAGONAL")
      call compute_sckh_offdiagonal(p)
    case("RESONANT")
      call calculate_SCKH_res_FC(p)

    end select

    ! utility to broaden spectra afterwards
    case("BROADEN_RIXS")
      call read_broaden_write_rixs(p)
      
  case default
    write(6,*) "runmode must be one of 'XAS', 'KH', 'KH_RESONANT', 'KH_RESONANT_EL',&
         'SCKH_PES', 'SCKH' ", upper(p % runmode)
    stop
    
  end select
  
end program sckh_main

