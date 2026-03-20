program check_E_mean
  use m_sckh_params_t, only: sckh_params_t, init_sckh_params_t
  use m_upper, only : upper
  use m_log, only : init_logs, ilog, log_timing_note
  use m_SCKH_utils, only: read_one_sckh_traj, read_one_sckh_res_traj
  use m_check_E_mean, only: check_e_mean_nonresonant, check_e_mean_resonant

  implicit none

  type(sckh_params_t):: p
  integer :: iv
  iv = 1

  call init_logs(-1, '.SCKH.txt', 1, iv, ilog)

  call init_sckh_params_t(p)

  select case(upper(p % runmode))
    
  case("SCKH")

    select case(upper(p%runmode_sckh))
    case( "NONRESONANT")
      call check_e_mean_nonresonant(p)

    case("RESONANT")
      call check_e_mean_resonant(p)

    end select

  end select

end program check_E_mean
