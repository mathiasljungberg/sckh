module m_broaden_rixs
  
#include "m_define_macro.F90"  
  implicit none

contains  

  subroutine read_broaden_write_rixs(p)
    use m_precision, only: wp
    use m_sckh_params_t, only: sckh_params_t 
    use m_io, only: get_free_handle
    use m_rixs_io, only: write_RIXS_spectra
    
    type(sckh_params_t), intent(inout):: p 

    character(80):: file,string, dummy
    integer:: ifile, n_omega_in, n_omega_out
    integer:: i,j
    real(wp),allocatable:: omega_in(:)
    real(wp),allocatable:: omega_out(:)
    real(wp), allocatable:: lambda_lp(:,:), lambda_ln(:,:), lambda_cp(:,:)
    
    ! read RIXS spectrum from file
    file="_sigma_all_nogrid_nobroad"
    file = trim(adjustl(p % outfile)) //  trim(adjustl(file)) // ".dat"
    
    ifile = get_free_handle()
    open(ifile,file=file,status='unknown')
    
    read(ifile,*) dummy, n_omega_in, n_omega_out
    
    allocate(omega_in(n_omega_in))
    allocate(omega_out(n_omega_out))
    allocate(lambda_lp(n_omega_in, n_omega_out))
    allocate(lambda_ln(n_omega_in, n_omega_out))
    allocate(lambda_cp(n_omega_in, n_omega_out))
    
    do j=1,  n_omega_in
      do i=1,  n_omega_out
        read(ifile,*) omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do
    read(ifile, *) 
    read(ifile, *) 
  end do
  
  close(ifile) 

  ! broaden
  call broaden_rixs(omega_in, omega_out, lambda_lp, lambda_ln, lambda_cp, &
       p % gamma_inc_FWHM / 2.0_wp, p % gamma_instr_FWHM / 2.0_wp, &
       p % broadening_func_inc, p % broadening_func_instr)

  ! write 
  call write_RIXS_spectra(omega_in, omega_out, lambda_lp, lambda_ln, &
       lambda_cp, p % outfile, p % kh_print_stride, .true.)
  
end subroutine read_broaden_write_rixs

  !
  ! broaden
  !
  subroutine broaden_rixs(omega_in, omega_out, lambda_lp, lambda_ln, lambda_cp, &
       gamma_inc, gamma_instr, broadening_func_inc, broadening_func_instr) 

    use m_precision, only: wp
    use m_spectrum_utils, only: convolution_lorentzian_grid_fft_many_freq
    use m_upper, only : upper
    
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(inout):: lambda_lp(:,:)
    real(wp), intent(inout):: lambda_ln(:,:)
    real(wp), intent(inout):: lambda_cp(:,:)
    real(wp), intent(in):: gamma_inc, gamma_instr
    character(*), intent(in):: broadening_func_inc
    character(*), intent(in):: broadening_func_instr
    
    real(wp), allocatable:: sigma_tmp(:,:)
    
    allocate(sigma_tmp(size(omega_in), size(omega_out)))
    
    if(gamma_inc .gt. 1d-5) then
      
      write(6,*) "Entering convolute_incoming, broadening ",  upper(broadening_func_inc)
      
      sigma_tmp = lambda_lp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_lp, 1,  upper(broadening_func_inc))
      sigma_tmp = lambda_ln
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_ln, 1, upper(broadening_func_inc))
      sigma_tmp = lambda_cp
      call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
           2.0_wp * gamma_inc, omega_out, lambda_cp, 1, upper(broadening_func_inc))
      
    write(6,*) "Done"
    
  end if
      
  if(gamma_instr .gt. 1d-5) then
        
    write(6,*) "Entering convolute_instrumental, broadening  ",  upper(broadening_func_instr) 
        
    sigma_tmp = lambda_lp
    call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
         2.0_wp * gamma_instr, omega_out, lambda_lp, 2, upper(broadening_func_instr))
    sigma_tmp = lambda_ln
    call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
         2.0_wp * gamma_instr, omega_out, lambda_ln, 2, upper(broadening_func_instr))
    sigma_tmp = lambda_cp
    call convolution_lorentzian_grid_fft_many_freq(omega_in, sigma_tmp, &
         2.0_wp * gamma_instr, omega_out, lambda_cp, 2, upper(broadening_func_instr))
  end if
  
  write(6,*) "Done"

end subroutine broaden_rixs


end module m_broaden_rixs

  


  
