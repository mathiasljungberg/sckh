module m_rixs_io
  
#include "m_define_macro.F90"  
  implicit none

contains  
  

  !
  ! to write RIXS spectra
  !

  subroutine write_RIXS_map(omega_in, omega_out, lambda_lp, lambda_ln, lambda_cp, outfile, kh_print_stride)
    use m_precision, only: wp
    use m_io, only: get_free_handle

    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: lambda_lp(:,:)
    real(wp), intent(in):: lambda_ln(:,:)
    real(wp), intent(in):: lambda_cp(:,:)
    character(*), intent(in):: outfile
    integer, intent(in):: kh_print_stride

    integer:: i,j
    integer:: ifile, n_omega_in, n_omega_out
    character(80):: string, file

    n_omega_in = size(omega_in)
    n_omega_out = size(omega_out)
    
    file="_sigma_all_nogrid_nobroad"
    file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
    
    ifile = get_free_handle()
    open(ifile,file=file,status='unknown')
    
    write(ifile,*) "# ",  n_omega_in, n_omega_out

     do j=1, n_omega_in, kh_print_stride
      do i=1, n_omega_out, kh_print_stride
         write(ifile,'(5ES20.10E3)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do
      write(ifile, *) 
      write(ifile, *) 
    end do
    
    close(ifile) 
    
  end subroutine write_RIXS_map
    
  
  subroutine write_RIXS_spectra(omega_in, omega_out, lambda_lp, lambda_ln, lambda_cp, outfile, &
       kh_print_stride, flag_norm_in)
    use m_precision, only: wp
    use m_io, only: get_free_handle
    
    real(wp), intent(in):: omega_in(:)
    real(wp), intent(in):: omega_out(:)
    real(wp), intent(in):: lambda_lp(:,:)
    real(wp), intent(in):: lambda_ln(:,:)
    real(wp), intent(in):: lambda_cp(:,:)
    character(*), intent(in):: outfile
    integer, intent(in):: kh_print_stride
    logical, intent(in), optional:: flag_norm_in
    
    integer:: i,j
    integer:: ifile, n_omega_in, n_omega_out
    character(80):: string, file
    logical:: flag_norm
    real(wp):: int_lp, int_ln, int_cp, d_om_out, d_om_in
    
    if(present(flag_norm_in)) then
      flag_norm = flag_norm_in
    else
      flag_norm = .false.
    end if
    
    n_omega_in = size(omega_in)
    n_omega_out = size(omega_out)

    
    ! write spectra to individual files
    do j=1, n_omega_in, kh_print_stride
      
      file="_sigma_"
      write(string,'(F6.2)') omega_in(j)   

      file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
      
      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')
      
      do i=1, n_omega_out, kh_print_stride
        
        write(ifile,'(4ES20.10E3)') omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
      end do
      
      close(ifile) 
      
   end do



   !if(.true.) then
   
    ! write spectra to file
   file="_sigma_all"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in, kh_print_stride
     do i=1, n_omega_out, kh_print_stride
       write(ifile,'(5ES20.10E3)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
   end do
   
   close(ifile) 

   ! write spectra to file
   file="_sigma_all_nogrid"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do j=1, n_omega_in, kh_print_stride
     do i=1, n_omega_out, kh_print_stride
       write(ifile,'(5ES20.10E3)') omega_in(j), omega_out(i), lambda_lp(j,i), lambda_ln(j,i), lambda_cp(j,i)
     end do
     write(ifile, *) 
     write(ifile, *) 
   end do
   
   close(ifile) 

   ! write spectra summed over incoming frequencies
   file="_nonres"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')

   do i=1, n_omega_out, kh_print_stride
     write(ifile,'(5ES20.10E3)') omega_out(i), sum(lambda_lp(:,i)), sum(lambda_ln(:,i)), sum(lambda_cp(:,i))
   end do
   
   close(ifile) 

   ! write spectra summed over outgoing frequencies
   file="_xas"
   !write(string,'(F6.2)') omega_in(j)   
   file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
   
   ifile = get_free_handle()
   open(ifile,file=file,status='unknown')
   
   do i=1, n_omega_in, kh_print_stride
     write(ifile,'(5ES20.10E3)') omega_in(i), sum(lambda_lp(i,:)), sum(lambda_ln(i,:)), sum(lambda_cp(i,:))
   end do
   
   close(ifile) 

   
   ! write normalized emission spectra to individual files
   if (flag_norm) then

     d_om_out=(omega_out(2)-omega_out(1)) !* kh_print_stride
     d_om_in=(omega_in(2)-omega_in(1)) !* kh_print_stride

     write(6,*) "d_om_out", d_om_out
     write(6,*) "d_om_in", d_om_in
     write(6,*) "kh_print_stride", kh_print_stride
     
     do j=1, n_omega_in, kh_print_stride
       
       file="_sigma_norm_"
       write(string,'(F6.2)') omega_in(j)   
       
       file = trim(adjustl(outfile)) //  trim(adjustl(file)) // trim(adjustl(string)) // ".dat"
       
       ifile = get_free_handle()
       open(ifile,file=file,status='unknown')
       
       ! normalize spectrum
       int_lp = sum(lambda_lp(j,:)) * d_om_out * kh_print_stride
       int_ln = sum(lambda_ln(j,:)) * d_om_out * kh_print_stride
       int_cp = sum(lambda_cp(j,:)) * d_om_out * kh_print_stride
       
       do i=1, n_omega_out, kh_print_stride
         write(ifile,'(4ES20.10E3)') omega_out(i), lambda_lp(j,i)/int_lp, &
              lambda_ln(j,i)/int_ln, lambda_cp(j,i)/int_cp
       end do
       
       close(ifile) 
       
     end do
     

     ! write spectra summed over incoming frequencies
     file="_nonres_norm"
     !write(string,'(F6.2)') omega_in(j)   
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
     
     ifile = get_free_handle()
     open(ifile,file=file,status='unknown')

     ! normalize spectrum
     int_lp = sum(lambda_lp(:,:)) * d_om_out  * kh_print_stride
     int_ln = sum(lambda_ln(:,:)) * d_om_out  * kh_print_stride
     int_cp = sum(lambda_cp(:,:)) * d_om_out  * kh_print_stride
          
     do i=1, n_omega_out, kh_print_stride
       write(ifile,'(5ES20.10E3)') omega_out(i), sum(lambda_lp(:,i)) / int_lp, &
            sum(lambda_ln(:,i))/int_ln, sum(lambda_cp(:,i))/int_cp
     end do
     
     close(ifile) 
     
     ! write spectra summed over outgoing frequencies
     file="_xas_norm"
     !write(string,'(F6.2)') omega_in(j)   
     file = trim(adjustl(outfile)) //  trim(adjustl(file)) // ".dat"
     
     ifile = get_free_handle()
     open(ifile,file=file,status='unknown')
     
     ! normalize spectrum
     int_lp = sum(lambda_lp(:,:)) * d_om_in  * kh_print_stride
     int_ln = sum(lambda_ln(:,:)) * d_om_in  * kh_print_stride
     int_cp = sum(lambda_cp(:,:)) * d_om_in  * kh_print_stride

     do i=1, n_omega_in, kh_print_stride
       write(ifile,'(5ES20.10E3)') omega_in(i), sum(lambda_lp(i,:))/int_lp, &
            sum(lambda_ln(i,:))/int_ln, sum(lambda_cp(i,:))/int_cp
     end do
     
     close(ifile) 
     
   end if ! if (flag_norm) then



 end subroutine write_RIXS_spectra


end module m_rixs_io
