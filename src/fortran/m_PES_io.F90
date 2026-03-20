module m_PES_io

  implicit none

contains

subroutine read_PES_file(filename, npoints_in, nstates, X_dvr, E, units_in)
  use m_precision, only: wp
  use m_constants, only: const
  use m_KH_functions, only: check_dvr_bounds 
  use m_splines, only: spline_easy
  use m_upper, only : upper
  
  character(80), intent(in):: filename
  integer, intent(in):: npoints_in, nstates
  real(kind=wp), intent(inout):: X_dvr(nstates)
  real(kind=wp), intent(out):: E(nstates)
  character(*), intent(in), optional:: units_in
  
  real(kind=wp):: x_in(npoints_in), E_in(npoints_in)
  integer:: i
  character(80):: units

  if(present(units_in)) then
    units = units_in
  else
    units = "ANG"
  end if
  
  !
  ! read pes from files (allow for PES:s with arbitrary points)
  !

  open(10,file= filename,status='unknown')

  do i=1,npoints_in
     read(10,*) x_in(i), E_in(i)
  end do

  close(10) 

  if(upper(units) .eq. "ANG") then
    x_in = x_in*1d-10
  else if(upper(units) .eq. "BOHR") then
    x_in = x_in * 1d-10 * const % bohr 
  else
    write(6,*) "read_PES_file: units = ",units," must be 'ANG' or 'BOHR'"
  end if
  
   E_in = E_in *const % hartree
  
  call check_dvr_bounds(X_dvr, x_in)
  call spline_easy(x_in, E_in, npoints_in, X_dvr, E, nstates)

end subroutine read_PES_file


subroutine read_dipole_file(filename, npoints_in, nstates, X_dvr, dipole, units_in)
  use m_precision, only: wp
  use m_constants, only: const
  use m_KH_functions, only: check_dvr_bounds 
  use m_splines, only: spline_easy
  use m_upper, only : upper
  
  character(80), intent(in):: filename
  integer, intent(in):: npoints_in, nstates
  real(kind=wp), intent(inout):: X_dvr(nstates)
  real(kind=wp), intent(out):: dipole(nstates,3)
  character(*), intent(in), optional:: units_in
  
  real(kind=wp):: x_in(npoints_in), dipole_in(npoints_in,3)
  integer:: i 
  character(80):: units

  if(present(units_in)) then
    units = units_in
  else
    units = "ANG"
  end if
  
  open(10,file=filename ,status='unknown')
  
  do i=1,npoints_in
     read(10,*) x_in(i), dipole_in(i,1),dipole_in(i,2),dipole_in(i,3)
  end do
  
  close(10) 

 if(upper(units) .eq. "ANG") then
    x_in = x_in*1d-10
  else if(upper(units) .eq. "BOHR") then
    x_in = x_in * 1d-10 * const % bohr 
  else
    write(6,*) "read_dipole_file: units = ",units," must be 'ANG' or 'BOHR'"
  end if
  
  call check_dvr_bounds(X_dvr, x_in)
  call spline_easy(x_in, dipole_in(:,1), npoints_in, X_dvr, dipole(:,1), nstates)
  call spline_easy(x_in, dipole_in(:,2), npoints_in, X_dvr, dipole(:,2), nstates)
  call spline_easy(x_in, dipole_in(:,3), npoints_in, X_dvr, dipole(:,3), nstates)

end subroutine read_dipole_file


subroutine read_nac_file(filename, npoints_in, nstates, X_dvr, npesfiles, nac)
  use m_precision, only: wp
  use m_constants, only: const
  !use m_KH_functions, only: check_dvr_bounds 
  use m_splines, only: spline_easy

  character(80), intent(in):: filename
  integer, intent(in):: npoints_in, nstates
  real(kind=wp), intent(inout):: X_dvr(nstates)
  integer, intent(in):: npesfiles
  real(kind=wp), intent(out):: nac(npesfiles,npesfiles,nstates, 2)  

  real(kind=wp):: nac_in(npesfiles,npesfiles,npoints_in, 2), x_in(npoints_in)  
  integer:: j,k,l 

  open(10,file=filename ,status='unknown')

  do j=1,npesfiles
     do k=1,npesfiles
        do l=1,npoints_in
           ! read < j(R) | d/dR | k(R) > and < j(R) | d2/dR2 | k(R) >    
           read(10,*) x_in(l), nac_in(j,k,l,1), nac_in(j,k,l,2)
        end do

        call spline_easy(x_in, nac_in(j,k,:,1), npoints_in, X_dvr, nac(j,k,:,1), nstates)
        call spline_easy(x_in, nac_in(j,k,:,2), npoints_in, X_dvr, nac(j,k,:,2), nstates)
     end do
  end do

end subroutine read_nac_file

subroutine get_projections(p)
  use m_sckh_params_t, only: sckh_params_t 
  use m_precision, only: wp
  use m_io, only: get_free_handle
  
  type(sckh_params_t), intent(inout):: p 

  integer:: ifile,i
  real(8):: dnrm2
  
  if(p % use_proj) then
    
    ifile = get_free_handle()
    open(ifile, file= p % proj_file, action='read')    
    read(ifile,*) p % nproj
    
    if(allocated(p % projvec)) deallocate(p % projvec)
    allocate(p % projvec(p % nproj,3))
    
    do i=1, p % nproj                                         
      read(ifile,*) p % projvec(i,1), p % projvec(i,2), p % projvec(i,3) 
      
      !normalize projvec                               
      p % projvec(i,:) = p % projvec(i,:) / dnrm2(3,p % projvec(i,:),1)
      write(6,*) "projvector", i,  p % projvec(i,:)            
    end do
    
    close(ifile)
  else
    ! the three cartesian directions
    p % nproj =3
    
    if(allocated(p % projvec)) deallocate(p % projvec)
    allocate(p % projvec(p % nproj,3))
    
    p % projvec =0.0_wp
    p % projvec(1,1) =1.0_wp
    p % projvec(2,2) =1.0_wp
    p % projvec(3,3) =1.0_wp
    
  end if
  
end subroutine get_projections

subroutine read_file_list(filename, nfiles, file_list)
  use m_io, only: get_free_handle
  
  character(*), intent(in):: filename
  integer, intent(in):: nfiles
  character(*), intent(out), allocatable:: file_list(:)

  integer:: ifile, i
  
  ifile = get_free_handle()
  open(ifile, file= filename, action='read')
  
  if(allocated(file_list)) deallocate(file_list)
  allocate(file_list(nfiles))

  do i=1, nfiles
    read(ifile,*) file_list(i)
    write(6,*) file_list(i)
  end do
  
  close(ifile)

end subroutine read_file_list

end module m_PES_io


