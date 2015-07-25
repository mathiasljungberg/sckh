module m_log_type

  !!
  !! Because of this use statements the subroutines die(...) and warn(...)
  !! could be "used" via the module m_log. 
  !!
  use m_die, only : die
  use m_warn, only : warn
  
  implicit none

#include "m_define_macro.F90"
  
  type log_t 
    integer :: iverb=0  ! verbosity level
    integer :: stdout=0  ! the file handle to write stdout messages

    integer :: times_handle = -999;
    character(1000) :: times_fname='';
   
    integer :: sizes_handle = -999;
    character(1000) :: sizes_fname='';

    character(1000) :: memory_fname=''
    integer :: mem_note_first = 1
    
  end type !log_t  
  
  
  contains

!
! Initialize file names and file units (numbers) for log- and times- files.
!
subroutine init_log(fname_suffix, mynode, iv, lb)
#define _sname "init_log"
  use m_io, only : get_free_handle
  use m_color, only : bc, color_on  
  implicit none
  character(len=*), intent(in) :: fname_suffix
  integer, intent(in) :: mynode, iv
  type(log_t), intent(inout) :: lb ! log_book
  
  ! internal
  integer :: ios

  lb%stdout = 6
  lb%iverb  = iv
  
  write(lb%times_fname, '(a,i0.4,a)') 'TIMES.', mynode, trim(fname_suffix);
  lb%times_handle=get_free_handle();
  open(lb%times_handle, file=lb%times_fname, action='write', iostat=ios);
  if(ios==0) then;
    if(iv>0) write(lb%stdout,*)_sname, trim(lb%times_fname), ' is open...';
  else
    write(0,*) _sname, trim(lb%times_fname), ' cannot be open.';
  endif

  write(lb%sizes_fname, '(a,i0.4,a)') 'SIZES.', mynode, trim(fname_suffix);
  lb%sizes_handle=get_free_handle();
  open(lb%sizes_handle, file=lb%sizes_fname, action='write', iostat=ios);
  if(ios==0) then;
    if(iv>0) write(lb%stdout,*)_sname, trim(lb%sizes_fname), ' is open...';
  else
    write(0,*) _sname, trim(lb%sizes_fname), ' cannot be open.';
  endif

  write(lb%memory_fname, '(a,i0.4,a)') 'MEMORY.', mynode, trim(fname_suffix);
  if(iv>0) write(lb%stdout,*)_sname, trim(lb%memory_fname), ' will be open...';

#undef _sname
end subroutine !init_log

end module !m_log_type
