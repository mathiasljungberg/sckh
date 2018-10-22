module m_check_E_mean

  implicit none

contains

  subroutine check_e_mean_nonresonant(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_sckh_params_t, only: sckh_params_t
    use m_io, only: get_free_handle
    use m_SCKH_utils, only: read_one_sckh_traj

    type(sckh_params_t), intent(inout):: p 

    character(80), allocatable:: traj_files(:)
    character(80):: string
    character(80):: file

    real(kind=wp), allocatable:: time(:)
    real(kind=wp), allocatable:: time_inp(:)
    real(kind=wp), allocatable:: time_inp2(:)
    real(kind=wp), allocatable:: E_n_inp(:) 
    real(kind=wp), allocatable:: traj_weights(:) 
    real(kind=wp), allocatable:: E_f_inp(:,:)
    real(kind=wp), allocatable:: E_trans(:,:)
    real(kind=wp), allocatable:: D_fn_inp(:,:,:)
    real(kind=wp), allocatable:: E_IP1s(:) 
    real(kind=wp), allocatable:: E_gs_inp(:) 
    real(kind=wp):: time_l, E, delta_t, time_l2, E_nf_mean

    integer:: ntraj, ntsteps, ntsteps_inp, nfinal, ifile
    integer:: i, j, traj

    logical:: check_time
    
    ! set some local variables
    ntraj = p % ntraj
    ntsteps_inp = p % ntsteps
    ntsteps = p % ntsteps2
    nfinal = p % nfinal 
    delta_t = p % delta_t

    allocate(traj_files(ntraj))
    allocate(traj_weights(ntraj))

    ! read names of trajectory files
    ifile = get_free_handle()
    open(ifile,file= p % traj_file_list,status='unknown')

    do i = 1 , ntraj
      read(ifile,*) traj_files(i), traj_weights(i)
    end do

    close(ifile)

    allocate(time_inp(ntsteps_inp))
    allocate(time_inp2(ntsteps_inp))
    allocate(E_gs_inp(ntsteps_inp))
    allocate(E_n_inp(ntsteps_inp))
    allocate(E_f_inp(nfinal,ntsteps_inp))
    allocate(D_fn_inp(nfinal,ntsteps_inp,3))
    allocate(E_IP1s(ntsteps_inp))
    allocate(E_trans(nfinal,ntsteps_inp))

    ! create time_inp, in s
    delta_t = delta_t * 1.d-15 ! femtoseconds
    time_l = (ntsteps_inp - 1) * delta_t
    do i = 1 , ntsteps_inp
      time_inp(i)= (i-1) * delta_t
    end do

    delta_t = time(2) - time(1)
    time_l2 = (ntsteps-1) * delta_t 

    !write(6,*) "ntsteps", ntsteps
    !write(6,*) "delta_t", delta_t
    !write(6,*) "time_l", time_l

    !
    ! Loop over trajectories
    !

    E_nf_mean = 0.0_wp

    do traj = 1 , ntraj

      call read_one_sckh_traj(ntsteps_inp, nfinal, traj_files(traj), time_inp, &
           time_inp2, E_gs_inp, E_n_inp, E_IP1s, E_trans, E_f_inp, D_fn_inp, check_time)

      E_nf_mean = E_nf_mean + sum(E_n_inp(:) - E_f_inp(nfinal,:)) / max(1,size(E_n_inp(:)))

    ! write a list of energy 

      write(string,'(i5)') traj
      file = "energy_" // trim(adjustl(string)) // ".input"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i = 1 , ntsteps_inp
        write(ifile,'(f10.3,200('', '',e13.6))') &
             time_inp(i) * 1.d15, E_gs_inp(i), &
             E_n_inp(i) * (const % eV / const % Hartree), &
             (E_f_inp(j,i) * (const % eV / const % Hartree), j=1,nfinal)
      end do

      close(ifile)

    end do !end traj

    E_nf_mean = E_nf_mean / ntraj

    write(6,*) 'E_nf_mean = ', E_nf_mean

  end subroutine check_e_mean_nonresonant


  subroutine check_e_mean_resonant(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_sckh_params_t, only: sckh_params_t 
    use m_io, only: get_free_handle
    use m_SCKH_utils, only: read_one_sckh_res_traj
    
    type(sckh_params_t), intent(inout):: p 

    character(80),allocatable:: traj_files(:)
    character(80):: string
    character(80):: file

    integer, dimension(1)::ind  
    integer:: ifile, ntraj, ntsteps, ntsteps_inp, ninter, nfinal
    integer:: i, j, traj

    real(wp), allocatable:: time(:)
    real(wp), allocatable:: time_inp(:)
    real(wp), allocatable:: time_inp2(:)
    real(wp), allocatable:: E_n0(:)
    real(wp), allocatable:: E_fn_corr(:,:)
    real(wp), allocatable:: E_i_inp(:)
    real(wp), allocatable:: E_n_inp(:,:)
    real(wp), allocatable:: E_f_inp(:,:)
    real(wp), allocatable:: D_fn_inp(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:)
    real(wp), allocatable:: traj_weights(:)
    real(wp), allocatable:: E_IP1s(:)
    real(wp), allocatable:: E_trans(:,:)
    real(wp):: E_nf_mean, E_fi_mean, E_ni_mean
    !real(wp),allocatable:: E_IP1s(:)
    
    real(wp):: delta_t, time_l, time_l2
    
!    logical:: check_time
    
    ! set some local variables
    ntraj = p % ntraj
    ntsteps_inp = p % ntsteps
    ntsteps = p % ntsteps2
    nfinal = p % nfinal
    ninter = p % ninter
    delta_t = p % delta_t

    allocate(traj_files(ntraj), traj_weights(ntraj))
    
    ! read names of trajectory files
    ifile = get_free_handle()
    open(ifile, file= p % traj_file_list,status='unknown')
    
    do i = 1 , ntraj
      read(ifile,*) traj_files(i), traj_weights(i)
    end do
    
    close(ifile)

    ! allocate variables
    allocate(time_inp(ntsteps_inp))
    allocate(time_inp2(ntsteps_inp))

    allocate(E_i_inp(ntsteps_inp))
    allocate(E_n_inp(ninter,ntsteps_inp))
    allocate(E_f_inp(nfinal,ntsteps_inp))
    allocate(E_n0(ntsteps_inp))
    allocate(D_ni_inp(ninter, ntsteps_inp,3))
    allocate(D_fn_inp(nfinal, ninter, ntsteps_inp,3))
    allocate(E_IP1s(ntsteps_inp))
    allocate(E_trans(nfinal,ntsteps_inp))

    !
    ! Here starts the program proper
    ! 

    delta_t = delta_t * 1.d-15 ! femtoseconds
    time_l = (ntsteps_inp-1) * delta_t
    do i = 1 , ntsteps_inp
      time_inp(i) = (i-1) * delta_t
    end do

    delta_t = time(2) - time(1)
    time_l2 = (ntsteps-1) * delta_t 
    
    !write(6,*) "ntsteps", ntsteps
    !write(6,*) "delta_t", delta_t
    !write(6,*) "time_l", time_l

    E_nf_mean = 0.0_wp
    E_ni_mean = 0.0_wp

    do traj = 1 , ntraj  
    !
    ! read trajectory file
    !  
      call read_one_sckh_res_traj(ntsteps_inp, nfinal, ninter, &
           traj_files(traj), time_inp, &
           time_inp2, E_i_inp,  E_n_inp, &
           E_IP1s, E_trans, &
           E_f_inp, D_fn_inp(:,1,:,:), &
           !E_XAS_inp, E_IP1s_XAS, E_trans_XAS, &
           D_ni_inp(:,1,:),&
           E_n0, &   ! reference state
           E_fn_corr, &
           !norbs_gs, nocc_gs, eps_gs, norbs_exc, nocc_exc, eps_exc,&
           .false.)
      !check_time)

      E_nf_mean = E_nf_mean + sum(E_n_inp(1,:) - E_f_inp(nfinal,:)) / max(1,size(E_n_inp(1,:)))
      E_ni_mean = E_ni_mean + sum(E_n_inp(1,:) - E_i_inp(:)) / max(1,size(E_n_inp(1,:)))
      !write(6,*) "E_nf_mean", E_nf_mean
      !write(6,*) "E_ni_mean", E_ni_mean

      !ind = minloc(E_i_inp)
      !E_nf_mean = E_nf_mean + (E_n_inp(1,ind(1)) - E_f_inp(1,ind(1)))
      !E_ni_mean = E_ni_mean + (E_n_inp(1,ind(1)) - E_i_inp(ind(1)))

      !write(6,*) "E_nf_mean", E_nf_mean
      !write(6,*) "E_ni_mean", E_ni_mean
        
    ! write a list of energy 

      write(string,'(i5)') traj
      file = "energy_" // trim(adjustl(string)) // ".input"

      ifile = get_free_handle()
      open(ifile,file=file,status='unknown')

      do i = 1 , ntsteps_inp
        write(ifile,'(f10.3,200('', '',e13.6))') &
             time_inp(i) * 1.d15, E_i_inp(i), &
             (E_n_inp(j,i) * (const % eV / const % Hartree), j=1,ninter), &
             (E_f_inp(j,i) * (const % eV / const % Hartree), j=1,nfinal)
      end do

      close(ifile)

    end do ! do traj = 1, ntraj

    E_nf_mean = E_nf_mean / ntraj
    E_ni_mean = E_ni_mean / ntraj
    E_fi_mean = E_ni_mean - E_nf_mean
    
    write(6,*) 'E_nf_mean = ', E_nf_mean
    write(6,*) 'E_ni_mean = ', E_ni_mean
    write(6,*) 'E_fi_mean = ', E_fi_mean

  end subroutine check_e_mean_resonant



end module m_check_E_mean
