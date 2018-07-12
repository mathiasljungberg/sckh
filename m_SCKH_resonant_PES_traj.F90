module m_SCKH_resonant_PES_traj

  implicit none

contains

  subroutine calculate_SCKH_res_PES_traj(p)
    use m_precision, only: wp
    use m_constants, only: const
    use m_SCKH_utils, only: sample_x_mom_modes
    use m_SCKH_utils, only: verlet_trajectory_xva
    use m_sckh_params_t, only: sckh_params_t 
    use m_io, only: get_free_handle
    use m_splines, only: spline_easy
    use m_PES_io, only: read_dipole_file
    use m_PES_io, only: read_PES_file
    use m_PES_io, only: get_projections
    use m_PES_io, only: read_file_list
    use m_KH_functions, only: solve_vib_problem
    use m_upper, only : upper
!    use m_SCKH_resonant_PES, only: compute_E_means_and_omegas_one
    
    type(sckh_params_t), intent(inout):: p 

    real(wp), allocatable:: time(:)
    real(wp), allocatable:: E_dyn_inp(:)
    real(wp), allocatable:: E_dyn2_inp(:)  
    real(wp), allocatable:: E_n0_inp(:)
!    real(wp), allocatable:: E_n0(:)
    real(wp), allocatable:: E_corr(:,:)
    real(wp), allocatable:: E_fn_corr(:,:)
    real(wp), allocatable:: E_lp_corr(:)
    real(wp), allocatable:: shift(:)
    real(wp), allocatable:: E_i_inp(:)
    real(wp), allocatable:: E_i1(:)
!    real(wp), allocatable:: E_i2(:)
    real(wp), allocatable:: E_n_inp(:,:)
    real(wp), allocatable:: E_n1(:,:)
!    real(wp), allocatable:: E_n2(:,:)
    real(wp), allocatable:: E_f(:,:)
    real(wp), allocatable:: E_f_inp(:,:)
    real(wp), allocatable:: E_f1(:,:)
!    real(wp), allocatable:: E_f2(:,:)
    real(wp), allocatable:: E_fc1(:,:,:)
!    real(wp), allocatable:: E_fc2(:,:,:)
    real(wp), allocatable:: D_fn_inp(:,:,:,:)
    real(wp), allocatable:: D_fn1(:,:,:,:)
!    real(wp), allocatable:: D_fn2(:,:,:,:)
    real(wp), allocatable:: D_ni_inp(:,:,:)
    real(wp), allocatable:: D_ni1(:,:,:)
!    real(wp), allocatable:: D_ni2(:,:,:) 
!    real(wp):: E_nf_mean
!    real(wp):: E_fi_mean
!    real(wp):: E_ni_mean
    
    character(80):: string
    character(80):: file
    real(wp), allocatable:: x_sampl(:)
    real(wp), allocatable:: mom_sampl(:)
    real(wp), allocatable:: x_mom_sampl(:,:)
    real(wp), allocatable:: x_new(:)
    real(wp), allocatable:: v_new(:)
    real(wp), allocatable:: a_new(:)
    
    integer:: npoints_x_sampl
    integer:: npoints_mom_sampl
    integer:: npoints_x_mom_sampl
    integer:: ninter
    integer:: nfinal
    integer:: ntsteps  
    integer:: npoints_in
    integer:: nfinal_tot

    real(wp), allocatable:: c_i(:,:)
    real(wp), allocatable:: eig_i(:)
    
    real(wp):: time_l
    real(wp):: delta_t
    real(wp):: norm
    real(wp):: mu_SI
    real(wp):: dx
    real(wp):: fac_t
    
    integer, dimension(1):: ind  
    integer:: ifile
    
    real(wp), allocatable:: X_r(:)
    real(wp):: dvr_start 
    integer:: npoints, ii
    real(8):: dnrm2
    integer:: i
    integer:: j
    integer:: m
    integer:: m1
    integer:: m2
    integer:: m3
    integer:: m4
    integer:: traj
    integer:: traj2
    integer:: n_e
    integer:: f_e
    integer:: fn_e

    real(wp):: E_ni
    real(wp):: E_nf
    real(wp):: E_i_min, E_i_av
    

    ! set some local variables
    ntsteps = p % ntsteps
    npoints_x_sampl = p % npoints_x_sampl
    npoints_mom_sampl = p % npoints_mom_sampl
    npoints_in = p % npoints_in ! this should be renamed, corresponds to p % nstates
    nfinal = p % npesfile_f
    ninter = p % npesfile_n
    mu_SI = p % mu * const % u

    dvr_start = p % dvr_start_in * 1.0d-10
    dx = p % dx_in * 1.0d-10
    delta_t = p % delta_t

    ! projections
    call  get_projections(p)

    allocate(E_i_inp(npoints_in))
    allocate(E_n_inp(ninter,npoints_in))
    allocate(E_dyn_inp(npoints_in))
    allocate(E_dyn2_inp(npoints_in))
    allocate(E_f_inp(nfinal,npoints_in))
    allocate(E_n0_inp(npoints_in))
    allocate(E_lp_corr(npoints_in))
    allocate(D_ni_inp(ninter,npoints_in,3))
    allocate(D_fn_inp(nfinal,ninter,npoints_in,3))
    allocate(time(ntsteps))
!    allocate(E_n0(ntsteps))
    allocate(E_i1(ntsteps))
!    allocate(E_i2(ntsteps))
    allocate(E_n1(ninter,ntsteps))
!    allocate(E_n2(ninter, ntsteps))
    allocate(D_ni1(ninter,ntsteps,3))
!    allocate(D_ni2(ninter, ntsteps,3))
!    allocate(D_fn2(nfinal, ninter, ntsteps,3))
    allocate(D_fn1(nfinal,ninter,ntsteps,3))
    !allocate(E_i_mean(ntsteps))
    allocate(shift(npoints_in))
    allocate(c_i(npoints_in,npoints_in))
    allocate(eig_i(npoints_in))
    allocate(x_new(ntsteps))
    allocate(v_new(ntsteps))
    allocate(a_new(ntsteps))
    allocate(X_r(npoints_in))
    
    ! in ORB mode there are nfinal * ninter final states
    write(6,*) "p % KH_states_mode = ", p % KH_states_mode
    write(6,*) "p % sckh_pes_dyn_mode  = ", p % sckh_pes_dyn_mode
    write(6,*) "p % sckh_alt_mode  = ", p % sckh_alt_mode
    if (upper(p % KH_states_mode) .eq. "STATES") then
      nfinal_tot = nfinal

      allocate(E_fc1(nfinal,1, ntsteps))
!      allocate(E_fc1(nfinal,1, ntsteps), &
!           E_fc2(nfinal,1,ntsteps))
      
    else if (upper(p % KH_states_mode) .eq. "ORBS") then
      nfinal_tot = nfinal * ninter

      allocate(E_fc1(nfinal,ninter, ntsteps))
      allocate(E_corr(ninter, ntsteps))
      allocate(E_f(nfinal, ntsteps))
!      allocate(E_fc1(nfinal,ninter, ntsteps), &
!           E_fc2(nfinal,ninter, ntsteps))

    else
      write(6,*) "p % KH_states_mode should be either 'STATES' or 'ORBS'"
    end if

    
    ! set up grid points
    do i = 1, npoints_in
      X_r(i) = (i-1) * dx + dvr_start
    end do
    
    ! read PES files
    call read_PES_file(p % pes_file_i, p % npoints_in, p % npoints_in, X_r, E_i_inp)

    ! intermediate state reference energy (lowest state) 
    if( p % use_n0_state ) then 
      call read_PES_file(p % pes_file_n, p % npoints_in, p % npoints_in, X_r, E_n0_inp)
    else
      E_n0_inp = 0.0d0
    end if

    ! read list of intermediate state pes_files and dipole_files
    call read_file_list(p % pes_file_list_n, p % npesfile_n, p % pes_files_n)
    call read_file_list(p % dipole_file_list_n, p % npesfile_n, p % dipolefile_n)

    do j = 1 , p % npesfile_n
      call read_PES_file(p % pes_files_n(j), p % npoints_in, p % npoints_in, X_r, E_n_inp(j,:))
      call read_dipole_file(p % dipolefile_n(j), p % npoints_in, p % npoints_in, X_r, &
           D_ni_inp(j,:,:))
    end do
    
    ! read PES file where the first and second dynamics are run 
    if ( p % use_dynamics_file ) then
      call read_PES_file(p % pes_file_dyn, p % npoints_in, p % npoints_in, X_r, E_dyn_inp)
      call read_PES_file(p % pes_file_dyn2, p % npoints_in, p % npoints_in, X_r, E_dyn2_inp)
    else
      E_dyn_inp = E_i_inp
      E_dyn2_inp = E_n_inp(1,:) + E_n0_inp ! add reference state to dynamics 
    end if

    ! read list of final state pes_files and dipole_files
    call read_file_list(p % pes_file_list_f, p % npesfile_f, p % pes_files_f)
    call read_file_list(p % dipole_file_list_f, p % npesfile_f, p % dipolefile_f)

    do j = 1 , p % npesfile_f
      call read_PES_file(p % pes_files_f(j), p % npoints_in, p % npoints_in, X_r, E_f_inp(j,:))
      ! temporary hack: only one intermediate state
      call read_dipole_file(p % dipolefile_f(j), p % npoints_in, p % npoints_in, X_r, &
           D_fn_inp(j,1,:,:))
    end do

    ! ugly hack again to play
    !E_dyn_inp = (E_i_inp +  E_f_inp(1,:))/2.0d0
    !E_dyn2_inp = (E_i_inp + E_n_inp(1,:) +  E_f_inp(1,:))/3.0d0

    ! read list of corrections to the final state files coming from the excited electron
    if ( upper(p % KH_states_mode) .eq. "ORBS" ) then    
      allocate(E_fn_corr(p % npesfile_n, p % npoints_in))

      call read_file_list(p % pes_file_list_fn_corr, p % npesfile_n, p % pes_files_fn_corr)
      
      do j = 1 , p % npesfile_n
        call read_PES_file(p % pes_files_fn_corr(j), p % npoints_in, &
             p % npoints_in, X_r, E_fn_corr(j,:))
      end do
!    else
      
    end if

    ! Shift orbital energies so that E_f(1,:) have energies E_lp_corr
    ! and the spacing between the intermediate and final states are preserved
    if( p % shift_PES .eq. 1 ) then
      call read_PES_file(p % pes_file_lp_corr, p % npoints_in, p % npoints_in, X_r, &
           E_lp_corr)

      shift = E_lp_corr - E_f_inp(1,:) 
      
      do j = 1 , p % npesfile_f
        E_f_inp(j,:) = E_f_inp(j,:) + shift
      end do
      write(6,*) "Shifted PES:s"
    end if

    ! Solve the vibrational problem for initial state to be able to sample initial distribution
    call solve_vib_problem(dx, E_i_inp, eig_i, c_i, mu_SI, p % vib_solver)
    write(6,*) "Calculated initial state eigenfunctions"
    write(6,*) "Initial state fundamental", (eig_i(2) -eig_i(1))*const % cm

    ! convert to eV units
    E_i_inp = E_i_inp / const % eV
    E_n_inp = E_n_inp / const % eV
    E_n0_inp = E_n0_inp / const % eV
    E_dyn_inp = E_dyn_inp / const % eV
    E_dyn2_inp = E_dyn2_inp / const % eV
    E_f_inp = E_f_inp / const % eV

    ind = minloc(E_i_inp)
    E_i_min = E_i_inp(ind(1))
    E_i_av = sum(E_i_inp * c_i(:,1)**2) - E_i_min! ZPE

    write(6,*) "Average potenital energy, from minimum", E_i_av
    
    if(allocated(E_fn_corr)) E_fn_corr = E_fn_corr / const % eV

    ifile = get_free_handle()
    open(ifile, file="inital_state_eigvec.txt", action='write')
    do i = 1 , npoints_in
      write(ifile,'(3ES16.6)') X_r(i), c_i(i,1), c_i(i,1) ** 2
    end do
    close(ifile)

    ! sample the positions and momenta
    call sample_x_mom_modes(npoints_x_sampl, npoints_mom_sampl, &
         p % samplemode, X_r, c_i(:,1), x_mom_sampl)
    npoints_x_mom_sampl = size(x_mom_sampl,1)

    delta_t = delta_t * 1.d-15 ! femtoseconds
    time_l = (ntsteps-1) * delta_t

    write(6,*) "outfile", p % outfile
    write(6,*) "mu_SI", mu_SI 
    write(6,*) "time_l", time_l
    write(6,*)
    write(6,*) "Fundamental frequency resolution", 2.0_wp * const % pi * const % hbar &
         / ( time_l * const % eV )
    write(6,*) "delta t", delta_t

    !
    ! Here starts the program proper
    ! 

    do i = 1 , ntsteps
      time(i) = (i-1) * delta_t
    end do
    
!    call compute_E_means_and_omegas_one(E_i_inp, E_i_inp, E_n_inp(1,:)+ E_n0, &
!         E_f_inp(1,:), &
!         E_ni_mean, E_fi_mean, E_nf_mean, time_l, omega_in, omega_out)    

!    if (p % use_E_mean) then
!      write(6,*) "p % use_E_mean", p % use_E_mean
!      E_nf_mean = p % E_nf_mean
!      E_ni_mean = p % E_ni_mean
!    else
!      write(6,*) "p % use_E_mean", p % use_E_mean
!      E_nf_mean =   E_n_inp(1,ind(1))+ E_n0(ind(1)) - E_f_inp(1,ind(1)) !E_n1(1,ind(1)) - E_fc1(1,1,ind(1))
!      E_ni_mean =   E_n_inp(1,ind(1))+ E_n0(ind(1)) - E_i_inp(ind(1))   !E_n1(1,ind(1)) - E_i1(ind(1))
!    end if
!    E_fi_mean = E_ni_mean -E_nf_mean
!
!    write(6,*) "E_nf_mean", E_nf_mean
!    write(6,*) "E_ni_mean", E_ni_mean
!    write(6,*) "E_fi_mean", E_fi_mean
    
    do traj = 1 , npoints_x_mom_sampl

      write(6,*) "traj ", traj, "out of ", npoints_x_mom_sampl
      write(6,*) "x_ini: ", x_mom_sampl(traj,1), "v_ini: ", x_mom_sampl(traj,2)/mu_SI

      ! new loop here over intermediate states i order to do separate dynamics
      do n_e = 1 , ninter

        if( upper(p % sckh_pes_dyn_mode) .eq. "SEPARATE" ) then

          write(6,*) "dynamics mode: SEPARATE"
          ! compute mean energies on E_dyn_inp that will match with the other routines
          call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
               (E_n_inp(n_e,:) + E_n0_inp) * const % eV, delta_t, mu_SI, x_new, v_new, a_new, &
               p % use_abs_bc, p % abs_bc_max_x * 1.0d-10)
          
        else if ( upper(p % sckh_pes_dyn_mode) .eq. "SINGLE" ) then

          ! only do dynamics at first n
          if( n_e .eq. 1 ) then

            write(6,*) "dynamics mode: SINGLE"
            
            ! compute mean energies on E_dyn_inp that will match with the other routines
            call verlet_trajectory_xva(x_mom_sampl(traj,1), x_mom_sampl(traj,2)/mu_SI, X_r, &
                 E_dyn2_inp * const % eV, delta_t, mu_SI, x_new, v_new, a_new, &
                 p % use_abs_bc, p % abs_bc_max_x * 1.0d-10)
          end if

        else

          write(6,*) "p % sckh_pes_dyn_mode must be 'SINGLE' or 'SEPARATE'"
          
        end if ! if(upper(p % dynamics_mode .eq. "SEPARATE")) then

        ! a bit unnecessary splines but well...
        call spline_easy(X_r, E_i_inp, npoints_in, x_new, E_i1, ntsteps)
      
        !do n_e=1,ninter
        call spline_easy(X_r, E_n_inp(n_e,:) + E_n0_inp, npoints_in, x_new, E_n1(n_e,:), ntsteps)  
!        call spline_easy(X_r, E_n0_inp, npoints_in, x_new, E_n0, ntsteps)  
        !  end do
        
        ! options for dipole moments in XAS
        if ( upper(p % dipole_mode) .eq. "DIPOLE" ) then
          !do n_e=1,ninter
          do m = 1 , 3
            call spline_easy(X_r, D_ni_inp(n_e,:,m), npoints_in, x_new, D_ni1(n_e,:,m), ntsteps)  
            ! end do
          end do
        else if( upper(p % dipole_mode) .eq. "FC" ) then
          D_ni1 = 1.0_wp
        else if( upper(p % dipole_mode) .eq. "DIPOLE_X0" ) then
          do i = 1 , npoints_in !p % nstates
            D_ni1(:,i,:) = D_ni_inp(:,ind(1),:) 
          end do
        else
          write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
          stop
        end if

        ! options for dipole moments in XES
        if ( upper(p % dipole_mode) .eq. "DIPOLE" ) then        
          do f_e = 1 , nfinal
            do m = 1 , 3
              call spline_easy(X_r, D_fn_inp(f_e,1,:,m), npoints_in, x_new, D_fn1(f_e,1,:,m), &
                   ntsteps)  
            end do
          end do
        else if( upper(p % dipole_mode) .eq. "FC" ) then
          D_fn1 = 1.0_wp
        else if( upper(p % dipole_mode) .eq. "DIPOLE_X0" ) then
          do i = 1 , npoints_in !p % nstates
            D_fn1(:,1,i,:) = D_fn_inp(:,1,ind(1),:) 
          end do
        else
          write(6,*) "p % dipole_mode must be DIPOLE, FC or DIPOLE_X0"
          stop
        end if

      !F_if_om_omp = 0.0_wp
      
        if ( upper(p % KH_states_mode) .eq. "STATES" ) then
        
          do f_e = 1 , nfinal
            call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new, E_fc1(f_e,1,:), ntsteps)  
          end do
        
          do f_e = 1 , nfinal

            !  do n_e=1,ninter
            
            !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV

            if( p % include_ZPE ) then
              E_ni = E_n1(n_e,1) - (E_i_min + E_i_av)
              !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            else
              E_ni = E_n1(n_e,1) - E_i1(1) 
              !E_ni = E_n1(n_e,1)-E_i1(1) + (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            end if
            
            E_nf = E_n1(n_e,1) - E_fc1(f_e,1,1)

          end do

        else if ( upper(p % KH_states_mode) .eq. "ORBS" ) then

          call spline_easy(X_r, E_fn_corr(n_e,:), npoints_in, x_new, E_corr(n_e,:), ntsteps)

          do f_e = 1 , nfinal
            !do n_e=1,ninter

            fn_e = ninter * (f_e - 1) + n_e  

            write(6,*) "f_e, n_e ", f_e, n_e, " out of ", nfinal, ninter
            
            call spline_easy(X_r, E_f_inp(f_e,:), npoints_in, x_new, E_f(f_e,:), ntsteps)  
            call spline_easy(X_r, E_f_inp(f_e,:) + E_fn_corr(n_e,:), npoints_in, x_new, &
                 E_fc1(f_e,n_e,:), ntsteps)  
            
            if( p % include_ZPE ) then
              E_ni = E_n1(n_e,1) - (E_i_min + E_i_av)
            else
              E_ni = E_n1(n_e,1) - E_i1(1) !+ (x_mom_sampl(traj,2) **2 / (2.0_wp*mu_SI)) / const % eV
            end if

            E_nf = E_n1(n_e,1) - E_fc1(f_e,n_e,1)
            
            !end do
          end do ! do f_e=1,nfinal
        
        end if ! if (upper(p % KH_states_mode) .eq. "STATES") then
      
      end do ! do n_e=1,ninter

      ! add by O.Takahashi 2018/06/29
      ! write trajectory data to file
      if ( index(upper(p % output), "TRAJECTORY") .ne. 0 ) then

        write(6,*) "write trajectory data to file"
        file = "_traj_"
        write(string,'(I4)') traj
        file = trim(adjustl(p % outfile)) // trim(adjustl(file)) // &
             trim(adjustl(string)) // ".dat"
        write(6,*) file

        if (upper(p % KH_states_mode) .eq. "STATES") then
        
          call write_one_sckh_res_PES_traj(ntsteps, nfinal, ninter, file, &
               time, E_i1, E_n1, E_fc1(:,1,:), D_ni1(:,1,:), D_fn1(:,1,:,:))

        else if (upper(p % KH_states_mode) .eq. "ORBS") then

          call write_one_sckh_res_PES_traj_orbs(ntsteps, nfinal, ninter, file, &
               time, E_i1, E_n1, E_f, E_corr, D_ni1(:,1,:), D_fn1(:,1,:,:))

        end if

      end if

    end do ! do traj=1, npoints_x_mom_sampl
    
  end subroutine calculate_SCKH_res_PES_traj

  subroutine write_one_sckh_res_PES_traj(ntsteps_inp, nfinal, ninter, traj_file, &
       time_inp, E_i1, E_n1, E_fc1, D_ni1, D_fn1)
    use m_precision, only: wp
    use m_constants, only: const
    use m_io, only: get_free_handle

    integer, intent(in):: ntsteps_inp, nfinal, ninter
    character(*), intent(in):: traj_file
    real(wp), intent(in):: time_inp(:), E_i1(:), E_n1(:,:)
    real(wp), intent(in):: E_fc1(:,:), D_fn1(:,:,:)
    real(wp), intent(in):: D_ni1(:,:)

    integer:: ifile, i, j, m, ntrans !, jj
    character(80):: dummy



    ifile = get_free_handle()
    open(ifile,file=traj_file,status='unknown')

    ! XAS for first time step
    dummy = 'XAS '
    write(ifile,'(a5,i5)') dummy, ninter ! number of x-ray transitions, should be the same number as the number of unocc states used

    do j = 1 , ninter
!      write(ifile,'(4e20.10)') E_n1(j,1), D_ni1(j,1), D_ni1(j,2), D_ni1(j,3)
      write(ifile,'(4e20.10)') E_n1(j,1), (D_ni1(j,m), m=1,3)
    end do

    ! now read trajectory with XES and more stuff
    do i = 1 , ntsteps_inp

      write(ifile,'(e20.10)') time_inp(i) * 1.0D15
      write(ifile,'(e20.10)') E_i1(i)
      do j = 1 , ninter
        write(ifile,'(e20.10)') E_n1(j,i)
      end do

      ! XES

      dummy = 'XES '
      write(ifile,'(a5,i5)') dummy, nfinal

      do j = 1 , nfinal
!        jj = j + ntrans - nfinal 
!        write(ifile,*) E_f_inp(jj,i), D_fn1(jj,i,1), D_fn1(jj,i,2), D_fn1(jj,i,3)
!        write(ifile,'(f15.8,3e17.8)') E_fc1(j,i), D_fn1(j,i,1), D_fn1(j,i,2), D_fn1(j,i,3)
        write(ifile,'(f15.8,3e17.8)') E_fc1(j,i), (D_fn1(j,i,m), m=1,3)
      end do

    end do !i

    close(ifile)

  end subroutine write_one_sckh_res_PES_traj

  subroutine write_one_sckh_res_PES_traj_orbs(ntsteps_inp, nfinal, ninter, traj_file, &
       time_inp, E_i1, E_n1, E_f, E_fn_corr, D_ni1, D_fn1)

    use m_precision, only: wp
    use m_constants, only: const
    use m_io, only: get_free_handle

    integer, intent(in):: ntsteps_inp, nfinal, ninter
    character(*), intent(in):: traj_file
    real(wp), intent(in):: time_inp(:), E_i1(:)
    real(wp), intent(in):: E_n1(:,:)
    real(wp), intent(in):: E_f(:,:), E_fn_corr(:,:), D_fn1(:,:,:)
    real(wp), intent(in):: D_ni1(:,:)

    real(wp):: rdummy
    integer:: ifile, i, j, k, m, ntrans !, jj
    integer:: norbs_gs, nocc_gs, norbs_exc, nocc_exc
    character(80):: dummy
    character(80):: traj_file_demon


! write PES format

    ifile = get_free_handle()
    open(ifile,file=traj_file,status='unknown')

    ! XAS for first time step
    dummy = 'XAS'
    write(ifile,'(a5,i5)') dummy, ninter ! number of x-ray transitions, should be the same number as the number of unocc states used

    do j = 1 , ninter
!      write(ifile,'(4e20.10)') E_n1(j,1), D_ni1(j,1), D_ni1(j,2), D_ni1(j,3)
      write(ifile,'(4e20.10)') E_n1(j,1), (D_ni1(j,m), m=1,3)
    end do

    ! now read trajectory with XES and more stuff
    do i = 1 , ntsteps_inp

      write(ifile,'(e20.10)') time_inp(i) * 1.0D15
      write(ifile,'(e20.10)') E_i1(i)
      do k = 1 , ninter
        write(ifile,'(e20.10)') E_n1(k,i)
      end do

      ! XES

      dummy = 'XES'
      write(ifile,'(a5,i5)') dummy, nfinal

        do j = 1 , nfinal
!        jj = j + ntrans - nfinal 
!        write(ifile,*) E_f_inp(jj,i), D_fn1(jj,i,1), D_fn1(jj,i,2), D_fn1(jj,i,3)
!          write(ifile,'(4e17.8)') E_f(j,i), D_fn1(j,i,1), D_fn1(j,i,2), D_fn1(j,i,3)
          write(ifile,'(4e17.8)') E_f(j,i), (D_fn1(j,i,m), m=1,3)
        end do
        
      do k = 1 , ninter
        write(ifile,'(4e17.8)') E_fn_corr(k,i)
      end do

    end do !i

    close(ifile)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write deMon format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ifile = get_free_handle()
    traj_file_demon = trim(adjustl(traj_file)) // "demon"
    open(ifile,file=traj_file_demon,status='unknown')

    ! XAS for first time step
    dummy = 'XAS'
    write(ifile,'(e20.10)') rdummy ! E_XAS_inp ! total energy
    write(ifile,'(e20.10)') rdummy ! E_IP1s_XAS ! 1s orbital energy 
    write(ifile,'(a5,i5)') dummy, ninter ! number of x-ray transitions, should be the same number as the number of unocc states used

    do j = 1 , ninter
!      write(ifile,'(4e20.10)') E_n1(j,1), D_ni1(j,1), D_ni1(j,2), D_ni1(j,3)
      write(ifile,'(f15.10,3e20.10)') &
           ( E_n1(j,1) - E_i1(1) ) * const % eV / const % hartree, &
           (D_ni1(j,m), m=1,3)
    end do

    ! now read trajectory with XES and more stuff
    do i = 1 , ntsteps_inp

      write(ifile,'(f10.5)') time_inp(i) * 1.0D15  ! time_inp(i)
      write(ifile,'(e20.10)') E_i1(i) * const % eV / const % hartree   ! E_gs_inp(i)
      write(ifile,'(e20.10)') E_n1(1,i) * const % eV / const % hartree ! E_n0(i) = E_n1(1,i)
      write(ifile,'(e20.10)') E_n1(1,i) - E_i1(i)  ! E_IP1s(i)

      ! XES

      dummy = 'XES'
      write(ifile,'(a5,i5)') dummy, nfinal

      do j = 1 , nfinal
!        write(ifile,'(4e17.8)') E_f(j,i), D_fn1(j,i,1), D_fn1(j,i,2), D_fn1(j,i,3)
        write(ifile,'(f15.10,3e20.10)') &
             ( E_n1(1,i) - E_f(j,i) ) * const % eV / const % hartree, &
             (D_fn1(j,i,m), m=1,3)
      end do

      ! orbital energies for the ground state

      nocc_gs = 0
      norbs_gs = ninter + nocc_gs
      write(ifile,*) norbs_gs, nocc_gs
      ! check
      if ( norbs_gs .lt. nocc_gs ) then
        write(6,*) "Error, norbs_gs should be larger than nocc_gs", &
             norbs_gs, nocc_gs
        stop
      end if
      if ( norbs_gs .lt. nocc_gs+ninter-1 ) then
        write(6,*) "Error, norbs_gs should be larger than nocc_gs+ninter-1", &
             norbs_gs, nocc_gs, ninter 
        stop
      end if
      do j = nocc_gs + 1 , norbs_gs
        write(ifile,'(f20.10)') E_fn_corr(j,i) * const % eV / const % hartree
      end do

      ! orbital energies for the excited state
      nocc_exc = 1
      norbs_exc = ninter + nocc_exc - 1
      write(ifile,*) norbs_exc, nocc_exc
      ! check
      if ( norbs_exc .lt. nocc_exc ) then
        write(6,*) "Error, norbs_exc should be larger than nocc_exc", &
             norbs_exc, nocc_exc
        stop
      end if
      if ( norbs_exc .lt. nocc_exc+ninter-1 ) then
        write(6,*) "Error, norbs_exc should be larger than nocc_exc+ninter-1", &
             norbs_exc, nocc_exc, ninter 
        stop
      end if
      do j = nocc_exc , norbs_exc
        write(ifile,'(f20.10)') (E_n1(j,i) - E_n1(1,i)) * const % eV / const % hartree
      end do

     end do !i

    close(ifile)

  end subroutine write_one_sckh_res_PES_traj_orbs

  subroutine read_one_sckh_res_PES_traj(ntsteps_inp, nfinal, ninter, traj_file, &
       time_inp, E_i_inp, E_n_inp, E_f_inp, D_ni1, D_fn1)

    use m_precision, only: wp
    use m_constants, only: const
    use m_io, only: get_free_handle

    integer, intent(in):: ntsteps_inp, nfinal, ninter
    character(*), intent(in):: traj_file

    real(wp), intent(in):: time_inp(:)
    real(wp), intent(out):: E_i_inp(:), E_n_inp(:,:)
    real(wp), intent(out):: E_f_inp(:,:), D_fn1(:,:,:)
    real(wp), intent(out):: D_ni1(:,:)

    integer:: ifile, i, j, k, m, ntrans !, jj
    character(80):: dummy
    real(wp):: rdummy

!    write(6,*) 'filename: ', traj_file

    ifile = get_free_handle()
    open(ifile,file=traj_file,status='old')  

    ! XAS for first time step
    read(ifile,*) dummy, ntrans ! number of x-ray transitions, should be the same number as the number of unocc states used

    ! check that
    if ( ntrans .ne. ninter ) then
      write(6,*) "Error, ntrans != ninter", ntrans, ninter
    end if

    do j = 1 , ntrans
!      read(ifile,*) E_n_inp(j,1), D_ni1(j,1), D_ni1(j,2), D_ni1(j,3)
      read(ifile,*) E_n_inp(j,1), (D_ni1(j,m), m=1,3)
    end do

    ! now read trajectory with XES and more stuff
    do i = 1 , ntsteps_inp

      read(ifile,*) rdummy ! time_inp2(i)
      read(ifile,*) E_i_inp(i)
      do k = 1 , ninter
        read(ifile,*) E_n_inp(k,i)
      end do
      
      ! XES

      read(ifile,*) dummy, ntrans

      ! check that
      if ( ntrans .ne. nfinal ) then
        write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
      end if

      do j = 1 , ntrans
!        jj = j + ntrans - nfinal 
!        read(ifile,*) E_f_inp(jj,i), D_fn1(jj,i,1), D_fn1(jj,i,2), D_fn1(jj,i,3)
!        read(ifile,*) E_f_inp(j,i), D_fn1(j,i,1), D_fn1(j,i,2), D_fn1(j,i,3)
        read(ifile,*) E_f_inp(j,i), (D_fn1(j,i,m), m=1,3)
      end do

!      write(6,*) 'filename: ', traj_file,i

    end do !i

!    write(6,*) 'final filename: ', traj_file

    close(ifile)

  end subroutine read_one_sckh_res_PES_traj

  subroutine read_one_sckh_res_PES_traj_orbs(ntsteps_inp, nfinal, ninter, &
       traj_file, time_inp, E_i_inp, E_n_inp, E_f_inp, E_fn_corr, D_ni_inp, D_fn_inp)

    use m_precision, only: wp
    use m_constants, only: const
    use m_io, only: get_free_handle

    integer, intent(in):: ntsteps_inp, nfinal, ninter
    character(*), intent(in):: traj_file

    real(wp), intent(in):: time_inp(:)
    real(wp), intent(out):: E_i_inp(:), E_n_inp(:,:)
    real(wp), intent(out):: E_f_inp(:,:), E_fn_corr(:,:), D_fn_inp(:,:,:)
    real(wp), intent(out):: D_ni_inp(:,:)

    integer:: ifile, i, j, k, m, ntrans !, jj
    character(80):: dummy
    real(wp):: rdummy

!    write(6,*) 'filename: ', traj_file

    ifile = get_free_handle()
    open(ifile,file=traj_file,status='old')  

    ! XAS for first time step
    read(ifile,*) dummy, ntrans ! number of x-ray transitions, should be the same number as the number of unocc states used

    ! check that
    if ( ntrans .ne. ninter ) then
      write(6,*) "Error, ntrans != ninter", ntrans, ninter
    end if

    do j = 1 , ntrans
!      read(ifile,*) E_n_inp(j,1), D_ni_inp(j,1), D_ni_inp(j,2), D_ni_inp(j,3)
      read(ifile,*) E_n_inp(j,1), (D_ni_inp(j,m), m=1,3)
    end do

!    write(6,*) 'filename: ', traj_file

    ! now read trajectory with XES and more stuff
    do i = 1 , ntsteps_inp

      read(ifile,*) rdummy ! time_inp2(i)
      read(ifile,*) E_i_inp(i)
      do k = 1 , ninter
        read(ifile,*) E_n_inp(k,i)
      end do
      
      ! XES

      read(ifile,*) dummy, ntrans

      ! check that
      if ( ntrans .ne. nfinal ) then
        write(6,*) "Error, ntrans != nfinal", ntrans, nfinal
      end if

      do j = 1 , nfinal
!        jj = j + ntrans - nfinal 
!        read(ifile,*) E_f_inp_inp(jj,i), D_fn_inp(jj,i,1), D_fn_inp(jj,i,2), D_fn_inp(jj,i,3)
!        read(ifile,*) E_f_inp(j,i), D_fn_inp(j,i,1), D_fn_inp(j,i,2), D_fn_inp(j,i,3)
        read(ifile,*) E_f_inp(j,i), (D_fn_inp(j,i,m), m=1,3)
      end do
      do k = 1 , ninter
        read(ifile,*) E_fn_corr(k,i)
      end do

    end do !i

    close(ifile)

  end subroutine read_one_sckh_res_PES_traj_orbs

 
end module m_SCKH_resonant_PES_traj
