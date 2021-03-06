module m_sckh_params_t
  use m_precision, only: wp
  implicit none


  ! this typs should contain all input parameters in the program
  type sckh_params_t

     ! runmode
     character(80):: runmode

     ! PES and dipole files (PES modes)
     character(80):: pes_file_i
     character(80):: pes_file_n
     character(80):: pes_file_dyn
     character(80):: pes_file_dyn2
     character(80):: pes_file_lp_corr
     character(80):: pes_file_ch_corr

     character(80):: pes_file_list_n
     character(80):: pes_file_list_f
     character(80):: pes_file_list_fn_corr
     character(80):: dipole_file_list_n
     character(80):: dipole_file_list_f

     character(80):: traj_file_list

     character(80), allocatable:: pes_files_n(:)
     character(80), allocatable:: pes_files_f(:)
     character(80), allocatable:: pes_files_fn_corr(:)
     character(80), allocatable:: dipolefile_n(:)
     character(80), allocatable:: dipolefile_f(:)

     !character(80), allocatable:: traj_files(:)

     character(80)::outfile

     character(80):: nac_file

     character(80):: proj_file

     ! 
     character(80):: PES_units
     
     ! 
     integer:: nstates
     integer:: npoints_in

     ! omega_in
     integer:: n_omega_in
     real(wp):: omega_in_start
     real(wp):: omega_in_end

     ! omega_out
     integer:: n_omega_out
     real(wp):: omega_out_start
     real(wp):: omega_out_end

     integer:: npesfile_n
     integer:: npesfile_f 
     integer:: shift_PES
     integer:: nonadiabatic
     real(wp):: mu
     real(wp):: dvr_start_in
     real(wp):: dx_in
     real(wp):: gamma_FWHM
     real(wp):: gamma_instr_FWHM
     real(wp):: gamma_inc_FWHM
     real(wp):: gamma_R_FWHM 

     character(80):: broadening_func_inc
     character(80):: broadening_func_instr
     
     character(80):: vib_solver
     character(80):: KH_amplitude_mode
     character(80):: KH_states_mode
     logical:: KH_bin_mode
     logical:: KH_nonres_mode
     
     logical:: use_n0_state
     
     ! sckh parameters
     ! for SCXAS
     character(80):: runmode_scxas
     logical:: include_ZPE
     character(80):: dipole_mode
     
     ! for SCKH_PES
     integer:: samplemode
     integer:: npoints_x_sampl 
     integer:: npoints_mom_sampl 
     integer:: ntsteps
     integer:: ntsteps2
     real(kind=wp):: delta_t
     real(kind=wp):: delta_t2
     logical:: use_dynamics_file
     character(80):: runmode_sckh_res
     logical:: sckh_dyn_separate
     character(80):: sckh_pes_dyn_mode
     character(80):: sckh_alt_mode

     logical:: use_E_mean
     real(wp):: E_nf_mean
     real(wp):: E_ni_mean

     logical:: use_abs_bc
     real(wp):: abs_bc_max_x
     
     ! for SCKH
     integer:: ntraj
     integer:: nfinal
     ! add by O.Takahashi 2018/07/01
     integer:: ninter
     character(80):: runmode_sckh
     logical:: use_mass_scaling
     real(wp):: mass_scaling
     
     ! projections
     logical:: use_proj
     integer:: nproj
     real(kind=wp), allocatable:: projvec(:,:)

     ! printing
     integer:: kh_print_stride

     ! add by O.Takahashi 2018/06/29
     character(80) :: output
     character(80) :: read_traj
     integer:: dyn_state_on_n
     
  end type sckh_params_t

contains

  subroutine init_sckh_params_t(p)
    use m_precision, only: wp
    use m_log, only : ilog
    use m_input, only : init_parameter, input_t, input_calc_params
    use m_upper, only : upper

    type(sckh_params_t), intent(out) :: p
    
    type(input_t) :: inp
    integer :: iv_in, iv
    
    iv_in = 1
    iv = 1 
    
    ! read input file
    call input_calc_params('input.txt', inp, iv_in, ilog)
    
    ! runmode: either "SCKH" or "KH"
    call init_parameter('runmode',inp, "KH", p % runmode ,iv)
    
    ! pes_file_i: PES file for the initial state
    call init_parameter('pes_file_i',inp, "pes_file_i.txt", p % pes_file_i  ,iv)
    ! pes_file_n: PES file for the intermediate state, XES
    call init_parameter('pes_file_n',inp, "pes_file_n.txt", p % pes_file_n  ,iv)
    ! pes_file_n: PES file for the dynamics, usually equal to pes_file_n, XES, or pes_file_i for RIXS
    call init_parameter('pes_file_dyn',inp, "pes_file_n.txt", p % pes_file_dyn  ,iv)
    ! pes_file_n: PES file for the second dynamics in RIXS, usually equal to pes_file_n
    call init_parameter('pes_file_dyn2',inp, "pes_file_n.txt", p % pes_file_dyn2  ,iv)
    ! pes_file_lp_corr: PES for first final state, computed more accurately
    call init_parameter('pes_file_lp_corr',inp, "pes_file_lp_corr.txt", p % pes_file_lp_corr ,iv)
    ! pes_file_ch_corr: PES for the core hole state, computed more accurately
    call init_parameter('pes_file_ch_corr',inp, "pes_file_ch_corr.txt", p % pes_file_ch_corr ,iv)
    
    ! pes_file_list_n: in case of many interemedite PES:es, list of these files
    call init_parameter('pes_file_list_n',inp, "pes_file_list_n.txt", p % pes_file_list_n ,iv)
    ! pes_file_list_n: list of the final PES:es 
    call init_parameter('pes_file_list_f',inp, "pes_file_list_f.txt", p % pes_file_list_f ,iv)
    ! pes_file_list_fn_corr: correction for final state PESes in the KH_resonant_orb program
    call init_parameter('pes_file_list_fn_corr',inp, "pes_file_list_fn_corr.txt", p % pes_file_list_fn_corr ,iv)
    ! dipole files for intermediate states
    call init_parameter('dipole_file_list_n',inp, "dipole_file_list_n.txt", p % dipole_file_list_n ,iv)
    ! dipole files for final states
    call init_parameter('dipole_file_list_f',inp, "dipole_file_list_f.txt", p % dipole_file_list_f ,iv)

    ! list of trajectroy files for SCKH
    call init_parameter('traj_file_list',inp, "traj_file_list.txt", p % traj_file_list ,iv)
    
    ! outfile: the base name of the output files
    call init_parameter('outfile',inp, "outfile,txt", p % outfile ,iv)
    
    ! nac_file: the file containing the non-adiabatic couplings
    call init_parameter('nac_file',inp, "nac_file.txt", p % nac_file ,iv)

    ! proj_file: the file containing the projections, if used
    call init_parameter('proj_file',inp, "proj_file.txt", p % proj_file ,iv)
    
    ! PES_units, for reading in PES and dipole files. "ANG" or "BOHR"
    call init_parameter('PES_units',inp, "ANG", p % PES_units ,iv)
    
    ! nstates: the number of dvr points, equally the number of vibrational eigenstates
    call init_parameter('nstates',inp, 100, p % nstates ,iv)
    ! npoints_in: the number of points in the supplied PES files
    call init_parameter('npoints_in',inp, 100, p % npoints_in ,iv)

    ! n_omega_in: number of absorption frequencies
    call init_parameter('n_omega_in',inp, 100, p % n_omega_in ,iv)
    ! omega_in_start: starting absorption frequency [eV]
    call init_parameter('omega_in_start',inp, 0.0_wp, p % omega_in_start ,iv)  
    ! omega_in_end: ending absorption frequency [eV]
    call init_parameter('omega_in_end',inp, 100.0_wp, p % omega_in_end ,iv)

    
    ! n_omega_out: number of emission frequencies
    call init_parameter('n_omega_out',inp, 100, p % n_omega_out ,iv)
    ! omega_start: starting emission frequency [eV]
    call init_parameter('omega_out_start',inp, 0.0_wp, p % omega_out_start ,iv)  
    ! omega_end: ending emission frequency [eV]
    call init_parameter('omega_out_end',inp, 100.0_wp, p % omega_out_end ,iv)


    ! npesfile_n: the number of intermediate state PES files
    call init_parameter('npesfile_n',inp, 1, p % npesfile_n ,iv)
    ! npesfile_f: the number of final state PES files
    call init_parameter('npesfile_f',inp, 1, p % npesfile_f ,iv)
    ! shift_PES: =1, then use pes_file_lp_correct to correct the first final state, =0, do nothing
    call init_parameter('shift_PES',inp, 0, p % shift_PES,iv)
    ! nonadiabatic: if 1 yes, if 0 no
    call init_parameter('nonadiabatic',inp, 0, p % nonadiabatic ,iv)
    ! mu, reduced mass
    call init_parameter('mu',inp, 1.0_wp, p % mu ,iv)
    ! dvr_start_in: x coordinate of the starting point of the dvr points [Angstrom]
    call init_parameter('dvr_start_in',inp, 0.0_wp, p % dvr_start_in ,iv)
    ! dx_in: spacing between dvr points [Angstroms]
    call init_parameter('dx_in',inp, 0.1_wp, p % dx_in ,iv)
    ! gamma_FWHM: lifetime broadening [eV]
    call init_parameter('gamma_FWHM',inp, 0.1_wp, p % gamma_FWHM ,iv)
    ! gamma_instr_FWHM: instrumental broadening [eV]
    call init_parameter('gamma_instr_FWHM',inp, 0.1_wp, p % gamma_instr_FWHM ,iv)
    ! gamma_inc_FWHM: incoming distribution broadening [eV]
    call init_parameter('gamma_inc_FWHM',inp, 0.1_wp, p % gamma_inc_FWHM ,iv)
    ! gamma_R_FWHM: broadening for resonance factor [eV]
    call init_parameter('gamma_R_FWHM',inp, 0.1_wp, p % gamma_R_FWHM ,iv)
    ! broadening function incoming  [eV]
    call init_parameter('broadening_func_inc',inp, "GAUSSIAN", p % broadening_func_inc ,iv)
    ! broadening function incoming  [eV]
    call init_parameter('broadening_func_instr',inp, "GAUSSIAN", p % broadening_func_instr ,iv)

    ! vib_solver: 1d vibrational solver, SINC_DVR or FOURIER_REAL
    call init_parameter('vib_solver',inp, "SINC_DVR", p % vib_solver ,iv)
    ! KH_amplitude_mode: OUTGOING or INGOING
    call init_parameter('KH_amplitude_mode',inp, "OUTGOING", p % KH_amplitude_mode ,iv)
    ! KH_states_mode: ORBS or STATES
    call init_parameter('KH_states_mode',inp, "STATES", p % KH_states_mode ,iv)
    ! KH_bin_mode: true or false 
    call init_parameter('KH_bin_mode',inp, .false., p % KH_bin_mode ,iv)
    ! KH_nonres_mode: true or false. Include nonresonant term or not  
    call init_parameter('KH_nonres_mode',inp, .true., p % KH_nonres_mode ,iv)
    ! use_n0_state: if .true., read reference state for |n>  
    call init_parameter('use_n0_state',inp, .false.,p % use_n0_state,iv)
    
    ! sckh parameters
    call init_parameter('runmode_scxas',inp, "full", p % runmode_scxas,iv)

    call init_parameter('include_ZPE',inp, .false., p % include_ZPE,iv)
    
    call init_parameter('dipole_mode',inp, "DIPOLE", p % dipole_mode,iv)
    
    call init_parameter('samplemode',inp, 1, p % samplemode,iv)

    call init_parameter('npoints_x_sampl',inp, 1, p % npoints_x_sampl,iv)

    call init_parameter('npoints_mom_sampl',inp,1 , p % npoints_mom_sampl,iv)

    call init_parameter('ntsteps',inp, 1, p % ntsteps,iv)

    call init_parameter('ntsteps2',inp, 1, p % ntsteps2,iv)

    call init_parameter('delta_t',inp, 0.1_wp, p % delta_t,iv)

    call init_parameter('delta_t2',inp, 0.1_wp, p % delta_t2,iv)

    call init_parameter('ntraj',inp, 1, p % ntraj, iv)

    call init_parameter('nfinal',inp, 1, p % nfinal, iv)

    ! add by O.Takahashi 2018/07/01
    call init_parameter('ninter',inp, 1, p % ninter, iv)

    call init_parameter('runmode_sckh',inp,'nonresonant',p%runmode_sckh,iv)

    call init_parameter('use_mass_scaling',inp, .false.,p % use_mass_scaling,iv)

    call init_parameter('mass_scaling',inp, 1.0d0 ,p % mass_scaling,iv)

    call init_parameter('runmode_sckh_res',inp,'FULL',p%runmode_sckh_res,iv)

    ! if .true. run D_in on E_dyn and D_nf on E_dyn2
    ! if .false, run t on E_dyn and t' on E_dyn2
    call init_parameter('sckh_dyn_separate',inp,.false.,p%sckh_dyn_separate,iv)

    ! option in SCKH_resonant_PES_FC 
    ! if "SINGLE", run dynamics on E_dyn2_inp
    ! if "SEPARATE", run dynamics on each intermediate state
    call init_parameter('sckh_pes_dyn_mode',inp, "SINGLE" ,p % sckh_pes_dyn_mode,iv)

    call init_parameter('sckh_alt_mode',inp, "NORMAL" ,p % sckh_alt_mode,iv)
    
    ! if use_E_mean= .true. then we use the input values, instead of computing them on the fly.
    ! E_nf_mean and E_ni_mean in eV
    call init_parameter('use_E_mean',inp, .false. ,p % use_E_mean,iv)
    call init_parameter('E_nf_mean',inp, 0.0_wp ,p % E_nf_mean,iv)
    call init_parameter('E_ni_mean',inp, 0.0_wp ,p % E_ni_mean,iv)


    ! Absorbing boundary conditions in verlet_trajectory_xva, if use_abs_bc = .true,
    ! then at potential will be absorbing at the value abs_bc_max_x [Angstrom] 
    call init_parameter('use_abs_bc',inp, .false. ,p % use_abs_bc,iv)
    call init_parameter('abs_bc_max_x',inp, 0.0_wp ,p % abs_bc_max_x,iv)
    
    ! stride in frequencies when printing KH cross sections
    call init_parameter('kh_print_stride',inp,1,p % kh_print_stride,iv)

    ! projections
    call init_parameter('use_proj',inp, .false., p % use_proj,iv)

    ! flag to use an additional dynamics file (pes_file_dyn)
    call init_parameter('use_dynamics_file',inp, .false., p % use_dynamics_file,iv)
    
    ! add by O.Takahashi 2018/06/29
    ! output things, trajectory, etc
    call init_parameter('output', inp, "NONE", p % output,iv)

    ! read trajectory from file, DEMON or PES
    call init_parameter('read_traj', inp, "DEMON", p % read_traj, iv)

    ! save dynamics state on intermediate states
    call init_parameter('dyn_state_on_n', inp, 1, p % dyn_state_on_n, iv)

  end subroutine init_sckh_params_t
  
end module m_sckh_params_t
