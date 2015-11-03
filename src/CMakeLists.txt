cmake_minimum_required(VERSION 2.4)

PROJECT(fortran-utilities Fortran)


IF(CMAKE_Fortran_COMPILER MATCHES ifort)

  MESSAGE(STATUS  "ifort compiler found")

  SET(FLAGS "-O0 -warn all -g -CB")
  SET(INCDIR "-I/usr/include")
  #SET(OMP_FLAGS "")
  #SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")

   SET(MKLROOT "/opt/intel/composer_xe_2013_sp1/mkl")
   SET(LAPACK_LIB " -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm")
   IF(USER_OS STREQUAL LinuxMint )
   SET(MKLROOT "/opt/intel/composer_xe_2013_sp1/mkl")
   SET(INCDIR " -I${MKLROOT}/include")
  #SET(OMP_FLAGS "")
  #SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
   SET(LAPACK_LIB " -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm")
   ENDIF(USER_OS STREQUAL LinuxMint)

ELSE(CMAKE_Fortran_COMPILER MATCHES gfortran)

  MESSAGE(STATUS  "gfortran compiler found")

  IF(USER_OS STREQUAL mac)

    MESSAGE(STATUS  "option -DUSER_OS=mac, gfortran for mac")

    # gfortran/mac
    #SET(FLAGS "-O0  -fbounds-check -ffree-line-length-none -framework vecLib -DUSE_OWN_CZDOT")
    SET(FLAGS "-O0  -fbounds-check -framework vecLib")  
    SET(INCDIR "-I/usr/local/include -I/usr/include")
    #SET(OMP_FLAGS "")
    #SET(FFTW_LIB "-L/usr/local/lib -lfftw3 -lm")
    SET(LAPACK_LIB "")
  ELSE(USER_OS STREQUAL mac)

    MESSAGE(STATUS  "gfortran for Linux")

    # gfortran/linux
    SET(FLAGS "-O0  -fbounds-check -ffree-line-length-none -fmax-errors=5")
    SET(INCDIR "-I/usr/include")
    SET(OMP_FLAGS "")
    SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
    SET(MKLPATH "")
    SET(LAPACK_LIB "-L/usr/lib -lblas -llapack")

  ENDIF(USER_OS STREQUAL mac)

ENDIF(CMAKE_Fortran_COMPILER MATCHES ifort)


SET(ALL_MODULES_F90 
m_precision.F90
m_constants.F90
m_vib_finite_diff.F90
m_fourier_grid.F90
m_splines.F90
m_algebra.F90
m_color.F90
m_define_macro.F90
m_die.F90
m_sort.F90
m_upper.F90
m_warn.F90
m_sckh_params_t.F90
m_log.F90
m_log_type.F90
m_io.F90
m_input.F90
m_rkf45_matrix.F90
m_func.F90
m_KH.F90
m_PES_io.F90
m_KH_functions.F90
m_KH_utils.F90
m_SCKH_PES.F90
m_SCKH_utils.F90
hist_class.F90
m_FFT.F90
m_SCKH.F90
)

SET(ALL_MODULES_F 
sffteu.f
)
#bubblesort.f    crossprod.f     factorials.f    rnd.f           sffteu.f        splib.f         splines.f

#ADD_LIBRARY(fortran-utilities STATIC ${ALL_MODULES_F90} ${ALL_MODULES_F})

#add_subdirectory(morsefitting morsefitting)

SET(EXECUTABLES
vib_finite_diff sckh_main
)

FOREACH(EXEC ${EXECUTABLES})
  ADD_EXECUTABLE(${EXEC} ${EXEC}.F90 ${ALL_MODULES_F90} ${ALL_MODULES_F})
  TARGET_LINK_LIBRARIES(${EXEC} ${FLAGS} ${LAPACK_LIB}
    ${LINK_ADDONS})
ENDFOREACH(EXEC)

