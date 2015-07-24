cmake_minimum_required(VERSION 2.4)

PROJECT(fortran-utilities Fortran)

IF(CMAKE_Fortran_COMPILER MATCHES ifort)

  MESSAGE(STATUS  "ifort compiler found")

  SET(FLAGS "-O0 -warn all -g -CB")
  SET(INCDIR "-I/usr/include")
  #SET(OMP_FLAGS "")
  #SET(FFTW_LIB "-L/usr/lib -lfftw3 -lm")
  #SET(MKLPATH "/opt/intelcompiler/mkl/lib/intel64")
  SET(LAPACK_LIB "-Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_sequential.a ${MKLPATH}/libmkl_core.a -Wl,--end-group")

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
)

#SET(ALL_MODULES_F 
#bubblesort.f    crossprod.f     factorials.f    rnd.f           sffteu.f        splib.f         splines.f
#)


#ADD_LIBRARY(fortran-utilities STATIC ${ALL_MODULES_F90} ${ALL_MODULES_F})

#add_subdirectory(morsefitting morsefitting)

SET(EXECUTABLES
vib_finite_diff
)

FOREACH(EXEC ${EXECUTABLES})
  ADD_EXECUTABLE(${EXEC} ${EXEC}.F90 ${ALL_MODULES_F90})
  TARGET_LINK_LIBRARIES(${EXEC} ${FLAGS} ${LAPACK_LIBRARIES}
    ${LINK_ADDONS})
ENDFOREACH(EXEC)

