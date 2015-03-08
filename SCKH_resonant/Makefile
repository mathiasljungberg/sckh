#F90 =/opt/intel/fc/9.1.036/bin/ifort # ny kompilator på qcllx07
#F90 =/opt/intel/fc/9.0/bin/ifort # kompilator på qcllx18
#F90 =/opt/intel/fc/10.1.015/bin/ifort # kompilator p qcllx21
#F90 =gfortran
#F90= /afs/physto.se/qcl/prog/g95/g95-install/bin/i686-suse-linux-gnu-g95
F90=gfortran
#F90=pgf90 #qcllx22
#F90FLAGS = -O3 -CB -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/   #intel compiler
#F90FLAGS = -O3 -fno-underscoring -fbounds-check -J/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #gfortran compiler
#F90FLAGS = -O3 -fbounds-check -ffree-line-length-huge -I/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #g95
#F90FLAGS = -O3 -C -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #pgf90
#F90FLAGS = -O3 -C -module .    #pgf90
#F90FLAGS = -Wall -fbounds-check -ffree-line-length-none -ffixed-line-length-none #-O3 -C -module .    #pgf90
#F77FLAGS = -Wall -fbounds-check #-ffree-line-length-none -ffixed-line-length-none #-O3 -C -module .    #pgf90
F90FLAGS= -I/Users/mathias/work/work_bcn/fortranfiler/LIB -J/Users/mathias/work/work_stockholm/fortranfiler/Wave_packet_dynamics/SCKH_resonant -fno-underscoring -fbounds-check -ffree-line-length-none -framework vecLib # gfortran/mac lion

F77FLAGS= -I/Users/mathias/work/work_bcn/fortranfiler/LIB -J/Users/mathias/work/work_stockholm/fortranfiler/Wave_packet_dynamics/SCKH_resonant -fno-underscoring -fbounds-check -ffree-line-length-none -framework vecLib # gfortran/mac lion


#LFLAGS = -static -static-libcxa #qcllx07
#LFLAGS = -static-libcxa #qcllx18

#LOPT= -L/usr/lib/  -L/usr/libexec/CERNLIB/2004/lib  -llapack3 -lblas  -lg2c # intel compiler
LOPT= -L/Developer/SDKs/MacOSX10.5.sdk/usr/lib -llapack -lblas 

# order important in the following sequence: -llapack3 -lblas -lg2c 


#LOBJS= linspace.o\
#splines.o\
#spline_easy.o\
#sffteu.o\
#SCKH_functions.o\
#KH_functions.o

LOBJS= splines.o\
sffteu.o\
spline_m.o\
FFT_m.o\
SCKH_functions.o\
KH_functions.o\
hist_class.o

all: parameters FFT_m spline_m KH_functions SCKH_functions hist_class SCKH_resonant 


FFT_m: FFT_m.f90
	$(F90) $(F90FLAGS) -c FFT_m.f90

hist_class: hist_class.f90
	$(F90) $(F90FLAGS) -c hist_class.f90

parameters: parameters.f90
	$(F90) $(F90FLAGS) -c parameters.f90

spline_m: spline_m.f90 parameters.f90
	$(F90) $(F90FLAGS) -c spline_m.f90

KH_functions: KH_functions.f90
	$(F90) $(F90FLAGS) -c KH_functions.f90

SCKH_functions: SCKH_functions.f90
	$(F90) $(F90FLAGS) -c SCKH_functions.f90

SCKH_resonant: $(LOBJS) SCKH_resonant.o 
	$(F90) $(F90FLAGS) -o SCKH_resonant SCKH_resonant.o \
	$(LOBJS) #$(LOPT)

%.o: %.f
	$(F90) $(F77FLAGS) -c $*.f

%.o: %.f90
	$(F90) $(F90FLAGS) -c $*.f90

clean: 
	rm -f *.o *.mod SCKH_resonant 

#depend:
#	makedepend -- $(FLAGS) -- *.cpp
