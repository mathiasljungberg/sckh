#F90 =/opt/intel/fc/9.1.036/bin/ifort # ny kompilator på qcllx07
#F90 =/opt/intel/fc/9.0/bin/ifort # kompilator på qcllx18
#F90 =/opt/intel/fc/10.1.015/bin/ifort # kompilator p qcllx21
#F90 =gfortran
#F90= /afs/physto.se/qcl/prog/g95/g95-install/bin/i686-suse-linux-gnu-g95
#F90=pgf90 #qcllx22
#F90FLAGS = -O3 -CB -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/   #intel compiler
#F90FLAGS = -O3 -fno-underscoring -fbounds-check -J/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #gfortran compiler
#F90FLAGS = -O3 -fbounds-check -ffree-line-length-huge -I/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #g95
#F90FLAGS = -O3 -C -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #pgf90
#F90FLAGS = -O3 -C -module .    #pgf90
F90=

#LFLAGS = -static -static-libcxa #qcllx07
#LFLAGS = -static-libcxa #qcllx18

#LOPT= -L/usr/lib/  -L/usr/libexec/CERNLIB/2004/lib  -llapack3 -lblas  -lg2c # intel compiler

# order important in the following sequence: -llapack3 -lblas -lg2c 

LOBJS= splines.o\
spline_m.o\
KH_functions.o\

all: parameters spline_m KH_functions sinc_DVR

parameters: parameters.f90
	$(F90) $(F90FLAGS) -c parameters.f90

spline_m: spline_m.f90 parameters.f90
	$(F90) $(F90FLAGS) -c spline_m.f90

KH_functions: KH_functions.f90
	$(F90) $(F90FLAGS) -c KH_functions.f90

sinc_DVR: sinc_DVR.o $(LOBJS) 
	$(F90) $(LFLAGS) -o sinc_DVR sinc_DVR.o\
	$(LOBJS) $(LOPT)

%.o: %.f
	$(F90) $(F90FLAGS) -c $*.f

%.o: %.f90
	$(F90) $(F90FLAGS) -c $*.f90

clean: 
	rm -f *.o SCKH

#depend:
#	makedepend -- $(FLAGS) -- *.cpp
