#F90 =/opt/intel/fc/9.1.036/bin/ifort # ny kompilator p� qcllx07
#F90 =/opt/intel/fc/9.0/bin/ifort # kompilator p� qcllx18
#F90 =/opt/intel/fc/10.1.015/bin/ifort # kompilator p qcllx21
#F90 =gfortran
#F90= /afs/physto.se/qcl/prog/g95/g95-install/bin/i686-suse-linux-gnu-g95
#F90=pgf90 #qcllx22
F90=gfortran
#F90FLAGS = -O3 -CB -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/   #intel compiler
#F90FLAGS = -O3 -fno-underscoring -fbounds-check -J/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #gfortran compiler
#F90FLAGS = -O3 -fbounds-check -ffree-line-length-huge -I/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #g95
#F90FLAGS = -O3 -C -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/    #pgf90
#F90FLAGS = -O3 -C -module .    #pgf90
F90FLAGS = -Wall -fbounds-check -ffree-line-length-none -ffixed-line-length-none # gfortran
#LFLAGS = -static -static-libcxa #qcllx07
#LFLAGS = -static-libcxa #qcllx18

# order important in the following sequence: -llapack3 -lblas -lg2c 
#LOPT= -L/usr/lib/  -L/usr/libexec/CERNLIB/2004/lib  -llapack3 -lblas  -lg2c # intel compiler

LOPT= -framework veclib #-L/Developer/SDKs/MacOSX10.5.sdk/usr/lib -llapack -lblas # gfortran/mac

LOBJS= splines.o \
spline_m.o \
KH_functions.o \

all: parameters spline_m KH_functions KH

parameters: parameters.f90
	$(F90) $(F90FLAGS) -c parameters.f90

spline_m: spline_m.f90 parameters.f90
	$(F90) $(F90FLAGS) -c spline_m.f90

KH_functions: KH_functions.f90
	$(F90) $(F90FLAGS) -c KH_functions.f90

KH: KH.o $(LOBJS) 
	$(F90) $(LFLAGS) -o KH KH.o \
	$(LOBJS) $(LOPT)

%.o: %.f
	$(F90) $(F90FLAGS) -c $*.f

%.o: %.f90
	$(F90) $(F90FLAGS) -c $*.f90

clean: 
	rm -f *.o SCKH

#depend:
#	makedepend -- $(FLAGS) -- *.cpp
