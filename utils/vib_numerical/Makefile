#F90 =/opt/intel/fc/9.1.036/bin/ifort # ny kompilator på qcllx07
#F90 =/opt/intel/fc/9.0/bin/ifort # kompilator på qcllx18
#F90FLAGS = -O3 -module /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/
#LFLAGS = -static -static-libcxa #qcllx07
#LFLAGS = -static-libcxa #qcllx18
#F90 = gfortran

#LOBJS8= -L/opt/intel/fc/9.1.036/lib -L/opt/intel/fc/9.1.036/lib -L/usr/libexec/CERNLIB/2004/lib/ -L/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/ -llapack3 -lblas 
#LOBJS8= /usr/libexec/CERNLIB/2004/lib/liblapack3.a  /usr/libexec/CERNLIB/2004/lib/libblas.a
#LOPT= -L/usr/lib/  -L/usr/libexec/CERNLIB/2004/lib  -lifport -lifcoremt -limf -lm \
#-lcxa -lpthread -lirc -lunwind -lc -lirc_s -llapack3 -lblas -lg2c 

#LOPT= -L/usr/lib/   -llapack -lblas  -lg2c 

#LOPT= -L/usr/lib/  -L/qcl/prog/ATLAS/lib/Linux_HAMMER32SSE2_2/ -llapack3 -lblas  -lg2c 

LOPT = -O2 -g -Wall -fbounds-check -ffree-line-length-none -framework vecLib

# order important in the following sequence: -llapack3 -lblas -lg2c 

#-L/opt/intel/fc/9.1.036/lib -L/usr/libexec/CERNLIB/2004/lib \
#-lifport -lifcoremt -limf -lm -lcxa -lpthread -lirc -lunwind -lc -lirc_s \
#-llapack3 -lblas -lg2c

LIBDIR = /Users/mathias/work/fortran-utilities_build/CMakeFiles/fortran-utilities.dir
INCDIR =  /Users/mathias/work/fortran-utilities_build
#LIBS = fortran-utilities

F90 = gfortran $(LOPT)  -I$(INCDIR) 

#-L$(LIBDIR) -l$(LIBS)

#LOBJS= /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/splines.o \
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/spline_easy.o \
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/linspace.o
#LOBJS2= /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/linspace.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/poly_n_fit.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/fcn3params.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/morsefitting/lmder.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/morsefitting/fitting_routines.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/splines.o \
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/spline_easy.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/factorials.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/splib.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/eigenfunctions.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/mat_elem_xn.o
#LOBJS3= /afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/mat_elem_xn.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/delta.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/svdcmp.o
#LOBJS4= /afs/physto.se/home/m/mathiasl/fortranfiler/LIB/crossprod.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/H2O.pjt2.o
#LOBJS5= /afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/mat_elem_xn.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/delta.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/svdcmp.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/splib.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/LIB/factorials.o\
#	/afs/physto.se/home/m/mathiasl/fortranfiler/vib_numerical/eigenfunctions.o


vib_finite_diff: vib_finite_diff.o
	$(F90) $(LFLAGS)  -o vib_finite_diff vib_finite_diff.o\
	 $(LIBDIR)/splines.o $(LIBDIR)/spline_easy.o $(LIBDIR)/linspace.o


#	 -L$(LIBDIR) -lfortran-utilities
vib_harm: vib_harm.o
	$(F90) $(LFLAGS)  -o vib_harm vib_harm.o\
	 $(LIBDIR)/linspace.o\
	 $(LIBDIR)/poly_n_fit.o\
	 $(LIBDIR)/fcn3params.o\
	 $(LIBDIR)/morsefitting/lmder.o\
	 $(LIBDIR)/morsefitting/fitting_routines.o\
	 $(LIBDIR)/splines.o \
	 $(LIBDIR)/spline_easy.o\
	 $(LIBDIR)/factorials.o\
	 $(LIBDIR)/splib.o\
	 $(LOBJS2) $(LOPT)
vib_normal3: vib_normal3.o
	$(F90) $(LFLAGS)  -o vib_normal3.x vib_normal3.o\
	 $(LOBJS3) $(LOPT)
nmode: nmode.o
	$(F90) $(LFLAGS)  -o nmode nmode.o\
	 $(LOBJS4) $(LOPT)
vib_normal3_int: vib_normal3_int.o fit_poly_surface.o
	$(F90) $(LFLAGS)  -o vib_normal3_int.x vib_normal3_int.o fit_poly_surface.o\
	 $(LOBJS5) $(LOPT)
convert_normal_to_internal: convert_normal_to_internal.o 
	$(F90) $(LFLAGS)  -o convert_normal_to_internal.x convert_normal_to_internal.o\
	 $(LOBJS4) $(LOPT)
vib_harm2: vib_harm2.o
	$(F90) $(LFLAGS)  -o vib_harm2 vib_harm2.o\
	 $(LOBJS2) $(LOPT)
vib_finite_diff_intensities: vib_finite_diff_intensities.o
	$(F90) $(LFLAGS) -o vib_finite_diff_intensities vib_finite_diff_intensities.o\
	 $(LOBJS) $(LOPT)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $*.f90



#clean: 
#	rm -f *.o watersort.out replaceSPC.out calculate_e-field.out \
#	 makeperiodic.out cutoutnearest.out makeinput.out frequencies.out\
#	mdtoxyz.out

#depend:
#	makedepend -- $(FLAGS) -- *.cpp
