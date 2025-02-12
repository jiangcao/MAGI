# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2019.0.117/linux/mkl
FC90 = /usr/sepp/bin/ifort-2020-af
MPIFC90 = /usr/pack/mpich-3.2.1-af/linux-x64/bin/mpif90
#F90_FLAGS =  -r8 -check bounds -traceback -fpp 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#F90_FLAGS =  -r8 -O2 -fpp -mkl -traceback  -qopenmp # -check bounds
F90_FLAGS = -O2 -march=native -mtune=native -mavx -ffree-line-length-none -fopenmp -fbacktrace -fbounds-check  -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -ffast-math #-Wall -Wextra 
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl 


# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: bse.x  gw.x  post.x 

bse.x : bse_abs.f90 parameters_mod.o sort.o xsf_mod.o utilities.o  linalg.o wannierHam3d.o setup_mod.o  bse_mod.o

	$(MPIFC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)

gw.x : gw.f90 parameters_mod.o sort.o xsf_mod.o utilities.o mkl_dfti.o linalg.o fft_mod.o wannierHam3d.o green_bse.o green.o green_rgf.o setup_mod.o

	$(MPIFC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)

post.x : post.f90 parameters_mod.o setup_mod.o utilities.o post_mod.o  linalg.o

	$(MPIFC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod $(SRCDIR)/*.mod

# implicit rules

%.o : %.f90
	$(MPIFC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS)  -I$(MODDIR) -J$(MODDIR)
