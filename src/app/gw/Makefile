# Fortran 90 compiler
# MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2019.0.117/linux/mkl/
# FC90 = /usr/sepp/bin/ifort-2020-af
# FC90 = /home/jiacao/openmpi/4.1.1-ifort/bin/mpif90
#FC90 = /usr/zupo/local/linux-local/mpich-4.0.2/gcc11/bin/mpifort

MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2020.0.166/linux/mkl/
FC90 = /usr/pack/mpich-3.2.1-af/linux-x64/bin/mpif90
F90_FLAGS = -Wall -Wextra -O2 -march=native -ffast-math -ffree-line-length-none -fopenmp -fbacktrace  -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function
#LIBS = -L ${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl 

MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#F90_FLAGS =  -I"${MKLROOT}/include" -r8 -O2 -fpp -mkl -traceback  -qopenmp -qopt-matmul # -check bounds  

#F90_FLAGS =  -g -O2 -ffree-line-length-none -fopenmp -fbacktrace

#LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl 

# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: gw_mpi.x 

gw_mpi.x : main3d.f90 mkl_dfti.o parameters.o fft_mod.o block_sparse.o wannierHam3d.o green.o gw_dense.o green_rgf.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)


liteQT_n.x : main.f90 mkl_dfti.o wannierHam.o green.o green_rgf.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)

bse.x : bse_abs.f90 wannierHam.o bse_mod.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)


.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod $(SRCDIR)/*.o $(SRCDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -I$(MODDIR) -J$(MODDIR) 

