# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2020.0.166/linux/mkl/
FC90 = nvfortran 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
F90_FLAGS = -Wall -Wextra  -traceback  -cudalib=cublas -lblas -mp=gpu -gpu=nvlamath -cudalib=nvlamath  -llapack -O3 -g
LIBS = -stdpar -gpu=nvlamath -cudalib=nvlamath   -llapack 

# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: qt.x 

qt.x : main.f90  static.o  mkl_dfti.o linalg.o matrix_c.o graph_partition.o wannierHam3d.o deviceHam_mod.o  rgf_mod.o cuda_rgf_mod.o output_mod.o negf_mod.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -module $(MODDIR) 

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -module $(MODDIR) 
