F90FLAGS=-march=native -ftree-parallelize-loops=28 -ffree-line-length-none -fbounds-check -fbacktrace -ffast-math -fopenmp -fexternal-blas
NVF90FLAGS=-g -traceback -fast -O3 -stdpar=gpu -mp -cudalib=cublas -lblas
MKLROOT="/usr/pack/intel_compiler-2020-af/x64/mkl/"
FC=gfortran
NVFC=nvfortran
F2PY=f2py
GCC=gcc
F2PYFLAGS=-c --fcompiler=gnu95 -lgomp -llapack -lblas --f90flags="${F90FLAGS}"
F2PYNVF90FLAGS=-c --fcompiler=nv -lgomp -llapack -lblas --f90flags="${NVF90FLAGS}"
OUTDIR=bin


nvhome=/usr/local/nvida_hpc_sdk
target=Linux_x86_64
version=23.9
nvcudadir="${nvhome}/${target}/${version}/cuda"
nvcompdir="${nvhome}/${target}/${version}/compilers"
nvmathdir="${nvhome}/${target}/${version}/math_libs"
nvcommdir="${nvhome}/${target}/${version}/comm_libs"

.PHONY: all clean

all: directories util wannier negf clean_compile

directories:
	mkdir -p ${OUTDIR}

util:
	F77=gfortran F90=gfortran CC=gcc ${F2PY} ${F2PYFLAGS} -m util src/magi/core/interface/util.f95
	mv util*.so ${OUTDIR}

wannier:
	F77=gfortran F90=gfortran CC=gcc ${F2PY} ${F2PYFLAGS} -m wannier src/magi/core/interface/wannier.f95
	mv wannier*.so ${OUTDIR}

negf: mkl_dfti.mod 
	F77=gfortran F90=gfortran CC=gcc ${F2PY} ${F2PYFLAGS} -m negf src/magi/core/interface/negf.f95 skip: trimul_c mul_c :
	mv negf*.so ${OUTDIR}

mkl_dfti.mod: 
	${FC} -c ${MKLROOT}/include/mkl_dfti.f90 

fortran.o:
	${GCC} -c src/magi/core/modules/gpu/fortran.c -I${nvcudadir}/include/ -DCUBLAS_GFORTRAN

clean_compile:
	rm -f *.mod *.o
	rm -f src/magi/core/interface/*.mod src/magi/core/interface/*.o

clean:
	rm -f ${OUTDIR}/*
