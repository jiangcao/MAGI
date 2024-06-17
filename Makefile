F90FLAGS=-march=native -ffree-line-length-none -fbounds-check -fbacktrace -ffast-math -fopenmp -fexternal-blas
MKLROOT="/usr/pack/intel_compiler-2020-af/x64/mkl/"
FC=gfortran
F2PY=f2py
F2PYFLAGS=-c --fcompiler=gnu95 -lgomp -llapack -lblas --f90flags="${F90FLAGS}"
OUTDIR=bin

.PHONY: all clean

all: directories util wannier negf clean_compile

directories:
	mkdir -p ${OUTDIR}

util:

	${F2PY} ${F2PYFLAGS} -m util src/magi/core/interface/util.f95
	mv util*.so ${OUTDIR}

wannier:
	${F2PY} ${F2PYFLAGS} -m wannier src/magi/core/interface/wannier.f95
	mv wannier*.so ${OUTDIR}

negf: mkl_dfti.mod
	${F2PY} ${F2PYFLAGS} -m negf src/magi/core/interface/negf.f95 skip: trimul_c mul_c :
	mv negf*.so ${OUTDIR}

mkl_dfti.mod: 
	${FC} -c ${MKLROOT}/include/mkl_dfti.f90 

clean_compile:
	rm -f *.mod
	rm -f src/magi/core/interface/*.mod

clean:
	rm -f ${OUTDIR}/*