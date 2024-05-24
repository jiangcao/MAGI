F90FLAGS="-ffree-line-length-none -fbounds-check -fbacktrace -fopenmp"
rm src/*.mod
rm src/*.o
rm *.mod *.o
gfortran -c mkl_dfti.f90
f2py -c --fcompiler=gnu95 -m util src/util.f95 --f90flags="${F90FLAGS}" 
f2py -c --fcompiler=gnu95 -m wannier src/wannier.f95 --f90flags="${F90FLAGS}"
f2py -c --fcompiler=gnu95 -m negf src/negf.f95  skip: trimul_c mul_c : --f90flags="${F90FLAGS}"
