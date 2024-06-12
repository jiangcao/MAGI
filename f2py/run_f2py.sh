
F90FLAGS="-g -march=native -O2 -ffree-line-length-none -fbounds-check -fbacktrace -ffast-math -fopenmp "
MKLROOT="/usr/pack/intel_compiler-2020-af/x64/mkl/"

rm src/*.mod
rm src/*.o
rm *.mod *.o
gfortran -c $MKLROOT/include/mkl_dfti.f90

python -m numpy.f2py  -c ../Modules/parameters.f95 --backend meson --fcompiler=gnu95 -m solver ../Modules/sinv.f95 --f90flags="${F90FLAGS}" --dep lapack  
python -m numpy.f2py  -c ../Modules/parameters.f95 --backend meson --fcompiler=gnu95 -m util src/util.f95 --f90flags="${F90FLAGS}" --dep lapack  
python -m numpy.f2py  -c ../Modules/parameters.f95  --backend meson --fcompiler=gnu95 -m wannier src/wannier.f95 --f90flags="${F90FLAGS}"  --dep lapack  
python -m numpy.f2py  -c mkl_dfti.f90 ../Modules/parameters.f95 ../Modules/fft.f95 ../Modules/matrix_c.f95  ../Modules/legendre.f95 \
 ../Modules/output.f95 ../Modules/observ.f95 ../Modules/open_boundary.f95  ../Modules/polarization.f95 \
 ../Modules/gw_dense.f95 ../Modules/bse_dense.f95  ../Modules/bse_sparse.f95 --fcompiler=gnu95 -m negf  skip: trimul_c mul_c : --f90flags="${F90FLAGS}"  --dep lapack 


