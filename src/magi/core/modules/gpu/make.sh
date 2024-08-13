# module load ~/modulefiles/nvhpc_23.9

gcc -c fortran.c -I/usr/local/nvida_hpc_sdk/Linux_x86_64/24.5/cuda/12.4/include/ -DCUBLAS_GFORTRAN
nvfortran -c gpu_polarization.f95 fortran.o -cudalib=cublas -lblas -mp -stdpar 
nvfortran -c gpu_selfenergy.f95 fortran.o -cudalib=cublas -lblas -mp -stdpar
nvfortran test_polarization.f90 fortran.o gpu_polarization.o -o gpu_polarizability.x -cudalib=cublas -lblas -mp -stdpar
nvfortran test_selfenergy.f90 fortran.o gpu_selfenergy.o -o gpu_sigma.x -cudalib=cublas -lblas -mp -stdpar