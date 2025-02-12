PROGRAM Poisson
USE fem_mod, only: fem_load_gmsh,type_mesh,fem_free_memory,fem_remove_duplicate_points

implicit none

type(type_mesh)::mesh 
! MPI variables
integer ( kind = 4 ) ierr

include "mpif.h"

call MPI_Init(ierr)

    
call fem_load_gmsh('struct.msh', mesh)

call fem_remove_duplicate_points(mesh)

call fem_free_memory(mesh)

call MPI_Finalize( ierr )

END PROGRAM Poisson


