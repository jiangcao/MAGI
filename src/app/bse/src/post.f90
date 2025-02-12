PROGRAM post
    use parameters_mod,only:dp
    use setup_mod
    use post_mod, only : W_to_real_space
    use wannierHam3d, only : nb,w90_load_from_file,w90_free_memory

    implicit none
    complex(kind=dp), allocatable :: wrspace(:,:,:,:)

    lread_input_bse = .true.
    lread_input_post = .true.
    lread_input_ham = .true.

    call begin_setup()
    
    call W_to_real_space(trim(dataset), number, nky, nkz, nen, nb, wrspace)

    deallocate(wrspace)
    call w90_free_memory()
    call MPI_Finalize( ierr )
end program post
