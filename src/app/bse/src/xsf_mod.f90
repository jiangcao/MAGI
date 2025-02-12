module xsf_mod

    use parameters_mod, only : dp

    implicit none

    private 

    public :: xsf_write_header, xsf_write_3ddatablock

    contains


    
    SUBROUTINE xsf_write_header (primvec,convvec,coor,iat, nat, ounit)
    !     writes the header for XSF structure file    
    integer :: i,j
    integer,intent(in)::iat, ounit
    character(len=2),intent(in):: nat(iat)
    real(dp),intent(in):: primvec(3,3)        ! primitive lattice vectors
    real(dp),intent(in):: convvec(3,3)        ! conventional lattice vectors
    real(dp),intent(in):: coor(3,iat)     ! atomic coordinates
    !    
    write(ounit,*) 'CRYSTAL'
    write(ounit,*) 'PRIMVEC'
    write(ounit,1000) primvec
    write(ounit,*) 'CONVVEC'
    write(ounit,1000) convvec
    write(ounit,*) 'PRIMCOORD'
    write(ounit,*) iat, 1
    do i=1,iat
        write(ounit,1001) nat(i), (coor(j,i), j=1,3)
    enddo
    !
    1000 format(2(3(F15.9,2X),/),3(F15.9,2X))
    1001 format(a2,3x,3(F15.9,2X))
    end subroutine xsf_write_header
    


    subroutine xsf_write_3ddatablock(value, nx, ny, nz,origin,primvec, ounit)
        real(dp), intent(in) :: value(:,:,:)
        integer, intent(in)  :: nx,ny,nz , ounit
        real(dp),intent(in):: origin(3)        ! origin
        real(dp),intent(in):: primvec(3,3)       ! primitive lattice vectors
        integer::ix,iy,iz
        write(ounit,*) 'BEGIN_BLOCK_DATAGRID_3D'
        write(ounit,*) '3D_field'
        write(ounit,*) 'BEGIN_DATAGRID_3D'
        write(ounit,1000) nx,ny,nz
        write(ounit,1001) origin
        write(ounit,1001) primvec
        write(ounit,'(6E15.5)') (((value(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
        write(ounit,*) 'END_DATAGRID_3D'
        write(ounit,*) 'END_BLOCK_DATAGRID_3D'
        1000 format(3(I15,2X))
        1001 format(3(F15.9,2X))
    end subroutine xsf_write_3ddatablock

end module xsf_mod