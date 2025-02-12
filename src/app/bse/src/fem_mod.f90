!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 12/2023

module fem_mod

implicit none

private ! all entities are module-private by default
public type_mesh, fem_load_gmsh, fem_free_memory, fem_remove_duplicate_points

type type_mesh
   integer :: NP  ! number of nodes
   integer :: NE  ! number of elements
   integer :: NP_D  ! number of nodes of Dirichlet
   integer, allocatable :: TYP(:), TYP_D(:), TYP_S(:), REG(:)
   integer, allocatable :: NBN(:)     ! nbr of points in an element
   integer, allocatable :: NRG(:)     ! Region ID for an element
   integer, allocatable :: geo_ent(:) !
   integer, allocatable :: IND(:,:)   ! connection table
   real(8), allocatable :: XN(:)  ! coordinates of nodes
   real(8), allocatable :: YN(:)  ! coordinates of nodes
   real(8), allocatable :: ZN(:)  ! coordinates of nodes
end type type_mesh

contains

subroutine fem_free_memory(mesh)
type(type_mesh),intent(inout) :: mesh 
	if (allocated(mesh%typ)) deallocate(mesh%typ)
	if (allocated(mesh%typ_d)) deallocate(mesh%typ_d)
	if (allocated(mesh%typ_s)) deallocate(mesh%typ_s)
	if (allocated(mesh%reg)) deallocate(mesh%reg)
	if (allocated(mesh%nbn)) deallocate(mesh%nbn)
	if (allocated(mesh%nrg)) deallocate(mesh%nrg)
	if (allocated(mesh%geo_ent)) deallocate(mesh%geo_ent)
	if (allocated(mesh%ind)) deallocate(mesh%ind)
	if (allocated(mesh%xn)) deallocate(mesh%xn)
	if (allocated(mesh%yn)) deallocate(mesh%yn)
	if (allocated(mesh%zn)) deallocate(mesh%zn)
end subroutine fem_free_memory

subroutine fem_remove_duplicate_points(mesh)
	type(type_mesh),intent(inout) :: mesh 
	! ----
	logical,allocatable :: pt_toDelete(:),elemt_toDelete(:)
	integer :: NP_deleted,jj,NE_deleted,i,j,elemt,hash_val,ii
	integer,allocatable :: pt_shift(:),ind_new(:,:),nrg_new(:),TYP_new(:),NBN_new(:)
	integer, allocatable :: hash_elemt(:,:)
	real(8), allocatable :: XN_new(:),YN_new(:),ZN_new(:)
	!
	write(*,*) 'Double-point checking ...'
	allocate(pt_toDelete(mesh%NP))
	allocate(pt_shift(mesh%NP))
	pt_toDelete = .false.
	NP_deleted = 0
	allocate(hash_elemt(1:1000,1:(mesh%NP/100))) !! Hash-table for checking Double-points
	hash_elemt = 0
	do i = 1, mesh%NP
		hash_val = MOD(floor((abs(mesh%XN(i))+abs(mesh%YN(i))+abs(mesh%ZN(i)))*1e9),1000)+1    
		if (hash_elemt(hash_val,1) .gt. 0) then
			do elemt=1, hash_elemt(hash_val,1)            
				j = hash_elemt(hash_val,elemt+1)
				!! CHeck is two points are equal
				if (abs(mesh%XN(i)-mesh%XN(j)).lt.1.0D-10 .and. &
					abs(mesh%YN(i)-mesh%YN(j)).lt.1.0D-10 .and. &
					abs(mesh%ZN(i)-mesh%ZN(j)).lt.1.0D-10 )then
					write(*,*) '-- Find double-point:',j,i,' combining them ...'
					pt_toDelete(i) = .true.  !! prepare to delete point-i
					NP_deleted = NP_deleted + 1
					do ii  = 1, mesh%NE   !! substitute point-j by i in the connectivity table
						do jj = 1, mesh%NBN(ii)
							if (mesh%ind(ii,jj)==i) then
								mesh%ind(ii,jj) = j
							endif
						enddo
					enddo
				endif
			enddo
		endif
		if (.not. pt_toDelete(i)) then
			!! add point-i into the hash-table
			hash_elemt(hash_val,1) = hash_elemt(hash_val,1)+1
			if (hash_elemt(hash_val,1) .gt. (size(hash_elemt,2)-1)) then 
				write(*,*) "Failed, hash-table full"
				stop
			else
				hash_elemt(hash_val,hash_elemt(hash_val,1)+1)=i
			endif
		endif
	enddo
	if (NP_deleted > 0) then
		!! Be careful also to update the connectivity table
		write(*,*) "  Compressing point table ..."
		allocate(XN_new(mesh%NP-NP_deleted))
		allocate(YN_new(mesh%NP-NP_deleted))
		allocate(ZN_new(mesh%NP-NP_deleted))
		NP_deleted = 0
		do i = 1, mesh%NP
			if (pt_toDelete(i)) then
				NP_deleted = NP_deleted + 1
			else
				XN_new(i-NP_deleted) = mesh%XN(i)
				YN_new(i-NP_deleted) = mesh%YN(i)
				ZN_new(i-NP_deleted) = mesh%ZN(i)
			endif
			pt_shift(i) = NP_deleted
		enddo
		write(*,*) "Update connectivity table ... "
		!! update connectivity table
		do ii = 1, mesh%NE
			do jj = 1, mesh%NBN(ii)
				if (pt_toDelete(mesh%ind(ii,jj))) then
				write(*,*) "  Something wrong, point should have been deleted: ", mesh%ind(ii,jj),&
				" in element: ", ii
				stop
			endif
			mesh%ind(ii,jj) = mesh%ind(ii,jj) - pt_shift(mesh%ind(ii,jj))
			enddo
		enddo
		mesh%NP = mesh%NP-NP_deleted
		deallocate(mesh%XN,mesh%YN,mesh%ZN)
		call move_alloc(XN_new,mesh%XN)
		call move_alloc(YN_new,mesh%YN)
		call move_alloc(ZN_new,mesh%ZN)	
	endif
deallocate(pt_toDelete,pt_shift,hash_elemt)

end subroutine fem_remove_duplicate_points



subroutine fem_load_gmsh(fname, mesh)
	character(len=*),intent(in) :: fname 
	type(type_mesh),intent(inout) :: mesh 
	! ---- 
	integer :: fu,rc
	character(len=500) :: line, comment, gmesh_version,subline
	integer :: i,j,k,nb_tag,ident
	!
	open (action='read', file=fname, iostat=rc, newunit=fu)
	if (rc /= 0) then 
		write(*,*) 'file open error, file=',fname,', error code=',rc
		call abort
	endif
	read(fu,*)comment         ! $MeshFormat
	read(fu,*)gmesh_version   ! 2.2 0 8
	read(fu,*)comment         ! $EndMeshFormat
	read(fu,*)comment         ! $Nodes
	read(fu,*)mesh%NP              ! Point number

	write(*,*)'number of nodes in the structure:',mesh%NP

	allocate(mesh%XN(1:mesh%NP))
	allocate(mesh%YN(1:mesh%NP))
	allocate(mesh%ZN(1:mesh%NP))

	do i=1,mesh%NP
	   read(fu,*)j,mesh%XN(i),mesh%YN(i),mesh%ZN(i)
	enddo

	mesh%XN=mesh%XN*1.0d-7   ! transforms nm -> cm 
	mesh%YN=mesh%YN*1.0d-7   ! transforms nm -> cm
	mesh%ZN=mesh%ZN*1.0d-7   ! transforms nm -> cm

	read(fu,*)comment   ! $EndNodes
	read(fu,*)comment   ! $Elements
	read(fu,*)mesh%NE
	write(*,*) 'number of elements=',mesh%NE
	allocate(mesh%TYP(1:mesh%NE))
	allocate(mesh%NBN(1:mesh%NE))
	allocate(mesh%NRG(1:mesh%NE))
	allocate(mesh%ind(1:mesh%NE,1:8))
	allocate(mesh%geo_ent(1:mesh%NE))

	mesh%TYP=0
	mesh%NBN=0
	mesh%NRG=0
	mesh%ind=0
	mesh%geo_ent=0

	do i=1,mesh%NE
	    read(fu,'(A)',iostat=rc) line
	    if (rc.ne.0) exit

	    k = INDEX(line,' ')
	    k = k + 1   

	    if (line(k:k).eq.'1') then   !!! Line
	        read(line,*) ident,mesh%TYP(i),nb_tag,mesh%NRG(i), mesh%geo_ent(i),mesh%ind(i,1:2)
	        mesh%NBN(i) = 2
	    else if(line(k:k).eq.'2') then   !!! Triangle
	        read(line,*) ident,mesh%TYP(i),nb_tag,mesh%NRG(i), mesh%geo_ent(i),mesh%ind(i,1:3)		
	        mesh%NBN(i) = 3
	    else if(line(k:k).eq.'4') then  !!! Tetraede
	        read(line,*) ident,mesh%TYP(i),nb_tag,mesh%NRG(i), mesh%geo_ent(i),mesh%ind(i,1:4)
	        mesh%NBN(i) = 4
	   else if(line(k:k).eq.'3') then   !!! Rectangle
	       read(line,*) ident,mesh%TYP(i),nb_tag,mesh%NRG(i), mesh%geo_ent(i),mesh%ind(i,1:4)		
	       mesh%NBN(i) = 4
	   else if(line(k:k).eq.'5') then   !!! Cube
	       read(line,*) ident,mesh%TYP(i),nb_tag,mesh%NRG(i), mesh%geo_ent(i),mesh%ind(i,1:8)		
	       mesh%NBN(i) = 8
	    else 
	       write(*,*) 'Read unsupported element type ! in line ',i
	       write(*,*) line
	       stop 
	    endif
	enddo

	close(fu)

end subroutine fem_load_gmsh




end module fem_mod