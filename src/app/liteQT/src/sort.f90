!===============================================================================
! Copyright (C) 2023 Jiang Cao
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: Jiang Cao <jiacao@ethz.ch>
! Comment:
!  
! Maintenance:
!===============================================================================
module sort
    use parameters_mod,only:dp
    implicit none 

    private 

    public :: dmerge_sort

    contains

    ! INPUT
    !     V       : vector to sort                                      : real(double), dimension(:)
    !     Z       : vector of the elements position                     : integer(:)
    !     p, r    : staring and ending index where perform the sorting  : integer, optional
    ! OUTPUT
    !     Z       : Z(p:q) contains the indices of the ordered elements : integer(:)
    ! INTERNAL
    !     pp, rr  : internal copy of p and r                            : integer
    !     q       : intemediate index                                   : integer
    pure recursive subroutine dmerge_sort(V,Z,p,r)
    implicit none
    integer, intent(in), optional :: p, r
    real(dp), dimension(:), intent(in) :: V
    integer, dimension(:), intent(inout) :: Z
    integer :: pp, rr, q
    if(present(p)) then
        pp=p
    else
        pp=1               ! if p is not present we are at the first iteration and p is set to 1
    end if
    if(present(r)) then
        rr=r
    else
        rr=size(V)         ! if r si not present we are at the first iteration and r is set to the length of V
    end if
    if(.not.present(p).and..not.present(r)) then
        do q=1,size(V)
            Z(q)=q          ! if p and r are not present we are at the first iteration and Z=1,2,3,... length of V
        end do
    end if
    if(pp<rr) then        ! if p<r split the calculation
        q=(pp+rr)/2
        call dmerge_sort(V=V,Z=Z,P=pp,R=q)        ! sorting from index p to index q
        call dmerge_sort(V=V,Z=Z,P=q+1,R=rr)      ! sorting from index q+1 to index r
        call dmerge(V=V,Z=Z,P=pp,Q=q,R=rr)        ! merging the two sortes parts
    end if
    end subroutine dmerge_sort


    pure subroutine dmerge(V,Z,p,q,r)
    implicit none
    integer, intent(in) :: p, q ,r
    real(dp), dimension(:), intent(in) :: V
    integer, dimension(:), intent(inout) :: Z
    integer, dimension(:), allocatable :: zleft, zright
    integer :: n1, n2 
    integer :: i, j, k
    n1=q-p+1
    n2=r-q
    allocate(zleft(n1),zright(n2))
    zleft(1:n1)=Z(p:q)
    zright(1:n2)=Z(q+1:r)
    i=1
    j=1
    do k=p,r
        if(i<=n1.and.j<=n2) then
            if(V(zleft(i))<V(zright(j))) then
            Z(k)=zleft(i)
            i=i+1
            else
            Z(k)=zright(j)
            j=j+1
            end if
        elseif(i>n1.and.j<=n2) then
            Z(k)=zright(j)
            j=j+1
        elseif(j>n2.and.i<=n1) then
            Z(k)=zleft(i)
            i=i+1
        end if
    end do
    deallocate(zleft,zright)  
    end subroutine dmerge

end module sort 