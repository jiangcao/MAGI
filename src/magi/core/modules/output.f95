!===============================================================================
! Copyright (C) 2023 Jiang Cao
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jiacao <jiacao@ethz.ch>
! Comment:
!  
! Maintenance:
!===============================================================================
module output 
    implicit none
    contains 

    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_per_kz(dataset,i,G,nen,en,nky,nkz,length,NB,Lx,coeff,at_ky,at_kz)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:) ! (m,m,e,k) kz is the fast-running index in k
        integer, intent(in)::i,nen,length,NB,nky,nkz
        real(8), intent(in)::Lx,en(nen),coeff(2)
        real(8), intent(in), optional::at_ky,at_kz
        integer:: ie,j,ib,ikz,iky,k,kb
        real(8):: kcenter(2),dky,dkz
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        kcenter=0.0d0
        if (present(at_ky)) then
            kcenter(1)=at_ky
        endif
        if (present(at_kz)) then
            kcenter(2)=at_kz
        endif
        dky=1.0d0/dble(nky)
        dkz=1.0d0/dble(nkz)
        open(unit=11,file=trim(dataset)//i_str//'_ky.dat',status='unknown')
        do ie = 1,nen
            ikz=max(min(nkz/2+ floor(kcenter(2)/dkz)+1 , nkz) , 1)            
            do iky=1,nky
                tr=0.0d0          
                do j = 1,length
                    do ib=1,nb
                        tr = tr+G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz+nkz*(iky-1))            
                    enddo
                end do
                write(11,'(4E20.6)') dble(iky)/dble(nky), en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            enddo
            write(11,*)    
        enddo
        close(11)

        open(unit=11,file=trim(dataset)//i_str//'_kz.dat',status='unknown')
        do ie = 1,nen
            iky=max(min(nky/2+ floor(kcenter(1)/dky)+1 , nkz) , 1)           
            do ikz=1,nkz
                tr=0.0d0          
                do j = 1,length
                    do ib=1,nb
                        tr = tr+G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz+nkz*(iky-1))            
                    enddo
                end do
                write(11,'(4E20.6)') dble(ikz)/dble(nkz), en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            enddo
            write(11,*)    
        enddo
        close(11)
    end subroutine write_spectrum_per_kz


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_summed_over_kz(dataset,i,G,nen,nsub,en,ensub,nkz,length,NB,Lx,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:,:)
        integer, intent(in)::i,nen,length,NB,nkz,nsub
        real(8), intent(in)::Lx,en(nen),coeff(2),ensub(nsub)
        integer:: ie,j,ib,ikz,isub
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do isub = 1,nsub
                do j = 1,length
                    tr=0.0d0          
                    do ib=1,nb
                    do ikz=1,nkz
                        tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,isub,ikz)            
                    enddo
                    enddo
                    tr=tr/dble(nkz)
                    write(11,'(4E20.6)') dble(j-1)*Lx, en(ie)+ensub(isub), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
                end do
                write(11,*)   
            enddo 
        enddo
        close(11)
    end subroutine write_spectrum_summed_over_kz


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_nosub(dataset,i,G,nen,en,length,NB,Lx,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:)
        integer, intent(in)::i,nen,length,NB
        real(8), intent(in)::Lx,en(nen),coeff(2)
        integer:: ie,j,ib
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do j = 1,length
                tr=0.0d0          
                do ib=1,nb                    
                    tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie)                                
                enddo                    
                write(11,'(4E20.6)') dble(j-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            end do
            write(11,*)   
        enddo
        close(11)
    end subroutine write_spectrum_nosub


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum(dataset,i,G,nen,nsub,en,ensub,length,NB,Lx,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:)
        integer, intent(in)::i,nen,length,NB,nsub
        real(8), intent(in)::Lx,en(nen),coeff(2),ensub(nsub)
        integer:: ie,j,ib,isub
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do isub = 1,nsub
                do j = 1,length
                    tr=0.0d0          
                    do ib=1,nb                    
                        tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,isub)                                
                    enddo                    
                    write(11,'(4E20.6)') dble(j-1)*Lx, en(ie)+ensub(isub), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
                end do
                write(11,*)   
            enddo 
        enddo
        close(11)
    end subroutine write_spectrum


    ! write current into file 
    subroutine write_current(dataset,i,cur,length,NB,NS,Lx)
        character(len=*), intent(in) :: dataset
        real(8), intent(in) :: cur(:,:)
        integer, intent(in)::i,length,NB,NS
        real(8), intent(in)::Lx
        integer:: j,ib,jb,ii
        real(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ii = 1,length-1
            tr=0.0d0          
            do ib=1,nb  
            do jb=1,nb       
                do j=ii,min(ii+NS-1,length-1)
                tr = tr+ cur((ii-1)*nb+ib,j*nb+jb)
                enddo
            enddo                        
            end do
            write(11,'(2E20.6)') dble(ii)*Lx, tr
        end do
    end subroutine write_current
    
    ! write current spectrum into file (pm3d map)
    subroutine write_current_spectrum(dataset,i,cur,nen,en,length,NB,Lx)
        character(len=*), intent(in) :: dataset
        real(8), intent(in) :: cur(:,:,:)
        integer, intent(in)::i,nen,length,NB
        real(8), intent(in)::Lx,en(nen)
        integer:: ie,j,ib,jb
        real(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do j = 1,length-1
                tr=0.0d0          
                do ib=1,nb  
                do jb=1,nb        
                    tr = tr+ cur((j-1)*nb+ib,j*nb+jb,ie)
                enddo                        
                end do
                write(11,'(3E20.6)') dble(j)*Lx, en(ie), tr
            end do
            write(11,*)    
        end do
        close(11)
    end subroutine write_current_spectrum

    ! write transmission spectrum into file
    subroutine write_transmission_spectrum(dataset,i,tr,nen,en)
        character(len=*), intent(in) :: dataset
        real(8), intent(in) :: tr(:)
        integer, intent(in)::i,nen
        real(8), intent(in)::en(nen)
        integer:: ie,j,ib
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen    
        write(11,'(2E20.6)') en(ie), dble(tr(ie))      
        end do
        close(11)
    end subroutine write_transmission_spectrum


    ! write trace of diagonal blocks
    subroutine write_trace(dataset,i,G,length,NB,Lx,coeff,E)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:)
        integer, intent(in)::i,length,NB
        real(8), intent(in)::Lx,coeff(2)
        real(8), intent(in),optional::E
        integer:: ie,j,ib
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown', position="append", action="write")
        do j = 1,length
            tr=0.0d0          
            do ib=1,nb
                tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib)            
            end do
            if (.not.(present(E))) then
            write(11,'(3E18.4)') (j-1)*Lx, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            else
            write(11,'(4E18.4)') (j-1)*Lx, E, dble(tr)*coeff(1), aimag(tr)*coeff(2)         
            endif
        end do
        write(11,*)
        close(11)
    end subroutine write_trace


    ! write a matrix for one energy index into a file
    subroutine write_matrix(dataset,i,G,en,length,NB,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:)
        integer, intent(in)::i,length,NB
        real(8), intent(in)::en,coeff(2)
        integer:: ie,j,ib,l
        complex(8)::tr
        logical :: lexist
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        inquire(file=trim(dataset)//i_str//'.dat', exist=lexist)
        if (lexist) then
            open(11, file=trim(dataset)//i_str//'.dat', status="old", position="append", action="write")
        else
            open(11, file=trim(dataset)//i_str//'.dat', status="new", action="write")
        end if
        fmt = '(E18.6,2I8,2E18.6)'
        do l=1,length*NB
            do j = 1,length*NB
                tr = G(l,j)            
                write(11,fmt) en,l,j, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            end do
        end do
        write(11,*)    
        close(11)
    end subroutine write_matrix

end module output


