
program test_polarization
    use gpu_polarization
    implicit none
    complex(8), dimension(:,:), allocatable :: A, B, C, ref
    complex(8), dimension(:,:), allocatable :: G_lesser,G_greater,G_retarded,G_advanced
    real(8) :: tstart, tstop, elapsed_time
    real(8) :: gflops, sum, L2
    real :: randnum
    integer, allocatable, dimension(:) :: jl,ki
    integer(8) :: devPtrA, devPtrB, devPtrC, devPtrGL,devPtrGG,devPtrGR,devPtrGA
    integer :: n=1400
    integer :: m=140
    integer :: nen=2048
    integer :: nop=10
    ! integer :: size_of_real=8 !4->single precision; 8->double precision
    ! integer :: size_of_complex=16 !4->single precision; 8->double precision
    integer :: i,j,ie
  
    integer,dimension(8) :: values
    complex(8) :: a1,a2
    integer :: seed
    integer :: index
  
    real(8),parameter :: pi = 4.0*atan(1.0)
    real(8),parameter :: alpha = 0.99

    a1 = 1.0 - alpha
    a2 = alpha * 0.5

    call date_and_time(VALUES=values) !values(8) = milisecs of the second
    seed = values(8) !using value in milisecs as seeder
    print *, seed
    call srand(seed) !not a std implementation, but i like it better.
    
    allocate(G_lesser(nen,m*m))
    allocate(G_greater(nen,m*m))
    allocate(G_retarded(nen,m*m))
    allocate(G_advanced(nen,m*m))

    allocate(jl(n))
    allocate(ki(n))
    allocate(C(n,n))
    allocate(ref(n,n))
    
    !solution matrix
    C = dcmplx(0.0,0.0)
    !reference matrix
    ref = dcmplx(0.0,0.0)
    !initiating matrix value
    
    do j=1,m*m
        do i=1,nen
            call random_number(randnum)
            G_greater(i,j) = dcmplx(randnum,0.0)
            call random_number(randnum)
            G_lesser(i,j) = dcmplx(randnum,0.0)
            call random_number(randnum)
            G_advanced(i,j) = dcmplx(randnum,0.0)
            call random_number(randnum)
            g_retarded(i,j) = dcmplx(randnum,0.0)
        enddo
    enddo

    do j = 1, n
        jl(j) = j
        ki(j) = j 
    enddo 

    call cpu_time(tstart)
    call gpu_polarization(a1=a1,a2=a2,nop=nop,nen=nen,nm=m,num_jl=n,jl=jl,num_ki=n,ki=ki,copy_to_gpu=.true.,&
        G_lesser=G_lesser,G_greater=G_greater,G_retarded=G_retarded,G_advanced=G_advanced,&
        devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,partial_P=C)
    call cpu_time(tstop)

    elapsed_time = tstop - tstart !in seconds
    write(*,20) 'Elapsed time : ',elapsed_time, 'secs'

    ! print *,C(1,1)
    
    20 format(A15,2X,1F0.8,2X,A4)  

    do j=1,200000
        ki = ki + 1
        if (ki(1) > m*m) then 
            ki = 1
        endif

        call cpu_time(tstart)

        call gpu_polarization(a1=a1,a2=a2,nop=nop,nen=nen,nm=m,num_jl=n,jl=jl,num_ki=n,ki=ki,copy_to_gpu=.false.,&
            G_lesser=G_lesser,G_greater=G_greater,G_retarded=G_retarded,G_advanced=G_advanced,&
            devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,partial_P=C)

        call cpu_time(tstop)

        elapsed_time = tstop - tstart !in seconds
        ! print *,C(1,1)
        write(*,20) 'Elapsed time : ',elapsed_time, 'secs'
    enddo
    
    print *,C(1,1)
    
    !Free GPU memory
    call cublas_free(devPtrGR)
    call cublas_free(devPtrGL)
    call cublas_free(devPtrGG)

!    !reference from CPU
!    call cpu_time(tstart)
!    do i=1,n 
!        do j=1,n 
!            ref(i,j) = a2* dot_product( G_greater(1+nop:nen, jl(i)) , G_lesser(1:nen-nop, ki(j)) ) &
!                     - a2* dot_product( G_lesser(1+nop:nen, jl(i)) , G_greater(1:nen-nop, ki(j)) ) &
!                     + a1* dot_product( G_lesser(1+nop:nen, jl(i)) , G_advanced(1:nen-nop, ki(j)) ) &
!                     + a1* dot_product( G_retarded(1+nop:nen, jl(i)) , G_lesser(1:nen-nop, ki(j)) ) 
!
!        enddo 
!    enddo 
!    call cpu_time(tstop)
!    
!    elapsed_time = tstop - tstart !in seconds
!    write(*,20) 'CPU Elapsed t : ',elapsed_time, 'secs'
!    
!    print *,C(1,1)
!    print *,ref(1,1)
!    index = n
!    sum = 0.0
!    !show result
!    do j = 1, index
!       do i = 1, index
!          sum = sum + abs( ref(i,j)- C(i,j) )**2
!       end do
!    end do
!    L2 = sqrt( sum / ( dble(index)*dble(index) ) )
!    write(*,30) 'L2-residual :', L2
!  30 format(A15,2X,1E18.8)
!    
    
  end program test_polarization
  
