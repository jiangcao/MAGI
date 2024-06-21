
program test_selfenergy
    use gpu_selfenergy
    implicit none
    complex(8), dimension(:,:), allocatable :: A, B, C, ref
    complex(8), dimension(:,:), allocatable :: G_lesser,G_greater,G_retarded
    complex(8), dimension(:), allocatable :: W_lesser,W_greater,W_retarded
    complex(8), dimension(:,:), allocatable :: sig_lesser,sig_greater,sig_retarded
    real(8) :: tstart, tstop, elapsed_time
    real(8) :: gflops, sum, L2
    real :: randnum
    integer, allocatable, dimension(:) :: ij
    integer(8) :: devPtrA, devPtrB, devPtrC, devPtrGL,devPtrGG,devPtrGR
    integer :: n=1400
    integer :: m=140
    integer :: nen=2048
    integer :: nop=10
    ! integer :: size_of_real=8 !4->single precision; 8->double precision
    ! integer :: size_of_complex=16 !4->single precision; 8->double precision
    integer :: i,j,ie,ie1,ie2
  
    integer,dimension(8) :: values
    complex(8) :: a1,a2
    integer :: seed
    integer :: index
  
    real(8),parameter :: pi = 4.0*atan(1.0)

    call date_and_time(VALUES=values) !values(8) = milisecs of the second
    seed = values(8) !using value in milisecs as seeder
    print *, seed
    call srand(seed) !not a std implementation, but i like it better.
    
    allocate(G_lesser(nen,m*m))
    allocate(G_greater(nen,m*m))
    allocate(G_retarded(nen,m*m))

    allocate(sig_lesser(nen,m*m))
    allocate(sig_greater(nen,m*m))
    allocate(sig_retarded(nen,m*m))

    allocate(w_lesser(m*m))
    allocate(w_greater(m*m))
    allocate(w_retarded(m*m))

    allocate(ij(m*m))
    ! allocate(C(n,n))
    allocate(ref(nen,m*m))
    
    !solution matrix
    ! C = dcmplx(0.0,0.0)
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
            g_retarded(i,j) = dcmplx(randnum,0.0)
        enddo
        call random_number(randnum)
        w_greater(j) = dcmplx(randnum,0.0)
        call random_number(randnum)
        w_retarded(j) = dcmplx(randnum,0.0)
        call random_number(randnum)
        w_lesser(j) = dcmplx(randnum,0.0)
    enddo

    do j = 1, m*m
        ij(j) = j
    enddo 
    sig_lesser = dcmplx(0.0,0.0)
    sig_greater= dcmplx(0.0,0.0)
    sig_retarded = dcmplx(0.0,0.0)
    call cpu_time(tstart)
    call gpu_selfenergy_GW(nen=nen,nop=nop,nm=m,num_ij=m*m,ij=ij,copy_to_gpu=.true.,&
                        G_lesser=g_lesser,G_greater=g_greater,G_retarded=g_retarded,&
                        W_lesser=W_lesser,W_greater=W_greater,W_retarded=W_retarded,&
                        devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,&
                        sig_lesser=sig_lesser,sig_greater=sig_greater,sig_retarded=sig_retarded)
    call cpu_time(tstop)

    elapsed_time = tstop - tstart !in seconds
    write(*,20) 'Elapsed time : ',elapsed_time, 'secs'
    print *,sig_retarded(nen/2,10)

    ! print *,C(1,1)
    
    20 format(A15,2X,1F0.8,2X,A4)  

    do j=1,10000
    sig_lesser = dcmplx(0.0,0.0)
    sig_greater= dcmplx(0.0,0.0)
    sig_retarded = dcmplx(0.0,0.0)        
        call cpu_time(tstart)

        call gpu_selfenergy_GW(nen=nen,nop=nop,nm=m,num_ij=m*m,ij=ij,copy_to_gpu=.false.,&
                            G_lesser=g_lesser,G_greater=g_greater,G_retarded=g_retarded,&
                            W_lesser=W_lesser,W_greater=W_greater,W_retarded=W_retarded,&
                            devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,&
                            sig_lesser=sig_lesser,sig_greater=sig_greater,sig_retarded=sig_retarded)

        call cpu_time(tstop)

        elapsed_time = tstop - tstart !in seconds
        ! print *,C(1,1)
        write(*,20) 'Elapsed time : ',elapsed_time, 'secs'
        print *,sig_retarded(nen/2,10)
    enddo

    !Free GPU memory
    call cublas_free(devPtrGR)
    call cublas_free(devPtrGL)
    call cublas_free(devPtrGG)

    sig_lesser = dcmplx(0.0,0.0)
    sig_greater= dcmplx(0.0,0.0)
    ref = dcmplx(0.0,0.0)
    !reference from CPU
    call cpu_time(tstart)
    !$omp parallel default(shared) private(i,ij,ie1,ie2) 
    !$omp do
    do i=1,m*m                                
        ie1 = max(nop,0) + 1
        ie2 = min(nen+nop,nen)
        Sig_lesser(ie1:ie2,ij(i))=Sig_lesser(ie1:ie2,ij(i)) + G_lesser((ie1-nop):(ie2-nop),ij(i)) * W_lesser(ij(i))                                
        Sig_greater(ie1:ie2,ij(i))=Sig_greater(ie1:ie2,ij(i)) + G_greater((ie1-nop):(ie2-nop),ij(i)) * W_greater(ij(i))   
        ref(ie1:ie2,ij(i))=ref(ie1:ie2,ij(i)) + &
                                G_lesser((ie1-nop):(ie2-nop),ij(i)) * W_retarded(ij(i)) + &                                      
                                G_retarded((ie1-nop):(ie2-nop),ij(i)) * W_lesser(ij(i)) + &
                                G_retarded((ie1-nop):(ie2-nop),ij(i)) * W_retarded(ij(i))                                                  
        !
        ie1 = max(-nop,0) + 1
        ie2 = min(nen-nop,nen)
        Sig_lesser(ie1:ie2,ij(i))=Sig_lesser(ie1:ie2,ij(i)) + G_lesser((ie1+nop):(ie2+nop),ij(i)) * W_greater(ij(i))   
        Sig_greater(ie1:ie2,ij(i))=Sig_greater(ie1:ie2,ij(i)) + G_greater((ie1+nop):(ie2+nop),ij(i)) * W_lesser(ij(i))   
        ref(ie1:ie2,ij(i))=ref(ie1:ie2,ij(i)) - &
                                G_lesser((ie1+nop):(ie2+nop),ij(i)) * conjg(W_retarded(ij(i))) - &                                      
                                G_retarded((ie1+nop):(ie2+nop),ij(i)) * conjg(W_greater(ij(i))) - &
                                G_retarded((ie1+nop):(ie2+nop),ij(i)) * conjg(W_retarded(ij(i)))     
    
    enddo
    !$omp end do
    !$omp end parallel 
    call cpu_time(tstop)
    
    elapsed_time = tstop - tstart !in seconds
    write(*,20) 'CPU Elapsed t : ',elapsed_time, 'secs'
    
    print *,sig_retarded(1,10)
    print *,ref(1,10)
    index = n
    sum = 0.0
    !show result
    do j = 1, m*m
       do i = 1, nen
          sum = sum + abs( ref(i,j)- sig_retarded(i,j) )**2
       end do
    end do
    L2 = sqrt( sum / ( dble(m*m)*dble(nen) ) )
    write(*,30) 'L2-residual :', L2
  30 format(A15,2X,1E18.8)
    
    
  end program test_selfenergy
  
