program CSR
    use diffusion_equ

    implicit none

    ! Arrays for the matrix sparse row (CSR) format
    integer*4, allocatable :: ia(:), ja(:)
    real*8, allocatable :: a(:), f(:), c(:)

    ! Maximum size of matrix and the maximum number of non-zero entries
    integer, parameter :: maxn = 100000 
    integer, parameter :: maxnz = 1000000

    ! Work arrays to keep ILU factors and 8 vectors 
    integer, parameter :: MaxWi = maxnz + 2*maxn + 1
    integer, parameter :: MaxWr = maxnz + 8*maxn
    real*8 :: rW(MaxWr)
    integer :: iW(MaxWi)

    ! BiCGStab data
    external matvec, prevec0
    integer :: ITER, INFO, NUNIT
    real*8 :: RESID

    ! ILU0 data
    integer :: ierr, ipaLU, ipjLU, ipjU, ipiw

    ! External routines from BLAS library
    real*8 :: ddot
    external ddot, dcopy    

    ! Local variables
    integer :: imatvec(1), iprevec(1), ipBCG
    real*8  :: resinit

    ! Local 
    integer*4 :: N, i, j, nz
    real*8 :: c_error_system

    N = 64
    
    ! 5 - diagonal matrix
    allocate( ia((N-1)*(N-1) + 1) )
    allocate( ja( 5*(N-1)*(N-1) ) )
    allocate( a( 5*(N-1)*(N-1) ) )
    allocate( f((N-1)*(N-1)), c((N-1)*(N-1)) )

    !ia = 0
    !ja = 0
    !a = 0

    ! ======================================================================
    !  Generate system in CSR format
    ! ======================================================================
    call generate_matrix(N, ia, ja, a, f)	

    !print *, IA
    !print *, A
    !print *, JA

    ! Print matrix
    call draw_matrix((N-1)*(N-1), IA, JA, "matrix.ps")
    
    nz = ia((N-1)*(N-1) + 1)
    print *, nz

    ! ======================================================================
    !  Initialization of the preconditioner
    ! ======================================================================
    ipaLU = 1
    ipBCG = ipaLU + nz
    ipjU  = 1
    ipjLU = ipjU + (N-1)*(N-1) + 1
    ipiw  = ipjLU + nz ! work array of length n

    call ilu0((N-1)*(N-1), a, ja, ia, rW(ipaLU),iW(ipjLU),iW(ipjU),iW(ipiw),ierr)

    if (ierr .ne. 0) then
        write(*, '(A,I7)') 'initialization of ilu0 failed, ierr =', ierr
        goto 1002
    endif


    ! ======================================================================
    !  Set initial guess and compute initial residual
    ! ======================================================================
    ! set initial guess to 0
    call dcopy((N-1)*(N-1), 0d0, 0, c, 1)

    ! compute initial residual norm
    resinit = ddot((N-1)*(N-1), f, 1, f, 1)
    resinit = dsqrt(resinit)
    if (resinit .eq. 0d0) then
        write(*, '(A)') 'rhs=0, nothing to solve!'
        goto 1002
    endif

    ! ======================================================================	
    !  Iterative solution
    ! ======================================================================
    ITER = 1000              ! max number of iterations
    RESID = 1d-6 * resinit   ! threshold for \|RESID\|
    INFO  = 0                ! no troubles on input
    NUNIT = 6                ! output channel
    iprevec(1) = (N-1)*(N-1) ! single entry required: system size 
    imatvec(1) = (N-1)*(N-1) ! single entry required: system size 

    call slpbcgs(                    &
         prevec0, iprevec, iW,rW,	 &
         matvec,  imatvec, ia,ja,a,  &
         rW(ipBCG), (N-1)*(N-1), 8,  &
         (N-1)*(N-1), f, c,		     &
         ITER, RESID,                &
         INFO, NUNIT)
    if (INFO .ne. 0) then
        write(*, '(A)') 'BiCGStab failed'
        goto 1002
    endif

    ! Check error
    c_error_system = c_check_err(N, c)
    print *, 'error is ', c_error_system


    ! Print solve
    print *, 'C is:'
    do i = 1, N-1
        do j = 1, N-1
            print *, c((j-1)*(N-1) + i)
        enddo
    enddo	
    ! ======================================================================
    !  Deallocate area
    ! ======================================================================
    deallocate( ia )
    deallocate( ja )
    deallocate( a )
    deallocate( f, c )

    stop

    ! ======================================================================
    !  ERROR SECTION
    ! ======================================================================
    1000 continue
    write(*, '(A)') 'Cannot open file CSRsystem'
    stop 911
    1001 continue
    write(*, '(A)') 'Corrupted data in CSRsystem'
    stop 911
    1002 continue
    stop 911

end



