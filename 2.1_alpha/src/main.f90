! Diffusion equation in square area.
! Finite difference method.
! 
!----------  AREA:  ----------!
! y ^
!   |      D1
!  1|-----------|
!   |           |
!   |           |
! N1|		    | D2
!   |           |
!   |-----------|-->
!  0      N2    1  x
! 
! Generate: 
! Dv = f,
! v = [c_1:1 c_2:1, ... , c_N-1:1, c_1:2, ... , c_N-1:2, c_1:3, ... ], 
! c_i:j = c(x_i, y_j)
!
module diffusion_equ
	implicit none

	real*8, parameter :: pi = 4.0d0 * datan(1.0d0)
	real*8, parameter :: pi2 = pi**2

	real*8, parameter :: eps = 1.0

	real*8, parameter :: dx = 1.0
	real*8, parameter :: dy = eps

	contains

	real*8 function f_diff(x, y)
		implicit none
		real*8 :: x, y

		f_diff = (1 + eps) * cos(pi * x) * cos(pi * y) * pi2
		return
	end function

	real*8 function g_D1(x)
		implicit none
		real*8 :: x

		g_D1 = -cos(pi * x)
		return
	end function


	real*8 function g_D2(y)
		implicit none
		real*8 :: y

		g_D2 = -cos(pi * y)
		return
	end function

	real*8 function g_N1(y)
		implicit none
		real*8 :: y

		g_N1 = 0.0
		return
	end function

	real*8 function g_N2(x)
		implicit none
		real*8 :: x

		g_N2 = 0.0
		return
	end function

	subroutine generate_matrix(N, IA, JA, A, F)
        implicit none
        integer*4 :: N
		integer*4, intent(in out) :: IA(:), JA(:)
        real*8, intent(in out) :: A(:), F(:)
        ! Local
        integer*4 :: i, j
        real*8 :: h, hh

		! 
		h = 1.0 / N
		hh = h**2

		! BOUNDARY POINTS (N2)
		! Edge point (N1-N2)
		i = 1	
		j = 1
		IA( (j-1)*(N-1) + i ) = (j-1)*(N-1) + i ! (i,j) point.

		JA(IA( (j-1)*(N-1) + i )) = (j-1)*(N-1) + i ! 1
		A(IA( (j-1)*(N-1) + i )) = (dx + dy) / hh

		JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i + 1 ! 2
		A(IA( (j-1)*(N-1) + i )  + 1) = -dx / hh

		JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1 + 1)*(N-1) + i ! N
		A(IA( (j-1)*(N-1) + i )  + 2) = -dy / hh 
		
		F((j-1)*(N-1) + i) = f_diff(i*h, j*h) - (dx * g_N1(j*h) / h) - (dy * g_N2(i*h) / h)
		IA( (j-1)*(N-1) + i + 1 ) = IA( (j-1)*(N-1) + i ) + 3

		print *, 'N1-N2 is OK'

		! Inner points (N2)
		j = 1
		do i = 2, N-2
			JA(IA( (j-1)*(N-1) + i ) ) = (j-1)*(N-1) + i - 1 ! i-1			
			A(IA( (j-1)*(N-1) + i )  ) = -dx / hh
			
			JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i ! i			
			A(IA( (j-1)*(N-1) + i )  + 1) = (2.0 * dx + dy) / hh	
		
			JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i + 1 ! i+1
			A(IA( (j-1)*(N-1) + i )  + 2) = -dx / hh			

			JA(IA( (j-1)*(N-1) + i ) + 3) = (j-1 + 1)*(N-1) + i ! (N-1) + i			
			A(IA( (j-1)*(N-1) + i )  + 3) = -dy / hh			

			F((j-1)*(N-1) + i) = f_diff(i*h, j*h) - (dy * g_N2(i*h) / h)
			IA( (j-1)*(N-1) + i + 1 ) = IA( (j-1)*(N-1) + i ) + 4

			print *, 'N2 is OK'
		enddo
		! Edge point (N2-D2)
		i = N-1
		j = 1
		
		JA(IA( (j-1)*(N-1) + i ) ) = (j-1)*(N-1) + i - 1 !N-2
		A(IA( (j-1)*(N-1) + i )  ) = -dx / hh

		JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i !N-1
		A(IA( (j-1)*(N-1) + i )  + 1) = (2.0*dx + dy) / hh

		JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1 + 1)*(N-1) + i !N-1 + N-1
		A(IA( (j-1)*(N-1) + i )  + 2) = -dy / hh
		
		F((j-1)*(N-1) + i) = f_diff(i*h, j*h) - (dy * g_N2(i*h) / h) + (dx * g_D2(j*h) / hh)
        IA( (j-1)*(N-1) + i + 1 ) = IA( (j-1)*(N-1) + i ) + 3

		print *, 'N2-D2 is OK'

		! INNER AREA
		do j = 2, N-2
			! Boundary points (N1)
			i = 1
			JA(IA( (j-1)*(N-1) + i ) ) = (j-1 - 1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  ) = -dy / hh

			JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  + 1) = (dx + 2.0*dy) / hh		

			JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i + 1
			A(IA( (j-1)*(N-1) + i )  + 2) = -dx / hh

			JA(IA( (j-1)*(N-1) + i ) + 3) = (j-1 + 1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  + 3) = -dy / hh

			F((j-1)*(N-1) + i) = f_diff(i*h, j*h) - (dx * g_N1(j*h) / h)
			IA( (j-1)*(N-1) + i + 1 ) = IA( (j-1)*(N-1) + i ) + 4

			print *, 'N1 is OK'

			! Inner points
			do i = 2, N-2
				JA(IA( (j-1)*(N-1) + i ) ) = (j-1 - 1)*(N-1) + i
				A(IA( (j-1)*(N-1) + i )  ) = -dy / hh

				JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i - 1
				A(IA( (j-1)*(N-1) + i )  + 1) = -dx / hh 

				JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i
				A(IA( (j-1)*(N-1) + i )  + 2) = (2.0*dx + 2.0*dy) / hh

				JA(IA( (j-1)*(N-1) + i ) + 3) = (j-1)*(N-1) + i + 1
				A(IA( (j-1)*(N-1) + i )  + 3) = -dx / hh 

				JA(IA( (j-1)*(N-1) + i ) + 4) = (j-1 + 1)*(N-1) + i
				A(IA( (j-1)*(N-1) + i )  + 4) = -dy / hh
 
				F((j-1)*(N-1) + i) = f_diff(i*h, j*h)
				IA( (j-1)*(N-1) + i + 1) = IA( (j-1)*(N-1) + i ) + 5

				print *, 'is OK'
			enddo

			! Boundary points (D2)
			i = N-1
			JA(IA( (j-1)*(N-1) + i ) ) = (j-1 - 1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  ) = -dy / hh

			JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i - 1
			A(IA( (j-1)*(N-1) + i )  + 1) = -dx / hh

			JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  + 2) = (2.0*dx + 2.0*dy) / hh

			JA(IA( (j-1)*(N-1) + i ) + 3) = (j-1 + 1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  + 3) = -dy / hh

			F((j-1)*(N-1) + i) = f_diff(i*h, j*h) + (dx * g_D2(j*h) / hh)
		    IA( (j-1)*(N-1) + i + 1) = IA( (j-1)*(N-1) + i ) + 4 

			print *, 'D2 is OK'
		enddo		
		
		! BOUNDARY POINTS (D1)
		! Edge point (N1-D1)
		i = 1
		j = N-1
		
		JA(IA( (j-1)*(N-1) + i ) ) = (j-1 - 1)*(N-1) + i
		A(IA( (j-1)*(N-1) + i )  ) = -dy / hh

		JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i
		A(IA( (j-1)*(N-1) + i )  + 1) = (dx + 2.0*dy) / hh

		JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i + 1
		A(IA( (j-1)*(N-1) + i )  + 2) = -dx / hh

		F((j-1)*(N-1) + i) = f_diff(i*h, j*h) - (dx * g_N1(j*h) / h) + (dy * g_D1(i*h) / hh)
		IA( (j-1)*(N-1) + i + 1) = IA( (j-1)*(N-1) + i ) + 3

		print *, 'N1-D1 is OK'

		! Inner points (D1)
		j = N-1
		do i = 2, N-2
			JA(IA( (j-1)*(N-1) + i ) ) = (j-1 - 1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  ) = -dy / hh
			
			JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i - 1			
			A(IA( (j-1)*(N-1) + i )  + 1) = -dx / hh	
		
			JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i
			A(IA( (j-1)*(N-1) + i )  + 2) = (2.0*dx + 2.0*dy) / hh			

			JA(IA( (j-1)*(N-1) + i ) + 3) = (j-1)*(N-1) + i + 1
			A(IA( (j-1)*(N-1) + i )  + 3) = -dx / hh			

			F((j-1)*(N-1) + i) = f_diff(i*h, j*h) + (dy * g_D1(i*h) / hh)
			IA( (j-1)*(N-1) + i + 1 ) = IA( (j-1)*(N-1) + i ) + 4

			print *, 'D1 is OK'
		enddo
		! Edge point (D1-D2)
		i = N-1
		j = N-1
		
		JA(IA( (j-1)*(N-1) + i ) ) = (j-1 - 1)*(N-1) + i
		A(IA( (j-1)*(N-1) + i )  ) = -dy / hh

		JA(IA( (j-1)*(N-1) + i ) + 1) = (j-1)*(N-1) + i - 1
		A(IA( (j-1)*(N-1) + i )  + 1) = -dx / hh

		JA(IA( (j-1)*(N-1) + i ) + 2) = (j-1)*(N-1) + i
		A(IA( (j-1)*(N-1) + i )  + 2) = (2.0*dx + 2.0*dy) / hh
		
		F((j-1)*(N-1) + i) = f_diff(i*h, j*h) + (dy * g_D1(i*h) / hh) + (dx * g_D2(j*h) / hh)
        IA( (j-1)*(N-1) + i + 1 ) = IA( (j-1)*(N-1) + i ) + 3
		
		print *, 'D1-D2 is OK'

		return
	end subroutine

	real*8 function solve(x, y)
		implicit none
		real*8 :: x, y

		solve = cos(pi * x) * cos(pi * y)
		return
	end function

	real*8 function simple_check_err(N, c)
		implicit none
		integer*4 :: N
		real*8, intent(in out) :: c(:)

		integer*4 :: i, j
		real*8 :: h, er, mx

		h = 1.0 / N
		er = 0.0d0
		mx = 0.0d0
		do j = 1, N-1
			do i = 1, N-1
				! print *, 'c = ', c((j-1)*(N-1) + i)
				! print *, 's = ', solve(i*h, j*h)
				er = abs( c((j-1)*(N-1) + i) - solve(i*h, j*h) ) 
				if (er > mx) then
					mx = er
				endif
			enddo
		enddo

		simple_check_err = er
		return
	end function

end module

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
	real*8 :: error_system

	N = 256
	
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

	error_system = simple_check_err(N, c)
	print *, 'error is ', error_system

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



