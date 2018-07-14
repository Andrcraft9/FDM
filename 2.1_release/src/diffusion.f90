! Diffusion equation in square area.
! Finite difference method.
! 
! 
! VERSION 1 
! 
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

    real*8, parameter :: eps = 100.0

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
        integer*4 :: i, j, k
        real*8 :: h, hh

        ! 
        h = 1.0 / N
        hh = h**2

        ! BOUNDARY POINTS (N2)
        ! Edge point (N1-N2)
        i = 1	
        j = 1
        k = (j-1)*(N-1) + i 

        IA(k) = (j-1)*(N-1) + i ! (i,j) point.

        JA(IA(k)) = (j-1)*(N-1) + i ! 1
        A(IA(k )) = (dx + dy) / hh

        JA(IA(k) + 1) = (j-1)*(N-1) + i + 1 ! 2
        A(IA(k)  + 1) = -dx / hh

        JA(IA(k) + 2) = (j-1 + 1)*(N-1) + i ! N
        A(IA(k)  + 2) = -dy / hh 
            
        F(k) = f_diff(i*h, j*h) + (g_N1(j*h) / h) + (g_N2(i*h) / h)			
        IA(k + 1) = IA(k) + 3			

        !print *, 'N1-N2 is OK'

        ! Inner points (N2)
        j = 1
        do i = 2, N-2
            k = (j-1)*(N-1) + i

            JA(IA(k)) = (j-1)*(N-1) + i - 1 ! i-1			
            A(IA(k) ) = -dx / hh
        
            JA(IA(k) + 1) = (j-1)*(N-1) + i ! i			
            A(IA(k)  + 1) = (2.0 * dx + dy) / hh	
    
            JA(IA(k) + 2) = (j-1)*(N-1) + i + 1 ! i+1
            A(IA(k)  + 2) = -dx / hh			

            JA(IA(k) + 3) = (j-1 + 1)*(N-1) + i ! (N-1) + i			
            A(IA(k)  + 3) = -dy / hh			

            F(k) = f_diff(i*h, j*h) + (g_N2(i*h) / h)
            IA(k + 1) = IA(k) + 4

            !print *, 'N2 is OK'
        enddo
        ! Edge point (N2-D2)
        i = N-1
        j = 1
        k = (j-1)*(N-1) + i

        JA(IA(k)) = (j-1)*(N-1) + i - 1 !N-2
        A(IA(k) ) = -dx / hh

        JA(IA(k) + 1) = (j-1)*(N-1) + i !N-1
        A(IA(k)  + 1) = (2.0*dx + dy) / hh

        JA(IA(k) + 2) = (j-1 + 1)*(N-1) + i !N-1 + N-1
        A(IA(k)  + 2) = -dy / hh
        
        F(k) = f_diff(i*h, j*h) + (g_N2(i*h) / h) + (dx * g_D2(j*h) / hh)
        IA(k + 1) = IA(k) + 3

        !print *, 'N2-D2 is OK'

        ! INNER AREA
        do j = 2, N-2
            ! Boundary points (N1)
            i = 1
            k = (j-1)*(N-1) + i

            JA(IA(k)) = (j-1 - 1)*(N-1) + i
            A(IA(k) ) = -dy / hh

            JA(IA(k) + 1) = (j-1)*(N-1) + i
            A(IA(k)  + 1) = (dx + 2.0*dy) / hh		

            JA(IA(k) + 2) = (j-1)*(N-1) + i + 1
            A(IA(k)  + 2) = -dx / hh

            JA(IA(k) + 3) = (j-1 + 1)*(N-1) + i
            A(IA(k)  + 3) = -dy / hh

            F(k) = f_diff(i*h, j*h) + (g_N1(j*h) / h)
            IA(k + 1) = IA(k) + 4

            !print *, 'N1 is OK'

            ! Inner points
            do i = 2, N-2
                k = (j-1)*(N-1) + i

                JA(IA(k)) = (j-1 - 1)*(N-1) + i
                A(IA(k) ) = -dy / hh

                JA(IA(k) + 1) = (j-1)*(N-1) + i - 1
                A(IA(k)  + 1) = -dx / hh 

                JA(IA(k) + 2) = (j-1)*(N-1) + i
                A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh

                JA(IA(k) + 3) = (j-1)*(N-1) + i + 1
                A(IA(k)  + 3) = -dx / hh 

                JA(IA(k) + 4) = (j-1 + 1)*(N-1) + i
                A(IA(k)  + 4) = -dy / hh
 
                F(k) = f_diff(i*h, j*h)
                IA(k + 1) = IA(k) + 5

                !print *, 'is OK'
            enddo

            ! Boundary points (D2)
            i = N-1
            k = (j-1)*(N-1) + i

            JA(IA(k)) = (j-1 - 1)*(N-1) + i
            A(IA(k) ) = -dy / hh

            JA(IA(k) + 1) = (j-1)*(N-1) + i - 1
            A(IA(k)  + 1) = -dx / hh

            JA(IA(k) + 2) = (j-1)*(N-1) + i
            A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh

            JA(IA(k) + 3) = (j-1 + 1)*(N-1) + i
            A(IA(k)  + 3) = -dy / hh

            F(k) = f_diff(i*h, j*h) + (dx * g_D2(j*h) / hh)
            IA(k + 1) = IA(k) + 4 

            !print *, 'D2 is OK'
        enddo		
        
        ! BOUNDARY POINTS (D1)
        ! Edge point (N1-D1)
        i = 1
        j = N-1
        k = (j-1)*(N-1) + i

        JA(IA(k)) = (j-1 - 1)*(N-1) + i
        A(IA(k) ) = -dy / hh

        JA(IA(k) + 1) = (j-1)*(N-1) + i
        A(IA(k)  + 1) = (dx + 2.0*dy) / hh

        JA(IA(k) + 2) = (j-1)*(N-1) + i + 1
        A(IA(k)  + 2) = -dx / hh

        F(k) = f_diff(i*h, j*h) + (g_N1(j*h) / h) + (dy * g_D1(i*h) / hh)
        IA(k + 1) = IA(k) + 3

        !print *, 'N1-D1 is OK'

        ! Inner points (D1)
        j = N-1
        do i = 2, N-2
            k = (j-1)*(N-1) + i

            JA(IA(k)) = (j-1 - 1)*(N-1) + i
            A(IA(k) ) = -dy / hh
            
            JA(IA(k) + 1) = (j-1)*(N-1) + i - 1			
            A(IA(k)  + 1) = -dx / hh
        
            JA(IA(k) + 2) = (j-1)*(N-1) + i
            A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh			

            JA(IA(k) + 3) = (j-1)*(N-1) + i + 1
            A(IA(k)  + 3) = -dx / hh			

            F(k) = f_diff(i*h, j*h) + (dy * g_D1(i*h) / hh)
            IA(k + 1) = IA(k) + 4

            !print *, 'D1 is OK'
        enddo
        ! Edge point (D1-D2)
        i = N-1
        j = N-1
        k = (j-1)*(N-1) + i

        JA(IA(k)) = (j-1 - 1)*(N-1) + i
        A(IA(k) ) = -dy / hh

        JA(IA(k) + 1) = (j-1)*(N-1) + i - 1
        A(IA(k)  + 1) = -dx / hh

        JA(IA(k) + 2) = (j-1)*(N-1) + i
        A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh
        
        F(k) = f_diff(i*h, j*h) + (dy * g_D1(i*h) / hh) + (dx * g_D2(j*h) / hh)
        IA(k + 1) = IA(k) + 3
        
        !print *, 'D1-D2 is OK'

        return
    end subroutine

    real*8 function solve(x, y)
        implicit none
        real*8 :: x, y

        solve = cos(pi * x) * cos(pi * y)
        return
    end function

    real*8 function c_check_err(N, c)
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

        c_check_err = mx
        return
    end function

end module
