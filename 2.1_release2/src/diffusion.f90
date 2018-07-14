! Diffusion equation in square area.
! Finite difference method.
! 
! 
! VERSION 2
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
! v = [c_0:0, c_1:0, ... , c_N:0, c_0:1, ... , c_N:1, c_0:2, ... ], 
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
        j = 0
        IA(1) = 1
        do i = 0, N-1
            k = N*j + i + 1
            JA(IA(k)) = k
            A(IA(k) ) = -dy / hh

            JA(IA(k) + 1) = k + N
            A(IA(k)  + 1) = dy / hh
            
            F(k) = g_N2(i*h) / h
            IA(k + 1) = IA(k) + 2 
        enddo	
        
        ! INNER AREA
        do j = 1, N-2
            ! Boundary points (N1)
            i = 0
            k = N*j + i + 1

            JA(IA(k)) = k
            A(IA(k) ) = -dx / hh
            
            JA(IA(k) + 1) = k + 1
            A(IA(k)  + 1) = dx / hh
            
            F(k) = g_N1(j*h) / h
            IA(k + 1) = IA(k) + 2 

            ! Inner points
            do i = 1, N-2
                k = N*j + i + 1

                JA(IA(k)) = k - N
                A(IA(k) ) = -dy / hh

                JA(IA(k) + 1) = k - 1
                A(IA(k)  + 1) = -dx / hh 

                JA(IA(k) + 2) = k
                A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh

                JA(IA(k) + 3) = k + 1
                A(IA(k)  + 3) = -dx / hh 

                JA(IA(k) + 4) = k + N
                A(IA(k)  + 4) = -dy / hh
 
                F(k) = f_diff(i*h, j*h)
                IA(k + 1) = IA(k) + 5
            enddo

            ! Boundary points (D2)
            i = N - 1
            k = N*j + i + 1

            JA(IA(k)) = k - N
            A(IA(k) ) = -dy / hh

            JA(IA(k) + 1) = k - 1
            A(IA(k)  + 1) = -dx / hh

            JA(IA(k) + 2) = k
            A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh

            JA(IA(k) + 3) = k + N
            A(IA(k)  + 3) = -dy / hh

            F(k) = f_diff(i*h, j*h) + (dx * g_D2(j*h) / hh)
            IA(k + 1) = IA(k) + 4 

            !print *, 'D2 is OK'
        enddo		

        ! Boundary points (D2)
        j = N-1

        ! Boundary points (N1)
        i = 0
        k = N*j + i + 1

        JA(IA(k)) = k
        A(IA(k) ) = -dx / hh
        
        JA(IA(k) + 1) = k + 1
        A(IA(k)  + 1) = dx / hh
        
        F(k) = g_N1(j*h) / h
        IA(k + 1) = IA(k) + 2 

        do i = 1, N-2
            k = N*j + i + 1

            JA(IA(k)) = k - N
            A(IA(k) ) = -dy / hh
            
            JA(IA(k) + 1) = k - 1
            A(IA(k)  + 1) = -dx / hh
        
            JA(IA(k) + 2) = k
            A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh			

            JA(IA(k) + 3) = k + 1
            A(IA(k)  + 3) = -dx / hh			

            F(k) = f_diff(i*h, j*h) + (dy * g_D1(i*h) / hh)
            IA(k + 1) = IA(k) + 4
        enddo

        ! Edge point (D1-D2)
        i = N-1
        k = N*j + i + 1

        JA(IA(k)) = k - N
        A(IA(k) ) = -dy / hh

        JA(IA(k) + 1) = k - 1
        A(IA(k)  + 1) = -dx / hh

        JA(IA(k) + 2) = k
        A(IA(k)  + 2) = (2.0*dx + 2.0*dy) / hh
        
        F(k) = f_diff(i*h, j*h) + (dy * g_D1(i*h) / hh) + (dx * g_D2(j*h) / hh)
        IA(k + 1) = IA(k) + 3

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
        do j = 0, N-1
            do i = 0, N-1
                ! print *, 'c = ', c((j-1)*(N-1) + i)
                ! print *, 's = ', solve(i*h, j*h)
                er = abs( c(N*j + i + 1) - solve(i*h, j*h) ) 
                if (er > mx) then
                    mx = er
                endif
            enddo
        enddo

        c_check_err = mx
        return
    end function

    real*8 function L2_check_err(N, c)
        implicit none
        integer*4 :: N
        real*8, intent(in out) :: c(:)

        real*8 :: A(3, 9)
        integer*4 :: i, j, k, m
        real*8 :: h, S, T
        real*8 :: f1, f2, f3, f4, c_m
        real*8 :: x1, x2, x3, y1, y2, y3, x, y

        real*8 :: a11, a12, a21, a22, a23
        real*8 :: w1, w2
        
        a11 = 0.12495d0
        a12 = 0.43752d0
        a21 = 0.79711d0

        a22 = 0.16541d0
        a23 = 0.03748d0

        w1 = 0.20595d0
        w2 = 0.06369d0

        ! Init A
        A(1,1) = a11 
        A(2,1) = a12
        A(3,1) = a12

        A(1,2) = a12 
        A(2,2) = a11
        A(3,2) = a12

        A(1,3) = a12 
        A(2,3) = a12
        A(3,3) = a11

        A(1,4) = a21 
        A(2,4) = a22
        A(3,4) = a23

        A(1,5) = a21 
        A(2,5) = a23
        A(3,5) = a22

        A(1,6) = a22 
        A(2,6) = a21
        A(3,6) = a23

        A(1,7) = a22 
        A(2,7) = a23
        A(3,7) = a21

        A(1,8) = a23 
        A(2,8) = a21
        A(3,8) = a22

        A(1,9) = a23 
        A(2,9) = a22
        A(3,9) = a21

        ! Program
        h = 1.0 / N
        S = 0.0d0

        ! Compute integral
        T = 0.5d0 * h * h
        do j = 0, N-1
            do i = 0, N-1
                k = N*j + i + 1												

                ! Edge point
                f1 = c(k)

                if (i < (N-1)) then
                    f2 = c(k+1)
                else
                    f2 = g_D2(j*h)
                endif

                if (j < (N-1)) then
                    f3 = c(k+N)
                else
                    f3 = g_D2(i*h)
                endif

                if (i < (N-1)) then
                    if (j < (N-1)) then
                        f4 = c(k+N+1)
                    else
                        f4 = g_D1((i+1)*h) 
                    endif
                else
                    f4 = g_D2((j+1)*h)
                endif			
                
                x1 = i*h
                y1 = j*h
                x3 = x1 + h
                y3 = y1 + h				

                ! Down-triangle
                x2 = x1 + h
                y2 = y1
                do m = 1, 3
                    x = x1*A(1, m) + x2*A(2, m) + x3*A(3, m)					
                    y = y1*A(1, m) + y2*A(2, m) + y3*A(3, m)

                    c_m = f1 + (f2 - f1)*(x - x1)/h + (f4 - f2)*(y - y1)/h
                    S = S + T*w1*( (c_m - solve(x,y))**2 )				
                enddo

                do m = 4, 9
                    x = x1*A(1, m) + x2*A(2, m) + x3*A(3, m)					
                    y = y1*A(1, m) + y2*A(2, m) + y3*A(3, m)
                    
                    c_m = f1 + (f2 - f1)*(x - x1)/h + (f4 - f2)*(y - y1)/h
                    S = S + T*w2*( (c_m - solve(x,y))**2 )
                enddo

                ! Up-triangle
                x2 = x1
                y2 = y1 + h
                do m = 1, 3
                    x = x1*A(1, m) + x2*A(2, m) + x3*A(3, m)					
                    y = y1*A(1, m) + y2*A(2, m) + y3*A(3, m)

                    c_m = f1 + (f4 - f3)*(x - x1)/h + (f3 - f1)*(y - y1)/h
                    S = S + T*w1*( (c_m - solve(x,y))**2 )				
                enddo
                do m = 4, 9
                    x = x1*A(1, m) + x2*A(2, m) + x3*A(3, m)					
                    y = y1*A(1, m) + y2*A(2, m) + y3*A(3, m)

                    c_m = f1 + (f4 - f3)*(x - x1)/h + (f3 - f1)*(y - y1)/h
                    S = S + T*w2*( (c_m - solve(x,y))**2 )
                enddo

            enddo
        enddo
            
        L2_check_err = sqrt(S)

        return
    end function

end module
