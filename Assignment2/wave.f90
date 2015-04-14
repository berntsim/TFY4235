MODULE wave
use waveparameters
IMPLICIT NONE
CONTAINS

SUBROUTINE initXY()
        do i = 1,N
                xPoints(i) = i*h
        end do
        yPoints = xPoints
END SUBROUTINE

SUBROUTINE calculate_sines()
        do i = 1,N
                sinXVector(i) = dsin(pi*xPoints(i))
                sinYVector(i) = dsin(2._wp*pi*yPoints(i))
        end do
END SUBROUTINE

SUBROUTINE initialize_system()
        do i = 2,N-1
                do j = 2,N-1
                        uPrev(i,j)=sinXVector(i)*sinYVector(i)                                                                  !
                end do        
        end do
        do i = 2,N-1
                do j = 2,N-1
                        u(i,j)=uPrev(i,j)*(1-2*alpha)-(alpha/2._wp)*(uPrev(i+1,j)+uPrev(i-1,j)+uPrev(i,j+1)+uPrev(i,j-1))!Using the reflective BC to calculate u(n=-1)
                end do
                u(i,1) = 0._wp           !Setting boundary conditions
                u(i,N) = 0._wp
                u(1,i) = 0._wp
                u(N,i) = 0._wp

                uNext(i,1) = 0._wp
                uNext(i,N) = 0._wp
                uNext(1,i) = 0._wp
                uNext(N,i) = 0._wp
        end do
END SUBROUTINE

SUBROUTINE iterate_time()
        do i = 2,N-1
                do j = 2,N-1
                        uNext(i,j) = u(i,j)*(2._wp-4._wp*alpha)-uPrev(i,j)+alpha*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)) 
                end do 
        end do
        do i = 2,N-1
                do j = 2,N-1
                        uPrev(i,j) = u(i,j)
                        u(i,j) = uNext(i,j)
                end do
        end do
END SUBROUTINE

SUBROUTINE wave_equation()
        call initXY()
        call calculate_sines()
        call initialize_system()
        call iterate_time()
        !print*, sinXvector
        !print*, sinYVector
        !print*,
        !print*, u
END SUBROUTINE

END MODULE
