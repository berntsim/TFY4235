MODULE wave
use waveparameters
IMPLICIT NONE
CONTAINS

SUBROUTINE initXY()
        integer :: i_init       
        do i_init = 1,N
                xPoints(i_init) = (i_init-1)*h
        end do
        yPoints = xPoints
END SUBROUTINE

SUBROUTINE initialize_system()
        integer :: i,j
        uPrev(1,:) = 0._wp
        uPrev(N,:) = 0._wp
        uPrev(:,N) = 0._wp
        uPrev(:,1) = 0._wp
        do i = 1,N
                do j = 1,N
                        if (water .eqv. .false.) then
                                uPrev(i,j)=dsin(xPoints(i)*pi)*dsin(yPoints(j)*2._wp*pi)                                                                  !
                        else
                                uPrev(i,j)=dexp(-((xPoints(i)-0.5_wp)**2+(yPoints(j)-0.5_wp)**2)/sigma)
                        end if
                end do        
        end do
        do i = 2,N-1
                do j = 2,N-1
                        u(i,j)=uPrev(i,j)*(1-2*alpha)+(alpha/2._wp)*(uPrev(i+1,j)+uPrev(i-1,j)+uPrev(i,j+1)+uPrev(i,j-1))!Using the reflective BC to calculate u(n=-1)
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
        integer:: i2,j2
        do i2 = 2,N-1
                do j2 = 2,N-1
                        uNext(i2,j2) = u(i2,j2)*(2._wp-4._wp*alpha)-uPrev(i2,j2)+alpha*(u(i2+1,j2)+u(i2-1,j2)+u(i2,j2+1)+u(i2,j2-1)) 
                end do 
        end do
        do i2 = 2,N-1
                do j2 = 2,N-1
                        uPrev(i2,j2) = u(i2,j2)
                        u(i2,j2) = uNext(i2,j2)
                end do
        end do
END SUBROUTINE

SUBROUTINE plot_gnuplot(counter)
        integer, parameter :: gnuplotter = 28
        integer, intent(in) :: counter  !counter for naming files
        open (unit=gnuplotter, file = "plot_test.gnu")
        write(gnuplotter,*) 'set terminal png size 600,500 enhanced font "Helvetica,10"'
        write(gnuplotter,'(A,i3.3,A)') 'set output "fig/test',counter,'.png"'
        write(gnuplotter,*) 'set xlabel "x-points"'
        write(gnuplotter,*) 'set xrange [0:1]'
        write(gnuplotter,*) 'set ylabel "y-points"'
        write(gnuplotter,*) 'set yrange [0:1]'
        write(gnuplotter,*) 'set zrange[-1:1]'
        write(gnuplotter,*) 'set pm3d'
        write(gnuplotter,*) 'set cbrange [-0.15:0.15]'
         write(gnuplotter,*) 'splot "wave.dat" with dots '
        Call SYSTEM('gnuplot -p "plot_test.gnu"') 
        Call SYSTEM('rm plot_test.gnu')
END SUBROUTINE




SUBROUTINE wave_equation()
        integer :: j_we,i1,j1, pngCounter = 0
        call initXY()
        call initialize_system()
        print*, Nt/plotInterval
        do i1 = 1,Nt
                call iterate_time()
                if (mod(i1,plotInterval) .eq. 0) then
                        open (unit=out_unit, file="wave.dat") 
                        do j_we = 1,N
                                do j1 = 1,N
                                        write(out_unit,*) xPoints(j_we), yPoints(j1), u(j_we,j1)
                                end do
                                write(out_unit,*)
                        end do 
                        call plot_gnuplot(pngCounter)
                        pngCounter = pngCounter + 1
                        Call SYSTEM('rm wave.dat')
                end if
        end do
END SUBROUTINE

END MODULE
