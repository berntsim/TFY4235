MODULE advection
use advectionparameters
IMPLICIT NONE
contains

SUBROUTINE initialize()
        integer :: i = 0
        do i = 1,N
                xPoints(i) = (i-1)*h
                u(i) = 1._wp +  dexp(-((xPoints(i)-xZero)**2)/(2._wp*sigma))
        end do
END SUBROUTINE

SUBROUTINE exact(exactTime)
        integer :: ie = 0
        real(kind=wp), intent(in) :: exactTime
        do ie = 1,N
                uExact(ie) =1._wp +  dexp(-((xPoints(ie) - xZero - c*exactTime)**2)/(2._wp*sigma))
        end do
END SUBROUTINE

SUBROUTINE LW()
        integer :: i_LW = 0
        do i_LW = 2,N-1
                uNext(i_LW) = u(i_LW) - beta*(u(i_LW+1) - u(i_LW-1)) + alpha*(u(i_LW+1) - 2*u(i_LW) + u(i_LW-1))
        end do
        u = uNext
END SUBROUTINE

SUBROUTINE LW_Burger()
        integer :: ib = 0
        do ib = 2,N-1
                uNext(ib) = u(ib) - kappa*(u(ib+1)*u(ib+1) - u(ib-1)*u(ib-1)) + iota*((u(ib+1) + u(ib))* &
                &(u(ib+1)*u(ib+1) - u(ib)*u(ib)) - (u(ib) + u(ib-1))*(u(ib)*u(ib) - u(ib-1)*u(ib-1)))
        end do
        u = uNext
        u(1) = 1._wp
END SUBROUTINE

SUBROUTINE plot_gnuplot(counter,timer)
        integer, parameter :: gnuplotter = 28
        integer, intent(in) :: counter  !counter for naming files
        real(kind=wp), intent(in) :: timer
        open (unit=gnuplotter, file = "plot_advection.gnu")
        write(gnuplotter,*) 'set terminal png size 600,500 enhanced font "Helvetica,10"'
        write(gnuplotter,'(A,i3.3,A)') 'set output "adfig/test',counter,'.png"'
        write(gnuplotter,*) 'set xlabel "Position"'
        write(gnuplotter,*) 'set title "Time = ',timer,' seconds"'
        write(gnuplotter,*) 'set xrange [0:1]'
        write(gnuplotter,*) 'set ylabel "Concentration"'
        write(gnuplotter,*) 'set yrange [0:1]'
        write(gnuplotter,*) 'plot "advection.dat" w l title "Numerical solution", "exact_advection.dat" w l title "Exact"'
        write(gnuplotter,*) 'plot "exact_advection.dat" with lines title "Exact solution"'
        Call SYSTEM('gnuplot -p "plot_advection.gnu"') 
        Call SYSTEM('rm plot_advection.gnu')
END SUBROUTINE

SUBROUTINE plot_burger(counter,timer)
        integer, parameter :: gnuplotter = 28
        integer, intent(in) :: counter  !counter for naming files
        real(kind=wp), intent(in) :: timer
        open (unit=gnuplotter, file = "plot_burger.gnu")
        write(gnuplotter,*) 'set terminal png size 600,500 enhanced font "Helvetica,10"'
        write(gnuplotter,'(A,i3.3,A)') 'set output "burgerfig/test',counter,'.png"'
        write(gnuplotter,*) 'set xlabel "Position"'
        write(gnuplotter,*) 'set title "Time = ',timer,' seconds"'
        write(gnuplotter,*) 'set xrange [0:1]'
        write(gnuplotter,*) 'set ylabel "Concentration"'
        write(gnuplotter,*) 'set yrange [0:2]'
        write(gnuplotter,*) 'plot "burger.dat" w l title "Hopf", "exact_advection.dat" w l title "Advection"'
        Call SYSTEM('gnuplot -p "plot_burger.gnu"') 
        Call SYSTEM('rm plot_burger.gnu')
END SUBROUTINE


SUBROUTINE advection_()
        integer :: i = 0
        integer :: j = 0
        integer :: pngcounter = 0
        real(kind=wp) :: time = 0._wp
        call initialize()
        print*, Nt/plotInterval        
        do i = 1,Nt
                call LW()
                call exact(time)
                if (mod(i,plotInterval) .eq. 0) then 
                        open (unit=out_unit, file="advection.dat")
                        open (unit=out_unit1, file="exact_advection.dat")
                        do j = 1,N
                                write(out_unit,*) xPoints(j), u(j)
                                write(out_unit1,*) xPoints(j), uExact(j)
                        end do
                        call plot_gnuplot(pngcounter,time)
                        pngcounter = pngcounter + 1
                        call SYSTEM('rm advection.dat')
                        call SYSTEM('rm exact_advection.dat')
                end if
                time = time + dt
        end do
END SUBROUTINE

SUBROUTINE inviscid_burger()
        integer :: i=0
        integer :: j=0
        integer :: pngcounter = 0
        real(kind=wp) :: time = 0._wp
        call initialize()
        print*, Nt/plotInterval        
        do i = 1,Nt
                call LW_Burger()
                call exact(time)
                if (mod(i,plotInterval) .eq. 0) then 
                        open (unit=out_unit, file="burger.dat")
                        open (unit=out_unit1, file="exact_advection.dat")
                        do j = 1,N
                                write(out_unit,*) xPoints(j), u(j)
                                write(out_unit1,*) xPoints(j), uExact(j)
                        end do
                        call plot_burger(pngcounter,time)
                        pngcounter = pngcounter + 1
                        call SYSTEM('rm burger.dat')
                        call SYSTEM('rm exact_advection.dat')
                end if
                time = time + dt
        end do    


END SUBROUTINE
END MODULE




