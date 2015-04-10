program main
use parameters
use random
implicit none
!--------------------------- Parameters ----------------------------------
real(kind=wp)	:: Ur, f, potential, force, euler, Pu, tmp = 0._wp
real(kind=wp),dimension(100000)	:: test_gaussian
integer	:: i = 0, z = 0, j = 0

!--------------------- Random number generator seeding --------------------------
integer :: values(1:8), k
integer, dimension(:), allocatable :: seed
call date_and_time(values=values)
call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)

!-------------------- Changing to dimensionless units ---------------------------
gamma_1 = 6._wp*pi*eta*r1							!The main point of this part is to initialize and make units reduced
gamma_2 = 6._wp*pi*eta*r2

allocate (xArray(N))
allocate (tArray(N))
allocate (potArray(N))
allocate (randNumbers(N))
allocate (velocities(particles))
allocate (errors(tauIterations))
allocate (particlePos(9*particles))						!This 9 is because we check the particle's position at 9 different times

randNumbers(1) = random_normal()
xArray(1) = startX
tArray(1) = startT
xArray = xArray/L
omega = (1.60217657_wp*10._wp**(-19))*deltaU/(gamma_1*L**2)
tArray = tArray*omega
D = KbT/deltaU
tau = tau*omega
const = sqrt(2._wp*D*deltaT)							!Just making sure that we don't spend computational power on calculating sqrt all the time.
!print*, omega
! ---------- Testing section -------------------------
!open (unit=out_unit2,file="potentialtest.dat",action="write",status="replace")
!open (unit=out_unit3,file="forcetest.dat",action="write",status="replace")
!do i = 1,N-1									!Testing if the potential looks like it should, by generating data and plotting it
!	tArray(i+1) = tArray(i) + deltaT					!using gnuplot for instance.
!	xArray(i+1) = potential(tArray(i), xArray(i))
!	if (MODULO(i,1000) .eq. 0) then
!		write(out_unit3,*) tArray(i), force(tArray(i))
!		write(out_unit2,*) tArray(i), xArray(i)
!	end if	
!end do

!open (unit=out_unit,file="test.dat",action="write",status="replace")		!Write random numbers to a file, so that one can check if they are are really random,
!do i = 1,100000								!normally distributed numbers. The actual plotting of this was done in python.
!test_gaussian(i) = random_normal()						!See random.f95 for module including random_normal(). This is not box-müller, but the other.		
!	write(out_unit,*) test_gaussian(i)
!end do
close(out_unit)


!------------------ Determining the step size ----------------------------------
testMaxForce(1) = 1._wp/alpha							!This is the requirement stated in the assignment that our time step will not jump through several 
testMaxForce(2) = 1._wp/(1._wp - alpha)						!variations of the potential. It has been modified to be dimensionless
leftSideStep = MAXVAL(testMaxForce)*deltaT + 4*sqrt(2*D*deltaT)

if (leftSideStep .gt. (alpha/10._wp)) then
	write(*,*) "Timestep is too large, choose another"
	return
end if

!----------------------Euler scheme and iterations ------------------------------
if (findAvgVelocity .eqv. .true.) then							!Checking wether or not we want to look at one particle, or many.
	open (unit=out_unit5,file="avgVelocity_r2_tauOP1.dat",action="write",status="replace")	
	write(out_unit5,*) 0.0_wp,0.0_wp, 0.0_wp	
	do j = 1,tauIterations
		tau = tauStep*j*omega							!Iterates over tau, so that we may find a tau_optimal
		do z = 1,particles							!Looks at #particles amount of particles, and find their average speed
			do i = 1,N-1							!Do the usual euler scheme
				randNumbers(i+1) = random_normal()							
				xArray(i+1) = euler(xArray(i),tArray(i),randNumbers(i))
				tArray(i+1) = tArray(i) + deltaT
			end do				
			velocities(z) = (L*xArray(N))/(tArray(N)/omega)
			veloSum = veloSum + velocities(z)				!Finding the sum of the velocities to average later
			xArray(1) = startX/L						!Reset where to start in time and space
			tArray(1) = startT*omega
			print*, z					
		end do

		do i = 1,particles
		errorSum = errorSum + (velocities(i) - (veloSum/particles))**2
		end do

		s = sqrt((1.0_wp/(particles-1.0_wp))*errorSum)
		errors(j) = s
		write(out_unit5,*)  tau/omega, (veloSum/particles), errors(j)			
		errorSum = 0.0_wp
		veloSum = 0.0_wp							!Reset the sum of the velocities
		write (*,*) "tau = ", tau/omega						!Print tau so user knows ~what stage the program is at
	end do
else if (ensemble .eqv. .true.) then
	open (unit=out_unit6,file="task12.dat",action="write",status="replace")
	write(out_unit6,*) deltaT							!Writing timesteps and iterations, so that one may change them here and the python 
	write(out_unit6,*) N								!code will still work.
	write(out_unit6,*) omega
	write(out_unit6,*) pickInterval
	write(out_unit6,*) particles
	do z = 1,particles								!Looks at #particles amount of particles for calculating ensembles
		do i = 1,N-1								!Do the usual euler scheme
			randNumbers(i+1) = random_normal()							
			xArray(i+1) = euler(xArray(i),tArray(i),randNumbers(i))
			tArray(i+1) = tArray(i) + deltaT
			if (MODULO(i,pickInterval) .eq. 0) then
				write(out_unit6,*) xArray(i)*L
			end if
		end do
		xArray(1) = startX/L							!Reset where to start in time and space
		tArray(1) = startT*omega
	print*, z					
	end do
else
	open (unit=out_unit1,file="walkerTask6.dat",action="write",status="replace")	!Open file to which the position data will be written
	open (unit=out_unit4,file="Boltzmanndata.dat",action="write",status="replace")	!Open file to which the potential data will be written
	do i = 1,N-1									
		randNumbers(i+1) = random_normal()					! The implementation of this euler scheme can be found in euler.f95	
		xArray(i+1) = euler(xArray(i),tArray(i),randNumbers(i))
		tArray(i+1) = tArray(i) + deltaT
		potArray(i) = potential(xArray(i), tArray(i))
		if (MODULO(i,1000) .eq. 0) then
			write(out_unit1,*) tArray(i)/omega, xArray(i)*L			!Write to file
			write(out_unit4,*) potArray(i) 
		end if
	end do
end if

end program



