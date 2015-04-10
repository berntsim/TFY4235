MODULE parameters
implicit none
integer, parameter :: wp=kind(0.d0)
real(kind=wp),parameter		:: pi = 3.141592653589793238_wp, alpha = 0.2_wp, deltaT = 0.0001_wp
real(kind=wp),parameter		:: L = 0.000020_wp, KbT = 0.026_wp, startX = 0._wp, startT = 0._wp, r1 = 0.000000012_wp, eta = 0.001_wp
real(kind=wp)			:: D = 0._wp, omega = 1._wp, leftSideStep = 0._wp, gamma_1 = 1._wp, gamma_2 = 1._wp, deltaU  = 80.0_wp
real(kind=wp)			:: dc = 0._wp, r2 = 3.0_wp*r1, veloSum = 0.0_wp, veloAvg = 0.0_wp, error = 0.0_wp, tauStep = 0.5_wp, tau = 0.5_wp
real(kind=wp)			:: s = 0.0_wp, errorSum = 0.0_wp, const = 0.0_wp
integer, parameter 		:: out_unit=29, N = 10000000, out_unit1 = 28, out_unit2 = 27, out_unit3 = 26, out_unit4 = 25, out_unit5 = 24
integer, parameter		:: out_unit6 = 23
integer, parameter		:: particles = 100, tauIterations = 1, pickInterval = 1000000
real(kind=wp),allocatable	:: xArray(:), tArray(:), randNumbers(:), potArray(:), velocities(:), errors(:), particlePos(:)
real(kind=wp),dimension(2)	:: testMaxForce
logical				:: flashing = .false.				!Used for task 8 and onwards.
logical				:: findAvgVelocity = .false.			!Used for task 9 and onwards.
logical				:: ensemble = .true.				!Used for task 12 and onwards.

END MODULE

!To check for Boltzmann distribution, make sure deltaU = 10*KbT, flashing = false, findAvgVelocity = false, N = 10000000, deltaT = 0.000001_wp
!Max N able to run using allocate is N = 90000000



!Med const: 50 partikler t = 1 min, 17 sekunder og 4 tideler
!Uten const: 50 partikler t = 1 min, 21 sekunder og 7 tideler
