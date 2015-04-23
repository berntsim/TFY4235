MODULE advectionparameters
integer, parameter                      :: wp=kind(0.d0)        !Set working precission
integer, parameter                      :: N = 1000               !Setting number of iterations for the program to run
real(kind=wp), parameter                :: r = 0.25_wp          !For the stability of euler
real(kind=wp), parameter                :: t = 0.1_wp           !Setting simulation time
real(kind=wp), parameter                :: c = 1._wp            !setting speed for the transport
real(kind=wp), dimension(N)             :: xPoints              !x-values used in the problem
real(kind=wp), parameter                :: xZero = 0.4_wp       !Setting the point where the initial concentration will be centered
real(kind=wp), dimension(N)             :: u                    !The vector in which we will store the concentration
real(kind=wp), dimension(N)             :: uNext                !For updating in the scheme
real(kind=wp), dimension(N)             :: uExact               !For calculating the exact solution
real(kind=wp), parameter                :: sigma = 0.001        !Setting the standard deviation of the gaussian used as initial concentration
real(kind=wp), parameter                :: L = 2.0_wp           !Setting length of the system
integer, parameter                      :: plotInterval = 1    !Setting how often the program should plot the data
integer, parameter                      :: out_unit = 29        !Out unit for numerical data
integer, parameter                      :: out_unit1 = 28       !Out unit for exact solution

!--------------------------Section for calculating dt and Nt----------------

real(kind=wp), parameter                :: h = L/(N-1)          !Setting space steps
real(kind=wp), parameter                :: dt = dsqrt(r*h*h)    !Calculating the time step size
real(kind=wp), parameter                :: alpha = c*c*r/2._wp  !Constant used in the scheme
real(kind=wp), parameter                :: beta = c*dt/(2._wp*h)!Setting the other constant used in the scheme
real(kind=wp), parameter                :: kappa = dt/(4._wp*h) !Setting the first constant used in the LW_Burger scheme
real(kind=wp), parameter                :: iota = dt*dt/(8._wp*h*h)      !Setting the second constant used in the LW_Burger scheme
integer, parameter                      :: Nt = nint(t/dt)      !Setting time iterations to run according to simulation time
END MODULE
