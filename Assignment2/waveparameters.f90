MODULE waveparameters
integer, parameter                      :: wp=kind(0.d0)        !Set working precission
real(kind=wp), parameter                :: c = 1._wp            !Set speed at which the wave travels
integer, parameter                      :: N = 5               !Set number of grid points in the x and y-direction
real(kind=wp), parameter                :: t = 1._wp            !Set simulation time
real(kind=wp), dimension(N,N)           :: u                    !Initializing current vector
real(kind=wp), dimension(N,N)           :: uPrev                !Initializing previous vector
real(kind=wp), dimension(N,N)           :: uNext                !Initializing next vector
real(kind=wp), dimension(N)             :: xPoints              !Initializing x-points
real(kind=wp), dimension(N)             :: yPoints              !Initializing y-points
real(kind=wp), dimension(N)             :: sinXVector           !Making a vector for calculating the sin values in x-direction
real(kind=wp), dimension(N)             :: sinYVector           !Same, but in y-direction
real(kind=wp), parameter                :: pi = 3.141592653589793238_wp
real(kind=wp), parameter                :: r = 0.25             !Setting parameter to maintain stability criteria
real(kind=wp), parameter                :: L = 1._wp            !Setting system lenght. It is such that Lx = Ly = L
integer                                 :: i,j                  !Iteration variables

!------------Calculating dt and NT-------------

real(kind=wp), parameter                :: h = L/(N+1)          !Setting step size in the xy-grid
real(kind=wp), parameter                :: dt = r*h*h           !Calculating the time step size
real(kind=wp), parameter                :: alpha = c*dt*dt/(h*h)!Constant used in the scheme
integer, parameter                      :: Nt = nint(t/dt)      !Setting time iterations to run according to simulation time

END MODULE